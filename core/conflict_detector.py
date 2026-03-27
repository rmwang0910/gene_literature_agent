"""
冲突检测模块 - 使用Sentence-BERT进行语义分析，检测文献结论之间的矛盾和冲突
"""
import os
import json
import logging
import re
from typing import List, Dict, Optional, Tuple, Set
from dataclasses import dataclass, field
from collections import defaultdict

# 设置离线模式，避免 HuggingFace 网络请求超时
os.environ["HF_HUB_OFFLINE"] = "1"
os.environ["TRANSFORMERS_OFFLINE"] = "1"

# 抑制 HuggingFace 警告
import warnings
warnings.filterwarnings("ignore", message=".*timed out.*")
warnings.filterwarnings("ignore", message=".*Retrying.*")

from ..config import default_config
from ..prompts import CONFLICT_DETECTION_PROMPT
from .conclusion_extractor import GeneConclusion, LLMProvider, OpenAIProvider, MockProvider
from ..constants import CONFLICT_DIMENSIONS, CORE_SEMANTIC_OPPOSITES, DISEASE_KEYWORDS

# 配置日志
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class Conflict:
    """冲突数据类"""
    topic: str  # 冲突主题
    conclusion_a: str  # 结论A
    pmid_a: str  # 结论A的PMID
    conclusion_b: str  # 结论B
    pmid_b: str  # 结论B的PMID
    resolution_suggestion: str = ""  # 解决建议
    severity: str = "medium"  # high, medium, low
    conflict_dimension: str = ""  # 冲突维度（如：疾病关系方向、预后）


@dataclass
class ConflictReport:
    """冲突报告"""
    gene: str
    conflicts: List[Conflict]
    consensus: List[str]  # 共识结论
    merged_conclusions: List[Dict]  # 合并后的结论
    needs_review: bool
    total_conclusions: int
    conflict_ratio: float  # 冲突比例
    similarity_matrix: Optional[List[List[float]]] = None  # 相似度矩阵


class SentenceEmbedder:
    """Sentence-BERT 向量化器"""

    # 模型优先级列表（按可用性排序，中文优先）
    MODEL_FALLBACKS = [
        "BAAI/bge-small-zh-v1.5",  # 中文模型
        "sentence-transformers/all-MiniLM-L6-v2",  # 英文模型
        "paraphrase-multilingual-MiniLM-L12-v2",  # 多语言模型
    ]

    def __init__(self, model_name: str = None):
        """
        初始化 Sentence-BERT 模型

        Args:
            model_name: 模型名称，默认尝试使用已缓存的模型
        """
        self.model = None
        self.model_name = model_name
        self._initialized = False

    def _lazy_init(self):
        """懒加载模型"""
        if self._initialized:
            return

        try:
            from sentence_transformers import SentenceTransformer

            # 如果指定了模型，尝试使用
            if self.model_name:
                try:
                    self.model = SentenceTransformer(self.model_name)
                    self._initialized = True
                    logger.info(f"Sentence-BERT 模型加载成功: {self.model_name}")
                    return
                except Exception as e:
                    logger.warning(f"加载指定模型失败 {self.model_name}: {e}")

            # 尝试回退模型列表
            for model_name in self.MODEL_FALLBACKS:
                try:
                    self.model = SentenceTransformer(model_name)
                    self.model_name = model_name
                    self._initialized = True
                    logger.info(f"Sentence-BERT 模型加载成功: {model_name}")
                    return
                except Exception as e:
                    logger.debug(f"尝试加载模型 {model_name} 失败: {e}")
                    continue

            # 所有模型都加载失败
            logger.warning("所有 Sentence-BERT 模型加载失败，将使用简单的文本匹配。")
            self.model = None
            self._initialized = True

        except ImportError:
            logger.warning("未安装 sentence-transformers，将使用简单的文本匹配。请运行: pip install sentence-transformers")
            self.model = None
            self._initialized = True

    def encode(self, texts: List[str]) -> List[List[float]]:
        """
        将文本转换为向量

        Args:
            texts: 文本列表

        Returns:
            向量列表
        """
        self._lazy_init()

        if self.model is None:
            # 降级方案：返回空向量
            return [[0.0] * 384 for _ in texts]

        embeddings = self.model.encode(texts, convert_to_numpy=True)
        return embeddings.tolist()

    def compute_similarity(self, vec_a: List[float], vec_b: List[float]) -> float:
        """
        计算两个向量的余弦相似度

        Args:
            vec_a: 向量A
            vec_b: 向量B

        Returns:
            相似度分数 (0-1)
        """
        import math

        dot_product = sum(a * b for a, b in zip(vec_a, vec_b))
        norm_a = math.sqrt(sum(a * a for a in vec_a))
        norm_b = math.sqrt(sum(b * b for b in vec_b))

        if norm_a == 0 or norm_b == 0:
            return 0.0

        return dot_product / (norm_a * norm_b)

    def compute_similarity_matrix(self, texts: List[str]) -> Tuple[List[List[float]], List[List[float]]]:
        """
        计算文本列表的相似度矩阵

        Args:
            texts: 文本列表

        Returns:
            (向量列表, 相似度矩阵)
        """
        if not texts:
            return [], []

        embeddings = self.encode(texts)
        n = len(texts)
        similarity_matrix = [[0.0] * n for _ in range(n)]

        for i in range(n):
            for j in range(n):
                if i == j:
                    similarity_matrix[i][j] = 1.0
                else:
                    similarity_matrix[i][j] = self.compute_similarity(embeddings[i], embeddings[j])

        return embeddings, similarity_matrix


class ConflictRules:
    """冲突规则定义"""

    # 使用常量中的冲突维度
    DIMENSION_OPPOSITES = CONFLICT_DIMENSIONS
    CORE_SEMANTIC_OPPOSITES = CORE_SEMANTIC_OPPOSITES

    # 非冲突的同义表述（仅疾病/组织不同，但结论方向一致）
    NON_CONFLICT_PATTERNS = [
        # 同一方向，不同癌症类型
        (r"促进.*?癌", r"促进.*?癌"),
        (r"抑制.*?癌", r"抑制.*?癌"),
        # 同一方向，不同组织
        (r"高表达于.*?组织", r"高表达于.*?组织"),
        (r"低表达于.*?组织", r"低表达于.*?组织"),
    ]

    @classmethod
    def extract_dimension_value(cls, text: str, dimension: str) -> Optional[str]:
        """
        提取文本中某个维度的值

        Args:
            text: 文本
            dimension: 维度名称

        Returns:
            维度值或None
        """
        text_lower = text.lower()
        opposites = cls.DIMENSION_OPPOSITES.get(dimension, [])

        for pair in opposites:
            if pair[0] in text_lower:
                return pair[0]
            if pair[1] in text_lower:
                return pair[1]

        return None

    @classmethod
    def check_dimension_conflict(cls, text_a: str, text_b: str, dimension: str) -> Tuple[bool, str]:
        """
        检查两个文本在特定维度上是否存在冲突

        Args:
            text_a: 文本A
            text_b: 文本B
            dimension: 维度名称

        Returns:
            (是否冲突, 维度值A, 维度值B)
        """
        text_a_lower = text_a.lower()
        text_b_lower = text_b.lower()

        opposites = cls.DIMENSION_OPPOSITES.get(dimension, [])

        for pair in opposites:
            a_has_first = pair[0] in text_a_lower
            a_has_second = pair[1] in text_a_lower
            b_has_first = pair[0] in text_b_lower
            b_has_second = pair[1] in text_b_lower

            # 检查是否对立
            if (a_has_first and b_has_second) or (a_has_second and b_has_first):
                return True, f"'{pair[0]}' vs '{pair[1]}'"

        return False, ""

    @classmethod
    def is_true_conflict(cls, text_a: str, text_b: str) -> Tuple[bool, str]:
        """
        判断两个结论是否真正冲突

        检测逻辑：
        1. 首先检查核心语义对立（强冲突，必须检测）
        2. 然后检查各维度的对立关系
        3. "促进肺癌" vs "促进乳腺癌" 不算冲突
        4. "促进肿瘤" vs "抑制肿瘤" 算冲突

        Args:
            text_a: 文本A
            text_b: 文本B

        Returns:
            (是否冲突, 冲突维度描述)
        """
        text_a_lower = text_a.lower()
        text_b_lower = text_b.lower()

        # 优先检查核心语义对立（强冲突）
        for phrase_a, phrase_b in cls.CORE_SEMANTIC_OPPOSITES:
            a_has_first = phrase_a in text_a_lower
            a_has_second = phrase_b in text_a_lower
            b_has_first = phrase_a in text_b_lower
            b_has_second = phrase_b in text_b_lower

            if (a_has_first and b_has_second) or (a_has_second and b_has_first):
                return True, f"核心语义冲突: '{phrase_a}' vs '{phrase_b}'"

        # 检查各维度
        for dimension in cls.DIMENSION_OPPOSITES.keys():
            is_conflict, detail = cls.check_dimension_conflict(text_a, text_b, dimension)
            if is_conflict:
                return True, f"{dimension}: {detail}"

        return False, ""

    @classmethod
    def is_same_direction_different_target(cls, text_a: str, text_b: str) -> bool:
        """
        检查是否是同方向、不同目标的情况

        例如：
        - "促进肺癌" vs "促进乳腺癌" -> True（不冲突）
        - "抑制肺癌" vs "抑制乳腺癌" -> True（不冲突）

        Args:
            text_a: 文本A
            text_b: 文本B

        Returns:
            是否为同方向不同目标
        """
        text_a_lower = text_a.lower()
        text_b_lower = text_b.lower()

        # 检查是否同时包含相同的方向词
        direction_words = ["促进", "抑制", "致癌", "抑癌", "高表达", "低表达"]

        for word in direction_words:
            if word in text_a_lower and word in text_b_lower:
                # 同方向，检查目标是否不同
                # 提取疾病/组织词（从配置加载）
                disease_patterns = DISEASE_KEYWORDS

                diseases_a = [d for d in disease_patterns if d in text_a_lower]
                diseases_b = [d for d in disease_patterns if d in text_b_lower]

                if diseases_a and diseases_b and diseases_a != diseases_b:
                    return True

        return False


class ConflictDetector:
    """冲突检测器 - 使用Sentence-BERT进行语义分析"""

    def __init__(
        self,
        provider: Optional[LLMProvider] = None,
        similarity_threshold: float = 0.8,
        conflict_similarity_threshold: float = 0.65,  # 新增：高于此相似度不视为冲突
        use_embedding: bool = True
    ):
        """
        初始化冲突检测器

        Args:
            provider: LLM 提供者（用于复杂冲突检测）
            similarity_threshold: 相似度阈值，高于此值的结论将被合并
            conflict_similarity_threshold: 冲突相似度阈值，高于此值不视为冲突
            use_embedding: 是否使用Sentence-BERT向量化
        """
        self.provider = provider or self._get_default_provider()
        self.similarity_threshold = similarity_threshold
        self.conflict_similarity_threshold = conflict_similarity_threshold
        self.use_embedding = use_embedding
        self.embedder = SentenceEmbedder() if use_embedding else None

    def _get_default_provider(self) -> LLMProvider:
        """获取默认 LLM 提供者"""
        if default_config.llm_api_key:
            return OpenAIProvider()
        return MockProvider()

    def _get_conclusion_text(self, c: GeneConclusion) -> str:
        """获取结论的完整文本表示"""
        parts = []
        if c.disease_relation:
            parts.append(f"疾病关系: {c.disease_relation}")
        if c.clinical_significance:
            parts.append(f"临床意义: {c.clinical_significance}")
        if c.mutation_effects and c.mutation_effects != "未提及":
            parts.append(f"突变效应: {c.mutation_effects}")
        return " | ".join(parts)

    def compute_similarities(self, conclusions: List[GeneConclusion]) -> Tuple[List[List[float]], List[str]]:
        """
        计算结论之间的语义相似度矩阵

        Args:
            conclusions: 结论列表

        Returns:
            (相似度矩阵, 文本列表)
        """
        if not conclusions or not self.use_embedding:
            return [], []

        texts = [self._get_conclusion_text(c) for c in conclusions]
        _, similarity_matrix = self.embedder.compute_similarity_matrix(texts)

        return similarity_matrix, texts

    def detect_conflicts_with_rules(
        self,
        conclusions: List[GeneConclusion],
        similarity_matrix: Optional[List[List[float]]] = None
    ) -> List[Conflict]:
        """
        使用规则检测真正的冲突

        仅当满足以下条件时才标记为冲突：
        1. 语义相似度 < conflict_similarity_threshold（避免相似结论被误判）
        2. 在关键维度上明显对立（如"促进"vs"抑制"）

        Args:
            conclusions: 结论列表
            similarity_matrix: 相似度矩阵（可选）

        Returns:
            检测到的冲突列表
        """
        conflicts = []
        n = len(conclusions)

        for i in range(n):
            for j in range(i + 1, n):
                conc_a = conclusions[i]
                conc_b = conclusions[j]

                # 检查语义相似度，如果太高则跳过冲突检测
                if similarity_matrix and i < len(similarity_matrix) and j < len(similarity_matrix[i]):
                    similarity = similarity_matrix[i][j]
                    if similarity >= self.conflict_similarity_threshold:
                        # 高相似度的结论不视为冲突，跳过
                        logger.debug(
                            f"跳过潜在冲突 (相似度 {similarity:.3f} >= {self.conflict_similarity_threshold}): "
                            f"PMID {conc_a.pmid} vs {conc_b.pmid}"
                        )
                        continue

                # 检查疾病关系冲突
                if conc_a.disease_relation and conc_b.disease_relation:
                    # 首先检查是否是"同方向不同目标"的情况
                    if ConflictRules.is_same_direction_different_target(
                        conc_a.disease_relation, conc_b.disease_relation
                    ):
                        # 不算冲突，跳过
                        continue

                    # 检查真正的维度冲突
                    is_conflict, detail = ConflictRules.is_true_conflict(
                        conc_a.disease_relation, conc_b.disease_relation
                    )
                    if is_conflict:
                        conflicts.append(Conflict(
                            topic="疾病关系",
                            conclusion_a=conc_a.disease_relation,
                            pmid_a=conc_a.pmid,
                            conclusion_b=conc_b.disease_relation,
                            pmid_b=conc_b.pmid,
                            severity="high",
                            conflict_dimension=detail
                        ))

                # 检查临床意义冲突
                if conc_a.clinical_significance and conc_b.clinical_significance:
                    is_conflict, detail = ConflictRules.is_true_conflict(
                        conc_a.clinical_significance, conc_b.clinical_significance
                    )
                    if is_conflict:
                        conflicts.append(Conflict(
                            topic="临床意义",
                            conclusion_a=conc_a.clinical_significance,
                            pmid_a=conc_a.pmid,
                            conclusion_b=conc_b.clinical_significance,
                            pmid_b=conc_b.pmid,
                            severity="high",
                            conflict_dimension=detail
                        ))

                # 检查突变效应冲突
                if conc_a.mutation_effects and conc_b.mutation_effects:
                    if conc_a.mutation_effects != "未提及" and conc_b.mutation_effects != "未提及":
                        is_conflict, detail = ConflictRules.is_true_conflict(
                            conc_a.mutation_effects, conc_b.mutation_effects
                        )
                        if is_conflict:
                            conflicts.append(Conflict(
                                topic="突变效应",
                                conclusion_a=conc_a.mutation_effects,
                                pmid_a=conc_a.pmid,
                                conclusion_b=conc_b.mutation_effects,
                                pmid_b=conc_b.pmid,
                                severity="medium",
                                conflict_dimension=detail
                            ))

        return conflicts

    def merge_similar_conclusions(
        self,
        conclusions: List[GeneConclusion],
        similarity_matrix: List[List[float]]
    ) -> List[Dict]:
        """
        合并相似的结论（相似度 > 阈值）

        注意：不会合并存在冲突的结论

        Args:
            conclusions: 结论列表
            similarity_matrix: 相似度矩阵

        Returns:
            合并后的结论列表
        """
        if not conclusions:
            return []

        n = len(conclusions)

        # 首先检测冲突，构建冲突图（传入相似度矩阵以过滤高相似度结论）
        conflict_pairs = set()
        conflicts = self.detect_conflicts_with_rules(conclusions, similarity_matrix)
        for c in conflicts:
            pair = tuple(sorted([c.pmid_a, c.pmid_b]))
            conflict_pairs.add(pair)

        # 构建 PMID 到索引的映射
        pmid_to_idx = {c.pmid: i for i, c in enumerate(conclusions)}

        merged = []
        merged_indices = set()

        for i in range(n):
            if i in merged_indices:
                continue

            # 找出与当前结论相似的所有结论（排除冲突的）
            similar_group = [i]
            for j in range(i + 1, n):
                if j in merged_indices:
                    continue

                # 检查是否存在冲突
                pair = tuple(sorted([conclusions[i].pmid, conclusions[j].pmid]))
                if pair in conflict_pairs:
                    # 存在冲突，不合并
                    continue

                if similarity_matrix and i < len(similarity_matrix) and j < len(similarity_matrix[i]):
                    if similarity_matrix[i][j] >= self.similarity_threshold:
                        similar_group.append(j)

            # 合并该组结论
            if len(similar_group) == 1:
                # 无相似结论
                merged.append({
                    "representative": conclusions[i].to_dict(),
                    "pmids": [conclusions[i].pmid],
                    "similarity_score": 1.0,
                    "is_merged": False
                })
            else:
                # 合并相似结论
                group_conclusions = [conclusions[idx] for idx in similar_group]
                pmids = [c.pmid for c in group_conclusions]

                # 选择最具代表性的结论（通常选择置信度最高的或最早的）
                representative = max(group_conclusions, key=lambda c: (
                    1 if c.confidence == "high" else (0.5 if c.confidence == "medium" else 0),
                    c.year if c.year else 0
                ))

                # 计算平均相似度
                avg_similarity = sum(
                    similarity_matrix[similar_group[0]][j] for j in similar_group[1:]
                ) / max(len(similar_group) - 1, 1) if len(similar_group) > 1 else 1.0

                merged.append({
                    "representative": representative.to_dict(),
                    "pmids": pmids,
                    "similarity_score": round(avg_similarity, 3),
                    "is_merged": True,
                    "merged_count": len(similar_group)
                })

                merged_indices.update(similar_group)

        logger.info(f"结论合并完成: 原始 {n} 条 -> 合并后 {len(merged)} 条 (排除 {len(conflict_pairs)} 对冲突)")
        return merged

    def detect_with_llm(
        self,
        gene: str,
        conclusions: List[GeneConclusion]
    ) -> Tuple[List[Conflict], List[str]]:
        """
        使用 LLM 进行更智能的冲突检测

        Args:
            gene: 基因名
            conclusions: 结论列表

        Returns:
            (冲突列表, 共识列表)
        """
        if not conclusions:
            return [], []

        # 构建提示词
        conclusions_text = "\n\n".join([
            f"PMID: {c.pmid}\n"
            f"疾病关系: {c.disease_relation}\n"
            f"信号通路: {', '.join(c.pathways)}\n"
            f"临床意义: {c.clinical_significance}\n"
            f"突变效应: {c.mutation_effects}\n"
            f"关键发现: {'; '.join(c.key_findings)}"
            for c in conclusions
        ])

        prompt = CONFLICT_DETECTION_PROMPT.format(
            gene=gene,
            conclusions=conclusions_text
        )

        try:
            response = self.provider.generate(prompt)

            # 解析响应
            json_match = re.search(r'\{[\s\S]*\}', response)
            if json_match:
                data = json.loads(json_match.group())
            else:
                data = json.loads(response)

            conflicts = []
            for cf in data.get("conflicts", []):
                conflicts.append(Conflict(
                    topic=cf.get("topic", "未知"),
                    conclusion_a=cf.get("conclusion_a", ""),
                    pmid_a=cf.get("pmid_a", ""),
                    conclusion_b=cf.get("conclusion_b", ""),
                    pmid_b=cf.get("pmid_b", ""),
                    resolution_suggestion=cf.get("resolution_suggestion", "")
                ))

            consensus = data.get("consensus", [])
            return conflicts, consensus

        except Exception as e:
            logger.error(f"LLM 冲突检测失败: {e}")
            return [], []

    def detect(
        self,
        gene: str,
        conclusions: List[GeneConclusion],
        use_llm: bool = True
    ) -> ConflictReport:
        """
        综合冲突检测

        Args:
            gene: 基因名
            conclusions: 结论列表
            use_llm: 是否使用 LLM 辅助检测

        Returns:
            ConflictReport 对象
        """
        if not conclusions:
            return ConflictReport(
                gene=gene,
                conflicts=[],
                consensus=[],
                merged_conclusions=[],
                needs_review=False,
                total_conclusions=0,
                conflict_ratio=0.0
            )

        # 1. 计算语义相似度矩阵
        similarity_matrix, texts = self.compute_similarities(conclusions)

        # 2. 基于规则检测冲突（仅检测真正的对立冲突，使用相似度过滤）
        rule_conflicts = self.detect_conflicts_with_rules(conclusions, similarity_matrix)
        logger.info(
            f"冲突检测完成: {len(conclusions)} 条结论中检测到 {len(rule_conflicts)} 个冲突 "
            f"(相似度阈值={self.conflict_similarity_threshold})"
        )

        # 3. 合并相似结论
        merged_conclusions = self.merge_similar_conclusions(conclusions, similarity_matrix)

        # 4. LLM 检测（可选）
        llm_conflicts = []
        consensus = []
        if use_llm and len(conclusions) > 1:
            try:
                llm_conflicts, consensus = self.detect_with_llm(gene, conclusions)
            except Exception as e:
                logger.warning(f"LLM 检测失败，仅使用规则检测: {e}")

        # 5. 合并冲突（去重）
        all_conflicts = self._merge_conflicts(rule_conflicts, llm_conflicts)

        # 6. 计算冲突比例
        conflict_ratio = len(all_conflicts) / max(len(conclusions), 1)

        # 7. 判断是否需要人工审核
        needs_review = len(all_conflicts) > 0 or conflict_ratio > 0.3

        return ConflictReport(
            gene=gene,
            conflicts=all_conflicts,
            consensus=consensus,
            merged_conclusions=merged_conclusions,
            needs_review=needs_review,
            total_conclusions=len(conclusions),
            conflict_ratio=conflict_ratio,
            similarity_matrix=similarity_matrix
        )

    def _merge_conflicts(
        self,
        rule_conflicts: List[Conflict],
        llm_conflicts: List[Conflict]
    ) -> List[Conflict]:
        """合并冲突列表，去除重复"""
        merged = list(rule_conflicts)
        existing_pairs = set()

        for c in rule_conflicts:
            pair = tuple(sorted([c.pmid_a, c.pmid_b]))
            existing_pairs.add(pair)

        for c in llm_conflicts:
            pair = tuple(sorted([c.pmid_a, c.pmid_b]))
            if pair not in existing_pairs:
                merged.append(c)
                existing_pairs.add(pair)

        return merged


class ConclusionMerger:
    """结论合并器"""

    def merge_conclusions(
        self,
        conclusions: List[GeneConclusion],
        conflict_report: ConflictReport
    ) -> Dict[str, any]:
        """
        合并结论

        Args:
            conclusions: 结论列表
            conflict_report: 冲突报告

        Returns:
            合并后的结论摘要
        """
        if not conclusions:
            return {
                "gene": conflict_report.gene,
                "total_papers": 0,
                "summary": "无文献数据",
                "disease_relations": [],
                "pathways": {},
                "clinical_significance": [],
                "mutations": [],
                "conflicts": [],
                "merged_conclusions": []
            }

        # 收集疾病关系
        disease_relations = defaultdict(list)
        for c in conclusions:
            if c.disease_relation and c.disease_relation not in ["解析失败", "提取失败"]:
                disease_relations[c.disease_relation].append(c.pmid)

        # 收集信号通路
        pathways = defaultdict(list)
        for c in conclusions:
            for pathway in c.pathways:
                pathways[pathway].append(c.pmid)

        # 收集临床意义
        clinical_sigs = defaultdict(list)
        for c in conclusions:
            if c.clinical_significance and c.clinical_significance not in ["解析失败", "提取失败"]:
                clinical_sigs[c.clinical_significance].append(c.pmid)

        # 收集突变信息
        mutations = []
        for c in conclusions:
            if c.mutation_effects and c.mutation_effects not in ["未提及", "解析失败", "提取失败"]:
                mutations.append({
                    "effect": c.mutation_effects,
                    "pmid": c.pmid
                })

        # 收集基因功能语境
        gene_contexts = defaultdict(list)
        for c in conclusions:
            if c.gene_function_context:
                gene_contexts[c.gene_function_context].append(c.pmid)

        return {
            "gene": conflict_report.gene,
            "total_papers": len(conclusions),
            "disease_relations": dict(disease_relations),
            "pathways": dict(pathways),
            "clinical_significance": dict(clinical_sigs),
            "mutations": mutations,
            "gene_contexts": dict(gene_contexts),
            "conflicts": [
                {
                    "topic": c.topic,
                    "conclusion_a": c.conclusion_a,
                    "pmid_a": c.pmid_a,
                    "conclusion_b": c.conclusion_b,
                    "pmid_b": c.pmid_b,
                    "severity": c.severity,
                    "conflict_dimension": c.conflict_dimension
                }
                for c in conflict_report.conflicts
            ],
            "consensus": conflict_report.consensus,
            "needs_review": conflict_report.needs_review,
            "merged_conclusions": conflict_report.merged_conclusions
        }


if __name__ == "__main__":
    # 测试
    from .conclusion_extractor import GeneConclusion

    # 创建测试结论
    conclusions = [
        GeneConclusion(
            gene="BRAF",
            pmid="12345678",
            disease_relation="突变型BRAF促进肺癌进展",
            pathways=["MAPK", "ERK"],
            clinical_significance="与不良预后相关",
            mutation_effects="V600E突变激活激酶活性",
            key_findings=["发现1"],
            confidence="high"
        ),
        GeneConclusion(
            gene="BRAF",
            pmid="23456789",
            disease_relation="突变型BRAF促进乳腺癌进展",
            pathways=["MAPK"],
            clinical_significance="与不良预后相关",
            mutation_effects="V600E突变激活激酶活性",
            key_findings=["发现2"],
            confidence="medium"
        ),
        GeneConclusion(
            gene="BRAF",
            pmid="34567890",
            disease_relation="突变型BRAF抑制肿瘤生长",
            pathways=["MAPK"],
            clinical_significance="与良好预后相关",
            mutation_effects="突变失活",
            key_findings=["发现3"],
            confidence="low"
        )
    ]

    detector = ConflictDetector(use_embedding=False)
    report = detector.detect("BRAF", conclusions, use_llm=False)

    print(f"发现 {len(report.conflicts)} 个冲突")
    for conflict in report.conflicts:
        print(f"  主题: {conflict.topic}")
        print(f"  冲突维度: {conflict.conflict_dimension}")
        print(f"  结论A (PMID {conflict.pmid_a}): {conflict.conclusion_a}")
        print(f"  结论B (PMID {conflict.pmid_b}): {conflict.conclusion_b}")

    print(f"\n合并后的结论: {len(report.merged_conclusions)} 条")
    for mc in report.merged_conclusions:
        if mc["is_merged"]:
            print(f"  合并 {mc['merged_count']} 条 -> PMIDs: {mc['pmids']}")
        else:
            print(f"  单条 -> PMID: {mc['pmids']}")