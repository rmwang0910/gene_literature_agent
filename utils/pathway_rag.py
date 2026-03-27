"""
通路RAG系统 - 基于向量检索的通路名称标准化
无需硬编码，通过向量相似度匹配实现模糊查询
"""
import os
import json
import logging
import pickle
import time
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass, field, asdict
from pathlib import Path

# 配置日志
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class PathwayDocument:
    """通路文档"""
    id: str  # 数据库ID (KEGG/Reactome)
    name: str  # 标准名称
    source: str  # 数据来源
    synonyms: List[str] = field(default_factory=list)  # 同义词
    description: str = ""  # 描述（可选）

    def to_dict(self) -> Dict:
        return asdict(self)


class PathwayRAGBuilder:
    """通路RAG索引构建器"""

    def __init__(self, data_dir: str = None):
        self.data_dir = data_dir or os.path.join(os.path.dirname(__file__), "..", "data", "pathway")
        self.documents: List[PathwayDocument] = []
        self.name_to_doc: Dict[str, PathwayDocument] = {}  # 名称到文档的映射
        self.embeddings = None
        self.embedder = None

    def load_from_downloader(self, mapping: Dict[str, str], source_info: Dict[str, str] = None):
        """
        从下载器的映射表加载

        Args:
            mapping: {同义词: 标准名称} 映射
            source_info: {标准名称: 数据库ID} 可选
        """
        logger.info("正在构建通路文档...")

        # 按标准名称分组，收集所有同义词
        grouped = {}
        for synonym, standard in mapping.items():
            if standard not in grouped:
                grouped[standard] = set()
            # 添加所有同义词（包括中英文）
            if synonym.lower() != standard.lower():
                grouped[standard].add(synonym)

        # 创建文档，包含所有同义词
        for standard, synonyms in grouped.items():
            source_id = source_info.get(standard, "") if source_info else ""
            source = "KEGG" if source_id.startswith("path:") else "Reactome" if source_id.startswith("R-") else "Unknown"

            doc = PathwayDocument(
                id=source_id,
                name=standard,
                source=source,
                synonyms=list(synonyms)
            )
            self.documents.append(doc)
            self.name_to_doc[standard.lower()] = doc

        logger.info(f"构建了 {len(self.documents)} 个通路文档")

    def load_from_kegg_reactome(self):
        """直接从KEGG和Reactome下载并构建"""
        import requests

        logger.info("正在从KEGG和Reactome下载通路数据...")

        # 下载KEGG通路
        try:
            resp = requests.get("https://rest.kegg.jp/list/pathway/hsa", timeout=30)
            resp.raise_for_status()

            for line in resp.text.strip().split('\n'):
                parts = line.split('\t')
                if len(parts) >= 2:
                    import re
                    pathway_id = parts[0]
                    name = re.sub(r'\s*-\s*Homo sapiens.*$', '', parts[1]).strip()

                    if name and name.lower() not in self.name_to_doc:
                        doc = PathwayDocument(
                            id=pathway_id,
                            name=name,
                            source="KEGG"
                        )
                        self.documents.append(doc)
                        self.name_to_doc[name.lower()] = doc

            logger.info(f"KEGG: {len([d for d in self.documents if d.source == 'KEGG'])} 条通路")
        except Exception as e:
            logger.warning(f"KEGG下载失败: {e}")

        # 下载Reactome通路
        try:
            resp = requests.get("https://reactome.org/download/current/ReactomePathways.txt", timeout=60)
            resp.raise_for_status()

            for line in resp.text.strip().split('\n'):
                parts = line.split('\t')
                if len(parts) >= 3 and parts[2].lower() == "homo sapiens":
                    pathway_id = parts[0]
                    name = parts[1].strip()

                    if name and name.lower() not in self.name_to_doc:
                        doc = PathwayDocument(
                            id=pathway_id,
                            name=name,
                            source="Reactome"
                        )
                        self.documents.append(doc)
                        self.name_to_doc[name.lower()] = doc

            logger.info(f"Reactome: {len([d for d in self.documents if d.source == 'Reactome'])} 条通路")
        except Exception as e:
            logger.warning(f"Reactome下载失败: {e}")

        logger.info(f"总计 {len(self.documents)} 条通路")

    def build_embeddings(self, model_name: str = "BAAI/bge-small-zh-v1.5"):
        """
        构建向量索引

        Args:
            model_name: Sentence-BERT模型名称
        """
        logger.info(f"正在构建向量索引 (模型: {model_name})...")

        try:
            from sentence_transformers import SentenceTransformer
            import numpy as np

            # 设置离线模式
            os.environ["HF_HUB_OFFLINE"] = "1"
            os.environ["TRANSFORMERS_OFFLINE"] = "1"

            self.embedder = SentenceTransformer(model_name)

            # 构建待嵌入文本列表
            texts = []
            for doc in self.documents:
                # 组合名称和同义词
                text = doc.name
                if doc.synonyms:
                    text += " " + " ".join(doc.synonyms[:5])  # 最多取5个同义词
                texts.append(text)

            # 批量嵌入
            self.embeddings = self.embedder.encode(texts, show_progress_bar=True)
            logger.info(f"向量维度: {self.embeddings.shape}")

        except ImportError:
            logger.warning("未安装 sentence-transformers，跳过向量索引构建")
        except Exception as e:
            logger.warning(f"向量索引构建失败: {e}")

    def save(self, filename: str = "pathway_rag_index.pkl"):
        """
        保存RAG索引到文件

        Args:
            filename: 输出文件名
        """
        output_path = os.path.join(self.data_dir, filename)

        data = {
            "documents": [doc.to_dict() for doc in self.documents],
            "name_to_doc": {k: v.to_dict() for k, v in self.name_to_doc.items()},
            "embeddings": self.embeddings.tolist() if self.embeddings is not None else None,
            "model_name": getattr(self.embedder, 'model_name', None) if self.embedder else None,
            "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
            "version": "1.0"
        }

        with open(output_path, 'wb') as f:
            pickle.dump(data, f)

        logger.info(f"RAG索引已保存: {output_path}")

        # 同时保存JSON格式（可读）
        json_path = output_path.replace('.pkl', '.json')
        json_data = {
            "created_at": data["created_at"],
            "version": data["version"],
            "total_documents": len(self.documents),
            "pathways": [doc.to_dict() for doc in self.documents]
        }
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(json_data, f, ensure_ascii=False, indent=2)

        logger.info(f"JSON文件已保存: {json_path}")

        return output_path


class PathwayRAGRetriever:
    """通路RAG检索器"""

    def __init__(self, index_path: str = None, data_dir: str = None):
        self.data_dir = data_dir or os.path.join(os.path.dirname(__file__), "..", "data", "pathway")
        self.index_path = index_path or os.path.join(self.data_dir, "pathway_rag_index.pkl")

        self.documents: List[PathwayDocument] = []
        self.name_to_doc: Dict[str, PathwayDocument] = {}
        self.embeddings = None
        self.embedder = None

        self._load()

    def _load(self):
        """加载RAG索引"""
        if not os.path.exists(self.index_path):
            logger.warning(f"RAG索引文件不存在: {self.index_path}")
            return

        with open(self.index_path, 'rb') as f:
            data = pickle.load(f)

        self.documents = [PathwayDocument(**d) for d in data.get("documents", [])]
        self.name_to_doc = {k: PathwayDocument(**v) for k, v in data.get("name_to_doc", {}).items()}

        if data.get("embeddings") is not None:
            import numpy as np
            self.embeddings = np.array(data["embeddings"])

        logger.info(f"加载了 {len(self.documents)} 个通路文档")

    def _init_embedder(self):
        """延迟初始化嵌入模型"""
        if self.embedder is not None:
            return

        try:
            from sentence_transformers import SentenceTransformer

            os.environ["HF_HUB_OFFLINE"] = "1"
            os.environ["TRANSFORMERS_OFFLINE"] = "1"

            model_name = "BAAI/bge-small-zh-v1.5"
            self.embedder = SentenceTransformer(model_name)
            logger.info(f"嵌入模型已加载: {model_name}")
        except Exception as e:
            logger.warning(f"嵌入模型加载失败: {e}")

    def exact_match(self, query: str) -> Optional[str]:
        """
        精确匹配

        Args:
            query: 查询字符串

        Returns:
            标准通路名称或None
        """
        query_lower = query.lower().strip()

        # 直接匹配
        if query_lower in self.name_to_doc:
            return self.name_to_doc[query_lower].name

        # 在同义词中匹配
        for doc in self.documents:
            if query_lower in [s.lower() for s in doc.synonyms]:
                return doc.name

        return None

    def fuzzy_match(self, query: str, top_k: int = 3, threshold: float = 0.7) -> List[Tuple[str, float]]:
        """
        模糊匹配（基于向量相似度）

        Args:
            query: 查询字符串
            top_k: 返回前K个结果
            threshold: 相似度阈值

        Returns:
            [(标准名称, 相似度), ...]
        """
        if self.embeddings is None or len(self.embeddings) == 0:
            # 降级到精确匹配
            result = self.exact_match(query)
            if result:
                return [(result, 1.0)]
            return []

        self._init_embedder()

        if self.embedder is None:
            return []

        try:
            import numpy as np

            # 嵌入查询
            query_embedding = self.embedder.encode([query])[0]

            # 计算余弦相似度
            similarities = np.dot(self.embeddings, query_embedding) / (
                np.linalg.norm(self.embeddings, axis=1) * np.linalg.norm(query_embedding)
            )

            # 获取top_k
            top_indices = np.argsort(similarities)[::-1][:top_k]

            results = []
            for idx in top_indices:
                if similarities[idx] >= threshold:
                    results.append((self.documents[idx].name, float(similarities[idx])))

            return results

        except Exception as e:
            logger.warning(f"向量检索失败: {e}")
            return []

    def normalize(self, query: str, use_fuzzy: bool = True) -> str:
        """
        规范化通路名称

        Args:
            query: 原始通路名称
            use_fuzzy: 是否使用模糊匹配

        Returns:
            规范化后的名称（未匹配则返回原值）
        """
        if not query or not query.strip():
            return query

        # 1. 先尝试精确匹配
        result = self.exact_match(query)
        if result:
            return result

        # 2. 模糊匹配
        if use_fuzzy:
            fuzzy_results = self.fuzzy_match(query, top_k=1, threshold=0.75)
            if fuzzy_results:
                return fuzzy_results[0][0]

        # 3. 未匹配，返回原值（做简单清理）
        return query.strip()

    def normalize_batch(self, queries: List[str], use_fuzzy: bool = True) -> List[str]:
        """
        批量规范化

        Args:
            queries: 通路名称列表
            use_fuzzy: 是否使用模糊匹配

        Returns:
            规范化后的列表
        """
        return [self.normalize(q, use_fuzzy) for q in queries]


def build_pathway_rag(output_dir: str = None):
    """
    构建通路RAG索引

    Args:
        output_dir: 输出目录
    """
    print("=" * 60)
    print("通路RAG索引构建器")
    print("=" * 60)

    builder = PathwayRAGBuilder(data_dir=output_dir)

    # 方法1: 从已有的mapping.json加载
    mapping_path = os.path.join(builder.data_dir, "pathway_mapping.json")
    if os.path.exists(mapping_path):
        print(f"从 {mapping_path} 加载映射...")
        with open(mapping_path, 'r', encoding='utf-8') as f:
            mapping = json.load(f)
        builder.load_from_downloader(mapping)
    else:
        # 方法2: 直接从数据库下载
        print("从KEGG和Reactome下载...")
        builder.load_from_kegg_reactome()

    # 构建向量索引
    builder.build_embeddings()

    # 保存
    index_path = builder.save()

    print()
    print("=" * 60)
    print(f"构建完成!")
    print(f"索引文件: {index_path}")
    print(f"文档数量: {len(builder.documents)}")
    print("=" * 60)

    return index_path


# 模块级单例检索器
_retriever: Optional[PathwayRAGRetriever] = None


def get_pathway_retriever() -> PathwayRAGRetriever:
    """获取全局检索器实例"""
    global _retriever
    if _retriever is None:
        _retriever = PathwayRAGRetriever()
    return _retriever


def normalize_pathway_rag(query: str) -> str:
    """
    使用RAG规范化通路名称（模块级便捷函数）

    Args:
        query: 原始通路名称

    Returns:
        规范化后的名称
    """
    retriever = get_pathway_retriever()
    if retriever.documents:
        return retriever.normalize(query)
    # RAG索引不存在时，返回原值
    return query


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="通路RAG系统")
    parser.add_argument("--build", action="store_true", help="构建RAG索引")
    parser.add_argument("--query", type=str, help="测试查询")
    parser.add_argument("--output-dir", default=None, help="输出目录")
    args = parser.parse_args()

    if args.build:
        build_pathway_rag(args.output_dir)

    if args.query:
        retriever = PathwayRAGRetriever(data_dir=args.output_dir)
        print(f"\n查询: {args.query}")
        print(f"精确匹配: {retriever.exact_match(args.query)}")
        print(f"模糊匹配: {retriever.fuzzy_match(args.query)}")
        print(f"规范化: {retriever.normalize(args.query)}")