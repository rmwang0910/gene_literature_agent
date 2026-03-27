"""
工具函数模块
"""
import os
import json
import re
import time
from typing import List, Dict, Optional, Any, Callable
from functools import wraps
from pathlib import Path

from ..config import default_config


# ============================================================================
# 信号通路名称规范化 - 使用RAG检索
# ============================================================================

# RAG检索器（延迟加载）
_pathway_rag = None


def _get_pathway_rag():
    """获取通路RAG检索器（延迟加载）"""
    global _pathway_rag
    if _pathway_rag is None:
        try:
            from .pathway_rag import PathwayRAGRetriever
            _pathway_rag = PathwayRAGRetriever()
        except Exception:
            _pathway_rag = False  # 标记为不可用
    return _pathway_rag if _pathway_rag else None


def normalize_pathway(pathway: str) -> str:
    """
    规范化单个通路名称（使用RAG检索）

    Args:
        pathway: 原始通路名称

    Returns:
        规范化后的通路名称
    """
    if not pathway:
        return pathway

    pathway = pathway.strip()

    # 尝试使用RAG检索
    rag = _get_pathway_rag()
    if rag:
        result = rag.normalize(pathway, use_fuzzy=False)  # 先不使用模糊匹配，提高速度
        if result:
            return result

    # 没有匹配，返回原始值
    return pathway


def normalize_pathways(pathways: List[str], lang: str = "zh") -> List[str]:
    """
    规范化通路名称列表，去重，并转换为显示名称

    使用RAG向量检索智能合并中英文同义词，无需硬编码映射

    Args:
        pathways: 通路名称列表
        lang: 显示语言 ("zh" 中文, "en" 英文)

    Returns:
        规范化并去重后的通路列表
    """
    if not pathways:
        return []

    rag = _get_pathway_rag()
    normalized = []
    seen_standards = set()  # 记录已处理的标准名称（小写）

    for p in pathways:
        if not p or p in ["未提及", "解析失败", ""]:
            continue

        p = p.strip()

        # 使用RAG查找标准名称
        standard_name = p
        if rag:
            # 先精确匹配
            result = rag.exact_match(p)
            if result:
                standard_name = result
            else:
                # 模糊匹配（向量相似度）- 仅对非中文输入使用
                if not any('\u4e00' <= c <= '\u9fff' for c in p):
                    fuzzy_results = rag.fuzzy_match(p, top_k=1, threshold=0.85)
                    if fuzzy_results:
                        standard_name = fuzzy_results[0][0]

        # 基于标准名称去重
        standard_lower = standard_name.lower()
        if standard_lower in seen_standards:
            continue
        seen_standards.add(standard_lower)

        # 转换为显示名称
        display_name = standard_name
        if lang == "zh" and rag:
            # 查找中文同义词作为显示名称
            for doc in rag.documents:
                if doc.name.lower() == standard_lower:
                    chinese_synonyms = [
                        syn for syn in doc.synonyms
                        if any('\u4e00' <= c <= '\u9fff' for c in syn)
                    ]
                    if chinese_synonyms:
                        # 优先选择包含"信号通路"且首字符是大写或中文的
                        for syn in chinese_synonyms:
                            if "信号通路" in syn:
                                first_char = syn[0]
                                if first_char.isupper() or '\u4e00' <= first_char <= '\u9fff':
                                    display_name = syn
                                    break
                        else:
                            # 次优先选择包含"通路"的
                            for syn in chinese_synonyms:
                                if "通路" in syn:
                                    display_name = syn
                                    break
                            else:
                                # 最后选择第一个中文同义词
                                display_name = chinese_synonyms[0]
                    break

        normalized.append(display_name)

    return normalized


def get_pathway_display_name(pathway: str, lang: str = "zh") -> str:
    """
    获取通路的显示名称

    Args:
        pathway: 规范化的通路名称（英文标准名称）
        lang: 语言 ("zh" 或 "en")

    Returns:
        显示名称
    """
    if lang == "zh":
        # 尝试通过RAG查找中文名称
        rag = _get_pathway_rag()
        if rag:
            # 在文档中查找对应的中文同义词
            for doc in rag.documents:
                if doc.name.lower() == pathway.lower():
                    # 找包含中文的同义词
                    for syn in doc.synonyms:
                        if any('\u4e00' <= c <= '\u9fff' for c in syn):
                            return syn
        return pathway
    else:
        return pathway


# ============================================================================
# 突变位点规范化
# ============================================================================

# 氨基酸三字母到单字母映射
AMINO_ACID_MAP = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
    'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
    'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
    'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
    # 小写版本
    'ala': 'A', 'arg': 'R', 'asn': 'N', 'asp': 'D', 'cys': 'C',
    'gln': 'Q', 'glu': 'E', 'gly': 'G', 'his': 'H', 'ile': 'I',
    'leu': 'L', 'lys': 'K', 'met': 'M', 'phe': 'F', 'pro': 'P',
    'ser': 'S', 'thr': 'T', 'trp': 'W', 'tyr': 'Y', 'val': 'V',
}


def normalize_mutation_site(site: str) -> str:
    """
    规范化单个突变位点

    将各种表述形式统一为标准格式，如 "R175H"、"V600E"
    不预设热点突变，所有热点识别基于文献频次

    Args:
        site: 原始突变位点表述

    Returns:
        规范化后的突变位点
    """
    if not site:
        return site

    site = site.strip()

    # 去除常见的后缀
    site_clean = re.sub(r'(突变|mutation|mutant|variant)$', '', site, flags=re.IGNORECASE).strip()

    # 尝试解析三字母格式（如 Arg175His → R175H）
    three_letter_pattern = r'^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$'
    match = re.match(three_letter_pattern, site_clean, re.IGNORECASE)
    if match:
        ref_aa, position, alt_aa = match.groups()
        ref_single = AMINO_ACID_MAP.get(ref_aa, ref_aa[0].upper())
        alt_single = AMINO_ACID_MAP.get(alt_aa, alt_aa[0].upper())
        return f"{ref_single}{position}{alt_single}"

    # 尝试解析单字母格式（如 R175H）
    single_letter_pattern = r'^([A-Z])(\d+)([A-Z])$'
    if re.match(single_letter_pattern, site_clean, re.IGNORECASE):
        return site_clean.upper()

    # 尝试解析带基因名前缀的格式（如 BRAF V600E → V600E）
    gene_prefix_pattern = r'^[A-Z0-9]+\s+([A-Z]\d+[A-Z])$'
    match = re.match(gene_prefix_pattern, site_clean, re.IGNORECASE)
    if match:
        return match.group(1).upper()

    # Exon 删除格式
    exon_del_pattern = r'(exon\s*(\d+)\s*del|E(\d+)del)'
    match = re.search(exon_del_pattern, site_clean, re.IGNORECASE)
    if match:
        exon_num = match.group(2) or match.group(3)
        return f"E{exon_num}del"

    # 融合基因格式（如 EML4-ALK）
    fusion_pattern = r'^([A-Z0-9]+)-([A-Z0-9]+)$'
    if re.match(fusion_pattern, site_clean):
        return site_clean.upper()

    # 如果都不匹配，返回清理后的原始值（大写）
    return site_clean.upper()


def normalize_mutation_sites(sites: List[str]) -> List[str]:
    """
    规范化突变位点列表，去重

    Args:
        sites: 突变位点列表

    Returns:
        规范化并去重后的突变位点列表
    """
    if not sites:
        return []

    normalized = []
    seen = set()

    for s in sites:
        if not s or s in ["未提及", "解析失败", "", "无"]:
            continue

        norm_s = normalize_mutation_site(s)

        # 去重
        if norm_s.lower() not in seen:
            seen.add(norm_s.lower())
            normalized.append(norm_s)

    return normalized


def setup_environment(
    ncbi_email: Optional[str] = None,
    ncbi_api_key: Optional[str] = None,
    llm_api_key: Optional[str] = None,
    llm_provider: str = "openai"
) -> Dict[str, str]:
    """
    设置环境变量

    Args:
        ncbi_email: NCBI 邮箱
        ncbi_api_key: NCBI API Key
        llm_api_key: LLM API Key
        llm_provider: LLM 提供者

    Returns:
        设置的环境变量字典
    """
    env_vars = {}

    if ncbi_email:
        os.environ["NCBI_EMAIL"] = ncbi_email
        env_vars["NCBI_EMAIL"] = ncbi_email

    if ncbi_api_key:
        os.environ["NCBI_API_KEY"] = ncbi_api_key
        env_vars["NCBI_API_KEY"] = "***"  # 隐藏实际值

    if llm_api_key:
        os.environ["LLM_API_KEY"] = llm_api_key
        env_vars["LLM_API_KEY"] = "***"

    os.environ["LLM_PROVIDER"] = llm_provider
    env_vars["LLM_PROVIDER"] = llm_provider

    return env_vars


def parse_gene_list(input_text: str) -> List[str]:
    """
    解析基因列表输入

    Args:
        input_text: 输入文本（可以是逗号分隔、换行分隔或空格分隔）

    Returns:
        基因名列表
    """
    # 替换各种分隔符为统一格式
    text = input_text.replace(',', '\n').replace(';', '\n').replace('\t', '\n')

    # 分割并清理
    genes = []
    for line in text.split('\n'):
        gene = line.strip()
        if gene:
            genes.append(gene)

    return genes


def format_pmid_link(pmid: str) -> str:
    """
    格式化 PMID 链接

    Args:
        pmid: PubMed ID

    Returns:
        格式化的链接
    """
    return f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"


def batch_process(
    items: List[Any],
    processor: Callable,
    batch_size: int = 10,
    delay: float = 0.5,
    progress_callback: Optional[Callable] = None
) -> List[Any]:
    """
    批量处理工具

    Args:
        items: 待处理项目列表
        processor: 处理函数
        batch_size: 批次大小
        delay: 批次间延迟
        progress_callback: 进度回调函数

    Returns:
        处理结果列表
    """
    results = []
    total = len(items)

    for i, item in enumerate(items):
        try:
            result = processor(item)
            results.append(result)
        except Exception as e:
            print(f"处理失败: {e}")
            results.append(None)

        if progress_callback:
            progress_callback(i + 1, total)

        if (i + 1) % batch_size == 0 and i + 1 < total:
            time.sleep(delay)

    return results


def retry_on_failure(
    func: Callable,
    max_retries: int = 3,
    delay: float = 1.0,
    backoff: float = 2.0
) -> Callable:
    """
    失败重试装饰器

    Args:
        func: 要包装的函数
        max_retries: 最大重试次数
        delay: 初始延迟
        backoff: 退避系数

    Returns:
        包装后的函数
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        last_error = None
        current_delay = delay

        for attempt in range(max_retries + 1):
            try:
                return func(*args, **kwargs)
            except Exception as e:
                last_error = e
                if attempt < max_retries:
                    time.sleep(current_delay)
                    current_delay *= backoff

        raise last_error

    return wrapper


def save_cache(cache_file: str, data: Any) -> None:
    """
    保存缓存

    Args:
        cache_file: 缓存文件路径
        data: 要缓存的数据
    """
    cache_path = Path(cache_file)
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    with open(cache_path, 'w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=2)


def load_cache(cache_file: str) -> Optional[Any]:
    """
    加载缓存

    Args:
        cache_file: 缓存文件路径

    Returns:
        缓存的数据，如果不存在返回 None
    """
    cache_path = Path(cache_file)

    if not cache_path.exists():
        return None

    try:
        with open(cache_path, 'r', encoding='utf-8') as f:
            return json.load(f)
    except (json.JSONDecodeError, IOError):
        return None


def format_duration(seconds: float) -> str:
    """
    格式化持续时间

    Args:
        seconds: 秒数

    Returns:
        格式化的时间字符串
    """
    if seconds < 60:
        return f"{seconds:.1f}秒"
    elif seconds < 3600:
        minutes = seconds / 60
        return f"{minutes:.1f}分钟"
    else:
        hours = seconds / 3600
        return f"{hours:.1f}小时"


def estimate_cost(
    num_abstracts: int,
    avg_tokens_per_abstract: int = 500,
    cost_per_1k_tokens: float = 0.00015
) -> Dict[str, float]:
    """
    估算 LLM 调用成本

    Args:
        num_abstracts: 摘要数量
        avg_tokens_per_abstract: 每个摘要的平均 token 数
        cost_per_1k_tokens: 每 1000 token 的成本

    Returns:
        成本估算字典
    """
    total_tokens = num_abstracts * avg_tokens_per_abstract * 2  # 输入+输出
    total_cost = (total_tokens / 1000) * cost_per_1k_tokens

    return {
        "num_abstracts": num_abstracts,
        "total_tokens": total_tokens,
        "estimated_cost_usd": total_cost,
        "estimated_cost_cny": total_cost * 7.2
    }


def validate_gene_name(gene: str) -> bool:
    """
    验证基因名称格式

    Args:
        gene: 基因名

    Returns:
        是否有效
    """
    if not gene:
        return False

    # 基因名通常只包含字母、数字和连字符
    import re
    pattern = r'^[A-Za-z][A-Za-z0-9\-]*$'
    return bool(re.match(pattern, gene))


def create_summary_statistics(
    conclusions: List[Dict],
    conflicts: List[Dict]
) -> Dict[str, Any]:
    """
    创建汇总统计

    Args:
        conclusions: 结论列表
        conflicts: 冲突列表

    Returns:
        统计字典
    """
    if not conclusions:
        return {
            "total_papers": 0,
            "total_conflicts": 0,
            "conflict_rate": 0,
            "high_confidence_count": 0,
            "pathways": {},
            "top_pathways": []
        }

    # 统计置信度
    confidence_counts = {"high": 0, "medium": 0, "low": 0}
    pathway_counts = {}

    for c in conclusions:
        conf = c.get("confidence", "medium")
        confidence_counts[conf] = confidence_counts.get(conf, 0) + 1

        for pathway in c.get("pathways", []):
            pathway_counts[pathway] = pathway_counts.get(pathway, 0) + 1

    # 排序通路
    top_pathways = sorted(
        pathway_counts.items(),
        key=lambda x: x[1],
        reverse=True
    )[:5]

    return {
        "total_papers": len(conclusions),
        "total_conflicts": len(conflicts),
        "conflict_rate": len(conflicts) / len(conclusions) if conclusions else 0,
        "confidence_distribution": confidence_counts,
        "pathways": pathway_counts,
        "top_pathways": top_pathways
    }


class ProgressTracker:
    """进度跟踪器"""

    def __init__(self, total: int, description: str = "处理中"):
        """
        初始化进度跟踪器

        Args:
            total: 总任务数
            description: 任务描述
        """
        self.total = total
        self.description = description
        self.current = 0
        self.start_time = time.time()

    def update(self, increment: int = 1) -> None:
        """
        更新进度

        Args:
            increment: 增量
        """
        self.current += increment

    def get_progress(self) -> Dict[str, Any]:
        """
        获取进度信息

        Returns:
            进度信息字典
        """
        elapsed = time.time() - self.start_time
        progress = self.current / self.total if self.total > 0 else 0

        if progress > 0:
            estimated_total = elapsed / progress
            estimated_remaining = estimated_total - elapsed
        else:
            estimated_remaining = 0

        return {
            "description": self.description,
            "current": self.current,
            "total": self.total,
            "progress": progress,
            "elapsed": elapsed,
            "estimated_remaining": estimated_remaining
        }

    def format_progress(self) -> str:
        """
        格式化进度显示

        Returns:
            格式化的进度字符串
        """
        info = self.get_progress()
        return (
            f"{info['description']}: {info['current']}/{info['total']} "
            f"({info['progress']:.0%}) - "
            f"已用时: {format_duration(info['elapsed'])}, "
            f"预计剩余: {format_duration(info['estimated_remaining'])}"
        )


if __name__ == "__main__":
    # 测试
    genes = parse_gene_list("BRAF, TP53\nEGFR MYC")
    print(f"解析基因: {genes}")

    cost = estimate_cost(100)
    print(f"成本估算: {cost}")

    tracker = ProgressTracker(10, "处理摘要")
    for i in range(10):
        tracker.update()
        print(tracker.format_progress())
        time.sleep(0.1)