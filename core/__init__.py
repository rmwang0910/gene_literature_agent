"""
核心业务模块
Core Business Modules
"""

from .gene_normalizer import GeneNormalizer, GeneInfo, normalize_gene
from .pubmed_searcher import PubMedSearcher, PubMedArticle, SearchResult
from .abstract_fetcher import AbstractProcessor, AbstractFetcher, ProcessedAbstract
from .citation_fetcher import (
    CitationEnricher,
    CitationInfo,
    get_citation_stars,
    get_impact_badge,
    calculate_confidence_score,
)
from .conclusion_extractor import ConclusionExtractor, GeneConclusion, create_provider
from .conflict_detector import ConflictDetector, ConflictReport, Conflict
from .report_generator import ReportGenerator, generate_report

__all__ = [
    # 基因规范化
    "GeneNormalizer",
    "GeneInfo",
    "normalize_gene",
    # PubMed 检索
    "PubMedSearcher",
    "PubMedArticle",
    "SearchResult",
    # 摘要处理
    "AbstractProcessor",
    "AbstractFetcher",
    "ProcessedAbstract",
    # 引用数据
    "CitationEnricher",
    "CitationInfo",
    "get_citation_stars",
    "get_impact_badge",
    "calculate_confidence_score",
    # 结论提取
    "ConclusionExtractor",
    "GeneConclusion",
    "create_provider",
    # 冲突检测
    "ConflictDetector",
    "ConflictReport",
    "Conflict",
    # 报告生成
    "ReportGenerator",
    "generate_report",
]
