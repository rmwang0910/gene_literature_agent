"""
基因文献智能检索助手
Gene Literature Intelligent Retrieval Assistant

一个智能的基因文献分析工具，支持：
- 基因名称标准化
- PubMed 文献检索
- 摘要自动分析
- 结构化结论提取
- 冲突检测与报告生成
- 外部数据源整合（Open Targets、ClinVar、COSMIC、KEGG/Reactome）
"""

__version__ = "1.1.0"
__author__ = "Gene Literature Agent Team"

# 从配置模块导入
from .config import Config, default_config

# 从核心模块导入
from .core import (
    # 基因规范化
    GeneNormalizer,
    GeneInfo,
    normalize_gene,
    # PubMed 检索
    PubMedSearcher,
    PubMedArticle,
    SearchResult,
    # 摘要处理
    AbstractProcessor,
    AbstractFetcher,
    ProcessedAbstract,
    # 引用数据
    CitationEnricher,
    # 结论提取
    ConclusionExtractor,
    GeneConclusion,
    create_provider,
    # 冲突检测
    ConflictDetector,
    ConflictReport,
    Conflict,
    # 报告生成
    ReportGenerator,
    generate_report,
)

# 从工具模块导入
from .utils import (
    parse_gene_list,
    batch_process,
    estimate_cost,
    ProgressTracker,
)

__all__ = [
    # 版本信息
    "__version__",
    "__author__",
    # 配置
    "Config",
    "default_config",
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
    # 工具函数
    "parse_gene_list",
    "batch_process",
    "estimate_cost",
    "ProgressTracker",
]
