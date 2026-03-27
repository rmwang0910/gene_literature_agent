"""
配置模块 - 管理API密钥和系统配置
支持从 config.yaml 文件或环境变量读取配置
"""
import os
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any
from pathlib import Path


def load_yaml_config(config_path: Optional[str] = None) -> dict:
    """
    从 YAML 文件加载配置

    Args:
        config_path: 配置文件路径，默认为当前目录的 config.yaml

    Returns:
        配置字典
    """
    if config_path is None:
        config_path = Path(__file__).parent / "config.yaml"
    else:
        config_path = Path(config_path)

    if not config_path.exists():
        return {}

    try:
        import yaml
        with open(config_path, 'r', encoding='utf-8') as f:
            return yaml.safe_load(f) or {}
    except ImportError:
        print("警告: 未安装 PyYAML，无法读取 config.yaml。请运行: pip install pyyaml")
        return {}
    except Exception as e:
        print(f"警告: 读取配置文件失败: {e}")
        return {}


# 加载 YAML 配置
_yaml_config = load_yaml_config()


def get_api_url(api_name: str) -> str:
    """
    获取 API URL

    Args:
        api_name: API 名称 (crossref, openalex, ncbi, etc.)

    Returns:
        API URL
    """
    default_urls = {
        "crossref": "https://api.crossref.org",
        "openalex": "https://api.openalex.org",
        "ncbi": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils",
        "mygene": "https://mygene.info/v3",
        "kegg": "https://rest.kegg.jp",
        "reactome": "https://reactome.org",
        "wikipathways": "https://sparql.wikipathways.org/sparql",
        "cbioportal": "https://www.cbioportal.org/api",
        "cosmic": "https://cancer.sanger.ac.uk/api",
        "oncokb": "https://www.oncokb.org/api/v1",
    }
    return _yaml_config.get("apis", {}).get(api_name, default_urls.get(api_name, ""))


def get_wildtype_keywords() -> List[str]:
    """获取野生型关键词列表"""
    return _yaml_config.get("wildtype_keywords", [
        "野生型", "wild-type", "wildtype", "wild type",
        "正常", "normal", "native", "生理",
        "功能正常", "functional", "非突变"
    ])


def get_mutant_keywords() -> List[str]:
    """获取突变型关键词列表"""
    return _yaml_config.get("mutant_keywords", [
        "突变", "mutation", "mutant", "mutated",
        "缺失", "deletion", "deleted",
        "失活", "inactivat", "loss-of-function", "功能丧失",
        "激活", "activat", "gain-of-function", "功能获得",
        "变异", "variant", "variation",
        "点突变", "point mutation", "移码", "frameshift",
        "插入", "insertion", "插入缺失", "indel",
        "融合", "fusion",
        "过表达", "overexpress", "高表达",
        "低表达", "underexpress", "down-regul",
        "扩增", "amplif", "amplification",
        "重排", "rearrangement", "translocation"
    ])


def get_disease_keywords() -> List[str]:
    """获取疾病关键词列表（按优先级排序）"""
    return _yaml_config.get("disease_keywords", [
        "非小细胞肺癌", "小细胞肺癌", "肺腺癌", "肺鳞癌", "肺癌",
        "三阴性乳腺癌", "HER2阳性乳腺癌", "乳腺癌",
        "肝细胞癌", "胆管细胞癌", "肝癌",
        "胃癌", "结直肠癌", "结肠癌", "直肠癌", "胰腺癌", "食管癌",
        "卵巢癌", "宫颈癌", "子宫内膜癌",
        "前列腺癌", "膀胱癌", "肾癌",
        "白血病", "淋巴瘤", "骨髓瘤",
        "胶质瘤", "胶质母细胞瘤", "脑膜瘤", "神经母细胞瘤",
        "黑色素瘤", "骨肉瘤", "甲状腺癌", "头颈癌",
        "恶性肿瘤", "肿瘤", "癌症", "癌"
    ])


def get_conflict_dimensions() -> Dict[str, List[List[str]]]:
    """获取冲突维度定义"""
    default = {
        "疾病关系方向": [
            ["促进", "抑制"], ["致癌", "抑癌"], ["驱动", "抑制"],
            ["激活", "失活"], ["上调", "下调"], ["高表达", "低表达"],
        ],
        "预后": [
            ["不良", "良好"], ["差", "好"], ["预后差", "预后好"],
        ],
        "治疗响应": [
            ["耐药", "敏感"], ["抵抗", "敏感"], ["无效", "有效"],
        ],
    }
    return _yaml_config.get("conflict_dimensions", default)


def get_core_semantic_opposites() -> List[List[str]]:
    """获取核心语义对立词组"""
    return _yaml_config.get("core_semantic_opposites", [
        ["预后好", "预后差"],
        ["预后良好", "预后不良"],
        ["生存率高", "生存率低"],
        ["促进肿瘤", "抑制肿瘤"],
        ["致癌", "抑癌"],
        ["功能获得", "功能丧失"],
        ["对治疗敏感", "对治疗耐药"],
        ["治疗有效", "治疗无效"],
    ])


def get_impact_factor_tiers() -> Dict[str, Dict[str, Any]]:
    """获取影响因子阈值配置"""
    return _yaml_config.get("impact_factor_tiers", {
        "top": {"threshold": 50, "label": "顶级期刊", "badge": "🏆"},
        "high": {"threshold": 20, "label": "高影响力", "badge": "⭐"},
        "medium_high": {"threshold": 10, "label": "中高影响力", "badge": "📈"},
        "medium": {"threshold": 5, "label": "中等影响力", "badge": "📄"},
        "low": {"threshold": 0, "label": "", "badge": ""},
    })


def get_confidence_levels() -> Dict[str, Dict[str, Any]]:
    """获取置信度等级配置"""
    return _yaml_config.get("confidence_levels", {
        "high": {"threshold": 80, "label": "高置信度", "stars": "★★★★★"},
        "medium": {"threshold": 50, "label": "中等置信度", "stars": "★★★☆☆"},
        "low": {"threshold": 0, "label": "低置信度", "stars": "★☆☆☆☆"},
    })


@dataclass
class Config:
    """系统配置类"""

    # NCBI 配置
    ncbi_email: str = field(default_factory=lambda: (
        _yaml_config.get("ncbi", {}).get("email", "") or
        os.getenv("NCBI_EMAIL", "")
    ))
    ncbi_api_key: str = field(default_factory=lambda: (
        _yaml_config.get("ncbi", {}).get("api_key", "") or
        os.getenv("NCBI_API_KEY", "")
    ))

    # LLM 配置
    llm_provider: str = field(default_factory=lambda: (
        _yaml_config.get("llm", {}).get("provider", "openai") or
        os.getenv("LLM_PROVIDER", "openai")
    ))
    llm_api_key: str = field(default_factory=lambda: (
        _yaml_config.get("llm", {}).get("api_key", "") or
        os.getenv("LLM_API_KEY", "")
    ))
    llm_model: str = field(default_factory=lambda: (
        _yaml_config.get("llm", {}).get("model", "qwen-plus") or
        os.getenv("LLM_MODEL", "qwen-plus")
    ))
    llm_base_url: Optional[str] = field(default_factory=lambda: (
        _yaml_config.get("llm", {}).get("base_url", "") or
        os.getenv("LLM_BASE_URL", "") or None
    ))
    llm_timeout: int = field(default_factory=lambda: (
        _yaml_config.get("llm", {}).get("timeout", 120) or
        int(os.getenv("LLM_TIMEOUT", "120"))
    ))
    llm_max_tokens: int = field(default_factory=lambda: (
        _yaml_config.get("llm", {}).get("max_tokens", 4096) or
        int(os.getenv("LLM_MAX_TOKENS", "4096"))
    ))
    llm_temperature: float = field(default_factory=lambda: (
        _yaml_config.get("llm", {}).get("temperature", 0.3) or
        float(os.getenv("LLM_TEMPERATURE", "0.3"))
    ))

    # 检索配置
    default_max_results: int = field(default_factory=lambda: (
        _yaml_config.get("search", {}).get("default_max_results", 10)
    ))
    request_delay: float = field(default_factory=lambda: (
        _yaml_config.get("search", {}).get("request_delay", 0.34)
    ))
    sort_by: str = field(default_factory=lambda: (
        _yaml_config.get("search", {}).get("sort_by", "relevance")
    ))
    article_types: List[str] = field(default_factory=lambda: (
        _yaml_config.get("search", {}).get("article_types", [])
    ))
    min_citations: int = field(default_factory=lambda: (
        _yaml_config.get("search", {}).get("min_citations", 0)
    ))
    min_impact_factor: float = field(default_factory=lambda: (
        _yaml_config.get("search", {}).get("min_impact_factor", 0.0)
    ))
    open_access_only: bool = field(default_factory=lambda: (
        _yaml_config.get("search", {}).get("open_access_only", False)
    ))

    # 缓存配置
    cache_dir: str = field(default_factory=lambda: os.path.expanduser(
        _yaml_config.get("cache", {}).get("dir", "~/.gene_literature_cache")
    ))
    citation_ttl_days: int = field(default_factory=lambda: (
        _yaml_config.get("cache", {}).get("citation_ttl_days", 30)
    ))
    conclusion_ttl_days: int = field(default_factory=lambda: (
        _yaml_config.get("cache", {}).get("conclusion_ttl_days", 7)
    ))

    # 报告配置
    report_output_dir: str = field(default_factory=lambda: os.path.expanduser(
        _yaml_config.get("report", {}).get("output_dir", "~/gene_reports")
    ))

    def __post_init__(self):
        """确保目录存在"""
        os.makedirs(self.cache_dir, exist_ok=True)
        os.makedirs(self.report_output_dir, exist_ok=True)

    @property
    def is_ncbi_configured(self) -> bool:
        """检查NCBI是否已配置"""
        return bool(self.ncbi_email)

    @property
    def is_llm_configured(self) -> bool:
        """检查LLM是否已配置"""
        return bool(self.llm_api_key)

    # API URL 属性
    @property
    def crossref_url(self) -> str:
        return get_api_url("crossref")

    @property
    def openalex_url(self) -> str:
        return get_api_url("openalex")

    @property
    def ncbi_url(self) -> str:
        return get_api_url("ncbi")

    @property
    def mygene_url(self) -> str:
        return get_api_url("mygene")

    @property
    def kegg_url(self) -> str:
        return get_api_url("kegg")

    @property
    def reactome_url(self) -> str:
        return get_api_url("reactome")


# 提示词模板已移至 prompts/ 模块
# 为保持向后兼容，从 prompts 模块重新导出
from .prompts import (
    EXTRACTION_PROMPT_TEMPLATE,
    CONFLICT_DETECTION_PROMPT,
    SUMMARY_PROMPT_TEMPLATE,
    GENE_SUMMARY_PROMPT_TEMPLATE,
)

# 默认配置实例
default_config = Config()
