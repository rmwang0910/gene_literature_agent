"""
常量定义模块 - 从配置文件加载常量值
"""
from .config import (
    get_wildtype_keywords,
    get_mutant_keywords,
    get_disease_keywords,
    get_conflict_dimensions,
    get_core_semantic_opposites,
    get_impact_factor_tiers,
    get_confidence_levels,
)

# ============================================================================
# 缺失值表示（不可配置，系统内部使用）
# ============================================================================
MISSING_VALUE = "-"
NOT_MENTIONED = "未提及"
PARSE_FAILED = "解析失败"

# 所有被视为"无效"的值列表
INVALID_VALUES = {"", "-", "未提及", "解析失败", "无", "None", "null"}

# ============================================================================
# 从配置文件加载的关键词（可配置）
# ============================================================================

# 野生型关键词
WILDTYPE_KEYWORDS = get_wildtype_keywords()

# 突变型关键词
MUTANT_KEYWORDS = get_mutant_keywords()

# 疾病关键词（按优先级排序）
DISEASE_KEYWORDS = get_disease_keywords()

# 冲突维度定义
CONFLICT_DIMENSIONS = get_conflict_dimensions()

# 核心语义对立词组
CORE_SEMANTIC_OPPOSITES = get_core_semantic_opposites()

# ============================================================================
# 置信度和影响力配置
# ============================================================================

# 置信度等级阈值
CONFIDENCE_LEVELS = get_confidence_levels()

# 影响因子阈值
IMPACT_FACTOR_TIERS = get_impact_factor_tiers()
