"""
文本处理工具模块 - 基因功能语境解析等
"""
import re
import logging
from typing import Tuple, List, Optional

from ..constants import (
    WILDTYPE_KEYWORDS, MUTANT_KEYWORDS, MISSING_VALUE,
    NOT_MENTIONED, INVALID_VALUES
)

logger = logging.getLogger(__name__)


def split_wildtype_mutant(context: str) -> Tuple[str, str]:
    """
    智能分割野生型和突变型描述

    支持多种格式：
    1. 标准格式：野生型：xxx；突变型：yyy
    2. 带基因名：野生型TP53：xxx；突变型TP53：yyy
    3. 隐式对比：正常功能xxx，突变后yyy
    4. 英文格式：wild-type: xxx; mutant: yyy
    5. 无明确分隔：野生型xxx突变型yyy

    Args:
        context: 原始语境文本

    Returns:
        (野生型描述, 突变型描述)
    """
    if not context or context.strip() in INVALID_VALUES:
        return MISSING_VALUE, MISSING_VALUE

    context = context.strip()
    wildtype_desc = ""
    mutant_desc = ""

    # 策略1: 标准格式匹配（带可选基因名）
    # 匹配 "野生型[基因名]：功能描述"
    wt_match = re.search(
        r'野生型[^：:\n]*[：:]\s*(.+?)(?=突变型[^：:\n]*[：:]|$)',
        context, re.DOTALL
    )
    mt_match = re.search(
        r'突变型[^：:\n]*[：:]\s*(.+?)(?=野生型[^：:\n]*[：:]|$)',
        context, re.DOTALL
    )

    if wt_match:
        wildtype_desc = _clean_description(wt_match.group(1))
    if mt_match:
        mutant_desc = _clean_description(mt_match.group(1))

    # 策略2: 英文格式匹配
    if not wildtype_desc and not mutant_desc:
        wt_match = re.search(
            r'wild-?type[^：:\n]*[：:]\s*(.+?)(?=mutant|$)',
            context, re.IGNORECASE | re.DOTALL
        )
        mt_match = re.search(
            r'mutant[^：:\n]*[：:]\s*(.+?)(?=wild-?type|$)',
            context, re.IGNORECASE | re.DOTALL
        )

        if wt_match:
            wildtype_desc = _clean_description(wt_match.group(1))
        if mt_match:
            mutant_desc = _clean_description(mt_match.group(1))

    # 策略3: 隐式对比（"正常...，突变后..."）
    if not wildtype_desc and not mutant_desc:
        # 检测隐式野生型
        implicit_wt = re.search(
            r'(?:正常|生理|野生型)[^，,；;。]+?(?=，|,|；|;|突变|$)',
            context
        )
        implicit_mt = re.search(
            r'(?:突变|变异)[^，,；;。]*?(?:后|导致|引起|使)?[^，,；;。]+',
            context
        )

        if implicit_wt:
            wildtype_desc = _clean_description(implicit_wt.group())
        if implicit_mt:
            mutant_desc = _clean_description(implicit_mt.group())

    # 策略4: 按分隔符切分并分类
    if not wildtype_desc and not mutant_desc:
        parts = re.split(r'[;；。，,]', context)

        for part in parts:
            part = part.strip()
            if not part:
                continue

            # 判断属于野生型还是突变型
            is_wildtype = any(kw in part.lower() for kw in [kw.lower() for kw in WILDTYPE_KEYWORDS])
            is_mutant = any(kw in part.lower() for kw in [kw.lower() for kw in MUTANT_KEYWORDS])

            # 排除同时包含两者的部分
            if is_wildtype and not is_mutant:
                wildtype_desc = part
            elif is_mutant and not is_wildtype:
                mutant_desc = part
            elif is_mutant and is_wildtype:
                # 同时包含，尝试进一步分割
                sub_parts = re.split(r'[，,]', part)
                for sub in sub_parts:
                    sub = sub.strip()
                    if any(kw in sub for kw in WILDTYPE_KEYWORDS):
                        wildtype_desc = sub
                    elif any(kw in sub for kw in MUTANT_KEYWORDS):
                        mutant_desc = sub

    # 策略5: 仅检测到一种类型时的推断
    if wildtype_desc and not mutant_desc:
        # 尝试从剩余文本提取突变型
        remaining = context.replace(wildtype_desc, "").strip()
        if remaining and any(kw in remaining for kw in MUTANT_KEYWORDS):
            # 清理分隔符
            remaining = re.sub(r'^[：:，,；;。]+', '', remaining)
            remaining = re.sub(r'[：:，,；;。]+$', '', remaining)
            if remaining:
                mutant_desc = _clean_description(remaining)

    if mutant_desc and not wildtype_desc:
        # 尝试从剩余文本提取野生型
        remaining = context.replace(mutant_desc, "").strip()
        if remaining and any(kw in remaining for kw in WILDTYPE_KEYWORDS):
            remaining = re.sub(r'^[：:，,；;。]+', '', remaining)
            remaining = re.sub(r'[：:，,；;。]+$', '', remaining)
            if remaining:
                wildtype_desc = _clean_description(remaining)

    return wildtype_desc or MISSING_VALUE, mutant_desc or MISSING_VALUE


def _clean_description(desc: str) -> str:
    """
    清理描述文本

    Args:
        desc: 原始描述

    Returns:
        清理后的描述
    """
    if not desc:
        return ""

    # 去除首尾空白和标点
    desc = desc.strip()
    desc = re.sub(r'^[：:，,；;。\s]+', '', desc)
    desc = re.sub(r'[：:，,；;。\s]+$', '', desc)

    # 如果过长，按分隔符截断
    if len(desc) > 200:
        parts = re.split(r'[;；。]', desc, maxsplit=1)
        desc = parts[0].strip()

    # 移除重复的关键词前缀（如开头的"野生型"或"突变型"）
    desc = re.sub(r'^(野生型|突变型)[^：:]*[：:]\s*', '', desc)

    return desc


def extract_disease_names(text: str, disease_keywords: List[str] = None) -> List[str]:
    """
    从文本中提取疾病名称

    Args:
        text: 输入文本
        disease_keywords: 疾病关键词列表（默认使用常量中的列表）

    Returns:
        提取到的疾病列表（去重，按优先级）
    """
    from ..constants import DISEASE_KEYWORDS

    keywords = disease_keywords or DISEASE_KEYWORDS
    found_diseases = set()

    for disease in keywords:
        if disease in text:
            # 检查是否已被更具体的疾病名称包含
            is_subsumed = False
            for existing in list(found_diseases):
                if disease in existing:
                    is_subsumed = True
                    break
                elif existing in disease:
                    found_diseases.remove(existing)

            if not is_subsumed:
                found_diseases.add(disease)

    # 按原始关键词顺序排序
    return [d for d in keywords if d in found_diseases]


def normalize_gene_context(context: str) -> str:
    """
    规范化基因功能语境

    将各种格式统一为："野生型：xxx；突变型：yyy"

    Args:
        context: 原始语境

    Returns:
        规范化后的语境
    """
    wt, mt = split_wildtype_mutant(context)

    if wt == MISSING_VALUE and mt == MISSING_VALUE:
        return context  # 无法解析，返回原值

    parts = []
    if wt != MISSING_VALUE:
        parts.append(f"野生型：{wt}")
    if mt != MISSING_VALUE:
        parts.append(f"突变型：{mt}")

    return "；".join(parts)


def detect_function_type(text: str) -> str:
    """
    检测文本描述的是野生型还是突变型功能

    Args:
        text: 文本描述

    Returns:
        "wildtype", "mutant", "both", 或 "unknown"
    """
    text_lower = text.lower()

    has_wildtype = any(kw.lower() in text_lower for kw in WILDTYPE_KEYWORDS)
    has_mutant = any(kw.lower() in text_lower for kw in MUTANT_KEYWORDS)

    if has_wildtype and has_mutant:
        return "both"
    elif has_wildtype:
        return "wildtype"
    elif has_mutant:
        return "mutant"
    else:
        return "unknown"


# ============================================================================
# 测试
# ============================================================================

if __name__ == "__main__":
    # 测试用例
    test_cases = [
        # 标准格式
        "野生型STAT3：调控炎症与免疫应答，维持免疫稳态；突变型STAT3：超激活或失活均导致免疫缺陷、自身免疫病或癌症。",
        # 简单格式
        "野生型：正常功能；突变型：致癌",
        # 带基因名
        "野生型TP53：作为抑癌基因，调控细胞周期和DNA修复；突变型TP53：丧失抑癌功能，获得促癌功能",
        # 仅突变型
        "突变型TP53：促进肿瘤进展",
        # 英文格式
        "Wild-type: tumor suppressor; Mutant: oncogenic function",
        # 隐式对比
        "正常BRCA1修复DNA损伤，突变后导致癌症风险增加",
        # 无冒号格式
        "野生型BRCA1修复DNA损伤，突变型导致癌症风险增加",
    ]

    print("=" * 70)
    print("野生型/突变型分割测试")
    print("=" * 70)

    for test in test_cases:
        wt, mt = split_wildtype_mutant(test)
        print(f"\n输入: {test[:50]}...")
        print(f"野生型: {wt}")
        print(f"突变型: {mt}")
