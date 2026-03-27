"""
报告生成模块 - 生成 Markdown、PDF、HTML 和 Excel 报告
"""
import os
import logging
from typing import List, Dict, Optional, Any
from datetime import datetime
from pathlib import Path
import json

from ..config import default_config
from ..prompts import SUMMARY_PROMPT_TEMPLATE, GENE_SUMMARY_PROMPT_TEMPLATE
from .conclusion_extractor import GeneConclusion, OpenAIProvider, LLMProvider
from .conflict_detector import ConflictReport, Conflict, ConclusionMerger
from .citation_fetcher import get_citation_stars, get_impact_badge, calculate_confidence_score
from ..utils.helpers import normalize_pathways
from ..utils.text_utils import split_wildtype_mutant, detect_function_type
from ..constants import MISSING_VALUE, NOT_MENTIONED, INVALID_VALUES
from collections import Counter

# 配置日志
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def build_mutation_disease_matrix(
    conclusions: List[GeneConclusion]
) -> Dict[str, Dict[str, List[str]]]:
    """
    构建突变-疾病关联矩阵

    Args:
        conclusions: 结论列表

    Returns:
        {突变位点: {疾病类型: [PMID列表]}}
    """
    matrix = {}

    for c in conclusions:
        for assoc in c.mutation_disease_associations:
            mutation = assoc.get("mutation", "")
            disease = assoc.get("disease", "")

            if mutation and disease:
                # 规范化突变位点
                from ..utils.helpers import normalize_mutation_site
                mutation = normalize_mutation_site(mutation)

                if mutation not in matrix:
                    matrix[mutation] = {}
                if disease not in matrix[mutation]:
                    matrix[mutation][disease] = []
                # 避免重复添加同一篇文献
                if c.pmid not in matrix[mutation][disease]:
                    matrix[mutation][disease].append(c.pmid)

    return matrix


def format_mutation_disease_matrix(
    matrix: Dict[str, Dict[str, List[str]]],
    min_count: int = 1
) -> str:
    """
    格式化突变-疾病关联矩阵为Markdown表格

    Args:
        matrix: 关联矩阵 {突变: {疾病: [PMID列表]}}
        min_count: 最小出现次数阈值

    Returns:
        Markdown表格字符串
    """
    if not matrix:
        return ""

    # 收集所有疾病类型
    all_diseases = set()
    for diseases in matrix.values():
        all_diseases.update(diseases.keys())

    if not all_diseases:
        return ""

    # 按出现次数排序疾病
    disease_counts = Counter()
    for diseases in matrix.values():
        for disease, pmids in diseases.items():
            disease_counts[disease] += len(pmids)
    sorted_diseases = [d for d, _ in disease_counts.most_common()]

    # 按总关联数排序突变
    mutation_totals = {
        mut: sum(len(pmids) for pmids in diseases.values())
        for mut, diseases in matrix.items()
    }
    sorted_mutations = sorted(
        matrix.keys(),
        key=lambda m: mutation_totals[m],
        reverse=True
    )

    # 构建表格
    lines = []
    header = "| 突变位点 | " + " | ".join(sorted_diseases) + " | 总计 |"
    separator = "| --- | " + " | ".join(["---"] * len(sorted_diseases)) + " | --- |"

    lines.append(header)
    lines.append(separator)

    for mutation in sorted_mutations:
        row = [f"**{mutation}**"]
        total = 0
        for disease in sorted_diseases:
            pmids = matrix[mutation].get(disease, [])
            count = len(pmids)
            if count >= min_count:
                # 生成带 PMID 链接的单元格
                pmid_links = ", ".join([
                    f"[{pmid}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/)"
                    for pmid in pmids
                ])
                row.append(pmid_links)
                total += count
            else:
                row.append("-")
        row.append(str(total))
        lines.append("| " + " | ".join(row) + " |")

    return "\n".join(lines)


def build_integrated_mutation_disease_heatmap(
    conclusions: List[GeneConclusion],
    external_evidence: Optional[Dict[str, Any]] = None
) -> str:
    """
    构建整合数据库与文献的突变-疾病关联热图（HTML表格）

    Args:
        conclusions: 文献结论列表
        external_evidence: 外部数据源证据

    Returns:
        HTML 热图表格字符串
    """
    # 收集所有突变-疾病关联
    # 结构: {突变: {疾病: {"count": int, "sources": [来源列表], "pmids": [PMID列表]}}}
    associations = {}

    # 1. 从文献结论中收集
    for c in conclusions:
        for assoc in c.mutation_disease_associations:
            mutation = assoc.get("mutation", "")
            disease = assoc.get("disease", "")

            if mutation and disease:
                from ..utils.helpers import normalize_mutation_site
                mutation = normalize_mutation_site(mutation)

                if mutation not in associations:
                    associations[mutation] = {}
                if disease not in associations[mutation]:
                    associations[mutation][disease] = {"count": 0, "sources": [], "pmids": []}

                associations[mutation][disease]["count"] += 1
                if "文献" not in associations[mutation][disease]["sources"]:
                    associations[mutation][disease]["sources"].append("文献")
                if c.pmid not in associations[mutation][disease]["pmids"]:
                    associations[mutation][disease]["pmids"].append(c.pmid)

    # 2. 从COSMIC热点突变收集
    if external_evidence:
        hotspots = external_evidence.get("hotspot_mutations", [])
        for m in hotspots:
            mutation = m.get("mutation", "")
            cancer_types = m.get("cancer_types", [])
            sample_count = m.get("sample_count", 0)

            if mutation and cancer_types:
                if mutation not in associations:
                    associations[mutation] = {}

                for cancer in cancer_types[:3]:  # 取前3个癌症类型
                    if cancer not in associations[mutation]:
                        associations[mutation][cancer] = {"count": 0, "sources": [], "pmids": []}

                    # 根据样本数估算关联强度 (log scale)
                    import math
                    strength = min(5, int(math.log10(sample_count + 1)) + 1) if sample_count > 0 else 1
                    associations[mutation][cancer]["count"] += strength
                    if "COSMIC" not in associations[mutation][cancer]["sources"]:
                        associations[mutation][cancer]["sources"].append("COSMIC")

        # 3. 从Open Targets疾病关联收集
        disease_assocs = external_evidence.get("disease_associations", [])
        related_cancers = external_evidence.get("related_cancer_types", [])

        # 将高关联度疾病与所有突变关联
        for d in disease_assocs[:5]:
            disease_name = d.get("disease_name", "")
            score = d.get("score", 0)
            if disease_name and score > 0.5:
                for mutation in associations.keys():
                    if disease_name not in associations[mutation]:
                        associations[mutation][disease_name] = {"count": 0, "sources": [], "pmids": []}
                    if "OpenTargets" not in associations[mutation][disease_name]["sources"]:
                        associations[mutation][disease_name]["sources"].append("OpenTargets")
                        # 根据评分添加关联强度
                        associations[mutation][disease_name]["count"] += int(score * 3)

    if not associations:
        return ""

    # 收集所有疾病，按关联总数排序
    all_diseases = set()
    for diseases in associations.values():
        all_diseases.update(diseases.keys())

    if not all_diseases:
        return ""

    disease_totals = Counter()
    for diseases in associations.values():
        for disease, data in diseases.items():
            disease_totals[disease] += data["count"]
    sorted_diseases = [d for d, _ in disease_totals.most_common()[:8]]  # 最多8列

    # 按关联总数排序突变
    mutation_totals = {}
    for mut, diseases in associations.items():
        mutation_totals[mut] = sum(d["count"] for d in diseases.values())
    sorted_mutations = sorted(
        associations.keys(),
        key=lambda m: mutation_totals[m],
        reverse=True
    )[:10]  # 最多10行

    if not sorted_mutations or not sorted_diseases:
        return ""

    # 计算最大值用于归一化
    max_count = max(
        associations[mut].get(dis, {}).get("count", 0)
        for mut in sorted_mutations
        for dis in sorted_diseases
    )
    if max_count == 0:
        max_count = 1

    # 构建HTML热图表格
    lines = []
    lines.append("<table style=\"border-collapse:collapse;width:100%;font-size:0.9em\">")

    # 表头
    lines.append("<tr>")
    lines.append("<th style=\"padding:8px;border:1px solid #ddd;background:#f5f5f5\">突变</th>")
    for disease in sorted_diseases:
        # 疾病名太长时截断
        short_disease = disease[:15] + "..." if len(disease) > 15 else disease
        lines.append(f"<th style=\"padding:8px;border:1px solid #ddd;background:#f5f5f5;writing-mode:vertical-rl;text-orientation:mixed;max-width:80px\">{short_disease}</th>")
    lines.append("</tr>")

    # 数据行
    for mutation in sorted_mutations:
        lines.append("<tr>")
        lines.append(f"<td style=\"padding:8px;border:1px solid #ddd;font-weight:bold\">{mutation}</td>")

        for disease in sorted_diseases:
            data = associations[mutation].get(disease, {"count": 0, "sources": [], "pmids": []})
            count = data["count"]
            sources = data["sources"]
            pmids = data["pmids"]

            if count > 0:
                # 计算热度颜色 (红色渐变)
                intensity = min(1.0, count / max_count)
                r, g, b = 255, int(255 * (1 - intensity * 0.7)), int(255 * (1 - intensity * 0.7))

                # 来源标记
                source_badges = []
                if "COSMIC" in sources:
                    source_badges.append("🔬")
                if "文献" in sources:
                    source_badges.append("📄")
                if "OpenTargets" in sources:
                    source_badges.append("🎯")

                badge_str = "".join(source_badges)

                # PMID链接（如果有）
                if pmids:
                    pmid_links = ", ".join([f"<a href='https://pubmed.ncbi.nlm.nih.gov/{p}/'>{p}</a>" for p in pmids[:2]])
                    cell_content = f"{badge_str}<br><small>{pmid_links}</small>"
                else:
                    cell_content = badge_str if badge_str else "●"

                lines.append(f"<td style=\"padding:8px;border:1px solid #ddd;background-color:rgb({r},{g},{b});text-align:center\">{cell_content}</td>")
            else:
                lines.append("<td style=\"padding:8px;border:1px solid #ddd;text-align:center;color:#ccc\">-</td>")

        lines.append("</tr>")

    lines.append("</table>")

    # 图例
    lines.append("")
    lines.append("<div style=\"margin-top:12px;font-size:0.85em;color:#666\">")
    lines.append("  <strong>图例:</strong> 🔬 COSMIC | 📄 文献 | 🎯 OpenTargets | 颜色深浅表示关联强度")
    lines.append("</div>")

    return "\n".join(lines)


class GeneSummaryGenerator:
    """基因总结生成器 - 使用LLM生成综合总结"""

    def __init__(self, provider: Optional[LLMProvider] = None):
        """
        初始化总结生成器

        Args:
            provider: LLM提供者
        """
        self.provider = provider or self._get_default_provider()

    def _get_default_provider(self) -> LLMProvider:
        """获取默认LLM提供者"""
        if default_config.llm_api_key:
            return OpenAIProvider()
        raise ValueError("未配置LLM API Key，无法生成基因总结")

    def _format_conclusions_for_prompt(
        self,
        conclusions: List[GeneConclusion]
    ) -> str:
        """格式化结论数据用于提示词"""
        formatted = []
        for i, c in enumerate(conclusions, 1):
            # 构建置信度显示
            conf_display = c.confidence
            if c.citation_count > 0 or c.impact_factor > 0:
                conf_result = calculate_confidence_score(
                    citation_count=c.citation_count,
                    impact_factor=c.impact_factor,
                    publication_year=c.year or 0,
                    llm_confidence=c.confidence
                )
                from .citation_fetcher import get_confidence_stars
                conf_display = f"{c.confidence} ({get_confidence_stars(conf_result['score'])})"

            # 基因功能语境标签
            func_tag = ""
            if c.gene_function_context:
                if "野生型" in c.gene_function_context and "突变型" not in c.gene_function_context:
                    func_tag = '<span class="tag tag-wildtype">野生型</span>'
                elif "突变型" in c.gene_function_context and "野生型" not in c.gene_function_context:
                    func_tag = '<span class="tag tag-mutant">突变型</span>'
                elif "野生型" in c.gene_function_context and "突变型" in c.gene_function_context:
                    func_tag = '<span class="tag tag-wildtype">野生型</span> <span class="tag tag-mutant">突变型</span>'
                else:
                    func_tag = '<span class="tag tag-unknown">未区分</span>'

            # PMID 链接
            pmid_link = f'<a href="https://pubmed.ncbi.nlm.nih.gov/{c.pmid}/" target="_blank" class="pmid-link">PMID: {c.pmid} ↗</a>'

            formatted.append(f"""
### 文献 {i} ({pmid_link})
{func_tag}
- 疾病关系: {c.disease_relation or '未提及'}
- 信号通路: {', '.join(c.pathways) if c.pathways else '未提及'}
- 临床意义: {c.clinical_significance or '未提及'}
- 突变效应: {c.mutation_effects or '未提及'}
- 基因功能语境: {c.gene_function_context or '未提及'}
- 置信度: {conf_display}
- 发表年份: {c.year or '未知'}
- 期刊: {c.journal or '未知'} {get_impact_badge(c.impact_factor)}
- 引用次数: {c.citation_count}
""")
        return "\n".join(formatted)

    def _format_conflicts_for_prompt(
        self,
        conflict_report: ConflictReport
    ) -> str:
        """格式化冲突数据用于提示词"""
        if not conflict_report.conflicts:
            return "未检测到明显冲突"

        formatted = []
        for i, c in enumerate(conflict_report.conflicts, 1):
            formatted.append(f"""
冲突 {i}: {c.topic}
- PMID {c.pmid_a}: {c.conclusion_a}
- PMID {c.pmid_b}: {c.conclusion_b}
- 冲突维度: {c.conflict_dimension or '未分类'}
- 建议: {c.resolution_suggestion or '需要进一步验证'}
""")
        return "\n".join(formatted)

    def generate_summary(
        self,
        gene: str,
        conclusions: List[GeneConclusion],
        conflict_report: ConflictReport
    ) -> str:
        """
        生成基因综合总结

        Args:
            gene: 基因名
            conclusions: 结论列表
            conflict_report: 冲突报告

        Returns:
            基因总结文本
        """
        if not conclusions:
            return f"暂无{gene}相关文献数据。"

        # 格式化数据
        conclusions_data = self._format_conclusions_for_prompt(conclusions)
        conflicts_data = self._format_conflicts_for_prompt(conflict_report)

        # 构建提示词
        prompt = GENE_SUMMARY_PROMPT_TEMPLATE.format(
            gene=gene,
            conclusions_data=conclusions_data,
            conflicts_data=conflicts_data
        )

        try:
            response = self.provider.generate(
                prompt,
                max_tokens=500,
                temperature=0.3
            )
            return response.strip()
        except Exception as e:
            logger.warning(f"LLM生成基因总结失败: {e}")
            # 降级：生成简单总结
            return self._generate_fallback_summary(gene, conclusions, conflict_report)

    def _generate_fallback_summary(
        self,
        gene: str,
        conclusions: List[GeneConclusion],
        conflict_report: ConflictReport
    ) -> str:
        """生成降级总结（当LLM不可用时）"""
        # 收集关键信息
        diseases = set()
        pathways = set()
        clinical_sigs = set()
        mutations = set()
        wildtype_funcs = []
        mutant_funcs = []

        for c in conclusions:
            if c.disease_relation:
                diseases.add(c.disease_relation)
            # 规范化通路名称后再添加
            pathways.update(normalize_pathways(c.pathways))
            if c.clinical_significance:
                clinical_sigs.add(c.clinical_significance)
            if c.mutation_effects and c.mutation_effects != "未提及":
                mutations.add(c.mutation_effects)
            if c.gene_function_context:
                if "野生型" in c.gene_function_context:
                    wildtype_funcs.append(c.gene_function_context)
                elif "突变型" in c.gene_function_context:
                    mutant_funcs.append(c.gene_function_context)

        # 构建总结
        parts = [f"**{gene}**是重要的肿瘤相关基因。"]

        if wildtype_funcs or mutant_funcs:
            func_desc = []
            if wildtype_funcs:
                func_desc.append("野生型具有正常生理功能")
            if mutant_funcs:
                func_desc.append("突变型与肿瘤发生发展相关")
            parts.append("，".join(func_desc) + "。")

        if diseases:
            disease_str = "、".join(list(diseases)[:3])
            parts.append(f"与{disease_str}等疾病相关。")

        if pathways:
            pathway_str = "、".join(list(pathways)[:5])
            parts.append(f"主要涉及{pathway_str}等信号通路。")

        if clinical_sigs:
            parts.append(f"临床意义：{list(clinical_sigs)[0]}。")

        if conflict_report.conflicts:
            parts.append(f"文献结论存在{len(conflict_report.conflicts)}处争议，需进一步验证。")

        return "".join(parts)


class ReportGenerator:
    """报告生成器"""

    def __init__(self, output_dir: Optional[str] = None, enable_external_sources: bool = True):
        """
        初始化报告生成器

        Args:
            output_dir: 输出目录
            enable_external_sources: 是否启用外部数据源（Open Targets, ClinVar, COSMIC, Pathway）
        """
        self.output_dir = Path(output_dir or default_config.report_output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.merger = ConclusionMerger()
        self.enable_external_sources = enable_external_sources
        self._evidence_enricher = None

    def _get_evidence_enricher(self):
        """获取或创建证据整合器"""
        if self._evidence_enricher is None and self.enable_external_sources:
            try:
                from ..datasources.evidence_enricher import EvidenceEnricher
                self._evidence_enricher = EvidenceEnricher()
            except ImportError as e:
                logger.warning(f"无法加载外部数据源模块: {e}")
        return self._evidence_enricher

    def enrich_with_external_sources(
        self,
        gene: str,
        variant: Optional[str] = None
    ) -> Optional[Dict[str, Any]]:
        """
        从外部数据源获取证据

        Args:
            gene: 基因符号
            variant: 可选的变异描述

        Returns:
            整合后的证据字典
        """
        enricher = self._get_evidence_enricher()
        if not enricher:
            return None

        try:
            logger.info(f"正在从外部数据源获取 {gene} 的证据...")
            evidence = enricher.enrich_gene(gene, variant)
            logger.info(f"外部数据源查询完成，数据来源: {', '.join(evidence.sources_queried)}")
            return {
                "gene": evidence.gene,
                "is_cancer_gene": evidence.is_cancer_gene,
                "cancer_gene_role": evidence.cancer_gene_role,
                "disease_associations": evidence.disease_associations,
                "drugs": evidence.drugs,
                "opentargets_pmids": evidence.opentargets_pmids,
                "pathogenic_variants": evidence.pathogenic_variants,
                "variant_significance": evidence.variant_significance,
                "hotspot_mutations": evidence.hotspot_mutations,
                "related_cancer_types": evidence.related_cancer_types,
                "pathways": evidence.pathways,
                "sources_queried": evidence.sources_queried
            }
        except Exception as e:
            logger.warning(f"获取外部数据源证据失败: {e}")
            return None

    # _format_external_evidence_section 已废弃，外部证据已整合到 generate_markdown 各章节中

    def enrich_conclusions_with_citations(
        self,
        conclusions: List[GeneConclusion],
        delay: float = 0.15
    ) -> List[GeneConclusion]:
        """
        为结论添加引用次数和期刊信息

        Args:
            conclusions: 结论列表
            delay: API 请求间隔

        Returns:
            丰富后的结论列表
        """
        from .citation_fetcher import CitationEnricher

        if not conclusions:
            return conclusions

        # 获取所有 PMID
        pmids = [c.pmid for c in conclusions]
        doi_map = {c.pmid: c.doi for c in conclusions if c.doi}

        logger.info(f"正在获取 {len(pmids)} 篇文献的引用数据...")

        # 批量获取引用信息
        enricher = CitationEnricher()
        citation_info = enricher.enrich_pmids(pmids, doi_map, delay)

        # 更新结论
        for c in conclusions:
            info = citation_info.get(c.pmid)
            if info:
                c.citation_count = info.citation_count
                c.journal = info.journal
                c.impact_factor = info.impact_factor
                if info.doi and not c.doi:
                    c.doi = info.doi
                if info.year and not c.year:
                    c.year = info.year
                if info.publication_type:
                    c.publication_type = info.publication_type

                # 计算综合置信度分数
                conf_result = calculate_confidence_score(
                    citation_count=info.citation_count,
                    impact_factor=info.impact_factor,
                    publication_year=info.year or c.year,
                    llm_confidence=c.confidence
                )
                # 更新置信度为计算后的等级
                c.confidence = conf_result["level"]

        logger.info("引用数据获取完成")
        return conclusions

    def _format_reference_card(
        self,
        c: GeneConclusion,
        compact: bool = False
    ) -> List[str]:
        """
        格式化单篇文献为卡片样式

        Args:
            c: 结论对象
            compact: 是否使用紧凑模式（普通文献）

        Returns:
            Markdown 行列表
        """
        lines = []
        from .citation_fetcher import get_normalized_score

        title = c.title if c.title else "无标题"

        # 文献类型标识
        type_tag = ""
        type_color = "#666"
        if c.publication_type == "review":
            type_tag = "综述"
            type_color = "#9c27b0"
        elif c.publication_type == "research":
            type_tag = "原创研究"
            type_color = "#1976d2"

        # 获取星级
        stars = get_citation_stars(c.citation_count, c.impact_factor, c.year)
        impact_badge = get_impact_badge(c.impact_factor)

        # 计算年均引用
        current_year = datetime.now().year
        if c.year and c.year > 0 and c.year <= current_year:
            years = max(1, current_year - c.year)
            annual_citations = c.citation_count / years
            annual_str = f"{annual_citations:.1f}"
        else:
            annual_str = "-"

        # 影响力评分
        normalized_score = get_normalized_score(c.citation_count, c.impact_factor, c.year)

        if compact:
            # 紧凑模式：单行卡片
            journal_info = c.journal if c.journal else "未知期刊"
            if c.impact_factor > 0:
                journal_info += f" (IF: {c.impact_factor:.1f})"

            lines.append(f"<div style=\"padding:12px;margin:8px 0;background:#f8f9fa;border-left:3px solid #dee2e6;border-radius:4px\">")
            lines.append(f"  <div style=\"display:flex;justify-content:space-between;align-items:center\">")
            lines.append(f"    <div><a href=\"https://pubmed.ncbi.nlm.nih.gov/{c.pmid}/\" target=\"_blank\"><strong>PMID: {c.pmid}</strong></a> {title[:60]}{'...' if len(title) > 60 else ''}</div>")
            lines.append(f"    <div style=\"color:#666;font-size:0.85em\">{stars}</div>")
            lines.append(f"  </div>")
            lines.append(f"  <div style=\"font-size:0.85em;color:#666;margin-top:4px\">{journal_info} | {c.year} | 引用: {c.citation_count}</div>")
            lines.append(f"</div>")
        else:
            # 完整卡片模式：高影响力文献
            # 金色边框表示高影响力，绿色表示中等影响力
            border_color = "#ffc107"  # 金色
            bg_color = "#fffde7"  # 浅黄背景

            lines.append(f"<div style=\"padding:16px;margin:12px 0;background:{bg_color};border-left:4px solid {border_color};border-radius:8px;box-shadow:0 1px 3px rgba(0,0,0,0.1)\">")

            # 标题行
            lines.append(f"  <div style=\"display:flex;justify-content:space-between;align-items:flex-start\">")
            lines.append(f"    <div style=\"flex:1\">")
            lines.append(f"      <a href=\"https://pubmed.ncbi.nlm.nih.gov/{c.pmid}/\" target=\"_blank\" style=\"font-weight:600;color:#1976d2;text-decoration:none\">PMID: {c.pmid}</a>")
            lines.append(f"      <div style=\"font-size:1.05em;margin:4px 0\">{title[:80]}{'...' if len(title) > 80 else ''}</div>")
            lines.append(f"    </div>")
            lines.append(f"    <div style=\"text-align:right;min-width:80px\">")
            lines.append(f"      <div style=\"font-size:1.2em\">{stars}</div>")
            lines.append(f"      <div style=\"font-size:0.75em;color:#888\">评分: {normalized_score:.0f}</div>")
            lines.append(f"    </div>")
            lines.append(f"  </div>")

            # 元信息行
            meta_parts = []
            if c.journal:
                meta_parts.append(f"<strong>{c.journal}</strong>")
            if c.year:
                meta_parts.append(f"{c.year}年")
            if c.impact_factor > 0:
                meta_parts.append(f"IF: {c.impact_factor:.1f}")
            meta_str = " | ".join(meta_parts) if meta_parts else ""

            lines.append(f"  <div style=\"font-size:0.9em;color:#555;margin:8px 0\">{meta_str}</div>")

            # 统计行
            stats_parts = [f"📖 引用: {c.citation_count}"]
            if annual_str != "-":
                stats_parts.append(f"📈 年均: {annual_str}")
            if type_tag:
                stats_parts.append(f"<span style=\"color:{type_color}\">📄 {type_tag}</span>")
            if impact_badge:
                stats_parts.append(impact_badge)

            lines.append(f"  <div style=\"font-size:0.85em;color:#666\">{' | '.join(stats_parts)}</div>")
            lines.append(f"</div>")

        return lines

    def generate_markdown(
        self,
        gene: str,
        conclusions: List[GeneConclusion],
        conflict_report: ConflictReport,
        articles: Optional[List[Dict]] = None,
        include_summary: bool = True,
        variant: Optional[str] = None
    ) -> str:
        """
        生成 Markdown 报告（重构版 - 整合相似内容）

        Args:
            gene: 基因名
            conclusions: 结论列表
            conflict_report: 冲突报告
            articles: 原始文章信息列表
            include_summary: 是否包含总结
            variant: 可选的变异描述（用于查询外部数据源）

        Returns:
            Markdown 文本
        """
        merged = self.merger.merge_conclusions(conclusions, conflict_report)

        # 获取外部数据源证据
        external_evidence = None
        if self.enable_external_sources:
            external_evidence = self.enrich_with_external_sources(gene, variant)

        lines = []

        # =====================================================================
        # 1. 标题和元信息
        # =====================================================================
        lines.append(f"# 基因文献分析报告: {gene}")
        lines.append("")
        lines.append(f"> 生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        lines.append(f"> 文献数量: {merged['total_papers']}")
        if external_evidence:
            lines.append(f"> 数据来源: 文献 + {', '.join(external_evidence.get('sources_queried', []))}")
        lines.append("")

        # =====================================================================
        # 2. 概述
        # =====================================================================
        lines.append("## 概述")
        lines.append("")

        # LLM生成的基因总结
        if include_summary and conclusions:
            try:
                summary_generator = GeneSummaryGenerator()
                gene_summary = summary_generator.generate_summary(gene, conclusions, conflict_report)
                lines.append(gene_summary)
                lines.append("")
            except Exception as e:
                logger.warning(f"生成基因总结失败: {e}")

        if conflict_report.needs_review:
            lines.append("> ⚠️ **注意**: 检测到文献结论存在冲突，建议仔细核实。")
            lines.append("")

        # =====================================================================
        # 3. 基因功能（整合野生型/突变型功能 + 突变效应）
        # =====================================================================
        lines.append("## 基因功能")
        lines.append("")

        # 3.1 功能语境（野生型 vs 突变型）
        all_contexts = []
        for c in conclusions:
            if c.gene_function_context:
                all_contexts.append((c.gene_function_context, c))

        if all_contexts:
            lines.append("### 野生型与突变型功能对比")
            lines.append("")
            lines.append("| 野生型功能 | 突变型功能 | 来源 |")
            lines.append("| --- | --- | --- |")
            for context, c in all_contexts:
                wildtype_desc, mutant_desc = split_wildtype_mutant(context)
                stars = get_citation_stars(c.citation_count, c.impact_factor, c.year)
                lines.append(f"| {wildtype_desc} | {mutant_desc} | [{c.pmid}](https://pubmed.ncbi.nlm.nih.gov/{c.pmid}/) {stars} |")
            lines.append("")

        # 3.2 突变效应详情
        if merged["mutations"]:
            lines.append("### 突变效应")
            lines.append("")
            lines.append("| 效应描述 | 来源 |")
            lines.append("| --- | --- |")
            for mut in merged["mutations"]:
                lines.append(f"| {mut['effect']} | [{mut['pmid']}](https://pubmed.ncbi.nlm.nih.gov/{mut['pmid']}/) |")
            lines.append("")

        if not all_contexts and not merged["mutations"]:
            lines.append("暂无相关基因功能信息")
            lines.append("")

        # =====================================================================
        # 4. 疾病关联（整合 Open Targets + 文献 + 突变-疾病矩阵）
        # =====================================================================
        lines.append("## 疾病关联")
        lines.append("")

        # 4.1 数据库来源
        disease_assocs = external_evidence.get("disease_associations", []) if external_evidence else []
        has_db_disease = bool(disease_assocs) or (external_evidence and external_evidence.get("is_cancer_gene"))

        if has_db_disease:
            lines.append("### 数据库来源")
            lines.append("")

            # Open Targets 基因-疾病关联
            if disease_assocs:
                lines.append("#### Open Targets 基因-疾病关联")
                lines.append("")
                lines.append("| 疾病 | 关联评分 | 证据类型 | 详情 |")
                lines.append("| --- | --- | --- | --- |")
                for d in disease_assocs[:10]:
                    score = d.get("score", 0)
                    score_bar = "█" * int(score * 10) + "░" * (10 - int(score * 10))
                    evidence_types = d.get("evidence_types", [])
                    evidence_str = ", ".join(evidence_types[:3]) if evidence_types else "-"
                    url = d.get("url", "")
                    lines.append(f"| {d.get('disease_name', '')} | {score_bar} {score:.2f} | {evidence_str} | [详情]({url}) |")
                lines.append("")

            # COSMIC 癌基因信息
            if external_evidence and external_evidence.get("is_cancer_gene"):
                lines.append("#### COSMIC 癌基因信息")
                lines.append("")
                role = external_evidence.get("cancer_gene_role", "")
                role_cn = {"oncogene": "癌基因", "TSG": "抑癌基因", "fusion": "融合基因"}.get(role, role)
                cancer_types = external_evidence.get("related_cancer_types", [])
                lines.append(f"- **角色**: {role_cn}")
                if cancer_types:
                    lines.append(f"- **相关癌症类型**: {', '.join(cancer_types[:8])}")
                lines.append("")

        # 4.2 文献来源
        all_relations = [(c.disease_relation, c) for c in conclusions if c.disease_relation]
        if all_relations:
            lines.append("### 文献来源")
            lines.append("")
            lines.append("| 疾病关系描述 | 来源 |")
            lines.append("| --- | --- |")
            for relation, c in all_relations:
                # 清理野生型/突变型前缀，直接展示核心内容
                clean_relation = relation
                for prefix in ["突变型", "野生型"]:
                    if prefix in clean_relation:
                        parts = clean_relation.split(prefix)
                        if len(parts) > 1:
                            clean_relation = parts[-1].strip()
                            if clean_relation.startswith("：") or clean_relation.startswith(":"):
                                clean_relation = clean_relation[1:].strip()
                stars = get_citation_stars(c.citation_count, c.impact_factor, c.year)
                lines.append(f"| {clean_relation} | [{c.pmid}](https://pubmed.ncbi.nlm.nih.gov/{c.pmid}/) {stars} |")
            lines.append("")

        # 4.3 综合分析：突变-疾病关联热图
        heatmap_html = build_integrated_mutation_disease_heatmap(conclusions, external_evidence)
        if heatmap_html:
            lines.append("### 综合分析：突变-疾病关联热图")
            lines.append("")
            lines.append("> 整合 COSMIC、OpenTargets 数据库与文献报道的突变-疾病关联")
            lines.append("")
            lines.append(heatmap_html)
            lines.append("")

        if not has_db_disease and not all_relations:
            lines.append("暂无疾病关联信息")
            lines.append("")

        # =====================================================================
        # 5. 变异与临床意义（整合 COSMIC热点 + ClinVar + 文献临床意义）
        # =====================================================================
        lines.append("## 变异与临床意义")
        lines.append("")

        # 5.1 数据库来源
        hotspots = external_evidence.get("hotspot_mutations", []) if external_evidence else []
        pathogenic = external_evidence.get("pathogenic_variants", []) if external_evidence else []
        has_db_variants = bool(hotspots) or bool(pathogenic)

        if has_db_variants:
            lines.append("### 数据库来源")
            lines.append("")

            # COSMIC 热点突变
            if hotspots:
                lines.append("#### COSMIC 热点突变")
                lines.append("")
                lines.append("| 突变 | 样本数 | 频率 | 相关癌症 |")
                lines.append("| --- | --- | --- | --- |")
                for m in hotspots[:8]:
                    freq = f"{m.get('frequency', 0) * 100:.1f}%" if m.get('frequency') else "-"
                    cancer_types = ", ".join(m.get("cancer_types", [])[:3]) if m.get("cancer_types") else "-"
                    hotspot_tag = "🔥 " if m.get("is_hotspot") else ""
                    lines.append(f"| {hotspot_tag}**{m.get('mutation', '')}** | {m.get('sample_count', 0)} | {freq} | {cancer_types} |")
                lines.append("")

            # ClinVar 致病变异
            if pathogenic:
                lines.append("#### ClinVar 致病变异")
                lines.append("")
                lines.append("| 变异 | 临床意义 | 星级 | 详情 |")
                lines.append("| --- | --- | --- | --- |")
                for v in pathogenic[:8]:
                    stars = "★" * v.get("stars", 0) if v.get("stars", 0) > 0 else "-"
                    sig = v.get('significance_cn', v.get('significance', '')) or '-'
                    lines.append(f"| {v.get('variant', '')[:45]} | {sig} | {stars} | [ClinVar]({v.get('url', '')}) |")
                lines.append("")

        # 5.2 文献来源
        all_clinical = [(c.clinical_significance, c) for c in conclusions if c.clinical_significance]
        if all_clinical:
            lines.append("### 文献来源")
            lines.append("")
            lines.append("| 临床意义 | 来源 |")
            lines.append("| --- | --- |")
            for sig, c in all_clinical:
                stars = get_citation_stars(c.citation_count, c.impact_factor, c.year)
                # 清理野生型/突变型前缀，保留核心内容
                clean_sig = sig.replace("突变型", "").replace("野生型", "").strip()
                if clean_sig.startswith("：") or clean_sig.startswith(":"):
                    clean_sig = clean_sig[1:].strip()
                lines.append(f"| {clean_sig} | [{c.pmid}](https://pubmed.ncbi.nlm.nih.gov/{c.pmid}/) {stars} |")
            lines.append("")

        if not has_db_variants and not all_clinical:
            lines.append("暂无变异与临床意义相关信息")
            lines.append("")

        # =====================================================================
        # 6. 信号通路（整合 KEGG/Reactome + 文献）
        # =====================================================================
        lines.append("## 信号通路")
        lines.append("")

        # 6.1 数据库来源 - 双栏卡片布局
        db_pathways = external_evidence.get("pathways", []) if external_evidence else []
        kegg_pathways = [p for p in db_pathways if p.get("source") == "KEGG"][:10]
        reactome_pathways = [p for p in db_pathways if p.get("source") == "Reactome"][:10]
        has_db_pathways = bool(kegg_pathways) or bool(reactome_pathways)

        if has_db_pathways:
            lines.append("### 数据库来源")
            lines.append("")

            # 双栏布局容器
            lines.append("<div style=\"display:flex;gap:24px;flex-wrap:wrap\">")
            lines.append("")

            # KEGG 通路（左栏）
            if kegg_pathways:
                lines.append("<div style=\"flex:1;min-width:280px\">")
                lines.append("<div style=\"background:#e3f2fd;padding:12px 16px;border-radius:8px 8px 0 0;border-bottom:2px solid #1976d2\">")
                lines.append(f"<strong style=\"color:#1565c0\">🧬 KEGG</strong> <span style=\"color:#666;font-size:0.85em\">({len(kegg_pathways)}个通路)</span>")
                lines.append("</div>")
                lines.append("<div style=\"background:#f8f9fa;padding:8px;border-radius:0 0 8px 8px\">")
                for p in kegg_pathways:
                    name = p.get('name', '')
                    url = p.get('url', '')
                    lines.append(f"<div style=\"padding:8px 12px;margin:4px 0;background:white;border-radius:6px;border-left:3px solid #1976d2\">")
                    lines.append(f"<a href=\"{url}\" target=\"_blank\" style=\"color:#333;text-decoration:none\">{name}</a>")
                    lines.append("</div>")
                lines.append("</div>")
                lines.append("</div>")
                lines.append("")

            # Reactome 通路（右栏）
            if reactome_pathways:
                lines.append("<div style=\"flex:1;min-width:280px\">")
                lines.append("<div style=\"background:#f3e5f5;padding:12px 16px;border-radius:8px 8px 0 0;border-bottom:2px solid #7b1fa2\">")
                lines.append(f"<strong style=\"color:#6a1b9a\">🔬 Reactome</strong> <span style=\"color:#666;font-size:0.85em\">({len(reactome_pathways)}个通路)</span>")
                lines.append("</div>")
                lines.append("<div style=\"background:#f8f9fa;padding:8px;border-radius:0 0 8px 8px\">")
                for p in reactome_pathways:
                    name = p.get('name', '')
                    url = p.get('url', '')
                    lines.append(f"<div style=\"padding:8px 12px;margin:4px 0;background:white;border-radius:6px;border-left:3px solid #7b1fa2\">")
                    lines.append(f"<a href=\"{url}\" target=\"_blank\" style=\"color:#333;text-decoration:none\">{name}</a>")
                    lines.append("</div>")
                lines.append("</div>")
                lines.append("</div>")
                lines.append("")

            lines.append("</div>")
            lines.append("")

        # 6.2 文献来源
        if merged["pathways"]:
            lines.append("### 文献来源")
            lines.append("")
            sorted_pathways = sorted(merged["pathways"].items(), key=lambda x: len(x[1]), reverse=True)

            # 使用小卡片展示文献报道的通路
            lines.append("<div style=\"display:flex;flex-wrap:wrap;gap:8px\">")
            for pathway, pmids in sorted_pathways:
                pmid_links = ", ".join([f"<a href='https://pubmed.ncbi.nlm.nih.gov/{p}/'>{p}</a>" for p in pmids])
                lines.append(f"<div style=\"padding:8px 14px;background:#fff8e1;border:1px solid #ffb300;border-radius:16px;font-size:0.9em\">")
                lines.append(f"<strong>{pathway}</strong> <span style=\"color:#666\">({len(pmids)}篇)</span>")
                lines.append(f"<div style=\"font-size:0.8em;color:#888;margin-top:4px\">{pmid_links}</div>")
                lines.append("</div>")
            lines.append("</div>")
            lines.append("")

        if not has_db_pathways and not merged["pathways"]:
            lines.append("暂无信号通路相关信息")
            lines.append("")

        # =====================================================================
        # 7. 靶向治疗（如果有药物数据）
        # =====================================================================
        drugs = external_evidence.get("drugs", []) if external_evidence else []
        if drugs:
            approved_drugs = [d for d in drugs if d.get("is_approved")]
            clinical_drugs = [d for d in drugs if not d.get("is_approved") and d.get("phase", 0) >= 2]

            if approved_drugs or clinical_drugs:
                lines.append("## 靶向治疗")
                lines.append("")
                lines.append("> 数据来源: Open Targets")
                lines.append("")

                if approved_drugs:
                    lines.append("### 已批准药物")
                    lines.append("")
                    lines.append("| 药物 | 类型 | 作用机制 | 适应症 |")
                    lines.append("| --- | --- | --- | --- |")
                    for d in approved_drugs[:5]:
                        lines.append(f"| **{d.get('drug_name', '')}** | {d.get('drug_type', '')} | {d.get('mechanism', '')[:35]} | {d.get('indication', '')[:25]} |")
                    lines.append("")

                if clinical_drugs:
                    lines.append("### 临床试验药物")
                    lines.append("")
                    for d in clinical_drugs[:5]:
                        lines.append(f"- {d.get('drug_name', '')} (Phase {d.get('phase', '')}) - {d.get('indication', '')}")
                    lines.append("")

        # 冲突与争议 - 只在有冲突或需要显示一致性分析时才显示
        has_conflicts = conflict_report.conflicts and len(conflict_report.conflicts) > 0
        has_similarity_data = conflict_report.similarity_matrix and len(conflict_report.similarity_matrix) > 1

        if has_conflicts or has_similarity_data:
            lines.append("## 冲突与争议")
            lines.append("")

            # 显示分析方法和一致性评估
            if has_similarity_data:
                # 计算平均相似度
                n = len(conflict_report.similarity_matrix)
                total_sim = sum(
                    conflict_report.similarity_matrix[i][j]
                    for i in range(n) for j in range(i + 1, n)
                )
                avg_sim = total_sim / (n * (n - 1) / 2)
                consistency_level = "高" if avg_sim > 0.7 else "中等" if avg_sim > 0.5 else "较低"
                lines.append(f"> 📊 **语义分析结果**: 结论一致性{consistency_level}（平均相似度 {avg_sim:.1%}）")
                lines.append("")

            if has_conflicts:
                lines.append(f"⚠️ 检测到 **{len(conflict_report.conflicts)}** 处潜在冲突：")
                lines.append("")

                for i, conflict in enumerate(conflict_report.conflicts, 1):
                    # 冲突卡片标题
                    lines.append(f"### 冲突 {i}：{conflict.topic}")
                    lines.append("")

                    # 冲突维度标签
                    if conflict.conflict_dimension:
                        lines.append(f"<span style=\"background:#fff3cd;padding:2px 8px;border-radius:4px;font-size:0.85em\">冲突维度：{conflict.conflict_dimension}</span>")
                        lines.append("")

                    # 对比式布局
                    lines.append("<div style=\"display:flex;gap:16px;margin:16px 0\">")
                    lines.append("")

                    # 结论A（左侧）
                    lines.append("<div style=\"flex:1;padding:16px;background:#ffebee;border-left:4px solid #e53935;border-radius:8px\">")
                    lines.append(f"<div style=\"font-weight:600;color:#c62828;margin-bottom:8px\">结论 A</div>")
                    lines.append(f"<div style=\"font-size:0.85em;color:#666;margin-bottom:8px\">[PMID:{conflict.pmid_a}](https://pubmed.ncbi.nlm.nih.gov/{conflict.pmid_a}/)</div>")
                    lines.append(f"<div>{conflict.conclusion_a}</div>")
                    lines.append("</div>")
                    lines.append("")

                    # VS分隔符
                    lines.append("<div style=\"display:flex;align-items:center;font-weight:bold;color:#757575\">VS</div>")

                    # 结论B（右侧）
                    lines.append("<div style=\"flex:1;padding:16px;background:#e3f2fd;border-left:4px solid #1976d2;border-radius:8px\">")
                    lines.append(f"<div style=\"font-weight:600;color:#1565c0;margin-bottom:8px\">结论 B</div>")
                    lines.append(f"<div style=\"font-size:0.85em;color:#666;margin-bottom:8px\">[PMID:{conflict.pmid_b}](https://pubmed.ncbi.nlm.nih.gov/{conflict.pmid_b}/)</div>")
                    lines.append(f"<div>{conflict.conclusion_b}</div>")
                    lines.append("</div>")
                    lines.append("")

                    lines.append("</div>")
                    lines.append("")

                    # 解决建议
                    if conflict.resolution_suggestion:
                        lines.append(f"<div style=\"background:#f5f5f5;padding:12px;border-radius:8px;margin-top:8px\">")
                        lines.append(f"<strong>💡 建议：</strong>{conflict.resolution_suggestion}")
                        lines.append("</div>")
                        lines.append("")
            else:
                lines.append("✅ **经语义分析，本报告结论一致性较高，未检测到明显冲突**")
                lines.append("")
                # 添加可展开的原文对比区域
                lines.append("<details>")
                lines.append("<summary><strong>📋 点击展开查看所有结论详情</strong></summary>")
                lines.append("")
                # 使用卡片式布局展示每个结论
                for i, c in enumerate(conclusions, 1):
                    lines.append(f"<div style=\"padding:12px;margin:8px 0;background:#f8f9fa;border-left:3px solid #4caf50;border-radius:4px\">")

                    # PMID 和标题
                    title_short = c.title[:50] + "..." if c.title and len(c.title) > 50 else (c.title or "无标题")
                    lines.append(f"  <div style=\"margin-bottom:6px\">")
                    lines.append(f"    <a href=\"https://pubmed.ncbi.nlm.nih.gov/{c.pmid}/\" target=\"_blank\" style=\"font-weight:600;color:#1976d2\">PMID: {c.pmid}</a>")
                    lines.append(f"    <span style=\"color:#666;font-size:0.9em\"> - {title_short}</span>")
                    lines.append(f"  </div>")

                    # 疾病关系
                    disease = c.disease_relation[:80] + "..." if len(c.disease_relation) > 80 else c.disease_relation
                    lines.append(f"  <div style=\"font-size:0.9em;margin:4px 0\">")
                    lines.append(f"    <span style=\"color:#e53935;font-weight:500\">疾病关系:</span> {disease}")
                    lines.append(f"  </div>")

                    # 临床意义
                    clinical = c.clinical_significance[:60] + "..." if len(c.clinical_significance) > 60 else c.clinical_significance
                    if clinical and clinical not in INVALID_VALUES:
                        lines.append(f"  <div style=\"font-size:0.9em;margin:4px 0\">")
                        lines.append(f"    <span style=\"color:#1976d2;font-weight:500\">临床意义:</span> {clinical}")
                        lines.append(f"  </div>")

                    # 突变效应
                    mutation = c.mutation_effects[:60] + "..." if len(c.mutation_effects) > 60 else c.mutation_effects
                    if mutation and mutation not in INVALID_VALUES:
                        lines.append(f"  <div style=\"font-size:0.9em;margin:4px 0\">")
                        lines.append(f"    <span style=\"color:#7b1fa2;font-weight:500\">突变效应:</span> {mutation}")
                        lines.append(f"  </div>")

                    lines.append(f"</div>")
                lines.append("")
                lines.append("</details>")
                lines.append("")

        # 参考文献
        lines.append("## 参考文献")
        lines.append("")
        lines.append("> 按影响力评分排序，综合考虑年均引用次数、期刊影响因子和发表年份")
        lines.append("")

        # 导入归一化评分函数
        from .citation_fetcher import get_normalized_score

        if conclusions:
            # 按影响力评分排序（综合年均引用和IF）
            sorted_conclusions = sorted(
                conclusions,
                key=lambda c: get_normalized_score(c.citation_count, c.impact_factor, c.year),
                reverse=True
            )

            # 分组：高影响力 vs 普通
            high_impact = []
            normal = []
            for c in sorted_conclusions:
                score = get_normalized_score(c.citation_count, c.impact_factor, c.year)
                if score >= 50:
                    high_impact.append(c)
                else:
                    normal.append(c)

            # 渲染高影响力文献
            if high_impact:
                lines.append("### 🏆 高影响力文献")
                lines.append("")
                for c in high_impact:
                    lines.extend(self._format_reference_card(c))
                lines.append("")

            # 渲染普通文献
            if normal:
                lines.append("### 📚 其他参考文献")
                lines.append("")
                for c in normal:
                    lines.extend(self._format_reference_card(c, compact=True))
                lines.append("")
        else:
            lines.append("暂无参考文献")
        lines.append("")

        # 页脚
        lines.append("---")
        lines.append("")
        lines.append("*本报告由基因文献智能检索助手自动生成*")

        return "\n".join(lines)

    def save_markdown(self, gene: str, content: str, filename: Optional[str] = None) -> str:
        """
        保存 Markdown 文件

        Args:
            gene: 基因名
            content: Markdown 内容
            filename: 文件名（可选）

        Returns:
            保存的文件路径
        """
        if filename is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = f"{gene}_report_{timestamp}.md"

        filepath = self.output_dir / filename
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(content)

        return str(filepath)

    def generate_excel(
        self,
        gene: str,
        conclusions: List[GeneConclusion],
        conflict_report: ConflictReport
    ) -> str:
        """
        生成 Excel 报告

        Args:
            gene: 基因名
            conclusions: 结论列表
            conflict_report: 冲突报告

        Returns:
            Excel 文件路径
        """
        try:
            import pandas as pd
        except ImportError:
            raise ImportError("请安装 pandas: pip install pandas openpyxl")

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"{gene}_report_{timestamp}.xlsx"
        filepath = self.output_dir / filename

        with pd.ExcelWriter(filepath, engine='openpyxl') as writer:
            # 结论表
            if conclusions:
                df_conclusions = pd.DataFrame([c.to_dict() for c in conclusions])
                df_conclusions.to_excel(writer, sheet_name='结论详情', index=False)
            else:
                pd.DataFrame({'提示': ['暂无数据']}).to_excel(writer, sheet_name='结论详情', index=False)

            # 冲突表
            if conflict_report.conflicts:
                df_conflicts = pd.DataFrame([
                    {
                        '冲突主题': c.topic,
                        '文献A PMID': c.pmid_a,
                        '结论A': c.conclusion_a,
                        '文献B PMID': c.pmid_b,
                        '结论B': c.conclusion_b,
                        '严重程度': c.severity,
                        '建议': c.resolution_suggestion
                    }
                    for c in conflict_report.conflicts
                ])
                df_conflicts.to_excel(writer, sheet_name='冲突检测', index=False)
            else:
                pd.DataFrame({'提示': ['未检测到冲突']}).to_excel(writer, sheet_name='冲突检测', index=False)

            # 汇总表
            merged = self.merger.merge_conclusions(conclusions, conflict_report)
            summary_data = {
                '指标': ['总文献数', '冲突数量', '冲突比例', '需要审核'],
                '值': [
                    merged['total_papers'],
                    len(conflict_report.conflicts),
                    f"{conflict_report.conflict_ratio:.1%}",
                    '是' if conflict_report.needs_review else '否'
                ]
            }
            pd.DataFrame(summary_data).to_excel(writer, sheet_name='汇总', index=False)

        return str(filepath)

    def generate_json(
        self,
        gene: str,
        conclusions: List[GeneConclusion],
        conflict_report: ConflictReport
    ) -> str:
        """
        生成 JSON 报告

        Args:
            gene: 基因名
            conclusions: 结论列表
            conflict_report: 冲突报告

        Returns:
            JSON 文件路径
        """
        merged = self.merger.merge_conclusions(conclusions, conflict_report)

        report = {
            "gene": gene,
            "generated_at": datetime.now().isoformat(),
            "summary": merged,
            "conclusions": [c.to_dict() for c in conclusions],
            "conflicts": [
                {
                    "topic": c.topic,
                    "pmid_a": c.pmid_a,
                    "conclusion_a": c.conclusion_a,
                    "pmid_b": c.pmid_b,
                    "conclusion_b": c.conclusion_b,
                    "severity": c.severity,
                    "resolution_suggestion": c.resolution_suggestion
                }
                for c in conflict_report.conflicts
            ]
        }

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"{gene}_report_{timestamp}.json"
        filepath = self.output_dir / filename

        with open(filepath, 'w', encoding='utf-8') as f:
            json.dump(report, f, ensure_ascii=False, indent=2)

        return str(filepath)

    def generate_pdf(self, markdown_content: str, gene: str) -> Optional[str]:
        """
        从 Markdown 生成 PDF

        Args:
            markdown_content: Markdown 内容
            gene: 基因名

        Returns:
            PDF 文件路径，如果失败返回 None
        """
        try:
            import markdown
            import pdfkit
        except ImportError:
            print("PDF 生成需要安装: pip install markdown pdfkit")
            return None

        # 转换 Markdown 为 HTML
        html_content = markdown.markdown(
            markdown_content,
            extensions=['tables', 'toc']
        )

        # 添加样式
        styled_html = f"""
        <html>
        <head>
            <meta charset="utf-8">
            <style>
                body {{
                    font-family: "Microsoft YaHei", "SimSun", sans-serif;
                    line-height: 1.6;
                    padding: 20px;
                    max-width: 800px;
                    margin: 0 auto;
                }}
                h1 {{ color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }}
                h2 {{ color: #34495e; margin-top: 30px; }}
                h3 {{ color: #7f8c8d; }}
                table {{
                    border-collapse: collapse;
                    width: 100%;
                    margin: 20px 0;
                }}
                th, td {{
                    border: 1px solid #bdc3c7;
                    padding: 8px 12px;
                    text-align: left;
                }}
                th {{ background-color: #ecf0f1; }}
                blockquote {{
                    border-left: 4px solid #3498db;
                    padding-left: 15px;
                    color: #7f8c8d;
                    margin: 20px 0;
                }}
                code {{
                    background-color: #f4f4f4;
                    padding: 2px 6px;
                    border-radius: 3px;
                }}
                a {{ color: #3498db; }}
            </style>
        </head>
        <body>
            {html_content}
        </body>
        </html>
        """

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"{gene}_report_{timestamp}.pdf"
        filepath = self.output_dir / filename

        try:
            pdfkit.from_string(styled_html, str(filepath))
            return str(filepath)
        except Exception as e:
            print(f"PDF 生成失败: {e}")
            return None

    def generate_html(
        self,
        gene: str,
        conclusions: List[GeneConclusion],
        conflict_report: ConflictReport,
        include_charts: bool = True
    ) -> str:
        """
        生成包含图表的完整 HTML 报告

        Args:
            gene: 基因名
            conclusions: 结论列表
            conflict_report: 冲突报告
            include_charts: 是否包含图表

        Returns:
            HTML 文件路径
        """
        try:
            from ..utils.visualizer import ReportVisualizer
            visualizer = ReportVisualizer()
        except ImportError:
            logger.warning("无法导入可视化模块，将生成不含图表的HTML报告")
            visualizer = None

        # 生成 Markdown 内容
        md_content = self.generate_markdown(gene, conclusions, conflict_report)

        if visualizer and include_charts:
            # 使用可视化器生成完整 HTML
            html_content = visualizer.generate_html_report(
                gene, conclusions, conflict_report, md_content
            )
        else:
            # 简单 HTML
            import markdown
            html_content = markdown.markdown(md_content, extensions=['tables', 'toc'])
            html_content = f"""
<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <title>基因文献分析报告: {gene}</title>
    <style>
        body {{ font-family: sans-serif; line-height: 1.6; max-width: 800px; margin: 0 auto; padding: 20px; }}
        h1 {{ color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }}
        h2 {{ color: #34495e; }}
        table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #3498db; color: white; }}
        blockquote {{ border-left: 4px solid #3498db; padding-left: 15px; color: #7f8c8d; }}
        a {{ color: #3498db; }}
    </style>
</head>
<body>
    {html_content}
    <footer style="text-align: center; color: #7f8c8d; margin-top: 40px;">
        <p>本报告由基因文献智能检索助手自动生成</p>
    </footer>
</body>
</html>
"""

        # 保存文件
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"{gene}_report_{timestamp}.html"
        filepath = self.output_dir / filename

        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(html_content)

        return str(filepath)


def generate_report(
    gene: str,
    conclusions: List[GeneConclusion],
    conflict_report: ConflictReport,
    formats: List[str] = ["markdown", "json", "excel", "html"],
    output_dir: Optional[str] = None,
    include_charts: bool = True,
    enrich_citations: bool = True
) -> Dict[str, str]:
    """
    生成报告的便捷函数

    Args:
        gene: 基因名
        conclusions: 结论列表
        conflict_report: 冲突报告
        formats: 输出格式列表 (markdown, json, excel, pdf, html)
        output_dir: 输出目录
        include_charts: 是否在 HTML 中包含图表
        enrich_citations: 是否获取引用数据

    Returns:
        生成的文件路径字典
    """
    generator = ReportGenerator(output_dir)
    results = {}

    # 获取引用数据
    if enrich_citations and conclusions:
        conclusions = generator.enrich_conclusions_with_citations(conclusions)

    if "markdown" in formats:
        md_content = generator.generate_markdown(gene, conclusions, conflict_report)
        results["markdown"] = generator.save_markdown(gene, md_content)

    if "json" in formats:
        results["json"] = generator.generate_json(gene, conclusions, conflict_report)

    if "excel" in formats:
        results["excel"] = generator.generate_excel(gene, conclusions, conflict_report)

    if "pdf" in formats:
        md_content = generator.generate_markdown(gene, conclusions, conflict_report)
        results["pdf"] = generator.generate_pdf(md_content, gene)

    if "html" in formats:
        results["html"] = generator.generate_html(
            gene, conclusions, conflict_report, include_charts=include_charts
        )

    return results


if __name__ == "__main__":
    # 测试
    from .conclusion_extractor import GeneConclusion
    from .conflict_detector import ConflictReport, Conflict

    conclusions = [
        GeneConclusion(
            gene="BRAF",
            pmid="12345678",
            title="BRAF mutations in lung cancer",
            disease_relation="促进肺癌进展",
            pathways=["MAPK", "ERK"],
            clinical_significance="与不良预后相关",
            mutation_effects="V600E突变激活激酶活性",
            key_findings=["发现1", "发现2"],
            confidence="high",
            year=2023
        )
    ]

    conflicts = [Conflict(
        topic="预后意义",
        conclusion_a="与不良预后相关",
        pmid_a="12345678",
        conclusion_b="与良好预后相关",
        pmid_b="23456789",
        severity="high"
    )]

    report = ConflictReport(
        gene="BRAF",
        conflicts=conflicts,
        consensus=["BRAF突变在肺癌中常见"],
        needs_review=True,
        total_conclusions=2,
        conflict_ratio=0.5
    )

    generator = ReportGenerator()
    md = generator.generate_markdown("BRAF", conclusions, report)
    print(md)