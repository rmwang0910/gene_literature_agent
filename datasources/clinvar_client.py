"""
ClinVar API 客户端

提供变异临床意义查询功能。
基于 NCBI E-utilities API: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/
"""

import requests
import logging
import time
import xml.etree.ElementTree as ET
from typing import Optional, List, Dict, Any
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)


@dataclass
class VariantClinicalSignificance:
    """变异临床意义数据"""
    variation_id: str
    variant_name: str  # 如 "NM_004333.6:c.1799T>A (p.Val600Glu)"
    gene: str
    clinical_significance: str  # Pathogenic, Likely pathogenic, VUS, etc.
    review_status: str  # 星级评审状态
    star_rating: int  # 0-4 星
    condition: str  # 相关疾病
    accession: str  # RCV/VCV accession
    last_evaluated: str = ""
    pmids: List[str] = field(default_factory=list)


@dataclass
class GeneVariantSummary:
    """基因变异汇总"""
    gene: str
    total_variants: int
    pathogenic_count: int
    likely_pathogenic_count: int
    vus_count: int
    benign_count: int
    likely_benign_count: int
    top_conditions: List[str] = field(default_factory=list)


class ClinVarClient:
    """ClinVar E-utilities API 客户端"""

    BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    # 临床意义映射
    CLINICAL_SIGNIFICANCE_MAP = {
        "pathogenic": "致病",
        "likely pathogenic": "可能致病",
        "uncertain significance": "意义不明",
        "likely benign": "可能良性",
        "benign": "良性",
        "conflicting interpretations of pathogenicity": "意义冲突",
        "drug response": "药物响应",
        "association": "关联",
        "risk factor": "风险因子",
        "protective": "保护性",
        "affects": "影响",
        "not provided": "未提供",
    }

    # 评审状态星级
    REVIEW_STATUS_STARS = {
        "practice guideline": 4,
        "reviewed by expert panel": 3,
        "criteria provided, multiple submitters, no conflicts": 2,
        "criteria provided, conflicting interpretations": 1,
        "criteria provided, single submitter": 1,
        "no assertion criteria provided": 0,
        "no assertion provided": 0,
    }

    def __init__(self, api_key: Optional[str] = None, email: Optional[str] = None, timeout: int = 30):
        """
        初始化 ClinVar 客户端

        Args:
            api_key: NCBI API key（可选，有 key 速率限制 10 req/s，无 key 3 req/s）
            email: 联系邮箱（推荐）
            timeout: 请求超时时间
        """
        self.api_key = api_key
        self.email = email
        self.timeout = timeout
        self._session = requests.Session()
        self._last_request_time = 0
        self._min_interval = 0.34 if api_key else 1.0  # 遵守速率限制

    def _rate_limit(self):
        """速率限制"""
        elapsed = time.time() - self._last_request_time
        if elapsed < self._min_interval:
            time.sleep(self._min_interval - elapsed)
        self._last_request_time = time.time()

    def _build_params(self, **kwargs) -> Dict[str, str]:
        """构建请求参数"""
        params = {k: v for k, v in kwargs.items() if v is not None}
        if self.api_key:
            params["api_key"] = self.api_key
        if self.email:
            params["email"] = self.email
        return params

    def search_variants(
        self,
        gene: str,
        clinical_significance: Optional[str] = None,
        max_results: int = 100
    ) -> List[str]:
        """
        搜索基因的变异

        Args:
            gene: 基因符号
            clinical_significance: 过滤临床意义（pathogenic, likely pathogenic, etc.）
            max_results: 最大返回数量

        Returns:
            ClinVar Variation ID 列表
        """
        self._rate_limit()

        # 构建查询
        query_parts = [f"{gene}[gene]"]
        if clinical_significance:
            query_parts.append(f"{clinical_significance}[CLNSIG]")

        query = " AND ".join(query_parts)

        params = self._build_params(
            db="clinvar",
            term=query,
            retmax=max_results,
            retmode="json"
        )

        try:
            response = self._session.get(
                f"{self.BASE_URL}/esearch.fcgi",
                params=params,
                timeout=self.timeout
            )
            response.raise_for_status()
            data = response.json()

            id_list = data.get("esearchresult", {}).get("idlist", [])
            logger.info(f"ClinVar 搜索 {gene}: 找到 {len(id_list)} 个变异")
            return id_list

        except requests.RequestException as e:
            logger.error(f"ClinVar 搜索失败: {e}")
            return []

    def get_variant_summary(self, variant_ids: List[str]) -> List[Dict[str, Any]]:
        """
        获取变异摘要信息

        Args:
            variant_ids: ClinVar Variation ID 列表

        Returns:
            变异摘要列表
        """
        if not variant_ids:
            return []

        self._rate_limit()

        params = self._build_params(
            db="clinvar",
            id=",".join(variant_ids[:100]),  # 限制每次最多 100 个
            retmode="json"
        )

        try:
            response = self._session.get(
                f"{self.BASE_URL}/esummary.fcgi",
                params=params,
                timeout=self.timeout
            )
            response.raise_for_status()
            data = response.json()

            results = []
            result_data = data.get("result", {})

            for vid in variant_ids:
                if vid in result_data:
                    entry = result_data[vid]

                    # 临床意义可能在 germline_classification 或 clinical_significance 下
                    germline = entry.get("germline_classification", {})
                    clinical = entry.get("clinical_significance", {})

                    # 优先使用 germline_classification，回退到 clinical_significance
                    significance = germline.get("description", "") or clinical.get("description", "")
                    review_status = germline.get("review_status", "") or clinical.get("review_status", "")

                    results.append({
                        "variation_id": vid,
                        "title": entry.get("title", ""),
                        "gene": self._extract_gene(entry),
                        "clinical_significance": significance,
                        "review_status": review_status,
                        "protein_change": entry.get("protein_change", ""),
                        "accession": entry.get("accession", ""),
                    })

            return results

        except requests.RequestException as e:
            logger.error(f"获取 ClinVar 变异摘要失败: {e}")
            return []

    def _extract_gene(self, entry: Dict) -> str:
        """从 esummary 结果中提取基因名"""
        genes = entry.get("genes", [])
        if genes and isinstance(genes, list):
            return genes[0].get("symbol", "") if genes[0] else ""
        return ""

    def get_pathogenic_variants(
        self,
        gene: str,
        include_likely: bool = True,
        max_results: int = 50
    ) -> List[VariantClinicalSignificance]:
        """
        获取基因的致病变异

        Args:
            gene: 基因符号
            include_likely: 是否包含可能致病
            max_results: 最大结果数

        Returns:
            致病变异列表
        """
        # 搜索致病变异
        pathogenic_ids = self.search_variants(gene, "pathogenic", max_results)

        if include_likely:
            likely_ids = self.search_variants(gene, "likely pathogenic", max_results)
            pathogenic_ids = list(set(pathogenic_ids + likely_ids))[:max_results]

        if not pathogenic_ids:
            return []

        # 获取详细信息
        summaries = self.get_variant_summary(pathogenic_ids)

        results = []
        for s in summaries:
            review_status = s.get("review_status", "").lower()
            star_rating = 0
            for status, stars in self.REVIEW_STATUS_STARS.items():
                if status in review_status:
                    star_rating = stars
                    break

            results.append(VariantClinicalSignificance(
                variation_id=s.get("variation_id", ""),
                variant_name=s.get("title", ""),
                gene=s.get("gene", gene),
                clinical_significance=s.get("clinical_significance", ""),
                review_status=s.get("review_status", ""),
                star_rating=star_rating,
                condition="",  # 需要额外查询获取
                accession=s.get("accession", ""),
            ))

        # 按星级排序
        results.sort(key=lambda x: x.star_rating, reverse=True)
        return results

    def get_gene_variant_summary(self, gene: str) -> GeneVariantSummary:
        """
        获取基因变异统计摘要

        Args:
            gene: 基因符号

        Returns:
            基因变异统计
        """
        # 分别查询各类意义的变异数量
        pathogenic = len(self.search_variants(gene, "pathogenic", 1000))
        likely_pathogenic = len(self.search_variants(gene, "likely pathogenic", 1000))
        vus = len(self.search_variants(gene, "uncertain significance", 1000))
        benign = len(self.search_variants(gene, "benign", 1000))
        likely_benign = len(self.search_variants(gene, "likely benign", 1000))

        total = pathogenic + likely_pathogenic + vus + benign + likely_benign

        return GeneVariantSummary(
            gene=gene,
            total_variants=total,
            pathogenic_count=pathogenic,
            likely_pathogenic_count=likely_pathogenic,
            vus_count=vus,
            benign_count=benign,
            likely_benign_count=likely_benign
        )

    def get_variant_clinical_significance(
        self,
        gene: str,
        variant: str
    ) -> Optional[VariantClinicalSignificance]:
        """
        查询特定变异的临床意义

        Args:
            gene: 基因符号
            variant: 变异描述（如 "V600E", "c.1799T>A"）

        Returns:
            变异临床意义，未找到返回 None
        """
        self._rate_limit()

        # 构建查询：基因 + 变异
        query = f"{gene}[gene] AND {variant}[variant name]"

        params = self._build_params(
            db="clinvar",
            term=query,
            retmax=10,
            retmode="json"
        )

        try:
            response = self._session.get(
                f"{self.BASE_URL}/esearch.fcgi",
                params=params,
                timeout=self.timeout
            )
            response.raise_for_status()
            data = response.json()

            id_list = data.get("esearchresult", {}).get("idlist", [])

            if not id_list:
                # 尝试更宽松的搜索
                query = f"{gene}[gene] AND {variant}"
                params["term"] = query
                response = self._session.get(
                    f"{self.BASE_URL}/esearch.fcgi",
                    params=params,
                    timeout=self.timeout
                )
                data = response.json()
                id_list = data.get("esearchresult", {}).get("idlist", [])

            if not id_list:
                return None

            summaries = self.get_variant_summary(id_list[:1])
            if not summaries:
                return None

            s = summaries[0]
            review_status = s.get("review_status", "").lower()
            star_rating = 0
            for status, stars in self.REVIEW_STATUS_STARS.items():
                if status in review_status:
                    star_rating = stars
                    break

            return VariantClinicalSignificance(
                variation_id=s.get("variation_id", ""),
                variant_name=s.get("title", ""),
                gene=s.get("gene", gene),
                clinical_significance=s.get("clinical_significance", ""),
                review_status=s.get("review_status", ""),
                star_rating=star_rating,
                condition="",
                accession=s.get("accession", ""),
            )

        except requests.RequestException as e:
            logger.error(f"查询 ClinVar 变异失败: {e}")
            return None

    def translate_significance(self, significance: str) -> str:
        """将临床意义翻译为中文"""
        sig_lower = significance.lower().strip()
        return self.CLINICAL_SIGNIFICANCE_MAP.get(sig_lower, significance)


# 便捷函数
def get_variant_significance(gene: str, variant: str) -> Optional[Dict[str, Any]]:
    """
    获取变异临床意义（便捷函数）

    Args:
        gene: 基因符号
        variant: 变异描述

    Returns:
        {"significance": "...", "significance_cn": "...", "stars": int, ...}
    """
    client = ClinVarClient()
    result = client.get_variant_clinical_significance(gene, variant)

    if not result:
        return None

    return {
        "gene": result.gene,
        "variant": result.variant_name,
        "significance": result.clinical_significance,
        "significance_cn": client.translate_significance(result.clinical_significance),
        "stars": result.star_rating,
        "review_status": result.review_status,
        "accession": result.accession,
        "source": "ClinVar",
        "source_url": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{result.variation_id}/"
    }


if __name__ == "__main__":
    # 测试
    client = ClinVarClient()

    print("=== 测试 BRAF V600E 临床意义 ===")
    result = client.get_variant_clinical_significance("BRAF", "V600E")
    if result:
        print(f"变异: {result.variant_name}")
        print(f"临床意义: {result.clinical_significance}")
        print(f"中文: {client.translate_significance(result.clinical_significance)}")
        print(f"星级: {'★' * result.star_rating}{'☆' * (4 - result.star_rating)}")
        print(f"评审状态: {result.review_status}")

    print("\n=== 测试 TP53 致病变异统计 ===")
    summary = client.get_gene_variant_summary("TP53")
    print(f"基因: {summary.gene}")
    print(f"总变异数: {summary.total_variants}")
    print(f"  致病: {summary.pathogenic_count}")
    print(f"  可能致病: {summary.likely_pathogenic_count}")
    print(f"  VUS: {summary.vus_count}")
    print(f"  良性: {summary.benign_count}")
    print(f"  可能良性: {summary.likely_benign_count}")

    print("\n=== 测试 BRCA1 致病变异列表 ===")
    variants = client.get_pathogenic_variants("BRCA1", max_results=5)
    for v in variants[:5]:
        print(f"  - {v.variant_name}")
        print(f"    意义: {v.clinical_significance} ({'★' * v.star_rating})")
