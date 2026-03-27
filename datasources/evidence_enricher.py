"""
证据整合模块

整合多个外部数据源的证据，为文献结论提供权威支撑：
- Open Targets: 基因-疾病关联
- ClinVar: 变异临床意义
- COSMIC: 热点突变、癌症基因
- KEGG/Reactome: 通路数据
"""

import logging
from typing import List, Dict, Any, Optional
from dataclasses import dataclass, field

from .opentargets_client import OpenTargetsClient, DiseaseAssociation
from .clinvar_client import ClinVarClient, VariantClinicalSignificance
from .cosmic_client import COSMICClient, CancerMutation
from .pathway_client import PathwayClient, Pathway

logger = logging.getLogger(__name__)


@dataclass
class EnrichedEvidence:
    """整合后的证据数据"""
    gene: str

    # Open Targets 数据
    disease_associations: List[Dict[str, Any]] = field(default_factory=list)
    drugs: List[Dict[str, Any]] = field(default_factory=list)
    opentargets_pmids: List[str] = field(default_factory=list)

    # ClinVar 数据
    pathogenic_variants: List[Dict[str, Any]] = field(default_factory=list)
    variant_significance: Optional[Dict[str, Any]] = None

    # COSMIC 数据
    is_cancer_gene: bool = False
    cancer_gene_role: str = ""
    hotspot_mutations: List[Dict[str, Any]] = field(default_factory=list)
    related_cancer_types: List[str] = field(default_factory=list)

    # Pathway 数据
    pathways: List[Dict[str, Any]] = field(default_factory=list)

    # 元数据
    sources_queried: List[str] = field(default_factory=list)


class EvidenceEnricher:
    """证据整合器"""

    def __init__(
        self,
        cache_dir: Optional[str] = None,
        clinvar_api_key: Optional[str] = None,
        cosmic_data_dir: Optional[str] = None
    ):
        """
        初始化证据整合器

        Args:
            cache_dir: 缓存目录
            clinvar_api_key: ClinVar API key
            cosmic_data_dir: COSMIC 本地数据目录
        """
        self.opentargets = OpenTargetsClient(cache_dir=cache_dir)
        self.clinvar = ClinVarClient(api_key=clinvar_api_key)
        self.cosmic = COSMICClient(data_dir=cosmic_data_dir)
        self.pathway = PathwayClient()

    def enrich_gene(
        self,
        gene: str,
        variant: Optional[str] = None,
        include_sources: Optional[List[str]] = None
    ) -> EnrichedEvidence:
        """
        整合基因的多源证据

        Args:
            gene: 基因符号
            variant: 可选的变异描述（如 "V600E"）
            include_sources: 要查询的数据源列表，默认全部

        Returns:
            整合后的证据
        """
        if include_sources is None:
            include_sources = ["OpenTargets", "ClinVar", "COSMIC", "Pathway"]

        evidence = EnrichedEvidence(gene=gene)

        # 1. Open Targets 数据
        if "OpenTargets" in include_sources:
            try:
                self._enrich_opentargets(evidence, gene)
                evidence.sources_queried.append("OpenTargets")
            except Exception as e:
                logger.warning(f"Open Targets 查询失败: {e}")

        # 2. ClinVar 数据
        if "ClinVar" in include_sources:
            try:
                self._enrich_clinvar(evidence, gene, variant)
                evidence.sources_queried.append("ClinVar")
            except Exception as e:
                logger.warning(f"ClinVar 查询失败: {e}")

        # 3. COSMIC 数据
        if "COSMIC" in include_sources:
            try:
                self._enrich_cosmic(evidence, gene, variant)
                evidence.sources_queried.append("COSMIC")
            except Exception as e:
                logger.warning(f"COSMIC 查询失败: {e}")

        # 4. Pathway 数据
        if "Pathway" in include_sources:
            try:
                self._enrich_pathway(evidence, gene)
                evidence.sources_queried.append("Pathway")
            except Exception as e:
                logger.warning(f"Pathway 查询失败: {e}")

        return evidence

    def _enrich_opentargets(self, evidence: EnrichedEvidence, gene: str):
        """从 Open Targets 获取数据"""
        # 搜索基因
        target = self.opentargets.search_target(gene)
        if not target:
            return

        # 获取疾病关联
        associations = self.opentargets.get_disease_associations(
            target.ensembl_id, min_score=0.1, size=20
        )
        evidence.disease_associations = [
            {
                "disease_id": a.disease_id,
                "disease_name": a.disease_name,
                "score": round(a.score, 3),
                "evidence_types": a.evidence_types,
                "source": "OpenTargets",
                "url": f"https://platform.opentargets.org/evidence/{target.ensembl_id}/{a.disease_id}"
            }
            for a in associations
        ]

        # 获取药物
        drugs = self.opentargets.get_associated_drugs(target.ensembl_id, size=20)
        evidence.drugs = [
            {
                "drug_id": d.drug_id,
                "drug_name": d.drug_name,
                "drug_type": d.drug_type,
                "mechanism": d.mechanism_of_action,
                "phase": d.phase,
                "is_approved": d.is_approved,
                "indication": d.indication,
                "source": "OpenTargets"
            }
            for d in drugs
        ]

        # 获取 PMIDs（从证据中提取）
        if associations:
            top_disease = associations[0]
            evidences = self.opentargets.get_target_disease_evidence(
                target.ensembl_id, top_disease.disease_id, size=50
            )
            pmids = set()
            for ev in evidences:
                pmids.update(ev.get("pmids", []))
            evidence.opentargets_pmids = list(pmids)[:20]

    def _enrich_clinvar(self, evidence: EnrichedEvidence, gene: str, variant: Optional[str]):
        """从 ClinVar 获取数据"""
        # 获取致病变异列表
        pathogenic = self.clinvar.get_pathogenic_variants(gene, max_results=10)
        evidence.pathogenic_variants = [
            {
                "variant": v.variant_name,
                "significance": v.clinical_significance,
                "significance_cn": self.clinvar.translate_significance(v.clinical_significance),
                "stars": v.star_rating,
                "review_status": v.review_status,
                "accession": v.accession,
                "source": "ClinVar",
                "url": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{v.variation_id}/"
            }
            for v in pathogenic
        ]

        # 如果指定了变异，查询其临床意义
        if variant:
            var_sig = self.clinvar.get_variant_clinical_significance(gene, variant)
            if var_sig:
                evidence.variant_significance = {
                    "variant": var_sig.variant_name,
                    "significance": var_sig.clinical_significance,
                    "significance_cn": self.clinvar.translate_significance(var_sig.clinical_significance),
                    "stars": var_sig.star_rating,
                    "review_status": var_sig.review_status,
                    "source": "ClinVar",
                    "url": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{var_sig.variation_id}/"
                }

    def _enrich_cosmic(self, evidence: EnrichedEvidence, gene: str, variant: Optional[str]):
        """从 COSMIC 获取数据"""
        # 检查是否为癌症基因
        cancer_gene = self.cosmic.get_cancer_gene_info(gene)
        evidence.is_cancer_gene = cancer_gene is not None
        evidence.cancer_gene_role = cancer_gene.role if cancer_gene else ""

        # 获取热点突变
        hotspots = self.cosmic.get_hotspot_mutations(gene)
        evidence.hotspot_mutations = [
            {
                "mutation": m.mutation,
                "sample_count": m.sample_count,
                "frequency": round(m.frequency, 4),
                "cancer_types": m.cancer_types[:5],
                "is_hotspot": True,
                "source": "COSMIC"
            }
            for m in hotspots[:10]
        ]

        # 如果指定了变异，标记是否为热点
        if variant:
            is_hotspot = self.cosmic.is_hotspot(gene, variant)
            mutation_info = self.cosmic.get_mutation_info(gene, variant)
            if mutation_info:
                evidence.hotspot_mutations.insert(0, {
                    "mutation": variant,
                    "sample_count": mutation_info.sample_count,
                    "frequency": round(mutation_info.frequency, 4),
                    "cancer_types": mutation_info.cancer_types[:5],
                    "is_hotspot": is_hotspot,
                    "is_query_variant": True,
                    "source": "COSMIC"
                })

        # 获取相关癌症类型
        evidence.related_cancer_types = self.cosmic.get_gene_cancer_types(gene)[:10]

    def _enrich_pathway(self, evidence: EnrichedEvidence, gene: str):
        """从通路数据库获取数据"""
        pathways = self.pathway.get_gene_pathways(gene)

        # 分别获取 KEGG 和 Reactome 通路，确保两者都能展示
        kegg_pathways = [p for p in pathways if p.source == "KEGG"][:10]
        reactome_pathways = [p for p in pathways if p.source == "Reactome"][:10]

        # 合并结果
        all_pathways = kegg_pathways + reactome_pathways
        evidence.pathways = [
            {
                "id": p.id,
                "name": p.name,
                "source": p.source,
                "url": p.url
            }
            for p in all_pathways
        ]

    def get_disease_association_with_pmids(
        self,
        gene: str,
        disease_name: str
    ) -> Dict[str, Any]:
        """
        获取特定基因-疾病关联的详细证据（含 PMIDs）

        Args:
            gene: 基因符号
            disease_name: 疾病名称

        Returns:
            关联证据
        """
        # 搜索基因
        target = self.opentargets.search_target(gene)
        if not target:
            return {"error": f"未找到基因 {gene}"}

        # 搜索疾病
        disease = self.opentargets.search_disease(disease_name)
        if not disease:
            return {"error": f"未找到疾病 {disease_name}"}

        # 获取详细证据
        evidences = self.opentargets.get_target_disease_evidence(
            target.ensembl_id, disease["id"], size=100
        )

        pmids = set()
        evidence_by_type = {}

        for ev in evidences:
            dtype = ev.get("datatype", "unknown")
            if dtype not in evidence_by_type:
                evidence_by_type[dtype] = []
            evidence_by_type[dtype].append(ev)
            pmids.update(ev.get("pmids", []))

        return {
            "gene": gene,
            "disease": disease["name"],
            "disease_id": disease["id"],
            "total_evidence": len(evidences),
            "evidence_by_type": {
                k: len(v) for k, v in evidence_by_type.items()
            },
            "pmids": list(pmids),
            "source": "OpenTargets",
            "url": f"https://platform.opentargets.org/evidence/{target.ensembl_id}/{disease['id']}"
        }

    def get_clinical_significance_with_source(
        self,
        gene: str,
        variant: str
    ) -> Dict[str, Any]:
        """
        获取变异的临床意义（带数据来源）

        Args:
            gene: 基因符号
            variant: 变异描述

        Returns:
            临床意义信息
        """
        result = {
            "gene": gene,
            "variant": variant,
            "sources": []
        }

        # ClinVar
        clinvar_sig = self.clinvar.get_variant_clinical_significance(gene, variant)
        if clinvar_sig:
            result["sources"].append({
                "source": "ClinVar",
                "significance": clinvar_sig.clinical_significance,
                "significance_cn": self.clinvar.translate_significance(clinvar_sig.clinical_significance),
                "stars": clinvar_sig.star_rating,
                "review_status": clinvar_sig.review_status,
                "url": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{clinvar_sig.variation_id}/"
            })

        # COSMIC 热点检查
        is_hotspot = self.cosmic.is_hotspot(gene, variant)
        mutation_info = self.cosmic.get_mutation_info(gene, variant)
        if is_hotspot or mutation_info:
            cosmic_data = {
                "source": "COSMIC",
                "is_hotspot": is_hotspot,
            }
            if mutation_info:
                cosmic_data.update({
                    "sample_count": mutation_info.sample_count,
                    "cancer_types": mutation_info.cancer_types[:5]
                })
            result["sources"].append(cosmic_data)

        # 综合评估
        if result["sources"]:
            # 优先使用 ClinVar 的结论
            clinvar_source = next((s for s in result["sources"] if s["source"] == "ClinVar"), None)
            if clinvar_source:
                result["primary_significance"] = clinvar_source["significance"]
                result["primary_significance_cn"] = clinvar_source["significance_cn"]
            elif is_hotspot:
                result["primary_significance"] = "Likely oncogenic (hotspot)"
                result["primary_significance_cn"] = "可能致癌（热点突变）"

        return result


# 便捷函数
def enrich_gene_evidence(gene: str, variant: Optional[str] = None) -> Dict[str, Any]:
    """
    获取基因的多源证据（便捷函数）

    Args:
        gene: 基因符号
        variant: 可选的变异描述

    Returns:
        整合后的证据字典
    """
    enricher = EvidenceEnricher()
    evidence = enricher.enrich_gene(gene, variant)

    return {
        "gene": evidence.gene,
        "is_cancer_gene": evidence.is_cancer_gene,
        "cancer_gene_role": evidence.cancer_gene_role,
        "disease_associations": evidence.disease_associations[:5],
        "drugs": [d for d in evidence.drugs if d["is_approved"]][:5],
        "pathogenic_variants": evidence.pathogenic_variants[:5],
        "variant_significance": evidence.variant_significance,
        "hotspot_mutations": evidence.hotspot_mutations[:5],
        "pathways": evidence.pathways[:10],
        "sources_queried": evidence.sources_queried
    }


if __name__ == "__main__":
    # 测试
    enricher = EvidenceEnricher()

    print("=== 测试 BRAF V600E 证据整合 ===")
    evidence = enricher.enrich_gene("BRAF", "V600E")

    print(f"\n基因: {evidence.gene}")
    print(f"是否癌症基因: {evidence.is_cancer_gene}")
    print(f"角色: {evidence.cancer_gene_role}")
    print(f"查询的数据源: {', '.join(evidence.sources_queried)}")

    print(f"\n疾病关联 ({len(evidence.disease_associations)}):")
    for d in evidence.disease_associations[:3]:
        print(f"  - {d['disease_name']} (评分: {d['score']})")

    print(f"\n药物 ({len(evidence.drugs)}):")
    for d in evidence.drugs[:3]:
        status = "已批准" if d['is_approved'] else f"Phase {d['phase']}"
        print(f"  - {d['drug_name']} [{status}]")

    print(f"\n热点突变 ({len(evidence.hotspot_mutations)}):")
    for m in evidence.hotspot_mutations[:3]:
        print(f"  - {m['mutation']}: {m['sample_count']} 样本")

    if evidence.variant_significance:
        print(f"\nV600E 临床意义:")
        print(f"  意义: {evidence.variant_significance['significance']}")
        print(f"  中文: {evidence.variant_significance['significance_cn']}")
        print(f"  星级: {'★' * evidence.variant_significance['stars']}")

    print(f"\n通路 ({len(evidence.pathways)}):")
    for p in evidence.pathways[:5]:
        print(f"  - [{p['source']}] {p['name']}")
