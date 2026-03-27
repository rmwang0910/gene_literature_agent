"""
Open Targets API 客户端

提供基因-疾病关联、药物靶点、证据数据的查询功能。
基于 GraphQL API: https://api.platform.opentargets.org/api/v4/graphql
"""

import requests
import logging
from typing import Optional, List, Dict, Any
from dataclasses import dataclass, field
import json
from pathlib import Path

logger = logging.getLogger(__name__)


@dataclass
class DiseaseAssociation:
    """疾病关联数据"""
    disease_id: str  # EFO ID
    disease_name: str
    score: float  # 关联评分 0-1
    evidence_count: int = 0
    evidence_types: List[str] = field(default_factory=list)  # genetic_association, somatic_mutation, etc.
    literature_pmids: List[str] = field(default_factory=list)


@dataclass
class DrugAssociation:
    """药物关联数据"""
    drug_id: str  # ChEMBL ID
    drug_name: str
    drug_type: str  # Small molecule, Antibody, etc.
    mechanism_of_action: str
    phase: int  # 临床试验阶段
    is_approved: bool
    indication: str = ""


@dataclass
class TargetInfo:
    """靶点（基因）信息"""
    ensembl_id: str
    symbol: str
    name: str
    description: str = ""
    tractability: Dict[str, Any] = field(default_factory=dict)
    safety_liabilities: List[str] = field(default_factory=list)


class OpenTargetsClient:
    """Open Targets Platform API 客户端"""

    API_ENDPOINT = "https://api.platform.opentargets.org/api/v4/graphql"

    def __init__(self, cache_dir: Optional[str] = None, timeout: int = 30):
        self.timeout = timeout
        self.cache_dir = Path(cache_dir) if cache_dir else None
        if self.cache_dir:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
        self._session = requests.Session()
        self._session.headers.update({
            'Content-Type': 'application/json',
            'Accept': 'application/json'
        })

    def _execute_query(self, query: str, variables: Optional[Dict] = None) -> Dict:
        """执行 GraphQL 查询"""
        payload = {"query": query}
        if variables:
            payload["variables"] = variables

        try:
            response = self._session.post(
                self.API_ENDPOINT,
                json=payload,
                timeout=self.timeout
            )
            response.raise_for_status()
            result = response.json()

            if "errors" in result:
                logger.warning(f"GraphQL errors: {result['errors']}")

            return result.get("data", {})
        except requests.RequestException as e:
            logger.error(f"Open Targets API 请求失败: {e}")
            return {}

    def search_target(self, gene_symbol: str) -> Optional[TargetInfo]:
        """
        通过基因符号搜索靶点，获取 Ensembl ID

        Args:
            gene_symbol: 基因符号，如 "BRAF", "TP53"

        Returns:
            TargetInfo 对象
        """
        query = """
        query searchTarget($queryString: String!) {
            search(queryString: $queryString, entityNames: ["target"], page: {size: 5, index: 0}) {
                hits {
                    id
                    name
                    description
                    entity
                }
            }
        }
        """

        result = self._execute_query(query, {"queryString": gene_symbol})
        hits = result.get("search", {}).get("hits", [])

        for hit in hits:
            if hit.get("entity") == "target":
                # 精确匹配基因符号
                if hit.get("name", "").upper() == gene_symbol.upper():
                    return TargetInfo(
                        ensembl_id=hit.get("id", ""),
                        symbol=hit.get("name", ""),
                        name=hit.get("name", ""),
                        description=hit.get("description", "")
                    )

        # 如果没有精确匹配，返回第一个
        if hits:
            hit = hits[0]
            return TargetInfo(
                ensembl_id=hit.get("id", ""),
                symbol=hit.get("name", ""),
                name=hit.get("name", ""),
                description=hit.get("description", "")
            )

        return None

    def get_disease_associations(
        self,
        ensembl_id: str,
        min_score: float = 0.1,
        size: int = 50
    ) -> List[DiseaseAssociation]:
        """
        获取基因的疾病关联

        Args:
            ensembl_id: Ensembl 基因 ID
            min_score: 最小关联评分
            size: 返回数量

        Returns:
            疾病关联列表
        """
        query = """
        query targetDiseases($ensemblId: String!, $size: Int!) {
            target(ensemblId: $ensemblId) {
                id
                approvedSymbol
                associatedDiseases(page: {size: $size, index: 0}) {
                    count
                    rows {
                        disease {
                            id
                            name
                        }
                        score
                        datatypeScores {
                            id
                            score
                        }
                    }
                }
            }
        }
        """

        result = self._execute_query(query, {"ensemblId": ensembl_id, "size": size})
        target_data = result.get("target", {})

        if not target_data:
            return []

        associations = []
        rows = target_data.get("associatedDiseases", {}).get("rows", [])

        for row in rows:
            score = row.get("score", 0)
            if score < min_score:
                continue

            disease = row.get("disease", {})
            datatype_scores = row.get("datatypeScores", [])
            evidence_types = [dt["id"] for dt in datatype_scores if dt.get("score", 0) > 0]

            associations.append(DiseaseAssociation(
                disease_id=disease.get("id", ""),
                disease_name=disease.get("name", ""),
                score=score,
                evidence_types=evidence_types
            ))

        return associations

    def get_target_disease_evidence(
        self,
        ensembl_id: str,
        efo_id: str,
        size: int = 100
    ) -> List[Dict[str, Any]]:
        """
        获取特定基因-疾病对的详细证据

        Args:
            ensembl_id: Ensembl 基因 ID
            efo_id: EFO 疾病 ID
            size: 返回数量

        Returns:
            证据列表，包含来源、评分、文献等
        """
        query = """
        query targetDiseaseEvidence($ensemblId: String!, $efoId: String!, $size: Int!) {
            target(ensemblId: $ensemblId) {
                evidences(efoIds: [$efoId], size: $size) {
                    count
                    rows {
                        id
                        score
                        datasourceId
                        datatypeId
                        literature
                        diseaseFromSource
                        targetFromSourceId
                    }
                }
            }
        }
        """

        result = self._execute_query(query, {
            "ensemblId": ensembl_id,
            "efoId": efo_id,
            "size": size
        })

        evidences = result.get("target", {}).get("evidences", {}).get("rows", [])

        processed = []
        for ev in evidences:
            pmids = []
            if ev.get("literature"):
                pmids = [str(lit) for lit in ev["literature"]]

            processed.append({
                "id": ev.get("id", ""),
                "score": ev.get("score", 0),
                "datasource": ev.get("datasourceId", ""),
                "datatype": ev.get("datatypeId", ""),
                "pmids": pmids,
                "disease_from_source": ev.get("diseaseFromSource", "")
            })

        return processed

    def get_associated_drugs(
        self,
        ensembl_id: str,
        size: int = 50
    ) -> List[DrugAssociation]:
        """
        获取靶向该基因的药物

        Args:
            ensembl_id: Ensembl 基因 ID
            size: 返回数量

        Returns:
            药物关联列表
        """
        query = """
        query targetDrugs($ensemblId: String!, $size: Int!) {
            target(ensemblId: $ensemblId) {
                id
                approvedSymbol
                knownDrugs(size: $size) {
                    count
                    rows {
                        drug {
                            id
                            name
                            drugType
                            maximumClinicalTrialPhase
                            isApproved
                            hasBeenWithdrawn
                        }
                        mechanismOfAction
                        phase
                        disease {
                            id
                            name
                        }
                    }
                }
            }
        }
        """

        result = self._execute_query(query, {"ensemblId": ensembl_id, "size": size})
        rows = result.get("target", {}).get("knownDrugs", {}).get("rows", [])

        drugs = []
        seen_drugs = set()  # 去重

        for row in rows:
            drug = row.get("drug", {})
            drug_id = drug.get("id", "")

            if drug_id in seen_drugs:
                continue
            seen_drugs.add(drug_id)

            disease = row.get("disease", {})

            drugs.append(DrugAssociation(
                drug_id=drug_id,
                drug_name=drug.get("name", ""),
                drug_type=drug.get("drugType", ""),
                mechanism_of_action=row.get("mechanismOfAction", ""),
                phase=row.get("phase", 0),
                is_approved=drug.get("isApproved", False),
                indication=disease.get("name", "")
            ))

        return drugs

    def search_disease(self, disease_name: str) -> Optional[Dict[str, str]]:
        """
        搜索疾病获取 EFO ID

        Args:
            disease_name: 疾病名称

        Returns:
            {"id": "EFO_xxx", "name": "..."}
        """
        query = """
        query searchDisease($queryString: String!) {
            search(queryString: $queryString, entityNames: ["disease"], page: {size: 5, index: 0}) {
                hits {
                    id
                    name
                    entity
                }
            }
        }
        """

        result = self._execute_query(query, {"queryString": disease_name})
        hits = result.get("search", {}).get("hits", [])

        for hit in hits:
            if hit.get("entity") == "disease":
                return {
                    "id": hit.get("id", ""),
                    "name": hit.get("name", "")
                }

        return None

    def get_gene_disease_summary(self, gene_symbol: str, top_n: int = 10) -> Dict[str, Any]:
        """
        获取基因的疾病关联摘要（便捷方法）

        Args:
            gene_symbol: 基因符号
            top_n: 返回 Top N 疾病

        Returns:
            {
                "gene": {...},
                "diseases": [...],
                "total_count": int
            }
        """
        # 1. 搜索基因获取 Ensembl ID
        target = self.search_target(gene_symbol)
        if not target:
            logger.warning(f"在 Open Targets 中未找到基因: {gene_symbol}")
            return {"gene": None, "diseases": [], "total_count": 0}

        # 2. 获取疾病关联
        associations = self.get_disease_associations(target.ensembl_id, size=top_n)

        return {
            "gene": {
                "symbol": target.symbol,
                "ensembl_id": target.ensembl_id,
                "name": target.name,
                "description": target.description
            },
            "diseases": [
                {
                    "id": a.disease_id,
                    "name": a.disease_name,
                    "score": a.score,
                    "evidence_types": a.evidence_types
                }
                for a in associations
            ],
            "total_count": len(associations)
        }


# 便捷函数
def get_gene_disease_associations(gene_symbol: str, min_score: float = 0.1) -> List[Dict]:
    """
    获取基因的疾病关联（便捷函数）

    Args:
        gene_symbol: 基因符号
        min_score: 最小关联评分

    Returns:
        疾病关联列表
    """
    client = OpenTargetsClient()
    summary = client.get_gene_disease_summary(gene_symbol)
    return [d for d in summary.get("diseases", []) if d.get("score", 0) >= min_score]


if __name__ == "__main__":
    # 测试
    client = OpenTargetsClient()

    print("=== 测试 BRAF 基因-疾病关联 ===")
    summary = client.get_gene_disease_summary("BRAF", top_n=5)

    if summary.get("gene"):
        print(f"基因: {summary['gene']['symbol']} ({summary['gene']['ensembl_id']})")
        print(f"描述: {summary['gene']['description'][:100]}...")
        print(f"\nTop {len(summary['diseases'])} 疾病关联:")
        for d in summary["diseases"]:
            print(f"  - {d['name']} (评分: {d['score']:.3f})")
            print(f"    证据类型: {', '.join(d['evidence_types'])}")

    print("\n=== 测试 BRAF 靶向药物 ===")
    target = client.search_target("BRAF")
    if target:
        drugs = client.get_associated_drugs(target.ensembl_id, size=5)
        for drug in drugs:
            status = "已批准" if drug.is_approved else f"Phase {drug.phase}"
            print(f"  - {drug.drug_name} ({drug.drug_type}) [{status}]")
            print(f"    作用机制: {drug.mechanism_of_action}")
            print(f"    适应症: {drug.indication}")
