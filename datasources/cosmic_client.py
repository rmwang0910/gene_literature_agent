"""
COSMIC 数据库客户端

提供癌症体细胞突变数据查询功能。
注意：COSMIC 需要认证才能下载数据，本模块提供：
1. 本地数据文件访问
2. Cancer Gene Census 查询
3. 热点突变数据
"""

import json
import logging
import os
from pathlib import Path
from typing import Optional, List, Dict, Any, Set
from dataclasses import dataclass, field
import requests
import base64

logger = logging.getLogger(__name__)


@dataclass
class CancerMutation:
    """癌症突变数据"""
    gene: str
    mutation: str  # 如 "V600E"
    aa_position: int = 0
    cancer_types: List[str] = field(default_factory=list)
    sample_count: int = 0
    frequency: float = 0.0  # 在该基因所有突变中的占比
    is_hotspot: bool = False
    cosmic_id: str = ""


@dataclass
class CancerGene:
    """Cancer Gene Census 基因"""
    gene: str
    name: str
    role: str  # oncogene, TSG, fusion
    cancer_types: List[str] = field(default_factory=list)
    mutation_types: List[str] = field(default_factory=list)  # Mis, N, F, etc.
    tier: int = 1  # 1 或 2


class COSMICClient:
    """COSMIC 数据库客户端"""

    # COSMIC 下载 API（需要认证）
    DOWNLOAD_API = "https://cancer.sanger.ac.uk/cosmic/file_download"

    def __init__(
        self,
        data_dir: Optional[str] = None,
        email: Optional[str] = None,
        password: Optional[str] = None
    ):
        """
        初始化 COSMIC 客户端

        Args:
            data_dir: 本地数据目录
            email: COSMIC 账号邮箱（用于下载）
            password: COSMIC 账号密码
        """
        self.data_dir = Path(data_dir) if data_dir else Path(__file__).parent.parent / "data" / "cosmic"
        self.data_dir.mkdir(parents=True, exist_ok=True)

        self.email = email or os.environ.get("COSMIC_EMAIL")
        self.password = password or os.environ.get("COSMIC_PASSWORD")

        # 加载本地数据
        self._hotspots: Dict[str, List[Dict]] = {}
        self._cancer_genes: Dict[str, CancerGene] = {}
        self._load_local_data()

    def _load_local_data(self):
        """加载本地数据文件"""
        # 加载热点突变数据
        hotspot_file = self.data_dir / "hotspot_mutations.json"
        if hotspot_file.exists():
            try:
                with open(hotspot_file, 'r', encoding='utf-8') as f:
                    self._hotspots = json.load(f)
                logger.info(f"已加载 {len(self._hotspots)} 个基因的热点突变数据")
            except Exception as e:
                logger.warning(f"加载热点突变数据失败: {e}")

        # 加载 Cancer Gene Census
        cgc_file = self.data_dir / "cancer_gene_census.json"
        if cgc_file.exists():
            try:
                with open(cgc_file, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                    for gene_data in data:
                        gene = gene_data.get("gene", "")
                        self._cancer_genes[gene.upper()] = CancerGene(
                            gene=gene,
                            name=gene_data.get("name", ""),
                            role=gene_data.get("role", ""),
                            cancer_types=gene_data.get("cancer_types", []),
                            mutation_types=gene_data.get("mutation_types", []),
                            tier=gene_data.get("tier", 1)
                        )
                logger.info(f"已加载 {len(self._cancer_genes)} 个癌症基因")
            except Exception as e:
                logger.warning(f"加载 Cancer Gene Census 失败: {e}")

    def is_cancer_gene(self, gene: str) -> bool:
        """
        检查是否为 Cancer Gene Census 中的癌症基因

        Args:
            gene: 基因符号

        Returns:
            是否为已知癌症基因
        """
        return gene.upper() in self._cancer_genes

    def get_cancer_gene_info(self, gene: str) -> Optional[CancerGene]:
        """
        获取癌症基因信息

        Args:
            gene: 基因符号

        Returns:
            癌症基因信息
        """
        return self._cancer_genes.get(gene.upper())

    def get_hotspot_mutations(self, gene: str) -> List[CancerMutation]:
        """
        获取基因的热点突变

        Args:
            gene: 基因符号

        Returns:
            热点突变列表
        """
        gene_upper = gene.upper()

        if gene_upper not in self._hotspots:
            return []

        mutations = []
        for m in self._hotspots[gene_upper]:
            mutations.append(CancerMutation(
                gene=gene_upper,
                mutation=m.get("mutation", ""),
                aa_position=m.get("position", 0),
                cancer_types=m.get("cancer_types", []),
                sample_count=m.get("sample_count", 0),
                frequency=m.get("frequency", 0.0),
                is_hotspot=True,
                cosmic_id=m.get("cosmic_id", "")
            ))

        # 按样本数排序
        mutations.sort(key=lambda x: x.sample_count, reverse=True)
        return mutations

    def is_hotspot(self, gene: str, mutation: str) -> bool:
        """
        检查突变是否为热点

        Args:
            gene: 基因符号
            mutation: 突变描述（如 "V600E"）

        Returns:
            是否为热点突变
        """
        gene_upper = gene.upper()
        mutation_upper = mutation.upper()

        if gene_upper not in self._hotspots:
            return False

        for m in self._hotspots[gene_upper]:
            if m.get("mutation", "").upper() == mutation_upper:
                return True

        return False

    def get_mutation_info(self, gene: str, mutation: str) -> Optional[CancerMutation]:
        """
        获取特定突变的信息

        Args:
            gene: 基因符号
            mutation: 突变描述

        Returns:
            突变信息
        """
        gene_upper = gene.upper()
        mutation_upper = mutation.upper()

        if gene_upper not in self._hotspots:
            return None

        for m in self._hotspots[gene_upper]:
            if m.get("mutation", "").upper() == mutation_upper:
                return CancerMutation(
                    gene=gene_upper,
                    mutation=m.get("mutation", ""),
                    aa_position=m.get("position", 0),
                    cancer_types=m.get("cancer_types", []),
                    sample_count=m.get("sample_count", 0),
                    frequency=m.get("frequency", 0.0),
                    is_hotspot=True,
                    cosmic_id=m.get("cosmic_id", "")
                )

        return None

    def get_gene_cancer_types(self, gene: str) -> List[str]:
        """
        获取与基因相关的癌症类型

        Args:
            gene: 基因符号

        Returns:
            相关癌症类型列表
        """
        cancer_types: Set[str] = set()

        # 从 Cancer Gene Census 获取
        cancer_gene = self.get_cancer_gene_info(gene)
        if cancer_gene:
            cancer_types.update(cancer_gene.cancer_types)

        # 从热点突变获取
        hotspots = self.get_hotspot_mutations(gene)
        for m in hotspots:
            cancer_types.update(m.cancer_types)

        return list(cancer_types)

    def get_gene_role(self, gene: str) -> str:
        """
        获取基因在癌症中的角色

        Args:
            gene: 基因符号

        Returns:
            角色描述（oncogene, TSG, fusion, unknown）
        """
        cancer_gene = self.get_cancer_gene_info(gene)
        if cancer_gene:
            return cancer_gene.role
        return "unknown"

    def download_file(
        self,
        filepath: str,
        output_filename: Optional[str] = None
    ) -> Optional[str]:
        """
        从 COSMIC 下载数据文件（需要认证）

        Args:
            filepath: COSMIC 文件路径（如 "GRCh38/cosmic/latest/cancer_gene_census.csv"）
            output_filename: 输出文件名

        Returns:
            下载的文件路径，失败返回 None
        """
        if not self.email or not self.password:
            logger.error("需要 COSMIC 账号才能下载数据。请设置 COSMIC_EMAIL 和 COSMIC_PASSWORD 环境变量")
            return None

        try:
            # Step 1: 获取下载 URL
            auth = base64.b64encode(f"{self.email}:{self.password}".encode()).decode()
            headers = {"Authorization": f"Basic {auth}"}

            response = requests.get(
                f"{self.DOWNLOAD_API}/{filepath}",
                headers=headers,
                timeout=30
            )
            response.raise_for_status()
            download_url = response.json().get("url")

            if not download_url:
                logger.error("无法获取下载 URL")
                return None

            # Step 2: 下载文件
            output_path = self.data_dir / (output_filename or filepath.split("/")[-1])
            response = requests.get(download_url, stream=True, timeout=300)
            response.raise_for_status()

            with open(output_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)

            logger.info(f"下载完成: {output_path}")
            return str(output_path)

        except requests.RequestException as e:
            logger.error(f"下载 COSMIC 数据失败: {e}")
            return None

    def get_mutation_summary(self, gene: str) -> Dict[str, Any]:
        """
        获取基因突变摘要

        Args:
            gene: 基因符号

        Returns:
            突变摘要信息
        """
        cancer_gene = self.get_cancer_gene_info(gene)
        hotspots = self.get_hotspot_mutations(gene)
        cancer_types = self.get_gene_cancer_types(gene)

        return {
            "gene": gene.upper(),
            "is_cancer_gene": cancer_gene is not None,
            "role": cancer_gene.role if cancer_gene else "unknown",
            "tier": cancer_gene.tier if cancer_gene else None,
            "hotspot_count": len(hotspots),
            "top_hotspots": [
                {
                    "mutation": m.mutation,
                    "sample_count": m.sample_count,
                    "frequency": m.frequency,
                    "cancer_types": m.cancer_types[:3]
                }
                for m in hotspots[:5]
            ],
            "cancer_types": cancer_types[:10],
            "source": "COSMIC"
        }


# 便捷函数
def is_hotspot_mutation(gene: str, mutation: str) -> bool:
    """检查是否为热点突变（便捷函数）"""
    client = COSMICClient()
    return client.is_hotspot(gene, mutation)


def get_cancer_gene_role(gene: str) -> str:
    """获取癌症基因角色（便捷函数）"""
    client = COSMICClient()
    return client.get_gene_role(gene)


if __name__ == "__main__":
    # 测试
    client = COSMICClient()

    print("=== 测试 BRAF 突变信息 ===")
    summary = client.get_mutation_summary("BRAF")
    print(f"基因: {summary['gene']}")
    print(f"是否癌症基因: {summary['is_cancer_gene']}")
    print(f"角色: {summary['role']}")
    print(f"热点突变数: {summary['hotspot_count']}")
    print("Top 热点:")
    for h in summary['top_hotspots']:
        print(f"  - {h['mutation']}: {h['sample_count']} 样本")
    print(f"相关癌症: {', '.join(summary['cancer_types'][:5])}")

    print("\n=== 测试 V600E 是否为热点 ===")
    print(f"BRAF V600E 是热点: {client.is_hotspot('BRAF', 'V600E')}")
    print(f"BRAF V600K 是热点: {client.is_hotspot('BRAF', 'V600K')}")
    print(f"BRAF V600X 是热点: {client.is_hotspot('BRAF', 'V600X')}")

    print("\n=== 测试 TP53 角色 ===")
    print(f"TP53 角色: {client.get_gene_role('TP53')}")
    print(f"EGFR 角色: {client.get_gene_role('EGFR')}")
