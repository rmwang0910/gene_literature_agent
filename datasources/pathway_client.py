"""
通路数据库客户端

提供 KEGG 和 Reactome 通路查询功能。
使用 REST API 实时查询，替代本地 RAG 索引。
"""

import requests
import logging
import re
from typing import Optional, List, Dict, Any
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)


@dataclass
class Pathway:
    """通路数据"""
    id: str  # 如 "hsa04010" (KEGG) 或 "R-HSA-109582" (Reactome)
    name: str
    description: str = ""
    source: str = ""  # "KEGG" 或 "Reactome"
    genes: List[str] = field(default_factory=list)
    url: str = ""


@dataclass
class PathwayEnrichment:
    """通路富集结果"""
    pathway: Pathway
    p_value: float = 1.0
    gene_ratio: str = ""  # 如 "5/20"
    overlap_genes: List[str] = field(default_factory=list)


class KEGGClient:
    """KEGG REST API 客户端"""

    BASE_URL = "https://rest.kegg.jp"

    def __init__(self, timeout: int = 30):
        self.timeout = timeout
        self._session = requests.Session()

    def get_gene_pathways(self, gene: str, organism: str = "hsa") -> List[Pathway]:
        """
        获取基因参与的 KEGG 通路

        Args:
            gene: 基因符号
            organism: KEGG 物种代码（hsa=人类）

        Returns:
            通路列表
        """
        try:
            # 首先通过基因符号查找 KEGG 基因 ID
            response = self._session.get(
                f"{self.BASE_URL}/find/genes/{gene}",
                timeout=self.timeout
            )
            response.raise_for_status()

            # 解析结果找到人类基因
            gene_id = None
            for line in response.text.strip().split('\n'):
                if line and f"{organism}:" in line:
                    parts = line.split('\t')
                    if parts:
                        gene_id = parts[0]
                        break

            if not gene_id:
                return []

            # 获取基因关联的通路
            response = self._session.get(
                f"{self.BASE_URL}/link/pathway/{gene_id}",
                timeout=self.timeout
            )
            response.raise_for_status()

            pathway_ids = []
            for line in response.text.strip().split('\n'):
                if line:
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        pathway_ids.append(parts[1].replace("path:", ""))

            # 获取通路详细信息
            pathways = []
            for pid in pathway_ids[:20]:  # 限制数量
                pathway = self.get_pathway_info(pid)
                if pathway:
                    pathways.append(pathway)

            return pathways

        except requests.RequestException as e:
            logger.error(f"KEGG API 请求失败: {e}")
            return []

    def get_pathway_info(self, pathway_id: str) -> Optional[Pathway]:
        """
        获取通路详细信息

        Args:
            pathway_id: KEGG 通路 ID

        Returns:
            通路信息
        """
        try:
            response = self._session.get(
                f"{self.BASE_URL}/get/{pathway_id}",
                timeout=self.timeout
            )
            response.raise_for_status()

            text = response.text
            name = ""
            description = ""

            for line in text.split('\n'):
                if line.startswith("NAME"):
                    name = line.replace("NAME", "").strip()
                    # 去掉物种后缀
                    name = re.sub(r' - Homo sapiens \(human\)$', '', name)
                elif line.startswith("DESCRIPTION"):
                    description = line.replace("DESCRIPTION", "").strip()

            if name:
                return Pathway(
                    id=pathway_id,
                    name=name,
                    description=description,
                    source="KEGG",
                    url=f"https://www.kegg.jp/kegg-bin/show_pathway?{pathway_id}"
                )

            return None

        except requests.RequestException as e:
            logger.debug(f"获取 KEGG 通路信息失败: {e}")
            return None

    def search_pathways(self, keyword: str, organism: str = "hsa") -> List[Pathway]:
        """
        搜索 KEGG 通路

        Args:
            keyword: 关键词
            organism: KEGG 物种代码

        Returns:
            通路列表
        """
        try:
            response = self._session.get(
                f"{self.BASE_URL}/find/pathway/{keyword}",
                timeout=self.timeout
            )
            response.raise_for_status()

            pathways = []
            for line in response.text.strip().split('\n'):
                if line:
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        pid = parts[0].replace("path:", "")
                        # 只保留指定物种的通路
                        if pid.startswith(organism):
                            name = parts[1] if len(parts) > 1 else ""
                            pathways.append(Pathway(
                                id=pid,
                                name=name,
                                source="KEGG",
                                url=f"https://www.kegg.jp/kegg-bin/show_pathway?{pid}"
                            ))

            return pathways[:20]

        except requests.RequestException as e:
            logger.error(f"搜索 KEGG 通路失败: {e}")
            return []


class ReactomeClient:
    """Reactome REST API 客户端"""

    BASE_URL = "https://reactome.org/ContentService"

    def __init__(self, timeout: int = 30):
        self.timeout = timeout
        self._session = requests.Session()
        self._session.headers.update({
            'Accept': 'application/json'
        })

    def get_gene_pathways(self, gene: str) -> List[Pathway]:
        """
        获取基因参与的 Reactome 通路

        Args:
            gene: 基因符号

        Returns:
            通路列表
        """
        try:
            # 1. 首先搜索基因获取 Reactome 实体 ID
            search_response = self._session.get(
                f"{self.BASE_URL}/search/query",
                params={
                    "query": gene,
                    "species": "Homo sapiens",
                    "types": "Protein",
                    "cluster": "true"
                },
                timeout=self.timeout
            )
            search_response.raise_for_status()
            search_data = search_response.json()

            # 提取实体 ID
            entity_id = None
            results = search_data.get("results", [])
            if results:
                entries = results[0].get("entries", [])
                if entries:
                    # 找到精确匹配的基因（不是突变体）
                    for entry in entries:
                        name = entry.get("name", "")
                        # 移除 HTML 标签
                        clean_name = re.sub(r'<[^>]+>', '', name)
                        if clean_name.upper() == gene.upper():
                            entity_id = entry.get("stId")
                            break
                    # 如果没有精确匹配，使用第一个
                    if not entity_id and entries:
                        entity_id = entries[0].get("stId")

            if not entity_id:
                return []

            # 2. 使用实体 ID 查询参与的通路
            response = self._session.get(
                f"{self.BASE_URL}/data/pathways/low/entity/{entity_id}",
                timeout=self.timeout
            )
            response.raise_for_status()
            data = response.json()

            pathways = []
            for item in data[:20]:
                pathways.append(Pathway(
                    id=item.get("stId", ""),
                    name=item.get("displayName", ""),
                    description="",
                    source="Reactome",
                    url=f"https://reactome.org/content/detail/{item.get('stId', '')}"
                ))

            return pathways

        except requests.RequestException as e:
            logger.debug(f"Reactome API 请求失败: {e}")
            return []

    def get_pathway_info(self, pathway_id: str) -> Optional[Pathway]:
        """
        获取通路详细信息

        Args:
            pathway_id: Reactome 通路 ID

        Returns:
            通路信息
        """
        try:
            response = self._session.get(
                f"{self.BASE_URL}/data/query/{pathway_id}",
                timeout=self.timeout
            )
            response.raise_for_status()
            data = response.json()

            return Pathway(
                id=pathway_id,
                name=data.get("displayName", ""),
                description=data.get("summation", [{}])[0].get("text", "") if data.get("summation") else "",
                source="Reactome",
                url=f"https://reactome.org/content/detail/{pathway_id}"
            )

        except requests.RequestException as e:
            logger.debug(f"获取 Reactome 通路信息失败: {e}")
            return None

    def search_pathways(self, keyword: str) -> List[Pathway]:
        """
        搜索 Reactome 通路

        Args:
            keyword: 关键词

        Returns:
            通路列表
        """
        try:
            response = self._session.get(
                f"{self.BASE_URL}/search/query",
                params={
                    "query": keyword,
                    "species": "Homo sapiens",
                    "types": "Pathway",
                    "cluster": "true"
                },
                timeout=self.timeout
            )
            response.raise_for_status()
            data = response.json()

            pathways = []
            results = data.get("results", [])
            for result in results:
                entries = result.get("entries", [])
                for entry in entries[:10]:
                    pathways.append(Pathway(
                        id=entry.get("stId", ""),
                        name=entry.get("name", ""),
                        description=entry.get("summation", ""),
                        source="Reactome",
                        url=f"https://reactome.org/content/detail/{entry.get('stId', '')}"
                    ))

            return pathways[:20]

        except requests.RequestException as e:
            logger.error(f"搜索 Reactome 通路失败: {e}")
            return []


class PathwayClient:
    """统一通路客户端（整合 KEGG 和 Reactome）"""

    def __init__(self, timeout: int = 30):
        self.kegg = KEGGClient(timeout=timeout)
        self.reactome = ReactomeClient(timeout=timeout)

    def get_gene_pathways(
        self,
        gene: str,
        sources: Optional[List[str]] = None
    ) -> List[Pathway]:
        """
        获取基因参与的通路（整合多数据源）

        Args:
            gene: 基因符号
            sources: 数据源列表，默认 ["KEGG", "Reactome"]

        Returns:
            通路列表（去重）
        """
        if sources is None:
            sources = ["KEGG", "Reactome"]

        pathways = []

        if "KEGG" in sources:
            kegg_pathways = self.kegg.get_gene_pathways(gene)
            pathways.extend(kegg_pathways)
            logger.info(f"KEGG: 找到 {len(kegg_pathways)} 个 {gene} 相关通路")

        if "Reactome" in sources:
            reactome_pathways = self.reactome.get_gene_pathways(gene)
            pathways.extend(reactome_pathways)
            logger.info(f"Reactome: 找到 {len(reactome_pathways)} 个 {gene} 相关通路")

        return pathways

    def search_pathways(
        self,
        keyword: str,
        sources: Optional[List[str]] = None
    ) -> List[Pathway]:
        """
        搜索通路

        Args:
            keyword: 关键词
            sources: 数据源列表

        Returns:
            通路列表
        """
        if sources is None:
            sources = ["KEGG", "Reactome"]

        pathways = []

        if "KEGG" in sources:
            pathways.extend(self.kegg.search_pathways(keyword))

        if "Reactome" in sources:
            pathways.extend(self.reactome.search_pathways(keyword))

        return pathways

    def get_pathway_summary(self, gene: str) -> Dict[str, Any]:
        """
        获取基因通路摘要

        Args:
            gene: 基因符号

        Returns:
            通路摘要
        """
        pathways = self.get_gene_pathways(gene)

        kegg_pathways = [p for p in pathways if p.source == "KEGG"]
        reactome_pathways = [p for p in pathways if p.source == "Reactome"]

        return {
            "gene": gene,
            "total_pathways": len(pathways),
            "kegg_count": len(kegg_pathways),
            "reactome_count": len(reactome_pathways),
            "pathways": [
                {
                    "id": p.id,
                    "name": p.name,
                    "source": p.source,
                    "url": p.url
                }
                for p in pathways
            ]
        }


# 便捷函数
def get_gene_pathways(gene: str) -> List[Dict[str, str]]:
    """
    获取基因参与的通路（便捷函数）

    Args:
        gene: 基因符号

    Returns:
        通路列表
    """
    client = PathwayClient()
    pathways = client.get_gene_pathways(gene)
    return [
        {
            "id": p.id,
            "name": p.name,
            "source": p.source,
            "url": p.url
        }
        for p in pathways
    ]


if __name__ == "__main__":
    # 测试
    client = PathwayClient()

    print("=== 测试 BRAF 通路 ===")
    summary = client.get_pathway_summary("BRAF")
    print(f"基因: {summary['gene']}")
    print(f"总通路数: {summary['total_pathways']}")
    print(f"  KEGG: {summary['kegg_count']}")
    print(f"  Reactome: {summary['reactome_count']}")
    print("\n通路列表:")
    for p in summary['pathways'][:10]:
        print(f"  - [{p['source']}] {p['name']}")
        print(f"    {p['url']}")

    print("\n=== 搜索 MAPK 通路 ===")
    mapk_pathways = client.search_pathways("MAPK")
    for p in mapk_pathways[:5]:
        print(f"  - [{p.source}] {p.name}")
