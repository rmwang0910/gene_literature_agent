"""
基因规范化模块 - 使用 MyGene.info API 标准化基因名称
"""
import requests
from typing import Tuple, List, Optional, Dict, Any
from dataclasses import dataclass
import json
from pathlib import Path

from ..config import default_config


@dataclass
class GeneInfo:
    """基因信息数据类"""
    symbol: str  # 官方基因符号
    name: str  # 基因全称
    aliases: List[str]  # 别名列表
    entrez_id: Optional[str] = None  # Entrez Gene ID
    ensembl_id: Optional[str] = None  # Ensembl ID
    taxid: int = 9606  # 物种ID，默认人类


class GeneNormalizer:
    """基因名称规范化类"""

    MYGENE_BASE_URL = "https://mygene.info/v3"

    def __init__(self, cache_dir: Optional[str] = None):
        """
        初始化基因规范化器

        Args:
            cache_dir: 缓存目录路径
        """
        self.cache_dir = Path(cache_dir or default_config.cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.cache_file = self.cache_dir / "gene_cache.json"
        self._cache = self._load_cache()

    def _load_cache(self) -> Dict[str, Any]:
        """加载缓存"""
        if self.cache_file.exists():
            try:
                with open(self.cache_file, 'r', encoding='utf-8') as f:
                    return json.load(f)
            except (json.JSONDecodeError, IOError):
                return {}
        return {}

    def _save_cache(self):
        """保存缓存"""
        with open(self.cache_file, 'w', encoding='utf-8') as f:
            json.dump(self._cache, f, ensure_ascii=False, indent=2)

    def normalize(self, gene_name: str, species: str = "human") -> Optional[GeneInfo]:
        """
        规范化基因名称

        Args:
            gene_name: 输入的基因名称（可以是符号、别名或全称）
            species: 物种名称，默认人类

        Returns:
            GeneInfo 对象，如果未找到返回 None
        """
        # 检查缓存
        cache_key = f"{gene_name.lower()}_{species}"
        if cache_key in self._cache:
            cached = self._cache[cache_key]
            return GeneInfo(**cached)

        # 调用 MyGene.info API
        try:
            # 构建查询
            taxid = self._get_taxid(species)
            url = f"{self.MYGENE_BASE_URL}/query"
            params = {
                "q": gene_name,
                "fields": "symbol,name,alias,entrezgene,ensembl.gene,taxid",
                "species": taxid,
                "size": 1
            }

            response = requests.get(url, params=params, timeout=10)
            response.raise_for_status()
            data = response.json()

            if not data.get("hits"):
                return None

            hit = data["hits"][0]

            # 解析别名
            aliases = []
            if "alias" in hit:
                alias_data = hit["alias"]
                if isinstance(alias_data, list):
                    aliases = [str(a) for a in alias_data]
                else:
                    aliases = [str(alias_data)]

            # 解析 Ensembl ID
            ensembl_id = None
            if "ensembl" in hit:
                ensembl_data = hit["ensembl"]
                if isinstance(ensembl_data, dict):
                    ensembl_id = ensembl_data.get("gene")
                elif isinstance(ensembl_data, list) and ensembl_data:
                    ensembl_id = ensembl_data[0].get("gene")

            # 构建 GeneInfo
            gene_info = GeneInfo(
                symbol=hit.get("symbol", gene_name.upper()),
                name=hit.get("name", ""),
                aliases=aliases,
                entrez_id=str(hit.get("entrezgene", "")) if hit.get("entrezgene") else None,
                ensembl_id=ensembl_id,
                taxid=hit.get("taxid", taxid)
            )

            # 缓存结果
            self._cache[cache_key] = {
                "symbol": gene_info.symbol,
                "name": gene_info.name,
                "aliases": gene_info.aliases,
                "entrez_id": gene_info.entrez_id,
                "ensembl_id": gene_info.ensembl_id,
                "taxid": gene_info.taxid
            }
            self._save_cache()

            return gene_info

        except requests.RequestException as e:
            print(f"查询 MyGene.info API 失败: {e}")
            return None

    def normalize_batch(self, gene_names: List[str], species: str = "human") -> Dict[str, Optional[GeneInfo]]:
        """
        批量规范化基因名称

        Args:
            gene_names: 基因名称列表
            species: 物种名称

        Returns:
            字典，键为原始输入，值为 GeneInfo 对象
        """
        results = {}
        for gene in gene_names:
            gene = gene.strip()
            if gene:
                results[gene] = self.normalize(gene, species)
        return results

    def build_search_terms(self, gene_info: GeneInfo, max_aliases: int = 3) -> List[str]:
        """
        构建 PubMed 检索词列表

        Args:
            gene_info: 基因信息
            max_aliases: 最大别名数量

        Returns:
            检索词列表
        """
        terms = [f'"{gene_info.symbol}"[Title/Abstract]']

        # 添加常用别名（限制数量避免检索式过长）
        for alias in gene_info.aliases[:max_aliases]:
            if alias.upper() != gene_info.symbol.upper():
                terms.append(f'"{alias}"[Title/Abstract]')

        return terms

    @staticmethod
    def _get_taxid(species: str) -> int:
        """获取物种的 Taxonomy ID"""
        species_map = {
            "human": 9606,
            "mouse": 10090,
            "rat": 10116,
            "zebrafish": 7955,
            "drosophila": 7227,
            "c. elegans": 6239,
            "yeast": 559292
        }
        return species_map.get(species.lower(), 9606)


def normalize_gene(gene_name: str) -> Tuple[str, List[str]]:
    """
    简单的基因规范化函数

    Args:
        gene_name: 基因名称

    Returns:
        (标准化符号, 别名列表)
    """
    normalizer = GeneNormalizer()
    result = normalizer.normalize(gene_name)
    if result:
        return result.symbol, result.aliases
    return gene_name.upper(), []


if __name__ == "__main__":
    # 测试
    normalizer = GeneNormalizer()
    test_genes = ["BRAF", "TP53", "EGFR", "MYC"]
    for gene in test_genes:
        info = normalizer.normalize(gene)
        if info:
            print(f"{gene} -> {info.symbol} ({info.name})")
            print(f"  别名: {info.aliases[:5]}")
            print(f"  Entrez: {info.entrez_id}")
            print()