"""
热点突变数据下载器 - 从 COSMIC、OncoKB 等数据库获取突变热点信息
"""
import os
import json
import logging
import requests
import time
import re
from typing import Dict, List, Set, Tuple, Optional
from dataclasses import dataclass, field
from pathlib import Path
from collections import Counter

# 配置日志
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class HotspotMutation:
    """热点突变数据"""
    gene: str
    mutation: str  # 如 R175H
    position: int  # 氨基酸位置
    ref_aa: str  # 参考氨基酸
    alt_aa: str  # 变异氨基酸
    cancer_types: List[str] = field(default_factory=list)  # 相关癌症类型
    frequency: int = 0  # 在数据库中的频次
    source: str = ""  # 数据来源
    clinical_significance: str = ""  # 临床意义


class COSMICDownloader:
    """COSMIC/TCGA 数据库下载器"""

    BASE_URL = "https://cancer.sanger.ac.uk/api"
    CBIOPORTAL_URL = "https://www.cbioportal.org/api"

    def __init__(self, email: str = None):
        self.email = email
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'GeneLiteratureAgent/1.0',
            'Accept': 'application/json'
        })

    def download_from_cbioportal(self) -> Dict[str, List[HotspotMutation]]:
        """
        从 cBioPortal 获取真实的TCGA突变数据（免费，无需API Key）

        Returns:
            {基因: [HotspotMutation列表]}
        """
        logger.info("正在从 cBioPortal 获取TCGA突变数据...")

        # 主要癌症类型的TCGA研究
        studies = {
            "brca_tcga": "乳腺癌",
            "luad_tcga": "肺腺癌",
            "lusc_tcga": "肺鳞癌",
            "coadread_tcga": "结直肠癌",
            "hnsc_tcga": "头颈癌",
            "gbm_tcga": "胶质母细胞瘤",
            "kirc_tcga": "肾透明细胞癌",
            "lihc_tcga": "肝癌",
            "paad_tcga": "胰腺癌",
            "ov_tcga": "卵巢癌",
        }

        # 关键癌症基因及其Entrez ID
        genes = {
            "TP53": 7157, "KRAS": 3845, "BRAF": 673, "EGFR": 1956,
            "PIK3CA": 5290, "NRAS": 4893, "IDH1": 3417, "IDH2": 3418,
            "PTEN": 5728, "APC": 324, "VHL": 7428, "JAK2": 3717,
            "KIT": 3815, "ALK": 238, "ROS1": 6098, "HRAS": 3265,
        }

        # 收集突变数据
        all_mutations = {gene: {} for gene in genes}

        total_studies = len(studies)
        for i, (study_id, cancer_name) in enumerate(studies.items()):
            logger.info(f"处理研究 [{i+1}/{total_studies}]: {cancer_name}")
            molecular_profile = f"{study_id}_mutations"

            for gene_name, gene_id in genes.items():
                try:
                    url = f"{self.CBIOPORTAL_URL}/mutations/fetch"
                    response = self.session.post(
                        url,
                        headers={"Content-Type": "application/json"},
                        json={
                            "entrezGeneIds": [gene_id],
                            "molecularProfileIds": [molecular_profile]
                        },
                        timeout=30
                    )

                    if response.status_code == 200:
                        mutations = response.json()
                        for m in mutations:
                            protein_change = m.get('proteinChange', '')
                            if protein_change and protein_change != 'null':
                                # 过滤掉splice突变
                                if 'splice' not in protein_change.lower():
                                    if protein_change not in all_mutations[gene_name]:
                                        all_mutations[gene_name][protein_change] = 0
                                    all_mutations[gene_name][protein_change] += 1

                    time.sleep(0.2)  # 避免请求过快

                except Exception as e:
                    logger.debug(f"获取 {gene_name} @ {study_id} 失败: {e}")

        # 构建HotspotMutation对象
        hotspots = {}
        for gene, mutations in all_mutations.items():
            if mutations:
                hotspots[gene] = []
                for mutation, count in mutations.items():
                    # 解析突变位点
                    match = re.match(r'([A-Z])(\d+)([A-Z]|\*)', mutation)
                    if match:
                        ref_aa, position, alt_aa = match.groups()
                        # 只保留频次>=3的作为热点
                        if count >= 3:
                            hotspots[gene].append(HotspotMutation(
                                gene=gene,
                                mutation=mutation,
                                position=int(position),
                                ref_aa=ref_aa,
                                alt_aa=alt_aa,
                                cancer_types=["多种癌症"],
                                frequency=count,
                                source="TCGA"
                            ))

                # 按频次排序
                hotspots[gene].sort(key=lambda x: x.frequency, reverse=True)

        logger.info(f"从TCGA获取了 {len(hotspots)} 个基因的热点突变")
        return hotspots

    def download_cancer_gene_census(self) -> Dict[str, dict]:
        """
        下载癌症基因普查列表

        Returns:
            基因信息字典
        """
        logger.info("正在下载 COSMIC 癌症基因普查...")

        # COSMIC 癌症基因普查可公开访问
        url = "https://cancer.sanger.ac.uk/census/Census_all.tsv"

        try:
            response = self.session.get(url, timeout=60)
            response.raise_for_status()

            genes = {}
            lines = response.text.strip().split('\n')
            if lines:
                header = lines[0].split('\t')
                for line in lines[1:]:
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        gene_symbol = parts[0]
                        genes[gene_symbol] = {
                            'name': gene_symbol,
                            'entrez': parts[1] if len(parts) > 1 else '',
                            'tier': parts[7] if len(parts) > 7 else '',
                            'hallmark': parts[8] if len(parts) > 8 else '',
                            'somatic': parts[9] if len(parts) > 9 else '',
                            'germline': parts[10] if len(parts) > 10 else '',
                        }

            logger.info(f"获取 {len(genes)} 个癌症相关基因")
            return genes

        except Exception as e:
            logger.warning(f"COSMIC 基因普查下载失败: {e}")
            return {}

    def download_hotspot_mutations(self, use_online: bool = True) -> Dict[str, List[HotspotMutation]]:
        """
        下载热点突变数据

        仅从 cBioPortal 获取真实TCGA数据，不使用硬编码数据。

        Args:
            use_online: 是否从在线API获取数据

        Returns:
            {基因: [HotspotMutation列表]}

        Raises:
            RuntimeError: 如果数据获取失败
        """
        logger.info("正在获取热点突变数据（来源：cBioPortal TCGA API）...")

        if use_online:
            # 从 cBioPortal 获取真实TCGA数据（免费，无需API Key）
            try:
                hotspots = self.download_from_cbioportal()
                if hotspots:
                    logger.info(f"成功从 cBioPortal 获取 {len(hotspots)} 个基因的热点突变")
                    return hotspots
                else:
                    logger.error("cBioPortal 返回空数据，请检查网络连接")
            except Exception as e:
                logger.error(f"cBioPortal 数据获取失败: {e}")

        # 无硬编码备份数据
        raise RuntimeError(
            "热点突变数据获取失败。请确保网络连接正常，然后运行：\n"
            "  python hotspot_downloader.py --build\n"
            "数据将自动从 cBioPortal TCGA API 下载（免费，无需API Key）"
        )

    def _get_literature_hotspots(self) -> Dict[str, List[HotspotMutation]]:
        """
        备用数据获取方法（已禁用硬编码数据）

        所有热点突变数据应从数据库API获取，不接受硬编码数据。
        如果在线数据获取失败，返回空字典。

        Returns:
            {基因: [HotspotMutation列表]}
        """
        # 不使用硬编码数据，所有数据必须来自数据库API
        logger.warning("在线数据获取失败，且不使用硬编码数据。请检查网络连接后重试。")
        logger.info("提示: 运行 python hotspot_downloader.py --build 重新从cBioPortal获取数据")
        return {}


class OncoKBDownloader:
    """OncoKB 数据下载器"""

    BASE_URL = "https://www.oncokb.org/api/v1"

    def __init__(self, api_key: str = None):
        self.api_key = api_key
        self.session = requests.Session()
        if api_key:
            self.session.headers.update({
                'Authorization': f'Bearer {api_key}'
            })

    def download_annotated_variants(self) -> Dict[str, List[dict]]:
        """
        下载 OncoKB 注释的变异

        Returns:
            {基因: [变异信息列表]}
        """
        if not self.api_key:
            logger.warning("OncoKB 需要 API Key，跳过")
            return {}

        logger.info("正在下载 OncoKB 变异注释...")

        try:
            url = f"{self.BASE_URL}/annotate/mutations"
            response = self.session.get(url, timeout=60)

            if response.status_code == 200:
                data = response.json()
                logger.info(f"获取 OncoKB 变异数据")
                return data
            else:
                logger.warning(f"OncoKB API 返回错误: {response.status_code}")
                return {}

        except Exception as e:
            logger.warning(f"OncoKB 下载失败: {e}")
            return {}


class HotspotBuilder:
    """热点突变数据构建器"""

    def __init__(self, output_dir: str = None):
        self.output_dir = output_dir or os.path.join(os.path.dirname(__file__), "..", "data", "mutation")
        self.cosmic = COSMICDownloader()
        self.hotspots: Dict[str, List[HotspotMutation]] = {}

    def build_hotspot_database(self) -> Dict[str, List[HotspotMutation]]:
        """
        构建热点突变数据库

        Returns:
            {基因: [HotspotMutation列表]}
        """
        logger.info("正在构建热点突变数据库...")

        # 1. 从 COSMIC 获取
        cosmic_hotspots = self.cosmic.download_hotspot_mutations()

        # 合并数据
        self.hotspots = cosmic_hotspots

        logger.info(f"热点突变数据库构建完成，共 {len(self.hotspots)} 个基因")

        return self.hotspots

    def get_position_hotspots(self, gene: str) -> Dict[int, List[str]]:
        """
        获取基因的氨基酸位置热点

        Args:
            gene: 基因名

        Returns:
            {位置: [突变列表]}
        """
        if gene not in self.hotspots:
            # 尝试加载已保存的数据
            self._load_hotspots()

        if gene not in self.hotspots:
            return {}

        position_hotspots = {}
        for mut in self.hotspots[gene]:
            if mut.position not in position_hotspots:
                position_hotspots[mut.position] = []
            position_hotspots[mut.position].append(mut.mutation)

        return position_hotspots

    def is_hotspot(self, gene: str, mutation: str) -> Tuple[bool, str]:
        """
        判断突变是否为热点突变

        Args:
            gene: 基因名
            mutation: 突变（如 R175H）

        Returns:
            (是否热点, 说明)
        """
        # 解析突变
        match = re.match(r'([A-Z])(\d+)([A-Z]|\*)', mutation, re.IGNORECASE)
        if not match:
            return False, ""

        ref_aa, position, alt_aa = match.groups()
        position = int(position)

        # 获取该基因的热点位置
        position_hotspots = self.get_position_hotspots(gene.upper())

        if position in position_hotspots:
            # 检查是否为已知热点突变
            if mutation.upper() in [m.upper() for m in position_hotspots[position]]:
                return True, f"热点突变({mutation})"
            else:
                # 该位置是热点位置，但具体突变类型不同
                known_mutations = position_hotspots[position]
                return True, f"热点位置({position}，已知突变: {', '.join(known_mutations[:3])})"

        return False, ""

    def save(self, filename: str = "hotspot_mutations.json"):
        """
        保存热点突变数据

        Args:
            filename: 输出文件名
        """
        output_path = os.path.join(self.output_dir, filename)

        # 统计数据来源
        sources = set()
        for gene, mutations in self.hotspots.items():
            for m in mutations:
                if m.source:
                    sources.add(m.source)
        source_str = ", ".join(sorted(sources)) if sources else "Unknown"

        data = {
            "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
            "source": f"cBioPortal TCGA API ({source_str})",
            "data_url": "https://www.cbioportal.org/api",
            "description": "数据来自TCGA项目，通过cBioPortal API免费获取，无需API Key",
            "genes": {}
        }

        for gene, mutations in self.hotspots.items():
            data["genes"][gene] = [
                {
                    "mutation": m.mutation,
                    "position": m.position,
                    "ref_aa": m.ref_aa,
                    "alt_aa": m.alt_aa,
                    "cancer_types": m.cancer_types,
                    "frequency": m.frequency,
                    "source": m.source
                }
                for m in mutations
            ]

        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(data, f, ensure_ascii=False, indent=2)

        logger.info(f"热点突变数据已保存: {output_path}")
        return output_path

    def _load_hotspots(self):
        """加载已保存的热点数据"""
        cache_path = os.path.join(self.output_dir, "hotspot_mutations.json")

        if not os.path.exists(cache_path):
            self.build_hotspot_database()
            return

        try:
            with open(cache_path, 'r', encoding='utf-8') as f:
                data = json.load(f)

            self.hotspots = {}
            for gene, mutations in data.get("genes", {}).items():
                self.hotspots[gene] = [
                    HotspotMutation(**m) for m in mutations
                ]

            logger.info(f"从缓存加载了 {len(self.hotspots)} 个基因的热点突变")

        except Exception as e:
            logger.warning(f"加载热点缓存失败: {e}")
            self.build_hotspot_database()


# 全局热点检索器
_hotspot_builder: Optional[HotspotBuilder] = None


def get_hotspot_builder() -> HotspotBuilder:
    """获取全局热点构建器实例"""
    global _hotspot_builder
    if _hotspot_builder is None:
        _hotspot_builder = HotspotBuilder()
    return _hotspot_builder


def is_hotspot_mutation(gene: str, mutation: str) -> Tuple[bool, str]:
    """
    判断突变是否为热点突变（便捷函数）

    Args:
        gene: 基因名
        mutation: 突变（如 R175H）

    Returns:
        (是否热点, 说明)
    """
    builder = get_hotspot_builder()
    return builder.is_hotspot(gene, mutation)


def annotate_mutations_with_hotspot(
    gene: str,
    mutations: List[str]
) -> List[Dict[str, str]]:
    """
    标注突变列表中的热点突变

    Args:
        gene: 基因名
        mutations: 突变列表

    Returns:
        标注后的突变列表 [{"mutation": "R175H", "is_hotspot": True, "label": "🔥 热点突变"}]
    """
    builder = get_hotspot_builder()
    result = []

    for mut in mutations:
        is_hs, description = builder.is_hotspot(gene, mut)
        result.append({
            "mutation": mut,
            "is_hotspot": is_hs,
            "label": f"🔥 {description}" if is_hs else mut,
            "description": description
        })

    return result


def main():
    """主函数"""
    import argparse

    parser = argparse.ArgumentParser(description="热点突变数据下载器")
    parser.add_argument("--build", action="store_true", help="构建热点突变数据库")
    parser.add_argument("--query", nargs=2, metavar=("GENE", "MUTATION"), help="查询突变是否为热点")
    parser.add_argument("--output-dir", default=None, help="输出目录")
    args = parser.parse_args()

    print("=" * 60)
    print("热点突变数据下载器")
    print("=" * 60)

    builder = HotspotBuilder(output_dir=args.output_dir)

    if args.build:
        builder.build_hotspot_database()
        builder.save()

        print()
        print("=" * 60)
        print("构建完成!")
        print(f"基因数量: {len(builder.hotspots)}")

        # 显示部分结果
        print()
        print("热点突变示例:")
        for gene in ["TP53", "KRAS", "BRAF"]:
            if gene in builder.hotspots:
                mutations = [m.mutation for m in builder.hotspots[gene][:5]]
                print(f"  {gene}: {', '.join(mutations)}")

        print("=" * 60)

    if args.query:
        gene, mutation = args.query
        is_hs, description = builder.is_hotspot(gene, mutation)
        if is_hs:
            print(f"✅ {gene} {mutation} 是热点突变")
            print(f"   {description}")
        else:
            print(f"❌ {gene} {mutation} 不是已知热点突变")


if __name__ == "__main__":
    main()