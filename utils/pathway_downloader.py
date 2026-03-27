"""
通路数据下载器 - 从 KEGG、Reactome、WikiPathways 获取标准通路名称
生成中英文映射表供系统使用
"""
import os
import re
import json
import logging
import requests
import time
from typing import Dict, List, Tuple, Set
from dataclasses import dataclass, field
from pathlib import Path

# 配置日志
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class PathwayInfo:
    """通路信息"""
    id: str  # 数据库ID
    name: str  # 标准名称
    synonyms: List[str] = field(default_factory=list)  # 同义词
    source: str = ""  # 数据来源


class PathwayDownloader:
    """通路数据下载器"""

    def __init__(self, output_dir: str = None):
        self.output_dir = output_dir or os.path.join(os.path.dirname(__file__), "..", "data", "pathway")
        self.pathways: Dict[str, PathwayInfo] = {}
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'Mozilla/5.0 (compatible; GeneLiteratureAgent/1.0)'
        })

    def download_kegg_pathways(self) -> Dict[str, PathwayInfo]:
        """
        从 KEGG 下载人类通路列表

        Returns:
            通路字典 {name: PathwayInfo}
        """
        logger.info("正在从 KEGG 下载通路数据...")
        pathways = {}

        try:
            # 获取人类通路列表
            url = "https://rest.kegg.jp/list/pathway/hsa"
            response = self.session.get(url, timeout=30)
            response.raise_for_status()

            for line in response.text.strip().split('\n'):
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) >= 2:
                    pathway_id = parts[0]  # path:hsa04010
                    name = parts[1]  # MAPK signaling pathway - Homo sapiens (human)

                    # 清理名称（去除物种后缀）
                    name = re.sub(r'\s*-\s*Homo sapiens.*$', '', name)
                    name = name.strip()

                    if name:
                        key = name.lower()
                        pathways[key] = PathwayInfo(
                            id=pathway_id,
                            name=name,
                            source="KEGG"
                        )

            logger.info(f"KEGG: 获取 {len(pathways)} 条通路")

        except Exception as e:
            logger.warning(f"KEGG 下载失败: {e}")

        return pathways

    def download_reactome_pathways(self) -> Dict[str, PathwayInfo]:
        """
        从 Reactome 下载通路列表

        Returns:
            通路字典 {name: PathwayInfo}
        """
        logger.info("正在从 Reactome 下载通路数据...")
        pathways = {}

        try:
            # 下载 Reactome 通路列表文件
            url = "https://reactome.org/download/current/ReactomePathways.txt"
            response = self.session.get(url, timeout=60)
            response.raise_for_status()

            for line in response.text.strip().split('\n'):
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) >= 3:
                    pathway_id = parts[0]  # R-HSA-162582
                    name = parts[1]  # Signaling by WNT
                    species = parts[2]  # Homo sapiens

                    # 只保留人类通路
                    if species.lower() == "homo sapiens":
                        key = name.lower()
                        if key not in pathways:
                            pathways[key] = PathwayInfo(
                                id=pathway_id,
                                name=name,
                                source="Reactome"
                            )

            logger.info(f"Reactome: 获取 {len(pathways)} 条通路")

        except Exception as e:
            logger.warning(f"Reactome 下载失败: {e}")

        return pathways

    def download_wikipathways(self) -> Dict[str, PathwayInfo]:
        """
        从 WikiPathways 下载人类通路列表

        Returns:
            通路字典 {name: PathwayInfo}
        """
        logger.info("正在从 WikiPathways 下载通路数据...")
        pathways = {}

        try:
            # 使用 SPARQL 查询
            sparql_url = "https://sparql.wikipathways.org/sparql"
            query = """
            SELECT ?pathway ?label WHERE {
                ?pathway a wp:Pathway ;
                         rdfs:label ?label ;
                         wp:organism <http://purl.org/twc/vocab/vann/Homo_sapiens> .
            }
            """

            headers = {
                'Accept': 'application/sparql-results+json'
            }

            response = self.session.get(
                sparql_url,
                params={'query': query, 'format': 'json'},
                headers=headers,
                timeout=60
            )
            response.raise_for_status()

            data = response.json()
            for binding in data.get('results', {}).get('bindings', []):
                pathway_uri = binding.get('pathway', {}).get('value', '')
                name = binding.get('label', {}).get('value', '')

                if name:
                    key = name.lower()
                    if key not in pathways:
                        pathways[key] = PathwayInfo(
                            id=pathway_uri,
                            name=name,
                            source="WikiPathways"
                        )

            logger.info(f"WikiPathways: 获取 {len(pathways)} 条通路")

        except Exception as e:
            logger.warning(f"WikiPathways 下载失败: {e}")

        return pathways

    def generate_chinese_synonyms(self, english_name: str) -> List[str]:
        """
        根据英文通路名称生成中文同义词

        从外部JSON文件加载中文翻译，用户可自行扩展该文件。
        对于包含pathway的名称，自动生成"XXX通路"格式的翻译。

        Args:
            english_name: 英文通路名称

        Returns:
            中文同义词列表
        """
        synonyms = []

        # 1. 从外部JSON文件加载预定义的中文翻译
        synonyms_file = os.path.join(self.output_dir, "pathway_chinese_synonyms.json")
        if os.path.exists(synonyms_file):
            try:
                with open(synonyms_file, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                mappings = data.get("mappings", {})
                name_lower = english_name.lower()
                for eng_name, cn_list in mappings.items():
                    if eng_name.lower() == name_lower:
                        synonyms.extend(cn_list)
                        break
            except Exception as e:
                logger.debug(f"加载中文同义词文件失败: {e}")

        # 2. 通用翻译规则（仅对包含pathway的简短通路名）
        name_lower = english_name.lower()
        if len(english_name) < 50 and "pathway" in name_lower:
            # 提取通路基础名称（移除"pathway"）
            base = re.sub(r'\s*(signaling\s+)?(pathway|pathways)\s*$', '', english_name, flags=re.IGNORECASE)
            # 只有当base不为空且不是纯英文单词组合时才生成
            if base and base.lower() != english_name.lower():
                synonyms.append(f"{base}通路")
                # 如果不包含"signaling"，也生成"信号通路"版本
                if "signaling" not in name_lower:
                    synonyms.append(f"{base}信号通路")

        return list(set(synonyms))

    def build_mapping_table(self) -> Dict[str, str]:
        """
        构建通路名称映射表

        Returns:
            映射字典 {同义词: 标准名称}
        """
        logger.info("正在构建通路映射表...")

        # 下载各数据库通路
        kegg = self.download_kegg_pathways()
        reactome = self.download_reactome_pathways()
        wikipathways = self.download_wikipathways()

        logger.info(f"合并后: KEGG {len(kegg)}, Reactome {len(reactome)}, WikiPathways {len(wikipathways)}")

        # 构建映射表
        mapping = {}

        # 中文同义词专用映射（优先级最高）
        chinese_mapping = {}

        # 常见通路名称变体映射（简写→标准名称）
        # 这些是基于命名规则的逻辑映射，而非硬编码的翻译
        pathway_variants = {
            # 信号通路简写形式 -> KEGG标准名称
            "p53 pathway": "p53 signaling pathway",
            "mapk pathway": "MAPK signaling pathway",
            "pi3k pathway": "PI3K-Akt signaling pathway",
            "pi3k/akt pathway": "PI3K-Akt signaling pathway",
            "pi3k-akt pathway": "PI3K-Akt signaling pathway",
            "akt pathway": "PI3K-Akt signaling pathway",
            "wnt pathway": "Wnt signaling pathway",
            "nf-kb pathway": "NF-kappa B signaling pathway",
            "jak-stat pathway": "JAK-STAT signaling pathway",
            "jak/stat pathway": "JAK-STAT signaling pathway",
            "notch pathway": "Notch signaling pathway",
            "hedgehog pathway": "Hedgehog signaling pathway",
            "hh pathway": "Hedgehog signaling pathway",
            "tgf-beta pathway": "TGF-beta signaling pathway",
            "mtor pathway": "mTOR signaling pathway",
            "hippo pathway": "Hippo signaling pathway",
            "vegf pathway": "VEGF signaling pathway",
            "ras pathway": "Ras signaling pathway",
            "erk pathway": "MAPK signaling pathway",
            "erk1/2 pathway": "MAPK signaling pathway",
            "ampk pathway": "AMPK signaling pathway",
            "dna damage response pathway": "DNA damage response",
        }

        # 1. 首先处理 KEGG（优先级最高，包含简短的标准通路名）
        for key, info in kegg.items():
            standard_name = info.name
            name_lower = standard_name.lower()

            # 添加标准名称
            mapping[standard_name] = standard_name
            mapping[name_lower] = standard_name

            # KEGG 通路的中文同义词（单独存储，避免被覆盖）
            chinese_synonyms = self.generate_chinese_synonyms(standard_name)
            for cn in chinese_synonyms:
                chinese_mapping[cn] = standard_name
                chinese_mapping[cn.lower()] = standard_name

            # signaling pathway 变体
            if "signaling pathway" in name_lower:
                short = re.sub(r'\s*signaling\s+', '', standard_name, flags=re.IGNORECASE).strip()
                if short != standard_name:
                    mapping[short] = standard_name
                    mapping[short.lower()] = standard_name

        # 2. 处理 Reactome（不覆盖 KEGG 和中文映射）
        for key, info in reactome.items():
            standard_name = info.name
            name_lower = standard_name.lower()

            # 仅添加未被 KEGG 覆盖的
            if name_lower not in mapping:
                mapping[standard_name] = standard_name
                mapping[name_lower] = standard_name

        # 3. 处理 WikiPathways（最低优先级）
        for key, info in wikipathways.items():
            standard_name = info.name
            name_lower = standard_name.lower()

            if name_lower not in mapping:
                mapping[standard_name] = standard_name
                mapping[name_lower] = standard_name

        # 4. 添加常见变体映射
        for variant, standard in pathway_variants.items():
            variant_lower = variant.lower()
            # 查找标准名称在映射中的目标
            if standard.lower() in mapping:
                target = mapping[standard.lower()]
            else:
                target = standard
            mapping[variant] = target
            mapping[variant_lower] = target

        # 5. 最后合并中文映射（最高优先级）
        mapping.update(chinese_mapping)

        # 6. 为所有标准通路名称生成中文同义词（包括不在KEGG中的）
        all_standards = set(mapping.values())  # 所有不重复的标准名称
        for standard in all_standards:
            chinese_synonyms = self.generate_chinese_synonyms(standard)
            for cn in chinese_synonyms:
                if cn not in mapping:  # 不覆盖已有映射
                    mapping[cn] = standard
                    chinese_mapping[cn] = standard

        # 再次更新映射
        mapping.update(chinese_mapping)

        logger.info(f"映射表共 {len(mapping)} 条记录（含 {len(chinese_mapping)} 条中文同义词）")

        return mapping

    def save_mapping(self, mapping: Dict[str, str], filename: str = "pathway_mapping.json"):
        """
        保存映射表到文件

        Args:
            mapping: 映射字典
            filename: 输出文件名
        """
        output_path = os.path.join(self.output_dir, filename)

        # 按键排序
        sorted_mapping = dict(sorted(mapping.items()))

        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(sorted_mapping, f, ensure_ascii=False, indent=2)

        logger.info(f"映射表已保存到: {output_path}")

        return output_path

    def generate_python_dict(self, mapping: Dict[str, str], filename: str = "pathway_synonyms_generated.py"):
        """
        生成 Python 字典格式的映射表

        Args:
            mapping: 映射字典
            filename: 输出文件名
        """
        output_path = os.path.join(self.output_dir, filename)

        # 按值分组
        grouped = {}
        for synonym, standard in mapping.items():
            if standard not in grouped:
                grouped[standard] = []
            grouped[standard].append(synonym)

        with open(output_path, 'w', encoding='utf-8') as f:
            f.write('"""\n自动生成的通路名称映射表\n数据来源: KEGG, Reactome, WikiPathways\n生成时间: {}\n"""\n\n'.format(
                time.strftime("%Y-%m-%d %H:%M:%S")
            ))
            f.write("PATHWAY_SYNONYMS = {\n")

            for standard in sorted(grouped.keys()):
                synonyms = sorted(set(grouped[standard]))
                for syn in synonyms:
                    if syn != standard:  # 不写入自映射
                        f.write(f'    "{syn}": "{standard}",\n')

            f.write("}\n")

        logger.info(f"Python字典已保存到: {output_path}")

        return output_path

    def update_utils_file(self, mapping: Dict[str, str], utils_path: str = None):
        """
        更新 utils.py 中的 PATHWAY_SYNONYMS 字典

        Args:
            mapping: 映射字典
            utils_path: utils.py 文件路径
        """
        if utils_path is None:
            utils_path = os.path.join(self.output_dir, "utils.py")

        if not os.path.exists(utils_path):
            logger.warning(f"utils.py 文件不存在: {utils_path}")
            return None

        # 读取现有文件
        with open(utils_path, 'r', encoding='utf-8') as f:
            content = f.read()

        # 查找 PATHWAY_SYNONYMS 字典位置
        import re
        pattern = r'PATHWAY_SYNONYMS\s*=\s*\{[^}]+\}'
        match = re.search(pattern, content, re.DOTALL)

        if not match:
            logger.warning("未找到 PATHWAY_SYNONYMS 字典")
            return None

        # 生成新的字典内容
        dict_lines = ["PATHWAY_SYNONYMS = {"]
        for syn in sorted(mapping.keys()):
            standard = mapping[syn]
            if syn.lower() != standard.lower():  # 不包含自映射
                dict_lines.append(f'    "{syn}": "{standard}",')
        dict_lines.append("}")

        new_dict = "\n".join(dict_lines)

        # 替换
        new_content = content[:match.start()] + new_dict + content[match.end():]

        # 写回文件
        with open(utils_path, 'w', encoding='utf-8') as f:
            f.write(new_content)

        logger.info(f"已更新 utils.py 中的 PATHWAY_SYNONYMS ({len(mapping)} 条记录)")
        return utils_path


def main():
    """主函数"""
    import argparse

    parser = argparse.ArgumentParser(description="通路数据下载器")
    parser.add_argument("--update-utils", action="store_true", help="更新 utils.py 中的映射表")
    parser.add_argument("--output-dir", default=None, help="输出目录")
    args = parser.parse_args()

    print("=" * 60)
    print("通路数据下载器")
    print("=" * 60)
    print()

    downloader = PathwayDownloader(output_dir=args.output_dir)

    # 下载并构建映射表
    mapping = downloader.build_mapping_table()

    # 保存结果
    json_path = downloader.save_mapping(mapping)
    py_path = downloader.generate_python_dict(mapping)

    print()
    print("=" * 60)
    print("下载完成!")
    print(f"JSON文件: {json_path}")
    print(f"Python文件: {py_path}")
    print(f"映射条目: {len(mapping)}")

    # 可选：更新 utils.py
    if args.update_utils:
        utils_path = downloader.update_utils_file(mapping)
        if utils_path:
            print(f"已更新: {utils_path}")

    print("=" * 60)

    # 显示部分结果
    print("\n关键通路映射示例:")
    examples = ["p53信号通路", "MAPK信号通路", "Wnt信号通路", "凋亡通路", "细胞周期通路"]
    for key in examples:
        if key in mapping:
            print(f"  '{key}' -> '{mapping[key]}'")


if __name__ == "__main__":
    main()