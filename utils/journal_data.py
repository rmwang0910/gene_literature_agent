"""
期刊影响因子模块 - 加载 SJR 数据
"""
import csv
import logging
from pathlib import Path
from typing import Dict, Optional

logger = logging.getLogger(__name__)

# 数据文件路径
DATA_FILE = Path(__file__).parent.parent / "data" / "journal" / "scimagojr_2024.csv"

# 全局缓存
_journal_data: Optional[Dict[str, dict]] = None


def _parse_sjr(value: str) -> float:
    """解析 SJR 值（逗号作为小数点）"""
    if not value:
        return 0.0
    try:
        # 欧洲格式：逗号作为小数点
        return float(value.replace(",", "."))
    except ValueError:
        return 0.0


def load_journal_data() -> Dict[str, dict]:
    """
    加载期刊数据

    CSV 列顺序：
    Rank(0);Sourceid(1);Title(2);Type(3);Issn(4);Publisher(5);Open Access(6);
    Open Access Diamond(7);SJR(8);SJR Best Quartile(9);H index(10);...

    Returns:
        {期刊名: {sjr, h_index, quartile, ...}} 字典
    """
    global _journal_data

    if _journal_data is not None:
        return _journal_data

    _journal_data = {}

    if not DATA_FILE.exists():
        logger.warning(f"期刊数据文件不存在: {DATA_FILE}")
        return _journal_data

    try:
        with open(DATA_FILE, "r", encoding="utf-8") as f:
            reader = csv.reader(f, delimiter=";")

            # 跳过标题行
            next(reader)

            for row in reader:
                if len(row) < 11:
                    continue

                # 提取字段（按正确的列索引）
                title = row[2].strip('"')    # Title (列2)
                sjr = _parse_sjr(row[8])      # SJR (列8)
                quartile = row[9]             # SJR Best Quartile (列9)
                h_index_str = row[10]         # H index (列10)

                if not title:
                    continue

                try:
                    h_index = int(h_index_str) if h_index_str else 0
                except ValueError:
                    h_index = 0

                # 存储数据
                _journal_data[title.lower()] = {
                    "title": title,
                    "sjr": sjr,
                    "quartile": quartile,
                    "h_index": h_index
                }

        logger.info(f"加载了 {len(_journal_data)} 种期刊数据")

    except Exception as e:
        logger.error(f"加载期刊数据失败: {e}")

    return _journal_data


def get_sjr(journal_name: str) -> float:
    """
    获取期刊的 SJR 值

    Args:
        journal_name: 期刊名称

    Returns:
        SJR 值，未找到返回 0.0
    """
    data = load_journal_data()
    journal_key = journal_name.lower().strip()

    if journal_key in data:
        return data[journal_key]["sjr"]

    return 0.0


def get_journal_info(journal_name: str) -> Optional[dict]:
    """
    获取期刊完整信息

    Args:
        journal_name: 期刊名称

    Returns:
        期刊信息字典，未找到返回 None
    """
    data = load_journal_data()
    journal_key = journal_name.lower().strip()

    return data.get(journal_key)


def get_impact_factor(journal_name: str) -> float:
    """
    获取期刊影响因子（使用 SJR 作为替代指标）

    SJR 和 IF 的换算参考：
    - SJR > 10 ≈ IF > 20 (顶级期刊)
    - SJR > 5 ≈ IF > 10 (高影响力)
    - SJR > 2 ≈ IF > 5 (中等影响力)

    Args:
        journal_name: 期刊名称

    Returns:
        估算的影响因子值
    """
    sjr = get_sjr(journal_name)

    if sjr <= 0:
        return 0.0

    # SJR 到 IF 的近似换算
    # 这是一个粗略估算，SJR 通常低于 IF
    # 参考：https://www.scimagojr.com/help.php
    return sjr


# 预加载常用期刊（用于快速匹配）
_COMMON_JOURNALS = [
    "nature", "science", "cell", "nature medicine", "nature genetics",
    "nature biotechnology", "cancer cell", "cancer discovery",
    "journal of clinical oncology", "new england journal of medicine",
    "the lancet", "lancet oncology", "jama", "bmj"
]


def preload_common_journals():
    """预加载常用期刊数据"""
    load_journal_data()


if __name__ == "__main__":
    # 测试
    data = load_journal_data()
    print(f"共加载 {len(data)} 种期刊")

    # 测试查找
    test_journals = ["Nature", "Cell", "Science", "Cancer Research", "Unknown Journal"]
    for j in test_journals:
        info = get_journal_info(j)
        if info:
            print(f"{j}: SJR={info['sjr']}, Q={info['quartile']}, H={info['h_index']}")
        else:
            print(f"{j}: 未找到")
