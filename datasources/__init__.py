"""
数据源集成模块

提供对外部数据库的统一访问接口：
- Open Targets: 基因-疾病关联
- ClinVar: 变异临床意义
- COSMIC: 癌症体细胞突变
- KEGG/Reactome: 通路数据
"""

from .opentargets_client import OpenTargetsClient
from .clinvar_client import ClinVarClient
from .cosmic_client import COSMICClient
from .pathway_client import PathwayClient

__all__ = [
    'OpenTargetsClient',
    'ClinVarClient',
    'COSMICClient',
    'PathwayClient',
]
