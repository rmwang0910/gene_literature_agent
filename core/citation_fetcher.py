"""
引用数据查询模块 - 通过 CrossRef 和 NCBI 获取引用次数和期刊信息
"""
import time
import logging
from typing import List, Dict, Optional, Tuple, Any
from dataclasses import dataclass
import requests
import xml.etree.ElementTree as ET
from abc import ABC, abstractmethod

from ..config import default_config, get_api_url
from ..utils.journal_data import get_sjr, get_journal_info

# 配置日志
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class CitationInfo:
    """引用信息数据类"""
    pmid: str
    doi: str = ""
    citation_count: int = 0
    journal: str = ""
    impact_factor: float = 0.0
    year: int = 0
    is_open_access: bool = False
    source: str = ""  # 数据来源
    publication_type: str = ""  # 文献类型: "review", "research", ""(未知)


class CitationProvider(ABC):
    """引用数据提供者抽象基类"""

    @abstractmethod
    def get_citation_info(self, pmid: str, doi: str = "") -> CitationInfo:
        """获取引用信息"""
        pass

    @abstractmethod
    def get_citation_info_batch(self, pmids: List[str]) -> Dict[str, CitationInfo]:
        """批量获取引用信息"""
        pass


class CrossRefProvider(CitationProvider):
    """CrossRef API 提供者"""

    def __init__(self, email: Optional[str] = None):
        """
        初始化 CrossRef 提供者

        Args:
            email: 用于礼貌池（polite pool）的邮箱
        """
        self.email = email or default_config.ncbi_email or "research@example.com"
        self.headers = {
            "User-Agent": f"GeneLiteratureAgent/1.0 (mailto:{self.email})"
        }

    def get_citation_info(self, pmid: str, doi: str = "") -> CitationInfo:
        """
        获取单篇文章的引用信息

        Args:
            pmid: PubMed ID
            doi: DOI（可选，优先使用）

        Returns:
            CitationInfo 对象
        """
        work_data = None

        # 1. 优先通过 DOI 查询
        if doi:
            work_data = self._get_work_by_doi(doi)

        # 2. 通过 PMID 查询（如果 DOI 查询失败）
        if not work_data and pmid:
            work_data = self._get_work_by_pmid(pmid)

        if work_data:
            return self._parse_work_data(pmid, work_data)

        return CitationInfo(pmid=pmid, source="crossref")

    def _get_work_by_doi(self, doi: str) -> Optional[Dict]:
        """通过 DOI 获取文章数据"""
        # 清理 DOI
        doi = doi.replace("doi:", "").strip()
        url = f"{get_api_url('crossref')}/works/{requests.utils.quote(doi)}"

        try:
            response = requests.get(url, headers=self.headers, timeout=15)
            if response.status_code == 200:
                data = response.json()
                return data.get("message", {})
        except Exception as e:
            logger.debug(f"CrossRef DOI 查询失败 ({doi}): {e}")

        return None

    def _get_work_by_pmid(self, pmid: str) -> Optional[Dict]:
        """通过 PMID 获取文章数据"""
        url = f"{get_api_url('crossref')}/works?filter=pmid:{pmid}&rows=1"

        try:
            response = requests.get(url, headers=self.headers, timeout=15)
            if response.status_code == 200:
                data = response.json()
                items = data.get("message", {}).get("items", [])
                if items:
                    return items[0]
        except Exception as e:
            logger.debug(f"CrossRef PMID 查询失败 ({pmid}): {e}")

        return None

    def _parse_work_data(self, pmid: str, work_data: Dict) -> CitationInfo:
        """解析 CrossRef 工作数据"""
        # 引用次数
        citation_count = work_data.get("is-referenced-by-count", 0)

        # DOI
        doi = work_data.get("DOI", "")

        # 期刊
        journal = ""
        if "container-title" in work_data and work_data["container-title"]:
            journal = work_data["container-title"][0]
        elif "short-container-title" in work_data and work_data["short-container-title"]:
            journal = work_data["short-container-title"][0]

        # 年份
        year = 0
        published = work_data.get("published", {})
        if published and "date-parts" in published:
            date_parts = published["date-parts"]
            if date_parts and date_parts[0]:
                year = date_parts[0][0]

        # 影响因子
        impact_factor = get_sjr(journal)

        # 开放获取
        is_oa = work_data.get("is-oa", False)

        return CitationInfo(
            pmid=pmid,
            doi=doi,
            citation_count=citation_count,
            journal=journal,
            impact_factor=impact_factor,
            year=year,
            is_open_access=is_oa,
            source="crossref"
        )

    def get_citation_info_batch(
        self,
        pmids: List[str],
        delay: float = 0.1
    ) -> Dict[str, CitationInfo]:
        """批量获取引用信息"""
        results = {}

        for pmid in pmids:
            try:
                results[pmid] = self.get_citation_info(pmid)
                time.sleep(delay)
            except Exception as e:
                logger.warning(f"获取 PMID {pmid} 引用信息失败: {e}")
                results[pmid] = CitationInfo(pmid=pmid, source="crossref")

        return results


class OpenAlexProvider(CitationProvider):
    """OpenAlex API 提供者"""

    def __init__(self, email: Optional[str] = None):
        self.email = email or default_config.ncbi_email
        self.headers = {}
        if self.email:
            self.headers["User-Agent"] = f"mailto:{self.email}"

    def get_citation_info(self, pmid: str, doi: str = "") -> CitationInfo:
        """获取单篇文章的引用信息"""
        work_data = None

        # 1. 通过 DOI 查询
        if doi:
            doi_clean = doi.replace("doi:", "").strip()
            work_data = self._get_work(f"doi:{doi_clean}")

        # 2. 通过 PMID 查询
        if not work_data and pmid:
            work_data = self._get_work(f"pmid:{pmid}")

        if work_data:
            return self._parse_work_data(pmid, work_data)

        return CitationInfo(pmid=pmid, source="openalex")

    def _get_work(self, openalex_id: str) -> Optional[Dict]:
        """获取 OpenAlex 工作数据"""
        url = f"{get_api_url('openalex')}/works/{openalex_id}"

        try:
            response = requests.get(url, headers=self.headers, timeout=15)
            if response.status_code == 200:
                return response.json()
        except Exception as e:
            logger.debug(f"OpenAlex 查询失败 ({openalex_id}): {e}")

        return None

    def _parse_work_data(self, pmid: str, work_data: Dict) -> CitationInfo:
        """解析 OpenAlex 工作数据"""
        # 引用次数
        citation_count = work_data.get("cited_by_count", 0)

        # DOI
        doi = work_data.get("doi", "")

        # 期刊
        journal = ""
        source_info = {}
        if work_data.get("primary_location"):
            source_info = work_data["primary_location"].get("source", {})
        if source_info:
            journal = source_info.get("display_name", "")

        # 年份
        year = work_data.get("publication_year", 0)

        # 影响因子
        impact_factor = get_sjr(journal)

        # 开放获取
        is_oa = False
        oa_info = work_data.get("open_access", {})
        if oa_info:
            is_oa = oa_info.get("is_oa", False)

        return CitationInfo(
            pmid=pmid,
            doi=doi,
            citation_count=citation_count,
            journal=journal,
            impact_factor=impact_factor,
            year=year,
            is_open_access=is_oa,
            source="openalex"
        )

    def get_citation_info_batch(
        self,
        pmids: List[str],
        delay: float = 0.1
    ) -> Dict[str, CitationInfo]:
        """批量获取引用信息"""
        results = {}

        for pmid in pmids:
            try:
                results[pmid] = self.get_citation_info(pmid)
                time.sleep(delay)
            except Exception as e:
                logger.warning(f"获取 PMID {pmid} 引用信息失败: {e}")
                results[pmid] = CitationInfo(pmid=pmid, source="openalex")

        return results


class NCBICitationProvider(CitationProvider):
    """NCBI E-utilities API 提供者"""

    def __init__(
        self,
        email: Optional[str] = None,
        api_key: Optional[str] = None
    ):
        self.email = email or default_config.ncbi_email
        self.api_key = api_key or default_config.ncbi_api_key
        self.request_delay = default_config.request_delay

    def get_citation_info(self, pmid: str, doi: str = "") -> CitationInfo:
        """通过 NCBI 获取引用信息"""
        try:
            # 使用 elink 获取引用数
            params = {
                "dbfrom": "pubmed",
                "db": "pubmed",
                "linkname": "pubmed_pubmed_citedin",
                "id": pmid,
                "retmode": "json"
            }

            if self.email:
                params["email"] = self.email
            if self.api_key:
                params["api_key"] = self.api_key

            url = f"{get_api_url('ncbi')}/elink.fcgi"
            response = requests.get(url, params=params, timeout=15)

            citation_count = 0
            if response.status_code == 200:
                data = response.json()
                linksets = data.get("linksets", [])
                if linksets and "linksetdbs" in linksets[0]:
                    for linksetdb in linksets[0]["linksetdbs"]:
                        if "links" in linksetdb:
                            citation_count = len(linksetdb["links"])

            # 获取文章详细信息
            article_info = self._get_article_info(pmid)

            return CitationInfo(
                pmid=pmid,
                doi=article_info.get("doi", ""),
                citation_count=citation_count,
                journal=article_info.get("journal", ""),
                impact_factor=get_sjr(article_info.get("journal", "")),
                year=article_info.get("year", 0),
                source="ncbi",
                publication_type=article_info.get("publication_type", "")
            )

        except Exception as e:
            logger.warning(f"NCBI 查询失败 ({pmid}): {e}")
            return CitationInfo(pmid=pmid, source="ncbi")

    def _get_article_info(self, pmid: str) -> Dict:
        """获取文章基本信息"""
        try:
            params = {
                "db": "pubmed",
                "id": pmid,
                "rettype": "xml",
                "retmode": "xml"
            }

            if self.email:
                params["email"] = self.email
            if self.api_key:
                params["api_key"] = self.api_key

            url = f"{get_api_url('ncbi')}/efetch.fcgi"
            response = requests.get(url, params=params, timeout=15)

            if response.status_code == 200:
                root = ET.fromstring(response.text)

                # 期刊
                journal_elem = root.find(".//Journal/Title")
                journal = journal_elem.text if journal_elem is not None else ""

                # 年份
                year = 0
                year_elem = root.find(".//PubDate/Year")
                if year_elem is not None and year_elem.text:
                    year = int(year_elem.text)

                # DOI
                doi = ""
                for article_id in root.findall(".//ArticleId"):
                    if article_id.get("IdType") == "doi":
                        doi = article_id.text or ""
                        break

                # 文献类型
                publication_type = ""
                pub_type_elems = root.findall(".//PublicationType")
                pub_types = [pt.text.lower() if pt.text else "" for pt in pub_type_elems]

                # 判断是否为综述
                is_review = any(
                    "review" in pt or "meta-analysis" in pt or "systematic" in pt
                    for pt in pub_types
                )
                if is_review:
                    publication_type = "review"
                else:
                    # 判断是否为原创研究
                    is_research = any(
                        "research" in pt or "article" in pt or "clinical trial" in pt
                        for pt in pub_types
                    )
                    if is_research:
                        publication_type = "research"

                return {
                    "journal": journal,
                    "year": year,
                    "doi": doi,
                    "publication_type": publication_type
                }

        except Exception as e:
            logger.debug(f"NCBI 文章信息获取失败 ({pmid}): {e}")

        return {}

    def get_citation_info_batch(
        self,
        pmids: List[str],
        delay: float = None
    ) -> Dict[str, CitationInfo]:
        """批量获取引用信息"""
        delay = delay or self.request_delay
        results = {}

        for pmid in pmids:
            try:
                results[pmid] = self.get_citation_info(pmid)
                time.sleep(delay)
            except Exception as e:
                logger.warning(f"获取 PMID {pmid} 引用信息失败: {e}")
                results[pmid] = CitationInfo(pmid=pmid, source="ncbi")

        return results


class CitationEnricher:
    """引用数据丰富器 - 整合多个数据源"""

    def __init__(
        self,
        providers: Optional[List[CitationProvider]] = None,
        cache_dir: Optional[str] = None,
        use_cache: bool = True
    ):
        """
        初始化引用数据丰富器

        Args:
            providers: 数据提供者列表（按优先级排序）
            cache_dir: 缓存目录
            use_cache: 是否使用缓存
        """
        if providers:
            self.providers = providers
        else:
            # 默认提供者优先级：CrossRef > OpenAlex > NCBI
            self.providers = [
                CrossRefProvider(),
                OpenAlexProvider(),
                NCBICitationProvider()
            ]

        self.cache_dir = cache_dir or default_config.cache_dir
        self.use_cache = use_cache

        # 初始化缓存
        if self.use_cache:
            try:
                from ..utils.cache_manager import get_citation_cache
                self._cache = get_citation_cache()
            except ImportError:
                self._cache = None
                self.use_cache = False
        else:
            self._cache = None

    def enrich_pmid(self, pmid: str, doi: str = "", use_cache: bool = True) -> CitationInfo:
        """
        丰富单篇文章的引用信息

        Args:
            pmid: PubMed ID
            doi: DOI（可选）
            use_cache: 是否使用缓存

        Returns:
            CitationInfo 对象
        """
        # 尝试从缓存获取
        if use_cache and self._cache:
            cached = self._cache.get(pmid)
            if cached:
                logger.debug(f"从缓存获取引用数据: {pmid}")
                return CitationInfo(**cached)

        for provider in self.providers:
            try:
                info = provider.get_citation_info(pmid, doi)
                # 如果获取到了有效数据，返回
                if info.citation_count > 0 or info.journal:
                    # 存入缓存
                    if self._cache and use_cache:
                        self._cache.set(pmid, {
                            'pmid': info.pmid,
                            'doi': info.doi,
                            'citation_count': info.citation_count,
                            'journal': info.journal,
                            'impact_factor': info.impact_factor,
                            'year': info.year,
                            'publication_type': info.publication_type,
                            'source': info.source
                        })
                    return info
            except Exception as e:
                logger.debug(f"{provider.__class__.__name__} 查询失败: {e}")
                continue

        # 所有提供者都失败，返回空信息
        return CitationInfo(pmid=pmid, source="none")

    def enrich_pmids(
        self,
        pmids: List[str],
        doi_map: Optional[Dict[str, str]] = None,
        delay: float = 0.15
    ) -> Dict[str, CitationInfo]:
        """
        批量丰富 PMID 列表的引用信息

        Args:
            pmids: PMID 列表
            doi_map: PMID 到 DOI 的映射
            delay: 请求间隔

        Returns:
            PMID 到 CitationInfo 的映射
        """
        doi_map = doi_map or {}
        results = {}

        # 获取缓存的PMID
        cached_pmids = set()
        if self._cache:
            cached_data = self._cache.get_batch(pmids)
            for pmid, data in cached_data.items():
                results[pmid] = CitationInfo(**data)
                cached_pmids.add(pmid)

        # 获取未缓存的PMID
        missing_pmids = [p for p in pmids if p not in cached_pmids]

        if missing_pmids:
            logger.info(f"从API获取 {len(missing_pmids)} 个PMID的引用数据 (缓存命中: {len(cached_pmids)})")

        for i, pmid in enumerate(missing_pmids):
            if i > 0:
                time.sleep(delay)

            doi = doi_map.get(pmid, "")
            results[pmid] = self.enrich_pmid(pmid, doi, use_cache=True)

            # 日志
            if (i + 1) % 10 == 0:
                logger.info(f"已处理 {i + 1}/{len(missing_pmids)} 个 PMID")

        return results


def get_citation_stars(
    citation_count: int,
    impact_factor: float = 0.0,
    year: int = 0,
    current_year: int = None
) -> str:
    """
    根据引用次数、影响因子和年份生成星级

    引入年均引用次数，使较新文献在相同引用次数下获得更高评分

    Args:
        citation_count: 引用次数
        impact_factor: 影响因子
        year: 发表年份
        current_year: 当前年份

    Returns:
        星级字符串（如 "★★★☆☆"）
    """
    import datetime
    import math

    if current_year is None:
        current_year = datetime.datetime.now().year

    # 计算年均引用次数
    if year > 0 and year <= current_year:
        years_since_pub = max(1, current_year - year)  # 至少1年
        annual_citations = citation_count / years_since_pub
    else:
        annual_citations = citation_count  # 未知年份，使用原始引用次数
        years_since_pub = 1

    # 综合评分：年均引用次数权重 40%，原始引用次数 20%，影响因子权重 40%
    # 年均引用分档：0-5=1星, 5-15=2星, 15-30=3星, 30-50=4星, 50+=5星
    # 原始引用分档：0-10=1星, 11-50=2星, 51-100=3星, 101-500=4星, 500+=5星
    # 影响因子分档：0-5=1星, 5-10=2星, 10-20=3星, 20-50=4星, 50+=5星

    def get_stars_for_annual_citations(annual: float) -> int:
        if annual >= 50:
            return 5
        elif annual >= 30:
            return 4
        elif annual >= 15:
            return 3
        elif annual >= 5:
            return 2
        else:
            return 1

    def get_stars_for_citations(citations: int) -> int:
        if citations >= 500:
            return 5
        elif citations >= 100:
            return 4
        elif citations >= 50:
            return 3
        elif citations >= 10:
            return 2
        else:
            return 1

    def get_stars_for_if(if_value: float) -> int:
        if if_value >= 50:
            return 5
        elif if_value >= 20:
            return 4
        elif if_value >= 10:
            return 3
        elif if_value >= 5:
            return 2
        else:
            return 1

    annual_stars = get_stars_for_annual_citations(annual_citations)
    citation_stars = get_stars_for_citations(citation_count)
    if_stars = get_stars_for_if(impact_factor)

    # 加权平均：年均引用 40%, 总引用 20%, 影响因子 40%
    total_stars = annual_stars * 0.4 + citation_stars * 0.2 + if_stars * 0.4
    final_stars = round(total_stars)

    # 返回星级字符串
    return "★" * final_stars + "☆" * (5 - final_stars)


def get_normalized_score(
    citation_count: int,
    impact_factor: float = 0.0,
    year: int = 0,
    current_year: int = None
) -> float:
    """
    计算归一化影响力评分 (0-100)

    综合考虑年均引用次数和影响因子

    Args:
        citation_count: 引用次数
        impact_factor: 影响因子
        year: 发表年份
        current_year: 当前年份

    Returns:
        归一化评分 (0-100)
    """
    import datetime
    import math

    if current_year is None:
        current_year = datetime.datetime.now().year

    # 计算年均引用次数
    if year > 0 and year <= current_year:
        years_since_pub = max(1, current_year - year)
        annual_citations = citation_count / years_since_pub
    else:
        annual_citations = citation_count
        years_since_pub = 1

    # 年均引用分数 (0-50分)，使用对数归一化
    if annual_citations > 0:
        annual_score = min(50, 15 * math.log10(annual_citations + 1))
    else:
        annual_score = 0

    # 影响因子分数 (0-50分)
    if impact_factor >= 50:
        if_score = 50
    elif impact_factor >= 20:
        if_score = 40
    elif impact_factor >= 10:
        if_score = 30
    elif impact_factor >= 5:
        if_score = 20
    elif impact_factor > 0:
        if_score = 10
    else:
        if_score = 5

    return round(annual_score + if_score, 1)


def get_impact_badge(impact_factor: float) -> str:
    """
    根据影响因子生成标识

    Args:
        impact_factor: 影响因子

    Returns:
        影响因子标识
    """
    if impact_factor >= 50:
        return "🏆 顶级期刊"
    elif impact_factor >= 20:
        return "⭐ 高影响力"
    elif impact_factor >= 10:
        return "📈 中高影响力"
    elif impact_factor >= 5:
        return "📄 中等影响力"
    else:
        return ""


def calculate_confidence_score(
    citation_count: int,
    impact_factor: float,
    publication_year: int,
    current_year: int = None,
    llm_confidence: str = "medium"
) -> Dict[str, Any]:
    """
    计算综合置信度分数

    考虑因素：
    1. 发表年份权重（近5年权重更高）
    2. 引用次数权重
    3. 期刊影响因子权重
    4. LLM 初始置信度

    Args:
        citation_count: 引用次数
        impact_factor: 影响因子
        publication_year: 发表年份
        current_year: 当前年份（默认使用系统时间）
        llm_confidence: LLM 初始置信度 (high/medium/low)

    Returns:
        包含分数和等级的字典
    """
    import datetime

    if current_year is None:
        current_year = datetime.datetime.now().year

    # 1. 年份权重（近5年权重更高）
    year_weight = 1.0
    if publication_year > 0:
        years_old = current_year - publication_year
        if years_old <= 1:
            year_weight = 1.3  # 1年内：30%加成
        elif years_old <= 3:
            year_weight = 1.2  # 1-3年：20%加成
        elif years_old <= 5:
            year_weight = 1.1  # 3-5年：10%加成
        elif years_old <= 10:
            year_weight = 1.0  # 5-10年：正常权重
        else:
            year_weight = 0.9  # 10年以上：10%降权

    # 2. 引用次数分数 (0-100)
    # 使用对数归一化，避免极端值影响
    import math
    if citation_count > 0:
        citation_score = min(100, 20 * math.log10(citation_count + 1))
    else:
        citation_score = 0

    # 3. 影响因子分数 (0-100)
    if impact_factor >= 50:
        if_score = 100
    elif impact_factor >= 20:
        if_score = 80
    elif impact_factor >= 10:
        if_score = 60
    elif impact_factor >= 5:
        if_score = 40
    elif impact_factor > 0:
        if_score = 20
    else:
        if_score = 10  # 未知期刊给基础分

    # 4. LLM 置信度分数
    llm_score_map = {"high": 100, "medium": 70, "low": 40}
    llm_score = llm_score_map.get(llm_confidence, 70)

    # 5. 综合计算（加权平均）
    # 权重：引用 25%, 影响因子 25%, 年份 20%, LLM 30%
    raw_score = (
        citation_score * 0.25 +
        if_score * 0.25 +
        llm_score * 0.30 +
        (year_weight * 70) * 0.20  # 年份基准分70，乘以权重系数
    )

    # 归一化到 0-100
    final_score = min(100, max(0, raw_score))

    # 6. 置信度等级
    if final_score >= 80:
        level = "high"
        level_desc = "高置信度"
    elif final_score >= 50:
        level = "medium"
        level_desc = "中等置信度"
    else:
        level = "low"
        level_desc = "低置信度"

    return {
        "score": round(final_score, 1),
        "level": level,
        "level_desc": level_desc,
        "year_weight": year_weight,
        "components": {
            "citation_score": round(citation_score, 1),
            "impact_factor_score": round(if_score, 1),
            "llm_score": llm_score,
            "year_score": round(year_weight * 70, 1)
        }
    }


def get_confidence_stars(score: float) -> str:
    """
    根据置信度分数生成星级

    Args:
        score: 置信度分数 (0-100)

    Returns:
        星级字符串
    """
    if score >= 80:
        return "★★★★★"
    elif score >= 60:
        return "★★★★☆"
    elif score >= 40:
        return "★★★☆☆"
    elif score >= 20:
        return "★★☆☆☆"
    else:
        return "★☆☆☆☆"


if __name__ == "__main__":
    # 测试
    enricher = CitationEnricher()

    # 测试单个 PMID
    pmid = "33208945"  # Nature 文章
    info = enricher.enrich_pmid(pmid)
    print(f"PMID: {info.pmid}")
    print(f"期刊: {info.journal}")
    print(f"影响因子: {info.impact_factor}")
    print(f"引用次数: {info.citation_count}")
    print(f"星级: {get_citation_stars(info.citation_count, info.impact_factor)}")
    print(f"标识: {get_impact_badge(info.impact_factor)}")