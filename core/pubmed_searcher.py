"""
PubMed检索模块 - 使用 NCBI E-utilities API 检索文献
支持多种筛选条件：日期排序、引用次数、文章类型、期刊质量
"""
import time
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass, field
from datetime import datetime
import xml.etree.ElementTree as ET
from pathlib import Path
import json
import requests
import urllib.parse

from ..config import default_config, get_api_url
from ..utils.journal_data import get_sjr


@dataclass
class PubMedArticle:
    """PubMed 文章数据类"""
    pmid: str
    title: str
    abstract: str = ""
    authors: List[str] = field(default_factory=list)
    journal: str = ""
    publication_date: str = ""
    year: int = 0
    doi: str = ""
    keywords: List[str] = field(default_factory=list)
    mesh_terms: List[str] = field(default_factory=list)
    citation_count: int = 0
    impact_factor: float = 0.0
    is_open_access: bool = False
    article_type: str = ""


@dataclass
class SearchResult:
    """检索结果数据类"""
    query: str
    total_count: int
    articles: List[PubMedArticle]
    search_time: str = field(default_factory=lambda: datetime.now().isoformat())
    filters_applied: Dict = field(default_factory=dict)


# 文章类型映射
ARTICLE_TYPE_MAP = {
    "review": "Review",
    "clinical_trial": "Clinical Trial",
    "randomized_controlled_trial": "Randomized Controlled Trial",
    "meta_analysis": "Meta-Analysis",
    "case_report": "Case Report",
    "observational_study": "Observational Study",
    "systematic_review": "Systematic Review",
    "editorial": "Editorial",
    "letter": "Letter",
    "comment": "Comment",
}


class OpenAlexEnricher:
    """OpenAlex 数据丰富器 - 获取引用次数和期刊信息"""

    def __init__(self, email: Optional[str] = None):
        self.email = email or default_config.ncbi_email
        self.headers = {}
        if self.email:
            self.headers["User-Agent"] = f"mailto:{self.email}"

    def enrich_article(self, article: PubMedArticle) -> PubMedArticle:
        """丰富单篇文章的引用数据"""
        work_data = None

        # 1. 优先通过 DOI 查询
        if article.doi:
            doi = article.doi.replace("doi:", "").strip()
            work_data = self._get_work_by_id(f"doi:{doi}")

        # 2. 通过 PMID 查询
        if not work_data and article.pmid:
            work_data = self._get_work_by_id(f"pmid:{article.pmid}")

        if work_data:
            # 更新引用次数
            if 'cited_by_count' in work_data:
                article.citation_count = work_data['cited_by_count']

            # 更新开放获取状态
            oa = work_data.get('open_access', {})
            if oa.get('is_oa'):
                article.is_open_access = True

            # 更新期刊影响因子（从预设列表获取）
            source = {}
            if work_data.get('primary_location'):
                source = work_data['primary_location'].get('source', {})
            if source and source.get('display_name'):
                journal_name = source['display_name']
                article.impact_factor = get_sjr(journal_name)

        return article

    def _get_work_by_id(self, openalex_id: str) -> Optional[Dict]:
        """获取 OpenAlex 工作数据"""
        url = f"{get_api_url('openalex')}/works/{openalex_id}"
        try:
            response = requests.get(url, headers=self.headers, timeout=10)
            if response.status_code == 200:
                return response.json()
        except Exception:
            pass
        return None

    def enrich_articles(self, articles: List[PubMedArticle]) -> List[PubMedArticle]:
        """批量丰富文章数据"""
        for article in articles:
            try:
                self.enrich_article(article)
                time.sleep(0.1)  # 避免 API 限速
            except Exception:
                pass
        return articles


class PubMedSearcher:
    """PubMed 检索器类"""

    def __init__(
        self,
        email: Optional[str] = None,
        api_key: Optional[str] = None,
        cache_dir: Optional[str] = None
    ):
        """
        初始化 PubMed 检索器

        Args:
            email: NCBI 账号邮箱
            api_key: NCBI API Key
            cache_dir: 缓存目录
        """
        self.email = email or default_config.ncbi_email or "your_email@example.com"
        self.api_key = api_key or default_config.ncbi_api_key

        self.cache_dir = Path(cache_dir or default_config.cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.request_delay = default_config.request_delay
        self._last_request_time = 0

        # OpenAlex 丰富器
        self._enricher = None

    @property
    def enricher(self):
        """懒加载 OpenAlex 丰富器"""
        if self._enricher is None:
            self._enricher = OpenAlexEnricher(self.email)
        return self._enricher

    def _rate_limit(self):
        """限制请求频率"""
        elapsed = time.time() - self._last_request_time
        if elapsed < self.request_delay:
            time.sleep(self.request_delay - elapsed)
        self._last_request_time = time.time()

    def _make_request(self, endpoint: str, params: Dict) -> str:
        """发送 HTTP 请求"""
        self._rate_limit()

        params["email"] = self.email
        if self.api_key:
            params["api_key"] = self.api_key

        url = f"{get_api_url('ncbi')}/{endpoint}"
        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()

        return response.text

    def build_query(
        self,
        gene_symbol: str,
        aliases: Optional[List[str]] = None,
        background: Optional[str] = None,
        date_range: Optional[Tuple[int, int]] = None,
        article_types: Optional[List[str]] = None
    ) -> str:
        """
        构建 PubMed 检索式

        Args:
            gene_symbol: 基因符号
            aliases: 基因别名列表
            background: 背景（疾病/组织）
            date_range: 日期范围 (起始年, 结束年)
            article_types: 文章类型限制

        Returns:
            检索式字符串
        """
        # 基因相关词
        gene_terms = [f'{gene_symbol}[Title/Abstract]']
        if aliases:
            for alias in aliases[:3]:
                if alias and not any(c in alias for c in '()[]'):
                    gene_terms.append(f'{alias}[Title/Abstract]')

        # 用括号包围OR条件
        if len(gene_terms) > 1:
            gene_query = "(" + " OR ".join(gene_terms[:4]) + ")"
        else:
            gene_query = gene_terms[0]

        # 组合检索式
        if background:
            query = f"({gene_query}) AND {background}[Title/Abstract]"
        else:
            query = gene_query

        # 添加日期限制
        if date_range:
            start_year, end_year = date_range
            query += f' AND ({start_year}[PDAT] : {end_year}[PDAT])'

        # 添加文章类型限制
        if article_types:
            type_filters = []
            for t in article_types:
                mapped = ARTICLE_TYPE_MAP.get(t.lower(), t)
                type_filters.append(f'{mapped}[PTYP]')
            if type_filters:
                query += f" AND ({' OR '.join(type_filters)})"

        return query

    def search(
        self,
        query: str,
        max_results: int = 10,
        sort: str = "relevance"
    ) -> List[str]:
        """
        执行检索，返回 PMID 列表

        Args:
            query: 检索式
            max_results: 最大结果数
            sort: 排序方式 (relevance, pub_date, first_author, last_author, journal)

        Returns:
            PMID 列表
        """
        try:
            params = {
                "db": "pubmed",
                "term": query,
                "retmax": max_results,
                "sort": sort,
                "retmode": "json"
            }

            response_text = self._make_request("esearch.fcgi", params)
            data = json.loads(response_text)

            return data.get("esearchresult", {}).get("idlist", [])

        except Exception as e:
            print(f"PubMed 检索失败: {e}")
            return []

    def fetch_details(self, pmids: List[str]) -> List[PubMedArticle]:
        """获取文章详细信息"""
        if not pmids:
            return []

        articles = []

        # 分批获取（每批最多200篇）
        batch_size = 200
        for i in range(0, len(pmids), batch_size):
            batch = pmids[i:i + batch_size]

            try:
                params = {
                    "db": "pubmed",
                    "id": ",".join(batch),
                    "rettype": "xml",
                    "retmode": "xml"
                }

                xml_data = self._make_request("efetch.fcgi", params)
                root = ET.fromstring(xml_data)

                for article_elem in root.findall(".//PubmedArticle"):
                    article = self._parse_article(article_elem)
                    if article:
                        articles.append(article)

            except Exception as e:
                print(f"获取文章详情失败: {e}")
                continue

        return articles

    def _parse_article(self, article_elem) -> Optional[PubMedArticle]:
        """解析单篇文章"""
        try:
            # PMID
            pmid_elem = article_elem.find(".//PMID")
            pmid = pmid_elem.text if pmid_elem is not None else ""

            # Title
            title_elem = article_elem.find(".//ArticleTitle")
            title = "".join(title_elem.itertext()) if title_elem is not None else ""

            # Abstract
            abstract_parts = []
            for abstract_elem in article_elem.findall(".//AbstractText"):
                label = abstract_elem.get("Label", "")
                text = "".join(abstract_elem.itertext())
                if label:
                    abstract_parts.append(f"{label}: {text}")
                else:
                    abstract_parts.append(text)
            abstract = " ".join(abstract_parts)

            # Authors
            authors = []
            for author_elem in article_elem.findall(".//Author"):
                lastname = author_elem.findtext("LastName", "")
                forename = author_elem.findtext("ForeName", "")
                if lastname:
                    authors.append(f"{lastname} {forename}".strip())

            # Journal
            journal_elem = article_elem.find(".//Journal/Title")
            journal = journal_elem.text if journal_elem is not None else ""

            # Publication Date
            pub_date = ""
            year = 0
            year_elem = article_elem.find(".//PubDate/Year")
            if year_elem is not None and year_elem.text:
                year = int(year_elem.text)
                pub_date = year_elem.text

            medline_date = article_elem.find(".//PubDate/MedlineDate")
            if medline_date is not None and medline_date.text:
                pub_date = medline_date.text
                try:
                    year = int(medline_date.text[:4])
                except ValueError:
                    pass

            # DOI
            doi = ""
            for article_id in article_elem.findall(".//ArticleId"):
                if article_id.get("IdType") == "doi":
                    doi = article_id.text or ""
                    break

            # Keywords
            keywords = []
            for keyword_elem in article_elem.findall(".//Keyword"):
                if keyword_elem.text:
                    keywords.append(keyword_elem.text)

            # MeSH Terms
            mesh_terms = []
            for mesh_elem in article_elem.findall(".//MeshHeading/DescriptorName"):
                if mesh_elem.text:
                    mesh_terms.append(mesh_elem.text)

            # Article Type
            article_type = ""
            pub_type_elem = article_elem.find(".//PublicationType")
            if pub_type_elem is not None and pub_type_elem.text:
                article_type = pub_type_elem.text

            # 从数据库获取影响因子(SJR)
            impact_factor = get_sjr(journal)

            return PubMedArticle(
                pmid=pmid,
                title=title,
                abstract=abstract,
                authors=authors,
                journal=journal,
                publication_date=pub_date,
                year=year,
                doi=doi,
                keywords=keywords,
                mesh_terms=mesh_terms,
                impact_factor=impact_factor,
                article_type=article_type
            )

        except Exception as e:
            print(f"解析文章失败: {e}")
            return None

    def search_and_fetch(
        self,
        query: str,
        max_results: int = 10,
        sort: str = "relevance",
        min_citations: int = 0,
        min_impact_factor: float = 0.0,
        open_access_only: bool = False,
        enrich: bool = True
    ) -> SearchResult:
        """
        检索并获取详细信息

        Args:
            query: 检索式
            max_results: 最大结果数
            sort: 排序方式
            min_citations: 最小引用次数
            min_impact_factor: 最小影响因子
            open_access_only: 仅开放获取
            enrich: 是否用 OpenAlex 丰富数据

        Returns:
            SearchResult 对象
        """
        pmids = self.search(query, max_results, sort)
        articles = self.fetch_details(pmids)

        # 记录应用的筛选条件
        filters_applied = {
            "sort": sort,
            "min_citations": min_citations,
            "min_impact_factor": min_impact_factor,
            "open_access_only": open_access_only
        }

        # 用 OpenAlex 丰富数据（获取引用次数）
        if enrich and (min_citations > 0 or open_access_only):
            print(f"  正在获取引用数据...")
            articles = self.enricher.enrich_articles(articles)

        # 应用筛选条件
        filtered_articles = []
        for article in articles:
            # 引用次数筛选
            if min_citations > 0 and article.citation_count < min_citations:
                continue
            # 影响因子筛选
            if min_impact_factor > 0 and article.impact_factor < min_impact_factor:
                continue
            # 开放获取筛选
            if open_access_only and not article.is_open_access:
                continue
            filtered_articles.append(article)

        # 获取总结果数
        try:
            params = {
                "db": "pubmed",
                "term": query,
                "retmax": 0,
                "retmode": "json"
            }
            response_text = self._make_request("esearch.fcgi", params)
            data = json.loads(response_text)
            total_count = int(data.get("esearchresult", {}).get("count", 0))
        except Exception:
            total_count = len(filtered_articles)

        return SearchResult(
            query=query,
            total_count=total_count,
            articles=filtered_articles,
            filters_applied=filters_applied
        )


def search_pubmed(
    query: str,
    max_results: int = 10,
    sort: str = "relevance",
    min_citations: int = 0,
    min_impact_factor: float = 0.0,
    article_types: Optional[List[str]] = None
) -> List[PubMedArticle]:
    """
    简单的检索函数

    Args:
        query: 检索式
        max_results: 最大结果数
        sort: 排序方式
        min_citations: 最小引用次数
        min_impact_factor: 最小影响因子
        article_types: 文章类型列表

    Returns:
        文章列表
    """
    searcher = PubMedSearcher()

    # 如果有文章类型限制，添加到查询中
    if article_types:
        type_filters = []
        for t in article_types:
            mapped = ARTICLE_TYPE_MAP.get(t.lower(), t)
            type_filters.append(f'{mapped}[PTYP]')
        if type_filters:
            query += f" AND ({' OR '.join(type_filters)})"

    result = searcher.search_and_fetch(
        query,
        max_results,
        sort,
        min_citations,
        min_impact_factor
    )
    return result.articles


if __name__ == "__main__":
    # 测试
    searcher = PubMedSearcher()

    # 构建检索式
    query = searcher.build_query(
        gene_symbol="BRAF",
        aliases=["BRAF1"],
        background="lung cancer",
        article_types=["review"]
    )
    print(f"检索式: {query}")

    # 检索（按日期排序，筛选高引用）
    result = searcher.search_and_fetch(
        query,
        max_results=10,
        sort="pub_date",
        min_citations=10,
        enrich=True
    )
    print(f"找到 {result.total_count} 篇文献，获取 {len(result.articles)} 篇")
    print(f"筛选条件: {result.filters_applied}")

    for article in result.articles:
        print(f"\nPMID: {article.pmid}")
        print(f"标题: {article.title}")
        print(f"期刊: {article.journal} (IF: {article.impact_factor})")
        print(f"引用: {article.citation_count}")
        print(f"年份: {article.year}")