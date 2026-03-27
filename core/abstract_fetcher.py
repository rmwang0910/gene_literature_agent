"""
摘要获取与预处理模块
"""
import re
from typing import List, Dict, Optional
from dataclasses import dataclass

from .pubmed_searcher import PubMedArticle, PubMedSearcher


@dataclass
class ProcessedAbstract:
    """处理后的摘要"""
    pmid: str
    title: str
    abstract: str
    word_count: int
    sections: Dict[str, str]  # 分段摘要（如 Background, Methods, Results 等）


class AbstractProcessor:
    """摘要处理器"""

    # 常见摘要分段标识
    SECTION_PATTERNS = {
        "background": r"(?i)(background|objective|aim|introduction)[:：]?\s*",
        "methods": r"(?i)(methods|methodology|materials and methods|design)[:：]?\s*",
        "results": r"(?i)(results|findings)[:：]?\s*",
        "conclusions": r"(?i)(conclusions?|summary|interpretation)[:：]?\s*",
        "significance": r"(?i)(significance|clinical relevance)[:：]?\s*"
    }

    def __init__(self, max_length: int = 4000):
        """
        初始化摘要处理器

        Args:
            max_length: 摘要最大长度（字符数），用于限制 LLM 输入
        """
        self.max_length = max_length

    def clean_text(self, text: str) -> str:
        """
        清洗文本

        Args:
            text: 原始文本

        Returns:
            清洗后的文本
        """
        if not text:
            return ""

        # 移除多余的空白字符
        text = re.sub(r'\s+', ' ', text)

        # 移除特殊字符（保留基本标点）
        text = re.sub(r'[^\w\s\.,;:!?()\[\]\-"\']', '', text)

        # 统一标点
        text = text.replace('。', '.').replace('，', ',').replace('：', ':')
        text = text.replace('；', ';').replace('！', '!').replace('？', '?')

        return text.strip()

    def extract_sections(self, abstract: str) -> Dict[str, str]:
        """
        提取摘要的结构化分段

        Args:
            abstract: 摘要文本

        Returns:
            分段字典
        """
        sections = {}

        # 尝试匹配分段标题
        remaining = abstract
        for section_name, pattern in self.SECTION_PATTERNS.items():
            match = re.search(pattern, remaining)
            if match:
                # 找到该段落的起始位置
                start = match.end()
                # 查找下一个分段标题
                next_match = None
                for other_pattern in self.SECTION_PATTERNS.values():
                    m = re.search(other_pattern, remaining[start:])
                    if m:
                        if next_match is None or m.start() < next_match.start():
                            next_match = m

                if next_match:
                    end = start + next_match.start()
                    sections[section_name] = remaining[start:end].strip()
                    remaining = remaining[start:]
                else:
                    sections[section_name] = remaining[start:].strip()
                    remaining = ""

        # 如果没有找到分段，返回整个摘要
        if not sections and abstract:
            sections["full"] = abstract

        return sections

    def truncate(self, text: str, max_length: Optional[int] = None) -> str:
        """
        截断文本到指定长度

        Args:
            text: 原始文本
            max_length: 最大长度

        Returns:
            截断后的文本
        """
        max_length = max_length or self.max_length
        if len(text) <= max_length:
            return text

        # 尝试在句子边界截断
        truncated = text[:max_length]
        last_period = truncated.rfind('.')
        last_question = truncated.rfind('?')
        last_exclaim = truncated.rfind('!')

        last_sentence_end = max(last_period, last_question, last_exclaim)

        if last_sentence_end > max_length * 0.7:
            return truncated[:last_sentence_end + 1]

        return truncated + "..."

    def process(self, article: PubMedArticle) -> ProcessedAbstract:
        """
        处理单篇文章

        Args:
            article: PubMed 文章

        Returns:
            ProcessedAbstract 对象
        """
        # 合并标题和摘要
        full_text = f"{article.title}. {article.abstract}"

        # 清洗文本
        cleaned = self.clean_text(full_text)

        # 提取分段
        sections = self.extract_sections(cleaned)

        # 截断
        truncated = self.truncate(cleaned)

        return ProcessedAbstract(
            pmid=article.pmid,
            title=article.title,
            abstract=truncated,
            word_count=len(truncated.split()),
            sections=sections
        )

    def process_batch(self, articles: List[PubMedArticle]) -> List[ProcessedAbstract]:
        """
        批量处理文章

        Args:
            articles: 文章列表

        Returns:
            处理后的摘要列表
        """
        return [self.process(article) for article in articles]

    def prepare_for_llm(
        self,
        article: PubMedArticle,
        include_title: bool = True,
        max_length: Optional[int] = None
    ) -> str:
        """
        准备 LLM 输入文本

        Args:
            article: 文章
            include_title: 是否包含标题
            max_length: 最大长度

        Returns:
            格式化的文本
        """
        parts = []

        if include_title:
            parts.append(f"Title: {article.title}")

        if article.journal and article.year:
            parts.append(f"Source: {article.journal}, {article.year}")

        if article.abstract:
            parts.append(f"Abstract: {article.abstract}")

        full_text = "\n".join(parts)
        return self.truncate(full_text, max_length)


class AbstractFetcher:
    """摘要获取器（用于单独获取摘要）"""

    def __init__(self, email: Optional[str] = None, api_key: Optional[str] = None):
        """
        初始化摘要获取器

        Args:
            email: NCBI 邮箱
            api_key: NCBI API Key
        """
        self.searcher = PubMedSearcher(email=email, api_key=api_key)
        self.processor = AbstractProcessor()

    def fetch_by_pmid(self, pmid: str) -> Optional[ProcessedAbstract]:
        """
        根据 PMID 获取摘要

        Args:
            pmid: PubMed ID

        Returns:
            处理后的摘要
        """
        articles = self.searcher.fetch_details([pmid])
        if articles:
            return self.processor.process(articles[0])
        return None

    def fetch_by_pmids(self, pmids: List[str]) -> List[ProcessedAbstract]:
        """
        批量获取摘要

        Args:
            pmids: PMID 列表

        Returns:
            处理后的摘要列表
        """
        articles = self.searcher.fetch_details(pmids)
        return self.processor.process_batch(articles)


if __name__ == "__main__":
    # 测试
    processor = AbstractProcessor()

    # 测试文本清洗
    test_abstract = """
    BACKGROUND: The BRAF gene is frequently mutated in various cancers.
    METHODS: We analyzed 100 patients with lung cancer.
    RESULTS: BRAF V600E mutation was found in 15% of patients.
    CONCLUSIONS: BRAF mutations are common in lung cancer.
    """

    sections = processor.extract_sections(test_abstract)
    for section, content in sections.items():
        print(f"{section}: {content[:50]}...")