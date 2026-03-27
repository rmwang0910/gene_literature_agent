"""
结论提取模块 - 使用 LLM 从摘要中提取结构化结论
"""
import json
import logging
import os
import re
import time
from typing import List, Dict, Optional, Any
from dataclasses import dataclass, field, asdict
from abc import ABC, abstractmethod

from ..config import default_config
from ..prompts import EXTRACTION_PROMPT_TEMPLATE
from .abstract_fetcher import ProcessedAbstract
from ..utils.helpers import normalize_pathways, normalize_mutation_sites
from ..constants import MUTANT_KEYWORDS, WILDTYPE_KEYWORDS

# 配置日志
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class GeneConclusion:
    """基因结论数据类"""
    gene: str
    pmid: str
    title: str = ""
    disease_relation: str = ""
    pathways: List[str] = field(default_factory=list)
    clinical_significance: str = ""
    mutation_effects: str = ""
    mutation_sites: List[str] = field(default_factory=list)  # 具体突变位点列表
    mutation_disease_associations: List[Dict[str, str]] = field(default_factory=list)  # 突变-疾病关联
    gene_function_context: str = ""  # 基因功能语境（野生型/突变型）
    key_findings: List[str] = field(default_factory=list)
    confidence: str = "medium"
    year: int = 0
    citation_count: int = 0
    journal: str = ""
    impact_factor: float = 0.0
    doi: str = ""
    publication_type: str = ""  # 文献类型: "review", "research", ""(未知)

    def to_dict(self) -> Dict[str, Any]:
        """转换为字典"""
        return asdict(self)


def postprocess_gene_context(conclusion: GeneConclusion) -> GeneConclusion:
    """
    后处理基因结论，统一格式

    确保所有条目都包含野生型和突变型描述，格式统一为：
    "野生型：[功能]；突变型：[效应]"

    缺失信息时明确标注"摘要未提及"

    Args:
        conclusion: 原始结论

    Returns:
        处理后的结论
    """
    gene = conclusion.gene or ""
    ctx = conclusion.gene_function_context or ""

    # 检查现有语境中是否已包含野生型和突变型
    has_wildtype_in_ctx = "野生型" in ctx
    has_mutant_in_ctx = "突变型" in ctx

    # 提取文本用于关键词检测
    text_to_check = " ".join([
        conclusion.disease_relation,
        conclusion.clinical_significance,
        conclusion.mutation_effects,
        " ".join(conclusion.key_findings)
    ]).lower()

    has_mutation_kw = any(kw.lower() in text_to_check for kw in MUTANT_KEYWORDS)
    has_wildtype_kw = any(kw.lower() in text_to_check for kw in WILDTYPE_KEYWORDS)

    mutation_detail = conclusion.mutation_effects if conclusion.mutation_effects and conclusion.mutation_effects != "未提及" else ""

    # 情况1：已有完整的对比格式
    if has_wildtype_in_ctx and has_mutant_in_ctx:
        return conclusion

    # 情况2：仅提及突变型，需补充野生型标注
    if has_mutant_in_ctx and not has_wildtype_in_ctx:
        conclusion.gene_function_context = f"野生型：摘要未提及；{ctx}"
        return conclusion

    # 情况3：仅提及野生型，需补充突变型标注
    if has_wildtype_in_ctx and not has_mutant_in_ctx:
        if has_mutation_kw and mutation_detail:
            conclusion.gene_function_context = f"{ctx}；突变型：{mutation_detail}"
        elif has_mutation_kw:
            conclusion.gene_function_context = f"{ctx}；突变型：功能改变"
        else:
            conclusion.gene_function_context = f"{ctx}；突变型：摘要未提及"
        return conclusion

    # 情况4：都未明确提及，根据关键词推断
    if has_mutation_kw and not has_wildtype_kw:
        if mutation_detail:
            conclusion.gene_function_context = f"野生型：摘要未提及；突变型：{mutation_detail}"
        else:
            conclusion.gene_function_context = "野生型：摘要未提及；突变型：功能改变"
    elif has_wildtype_kw and not has_mutation_kw:
        conclusion.gene_function_context = "野生型：功能正常；突变型：摘要未提及"
    elif has_mutation_kw and has_wildtype_kw:
        if mutation_detail:
            conclusion.gene_function_context = f"野生型：正常功能；突变型：{mutation_detail}"
        else:
            conclusion.gene_function_context = "野生型：正常功能；突变型：功能改变"
    else:
        conclusion.gene_function_context = "野生型：摘要未提及；突变型：摘要未提及"

    return conclusion


class LLMProvider(ABC):
    """LLM 提供者抽象基类"""

    @abstractmethod
    def generate(self, prompt: str, max_tokens: int = 1000, temperature: float = 0.3) -> str:
        """生成文本"""
        pass

    def generate_structured(
        self,
        prompt: str,
        schema: Dict[str, Any],
        max_tokens: int = 1000,
        temperature: float = 0.3
    ) -> Dict[str, Any]:
        """生成结构化输出"""
        response = self.generate(prompt, max_tokens, temperature)

        # 尝试提取 JSON
        content = response.strip()
        if content.startswith("```json"):
            content = content[7:]
        if content.startswith("```"):
            content = content[3:]
        if content.endswith("```"):
            content = content[:-3]

        try:
            return json.loads(content.strip())
        except json.JSONDecodeError:
            # 尝试正则提取
            json_match = re.search(r'\{[\s\S]*\}', content)
            if json_match:
                return json.loads(json_match.group())
            return {"raw_response": response, "parse_error": True}


class OpenAIProvider(LLMProvider):
    """OpenAI API 提供者（参考 bioinfo_ai_literature_daily）"""

    def __init__(
        self,
        api_key: Optional[str] = None,
        model: Optional[str] = None,
        base_url: Optional[str] = None,
        timeout: Optional[float] = None,
        max_retries: int = 3,
        retry_backoff: float = 1.5
    ):
        self.api_key = api_key or default_config.llm_api_key
        self.model = model or default_config.llm_model or "qwen-plus"
        self.base_url = base_url or default_config.llm_base_url or "https://dashscope.aliyuncs.com/compatible-mode/v1"
        self.timeout = timeout or default_config.llm_timeout or 120
        self.max_retries = max(1, int(os.getenv("LLM_MAX_RETRIES", str(max_retries))))
        self.retry_backoff = max(0.0, float(os.getenv("LLM_RETRY_BACKOFF", str(retry_backoff))))

        # 可重试的状态码
        retry_status_codes_raw = os.getenv("LLM_RETRY_STATUS_CODES", "408,409,429,500,502,503,504")
        self.retry_status_codes = {
            int(code.strip())
            for code in retry_status_codes_raw.split(",")
            if code.strip().isdigit()
        }

        # 初始化 OpenAI 客户端
        try:
            from openai import OpenAI, APIConnectionError, APITimeoutError, RateLimitError, APIStatusError
            self.OpenAI = OpenAI
            self.APIConnectionError = APIConnectionError
            self.APITimeoutError = APITimeoutError
            self.RateLimitError = RateLimitError
            self.APIStatusError = APIStatusError
            self.client = OpenAI(
                api_key=self.api_key,
                base_url=self.base_url,
                timeout=self.timeout
            )
        except ImportError:
            raise ImportError("请安装 openai: pip install openai")

        logger.info(f"OpenAI Provider 初始化: base_url={self.base_url}, model={self.model}, timeout={self.timeout}")

    def generate(
        self,
        prompt: str,
        max_tokens: int = 1000,
        temperature: float = 0.3
    ) -> str:
        """调用 OpenAI API，带重试机制"""
        last_exception = None

        for attempt in range(1, self.max_retries + 1):
            try:
                response = self.client.chat.completions.create(
                    model=self.model,
                    messages=[
                        {"role": "system", "content": "你是一个专业的生物医学信息提取助手，擅长从文献摘要中提取结构化信息。请始终以JSON格式输出结果。"},
                        {"role": "user", "content": prompt}
                    ],
                    max_tokens=max_tokens,
                    temperature=temperature
                )
                return response.choices[0].message.content or ""

            except Exception as e:
                last_exception = e
                can_retry = self._can_retry(e) and attempt < self.max_retries

                if can_retry:
                    wait_seconds = 0.0 if self.retry_backoff <= 0 else self.retry_backoff ** (attempt - 1)
                    logger.warning(
                        f"LLM请求失败（第{attempt}/{self.max_retries}次）: {e}; "
                        f"{wait_seconds:.1f}s 后重试"
                    )
                    if wait_seconds > 0:
                        time.sleep(wait_seconds)
                    continue

                logger.error(f"LLM生成失败: {e}")
                raise

        if last_exception:
            raise last_exception
        raise RuntimeError("LLM生成失败：未知错误")

    def _can_retry(self, error: Exception) -> bool:
        """判断是否可以重试"""
        # 连接错误、超时、限流可以重试
        if isinstance(error, (self.APIConnectionError, self.APITimeoutError, self.RateLimitError)):
            return True

        # API状态错误，检查状态码
        if isinstance(error, self.APIStatusError):
            status_code = getattr(error, "status_code", None)
            return status_code in self.retry_status_codes

        # 检查错误消息中的可重试信号
        message = str(error).lower()
        transient_signals = [
            "connection error",
            "timed out",
            "timeout",
            "temporarily unavailable",
            "try again",
            "rate limit"
        ]
        return any(signal in message for signal in transient_signals)


class GroqProvider(LLMProvider):
    """Groq API 提供者（使用 Llama 模型）"""

    def __init__(self, api_key: Optional[str] = None, model: str = "llama-3.3-70b-versatile"):
        self.api_key = api_key
        self.model = model

    def generate(self, prompt: str, max_tokens: int = 1000, temperature: float = 0.3) -> str:
        """调用 Groq API"""
        try:
            import groq

            client = groq.Groq(api_key=self.api_key)

            response = client.chat.completions.create(
                model=self.model,
                messages=[
                    {"role": "system", "content": "你是一个专业的生物医学信息提取助手。请始终以JSON格式输出结果。"},
                    {"role": "user", "content": prompt}
                ],
                temperature=temperature,
                max_tokens=max_tokens
            )

            return response.choices[0].message.content

        except ImportError:
            raise ImportError("请安装 groq: pip install groq")
        except Exception as e:
            raise RuntimeError(f"Groq API 调用失败: {e}")


class MockProvider(LLMProvider):
    """Mock 提供者（用于测试）"""

    def generate(self, prompt: str, max_tokens: int = 1000, temperature: float = 0.3) -> str:
        """返回模拟结果"""
        return json.dumps({
            "gene": "UNKNOWN",
            "disease_relation": "模拟结果 - 促进肿瘤生长",
            "pathways": ["MAPK", "ERK"],
            "clinical_significance": "与不良预后相关",
            "mutation_effects": "V600E突变激活激酶活性",
            "key_findings": ["发现1", "发现2", "发现3"],
            "confidence": "low",
            "pmid": "00000000"
        })


class ConclusionExtractor:
    """结论提取器"""

    def __init__(
        self,
        provider: Optional[LLMProvider] = None,
        cache_dir: Optional[str] = None
    ):
        """
        初始化结论提取器

        Args:
            provider: LLM 提供者
            cache_dir: 缓存目录
        """
        self.provider = provider or self._get_default_provider()
        self.cache_dir = cache_dir or default_config.cache_dir

    def _get_default_provider(self) -> LLMProvider:
        """获取默认 LLM 提供者"""
        llm_provider = default_config.llm_provider.lower()

        if llm_provider == "openai":
            return OpenAIProvider()
        elif llm_provider == "groq":
            return GroqProvider()
        else:
            # 默认尝试 OpenAI
            if default_config.llm_api_key:
                return OpenAIProvider()
            return MockProvider()

    def _build_prompt(self, gene: str, abstract: ProcessedAbstract) -> str:
        """构建提示词"""
        return EXTRACTION_PROMPT_TEMPLATE.format(
            gene=gene,
            abstract=abstract.abstract,
            pmid=abstract.pmid
        )

    def _parse_response(self, response: str, gene: str, pmid: str, title: str = "") -> GeneConclusion:
        """解析 LLM 响应"""
        try:
            # 尝试提取 JSON
            json_match = re.search(r'\{[\s\S]*\}', response)
            if json_match:
                data = json.loads(json_match.group())
            else:
                data = json.loads(response)

            conclusion = GeneConclusion(
                gene=gene,
                pmid=pmid,
                title=title,
                disease_relation=data.get("disease_relation", ""),
                pathways=normalize_pathways(data.get("pathways", [])),
                clinical_significance=data.get("clinical_significance", ""),
                mutation_effects=data.get("mutation_effects", ""),
                mutation_sites=normalize_mutation_sites(data.get("mutation_sites", [])),
                mutation_disease_associations=data.get("mutation_disease_associations", []),
                gene_function_context=data.get("gene_function_context", ""),
                key_findings=data.get("key_findings", []),
                confidence=data.get("confidence", "medium")
            )

            # 后处理：自动补全基因功能语境
            conclusion = postprocess_gene_context(conclusion)

            return conclusion

        except (json.JSONDecodeError, KeyError) as e:
            logger.warning(f"解析响应失败: {e}")
            return GeneConclusion(
                gene=gene,
                pmid=pmid,
                title=title,
                disease_relation="解析失败",
                confidence="low"
            )

    def extract(self, gene: str, abstract: ProcessedAbstract) -> GeneConclusion:
        """
        从单篇摘要中提取结论

        Args:
            gene: 基因名
            abstract: 处理后的摘要

        Returns:
            GeneConclusion 对象
        """
        prompt = self._build_prompt(gene, abstract)

        try:
            response = self.provider.generate(prompt)
            return self._parse_response(response, gene, abstract.pmid, abstract.title)
        except Exception as e:
            logger.error(f"提取结论失败 (PMID: {abstract.pmid}): {e}")
            return GeneConclusion(
                gene=gene,
                pmid=abstract.pmid,
                title=abstract.title,
                disease_relation=f"提取失败: {str(e)}",
                confidence="low"
            )

    def extract_batch(
        self,
        gene: str,
        abstracts: List[ProcessedAbstract],
        delay: float = 0.5
    ) -> List[GeneConclusion]:
        """
        批量提取结论

        Args:
            gene: 基因名
            abstracts: 摘要列表
            delay: 请求间隔（秒）

        Returns:
            结论列表
        """
        conclusions = []

        for i, abstract in enumerate(abstracts):
            if i > 0:
                time.sleep(delay)

            logger.info(f"正在处理 {gene} - PMID {abstract.pmid} ({i+1}/{len(abstracts)})")
            conclusion = self.extract(gene, abstract)
            conclusions.append(conclusion)

        return conclusions

    def extract_with_retry(
        self,
        gene: str,
        abstract: ProcessedAbstract,
        max_retries: int = 3
    ) -> GeneConclusion:
        """
        带重试的提取

        Args:
            gene: 基因名
            abstract: 摘要
            max_retries: 最大重试次数

        Returns:
            GeneConclusion 对象
        """
        last_error = None

        for attempt in range(max_retries):
            try:
                return self.extract(gene, abstract)
            except Exception as e:
                last_error = e
                if attempt < max_retries - 1:
                    time.sleep(2 ** attempt)  # 指数退避

        return GeneConclusion(
            gene=gene,
            pmid=abstract.pmid,
            title=abstract.title,
            disease_relation=f"提取失败: {str(last_error)}",
            confidence="low"
        )


def create_provider(
    provider_type: str,
    api_key: Optional[str] = None,
    model: Optional[str] = None,
    base_url: Optional[str] = None
) -> LLMProvider:
    """
    创建 LLM 提供者

    Args:
        provider_type: 提供者类型 (openai, groq, mock)
        api_key: API Key
        model: 模型名称
        base_url: 基础 URL

    Returns:
        LLMProvider 对象
    """
    provider_type = provider_type.lower()

    if provider_type == "openai":
        return OpenAIProvider(
            api_key=api_key,
            model=model or "qwen-plus",
            base_url=base_url
        )
    elif provider_type == "groq":
        return GroqProvider(
            api_key=api_key,
            model=model or "llama-3.3-70b-versatile"
        )
    else:
        return MockProvider()


if __name__ == "__main__":
    # 测试
    from .abstract_fetcher import ProcessedAbstract

    # 创建测试摘要
    test_abstract = ProcessedAbstract(
        pmid="12345678",
        title="BRAF mutations in lung cancer",
        abstract="Background: BRAF gene mutations are common in various cancers. Methods: We analyzed 100 lung cancer patients. Results: BRAF V600E mutation was found in 15% of patients and was associated with poor prognosis. Conclusions: BRAF mutations are important therapeutic targets in lung cancer.",
        word_count=50,
        sections={}
    )

    # 使用 Mock 提供者测试
    extractor = ConclusionExtractor(provider=MockProvider())
    conclusion = extractor.extract("BRAF", test_abstract)
    print(json.dumps(conclusion.to_dict(), indent=2, ensure_ascii=False))