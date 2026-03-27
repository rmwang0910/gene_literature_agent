"""
缓存管理模块 - 持久化存储API响应数据
"""
import os
import json
import time
import hashlib
import logging
from typing import Dict, Optional, Any, List
from dataclasses import dataclass, asdict
from pathlib import Path
from datetime import datetime, timedelta

logger = logging.getLogger(__name__)


@dataclass
class CacheEntry:
    """缓存条目"""
    key: str
    value: Any
    created_at: float
    expires_at: float
    source: str = ""

    def is_expired(self) -> bool:
        return time.time() > self.expires_at

    def to_dict(self) -> Dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Dict) -> 'CacheEntry':
        return cls(**data)


class CitationCache:
    """引用数据缓存管理器"""

    DEFAULT_TTL_DAYS = 30  # 默认缓存30天

    def __init__(self, cache_dir: str = None, ttl_days: int = None):
        """
        初始化缓存管理器

        Args:
            cache_dir: 缓存目录
            ttl_days: 缓存有效期（天）
        """
        self.cache_dir = Path(cache_dir or os.path.expanduser("~/.gene_literature_cache"))
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.ttl_days = ttl_days or self.DEFAULT_TTL_DAYS
        self.cache_file = self.cache_dir / "citation_cache.json"
        self._cache: Dict[str, CacheEntry] = {}
        self._load()

    def _load(self):
        """从文件加载缓存"""
        if not self.cache_file.exists():
            return

        try:
            with open(self.cache_file, 'r', encoding='utf-8') as f:
                data = json.load(f)

            for key, entry_data in data.items():
                entry = CacheEntry.from_dict(entry_data)
                # 移除过期条目
                if not entry.is_expired():
                    self._cache[key] = entry

            logger.info(f"加载了 {len(self._cache)} 条缓存记录")

        except Exception as e:
            logger.warning(f"加载缓存失败: {e}")
            self._cache = {}

    def _save(self):
        """保存缓存到文件"""
        try:
            data = {k: v.to_dict() for k, v in self._cache.items()}
            with open(self.cache_file, 'w', encoding='utf-8') as f:
                json.dump(data, f, ensure_ascii=False, indent=2)
        except Exception as e:
            logger.warning(f"保存缓存失败: {e}")

    def _make_key(self, pmid: str) -> str:
        """生成缓存键"""
        return f"pmid:{pmid}"

    def get(self, pmid: str) -> Optional[Dict]:
        """
        获取缓存的引用数据

        Args:
            pmid: PubMed ID

        Returns:
            缓存的数据或None
        """
        key = self._make_key(pmid)
        entry = self._cache.get(key)

        if entry is None:
            return None

        if entry.is_expired():
            del self._cache[key]
            self._save()
            return None

        return entry.value

    def set(self, pmid: str, data: Dict, source: str = ""):
        """
        设置缓存

        Args:
            pmid: PubMed ID
            data: 要缓存的数据
            source: 数据来源
        """
        key = self._make_key(pmid)
        now = time.time()
        expires_at = now + (self.ttl_days * 24 * 3600)

        entry = CacheEntry(
            key=key,
            value=data,
            created_at=now,
            expires_at=expires_at,
            source=source
        )

        self._cache[key] = entry
        self._save()

    def get_batch(self, pmids: List[str]) -> Dict[str, Dict]:
        """
        批量获取缓存

        Args:
            pmids: PubMed ID列表

        Returns:
            {pmid: data} 字典（仅包含命中的）
        """
        results = {}
        for pmid in pmids:
            data = self.get(pmid)
            if data is not None:
                results[pmid] = data
        return results

    def set_batch(self, data: Dict[str, Dict], source: str = ""):
        """
        批量设置缓存

        Args:
            data: {pmid: data} 字典
            source: 数据来源
        """
        for pmid, value in data.items():
            self.set(pmid, value, source)

    def get_missing_pmids(self, pmids: List[str]) -> List[str]:
        """
        获取未缓存的PMID列表

        Args:
            pmids: 请求的PMID列表

        Returns:
            未缓存的PMID列表
        """
        return [p for p in pmids if self.get(p) is None]

    def clear_expired(self):
        """清除过期缓存"""
        expired_keys = [k for k, v in self._cache.items() if v.is_expired()]
        for key in expired_keys:
            del self._cache[key]

        if expired_keys:
            self._save()
            logger.info(f"清除了 {len(expired_keys)} 条过期缓存")

    def clear_all(self):
        """清除所有缓存"""
        self._cache = {}
        if self.cache_file.exists():
            self.cache_file.unlink()
        logger.info("已清除所有缓存")

    def get_stats(self) -> Dict[str, Any]:
        """获取缓存统计信息"""
        total = len(self._cache)
        expired = sum(1 for v in self._cache.values() if v.is_expired())

        return {
            "total_entries": total,
            "valid_entries": total - expired,
            "expired_entries": expired,
            "cache_file": str(self.cache_file),
            "ttl_days": self.ttl_days
        }


class ConclusionCache:
    """结论缓存管理器（缓存LLM提取结果）"""

    DEFAULT_TTL_DAYS = 7  # 结论缓存7天（LLM结果相对稳定）

    def __init__(self, cache_dir: str = None, ttl_days: int = None):
        self.cache_dir = Path(cache_dir or os.path.expanduser("~/.gene_literature_cache"))
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.ttl_days = ttl_days or self.DEFAULT_TTL_DAYS
        self.cache_file = self.cache_dir / "conclusion_cache.json"
        self._cache: Dict[str, CacheEntry] = {}
        self._load()

    def _load(self):
        """加载缓存"""
        if not self.cache_file.exists():
            return

        try:
            with open(self.cache_file, 'r', encoding='utf-8') as f:
                data = json.load(f)

            for key, entry_data in data.items():
                entry = CacheEntry.from_dict(entry_data)
                if not entry.is_expired():
                    self._cache[key] = entry

        except Exception as e:
            logger.warning(f"加载结论缓存失败: {e}")

    def _save(self):
        """保存缓存"""
        try:
            data = {k: v.to_dict() for k, v in self._cache.items()}
            with open(self.cache_file, 'w', encoding='utf-8') as f:
                json.dump(data, f, ensure_ascii=False, indent=2)
        except Exception as e:
            logger.warning(f"保存结论缓存失败: {e}")

    def _make_key(self, gene: str, pmid: str, abstract_hash: str = None) -> str:
        """生成缓存键"""
        if abstract_hash:
            return f"conclusion:{gene}:{pmid}:{abstract_hash}"
        return f"conclusion:{gene}:{pmid}"

    def get(self, gene: str, pmid: str, abstract_hash: str = None) -> Optional[Dict]:
        """获取缓存的结论"""
        key = self._make_key(gene, pmid, abstract_hash)
        entry = self._cache.get(key)

        if entry is None or entry.is_expired():
            return None

        return entry.value

    def set(self, gene: str, pmid: str, conclusion: Dict, abstract_hash: str = None):
        """设置结论缓存"""
        key = self._make_key(gene, pmid, abstract_hash)
        now = time.time()
        expires_at = now + (self.ttl_days * 24 * 3600)

        entry = CacheEntry(
            key=key,
            value=conclusion,
            created_at=now,
            expires_at=expires_at,
            source="llm"
        )

        self._cache[key] = entry
        self._save()

    def clear_all(self):
        """清除所有缓存"""
        self._cache = {}
        if self.cache_file.exists():
            self.cache_file.unlink()


def compute_text_hash(text: str) -> str:
    """计算文本哈希值"""
    return hashlib.md5(text.encode('utf-8')).hexdigest()[:8]


# 全局缓存实例（延迟初始化）
_citation_cache: Optional[CitationCache] = None
_conclusion_cache: Optional[ConclusionCache] = None


def get_citation_cache() -> CitationCache:
    """获取全局引用缓存实例"""
    global _citation_cache
    if _citation_cache is None:
        _citation_cache = CitationCache()
    return _citation_cache


def get_conclusion_cache() -> ConclusionCache:
    """获取全局结论缓存实例"""
    global _conclusion_cache
    if _conclusion_cache is None:
        _conclusion_cache = ConclusionCache()
    return _conclusion_cache


if __name__ == "__main__":
    # 测试缓存
    cache = CitationCache()
    print(f"缓存统计: {cache.get_stats()}")

    # 测试存储
    cache.set("12345678", {"citation_count": 100, "journal": "Nature"})
    print(f"获取缓存: {cache.get('12345678')}")

    # 测试批量
    missing = cache.get_missing_pmids(["12345678", "87654321"])
    print(f"缺失的PMID: {missing}")
