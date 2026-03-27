#!/usr/bin/env python
"""
基因文献智能检索助手 - 命令行入口
"""
import argparse
import sys
from typing import List, Optional

# 导入模块
from .config import default_config
from .core import (
    GeneNormalizer,
    GeneInfo,
    PubMedSearcher,
    AbstractProcessor,
    ConclusionExtractor,
    ConflictDetector,
    ReportGenerator,
)
from .utils import estimate_cost


def analyze_gene(
    gene: str,
    background: Optional[str] = None,
    max_results: int = 10,
    output_format: str = "markdown",
    sort: str = "relevance",
    article_types: Optional[List[str]] = None,
    min_citations: int = 0,
    min_impact_factor: float = 0.0,
    open_access_only: bool = False
) -> dict:
    """
    分析单个基因

    Args:
        gene: 基因名称
        background: 疾病背景
        max_results: 最大检索文献数
        output_format: 输出格式
        sort: 排序方式 (relevance, pub_date, first_author, last_author, journal)
        article_types: 文章类型列表
        min_citations: 最小引用次数
        min_impact_factor: 最小影响因子
        open_access_only: 仅开放获取

    Returns:
        分析结果字典
    """
    print(f"\n{'='*60}")
    print(f"正在分析基因: {gene}")
    print(f"{'='*60}")

    # 初始化组件
    normalizer = GeneNormalizer()
    searcher = PubMedSearcher()
    processor = AbstractProcessor()
    extractor = ConclusionExtractor()
    detector = ConflictDetector()
    generator = ReportGenerator()

    # 1. 基因规范化
    print("\n[1/5] 规范化基因名称...")
    gene_info = normalizer.normalize(gene)
    if not gene_info:
        print(f"  ❌ 无法识别基因: {gene}")
        return {"gene": gene, "status": "error", "error": "无法识别基因名称"}

    print(f"  ✓ 官方符号: {gene_info.symbol}")
    print(f"  ✓ 全称: {gene_info.name}")
    if gene_info.aliases:
        print(f"  ✓ 别名: {', '.join(gene_info.aliases[:5])}")

    # 2. PubMed 检索
    print(f"\n[2/5] 检索 PubMed (最多 {max_results} 篇)...")
    query = searcher.build_query(
        gene_symbol=gene_info.symbol,
        aliases=gene_info.aliases,
        background=background,
        article_types=article_types
    )
    print(f"  检索式: {query}")

    search_result = searcher.search_and_fetch(
        query,
        max_results,
        sort=sort,
        min_citations=min_citations,
        min_impact_factor=min_impact_factor,
        open_access_only=open_access_only
    )
    print(f"  ✓ 找到 {search_result.total_count} 篇文献，获取 {len(search_result.articles)} 篇")

    if not search_result.articles:
        return {"gene": gene, "status": "no_results", "error": "未找到相关文献"}

    # 3. 处理摘要
    print(f"\n[3/5] 处理摘要...")
    abstracts = processor.process_batch(search_result.articles)
    print(f"  ✓ 处理完成 {len(abstracts)} 篇摘要")

    # 4. 提取结论
    print(f"\n[4/5] 提取结论...")
    conclusions = []
    for i, abstract in enumerate(abstracts):
        print(f"  处理 {i+1}/{len(abstracts)}: PMID {abstract.pmid}")
        conclusion = extractor.extract_with_retry(gene_info.symbol, abstract)
        conclusion.year = search_result.articles[i].year
        conclusions.append(conclusion)

    print(f"  ✓ 提取完成 {len(conclusions)} 条结论")

    # 5. 冲突检测
    print(f"\n[5/5] 检测冲突...")
    conflict_report = detector.detect(gene_info.symbol, conclusions)
    if conflict_report.conflicts:
        print(f"  ⚠️ 发现 {len(conflict_report.conflicts)} 处潜在冲突")
    else:
        print(f"  ✓ 未检测到明显冲突")

    # 生成报告
    print(f"\n生成报告...")
    if output_format == "markdown" or output_format == "all":
        md_content = generator.generate_markdown(
            gene_info.symbol,
            conclusions,
            conflict_report,
            [{"pmid": a.pmid, "title": a.title, "year": a.year} for a in search_result.articles]
        )
        filepath = generator.save_markdown(gene_info.symbol, md_content)
        print(f"  ✓ Markdown: {filepath}")

    if output_format == "json" or output_format == "all":
        filepath = generator.generate_json(gene_info.symbol, conclusions, conflict_report)
        print(f"  ✓ JSON: {filepath}")

    if output_format == "excel" or output_format == "all":
        try:
            filepath = generator.generate_excel(gene_info.symbol, conclusions, conflict_report)
            print(f"  ✓ Excel: {filepath}")
        except Exception as e:
            print(f"  ⚠️ Excel 生成失败: {e}")

    if output_format == "html" or output_format == "all":
        try:
            # 获取引用数据
            generator.enrich_conclusions_with_citations(conclusions)
            filepath = generator.generate_html(gene_info.symbol, conclusions, conflict_report)
            print(f"  ✓ HTML: {filepath}")
        except Exception as e:
            print(f"  ⚠️ HTML 生成失败: {e}")

    return {
        "gene": gene,
        "status": "completed",
        "gene_info": gene_info,
        "articles": search_result.articles,
        "conclusions": conclusions,
        "conflict_report": conflict_report
    }


def main():
    """命令行主函数"""
    parser = argparse.ArgumentParser(
        description="基因文献智能检索助手 - 自动检索PubMed文献并提取结构化结论",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 分析单个基因
  python cli.py BRAF

  # 分析多个基因，指定疾病背景
  python cli.py BRAF TP53 EGFR --background "lung cancer"

  # 按出版日期排序
  python cli.py BRAF --sort pub_date

  # 筛选综述文章
  python cli.py BRAF --article-types review

  # 筛选高引用文章（引用>=50）
  python cli.py BRAF --min-citations 50

  # 筛选高影响因子期刊（IF>=10）
  python cli.py BRAF --min-impact-factor 10

  # 仅开放获取文献
  python cli.py BRAF --open-access-only

  # 组合筛选：综述+高引用+高影响因子
  python cli.py BRAF --article-types review --min-citations 20 --min-impact-factor 5

  # 指定检索数量和输出格式
  python cli.py BRAF --max-results 20 --format json

  # 从文件读取基因列表
  python cli.py --file genes.txt --background "breast cancer"
        """
    )

    parser.add_argument(
        "genes",
        nargs="*",
        help="基因名称列表"
    )
    parser.add_argument(
        "-f", "--file",
        help="基因列表文件（每行一个基因）"
    )
    parser.add_argument(
        "-b", "--background",
        help="疾病/组织背景（如 'lung cancer'）"
    )
    parser.add_argument(
        "-m", "--max-results",
        type=int,
        default=10,
        help="每个基因最大检索文献数（默认: 10）"
    )
    parser.add_argument(
        "--format",
        choices=["markdown", "json", "excel", "html", "all"],
        default="all",
        help="报告输出格式（默认: all）"
    )
    parser.add_argument(
        "--sort",
        choices=["relevance", "pub_date", "first_author", "last_author", "journal"],
        default="relevance",
        help="排序方式（默认: relevance）"
    )
    parser.add_argument(
        "--article-types",
        nargs="+",
        choices=["review", "clinical_trial", "randomized_controlled_trial",
                 "meta_analysis", "case_report", "observational_study",
                 "systematic_review", "editorial", "letter", "comment"],
        help="文章类型筛选（可多选）"
    )
    parser.add_argument(
        "--min-citations",
        type=int,
        default=0,
        help="最小引用次数筛选（默认: 0，不限制）"
    )
    parser.add_argument(
        "--min-impact-factor",
        type=float,
        default=0.0,
        help="最小期刊影响因子筛选（默认: 0.0，不限制）"
    )
    parser.add_argument(
        "--open-access-only",
        action="store_true",
        help="仅检索开放获取文献"
    )
    parser.add_argument(
        "--estimate-cost",
        action="store_true",
        help="仅估算成本，不执行分析"
    )

    args = parser.parse_args()

    # 收集基因列表
    genes = list(args.genes)
    if args.file:
        try:
            with open(args.file, 'r', encoding='utf-8') as f:
                file_genes = [line.strip() for line in f if line.strip()]
                genes.extend(file_genes)
        except FileNotFoundError:
            print(f"错误: 文件不存在 - {args.file}")
            sys.exit(1)

    if not genes:
        parser.print_help()
        sys.exit(1)

    # 去重
    genes = list(dict.fromkeys(genes))

    # 估算成本
    if args.estimate_cost:
        cost = estimate_cost(len(genes) * args.max_results)
        print(f"\n成本估算:")
        print(f"  基因数量: {len(genes)}")
        print(f"  每基因文献数: {args.max_results}")
        print(f"  总摘要数: {len(genes) * args.max_results}")
        print(f"  预估 tokens: {cost['total_tokens']}")
        print(f"  预估成本: ${cost['estimated_cost_usd']:.4f}")
        sys.exit(0)

    # 显示分析信息
    print(f"\n{'#'*60}")
    print(f"# 基因文献智能检索助手")
    print(f"{'#'*60}")
    print(f"基因列表: {', '.join(genes)}")
    print(f"疾病背景: {args.background or '无'}")
    print(f"检索数量: {args.max_results} 篇/基因")
    print(f"排序方式: {args.sort}")
    if args.article_types:
        print(f"文章类型: {', '.join(args.article_types)}")
    if args.min_citations > 0:
        print(f"最小引用次数: {args.min_citations}")
    if args.min_impact_factor > 0:
        print(f"最小影响因子: {args.min_impact_factor}")
    if args.open_access_only:
        print(f"仅开放获取: 是")
    print(f"输出格式: {args.format}")

    # 执行分析
    results = {}
    for i, gene in enumerate(genes):
        print(f"\n进度: {i+1}/{len(genes)}")
        result = analyze_gene(
            gene=gene,
            background=args.background,
            max_results=args.max_results,
            output_format=args.format,
            sort=args.sort,
            article_types=args.article_types,
            min_citations=args.min_citations,
            min_impact_factor=args.min_impact_factor,
            open_access_only=args.open_access_only
        )
        results[gene] = result

    # 汇总
    print(f"\n{'='*60}")
    print("分析完成!")
    print(f"{'='*60}")

    success = sum(1 for r in results.values() if r["status"] == "completed")
    no_results = sum(1 for r in results.values() if r["status"] == "no_results")
    errors = sum(1 for r in results.values() if r["status"] == "error")

    print(f"成功: {success}, 无结果: {no_results}, 失败: {errors}")


if __name__ == "__main__":
    main()