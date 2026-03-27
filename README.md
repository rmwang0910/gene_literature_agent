# Gene Literature Agent

基因文献智能检索助手 - 一个智能的基因文献分析工具，自动从 PubMed 检索相关文献，提取核心结论并整合权威数据库证据。

## Features

### 核心功能
- **基因名称规范化** - 自动识别基因别名，转换为标准符号 (via MyGene.info)
- **PubMed 智能检索** - 支持背景词、时间范围、文章类型过滤
- **LLM 结论提取** - 利用大语言模型从摘要中提取结构化结论
- **野生型/突变型区分** - 智能区分野生型和突变型基因的功能描述
- **冲突检测** - 自动检测文献结论之间的矛盾与争议
- **多格式报告** - 支持 Markdown、Excel、JSON 输出
- **引用数据增强** - 通过 OpenAlex/CrossRef 获取引用信息和影响因子

### 外部数据源集成 (NEW)
- **Open Targets** - 基因-疾病关联证据、靶向药物信息
- **ClinVar** - 变异临床意义、ACMG 分类
- **COSMIC** - 癌症体细胞突变、热点突变、Cancer Gene Census
- **KEGG/Reactome** - 信号通路实时查询

## Installation

```bash
# Clone the repository
git clone https://github.com/your-username/gene_literature_agent.git
cd gene_literature_agent

# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # Linux/macOS
# or: venv\Scripts\activate  # Windows

# Install dependencies
pip install -r requirements.txt
```

## Configuration

**Important: Never commit your API keys to Git!**

1. Copy the example config file:
   ```bash
   cp config.yaml.example config.yaml
   ```

2. Edit `config.yaml` and fill in your credentials:
   ```yaml
   ncbi:
     email: "your_email@example.com"  # Required for NCBI API
     api_key: ""  # Optional, get from https://www.ncbi.nlm.nih.gov/account/settings/

   llm:
     api_key: "your_api_key_here"  # Required
     model: "qwen-plus"  # or gpt-4o-mini, etc.
     base_url: "https://dashscope.aliyuncs.com/compatible-mode/v1"  # Optional
   ```

3. Alternatively, use environment variables:
   ```bash
   export NCBI_EMAIL="your_email@example.com"
   export LLM_API_KEY="your_api_key"
   export LLM_MODEL="qwen-plus"
   export LLM_BASE_URL="https://dashscope.aliyuncs.com/compatible-mode/v1"
   ```

### Supported LLM Providers

| Provider | Model Examples | Base URL |
|----------|---------------|----------|
| OpenAI | gpt-4o-mini, gpt-4o | (default) |
| Alibaba Qwen | qwen-plus, qwen-turbo, qwen-max | https://dashscope.aliyuncs.com/compatible-mode/v1 |
| Other OpenAI-compatible | varies | Your custom endpoint |

## Usage

### Command Line

```bash
# Analyze a single gene
python cli.py BRAF

# Analyze multiple genes with disease context
python cli.py BRAF TP53 EGFR --background "lung cancer"

# Specify result count and output format
python cli.py BRAF --max-results 20 --format json

# Read genes from file
python cli.py --file genes.txt --background "breast cancer"

# Estimate cost before running
python cli.py BRAF TP53 --estimate-cost
```

### Command Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `genes` | Gene names to analyze | Required |
| `-f, --file` | Gene list file (one per line) | - |
| `-b, --background` | Disease/tissue context | None |
| `-m, --max-results` | Max papers per gene | 10 |
| `--format` | Output format (markdown/json/excel/all) | all |
| `--estimate-cost` | Estimate API cost only | - |

### Python API

```python
from gene_literature_agent import (
    GeneNormalizer,
    PubMedSearcher,
    ConclusionExtractor,
    ConflictDetector,
    ReportGenerator
)

# Normalize gene name
normalizer = GeneNormalizer()
gene_info = normalizer.normalize("BRAF")
print(f"Symbol: {gene_info.symbol}, Aliases: {gene_info.aliases}")

# Search PubMed
searcher = PubMedSearcher()
query = searcher.build_query("BRAF", background="lung cancer")
results = searcher.search_and_fetch(query, max_results=10)

# Extract conclusions
extractor = ConclusionExtractor()
conclusions = []
for article in results.articles:
    conclusion = extractor.extract("BRAF", article.abstract, article.pmid)
    conclusions.append(conclusion)

# Detect conflicts
detector = ConflictDetector()
conflict_report = detector.detect("BRAF", conclusions)

# Generate report with external database evidence
generator = ReportGenerator(enable_external_sources=True)
markdown = generator.generate_markdown("BRAF", conclusions, conflict_report, variant="V600E")
```

### External Database APIs (NEW)

```python
from datasources import (
    OpenTargetsClient,
    ClinVarClient,
    COSMICClient,
    PathwayClient
)
from datasources.evidence_enricher import EvidenceEnricher

# Query Open Targets for gene-disease associations
ot_client = OpenTargetsClient()
summary = ot_client.get_gene_disease_summary("BRAF", top_n=10)
print(f"Top diseases: {[d['disease_name'] for d in summary['diseases']]}")

# Query ClinVar for variant clinical significance
cv_client = ClinVarClient()
result = cv_client.get_variant_clinical_significance("BRAF", "V600E")
print(f"Clinical significance: {result.clinical_significance} ({result.star_rating} stars)")

# Query COSMIC for hotspot mutations
cosmic_client = COSMICClient()
hotspots = cosmic_client.get_hotspot_mutations("BRAF")
print(f"Top hotspots: {[h.mutation for h in hotspots[:5]]}")

# Query KEGG/Reactome for pathways
pathway_client = PathwayClient()
pathways = pathway_client.get_gene_pathways("BRAF")
print(f"Pathways: {[p.name for p in pathways[:5]]}")

# Integrated evidence enrichment
enricher = EvidenceEnricher()
evidence = enricher.enrich_gene("BRAF", variant="V600E")
print(f"Is cancer gene: {evidence.is_cancer_gene}")
print(f"Cancer gene role: {evidence.cancer_gene_role}")
print(f"Sources queried: {evidence.sources_queried}")
```

## Project Structure

```
gene_literature_agent/
├── __init__.py              # Package entry point
├── cli.py                   # Command line interface
├── config.py                # Configuration loader
├── config.yaml.example      # Configuration template (safe to commit)
├── constants.py             # Constants and enums
├── requirements.txt         # Python dependencies
├── README.md                # Documentation
├── .gitignore               # Git ignore rules
│
├── prompts/                 # LLM prompt templates
│   ├── __init__.py
│   ├── extraction.py        # Conclusion extraction prompts
│   ├── conflict.py          # Conflict detection prompts
│   └── summary.py           # Summary generation prompts
│
├── core/                    # Core business modules
│   ├── __init__.py
│   ├── gene_normalizer.py   # Gene name normalization (MyGene.info)
│   ├── pubmed_searcher.py   # PubMed search and retrieval
│   ├── abstract_fetcher.py  # Abstract processing
│   ├── citation_fetcher.py  # Citation data (OpenAlex/CrossRef)
│   ├── conclusion_extractor.py  # LLM-based conclusion extraction
│   ├── conflict_detector.py # Conflict detection between papers
│   └── report_generator.py  # Report generation
│
├── datasources/             # External database integrations (NEW)
│   ├── __init__.py
│   ├── opentargets_client.py   # Open Targets Platform API
│   ├── clinvar_client.py       # ClinVar E-utilities API
│   ├── cosmic_client.py        # COSMIC cancer mutations
│   ├── pathway_client.py       # KEGG/Reactome REST APIs
│   └── evidence_enricher.py    # Multi-source evidence integration
│
├── skills/                  # Reference skill definitions (from OpenClaw)
│   ├── opentargets-database/
│   ├── clinvar-database/
│   ├── cosmic-database/
│   ├── bio-pathway-kegg-pathways/
│   ├── bio-pathway-reactome/
│   ├── tooluniverse-cancer-variant-interpretation/
│   └── tooluniverse-literature-deep-research/
│
├── utils/                   # Utility modules
│   ├── __init__.py
│   ├── helpers.py           # General utility functions
│   ├── text_utils.py        # Text processing utilities
│   ├── cache_manager.py     # Cache management
│   ├── journal_data.py      # Journal impact factors
│   ├── visualizer.py        # Visualization generation
│   ├── pathway_rag.py       # Pathway RAG system (legacy)
│   ├── pathway_downloader.py # Pathway data downloader
│   └── hotspot_downloader.py # Cancer hotspot mutations
│
└── data/                    # Data files
    ├── journal/
    │   └── scimagojr_2024.csv    # Journal rankings
    ├── pathway/
    │   └── pathway_mapping.json  # Pathway name mappings
    ├── mutation/
    │   ├── hotspot_mutations.json     # Cancer hotspot data
    │   └── hotspot_mutations_tcga.json
    └── cosmic/              # COSMIC local data (if downloaded)
        ├── cancer_gene_census.json
        └── hotspot_mutations.json
```

## Data Sources

### Literature & Citations
- [NCBI PubMed](https://pubmed.ncbi.nlm.nih.gov/) - Literature search
- [OpenAlex](https://openalex.org/) - Citation data
- [CrossRef](https://www.crossref.org/) - DOI and citation data
- [SCImago Journal Rank](https://www.scimagojr.com/) - Journal impact factors

### Gene & Variant Annotation
- [MyGene.info](https://mygene.info/) - Gene name normalization
- [Open Targets Platform](https://platform.opentargets.org/) - Gene-disease associations, drug targets
- [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) - Variant clinical significance (ACMG classification)
- [COSMIC](https://cancer.sanger.ac.uk/cosmic) - Cancer somatic mutations, Cancer Gene Census

### Pathway Databases
- [KEGG](https://www.kegg.jp/) - Metabolic and signaling pathways
- [Reactome](https://reactome.org/) - Biological pathway database

### Cancer Genomics
- [cBioPortal](https://www.cbioportal.org/) - Cancer genomics (mutation prevalence)

## Output Example

The tool generates comprehensive reports including:

- Gene function summary (wild-type vs mutant)
- Literature evidence with citations
- Pathway associations
- Conflict analysis between papers
- Clinical significance
- Mutation hotspots (if applicable)

## License

MIT License

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Acknowledgments

- NCBI E-utilities for PubMed access
- MyGene.info for gene information
- OpenAlex for open citation data
