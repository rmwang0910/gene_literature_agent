"""
可视化模块 - 生成基因文献分析图表
"""
from __future__ import annotations

import io
import logging
import base64
from datetime import datetime
from typing import List, Dict, Optional, Tuple, TYPE_CHECKING
from collections import Counter
from dataclasses import dataclass

try:
    import matplotlib
    matplotlib.use('Agg')  # 非交互式后端
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    import matplotlib.font_manager as fm
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

try:
    import networkx as nx
    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False

try:
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False

from .helpers import normalize_pathways

if TYPE_CHECKING:
    from ..core.conclusion_extractor import GeneConclusion
    from ..core.conflict_detector import ConflictReport
from ..constants import DISEASE_KEYWORDS

# 配置日志
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# 中文字体配置
def get_chinese_font():
    """获取可用的中文字体"""
    if not HAS_MATPLOTLIB:
        return 'sans-serif'

    chinese_fonts = [
        'Arial Unicode MS',
        'SimHei',
        'Heiti SC',
        'STHeiti',
        'Microsoft YaHei',
        'PingFang SC',
        'Hiragino Sans GB',
    ]

    try:
        available_fonts = [f.name for f in fm.fontManager.ttflist]
    except Exception:
        return 'sans-serif'

    for font in chinese_fonts:
        if font in available_fonts:
            return font

    return 'sans-serif'


def setup_chinese_font():
    """设置中文字体（全局）"""
    font = get_chinese_font()
    plt.rcParams['font.family'] = font
    plt.rcParams['axes.unicode_minus'] = False
    return font


@dataclass
class ChartConfig:
    """图表配置"""
    width: int = 800
    height: int = 500
    dpi: int = 100
    font_family: str = None
    title_fontsize: int = 14
    label_fontsize: int = 10

    def __post_init__(self):
        if self.font_family is None:
            self.font_family = get_chinese_font()


class GeneDiseaseNetwork:
    """基因-疾病关联网络图"""

    def __init__(self, config: Optional[ChartConfig] = None):
        self.config = config or ChartConfig()

    def build_network(
        self,
        conclusions: List[GeneConclusion]
    ) -> Tuple['nx.Graph', Dict[str, List[str]], Dict[str, int]]:
        """
        构建基因-疾病关联网络

        Args:
            conclusions: 结论列表

        Returns:
            (NetworkX图, 节点分类, 节点频次)
        """
        if not HAS_NETWORKX:
            raise ImportError("请安装 networkx: pip install networkx")

        G = nx.Graph()
        node_types = {"gene": [], "disease": []}
        node_counts = {}  # 记录每个节点的关联文献数

        # 疾病关键词（从配置加载，按优先级排序）
        disease_keywords = DISEASE_KEYWORDS

        for c in conclusions:
            # 添加基因节点
            if c.gene:
                if c.gene not in G.nodes:
                    G.add_node(c.gene, type="gene")
                    node_types["gene"].append(c.gene)
                    node_counts[c.gene] = 1
                else:
                    node_counts[c.gene] += 1

            # 从疾病关系中提取疾病
            if c.disease_relation:
                matched_diseases = set()
                for disease in disease_keywords:
                    if disease in c.disease_relation:
                        # 避免重复匹配（如匹配了"肺癌"就不再匹配"癌"）
                        skip = False
                        for matched in matched_diseases:
                            if disease in matched or matched in disease:
                                skip = True
                                break
                        if skip:
                            continue

                        # 添加疾病节点
                        if disease not in G.nodes:
                            G.add_node(disease, type="disease")
                            node_types["disease"].append(disease)
                            node_counts[disease] = 1
                        else:
                            node_counts[disease] += 1

                        matched_diseases.add(disease)

                        # 添加边
                        if c.gene and disease:
                            if G.has_edge(c.gene, disease):
                                G[c.gene][disease]["weight"] += 1
                            else:
                                G.add_edge(c.gene, disease, weight=1)

        return G, node_types, node_counts

    def plot_matplotlib(
        self,
        conclusions: List[GeneConclusion],
        output_path: Optional[str] = None
    ) -> Optional[str]:
        """
        使用 Matplotlib 绘制网络图

        Args:
            conclusions: 结论列表
            output_path: 输出路径

        Returns:
            Base64编码的图片或文件路径
        """
        if not HAS_MATPLOTLIB or not HAS_NETWORKX:
            logger.warning("缺少依赖库，无法生成网络图")
            return None

        # 设置中文字体
        setup_chinese_font()

        G, node_types, node_counts = self.build_network(conclusions)

        if len(G.nodes) == 0:
            logger.warning("没有足够的数据生成网络图")
            return None

        # 节点数量过少时返回None（后续用表格替代）
        if len(G.nodes) <= 2:
            logger.info("节点数量过少，建议使用表格展示")
            return None

        # 创建图形
        fig, ax = plt.subplots(figsize=(self.config.width / self.config.dpi,
                                          self.config.height / self.config.dpi))

        # 优化布局参数
        k_value = max(1.5, 3.0 / (len(G.nodes) ** 0.3))  # 根据节点数调整
        pos = nx.spring_layout(G, k=k_value, iterations=100, seed=42)

        # 分离基因和疾病节点
        gene_nodes = node_types["gene"]
        disease_nodes = node_types["disease"]

        # 根据频次计算节点大小
        def get_node_size(node, base_size=800, max_size=2500):
            count = node_counts.get(node, 1)
            # 使用对数缩放避免极端值
            import math
            size = base_size + (max_size - base_size) * min(1, math.log(count + 1) / 3)
            return size

        # 绘制基因节点
        if gene_nodes:
            gene_sizes = [get_node_size(n, 1000, 3000) for n in gene_nodes]
            nx.draw_networkx_nodes(G, pos, nodelist=gene_nodes,
                                   node_color='#2ecc71', node_size=gene_sizes,
                                   alpha=0.9, ax=ax, edgecolors='#27ae60', linewidths=2)

        # 绘制疾病节点
        if disease_nodes:
            disease_sizes = [get_node_size(n, 600, 2000) for n in disease_nodes]
            nx.draw_networkx_nodes(G, pos, nodelist=disease_nodes,
                                   node_color='#e74c3c', node_size=disease_sizes,
                                   alpha=0.9, ax=ax, edgecolors='#c0392b', linewidths=2)

        # 绘制边（粗细根据权重）
        edges = G.edges(data=True)
        if edges:
            edge_weights = [max(1, d['weight'] * 1.5) for (_, _, d) in edges]
            edge_colors = ['#95a5a6' if d['weight'] == 1 else '#34495e' for (_, _, d) in edges]
            nx.draw_networkx_edges(G, pos, width=edge_weights,
                                   alpha=0.6, edge_color=edge_colors, ax=ax)

        # 绘制标签
        labels = {n: f"{n}\n({node_counts[n]}篇)" if node_counts[n] > 1 else n for n in G.nodes()}
        nx.draw_networkx_labels(G, pos, labels=labels, font_size=self.config.label_fontsize - 1,
                                font_family=self.config.font_family, ax=ax)

        ax.set_title('基因-疾病关联网络图', fontsize=self.config.title_fontsize, fontweight='bold')

        # 添加图例
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor='#2ecc71',
                   markersize=15, label='基因'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='#e74c3c',
                   markersize=12, label='疾病'),
            Line2D([0], [0], color='#34495e', linewidth=3, label='强关联'),
            Line2D([0], [0], color='#95a5a6', linewidth=1, label='关联')
        ]
        ax.legend(handles=legend_elements, loc='upper left', framealpha=0.9)
        ax.axis('off')

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=self.config.dpi, bbox_inches='tight')
            plt.close()
            return output_path
        else:
            # 返回Base64编码
            buf = io.BytesIO()
            plt.savefig(buf, format='png', dpi=self.config.dpi, bbox_inches='tight')
            plt.close()
            buf.seek(0)
            return base64.b64encode(buf.read()).decode('utf-8')

    def plot_plotly(
        self,
        conclusions: List[GeneConclusion]
    ) -> Optional[str]:
        """
        使用 Plotly 绘制交互式网络图

        Args:
            conclusions: 结论列表

        Returns:
            HTML字符串
        """
        if not HAS_PLOTLY or not HAS_NETWORKX:
            logger.warning("缺少依赖库，无法生成交互式网络图")
            return None

        G, node_types, node_counts = self.build_network(conclusions)

        if len(G.nodes) == 0:
            return None

        # 节点数量过少时返回None
        if len(G.nodes) <= 2:
            return None

        # 优化布局参数
        k_value = max(1.5, 3.0 / (len(G.nodes) ** 0.3))
        pos = nx.spring_layout(G, k=k_value, iterations=100, seed=42)

        # 创建边追踪
        edge_x = []
        edge_y = []
        edge_weights_list = []
        for edge in G.edges(data=True):
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])
            edge_weights_list.append(edge[2]['weight'])

        edge_trace = go.Scatter(
            x=edge_x, y=edge_y,
            line=dict(width=2, color='#888'),
            hoverinfo='none',
            mode='lines'
        )

        # 创建节点追踪（分别处理基因和疾病）
        gene_x, gene_y, gene_text, gene_sizes = [], [], [], []
        disease_x, disease_y, disease_text, disease_sizes = [], [], [], []

        for node in G.nodes():
            x, y = pos[node]
            count = node_counts.get(node, 1)
            # 节点大小根据频次缩放
            size = 25 + min(30, count * 8)
            hover_text = f"{node}<br>{count}篇文献支持"

            if node in node_types["gene"]:
                gene_x.append(x)
                gene_y.append(y)
                gene_text.append(hover_text)
                gene_sizes.append(size)
            else:
                disease_x.append(x)
                disease_y.append(y)
                disease_text.append(hover_text)
                disease_sizes.append(size)

        # 基因节点
        gene_trace = go.Scatter(
            x=gene_x, y=gene_y,
            mode='markers+text',
            text=[t.split('<br>')[0] for t in gene_text],
            textposition="top center",
            textfont=dict(size=12, color='#2ecc71'),
            marker=dict(
                size=gene_sizes,
                color='#2ecc71',
                line=dict(width=2, color='#27ae60'),
                opacity=0.9
            ),
            hoverinfo='text',
            hovertext=gene_text,
            name='基因'
        )

        # 疾病节点
        disease_trace = go.Scatter(
            x=disease_x, y=disease_y,
            mode='markers+text',
            text=[t.split('<br>')[0] for t in disease_text],
            textposition="top center",
            textfont=dict(size=11, color='#e74c3c'),
            marker=dict(
                size=disease_sizes,
                color='#e74c3c',
                line=dict(width=2, color='#c0392b'),
                opacity=0.9
            ),
            hoverinfo='text',
            hovertext=disease_text,
            name='疾病'
        )

        fig = go.Figure(data=[edge_trace, gene_trace, disease_trace],
                        layout=go.Layout(
                            title=dict(text='基因-疾病关联网络图', x=0.5, font=dict(size=16)),
                            showlegend=True,
                            legend=dict(x=0, y=1),
                            hovermode='closest',
                            margin=dict(b=20, l=5, r=5, t=40),
                            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                            width=None,  # 响应式宽度
                            height=self.config.height,
                            autosize=True,
                            plot_bgcolor='white'
                        ))

        return fig.to_html(include_plotlyjs='cdn', full_html=False)

    def generate_table(
        self,
        conclusions: List[GeneConclusion]
    ) -> Optional[str]:
        """
        生成基因-疾病关联表格（数据量小时使用）

        Args:
            conclusions: 结论列表

        Returns:
            Markdown表格字符串
        """
        G, node_types, node_counts = self.build_network(conclusions)

        if len(G.nodes) == 0:
            return None

        # 收集关联数据
        relations = []
        for gene in node_types["gene"]:
            for disease in G.neighbors(gene):
                weight = G[gene][disease]['weight']
                relations.append({
                    'gene': gene,
                    'disease': disease,
                    'count': weight
                })

        if not relations:
            return None

        # 按关联文献数排序
        relations.sort(key=lambda x: x['count'], reverse=True)

        # 生成表格
        lines = ["| 基因 | 疾病 | 支持文献数 |", "| --- | --- | --- |"]
        for r in relations:
            lines.append(f"| {r['gene']} | {r['disease']} | {r['count']}篇 |")

        return "\n".join(lines)


class PathwayFrequencyChart:
    """信号通路频率柱状图"""

    def __init__(self, config: Optional[ChartConfig] = None):
        self.config = config or ChartConfig()

    def get_pathway_counts(
        self,
        conclusions: List[GeneConclusion]
    ) -> Counter:
        """
        统计信号通路出现频率（自动规范化通路名称）

        Args:
            conclusions: 结论列表

        Returns:
            通路计数器
        """
        pathway_counts = Counter()
        for c in conclusions:
            # 规范化通路名称并去重
            normalized_pathways = normalize_pathways(c.pathways)
            for pathway in normalized_pathways:
                if pathway and pathway not in ["未提及", "解析失败"]:
                    pathway_counts[pathway] += 1
        return pathway_counts

    def plot_matplotlib(
        self,
        conclusions: List[GeneConclusion],
        output_path: Optional[str] = None,
        top_n: int = 10
    ) -> Optional[str]:
        """
        使用 Matplotlib 绘制柱状图

        Args:
            conclusions: 结论列表
            output_path: 输出路径
            top_n: 显示前N个通路

        Returns:
            Base64编码的图片或文件路径
        """
        if not HAS_MATPLOTLIB:
            logger.warning("缺少 matplotlib，无法生成柱状图")
            return None

        # 设置中文字体
        setup_chinese_font()

        pathway_counts = self.get_pathway_counts(conclusions)

        if not pathway_counts:
            return None

        # 取前N个
        top_pathways = pathway_counts.most_common(top_n)
        if not top_pathways:
            return None

        pathways = [p[0] for p in top_pathways]
        counts = [p[1] for p in top_pathways]

        # 创建图形
        fig, ax = plt.subplots(figsize=(self.config.width / self.config.dpi,
                                          self.config.height / self.config.dpi))

        # 绘制水平柱状图
        colors = plt.cm.Blues([0.4 + 0.5 * i / len(pathways) for i in range(len(pathways))])
        bars = ax.barh(pathways, counts, color=colors[::-1])

        # 添加数值标签
        for bar, count in zip(bars, counts):
            ax.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height()/2,
                   f'{count}', va='center', fontsize=self.config.label_fontsize)

        ax.set_xlabel('文献数量', fontsize=self.config.label_fontsize + 2)
        ax.set_title('信号通路频率分布', fontsize=self.config.title_fontsize)
        ax.invert_yaxis()  # 最高的在上面

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=self.config.dpi, bbox_inches='tight')
            plt.close()
            return output_path
        else:
            buf = io.BytesIO()
            plt.savefig(buf, format='png', dpi=self.config.dpi, bbox_inches='tight')
            plt.close()
            buf.seek(0)
            return base64.b64encode(buf.read()).decode('utf-8')

    def plot_plotly(
        self,
        conclusions: List[GeneConclusion],
        top_n: int = 10
    ) -> Optional[str]:
        """
        使用 Plotly 绘制交互式柱状图

        Args:
            conclusions: 结论列表
            top_n: 显示前N个通路

        Returns:
            HTML字符串
        """
        if not HAS_PLOTLY:
            logger.warning("缺少 plotly，无法生成交互式柱状图")
            return None

        pathway_counts = self.get_pathway_counts(conclusions)

        if not pathway_counts:
            return None

        top_pathways = pathway_counts.most_common(top_n)
        if not top_pathways:
            return None

        pathways = [p[0] for p in top_pathways[::-1]]  # 反转使最高的在上面
        counts = [p[1] for p in top_pathways[::-1]]

        fig = go.Figure(data=[
            go.Bar(
                x=counts,
                y=pathways,
                orientation='h',
                marker=dict(
                    color=counts,
                    colorscale='Blues',
                    showscale=False
                ),
                text=counts,
                textposition='outside'
            )
        ])

        # 根据名称长度动态调整左边距
        max_name_len = max(len(p) for p in pathways) if pathways else 20
        left_margin = max(120, min(250, max_name_len * 8))

        fig.update_layout(
            xaxis_title='文献数量',
            width=None,  # 响应式宽度
            height=self.config.height,
            margin=dict(l=left_margin, r=30, t=30, b=30),
            yaxis=dict(tickfont=dict(size=11)),
            autosize=True  # 自动调整大小
        )

        return fig.to_html(include_plotlyjs='cdn', full_html=False)


class PublicationTimeChart:
    """文献发表时间分布图"""

    def __init__(self, config: Optional[ChartConfig] = None):
        self.config = config or ChartConfig()

    def get_year_counts(
        self,
        conclusions: List[GeneConclusion]
    ) -> Dict[int, int]:
        """
        统计发表年份分布

        Args:
            conclusions: 结论列表

        Returns:
            年份计数字典
        """
        year_counts = Counter()
        for c in conclusions:
            if c.year and c.year > 1900 and c.year <= 2030:
                year_counts[c.year] += 1
        return dict(sorted(year_counts.items()))

    def plot_matplotlib(
        self,
        conclusions: List[GeneConclusion],
        output_path: Optional[str] = None
    ) -> Optional[str]:
        """
        使用 Matplotlib 绘制时间分布图

        Args:
            conclusions: 结论列表
            output_path: 输出路径

        Returns:
            Base64编码的图片或文件路径
        """
        if not HAS_MATPLOTLIB:
            logger.warning("缺少 matplotlib，无法生成时间分布图")
            return None

        # 设置中文字体
        setup_chinese_font()

        year_counts = self.get_year_counts(conclusions)

        if not year_counts:
            return None

        years = list(year_counts.keys())
        counts = list(year_counts.values())

        # 创建图形
        fig, ax = plt.subplots(figsize=(self.config.width / self.config.dpi,
                                          self.config.height / self.config.dpi))

        # 绘制柱状图
        bars = ax.bar(years, counts, color='#2ecc71', alpha=0.8, edgecolor='#27ae60')

        # 添加趋势线
        if len(years) > 1:
            ax.plot(years, counts, 'o-', color='#e74c3c', linewidth=2, markersize=6)

        # 添加数值标签
        for bar, count in zip(bars, counts):
            if count > 0:
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                       f'{count}', ha='center', fontsize=self.config.label_fontsize - 1)

        ax.set_xlabel('发表年份', fontsize=self.config.label_fontsize + 2)
        ax.set_ylabel('文献数量', fontsize=self.config.label_fontsize + 2)
        ax.set_title('文献发表时间分布', fontsize=self.config.title_fontsize)

        # 设置x轴刻度
        if len(years) > 10:
            ax.set_xticks(years[::2])
        else:
            ax.set_xticks(years)
        ax.set_xticklabels([str(y) for y in (years[::2] if len(years) > 10 else years)],
                          rotation=45, ha='right')

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=self.config.dpi, bbox_inches='tight')
            plt.close()
            return output_path
        else:
            buf = io.BytesIO()
            plt.savefig(buf, format='png', dpi=self.config.dpi, bbox_inches='tight')
            plt.close()
            buf.seek(0)
            return base64.b64encode(buf.read()).decode('utf-8')

    def plot_plotly(
        self,
        conclusions: List[GeneConclusion]
    ) -> Optional[str]:
        """
        使用 Plotly 绘制交互式时间分布图

        Args:
            conclusions: 结论列表

        Returns:
            HTML字符串
        """
        if not HAS_PLOTLY:
            logger.warning("缺少 plotly，无法生成交互式时间分布图")
            return None

        year_counts = self.get_year_counts(conclusions)

        if not year_counts:
            return None

        years = list(year_counts.keys())
        counts = list(year_counts.values())

        fig = go.Figure()

        # 柱状图
        fig.add_trace(go.Bar(
            x=years,
            y=counts,
            marker=dict(color='#2ecc71'),
            name='文献数量'
        ))

        # 趋势线
        if len(years) > 1:
            fig.add_trace(go.Scatter(
                x=years,
                y=counts,
                mode='lines+markers',
                marker=dict(color='#e74c3c', size=8),
                line=dict(width=2),
                name='趋势'
            ))

        fig.update_layout(
            title='文献发表时间分布',
            xaxis_title='发表年份',
            yaxis_title='文献数量',
            width=None,  # 响应式宽度
            height=self.config.height,
            autosize=True,
            showlegend=True,
            legend=dict(orientation='h', yanchor='bottom', y=1.02)
        )

        return fig.to_html(include_plotlyjs='cdn', full_html=False)


class ConflictSummaryChart:
    """冲突统计图表"""

    def __init__(self, config: Optional[ChartConfig] = None):
        self.config = config or ChartConfig()

    def plot_matplotlib(
        self,
        conflict_report: ConflictReport,
        output_path: Optional[str] = None
    ) -> Optional[str]:
        """
        绘制冲突统计饼图

        Args:
            conflict_report: 冲突报告
            output_path: 输出路径

        Returns:
            Base64编码的图片或文件路径
        """
        if not HAS_MATPLOTLIB:
            return None

        # 设置中文字体
        setup_chinese_font()

        # 统计冲突类型
        conflict_types = Counter()
        for c in conflict_report.conflicts:
            conflict_types[c.topic] += 1

        if not conflict_types:
            # 无冲突时显示单一结果
            labels = ['无冲突']
            sizes = [1]
            colors = ['#2ecc71']
            explode = (0,)
        else:
            labels = list(conflict_types.keys())
            sizes = list(conflict_types.values())
            colors = ['#e74c3c', '#f39c12', '#3498db', '#9b59b6'][:len(labels)]
            explode = tuple(0.05 for _ in labels)

        fig, ax = plt.subplots(figsize=(self.config.width / self.config.dpi,
                                          self.config.height / self.config.dpi))

        wedges, texts, autotexts = ax.pie(
            sizes, explode=explode, labels=labels, colors=colors,
            autopct='%1.1f%%', startangle=90,
            textprops={'fontsize': self.config.label_fontsize}
        )

        ax.set_title('冲突类型分布', fontsize=self.config.title_fontsize)
        ax.axis('equal')

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=self.config.dpi, bbox_inches='tight')
            plt.close()
            return output_path
        else:
            buf = io.BytesIO()
            plt.savefig(buf, format='png', dpi=self.config.dpi, bbox_inches='tight')
            plt.close()
            buf.seek(0)
            return base64.b64encode(buf.read()).decode('utf-8')

    def plot_plotly(
        self,
        conflict_report: ConflictReport
    ) -> Optional[str]:
        """
        使用 Plotly 绘制交互式冲突统计图

        Args:
            conflict_report: 冲突报告

        Returns:
            HTML字符串
        """
        if not HAS_PLOTLY:
            logger.warning("缺少 plotly，无法生成交互式冲突图")
            return None

        # 统计冲突类型
        conflict_types = Counter()
        for c in conflict_report.conflicts:
            conflict_types[c.topic] += 1

        if not conflict_types:
            # 无冲突时显示单一结果
            labels = ['无冲突']
            values = [1]
            colors = ['#2ecc71']
        else:
            labels = list(conflict_types.keys())
            values = list(conflict_types.values())
            colors = ['#e74c3c', '#f39c12', '#3498db', '#9b59b6', '#1abc9c', '#e67e22'][:len(labels)]

        fig = go.Figure(data=[go.Pie(
            labels=labels,
            values=values,
            marker=dict(colors=colors),
            hole=0.3,
            textinfo='label+percent',
            hovertemplate='%{label}<br>数量: %{value}<br>占比: %{percent}<extra></extra>'
        )])

        title = '冲突类型分布' if conflict_types else '结论一致性分析'
        fig.update_layout(
            title=title,
            width=None,  # 响应式宽度
            height=self.config.height,
            autosize=True,
            showlegend=True
        )

        return fig.to_html(include_plotlyjs='cdn', full_html=False)


class MutationDistributionChart:
    """突变类型分布饼图 - 基于文献频次动态识别热点突变"""

    # 常见突变类型关键词映射
    MUTATION_TYPE_KEYWORDS = {
        "错义突变": ["missense", "错义突变", "错义"],
        "无义突变": ["nonsense", "无义突变", "无义"],
        "移码突变": ["frameshift", "移码", "移码突变"],
        "缺失突变": ["deletion", "del", "缺失", "缺失突变"],
        "插入突变": ["insertion", "ins", "插入", "插入突变"],
        "剪接突变": ["splice", "剪接", "剪接突变"],
        "融合基因": ["fusion", "融合", "chimeric", "融合基因"],
        "基因扩增": ["amplification", "amplif", "扩增", "基因扩增"],
        "过表达": ["overexpress", "过表达", "高表达"],
        "低表达": ["underexpress", "低表达", "down-regul"],
    }

    def __init__(self, config: Optional[ChartConfig] = None):
        self.config = config or ChartConfig()

    def get_mutation_counts(
        self,
        conclusions: List[GeneConclusion],
        hotspot_min_count: int = 2
    ) -> Tuple[Counter, Dict[str, List[str]], Dict[str, int]]:
        """
        统计突变类型分布，基于文献频次动态识别热点突变

        Args:
            conclusions: 结论列表
            hotspot_min_count: 热点突变最少出现次数阈值（建议2-3篇）

        Returns:
            (突变类型计数器, 热点突变详情, 具体位点频次统计)
        """
        mutation_counts = Counter()
        hotspot_details = {}  # 记录每个热点突变的来源
        site_frequency = Counter()  # 具体位点出现频次

        # 第一步：收集所有具体突变位点及其来源
        site_sources = {}  # site -> list of pmids

        for c in conclusions:
            # 优先使用结构化的 mutation_sites 字段
            if hasattr(c, 'mutation_sites') and c.mutation_sites:
                for site in c.mutation_sites:
                    site_frequency[site] += 1
                    if site not in site_sources:
                        site_sources[site] = []
                    site_sources[site].append(c.pmid)

            # 同时从 mutation_effects 文本中提取补充
            mutation_text = (c.mutation_effects or "")
            if hasattr(c, 'key_findings'):
                mutation_text += " " + " ".join(c.key_findings or [])

            # 检测一般突变类型（仅在没有具体位点时）
            if not (hasattr(c, 'mutation_sites') and c.mutation_sites):
                if mutation_text.strip() and mutation_text.strip() != "未提及":
                    matched = False
                    for mut_type, keywords in self.MUTATION_TYPE_KEYWORDS.items():
                        for kw in keywords:
                            if kw.lower() in mutation_text.lower():
                                mutation_counts[mut_type] += 1
                                matched = True
                                break
                        if matched:
                            break

                    if not matched:
                        mutation_counts["其他突变"] += 1

        # 第二步：基于频次识别热点突变
        # 按频次排序
        sorted_sites = site_frequency.most_common()

        for site, count in sorted_sites:
            if count >= hotspot_min_count:
                # 达到阈值，标记为热点突变
                label = f"热点突变({site})"
                mutation_counts[label] = count
                hotspot_details[label] = [f"PMID:{pmid}" for pmid in site_sources[site]]
            else:
                # 未达到阈值，归类为"其他错义突变"
                if "其他错义突变" not in mutation_counts:
                    mutation_counts["其他错义突变"] = 0
                mutation_counts["其他错义突变"] += count

        return mutation_counts, hotspot_details, dict(site_frequency)

    def plot_matplotlib(
        self,
        conclusions: List[GeneConclusion],
        output_path: Optional[str] = None,
        source_info: Optional[str] = None
    ) -> Optional[str]:
        """
        绘制突变类型分布饼图

        Args:
            conclusions: 结论列表
            output_path: 输出路径
            source_info: 数据来源信息

        Returns:
            Base64编码的图片或文件路径
        """
        if not HAS_MATPLOTLIB:
            logger.warning("缺少 matplotlib，无法生成饼图")
            return None

        mutation_counts, hotspot_details, site_frequency = self.get_mutation_counts(conclusions)

        # 如果没有突变数据，返回 None
        if not mutation_counts:
            return None

        # 设置中文字体
        setup_chinese_font()

        labels = list(mutation_counts.keys())
        sizes = list(mutation_counts.values())
        total = sum(sizes)

        # 颜色方案：热点突变用暖色调，其他用冷色调
        base_colors = [
            '#e74c3c', '#c0392b', '#a93226', '#922b21',  # 红色系（热点突变）
            '#3498db', '#2ecc71', '#f39c12', '#9b59b6',  # 其他颜色
            '#1abc9c', '#e67e22', '#34495e', '#95a5a6', '#d35400'
        ]
        colors = base_colors[:len(labels)]

        # 热点突变稍微突出
        explode = []
        for label in labels:
            if "热点突变" in label:
                explode.append(0.08)
            else:
                explode.append(0.02)

        fig, ax = plt.subplots(figsize=(self.config.width / self.config.dpi,
                                          self.config.height / self.config.dpi + 0.5))

        # 绘制饼图，添加百分比标签
        wedges, texts, autotexts = ax.pie(
            sizes, explode=explode, labels=labels, colors=colors,
            autopct=lambda pct: f'{pct:.1f}%\n({int(round(pct/100*total))})' if pct > 3 else '',
            startangle=90,
            textprops={'fontsize': self.config.label_fontsize},
            pctdistance=0.75
        )

        # 设置百分比文字样式
        for autotext in autotexts:
            autotext.set_fontsize(9)
            autotext.set_color('white')
            autotext.set_fontweight('bold')

        ax.set_title(f'突变类型分布 (共{total}例)', fontsize=self.config.title_fontsize, fontweight='bold')
        ax.axis('equal')

        # 添加数据来源注释
        if source_info is None:
            num_articles = len(set(c.pmid for c in conclusions if c.pmid))
            hotspot_threshold = 2
            source_info = f"数据来源：基于本次检索的 {num_articles} 篇文献统计 | 热点突变阈值：≥{hotspot_threshold}篇"

        fig.text(0.5, 0.02, source_info, ha='center', fontsize=9,
                 style='italic', color='gray')

        plt.tight_layout(rect=[0, 0.05, 1, 1])

        if output_path:
            plt.savefig(output_path, dpi=self.config.dpi, bbox_inches='tight')
            plt.close()
            return output_path
        else:
            buf = io.BytesIO()
            plt.savefig(buf, format='png', dpi=self.config.dpi, bbox_inches='tight')
            plt.close()
            buf.seek(0)
            return base64.b64encode(buf.read()).decode('utf-8')

    def plot_plotly(
        self,
        conclusions: List[GeneConclusion]
    ) -> Optional[str]:
        """
        使用 Plotly 绘制交互式饼图

        Args:
            conclusions: 结论列表

        Returns:
            HTML字符串
        """
        if not HAS_PLOTLY:
            return None

        mutation_counts, hotspot_details, site_frequency = self.get_mutation_counts(conclusions)

        if not mutation_counts:
            return None

        labels = list(mutation_counts.keys())
        values = list(mutation_counts.values())
        total = sum(values)

        # 颜色方案
        colors = [
            '#e74c3c', '#c0392b', '#a93226', '#922b21',
            '#3498db', '#2ecc71', '#f39c12', '#9b59b6',
            '#1abc9c', '#e67e22', '#34495e', '#95a5a6', '#d35400'
        ][:len(labels)]

        # 构建悬停文本（包含热点突变详情）
        hover_texts = []
        for label, value in zip(labels, values):
            pct = value / total * 100
            hover_text = f"{label}<br>数量: {value}<br>占比: {pct:.1f}%"
            if label in hotspot_details:
                pmids = hotspot_details[label][:3]  # 最多显示3个来源
                hover_text += f"<br>来源: {', '.join(pmids)}"
                if len(hotspot_details[label]) > 3:
                    hover_text += f" 等{len(hotspot_details[label])}篇"
            hover_texts.append(hover_text)

        fig = go.Figure(data=[
            go.Pie(
                labels=labels,
                values=values,
                marker=dict(colors=colors),
                textinfo='label+percent+value',
                texttemplate='%{label}<br>%{percent:.1%} (%{value})',
                hole=0.35,
                hovertext=hover_texts,
                hoverinfo='text'
            )
        ])

        # 数据来源注释
        num_articles = len(set(c.pmid for c in conclusions if c.pmid))
        hotspot_threshold = 2
        source_text = f"数据来源：基于本次检索的 {num_articles} 篇文献统计 | 热点突变阈值：≥{hotspot_threshold}篇"

        # 如果有具体位点频次数据，添加说明
        if site_frequency:
            top_sites = sorted(site_frequency.items(), key=lambda x: x[1], reverse=True)[:3]
            if top_sites:
                sites_str = ", ".join([f"{s}({c})" for s, c in top_sites])
                source_text += f"<br>高频位点: {sites_str}"

        fig.update_layout(
            title=dict(
                text=f'<b>突变类型分布</b><br><span style="font-size:12px;color:gray">{source_text}</span>',
                x=0.5
            ),
            font=dict(size=12),
            showlegend=True,
            legend=dict(
                orientation="h",
                yanchor="bottom",
                y=-0.2,
                xanchor="center",
                x=0.5
            )
        )

        return fig.to_html(include_plotlyjs='cdn', full_html=False)


class MutationDiseaseHeatmap:
    """突变-疾病关联热图"""

    def __init__(self, config: Optional[ChartConfig] = None):
        self.config = config or ChartConfig()

    def build_matrix(
        self,
        conclusions: List[GeneConclusion]
    ) -> Tuple[List[str], List[str], List[List[int]]]:
        """
        构建突变-疾病关联矩阵

        Args:
            conclusions: 结论列表

        Returns:
            (突变列表, 疾病列表, 矩阵数据)
        """
        matrix = {}  # {mutation: {disease: count}}

        for c in conclusions:
            if hasattr(c, 'mutation_disease_associations') and c.mutation_disease_associations:
                for assoc in c.mutation_disease_associations:
                    mutation = assoc.get("mutation", "")
                    disease = assoc.get("disease", "")

                    if mutation and disease:
                        # 规范化突变位点
                        from .helpers import normalize_mutation_site
                        mutation = normalize_mutation_site(mutation)

                        if mutation not in matrix:
                            matrix[mutation] = {}
                        if disease not in matrix[mutation]:
                            matrix[mutation][disease] = 0
                        matrix[mutation][disease] += 1

        if not matrix:
            return [], [], []

        # 收集所有疾病类型
        all_diseases = set()
        for diseases in matrix.values():
            all_diseases.update(diseases.keys())

        # 按出现次数排序
        disease_counts = Counter()
        for diseases in matrix.values():
            for disease, count in diseases.items():
                disease_counts[disease] += count
        sorted_diseases = [d for d, _ in disease_counts.most_common()]

        # 按总关联数排序突变
        mutation_totals = {
            mut: sum(diseases.values())
            for mut, diseases in matrix.items()
        }
        sorted_mutations = sorted(
            matrix.keys(),
            key=lambda m: mutation_totals[m],
            reverse=True
        )

        # 构建矩阵数据
        data = []
        for mutation in sorted_mutations:
            row = []
            for disease in sorted_diseases:
                row.append(matrix[mutation].get(disease, 0))
            data.append(row)

        return sorted_mutations, sorted_diseases, data

    def plot_plotly(
        self,
        conclusions: List[GeneConclusion]
    ) -> Optional[str]:
        """
        使用 Plotly 绘制热图

        Args:
            conclusions: 结论列表

        Returns:
            HTML字符串
        """
        if not HAS_PLOTLY:
            logger.warning("缺少 plotly，无法生成热图")
            return None

        mutations, diseases, data = self.build_matrix(conclusions)

        if not mutations or not diseases:
            return None

        # 数据量过大时限制显示
        max_mutations = 20
        max_diseases = 15

        if len(mutations) > max_mutations:
            mutations = mutations[:max_mutations]
            data = data[:max_mutations]

        if len(diseases) > max_diseases:
            diseases = diseases[:max_diseases]
            data = [row[:max_diseases] for row in data]

        fig = go.Figure(data=go.Heatmap(
            z=data,
            x=diseases,
            y=mutations,
            colorscale=[
                [0, '#f8f9fa'],
                [0.3, '#d4edda'],
                [0.6, '#ffc107'],
                [1, '#dc3545']
            ],
            hovertemplate='%{y} - %{x}<br>关联次数: %{z}<extra></extra>',
            showscale=True,
            colorbar=dict(title="关联次数")
        ))

        fig.update_layout(
            title=dict(
                text='突变-疾病关联热图',
                x=0.5,
                font=dict(size=16)
            ),
            xaxis_title='疾病/癌症类型',
            yaxis_title='突变位点',
            xaxis=dict(tickangle=-45),
            width=None,  # 响应式宽度
            height=max(400, len(mutations) * 35 + 100),
            autosize=True,
            margin=dict(l=80, r=20, t=60, b=100)
        )

        return fig.to_html(include_plotlyjs='cdn', full_html=False)

    def generate_table(
        self,
        conclusions: List[GeneConclusion]
    ) -> Optional[str]:
        """
        生成突变-疾病关联表格（数据量小时使用）

        Args:
            conclusions: 结论列表

        Returns:
            Markdown表格字符串
        """
        from ..core.report_generator import build_mutation_disease_matrix, format_mutation_disease_matrix

        matrix = build_mutation_disease_matrix(conclusions)

        if not matrix:
            return None

        return format_mutation_disease_matrix(matrix)


class ReportVisualizer:
    """报告可视化整合器"""

    def __init__(self, config: Optional[ChartConfig] = None):
        self.config = config or ChartConfig()
        self.network_chart = GeneDiseaseNetwork(config)
        self.pathway_chart = PathwayFrequencyChart(config)
        self.time_chart = PublicationTimeChart(config)
        self.conflict_chart = ConflictSummaryChart(config)
        self.mutation_chart = MutationDistributionChart(config)
        self.mutation_disease_heatmap = MutationDiseaseHeatmap(config)

    def generate_all_charts(
        self,
        conclusions: List[GeneConclusion],
        conflict_report: ConflictReport,
        output_format: str = "html"
    ) -> Dict[str, str]:
        """
        生成所有图表（默认使用交互式Plotly图表）

        Args:
            conclusions: 结论列表
            conflict_report: 冲突报告
            output_format: 输出格式 ("html" 交互式, "base64" 静态图片)

        Returns:
            图表字典
        """
        charts = {}

        # 1. 基因-疾病网络图
        try:
            network_chart = self.network_chart.plot_plotly(conclusions)
            if network_chart:
                charts["network"] = network_chart
            else:
                # 节点过少时使用表格
                network_table = self.network_chart.generate_table(conclusions)
                if network_table:
                    charts["network_table"] = network_table
        except Exception as e:
            logger.warning(f"生成网络图失败: {e}")
            try:
                network_table = self.network_chart.generate_table(conclusions)
                if network_table:
                    charts["network_table"] = network_table
            except Exception:
                pass

        # 2. 信号通路频率图
        try:
            charts["pathway"] = self.pathway_chart.plot_plotly(conclusions)
        except Exception as e:
            logger.warning(f"生成通路图失败: {e}")

        # 3. 发表时间分布图
        try:
            charts["timeline"] = self.time_chart.plot_plotly(conclusions)
        except Exception as e:
            logger.warning(f"生成时间图失败: {e}")

        # 4. 冲突统计图
        try:
            charts["conflict"] = self.conflict_chart.plot_plotly(conflict_report)
        except Exception as e:
            logger.warning(f"生成冲突图失败: {e}")

        # 5. 突变类型分布图
        try:
            charts["mutation"] = self.mutation_chart.plot_plotly(conclusions)
        except Exception as e:
            logger.warning(f"生成突变分布图失败: {e}")

        # 6. 突变-疾病关联热图
        try:
            heatmap_chart = self.mutation_disease_heatmap.plot_plotly(conclusions)
            if heatmap_chart:
                charts["mutation_disease_heatmap"] = heatmap_chart
            else:
                # 数据量小时使用表格
                heatmap_table = self.mutation_disease_heatmap.generate_table(conclusions)
                if heatmap_table:
                    charts["mutation_disease_table"] = heatmap_table
        except Exception as e:
            logger.warning(f"生成突变-疾病热图失败: {e}")

        return {k: v for k, v in charts.items() if v is not None}

    def generate_html_report(
        self,
        gene: str,
        conclusions: List[GeneConclusion],
        conflict_report: ConflictReport,
        markdown_content: str
    ) -> str:
        """
        生成包含图表的完整HTML报告

        Args:
            gene: 基因名
            conclusions: 结论列表
            conflict_report: 冲突报告
            markdown_content: Markdown内容

        Returns:
            HTML字符串
        """
        import markdown

        # 生成交互式图表
        charts = self.generate_all_charts(conclusions, conflict_report, output_format="html")

        # 转换Markdown
        html_content = markdown.markdown(
            markdown_content,
            extensions=['tables', 'toc']
        )

        # 构建图表HTML（网格布局）
        charts_html = ""

        # 信号通路
        if "pathway" in charts:
            charts_html += f"""
            <div class="chart-card">
                <h3>信号通路</h3>
                <div class="chart-wrapper">{charts['pathway']}</div>
            </div>
            """

        # 基因-疾病关联
        if "network" in charts:
            charts_html += f"""
            <div class="chart-card">
                <h3>基因-疾病关联</h3>
                <div class="chart-wrapper">{charts['network']}</div>
            </div>
            """
        elif "network_table" in charts:
            charts_html += f"""
            <div class="chart-card">
                <h3>基因-疾病关联</h3>
                {markdown.markdown(charts['network_table'], extensions=['tables'])}
            </div>
            """

        # 突变-疾病关联
        if "mutation_disease_heatmap" in charts:
            charts_html += f"""
            <div class="chart-card">
                <h3>突变-疾病关联</h3>
                <div class="chart-wrapper">{charts['mutation_disease_heatmap']}</div>
            </div>
            """
        elif "mutation_disease_table" in charts:
            charts_html += f"""
            <div class="chart-card">
                <h3>突变-疾病关联</h3>
                {markdown.markdown(charts['mutation_disease_table'], extensions=['tables'])}
            </div>
            """

        # 突变类型
        if "mutation" in charts:
            charts_html += f"""
            <div class="chart-card">
                <h3>突变类型分布</h3>
                <div class="chart-wrapper">{charts['mutation']}</div>
            </div>
            """

        # 发表时间
        if "timeline" in charts:
            charts_html += f"""
            <div class="chart-card">
                <h3>发表时间分布</h3>
                <div class="chart-wrapper">{charts['timeline']}</div>
            </div>
            """

        # 一致性分析
        if "conflict" in charts:
            charts_html += f"""
            <div class="chart-card">
                <h3>结论一致性</h3>
                <div class="chart-wrapper">{charts['conflict']}</div>
            </div>
            """

        if not charts:
            charts_html = """<div class="chart-card"><p style="text-align:center;color:#999;">暂无数据</p></div>"""

        # 完整HTML模板
        full_html = f"""
<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>基因文献分析报告: {gene}</title>
    <style>
        * {{ box-sizing: border-box; }}
        html {{ scroll-behavior: smooth; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            line-height: 1.6;
            color: #333;
            margin: 0;
            padding: 0;
            background: #f5f7fa;
            display: flex;
            min-height: 100vh;
        }}
        /* 左侧导航栏 */
        .sidebar {{
            position: fixed;
            left: 0;
            top: 0;
            width: 220px;
            height: 100vh;
            background: linear-gradient(180deg, #667eea 0%, #764ba2 100%);
            padding: 20px 0;
            overflow-y: auto;
            z-index: 100;
            flex-shrink: 0;
        }}
        .sidebar-header {{
            padding: 0 20px 20px;
            border-bottom: 1px solid rgba(255,255,255,0.2);
            margin-bottom: 20px;
        }}
        .sidebar-header h1 {{
            color: white;
            font-size: 1.3em;
            margin: 0 0 8px 0;
        }}
        .sidebar-header .meta {{
            color: rgba(255,255,255,0.8);
            font-size: 0.8em;
        }}
        .nav-title {{
            color: rgba(255,255,255,0.6);
            font-size: 0.75em;
            text-transform: uppercase;
            letter-spacing: 1px;
            padding: 10px 20px;
            margin: 0;
        }}
        .nav-list {{
            list-style: none;
            margin: 0;
            padding: 0;
        }}
        .nav-list li {{ margin: 2px 0; }}
        .nav-list a {{
            display: flex;
            align-items: center;
            padding: 12px 20px;
            color: rgba(255,255,255,0.85);
            text-decoration: none;
            font-size: 0.9em;
            transition: all 0.2s;
            border-left: 3px solid transparent;
        }}
        .nav-list a:hover {{
            background: rgba(255,255,255,0.1);
            color: white;
        }}
        .nav-list a.active {{
            background: rgba(255,255,255,0.15);
            color: white;
            border-left-color: #fff;
            font-weight: 500;
        }}
        .nav-list a .dot {{
            width: 6px;
            height: 6px;
            border-radius: 50%;
            background: rgba(255,255,255,0.5);
            margin-right: 12px;
            transition: all 0.2s;
        }}
        .nav-list a.active .dot {{
            background: white;
            box-shadow: 0 0 8px rgba(255,255,255,0.8);
        }}
        /* 主内容区 */
        .main-content {{
            margin-left: 220px;
            flex: 1;
            padding: 30px 40px;
            min-width: 0;
        }}
        /* 内容卡片 */
        .content {{
            background: white;
            border-radius: 12px;
            padding: 30px;
            margin-bottom: 30px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.08);
        }}
        h2 {{
            color: #333;
            border-left: 4px solid #667eea;
            padding-left: 15px;
            margin-top: 30px;
            margin-bottom: 20px;
        }}
        h3 {{ color: #555; margin-top: 25px; margin-bottom: 15px; }}
        /* 表格 */
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 15px 0;
            font-size: 0.95em;
        }}
        th, td {{
            padding: 12px 15px;
            text-align: left;
            border-bottom: 1px solid #eee;
        }}
        th {{ background: #f8f9fa; font-weight: 600; color: #333; }}
        tr:hover {{ background: #f8f9fa; }}
        /* 图表网格 */
        .charts-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(min(100%, 500px), 1fr));
            gap: 20px;
            margin-top: 20px;
        }}
        .chart-card {{
            background: white;
            border-radius: 12px;
            padding: 20px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.08);
            overflow: hidden;
        }}
        .chart-card h3 {{
            margin: 0 0 15px 0;
            padding-bottom: 10px;
            border-bottom: 1px solid #eee;
            font-size: 1em;
            color: #333;
        }}
        .chart-wrapper {{
            min-height: 350px;
            width: 100%;
            overflow: hidden;
        }}
        .chart-wrapper > div {{ width: 100% !important; }}
        .chart-wrapper svg {{ max-width: 100%; height: auto; }}
        /* 链接 */
        a {{ color: #667eea; text-decoration: none; }}
        a:hover {{ text-decoration: underline; }}
        .pmid-link {{
            display: inline-block;
            padding: 2px 10px;
            background: #667eea;
            color: white !important;
            border-radius: 12px;
            font-size: 0.85em;
        }}
        .pmid-link:hover {{ background: #5a6fd6; text-decoration: none; }}
        /* 标签 */
        .tag {{
            display: inline-block;
            padding: 2px 10px;
            border-radius: 12px;
            font-size: 0.85em;
            margin-right: 5px;
        }}
        .tag-wildtype {{ background: #d4edda; color: #155724; }}
        .tag-mutant {{ background: #f8d7da; color: #721c24; }}
        /* 警告框 */
        .warning {{
            background: #fff3cd;
            border-left: 4px solid #ffc107;
            padding: 15px;
            margin: 15px 0;
            border-radius: 0 8px 8px 0;
        }}
        .success {{
            background: #d4edda;
            border-left: 4px solid #28a745;
            padding: 15px;
            margin: 15px 0;
            border-radius: 0 8px 8px 0;
        }}
        /* 返回顶部 */
        .back-top {{
            position: fixed;
            bottom: 30px;
            right: 30px;
            width: 50px;
            height: 50px;
            background: #667eea;
            color: white;
            border: none;
            border-radius: 50%;
            cursor: pointer;
            font-size: 20px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.2);
            opacity: 0;
            transition: opacity 0.3s;
            z-index: 101;
        }}
        .back-top.show {{ opacity: 1; }}
        .back-top:hover {{ background: #5a6fd6; }}
        /* 页脚 */
        footer {{
            text-align: center;
            color: #999;
            padding: 30px;
            font-size: 0.9em;
        }}
        /* 响应式 */
        @media (max-width: 1024px) {{
            .sidebar {{
                width: 180px;
            }}
            .main-content {{
                margin-left: 180px;
                padding: 20px;
            }}
        }}
        @media (max-width: 768px) {{
            .sidebar {{
                transform: translateX(-100%);
                transition: transform 0.3s;
            }}
            .sidebar.open {{
                transform: translateX(0);
            }}
            .main-content {{
                margin-left: 0;
                padding: 15px;
            }}
            .content {{ padding: 20px; }}
            .charts-grid {{ grid-template-columns: 1fr; }}
            .back-top {{ bottom: 20px; right: 20px; width: 44px; height: 44px; }}
            /* 移动端菜单按钮 */
            .menu-toggle {{
                display: block;
                position: fixed;
                top: 15px;
                left: 15px;
                z-index: 102;
                background: #667eea;
                color: white;
                border: none;
                padding: 10px 12px;
                border-radius: 8px;
                cursor: pointer;
            }}
        }}
        @media (min-width: 769px) {{
            .menu-toggle {{ display: none; }}
        }}
        @media print {{
            .sidebar, .back-top, .menu-toggle {{ display: none; }}
            .main-content {{ margin-left: 0; }}
            .content {{ box-shadow: none; }}
        }}
    </style>
</head>
<body>
    <!-- 移动端菜单按钮 -->
    <button class="menu-toggle" id="menuToggle">☰</button>

    <!-- 左侧导航栏 -->
    <nav class="sidebar" id="sidebar">
        <div class="sidebar-header">
            <h1>{gene}</h1>
            <div class="meta">{datetime.now().strftime('%Y-%m-%d %H:%M')}</div>
        </div>
        <p class="nav-title">快速导航</p>
        <ul class="nav-list">
            <li><a href="#section-概述"><span class="dot"></span>概述</a></li>
            <li><a href="#section-基因功能语境"><span class="dot"></span>基因功能语境</a></li>
            <li><a href="#section-疾病与突变关联"><span class="dot"></span>疾病与突变关联</a></li>
            <li><a href="#section-信号通路"><span class="dot"></span>信号通路</a></li>
            <li><a href="#section-临床意义"><span class="dot"></span>临床意义</a></li>
            <li><a href="#section-突变效应"><span class="dot"></span>突变效应</a></li>
            <li><a href="#section-数据可视化"><span class="dot"></span>数据可视化</a></li>
            <li><a href="#section-参考文献"><span class="dot"></span>参考文献</a></li>
        </ul>
    </nav>

    <!-- 主内容区 -->
    <main class="main-content">
        <div class="content">
            {html_content}
        </div>

        <h2 id="section-数据可视化">数据可视化</h2>
        <div class="charts-grid">
            {charts_html}
        </div>

        <footer>本报告由基因文献智能检索助手自动生成</footer>
    </main>

    <button class="back-top" id="backTop" title="返回顶部">↑</button>

    <script>
        // 为h2标题自动添加id锚点
        document.querySelectorAll('.content h2').forEach(h2 => {{
            const text = h2.textContent.trim();
            h2.id = 'section-' + text;
        }});

        // 返回顶部
        const btn = document.getElementById('backTop');
        window.addEventListener('scroll', () => {{
            btn.classList.toggle('show', window.scrollY > 300);
        }});
        btn.addEventListener('click', () => {{
            window.scrollTo({{ top: 0, behavior: 'smooth' }});
        }});

        // 平滑滚动 + 导航高亮
        const navLinks = document.querySelectorAll('.nav-list a');

        navLinks.forEach(a => {{
            a.addEventListener('click', e => {{
                e.preventDefault();
                const target = document.querySelector(a.getAttribute('href'));
                if (target) {{
                    target.scrollIntoView({{ behavior: 'smooth', block: 'start' }});
                    // 移动端点击后关闭菜单
                    document.getElementById('sidebar').classList.remove('open');
                }}
            }});
        }});

        // 滚动时高亮当前section
        const sections = document.querySelectorAll('h2[id^="section-"]');

        function highlightNav() {{
            let current = '';
            const scrollPos = window.scrollY + 150;

            sections.forEach(section => {{
                const sectionTop = section.offsetTop;
                const sectionHeight = section.parentElement.offsetHeight;

                if (scrollPos >= sectionTop && scrollPos < sectionTop + sectionHeight + 200) {{
                    current = section.getAttribute('id');
                }}
            }});

            navLinks.forEach(link => {{
                link.classList.remove('active');
                if (link.getAttribute('href') === '#' + current) {{
                    link.classList.add('active');
                }}
            }});
        }}

        window.addEventListener('scroll', highlightNav);
        highlightNav(); // 初始化

        // 移动端菜单切换
        document.getElementById('menuToggle').addEventListener('click', () => {{
            document.getElementById('sidebar').classList.toggle('open');
        }});

        // 点击内容区关闭移动端菜单
        document.querySelector('.main-content').addEventListener('click', () => {{
            document.getElementById('sidebar').classList.remove('open');
        }});

        // PMID新窗口
        document.querySelectorAll('a[href*="pubmed"]').forEach(a => {{
            a.setAttribute('target', '_blank');
            a.setAttribute('rel', 'noopener');
        }});

        // Plotly图表响应式调整
        function resizePlotlyCharts() {{
            document.querySelectorAll('.chart-wrapper .plotly-graph-div').forEach(chart => {{
                if (window.Plotly) {{
                    Plotly.Plots.resize(chart);
                }}
            }});
        }}
        window.addEventListener('resize', resizePlotlyCharts);
        setTimeout(resizePlotlyCharts, 100);
    </script>
</body>
</html>
"""

        return full_html


if __name__ == "__main__":
    # 测试
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).parent.parent.parent))
    from gene_literature_agent.core.conclusion_extractor import GeneConclusion
    from gene_literature_agent.core.conflict_detector import ConflictReport, Conflict

    conclusions = [
        GeneConclusion(
            gene="TP53",
            pmid="12345678",
            title="TP53 mutations in lung cancer",
            disease_relation="突变型TP53促进肺癌进展",
            pathways=["DNA损伤应答", "细胞周期", "凋亡"],
            clinical_significance="与不良预后相关",
            mutation_effects="功能丧失性突变",
            key_findings=["发现1"],
            confidence="high",
            year=2023
        ),
        GeneConclusion(
            gene="TP53",
            pmid="23456789",
            title="TP53 in breast cancer",
            disease_relation="突变型TP53促进乳腺癌进展",
            pathways=["DNA损伤应答", "细胞周期"],
            clinical_significance="与不良预后相关",
            mutation_effects="功能丧失性突变",
            key_findings=["发现2"],
            confidence="medium",
            year=2022
        ),
        GeneConclusion(
            gene="TP53",
            pmid="34567890",
            title="TP53 pathway analysis",
            disease_relation="突变型TP53与多种癌症相关",
            pathways=["DNA损伤应答", "凋亡", "MAPK"],
            clinical_significance="预测化疗敏感性",
            mutation_effects="多种突变类型",
            key_findings=["发现3"],
            confidence="high",
            year=2024
        ),
    ]

    conflicts = [
        Conflict(
            topic="疾病关系",
            conclusion_a="促进肿瘤",
            pmid_a="12345678",
            conclusion_b="抑制肿瘤",
            pmid_b="45678901",
            severity="high"
        )
    ]

    report = ConflictReport(
        gene="TP53",
        conflicts=conflicts,
        consensus=["TP53是重要的抑癌基因"],
        needs_review=True,
        total_conclusions=3,
        conflict_ratio=0.33
    )

    visualizer = ReportVisualizer()
    charts = visualizer.generate_all_charts(conclusions, report)

    print("生成的图表:")
    for name, data in charts.items():
        if data:
            print(f"  - {name}: {len(data)} bytes")

    # 测试HTML报告生成
    md_content = "# 基因文献分析报告: TP53\n\n## 概述\n\nTP53是重要的抑癌基因。"
    html_report = visualizer.generate_html_report("TP53", conclusions, report, md_content)

    # 保存测试
    output_path = "/tmp/test_report.html"
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html_report)
    print(f"\nHTML报告已保存到: {output_path}")