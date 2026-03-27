"""
冲突检测提示词模板
Conflict Detection Prompt Template
"""

CONFLICT_DETECTION_PROMPT = """你是一个专业的文献分析助手。请分析以下关于基因 {gene} 的多篇文献结论，找出其中的矛盾或冲突之处。

文献结论列表：
{conclusions}

请输出JSON格式：
{{
  "conflicts": [
    {{
      "topic": "冲突主题",
      "conclusion_a": "结论A内容",
      "pmid_a": "PMID",
      "conclusion_b": "结论B内容",
      "pmid_b": "PMID",
      "resolution_suggestion": "可能的解释或建议"
    }}
  ],
  "consensus": ["共识结论列表"],
  "needs_review": true/false
}}"""
