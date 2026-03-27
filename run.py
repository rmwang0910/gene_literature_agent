#!/usr/bin/env python
"""
基因文献智能检索助手 - 独立运行脚本
用法: python run.py BRAF
"""
import os
import sys
from pathlib import Path

# 确保从正确的位置导入包
# 添加父目录到路径的最前面
parent_dir = str(Path(__file__).parent.parent)
if parent_dir in sys.path:
    sys.path.remove(parent_dir)
sys.path.insert(0, parent_dir)

# 移除可能存在的当前目录，避免相对导入冲突
current_dir = str(Path(__file__).parent)
if current_dir in sys.path:
    sys.path.remove(current_dir)

# 清理已导入的模块缓存（如果有的话）
modules_to_remove = [key for key in sys.modules.keys()
                     if key.startswith('gene_literature_agent')
                     or key in ('config', 'core', 'utils', 'prompts', 'constants', 'datasources')]
for mod in modules_to_remove:
    del sys.modules[mod]

from gene_literature_agent.cli import main

if __name__ == "__main__":
    main()
