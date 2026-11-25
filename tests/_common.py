"""tests 公共入口:真实 config 加载(测试用真参数, 不另造第二事实来源)。"""
from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SCRIPTS = os.path.join(ROOT, "scripts")
for p in (HERE, SCRIPTS):
    if p not in sys.path:
        sys.path.insert(0, p)

from atx import config as atxconfig  # noqa: E402

CFG_PATH = os.path.join(ROOT, "config", "config.yaml")
CURATION_LOG = os.path.join(ROOT, "genes", "curation_log.md")


def cfg():
    return atxconfig.load(CFG_PATH)
