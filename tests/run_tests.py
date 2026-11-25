#!/usr/bin/env python
"""tests 运行器 — 环境无 pytest(实测 2026-11-03), 全部测试为纯 assert 脚本(spec §1/§3.13)。

用法:  python -I tests/run_tests.py [test_module_filter ...]
发现:  tests/test_*.py 内所有 test_* 函数;逐个执行, 汇总 PASS/FAIL, 有 FAIL 退出码 1。
"""
from __future__ import annotations

import importlib
import os
import sys
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for p in (HERE, os.path.join(ROOT, "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)


def main(argv):
    mods = sorted(f[:-3] for f in os.listdir(HERE)
                  if f.startswith("test_") and f.endswith(".py"))
    if argv:
        mods = [m for m in mods if any(a in m for a in argv)]
    n_pass = n_fail = 0
    t0 = time.time()
    for mn in mods:
        m = importlib.import_module(mn)
        fns = sorted(f for f in dir(m) if f.startswith("test_") and callable(getattr(m, f)))
        for fn in fns:
            try:
                getattr(m, fn)()
                n_pass += 1
                print(f"PASS  {mn}.{fn}")
            except Exception:
                n_fail += 1
                print(f"FAIL  {mn}.{fn}")
                traceback.print_exc()
    print(f"\n{n_pass}/{n_pass + n_fail} 通过, {n_fail} 失败  ({time.time() - t0:.1f}s)")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
