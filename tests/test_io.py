"""io.py 测试(spec §3.13 test_io):read_rows_csr 与全载切行一致;顺序还原;整数断言;
row_sums_csr;read_obs_var。合成小文件写临时 h5ad(不在数据目录写任何文件, 硬约束 1)。"""
from __future__ import annotations

import os
import tempfile

import numpy as np
import pandas as pd
import scipy.sparse as sp

from atx import io as atxio

N, G = 50, 10


def _make_tmp() -> tuple[str, np.ndarray]:
    import anndata
    rng = np.random.default_rng(20261013)
    counts = rng.integers(0, 20, size=(N, G)).astype(np.float32)
    tot = counts.sum(1, keepdims=True)
    logx = np.log1p(counts / np.where(tot == 0, 1, tot) * 1e4).astype(np.float32)
    obs = pd.DataFrame({"donor_id": [f"d{i%3}" for i in range(N)]}, index=[f"c{i}" for i in range(N)])
    var = pd.DataFrame({"feature_name": [f"G{i}" for i in range(G)]},
                       index=[f"ENSG{i:04d}" for i in range(G)])
    ad = anndata.AnnData(X=sp.csr_matrix(logx), obs=obs, var=var)
    ad.raw = anndata.AnnData(X=sp.csr_matrix(counts), obs=obs.copy(), var=var.copy())
    fd, tmp = tempfile.mkstemp(suffix=".h5ad")
    os.close(fd)
    ad.write_h5ad(tmp)
    return tmp, counts


def test_read_rows_csr_matches_full():
    tmp, counts = _make_tmp()
    try:
        rows = [7, 2, 49, 0, 3, 3, 1]        # 乱序 + 重复
        m, var_index = atxio.read_rows_csr(tmp, rows, layer="raw/X")
        assert sp.issparse(m) and m.shape == (len(rows), G)
        got = m.toarray()
        exp = counts[np.asarray(rows)]
        assert np.allclose(got, exp), "read_rows_csr != 全载后按同序切行"
        assert list(var_index) == [f"ENSG{i:04d}" for i in range(G)]
        # 大 chunk 跨越内部块的连续段读取:全行乱序
        rng = np.random.default_rng(1)
        all_rows = rng.permutation(N)
        m2, _ = atxio.read_rows_csr(tmp, all_rows, layer="raw/X", chunk=7)
        assert np.allclose(m2.toarray(), counts[all_rows])
    finally:
        os.unlink(tmp)


def test_read_rows_csr_integer_assert():
    tmp, counts = _make_tmp()
    try:
        # X 是 log 浮点层 → 显式 check_integer 必须拒绝(硬约束 4:计数在 raw/X)
        try:
            atxio.read_rows_csr(tmp, [0, 1, 2], layer="X", check_integer=True)
            raise AssertionError("log 层未触发整数断言")
        except AssertionError as e:
            assert "非整数" in str(e)
        # raw/X 默认启用整数断言 → 应通过
        m, _ = atxio.read_rows_csr(tmp, [0, 1, 2], layer="raw/X")
        assert np.allclose(m.toarray(), counts[:3])
        # 非整数行号 → 拒绝
        try:
            atxio.read_rows_csr(tmp, np.array([0.0, 1.0]), layer="raw/X")
            raise AssertionError("浮点行号未被拒绝")
        except AssertionError as e:
            assert "整数行号" in str(e)
        # 越界 → 拒绝
        try:
            atxio.read_rows_csr(tmp, [0, N], layer="raw/X")
            raise AssertionError("越界行号未被拒绝")
        except AssertionError:
            pass
    finally:
        os.unlink(tmp)


def test_row_sums_and_obs_var():
    tmp, counts = _make_tmp()
    try:
        import anndata
        ad = anndata.read_h5ad(tmp)
        rs = atxio.row_sums_csr(tmp, "X", block=37)   # 小 block 强制多轮分块
        exp = np.asarray(ad.X.sum(axis=1)).ravel()
        assert np.allclose(rs, exp, rtol=1e-5), "row_sums_csr != 全载行和"
        assert atxio.h5_shape(tmp, "X") == (N, G)
        assert atxio.h5_shape(tmp, "raw/X") == (N, G)
        obs, var = atxio.read_obs_var(tmp)
        assert list(obs.index) == [f"c{i}" for i in range(N)]
        assert list(var["feature_name"]) == [f"G{i}" for i in range(G)]
    finally:
        os.unlink(tmp)
