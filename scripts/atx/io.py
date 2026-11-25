"""atx.io — 三种读取模式(implementation_spec §3.3;硬约束 4)。

(a) read_full      anndata 全载(Sepp 人/鼠、Aldinger);
(b) read_obs_var   只读 obs/var(h5py + anndata.io.read_elem), 不触矩阵(Castonguay 100 万核);
(c) read_rows_csr  按**行号**从 CSR 的 indptr/indices/data 切连续块读取(Castonguay raw/X 专用)。

借鉴(Depression_CB 已验证陷阱):backed 模式的**列**切片是 O(nnz) 全扫会卡死 → 永远先切行/
全载, 再在内存中切列。行切片读取器对 raw/X 做整数性断言(Aldinger/Castonguay 的 X 是
log-normalized 浮点, 计数在 raw/X —— 硬约束 4)。
"""
from __future__ import annotations

import h5py
import numpy as np
import scipy.sparse as sp


def read_full(path: str):
    """anndata 全载(Sepp/Aldinger;内存见 spec §10)。"""
    import anndata
    return anndata.read_h5ad(path)


def read_obs_var(path: str):
    """只读 obs/var(不触矩阵;百万核文件的元数据通道)。返回 (obs_df, var_df)。"""
    from anndata.io import read_elem
    with h5py.File(path, "r") as f:
        obs = read_elem(f["obs"])
        var = read_elem(f["var"])
    return obs, var


def h5_shape(path: str, layer: str = "X") -> tuple[int, int]:
    """CSR 数据集的形状(不读矩阵内容)。"""
    with h5py.File(path, "r") as f:
        g = f[layer]
        if isinstance(g, h5py.Group) and "shape" in g.attrs:
            shp = g.attrs["shape"]
        else:
            shp = g.attrs.get("shape", g.shape)
        return int(shp[0]), int(shp[1])


def assert_integer_counts(data: np.ndarray, layer: str, sample: int = 1_000_000) -> None:
    """raw/X 整数性断言(spec §3.3):计数层必须全为整数值。"""
    d = data[:sample] if data.size > sample else data
    if d.size and not np.all(d == np.floor(d)):
        bad = int(np.sum(d != np.floor(d)))
        raise AssertionError(f"{layer} 非整数计数({bad}/{len(d)} 抽样值非整)——检查是否误读了 log 层")


def read_rows_csr(path: str, rows, layer: str = "raw/X", chunk: int = 5000,
                  check_integer: bool | None = None):
    """按行号切片读 CSR(spec §3.3 签名)。返回 (csr_matrix[len(rows) × n_vars], var_index)。

    内部:rows 排序后分块, 每块读 [indptr[first], indptr[last]+1) 的**连续** indices/data 大段
    (避免逐行小读的 h5py 开销), 再按原 rows 顺序还原。check_integer 缺省 = layer 含 'raw' 时启用。
    """
    rows = np.asarray(rows)
    if not np.issubdtype(rows.dtype, np.integer):
        raise AssertionError(f"rows 必须是整数行号(dtype={rows.dtype});不接受布尔/浮点索引")
    if check_integer is None:
        check_integer = "raw" in layer
    from anndata.io import read_elem
    n_out = len(rows)
    idx_parts: list[np.ndarray] = []
    dat_parts: list[np.ndarray] = []
    out_indptr = np.zeros(n_out + 1, dtype=np.int64)
    with h5py.File(path, "r") as f:
        g = f[layer]
        enc = g.attrs.get("encoding-type")
        if enc is not None and str(enc) != "csr_matrix":
            raise AssertionError(f"{path}:{layer} encoding-type={enc}, 本读取器只支持 csr_matrix")
        n_vars = int(g.attrs["shape"][1])
        indptr_d, indices_d, data_d = g["indptr"], g["indices"], g["data"]
        n_total = indptr_d.shape[0] - 1
        if n_out and (rows.min() < 0 or rows.max() >= n_total):
            raise AssertionError(f"rows 越界: [{rows.min()},{rows.max()}] vs n_total={n_total}")
        order = np.argsort(rows, kind="stable")
        srows = rows[order]
        var = read_elem(f["var"])
        var_index = var.index
        pos = 0
        while pos < n_out:
            sel = srows[pos:pos + chunk]
            # indptr 也按连续段读(h5py 花式索引要求严格递增, 重复行号会拒绝 → 读回内存再取)
            r_lo, r_hi = int(sel[0]), int(sel[-1])
            ip = indptr_d[r_lo:r_hi + 2].astype(np.int64)
            lo, hi = int(ip[0]), int(ip[-1])
            slab_i = indices_d[lo:hi]
            slab_d = data_d[lo:hi]
            if check_integer:
                assert_integer_counts(slab_d, layer)
            off = sel - r_lo
            p0 = ip[off] - lo
            p1 = ip[off + 1] - lo
            for j in range(len(sel)):
                a, b = int(p0[j]), int(p1[j])
                idx_parts.append(slab_i[a:b])
                dat_parts.append(slab_d[a:b])
                out_indptr[pos + j + 1] = out_indptr[pos + j] + (b - a)
            pos += len(sel)
    if idx_parts:
        indices = np.concatenate(idx_parts)
        data = np.concatenate(dat_parts)
    else:
        indices = np.zeros(0, dtype=np.int64)
        data = np.zeros(0, dtype=np.float32)
    m_sorted = sp.csr_matrix((data, indices, out_indptr), shape=(n_out, n_vars))
    inv = np.empty(n_out, dtype=np.int64)
    inv[order] = np.arange(n_out)
    return m_sorted[inv], var_index


def row_sums_csr(path: str, layer: str = "X", block: int = 10_000_000) -> np.ndarray:
    """流式计算 CSR 每行和(内存 O(block);s00 的 size_factor 复验用, 不必全载矩阵)。"""
    with h5py.File(path, "r") as f:
        g = f[layer]
        n = int(g.attrs["shape"][0])
        indptr = g["indptr"][:].astype(np.int64)
        data_d = g["data"]
        counts = np.diff(indptr)
        out = np.zeros(n, dtype=np.float64)
        r0 = 0
        while r0 < n:
            cum = np.cumsum(counts[r0:])
            take = int(np.searchsorted(cum, block, side="left")) + 1
            r1 = min(n, r0 + take)
            lo, hi = int(indptr[r0]), int(indptr[r1])
            slab = data_d[lo:hi]
            rid = np.repeat(np.arange(r0, r1, dtype=np.int64), counts[r0:r1])
            out[r0:r1] = np.bincount(rid - r0, weights=slab, minlength=r1 - r0)
            r0 = r1
    return out
