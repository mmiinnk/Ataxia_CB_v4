"""atx.config — config.yaml 唯一入口(硬约束 6)与物种派生配置(implementation_spec §3.1)。

职责:加载 YAML;派生每物种 bin_order(1..n)/period_of/短标签;组三层(main/descriptive/reference);
results/v{version}/ 路径解析;gene_set.version 与 curation_log.md 冻结版本的一致性门。

借鉴(Depression_CB 已验证口径,不拷代码):参数集中 YAML;时间顺序显式声明(h5ad categories 乱序,
存储顺序不可信——人 development_stage 的 categories 实测乱序,一切以 time_order 为准)。

pandas 3 注意:本模块不触 DataFrame,纯 dict/标量。
"""
from __future__ import annotations

import copy
import os
import re

import yaml

# curation_log.md 头部形如: > 版本：**v0.9.1-draft**（2026-10-29；…）
_FROZEN_VERSION_RE = re.compile(r"版本：\*\*v?([0-9][^*（(]*)\*\*")


def frozen_version(curation_log_path: str) -> str:
    """从 curation_log.md 头部解析冻结版本号(不含前导 v)。"""
    with open(curation_log_path, encoding="utf-8") as fh:
        head = fh.read(4000)
    m = _FROZEN_VERSION_RE.search(head)
    if m is None:
        raise RuntimeError(f"无法从 {curation_log_path} 头部解析 `版本：**v…**`")
    return m.group(1).strip()


class SpeciesCfg:
    """单物种派生配置。属性按 spec §3.1:time_order/bin_order/period_of/labels/donor_col/time_col/pooled。

    label_of: 原始时间值(h5ad 存储值, 如 'adult (9 weeks)') → 短标签('9wk', 文件名安全, §11 未决 8)。
    bin_order: 短标签 → 1..n(h5ad categories 乱序, 顺序只认 time_order)。
    period_of: 短标签 → period 名(config periods 以原始值书写, 此处换算)。
    """

    def __init__(self, name: str, d: dict):
        self.name = name
        self.d = d
        self.time_col = d["time_column"]
        self.donor_col = d["donor_column"]
        self.disease_col = d.get("disease_column")
        self.exclude_disease_donors = bool(d.get("exclude_disease_donors", False))
        self.n_cells_expected = int(d["n_cells_expected"])
        self.n_vars_expected = int(d["n_vars_expected"])
        self.time_order = [str(x) for x in d["time_order"]]
        self.labels = [str(x) for x in d["time_labels"]]
        if len(self.time_order) != len(self.labels):
            raise ValueError(f"{name}: time_order 与 time_labels 长度不一致")
        self.label_of = dict(zip(self.time_order, self.labels))
        self.bin_order = {lab: i + 1 for i, lab in enumerate(self.labels)}
        self.period_of: dict[str, str] = {}
        for period, stages in d["periods"].items():
            for s in stages:
                s = str(s)
                if s not in self.label_of:
                    raise ValueError(f"{name}: periods.{period} 含 time_order 外的值 {s!r}")
                self.period_of[self.label_of[s]] = period
        if set(self.period_of) != set(self.labels):
            raise ValueError(f"{name}: periods 未覆盖全部 bin: 缺 {set(self.labels) - set(self.period_of)}")
        # pooled_timepoints 以原始值书写(鼠 E10.5/E11.5 原始值==标签); single_donor_bins 以标签书写
        self.pooled = {self.label_of.get(str(x), str(x)) for x in (d.get("pooled_timepoints") or [])}
        self.single_donor_bins = {str(x) for x in (d.get("single_donor_bins") or [])}

    def collapsed(self, mapping: dict[str, str]) -> "SpeciesCfg":
        """敏感性开关 collapse_adult_bins(§11 未决 7):把若干标签合并为一个 bin。

        合并 bin 的供体数 = 各来源 bin 供体数之和(人 42/44/46/52y 各 1 供体 → adult 4 供体),
        故合并产物从 single_donor_bins 移除;bin_order 按首次出现序重排 1..n。
        """
        c = copy.copy(self)
        new_labels: list[str] = []
        for lab in self.labels:
            nl = mapping.get(lab, lab)
            if nl not in new_labels:
                new_labels.append(nl)
        c.labels = new_labels
        c.bin_order = {lab: i + 1 for i, lab in enumerate(new_labels)}
        c.label_of = {raw: mapping.get(lab, lab) for raw, lab in self.label_of.items()}
        c.period_of = {}
        for lab, period in self.period_of.items():
            c.period_of[mapping.get(lab, lab)] = period
        c.pooled = {mapping.get(lab, lab) for lab in self.pooled}
        c.single_donor_bins = {mapping.get(lab, lab) for lab in self.single_donor_bins} - set(mapping.values())
        return c

    def period(self, label: str) -> str:
        return self.period_of[label]


class Cfg:
    """全局配置。raw = YAML dict;所有取值走属性/方法, 脚本不得绕过本类硬编码(硬约束 6)。"""

    def __init__(self, raw: dict, path: str | None = None):
        self.raw = raw
        self.path = path
        self.meta = raw["meta"]
        self.seed = int(raw["meta"]["seed"])
        self.paths = raw["paths"]
        self.gene_set = raw["gene_set"]
        self.version = str(self.gene_set["version"])
        # draft = 非 1.0 冻结版本 → 探索性口径, 产物目录与文件均带版本标(spec §2 注释, 硬约束 6)
        self.draft = self.version != "1.0"
        self.bins_cfg = raw["bins"]
        self.stats_cfg = raw["stats"]
        self.sensitivity = raw.get("sensitivity", {}) or {}
        self._sp: dict[str, SpeciesCfg] = {}
        # 敏感性变体的产物隔离目录(方案A): 非空时 results() 在其下再加一层,
        # 使 results/v{ver}/<subdir>/{logs,intermediate,tables,figures,reports} 与正式产物完全隔离,
        # 避免翻开关重跑覆盖 results/v{ver}/ 的冻结产物(spec §9 sensitivity 目标的设计意图)。
        self.results_subdir = ""

    # ---- 路径 ----
    @property
    def root(self) -> str:
        return self.paths["project_root"]

    def path_for(self, key: str) -> str:
        return self.paths[key]

    def h5ad(self, key: str) -> str:
        """species.h5ad_key → 绝对路径。"""
        return self.paths[key]

    def results(self, *parts: str, mkdir: bool = False) -> str:
        """results/v{version}/… 拼接;draft 版本目录名即带 -draft 标。

        mkdir=True 时:末段含 '.'(视为文件名)→ 只创建其父目录;否则整路径当目录创建。
        这样 results("logs", "s00_x.tsv", mkdir=True) 可直接给 to_csv 用。
        """
        sub = (self.results_subdir,) if self.results_subdir else ()
        p = os.path.join(self.root, "results", f"v{self.version}", *sub, *parts)
        if mkdir:
            target = os.path.dirname(p) if "." in os.path.basename(p) else p
            os.makedirs(target, exist_ok=True)
        return p

    def write_draft_marker(self) -> None:
        """draft 模式下在 results/v{version}/ 放置 DRAFT 标记文件(探索性口径声明)。"""
        r = self.results()
        if self.draft:
            # 版本目录名含 '.'(如 v0.9.1-draft)时 results(mkdir=True) 会误判为文件只建父目录,
            # 故此处显式建目录(冻结 v1.0 时本方法为 no-op, 不建任何目录)。
            os.makedirs(r, exist_ok=True)
            with open(os.path.join(r, "DRAFT"), "w", encoding="utf-8") as fh:
                fh.write(
                    f"gene_set.version={self.version}(未冻结)。本目录全部产物为探索性口径,\n"
                    "不得作为论文数值引用;冻结(v1.0)后重跑。见 genes/curation_log.md 与 config.yaml 头注。\n"
                )

    # ---- 物种 ----
    def species(self, sp: str) -> SpeciesCfg:
        if sp not in self._sp:
            self._sp[sp] = SpeciesCfg(sp, self.raw["species"][sp])
        return self._sp[sp]

    def effective_species(self, sp: str) -> SpeciesCfg:
        """应用敏感性开关后的物种配置(现仅 collapse_adult_bins; 人 42/44/46/52y → 'AD')。"""
        spcfg = self.species(sp)
        if sp == "human" and self.sens("collapse_adult_bins", False):
            spcfg = spcfg.collapsed({"42y": "AD", "44y": "AD", "46y": "AD", "52y": "AD"})
        return spcfg

    # ---- 组层级(spec §3.1:main/descriptive/reference 三层) ----
    def sens(self, key: str, default=None):
        return self.sensitivity.get(key, default)

    def main_groups(self) -> list[str]:
        """统计口径分母组(tau/fraction/富集)。include_ubc_in_main 开关在此生效。"""
        tiers = self.raw["group_tiers"]
        m = list(tiers["main"])
        if self.sens("include_ubc_in_main", False) and "UBC" not in m:
            m = m + ["UBC"]
        return m

    def descriptive_groups(self) -> list[str]:
        t = self.raw["group_tiers"]
        d = list(t.get("descriptive", []))
        if self.sens("include_ubc_in_main", False):
            d = [g for g in d if g != "UBC"]
        return d

    def reference_groups(self) -> list[str]:
        return list(self.raw["group_tiers"].get("reference", []))

    def plot_groups(self) -> list[str]:
        """出图组 = main + descriptive + reference(默认 8+UBC+ALL=10, 即 cube_bin 全骨架行组)。"""
        return self.main_groups() + self.descriptive_groups() + self.reference_groups()

    def group_exprs(self) -> dict[str, str | None]:
        return {g: v["expr"] if isinstance(v, dict) else v for g, v in self.raw["cell_groups"].items()}

    def trend_anchors(self, sp: str) -> list[dict]:
        """spec §9(4) 门清单条目四的趋势方向锚点(config species.<sp>.trend_anchors)。

        返回 [{gene, group, sign, source}];sign ∈ {'+','−'→'-' 归一};未设 → []。
        2026-11-05 用户裁决: 人鼠两侧均设;锚点自 s01 常量迁入 config(硬约束 6)。
        """
        raw = self.raw["species"][sp].get("trend_anchors") or []
        out = []
        for a in raw:
            sign = str(a["sign"]).strip().replace("−", "-")
            if sign not in ("+", "-"):
                raise ValueError(f"trend_anchors sign={a['sign']!r} 须为 '+'/'-'({sp})")
            out.append({"gene": str(a["gene"]), "group": str(a["group"]), "sign": sign,
                        "source": str(a.get("source", ""))})
        return out

    # ---- 版本门(硬约束 6) ----
    def check_version_gate(self) -> str:
        """config gene_set.version 必须等于 curation_log.md 冻结版本;不等即 RuntimeError。"""
        fv = frozen_version(os.path.join(self.root, "genes", "curation_log.md"))
        if fv != self.version:
            raise RuntimeError(
                f"版本门失败: config gene_set.version={self.version!r} != curation_log 冻结版本 v{fv!r}"
                "(硬约束 6;先对齐再跑)")
        return fv


def apply_sens_overrides(cfg: "Cfg", pairs: list[str] | None) -> "Cfg":
    """CLI `--sens KEY=VALUE`(可多次) → 覆盖 cfg.sensitivity[KEY]。

    敏感性变体不改 config.yaml(改了会污染正式口径且易忘回滚), 而是在进程内覆盖;
    值按 true/false → bool、整数 → int、浮点 → float、其余 → str 解析。
    必须在 cfg.effective_species()/cfg.results()/sens_on 计算**之前**调用。
    """
    for kv in pairs or []:
        if "=" not in kv:
            raise ValueError(f"--sens 需 KEY=VALUE 形式, 得 {kv!r}")
        k, v = (s.strip() for s in kv.split("=", 1))
        if not k:
            raise ValueError(f"--sens 键名为空: {kv!r}")
        lv = v.lower()
        if lv in ("true", "false"):
            cfg.sensitivity[k] = (lv == "true")
        else:
            try:
                cfg.sensitivity[k] = int(v)
            except ValueError:
                try:
                    cfg.sensitivity[k] = float(v)
                except ValueError:
                    cfg.sensitivity[k] = v
        if k not in cfg.raw.get("sensitivity", {}):
            print(f"[sens] 注意: {k!r} 不在 config sensitivity 段内(新增键, 仅本次运行生效)")
    return cfg


def load(path: str) -> Cfg:
    """加载 config.yaml 并做最小结构校验(spec §3.1 签名: load(path) -> Cfg)。"""
    with open(path, encoding="utf-8") as fh:
        raw = yaml.safe_load(fh)
    for key in ("meta", "paths", "gene_set", "species", "cell_groups", "group_tiers",
                "bins", "stats", "state_alignment", "external", "figures"):
        if key not in raw:
            raise ValueError(f"config.yaml 缺顶层键: {key}")
    for sp in ("human", "mouse"):
        if sp not in raw["species"]:
            raise ValueError(f"config.yaml species 缺 {sp}")
    cfg = Cfg(raw, path)
    # 结构自检:组三层无重叠、plot_groups 含 ALL
    tiers = raw["group_tiers"]
    for a, b in (("main", "descriptive"), ("main", "reference"), ("descriptive", "reference")):
        ov = set(tiers[a]) & set(tiers[b])
        if ov:
            raise ValueError(f"group_tiers {a}/{b} 重叠: {ov}")
    if "ALL" not in cfg.reference_groups():
        raise ValueError("group_tiers.reference 必须含 ALL")
    return cfg
