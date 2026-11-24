"""Stage 1: load ST000390 (GC-TOF), build tidy tables, run integrity checks and a QC report.

No values are transformed, filtered or corrected here. Log2 + autoscaling is used
only for the QC PCA plot.

Outputs
  data/processed/intensities.csv   87 samples x 183 metabolites, raw counts as deposited
  data/processed/samples.csv       sample metadata, inferred patient, pair status, CV group
  data/processed/features.csv      metabolite annotation
  figures/01_*.png                 QC figures
  reports/01_qc.txt                the printed QC report
"""
import contextlib
import io
import json
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT = DATA / "processed"
FIG = ROOT / "figures"
REPORT = ROOT / "reports" / "01_qc.txt"
MWTAB = DATA / "ST000390_AN000626_mwtab.txt"
REST_DATA = DATA / "ST000390_data.json"

# Reference palette (dataviz skill), light mode; first three slots validate all-pairs.
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
COLORS = {"normal": "#2a78d6", "tumor": "#eb6834", "pool": "#1baf7a"}
MARKERS = {"normal": "o", "tumor": "s", "pool": "D"}


def block(lines, start, end):
    s = next(i for i, l in enumerate(lines) if l.startswith(start)) + 1
    e = lines.index(end)
    return pd.read_csv(io.StringIO("\n".join(lines[s:e])), sep="\t", dtype=str)


def norm_text(s):
    return re.sub(r"\s+", "", s.lower())


def load_mwtab(path):
    lines = path.read_text().splitlines()

    rows = []
    for ln in lines:
        if not ln.startswith("SUBJECT_SAMPLE_FACTORS"):
            continue
        p = ln.split("\t")
        r = {"sample": p[2].strip(), "subject": p[1].strip()}
        for kv in p[3].split("|"):
            k, _, v = kv.partition(":")
            r[k.strip()] = v.strip()
        for kv in p[4].split(";"):
            k, _, v = kv.partition("=")
            r[k.strip()] = v.strip()
        rows.append(r)
    samples = pd.DataFrame(rows).rename(columns={
        "Sex": "sex", "Smoking Status": "smoking", "Diagnosis": "diagnosis",
        "Tumor Present": "tumor_present", "Stage": "stage", "Age": "age", "Race": "race"})
    samples = samples.replace("-", np.nan)

    m = block(lines, "MS_METABOLITE_DATA_START", "MS_METABOLITE_DATA_END")
    m = m[m.iloc[:, 0] != "Factors"].set_index(m.columns[0])
    intens = m.apply(pd.to_numeric).T
    intens.index.name, intens.columns.name = "sample", None

    feats = block(lines, "METABOLITES_START", "METABOLITES_END").rename(columns={
        "metabolite_name": "name", "Retention Index": "retention_index",
        "Quantified m/z": "quant_mz", "KEGG ID": "kegg_id", "PubChem ID": "pubchem_id"})
    return samples, intens, feats


def infer_pairs(samples):
    """Group biological samples on all clinical fields; 1 tumor + 1 normal = one patient."""
    s = samples.copy()
    s["group"] = s.tumor_present.map({"Yes": "tumor", "No": "normal"}).fillna("pool")
    bio = s.group != "pool"
    fields = ["sex", "smoking", "diagnosis", "stage", "age", "race"]
    s.loc[bio, "_key"] = s.loc[bio, fields].astype(str).apply(lambda r: "|".join(map(norm_text, r)), axis=1)

    s["patient_id"], s["pair_status"], s["cv_group"] = pd.NA, "none", pd.NA
    n = 0
    # Number patients in order of their earliest sample ID so IDs are reproducible.
    for key, g in sorted(s[bio].groupby("_key"), key=lambda kg: kg[1]["sample"].min()):
        if sorted(g.group) == ["normal", "tumor"]:
            n += 1
            s.loc[g.index, ["patient_id", "pair_status", "cv_group"]] = [f"P{n:02d}", "unique", f"P{n:02d}"]
        else:
            s.loc[g.index, ["pair_status", "cv_group"]] = ["ambiguous", "AMB"]
    return s.drop(columns="_key")


def run_order(sample_ids):
    """Inferred from IDs like 121129bbdsa05_1: run date + position. Not confirmed by metadata."""
    parts = sample_ids.str.extract(r"^(\d{6})bbdsa(\d+)_")
    day = parts[0].map({"121129": "day1", "121130": "day2"})
    pos = parts[1].astype(int)
    order = (day == "day2") * 1000 + pos
    return day, pos, order.rank(method="first").astype(int)


class Tee:
    """Write to several streams at once (terminal + report file)."""
    def __init__(self, *streams):
        self.streams = streams

    def write(self, s):
        for st in self.streams:
            st.write(s)

    def flush(self):
        for st in self.streams:
            st.flush()


def check(cond, msg):
    if not cond:
        raise SystemExit(f"Integrity check failed: {msg}")
    print(f"  ok  {msg}")


def style(ax):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.grid(color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(exist_ok=True)
    plt.rcParams.update({"figure.facecolor": SURFACE, "text.color": INK, "axes.labelcolor": INK2,
                         "axes.titlesize": 11, "axes.labelsize": 9, "font.size": 9})

    print("ST000390 (AN000626, GC-TOF) stage 1 QC report")
    print("Generated by scripts/01_load_qc.py from data/ST000390_AN000626_mwtab.txt\n")

    samples, intens, feats = load_mwtab(MWTAB)
    samples = infer_pairs(samples)
    samples["run_day"], samples["run_pos"], samples["run_order"] = run_order(samples["sample"])
    samples = samples.set_index("sample")
    intens = intens.loc[samples.index]

    rest = json.loads(REST_DATA.read_text())
    rest_meta = pd.DataFrame([{"name": r["metabolite_name"], "metabolite_id": r["metabolite_id"],
                               "refmet_name": r["refmet_name"]} for r in rest.values()])
    feats = feats.merge(rest_meta, on="name", how="left")
    rest_mat = pd.DataFrame({r["metabolite_name"]: r["DATA"] for r in rest.values()}).apply(pd.to_numeric)

    # --- 1. Integrity checks -------------------------------------------------
    print("Integrity checks")
    counts = samples.group.value_counts().to_dict()
    check(intens.shape == (87, 183), f"matrix is 87 x 183 (got {intens.shape})")
    check(counts == {"tumor": 39, "normal": 39, "pool": 9}, f"39 tumor / 39 normal / 9 pool (got {counts})")
    check(samples.patient_id.nunique() == 36 and (samples.pair_status == "unique").sum() == 72,
          "36 uniquely matched pairs (72 samples)")
    check((samples.pair_status == "ambiguous").sum() == 6 and samples.loc[samples.cv_group == "AMB", "group"]
          .value_counts().to_dict() == {"tumor": 3, "normal": 3}, "one ambiguous group of 3 tumor + 3 normal")
    check(intens.notna().all().all() and (intens > 0).all().all(), "no missing values, zeros or negatives")
    check(feats.name.is_unique and set(feats.name) == set(intens.columns), "feature table matches matrix columns")
    check(rest_mat.reindex(index=intens.index, columns=intens.columns).equals(intens),
          "mwTab values identical to REST data endpoint")
    check(samples.run_day.notna().all() and samples.run_order.is_unique, "run day/position parsed for all samples")

    samples.to_csv(OUT / "samples.csv")
    intens.to_csv(OUT / "intensities.csv")
    feats.to_csv(OUT / "features.csv", index=False)
    print(f"  wrote {OUT.relative_to(ROOT)}/{{samples,intensities,features}}.csv")

    grp = samples.group
    bio = grp != "pool"

    # --- 2. Normalisation ----------------------------------------------------
    tot = intens.sum(axis=1)
    print("\nTotal signal per sample")
    print(f"  min {tot.min():,.0f}  median {tot.median():,.0f}  max {tot.max():,.0f}  "
          f"(range {100 * (tot.max() - tot.min()) / tot.median():.4f}% of median)")

    # --- 3. Pool variability -------------------------------------------------
    cv = lambda d: 100 * d.std(ddof=1) / d.mean()
    pool_cv = cv(intens[grp == "pool"])
    bio_cv = cv(intens[bio])
    feats_qc = pd.DataFrame({"pool_cv": pool_cv, "bio_cv": bio_cv, "cv_ratio": pool_cv / bio_cv})
    print("\nPool CV per metabolite (%)")
    print(f"  median {pool_cv.median():.1f}  IQR {pool_cv.quantile(.25):.1f}-{pool_cv.quantile(.75):.1f}  "
          f"max {pool_cv.max():.1f}  | biological median {bio_cv.median():.1f}")
    for t in (20, 30, 50):
        print(f"  pool CV > {t}%: {(pool_cv > t).sum()} metabolites")
    print(f"  pool CV >= biological CV: {(feats_qc.cv_ratio >= 1).sum()} metabolites")
    print("  highest pool CV:", ", ".join(f"{k} ({v:.0f}%)" for k, v in pool_cv.nlargest(8).items()))

    # --- 4. Run batches ------------------------------------------------------
    print("\nSamples per run day (day inferred from sample ID)")
    print(pd.crosstab(samples.run_day, grp).to_string())
    split = samples[samples.pair_status == "unique"].groupby("patient_id").run_day.nunique()
    print(f"  unique pairs split across days: {(split == 2).sum()} of {len(split)}")
    tab = pd.crosstab(samples.run_day[bio], grp[bio])
    print(f"  tumor/normal vs day, Fisher exact p = {stats.fisher_exact(tab.values)[1]:.3f}")
    rho, p = stats.spearmanr(samples.run_order[bio], (grp[bio] == "tumor").astype(int))
    print(f"  tumor status vs run order, Spearman rho = {rho:.2f} (p = {p:.2f})")

    # --- 5. Drift in pools ---------------------------------------------------
    pools = samples[grp == "pool"].sort_values("run_order")
    print("\nPools in run order:", ", ".join(f"{d}:{p}" for d, p in zip(pools.run_day, pools.run_pos)))
    lp = np.log2(intens.loc[pools.index])
    drift = lp.apply(lambda c: stats.spearmanr(pools.run_order, c)[0])
    print(f"  per-metabolite Spearman rho (pool intensity vs run order): median {drift.median():.2f}, "
          f"|rho| >= 0.8 in {(drift.abs() >= 0.8).sum()} metabolites")
    day_shift = (lp[pools.run_day == "day2"].mean() - lp[pools.run_day == "day1"].mean())
    print(f"  pool day2 - day1 log2 shift: median {day_shift.median():.2f}, "
          f"|shift| > 0.5 (1.4-fold) in {(day_shift.abs() > 0.5).sum()} metabolites")

    # --- 6. PCA (QC view only: log2 + autoscale) -----------------------------
    X = np.log2(intens)
    X = (X - X.mean()) / X.std(ddof=1)
    U, S, Vt = np.linalg.svd(X.values, full_matrices=False)
    scores = pd.DataFrame(U[:, :5] * S[:5], index=X.index, columns=[f"PC{i}" for i in range(1, 6)])
    var = 100 * S ** 2 / (S ** 2).sum()
    print(f"\nPCA (log2, autoscaled, all 87): PC1 {var[0]:.1f}%, PC2 {var[1]:.1f}%, PC3 {var[2]:.1f}%")
    pool_spread = scores.loc[grp == "pool", ["PC1", "PC2"]].std().mean()
    bio_spread = scores.loc[bio, ["PC1", "PC2"]].std().mean()
    print(f"  PC1/PC2 spread: pools {pool_spread:.2f} vs biological {bio_spread:.2f}")
    for pc in ("PC1", "PC2"):
        r, _ = stats.pearsonr(scores.loc[bio, pc], (grp[bio] == "tumor").astype(int))
        print(f"  {pc} correlation with tumor status: r = {r:.2f}")

    # Outliers: Hotelling T2 on the first k PCs (biological samples), F-distribution limit.
    k = 5
    Tb = scores.loc[bio].iloc[:, :k]
    n = len(Tb)
    t2 = ((Tb / Tb.std(ddof=1)) ** 2).sum(axis=1)
    lim95, lim99 = [k * (n - 1) * (n + 1) / (n * (n - k)) * stats.f.ppf(q, k, n - k) for q in (.95, .99)]
    print(f"  Hotelling T2 (5 PCs) limits: 95% {lim95:.1f}, 99% {lim99:.1f}")
    out = t2[t2 > lim95].sort_values(ascending=False)
    for sid, v in out.items():
        print(f"    {sid}  T2 = {v:.1f}  ({samples.loc[sid, 'group']}, {samples.loc[sid, 'cv_group']})"
              f"{'  > 99%' if v > lim99 else ''}")
    if out.empty:
        print("    none above 95%")
    samples.assign(t2=t2).to_csv(OUT / "samples.csv")
    feats.set_index("name").join(feats_qc).to_csv(OUT / "features.csv")

    # --- Figures -------------------------------------------------------------
    def legend(ax):
        ax.legend(frameon=False, fontsize=9, labelcolor=INK2, loc="best")

    fig, ax = plt.subplots(figsize=(6.4, 5))
    style(ax)
    for g in ("normal", "tumor", "pool"):
        d = scores[grp == g]
        ax.scatter(d.PC1, d.PC2, s=42, marker=MARKERS[g], color=COLORS[g], edgecolor=SURFACE,
                   linewidth=1.5, label=f"{g} (n={len(d)})", zorder=3)
    for sid in out.index:
        if t2[sid] > lim99:
            ax.annotate(sid, scores.loc[sid, ["PC1", "PC2"]], xytext=(5, 5), textcoords="offset points",
                        fontsize=8, color=INK2)
    ax.set(xlabel=f"PC1 ({var[0]:.1f}%)", ylabel=f"PC2 ({var[1]:.1f}%)")
    ax.set_title("ST000390 PCA, log2 autoscaled (QC view)", loc="left", color=INK)
    legend(ax)
    fig.tight_layout()
    fig.savefig(FIG / "01_pca_group.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 5))
    style(ax)
    for day, col, mk in (("day1", "#2a78d6", "o"), ("day2", "#eb6834", "s")):
        d = scores[samples.run_day == day]
        ax.scatter(d.PC1, d.PC2, s=42, marker=mk, color=col, edgecolor=SURFACE, linewidth=1.5,
                   label=f"{day} (n={len(d)})", zorder=3)
    ax.set(xlabel=f"PC1 ({var[0]:.1f}%)", ylabel=f"PC2 ({var[1]:.1f}%)")
    ax.set_title("ST000390 PCA by inferred run day", loc="left", color=INK)
    legend(ax)
    fig.tight_layout()
    fig.savefig(FIG / "01_pca_runday.png", dpi=150)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(7.5, 5.6), sharex=True)
    for ax, pc in zip(axes, ("PC1", "PC2")):
        style(ax)
        for g in ("normal", "tumor", "pool"):
            d = samples[grp == g]
            ax.scatter(d.run_order, scores.loc[d.index, pc], s=34, marker=MARKERS[g], color=COLORS[g],
                       edgecolor=SURFACE, linewidth=1.2, label=g, zorder=3)
        ax.axvline(samples.loc[samples.run_day == "day1", "run_order"].max() + 0.5, color=INK2,
                   linewidth=1, linestyle=(0, (3, 3)))
        ax.set_ylabel(f"{pc} score")
    axes[0].set_title("PCA scores by inferred run order (dashed line: day 1 | day 2)", loc="left", color=INK)
    axes[1].set_xlabel("Inferred run order")
    legend(axes[0])
    fig.tight_layout()
    fig.savefig(FIG / "01_run_order.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    style(ax)
    # Long biological tail is folded into a final ">=150%" bin so the 0-100% range stays readable.
    cap = 150
    bins = np.arange(0, cap + 10, 5)
    ax.hist(np.clip(bio_cv, 0, cap + 2.5), bins=bins, color=COLORS["normal"], alpha=0.85,
            label=f"biological samples (n=78); {(bio_cv >= cap).sum()} metabolites ≥{cap}%",
            edgecolor=SURFACE, linewidth=1)
    ax.hist(np.clip(pool_cv, 0, cap + 2.5), bins=bins, color=COLORS["pool"], alpha=0.85,
            label="QC pools (n=9)", edgecolor=SURFACE, linewidth=1)
    top = ax.get_ylim()[1]
    for t, y in ((20, 0.95), (30, 0.86)):
        ax.axvline(t, color=INK2, linewidth=1, linestyle=(0, (3, 3)))
        ax.text(t + 1.5, top * y, f"{t}%", color=INK2, fontsize=8)
    ax.set_xlim(0, cap + 10)
    ax.set_xticks(list(range(0, cap, 25)) + [cap])
    ax.set_xticklabels([str(x) for x in range(0, cap, 25)] + [f"≥{cap}"])
    ax.set(xlabel="Coefficient of variation per metabolite (%)", ylabel="Metabolites")
    ax.set_title("Metabolite variability: QC pools vs biological samples", loc="left", color=INK)
    legend(ax)
    fig.tight_layout()
    fig.savefig(FIG / "01_cv_pools.png", dpi=150)
    plt.close(fig)
    print(f"\nwrote figures/01_pca_group.png, 01_pca_runday.png, 01_run_order.png, 01_cv_pools.png")


if __name__ == "__main__":
    REPORT.parent.mkdir(exist_ok=True)
    with open(REPORT, "w") as fh, contextlib.redirect_stdout(Tee(sys.stdout, fh)):
        main()
