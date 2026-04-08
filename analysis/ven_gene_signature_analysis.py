"""
VEN Gene Signature Analysis - Frontal Cortex Spaceflight Data
NASA OSD-698 / GEO GSE239336

Question: Are genes associated with Von Economo Neuron (VEN) function and
fast-projection signalling specifically dysregulated in mouse frontal cortex
following real ISS spaceflight?

Data: GeoMx Digital Spatial Profiling of frontal cortex tissue from mice
flown on SpaceX-24 (35 days ISS). Ground Control vs Spaceflight (saline).

Usage:
    python analysis/ven_gene_signature_analysis.py

Outputs:
    figures/fig_ven_volcano.pdf - volcano plot with VEN genes highlighted
    figures/fig_ven_heatmap.pdf - VEN gene expression heatmap
    figures/fig_ven_significant.pdf - bar chart of significant VEN genes
    figures/fig_ven_categories.pdf - category-level summary with statistics
    results/ven_gene_results.json - full results table with category analysis
"""

import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from scipy import stats
from pathlib import Path

ROOT = Path(__file__).parent.parent
DATA = ROOT / "data"
FIGS = ROOT / "figures"
RES  = ROOT / "results"
FIGS.mkdir(exist_ok=True)
RES.mkdir(exist_ok=True)

#  VEN gene signature 
# Curated from: Allman et al. 2010, Stimpson et al. 2011, Hodge et al. 2020 (Allen Brain Atlas), von Economo 1925
VEN_GENES = {
    # Core VEN identity markers
    "VEN Identity": {
        "Fxyd1": "VEN marker (ion transport)",
        "Gadd45g": "VEN marker (stress response)",
        "Adra1a": "VEN marker (adrenergic receptor)",
        "Nos1": "nNOS  VEN biochemical marker",
        "Disc1": "VEN-associated (schizophrenia risk)",
        "Bcl11b": "Layer V projection neuron identity",
    },
    # Fast signalling - large neuron markers
    "Fast Signalling": {
        "Nefh": "Heavy neurofilament - large axon calibre",
        "Nefm": "Medium neurofilament - axon structure",
        "Nefl": "Light neurofilament - axon structure",
        "Syt1": "Synaptotagmin - fast synaptic release",
        "Snap25": "SNARE protein - fast vesicle fusion",
        "Vamp2":  "Vesicle-associated - fast transmission",
    },
    # Myelination - fast conduction speed
    "Myelination": {
        "Mbp": "Myelin basic protein",
        "Mog": "Myelin oligodendrocyte glycoprotein",
        "Plp1": "Proteolipid protein - myelin sheath",
        "Mag": "Myelin-associated glycoprotein",
        "Cnp": "2,3-cyclic nucleotide phosphodiesterase",
    },
    # Oxidative stress - VENs are large (high metabolic demand)
    "Oxidative Stress": {
        "Sod1": "Superoxide dismutase 1 - antioxidant",
        "Sod2": "Superoxide dismutase 2 - mitochondrial",
        "Cat":  "Catalase - H2O2 clearance",
        "Gpx1": "Glutathione peroxidase",
        "Nox4": "NADPH oxidase - ROS generation",
    },
    # Synaptic integration
    "Synaptic": {
        "Grin1":  "NMDA receptor subunit 1",
        "Grin2a": "NMDA receptor subunit 2A",
        "Grin2b": "NMDA receptor subunit 2B",
        "Dlg4":   "PSD-95 - postsynaptic scaffold",
        "Shank3": "Postsynaptic density protein",
    },
}

ALL_VEN_GENES = {g: cat for cat, genes in VEN_GENES.items() for g in genes}
GENE_DESC = {g: d for cat, genes in VEN_GENES.items() for g, d in genes.items()}

CAT_COLORS = {
    "VEN Identity":    "#E05252",
    "Fast Signalling": "#3A7FC1",
    "Myelination":     "#52C869",
    "Oxidative Stress":"#E8A020",
    "Synaptic":        "#9B59B6",
}


def load_de():
    path = DATA / "GSE239336_FCT_GCvsFLT-SAL_DEanalysis.txt"
    df = pd.read_csv(path, sep="\t", skiprows=6)
    df.columns = ["tag", "group", "gene", "log2fc",
                  "pvalue", "adj_pvalue", "neg_log10_p", "neg_log10_adjp"]
    df["gene"] = df["gene"].astype(str).str.strip()
    df["log2fc"] = pd.to_numeric(df["log2fc"], errors="coerce")
    df["pvalue"] = pd.to_numeric(df["pvalue"], errors="coerce")
    df["adj_pvalue"] = pd.to_numeric(df["adj_pvalue"], errors="coerce")
    df = df.dropna(subset=["gene", "log2fc", "pvalue"])
    df = df[df["gene"] != "nan"].reset_index(drop=True)
    return df


def load_expr():
    path = DATA / "GSE239336_FCT_GeneExpression_Q3norm.txt"
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str).str.strip()
    return df


def category_level_analysis(ven_de):
    """
    Category-level statistical analysis using one-sample t-tests.
    Tests whether each functional category shows coordinated expression change.
    More powerful than individual gene tests with small sample sizes.
    """
    print("CATEGORY-LEVEL ANALYSIS (Primary Statistical Test)")
    print(f"{'Category':<20} {'n':>3} {'Mean log2FC':>12} {'SEM':>8} {'t':>8} {'p-value':>10} {'Sig':>5}")
    
    results = []
    for cat in VEN_GENES.keys():
        sub = ven_de[ven_de.category == cat]
        log2fcs = sub.log2fc.values
        
        # One-sample t-test: is mean log2FC significantly different from 0?
        t_stat, p_val = stats.ttest_1samp(log2fcs, 0)
        sem = stats.sem(log2fcs)
        
        sig_marker = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else "ns"
        
        print(f"{cat:<20} {len(log2fcs):>3} {log2fcs.mean():>12.3f} {sem:>8.3f} "
              f"{t_stat:>8.3f} {p_val:>10.4f} {sig_marker:>5}")
        
        results.append({
            "category": cat,
            "n_genes": int(len(log2fcs)),
            "mean_log2fc": float(log2fcs.mean()),
            "sem": float(sem),
            "t_statistic": float(t_stat),
            "p_value": float(p_val),
            "df": int(len(log2fcs) - 1),
            "significant": bool(p_val < 0.05),
            "direction": "up" if log2fcs.mean() > 0 else "down"
        })
    
    print("Interpretation: One-sample t-tests vs. null hypothesis of no change (log2FC=0)")
    print("* p<0.05, ** p<0.01, *** p<0.001")
    
    return results


def run_analysis(de, expr):
    print("VEN GENE SIGNATURE - FRONTAL CORTEX SPACEFLIGHT (OSD-698)")
    print(f"\nTotal genes in dataset: {len(de)}")
    print(f"Significant overall (p<0.05): {(de.pvalue < 0.05).sum()}")
    print(f"Up in spaceflight (log2fc>0, p<0.05): "
          f"{((de.log2fc > 0) & (de.pvalue < 0.05)).sum()}")
    print(f"Down in spaceflight (log2fc<0, p<0.05): "
          f"{((de.log2fc < 0) & (de.pvalue < 0.05)).sum()}")

    # Extract VEN genes
    ven_de = de[de["gene"].isin(ALL_VEN_GENES)].copy()
    ven_de["category"] = ven_de["gene"].map(ALL_VEN_GENES)
    ven_de["description"] = ven_de["gene"].map(GENE_DESC)
    ven_de = ven_de.sort_values("pvalue")

    print(f"\nVEN genes found in dataset: {len(ven_de)}/{len(ALL_VEN_GENES)}")
    print(f"VEN genes significant (p<0.05): {(ven_de.pvalue < 0.05).sum()}")

    # CATEGORY-LEVEL ANALYSIS (PRIMARY)
    cat_results = category_level_analysis(ven_de)

    # INDIVIDUAL GENE ANALYSIS (SECONDARY/EXPLORATORY)
    print("INDIVIDUAL GENE RESULTS (Secondary/Exploratory)")
    print(f"{'Gene':<12} {'Category':<18} {'Log2FC':>8} {'p-value':>10} "
          f"{'adj.p':>8} {'Direction':<12}")
    for _, row in ven_de.iterrows():
        direction = "↑ spaceflight" if row.log2fc > 0 else "↓ spaceflight"
        sig = "***" if row.pvalue < 0.001 else "**" if row.pvalue < 0.01 \
              else "*" if row.pvalue < 0.05 else ""
        print(f"{row.gene:<12} {row.category:<18} {row.log2fc:>8.3f} "
              f"{row.pvalue:>10.4f} {str(row.adj_pvalue):>8} "
              f"{direction:<12} {sig}")

    # Enrichment test: are VEN genes more often significant than background?
    n_total = len(de)
    n_total_sig = (de.pvalue < 0.05).sum()
    n_ven = len(ven_de)
    n_ven_sig = (ven_de.pvalue < 0.05).sum()
    expected = n_total_sig / n_total * n_ven
    
    print("ENRICHMENT ANALYSIS")
    print(f"Background sig rate: {n_total_sig}/{n_total} = "
          f"{n_total_sig/n_total*100:.1f}%")
    print(f"VEN gene sig rate:   {n_ven_sig}/{n_ven} = "
          f"{n_ven_sig/n_ven*100:.1f}%")
    print(f"Expected by chance:  {expected:.1f}")
    print(f"Observed:            {n_ven_sig}")

    # Fisher exact test
    contingency = [[n_ven_sig, n_ven - n_ven_sig],
                   [n_total_sig - n_ven_sig,
                    n_total - n_total_sig - (n_ven - n_ven_sig)]]
    odds, p_fisher = stats.fisher_exact(contingency)
    print(f"Fisher exact: OR={odds:.2f}, p={p_fisher:.4f}")
    print("="*80)

    return ven_de, {
        "category_analysis": cat_results,
        "enrichment_p": float(p_fisher), 
        "odds_ratio": float(odds),
        "n_ven_sig": int(n_ven_sig), 
        "n_ven": int(n_ven),
        "background_rate": float(n_total_sig/n_total)
    }


def plot_volcano(de, ven_de):
    fig, ax = plt.subplots(figsize=(10, 7))

    non_sig = de[de.pvalue >= 0.05]
    sig_bg  = de[(de.pvalue < 0.05) & ~de.gene.isin(ALL_VEN_GENES)]
    ax.scatter(non_sig.log2fc, -np.log10(non_sig.pvalue + 1e-10),
               s=4, c="#cccccc", alpha=0.4, linewidths=0, label="All genes (ns)")
    ax.scatter(sig_bg.log2fc, -np.log10(sig_bg.pvalue + 1e-10),
               s=6, c="#aaaaaa", alpha=0.5, linewidths=0, label="Sig. (non-VEN)")

    for cat, genes in VEN_GENES.items():
        sub = ven_de[ven_de.category == cat]
        ax.scatter(sub.log2fc, -np.log10(sub.pvalue + 1e-10),
                   s=80, c=CAT_COLORS[cat], alpha=0.9, linewidths=0.5,
                   edgecolors="white", label=cat, zorder=5)

    # Label significant VEN genes
    sig_ven = ven_de[ven_de.pvalue < 0.1]
    for _, row in sig_ven.iterrows():
        ax.annotate(
            row.gene,
            xy=(row.log2fc, -np.log10(row.pvalue + 1e-10)),
            xytext=(8, 4), textcoords="offset points",
            fontsize=9, fontweight="bold",
            color=CAT_COLORS[row.category],
            path_effects=[pe.withStroke(linewidth=2, foreground="white")]
        )

    # Reference lines
    ax.axhline(-np.log10(0.05), color="gray", linestyle="--",
               alpha=0.6, linewidth=1, label="p=0.05")
    ax.axvline(0, color="gray", linestyle="-", alpha=0.3, linewidth=0.8)

    ax.set_xlabel("Log₂ Fold Change (Spaceflight vs Ground Control)", fontsize=12)
    ax.set_ylabel("-log₁₀(p-value)", fontsize=12)
    ax.set_title(
        "VEN Gene Signature in Frontal Cortex - Real ISS Spaceflight\n"
        "NASA OSD-698 (SpaceX-24, 35 days ISS, GeoMx DSP)",
        fontweight="bold", fontsize=12)
    ax.legend(fontsize=9, framealpha=0.9, loc="upper left",
              markerscale=1.5)
    ax.grid(alpha=0.2)
    ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    out = FIGS / "fig_ven_volcano.pdf"
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"\nSaved: {out}")


def plot_significant_bars(ven_de):
    # Show all VEN genes sorted by log2fc, coloured by significance
    df = ven_de.sort_values("log2fc")
    colors = []
    for _, row in df.iterrows():
        if row.pvalue < 0.05:
            colors.append("#E05252" if row.log2fc > 0 else "#3A7FC1")
        else:
            colors.append("#cccccc")

    fig, ax = plt.subplots(figsize=(10, 8))
    bars = ax.barh(df.gene, df.log2fc, color=colors, edgecolor="white",
                   linewidth=0.5, height=0.7)

    # Add p-value labels
    for i, (_, row) in enumerate(df.iterrows()):
        sig = "***" if row.pvalue < 0.001 else "**" if row.pvalue < 0.01 \
              else "*" if row.pvalue < 0.05 else f"p={row.pvalue:.2f}"
        x = row.log2fc + (0.02 if row.log2fc >= 0 else -0.02)
        ha = "left" if row.log2fc >= 0 else "right"
        ax.text(x, i, sig, va="center", ha=ha, fontsize=8,
                color="#333333" if row.pvalue >= 0.05 else "black",
                fontweight="bold" if row.pvalue < 0.05 else "normal")

    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Log₂ Fold Change (Spaceflight vs Ground Control)", fontsize=11)
    ax.set_title(
        "VEN-Associated Gene Expression Changes in Frontal Cortex\n"
        "Real ISS Spaceflight - NASA OSD-698 | "
        "Red=up, Blue=down (p<0.05) | Grey=ns",
        fontweight="bold", fontsize=11)
    ax.grid(axis="x", alpha=0.3)
    ax.spines[["top", "right"]].set_visible(False)

    # Category labels on y-axis
    ytick_colors = [CAT_COLORS.get(ALL_VEN_GENES.get(g, ""), "#333")
                    for g in df.gene]
    for tick, color in zip(ax.get_yticklabels(), ytick_colors):
        tick.set_color(color)
        tick.set_fontweight("bold")

    plt.tight_layout()
    out = FIGS / "fig_ven_significant.pdf"
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out}")


def plot_category_summary(ven_de, cat_results):
    """Mean log2fc per category with individual gene points and significance markers."""
    fig, ax = plt.subplots(figsize=(9, 5))

    cats = list(VEN_GENES.keys())
    x = np.arange(len(cats))
    width = 0.5

    # Create lookup for category p-values
    cat_pvals = {r["category"]: r["p_value"] for r in cat_results}

    for i, cat in enumerate(cats):
        sub = ven_de[ven_de.category == cat]
        means = sub.log2fc.mean()
        color = CAT_COLORS[cat]

        ax.bar(i, means, width=width, color=color, alpha=0.7,
               edgecolor="white")
        ax.scatter([i + np.random.uniform(-0.15, 0.15)
                    for _ in range(len(sub))],
                   sub.log2fc.values, s=60, c=color,
                   edgecolors="white", linewidths=0.5, zorder=5)

        # Mark significant genes
        sig = sub[sub.pvalue < 0.05]
        for _, row in sig.iterrows():
            ax.annotate(row.gene,
                        xy=(i, row.log2fc),
                        xytext=(5, 0), textcoords="offset points",
                        fontsize=8, fontweight="bold", color=color)
        
        # Add category-level significance marker
        p_val = cat_pvals.get(cat, 1.0)
        if p_val < 0.05:
            sig_marker = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*"
            # Position marker above the bar
            y_pos = means + 0.05 if means > 0 else means - 0.05
            ax.text(i, y_pos, sig_marker, ha="center", va="bottom" if means > 0 else "top",
                    fontsize=16, fontweight="bold", color=color)

    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(cats, rotation=20, ha="right", fontsize=10)
    ax.set_ylabel("Log₂ Fold Change (Spaceflight vs Ground)", fontsize=11)
    ax.set_title(
        "VEN Gene Categories - Mean Expression Change in Spaceflight\n"
        "NASA OSD-698 Frontal Cortex | Gene labels = individual p<0.05 | */*** = category p<0.05/0.001",
        fontweight="bold", fontsize=11)
    ax.grid(axis="y", alpha=0.3)
    ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    out = FIGS / "fig_ven_categories.pdf"
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out}")


def main():
    print("Loading data...")
    de   = load_de()
    expr = load_expr()

    ven_de, stats_res = run_analysis(de, expr)

    print("\nGenerating figures...")
    plot_volcano(de, ven_de)
    plot_significant_bars(ven_de)
    plot_category_summary(ven_de, stats_res["category_analysis"])

    # Save results JSON
    records = []
    for _, row in ven_de.iterrows():
        records.append({
            "gene": row.gene,
            "category":    row.category,
            "description": row.description,
            "log2fc": round(float(row.log2fc), 4),
            "pvalue": round(float(row.pvalue), 6),
            "adj_pvalue":  str(row.adj_pvalue),
            "significant": bool(row.pvalue < 0.05),
            "direction":   "up" if row.log2fc > 0 else "down",
        })

    out_json = RES / "ven_gene_results.json"
    with open(out_json, "w") as f:
        json.dump({
            "category_analysis": stats_res["category_analysis"],
            "enrichment": {
                "enrichment_p": stats_res["enrichment_p"],
                "odds_ratio": stats_res["odds_ratio"],
                "n_ven_sig": stats_res["n_ven_sig"],
                "n_ven": stats_res["n_ven"],
                "background_rate": stats_res["background_rate"]
            },
            "genes": records
        }, f, indent=2)
    print(f"Saved: {out_json}")
    print("ANALYSIS COMPLETE")
    print("\nKEY FINDINGS:")
    for cat_res in stats_res["category_analysis"]:
        if cat_res["significant"]:
            print(f"  ✓ {cat_res['category']}: mean log2FC = {cat_res['mean_log2fc']:.3f}, "
                  f"t({cat_res['df']}) = {cat_res['t_statistic']:.2f}, p = {cat_res['p_value']:.4f}")
    print("="*80)


if __name__ == "__main__":
    main()
