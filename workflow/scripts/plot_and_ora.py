#!/usr/bin/env python3
import os; os.environ["MPLBACKEND"] = "Agg"
import argparse
import os
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
import matplotlib.pyplot as plt

def bh_adjust(pvals):
    """Benjamini-Hochberg FDR correction."""
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    if n == 0:
        return pvals
    order = np.argsort(pvals)
    ranked = pvals[order]
    adj = np.empty(n, dtype=float)
    min_val = 1.0
    for i in range(n - 1, -1, -1):
        val = ranked[i] * n / (i + 1)
        min_val = min(min_val, val)
        adj[i] = min_val
    out = np.empty(n, dtype=float)
    out[order] = np.clip(adj, 0.0, 1.0)
    return out

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--deseq", required=True)
    parser.add_argument("--counts", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--gmt", required=True)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--lfc", type=float, default=1.0)
    parser.add_argument("--out-volcano", required=True)
    parser.add_argument("--out-pca", required=True)
    parser.add_argument("--out-ora", required=True)
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.out_volcano), exist_ok=True)
    os.makedirs(os.path.dirname(args.out_ora), exist_ok=True)

    # 1. Load DESeq2 results & classify DEGs
    res = pd.read_csv(args.deseq, sep="\t")
    tested = res.dropna(subset=["padj", "log2FoldChange"]).copy()
    tested["neg_log10_padj"] = -np.log10(np.clip(tested["padj"], 1e-300, 1.0))

    tested["status"] = "Not Significant"
    up_mask = (tested["padj"] < args.fdr) & (tested["log2FoldChange"] >= args.lfc)
    down_mask = (tested["padj"] < args.fdr) & (tested["log2FoldChange"] <= -args.lfc)
    tested.loc[up_mask, "status"] = "Up-regulated"
    tested.loc[down_mask, "status"] = "Down-regulated"

    # 2. Volcano Plot
    fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
    colors = {"Not Significant": "#9ca3af", "Up-regulated": "#dc2626", "Down-regulated": "#2563eb"}
    for status, group in tested.groupby("status"):
        ax.scatter(
            group["log2FoldChange"], group["neg_log10_padj"],
            c=colors[status], label=f"{status} (n={len(group)})",
            alpha=0.75, s=22, edgecolors="none"
        )
    ax.axvline(args.lfc, color="#374151", linestyle="--", linewidth=0.9)
    ax.axvline(-args.lfc, color="#374151", linestyle="--", linewidth=0.9)
    ax.axhline(-np.log10(args.fdr), color="#374151", linestyle="--", linewidth=0.9)

    # Annotate top 8 genes by adjusted p-value
    top_genes = tested.nsmallest(8, "padj")
    for _, row in top_genes.iterrows():
        ax.annotate(
            row["gene_id"],
            (row["log2FoldChange"], row["neg_log10_padj"]),
            xytext=(4, 4), textcoords="offset points", fontsize=8, fontweight="bold"
        )

    ax.set_xlabel("log2 Fold Change (Treated vs. Control)")
    ax.set_ylabel("-log10 Adjusted P-value (FDR)")
    ax.set_title("Differential Expression Volcano Plot (DESeq2)")
    ax.legend(frameon=True)
    fig.tight_layout()
    fig.savefig(args.out_volcano)
    fig.savefig(args.out_volcano.replace('.pdf', '.png'), dpi=200)
    plt.close(fig)

    # 3. Centered PCA on log2(normalized_counts + 1)
    counts_df = pd.read_csv(args.counts, sep="\t").set_index("gene_id")
    samples_df = pd.read_csv(args.samples, sep="\t").set_index("sample_id")
    log_counts = np.log2(counts_df[samples_df.index] + 1.0)

    # Select top 500 most variable genes
    top_var_genes = log_counts.var(axis=1).nlargest(min(500, len(log_counts))).index
    X = log_counts.loc[top_var_genes].T.values
    X_centered = X - X.mean(axis=0, keepdims=True)
    U, S, Vt = np.linalg.svd(X_centered, full_matrices=False)
    scores = U * S
    var_exp = (S ** 2) / np.sum(S ** 2) * 100.0

    fig, ax = plt.subplots(figsize=(7, 5.5), dpi=150)
    cond_colors = {"control": "#2563eb", "treated": "#dc2626"}
    for i, sample_id in enumerate(samples_df.index):
        cond = samples_df.loc[sample_id, "condition"]
        ax.scatter(scores[i, 0], scores[i, 1], c=cond_colors.get(cond, "#000000"), s=90, label=cond)
        ax.annotate(sample_id, (scores[i, 0], scores[i, 1]), xytext=(6, 4), textcoords="offset points", fontsize=9)

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), title="Condition")
    ax.set_xlabel(f"PC1 ({var_exp[0]:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({var_exp[1]:.1f}% variance)")
    ax.set_title("Sample PCA on Top Variable Genes (log2 Normalized Counts)")
    ax.grid(True, linestyle=":", alpha=0.5)
    fig.tight_layout()
    fig.savefig(args.out_pca)
    fig.savefig(args.out_pca.replace('.pdf', '.png'), dpi=200)
    plt.close(fig)

    # 4. One-sided Hypergeometric Over-Representation Analysis (ORA)
    universe = set(tested["gene_id"].unique())
    sig_genes = set(tested.loc[tested["padj"] < args.fdr, "gene_id"].unique())
    # Fallback to nominal p < 0.05 if subsampled read depth yields <5 FDR hits
    if len(sig_genes) < 5:
        sig_genes = set(tested.loc[tested["pvalue"] < 0.05, "gene_id"].unique())

    M = len(universe)
    N = len(sig_genes)

    ora_rows = []
    with open(args.gmt) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            pathway, desc, members = parts[0], parts[1], set(parts[2:])
            in_univ = members & universe
            n = len(in_univ)
            if n == 0:
                continue
            overlap = sig_genes & in_univ
            k = len(overlap)
            # P(X >= k) = sf(k - 1, M, n, N)
            pval = float(hypergeom.sf(k - 1, M, n, N)) if k > 0 else 1.0
            ora_rows.append({
                "pathway": pathway,
                "description": desc,
                "universe_size_M": M,
                "pathway_in_universe_n": n,
                "query_size_N": N,
                "overlap_k": k,
                "p_value": pval,
                "overlap_genes": ",".join(sorted(overlap))
            })

    ora_df = pd.DataFrame(ora_rows)
    if not ora_df.empty:
        ora_df["padj_bh"] = bh_adjust(ora_df["p_value"].values)
        ora_df = ora_df.sort_values("p_value")
    ora_df.to_csv(args.out_ora, sep="\t", index=False)
    print(f"Downstream plots and ORA saved. Top pathway: {ora_df.iloc[0]['pathway'] if not ora_df.empty else 'None'}")

if __name__ == "__main__":
    main()
