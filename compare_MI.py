#!/usr/bin/env python3

import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import warnings
from sklearn.metrics import mutual_info_score

warnings.filterwarnings('ignore')

def calculate_mi(X_ref, X_query, ref_labels, query_labels, n_bins=20):
    ref_celltypes = np.unique(ref_labels)
    query_celltypes = np.unique(query_labels)
    similarity_matrix = np.zeros((len(ref_celltypes), len(query_celltypes)))

    # Combine data once to define consistent bins per gene
    combined = np.vstack([X_ref, X_query])
    gene_bins = []
    for gene_idx in range(combined.shape[1]):
        gene_min = combined[:, gene_idx].min()
        gene_max = combined[:, gene_idx].max()
        bins = np.linspace(gene_min, gene_max, n_bins + 1)
        gene_bins.append(bins)

    for i, ref_ct in enumerate(ref_celltypes):
        ref_mask = (ref_labels == ref_ct)
        if np.sum(ref_mask) == 0:
            continue
        ref_data = X_ref[ref_mask]

        for j, query_ct in enumerate(query_celltypes):
            query_mask = (query_labels == query_ct)
            if np.sum(query_mask) == 0:
                continue
            query_data = X_query[query_mask]

            mi_scores = []
            for gene_idx in range(ref_data.shape[1]):
                bins = gene_bins[gene_idx]

                # Histogram counts per gene per cell type
                ref_hist, _ = np.histogram(ref_data[:, gene_idx], bins=bins)
                query_hist, _ = np.histogram(query_data[:, gene_idx], bins=bins)

                # Add pseudocount to avoid zero count issues
                ref_hist = ref_hist + 1
                query_hist = query_hist + 1

                mi = mutual_info_score(ref_hist, query_hist)
                mi_scores.append(mi)

            similarity_matrix[i, j] = np.mean(mi_scores)

    return ref_celltypes, query_celltypes, similarity_matrix


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input h5ad file")
    parser.add_argument("--sample1", required=True, help="Reference sample name")
    parser.add_argument("--sample2", required=True, help="Query sample name")
    parser.add_argument("--out", required=True, help="Output PNG file for heatmap")
    args = parser.parse_args()

    print("Loading data...")
    adata = sc.read_h5ad(args.input)

    celltype_col = "celltype"
    sample_col = "renamed_samples"

    # Filter samples
    adata_ref = adata[adata.obs[sample_col] == args.sample1].copy()
    adata_query = adata[adata.obs[sample_col] == args.sample2].copy()

    if adata_ref.n_obs == 0 or adata_query.n_obs == 0:
        print("ERROR: empty reference or query sample")
        sys.exit(1)

    print(f"Reference ({args.sample1}): {adata_ref.shape}")
    print(f"Query ({args.sample2}): {adata_query.shape}")

    # Merge for HVG selection
    print("Merging datasets for HVG selection...")
    adata_combined = adata_ref.concatenate(adata_query, batch_key="batch")

    # Select highly variable genes - use cell_ranger flavor to avoid loess errors
    print("Selecting highly variable genes...")
    sc.pp.highly_variable_genes(adata_combined, flavor='cell_ranger', n_top_genes=2000, subset=True, inplace=True)

    # Subset original adatas by HVGs (matching var_names from combined)
    var_names = adata_combined.var_names
    adata_ref = adata_ref[:, var_names].copy()
    adata_query = adata_query[:, var_names].copy()

    print(f"After HVG selection - Ref: {adata_ref.shape}, Query: {adata_query.shape}")

    # Convert to dense np arrays
    X_ref = adata_ref.X.toarray() if hasattr(adata_ref.X, "toarray") else adata_ref.X.copy()
    X_query = adata_query.X.toarray() if hasattr(adata_query.X, "toarray") else adata_query.X.copy()

    # Convert to float32 to save memory
    X_ref = X_ref.astype(np.float32)
    X_query = X_query.astype(np.float32)

    # Get labels
    ref_labels = adata_ref.obs[celltype_col].values
    query_labels = adata_query.obs[celltype_col].values

    print("Calculating mutual information similarity matrix...")
    ref_cts, query_cts, sim_matrix = calculate_mi(X_ref, X_query, ref_labels, query_labels)

    sim_df = pd.DataFrame(sim_matrix, index=ref_cts, columns=query_cts)

    # Plot heatmap
    print("Plotting similarity matrix...")
    plt.figure(figsize=(max(12, len(sim_df.columns) * 0.8),
                        max(10, len(sim_df.index) * 0.7)))

    im = plt.imshow(sim_df.values, aspect='auto', cmap='viridis')
    plt.colorbar(im, label="Mutual Information Score", shrink=0.8)

    plt.xticks(range(len(sim_df.columns)), sim_df.columns, rotation=45, ha='right', fontsize=10)
    plt.yticks(range(len(sim_df.index)), sim_df.index, fontsize=10)

    plt.xlabel(f"Query: {args.sample2}", fontsize=12, fontweight='bold')
    plt.ylabel(f"Reference: {args.sample1}", fontsize=12, fontweight='bold')
    plt.title(f"Cell Type Similarity: {args.sample1} → {args.sample2}\n(Mutual Information)", fontsize=14, fontweight='bold', pad=20)

    # Add text inside cells
    for i in range(len(sim_df.index)):
        for j in range(len(sim_df.columns)):
            val = sim_df.iloc[i, j]
            if val > 0.01:
                plt.text(j, i, f"{val:.3f}",
                         ha="center", va="center",
                         color="white" if val > 0.5 else "black",
                         fontsize=8 if len(sim_df.index) < 20 else 6,
                         fontweight='bold' if val > 0.7 else 'normal')

    plt.tight_layout()
    plt.savefig(args.out, dpi=300, bbox_inches='tight')
    plt.close()

    # Save similarity matrix CSV
    csv_path = args.out.rsplit('.', 1)[0] + '_MI_similarity.csv'
    sim_df.to_csv(csv_path)
    print(f"✓ Saved similarity matrix to: {csv_path}")

    # Save summary report
    report_path = args.out.rsplit('.', 1)[0] + '_MI_report.txt'
    with open(report_path, 'w') as f:
        f.write("Mutual Information Similarity Report\n")
        f.write("===============================\n")
        f.write(f"Reference sample: {args.sample1} ({adata_ref.shape[0]} cells)\n")
        f.write(f"Query sample: {args.sample2} ({adata_query.shape[0]} cells)\n")
        f.write("\nCell Types:\n")
        f.write(f"  Reference: {list(ref_cts)}\n")
        f.write(f"  Query: {list(query_cts)}\n\n")
        f.write("Top similarities (MI > 0.5):\n")
        for i in range(len(sim_df.index)):
            for j in range(len(sim_df.columns)):
                val = sim_df.iloc[i, j]
                if val > 0.5:
                    f.write(f"  {sim_df.index[i]} -> {sim_df.columns[j]}: {val:.4f}\n")

    print(f"✓ Saved report to: {report_path}")
    print(f"\n✅ DONE: {args.out}")


if __name__ == "__main__":
    main()

