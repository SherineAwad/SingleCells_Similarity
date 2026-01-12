#!/usr/bin/env python3

import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import warnings
from sklearn.metrics.pairwise import cosine_similarity
from scipy import sparse

warnings.filterwarnings('ignore')

def calculate_cosine_similarity_sparse(X_ref, X_query, ref_labels, query_labels):
    ref_celltypes = np.unique(ref_labels)
    query_celltypes = np.unique(query_labels)
    similarity_matrix = np.zeros((len(ref_celltypes), len(query_celltypes)))

    # Calculate mean expression profiles for each cell type using sparse operations
    ref_profiles = {}
    for ct in ref_celltypes:
        ref_mask = (ref_labels == ct)
        if np.sum(ref_mask) > 0:
            # Use sparse mean calculation
            if sparse.issparse(X_ref):
                ref_profiles[ct] = np.array(X_ref[ref_mask].mean(axis=0)).flatten()
            else:
                ref_profiles[ct] = X_ref[ref_mask].mean(axis=0)

    query_profiles = {}
    for ct in query_celltypes:
        query_mask = (query_labels == ct)
        if np.sum(query_mask) > 0:
            if sparse.issparse(X_query):
                query_profiles[ct] = np.array(X_query[query_mask].mean(axis=0)).flatten()
            else:
                query_profiles[ct] = X_query[query_mask].mean(axis=0)

    # Calculate cosine similarity between profiles
    for i, ref_ct in enumerate(ref_celltypes):
        if ref_ct not in ref_profiles:
            continue

        for j, query_ct in enumerate(query_celltypes):
            if query_ct not in query_profiles:
                continue

            # Reshape for cosine similarity
            ref_vec = ref_profiles[ref_ct].reshape(1, -1)
            query_vec = query_profiles[query_ct].reshape(1, -1)

            # Calculate cosine similarity
            similarity = cosine_similarity(ref_vec, query_vec)[0, 0]
            similarity_matrix[i, j] = similarity

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

    # USE ALL GENES - Remove HVG selection
    print("Using ALL genes for cosine similarity analysis...")

    # Ensure both datasets have the same genes in same order
    common_genes = adata_ref.var_names.intersection(adata_query.var_names)
    print(f"Number of common genes: {len(common_genes)}")

    # Subset to common genes
    adata_ref = adata_ref[:, common_genes].copy()
    adata_query = adata_query[:, common_genes].copy()

    print(f"After gene intersection - Ref: {adata_ref.shape}, Query: {adata_query.shape}")

    # KEEP AS SPARSE - DO NOT CONVERT TO DENSE
    print("Keeping data in sparse format for memory efficiency...")
    X_ref = adata_ref.X
    X_query = adata_query.X

    # Only convert small profiles to dense, not entire matrices
    if sparse.issparse(X_ref):
        print(f"Reference data is sparse: {X_ref.shape}, density: {X_ref.nnz / (X_ref.shape[0] * X_ref.shape[1]):.4f}")
    if sparse.issparse(X_query):
        print(f"Query data is sparse: {X_query.shape}, density: {X_query.nnz / (X_query.shape[0] * X_ref.shape[1]):.4f}")

    # Get labels
    ref_labels = adata_ref.obs[celltype_col].values
    query_labels = adata_query.obs[celltype_col].values

    print("Calculating cosine similarity matrix using sparse operations...")
    ref_cts, query_cts, sim_matrix = calculate_cosine_similarity_sparse(X_ref, X_query, ref_labels, query_labels)

    # Create DataFrame from similarity matrix
    sim_df = pd.DataFrame(sim_matrix, index=ref_cts, columns=query_cts)
    
    # REORDER CELL TYPES to match desired order
    desired_order = ['MG', 'MGPC', 'PR precursors', 'Rod', 'Cones', 'BC', 'AC', 'HC', 'RGC', 'Microglia_Immunecells', 'RPE', 'Melanocyte', 'Endothelial', 'Pericytes', 'Oligodendrocyte']
    
    # Filter to include only cell types present in the data but maintain desired order
    ref_order = [ct for ct in desired_order if ct in ref_cts]
    query_order = [ct for ct in desired_order if ct in query_cts]
    
    # Reorder the similarity matrix
    sim_df = sim_df.loc[ref_order, query_order]

    # Plot heatmap
    print("Plotting similarity matrix...")
    plt.figure(figsize=(max(12, len(sim_df.columns) * 0.8),
                        max(10, len(sim_df.index) * 0.7)))

    im = plt.imshow(sim_df.values, aspect='auto', cmap='viridis')
    plt.colorbar(im, label="Cosine Similarity Score", shrink=0.8)

    plt.xticks(range(len(sim_df.columns)), sim_df.columns, rotation=45, ha='right', fontsize=10)
    plt.yticks(range(len(sim_df.index)), sim_df.index, fontsize=10)

    plt.xlabel(f"Query: {args.sample2}", fontsize=12, fontweight='bold')
    plt.ylabel(f"Reference: {args.sample1}", fontsize=12, fontweight='bold')
    plt.title(f"Cell Type Similarity: {args.sample1} → {args.sample2}\n(Cosine Similarity - ALL Genes)", fontsize=14, fontweight='bold', pad=20)

    # Display ALL values - ALL BLACK FONT
    for i in range(len(sim_df.index)):
        for j in range(len(sim_df.columns)):
            val = sim_df.iloc[i, j]
            plt.text(j, i, f"{val:.3f}",
                     ha="center", va="center",
                     color="black",  # CHANGED: Always use black font
                     fontsize=8 if len(sim_df.index) < 20 else 6,
                     fontweight='bold' if val > 0.7 else 'normal')

    plt.tight_layout()
    plt.savefig(args.out, dpi=300, bbox_inches='tight')
    plt.close()

    # Save similarity matrix CSV
    csv_path = args.out.rsplit('.', 1)[0] + '_cosine_ALL_genes_similarity.csv'
    sim_df.to_csv(csv_path)
    print(f"✓ Saved similarity matrix to: {csv_path}")

    # Save summary report
    report_path = args.out.rsplit('.', 1)[0] + '_cosine_ALL_genes_report.txt'
    with open(report_path, 'w') as f:
        f.write("Cosine Similarity Report (ALL Genes)\n")
        f.write("====================================\n")
        f.write(f"Reference sample: {args.sample1} ({adata_ref.shape[0]} cells)\n")
        f.write(f"Query sample: {args.sample2} ({adata_query.shape[0]} cells)\n")
        f.write(f"Number of genes analyzed: {adata_ref.shape[1]}\n")
        f.write("\nCell Types:\n")
        f.write(f"  Reference: {list(sim_df.index)}\n")  # Updated to show reordered list
        f.write(f"  Query: {list(sim_df.columns)}\n\n")  # Updated to show reordered list
        f.write("Top similarities (cosine > 0.7):\n")
        count_high = 0
        for i in range(len(sim_df.index)):
            for j in range(len(sim_df.columns)):
                val = sim_df.iloc[i, j]
                if val > 0.7:
                    f.write(f"  {sim_df.index[i]} -> {sim_df.columns[j]}: {val:.4f}\n")
                    count_high += 1

        f.write(f"\nStatistics:\n")
        f.write(f"  Number of strong similarities (>0.7): {count_high}\n")
        f.write(f"  Average similarity: {np.nanmean(sim_df.values):.4f}\n")
        f.write(f"  Median similarity: {np.nanmedian(sim_df.values):.4f}\n")
        f.write(f"  Minimum similarity: {np.nanmin(sim_df.values):.4f}\n")
        f.write(f"  Maximum similarity: {np.nanmax(sim_df.values):.4f}\n")

    print(f"✓ Saved report to: {report_path}")
    print(f"\n✅ DONE: {args.out}")


if __name__ == "__main__":
    main()
