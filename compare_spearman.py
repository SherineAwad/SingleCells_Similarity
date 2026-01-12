#!/usr/bin/env python3

import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import warnings
from scipy.stats import spearmanr

warnings.filterwarnings('ignore')

def calculate_spearman_correlation(X_ref, X_query, ref_labels, query_labels):
    ref_celltypes = np.unique(ref_labels)
    query_celltypes = np.unique(query_labels)
    correlation_matrix = np.zeros((len(ref_celltypes), len(query_celltypes)))
    
    # Calculate median expression profiles for each cell type
    ref_profiles = {}
    for ct in ref_celltypes:
        ref_mask = (ref_labels == ct)
        if np.sum(ref_mask) > 0:
            # Use median expression across cells (more robust for correlation)
            ref_profiles[ct] = np.median(X_ref[ref_mask], axis=0)
    
    query_profiles = {}
    for ct in query_celltypes:
        query_mask = (query_labels == ct)
        if np.sum(query_mask) > 0:
            query_profiles[ct] = np.median(X_query[query_mask], axis=0)
    
    # Calculate Spearman correlation between profiles
    for i, ref_ct in enumerate(ref_celltypes):
        if ref_ct not in ref_profiles:
            continue
            
        for j, query_ct in enumerate(query_celltypes):
            if query_ct not in query_profiles:
                continue
                
            # Calculate Spearman correlation using ALL genes
            corr, p_value = spearmanr(ref_profiles[ref_ct], query_profiles[query_ct])
            correlation_matrix[i, j] = corr
    
    return ref_celltypes, query_celltypes, correlation_matrix


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

    # Use ALL genes (no HVG selection)
    print("Using ALL genes for correlation analysis...")
    
    # Ensure both datasets have the same genes in same order
    # Find intersection of genes
    common_genes = adata_ref.var_names.intersection(adata_query.var_names)
    print(f"Number of common genes: {len(common_genes)}")
    
    # Subset to common genes
    adata_ref = adata_ref[:, common_genes].copy()
    adata_query = adata_query[:, common_genes].copy()

    print(f"After gene intersection - Ref: {adata_ref.shape}, Query: {adata_query.shape}")

    # Convert to dense np arrays
    X_ref = adata_ref.X.toarray() if hasattr(adata_ref.X, "toarray") else adata_ref.X.copy()
    X_query = adata_query.X.toarray() if hasattr(adata_query.X, "toarray") else adata_query.X.copy()

    # Convert to float32
    X_ref = X_ref.astype(np.float32)
    X_query = X_query.astype(np.float32)

    # Get labels
    ref_labels = adata_ref.obs[celltype_col].values
    query_labels = adata_query.obs[celltype_col].values

    print("Calculating Spearman correlation matrix using ALL genes...")
    ref_cts, query_cts, corr_matrix = calculate_spearman_correlation(X_ref, X_query, ref_labels, query_labels)

    corr_df = pd.DataFrame(corr_matrix, index=ref_cts, columns=query_cts)

    # Plot heatmap
    print("Plotting correlation matrix...")
    plt.figure(figsize=(max(12, len(corr_df.columns) * 0.8),
                        max(10, len(corr_df.index) * 0.7)))

    im = plt.imshow(corr_df.values, aspect='auto', cmap='coolwarm', vmin=-1, vmax=1)
    plt.colorbar(im, label="Spearman Correlation (ρ)", shrink=0.8)

    plt.xticks(range(len(corr_df.columns)), corr_df.columns, rotation=45, ha='right', fontsize=10)
    plt.yticks(range(len(corr_df.index)), corr_df.index, fontsize=10)

    plt.xlabel(f"Query: {args.sample2}", fontsize=12, fontweight='bold')
    plt.ylabel(f"Reference: {args.sample1}", fontsize=12, fontweight='bold')
    plt.title(f"Cell Type Similarity: {args.sample1} → {args.sample2}\n(Spearman Correlation - ALL Genes)", 
              fontsize=14, fontweight='bold', pad=20)

    # Add text inside cells
    for i in range(len(corr_df.index)):
        for j in range(len(corr_df.columns)):
            val = corr_df.iloc[i, j]
            if abs(val) > 0.01:  # Only show non-negligible correlations
                plt.text(j, i, f"{val:.3f}",
                         ha="center", va="center",
                         color="white" if abs(val) > 0.5 else "black",
                         fontsize=8 if len(corr_df.index) < 20 else 6,
                         fontweight='bold' if abs(val) > 0.7 else 'normal')

    plt.tight_layout()
    plt.savefig(args.out, dpi=300, bbox_inches='tight')
    plt.close()

    # Save correlation matrix CSV
    csv_path = args.out.rsplit('.', 1)[0] + '_spearman_ALL_genes_correlation.csv'
    corr_df.to_csv(csv_path)
    print(f"✓ Saved correlation matrix to: {csv_path}")

    # Save summary report
    report_path = args.out.rsplit('.', 1)[0] + '_spearman_ALL_genes_report.txt'
    with open(report_path, 'w') as f:
        f.write("Spearman Correlation Report (ALL Genes)\n")
        f.write("=========================================\n")
        f.write(f"Reference sample: {args.sample1} ({adata_ref.shape[0]} cells)\n")
        f.write(f"Query sample: {args.sample2} ({adata_query.shape[0]} cells)\n")
        f.write(f"Number of genes analyzed: {adata_ref.shape[1]}\n")
        f.write("\nCell Types:\n")
        f.write(f"  Reference: {list(ref_cts)}\n")
        f.write(f"  Query: {list(query_cts)}\n\n")
        
        # Calculate and report average correlation
        avg_corr = np.nanmean(corr_matrix)
        median_corr = np.nanmedian(corr_matrix)
        f.write(f"Overall correlation statistics:\n")
        f.write(f"  Average correlation: {avg_corr:.4f}\n")
        f.write(f"  Median correlation: {median_corr:.4f}\n")
        
        f.write("\nTop positive correlations (ρ > 0.7):\n")
        positive_count = 0
        for i in range(len(corr_df.index)):
            for j in range(len(corr_df.columns)):
                val = corr_df.iloc[i, j]
                if val > 0.7:
                    f.write(f"  {corr_df.index[i]} -> {corr_df.columns[j]}: {val:.4f}\n")
                    positive_count += 1
        
        f.write(f"\nNumber of strong positive correlations (ρ > 0.7): {positive_count}\n")
        
        f.write("\nTop negative correlations (ρ < -0.5):\n")
        negative_count = 0
        for i in range(len(corr_df.index)):
            for j in range(len(corr_df.columns)):
                val = corr_df.iloc[i, j]
                if val < -0.5:
                    f.write(f"  {corr_df.index[i]} -> {corr_df.columns[j]}: {val:.4f}\n")
                    negative_count += 1
        
        f.write(f"\nNumber of strong negative correlations (ρ < -0.5): {negative_count}\n")

    print(f"✓ Saved report to: {report_path}")
    print(f"\n✅ DONE: {args.out}")


if __name__ == "__main__":
    main()
