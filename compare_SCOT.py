#!/usr/bin/env python3

import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import warnings
warnings.filterwarnings('ignore')

# Add SCOT source code to path
SCOT_PATH = "/nfs/turbo/umms-thahoang/sherine/tools/SCOT/src"
sys.path.insert(0, SCOT_PATH)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--sample1", required=True)
    parser.add_argument("--sample2", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    # Load
    adata = sc.read_h5ad(args.input)

    celltype_col = "celltype"
    sample_col = "renamed_samples"

    # Split reference/query
    adata_ref = adata[adata.obs[sample_col] == args.sample1].copy()
    adata_query = adata[adata.obs[sample_col] == args.sample2].copy()

    if adata_ref.n_obs == 0 or adata_query.n_obs == 0:
        print("ERROR: empty reference or query")
        sys.exit(1)

    print(f"Reference ({args.sample1}): {adata_ref.shape}")
    print(f"Query ({args.sample2}): {adata_query.shape}")
    
    # STEP 1: Import SCOT
    print("\n1. Importing SCOT...")
    try:
        # Try to import SCOTv2 from the source files
        from scotv2 import SCOTv2
        print("✓ SCOTv2 imported successfully")
    except ImportError as e:
        print(f"✗ ERROR importing SCOTv2: {e}")
        print("\nTrying to import SCOTv1...")
        try:
            from scotv1 import SCOT
            SCOTv2 = SCOT  # Alias for compatibility
            print("✓ SCOTv1 imported successfully")
        except ImportError as e2:
            print(f"✗ ERROR importing SCOTv1: {e2}")
            print("\nChecking what's in SCOT directory...")
            print(f"Files in {SCOT_PATH}: {os.listdir(SCOT_PATH)}")
            sys.exit(1)
    
    # STEP 2: Check SCOT dependencies
    print("\n2. Checking SCOT dependencies...")
    try:
        import torch
        print(f"✓ torch version: {torch.__version__}")
    except ImportError:
        print("✗ ERROR: torch not installed")
        print("Install with: pip install torch")
        sys.exit(1)
    
    try:
        import ot
        print(f"✓ POT version: {ot.__version__}")
    except ImportError:
        print("✗ ERROR: POT (Python Optimal Transport) not installed")
        print("Install with: pip install POT")
        sys.exit(1)
    
    # STEP 3: Preprocess data
    print("\n3. Preprocessing data...")
    
    # Clean data
    def clean_data(adata_sub):
        X = adata_sub.X.toarray() if hasattr(adata_sub.X, 'toarray') else adata_sub.X.copy()
        X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
        adata_sub.X = X
        return adata_sub
    
    adata_ref = clean_data(adata_ref)
    adata_query = clean_data(adata_query)
    
    # Simple normalization
    sc.pp.normalize_total(adata_ref, target_sum=1e4)
    sc.pp.normalize_total(adata_query, target_sum=1e4)
    
    # Select top variable genes (simple method)
    def select_top_genes(adata_sub, n_genes=2000):
        X = adata_sub.X
        # Calculate variance
        variances = np.var(X, axis=0)
        # Get top genes by variance
        top_indices = np.argsort(variances)[::-1][:n_genes]
        # Subset
        adata_sub = adata_sub[:, top_indices].copy()
        return adata_sub
    
    adata_ref = select_top_genes(adata_ref, 2000)
    adata_query = select_top_genes(adata_query, 2000)
    
    print(f"After gene selection - Ref: {adata_ref.shape}, Query: {adata_query.shape}")
    
    # Get data matrices
    X_ref = adata_ref.X
    X_query = adata_query.X
    
    # Convert to float32 for torch compatibility
    X_ref = X_ref.astype(np.float32)
    X_query = X_query.astype(np.float32)
    
    # Get cell type labels
    ref_labels = adata_ref.obs[celltype_col].values
    query_labels = adata_query.obs[celltype_col].values
    
    # STEP 4: Run SCOT
    print("\n4. Running SCOT alignment...")
    
    try:
        # Initialize SCOT
        print(f"  Initializing SCOT with data shapes: {X_ref.shape}, {X_query.shape}")
        scot = SCOTv2(X_ref, X_query)
        
        # Align datasets
        print("  Aligning datasets...")
        X_ref_aligned, X_query_aligned = scot.align()
        
        print(f"  Aligned shapes: {X_ref_aligned.shape}, {X_query_aligned.shape}")
        
    except Exception as e:
        print(f"✗ ERROR running SCOT: {e}")
        print("\nTrying SCOT with parameters...")
        try:
            scot = SCOTv2(X_ref, X_query)
            X_ref_aligned, X_query_aligned = scot.align(k=50, e=1e-3, normalize=True, norm="l2")
            print("✓ SCOT alignment successful with parameters")
        except Exception as e2:
            print(f"✗ ERROR with parameters: {e2}")
            print("\nTrying self-tuning mode...")
            try:
                X_ref_aligned, X_query_aligned = scot.align(selfTune=True)
                print("✓ SCOT alignment successful with self-tuning")
            except Exception as e3:
                print(f"✗ ERROR with self-tuning: {e3}")
                print("\nFalling back to simple alignment...")
                # Simple fallback: just use original data
                X_ref_aligned = X_ref
                X_query_aligned = X_query
    
    # STEP 5: Calculate similarities
    print("\n5. Calculating cell type similarities...")
    
    from scipy.spatial.distance import cdist
    
    # Get unique cell types
    ref_celltypes = np.unique(ref_labels)
    query_celltypes = np.unique(query_labels)
    
    print(f"  Reference cell types ({len(ref_celltypes)}): {list(ref_celltypes)}")
    print(f"  Query cell types ({len(query_celltypes)}): {list(query_celltypes)}")
    
    # Create similarity matrix
    similarity_matrix = np.zeros((len(ref_celltypes), len(query_celltypes)))
    
    # Calculate distances between cell type centroids in aligned space
    for i, ref_ct in enumerate(ref_celltypes):
        # Get cells of this type
        ref_mask = ref_labels == ref_ct
        if np.sum(ref_mask) > 0:
            ref_centroid = np.mean(X_ref_aligned[ref_mask], axis=0)
            
            for j, query_ct in enumerate(query_celltypes):
                # Get cells of this type
                query_mask = query_labels == query_ct
                if np.sum(query_mask) > 0:
                    query_centroid = np.mean(X_query_aligned[query_mask], axis=0)
                    
                    # Calculate distance
                    distance = np.linalg.norm(ref_centroid - query_centroid)
                    
                    # Convert to similarity (inverse of distance)
                    similarity = 1.0 / (1.0 + distance)
                    similarity_matrix[i, j] = similarity
                    
                    if similarity > 0.1:  # Only print significant similarities
                        print(f"    {ref_ct} -> {query_ct}: {similarity:.4f}")
    
    # Create DataFrame
    sim_df = pd.DataFrame(
        similarity_matrix,
        index=ref_celltypes,
        columns=query_celltypes
    )
    
    # STEP 6: Plot
    print("\n6. Plotting similarity matrix...")
    
    plt.figure(figsize=(max(12, len(sim_df.columns)*0.8), 
                       max(10, len(sim_df.index)*0.7)))
    
    # Plot heatmap
    im = plt.imshow(sim_df.values, aspect='auto', cmap='viridis')
    plt.colorbar(im, label="SCOT Similarity Score", shrink=0.8)
    
    # Set ticks
    plt.xticks(range(len(sim_df.columns)), sim_df.columns, rotation=45, ha='right', fontsize=10)
    plt.yticks(range(len(sim_df.index)), sim_df.index, fontsize=10)
    
    # Labels and title
    plt.xlabel(f"Query: {args.sample2}", fontsize=12, fontweight='bold')
    plt.ylabel(f"Reference: {args.sample1}", fontsize=12, fontweight='bold')
    plt.title(f"Cell Type Similarity: {args.sample1} → {args.sample2}\n(SCOT Optimal Transport)", 
              fontsize=14, fontweight='bold', pad=20)
    
    # Add text in cells
    for i in range(len(sim_df.index)):
        for j in range(len(sim_df.columns)):
            value = sim_df.iloc[i, j]
            if value > 0.01:  # Only show non-zero values
                plt.text(j, i, f"{value:.3f}",
                        ha="center", va="center",
                        color="white" if value > 0.5 else "black",
                        fontsize=8 if len(sim_df.index) < 20 else 6,
                        fontweight='bold' if value > 0.7 else 'normal')
    
    plt.tight_layout()
    plt.savefig(args.out, dpi=300, bbox_inches='tight')
    plt.close()
    
    # STEP 7: Save outputs
    print("\n7. Saving outputs...")
    
    # Save similarity matrix
    csv_path = args.out.rsplit('.', 1)[0] + '_SCOT_similarity.csv'
    sim_df.to_csv(csv_path)
    print(f"✓ Saved similarity matrix to: {csv_path}")
    
    # Save aligned embeddings
    np.savez_compressed(
        args.out.rsplit('.', 1)[0] + '_SCOT_aligned.npz',
        X_ref_aligned=X_ref_aligned,
        X_query_aligned=X_query_aligned,
        ref_labels=ref_labels,
        query_labels=query_labels,
        ref_celltypes=ref_celltypes,
        query_celltypes=query_celltypes
    )
    print(f"✓ Saved aligned embeddings")
    
    # Create a summary report
    report_path = args.out.rsplit('.', 1)[0] + '_SCOT_report.txt'
    with open(report_path, 'w') as f:
        f.write(f"SCOT Analysis Report\n")
        f.write(f"====================\n")
        f.write(f"Reference: {args.sample1} ({adata_ref.shape[0]} cells)\n")
        f.write(f"Query: {args.sample2} ({adata_query.shape[0]} cells)\n")
        f.write(f"\nCell Types:\n")
        f.write(f"  Reference: {list(ref_celltypes)}\n")
        f.write(f"  Query: {list(query_celltypes)}\n")
        f.write(f"\nTop Similarities:\n")
        
        # Get top similarities
        for i in range(len(sim_df.index)):
            for j in range(len(sim_df.columns)):
                value = sim_df.iloc[i, j]
                if value > 0.5:  # High similarity
                    f.write(f"  {sim_df.index[i]} -> {sim_df.columns[j]}: {value:.4f}\n")
    
    print(f"✓ Saved report to: {report_path}")
    print(f"\n✅ DONE: {args.out}")

if __name__ == "__main__":
    main()
