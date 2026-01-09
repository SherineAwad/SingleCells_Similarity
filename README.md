# Measuring similarities between cell types: other methods than correlations 


**Pearson correlation** measures global linear similarity across all genes and assumes comparable variance structure and direct gene-wise correspondence. In scRNA-seq data, biologically related cell types can violate these assumptions due to sparsity, nonlinear expression changes, or differences driven by a limited number of genes.


## Method 1:  OT using SCOT 


To relax these assumptions, we applied **Optimal Transport (OT)**–based alignment (SCOT). OT is more permissive than correlation, as it allows redistribution of expression mass and does not require one-to-one gene or cell matching. Conceptually, OT tests whether one population can be smoothly transformed into another in expression space.


SCOT applies Optimal Transport to align and compare single-cell datasets by preserving the local structure of cell populations. It uses OT to find the best way to match cells across conditions or samples based on their neighborhood relationships, enabling meaningful comparison of cell types without relying on shared labels or markers.


## What is Optimal Transport (OT)

Optimal Transport compares cell types by treating each as a cloud of cells and measuring how much “effort” it takes to move one cloud to match the other. Small effort means similar cell types; large effort means different.


## How OT Measures Similarity Between Cell Types in Single-Cell Data

In single-cell data, OT calculates how much cell “mass” from one type needs to be moved and reshaped to match another type’s distribution. The less movement required, the more similar the cell types are.


## How SCOT Uses OT

SCOT applies Optimal Transport to align and compare single-cell datasets by preserving the local structure of cell populations. It uses OT to find the best way to match cells across conditions or samples based on their neighborhood relationships, enabling meaningful comparison of cell types without relying on shared labels or markers.



## Using Zebrafish as example 


![Scot Control vs LD](SCOT_Ctrl_LD.png?v=1)


## SCOT Version

The similarity analysis was performed using **SCOT version 2** (SCOTv2), a Python tool for unsupervised alignment of single-cell multi-omics datasets.

### Citation



> Demetci, P., Santorella, R., Sandstede, B., Noble, W. S., and Singh, R. (2021).  
> Unsupervised integration of single-cell multi-omics datasets with disparities in cell-type representation.  
> bioRxiv 2021.11.09.467903.





2. Method 2: Mutual Information 

# Mutual Information (MI) for Measuring Similarity

Mutual Information (MI) is an information-theoretic measure that quantifies the amount of shared information between two variables. Unlike correlation measures, which often assume linear relationships, MI captures **both linear and nonlinear dependencies** between datasets.

## How MI Captures Similarity in Single-Cell Data

- MI evaluates the dependency between gene expression distributions across two samples or conditions.
- It measures how much knowledge of one gene’s expression in one dataset reduces the uncertainty about its expression in another dataset.
- This makes MI well-suited for detecting complex, nonlinear relationships and shared patterns in single-cell RNA-seq data, where gene expression changes may not be strictly linear.
- By aggregating MI scores across genes and cell types, one can build a similarity matrix that reflects the overall biological similarity between cell populations.

## Advantages of Using MI

- **Nonlinear relationships:** Captures complex patterns missed by Pearson or Spearman correlations.
- **Distribution-level similarity:** Goes beyond pointwise comparison to consider the overall information shared between datasets.
- **Robust to data transformations:** Works well with normalized or discretized data.


![MI](MI_Ctrl_LD.png?v=1)
