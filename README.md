# Measuring similarities between cell types: other methods than correlations 


## Method 1:  OT using SCOT 

SCOT applies Optimal Transport to align and compare single-cell datasets by preserving the local structure of cell populations. It uses OT to find the best way to match cells across conditions or samples based on their neighborhood relationships, enabling meaningful comparison of cell types without relying on shared labels or markers.


## What is Optimal Transport (OT)

Optimal Transport compares cell types by treating each as a cloud of cells and measuring how much “effort” it takes to move one cloud to match the other. Small effort means similar cell types; large effort means different.


## How OT Measures Similarity Between Cell Types in Single-Cell Data

In single-cell data, OT calculates how much cell “mass” from one type needs to be moved and reshaped to match another type’s distribution. The less movement required, the more similar the cell types are.


## How SCOT Uses OT

SCOT applies Optimal Transport to align and compare single-cell datasets by preserving the local structure of cell populations. It uses OT to find the best way to match cells across conditions or samples based on their neighborhood relationships, enabling meaningful comparison of cell types without relying on shared labels or markers.



## Using Zebrafish as example 


![Scot Control vs LD](SCOT_Ctrl_LD.png?v=1)

