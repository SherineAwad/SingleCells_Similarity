# Measuring similarities between cell types: other methods than correlations 


## Spearman correlation 

**Pearson correlation** measures global linear similarity across all genes and assumes comparable variance structure and direct gene-wise correspondence. In scRNA-seq data, biologically related cell types can violate these assumptions due to sparsity, nonlinear expression changes, or differences driven by a limited number of genes.

**Spearman correlation measures whether genes maintain their relative importance rankings between conditions.**

- Pearson: "Do they scale together linearly?"

- Spearman: "Do they keep the same ranking order?"

- Different cutoffs → Different actual values → Different Pearson correlation

- Different cutoffs → Same ranking order → Same Spearman correlation

- Spearman is robust to cutoff choices because it only cares about rankings, not exact values.


- A cell type could have high Spearman but low Pearson if genes change by different amounts but keep their importance rankings.

## **Spearman Correlation Key Points:**
- **Doesn't care** about exact expression amounts
- **Only cares** about ranking order (which genes are #1 most important, #2, #3, etc.)
- **High Spearman** = Cell type preserves its gene hierarchy despite changes
- **Low Spearman** = Gene importance rankings got reshuffled by stress

## **Biological Meaning:**
If Spearman is high between Control and LD cells, it means:  
"**This cell type kept its core identity priorities** - the same genes remained most important and least important, even if their actual expression levels changed."

![Spearmnan Control vs LD](Spearman_ALL_Ctrl_LD.png?v=2)

![Spearmnan Control vs NMDA](Spearman_ALL_Ctrl_NMDA.png?v=2)

## Cosine Similarity  
**Measures:** Similarity of *relative gene expression profiles* between cell types or conditions  
**Mechanism:** Computes the cosine of the angle between mean expression vectors in high-dimensional gene space (magnitude-independent)  
**Biological meaning:** How similar are the *relative priorities* of genes, regardless of total expression level?  
**For Control vs LD:** High cosine = Cells preserve transcriptional identity despite global shifts  

**Example:**  
Control macrophages:  
- Inflammatory genes: 80%  
- Housekeeping: 15%  
- Other genes: 5%  

LD macrophages:  
- Inflammatory genes: 75%  
- Housekeeping: 20%  
- Other genes: 5%  

Cosine similarity = **0.94** → Very similar gene ranking and proportional usage despite stress  

![Cosine Ctrl vs LD](Cosine_Ctrl_LD.png?v=3)

![Cosine Ctrl vs NMDA](Cosine_Ctrl_NMDA.png?v=3)

![Cosine LD vs NMDA](Cosine_LD_NMDA.png?v=2)
#### Citation 

Watson, E. R., Mora, A., Taherian Fard, A., & Mar, J. C. (2022). How does the structure of data impact cell–cell similarity? Evaluating how structural properties influence the performance of proximity metrics in single cell RNA-seq data. Briefings in bioinformatics, 23(6), bbac387.



---

## Mutual Information  
**Measures:** Similarity of *statistical dependencies* between genes or gene modules  
**Mechanism:** Quantifies how much knowing the expression of one gene reduces uncertainty about another  
**Biological meaning:** How well are *regulatory or signaling relationships* preserved?  
**For Control vs LD:** High MI = Cells maintain similar regulatory logic  

**Example:**  
Control T cells:  
- IFNγ ↑ → STAT1 consistently ↑ (strong dependency)  

LD T cells:  
- IFNγ ↑ → STAT1 sometimes ↑, sometimes unchanged (weakened dependency)  

Mutual information = **0.35** → Stress disrupts key gene–gene regulatory relationships  


![MI Control vs LD](MI_Ctrl_LD.png?v=3)

![MI Control vs NMDA](MI_Ctrl_NMDA.png?v=1)

#### Citation 
Chang, L. Y., Hao, T. Y., Wang, W. J., & Lin, C. Y. (2024). Inference of single-cell network using mutual information for scRNA-seq data analysis. BMC bioinformatics, 25(Suppl 2), 292.


---

## Optimal Transport (e.g., SCOT)  
**Measures:** Similarity of *cellular state distributions* after dataset alignment  
**Mechanism:** Computes the minimal “cost” required to transport one cell distribution onto another in expression or embedding space  
**Biological meaning:** How geometrically close are cellular states across conditions?  
**For Control vs LD:** Low transport cost / high similarity = Cells occupy similar phenotypic states  

**Example:**  
Control neurons cluster in the “mature excitatory neuron” region  
LD neurons cluster in nearly the same region after alignment  

OT similarity = **0.85** → Minimal shift in cellular state space under stress  


![Scot Control vs LD](SCOT_Ctrl_LD.png?v=3)

![Scot Control vs NMDA](SCOT_Ctrl_NMDA.png?v=1)

#### Citation

Demetci, P., Santorella, R., Sandstede, B., & Singh, R. (2021). Unsupervised integration of single-cell multi-omics datasets with disparities in cell-type representation. BioRxiv, 2021-11.








