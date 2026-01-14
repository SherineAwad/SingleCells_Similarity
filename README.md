# Interpreting Similarity Analyses Between Cell Types

# 🔬 🔬 🔬Before Assessing Whether Results “Make Sense”

Scientific interpretation requires separating **hypotheses**, **methodological assumptions**, and **observed outcomes**.  
A result should not be judged solely on whether it aligns with prior expectations.

---

## 1. What Is the Hypothesis?

Before applying any similarity or distance metric, we explicitly define:

- What biological notion of *similarity* is being tested  
  (e.g. shared gene expression patterns, co-variation structure, distributional similarity, or relative abundance shifts)

Different hypotheses imply **different mathematical representations** of similarity.

---

## 2. Which Methods Are Appropriate — Based on How They Work (Not on the Results)

Each method probes a distinct aspect of similarity, independent of the observed result:

| Method | What it measures | Key assumptions |
|------|------------------|----------------|
| Pearson correlation | Linear co-variation | Linear relationships, scale sensitivity |
| Spearman correlation | Monotonic rank relationships | Order-based similarity |
| Mutual Information (MI) | Statistical dependence | Captures non-linear relationships |
| Cosine similarity | Directional similarity | Scale-invariant, angle-based |
| Optimal Transport (OT) | Distributional alignment | Geometry and mass conservation |

A method is inappropriate **only if its assumptions conflict with the biological question**, not because the result is unexpected.

---

## 3. Can a Method Falsify the Hypothesis?

For each metric, we ask:

- Given how this method operates, could it detect the hypothesized relationship if it exists?
- If the hypothesis were true, would this metric be sensitive to it?

If the answer is *yes*, then unexpected results are informative rather than invalid.

---

## 4. Interpreting Unexpected Results

Unexpected outcomes may indicate:

1. The hypothesis is incomplete or incorrect  
2. The biological system violates prior assumptions  
3. The metric captures a different but valid notion of similarity  

Unexpected results do **not** imply methodological error by default.

---

## 5. Summary

Methods are tools that test hypotheses under defined assumptions.  
Scientific critique should therefore address:

- the **hypothesis**,  
- the **assumptions of the method**, or  
- the **interpretation of results**,  

rather than dismissing outcomes based on intuition alone.

# Methods in details 


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
**Biological meaning:** How similar are the *relative priorities* of genes?  
**Key point:** Cosine similarity reflects similarity **within the selected gene space**. Changing the gene set can change the biological interpretation ⚠️  

---

### Example 1: Using all genes (global transcriptional similarity)

Control macrophages:  
- Core inflammatory program (stable genes): 60%  
- Stress-induced inflammatory genes: 20%  
- Housekeeping genes: 15%  
- Other genes: 5%  

LD macrophages:  
- Core inflammatory program (stable genes): 55%  
- Stress-induced inflammatory genes: 25%  
- Housekeeping genes: 15%  
- Other genes: 5%  

Cosine similarity ≈ **high** → Overall inflammatory identity and global transcriptional structure are preserved  

**Interpretation:**  
> Macrophages maintain their core transcriptional identity under LD, with moderate stress-related shifts.

---

### ⚠️ Example 2: Using HVGs only (variable-program-focused similarity)

HVG selection preferentially retains **stress-responsive inflammatory genes** and removes:
- Stable inflammatory genes  
- Housekeeping genes  

Remaining gene space (HVGs):

Control macrophages:  
- Stress-induced inflammatory genes: 90%  
- Other variable genes: 10%  

LD macrophages:  
- Stress-induced inflammatory genes: 60%  
- Other variable genes: 40%  

Cosine similarity ≈ **lower** ⚠️  

**Interpretation:**  
> Stress-responsive programs differ substantially between Control and LD macrophages.

---

### 🚨 Why the story changes

- **All genes:** similarity reflects *global cell identity*  
- **HVGs only:** similarity reflects *how variable or stress-driven programs differ*  

> HVG-based cosine similarity emphasizes **change**, not **stability**.

---

### Warning 🚨

> **Cosine similarity is not invariant to gene selection.  
> HVG-based analyses answer a different biological question than all-gene analyses.**  

Both approaches are valid — but they tell **different stories**, and that choice must be stated explicitly.


![Cosine Ctrl vs LD](Cosine_Ctrl_LD.png?v=6)

![Cosine Ctrl vs NMDA](Cosine_Ctrl_NMDA.png?v=6)

![Cosine LD vs NMDA](Cosine_LD_NMDA.png?v=6)
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


# In Plain English 


### Pearson Correlation

**Asks?**

> Do genes increase and decrease together in a straight-line way?

Pearson correlation checks whether **highly expressed genes in one condition are also highly expressed in the other**, assuming a linear relationship.

* High Pearson → Expression changes are proportional
* Low Pearson → Expression changes are inconsistent or nonlinear

**Key intuition:**
Pearson cares about **exact values**, not just ordering.


### Spearman Correlation

**Asks?** 

> Are genes ranked in the same order?

Spearman correlation checks whether **genes keep the same relative ranking**, even if the actual expression values change.

* High Spearman → Gene importance order is preserved
* Low Spearman → Gene ranking is reshuffled

**Key intuition:**
Spearman ignores magnitude and focuses on **rank**.


## Cosine Similarity

**Asks?** 

> Are the same genes important in both conditions?

Cosine similarity looks at whether cells **prioritize the same genes**, regardless of overall expression level.

* High cosine → Same gene priorities → Same identity
* Low cosine → Different genes dominate → State or identity change


## Mutual Information (MI)

**Asks?**

> Do genes still relate to each other in the same way?

MI measures whether **gene–gene relationships are preserved** — when one gene changes, do the same other genes respond?

* High MI → Regulatory logic is intact
* Low MI → Gene coordination is disrupted


### Optimal Transport (OT)

**Asks?**

> *Did the cells move to a different biological state?*

OT measures how far cells must be “moved” to match another condition in biological state space.

* Low transport cost / high similarity → Same state
* High cost / low similarity → Different state








