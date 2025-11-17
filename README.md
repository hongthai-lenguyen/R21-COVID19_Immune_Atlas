# A High-Resolution Computational Atlas of Immune Plasticity to Define Immuno-Signatures of COVID-19 Severity and Convalescence

---
## Project overview

**Goal:** Build a high-resolution, reproducible computational atlas of circulating immune cells in acute and convalescent COVID-19 patients to (1) discover rare and transitional cellular states, (2) derive predictive multi-dimensional immuno-signatures of severity, and (3) infer pseudotemporal trajectories that expose checkpoints of incomplete immune recovery (PASC hypothesis). This is a hypothesis-driven, exploratory R21 style re-analysis of a publicly available, high-value dataset. 

**Why re-analyse?** The original analysis used visualization-first methods (viSNE/t-SNE) and clustered on 2D maps, which can distort global topology and obscure continuous phenotypic trajectories. Modern workflows (high-resolution clustering + topology-aware manifold/graph learning) can reveal continuous differentiation paths and rare subpopulations missed previously. (See Methods & citations below.)

**Note:** Raw FCS files are available at FlowRepository (do not re-host without permission).

---

## Key references & data access
- Mathew D. _et al._, **Deep immune profiling of COVID-19 patients**. _Science_ 2020 — foundational dataset and cohort description.
- Raw 28-color FCS files and clinical metadata: **FlowRepository ID: FR-FCM-Z6KL** (downloadable repository).
- `scCODA` (Bayesian compositional model for single-cell compositional analysis).
- `diffcyt` (high-resolution differential discovery framework for cytometry: DS/DA frameworks).
- Trajectory / graph tools: **Monocle 3** and **PAGA** (trajectory and graph abstraction approaches).

---
## Specific aims (concise)
1. **Aim 1 — Atlas:** Build an unsupervised, high-resolution atlas using FlowSOM directly on compensated 28-dimensional data + UMAP for visualization and topology preservation.
2. **Aim 2 — Signatures:** Derive predictive immuno-signatures using compositional DA (`scCODA`) and DS (`diffcyt`) analyses versus NIH/WHO 8-point severity.
3. **Aim 3 — Trajectories:** Infer pseudotemporal differentiation/activation trajectories (Monocle 3 / PAGA) to detect "stalled" convalescent states consistent with PASC-like immune dysregulation.

A comparative table (short) is provided in the repository `docs/` for grant-review readability.

---

## Pipeline summary (step-by-step)
1. **Data acquisition**
    - Download FCS files + metadata from FlowRepository (FR-FCM-Z6KL).
2. **Pre-processing**
    - Compensation (if needed), debris/doublet removal, viability gating, channel QC.
    - Transformation: apply $\mathrm{arcsinh}$ transform, e.g. $\mathrm{arcsinh}(x / 5)$ for each protein channel.
3. **High-resolution clustering (Aim 1)**
    - Run **FlowSOM** on the full 28-dimensional compensated+transformed matrix (no 2D pre-clustering).
    - Over-cluster (e.g., 20×20 SOM → 400 micro-clusters), then meta-cluster (hierarchical consensus to ~40–50 meta-clusters).
4. **Manifold & visualization**
    - Compute **UMAP** coordinates for visualization (maintain reproducible seeds + multiple initializations). (Use UMAP for interpretability; report limitations and compare with t-SNE results).
5. **Differential analyses (Aim 2)**
    - **DA:** use `scCODA` (Bayesian compositional) on the subject-level cell-type composition matrix; enforce compositional constraint $\sum_i p_i = 1$.
    - **DS:** within meta-clusters, use `diffcyt` to test marker expression changes vs severity and covariates.
6. **Trajectory inference (Aim 3)**
    - Subset T-cell and myeloid compartments; run **Monocle 3** `learn_graph()` and/or **PAGA** to infer connected graph structure; root trajectories using biologically defined naive populations.
7. **Validation & robustness**
    - Internal cross-validation, bootstrapping, sensitivity to clustering resolution, alternative manifold initializations, and multiple pseudotime initial roots.
8. **Packaging & release**
    - Release processed atlas (UMAP coordinates, FlowSOM labels), statistical results, reproducible R/Python notebooks, and an interactive HTML report.

---
## Methods & rationale — per aim (detailed)
### Aim 1 — Atlas (why and how)
- **Why:** Overcome topological distortions introduced by clustering on 2D viSNE maps; preserve continuous phenotypic structure (important for plasticity).
- **How (summary):** FlowSOM on the original 28 channels → over-clustering → meta-clusters → annotate via median marker profiles → visualize on UMAP; produce per-sample cell-type frequency matrix.

**Key parameters (suggested starting values):**
- FlowSOM: grid 20×20 (400 micro-clusters), metaclustering k ≈ 40–50.
- Transformation cofactor: 5 (standard for protein cytometry). Use $\mathrm{arcsinh}(x/5)$.
- UMAP: n_neighbors = 15–30 (sensitivity analysis), metric = euclidean, seed for reproducibility.

---
### Aim 2 — Immuno-signatures (DA + DS)
- **Compositional DA (scCODA):** model subject-level composition vector $\mathbf{p}=(p_1,\dots,p_K)$ with $\sum p_i=1$; use hierarchical Bayesian inference to identify cell types with credible changes when comparing Severe vs Moderate while controlling for covariates (age, sex, batch). `scCODA` was explicitly developed for this problem.
- **Differential State (diffcyt):** for each meta-cluster, test marker expression (median or robust summaries) versus severity using empirical Bayes moderated tests; controls false positive rates for rare populations.

**Important statistical notes:** treat counts as compositional, include covariates, correct for multiple hypothesis testing (Bayesian credible intervals + FDR where applicable), and perform permutation controls.

---

### Aim 3 — Trajectories (Monocle 3 / PAGA)
- **Rationale:** Cross-sectional data can be used to infer pseudotemporal progress through a biological process by arranging cells along learned graph structures; this can highlight whether convalescent cells "return" to healthy regions or remain stalled.
- **Workflow:** subset compartments → compute neighborhood graph → `learn_graph()` / PAGA → assign pseudotime via rooted graph → compare pseudotime distributions by clinical group (Healthy / Acute-Moderate / Acute-Severe / Convalescent).

**Caveat:** Pseudotime inferred from cross-sectional data is hypothesis-generating and requires careful robustness checks (alternate roots, bootstrapped graphs, independent validation cohorts where possible).

---
## Statistical considerations & notes
- **Compositional constraint:** cell type proportions sum to 1 per sample: $\sum_{i=1}^K p_i = 1$. Use compositional methods (scCODA) rather than naive tests on proportions.
- **Multiple testing:** use hierarchical Bayesian shrinkage (scCODA) and empirical Bayes (diffcyt).
- **Batch effects & confounders:** include sample collection batch, processing date, and donor metadata as covariates.
- **Rare populations:** over-clustering then meta-clustering helps preserve rare subsets; keep resolution sensitivity analysis.
- **Visualization pitfalls:** UMAP and t-SNE are visualization tools—interpret topology cautiously and report reproducibility across seeds and initializations.

---
## 🧑‍💻 Author and Contact
**Maintainer:** [HongThai LeNguyen]  
**Role:** Data Scientist & Bioinformatician  
**Institution:** *OUCRU*  
**Email:** *hongthai.lenguyen@outlook.com/thailnh@oucru.org*  

---
## 🧾 License
This project is licensed under the **[APACHE License](LICENSE)**.
  
---
## 🔬 Keywords

`COVID-19` · `SARS-CoV-2` · `immune atlas` · `immune plasticity` · `high-dimensional cytometry` · `flow cytometry` · `28-color panel` · `single-cell analysis` · `immuno-signatures` · `acute severity` · `convalescence` · `PASC` · `FlowSOM` · `UMAP` · `viSNE` · `differential abundance` · `differential state` · `scCODA` · `diffcyt` · `Monocle 3` · `PAGA` · `pseudotime trajectory` · `topology-aware analysis` · `high-resolution clustering` · `compositional data analysis` · `Mathew et al. 2020` · `FlowRepository FR-FCM-Z6KL` · `immunophenotyping` · `rare cell populations` · `meta-clustering`
