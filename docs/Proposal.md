# A High-Resolution Computational Atlas of Immune Plasticity to Define Immuno-Signatures of COVID-19 Severity and Convalescence

## Introduction & Rationale: The Unresolved Granularity of COVID-19 Immune Dysregulation

The SARS-CoV-2 pandemic has been characterized by profound clinical heterogeneity, with patient outcomes ranging from asymptomatic infection to severe acute respiratory distress syndrome (ARDS), multi-organ failure, and death. A central, unresolved question in infectious disease is to define the immunological mechanisms that determine this severe, life-threatening trajectory.

While acute disease remains a challenge, the most pressing, long-term clinical and economic crisis is the high prevalence of Post-Acute Sequelae of COVID-19 (PASC), or Long COVID. PASC is increasingly understood as a syndrome of persistent, systemic immune dysregulation. Seminal 2024 studies have provided a mechanistic framework for PASC, identifying "systemic inflammation and immune dysregulation" characterized by "exhausted SARS-CoV-2-specific $\text{CD8}^+\text{ T}$ cells," "T cell dysregulation," and a "mis-coordination between their $\text{SARS-CoV-2}$-specific T and B cell responses".1 This is compounded by "highly activated myeloid cells" and "altered state[s] of mitochondrial function in myeloid and lymphoid cells", which persist long after viral clearance. The factors distinguishing complete immune resolution from this chronic, unresolved inflammatory state remain dangerously ill-defined.2

In 2020, a foundational study by Mathew et al. (_Science_, 2020) provided the first deep immune atlas of acute COVID-19. Using 28-color flow cytometry on a large cohort of 125 acute hospitalized patients, 36 convalescent donors, and 60 healthy donors, the authors defined three broad "immunotypes" associated with clinical outcomes. This study established a critical link between T-cell activation, plasmablast responses, and disease severity. The full dataset was made publicly available (FlowRepository ID: FR-FCM-Z6KL), creating an invaluable resource for the scientific community.

However, this foundational analysis, while groundbreaking for its time, was limited by its computational methodology and the nascent understanding of COVID-19 immunopathology. The original analysis relied on viSNE (a t-SNE variant) for dimensionality reduction.3 It is now well-established that t-SNE and its variants are visualization algorithms that, while excellent at preserving local neighborhoods, fundamentally "struggle to preserve global structure".4 These methods can be "misleading" by creating "spurious clusters" and artifactually separating continuous biological processes into discrete, disconnected "islands".4 The original analysis appears to have run clustering _on the 2D viSNE map_ 3, thereby baking these topological distortions into their final model. This approach, by its very nature, obscures the precise biological phenomena we hypothesize are critical: "phenotypic plasticity" and "co-expression dynamics." The continuous differentiation pathways of T-cell exhaustion and myeloid activation—which do not exist as discrete islands but as phenotypic continua—are lost.

This proposal is founded on an R21-level opportunity: a novel, exploratory re-analysis of this high-value public dataset using a state-of-the-art computational pipeline that _corrects_ for these original limitations. This project aligns perfectly with NIH priorities for secondary analysis of existing human datasets 5 and the development of "innovative analytical methodologies" to address the COVID-19 crisis. We hypothesize that the original "immunotypes" are coarse bins that average out the most critical, predictive biological signals. We posit that **specific, rare cell subsets and complex, continuous marker co-expression patterns (immuno-signatures) within the myeloid and exhausted T-cell compartments—visible only through a high-resolution, topologically-aware computational framework—are predictive of acute disease severity and the incomplete immune resolution characteristic of PASC.**

The goal of this R21 project is to develop and apply this novel computational pipeline to the Mathew et al. dataset. We will re-map the entire immune atlas, define new granular immuno-signatures of severity, and, for the first time, model the pseudotime-based differentiation trajectories that distinguish acute, healthy, and PASC-like convalescent states.

## Specific Aims: A Novel Computational Workflow for Immune Re-stratification

This project will acquire the raw FCS files and clinical metadata for all 221 subjects from the Mathew et al. _Science_ 2020 study (FlowRepository ID: FR-FCM-Z6KL). We will deploy an integrated computational workflow designed to test our central hypothesis. The novelty of this R21 project lies in the application of this modern pipeline to extract latent biological information invisible to the original study's methods.5

**Table 1: Novel Computational Strategy and Rationale for Re-Analysis**

| **Aim**                  | **Core Question**                                                              | **Proposed Method(s)**             | **Original Method (Mathew et al.)**          | **Rationale & Novelty (with Citation Justification)**                                                                                                                                                                                                                           |
| ------------------------ | ------------------------------------------------------------------------------ | ---------------------------------- | -------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Aim 1** (Atlas)        | Can we identify previously unclassified rare and transitional cell states?     | **FlowSOM + UMAP**                 | viSNE + FlowSOM-on-viSNE 3                   | **Novelty:** Corrects a flawed workflow. UMAP preserves _global_ topology lost by viSNE 4, which is essential for mapping lineage continua and the _true_ relationships between clusters.                                                                                       |
| **Aim 2** (Signatures)   | What specific cell subsets and co-expression profiles define "Severe" disease? | **`diffcyt` (DS) + `scCODA` (DA)** | UMAP + Correlation of coarse "Immunotypes" 3 | **Novelty:** Moves from coarse "immunotypes" to granular, _predictive_ "immuno-signatures." `scCODA` is a Bayesian compositional model that avoids false positives. `diffcyt` enables Differential State (DS) analysis 6, testing the PI's "co-expression" hypothesis directly. |
| **Aim 3** (Trajectories) | Do "Severe" patients fail to resolve immune activation, leading to PASC?       | **Monocle 3 / PAGA**               | None (static analysis only)                  | **Novelty:** _The most exploratory R21 component._ We will infer _dynamics_ (pseudotime) from _static_ cross-sectional data 7 to model the _process_ of exhaustion and myeloid activation, providing a mechanistic model for PASC-related immune dysregulation.                 |

### Aim 1: To construct a high-resolution, unsupervised atlas of the circulating immune system in acute and convalescent COVID-19.

**Rationale:** To test our hypothesis, we must first generate a new, high-fidelity immune atlas that is free from the topological distortions of the original viSNE-based analysis.4 The original method, by clustering on a 2D visualization 3, aliased the high-dimensional data and obscured the very "phenotypic plasticity" we aim to discover. Our Aim 1 atlas will serve as the foundational, non-parametric map for the quantitative analyses in Aims 2 and 3.

**Methodology:**

1. **Data Acquisition and Pre-processing:** All 28-color FCS files from the FR-FCM-Z6KL dataset will be downloaded. Data will be rigorously pre-processed, including compensation spillover-correction (if not already applied), removal of debris, doublets, and non-viable cells (based on viability dye), and a standard $arcsinh$ transformation (e.g., cofactor 5) applied to all 28 protein channels.
    
2. **Unsupervised Clustering:** We will utilize **FlowSOM** (Self-Organizing Map). Critically, and in direct contrast to the original workflow 3, FlowSOM will be run _directly on the 28-dimensional compensated data_. This is the correct, robust application of the algorithm, ensuring we are clustering the true, high-dimensional cellular phenotypes. We will intentionally "over-cluster" the data by generating a high-resolution SOM (e.g., 20x20 grid, 400 micro-clusters), which will then be meta-clustered using hierarchical consensus clustering into approximately 40-50 high-fidelity meta-clusters. This over-clustering is essential for capturing the _rare cell subsets_ at the core of our hypothesis, as these populations are often averaged-out by lower-resolution methods.
    
3. **Visualization and Manifold Generation:** We will use **UMAP (Uniform Manifold Approximation and Projection)** to generate a new 2D visualization manifold. Unlike t-SNE/viSNE, UMAP is "generally better at preserving the global structure" of data. This property is essential for our hypothesis, as it will spatially organize clusters on the 2D map in a manner that reflects their true high-dimensional phenotypic relationships. For example, a continuous differentiation path from a naive to an effector cell will appear as a "bridge" of cells, rather than two disconnected "islands".4
    

**Expected Outcome:** A rich, high-resolution atlas that visualizes all 221 patient samples (125 acute, 36 convalescent, 60 healthy) in a single, unified topological space. This map, annotated with our 40-50 FlowSOM meta-clusters, will reveal the _true_ phenotypic landscape and the _relationships_ between cell states, providing the high-fidelity foundation for all subsequent Aims.

### Aim 2: To define robust, predictive "immuno-signatures" of acute disease severity.

**Rationale:** The original study's three "immunotypes" were descriptive, coarse bins. We will replace this model with a granular, quantitative, and _predictive_ "immuno-signature." We will test the association between the high-resolution immune features from Aim 1 and the clinical severity metadata collected in the original study. Our hypothesis is two-fold: that severity is driven by changes in the _abundance_ of rare cell subsets (Differential Abundance, DA) and by changes in the _phenotypic state_ (marker co-expression) of common cell subsets (Differential State, DS).

**Clinical Covariate:** We will utilize the **NIH/WHO 8-point ordinal scale** for disease severity, which was collected in the original study.3 This is a robust, standardized metric used in major COVID-19 clinical trials. We will stratify the 125 acute patients into "Moderate" (e.g., Scale 3-4: Hospitalized, no supplemental oxygen or oxygen by mask/nasal prongs) and "Severe" (e.g., Scale 5-7: Hospitalized, requiring non-invasive ventilation, high-flow oxygen, intubation, or ECMO).

**Methodology: A Two-Pronged Statistical Approach:**

1. **Differential Abundance (DA) Analysis with `scCODA`:**
    
    - **The Challenge:** Flow cytometry data is _compositional_. A massive increase in one population (e.g., the >30% plasmablast response noted in some patients) causes _relative_ decreases in all other populations, even if their absolute cell counts are stable. Standard statistical tests (e.g., t-test, Mann-Whitney) on these relative frequencies are statistically invalid and will generate a storm of false positives.
        
    - **Our Solution:** We will employ **`scCODA` (Compositional Data Analysis)**. `scCODA` is a "Bayesian model for compositional single-cell data analysis" that builds a hierarchical model to identify credible changes in cell-type abundance, all while properly accounting for the compositional nature of the data. Its power is established: `scCODA` has been shown to "detect cell-type changes in COVID-19 patients that were not detected with non-compositional tests". This is precisely the tool required to test our hypothesis about "rare cell subsets" [User Query].
        
2. **Differential State (DS) Analysis with `diffcyt`:**
    
    - **The Rationale:** This arm tests the "co-expression dynamics" portion of our hypothesis. We posit that severity is not just defined by _which cells_ are present, but by _how they are behaving_.
        
    - **Our Solution:** We will use the **`diffcyt`** framework. `diffcyt` is a "new computational framework for differential discovery... based on a combination of high-resolution clustering and empirical Bayes moderated tests".6 Its major advantage is the formal implementation of Differential State (DS) analysis. For each meta-cluster identified in Aim 1 (e.g., "Non-naive $\text{CD8}^+\text{ T}$-cells"), `diffcyt` will build a robust statistical model to test if the median expression of _state markers_ (e.g., PD-1, CD39, KI67, HLA-DR, as defined in the 28-color panel) within that cluster is significantly correlated with the NIH severity score. This method has "improved statistical performance, including for rare cell populations".6
        

**Expected Outcome:** A validated "immuno-signature" composed of both DA cell subsets (e.g., "a 1.5-fold credible increase in $\text{CD14}^+\text{CD16}^+$ myeloid subset X in Severe vs. Moderate") and DS marker profiles (e.g., "a significant upregulation of PD-1 and CD39 on $\text{CD8}^+\text{T}_{\text{EMRA}}$ cells in Severe vs. Moderate"). This signature will be a multi-dimensional, predictive model of severity, far surpassing the original coarse "immunotypes."

### Aim 3: To map the longitudinal trajectories of severity-associated immune populations from the acute to the convalescent phase.

**Rationale:** This is the project's most exploratory, innovative, and high-impact Aim, fully embracing the R21 mechanism. It directly addresses the PASC hypothesis by modeling the "phenotypic plasticity" of the immune system. The Mathew et al. dataset serendipitously contains a cohort of 36 "convalescent donors" sampled in 2020. We can now re-interrogate this cohort through the modern (2024/2025) lens of PASC immunology. We hypothesize that "Severe" acute disease permanently alters immune differentiation pathways, and this perturbation is visible as "incomplete recovery" in the convalescent cohort.

The Challenge & Solution: Trajectory Inference (TI) from Cross-Sectional Data:

We do not have longitudinal time-series data for individual patients. We have cross-sectional snapshots from three distinct groups (Healthy, Acute, Convalescent). To model dynamics from this static data, we will use Trajectory Inference (TI) algorithms.

- **Methodology:** We will employ **Monocle 3** and/or **PAGA**. These are state-of-the-art algorithms that "infer lineage relationships between clusters, _without a time series_".7 Monocle introduced the concept of **"pseudotime"**—a measure of progress through a biological process (e.g., differentiation, activation). It works by learning a "principal graph" that connects phenotypically similar clusters (from Aim 1) to form a continuous trajectory.7
    
- **Application:**
    
    1. **Subset Data:** We will create data subsets from our Aim 1 atlas focused on the key populations of interest: (i) the T-cell compartment (to model exhaustion) and (ii) the myeloid/monocyte compartment (to model activation/plasticity).
        
    2. **Learn Graph:** We will run Monocle 3's `learn_graph()` function 7 on these subsets, _simultaneously_ using all cells from the Healthy (n=60), Acute (n=125), and Convalescent (n=36) groups.
        
    3. **Root Trajectory:** We will use biological knowledge to "root" the trajectory, or define its starting point. For example, the "Naive $\text{CD8}^+\text{ T}$-cell" meta-cluster will be defined as the start ($pseudotime=0$) of the T-cell exhaustion/activation pathway.
        
    4. **Analyze Trajectories:** The algorithm will assign a pseudotime value to every cell. We will then de-convolve the graph and visualize the distribution of our 3 clinical groups (Healthy, Acute, Convalescent) along the trajectory.
        

**Hypothesized Outcome (The "PASC Map"):** We predict the resulting trajectories will be a "smoking gun" for PASC mechanisms. For the T-cell trajectory, we expect:

- **Healthy** cells will map almost exclusively near the naive root.
    
- **Acute-Severe** cells will map predominantly to the _end_ of the trajectory, a state defined by terminal exhaustion (e.g., high PD-1, CD39, T-bet expression).
    
- **Convalescent (PASC-proxy)** cells will be the critical test. We predict they will _not_ map back to the "Healthy" root. Instead, they will be **"stalled"** in intermediate or late-stage regions of the trajectory, co-localizing with the "exhausted" $\text{CD8}^+\text{ T}$ cell phenotypes 1 and "activated myeloid" phenotypes that define PASC.
    

**Expected Outcome:** A set of trajectory-based models that mechanistically visualize _how_ severe COVID-19 perturbs immune differentiation and _where_ PASC patients "get stuck." This provides the first map of the "checkpoints of incomplete immune recovery" [User Query], directly identifying the cellular states that define persistent immune dysregulation.

## Expected Outcomes & Impact

Expected Outcomes:

This 2-year R21 project is exploratory and will produce three major, tangible outcomes:

1. **A New, High-Resolution Immune Atlas:** We will deliver a new, definitive, high-resolution atlas (e.g., UMAP coordinates, FlowSOM cluster definitions, annotated cell populations) of the 221-subject Mathew et al. cohort. This atlas will correct the topological limitations of the original analysis 4 and will be shared as a new public resource, enabling further discovery by the community.
    
2. **Validated "Immuno-Signatures" of Severity:** We will deliver a validated, quantitative "immuno-signature"—a panel of DA cell subsets (from `scCODA`) and DS marker profiles (from `diffcyt`)—that is strongly predictive of the NIH 8-point ordinal severity score. This moves beyond the original _descriptive_ immunotypes to a _predictive_ prognostic tool.
    
3. **Mechanistic Trajectory Models of PASC:** We will generate the first _pseudotime-based trajectory models_ 7 of T-cell exhaustion and myeloid differentiation using cross-sectional acute, healthy, and convalescent data. This directly tests the "phenotypic plasticity" hypothesis and provides a powerful new model for understanding PASC.
    

Impact and Significance:

This project will, if successful, "challenge... and seek to shift current research... paradigms", fulfilling the core mission of the R21 mechanism.

- **Paradigm Shift:** This work will shift the field's understanding of COVID-19 immunopathology from _static, coarse "immunotypes"_ to a _dynamic, continuous, and high-resolution model_ of immune dysregulation. By demonstrating the power of this modern computational pipeline, this project will establish a new analytical standard for high-dimensional cytometry in infectious disease.
    
- **Prognostic Impact:** The "immuno-signature" from Aim 2 has clear translational potential. It represents a new, high-dimensional biomarker set that could, with further validation, be adapted to a clinical setting to risk-stratify acute COVID-19 patients, identifying those who would benefit most from immunomodulatory therapy.
    
- **Therapeutic Impact for PASC:** This is the project's most significant potential impact. The "stalled" cellular states and "checkpoints of incomplete immune recovery" identified in Aim 3 are _novel therapeutic targets_ for PASC. By pinpointing the _exact_ phenotypic state where the immune system of convalescent patients "gets stuck" 1, this work will provide a rational, mechanistic basis for designing new therapies (e.g., checkpoint inhibitors, metabolic modulators) to "un-stick" these cells and restore immune homeostasis, directly addressing the urgent need for targeted PASC interventions.2
    
- **Future Directions:** This 2-year R21 project will establish the feasibility of this computational approach and generate the critical preliminary data required for a subsequent, larger R01 application. That R01 will propose to validate these findings in new, prospective PASC cohorts and functionally test the therapeutic pathways identified in this exploratory analysis.