# Precision Vascular Disruption: A Network Engineering Approach Based on Rs

This repository contains the official source code and analysis scripts for the paper: **"Precision Vascular Disruption: A Network Engineering Approach Based on Network Resilience ($R_s$) Improves Melanoma Treatment Outcomes by 25%."**

*Authors: Madjid Eshaghi Gordji, et al.*
*Affiliation: Semnan University*
*Preprint/Publication Link: [Link to be added upon publication/preprint]*

## Overview

This research introduces a novel network resilience metric, **Structurally-Weighted Resistance Entropy ($R_s$)**, which unifies spectral and structural properties of complex networks. We demonstrate its power in two ways:
1.  **Theoretical:** We show that maximizing $R_s$ resolves the long-standing Schneider-Gao paradox in network science, producing networks with universal resilience against both targeted attacks and random failures.
2.  **Clinical:** We apply an $R_s$-guided greedy algorithm to personalize laser ablation therapy for melanoma, targeting the tumor vasculature. Our N=20 clinical trial shows that this method achieves a **25% relative improvement in Objective Response Rate (ORR)** while simultaneously reducing collateral tissue damage by **75%** and procedure time by **40%**.

## Repository Structure

This repository is organized to ensure clarity and reproducibility:

-   `/src`: Contains the core Python source code for:
    -   `calculate_Rs.py`: Functions to compute the $R_s$ metric.
    -   `rego_optimizer.py`: Implementation of the REGO optimization algorithm.
    -   `clinical_analysis.py`: Scripts for statistical analysis of clinical data.

-   `/data`: Contains sample data to demonstrate script functionality.
    -   **Note: Clinical data is anonymized and synthetic to protect patient privacy and comply with ethical standards.**

-   `/notebooks`: Jupyter notebooks that reproduce the main figures and results of the paper.
    -   `01_In_Silico_Simulations.ipynb`: Generates Figure 1 (Theoretical Results).
    -   `02_Clinical_Results_Analysis.ipynb`: Generates Figure 2 (Clinical Results).

## How to Use This Repository

The code is provided for transparency and reproducibility. The analysis can be fully replicated by running the Jupyter Notebooks in the `/notebooks` directory, which depend on the functions defined in `/src`.

A `requirements.txt` file is included, listing all necessary Python packages for execution.

## Citation

If you use this code or our findings in your research, please cite our paper. A BibTeX entry will be provided upon publication.
```bibtex
@article{Gordji2025Precision,
  title   = {Precision Vascular Disruption: A Network Engineering Approach Based on Network Resilience ($R_s$) Improves Melanoma Treatment Outcomes by 25\%},
  author  = {Gordji, Madjid Eshaghi and [Other Authors]},
  journal = {[Journal Name to be added]},
  year    = {2025},
  volume  = {[Volume to be added]},
  pages   = {[Pages to be added]}
}

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
`
