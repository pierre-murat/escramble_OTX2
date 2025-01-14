# Enhancer SCRaMbLE - OTX2 Super-Enhancer

This repository contains the source code used to reproduce the figures from the paper:

**Resolution of a Human Super-Enhancer by Targeted Genome Randomization**  
Authors: J. Koeppel, P. Murat, G. Girling, E. Madli Peets, M. Gouley, V. Rebernig, A. Maheshwari, J. Hepkema, J. Weller, J. H. Johnkingsly Jebaraj, R. Crawford, F. G. Liberante, L. Parts  

Developers:  
- J. Koeppel ([jk24@sanger.ac.uk](mailto:jk24@sanger.ac.uk))  
- P. Murat ([pm23@sanger.ac.uk](mailto:pm23@sanger.ac.uk))  
- A. Maheshwari ([am86@sanger.ac.uk](mailto:am86@sanger.ac.uk))  
- J. Hepkema ([jh47@sanger.ac.uk](mailto:jh47@sanger.ac.uk))  

## Repository Content

This repository includes bash scripts, Python codes, and R notebooks for the reproduction of the main and supplementary figures of the study. Below is an outline of the folder structure and corresponding analyses:

- **`001_cas9_sequencing_variants_call`**  
  Structural variant identification from Cas9 sequencing for OTX2-loxp3 and OTX2-loxp6 cell lines.  
  _Figures: Main Fig. 2 and Supplementary Fig. 4_

- **`002_cas9_sequencing_variants_analysis`**  
  Analysis of structural architectures generated from OTX2-loxp3 and OTX2-loxp6 cell lines identified via Cas9 sequencing.  
  _Figures: Main Fig. 2 and Supplementary Figs. 4, 5_

- **`003_Fiber-seq`**  
  Analysis of m6A (DNA accessibility) and meCpG modifications from long-read sequencing data.  
  _Figures: Main Fig. 3 and Supplementary Fig. 6_

- **`004_loxP6_pool_vs_clonal_comparison`**  
  Comparative analysis of architectures from the OTX2-loxp6 cell line and a heterogeneous cell population with OTX2-loxp6 design.  
  _Figures: Supplementary Figs. 7, 8_

- **`007_loxP7_architectures_analysis`**  
  Analysis of architectures generated in a heterogeneous cell population following the OTX2-loxp7 design.  
  _Figures: Main Figs. 4, 5 and Supplementary Figs. 11, 12, 13_

- **`005_PCR_sequencing_variants_call_analysis`**  
  Analysis of structural architectures from OTX2-loxp3 cell line identified via PCR sequencing.  
  _Figures: Supplementary Fig. 9_

- **`006_Custom_SV_caller`**  
  Custom structural variant caller for architecture analysis.  
  _Figures: Main Figs. 4, 5 and Supplementary Fig. 11_

- **`008_Prediction`**  
  Super-enhancer OTX2 interaction analysis using ABC and rE2G models, and gene expression prediction via the Borzoi model.  
  _Figures: Supplementary Fig. 14_

- **`009_TFBS_deletion_screen`**  
  Data processing and analysis for transcription factor binding site deletions.  
  _Figures: Main Fig. 6 and Supplementary Fig. 15_

- **`Source`**  
  Additional supporting scripts
---

For detailed usage instructions and further information, please refer to the scripts and documentation within each folder.