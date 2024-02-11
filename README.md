# CRISPRi in Zymomonas mobilis

This repository encompasses the core work presented in [Genetics of Aerotolerant Growth in an Alphaproteobacterium with a Reduced Genome](https://journals.asm.org/doi/10.1128/mbio.01487-23), and it is now evolving to include our latest research efforts. It serves as a foundation for both reproducing the published findings and exploring new directions in the study of *Zymomonas mobilis* through CRISPR interference (CRISPRi).

## Ongoing Research and Novel Techniques
We are developing and applying a novel analytical framework that integrates Limma-voom with STRING database insights and machine learning techniques. This innovative approach is designed to enhance gene set enrichment analysis by examining the effects of gene knockdowns in competition assays. Our aim is to refine the understanding of gene function and interaction networks within *Zymomonas mobilis*.

- **Limma-voom Integration**: Utilizes for robust normalization and differential expression analysis, preparing the data for deeper exploration.
- **STRING Database Analysis**: Employs for mapping gene interactions and networks, providing a context for the observed changes.
- **Machine Learning Algorithms**: Applied to identify patterns and predict outcomes from complex datasets, focusing on gene set enrichment and interaction effects.

This repository not only supports the reproducibility of our published work but also acts as a dynamic platform for sharing our ongoing research advancements.

## New additions:
- Leveraging Camera (from Limma) and Voom for normalization, alongside ggplot2 for violin and sina plots, to assess relative shifts in *Zymomonas mobilis* CRISPRi library composition, integrating [STRINGdb](https://string-db.org) for enriched gene interaction analysis. [Limma-Voom](https://bioconductor.org/packages/release/bioc/html/limma.html). *See* [this script](https://github.com/ryandward/zymomonas_mobilis_CRISPRi_seq/blob/main/enrichment_tests.r)
  
  ![image](https://github.com/ryandward/zymomonas_mobilis_CRISPRi_seq/assets/6970996/b51c7c8a-793d-4c98-ab31-bcba984f8215)
