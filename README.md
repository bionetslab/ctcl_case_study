# CTCL case study for the SHouT python package

## Overview of analyses

![ctcl_shout_overview](/readme_images/fig0_summary.jpeg)
*Fig: Overview of our analyses. (A) We generated multi-antigen images for 69 skin tissue samples from 21 CTCL, 23 AD, and 25 PSO patients. Subsequently, images were pre-processed via cell segmentation, cell-level protein abundance quantification, and cell type assignment. (B) We then computed spatial graph representations for all samples, which we analyzed using SquidPy as well as different heterogeneity scores implemented in our Python package SHouT: local and global entropy, local and global homophily, and egophily.*


## Cell type-wise SHouT scores

### Local entropy

![local_entropy_r=5](/readme_images/local_entropy_r_5.png)
*Fig: Distributions of local entropy with radius r = 5, across all samples, annotated with Bonferonni-corrected MWU P-values per condition pair.*


### Local homophily

![local_homophily_r=5](/readme_images/local_homophily_r_5.png)
*Fig: Distributions of local homophily with radius r = 5, across all samples, annotated with Bonferonni-corrected MWU P-values per condition pair.*


### Local Egophily

![egophily_r=5](/readme_images/egophily_r_5.png)
*Fig: Distributions of egophily with radius r = 5, across all samples, annotated with Bonferonni-corrected MWU P-values per condition pair.*


## Citing the work

Please cite the paper as follows:
- Will be updated once available.