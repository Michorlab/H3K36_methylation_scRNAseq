H3K36 methylation maintains cell identity by regulating opposing lineage programs
--------

This repository contains code, data, tables and plots to support single cell RNA sequencing analyses and reproduce results from the paper [H3K36 methylation maintains cell identity by regulating opposing lineage programs](https://), currently under revision as of February 2023.


Abstract
--------
While epigenetic regulation in early development and differentiation is relatively well-characterized, the precise mechanisms that subsequently maintain specialized cell states remain largely unexplored. Here, we employed histone H3.3 mutants to uncover a crucial role for H3K36 methylation in the maintenance of cell identities across diverse developmental paradigms including induced pluripotency, somatic lineage conversion and stem cell differentiation. Focusing on the experimentally induced conversion of fibroblasts to pluripotent stem cells, we show that H3K36M-mediated disruption of H3K36 methylation endows reprogramming intermediates with a plastic state poised to acquire pluripotency in nearly all cells. At a cellular level, H3K36M renders mesenchymal cells insensitive to TGFB signals and thus facilitates epithelial plasticity. At a molecular level, H3K36M leads to the downregulation of mesenchymal genes by depleting H3K36me2 and H3K27ac at associated enhancers. In parallel, H3K36M leads to the upregulation of epithelial and stem cell genes by enabling enhancer accessibility and hypomethylation. Mechanistically, we found that Tet-dependent demethylation uncouples pluripotency enhancer activation from somatic enhancer decommissioning in K36M cells. This dynamic switch of enhancer activities redirects Sox2 binding from promiscuous somatic targets to bona fide stem cell targets crucial for the establishment of the pluripotency network. Together, our findings reveal a dual role for H3K36 methylation in the maintenance of cell identity by integrating a key developmental pathway into sustained expression of cell type-specific programs, and by opposing the expression of alternative lineage programs through continual enhancer methylation. Our results provide crucial insight into the impact of dynamic H3K36 methylation patterns on physiological and pathological cell fate transitions, including development, tissue regeneration and cancer. 


Content
-------
* `/code/`: R scripts to reproduce analyses of single cell RNA sequencing data
* `/data/`: path folder to read in cellranger outputs, as well as preprocessed RDS data objects to be used as input for the scripts in this repo. This data needs to be downloaded from [Zenodo](https://zenodo.org/record/7454907#.Y97iXS-B2NF), and placed in the folder 'data'
* `/plots/`: output plots from the single cell RNA sequencing analyses, organized by self-explanatory folder names.
* `/tables/`: relevant gene lists (e.g. signatures) used throughout the analyses, together with coordinates of cells after dimensionality reduction


Data
-------
The scRNA-seq data generated in this study have been deposited in the NCBI GEO database under accession number [GSE203536](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203536).

Preprocessed RDS objects can be downloaded from [zenodo](https://zenodo.org/record/7454907#.Y97iXS-B2NF).


## Contact
For any questions, comments or suggestions regarding the analysis or interpretation of the single cell sequencing data in this manuscript, please feel to reach out to me via [GitHub](https://github.com/csimona), [Twitter](https://twitter.com/simocristea) or email (scristea@jimmy.harvard.edu).