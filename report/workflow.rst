This workflow performs differential expression analysis on single-end RNA-seq data of the diploid strains (NE/CA) and wild type strain (Sc/Sc).

Materials and methods
---------------------

We follow the data analysis procedure performed in Ghosh et al. (doi: 10.1016/j.compbiolchem.2020.107239). Version of tools and genomes have been updated and parameter used have been adapted according to our data.

The quality of the raw reads were assessed using FastQC v0.11.9 toolkit (Andrews, 2010). Adapters and low quality reads were trimmed using Trimmomatic v0.39 (Bolger et al., 2014).

HiSat2 v2.2.1 (Kim et al., 2015) was then used for mapping the raw reads to the reference genome Saccharomyces cerevisiae S288C (assembly R64). The expression for each gene were evaluated using featureCounts from the Subread v2.0.1 package (Liao et al., 2019) based on annotation of Saccharomyces cerevisiae assembly R64.

The low expressed genes which did not have more than one count per million reads (1CPM) in at least three samples within each datasets were removed from further analysis. The raw counts were then normalized and used for differential expression testing using DESeq2 v1.28.0 (Love et al., 2014). 

The genes with fold change > 2 and adjusted p-value (padj) < 0.01 were considered significantly upregulated and overexpressed while those with fold change < -2 and padj < 0.01 were considered significantly downregulated or underexpressed.
