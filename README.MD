

# Blocking of the CD80/86 axis with abatacept as a therapeutic approach to prevent progression to more severe forms of COVID-19


## Abstract

We have recently provided suggestive epidemiological evidence that rheumatic patients treated with abatacept -a CTLA4-Ig protein construct that inhibits co-stimulation through the CD80/86 axis- is associated with a lower incidence of COVID-19 associated symptoms. Following this observation, we hypothesized that blocking of the CD80/86 axis could be an effective therapeutic strategy against progression to more severe proinflammatory states. To support this hypothesis, we analyzed transcriptional data associated with the response to abatacept with pathological features associated with COVID-19, particularly in its more severe forms. We have found that most of the biological processes that have been consistently associated with COVID-19 pathology are reverted by treatment with abatacept, including downregulation of IL6 production, macrophage activation and the complement cascade. Comparing also to transcriptional data from SARS-CoV2 infected patients, we have found that the transcriptional response to abatacept has a very high level of antagonism to that elicited by COVID-19. Single cell RNA-Seq data from bronchoalveolar lavage fluid cells from COVID-19 patients and controls, showed a significant correlation along the main elements of the C80/86 axis: CD86+/80+ antigen presenting cells, activated CD4+(CTLA4+) T cells and IL6 production. Our in-silico study provides additional support to the hypothesis that blocking of the CD80/CD86 signaling axis may be protective of the excessive proinflammatory state associated with COVID-19 severe outcome. 

## Important Links

The manuscript can be found here **(link biorxiv)**

Raw RNA-sequencing data from the abatacept-treated patients will be deposited in the GEO database.

RNA-Sequencing pipeline for PBMCs from COVID-19 patients is available at the original authors repository: https://github.com/zhouyulab/ncov/ 

scRNA-Sequencing data from BALF samples from COVID-19 patients and healthy controls was downloadded from GEO accession number [GSE145926](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145926)


## Files

The main directory contains the scripts, input and output folders for the different analysis in the manuscript: 

* BALF_COVID19_scRNASeq_analysis.R: data processing from GEO download of the single cell RNA-seq data of the 12 BALF samples to figure generation. 
* DEG_GSEA_COVID19.R: differential expression analysis and GSEA analysis of the COVID-19 PBMC data.
* DEG_GSEA_RA_Abatacept.R: differential expression and GSEA analysis of the abatacept RNA-Seq longitudinal dataset.
* DEG_GSEA_RA_Abatacept_GranulocytesAdj.R: differential expression analysis and GSEA analysis of the abatacept RNA-seq longitudinal dataset, with adjustment for granulocyte percentage.
* Covid_BPs_in_Abatacept.R: analysis of the biological processes consistently associated with COVID-19 pathology and evaluation of antagonism by abatacept. 
* COVID_Abatacept_transcriptome_comparison.R: comparison of the BPs modified by COVID-19 and abatacept exposure and identification of antagonism.

## Session info
The analyses comparing abatacept transcriptome to COVID-19 curated processes and SARS-CoV-2 induced transcriptome were run under R version 3.6.3 (2020-02-29), platform x86_64-w64-mingw32/x64 (64-bit), Windows 10 x64 (build 18363).

The following packages were used:
shiny_1.4.0.2, edgeR_3.28.1, limma_3.42.2, fgsea_1.12.0, Rcpp_1.0.3, gridExtra_2.3, ggplot2_3.3.0, prabclus_2.3-2, mclust_5.4.6, MASS_7.3-51.6, eulerr_6.1.0, corrplot_0.84, metap_1.3, GO.db_3.10.0, org.Hs.eg.db_3.10.0, AnnotationDbi_1.48.0, IRanges_2.20.2, S4Vectors_0.24.4, Biobase_2.46.0, BiocGenerics_0.32.0, DeconCell_0.1.0, Seurat_3.1.5, sctransform_0.2.1, hdf5r_1.3.2 



