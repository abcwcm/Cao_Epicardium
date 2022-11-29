[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7374559.svg)](https://doi.org/10.5281/zenodo.7374559)      [**GEO: GSE202836**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202836)

# scRNA-seq of epicardial cells during zebrafish heart regeneration

>Xia et al.: *Activation of a transient progenitor state in the epicardium is required for zebrafish heart regeneration*. Nat.Commun.

The Applied Bioinformatics Core analyzed single-cell RNA sequencing data generated from epicardial cells of zebrafish models of cardiac injury.
The main goal was to understand whether different sub-populations of epicardial cells contribute to specific aspects of myocardial regeneration.

The raw reads were aligned and processed with the `CellRanger` pipeline (v. 3.0.2) using the zebrafish transcriptome version GRCz10.
Subsequent analyses were performed in R following the recommendations of [Amezquita et al., 2019](https://osca.bioconductor.org/), using numerous functions provided in the R packages `scater` and `scran`.

The analyses were also performed at a time when Seurat's `scTransform` method was en vogue.
For that, the [tutorials of the Satija Lab](https://satijalab.org/seurat/) were followed.
`monocle3` was used for pseudotime and trajectory analyses.

![](WCM_MB_LOGO_HZSS1L_CLR_RGB.png)

## DATA AVAILABILITY

* raw data (sequencing reads) and CellRanger output (matrices of read counts) can be downloaded from [**GEO: GSE202836**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202836)
* the `SingleCellExperiment` object named "sce_CellGeneFiltWithScTransformBatchCorrect_2019-06.rds" in the scripts here, can be downloaded [here](https://wcm.box.com/shared/static/uqv3zlp17txnb2548d7lvw43moka96yl.rds) -- caution, this is the direct link to the RDS file! Only open in R! That file is the result of the QC and filtering steps described in the [corresponding .Rmd](https://github.com/abcwcm/Cao_Epicardium/blob/main/Rmd/01_QC_and_filtering.Rmd) file above

For additional objects including the `SingleCellExperiment` objets that we generated from other publicly available scRNA-seq data, please get in touch with `abc at med.cornell.edu`, we're happy to share then with you via whatever route works best for you.

## CUSTOM PLOTTING FUNCTIONS

Most of the UMAPs were generated with internal functions for plotting, e.g. `scABC2::plot_reducedDims()`. These functions were originally based off of the plotting functions from the `scater` package, i.e. all plots can be reproduced with `scater`'s plotting functions, but we strongly recommend to rather check out the [`dittoSeq`](https://bioconductor.org/packages/release/bioc/html/dittoSeq.html) package instead. If you prefer to use our internal functions, do reach out, we're happy to share the two packages, scABC2 and ABCUtilities, but since they are unlikely to be maintained, those would be 'use at your own risk'.

## Filtering

Following CellRanger's generation of read count matrices, cells were required to have a minimum of 10e2.5 genes and max. 5% mitochondrial reads.
Genes were removed if they were detected in either fewer than 0.2% of the cells or in fewer than 5 cells per sample.

Read counts were normalized using `SCTransform` as implemented in `Seurat v.3.1` correcting for batch effects between the samples.
For visualizations and additional downstream analyses, the SCTransform-normalized (log-transformed) expression values were used unless noted otherwise.

For identifying clusters of cells with similar global transcriptomes, a shared nearest neighbor graph was constructed using `Seurat's` `FindNeighbors()` function with default settings (e.g. k = 20) and using the first 20 principal components following PCA.
Clusters were identified with `Seurat`'s `FindClusters()` function with the resolution parameter set to 0.3.
In addition, UMAP coordinates were calculated.
We assessed the identity of the cells within the resulting clusters using marker genes detected by `Seurat`'s `FindAllMarkers()` function with default settings as well as marker genes identified by `SC3` . 
Clusters with relatively low levels of *tcf21* and high levels of marker genes that indicate non-epicardial cell identities (e.g. immune cells with high levels of *cd74a* and mpeg* were removed from subsequent analyses.

After excluding the cells that belonged to the above mentioned clusters, we re-did all processing steps including (i) removal of genes that were expressed in fewer than 5 cells; (ii) calculation of normalized expression values correcting for the batch effect of the different conditions; (iii) PCA; (iv) clustering (resolution = 0.3) and (v) UMAP calculation.

## Marker gene detection

We used two complementary approaches to identify genes that may be representative of a given cell cluster: Seurat's `FindAllMarkers()` approach and the method implemented in `SC3` (Kiselev2017).
To determine the most robust candidate genes, we compared the results of both methods, using an adjusted p-value cutoff of 1\% and minimum log fold change of 0.8 for Seurat markers as well as a minimum AUROC value for SC3 markers of 0.7.

## Trajectory Inference

To infer the developmental order of certain subpopulations, particularly within the regenerating samples (d3, d7), we applied `monocle3`. 
Trajectory reconstruction methods typically need to achieve two key functions:

1. identification of the number and relationships of lineages (cells of similar transcription states), and 
2. calculation of a pseudotime measure that represents each cell's position along the trajectory in relationship to the terminal lineage.

We generally followed the instructions of the Trapnell Lab for applying the latest installment of the Monocle suite (<https://cole-trapnell-lab.github.io/monocle3/monocle3_docs/>). 
In brief, we converted the read counts into a monocle object, and re-processed the data using `preprocess_cds()` with the number of dimensions set to 100.
Cells were clustered with `cluster_cells()`, which relies on Louvain community detection similar to Seurat, but using the UMAP dimensions.
To identify the trajectories of individual cells through the UMAP space, `learn_graph()` was used.
To determine pseudotime values, root nodes were identified for each partition (as determined in the previous step) and pseudotime values were calculated based on each cell's projection on the principal graph.


## GO Term Overrepresentation

To identify functionally meaningful gene sets that are over-represented among gene lists of interest, we used the enrichment tests implemented in the functions `enrichKegg()`, `enrichPathway()` and `enrichGO()` of the R packages `clusterProfiler` and `ReactomePA`.
We tested different gene lists: (i)~marker genes determined by `Seurat` (logFC > 0 and adjusted p-value < 0.05) for each cluster, and (ii)~genes associated with temporal expression along the pseudotime trajectory determined by `Monocle~3` (q-value < 0.05 and Morans I statistic > 0.4).

Visualizations of gene set enrichments were done with the help of the `dotplot()`, `cnetplot()` and `heatplot()` functions of the `clusterProfiler()` package.

# References

* Amezquita, Robert, Lun, Aaron, Hicks, Stephanie, and Gottardo, Raphael (version 1.0.6). "Orchestrating Single-Cell Analysis with Bioconductor". <https://bioconductor.org/books/release/OSCA/>. Nat Methods. 2020 Feb;17(2):137-145. doi: 10.1038/s41592-019-0654-x. PMID: 31792435.
* **scater**: McCarthy DJ, Campbell KR, Lun ATL, Willis QF (2017). "Scater: pre-processing, quality control, normalisation and
visualisation of single-cell RNA-seq data in R." _Bioinformatics_, *33*, 1179-1186. <https://doi.org/10.1093/bioinformatics/btw777>.
* **scran**: L. Lun, A. T., Bach, K., & Marioni, J. C. (2016). Pooling across cells to normalize single-cell RNA sequencing data with many zero counts. Genome Biology, 17(1), 75. <https://doi.org/10.1186/s13059-016-0947-7>
* **scTransform**: Hafemeister, C., & Satija, R. (2019). Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. Genome Biology, 20(1), 1–15. <https://doi.org/10.1186/s13059-019-1874-1>
* **SC3** Kiselev, V. Y., Kirschner, K., Schaub, M. T., Andrews, T., Yiu, A., Chandra, T., … Hemberg, M. (2017). SC3: Consensus clustering of single-cell RNA-seq data. Nature Methods, 14(5), 483–486. <https://doi.org/10.1038/nmeth.4236Kiselev2017>
* **clusterProfiler**: Yu, G., Wang, L. G., Han, Y., & He, Q. Y. (2012). ClusterProfiler: An R package for comparing biological themes among gene clusters. OMICS A Journal of Integrative Biology, 5(16), 284–287. <https://doi.org/10.1089/omi.2011.0118>
* **monocle3**: Doesn't have its own publication yet, but see [here](https://cole-trapnell-lab.github.io/monocle3/papers/) for manuscripts related to previous versions of monocle.
