# Cao_Epicardium
Scripts related to scRNA-seq data analysis 

The Applied Bioinformatics Core analyzed single-cell RNA sequencing data generated from epicardial cells of zebrafish models of cardiac injury.
The main goal was to understand whether different sub-populations of epicardial cells contribute to specific aspects of myocardial regeneration.

The raw reads were aligned and processed with the CellRanger pipeline (v. 3.0.2) using the zebrafish transcriptome version GRCz10.
Subsequent analyses were performed in R following the recommendations of \citet{Amezquita2019} (\url{https://osca.bioconductor.org/}) using numerous functions provided in the R packages \texttt{scater} and \texttt{scran} \citep{McCarthy2017,scran} as well as \texttt{Seurat} following the tutorials of the Satija Lab (\url{https://satijalab.org/seurat/}).

Cells were required to have a minimum of 10e2.5 genes and max. 5% mitochondrial reads
Genes were removed if they were detected in either fewer than 0.2% of the cells or in fewer than 5cells per sample.

Read counts were ormalized using `SCTransform` as implemented in `Seurat~v.3.1` correcting for batch effects between the samples.
For visualizations and additional downstream analyses, the SCTransform-normalized (log-transformed) expression values were used unless noted otherwise.

For identifying clusters of cells with similar global transcriptomes, a shared nearest neighbor graph was constructed using `Seurat's` `FindNeighbors()` function with default settings (e.g. k = 20) and using the first 20 principal components following PCA.
Clusters were identified with `Seurat`'s `FindClusters()` function with the resolution parameter set to 0.3 .
In addition, UMAP coordinates were calculated.
We assessed the identity of the cells within the resulting clusters using marker genes detected by `Seurat`'s `FindAllMarkers()` function with default settings as well as marker genes identified by `SC3` . 
Clusters with relatively low levels of *tcf21* and high levels of marker genes that indicate non-epicardial cell identities (e.g. immune cells with high levels of *cd74a* and mpeg* were removed from subsequent analyses.

After excluding the cells that belonged to the above mentioned clusters, we re-did all processing steps including (i) removal of genes that were expressed in fewer than 5 cells; (ii) calculation of normalized expression values correcting for the batch effect of the different conditions; (iii) PCA; (iv) clustering (resolution = 0.3) and (v) UMAP calculation.

## Marker gene detection

We used two complementary approaches to identify genes that may be representative of a given cell cluster: Seurat's \texttt{FindAllMarkers} approach and the method implemented in \texttt{SC3} \citep{Kiselev2017}.
To determine the most robust candidate genes, we compared the results of both methods, using an adjusted p-value cutoff of 1\% and minimum log fold change of 0.8 for Seurat markers as well as a minimum AUROC value for SC3 markers of 0.7.

## Trajectory Inference

To infer the developmental order of certain subpopulations, particularly within the regenerating samples (d3, d7), we applied two different trajectory reconstruction algorithms: \texttt{Monocle~3} and \texttt{slingshot} \citep{monocle3,slingshot}.
Trajectory reconstruction methods typically need to achieve two key functions: (i)~identification of the number and relationships of lineages (cells of similar transcription states), and (ii)~calculation of a pseudotime measure that represents each cell's position along the trajectory in relationship to the terminal lineage.

We generally followed the instructions of the Trapnell Lab for applying the latest installment of the Monocle suite (\url{https://cole-trapnell-lab.github.io/monocle3/monocle3_docs/}).
In brief, we converted the read counts into a monocle object, and re-processed the data using \texttt{preprocess\_cds()} with the number of dimensions set to 100.
Cells were clustered with \texttt{cluster\_cells()}, which relies on Louvain community detection similar to Seurat, but using the UMAP dimensions.
To identify the trajectories of individual cells through the UMAP space, \texttt{learn\_graph()} was used.
To determine pseudotime values, root nodes were identified for each partition (as determined in the previous step) and pseudotime values were calculated based on each cell's projection on the principal graph.

## GO Term Overrepresentation

To identify functionally meaningful gene sets that are over-represented among gene lists of interest, we used the enrichment tests implemented in the functions \texttt{enrichKegg}, \texttt{enrichPathway} and \texttt{enrichGO} of the R packages \texttt{clusterProfiler} and \texttt{ReactomePA} \citep{clusterprofiler,reactomePA}.
We tested different gene lists: (i)~marker genes determined by \texttt{Seurat} (logFC \textgreater 0 and adjusted p-value \textless 0.05) for each cluster, and (ii)~genes associated with temporal expression along the pseudotime trajectory determined by \texttt{Monocle~3} (q-value \textless 0.05 and Morans I statistic \textgreater 0.4).

Visualizations of gene set enrichments were done with the help of the \texttt{dotplot}, \texttt{cnetplot} and \texttt{heatplot} functions of the \texttt{clusterProfiler} package \citep{clusterprofiler}.

# References

* Hafemeister2019
* Kiselev2017
* Waltman2013
* UMAP2019
