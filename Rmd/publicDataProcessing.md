# Revision August 2022

The reviewers asked to compare our findings to published data sets.
To this end, we reached out to multiple PIs/authors to obtain their processed data sets.

* [Hesse: Mouse EpiSC and aSCS](#hesse)
* [Kapuria: ZF unij. vs. regen. heart](#kapuria)
* [Sun: ZF unij. vs. regen. heart](#sun)

--------------------------------------

## DeBakkers et al.

[Prrx1b restricts fibrosis and promotes Nrg1-dependent cardiomyocyte proliferation during zebrafish heart regeneration](http://dx.doi.org/10.1242/DEV.198937)

Followed DeBakkers' scripts to generate an SCE object: `DeBakkers_sce.RDS`


<a name="hesse" />
## Hesse et al (2021): Mouse injured and regenerating

Data shared by Julia Hesse (Uni Duesseldorf).

>processed scRNAseq data of Figure 1 as R Object: `EPDC.SCTransform.integrated.rds` (renamed to `Hesse2021_EPDC.SCTransform.integrated.rds`)

>Please note that the naming of the clustering does not match that in the paper, as the cluster names were adjusted after data analysis, but you can easily assign labels using Figure 1 as a template.

```
library(Seurat); library(SingleCellExperiment)
seur <- readRDS("Hesse2021_EPDC.SCTransform.integrated.rds")
#An object of class Seurat
#38488 features across 13796 samples within 3 assays
#Active assay: SCT (17394 features, 0 variable features)
# 2 other assays present: RNA, integrated
# 2 dimensional reductions calculated: pca, umap
sce <- as.SingleCellExperiment(seur)
sce
#class: SingleCellExperiment 
#dim: 17394 13796 
#metadata(0):
#assays(2): counts logcounts
#rownames(17394): Xkr4 Mrpl15 ... Scgb3a2 Sall3
#rowData names(0):
#colnames(13796): AAACCTGAGACAGAGA-1 AAACCTGAGTGCGATG-1 ...
#  TTTGGTTCATGTTGAC-3 TTTGTCACACCAGATT-3
#colData names(29): orig.ident nCount_RNA ... migration_markers1 ident
#reducedDimNames(2): PCA UMAP
#altExpNames(0):

saveRDS(sce, file = "Hesse2021_EPDC.SCTransform.integrated_sce.rds")

```

<a name="sun" />
## Sun et al (2022): Zebrafish, uninjured vs regenerating heart (dpa7)

*hapln1 Defines an Epicardial Cell Subpopulation Required for Cardiomyocyte Expansion During Heart Morphogenesis and Regeneration*. Circulation (2022). 
Jisheng Sun, Elizabeth A. Peterson, Annabel Z. Wang, Jianhong Ou, Kieko E. Smith, Kenneth D. Poss and Jinhu Wang

<https://www.ahajournals.org/doi/full/10.1161/CIRCULATIONAHA.121.055468>

![](https://www.ahajournals.org/cms/asset/c4026df2-2c0a-42e7-bd31-fa75d94d5f9b/circulationaha.121.055468.fig01.jpg)

- Seurat object: `Sun2022_Epicard_uninured_and_injured.rds`

```
library(Seurat); library(SingleCellExperiment)
seur <- readRDS("Sun2022_Epicard_uninured_and_injured.rds")
#An object of class Seurat 
#21397 features across 7869 samples within 2 assays 
#Active assay: RNA (19397 features, 0 variable features)
# 1 other assay present: integrated
# 2 dimensional reductions calculated: pca, umap


sce <- as.SingleCellExperiment(seur)
sce
#class: SingleCellExperiment 
#dim: 19397 7869 
#metadata(0):
#assays(2): counts logcounts
#rownames(19397): ptpn12 phtf2.1 ... NC-002333.14 mCherry
#rowData names(0):
#colnames(7869): AAACCCAAGGTAGGCT-1_1 AAACCCAGTTGCCATA-1_1 ...
#  TTTGTTGCATTGCTTT-1_2 TTTGTTGTCCTACGGG-1_2
#colData names(7): nCount_RNA nFeature_RNA ... seurat_clusters ident
#reducedDimNames(2): PCA UMAP
#altExpNames(0):

saveRDS(sce, "Sun2022_Epicard_uninured_and_injured_sce.rds")

```

---------------------------------------

<a name="kapuria" />
## Kapuria et al (2022): Zebrafish, uninjured vs regenerating heart (dpa7)

* [publication: "Heterogeneous pdgfrb+ cells regulate coronary vessel development and revascularization during heart regeneration"](http://dx.doi.org/10.1242/dev.199752) 
* obtained Seurat RDS file from Kapuria 
    - the .rds file used in the paper for uninj_inj heart
    -`DefaultAssay(R-object name) <- "RNA"`, for comparing gene expression across conditions, cell-clusters
    - For cell labeling please follow the supplementary figure 5.1 from the paper.

Fig legend SupFig 5.1 refers to the MUT/WT comparison, I think he meant Fig S9:

**Characterization of *pdgfrb:EGFP* FACS sorted cells in uninjured and 7 dpa hearts.**
(A-A’) *pdgfrb:EGFP* cells FACS sorted from *Tg(pdgfrb:EGFP; cxcl12b:Citrine)* fish.

(A) UMAP plot of all FACS sorted cells (the uninjured and the 7 days post amputated samples integrated together), which form 15 clusters.

(A’) Relative abundance of uninjured and 7 dpa hearts in each cluster.

(B) Dot plot of the 3 top differentially expressed cluster marker genes of FACS sorted *pdgfrb:EGFP* expressing cells. Blue dot, uninjured. Red dot, 7 dpa. Expression level across cells within the cluster is shown in intensity of the color while the percentage of the cells expressing the marker gene is shown by the size of the dot (0- 100%). Differentially expressed genes were determined with minimum percent expression cut-off = 0.1 and minimum average log fold change = 0.25.

(C) Dot plot of differentially expressed cluster marker genes of smooth muscle, mural cells, epicardial cells and fibroblasts. Pdgfrb mRNA is specifically expressed in the cluster 5, 6, 10 and 12 of which, cluster 5, 6, and 10 express the epicardial markers (tcf21, tbx18, wt1a) and cluster 12 express the smooth muscle cell markers (tagln, acta2, myh11a) and mural cell markers (rgs5a, cd248a, kcne4, ndufa4l2a). Fibroblast marker postnb predominantly express in cluster 10 in uninjured heart but upregulated after injury in all pdgfrb + clusters (cluster 5, 6, 10, 12). Fibroblast marker pdgfra, specifically express in cluster 5, 6. Other fibroblast markers (lum, col1a2, col5a1, vim) express in cluster 5, 6 and 10, where cluster 6 has most fibroblast marker expression. lum express rarely but induced after injury in cluster 5, 6. col1a2, col5a1 induced after injury in cluster 5, 10 and 12.

(D) Cell identity for all clusters determined by the cluster marker genes. Cluster 1, 2, 3, 4 cells are identified as cardiomyocytes. cluster 5, 6, 10 cells are epicardial/EPDC/fibroblasts. Cluster 12 cells are mural cells. Cluster 0, 11, 13 are endothelial/endocardial cells. Cluster 8 cells are red blood cells (RBC). Cluster 9, 14 are macrophage/lymphocytes and cluster 7 cells are T cells/lymphocytes.

| Cluster | label |
|----------|------|
| 1, 2, 3, 4 |cardiomyocytes |
| 5, 6, 10 | epicardial/EPDC/fibroblasts |
| 12 | mural cells |
| 0, 11, 13 | endothelial/endocardial |
| 8 | red blood cells (RBC) |
| 9, 14 | macrophage/lymphocytes |
| 7  | T cells/lymphocytes |


```
> library(Seurat); library(SingleCellExperiment)
> seur <- readRDS("Kapuria2022_integrated_cells_uninj_7dpa_hrt.rds")
> seur
An object of class Seurat 
17216 features across 2728 samples within 2 assays 
Active assay: integrated (2000 features, 2000 variable features)
 1 other assay present: RNA
  3 dimensional reductions calculated: pca, umap, tsne

> DefaultAssay(seur) <- "RNA"

> head(seur@meta.data)
                           orig.ident nCount_RNA nFeature_RNA percent.mt
AAACCCACATGACAGG-1_1 pdgfrb_cell_ctrl      12712         1476 0.76305853
AAACCCATCATGCCAA-1_1 pdgfrb_cell_ctrl        770          258 2.98701299
AAACGAAGTCCCTGAG-1_1 pdgfrb_cell_ctrl       1945          407 0.05141388
AAACGAATCCTAGAGT-1_1 pdgfrb_cell_ctrl       5274         1343 1.49791430
AAAGAACAGACCATTC-1_1 pdgfrb_cell_ctrl       1066          271 3.37711069
AAAGAACGTAGATTGA-1_1 pdgfrb_cell_ctrl       1852          862 2.05183585
                     RNA_snn_res.0.5 seurat_clusters integrated_snn_res.0.5
AAACCCACATGACAGG-1_1               2               3                      3
AAACCCATCATGCCAA-1_1               3               4                      4
AAACGAAGTCCCTGAG-1_1               1               4                      4
AAACGAATCCTAGAGT-1_1               0               0                      0
AAAGAACAGACCATTC-1_1               3               4                      4
AAAGAACGTAGATTGA-1_1               0              13                     13

> unique(seur@meta.data$orig.ident)
[1] "pdgfrb_cell_ctrl" "pdgfrb_cell_7dpa"

> sce <- as.SingleCellExperiment(seur)
> sce
class: SingleCellExperiment 
dim: 15216 2728 
metadata(0):
assays(2): counts logcounts
rownames(15216): phtf2.1 dusp16 ... zgc:173552.7 igl4v8
rowData names(0):
colnames(2728): AAACCCACATGACAGG-1_1 AAACCCATCATGCCAA-1_1 ...
  TTTGGTTCATACAGAA-1_2 TTTGGTTGTCCTTAAG-1_2
colData names(8): orig.ident nCount_RNA ... integrated_snn_res.0.5
  ident
reducedDimNames(3): PCA UMAP TSNE
altExpNames(0):

> sce$ident <- as.character(sce$ident)
> sce$label <- ifelse(sce$ident %in% c("1","2","3","4"), "CM",
    ifelse(sce$ident %in% c("5","6","10"), "EPDC_fibroblasts",
        ifelse(sce$ident %in% c("12"), "mural",
            ifelse(sce$ident %in% c("0","11","13"), "endothelial_endocardial",
                ifelse(sce$ident == "8", "RBC",
                    ifelse(sce$ident %in% c("9","14"), "M0_lympho",
                        ifelse(sce$ident == "7", "T_lympho", "unknown")))))))

> saveRDS(sce, file = "Kapuria2022_integrated_cells_uninj_7dpa_hrt_sce.rds")
```
