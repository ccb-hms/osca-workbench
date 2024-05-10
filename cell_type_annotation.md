---
title: Cell type annotation
teaching: 30 # Minutes of teaching in the lesson
exercises: 15 # Minutes of exercises in the lesson
editor_options: 
  markdown: 
    wrap: 72
---

::: questions
-   ScRNA-seq: Can clustering group cells by type?
-   Marker genes: Key to identifying cell types in scRNA-seq clusters?
-   What's tricky about labeling cell types in scRNA-seq data?
-   How to confirm and improve cell type labels from scRNA-seq?
:::

::: objectives
-   Identify cell types in scRNA-seq data by clustering cells based on gene expression patterns.
-   Define cell types within scRNA-seq clusters using marker genes.
-   Compare methods for assigning cell type labels in scRNA-seq data: reference datasets vs. marker genes.
-   Evaluate scRNA-seq cell type labels by analyzing marker gene distribution within clusters.
:::

## Setup


```r
library(AUCell)
library(BiocStyle)
library(MouseGastrulationData)
library(SingleR)
library(bluster)
library(scater)
library(scran)
```

## Data retrieval


```r
sce <- WTChimeraData(samples = 5, type = "processed")
sce
```

```{.output}
class: SingleCellExperiment 
dim: 29453 2411 
metadata(0):
assays(1): counts
rownames(29453): ENSMUSG00000051951 ENSMUSG00000089699 ...
  ENSMUSG00000095742 tomato-td
rowData names(2): ENSEMBL SYMBOL
colnames(2411): cell_9769 cell_9770 ... cell_12178 cell_12179
colData names(11): cell barcode ... doub.density sizeFactor
reducedDimNames(2): pca.corrected.E7.5 pca.corrected.E8.5
mainExpName: NULL
altExpNames(0):
```

To speed up the computations, we subsample the dataset to 1,000 cells.


```r
set.seed(123)
ind <- sample(ncol(sce), 1000)
sce <- sce[,ind]
```

## Preprocessing


```r
sce <- logNormCounts(sce)
sce <- runPCA(sce)
```

## Clustering

Clustering is an unsupervised learning procedure that is used to
empirically define groups of cells with similar expression profiles. Its
primary purpose is to summarize complex scRNA-seq data into a digestible
format for human interpretation. This allows us to describe population
heterogeneity in terms of discrete labels that are easily understood,
rather than attempting to comprehend the high-dimensional manifold on
which the cells truly reside. After annotation based on marker genes,
the clusters can be treated as proxies for more abstract biological
concepts such as cell types or states.

Popularized by its use in
[Seurat](https://cran.r-project.org/web/packages/Seurat/index.html),
graph-based clustering is a flexible and scalable technique for
clustering large scRNA-seq datasets. We first build a graph where each
node is a cell that is connected to its nearest neighbors in the
high-dimensional space. Edges are weighted based on the similarity
between the cells involved, with higher weight given to cells that are
more closely related. We then apply algorithms to identify "communities"
of cells that are more connected to cells in the same community than
they are to cells of different communities. Each community represents a
cluster that we can use for downstream interpretation.

Here, we use the `clusterCells()` function from the
[scran](https://bioconductor.org/packages/scran) package to perform
graph-based clustering using the [Louvain
algorithm](https://doi.org/10.1088/1742-5468/2008/10/P10008) for
community detection. All calculations are performed using the top PCs to
take advantage of data compression and denoising. This function returns
a vector containing cluster assignments for each cell in our
`SingleCellExperiment` object.


```r
colLabels(sce) <- clusterCells(sce, use.dimred = "PCA",
                               BLUSPARAM = NNGraphParam(cluster.fun = "louvain"))
table(colLabels(sce))
```

```{.output}

  1   2   3   4   5   6   7   8   9  10  11 
100 160  99 141  63  93  60 108  44  91  41 
```

We assign the cluster assignments back into our `SingleCellExperiment`
object as a `factor` in the column metadata. This allows us to
conveniently visualize the distribution of clusters in eg. a *t*-SNE or
a UMAP.


```r
sce <- runUMAP(sce, dimred = "PCA")
plotReducedDim(sce, "UMAP", color_by = "label")
```

<img src="fig/cell_type_annotation-rendered-cluster-viz-1.png" style="display: block; margin: auto;" />

## Marker gene detection

To interpret clustering results as obtained in the previous section, we
identify the genes that drive separation between clusters. These marker
genes allow us to assign biological meaning to each cluster based on
their functional annotation. In the simplest case, we have *a priori*
knowledge of the marker genes associated with particular cell types,
allowing us to treat the clustering as a proxy for cell type identity.

The most straightforward approach to marker gene detection involves
testing for differential expression between clusters. If a gene is
strongly DE between clusters, it is likely to have driven the separation
of cells in the clustering algorithm.

Here, we perform a Wilcoxon rank sum test against a log2 fold change
threshold of 1, focusing on up-regulated (positive) markers in one
cluster when compared to another cluster.


```r
rownames(sce) <- rowData(sce)$SYMBOL
markers <- findMarkers(sce, test.type = "wilcox", direction = "up", lfc = 1)
markers
```

```{.output}
List of length 11
names(11): 1 2 3 4 5 6 7 8 9 10 11
```

The resulting object contains a sorted marker gene list for each
cluster, in which the top genes are those that contribute the most to
the separation of that cluster from mall other clusters.

Here, we inspect the ranked marker gene list for the first cluster.


```r
markers[[1]]
```

```{.output}
DataFrame with 29453 rows and 14 columns
                 Top     p.value         FDR summary.AUC     AUC.2     AUC.3
           <integer>   <numeric>   <numeric>   <numeric> <numeric> <numeric>
Crabp2             1 7.31206e-35 1.95784e-31    0.938625  0.938625 0.9159596
Ptn                1 2.26190e-43 6.66197e-39    0.983313  0.983313 0.9812121
Crabp1             1 1.13915e-32 2.09696e-29    0.926687  0.926687 0.7927273
Zeb2               2 3.98744e-10 1.38167e-07    0.801500  0.553125 0.3943434
Mest               2 1.10835e-24 9.60122e-22    0.883000  0.883000 0.0249495
...              ...         ...         ...         ...       ...       ...
AC125149.2     29448           1           1           0         0         0
AC125149.4     29449           1           1           0         0         0
AC234645.1     29450           1           1           0         0         0
AC168977.2     29451           1           1           0         0         0
Vmn2r122       29453           1           1           0         0         0
               AUC.4     AUC.5     AUC.6     AUC.7     AUC.8     AUC.9
           <numeric> <numeric> <numeric> <numeric> <numeric> <numeric>
Crabp2      0.861986  0.494444  0.690000  0.763833  0.845370  0.769091
Ptn         0.911773  0.771429  0.198710  0.830667  0.563333  0.898409
Crabp1      0.540071  0.714127  0.695591  0.767167  0.620556  0.524773
Zeb2        0.604184  0.390317  0.337849  0.801500  0.447963  0.341364
Mest        0.136028  0.476190  0.091828  0.386833  0.314444  0.773864
...              ...       ...       ...       ...       ...       ...
AC125149.2         0         0         0         0         0         0
AC125149.4         0         0         0         0         0         0
AC234645.1         0         0         0         0         0         0
AC168977.2         0         0         0         0         0         0
Vmn2r122           0         0         0         0         0         0
              AUC.10    AUC.11
           <numeric> <numeric>
Crabp2     0.9416484  0.882439
Ptn        0.9463736  0.967317
Crabp1     0.8584615  0.714878
Zeb2       0.6175824  0.501463
Mest       0.0154945  0.290976
...              ...       ...
AC125149.2         0         0
AC125149.4         0         0
AC234645.1         0         0
AC168977.2         0         0
Vmn2r122           0         0
```

The `Top` field provides the the minimum rank across all pairwise
comparisons. The `p.value` field provides the combined *p*-value across
all comparisons, and the `FDR` field the BH-adjusted *p*-value for each
gene. The `summary.AUC` provides area under the curve (here the
concordance probability) from the comparison with the lowest *p*-value,
the `AUC.n` fields provide the AUC for each pairwise comparison. The AUC
is the probability that a randomly selected cell in cluster *A* has a
greater expression of gene *X* than a randomly selected cell in *B*.

We can then inspect the top marker genes for the first cluster using the
`plotExpression` function from the
[scater](https://bioconductor.org/packages/scater) package.


```r
top.markers <- head(rownames(markers[[1]]))
plotExpression(sce, features = top.markers, x = "label", color_by = "label")
```

<img src="fig/cell_type_annotation-rendered-plot-markers-1.png" style="display: block; margin: auto;" />

## Cell type annotation

The most challenging task in scRNA-seq data analysis is arguably the
interpretation of the results. Obtaining clusters of cells is fairly
straightforward, but it is more difficult to determine what biological
state is represented by each of those clusters. Doing so requires us to
bridge the gap between the current dataset and prior biological
knowledge, and the latter is not always available in a consistent and
quantitative manner. Indeed, even the concept of a "cell type" is [not
clearly defined](https://doi.org/10.1016/j.cels.2017.03.006), with most
practitioners possessing a "I'll know it when I see it" intuition that
is not amenable to computational analysis. As such, interpretation of
scRNA-seq data is often manual and a common bottleneck in the analysis
workflow.

To expedite this step, we can use various computational approaches that
exploit prior information to assign meaning to an uncharacterized
scRNA-seq dataset. The most obvious sources of prior information are the
curated gene sets associated with particular biological processes, e.g.,
from the Gene Ontology (GO) or the Kyoto Encyclopedia of Genes and
Genomes (KEGG) collections. Alternatively, we can directly compare our
expression profiles to published reference datasets where each sample or
cell has already been annotated with its putative biological state by
domain experts. Here, we will demonstrate both approaches on the
wild-type chimera dataset.

### Assigning cell labels from reference data

A conceptually straightforward annotation approach is to compare the
single-cell expression profiles with previously annotated reference
datasets. Labels can then be assigned to each cell in our
uncharacterized test dataset based on the most similar reference
sample(s), for some definition of "similar". This is a standard
classification challenge that can be tackled by standard machine
learning techniques such as random forests and support vector machines.
Any published and labelled RNA-seq dataset (bulk or single-cell) can be
used as a reference, though its reliability depends greatly on the
expertise of the original authors who assigned the labels in the first
place.

In this section, we will demonstrate the use of the
*[SingleR](https://bioconductor.org/packages/3.18/SingleR)* method for cell type annotation [Aran et al.,
2019](https://www.nature.com/articles/s41590-018-0276-y). This method
assigns labels to cells based on the reference samples with the highest
Spearman rank correlations, using only the marker genes between pairs of
labels to focus on the relevant differences between cell types. It also
performs a fine-tuning step for each cell where the correlations are
recomputed with just the marker genes for the top-scoring labels. This
aims to resolve any ambiguity between those labels by removing noise
from irrelevant markers for other labels. Further details can be found
in the [*SingleR*
book](https://bioconductor.org/books/release/SingleRBook) from which
most of the examples here are derived.


```r
ref <- EmbryoAtlasData(samples = 29)
ref
```

```{.output}
class: SingleCellExperiment 
dim: 29452 7569 
metadata(0):
assays(1): counts
rownames(29452): ENSMUSG00000051951 ENSMUSG00000089699 ...
  ENSMUSG00000096730 ENSMUSG00000095742
rowData names(2): ENSEMBL SYMBOL
colnames(7569): cell_95727 cell_95728 ... cell_103294 cell_103295
colData names(17): cell barcode ... colour sizeFactor
reducedDimNames(2): pca.corrected umap
mainExpName: NULL
altExpNames(0):
```

```r
table(ref$stage)
```

```{.output}

E8.5 
7569 
```

To speed up the computations, we subsample the dataset to 1,000 cells.


```r
set.seed(123)
ind <- sample(ncol(ref), 1000)
ref <- ref[,ind]
```


```r
tab <- sort(table(ref$celltype), decreasing = TRUE)
tab
```

```{.output}

  Forebrain/Midbrain/Hindbrain                     Erythroid3 
                           131                             75 
             Paraxial mesoderm                            NMP 
                            69                             51 
                  ExE mesoderm               Surface ectoderm 
                            49                             47 
                     Allantois                     Mesenchyme 
                            46                             45 
                   Spinal cord            Pharyngeal mesoderm 
                            45                             41 
                  ExE endoderm                   Neural crest 
                            38                             35 
                           Gut Haematoendothelial progenitors 
                            30                             27 
         Intermediate mesoderm                 Cardiomyocytes 
                            27                             26 
              Somitic mesoderm                    Endothelium 
                            25                             23 
                    Erythroid2                  Def. endoderm 
                            11                              3 
                    Erythroid1            Blood progenitors 1 
                             2                              1 
           Blood progenitors 2                Caudal Mesoderm 
                             1                              1 
                           PGC 
                             1 
```


```r
ref <- logNormCounts(ref)
```

Some cleaning - remove cells of the reference dataset for which the cell
type annotation is missing.


```r
nna <- !is.na(ref$celltype)
ref <- ref[,nna]
```

Also remove cell types of very low abundance (here less than 10 cells)
to remove noise prior to subsequent annotation tasks.


```r
abu.ct <- names(tab)[tab >= 10]
ind <- ref$celltype %in% abu.ct
ref <- ref[,ind] 
```

Restrict to genes shared between query and reference dataset.


```r
rownames(ref) <- rowData(ref)$SYMBOL
isect <- intersect(rownames(sce), rownames(ref))
sce <- sce[isect,]
ref <- ref[isect,]
```

Convert sparse assay matrices to regular dense matrices for input to
SingleR.


```r
sce.mat <- as.matrix(assay(sce, "logcounts"))
ref.mat <- as.matrix(assay(ref, "logcounts"))
```


```r
res <- SingleR(test = sce.mat, ref = ref.mat, labels = ref$celltype)
res
```

```{.output}
DataFrame with 1000 rows and 4 columns
                                   scores                 labels delta.next
                                 <matrix>            <character>  <numeric>
cell_11995 0.348586:0.335451:0.314515:... Forebrain/Midbrain/H..  0.1285110
cell_10294 0.273570:0.260013:0.298932:...             Erythroid3  0.1381951
cell_9963  0.328538:0.291288:0.475611:...            Endothelium  0.2193295
cell_11610 0.281161:0.269245:0.299961:...             Erythroid3  0.0359215
cell_10910 0.422454:0.346897:0.355947:...           ExE mesoderm  0.0984285
...                                   ...                    ...        ...
cell_11597 0.323805:0.292967:0.300485:...                    NMP  0.1663369
cell_9807  0.464466:0.374189:0.381698:...             Mesenchyme  0.0833019
cell_10095 0.341721:0.288215:0.485324:...            Endothelium  0.0889931
cell_11706 0.267487:0.240215:0.286012:...             Erythroid2  0.0350557
cell_11860 0.345786:0.343437:0.313994:... Forebrain/Midbrain/H..  0.0117001
                    pruned.labels
                      <character>
cell_11995 Forebrain/Midbrain/H..
cell_10294             Erythroid3
cell_9963             Endothelium
cell_11610             Erythroid3
cell_10910           ExE mesoderm
...                           ...
cell_11597                    NMP
cell_9807              Mesenchyme
cell_10095            Endothelium
cell_11706             Erythroid2
cell_11860 Forebrain/Midbrain/H..
```

We inspect the results using a heatmap of the per-cell and label scores.
Ideally, each cell should exhibit a high score in one label relative to
all of the others, indicating that the assignment to that label was
unambiguous. This is largely the case for mesenchyme and endothelial
cells, whereas we see expectedly more ambiguity between the two
erythroid cell populations.


```r
plotScoreHeatmap(res)
```

<img src="fig/cell_type_annotation-rendered-score-heat-1.png" style="display: block; margin: auto;" />

We also compare the cell type assignments with the clustering results to
determine the identity of each cluster. Here, several cell type classes
are nested within the same cluster, indicating that these clusters are
composed of several transcriptomically similar cell populations (such as
cluster 4 and 6). On the other hand, there are also instances where we
have several clusters for the same cell type, indicating that the
clustering represents finer subdivisions within these cell types.


```r
library(pheatmap)
tab <- table(anno = res$pruned.labels, cluster = colLabels(sce))
pheatmap(log2(tab + 10), color = colorRampPalette(c("white", "blue"))(101))
```

<img src="fig/cell_type_annotation-rendered-unnamed-chunk-3-1.png" style="display: block; margin: auto;" />

As it so happens, we are in the fortunate position where our test
dataset also contains independently defined labels. We see strong
consistency between the two sets of labels, indicating that our
automatic annotation is comparable to that generated manually by domain
experts.


```r
tab <- table(res$pruned.labels, sce$celltype.mapped)
pheatmap(log2(tab + 10), color = colorRampPalette(c("white", "blue"))(101))
```

<img src="fig/cell_type_annotation-rendered-anno-vs-preanno-1.png" style="display: block; margin: auto;" />

### Assigning cell labels from gene sets

A related strategy is to explicitly identify sets of marker genes that
are highly expressed in each individual cell. This does not require
matching of individual cells to the expression values of the reference
dataset, which is faster and more convenient when only the identities of
the markers are available. We demonstrate this approach using cell type
markers derived from the mouse embryo atlas dataset.


```r
library(scran)
wilcox.z <- pairwiseWilcox(ref, ref$celltype, lfc = 1, direction = "up")
markers.z <- getTopMarkers(wilcox.z$statistics, wilcox.z$pairs, 
                           pairwise = FALSE, n = 50)
lengths(markers.z)
```

```{.output}
                     Allantois                 Cardiomyocytes 
                           106                            106 
                   Endothelium                     Erythroid2 
                           103                             54 
                    Erythroid3                   ExE endoderm 
                            84                            102 
                  ExE mesoderm   Forebrain/Midbrain/Hindbrain 
                            97                             97 
                           Gut Haematoendothelial progenitors 
                            90                             71 
         Intermediate mesoderm                     Mesenchyme 
                            70                            118 
                  Neural crest                            NMP 
                            66                             91 
             Paraxial mesoderm            Pharyngeal mesoderm 
                            88                             85 
              Somitic mesoderm                    Spinal cord 
                            86                             91 
              Surface ectoderm 
                            92 
```

Our test dataset will be as before the wild-type chimera dataset.


```r
sce
```

```{.output}
class: SingleCellExperiment 
dim: 29411 1000 
metadata(0):
assays(2): counts logcounts
rownames(29411): Xkr4 Gm1992 ... Vmn2r122 CAAA01147332.1
rowData names(2): ENSEMBL SYMBOL
colnames(1000): cell_11995 cell_10294 ... cell_11706 cell_11860
colData names(12): cell barcode ... sizeFactor label
reducedDimNames(4): pca.corrected.E7.5 pca.corrected.E8.5 PCA UMAP
mainExpName: NULL
altExpNames(0):
```

We use the *[AUCell](https://bioconductor.org/packages/3.18/AUCell)* package to identify marker sets that
are highly expressed in each cell. This method ranks genes by their
expression values within each cell and constructs a response curve of
the number of genes from each marker set that are present with
increasing rank. It then computes the area under the curve (AUC) for
each marker set, quantifying the enrichment of those markers among the
most highly expressed genes in that cell. This is roughly similar to
performing a Wilcoxon rank sum test between genes in and outside of the
set, but involving only the top ranking genes by expression in each
cell.


```r
library(GSEABase)
all.sets <- lapply(names(markers.z), 
                   function(x) GeneSet(markers.z[[x]], setName = x))
all.sets <- GeneSetCollection(all.sets)
all.sets
```

```{.output}
GeneSetCollection
  names: Allantois, Cardiomyocytes, ..., Surface ectoderm (19 total)
  unique identifiers: Prrx2, Spin2c, ..., Igf2bp3 (560 total)
  types in collection:
    geneIdType: NullIdentifier (1 total)
    collectionType: NullCollection (1 total)
```


```r
library(AUCell)
rankings <- AUCell_buildRankings(as.matrix(counts(sce)),
                                 plotStats = FALSE, verbose = FALSE)
cell.aucs <- AUCell_calcAUC(all.sets, rankings)
results <- t(assay(cell.aucs))
head(results)
```

```{.output}
            gene sets
cells         Allantois Cardiomyocytes Endothelium Erythroid2 Erythroid3
  cell_11995 0.21994609      0.2336827   0.1947563 0.12228508  0.1639749
  cell_10294 0.10066221      0.1246414   0.1104634 0.60732017  0.6081721
  cell_9963  0.34273069      0.3158431   0.5407815 0.13588372  0.1797007
  cell_11610 0.08244651      0.1340388   0.1048735 0.57892981  0.5706619
  cell_10910 0.31763336      0.2784001   0.2445795 0.08175859  0.1354818
  cell_11021 0.20549732      0.2320122   0.1869839 0.11351012  0.1769588
            gene sets
cells        ExE endoderm ExE mesoderm Forebrain/Midbrain/Hindbrain       Gut
  cell_11995   0.06053637    0.4371076                    0.6319109 0.3779960
  cell_10294   0.06369959    0.1900184                    0.2670988 0.1482988
  cell_9963    0.09466189    0.4165978                    0.4972886 0.3805059
  cell_11610   0.05720738    0.2060333                    0.2785028 0.1466698
  cell_10910   0.11787498    0.5770948                    0.5721069 0.4450836
  cell_11021   0.09285926    0.4736976                    0.6011933 0.3993608
            gene sets
cells        Haematoendothelial progenitors Intermediate mesoderm Mesenchyme
  cell_11995                      0.3329538             0.5396029  0.2456757
  cell_10294                      0.1621829             0.1781659  0.1078670
  cell_9963                       0.6010306             0.4623277  0.3566887
  cell_11610                      0.1526329             0.1968055  0.1118536
  cell_10910                      0.4333710             0.5481614  0.3512731
  cell_11021                      0.3336212             0.5197592  0.2403682
            gene sets
cells        Neural crest       NMP Paraxial mesoderm Pharyngeal mesoderm
  cell_11995    0.6535231 0.4841026         0.5510547           0.5076289
  cell_10294    0.2840264 0.2113283         0.2192031           0.2049349
  cell_9963     0.5374756 0.3740042         0.5007010           0.4618718
  cell_11610    0.2949987 0.2235049         0.2348724           0.2134701
  cell_10910    0.5472885 0.5225294         0.5506325           0.5564014
  cell_11021    0.5863399 0.5731637         0.4992671           0.4707530
            gene sets
cells        Somitic mesoderm Spinal cord Surface ectoderm
  cell_11995        0.4654095   0.5807133        0.5251957
  cell_10294        0.2141816   0.2290110        0.1874113
  cell_9963         0.4307335   0.4373858        0.4173470
  cell_11610        0.2339673   0.2564642        0.2054802
  cell_10910        0.5179489   0.5167149        0.5049140
  cell_11021        0.5186006   0.5526277        0.5052879
```

We assign cell type identity to each cell in the test dataset by taking
the marker set with the top AUC as the label for that cell. Our new
labels mostly agree with the original annotation (and, thus, also with
the reference-based annotation). Instances where the original annotation
is divided into several new label groups typically points to large
overlaps in their marker sets. In the absence of prior annotation, a
more general diagnostic check is to compare the assigned labels to
cluster identities, under the expectation that most cells of a single
cluster would have the same label (or, if multiple labels are present,
they should at least represent closely related cell states).


```r
new.labels <- colnames(results)[max.col(results)]
tab <- table(new.labels, sce$celltype.mapped)
tab
```

```{.output}
                                
new.labels                       Allantois Blood progenitors 1
  Allantois                             28                   0
  Cardiomyocytes                         0                   0
  Endothelium                            0                   0
  Erythroid2                             0                   0
  Erythroid3                             0                   0
  ExE mesoderm                           0                   0
  Forebrain/Midbrain/Hindbrain           1                   0
  Gut                                    0                   0
  Haematoendothelial progenitors         1                   0
  Intermediate mesoderm                  0                   0
  Mesenchyme                             0                   0
  Neural crest                          11                   4
  NMP                                    0                   0
  Paraxial mesoderm                      5                   0
  Pharyngeal mesoderm                    0                   0
  Somitic mesoderm                       0                   0
  Spinal cord                            0                   0
  Surface ectoderm                       0                   0
                                
new.labels                       Blood progenitors 2 Cardiomyocytes
  Allantois                                        0              0
  Cardiomyocytes                                   0             27
  Endothelium                                      0              0
  Erythroid2                                       1              0
  Erythroid3                                       0              0
  ExE mesoderm                                     0              0
  Forebrain/Midbrain/Hindbrain                     0              0
  Gut                                              0              0
  Haematoendothelial progenitors                   0              0
  Intermediate mesoderm                            0              0
  Mesenchyme                                       0              0
  Neural crest                                     6              1
  NMP                                              0              0
  Paraxial mesoderm                                0              4
  Pharyngeal mesoderm                              0              6
  Somitic mesoderm                                 0              0
  Spinal cord                                      0              0
  Surface ectoderm                                 0              0
                                
new.labels                       Caudal epiblast Caudal Mesoderm Def. endoderm
  Allantois                                    0               0             0
  Cardiomyocytes                               0               0             0
  Endothelium                                  0               0             0
  Erythroid2                                   0               0             0
  Erythroid3                                   0               0             0
  ExE mesoderm                                 0               0             0
  Forebrain/Midbrain/Hindbrain                 1               0             2
  Gut                                          0               0             0
  Haematoendothelial progenitors               0               0             0
  Intermediate mesoderm                        0               0             0
  Mesenchyme                                   0               0             0
  Neural crest                                 0               0             1
  NMP                                          0               0             0
  Paraxial mesoderm                            0               0             1
  Pharyngeal mesoderm                          0               0             0
  Somitic mesoderm                             0               1             1
  Spinal cord                                  0               0             0
  Surface ectoderm                             0               0             0
                                
new.labels                       Doublet Endothelium Erythroid1 Erythroid2
  Allantois                            0           0          0          0
  Cardiomyocytes                       0           0          0          0
  Endothelium                          0          13          0          0
  Erythroid2                           0           0         14         21
  Erythroid3                           0           0          6          5
  ExE mesoderm                         0           0          0          0
  Forebrain/Midbrain/Hindbrain         6           0          0          0
  Gut                                  0           0          0          0
  Haematoendothelial progenitors       0          19          0          0
  Intermediate mesoderm               13           0          0          0
  Mesenchyme                           0           0          0          0
  Neural crest                        20           1          3          0
  NMP                                  0           0          0          0
  Paraxial mesoderm                    5           0          0          0
  Pharyngeal mesoderm                  1           0          0          0
  Somitic mesoderm                     1           0          0          0
  Spinal cord                          1           0          0          0
  Surface ectoderm                     0           0          0          0
                                
new.labels                       Erythroid3 ExE mesoderm
  Allantois                               0            0
  Cardiomyocytes                          0            0
  Endothelium                             0            0
  Erythroid2                             51            0
  Erythroid3                             49            0
  ExE mesoderm                            0           23
  Forebrain/Midbrain/Hindbrain            0            1
  Gut                                     0            0
  Haematoendothelial progenitors          0            0
  Intermediate mesoderm                   0           16
  Mesenchyme                              0            0
  Neural crest                            0           17
  NMP                                     0            0
  Paraxial mesoderm                       0            0
  Pharyngeal mesoderm                     0            0
  Somitic mesoderm                        0            0
  Spinal cord                             0            0
  Surface ectoderm                        0            0
                                
new.labels                       Forebrain/Midbrain/Hindbrain Gut
  Allantois                                                 0   0
  Cardiomyocytes                                            0   0
  Endothelium                                               0   0
  Erythroid2                                                0   0
  Erythroid3                                                0   0
  ExE mesoderm                                              0   0
  Forebrain/Midbrain/Hindbrain                             83   4
  Gut                                                       0   5
  Haematoendothelial progenitors                            0   0
  Intermediate mesoderm                                     0   0
  Mesenchyme                                                0   0
  Neural crest                                             20   0
  NMP                                                       0   0
  Paraxial mesoderm                                         0   0
  Pharyngeal mesoderm                                       0   0
  Somitic mesoderm                                          0   0
  Spinal cord                                               0   0
  Surface ectoderm                                          0  12
                                
new.labels                       Haematoendothelial progenitors
  Allantois                                                   0
  Cardiomyocytes                                              0
  Endothelium                                                 0
  Erythroid2                                                  0
  Erythroid3                                                  0
  ExE mesoderm                                                0
  Forebrain/Midbrain/Hindbrain                                0
  Gut                                                         0
  Haematoendothelial progenitors                             16
  Intermediate mesoderm                                       0
  Mesenchyme                                                  0
  Neural crest                                                8
  NMP                                                         0
  Paraxial mesoderm                                           2
  Pharyngeal mesoderm                                         0
  Somitic mesoderm                                            0
  Spinal cord                                                 0
  Surface ectoderm                                            0
                                
new.labels                       Intermediate mesoderm Mesenchyme Neural crest
  Allantois                                          0          0            0
  Cardiomyocytes                                     0          0            0
  Endothelium                                        0          0            0
  Erythroid2                                         0          0            0
  Erythroid3                                         0          0            0
  ExE mesoderm                                       1          7            0
  Forebrain/Midbrain/Hindbrain                       0          3            0
  Gut                                                0          0            0
  Haematoendothelial progenitors                     0          0            0
  Intermediate mesoderm                              9         15            0
  Mesenchyme                                         0         69            0
  Neural crest                                       4          4           10
  NMP                                                0          0            0
  Paraxial mesoderm                                  0          7            0
  Pharyngeal mesoderm                                0         13            0
  Somitic mesoderm                                   2          0            0
  Spinal cord                                        0          0            0
  Surface ectoderm                                   0          0            0
                                
new.labels                       NMP Paraxial mesoderm Pharyngeal mesoderm
  Allantois                        0                 0                   0
  Cardiomyocytes                   0                 0                   0
  Endothelium                      0                 0                   0
  Erythroid2                       0                 0                   0
  Erythroid3                       0                 0                   0
  ExE mesoderm                     0                 0                   0
  Forebrain/Midbrain/Hindbrain    41                 2                   4
  Gut                              0                 0                   0
  Haematoendothelial progenitors   0                 0                   0
  Intermediate mesoderm            0                 1                   5
  Mesenchyme                       0                 0                   0
  Neural crest                     5                20                  10
  NMP                             16                 0                   0
  Paraxial mesoderm                0                43                   3
  Pharyngeal mesoderm              0                 0                  18
  Somitic mesoderm                 0                 0                   0
  Spinal cord                      0                 0                   0
  Surface ectoderm                 0                 0                   0
                                
new.labels                       Rostral neurectoderm Somitic mesoderm
  Allantois                                         0                0
  Cardiomyocytes                                    0                0
  Endothelium                                       0                0
  Erythroid2                                        0                0
  Erythroid3                                        0                0
  ExE mesoderm                                      0                0
  Forebrain/Midbrain/Hindbrain                     11                1
  Gut                                               0                0
  Haematoendothelial progenitors                    0                0
  Intermediate mesoderm                             0                2
  Mesenchyme                                        0                0
  Neural crest                                      0                7
  NMP                                               0                0
  Paraxial mesoderm                                 0                0
  Pharyngeal mesoderm                               0                0
  Somitic mesoderm                                  0               25
  Spinal cord                                       0                0
  Surface ectoderm                                  0                0
                                
new.labels                       Spinal cord Stripped Surface ectoderm
  Allantois                                0        0                0
  Cardiomyocytes                           0        0                0
  Endothelium                              0        1                0
  Erythroid2                               0        1                0
  Erythroid3                               0        1                0
  ExE mesoderm                             0        0                0
  Forebrain/Midbrain/Hindbrain            37        0               17
  Gut                                      0        0                2
  Haematoendothelial progenitors           0        0                0
  Intermediate mesoderm                    0        0                0
  Mesenchyme                               0        0                2
  Neural crest                             3        5                6
  NMP                                      0        0                0
  Paraxial mesoderm                        0        0                0
  Pharyngeal mesoderm                      0        0                0
  Somitic mesoderm                         0        0                0
  Spinal cord                              9        0                0
  Surface ectoderm                         0        0               20
```

As a diagnostic measure, we examine the distribution of AUCs across
cells for each label. In heterogeneous populations, the distribution for
each label should be bimodal with one high-scoring peak containing cells
of that cell type and a low-scoring peak containing cells of other
types. The gap between these two peaks can be used to derive a threshold
for whether a label is "active" for a particular cell. (In this case, we
simply take the single highest-scoring label per cell as the labels
should be mutually exclusive.) In populations where a particular cell
type is expected, lack of clear bimodality for the corresponding label
may indicate that its gene set is not sufficiently informative.


```r
par(mfrow = c(3,3))
AUCell_exploreThresholds(cell.aucs[1:9], plotHist = TRUE, assign = TRUE) 
```

<img src="fig/cell_type_annotation-rendered-auc-dist-1.png" style="display: block; margin: auto;" />

Shown is the distribution of AUCs in the wild-type chimera dataset for
each label in the embryo atlas dataset. The blue curve represents the
density estimate, the red curve represents a fitted two-component
mixture of normals, the pink curve represents a fitted three-component
mixture, and the grey curve represents a fitted normal distribution.
Vertical lines represent threshold estimates corresponding to each
estimate of the distribution.

## Session Info


```r
sessionInfo()
```

```{.output}
R version 4.3.3 (2024-02-29)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 22.04.4 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0

locale:
 [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
 [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
 [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
[10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   

time zone: UTC
tzcode source: system (glibc)

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] GSEABase_1.64.0              graph_1.80.0                
 [3] annotate_1.80.0              XML_3.99-0.16.1             
 [5] AnnotationDbi_1.64.1         pheatmap_1.0.12             
 [7] scran_1.30.2                 scater_1.30.1               
 [9] ggplot2_3.5.0                scuttle_1.12.0              
[11] bluster_1.12.0               SingleR_2.4.1               
[13] MouseGastrulationData_1.16.0 SpatialExperiment_1.12.0    
[15] SingleCellExperiment_1.24.0  SummarizedExperiment_1.32.0 
[17] Biobase_2.62.0               GenomicRanges_1.54.1        
[19] GenomeInfoDb_1.38.8          IRanges_2.36.0              
[21] S4Vectors_0.40.2             BiocGenerics_0.48.1         
[23] MatrixGenerics_1.14.0        matrixStats_1.2.0           
[25] BiocStyle_2.30.0             AUCell_1.24.0               

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3            jsonlite_1.8.8               
  [3] magrittr_2.0.3                ggbeeswarm_0.7.2             
  [5] magick_2.8.3                  farver_2.1.1                 
  [7] rmarkdown_2.26                zlibbioc_1.48.2              
  [9] vctrs_0.6.5                   memoise_2.0.1                
 [11] DelayedMatrixStats_1.24.0     RCurl_1.98-1.14              
 [13] htmltools_0.5.7               S4Arrays_1.2.1               
 [15] AnnotationHub_3.10.0          curl_5.2.1                   
 [17] BiocNeighbors_1.20.2          SparseArray_1.2.4            
 [19] htmlwidgets_1.6.4             plotly_4.10.4                
 [21] cachem_1.0.8                  igraph_2.0.3                 
 [23] mime_0.12                     lifecycle_1.0.4              
 [25] pkgconfig_2.0.3               rsvd_1.0.5                   
 [27] Matrix_1.6-5                  R6_2.5.1                     
 [29] fastmap_1.1.1                 GenomeInfoDbData_1.2.11      
 [31] shiny_1.8.0                   digest_0.6.35                
 [33] colorspace_2.1-0              dqrng_0.3.2                  
 [35] irlba_2.3.5.1                 ExperimentHub_2.10.0         
 [37] RSQLite_2.3.5                 beachmat_2.18.1              
 [39] labeling_0.4.3                filelock_1.0.3               
 [41] fansi_1.0.6                   httr_1.4.7                   
 [43] abind_1.4-5                   compiler_4.3.3               
 [45] bit64_4.0.5                   withr_3.0.0                  
 [47] BiocParallel_1.36.0           viridis_0.6.5                
 [49] DBI_1.2.2                     highr_0.10                   
 [51] R.utils_2.12.3                MASS_7.3-60.0.1              
 [53] rappdirs_0.3.3                DelayedArray_0.28.0          
 [55] rjson_0.2.21                  tools_4.3.3                  
 [57] vipor_0.4.7                   beeswarm_0.4.0               
 [59] interactiveDisplayBase_1.40.0 httpuv_1.6.14                
 [61] R.oo_1.26.0                   glue_1.7.0                   
 [63] nlme_3.1-164                  promises_1.2.1               
 [65] grid_4.3.3                    cluster_2.1.6                
 [67] generics_0.1.3                gtable_0.3.4                 
 [69] R.methodsS3_1.8.2             tidyr_1.3.1                  
 [71] data.table_1.15.2             metapod_1.10.1               
 [73] BiocSingular_1.18.0           ScaledMatrix_1.10.0          
 [75] utf8_1.2.4                    XVector_0.42.0               
 [77] ggrepel_0.9.5                 BiocVersion_3.18.1           
 [79] pillar_1.9.0                  limma_3.58.1                 
 [81] BumpyMatrix_1.10.0            later_1.3.2                  
 [83] splines_4.3.3                 dplyr_1.1.4                  
 [85] BiocFileCache_2.10.1          lattice_0.22-6               
 [87] survival_3.5-8                FNN_1.1.4                    
 [89] renv_1.0.7                    bit_4.0.5                    
 [91] tidyselect_1.2.1              locfit_1.5-9.9               
 [93] Biostrings_2.70.3             knitr_1.45                   
 [95] gridExtra_2.3                 edgeR_4.0.16                 
 [97] xfun_0.42                     mixtools_2.0.0               
 [99] statmod_1.5.0                 lazyeval_0.2.2               
[101] yaml_2.3.8                    evaluate_0.23                
[103] codetools_0.2-19              kernlab_0.9-32               
[105] tibble_3.2.1                  BiocManager_1.30.22          
[107] cli_3.6.2                     uwot_0.1.16                  
[109] xtable_1.8-4                  segmented_2.0-3              
[111] munsell_0.5.0                 Rcpp_1.0.12                  
[113] dbplyr_2.5.0                  png_0.1-8                    
[115] parallel_4.3.3                ellipsis_0.3.2               
[117] blob_1.2.4                    sparseMatrixStats_1.14.0     
[119] bitops_1.0-7                  viridisLite_0.4.2            
[121] scales_1.3.0                  purrr_1.0.2                  
[123] crayon_1.5.2                  rlang_1.1.3                  
[125] cowplot_1.1.3                 KEGGREST_1.42.0              
```

## Exercises

::: challenge
#### Exercise 1: Clustering

The [Leiden
algorithm](https://www.nature.com/articles/s41598-019-41695-z) is
similar to the Louvain algorithm, but it is faster and has been shown to
result in better connected communities. Modify the above call to
`clusterCells` to carry out the community detection with the Leiden
algorithm instead. Visualize the results in a UMAP plot.

::: hint
The `NNGraphParam` constructor has an argument `cluster.args`. This
allows to specify arguments passed on to the `cluster_leiden` function
from the
[igraph](https://cran.r-project.org/web/packages/igraph/index.html)
package. Use the `cluster.args` argument to parameterize the clustering
to use modularity as the objective function and a resolution parameter
of 0.5.
:::

::: solution
TODO
:::
:::

::: challenge
#### Exercise 2: Cluster annotation

Another strategy for annotating the clusters is to perform a gene set
enrichment analysis on the marker genes defining each cluster. This
identifies the pathways and processes that are (relatively) active in
each cluster based on upregulation of the associated genes compared to
other clusters. Focus on the top 100 up-regulated genes in a cluster of
your choice and perform a gene set enrichment analysis of biological
process (BP) gene sets from the Gene Ontology (GO).

::: hint
Use the `goana()` function from the *[limma](https://bioconductor.org/packages/3.18/limma)* package to
identify GO BP terms that are overrepresented in the list of marker
genes.
:::

::: solution
TODO
:::
:::

::: challenge
#### Exercise 3: Workflow

The [scRNAseq](https://bioconductor.org/packages/scRNAseq) package
provides gene-level counts for a collection of public scRNA-seq
datasets, stored as `SingleCellExperiment` objects with annotated cell-
and gene-level metadata. Consult the vignette of the
[scRNAseq](https://bioconductor.org/packages/scRNAseq) package to
inspect all available datasets and select a dataset of your choice.
Perform a typical scRNA-seq analysis on this dataset including QC,
normalization, feature selection, dimensionality reduction, clustering,
and marker gene detection.

::: solution
TODO
:::
:::

::: checklist
## Further Reading

-   OSCA book, [Chapters
    5-7](https://bioconductor.org/books/release/OSCA.basic/clustering.html)
-   Assigning cell types with SingleR ([the
    book](https://bioconductor.org/books/release/SingleRBook/)).
-   The [AUCell](https://bioconductor.org/packages/AUCell) package
    vignette.
:::

::: keypoints
-   TODO
:::