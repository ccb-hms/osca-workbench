---
title: Cell type annotation
teaching: 30 # Minutes of teaching in the lesson
exercises: 15 # Minutes of exercises in the lesson
editor_options: 
  markdown: 
    wrap: 72
---

::: questions
-   How can we identify groups of cells with similar expression profiles?
-   How can we identify genes that drive separation between these groups of cells?
-   How to leverage reference datasets and known marker genes for the cell type annotation of new datasets?
:::

::: objectives
-   Identify groups of cells by clustering cells based on gene expression patterns.
-   Identify marker genes through testing for differential expression between clusters.
-   Annotate cell types through annotation transfer from reference datasets.
-   Annotate cell types through marker gene set enrichment testing.
:::

## Setup



Again we'll start by loading the libraries we'll be using:


``` r
library(AUCell)
library(MouseGastrulationData)
library(SingleR)
library(bluster)
library(scater)
library(scran)
library(pheatmap)
library(GSEABase)
```

## Data retrieval

We'll be using the fifth processed sample from the WT chimeric mouse embryo data: 


``` r
sce <- WTChimeraData(samples = 5, type = "processed")
sce
```

``` output
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

To speed up the computations, we take a random subset of 1,000 cells.


``` r
set.seed(123)

ind <- sample(ncol(sce), 1000)

sce <- sce[,ind]
```

## Preprocessing

The SCE object needs to contain log-normalized expression counts as well as PCA coordinates in the reduced dimensions, so we compute those here: 


``` r
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
`SingleCellExperiment` object. We use the `colLabels()` function to assign the
cluster labels as a factor in the column data.


``` r
colLabels(sce) <- clusterCells(sce, use.dimred = "PCA",
                               BLUSPARAM = NNGraphParam(cluster.fun = "louvain"))

table(colLabels(sce))
```

``` output

  1   2   3   4   5   6   7   8   9  10  11 
100 160  99 141  63  93  60 108  44  91  41 
```
You can see we ended up with 11 clusters of varying sizes.

We can now overlay the cluster labels as color on a UMAP plot:


``` r
sce <- runUMAP(sce, dimred = "PCA")

plotReducedDim(sce, "UMAP", color_by = "label")
```

<img src="fig/cell_type_annotation-rendered-cluster-viz-1.png" style="display: block; margin: auto;" />

:::: challenge

Our clusters look semi-reasonable, but what if we wanted to make them less granular? Look at the help documentation for `?clusterCells` and `?NNGraphParam` to find out what we'd need to change to get fewer, larger clusters.

::: solution

We see in the help documentation for `?clusterCells` that all of the clustering algorithm details are handled through the `BLUSPARAM` argument, which needs to provide a `BlusterParam` object (of which `NNGraphParam` is a sub-class). Each type of clustering algorithm will have some sort of hyper-parameter that controls the granularity of the output clusters. Looking at `?NNGraphParam` specifically, we see an argument called `k` which is described as "An integer scalar specifying the number of nearest neighbors to consider during graph construction." If the clustering process has to connect larger sets of neighbors, the graph will tend to be cut into larger groups, resulting in less granular clusters. Try the two code blocks above once more with `k = 20`. Given their visual differences, do you think one set of clusters is "right" and the other is "wrong"?  
:::

::::


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

Here, we use `findMarkers()` to perform a Wilcoxon rank sum test against a log2 fold change
threshold of 1, focusing on up-regulated (positive) markers in one
cluster when compared to another cluster.


``` r
rownames(sce) <- rowData(sce)$SYMBOL

markers <- findMarkers(sce, test.type = "wilcox", direction = "up", lfc = 1)

markers
```

``` output
List of length 11
names(11): 1 2 3 4 5 6 7 8 9 10 11
```

<!-- TODO change this ^ to scoreMarkers() -->

The resulting object contains a sorted marker gene list for each
cluster, in which the top genes are those that contribute the most to
the separation of that cluster from mall other clusters.

Here, we inspect the ranked marker gene list for the first cluster.


``` r
markers[[1]]
```

``` output
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


``` r
top.markers <- head(rownames(markers[[1]]))
plotExpression(sce, features = top.markers, x = "label", color_by = "label")
```

<img src="fig/cell_type_annotation-rendered-plot-markers-1.png" style="display: block; margin: auto;" />

Clearly, not every marker gene distinguishes cluster 1 from every other cluster. However, with a combination of multiple marker genes it's possible to clearly identify gene patterns that are unique to cluster 1. It's sort of like the 20 questions game - with answers to the right questions about a cell (e.g. "Do you highly express Ptn?"), you can clearly identify what cluster it falls in.

:::: challenge 

Why do you think marker genes are found by aggregating pairwise comparisons rather than iteratively comparing each cluster to all other clusters? 

::: solution

One important reason why is because averages over all other clusters can be sensitive to the cell type composition. If a rare cell type shows up in one sample, the most discriminative marker genes found in this way could be very different from those found in another sample where the rare cell type is absent. 

Generally, it's good to keep in mind that the concept of "everything else" is not a stable basis for comparison. Read that sentence again, because its a subtle but broadly applicable point. Think about it and you can probably identify analogous issues in fields outside of single-cell analysis. It frequently comes up when comparisons between multiple categories are involved.

:::
::::


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
*[SingleR](https://bioconductor.org/packages/3.19/SingleR)* method for cell type annotation [Aran et al.,
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


``` r
ref <- EmbryoAtlasData(samples = 29)

ref
```

``` output
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

In order to reduce the computational load, we subsample the dataset to 1,000 cells.


``` r
set.seed(123)

ind <- sample(ncol(ref), 1000)

ref <- ref[,ind]
```

You can see we have an assortment of different cell types in the reference (with varying frequency):


``` r
tab <- sort(table(ref$celltype), decreasing = TRUE)
tab
```

``` output

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

We need the normalized log counts, so we add those on: 


``` r
ref <- logNormCounts(ref)
```

Some cleaning - remove cells of the reference dataset for which the cell
type annotation is missing:


``` r
nna <- !is.na(ref$celltype)

ref <- ref[,nna]
```

Also remove cell types of very low abundance (here less than 10 cells)
to remove noise prior to subsequent annotation tasks.


``` r
abu.ct <- names(tab)[tab >= 10]

ind <- ref$celltype %in% abu.ct

ref <- ref[,ind] 
```

Restrict to genes shared between query and reference dataset.


``` r
rownames(ref) <- rowData(ref)$SYMBOL

shared_genes <- intersect(rownames(sce), rownames(ref))

sce <- sce[shared_genes,]

ref <- ref[shared_genes,]
```

Convert sparse assay matrices to regular dense matrices for input to
SingleR:


``` r
sce.mat <- as.matrix(assay(sce, "logcounts"))

ref.mat <- as.matrix(assay(ref, "logcounts"))
```

Finally, run SingleR with the query and reference datasets:


``` r
res <- SingleR(test = sce.mat, 
               ref = ref.mat,
               labels = ref$celltype)
res
```

``` output
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
unambiguous. 


``` r
plotScoreHeatmap(res)
```

<img src="fig/cell_type_annotation-rendered-score-heat-1.png" style="display: block; margin: auto;" />

We obtained fairly unambiguous predictions for mesenchyme and endothelial
cells, whereas we see expectedly more ambiguity between the two
erythroid cell populations.

We can also compare the cell type assignments with the unsupervised clustering
results to determine the identity of each cluster. Here, several cell type
classes are nested within the same cluster, indicating that these clusters are
composed of several transcriptomically similar cell populations. On the other
hand, there are also instances where we have several clusters for the same cell
type, indicating that the clustering represents finer subdivisions within these
cell types.


``` r
tab <- table(anno = res$pruned.labels, cluster = colLabels(sce))

pheatmap(log2(tab + 10), color = colorRampPalette(c("white", "blue"))(101))
```

<img src="fig/cell_type_annotation-rendered-unnamed-chunk-3-1.png" style="display: block; margin: auto;" />

As it so happens, we are in the fortunate position where our test
dataset also contains independently defined labels. We see strong
consistency between the two sets of labels, indicating that our
automatic annotation is comparable to that generated manually by domain
experts.


``` r
tab <- table(res$pruned.labels, sce$celltype.mapped)

pheatmap(log2(tab + 10), color = colorRampPalette(c("white", "blue"))(101))
```

<img src="fig/cell_type_annotation-rendered-anno-vs-preanno-1.png" style="display: block; margin: auto;" />

:::: challenge

SingleR can be computationally expensive. How do you set it to run in parallel?

::: solution

Use `BiocParallel` and the `BPPARAM` argument! This example will set it to use four cores on your laptop, but you can also configure BiocParallel to use cluster jobs.


``` r
library(BiocParallel)

my_bpparam = MulticoreParam(workers = 4)

res <- SingleR(test = sce.mat, 
               ref = ref.mat,
               labels = ref$celltype,
               BPPARAM = my_bpparam)
```

`BiocParallel` is the most common way to enable parallel computation in Bioconductor packages, so you can expect to see it elsewhere outside of SingleR.

:::

::::

### Assigning cell labels from gene sets

A related strategy is to explicitly identify sets of marker genes that
are highly expressed in each individual cell. This does not require
matching of individual cells to the expression values of the reference
dataset, which is faster and more convenient when only the identities of
the markers are available. We demonstrate this approach using cell type
markers derived from the mouse embryo atlas dataset.


``` r
wilcox.z <- pairwiseWilcox(ref, ref$celltype, lfc = 1, direction = "up")

markers.z <- getTopMarkers(wilcox.z$statistics, wilcox.z$pairs, 
                           pairwise = FALSE, n = 50)

lengths(markers.z)
```

``` output
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


``` r
sce
```

``` output
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

We use the *[AUCell](https://bioconductor.org/packages/3.19/AUCell)* package to identify marker sets that
are highly expressed in each cell. This method ranks genes by their
expression values within each cell and constructs a response curve of
the number of genes from each marker set that are present with
increasing rank. It then computes the area under the curve (AUC) for
each marker set, quantifying the enrichment of those markers among the
most highly expressed genes in that cell. This is roughly similar to
performing a Wilcoxon rank sum test between genes in and outside of the
set, but involving only the top ranking genes by expression in each
cell.


``` r
all.sets <- lapply(names(markers.z), 
                   function(x) GeneSet(markers.z[[x]], setName = x))

all.sets <- GeneSetCollection(all.sets)

all.sets
```

``` output
GeneSetCollection
  names: Allantois, Cardiomyocytes, ..., Surface ectoderm (19 total)
  unique identifiers: Prrx2, Spin2c, ..., Igf2bp3 (560 total)
  types in collection:
    geneIdType: NullIdentifier (1 total)
    collectionType: NullCollection (1 total)
```


``` r
rankings <- AUCell_buildRankings(as.matrix(counts(sce)),
                                 plotStats = FALSE, verbose = FALSE)

cell.aucs <- AUCell_calcAUC(all.sets, rankings)

results <- t(assay(cell.aucs))

head(results)
```

``` output
            gene sets
cells        Allantois Cardiomyocytes Endothelium Erythroid2 Erythroid3
  cell_11995    0.2199          0.234       0.195     0.1223      0.164
  cell_10294    0.1007          0.125       0.110     0.6073      0.608
  cell_9963     0.3427          0.316       0.541     0.1359      0.180
  cell_11610    0.0824          0.134       0.105     0.5789      0.571
  cell_10910    0.3176          0.278       0.245     0.0818      0.135
  cell_11021    0.2055          0.232       0.187     0.1135      0.177
            gene sets
cells        ExE endoderm ExE mesoderm Forebrain/Midbrain/Hindbrain   Gut
  cell_11995       0.0605        0.437                        0.632 0.378
  cell_10294       0.0637        0.190                        0.267 0.148
  cell_9963        0.0947        0.417                        0.497 0.381
  cell_11610       0.0572        0.206                        0.279 0.147
  cell_10910       0.1179        0.577                        0.572 0.445
  cell_11021       0.0929        0.474                        0.601 0.399
            gene sets
cells        Haematoendothelial progenitors Intermediate mesoderm Mesenchyme
  cell_11995                          0.333                 0.540      0.246
  cell_10294                          0.162                 0.178      0.108
  cell_9963                           0.601                 0.462      0.357
  cell_11610                          0.153                 0.197      0.112
  cell_10910                          0.433                 0.548      0.351
  cell_11021                          0.334                 0.520      0.240
            gene sets
cells        Neural crest   NMP Paraxial mesoderm Pharyngeal mesoderm
  cell_11995        0.654 0.484             0.551               0.508
  cell_10294        0.284 0.211             0.219               0.205
  cell_9963         0.537 0.374             0.501               0.462
  cell_11610        0.295 0.224             0.235               0.213
  cell_10910        0.547 0.523             0.551               0.556
  cell_11021        0.586 0.573             0.499               0.471
            gene sets
cells        Somitic mesoderm Spinal cord Surface ectoderm
  cell_11995            0.465       0.581            0.525
  cell_10294            0.214       0.229            0.187
  cell_9963             0.431       0.437            0.417
  cell_11610            0.234       0.256            0.205
  cell_10910            0.518       0.517            0.505
  cell_11021            0.519       0.553            0.505
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


``` r
new.labels <- colnames(results)[max.col(results)]

tab <- table(new.labels, sce$celltype.mapped)

tab
```

``` output
                                
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


``` r
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

:::: challenge

The diagnostics don't look so good for some of the examples here. Which ones? Why?

::: solution

The example that jumps out most strongly to the eye is ExE endoderm, which doesn't show clear separate modes. Simultaneously, Endothelium seems to have three or four modes. 

Remember, this is an exploratory diagnostic, not the final word! At this point it'd be good to engage in some critical inspection of the results. Maybe we don't have enough / the best marker genes. In this particular case, the fact that we subsetted the reference set to 1000 cells probably didn't help.
:::

::::

## Session Info


``` r
sessionInfo()
```

``` output
R version 4.4.1 (2024-06-14)
Platform: x86_64-pc-linux-gnu
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
 [1] GSEABase_1.66.0              graph_1.82.0                
 [3] annotate_1.82.0              XML_3.99-0.16.1             
 [5] AnnotationDbi_1.66.0         pheatmap_1.0.12             
 [7] scran_1.32.0                 scater_1.32.0               
 [9] ggplot2_3.5.1                scuttle_1.14.0              
[11] bluster_1.14.0               SingleR_2.6.0               
[13] MouseGastrulationData_1.18.0 SpatialExperiment_1.14.0    
[15] SingleCellExperiment_1.26.0  SummarizedExperiment_1.34.0 
[17] Biobase_2.64.0               GenomicRanges_1.56.0        
[19] GenomeInfoDb_1.40.1          IRanges_2.38.0              
[21] S4Vectors_0.42.0             BiocGenerics_0.50.0         
[23] MatrixGenerics_1.16.0        matrixStats_1.3.0           
[25] AUCell_1.26.0                BiocStyle_2.32.0            

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3        jsonlite_1.8.8           
  [3] magrittr_2.0.3            ggbeeswarm_0.7.2         
  [5] magick_2.8.3              farver_2.1.2             
  [7] rmarkdown_2.27            zlibbioc_1.50.0          
  [9] vctrs_0.6.5               memoise_2.0.1            
 [11] DelayedMatrixStats_1.26.0 htmltools_0.5.8.1        
 [13] S4Arrays_1.4.1            AnnotationHub_3.12.0     
 [15] curl_5.2.1                BiocNeighbors_1.22.0     
 [17] SparseArray_1.4.8         htmlwidgets_1.6.4        
 [19] plotly_4.10.4             cachem_1.1.0             
 [21] igraph_2.0.3              mime_0.12                
 [23] lifecycle_1.0.4           pkgconfig_2.0.3          
 [25] rsvd_1.0.5                Matrix_1.7-0             
 [27] R6_2.5.1                  fastmap_1.2.0            
 [29] GenomeInfoDbData_1.2.12   digest_0.6.35            
 [31] colorspace_2.1-0          dqrng_0.4.1              
 [33] irlba_2.3.5.1             ExperimentHub_2.12.0     
 [35] RSQLite_2.3.7             beachmat_2.20.0          
 [37] labeling_0.4.3            filelock_1.0.3           
 [39] fansi_1.0.6               httr_1.4.7               
 [41] abind_1.4-5               compiler_4.4.1           
 [43] bit64_4.0.5               withr_3.0.0              
 [45] BiocParallel_1.38.0       viridis_0.6.5            
 [47] DBI_1.2.3                 highr_0.11               
 [49] R.utils_2.12.3            MASS_7.3-60.2            
 [51] rappdirs_0.3.3            DelayedArray_0.30.1      
 [53] rjson_0.2.21              tools_4.4.1              
 [55] vipor_0.4.7               beeswarm_0.4.0           
 [57] R.oo_1.26.0               glue_1.7.0               
 [59] nlme_3.1-164              grid_4.4.1               
 [61] cluster_2.1.6             generics_0.1.3           
 [63] gtable_0.3.5              R.methodsS3_1.8.2        
 [65] tidyr_1.3.1               data.table_1.15.4        
 [67] BiocSingular_1.20.0       ScaledMatrix_1.12.0      
 [69] metapod_1.12.0            utf8_1.2.4               
 [71] XVector_0.44.0            ggrepel_0.9.5            
 [73] BiocVersion_3.19.1        pillar_1.9.0             
 [75] limma_3.60.2              BumpyMatrix_1.12.0       
 [77] splines_4.4.1             dplyr_1.1.4              
 [79] BiocFileCache_2.12.0      lattice_0.22-6           
 [81] survival_3.6-4            FNN_1.1.4                
 [83] renv_1.0.7                bit_4.0.5                
 [85] tidyselect_1.2.1          locfit_1.5-9.9           
 [87] Biostrings_2.72.1         knitr_1.47               
 [89] gridExtra_2.3             edgeR_4.2.0              
 [91] xfun_0.44                 mixtools_2.0.0           
 [93] statmod_1.5.0             UCSC.utils_1.0.0         
 [95] lazyeval_0.2.2            yaml_2.3.8               
 [97] evaluate_0.23             codetools_0.2-20         
 [99] kernlab_0.9-32            tibble_3.2.1             
[101] BiocManager_1.30.23       cli_3.6.2                
[103] uwot_0.2.2                xtable_1.8-4             
[105] segmented_2.1-0           munsell_0.5.1            
[107] Rcpp_1.0.12               dbplyr_2.5.0             
[109] png_0.1-8                 parallel_4.4.1           
[111] blob_1.2.4                sparseMatrixStats_1.16.0 
[113] viridisLite_0.4.2         scales_1.3.0             
[115] purrr_1.0.2               crayon_1.5.2             
[117] rlang_1.1.3               cowplot_1.1.3            
[119] KEGGREST_1.44.0          
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
Use the `goana()` function from the *[limma](https://bioconductor.org/packages/3.19/limma)* package to
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
-   The two main approaches for cell type annotation are 1) manual annotation
    of clusters based on marker gene expression, and 2) computational annotation
    based on annotation transfer from reference datasets or marker gene set enrichment testing.
-   For manual annotation, cells are first clustered with unsupervised methods
    such as graph-based clustering followed by community detection algorithms such
    as Louvain or Leiden.
-   The `clusterCells` function from the *[scran](https://bioconductor.org/packages/3.19/scran)* package provides different
    algorithms that are commonly used for the clustering of scRNA-seq data.
-   Once clusters have been obtained, cell type labels are then manually
    assigned to cell clusters by matching cluster-specific upregulated marker
    genes with prior knowledge of cell-type markers.
-   The `findMarkers` function from the *[scran](https://bioconductor.org/packages/3.19/scran)* package 
    package can be used to find candidate marker genes for clusters of cells by
    testing for differential expression between pairs of clusters.
-   Computational annotation using published reference datasets or curated gene sets
    provides a fast, automated, and reproducible alternative to the manual
    annotation of cell clusters based on marker gene expression.
-   The *[SingleR](https://bioconductor.org/packages/3.19/SingleR)*
    package is a popular choice for reference-based annotation and assigns labels
    to cells based on the reference samples with the highest Spearman rank correlations.
-   The *[AUCell](https://bioconductor.org/packages/3.19/AUCell)* package provides an enrichment
    test to identify curated marker sets that are highly expressed in each cell. 
:::
