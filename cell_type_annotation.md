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

We see in the help documentation for `?clusterCells` that all of the clustering algorithm details are handled through the `BLUSPARAM` argument, which needs to provide a `BlusterParam` object (of which `NNGraphParam` is a sub-class). Each type of clustering algorithm will have some sort of hyper-parameter that controls the granularity of the output clusters. Looking at `?NNGraphParam` specifically, we see an argument called `k` which is described as "An integer scalar specifying the number of nearest neighbors to consider during graph construction." If the clustering process has to connect larger sets of neighbors, the graph will tend to be cut into larger groups, resulting in less granular clusters. Try the two code blocks above once more with `k = 30`. Given their visual differences, do you think one set of clusters is "right" and the other is "wrong"?


``` r
sce$clust2 <- clusterCells(sce, use.dimred = "PCA",
                           BLUSPARAM = NNGraphParam(cluster.fun = "louvain",
                                                    k = 30))

plotReducedDim(sce, "UMAP", color_by = "clust2")
```

<img src="fig/cell_type_annotation-rendered-unnamed-chunk-2-1.png" style="display: block; margin: auto;" />

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

Here, we use `scoreMarkers()` to perform pairwise comparisons of gene
expression, focusing on up-regulated (positive) markers in one cluster when
compared to another cluster.


``` r
rownames(sce) <- rowData(sce)$SYMBOL

markers <- scoreMarkers(sce)

markers
```

``` output
List of length 11
names(11): 1 2 3 4 5 6 7 8 9 10 11
```

The resulting object contains a sorted marker gene list for each
cluster, in which the top genes are those that contribute the most to
the separation of that cluster from all other clusters.

Here, we inspect the ranked marker gene list for the first cluster.


``` r
markers[[1]]
```

``` output
DataFrame with 29453 rows and 19 columns
               self.average other.average self.detected other.detected
                  <numeric>     <numeric>     <numeric>      <numeric>
Xkr4              0.0000000    0.00366101          0.00     0.00373784
Gm1992            0.0000000    0.00000000          0.00     0.00000000
Gm37381           0.0000000    0.00000000          0.00     0.00000000
Rp1               0.0000000    0.00000000          0.00     0.00000000
Sox17             0.0279547    0.18822927          0.02     0.09348375
...                     ...           ...           ...            ...
AC149090.1        0.3852624   0.352021067          0.33      0.2844935
DHRSX             0.4108022   0.491424091          0.35      0.3882325
Vmn2r122          0.0000000   0.000000000          0.00      0.0000000
CAAA01147332.1    0.0164546   0.000802687          0.01      0.0010989
tomato-td         0.6341678   0.624350570          0.51      0.4808379
               mean.logFC.cohen min.logFC.cohen median.logFC.cohen
                      <numeric>       <numeric>          <numeric>
Xkr4                 -0.0386672       -0.208498          0.0000000
Gm1992                0.0000000        0.000000          0.0000000
Gm37381               0.0000000        0.000000          0.0000000
Rp1                   0.0000000        0.000000          0.0000000
Sox17                -0.1383820       -1.292067          0.0324795
...                         ...             ...                ...
AC149090.1            0.0644403      -0.1263241          0.0366957
DHRSX                -0.1154163      -0.4619613         -0.1202781
Vmn2r122              0.0000000       0.0000000          0.0000000
CAAA01147332.1        0.1338463       0.0656709          0.1414214
tomato-td             0.0220121      -0.2535145          0.0196130
               max.logFC.cohen rank.logFC.cohen  mean.AUC   min.AUC median.AUC
                     <numeric>        <integer> <numeric> <numeric>  <numeric>
Xkr4                  0.000000             6949  0.498131  0.489247   0.500000
Gm1992                0.000000             6554  0.500000  0.500000   0.500000
Gm37381               0.000000             6554  0.500000  0.500000   0.500000
Rp1                   0.000000             6554  0.500000  0.500000   0.500000
Sox17                 0.200319             1482  0.462912  0.228889   0.499575
...                        ...              ...       ...       ...        ...
AC149090.1            0.427051             1685  0.518060  0.475000   0.508779
DHRSX                 0.130189             3431  0.474750  0.389878   0.471319
Vmn2r122              0.000000             6554  0.500000  0.500000   0.500000
CAAA01147332.1        0.141421             2438  0.504456  0.499560   0.505000
tomato-td             0.318068             2675  0.502868  0.427083   0.501668
                 max.AUC  rank.AUC mean.logFC.detected min.logFC.detected
               <numeric> <integer>           <numeric>          <numeric>
Xkr4                0.50      6882        -2.58496e-01       -1.58496e+00
Gm1992              0.50      6513        -8.00857e-17       -3.20343e-16
Gm37381             0.50      6513        -8.00857e-17       -3.20343e-16
Rp1                 0.50      6513        -8.00857e-17       -3.20343e-16
Sox17               0.51      3957        -4.48729e-01       -4.23204e+00
...                  ...       ...                 ...                ...
AC149090.1      0.588462      1932         2.34565e-01       -4.59278e-02
DHRSX           0.530054      2050        -1.27333e-01       -4.52151e-01
Vmn2r122        0.500000      6513        -8.00857e-17       -3.20343e-16
CAAA01147332.1  0.505000      4893         7.27965e-01       -6.64274e-02
tomato-td       0.576875      1840         9.80090e-02       -2.20670e-01
               median.logFC.detected max.logFC.detected rank.logFC.detected
                           <numeric>          <numeric>           <integer>
Xkr4                      0.00000000        3.20343e-16                5560
Gm1992                    0.00000000        3.20343e-16                5560
Gm37381                   0.00000000        3.20343e-16                5560
Rp1                       0.00000000        3.20343e-16                5560
Sox17                    -0.00810194        1.51602e+00                 341
...                              ...                ...                 ...
AC149090.1                 0.0821121        9.55592e-01                2039
DHRSX                     -0.1774204        2.28269e-01                3943
Vmn2r122                   0.0000000        3.20343e-16                5560
CAAA01147332.1             0.8267364        1.00000e+00                 898
tomato-td                  0.0805999        4.63438e-01                3705
```

Each column contains summary statistics for each gene in the given cluster.
These are usually the mean/median/min/max of statistics like Cohen's *d* and AUC
when comparing this cluster (cluster 1 in this case) to all other clusters.
`mean.AUC` is usually the most important to check. AUC is the probability that a
randomly selected cell in cluster *A* has a greater expression of gene
*X* than a randomly selected cell in cluster *B*. You can set `full.stats=TRUE` if you'd like the marker data frames to retain list columns containing each statistic for each pairwise comparison.

We can then inspect the top marker genes for the first cluster using the
`plotExpression` function from the
[scater](https://bioconductor.org/packages/scater) package.


``` r
c1_markers <- markers[[1]]

ord <- order(c1_markers$mean.AUC, 
             decreasing = TRUE)

top.markers <- head(rownames(c1_markers[ord,]))

plotExpression(sce, 
               features = top.markers, 
               x        = "label", 
               color_by = "label")
```

<img src="fig/cell_type_annotation-rendered-plot-markers-1.png" style="display: block; margin: auto;" />

Clearly, not every marker gene distinguishes cluster 1 from every other cluster. However, with a combination of multiple marker genes it's possible to clearly identify gene patterns that are unique to cluster 1. It's sort of like the 20 questions game - with answers to the right questions about a cell (e.g. "Do you highly express Ptn?"), you can clearly identify what cluster it falls in.

:::: challenge

Looking at the last plot, what clusters are most difficult to distinguish from cluster 1? Now re-run the UMAP plot from the previous section. Do the difficult-to-distinguish clusters make sense?

::: solution

You can see that at least among the top markers, cluster 6 (pale green) tends to have the least separation from cluster 1. 


``` r
plotReducedDim(sce, "UMAP", color_by = "label")
```

<img src="fig/cell_type_annotation-rendered-unnamed-chunk-3-1.png" style="display: block; margin: auto;" />

Looking at the UMAP again, we can see that the marker gene overlap of clusters 1 and 6 makes sense. They're right next to each other on the UMAP. They're probably closely related cell types, and a less granular clustering would probably lump them together.

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

Here we take a single sample from `EmbryoAtlasData` as our reference dataset. In practice you would want to take more/all samples, possibly with batch-effect correction (see the next episode).


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

<img src="fig/cell_type_annotation-rendered-unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

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

Assign the SingleR annotations as a column in the colData for the query object `sce`.

::: solution


``` r
sce$SingleR_label = res$pruned.labels
```

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

<!--- 

This version with scoreMarkers() produces worse looking diagnostics, so let's leave it with the pairwise Wilcox version.

``` r
ref_markers <- scoreMarkers(ref, groups = ref$celltype, lfc = 1)

get_top_markers <- function(marker_df, n = 100) {
  ord <- order(marker_df$mean.AUC, decreasing = TRUE)
  
  rownames(marker_df[ord,])[1:n]
}

markers.z <- lapply(ref_markers, get_top_markers)
```

-->

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
colData names(14): cell barcode ... clust2 SingleR_label
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
  unique identifiers: Phlda2, Spin2c, ..., Akr7a5 (991 total)
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
  cell_11995    0.0691         0.0536      0.0459     0.1056     0.0948
  cell_10294    0.0384         0.0370      0.0451     0.4910     0.5153
  cell_9963     0.2177         0.1011      0.4380     0.1070     0.1087
  cell_11610    0.0108         0.0428      0.0347     0.4687     0.4555
  cell_10910    0.1868         0.0766      0.0869     0.0766     0.0682
  cell_11021    0.0695         0.0596      0.0465     0.0991     0.0993
            gene sets
cells        ExE endoderm ExE mesoderm Forebrain/Midbrain/Hindbrain    Gut
  cell_11995       0.0342       0.1324                       0.3173 0.1011
  cell_10294       0.0680       0.0172                       0.0439 0.0458
  cell_9963        0.0686       0.1147                       0.1116 0.1237
  cell_11610       0.0573       0.0202                       0.0460 0.0268
  cell_10910       0.0764       0.3255                       0.1747 0.1816
  cell_11021       0.0680       0.2029                       0.2500 0.1277
            gene sets
cells        Haematoendothelial progenitors Intermediate mesoderm Mesenchyme
  cell_11995                         0.0396                0.1627     0.0869
  cell_10294                         0.0371                0.0387     0.0369
  cell_9963                          0.3906                0.1157     0.2338
  cell_11610                         0.0258                0.0444     0.0324
  cell_10910                         0.1477                0.2265     0.2042
  cell_11021                         0.0473                0.2016     0.0814
            gene sets
cells        Neural crest    NMP Paraxial mesoderm Pharyngeal mesoderm
  cell_11995        0.208 0.1805            0.1889              0.1844
  cell_10294        0.109 0.0445            0.0226              0.0263
  cell_9963         0.145 0.1045            0.1706              0.1561
  cell_11610        0.107 0.0640            0.0195              0.0294
  cell_10910        0.135 0.2025            0.1437              0.1686
  cell_11021        0.168 0.3432            0.1255              0.1515
            gene sets
cells        Somitic mesoderm Spinal cord Surface ectoderm
  cell_11995           0.1279      0.2592           0.1611
  cell_10294           0.0203      0.0633           0.0468
  cell_9963            0.1197      0.1014           0.0850
  cell_11610           0.0424      0.0692           0.0332
  cell_10910           0.1787      0.1521           0.1043
  cell_11021           0.2259      0.1974           0.1509
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
they should at least represent closely related cell states). We only print out the top-left corner of the table here, but you should try looking at the whole thing:


``` r
new.labels <- colnames(results)[max.col(results)]

tab <- table(new.labels, sce$celltype.mapped)

tab[1:4,1:4]
```

``` output
                
new.labels       Allantois Blood progenitors 1 Blood progenitors 2
  Allantois             44                   0                   0
  Cardiomyocytes         0                   0                   0
  Endothelium            0                   3                   0
  Erythroid2             0                   1                   7
                
new.labels       Cardiomyocytes
  Allantois                   0
  Cardiomyocytes             32
  Endothelium                 0
  Erythroid2                  0
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

Inspect the diagnostics for the next nine cell types. Do they look okay?

::: solution

``` r
par(mfrow = c(3,3))

AUCell_exploreThresholds(cell.aucs[10:18], plotHist = TRUE, assign = TRUE) 
```

<img src="fig/cell_type_annotation-rendered-unnamed-chunk-7-1.png" style="display: block; margin: auto;" />

``` output
$`Haematoendothelial progenitors`
$`Haematoendothelial progenitors`$aucThr
$`Haematoendothelial progenitors`$aucThr$selected
 R_k3 
0.273 

$`Haematoendothelial progenitors`$aucThr$thresholds
          threshold nCells
Global_k1     0.173    196
L_k2          0.199    138
R_k3          0.273     48

$`Haematoendothelial progenitors`$aucThr$comment
[1] "The global distribution overlaps the partial distributions. "


$`Haematoendothelial progenitors`$assignment
 [1] "cell_9963"  "cell_11842" "cell_9934"  "cell_11154" "cell_10322"
 [6] "cell_11337" "cell_11770" "cell_12021" "cell_11667" "cell_10047"
[11] "cell_11290" "cell_11733" "cell_11827" "cell_11404" "cell_11511"
[16] "cell_10146" "cell_10718" "cell_12090" "cell_10715" "cell_10458"
[21] "cell_10216" "cell_10010" "cell_11760" "cell_10126" "cell_10654"
[26] "cell_9813"  "cell_11916" "cell_12082" "cell_11769" "cell_11484"
[31] "cell_11414" "cell_11932" "cell_11791" "cell_9974"  "cell_11199"
[36] "cell_10256" "cell_10319" "cell_10977" "cell_10945" "cell_11077"
[41] "cell_10895" "cell_10168" "cell_9770"  "cell_11772" "cell_11771"
[46] "cell_11611" "cell_11287" "cell_10095"


$`Intermediate mesoderm`
$`Intermediate mesoderm`$aucThr
$`Intermediate mesoderm`$aucThr$selected
minimumDens 
     0.0723 

$`Intermediate mesoderm`$aucThr$thresholds
                threshold nCells
outlierOfGlobal    0.4388      0
Global_k1          0.1984    262
L_k2               0.2203    168
R_k3               0.1547    516
minimumDens        0.0723    841

$`Intermediate mesoderm`$aucThr$comment
[1] "The AUC might follow a normal distribution (random gene-set?). The global distribution overlaps the partial distributions. "


$`Intermediate mesoderm`$assignment
  [1] "cell_11995" "cell_9963"  "cell_10910" "cell_11021" "cell_10433"
  [6] "cell_11395" "cell_10779" "cell_10883" "cell_10721" "cell_10116"
 [11] "cell_10785" "cell_11781" "cell_11147" "cell_11842" "cell_11218"
 [16] "cell_11441" "cell_11558" "cell_9979"  "cell_10323" "cell_11663"
 [21] "cell_10935" "cell_11567" "cell_10673" "cell_11634" "cell_10815"
 [26] "cell_10077" "cell_9934"  "cell_11888" "cell_10356" "cell_10490"
 [31] "cell_9921"  "cell_10062" "cell_9809"  "cell_11906" "cell_10951"
 [36] "cell_10520" "cell_10142" "cell_12176" "cell_11850" "cell_10284"
 [41] "cell_9781"  "cell_11885" "cell_11687" "cell_10523" "cell_11100"
 [46] "cell_10881" "cell_10774" "cell_12107" "cell_11216" "cell_10439"
 [51] "cell_10902" "cell_10856" "cell_11702" "cell_11269" "cell_10248"
 [56] "cell_11883" "cell_11682" "cell_10877" "cell_9933"  "cell_10843"
 [61] "cell_11154" "cell_12052" "cell_9895"  "cell_12028" "cell_10454"
 [66] "cell_11750" "cell_12059" "cell_10615" "cell_10751" "cell_10387"
 [71] "cell_10245" "cell_11767" "cell_9919"  "cell_11976" "cell_11971"
 [76] "cell_10797" "cell_10094" "cell_11724" "cell_11030" "cell_11855"
 [81] "cell_10754" "cell_9905"  "cell_11530" "cell_11352" "cell_10357"
 [86] "cell_12074" "cell_12012" "cell_11561" "cell_11398" "cell_12075"
 [91] "cell_12177" "cell_11653" "cell_10644" "cell_11993" "cell_10322"
 [96] "cell_11900" "cell_11315" "cell_10474" "cell_10480" "cell_11220"
[101] "cell_10440" "cell_11596" "cell_10296" "cell_11337" "cell_11856"
[106] "cell_11728" "cell_10241" "cell_12016" "cell_11770" "cell_10917"
[111] "cell_11805" "cell_12002" "cell_10599" "cell_10846" "cell_10687"
[116] "cell_11330" "cell_11081" "cell_9953"  "cell_10181" "cell_10395"
[121] "cell_12021" "cell_11667" "cell_10547" "cell_11329" "cell_10332"
[126] "cell_10562" "cell_11183" "cell_11138" "cell_10385" "cell_10125"
[131] "cell_10047" "cell_10038" "cell_10386" "cell_11673" "cell_10105"
[136] "cell_12162" "cell_10749" "cell_11549" "cell_11458" "cell_11290"
[141] "cell_11213" "cell_10849" "cell_11465" "cell_9885"  "cell_10416"
[146] "cell_10847" "cell_11009" "cell_11733" "cell_10756" "cell_10479"
[151] "cell_11260" "cell_11002" "cell_10117" "cell_10933" "cell_11840"
[156] "cell_10922" "cell_11495" "cell_11380" "cell_10478" "cell_11026"
[161] "cell_11880" "cell_11896" "cell_12110" "cell_11115" "cell_10559"
[166] "cell_12151" "cell_11109" "cell_11014" "cell_11526" "cell_11018"
[171] "cell_11582" "cell_11476" "cell_10065" "cell_10628" "cell_11603"
[176] "cell_10825" "cell_9851"  "cell_10634" "cell_10791" "cell_9844" 
[181] "cell_10886" "cell_11931" "cell_11975" "cell_11150" "cell_10450"
[186] "cell_12115" "cell_11506" "cell_11239" "cell_11130" "cell_10011"
[191] "cell_11922" "cell_11827" "cell_10393" "cell_11982" "cell_10229"
[196] "cell_12140" "cell_9799"  "cell_10907" "cell_11957" "cell_9784" 
[201] "cell_11204" "cell_10292" "cell_11229" "cell_10692" "cell_11370"
[206] "cell_9972"  "cell_10914" "cell_10027" "cell_10262" "cell_11657"
[211] "cell_11552" "cell_11123" "cell_10900" "cell_11093" "cell_9778" 
[216] "cell_9938"  "cell_11001" "cell_11619" "cell_12160" "cell_9876" 
[221] "cell_9776"  "cell_10394" "cell_10029" "cell_11377" "cell_11679"
[226] "cell_11054" "cell_11528" "cell_11913" "cell_11515" "cell_10435"
[231] "cell_10758" "cell_10220" "cell_10624" "cell_11716" "cell_11371"
[236] "cell_12042" "cell_10414" "cell_10828" "cell_10659" "cell_11424"
[241] "cell_11096" "cell_12171" "cell_11803" "cell_11404" "cell_10561"
[246] "cell_11180" "cell_10046" "cell_10009" "cell_9792"  "cell_11853"
[251] "cell_11095" "cell_11511" "cell_11835" "cell_11463" "cell_10146"
[256] "cell_10317" "cell_11407" "cell_12173" "cell_10161" "cell_10438"
[261] "cell_10981" "cell_10079" "cell_11210" "cell_11174" "cell_10876"
[266] "cell_10087" "cell_11630" "cell_11434" "cell_9888"  "cell_10367"
[271] "cell_11107" "cell_11994" "cell_10192" "cell_11822" "cell_10647"
[276] "cell_10436" "cell_11592" "cell_9817"  "cell_11268" "cell_9961" 
[281] "cell_10071" "cell_9959"  "cell_11419" "cell_10585" "cell_11214"
[286] "cell_11084" "cell_11165" "cell_11327" "cell_11392" "cell_10718"
[291] "cell_10466" "cell_12090" "cell_11555" "cell_10995" "cell_10526"
[296] "cell_10761" "cell_11038" "cell_11232" "cell_10715" "cell_10458"
[301] "cell_10019" "cell_10328" "cell_10411" "cell_10313" "cell_9930" 
[306] "cell_11114" "cell_9936"  "cell_11234" "cell_10556" "cell_9846" 
[311] "cell_11871" "cell_11285" "cell_10991" "cell_11939" "cell_10213"
[316] "cell_9863"  "cell_11710" "cell_11256" "cell_10147" "cell_9989" 
[321] "cell_10464" "cell_10953" "cell_10388" "cell_10216" "cell_10010"
[326] "cell_11485" "cell_10695" "cell_11760" "cell_11328" "cell_11620"
[331] "cell_10694" "cell_10175" "cell_9997"  "cell_12147" "cell_10553"
[336] "cell_11266" "cell_10467" "cell_11233" "cell_10963" "cell_12149"
[341] "cell_9986"  "cell_11276" "cell_11904" "cell_12179" "cell_9847" 
[346] "cell_11303" "cell_11299" "cell_10956" "cell_10347" "cell_11599"
[351] "cell_10697" "cell_10190" "cell_10261" "cell_10498" "cell_11588"
[356] "cell_9977"  "cell_11839" "cell_11777" "cell_11941" "cell_10645"
[361] "cell_10126" "cell_10833" "cell_10672" "cell_10654" "cell_11024"
[366] "cell_10164" "cell_11522" "cell_10407" "cell_12134" "cell_11056"
[371] "cell_11489" "cell_9969"  "cell_9820"  "cell_11017" "cell_9835" 
[376] "cell_12064" "cell_11562" "cell_10345" "cell_10225" "cell_10783"
[381] "cell_12036" "cell_11333" "cell_10812" "cell_10998" "cell_9892" 
[386] "cell_10360" "cell_12152" "cell_10508" "cell_9813"  "cell_11124"
[391] "cell_11073" "cell_9859"  "cell_10421" "cell_10930" "cell_10374"
[396] "cell_12097" "cell_11217" "cell_10800" "cell_11063" "cell_10363"
[401] "cell_10515" "cell_10267" "cell_11047" "cell_11744" "cell_11564"
[406] "cell_9849"  "cell_11916" "cell_12082" "cell_12057" "cell_10664"
[411] "cell_11680" "cell_11542" "cell_10765" "cell_12066" "cell_11438"
[416] "cell_12072" "cell_10441" "cell_10396" "cell_10838" "cell_9832" 
[421] "cell_10299" "cell_11769" "cell_11402" "cell_11736" "cell_11012"
[426] "cell_11484" "cell_11962" "cell_10604" "cell_11966" "cell_10370"
[431] "cell_10099" "cell_10669" "cell_11774" "cell_11579" "cell_10059"
[436] "cell_10511" "cell_10597" "cell_12031" "cell_10961" "cell_10497"
[441] "cell_11934" "cell_11353" "cell_11133" "cell_9837"  "cell_11112"
[446] "cell_10272" "cell_10868" "cell_11734" "cell_10554" "cell_11683"
[451] "cell_10898" "cell_10623" "cell_10864" "cell_10217" "cell_10739"
[456] "cell_12174" "cell_12053" "cell_11859" "cell_10103" "cell_12167"
[461] "cell_10848" "cell_9975"  "cell_10204" "cell_11678" "cell_10769"
[466] "cell_10369" "cell_10155" "cell_10860" "cell_10512" "cell_10333"
[471] "cell_11149" "cell_10560" "cell_10510" "cell_11658" "cell_10354"
[476] "cell_11141" "cell_11263" "cell_12058" "cell_9856"  "cell_11208"
[481] "cell_9779"  "cell_11598" "cell_10231" "cell_11547" "cell_10468"
[486] "cell_11961" "cell_10610" "cell_10453" "cell_10048" "cell_11304"
[491] "cell_10571" "cell_11778" "cell_10546" "cell_12111" "cell_11120"
[496] "cell_9940"  "cell_12040" "cell_11471" "cell_11470" "cell_11321"
[501] "cell_10236" "cell_11176" "cell_12094" "cell_10509" "cell_10133"
[506] "cell_10359" "cell_11595" "cell_11019" "cell_10602" "cell_11414"
[511] "cell_11838" "cell_11932" "cell_11057" "cell_10303" "cell_10750"
[516] "cell_9880"  "cell_10224" "cell_10885" "cell_11281" "cell_10399"
[521] "cell_9816"  "cell_10737" "cell_10250" "cell_11040" "cell_9873" 
[526] "cell_11782" "cell_10227" "cell_12125" "cell_11714" "cell_11104"
[531] "cell_11354" "cell_10660" "cell_9907"  "cell_10249" "cell_10611"
[536] "cell_10018" "cell_11389" "cell_10098" "cell_10401" "cell_11418"
[541] "cell_12092" "cell_10890" "cell_10850" "cell_11791" "cell_11988"
[546] "cell_11897" "cell_10026" "cell_10403" "cell_10030" "cell_11688"
[551] "cell_11571" "cell_10696" "cell_9974"  "cell_9967"  "cell_11381"
[556] "cell_11972" "cell_12143" "cell_11837" "cell_10521" "cell_11280"
[561] "cell_11694" "cell_10166" "cell_11010" "cell_10947" "cell_11901"
[566] "cell_12044" "cell_10586" "cell_10428" "cell_10115" "cell_11357"
[571] "cell_11121" "cell_12114" "cell_12076" "cell_9823"  "cell_11162"
[576] "cell_10836" "cell_11954" "cell_11713" "cell_11168" "cell_10700"
[581] "cell_11845" "cell_11557" "cell_11346" "cell_9962"  "cell_11386"
[586] "cell_11799" "cell_10869" "cell_11411" "cell_10766" "cell_11450"
[591] "cell_11058" "cell_10481" "cell_10123" "cell_10238" "cell_10107"
[596] "cell_9968"  "cell_11951" "cell_11811" "cell_11519" "cell_11172"
[601] "cell_10661" "cell_10534" "cell_10904" "cell_12081" "cell_10237"
[606] "cell_10247" "cell_10380" "cell_11199" "cell_11447" "cell_11854"
[611] "cell_11524" "cell_10058" "cell_11586" "cell_10234" "cell_11743"
[616] "cell_10256" "cell_11775" "cell_11270" "cell_10089" "cell_9777" 
[621] "cell_10635" "cell_10049" "cell_10410" "cell_10675" "cell_10476"
[626] "cell_10832" "cell_10940" "cell_10788" "cell_11818" "cell_11502"
[631] "cell_10731" "cell_12166" "cell_10646" "cell_11397" "cell_10151"
[636] "cell_12085" "cell_11435" "cell_10239" "cell_11715" "cell_12060"
[641] "cell_11798" "cell_9929"  "cell_10189" "cell_10212" "cell_9797" 
[646] "cell_11529" "cell_11033" "cell_11794" "cell_9812"  "cell_9992" 
[651] "cell_11696" "cell_10680" "cell_10319" "cell_10120" "cell_10601"
[656] "cell_12145" "cell_10451" "cell_11243" "cell_9824"  "cell_11070"
[661] "cell_11983" "cell_11986" "cell_10422" "cell_11664" "cell_10143"
[666] "cell_11368" "cell_11580" "cell_11583" "cell_11312" "cell_11862"
[671] "cell_11493" "cell_11226" "cell_9785"  "cell_10417" "cell_9931" 
[676] "cell_11277" "cell_11261" "cell_11642" "cell_10977" "cell_10424"
[681] "cell_10588" "cell_9803"  "cell_11938" "cell_11877" "cell_10271"
[686] "cell_10177" "cell_10437" "cell_10741" "cell_10580" "cell_11102"
[691] "cell_10945" "cell_11077" "cell_9916"  "cell_9861"  "cell_10037"
[696] "cell_10485" "cell_10443" "cell_10140" "cell_10285" "cell_11618"
[701] "cell_9845"  "cell_11469" "cell_11140" "cell_10452" "cell_10352"
[706] "cell_9894"  "cell_11823" "cell_10944" "cell_10488" "cell_10867"
[711] "cell_11516" "cell_12116" "cell_10755" "cell_10072" "cell_11393"
[716] "cell_11650" "cell_10878" "cell_11400" "cell_9868"  "cell_10817"
[721] "cell_11700" "cell_11085" "cell_12007" "cell_11998" "cell_12158"
[726] "cell_10315" "cell_10698" "cell_10471" "cell_10111" "cell_11945"
[731] "cell_10999" "cell_11225" "cell_11083" "cell_11539" "cell_12113"
[736] "cell_11158" "cell_9802"  "cell_11173" "cell_10034" "cell_10193"
[741] "cell_12083" "cell_11296" "cell_11007" "cell_10895" "cell_10725"
[746] "cell_11709" "cell_10626" "cell_10168" "cell_9770"  "cell_10870"
[751] "cell_10950" "cell_9924"  "cell_12006" "cell_9983"  "cell_10787"
[756] "cell_10054" "cell_11118" "cell_10732" "cell_10198" "cell_11923"
[761] "cell_10518" "cell_10418" "cell_10484" "cell_10522" "cell_10153"
[766] "cell_11483" "cell_11309" "cell_11445" "cell_11754" "cell_11914"
[771] "cell_10677" "cell_11772" "cell_11697" "cell_9942"  "cell_10202"
[776] "cell_10383" "cell_11074" "cell_11212" "cell_10572" "cell_11342"
[781] "cell_10593" "cell_10283" "cell_12056" "cell_9819"  "cell_11771"
[786] "cell_11591" "cell_9998"  "cell_11144" "cell_11015" "cell_10381"
[791] "cell_10937" "cell_11265" "cell_10086" "cell_10187" "cell_10215"
[796] "cell_11611" "cell_9911"  "cell_11308" "cell_10735" "cell_11978"
[801] "cell_11067" "cell_10491" "cell_11003" "cell_10989" "cell_11480"
[806] "cell_10255" "cell_10312" "cell_11548" "cell_11207" "cell_9920" 
[811] "cell_10209" "cell_11211" "cell_11422" "cell_11360" "cell_10305"
[816] "cell_10821" "cell_10461" "cell_10008" "cell_10801" "cell_11319"
[821] "cell_10078" "cell_11651" "cell_11287" "cell_10389" "cell_10135"
[826] "cell_10773" "cell_11246" "cell_11623" "cell_11627" "cell_11031"
[831] "cell_10200" "cell_11113" "cell_9794"  "cell_12023" "cell_11197"
[836] "cell_11313" "cell_11129" "cell_11597" "cell_9807"  "cell_10095"
[841] "cell_11860"


$Mesenchyme
$Mesenchyme$aucThr
$Mesenchyme$aucThr$selected
minimumDens 
      0.418 

$Mesenchyme$aucThr$thresholds
            threshold nCells
Global_k1       0.261    168
L_k2            0.374     95
R_k3            0.171    346
minimumDens     0.418     81

$Mesenchyme$aucThr$comment
[1] ""


$Mesenchyme$assignment
 [1] "cell_10062" "cell_11100" "cell_12059" "cell_10387" "cell_9919" 
 [6] "cell_10754" "cell_12075" "cell_10474" "cell_10480" "cell_11081"
[11] "cell_10125" "cell_12162" "cell_11465" "cell_11840" "cell_9784" 
[16] "cell_11204" "cell_10900" "cell_9938"  "cell_10394" "cell_11377"
[21] "cell_11407" "cell_10367" "cell_9961"  "cell_9930"  "cell_11710"
[26] "cell_10388" "cell_12147" "cell_10963" "cell_11276" "cell_10498"
[31] "cell_10645" "cell_10164" "cell_11489" "cell_9969"  "cell_12064"
[36] "cell_10363" "cell_12072" "cell_10669" "cell_10511" "cell_11934"
[41] "cell_10864" "cell_10204" "cell_10155" "cell_11208" "cell_11120"
[46] "cell_9880"  "cell_10250" "cell_12092" "cell_11897" "cell_11571"
[51] "cell_11837" "cell_11010" "cell_10947" "cell_12044" "cell_11954"
[56] "cell_9962"  "cell_10237" "cell_11502" "cell_10151" "cell_11435"
[61] "cell_12060" "cell_10451" "cell_11243" "cell_11070" "cell_11983"
[66] "cell_11580" "cell_11583" "cell_9785"  "cell_10177" "cell_10485"
[71] "cell_11469" "cell_12113" "cell_10677" "cell_10383" "cell_11015"
[76] "cell_10491" "cell_11480" "cell_10461" "cell_10389" "cell_9794" 
[81] "cell_9807" 


$`Neural crest`
$`Neural crest`$aucThr
$`Neural crest`$aucThr$selected
outlierOfGlobal 
          0.331 

$`Neural crest`$aucThr$thresholds
                threshold nCells
outlierOfGlobal     0.331      0
Global_k1           0.181    264
L_k2                0.168    373
R_k3                0.254     12

$`Neural crest`$aucThr$comment
[1] "The AUC might follow a normal distribution (random gene-set?). The right distribution is taller. "


$`Neural crest`$assignment
character(0)


$NMP
$NMP$aucThr
$NMP$aucThr$selected
 L_k2 
0.314 

$NMP$aucThr$thresholds
          threshold nCells
Global_k1    0.1873    195
L_k2         0.3143     42
R_k3         0.0961    677

$NMP$aucThr$comment
[1] ""


$NMP$assignment
 [1] "cell_11021" "cell_10721" "cell_10116" "cell_10248" "cell_12052"
 [6] "cell_10440" "cell_11380" "cell_9851"  "cell_10886" "cell_11239"
[11] "cell_11180" "cell_10647" "cell_10411" "cell_11234" "cell_12149"
[16] "cell_9977"  "cell_9820"  "cell_12097" "cell_10597" "cell_10868"
[21] "cell_11149" "cell_10231" "cell_10453" "cell_11040" "cell_11104"
[26] "cell_10696" "cell_10380" "cell_10234" "cell_10832" "cell_9929" 
[31] "cell_9824"  "cell_11368" "cell_10140" "cell_11393" "cell_11998"
[36] "cell_10698" "cell_11083" "cell_12006" "cell_11207" "cell_10821"
[41] "cell_11319" "cell_11597"


$`Paraxial mesoderm`
$`Paraxial mesoderm`$aucThr
$`Paraxial mesoderm`$aucThr$selected
minimumDens 
     0.0736 

$`Paraxial mesoderm`$aucThr$thresholds
            threshold nCells
Global_k1      0.1970    181
L_k2           0.5372      0
R_k3           0.0667    841
minimumDens    0.0736    840

$`Paraxial mesoderm`$aucThr$comment
[1] "The global distribution overlaps the partial distributions. "


$`Paraxial mesoderm`$assignment
  [1] "cell_11995" "cell_9963"  "cell_10910" "cell_11021" "cell_10433"
  [6] "cell_11395" "cell_10779" "cell_10883" "cell_10721" "cell_10116"
 [11] "cell_10785" "cell_11781" "cell_11147" "cell_11842" "cell_11218"
 [16] "cell_11441" "cell_11558" "cell_9979"  "cell_10323" "cell_11663"
 [21] "cell_10935" "cell_11567" "cell_10673" "cell_11634" "cell_10815"
 [26] "cell_10077" "cell_9934"  "cell_11888" "cell_10356" "cell_10490"
 [31] "cell_9921"  "cell_10062" "cell_9809"  "cell_11906" "cell_10951"
 [36] "cell_10520" "cell_12176" "cell_11850" "cell_10284" "cell_9781" 
 [41] "cell_11885" "cell_11687" "cell_10523" "cell_11100" "cell_10881"
 [46] "cell_10774" "cell_12107" "cell_11216" "cell_10439" "cell_10902"
 [51] "cell_10856" "cell_11702" "cell_11269" "cell_10248" "cell_11883"
 [56] "cell_11682" "cell_10877" "cell_9933"  "cell_10843" "cell_11154"
 [61] "cell_12052" "cell_9895"  "cell_12028" "cell_10454" "cell_11750"
 [66] "cell_12059" "cell_10615" "cell_10751" "cell_10387" "cell_10245"
 [71] "cell_11767" "cell_9919"  "cell_11976" "cell_11971" "cell_10797"
 [76] "cell_10094" "cell_11724" "cell_11030" "cell_11855" "cell_10754"
 [81] "cell_9905"  "cell_11530" "cell_11352" "cell_10357" "cell_12074"
 [86] "cell_12012" "cell_11561" "cell_11398" "cell_12075" "cell_12177"
 [91] "cell_11653" "cell_10644" "cell_11993" "cell_10322" "cell_11900"
 [96] "cell_11315" "cell_10474" "cell_10480" "cell_11220" "cell_10440"
[101] "cell_11596" "cell_10296" "cell_11337" "cell_11856" "cell_11728"
[106] "cell_10241" "cell_12016" "cell_11770" "cell_10917" "cell_11805"
[111] "cell_12002" "cell_10599" "cell_10846" "cell_10687" "cell_11330"
[116] "cell_11081" "cell_9953"  "cell_10181" "cell_10395" "cell_12021"
[121] "cell_11667" "cell_10547" "cell_11329" "cell_10332" "cell_10562"
[126] "cell_11183" "cell_11138" "cell_10385" "cell_10125" "cell_10047"
[131] "cell_10038" "cell_10386" "cell_11673" "cell_10105" "cell_12162"
[136] "cell_10749" "cell_11549" "cell_11458" "cell_11290" "cell_11213"
[141] "cell_10849" "cell_11465" "cell_9885"  "cell_10416" "cell_10847"
[146] "cell_11009" "cell_11733" "cell_10756" "cell_10479" "cell_11260"
[151] "cell_11002" "cell_10117" "cell_10933" "cell_11840" "cell_10922"
[156] "cell_11495" "cell_11380" "cell_10478" "cell_11026" "cell_11880"
[161] "cell_11896" "cell_12110" "cell_11115" "cell_10559" "cell_12151"
[166] "cell_11109" "cell_11014" "cell_11526" "cell_11018" "cell_11582"
[171] "cell_11476" "cell_10065" "cell_10628" "cell_11603" "cell_10825"
[176] "cell_9851"  "cell_10634" "cell_10791" "cell_9844"  "cell_10886"
[181] "cell_11931" "cell_11975" "cell_11150" "cell_10450" "cell_12115"
[186] "cell_11506" "cell_11239" "cell_11130" "cell_10011" "cell_11922"
[191] "cell_11827" "cell_10393" "cell_11982" "cell_10229" "cell_12140"
[196] "cell_9799"  "cell_10907" "cell_11957" "cell_9784"  "cell_11204"
[201] "cell_10292" "cell_11229" "cell_10692" "cell_11370" "cell_9972" 
[206] "cell_10914" "cell_10027" "cell_10262" "cell_11657" "cell_11552"
[211] "cell_12029" "cell_11123" "cell_10900" "cell_11093" "cell_9778" 
[216] "cell_9938"  "cell_11001" "cell_11619" "cell_12160" "cell_9876" 
[221] "cell_9776"  "cell_10394" "cell_10029" "cell_11377" "cell_11679"
[226] "cell_11054" "cell_11528" "cell_11913" "cell_11515" "cell_10435"
[231] "cell_10758" "cell_10220" "cell_10624" "cell_11716" "cell_11371"
[236] "cell_12042" "cell_10414" "cell_10828" "cell_10659" "cell_11424"
[241] "cell_11096" "cell_12171" "cell_11803" "cell_11404" "cell_11180"
[246] "cell_10046" "cell_10009" "cell_9792"  "cell_11853" "cell_11095"
[251] "cell_11511" "cell_11835" "cell_11463" "cell_10146" "cell_10317"
[256] "cell_11407" "cell_12173" "cell_10161" "cell_10438" "cell_10981"
[261] "cell_10079" "cell_11210" "cell_11174" "cell_10876" "cell_10087"
[266] "cell_11630" "cell_11434" "cell_9888"  "cell_10367" "cell_11107"
[271] "cell_11994" "cell_10192" "cell_11822" "cell_10647" "cell_10436"
[276] "cell_11592" "cell_9817"  "cell_11268" "cell_9961"  "cell_10071"
[281] "cell_9959"  "cell_11419" "cell_10585" "cell_11214" "cell_11084"
[286] "cell_11165" "cell_11327" "cell_11392" "cell_10718" "cell_10466"
[291] "cell_12090" "cell_11555" "cell_10995" "cell_10526" "cell_10761"
[296] "cell_11038" "cell_11232" "cell_10715" "cell_10458" "cell_10019"
[301] "cell_10328" "cell_10411" "cell_10313" "cell_9930"  "cell_11114"
[306] "cell_9936"  "cell_11234" "cell_10556" "cell_9846"  "cell_11871"
[311] "cell_11285" "cell_10991" "cell_11939" "cell_10213" "cell_9863" 
[316] "cell_11710" "cell_11256" "cell_10147" "cell_9989"  "cell_10464"
[321] "cell_10953" "cell_10388" "cell_10216" "cell_10010" "cell_11485"
[326] "cell_10695" "cell_11760" "cell_11328" "cell_11620" "cell_10694"
[331] "cell_10175" "cell_9997"  "cell_12147" "cell_10553" "cell_11266"
[336] "cell_10467" "cell_11233" "cell_10963" "cell_12149" "cell_9986" 
[341] "cell_11276" "cell_11904" "cell_12179" "cell_9847"  "cell_11303"
[346] "cell_11299" "cell_10956" "cell_10347" "cell_11599" "cell_10697"
[351] "cell_10190" "cell_10261" "cell_10498" "cell_11588" "cell_9977" 
[356] "cell_11839" "cell_11777" "cell_11941" "cell_10645" "cell_10126"
[361] "cell_10833" "cell_10672" "cell_10654" "cell_11024" "cell_10164"
[366] "cell_11522" "cell_10407" "cell_12134" "cell_11056" "cell_11489"
[371] "cell_9969"  "cell_9820"  "cell_11017" "cell_9835"  "cell_12064"
[376] "cell_11562" "cell_10345" "cell_10225" "cell_10783" "cell_12036"
[381] "cell_11333" "cell_10812" "cell_10998" "cell_9892"  "cell_10360"
[386] "cell_12152" "cell_10508" "cell_9813"  "cell_11124" "cell_11073"
[391] "cell_9859"  "cell_10421" "cell_10930" "cell_10374" "cell_12097"
[396] "cell_11217" "cell_10800" "cell_11063" "cell_10363" "cell_10515"
[401] "cell_10267" "cell_11047" "cell_11744" "cell_11564" "cell_9849" 
[406] "cell_11916" "cell_12082" "cell_12057" "cell_10664" "cell_11680"
[411] "cell_11542" "cell_10765" "cell_12066" "cell_11438" "cell_12072"
[416] "cell_10441" "cell_10396" "cell_10838" "cell_9832"  "cell_10299"
[421] "cell_11769" "cell_11402" "cell_11736" "cell_11012" "cell_11484"
[426] "cell_11962" "cell_10604" "cell_11966" "cell_10370" "cell_10099"
[431] "cell_10669" "cell_11774" "cell_11579" "cell_10059" "cell_10511"
[436] "cell_10597" "cell_12031" "cell_10961" "cell_10497" "cell_11934"
[441] "cell_11353" "cell_11133" "cell_9837"  "cell_11112" "cell_10272"
[446] "cell_10868" "cell_11734" "cell_10554" "cell_11683" "cell_10898"
[451] "cell_10623" "cell_10864" "cell_10217" "cell_10739" "cell_12174"
[456] "cell_12053" "cell_11859" "cell_10103" "cell_12167" "cell_10848"
[461] "cell_9975"  "cell_10204" "cell_11678" "cell_10769" "cell_10369"
[466] "cell_10155" "cell_10860" "cell_10512" "cell_10333" "cell_11149"
[471] "cell_10560" "cell_10510" "cell_11658" "cell_10354" "cell_11141"
[476] "cell_11263" "cell_12058" "cell_9856"  "cell_11208" "cell_9779" 
[481] "cell_11598" "cell_10231" "cell_11547" "cell_10468" "cell_11961"
[486] "cell_10610" "cell_10453" "cell_10048" "cell_11304" "cell_10571"
[491] "cell_11778" "cell_10546" "cell_12111" "cell_11120" "cell_9940" 
[496] "cell_12040" "cell_11471" "cell_11470" "cell_11321" "cell_10236"
[501] "cell_11176" "cell_12094" "cell_10509" "cell_10133" "cell_10359"
[506] "cell_11595" "cell_11019" "cell_10602" "cell_11414" "cell_11838"
[511] "cell_11932" "cell_11057" "cell_10303" "cell_10750" "cell_9880" 
[516] "cell_10224" "cell_10885" "cell_11281" "cell_10399" "cell_9816" 
[521] "cell_10737" "cell_10250" "cell_11040" "cell_9873"  "cell_11782"
[526] "cell_10227" "cell_12125" "cell_11714" "cell_11104" "cell_11354"
[531] "cell_10660" "cell_9907"  "cell_10249" "cell_10611" "cell_10018"
[536] "cell_11389" "cell_10098" "cell_10401" "cell_11418" "cell_12092"
[541] "cell_10890" "cell_10850" "cell_11791" "cell_11988" "cell_11897"
[546] "cell_10026" "cell_10403" "cell_10030" "cell_11688" "cell_11571"
[551] "cell_10696" "cell_9974"  "cell_9967"  "cell_11381" "cell_11972"
[556] "cell_12143" "cell_11837" "cell_10521" "cell_11280" "cell_11694"
[561] "cell_10166" "cell_11010" "cell_10947" "cell_11901" "cell_12044"
[566] "cell_10586" "cell_10428" "cell_10115" "cell_11357" "cell_11121"
[571] "cell_12114" "cell_12076" "cell_9823"  "cell_11162" "cell_10836"
[576] "cell_11954" "cell_11713" "cell_11168" "cell_10700" "cell_11845"
[581] "cell_11557" "cell_11346" "cell_9962"  "cell_11386" "cell_11799"
[586] "cell_10869" "cell_11411" "cell_10766" "cell_11450" "cell_11058"
[591] "cell_10481" "cell_10123" "cell_10238" "cell_10107" "cell_9968" 
[596] "cell_11951" "cell_11811" "cell_11519" "cell_11172" "cell_10661"
[601] "cell_10534" "cell_10904" "cell_12081" "cell_10237" "cell_10247"
[606] "cell_10380" "cell_11199" "cell_11447" "cell_11854" "cell_11524"
[611] "cell_10058" "cell_11586" "cell_10234" "cell_11743" "cell_10256"
[616] "cell_11775" "cell_11270" "cell_10089" "cell_9777"  "cell_10635"
[621] "cell_10049" "cell_10410" "cell_10675" "cell_10476" "cell_10832"
[626] "cell_10940" "cell_10788" "cell_11818" "cell_11502" "cell_10731"
[631] "cell_12166" "cell_10646" "cell_11397" "cell_10151" "cell_12085"
[636] "cell_11435" "cell_10239" "cell_11715" "cell_12060" "cell_11798"
[641] "cell_9929"  "cell_10189" "cell_10212" "cell_9797"  "cell_11529"
[646] "cell_11033" "cell_11794" "cell_9812"  "cell_9992"  "cell_11696"
[651] "cell_10680" "cell_10319" "cell_10120" "cell_10601" "cell_12145"
[656] "cell_10451" "cell_11243" "cell_9824"  "cell_11070" "cell_11983"
[661] "cell_11986" "cell_10422" "cell_11664" "cell_10143" "cell_11368"
[666] "cell_11580" "cell_11583" "cell_11312" "cell_11862" "cell_11493"
[671] "cell_11226" "cell_9785"  "cell_10417" "cell_9931"  "cell_11277"
[676] "cell_11261" "cell_11642" "cell_10977" "cell_10424" "cell_10588"
[681] "cell_9803"  "cell_11938" "cell_11877" "cell_10271" "cell_10177"
[686] "cell_10437" "cell_10741" "cell_10580" "cell_11102" "cell_10945"
[691] "cell_11077" "cell_9916"  "cell_9861"  "cell_10037" "cell_10485"
[696] "cell_10443" "cell_10140" "cell_10285" "cell_11618" "cell_9845" 
[701] "cell_11469" "cell_11140" "cell_10452" "cell_10352" "cell_9894" 
[706] "cell_11823" "cell_10944" "cell_10488" "cell_10867" "cell_11516"
[711] "cell_12116" "cell_10755" "cell_10072" "cell_11393" "cell_11650"
[716] "cell_10878" "cell_11400" "cell_9868"  "cell_10817" "cell_11700"
[721] "cell_11085" "cell_12007" "cell_11998" "cell_12158" "cell_10315"
[726] "cell_10698" "cell_10471" "cell_10111" "cell_11945" "cell_10999"
[731] "cell_11225" "cell_11083" "cell_11539" "cell_12113" "cell_11158"
[736] "cell_9802"  "cell_11173" "cell_10034" "cell_10193" "cell_12083"
[741] "cell_11296" "cell_11007" "cell_10895" "cell_10725" "cell_11709"
[746] "cell_10626" "cell_10168" "cell_9770"  "cell_10870" "cell_10950"
[751] "cell_9924"  "cell_12006" "cell_9983"  "cell_10787" "cell_10054"
[756] "cell_11118" "cell_10732" "cell_10198" "cell_11923" "cell_10518"
[761] "cell_10418" "cell_10484" "cell_10522" "cell_10153" "cell_11483"
[766] "cell_11309" "cell_11445" "cell_11754" "cell_11914" "cell_10677"
[771] "cell_11772" "cell_11697" "cell_9942"  "cell_10202" "cell_10383"
[776] "cell_11074" "cell_11212" "cell_10572" "cell_11342" "cell_10593"
[781] "cell_10283" "cell_12056" "cell_9819"  "cell_11771" "cell_11591"
[786] "cell_9998"  "cell_11144" "cell_11015" "cell_10381" "cell_10937"
[791] "cell_11265" "cell_10086" "cell_10187" "cell_10215" "cell_11611"
[796] "cell_9911"  "cell_11308" "cell_10735" "cell_11978" "cell_11067"
[801] "cell_10491" "cell_11003" "cell_10989" "cell_11480" "cell_10255"
[806] "cell_10312" "cell_11548" "cell_11207" "cell_9920"  "cell_10209"
[811] "cell_11211" "cell_11422" "cell_11360" "cell_10305" "cell_10821"
[816] "cell_10461" "cell_10008" "cell_10801" "cell_11319" "cell_10078"
[821] "cell_11651" "cell_11287" "cell_10389" "cell_10135" "cell_10773"
[826] "cell_11246" "cell_11623" "cell_11627" "cell_11031" "cell_10200"
[831] "cell_11113" "cell_9794"  "cell_12023" "cell_11197" "cell_11313"
[836] "cell_11129" "cell_11597" "cell_9807"  "cell_10095" "cell_11860"


$`Pharyngeal mesoderm`
$`Pharyngeal mesoderm`$aucThr
$`Pharyngeal mesoderm`$aucThr$selected
minimumDens 
     0.0816 

$`Pharyngeal mesoderm`$aucThr$thresholds
            threshold nCells
Global_k1      0.2244    228
L_k2           0.5680      0
R_k3           0.0893    837
minimumDens    0.0816    838

$`Pharyngeal mesoderm`$aucThr$comment
[1] "The global distribution overlaps the partial distributions. "


$`Pharyngeal mesoderm`$assignment
  [1] "cell_11995" "cell_9963"  "cell_10910" "cell_11021" "cell_10433"
  [6] "cell_11395" "cell_10779" "cell_10883" "cell_10721" "cell_10116"
 [11] "cell_10785" "cell_11781" "cell_11147" "cell_11842" "cell_11218"
 [16] "cell_11441" "cell_11558" "cell_9979"  "cell_10323" "cell_11663"
 [21] "cell_10935" "cell_11567" "cell_10673" "cell_11634" "cell_10815"
 [26] "cell_10077" "cell_9934"  "cell_11888" "cell_10356" "cell_10490"
 [31] "cell_9921"  "cell_10062" "cell_9809"  "cell_11906" "cell_10951"
 [36] "cell_10520" "cell_12176" "cell_11850" "cell_10284" "cell_9781" 
 [41] "cell_11885" "cell_11687" "cell_10523" "cell_11100" "cell_10881"
 [46] "cell_10774" "cell_12107" "cell_11216" "cell_10439" "cell_10902"
 [51] "cell_10856" "cell_11702" "cell_11269" "cell_10248" "cell_11883"
 [56] "cell_11682" "cell_10877" "cell_9933"  "cell_10843" "cell_11154"
 [61] "cell_12052" "cell_9895"  "cell_12028" "cell_10454" "cell_11750"
 [66] "cell_12059" "cell_10615" "cell_10751" "cell_10387" "cell_10245"
 [71] "cell_11767" "cell_9919"  "cell_11976" "cell_11971" "cell_10797"
 [76] "cell_10094" "cell_11724" "cell_11030" "cell_11855" "cell_10754"
 [81] "cell_9905"  "cell_11530" "cell_11352" "cell_10357" "cell_12074"
 [86] "cell_12012" "cell_11561" "cell_11398" "cell_12075" "cell_12177"
 [91] "cell_11653" "cell_10644" "cell_11993" "cell_10322" "cell_11900"
 [96] "cell_11315" "cell_10474" "cell_10480" "cell_11220" "cell_10440"
[101] "cell_11596" "cell_10296" "cell_11337" "cell_11856" "cell_11728"
[106] "cell_10241" "cell_12016" "cell_11770" "cell_10917" "cell_11805"
[111] "cell_12002" "cell_10599" "cell_10846" "cell_10687" "cell_11330"
[116] "cell_11081" "cell_9953"  "cell_10181" "cell_10395" "cell_12021"
[121] "cell_11667" "cell_10547" "cell_11329" "cell_10332" "cell_10562"
[126] "cell_11183" "cell_11138" "cell_10385" "cell_10125" "cell_10047"
[131] "cell_10038" "cell_10386" "cell_11673" "cell_10105" "cell_12162"
[136] "cell_10749" "cell_11549" "cell_11458" "cell_11290" "cell_11213"
[141] "cell_10849" "cell_11465" "cell_9885"  "cell_10416" "cell_10847"
[146] "cell_11009" "cell_11733" "cell_10756" "cell_10479" "cell_11260"
[151] "cell_11002" "cell_10117" "cell_10933" "cell_11840" "cell_10922"
[156] "cell_11495" "cell_11380" "cell_10478" "cell_11026" "cell_11880"
[161] "cell_11896" "cell_12110" "cell_11115" "cell_10559" "cell_12151"
[166] "cell_11109" "cell_11014" "cell_11526" "cell_11018" "cell_11582"
[171] "cell_11476" "cell_10065" "cell_10628" "cell_11603" "cell_10825"
[176] "cell_9851"  "cell_10634" "cell_10791" "cell_9844"  "cell_10886"
[181] "cell_11931" "cell_11975" "cell_11150" "cell_10450" "cell_12115"
[186] "cell_11506" "cell_11239" "cell_11130" "cell_10011" "cell_11922"
[191] "cell_11827" "cell_10393" "cell_11982" "cell_10229" "cell_12140"
[196] "cell_9799"  "cell_10907" "cell_11957" "cell_9784"  "cell_11204"
[201] "cell_10292" "cell_11229" "cell_10692" "cell_11370" "cell_9972" 
[206] "cell_10914" "cell_10027" "cell_10262" "cell_11657" "cell_11552"
[211] "cell_11123" "cell_10900" "cell_11093" "cell_9778"  "cell_9938" 
[216] "cell_11001" "cell_11619" "cell_12160" "cell_9876"  "cell_9776" 
[221] "cell_10394" "cell_10029" "cell_11377" "cell_11679" "cell_11054"
[226] "cell_11528" "cell_11913" "cell_11515" "cell_10435" "cell_10758"
[231] "cell_10220" "cell_10624" "cell_11716" "cell_11371" "cell_12042"
[236] "cell_10414" "cell_10828" "cell_10659" "cell_11424" "cell_11096"
[241] "cell_12171" "cell_11803" "cell_11404" "cell_11180" "cell_10046"
[246] "cell_10009" "cell_9792"  "cell_11853" "cell_11095" "cell_11511"
[251] "cell_11835" "cell_11463" "cell_10146" "cell_10317" "cell_11407"
[256] "cell_12173" "cell_10161" "cell_10438" "cell_10981" "cell_10079"
[261] "cell_11210" "cell_11174" "cell_10876" "cell_10087" "cell_11630"
[266] "cell_11434" "cell_9888"  "cell_10367" "cell_11107" "cell_11994"
[271] "cell_10192" "cell_11822" "cell_10647" "cell_10436" "cell_11592"
[276] "cell_9817"  "cell_11268" "cell_9961"  "cell_10071" "cell_9959" 
[281] "cell_11419" "cell_10585" "cell_11214" "cell_11084" "cell_11165"
[286] "cell_11327" "cell_11392" "cell_10718" "cell_10466" "cell_12090"
[291] "cell_11555" "cell_10995" "cell_10526" "cell_10761" "cell_11038"
[296] "cell_11232" "cell_10715" "cell_10458" "cell_10019" "cell_10328"
[301] "cell_10411" "cell_10313" "cell_9930"  "cell_11114" "cell_9936" 
[306] "cell_11234" "cell_10556" "cell_9846"  "cell_11871" "cell_11285"
[311] "cell_10991" "cell_11939" "cell_10213" "cell_9863"  "cell_11710"
[316] "cell_11256" "cell_10147" "cell_9989"  "cell_10464" "cell_10953"
[321] "cell_10388" "cell_10216" "cell_10010" "cell_11485" "cell_10695"
[326] "cell_11760" "cell_11328" "cell_11620" "cell_10694" "cell_10175"
[331] "cell_9997"  "cell_12147" "cell_10553" "cell_11266" "cell_10467"
[336] "cell_11233" "cell_10963" "cell_12149" "cell_9986"  "cell_11276"
[341] "cell_11904" "cell_12179" "cell_9847"  "cell_11303" "cell_11299"
[346] "cell_10956" "cell_10347" "cell_11599" "cell_10697" "cell_10190"
[351] "cell_10261" "cell_10498" "cell_11588" "cell_9977"  "cell_11839"
[356] "cell_11777" "cell_11941" "cell_10645" "cell_10126" "cell_10833"
[361] "cell_10672" "cell_10654" "cell_11024" "cell_10164" "cell_11522"
[366] "cell_10407" "cell_12134" "cell_11056" "cell_11489" "cell_9969" 
[371] "cell_9820"  "cell_11017" "cell_9835"  "cell_12064" "cell_11562"
[376] "cell_10345" "cell_10225" "cell_10783" "cell_12036" "cell_11333"
[381] "cell_10812" "cell_10998" "cell_9892"  "cell_10360" "cell_12152"
[386] "cell_10508" "cell_9813"  "cell_11124" "cell_11073" "cell_9859" 
[391] "cell_10421" "cell_10930" "cell_10374" "cell_12097" "cell_11217"
[396] "cell_10800" "cell_11063" "cell_10363" "cell_10515" "cell_10267"
[401] "cell_11047" "cell_11744" "cell_11564" "cell_9849"  "cell_11916"
[406] "cell_12082" "cell_12057" "cell_10664" "cell_11680" "cell_11542"
[411] "cell_10765" "cell_12066" "cell_11438" "cell_12072" "cell_10441"
[416] "cell_10396" "cell_10838" "cell_9832"  "cell_10299" "cell_11769"
[421] "cell_11402" "cell_11736" "cell_11012" "cell_11484" "cell_11962"
[426] "cell_10604" "cell_11966" "cell_10370" "cell_10099" "cell_10669"
[431] "cell_11774" "cell_11579" "cell_10059" "cell_10511" "cell_10597"
[436] "cell_12031" "cell_10961" "cell_10497" "cell_11934" "cell_11353"
[441] "cell_11133" "cell_9837"  "cell_11112" "cell_10272" "cell_10868"
[446] "cell_11734" "cell_10554" "cell_11683" "cell_10898" "cell_10623"
[451] "cell_10864" "cell_10217" "cell_10739" "cell_12174" "cell_12053"
[456] "cell_11859" "cell_10103" "cell_12167" "cell_10848" "cell_9975" 
[461] "cell_10204" "cell_11678" "cell_10769" "cell_10369" "cell_10155"
[466] "cell_10860" "cell_10512" "cell_10333" "cell_11149" "cell_10560"
[471] "cell_10510" "cell_11658" "cell_10354" "cell_11141" "cell_11263"
[476] "cell_12058" "cell_9856"  "cell_11208" "cell_9779"  "cell_11598"
[481] "cell_10231" "cell_11547" "cell_10468" "cell_11961" "cell_10610"
[486] "cell_10453" "cell_10048" "cell_11304" "cell_10571" "cell_11778"
[491] "cell_10546" "cell_12111" "cell_11120" "cell_9940"  "cell_12040"
[496] "cell_11471" "cell_11470" "cell_11321" "cell_10236" "cell_11176"
[501] "cell_12094" "cell_10509" "cell_10133" "cell_10359" "cell_11595"
[506] "cell_11019" "cell_10602" "cell_11414" "cell_11838" "cell_11932"
[511] "cell_11057" "cell_10303" "cell_10750" "cell_9880"  "cell_10224"
[516] "cell_10885" "cell_11281" "cell_10399" "cell_9816"  "cell_10737"
[521] "cell_10250" "cell_11040" "cell_9873"  "cell_11782" "cell_10227"
[526] "cell_12125" "cell_11714" "cell_11104" "cell_11354" "cell_10660"
[531] "cell_9907"  "cell_10249" "cell_10611" "cell_10018" "cell_11389"
[536] "cell_10098" "cell_10401" "cell_11418" "cell_12092" "cell_10890"
[541] "cell_10850" "cell_11791" "cell_11988" "cell_11897" "cell_10026"
[546] "cell_10403" "cell_10030" "cell_11688" "cell_11571" "cell_10696"
[551] "cell_9974"  "cell_9967"  "cell_11381" "cell_11972" "cell_12143"
[556] "cell_11837" "cell_10521" "cell_11280" "cell_11694" "cell_10166"
[561] "cell_11010" "cell_10947" "cell_11901" "cell_12044" "cell_10586"
[566] "cell_10428" "cell_10115" "cell_11357" "cell_11121" "cell_12114"
[571] "cell_12076" "cell_9823"  "cell_11162" "cell_10836" "cell_11954"
[576] "cell_11713" "cell_11168" "cell_10700" "cell_11845" "cell_11557"
[581] "cell_11346" "cell_9962"  "cell_11799" "cell_10869" "cell_11411"
[586] "cell_10766" "cell_11450" "cell_11058" "cell_10481" "cell_10123"
[591] "cell_10238" "cell_10107" "cell_9968"  "cell_11951" "cell_11811"
[596] "cell_11519" "cell_11172" "cell_10661" "cell_10534" "cell_10904"
[601] "cell_12081" "cell_10237" "cell_10247" "cell_10380" "cell_11199"
[606] "cell_11447" "cell_11854" "cell_11524" "cell_10058" "cell_11586"
[611] "cell_10234" "cell_11743" "cell_10256" "cell_11775" "cell_11270"
[616] "cell_10089" "cell_9777"  "cell_10635" "cell_10049" "cell_10410"
[621] "cell_10675" "cell_10476" "cell_10832" "cell_10940" "cell_10788"
[626] "cell_11818" "cell_11502" "cell_10731" "cell_12166" "cell_10646"
[631] "cell_11397" "cell_10151" "cell_12085" "cell_11435" "cell_10239"
[636] "cell_11715" "cell_12060" "cell_11798" "cell_9929"  "cell_10189"
[641] "cell_10212" "cell_9797"  "cell_11529" "cell_11033" "cell_11794"
[646] "cell_9812"  "cell_9992"  "cell_11696" "cell_10680" "cell_10319"
[651] "cell_10120" "cell_10601" "cell_12145" "cell_10451" "cell_11243"
[656] "cell_9824"  "cell_11070" "cell_11983" "cell_11986" "cell_10422"
[661] "cell_11664" "cell_10143" "cell_11368" "cell_11580" "cell_11583"
[666] "cell_11312" "cell_11862" "cell_11493" "cell_11226" "cell_9785" 
[671] "cell_10417" "cell_9931"  "cell_11277" "cell_11261" "cell_11642"
[676] "cell_10977" "cell_10424" "cell_10588" "cell_9803"  "cell_11938"
[681] "cell_11877" "cell_10271" "cell_10177" "cell_10437" "cell_10741"
[686] "cell_10580" "cell_11102" "cell_10945" "cell_11077" "cell_9916" 
[691] "cell_9861"  "cell_10037" "cell_10485" "cell_10443" "cell_10140"
[696] "cell_10285" "cell_11618" "cell_9845"  "cell_11469" "cell_11140"
[701] "cell_10452" "cell_10352" "cell_9894"  "cell_11823" "cell_10944"
[706] "cell_10488" "cell_10867" "cell_11516" "cell_12116" "cell_10755"
[711] "cell_10072" "cell_11393" "cell_11650" "cell_10878" "cell_11400"
[716] "cell_9868"  "cell_10817" "cell_11700" "cell_11085" "cell_12007"
[721] "cell_11998" "cell_12158" "cell_10315" "cell_10698" "cell_10471"
[726] "cell_10111" "cell_11945" "cell_10999" "cell_11225" "cell_11083"
[731] "cell_11539" "cell_12113" "cell_11158" "cell_9802"  "cell_11173"
[736] "cell_10034" "cell_10193" "cell_12083" "cell_11296" "cell_11007"
[741] "cell_10895" "cell_10725" "cell_11709" "cell_10626" "cell_10168"
[746] "cell_9770"  "cell_10870" "cell_10950" "cell_9924"  "cell_12006"
[751] "cell_9983"  "cell_10787" "cell_10054" "cell_11118" "cell_10732"
[756] "cell_10198" "cell_11923" "cell_10518" "cell_10418" "cell_10484"
[761] "cell_10522" "cell_10153" "cell_11483" "cell_11309" "cell_11445"
[766] "cell_11754" "cell_11914" "cell_10677" "cell_11772" "cell_11697"
[771] "cell_9942"  "cell_10202" "cell_10383" "cell_11074" "cell_11212"
[776] "cell_10572" "cell_11342" "cell_10593" "cell_10283" "cell_12056"
[781] "cell_9819"  "cell_11771" "cell_11591" "cell_9998"  "cell_11144"
[786] "cell_11015" "cell_10381" "cell_10937" "cell_11265" "cell_10086"
[791] "cell_10187" "cell_10215" "cell_11611" "cell_9911"  "cell_11308"
[796] "cell_10735" "cell_11978" "cell_11067" "cell_10491" "cell_11003"
[801] "cell_10989" "cell_11480" "cell_10255" "cell_10312" "cell_11548"
[806] "cell_11207" "cell_9920"  "cell_10209" "cell_11211" "cell_11422"
[811] "cell_11360" "cell_10305" "cell_10821" "cell_10461" "cell_10008"
[816] "cell_10801" "cell_11319" "cell_10078" "cell_11651" "cell_11287"
[821] "cell_10389" "cell_10135" "cell_10773" "cell_11246" "cell_11623"
[826] "cell_11627" "cell_11031" "cell_10200" "cell_11113" "cell_9794" 
[831] "cell_12023" "cell_11197" "cell_11313" "cell_11129" "cell_11597"
[836] "cell_9807"  "cell_10095" "cell_11860"


$`Somitic mesoderm`
$`Somitic mesoderm`$aucThr
$`Somitic mesoderm`$aucThr$selected
minimumDens 
     0.0567 

$`Somitic mesoderm`$aucThr$thresholds
            threshold nCells
Global_k1      0.1638    171
L_k2           0.1739    159
minimumDens    0.0567    839

$`Somitic mesoderm`$aucThr$comment
[1] "The global distribution overlaps the partial distributions. "


$`Somitic mesoderm`$assignment
  [1] "cell_11995" "cell_9963"  "cell_10910" "cell_11021" "cell_10433"
  [6] "cell_11395" "cell_10779" "cell_10883" "cell_10721" "cell_10116"
 [11] "cell_10785" "cell_11781" "cell_11147" "cell_11842" "cell_11218"
 [16] "cell_11441" "cell_11558" "cell_9979"  "cell_10323" "cell_11663"
 [21] "cell_10935" "cell_11567" "cell_10673" "cell_11634" "cell_10815"
 [26] "cell_10077" "cell_9934"  "cell_11888" "cell_10356" "cell_10490"
 [31] "cell_9921"  "cell_10062" "cell_9809"  "cell_11906" "cell_10084"
 [36] "cell_10951" "cell_10520" "cell_10142" "cell_12176" "cell_11850"
 [41] "cell_10284" "cell_9781"  "cell_11885" "cell_11687" "cell_10523"
 [46] "cell_11100" "cell_10881" "cell_10774" "cell_12107" "cell_11216"
 [51] "cell_10439" "cell_10902" "cell_10856" "cell_11702" "cell_11269"
 [56] "cell_10248" "cell_11883" "cell_11682" "cell_10877" "cell_9933" 
 [61] "cell_10843" "cell_11154" "cell_12052" "cell_9895"  "cell_12028"
 [66] "cell_10454" "cell_11750" "cell_12059" "cell_10615" "cell_10751"
 [71] "cell_10387" "cell_10245" "cell_11767" "cell_9919"  "cell_11976"
 [76] "cell_11971" "cell_10797" "cell_10094" "cell_11724" "cell_11030"
 [81] "cell_11855" "cell_10754" "cell_9905"  "cell_11530" "cell_11352"
 [86] "cell_10357" "cell_12074" "cell_12012" "cell_11561" "cell_11398"
 [91] "cell_12075" "cell_12177" "cell_11653" "cell_10644" "cell_11993"
 [96] "cell_10322" "cell_11900" "cell_11315" "cell_10474" "cell_10480"
[101] "cell_11220" "cell_10440" "cell_11596" "cell_10296" "cell_11856"
[106] "cell_11728" "cell_10241" "cell_12016" "cell_11770" "cell_10917"
[111] "cell_11805" "cell_12002" "cell_10599" "cell_10846" "cell_10687"
[116] "cell_11330" "cell_11081" "cell_9953"  "cell_10181" "cell_10395"
[121] "cell_12021" "cell_11667" "cell_10547" "cell_11329" "cell_10332"
[126] "cell_10562" "cell_11183" "cell_11138" "cell_10385" "cell_10125"
[131] "cell_10047" "cell_10038" "cell_10386" "cell_11673" "cell_10105"
[136] "cell_12162" "cell_10749" "cell_11549" "cell_11458" "cell_11290"
[141] "cell_11213" "cell_10849" "cell_11465" "cell_9885"  "cell_10416"
[146] "cell_10847" "cell_11009" "cell_11733" "cell_10756" "cell_10479"
[151] "cell_11260" "cell_11002" "cell_10117" "cell_10933" "cell_11840"
[156] "cell_10922" "cell_11495" "cell_11380" "cell_10478" "cell_11026"
[161] "cell_11880" "cell_11896" "cell_12110" "cell_11115" "cell_10559"
[166] "cell_12151" "cell_11109" "cell_11014" "cell_11526" "cell_11018"
[171] "cell_11582" "cell_11476" "cell_10065" "cell_10628" "cell_11603"
[176] "cell_10825" "cell_9851"  "cell_10634" "cell_10791" "cell_9844" 
[181] "cell_10886" "cell_11931" "cell_11975" "cell_11150" "cell_10450"
[186] "cell_12115" "cell_11506" "cell_11239" "cell_11130" "cell_10011"
[191] "cell_11922" "cell_11827" "cell_10393" "cell_11982" "cell_10229"
[196] "cell_12140" "cell_9799"  "cell_10907" "cell_11957" "cell_9784" 
[201] "cell_11204" "cell_10292" "cell_11229" "cell_10692" "cell_11370"
[206] "cell_9972"  "cell_10914" "cell_10027" "cell_10262" "cell_11657"
[211] "cell_11552" "cell_11123" "cell_10900" "cell_11093" "cell_9778" 
[216] "cell_9938"  "cell_11001" "cell_11619" "cell_12160" "cell_9876" 
[221] "cell_9776"  "cell_10394" "cell_10029" "cell_11377" "cell_11679"
[226] "cell_11054" "cell_11528" "cell_11913" "cell_11515" "cell_10435"
[231] "cell_10758" "cell_10220" "cell_10624" "cell_11716" "cell_11371"
[236] "cell_12042" "cell_10414" "cell_10828" "cell_10659" "cell_11424"
[241] "cell_11096" "cell_12171" "cell_11803" "cell_11404" "cell_10561"
[246] "cell_11180" "cell_10046" "cell_10009" "cell_9792"  "cell_11853"
[251] "cell_11095" "cell_11511" "cell_11835" "cell_11463" "cell_10146"
[256] "cell_10317" "cell_11407" "cell_12173" "cell_10161" "cell_10438"
[261] "cell_10981" "cell_10079" "cell_11210" "cell_11174" "cell_10876"
[266] "cell_10087" "cell_11630" "cell_11434" "cell_9888"  "cell_10367"
[271] "cell_11107" "cell_11994" "cell_10192" "cell_11822" "cell_10647"
[276] "cell_10436" "cell_11592" "cell_9817"  "cell_11268" "cell_9961" 
[281] "cell_10071" "cell_9959"  "cell_11419" "cell_10585" "cell_11214"
[286] "cell_11084" "cell_11165" "cell_11327" "cell_11392" "cell_10718"
[291] "cell_10466" "cell_12090" "cell_11555" "cell_10995" "cell_10526"
[296] "cell_10761" "cell_11038" "cell_11232" "cell_10715" "cell_10458"
[301] "cell_10019" "cell_10328" "cell_10411" "cell_10313" "cell_9930" 
[306] "cell_11114" "cell_9936"  "cell_11234" "cell_10556" "cell_9846" 
[311] "cell_11871" "cell_11285" "cell_10991" "cell_11939" "cell_10213"
[316] "cell_9863"  "cell_11710" "cell_11256" "cell_10147" "cell_9989" 
[321] "cell_10464" "cell_10953" "cell_10388" "cell_10010" "cell_11485"
[326] "cell_10695" "cell_11760" "cell_11328" "cell_11620" "cell_10694"
[331] "cell_10175" "cell_9997"  "cell_12147" "cell_10553" "cell_11266"
[336] "cell_10467" "cell_11233" "cell_10963" "cell_12149" "cell_9986" 
[341] "cell_11276" "cell_11904" "cell_12179" "cell_9847"  "cell_11303"
[346] "cell_11299" "cell_10956" "cell_10347" "cell_11599" "cell_10697"
[351] "cell_10190" "cell_10261" "cell_10498" "cell_11588" "cell_9977" 
[356] "cell_11839" "cell_11777" "cell_11941" "cell_10645" "cell_10833"
[361] "cell_10672" "cell_10654" "cell_11024" "cell_10164" "cell_11522"
[366] "cell_10407" "cell_12134" "cell_11056" "cell_11489" "cell_9969" 
[371] "cell_9820"  "cell_11017" "cell_9835"  "cell_12064" "cell_11562"
[376] "cell_10345" "cell_10225" "cell_10783" "cell_12036" "cell_11333"
[381] "cell_10812" "cell_10998" "cell_9892"  "cell_10360" "cell_12152"
[386] "cell_10508" "cell_9813"  "cell_11124" "cell_11073" "cell_9859" 
[391] "cell_10421" "cell_10930" "cell_10374" "cell_12097" "cell_11217"
[396] "cell_10800" "cell_11063" "cell_10363" "cell_10515" "cell_10267"
[401] "cell_11047" "cell_11744" "cell_11564" "cell_9849"  "cell_11916"
[406] "cell_12082" "cell_12057" "cell_10664" "cell_11680" "cell_11542"
[411] "cell_10765" "cell_12066" "cell_11438" "cell_12072" "cell_10441"
[416] "cell_10396" "cell_10838" "cell_9832"  "cell_10299" "cell_11402"
[421] "cell_11736" "cell_11012" "cell_11484" "cell_11962" "cell_10604"
[426] "cell_11966" "cell_10370" "cell_10099" "cell_10669" "cell_11774"
[431] "cell_11579" "cell_10059" "cell_10511" "cell_10597" "cell_12031"
[436] "cell_10961" "cell_10497" "cell_11934" "cell_11353" "cell_11133"
[441] "cell_9837"  "cell_11112" "cell_10272" "cell_10868" "cell_11734"
[446] "cell_10554" "cell_11683" "cell_10898" "cell_10623" "cell_10864"
[451] "cell_10217" "cell_10739" "cell_12174" "cell_12053" "cell_11859"
[456] "cell_10103" "cell_12167" "cell_10848" "cell_9975"  "cell_10204"
[461] "cell_11678" "cell_10769" "cell_10369" "cell_10155" "cell_10860"
[466] "cell_10512" "cell_10333" "cell_11149" "cell_10560" "cell_10510"
[471] "cell_11658" "cell_10354" "cell_11141" "cell_11263" "cell_12058"
[476] "cell_9856"  "cell_11208" "cell_9779"  "cell_11598" "cell_10231"
[481] "cell_11547" "cell_10468" "cell_11961" "cell_10610" "cell_10453"
[486] "cell_10048" "cell_11304" "cell_10571" "cell_11778" "cell_10546"
[491] "cell_12111" "cell_11120" "cell_9940"  "cell_12040" "cell_11471"
[496] "cell_11470" "cell_11321" "cell_10236" "cell_11176" "cell_12094"
[501] "cell_10509" "cell_10133" "cell_10359" "cell_11595" "cell_11019"
[506] "cell_10602" "cell_11414" "cell_11838" "cell_11932" "cell_11057"
[511] "cell_10303" "cell_10750" "cell_9880"  "cell_10224" "cell_10885"
[516] "cell_11281" "cell_10399" "cell_9816"  "cell_10737" "cell_10250"
[521] "cell_11040" "cell_9873"  "cell_11782" "cell_10227" "cell_12125"
[526] "cell_11714" "cell_11104" "cell_11354" "cell_10660" "cell_9907" 
[531] "cell_10249" "cell_10611" "cell_10018" "cell_11389" "cell_10098"
[536] "cell_10401" "cell_11418" "cell_12092" "cell_10890" "cell_10850"
[541] "cell_11791" "cell_11988" "cell_11897" "cell_10026" "cell_10403"
[546] "cell_10030" "cell_11688" "cell_11571" "cell_10696" "cell_9974" 
[551] "cell_9967"  "cell_11381" "cell_11972" "cell_12143" "cell_11837"
[556] "cell_10521" "cell_11280" "cell_11694" "cell_10166" "cell_11010"
[561] "cell_10947" "cell_11901" "cell_12044" "cell_10586" "cell_10428"
[566] "cell_10115" "cell_11357" "cell_11121" "cell_12114" "cell_12076"
[571] "cell_9823"  "cell_11162" "cell_10836" "cell_11954" "cell_11713"
[576] "cell_11168" "cell_10700" "cell_11845" "cell_11557" "cell_11346"
[581] "cell_9962"  "cell_11386" "cell_11799" "cell_10869" "cell_11411"
[586] "cell_10766" "cell_11450" "cell_11058" "cell_10481" "cell_10123"
[591] "cell_10238" "cell_10107" "cell_9968"  "cell_11951" "cell_11811"
[596] "cell_11519" "cell_11172" "cell_10661" "cell_10534" "cell_10904"
[601] "cell_12081" "cell_10237" "cell_10247" "cell_10380" "cell_11199"
[606] "cell_11447" "cell_11854" "cell_11524" "cell_10058" "cell_11586"
[611] "cell_10234" "cell_11743" "cell_10256" "cell_11775" "cell_11270"
[616] "cell_10089" "cell_9777"  "cell_10635" "cell_10049" "cell_10410"
[621] "cell_10675" "cell_10476" "cell_10832" "cell_10940" "cell_10788"
[626] "cell_11818" "cell_11502" "cell_10731" "cell_12166" "cell_10646"
[631] "cell_11397" "cell_10151" "cell_12085" "cell_11435" "cell_10239"
[636] "cell_11715" "cell_12060" "cell_11798" "cell_9929"  "cell_10189"
[641] "cell_10212" "cell_9797"  "cell_11529" "cell_11033" "cell_11794"
[646] "cell_9812"  "cell_9992"  "cell_11696" "cell_10680" "cell_10319"
[651] "cell_10120" "cell_10601" "cell_12145" "cell_10451" "cell_11243"
[656] "cell_9824"  "cell_11070" "cell_11983" "cell_11986" "cell_10422"
[661] "cell_11664" "cell_10143" "cell_11368" "cell_11580" "cell_11583"
[666] "cell_11312" "cell_11862" "cell_11493" "cell_11226" "cell_9785" 
[671] "cell_10417" "cell_9931"  "cell_11277" "cell_11261" "cell_11642"
[676] "cell_10977" "cell_10424" "cell_10588" "cell_9803"  "cell_11938"
[681] "cell_11877" "cell_10271" "cell_10177" "cell_10437" "cell_10741"
[686] "cell_10580" "cell_11102" "cell_10945" "cell_11077" "cell_9916" 
[691] "cell_9861"  "cell_10037" "cell_10485" "cell_10443" "cell_10140"
[696] "cell_10285" "cell_11618" "cell_9845"  "cell_11469" "cell_11140"
[701] "cell_10452" "cell_10352" "cell_9894"  "cell_11823" "cell_10944"
[706] "cell_10488" "cell_10867" "cell_11516" "cell_12116" "cell_10755"
[711] "cell_10072" "cell_11393" "cell_11650" "cell_10878" "cell_11400"
[716] "cell_9868"  "cell_10817" "cell_11700" "cell_11085" "cell_12007"
[721] "cell_11998" "cell_12158" "cell_10315" "cell_10698" "cell_10471"
[726] "cell_10111" "cell_11945" "cell_10999" "cell_11225" "cell_11083"
[731] "cell_11539" "cell_12113" "cell_11158" "cell_9802"  "cell_11173"
[736] "cell_10034" "cell_10193" "cell_10831" "cell_12083" "cell_11296"
[741] "cell_11007" "cell_10895" "cell_10725" "cell_11709" "cell_10626"
[746] "cell_10168" "cell_9770"  "cell_10870" "cell_10950" "cell_9924" 
[751] "cell_12006" "cell_9983"  "cell_10787" "cell_10054" "cell_11118"
[756] "cell_10732" "cell_10198" "cell_11923" "cell_10518" "cell_10418"
[761] "cell_10484" "cell_10522" "cell_10153" "cell_11483" "cell_11309"
[766] "cell_11445" "cell_11754" "cell_11914" "cell_10677" "cell_11772"
[771] "cell_11697" "cell_9942"  "cell_10202" "cell_10383" "cell_11074"
[776] "cell_11212" "cell_10572" "cell_11342" "cell_10593" "cell_10283"
[781] "cell_12056" "cell_9819"  "cell_11771" "cell_11591" "cell_9998" 
[786] "cell_11144" "cell_11015" "cell_10381" "cell_10937" "cell_11265"
[791] "cell_10086" "cell_10187" "cell_10215" "cell_11611" "cell_9911" 
[796] "cell_11308" "cell_10735" "cell_11978" "cell_11067" "cell_10491"
[801] "cell_11003" "cell_10989" "cell_11480" "cell_10255" "cell_10312"
[806] "cell_11548" "cell_11207" "cell_9920"  "cell_10209" "cell_11211"
[811] "cell_11422" "cell_11360" "cell_10305" "cell_10821" "cell_10461"
[816] "cell_10008" "cell_10801" "cell_11319" "cell_10078" "cell_11651"
[821] "cell_11287" "cell_10389" "cell_10135" "cell_10773" "cell_11246"
[826] "cell_11623" "cell_11627" "cell_11031" "cell_10200" "cell_11113"
[831] "cell_9794"  "cell_12023" "cell_11197" "cell_11313" "cell_11129"
[836] "cell_11597" "cell_9807"  "cell_10095" "cell_11860"


$`Spinal cord`
$`Spinal cord`$aucThr
$`Spinal cord`$aucThr$selected
outlierOfGlobal 
          0.433 

$`Spinal cord`$aucThr$thresholds
                threshold nCells
outlierOfGlobal    0.4333      0
Global_k1          0.1939    239
L_k2               0.1241    632
R_k3               0.0913    754

$`Spinal cord`$aucThr$comment
[1] "The AUC might follow a normal distribution (random gene-set?). The right distribution is taller. "


$`Spinal cord`$assignment
character(0)
```

:::

::::

## Session Info


``` r
sessionInfo()
```

``` output
R version 4.4.1 (2024-06-14)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 22.04.5 LTS

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
 [83] renv_1.0.9                bit_4.0.5                
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

:::: challenge 

#### Extension Challenge 1: Group pair comparisons

Why do you think marker genes are found by aggregating pairwise comparisons rather than iteratively comparing each cluster to all other clusters? 

::: solution

One important reason why is because averages over all other clusters can be sensitive to the cell type composition. If a rare cell type shows up in one sample, the most discriminative marker genes found in this way could be very different from those found in another sample where the rare cell type is absent. 

Generally, it's good to keep in mind that the concept of "everything else" is not a stable basis for comparison. Read that sentence again, because its a subtle but broadly applicable point. Think about it and you can probably identify analogous issues in fields outside of single-cell analysis. It frequently comes up when comparisons between multiple categories are involved.

:::
::::

:::: challenge

#### Extension Challenge 2: Parallelizing SingleR

SingleR can be computationally expensive. How do you set it to run in parallel?

::: solution

Use `BiocParallel` and the `BPPARAM` argument! This example will set it to use four cores on your laptop, but you can also configure BiocParallel to use cluster jobs.


``` r
library(BiocParallel)

my_bpparam = MulticoreParam(workers = 4)

res2 <- SingleR(test = sce.mat, 
                ref = ref.mat,
                labels = ref$celltype,
                BPPARAM = my_bpparam)
```

`BiocParallel` is the most common way to enable parallel computation in Bioconductor packages, so you can expect to see it elsewhere outside of SingleR.

:::

::::

:::: challenge

#### Extension Challenge 3: Critical inspection of diagnostics

The first set of AUCell diagnostics don't look so good for some of the examples here. Which ones? Why?

::: solution

The example that jumps out most strongly to the eye is ExE endoderm, which doesn't show clear separate modes. Simultaneously, Endothelium seems to have three or four modes. 

Remember, this is an exploratory diagnostic, not the final word! At this point it'd be good to engage in some critical inspection of the results. Maybe we don't have enough / the best marker genes. In this particular case, the fact that we subsetted the reference set to 1000 cells probably didn't help.
:::

::::

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
-   The `scoreMarkers` function from the *[scran](https://bioconductor.org/packages/3.19/scran)* package 
    package can be used to find candidate marker genes for clusters of cells by
    ranking differential expression between pairs of clusters.
-   Computational annotation using published reference datasets or curated gene sets
    provides a fast, automated, and reproducible alternative to the manual
    annotation of cell clusters based on marker gene expression.
-   The *[SingleR](https://bioconductor.org/packages/3.19/SingleR)*
    package is a popular choice for reference-based annotation and assigns labels
    to cells based on the reference samples with the highest Spearman rank correlations.
-   The *[AUCell](https://bioconductor.org/packages/3.19/AUCell)* package provides an enrichment
    test to identify curated marker sets that are highly expressed in each cell. 
:::