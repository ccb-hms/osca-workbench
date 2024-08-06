---
title: Introduction to Bioconductor and the SingleCellExperiment class
teaching: 20 # Minutes of teaching in the lesson
exercises: 10 # Minutes of exercises in the lesson
---

:::::::::::::::::::::::::::::::::::::: questions 

- What is Bioconductor?
- How is single-cell data stored in the Bioconductor ecosystem?
- What is a `SingleCellExperiment` object?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Install and update Bioconductor packages. 
- Load data generated with common single-cell technologies as `SingleCellExperiment` objects. 
- Inspect and manipulate `SingleCellExperiment` objects. 

::::::::::::::::::::::::::::::::::::::::::::::::

## Setup




``` r
library(SingleCellExperiment)
library(MouseGastrulationData)
```

## Bioconductor

### Overview 

Within the R ecosystem, the Bioconductor project provides tools for the analysis and comprehension of high-throughput genomics data.
The scope of the project covers microarray data, various forms of sequencing (RNA-seq, ChIP-seq, bisulfite, genotyping, etc.), proteomics, flow cytometry and more.
One of Bioconductor's main selling points is the use of common data structures to promote interoperability between packages,
allowing code written by different people (from different organizations, in different countries) to work together seamlessly in complex analyses. 
By extending R to genomics, Bioconductor serves as a powerful addition to the computational biologist's toolkit.

### Installing Bioconductor Packages

The default repository for R packages is the [Comprehensive R Archive Network](https://cran.r-project.org/mirrors.html) (CRAN), which is home to over 13,000 different R packages. 
We can easily install packages from CRAN - say, the popular *[ggplot2](https://CRAN.R-project.org/package=ggplot2)* package for data visualization - by opening up R and typing in:


``` r
install.packages("ggplot2")
```

In our case, however, we want to install Bioconductor packages.
These packages are located in a separate repository (see comments below) so we first install the *[BiocManager](https://CRAN.R-project.org/package=BiocManager)* package to easily connect to the Bioconductor servers.


``` r
install.packages("BiocManager")
```

After that, we can use *[BiocManager](https://CRAN.R-project.org/package=BiocManager)*'s `install()` function to install any package from Bioconductor.
For example, the code chunk below uses this approach to install the *[SingleCellExperiment](https://bioconductor.org/packages/3.19/SingleCellExperiment)* package.


``` r
## The command below is a one-line shortcut for:
## library(BiocManager)
## install("SingleCellExperiment")
BiocManager::install("SingleCellExperiment")
```

Should we forget, the same instructions are present on the landing page of any Bioconductor package.
For example, looking at the [`scater`](https://bioconductor.org/packages/release/bioc/html/scater.html) package page on Bioconductor, we can see the following copy-pasteable instructions:


``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("scater")
```

Packages only need to be installed once, and then they are available for all subsequent uses of a particular R installation.
There is no need to repeat the installation every time we start R.

### Finding relevant packages

To find relevant Bioconductor packages, one useful resource is the [BiocViews](https://bioconductor.org/packages/release/BiocViews.html) page.
This provides a hierarchically organized view of annotations associated with each Bioconductor package.
For example, under the ["Software"](https://bioconductor.org/packages/release/BiocViews.html#___Software) label, we might be interested in a particular ["Technology"](https://bioconductor.org/packages/release/BiocViews.html#___Technology) such as... say, ["SingleCell"](https://bioconductor.org/packages/release/BiocViews.html#___SingleCell).
This gives us a listing of all Bioconductor packages that might be useful for our single-cell data analyses. 
CRAN uses the similar concept of ["Task views"](https://cran.r-project.org/web/views/), though this is understandably more general than genomics.
For example, the [Cluster task view page](https://cran.r-project.org/web/views/Cluster.html) lists an assortment of packages that are relevant to cluster analyses.

### Staying up to date

Updating all R/Bioconductor packages is as simple as running `BiocManager::install()` without any arguments.
This will check for more recent versions of each package (within a Bioconductor release) and prompt the user to update if any are available.


``` r
BiocManager::install()
```

## The `SingleCellExperiment` class

One of the main strengths of the Bioconductor project lies in the use of a common data infrastructure that powers interoperability across packages. 

Users should be able to analyze their data using functions from different Bioconductor packages without the need to convert between formats. To this end, the `SingleCellExperiment` class (from the _SingleCellExperiment_ package) serves as the common currency for data exchange across 70+ single-cell-related Bioconductor packages.

This class implements a data structure that stores all aspects of our single-cell data - gene-by-cell expression data, per-cell metadata and per-gene annotation - and manipulate them in a synchronized manner.

<img src="http://bioconductor.org/books/release/OSCA.intro/images/SingleCellExperiment.png" style="display: block; margin: auto;" />

Let's start with an example dataset.


``` r
sce <- WTChimeraData(samples=5)
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

We can think of this (and other) class as a _container_, that contains several different pieces of data in so-called _slots_.

The _getter_ methods are used to extract information from the slots and the _setter_ methods are used to add information into the slots. These are the only ways to interact with the objects (rather than directly accessing the slots).

Depending on the object, slots can contain different types of data (e.g., numeric matrices, lists, etc.). We will here review the main slots of the SingleCellExperiment class as well as their getter/setter methods.

:::: challenge

Before SingleCellExperiments, coders working with single cell data would sometimes keep all of these components in separate objects e.g. a matrix of counts, a data.frame of sample metadata, a data.frame of gene annotations and so on. What are the main disadvantages of this sort of "from scratch" approach?

::: solution

1. You have to do tons of manual book-keeping! If you perform a QC step that removes dead cells, now you also have to remember to remove that same set of cells from the cell-wise metadata. Dropped un-expressed genes? Don't forget to filter the gene metadata table too. 

2. All the downstream steps have to be "from scratch" as well! If your tables have some slight format difference from those of your lab-mate, suddenly the plotting code you're trying to re-use doesn't work! Agh!

:::

::::

### The `assays`

This is arguably the most fundamental part of the object that contains the count matrix, and potentially other matrices with transformed data. We can access the _list_ of matrices with the `assays` function and individual matrices with the `assay` function. If one of these matrices is called "counts", we can use the special `counts` getter (and the analogous `logcounts`).


``` r
identical(assay(sce), counts(sce))
```

``` output
[1] TRUE
```

``` r
counts(sce)[1:3, 1:3]
```

``` output
3 x 3 sparse Matrix of class "dgCMatrix"
                   cell_9769 cell_9770 cell_9771
ENSMUSG00000051951         .         .         .
ENSMUSG00000089699         .         .         .
ENSMUSG00000102343         .         .         .
```

You will notice that in this case we have a sparse matrix of class "dgTMatrix" inside the object. More generally, any "matrix-like" object can be used, e.g., dense matrices or HDF5-backed matrices (see "Working with large data").

### The `colData` and `rowData`

Conceptually, these are two data frames that annotate the columns and the rows of your assay, respectively.

One can interact with them as usual, e.g., by extracting columns or adding additional variables as columns.


``` r
colData(sce)
```

``` output
DataFrame with 2411 rows and 11 columns
                  cell          barcode    sample       stage    tomato
           <character>      <character> <integer> <character> <logical>
cell_9769    cell_9769 AAACCTGAGACTGTAA         5        E8.5      TRUE
cell_9770    cell_9770 AAACCTGAGATGCCTT         5        E8.5      TRUE
cell_9771    cell_9771 AAACCTGAGCAGCCTC         5        E8.5      TRUE
cell_9772    cell_9772 AAACCTGCATACTCTT         5        E8.5      TRUE
cell_9773    cell_9773 AAACGGGTCAACACCA         5        E8.5      TRUE
...                ...              ...       ...         ...       ...
cell_12175  cell_12175 TTTGGTTAGTCCGTAT         5        E8.5      TRUE
cell_12176  cell_12176 TTTGGTTAGTGTTGAA         5        E8.5      TRUE
cell_12177  cell_12177 TTTGGTTGTTAAAGAC         5        E8.5      TRUE
cell_12178  cell_12178 TTTGGTTTCAGTCAGT         5        E8.5      TRUE
cell_12179  cell_12179 TTTGGTTTCGCCATAA         5        E8.5      TRUE
                pool stage.mapped        celltype.mapped closest.cell
           <integer>  <character>            <character>  <character>
cell_9769          3        E8.25             Mesenchyme   cell_24159
cell_9770          3         E8.5            Endothelium   cell_96660
cell_9771          3         E8.5              Allantois  cell_134982
cell_9772          3         E8.5             Erythroid3  cell_133892
cell_9773          3        E8.25             Erythroid1   cell_76296
...              ...          ...                    ...          ...
cell_12175         3         E8.5             Erythroid3  cell_138060
cell_12176         3         E8.5 Forebrain/Midbrain/H..   cell_72709
cell_12177         3        E8.25       Surface ectoderm  cell_100275
cell_12178         3        E8.25             Erythroid2   cell_70906
cell_12179         3         E8.5            Spinal cord  cell_102334
           doub.density sizeFactor
              <numeric>  <numeric>
cell_9769    0.02985045    1.41243
cell_9770    0.00172753    1.22757
cell_9771    0.01338013    1.15439
cell_9772    0.00218402    1.28676
cell_9773    0.00211723    1.78719
...                 ...        ...
cell_12175   0.00129403   1.219506
cell_12176   0.01833074   1.095753
cell_12177   0.03104037   0.910728
cell_12178   0.00169483   2.061701
cell_12179   0.03767894   1.798687
```

``` r
rowData(sce)
```

``` output
DataFrame with 29453 rows and 2 columns
                              ENSEMBL         SYMBOL
                          <character>    <character>
ENSMUSG00000051951 ENSMUSG00000051951           Xkr4
ENSMUSG00000089699 ENSMUSG00000089699         Gm1992
ENSMUSG00000102343 ENSMUSG00000102343        Gm37381
ENSMUSG00000025900 ENSMUSG00000025900            Rp1
ENSMUSG00000025902 ENSMUSG00000025902          Sox17
...                               ...            ...
ENSMUSG00000095041 ENSMUSG00000095041     AC149090.1
ENSMUSG00000063897 ENSMUSG00000063897          DHRSX
ENSMUSG00000096730 ENSMUSG00000096730       Vmn2r122
ENSMUSG00000095742 ENSMUSG00000095742 CAAA01147332.1
tomato-td                   tomato-td      tomato-td
```

Note the `$` short cut.


``` r
identical(colData(sce)$sum, sce$sum)
```

``` output
[1] TRUE
```

``` r
sce$my_sum <- colSums(counts(sce))
colData(sce)
```

``` output
DataFrame with 2411 rows and 12 columns
                  cell          barcode    sample       stage    tomato
           <character>      <character> <integer> <character> <logical>
cell_9769    cell_9769 AAACCTGAGACTGTAA         5        E8.5      TRUE
cell_9770    cell_9770 AAACCTGAGATGCCTT         5        E8.5      TRUE
cell_9771    cell_9771 AAACCTGAGCAGCCTC         5        E8.5      TRUE
cell_9772    cell_9772 AAACCTGCATACTCTT         5        E8.5      TRUE
cell_9773    cell_9773 AAACGGGTCAACACCA         5        E8.5      TRUE
...                ...              ...       ...         ...       ...
cell_12175  cell_12175 TTTGGTTAGTCCGTAT         5        E8.5      TRUE
cell_12176  cell_12176 TTTGGTTAGTGTTGAA         5        E8.5      TRUE
cell_12177  cell_12177 TTTGGTTGTTAAAGAC         5        E8.5      TRUE
cell_12178  cell_12178 TTTGGTTTCAGTCAGT         5        E8.5      TRUE
cell_12179  cell_12179 TTTGGTTTCGCCATAA         5        E8.5      TRUE
                pool stage.mapped        celltype.mapped closest.cell
           <integer>  <character>            <character>  <character>
cell_9769          3        E8.25             Mesenchyme   cell_24159
cell_9770          3         E8.5            Endothelium   cell_96660
cell_9771          3         E8.5              Allantois  cell_134982
cell_9772          3         E8.5             Erythroid3  cell_133892
cell_9773          3        E8.25             Erythroid1   cell_76296
...              ...          ...                    ...          ...
cell_12175         3         E8.5             Erythroid3  cell_138060
cell_12176         3         E8.5 Forebrain/Midbrain/H..   cell_72709
cell_12177         3        E8.25       Surface ectoderm  cell_100275
cell_12178         3        E8.25             Erythroid2   cell_70906
cell_12179         3         E8.5            Spinal cord  cell_102334
           doub.density sizeFactor    my_sum
              <numeric>  <numeric> <numeric>
cell_9769    0.02985045    1.41243     27577
cell_9770    0.00172753    1.22757     29309
cell_9771    0.01338013    1.15439     28795
cell_9772    0.00218402    1.28676     34794
cell_9773    0.00211723    1.78719     38300
...                 ...        ...       ...
cell_12175   0.00129403   1.219506     26680
cell_12176   0.01833074   1.095753     19013
cell_12177   0.03104037   0.910728     24627
cell_12178   0.00169483   2.061701     46162
cell_12179   0.03767894   1.798687     38398
```

### The `reducedDims`

Everything that we have described so far (except for the `counts` getter) is part of the `SummarizedExperiment` class that SingleCellExperiment extends. You can find a complete lesson on the `SummarizedExperiment` class in [Introduction to data analysis with R and Bioconductor](https://carpentries-incubator.github.io/bioc-intro/60-next-steps.html) course.

One of the peculiarity of SingleCellExperiment is its ability to store reduced dimension matrices within the object. These may include PCA, t-SNE, UMAP, etc.


``` r
reducedDims(sce)
```

``` output
List of length 2
names(2): pca.corrected.E7.5 pca.corrected.E8.5
```

As for the other slots, we have the usual setter/getter, but it is somewhat rare to interact directly with these functions.

It is more common for other functions to _store_ this information in the object, e.g., the `runPCA` function from the `scater` package.

Here, we use `scater`'s `plotReducedDim` function as an example of how to extract this information _indirectly_ from the objects. Note that one could obtain the same results (somewhat less efficiently) by extracting the corresponding `reducedDim` matrix and `ggplot`.


``` r
library(scater)
plotReducedDim(sce, "pca.corrected.E8.5", colour_by = "celltype.mapped")
```

<img src="fig/intro-sce-rendered-unnamed-chunk-12-1.png" style="display: block; margin: auto;" />

:::::::::::::::::::::::::::::::::: challenge

#### Exercise 1

Create a `SingleCellExperiment` object: Try and create a `SingleCellExperiment` object "from scratch". Start from a `matrix` (either randomly generated or with some fake data in it) and add one or more columns as `colData`.

:::::::::::::: hint

The `SingleCellExperiment` constructor function can be used to create a new `SingleCellExperiment` object.

:::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::: challenge

#### Exercise 2

Combining two objects: The `MouseGastrulationData` package contains several datasets. Download sample 6 of the chimera experiment by running `sce6 <- WTChimeraData(sample=6)`. Use the `cbind` function to combine the new data with the `sce` object created before. 

:::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::: checklist

## Further Reading

* OSCA book, [Introduction](https://bioconductor.org/books/release/OSCA.intro)

::::::::::::::

::::::::::::::::::::::::::::::::::::: keypoints 

- The Bioconductor project provides open-source software packages for the comprehension of high-throughput biological data.
- A `SingleCellExperiment` object is an extension of the `SummarizedExperiment` object.
- `SingleCellExperiment` objects contain specialized data fields for storing data unique to single-cell analyses, such as the `reducedDims` field. 

::::::::::::::::::::::::::::::::::::::::::::::::


