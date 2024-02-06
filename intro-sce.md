---
title: Introduction to Bioconductor and the SingleCellExperiment class
teaching: 10 # Minutes of teaching in the lesson
exercises: 2 # Minutes of exercises in the lesson
---

:::::::::::::::::::::::::::::::::::::: questions 

- What is Bioconductor?
- How is single-cell data store in the Bioconductor ecosystem?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- TODO

::::::::::::::::::::::::::::::::::::::::::::::::

# Setup


```r
library("SummarizedExperiment")
library("SingleCellExperiment")
library("MouseGastrulationData")
```

# The `SingleCellExperiment` class

One of the main strengths of the Bioconductor project lies in the use of a common data infrastructure that powers interoperability across packages. 

Users should be able to analyze their data using functions from different Bioconductor packages without the need to convert between formats. To this end, the `SingleCellExperiment` class (from the _SingleCellExperiment_ package) serves as the common currency for data exchange across 70+ single-cell-related Bioconductor packages.

This class implements a data structure that stores all aspects of our single-cell data - gene-by-cell expression data, per-cell metadata and per-gene annotation - and manipulate them in a synchronized manner.

<img src="http://bioconductor.org/books/3.17/OSCA.intro/images/SingleCellExperiment.png" style="display: block; margin: auto;" />

Let's start with an example dataset.


```r
sce <- WTChimeraData(samples=5)
```

```{.error}
Error in `collect()`:
! Failed to collect lazy table.
Caused by error in `db_collect()` at dbplyr/R/verb-compute.R:131:2:
! Arguments in `...` must be used.
✖ Problematic argument:
• ..1 = Inf
ℹ Did you misspell an argument name?
```

```r
sce
```

```{.error}
Error in eval(expr, envir, enclos): object 'sce' not found
```

We can think of this (and other) class as a _container_, that contains several different pieces of data in so-called _slots_.

The _getter_ methods are used to extract information from the slots and the _setter_ methods are used to add information into the slots. These are the only ways to interact with the objects (rather than directly accessing the slots).

Depending on the object, slots can contain different types of data (e.g., numeric matrices, lists, etc.). We will here review the main slots of the SingleCellExperiment class as well as their getter/setter methods.

## The `assays`

This is arguably the most fundamental part of the object that contains the count matrix, and potentially other matrices with transformed data. We can access the _list_ of matrices with the `assays` function and individual matrices with the `assay` function. If one of these matrices is called "counts", we can use the special `counts` getter (and the analogous `logcounts`).


```r
identical(assay(sce), counts(sce))
```

```{.error}
Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'assay': object 'sce' not found
```

```r
counts(sce)[1:3, 1:3]
```

```{.error}
Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'counts': object 'sce' not found
```

You will notice that in this case we have a sparse matrix of class "dgTMatrix" inside the object. More generally, any "matrix-like" object can be used, e.g., dense matrices or HDF5-backed matrices (see "Working with large data").

## The `colData` and `rowData`

Conceptually, these are two data frames that annotate the columns and the rows of your assay, respectively.

One can interact with them as usual, e.g., by extracting columns or adding additional variables as columns.


```r
colData(sce)
```

```{.error}
Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'colData': object 'sce' not found
```

```r
rowData(sce)
```

```{.error}
Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'rowData': object 'sce' not found
```

Note the `$` short cut.


```r
identical(colData(sce)$sum, sce$sum)
```

```{.error}
Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'colData': object 'sce' not found
```

```r
sce$my_sum <- colSums(counts(sce))
```

```{.error}
Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'colSums': error in evaluating the argument 'object' in selecting a method for function 'counts': object 'sce' not found
```

```r
colData(sce)
```

```{.error}
Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'colData': object 'sce' not found
```

## The `reducedDims`

Everything that we have described so far (except for the `counts` getter) is part of the `SummarizedExperiment` class that SingleCellExperiment extends.

One of the peculiarity of SingleCellExperiment is its ability to store reduced dimension matrices within the object. These may include PCA, t-SNE, UMAP, etc.


```r
reducedDims(sce)
```

```{.error}
Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'reducedDims': object 'sce' not found
```

As for the other slots, we have the usual setter/getter, but it is somewhat rare to interact directly with these functions.

It is more common for other functions to _store_ this information in the object, e.g., the `runPCA` function from the `scater` package.

Here, we use `scater`'s `plotReducedDim` function as an example of how to extract this information _indirectly_ from the objects. Note that one could obtain the same results (somewhat less efficiently) by extracting the corresponding `reducedDim` matrix and `ggplot`.


```r
library(scater)
```

```{.output}
Loading required package: scuttle
```

```{.output}
Loading required package: ggplot2
```

```r
plotReducedDim(sce, "pca.corrected.E8.5", colour_by = "celltype.mapped")
```

```{.error}
Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'reducedDim': object 'sce' not found
```

:::::::::::::::::::::::::::::::::: challenge

#### Exercise 1

Create a `SingleCellExperiment` object: Try and create a SingleCellExperiment object "from scratch". Start from a matrix (either randomly generated or with some fake data in it) and add one or more columns as colData.

:::::::::::::: hint

The `SingleCellExperiment` function can be used to create a new SingleCellExperiment object

:::::::::::::::::::::::

:::::::::::::: solution

TODO
:::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::: challenge

#### Exercise 2

Combining two objects: The `MouseGastrulationData` package contains several datasets. Download sample 6 of the chimera experiment by running `sce6 <- WTChimeraData(sample=6)`. Use the `cbind` function to combine the new data with the `sce` object created before. 

:::::::::::::: solution

TODO
:::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::::::


::::::::::::::::::::::::::::::::::::: keypoints 

- TODO

::::::::::::::::::::::::::::::::::::::::::::::::

# Further Reading

* OSCA book, [Introduction](https://bioconductor.org/books/release/OSCA.intro)
