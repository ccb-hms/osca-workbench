# Orchestrating Single-Cell Analysis with Bioconductor

## Description

In the last few years, the profiling of a large number of genome-wide features
in individual cells has become routine. Consequently, a plethora of tools for
the analysis of single-cell data has been developed, making it hard to understand
the critical steps in the analysis workflow and the best methods for each objective
of one’s study.

This [Carpentries-style](https://carpentries.org/) tutorial aims to provide a
solid foundation in using [Bioconductor](https://bioconductor.org) tools for
single-cell RNA-seq analysis by walking through various steps of typical workflows
using example datasets.

This tutorial uses as a "text-book" the online book "Orchestrating Single-Cell
Analysis with Bioconductor" ([OSCA](https://bioconductor.org/books/release/OSCA/)),
[published in 2020](https://doi.org/10.1038%2Fs41592-019-0654-x), 
and continuously updated by many contributors from the Bioconductor
community. Like the book, this tutorial strives to be of interest to the
experimental biologists wanting to analyze their data and to the bioinformaticians
approaching single-cell data.

## Learning objectives

Attendees will learn how to analyze multi-condition single-cell RNA-seq from
raw data to statistical analyses and result interpretation. Students will learn
where the critical steps and methods choices are and will be able to leverage
large-data resources to analyze datasets comprising millions of cells.

In particular, participants will learn:

* How to access publicly available data, such as those from the Human Cell Atlas.
* How to perform data exploration, normalization, and dimensionality reduction.
* How to identify cell types/states and marker genes.
* How to correct for batch effects and integrate multiple samples.
* How to perform differential expression and differential abundance analysis between conditions.
* How to work with large out-of-memory datasets.
* How to interoperate with other popular single-cell analysis ecosystems.

## Other tools and tutorials for single-cell analysis

The focus of this tutorial is on single-cell analysis with R packages from the 
[Bioconductor](https://bioconductor.org) repository. Bioconductor packages are
collaboratively developed by an international community of developers that agree
on data and software standards to promote interoperability between packages,
extensibility of analysis workflows, and reproducibility of published research.

Other popular tools for single-cell analysis include:

* [Seurat](https://satijalab.org/seurat/), a stand-alone R package that has
pioneered elementary steps of typical single-cell analysis workflows, and 
* [scverse](https://scverse.org/), a collection of Python packages for single-cell
omics data analysis including [scanpy](https://scanpy.readthedocs.io) and
[scvi-tools](https://scvi-tools.org/).

Tutorials for working with these tools are available elsewhere and are not covered
in this tutorial. A demonstration of how to interoperate with `Seurat` and packages
from the `scverse` is given in [Session 5](https://ccb-hms.github.io/osca-workbench/large_data.html)
of this tutorial.

Other Carpentries-style tutorials for single-cell analysis with a different scope include:

- a [community-developed lesson](https://carpentries-incubator.github.io/scrna-seq-analysis/)
 that makes use of command-line utilities and `scanpy` for basic preprocessing steps, 
- and a [tutorial proposal](https://github.com/carpentries-incubator/proposals/issues/178)
based on `Seurat`. 

## Source

This lesson uses [The Carpentries Workbench](https://carpentries.github.io/sandpaper-docs/)
and is based on materials from the [OSCA tutorial at the ISMB 2023](https://bioconductor.github.io/ISMB.OSCA/).

## Citation

Amezquita RA, Lun ATL, Becht E, Carey VJ, Carpp LN, Geistlinger L, Marini F,
Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pagès H, Smith ML, Huber W,
Morgan M, Gottardo R, Hicks SC. Orchestrating single-cell analysis with
Bioconductor. *Nature Methods*, 2020.
doi: [10.1038/s41592-019-0654-x](https://doi.org/10.1038/s41592-019-0654-x)
