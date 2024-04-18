# Orchestrating Single-Cell Analysis with Bioconductor

This is a pre-alpha adaptation of the 2023 ISMB OSCA tutorial using The Carpentries Workbench lesson template. 

### Vignettes converted

As individual vignettes are converted into lessons, they can be added to `config.yaml` to be rendered and shown in the final Github Pages lesson.  

## Description

In the last few years, the profiling of a large number of genome-wide features
in individual cells has become routine. Consequently, a plethora of tools for
the analysis of single-cell data has been developed, making it hard to understand
the critical steps in the analysis workflow and the best methods for each objective
of one’s study.

This tutorial aims to provide a solid foundation in using Bioconductor tools
for single-cell RNA-seq analysis by walking through various steps of typical
workflows using example datasets.

This tutorial uses as a "text-book" the online book "Orchestrating Single-Cell
Analysis with Bioconductor"
([OSCA](https://bioconductor.org/books/release/OSCA/)), 
started in 2018 and continuously updated by many contributors from the Bioconductor
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

## Source

This lesson is a template lesson that uses [The Carpentries Workbench](https://carpentries.github.io/sandpaper-docs/) and is based on materials from the [OSCA tutorial at the ISMB 2023](https://bioconductor.github.io/ISMB.OSCA/).

