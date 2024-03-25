---
title: Accessing data from the Human Cell Atlas (HCA)
teaching: 20 # Minutes of teaching in the lesson
exercises: 10 # Minutes of exercises in the lesson
---

:::::::::::::::::::::::::::::::::::::: questions 

- TODO

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- TODO

::::::::::::::::::::::::::::::::::::::::::::::::

# HCA Project

The Human Cell Atlas (HCA) is a large project that aims to learn from and map
every cell type in the human body. The project extracts spatial and molecular
characteristics in order to understand cellular function and networks. It is an
international collaborative that charts healthy cells in the human body at all
ages. There are about 37.2 trillion cells in the human body. To read more about
the project, head over to their website at https://www.humancellatlas.org.

# CELLxGENE

CELLxGENE is a database and a suite of tools that help scientists to find,
download, explore, analyze, annotate, and publish single cell data. It includes
several analytic and visualization tools to help you to discover single cell
data patterns. To see the list of tools, browse to
https://cellxgene.cziscience.com/.

# CELLxGENE | Census

The Census provides efficient computational tooling to access, query, and
analyze all single-cell RNA data from CZ CELLxGENE Discover. Using a new access
paradigm of cell-based slicing and querying, you can interact with the data
through TileDB-SOMA, or get slices in AnnData or Seurat objects, thus
accelerating your research by significantly minimizing data harmonization at
https://chanzuckerberg.github.io/cellxgene-census/.

# The CuratedAtlasQueryR Project

To systematically characterize the immune system across tissues, demographics
and multiple studies, single cell transcriptomics data was harmonized from the
CELLxGENE database. Data from 28,975,366 cells that cover 156 tissues (excluding
cell cultures), 12,981 samples, and 324 studies were collected. The metadata was
standardized, including sample identifiers, tissue labels (based on anatomy) and
age. Also, the gene-transcript abundance of all samples was harmonized by
putting values on the positive natural scale (i.e. non-logarithmic).

To model the immune system across studies, we adopted a consistent immune
cell-type ontology appropriate for lymphoid and non-lymphoid tissues. We applied
a consensus cell labeling strategy between the Seurat blueprint and Monaco
[-@Monaco2019] to minimize biases in immune cell classification from
study-specific standards.

`CuratedAtlasQueryR` supports data access and programmatic exploration of the
harmonized atlas. Cells of interest can be selected based on ontology, tissue of
origin, demographics, and disease. For example, the user can select CD4 T helper
cells across healthy and diseased lymphoid tissue. The data for the selected
cells can be downloaded locally into popular single-cell data containers. Pseudo
bulk counts are also available to facilitate large-scale, summary analyses of
transcriptional profiles. This platform offers a standardized workflow for
accessing atlas-level datasets programmatically and reproducibly.


```{.error}
Error in knitr::include_graphics("figures/HCA_sccomp_SUPPLEMENTARY_technical_cartoon_curatedAtlasQuery.png"): Cannot find the file(s): "figures/HCA_sccomp_SUPPLEMENTARY_technical_cartoon_curatedAtlasQuery.png"
```

# Data Sources in R / Bioconductor

There are a few options to access single cell data with R / Bioconductor.

| Package | Target | Description |
|---------|-------------|---------|
| [hca](https://bioconductor.org/packages/hca) | [HCA Data Portal API](https://www.humancellatlas.org/data-portal/) | Project, Sample, and File level HCA data |
| [cellxgenedp](https://bioconductor.org/packages/cellxgenedp) | [CellxGene](https://cellxgene.cziscience.com/) | Human and mouse SC data including HCA |
| [CuratedAtlasQueryR](https://stemangiola.github.io/CuratedAtlasQueryR/) | [CellxGene](https://cellxgene.cziscience.com/) | fine-grained query capable CELLxGENE data including HCA |

# Installation


```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("stemangiola/CuratedAtlasQueryR")
```

# Package load 


```r
library(CuratedAtlasQueryR)
library(dplyr)
```

# HCA Metadata

The metadata allows the user to get a lay of the land of what is available
via the package. In this example, we are using the sample database URL which
allows us to get a small and quick subset of the available metadata.


```r
metadata <- get_metadata(remote_url = CuratedAtlasQueryR::SAMPLE_DATABASE_URL)
```

```{.output}
ℹ Downloading 1 file, totalling 0 GB
```

```{.output}
ℹ Downloading https://object-store.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/metadata/sample_metadata.0.2.3.parquet to /home/runner/.cache/R/CuratedAtlasQueryR/metadata.0.2.3.parquet
```

Get a view of the first 10 columns in the metadata with `glimpse`


```r
metadata |>
  select(1:10) |>
  glimpse()
```

```{.output}
Rows: ??
Columns: 10
Database: DuckDB v0.10.0 [unknown@Linux 6.5.0-1016-azure:R 4.3.3/:memory:]
$ cell_                             <chr> "TTATGCTAGGGTGTTG_12", "GCTTGAACATGG…
$ sample_                           <chr> "039c558ca1c43dc74c563b58fe0d6289", …
$ cell_type                         <chr> "mature NK T cell", "mature NK T cel…
$ cell_type_harmonised              <chr> "immune_unclassified", "cd8 tem", "i…
$ confidence_class                  <dbl> 5, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1, …
$ cell_annotation_azimuth_l2        <chr> "gdt", "cd8 tem", "cd8 tem", "cd8 te…
$ cell_annotation_blueprint_singler <chr> "cd4 tem", "cd8 tem", "cd8 tcm", "cl…
$ cell_annotation_monaco_singler    <chr> "natural killer", "effector memory c…
$ sample_id_db                      <chr> "0c1d320a7d0cbbc281a535912722d272", …
$ `_sample_name`                    <chr> "BPH340PrSF_Via___transition zone of…
```

# A note on the piping operator

The vignette materials provided by `CuratedAtlasQueryR` show the use of the
'native' R pipe (implemented after R version `4.1.0`). For those not familiar
with the pipe operator (`|>`), it allows you to chain functions by passing the
left-hand side (LHS) to the first input (typically) on the right-hand side
(RHS). 

In this example, we are extracting the `iris` data set from the `datasets`
package and 'then' taking a subset where the sepal lengths are greater than 5
and 'then' summarizing the data for each level in the `Species` variable with a
`mean`. The pipe operator can be read as 'then'.


```r
data("iris", package = "datasets")

iris |>
  subset(Sepal.Length > 5) |>
  aggregate(. ~ Species, data = _, mean)
```

```{.output}
     Species Sepal.Length Sepal.Width Petal.Length Petal.Width
1     setosa     5.313636    3.713636     1.509091   0.2772727
2 versicolor     5.997872    2.804255     4.317021   1.3468085
3  virginica     6.622449    2.983673     5.573469   2.0326531
```

# Summarizing the metadata

For each distinct tissue and dataset combination, count the number of datasets
by tissue type. 


```r
metadata |>
  distinct(tissue, dataset_id) |> 
  count(tissue)
```

```{.output}
# Source:   SQL [?? x 2]
# Database: DuckDB v0.10.0 [unknown@Linux 6.5.0-1016-azure:R 4.3.3/:memory:]
   tissue                          n
   <chr>                       <dbl>
 1 renal medulla                   6
 2 caecum                          1
 3 ileum                           1
 4 lymph node                      2
 5 transition zone of prostate     2
 6 peripheral zone of prostate     2
 7 fovea centralis                 1
 8 adrenal gland                   1
 9 heart left ventricle            7
10 bone marrow                     4
# ℹ more rows
```

# Columns available in the metadata


```r
head(names(metadata), 10)
```

```{.output}
! The `names()` method of <tbl_lazy> is for internal use only.
ℹ Did you mean `colnames()`?
```

```{.output}
[1] "src"        "lazy_query"
```

# Available assays


```r
metadata |>
    distinct(assay, dataset_id) |>
    count(assay)
```

```{.output}
# Source:   SQL [?? x 2]
# Database: DuckDB v0.10.0 [unknown@Linux 6.5.0-1016-azure:R 4.3.3/:memory:]
   assay                              n
   <chr>                          <dbl>
 1 10x 3' v3                         21
 2 Slide-seq                          4
 3 10x 3' v2                         27
 4 Visium Spatial Gene Expression     7
 5 10x 5' v1                          7
 6 scRNA-seq                          4
 7 Seq-Well                           2
 8 10x 5' v2                          2
 9 10x 3' v1                          1
10 Smart-seq2                         1
# ℹ more rows
```

# Available organisms


```r
metadata |>
    distinct(organism, dataset_id) |>
    count(organism)
```

```{.output}
# Source:   SQL [1 x 2]
# Database: DuckDB v0.10.0 [unknown@Linux 6.5.0-1016-azure:R 4.3.3/:memory:]
  organism         n
  <chr>        <dbl>
1 Homo sapiens    63
```

## Download single-cell RNA sequencing counts 

The data can be provided as either "counts" or counts per million "cpm" as given
by the `assays` argument in the `get_single_cell_experiment()` function. By
default, the `SingleCellExperiment` provided will contain only the 'counts'
data.

### Query raw counts


```r
single_cell_counts <- 
    metadata |>
    dplyr::filter(
        ethnicity == "African" &
        stringr::str_like(assay, "%10x%") &
        tissue == "lung parenchyma" &
        stringr::str_like(cell_type, "%CD4%")
    ) |>
    get_single_cell_experiment()
```

```{.output}
ℹ Realising metadata.
```

```{.output}
ℹ Synchronising files
```

```{.output}
ℹ Downloading 2 files, totalling 0.17 GB
```

```{.output}
ℹ Downloading https://object-store.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/cellxgene-0.2.1-hdf5/original/bc380dae8b14313a870973697842878b/assays.h5 to /home/runner/.cache/R/CuratedAtlasQueryR/0.2.1/original/bc380dae8b14313a870973697842878b/assays.h5
```

```{.output}
Downloading files ■■■■■■■■■■■■■■■■                  50% |  ETA: 13s
```

```{.output}
ℹ Downloading https://object-store.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/cellxgene-0.2.1-hdf5/original/bc380dae8b14313a870973697842878b/se.rds to /home/runner/.cache/R/CuratedAtlasQueryR/0.2.1/original/bc380dae8b14313a870973697842878b/se.rds
```

```{.output}
Downloading files ■■■■■■■■■■■■■■■■                  50% |  ETA: 13sℹ Reading files.
ℹ Compiling Single Cell Experiment.
```

```r
single_cell_counts
```

```{.output}
class: SingleCellExperiment 
dim: 36229 1571 
metadata(0):
assays(1): counts
rownames(36229): A1BG A1BG-AS1 ... ZZEF1 ZZZ3
rowData names(0):
colnames(1571): ACACCAAAGCCACCTG_SC18_1 TCAGCTCCAGACAAGC_SC18_1 ...
  CAGCATAAGCTAACAA_F02607_1 AAGGAGCGTATAATGG_F02607_1
colData names(56): sample_ cell_type ... updated_at_y original_cell_id
reducedDimNames(0):
mainExpName: NULL
altExpNames(0):
```

### Query counts scaled per million

This is helpful if just few genes are of interest, as they can be compared
across samples.


```r
metadata |>
  dplyr::filter(
      ethnicity == "African" &
      stringr::str_like(assay, "%10x%") &
      tissue == "lung parenchyma" &
      stringr::str_like(cell_type, "%CD4%")
  ) |>
  get_single_cell_experiment(assays = "cpm")
```

```{.output}
ℹ Realising metadata.
```

```{.output}
ℹ Synchronising files
```

```{.output}
ℹ Downloading 2 files, totalling 0.29 GB
```

```{.output}
ℹ Downloading https://object-store.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/cellxgene-0.2.1-hdf5/cpm/bc380dae8b14313a870973697842878b/assays.h5 to /home/runner/.cache/R/CuratedAtlasQueryR/0.2.1/cpm/bc380dae8b14313a870973697842878b/assays.h5
```

```{.output}
Downloading files ■■■■■■■■■■■■■■■■                  50% |  ETA: 21s
```

```{.output}
ℹ Downloading https://object-store.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/cellxgene-0.2.1-hdf5/cpm/bc380dae8b14313a870973697842878b/se.rds to /home/runner/.cache/R/CuratedAtlasQueryR/0.2.1/cpm/bc380dae8b14313a870973697842878b/se.rds
```

```{.output}
Downloading files ■■■■■■■■■■■■■■■■                  50% |  ETA: 21sℹ Reading files.
ℹ Compiling Single Cell Experiment.
```

```{.output}
class: SingleCellExperiment 
dim: 36229 1571 
metadata(0):
assays(1): cpm
rownames(36229): A1BG A1BG-AS1 ... ZZEF1 ZZZ3
rowData names(0):
colnames(1571): ACACCAAAGCCACCTG_SC18_1 TCAGCTCCAGACAAGC_SC18_1 ...
  CAGCATAAGCTAACAA_F02607_1 AAGGAGCGTATAATGG_F02607_1
colData names(56): sample_ cell_type ... updated_at_y original_cell_id
reducedDimNames(0):
mainExpName: NULL
altExpNames(0):
```

### Extract only a subset of genes


```r
single_cell_counts <-
    metadata |>
    dplyr::filter(
        ethnicity == "African" &
        stringr::str_like(assay, "%10x%") &
        tissue == "lung parenchyma" &
        stringr::str_like(cell_type, "%CD4%")
    ) |>
    get_single_cell_experiment(assays = "cpm", features = "PUM1")
```

```{.output}
ℹ Realising metadata.
```

```{.output}
ℹ Synchronising files
```

```{.output}
ℹ Downloading 0 files, totalling 0 GB
```

```{.output}
ℹ Reading files.
```

```{.output}
ℹ Compiling Single Cell Experiment.
```

```r
single_cell_counts
```

```{.output}
class: SingleCellExperiment 
dim: 1 1571 
metadata(0):
assays(1): cpm
rownames(1): PUM1
rowData names(0):
colnames(1571): ACACCAAAGCCACCTG_SC18_1 TCAGCTCCAGACAAGC_SC18_1 ...
  CAGCATAAGCTAACAA_F02607_1 AAGGAGCGTATAATGG_F02607_1
colData names(56): sample_ cell_type ... updated_at_y original_cell_id
reducedDimNames(0):
mainExpName: NULL
altExpNames(0):
```

### Extracting counts as a Seurat object

If needed, the H5 `SingleCellExperiment` can be converted into a Seurat object.
Note that it may take a long time and use a lot of memory depending on how many
cells you are requesting.


```r
single_cell_counts <-
    metadata |>
    dplyr::filter(
        ethnicity == "African" &
        stringr::str_like(assay, "%10x%") &
        tissue == "lung parenchyma" &
        stringr::str_like(cell_type, "%CD4%")
    ) |>
    get_seurat()

single_cell_counts
```

## Save your `SingleCellExperiment`

### Saving as HDF5 

The recommended way of saving these `SingleCellExperiment` objects, if
necessary, is to use `saveHDF5SummarizedExperiment` from the `HDF5Array`
package.


```r
single_cell_counts |> saveHDF5SummarizedExperiment("single_cell_counts")
```

# Exercises

:::::::::::::::::::::::::::::::::: challenge

#### Exercise 1

Use `count` and `arrange` to get the number of cells per tissue in descending
order.

:::::::::::::: solution


```r
metadata |>
    count(tissue) |>
    arrange(-n)
```
:::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::: challenge

#### Exercise 2

Use `dplyr`-isms to group by `tissue` and `cell_type` and get a tally of the
highest number of cell types per tissue combination. What tissue has the most
numerous type of cells? 

:::::::::::::: solution


```r
metadata |>
    group_by(tissue, cell_type) |>
    count() |>
    arrange(-n)
```
:::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::: challenge

#### Exercise 3

Spot some differences between the `tissue` and `tissue_harmonised` columns.
Use `count` to summarise.

:::::::::::::: solution


```r
metadata |>
    count(tissue) |>
    arrange(-n)

metadata |>
    count(tissue_harmonised) |>
    arrange(-n)
```
To see the full list of curated columns in the metadata, see the Details section
in the `?get_metadata` documentation page.
    
:::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::: challenge

#### Exercise 4

Now that we are a little familiar with navigating the metadata, let's obtain
a `SingleCellExperiment` of 10X scRNA-seq counts of `cd8 tem` `lung` cells for
females older than `80` with `COVID-19`. Note: Use the harmonized columns, where
possible. 

:::::::::::::: solution


```r
metadata |> 
    dplyr::filter(
        sex == "female" &
        age_days > 80 * 365 &
        stringr::str_like(assay, "%10x%") &
        disease == "COVID-19" &  
        tissue_harmonised == "lung" & 
        cell_type_harmonised == "cd8 tem"
    ) |>
    get_single_cell_experiment()
```

```{.output}
ℹ Realising metadata.
```

```{.output}
ℹ Synchronising files
```

```{.output}
ℹ Downloading 2 files, totalling 0 GB
```

```{.output}
ℹ Downloading https://object-store.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/cellxgene-0.2.1-hdf5/original/893d8537e318769108b4962020ddd846/assays.h5 to /home/runner/.cache/R/CuratedAtlasQueryR/0.2.1/original/893d8537e318769108b4962020ddd846/assays.h5
```

```{.output}
ℹ Downloading https://object-store.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/cellxgene-0.2.1-hdf5/original/893d8537e318769108b4962020ddd846/se.rds to /home/runner/.cache/R/CuratedAtlasQueryR/0.2.1/original/893d8537e318769108b4962020ddd846/se.rds
```

```{.output}
ℹ Reading files.
```

```{.output}
ℹ Compiling Single Cell Experiment.
```

```{.output}
class: SingleCellExperiment 
dim: 36229 12 
metadata(0):
assays(1): counts
rownames(36229): A1BG A1BG-AS1 ... ZZEF1 ZZZ3
rowData names(0):
colnames(12): TCATCATCATAACCCA_1 TATCTGTCAGAACCGA_1 ...
  CCCTTAGCATGACTTG_1 CAGTTCCGTAGCGTAG_1
colData names(56): sample_ cell_type ... updated_at_y original_cell_id
reducedDimNames(0):
mainExpName: NULL
altExpNames(0):
```
:::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: keypoints 

- TODO

::::::::::::::::::::::::::::::::::::::::::::::::

# Acknowledgements

Thank you to [Stefano Mangiola](https://github.com/stemangiola) and his team for
developing
[CuratedAtlasQueryR](https://github.com/stemangiola/CuratedAtlasQueryR) and
graciously providing the content from their vignette. Make sure to keep an eye
out for their publication for proper citation. Their bioRxiv paper can be found
at <https://www.biorxiv.org/content/10.1101/2023.06.08.542671v1>.
