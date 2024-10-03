---
title: Accessing data from the Human Cell Atlas (HCA)
teaching: 20 # Minutes of teaching in the lesson
exercises: 10 # Minutes of exercises in the lesson
---

:::::::::::::::::::::::::::::::::::::: questions 

- How to obtain single-cell reference maps from the Human Cell Atlas?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Learn about different resources for public single-cell RNA-seq data.
- Access data from the Human Cell Atlas using the `CuratedAtlasQueryR` package.
- Query for cells of interest and download them into a `SingleCellExperiment` object. 

::::::::::::::::::::::::::::::::::::::::::::::::


# Single Cell data sources

## HCA Project

The Human Cell Atlas (HCA) is a large project that aims to learn from and map
every cell type in the human body. The project extracts spatial and molecular
characteristics in order to understand cellular function and networks. It is an
international collaborative that charts healthy cells in the human body at all
ages. There are about 37.2 trillion cells in the human body. To read more about
the project, head over to their website at https://www.humancellatlas.org.

## CELLxGENE

CELLxGENE is a database and a suite of tools that help scientists to find,
download, explore, analyze, annotate, and publish single cell data. It includes
several analytic and visualization tools to help you to discover single cell
data patterns. To see the list of tools, browse to
https://cellxgene.cziscience.com/.

## CELLxGENE | Census

The Census provides efficient computational tooling to access, query, and
analyze all single-cell RNA data from CZ CELLxGENE Discover. Using a new access
paradigm of cell-based slicing and querying, you can interact with the data
through TileDB-SOMA, or get slices in AnnData or Seurat objects, thus
accelerating your research by significantly minimizing data harmonization at
https://chanzuckerberg.github.io/cellxgene-census/.

## The CuratedAtlasQueryR Project

The `CuratedAtlasQueryR` is an alternative package that can also be used to access the CELLxGENE data from R through a tidy API. The data has also been harmonized, curated, and re-annotated across studies.

`CuratedAtlasQueryR` supports data access and programmatic exploration of the
harmonized atlas. Cells of interest can be selected based on ontology, tissue of
origin, demographics, and disease. For example, the user can select CD4 T helper
cells across healthy and diseased lymphoid tissue. The data for the selected
cells can be downloaded locally into SingleCellExperiment objects. Pseudo
bulk counts are also available to facilitate large-scale, summary analyses of
transcriptional profiles. 

<img src="https://raw.githubusercontent.com/ccb-hms/osca-workbench/main/episodes/figures/curatedAtlasQuery.png" style="display: block; margin: auto;" />

## Data Sources in R / Bioconductor

There are a few options to access single cell data with R / Bioconductor.

| Package | Target | Description |
|---------|-------------|---------|
| [hca](https://bioconductor.org/packages/hca) | [HCA Data Portal API](https://www.humancellatlas.org/data-portal/) | Project, Sample, and File level HCA data |
| [cellxgenedp](https://bioconductor.org/packages/cellxgenedp) | [CellxGene](https://cellxgene.cziscience.com/) | Human and mouse SC data including HCA |
| [CuratedAtlasQueryR](https://stemangiola.github.io/CuratedAtlasQueryR/) | [CellxGene](https://cellxgene.cziscience.com/) | fine-grained query capable CELLxGENE data including HCA |

## Installation


``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("CuratedAtlasQueryR")
```

## Package load 




``` r
library(CuratedAtlasQueryR)
library(dplyr)
```

## HCA Metadata

The metadata allows the user to get a lay of the land of what is available
via the package. In this example, we are using the sample database URL which
allows us to get a small and quick subset of the available metadata.


``` r
metadata <- get_metadata(remote_url = CuratedAtlasQueryR::SAMPLE_DATABASE_URL) |> 
  collect()
```

Get a view of the first 10 columns in the metadata with `glimpse()`


``` r
metadata |>
  select(1:10) |>
  glimpse()
```

``` output
Rows: ??
Columns: 10
Database: DuckDB v0.10.2 [unknown@Linux 6.5.0-1021-azure:R 4.4.0/:memory:]
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

## A tangent on the pipe operator

The vignette materials provided by `CuratedAtlasQueryR` show the use of the
'native' R pipe (implemented after R version `4.1.0`). For those not familiar
with the pipe operator (`|>`), it allows you to chain functions by passing the
left-hand side as the first argument to the function on the right-hand side. It is used extensively in the `tidyverse` dialect of R, especially within the [`dplyr` package](https://dplyr.tidyverse.org/).

The pipe operator can be read as "and then". Thankfully, R doesn't care about whitespace, so it's common to start a new line after a pipe. Together these points enable users to "chain" complex sequences of commands into readable blocks.

In this example, we start with the built-in `mtcars` dataset and then filter to rows where `cyl` is not equal to 4, and then compute the mean `disp` value by each unique `cyl` value.


``` r
mtcars |> 
  filter(cyl != 4) |> 
  summarise(avg_disp = mean(disp),
            .by = cyl)
```

``` output
  cyl avg_disp
1   6 183.3143
2   8 353.1000
```

This command is equivalent to the following:


``` r
summarise(filter(mtcars, cyl != 4), mean_disp = mean(disp), .by = cyl)
```

## Exploring the metadata

Let's examine the metadata to understand what information it contains.

We can tally the tissue types across datasets to see what tissues the experimental data come from:


``` r
metadata |>
  distinct(tissue, dataset_id) |> 
  count(tissue) |> 
  arrange(-n)
```

``` output
# A tibble: 33 × 2
   tissue                   n
   <chr>                <int>
 1 blood                   17
 2 kidney                   8
 3 cortex of kidney         7
 4 heart left ventricle     7
 5 renal medulla            6
 6 respiratory airway       6
 7 bone marrow              4
 8 kidney blood vessel      4
 9 lung                     4
10 renal pelvis             4
# ℹ 23 more rows
```

We can do the same for the assay types:


``` r
metadata |>
    distinct(assay, dataset_id) |>
    count(assay)
```

``` output
# A tibble: 12 × 2
   assay                              n
   <chr>                          <int>
 1 10x 3' v1                          1
 2 10x 3' v2                         27
 3 10x 3' v3                         21
 4 10x 5' v1                          7
 5 10x 5' v2                          2
 6 Drop-seq                           1
 7 Seq-Well                           2
 8 Slide-seq                          4
 9 Smart-seq2                         1
10 Visium Spatial Gene Expression     7
11 scRNA-seq                          4
12 sci-RNA-seq                        1
```

:::: challenge

Look through the full list of metadata column names. Do any other metadata
columns jump out as interesting to you for your work?


``` r
names(metadata)
```

::::

## Downloading single cell data 

The data can be provided as either "counts" or counts per million "cpm" as given
by the `assays` argument in the `get_single_cell_experiment()` function. By
default, the `SingleCellExperiment` provided will contain only the 'counts'
data.

For the sake of demonstration, we'll focus this small subset of samples. We use the `filter()` function from the `dplyr` package to identify cells meeting the following criteria:

* African ethnicity
* 10x assay
* lung parenchyma tissue
* CD4 cells


``` r
sample_subset <- metadata |>
    filter(
        ethnicity == "African" &
        grepl("10x", assay) &
        tissue == "lung parenchyma" &
        grepl("CD4", cell_type)
    )
```

Out of the 111355 cells in the sample database, 1571 cells meet this criteria.

Now we can use `get_single_cell_experiment()`:


``` r
single_cell_counts <- sample_subset |>
    get_single_cell_experiment()

single_cell_counts
```

``` output
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

You can provide different arguments to `get_single_cell_experiment()` to get different formats or subsets of the data, like data scaled to counts per million:


``` r
sample_subset |>
  get_single_cell_experiment(assays = "cpm")
```

``` output
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

or data on only specific genes:


``` r
single_cell_counts <- sample_subset |>
    get_single_cell_experiment(assays = "cpm", features = "PUM1")

single_cell_counts
```

``` output
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

Or if needed, the H5 `SingleCellExperiment` can be returned a Seurat
object (note that this may take a long time and use a lot of memory depending on
how many cells you are requesting).


``` r
single_cell_counts <- sample_subset |>
    get_seurat()

single_cell_counts
```

## Save your `SingleCellExperiment`

Once you have a dataset you're happy with, you'll probably want to save it. The recommended way of saving these `SingleCellExperiment` objects is to use
`saveHDF5SummarizedExperiment` from the `HDF5Array` package.


``` r
single_cell_counts |> saveHDF5SummarizedExperiment("single_cell_counts")
```

## Exercises

:::::::::::::::::::::::::::::::::: challenge

#### Exercise 1: Basic counting + piping

Use `count` and `arrange` to get the number of cells per tissue in descending
order.

:::::::::::::: solution


``` r
metadata |>
    count(tissue) |>
    arrange(-n)
```
:::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::: challenge

#### Exercise 2: Tissue & type counting

`count()` can group by multiple factors by simply adding another grouping column
as an additional argument. Get a tally of the highest number of cell types per
tissue combination. What tissue has the most numerous type of cells?

:::::::::::::: solution


``` r
metadata |>
    count(tissue, cell_type) |>
    arrange(-n) |> 
    head(n = 1)
```
:::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::: challenge

#### Exercise 3: Comparing metadata categories

Spot some differences between the `tissue` and `tissue_harmonised` columns.
Use `count` to summarise.

:::::::::::::: solution


``` r
metadata |>
    count(tissue) |>
    arrange(-n)

metadata |>
    count(tissue_harmonised) |>
    arrange(-n)
```

For example you can see that `tissue_harmonised` merges the `cortex of kidney`
and `kidney` groups in `tissue`.

To see the full list of curated columns in the metadata, see the Details section
in the `?get_metadata` documentation page.
    
:::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::: challenge

#### Exercise 4: Highly specific cell groups

Now that we are a little familiar with navigating the metadata, let's obtain
a `SingleCellExperiment` of 10X scRNA-seq counts of `cd8 tem` `lung` cells for
females older than `80` with `COVID-19`. Note: Use the harmonized columns, where
possible. 

:::::::::::::: solution


``` r
metadata |> 
    filter(
        sex == "female" &
        age_days > 80 * 365 &
        grepl("10x", assay) &
        disease == "COVID-19" &  
        tissue_harmonised == "lung" & 
        cell_type_harmonised == "cd8 tem"
    ) |>
    get_single_cell_experiment()
```

``` output
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

You can see we don't get very many cells given the strict set of conditions we used.
:::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: keypoints 

- The `CuratedAtlasQueryR` package provides programmatic access to single-cell reference maps from the Human Cell Atlas.
- The package provides functionality to query for cells of interest and to download them into a `SingleCellExperiment` object.

::::::::::::::::::::::::::::::::::::::::::::::::


