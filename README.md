# virolink

`virolink` is an R package for visualizing viral integration events in host
genomes. It helps users explore where viral sequences integrate, how integration
patterns differ across samples or user-defined groups, and how host breakpoints connect to
viral genomic positions.

The package focuses on visualization and exploratory evidence views. It does
not perform read alignment, integration calling, or automatic biological
interpretation.

## Installation

Install the development version from GitHub:

```r
remotes::install_github("YuLab-SMU/virolink")
```

Load the package as:

```r
library(virolink)
```

## What virolink provides

`virolink` organizes viral-integration visualization around a practical
analysis workflow:

- standardize integration records into a consistent table structure;
- define host and virus genome coordinates explicitly;
- summarize integration burden across samples, cohorts, or user-defined groups;
- survey integration events across the host genome;
- compare host breakpoints with viral breakpoint positions;
- inspect local host genomic context around selected integration sites;
- visualize host-virus links in circular or linear layouts;
- build lightweight interactive HTML views for exploration.

## Input data

The core input is a tabular integration table with one row per integration
event. The standard columns are:

| Column | Meaning |
| --- | --- |
| `host_chr` | host chromosome |
| `host_pos` | host breakpoint position |
| `host_strand` | host-side breakpoint strand |
| `virus_chr` | virus sequence |
| `virus_pos` | viral breakpoint position |
| `virus_strand` | viral strand orientation |
| `support_reads` | number of reads supporting the integration event |
| `sample` | optional sample or group identifier |

Caller-specific importer functions can convert supported upstream result files
into this standard structure. Optional metadata such as `sample` is added only
when it is present in the input or supplied explicitly; it is not inferred from
file names or paths. Additional user-supplied columns are preserved as ordinary
metadata and can be used explicitly for grouping or coloring.

## Main user-facing functions

| Task | Functions |
| --- | --- |
| Standardize records | `read_integrations()`, `as_integrations()`, `validate_integrations()` |
| Import caller outputs | `read_caller_integrations()`, `read_ctat_vif_integrations()`, `read_virusfinder2_integrations()`, `read_virus_clip_integrations()`, `read_bsvf_integrations()` |
| Define genomes and annotations | `host_genome()`, `virus_genome()`, `virus_features()`, `hbv_features()`, `host_features()` |
| Compare samples or cohorts | `plot_cohort_comparison()` |
| Survey host-genome distribution | `plot_integration_landscape()` |
| Compare host and virus breakpoints | `plot_breakpoint_map()` |
| Inspect local host context | `plot_locus_context()` |
| Summarize signature or clonality evidence | `plot_integration_signature()` |
| Build circular linked views | `plot_integrations()` with `track_*()` functions |
| Build linear linked views | `plot_linear_integrations()` |
| Explore interactively | `plot_interactive_explorer()` |

## Documentation

The runnable workflow and demonstration plots are provided in the package
vignette:

```r
vignette("virolink", package = "virolink")
```

The vignette uses the small example files bundled under `inst/extdata`, so it
can be run without downloading external genome resources.

## Development status

`virolink` is under active development. The current API is centered on
user-friendly visualization workflows for viral integration data, and additional
views may be added as the package evolves.
