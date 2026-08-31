# virolink 0.1.7

- Added caller import helpers for converting CTAT-VIF and VirusFinder2 outputs
  into package-standard `vi_integrations` tables.
- Updated the development `vi_integrations` schema to use explicit host and
  virus breakpoint sides with `support_reads` as the read-support field.
- Removed `method` and `caller` as standardized integration-record fields;
  user-supplied grouping metadata is preserved as ordinary extra columns and
  must be referenced explicitly.
- Renamed circular track color mappings from method-specific `method_col` to
  generic `color_map`.
- Documented that importer helpers do not infer optional metadata such as
  `sample` from file names or paths.
- Added a top-level `README.md` describing the package scope, input format,
  main workflow, documentation entry point, and function map.

# virolink 0.1.6

- Fixed duplicate circular plots from `plot_integrations()` in Quarto/R
  Markdown output by drawing returned `vi_integration_plot` objects only once.
- Reused the same drawing path for direct plotting and printing so method
  legends remain consistent.
- Improved circular track labels and virus feature labels with smaller default
  text, softer colors, and optional width filtering for dense virus annotations.
- Added ideogram text-size controls and an `axis` option for circular plots;
  by default, genomic axis labels are now shown only on the virus sector to
  reduce host-chromosome label overlap.
- Updated the circular linked-view vignette example with shorter track labels.

# virolink 0.1.5

- Improved `plot_locus_context()` host feature annotations with softer
  background-like fills and labels positioned above annotation intervals.
- Added a separate host feature type legend so annotation colors are not
  confused with integration point colors.
- Added `show_feature_legend` to optionally hide the host feature type legend.

# virolink 0.1.4

- Added `plot_interactive_explorer()` for lightweight interactive exploration
  with a `ggiraph` host-virus breakpoint map and a searchable `DT` integration
  table.
- Added hover tooltips, point selection IDs, optional color and size mappings,
  and configurable table columns for interactive reports.
- Updated the vignette with a runnable interactive exploration example.

# virolink 0.1.3

- Updated oncoprint assembly to use the exported
  `aplotExtra::oncoplot_components()` interface instead of unexported
  `aplotExtra` internals.

# virolink 0.1.2

- Added `plot_integration_signature()` for visualizing candidate clone-like
  breakpoint clusters and viral-strand orientation signatures.
- Added host-only and host-virus proximity rules for grouping recurrent
  breakpoint evidence in the clonality view.
- Updated the vignette with a runnable mechanism-and-clonality evidence view.

# virolink 0.1.1

- Added `plot_cohort_comparison()` for comparing integration burden across
  samples, cohorts, methods, or other user-supplied grouping columns.
- Added support for event counts, total supporting reads, stacked composition,
  within-group proportions, and optional collapsing of low-frequency
  categories.
- Updated the vignette with a runnable cohort-comparison example.

# virolink 0.1.0

- Added `host_features()` to standardize local host annotation intervals for
  locus-level visualization.
- Added `plot_locus_context()` to inspect integration evidence around a
  selected host breakpoint with optional host feature annotations.
- Added a small demo host-feature table and a runnable locus-context example
  to the vignette workflow.

# virolink 0.0.9

- Added `plot_breakpoint_map()` to compare host-genome breakpoint positions
  with viral-genome breakpoint positions in a single two-dimensional view.
- Added optional virus feature bands, group-based point colors, and numeric
  point-size mapping for breakpoint maps.
- Updated the vignette with a runnable breakpoint-map example in the main
  viral-integration visualization workflow.

# virolink 0.0.8

- Added `plot_integration_landscape()` for a host-genome overview of viral
  integration sites or binned event density.
- Added optional grouping and color mapping to compare integration landscapes
  across samples, methods, or other columns.
- Reorganized the main vignette around the viral-integration analysis workflow
  and added a roadmap for planned visualization capabilities.
- Corrected `Authors@R` to identify a single package maintainer.

# virolink 0.0.7

- Added generic virus feature support for circos and linear plots:
  - `virus_features()` to standardize viral feature interval tables.
  - `hbv_features()` as a built-in helper for HBV genome annotation.
- Added two new circos tracks for virus-specific layers:
  - `track_virus_genes()` draws viral feature intervals on the virus sector.
  - `track_virus_density()` draws event-count density over `virus_pos`.
- Added `plot_linear_integrations()` for an independent non-circos view
  of host-virus integration positions with optional virus feature tracks.
- Fixed `R CMD check` NOTEs: double-escaped `\code` in roxygen, missing
  `importFrom(utils, tail)`.
- Changed `oncoplot()` example from `\donttest` to `\dontrun` to avoid
  optional-dependency failures in CI.

# ViHostPlot 0.0.6

- Added a new user-facing API for integration visualization:
  - `read_integrations()`, `as_integrations()`, and `validate_integrations()` for standardized input handling
  - `host_genome()` and `virus_genome()` for explicit genome objects
  - `track_ideogram()`, `track_sites()`, `track_density()`, and `track_links()` for declarative track construction
  - `plot_integrations()` as the main integration plotting entry point
  - `plot_oncoprint()` as the renamed oncoprint entry point
- Kept the legacy `visualize_viral_integration()` and `oncoplot()` interfaces as compatibility wrappers.
- Updated the example integration tables in `inst/extdata/` to the new column convention:
  - `host_chr`, `host_pos`, `virus_pos`, `support`, `sample`, `virus_strand`, `method`
- Standardized the package documentation and namespace around the new API.

# ViHostPlot 0.0.5

- Refined clinical oncoplot behavior and chromosome handling.
- Improved host genome resolution and circos plotting stability.

# ViHostPlot 0.0.4

- Refactored host genome handling for circos plots.
  - `host` now supports common built-in host names.
  - `chrom_file` supports custom chromosome-size tables.
- Moved host genome utilities into `host_chrom_utils.R` so genome resolution and plotting logic are easier to maintain separately.
- Changed virus handling to require explicit `virus_name` and `virus_length` instead of embedding the virus sequence in the chromosome-size table.
- Added UCSC-based host chromosome-size retrieval, so names such as `human` can resolve automatically to `hg38`.
- Removed `color_file` and `palettes`.
- `chrom_file` now only needs host chromosome information.
- Chromosome names are normalized by removing a `chr` prefix when needed, for example `chr1` becomes `1`.

# ViHostPlot 0.0.3

- Fixed a bug where some oncoplot titles were not shown during rendering.
- Updated the virus sequence lookup logic in the circos workflow so the virus name must be supplied explicitly.
- Added `histogram`-style genomic tracks to the circos plot workflow, allowing density views in addition to scatter views.

# ViHostPlot 0.0.2

- Removed the intermediate VCF generation step. `visualize_viral_integration()` now converts the input table directly into a `GInteractions` object through the new `create_gi_from_table()` function.
- Moved the core plotting engine to `circlize` genomic helper functions, including `circos.genomicInitialize()`, `circos.genomicTrackPlotRegion()`, `circos.genomicPoints()`, and `circos.genomicLink()`, for more stable genomic coordinate handling.
- Reduced heavy dependencies by removing several Bioconductor packages.
- Reduced tidyverse dependencies by removing `dplyr` and `stringr`.
- Added `create_host_virus_layout()` to generate a standard track layout list quickly.
- Added `match_chr_style()` to normalize chromosome naming styles, for example by matching `chr1` and `1`.
- Added `make_axis_labels()` to auto-format ideogram axis units, using Mb for host chromosomes and kb for the virus sequence.
- Added `validate_layout_list()` to catch invalid track definitions before plotting starts.

# ViHostPlot 0.0.1

- Reworked the `oncoplot()` layout and strengthened clinical annotation tracks.
- Removed `patchwork` / `ggplotify` dependencies and used native `aplot` layout composition.
- Simplified the mutation classification logic by replacing the old `discrete` and `continuous` labels with `binning_numeric()` for automatic binning and labeling of continuous variables.
- Exposed `n_breaks` and `colors` to improve customization of clinical annotation tracks.
