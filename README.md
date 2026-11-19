# methylomeTools

> **Development status:** Function testing, unit testing, and end-to-end testing are ongoing. This package is under active development.

R tools for methylation matrix preparation, surrogate variable analysis,
parallel DSS modeling, and CpG and region visualization.

The package brings these steps into a reusable workflow for sequencing-based
methylation studies. Array-to-CpG mapping and coverage filtering are also included.

## Installation

From the extracted package directory:

```sh
Rscript inst/examples/install_dependencies.R
R CMD INSTALL .
```

For a shared HPC library, see [HPC setup](docs/HPC.md). R >= 4.1.0 is the syntax
minimum; dependencies must also be compatible with your R and Bioconductor
versions. Cross-version support has not been established.

```r
library(methylomeTools)
help(package = "methylomeTools")
```

## Main functions

| Task | Functions | Input |
|---|---|---|
| Coverage filtering | `mask_matrices_by_coverage()`, `make_coverage_threshold_matrices()` | Coverage and optional methylation matrices |
| SVA | `prepare_methyl_matrix()`, `run_sva_on_matrix()`, `build_methyl_svs()` | Feature-by-sample beta or M-value matrix and aligned designs |
| DSS | `run_dss_parallel()` | BSseq read counts and sample metadata |
| Site plots | `plot_dss_manhattan()`, `plot_top_dml_betas_scatter()`, `plot_top_dml_betas_boxplot()` | DSS results; CpG plots also require BSseq and metadata |
| Region calling | `run_combp_on_dss()` | Genomic positions and site-level P values |
| Region plots | `make_dmr_site_table()`, `make_dmr_methylation_long()`, `plot_dmr_browser()` | Region coordinates, BSseq, and metadata |
| Annotation | `annotate_hg38_dmrs()`, CpG database helpers | Human hg38 regions or a compatible user-supplied SQLite database |

Matrices use features in rows and samples in columns. DSS requires methylated
and total read counts, not beta values alone. The BSseq assay named `M` contains
methylated counts; this differs from the log2 methylation M values used for SVA.
No patient data or lab databases are distributed.

## SVA example

Both model matrices use the same sample order as the methylation matrix.
Missing covariates can cause `model.matrix()` to omit rows.

```r
metadata <- metadata[match(colnames(beta_matrix), metadata$sample_id), , drop = FALSE]
stopifnot(identical(metadata$sample_id, colnames(beta_matrix)))
rownames(metadata) <- metadata$sample_id

mod <- model.matrix(~ exposure + age + sex + batch, data = metadata)
mod0 <- model.matrix(~ age + sex + batch, data = metadata)

sv <- build_methyl_svs(
  meth_mat = beta_matrix,
  mod = mod,
  mod0 = mod0,
  value_type = "beta",
  feature_strategy = "hvf_only",
  n_hvf = 100000,
  num_sv_method = "be",
  sva_method = "irw",
  workers = 1,
  blas_threads = 4
)

# Join by sample ID rather than assuming the returned row order.
idx <- match(metadata$sample_id, sv$sv$sample_id)
stopifnot(!anyNA(idx))
sv_names <- grep("^SV[0-9]+$", names(sv$sv), value = TRUE)
metadata[sv_names] <- sv$sv[idx, sv_names, drop = FALSE]
```

The feature count and thread count above are examples, not universal settings.
Selected beta values are converted to `log2(beta / (1 - beta))`, clipping endpoints
at 1e-6. Feature ranking uses the supplied scale. Supply `value_type = "M"` for
already transformed data. With `impute = FALSE`, selected features with missing
values are removed. Imputation is not implemented.

**Scale handling:** Older snapshots accepted `value_type = "beta"` without doing
the conversion. Beta-input results therefore differ from those snapshots, and
older beta-input caches are rejected. Already transformed inputs supplied as
`value_type = "M"` are used directly.

### Parallelism

The [BLAS setup guide](docs/BLAS.md) covers library discovery, thread settings,
external OpenBLAS configuration, and verification on other systems.

SVA runs in one R process. `workers` is retained for compatibility and does not
start preprocessing workers; values above one now warn. `blas_threads` requests
threads for matrix operations when R uses a suitable BLAS library. It does not
parallelize the estimator's iterations or guarantee a speedup. Setting environment
variables after R starts may be insufficient without `RhpcBLASctl`.

For batch jobs, the SLURM helper sets thread variables before launching R. It
honors an explicit `blas_threads` value and otherwise uses `cpus_per_task`.
Actual thread use depends on the installed BLAS implementation and CPU
allocation. Splitting chromosomes produces separate latent-variable estimates
rather than one joint SVA fit.

Prepare features once with `prepare_methyl_matrix()` and reuse the result with
`run_sva_on_matrix()` for models on the same samples. This avoids repeating
preprocessing. Independent models can be scheduled as separate jobs, with enough
memory for each matrix copy.

DSS supports process parallelism separately, using `n_workers` and
`parallel_backend`. Start with a small worker count because PSOCK workers can
copy count matrices.

## Examples and documentation

- [Analysis workflow](docs/WORKFLOW.md)
- [HPC and SLURM setup](docs/HPC.md)
- [BLAS libraries and threading](docs/BLAS.md)
- [Function signatures](docs/FUNCTION_REFERENCE.md)
- [Changes](docs/CHANGES.md)

A synthetic example is installed with the package:

```r
source(system.file("examples", "quickstart.R", package = "methylomeTools"))
```

## Current limitations

- The legacy `sva_method = "two_step"` setting is rejected. It previously called
  the sequencing estimator `svaseq`, rather than the named two-step estimator.
  Use `irw` for methylation M values.
- Plotting functions show observed beta values. They do not adjust those values
  for the covariates used in DSS.
- comb-p depends on ENmix, and generated SLURM scripts depend on the target
  cluster configuration.
- Annotation is human hg38-specific. Optional DMR tracks require compatible RDS
  annotation tables, and database helpers require a compatible SQLite schema.

## Citation and license

MIT license; see [LICENSE.md](LICENSE.md). Maintainer: Cam Lunn
<cameron.lunn@nyulangone.org>.

Use `citation("methylomeTools")` for software citation information and cite the
underlying methods used in your analysis, including sva, DSS, and ENmix where
applicable. See [citation guidance](docs/CITATION_GUIDE.md).
