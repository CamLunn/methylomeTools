# Changes

## 0.1.0.9002

- Removed extra empty arguments from help usage and function-reference signatures.
- Initialized PSOCK worker library paths and namespaces explicitly.
- Included unit-specific errors and log paths in failed DSS run messages.
- Forwarded live DSS output-compression and memory-log settings to each unit.
- Used portable character escapes in R source and registered data-mask columns
  for static code analysis.
- Removed an empty argument section from the no-argument database helper.

## 0.1.0.9001

- Converted selected beta values to log2 M values when `value_type = "beta"`.
  Results differ from older beta-input runs; older beta-input caches are rejected.
- Corrected SVA worker reporting. Preprocessing uses one R process; `workers`
  values above one warn. BLAS threading remains separate.
- Honored explicit BLAS thread counts in SLURM arguments, the startup environment,
  and the generated R driver. Thread counts above the CPU allocation are rejected.
- Added errors for rank-deficient or nonnested SVA designs, fractional manual SV
  counts, and SV counts that exhaust residual degrees of freedom.
- Disabled the legacy `two_step` option, which called `svaseq` on methylation
  values. The supported estimator is `sva::sva()` with `sva_method = "irw"`.
- Added regression tests for SVA scale handling, designs, and SLURM thread settings.
- Repaired help-page color defaults and synchronized documented signatures with
  function definitions.

## 0.1.0.9000

- Organized methylation preparation, SVA, DSS, comb-p, annotation, database access,
  and plotting functions into an R package.
- Added DSS PSOCK and fork backends, output configuration guards, and explicit
  resume behavior.
- Added sample-order checks and support for arbitrary unique feature IDs in SVA.
- Added `make_dmr_site_table()` to connect region results to browser inputs.
- Added array-to-CpG mapping with multiple probes, optional sample-name mapping,
  and explicit missing-value handling.
- Added optional installation-library support, MIT licensing, and citation files.

## Compatibility aliases

| Previous name | Current entry point |
|---|---|
| `run_dss_multifactor_by_chr()` | `run_dss_parallel()` |
| `make_truseq_mats_by_mincov_from_base()` | `make_coverage_threshold_matrices()` |
| `process_beta_EPIC()` | `align_array_beta_to_cpgs()` |

The historical `process_beta_EPIC()` wrapper retains its argument names and
`na.rm = FALSE` behavior. Available mapped probes are averaged; unmatched
mappings remain missing.
