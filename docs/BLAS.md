# BLAS libraries and SVA threading

BLAS and LAPACK provide matrix operations used by R and SVA. A threaded
implementation can accelerate some operations, but the benefit depends on matrix
size, the BLAS build, and the time spent in other parts of the estimator.

In methylomeTools, `workers` does not create SVA preprocessing processes.
`blas_threads` requests threads within one R process. DSS process parallelism is
controlled separately with `n_workers`.

## Identify the library used by R

Run in the R session that will perform the analysis:

```r
sessionInfo()[c("BLAS", "LAPACK")]

blas_path <- sessionInfo()$BLAS
if (!is.null(blas_path) && nzchar(blas_path)) {
  normalizePath(blas_path, mustWork = FALSE)
}

Sys.getenv(c("SLURM_CPUS_PER_TASK", "OPENBLAS_NUM_THREADS",
             "OMP_NUM_THREADS", "MKL_NUM_THREADS"))

if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  list(
    blas_available_cores = RhpcBLASctl::blas_get_num_procs(),
    omp_max_threads = RhpcBLASctl::omp_get_max_threads()
  )
}
```

`libRblas.so` is a common R library name. It may be R's reference BLAS, a symlink,
or another build, so the name alone does not establish threading capability.
`normalizePath()` resolves a symlink when the file is accessible. BLAS and LAPACK
may resolve to different files or to the same external library.

An empty thread environment variable means no value was explicitly set there;
libraries may still use defaults. `blas_get_num_procs()` reports available cores
for a recognized BLAS and can return one when unknown. `omp_get_max_threads()`
reports an OpenMP setting, not observed BLAS activity. Requested counts and these
queries do not replace elapsed-time measurements or CPU-usage observations.

## Find available libraries on an HPC

Environment module names and installation locations differ across institutions.
In a terminal with the module system initialized:

```sh
module -t avail 2>&1 | grep -Ei 'openblas|mkl|flexiblas|(^|/)r/'
module show OPENBLAS_MODULE 2>&1
module show R_MODULE 2>&1
```

`OPENBLAS_MODULE` and `R_MODULE` stand for names listed by the local module system.
Module definitions often expose library directories through `LD_LIBRARY_PATH`,
`LIBRARY_PATH`, `BLASDIR`, or similar variables. Loading a module usually changes
environment variables; it does not necessarily switch a separately built R
installation to that BLAS.

For a directory reported by the module definition, Linux commands such as these
list installed libraries and their dependencies:

```sh
ls -ld /path/to/openblas
find /path/to/openblas -maxdepth 5 \( -name '*.so*' -o -name '*.a' \) -print
ldd /absolute/path/to/libopenblas.so
```

A `.so` file is a shared library. A `.a` file is a static archive and cannot be
used as an `LD_PRELOAD` target. `not found` entries from `ldd` indicate unresolved
dependencies. Module definitions can contain stale paths, and login and compute
nodes can expose different filesystems. The actual library and its dependencies
must be accessible on the node running R.

## Inspect an OpenBLAS shared library

The package includes a Python 3 utility that loads one library in a separate
process and reports its OpenBLAS build configuration:

```r
system.file("examples", "inspect_openblas.py", package = "methylomeTools")
```

Use the returned script path in a terminal:

```sh
python3 /path/to/inspect_openblas.py /absolute/path/to/libopenblas.so
```

The report includes:

| Field | Meaning |
|---|---|
| `threading_mode = 0` | Serial OpenBLAS build |
| `threading_mode = 1` | Pthread-based multithreading |
| `threading_mode = 2` | OpenMP-based multithreading |
| `USE64BITINT` in configuration | 64-bit BLAS integer interface; requires an explicitly compatible R build |
| `MAX_THREADS` when reported | Build limit, not a scheduler allocation or observed thread count |
| CPU name or `DYNAMIC_ARCH` | Target architecture or runtime CPU selection |

The processor target and integer interface must match the R environment. A file
copied from a different machine is not interchangeable solely because its name
contains OpenBLAS. Cluster-supported R/BLAS combinations provide the preferred
starting point. A provider's Intel MKL configuration is best used through its
supported R installation or module rather than through the OpenBLAS argument.

## Activate threading

When R already links to a suitable threaded BLAS, `openblas_so = NULL` is
appropriate. Set the thread environment before R starts; the package also uses
RhpcBLASctl when available to request thread counts at runtime.

For OpenBLAS, pthread builds use `OPENBLAS_NUM_THREADS`; OpenMP builds use
`OMP_NUM_THREADS`. The latter also affects other OpenMP code. Thread counts must
fit within the allocated CPUs. A SLURM task with four CPUs supports a request of
four threads, not four processes each requesting four threads.

For a compatible external OpenBLAS library on Linux, a process-local preload can
select the shared library before R starts:

```sh
env LD_PRELOAD=/absolute/path/to/libopenblas.so \
    OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 \
    Rscript analysis.R
```

This applies only to that command and does not edit a shared R installation.
It cannot change BLAS in an already-running RStudio session. A library directory
is not a preload target; `LD_PRELOAD` takes the shared-library filename. Any
additional dependencies also need to be resolvable in that process.

macOS and Windows have different loader and R-distribution configurations. The
Linux preload example does not apply to them; an R distribution configured with
a suitable BLAS provides the appropriate route on those platforms.

## A small numerical check

With the same startup environment as the intended analysis, this R code checks
matrix multiplication against scalar dot products and reconstructs an SVD:

```r
print(sessionInfo()[c("BLAS", "LAPACK")])

set.seed(1)
x <- matrix(rnorm(6000), nrow = 200)
reference <- vapply(seq_len(ncol(x)), function(j) {
  vapply(seq_len(ncol(x)), function(i) {
    sum(x[, i] * x[, j])
  }, numeric(1))
}, numeric(ncol(x)))
stopifnot(max(abs(crossprod(x) - reference)) < 1e-10)

s <- svd(x)
reconstructed <- sweep(s$u, 2, s$d, "*") %*% t(s$v)
stopifnot(max(abs(x - reconstructed)) < 1e-10)
```

Successful execution establishes these small numerical checks in that process;
it does not establish a performance improvement or validate every analysis.

## SLURM submission with methylomeTools

The helper can preload an external OpenBLAS file and set thread variables before
launching the batch R process:

```r
blas_library <- "/absolute/path/to/libopenblas.so"

job <- submit_methyl_svs_slurm(
  meth_mat = m_matrix,
  mod = mod,
  mod0 = mod0,
  value_type = "M",
  feature_strategy = "hvf_only",
  n_hvf = 100000,
  workers = 1,
  cpus_per_task = 4,
  blas_threads = 4,
  openblas_so = blas_library,
  r_module = NULL,
  outdir = file.path(getwd(), "sva_run"),
  submit = TRUE
)
```

`r_module = NULL` uses the R executable available in the batch environment. A
cluster requiring module loading can supply its R module name. The package R
library, OpenBLAS file, and output directory need to be visible to compute nodes.
An explicit `blas_threads` count cannot exceed `cpus_per_task`; otherwise the
thread count defaults to the CPU request. Runtime logs report the BLAS path and
requested thread count.

## Performance comparisons

A useful comparison separates three effects:

- R's default BLAS with one requested thread;
- the external OpenBLAS library with one requested thread;
- the same external library with several requested thread counts.

The second comparison separates a change of implementation from a change in
thread count. Repeated trials in fresh R processes, using the same input matrix,
models, selected features, seed, and SV count, give comparable measurements.
SVA preparation and fitting times describe different work and can scale
differently. Automatic SV-count estimation adds a separate workload; a fixed
manual count isolates comparisons with the same number of fitted SVs.

Small matrices may be slower with additional threads. CPU allocation and thread
settings also affect resource use, so the largest request is not necessarily the
fastest. Numerical agreement can be assessed using the SV subspace, allowing
sign changes and rotations among equivalent vectors.

## References

- [R Installation and Administration: linear algebra](https://stat.ethz.ch/R-manual/R-patched/doc/manual/R-admin.html#Linear-algebra)
- [OpenBLAS build and runtime diagnostics](https://www.openmathlib.org/OpenBLAS/docs/faq/)
- [OpenBLAS runtime variables](https://www.openmathlib.org/OpenBLAS/docs/runtime_variables/)
- [RhpcBLASctl reference](https://cran.r-project.org/web/packages/RhpcBLASctl/refman/RhpcBLASctl.html)
