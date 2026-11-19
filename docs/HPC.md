# HPC, OnDemand and SLURM

See [BLAS libraries and threading](BLAS.md) for R library identification,
module discovery, external OpenBLAS inspection, and process-local startup settings.

OnDemand provides the session; the allocated CPUs and memory constrain the job.
Install the package in an R library available to both the interactive session and
compute nodes. No cluster path, R module, BLAS shared library or node exclusion is
hard-coded. SVA defaults to the installed package on compute nodes. A legacy
standalone script can still be specified with sva_lib.

DSS defaults to PSOCK processes, with fork available explicitly on Unix. PSOCK
may copy count data into each worker and use substantially more memory. Start with
n_workers=1, then a small worker count on representative data. Use only allocated
cores; the default uses parallelly::availableCores when installed and otherwise
falls back conservatively to SLURM_CPUS_PER_TASK or one core. Fork requires an
execution environment that supports forked R processes.

For SVA, use inst/examples/slurm_sva.R. Set r_module to your own module name only
when your cluster needs module loading. Rscript must resolve in the batch shell.
openblas_so is optional and has no default: use the cluster-supported R/BLAS
configuration rather than assuming a shared library from another machine applies.
More requested CPUs do not guarantee more BLAS threads. SVA workers does not
create a worker pool. Use workers=1 and request BLAS threads explicitly. For
example, cpus_per_task=8 with blas_threads=4 launches one R process requesting
four BLAS threads; a thread request above the CPU allocation is rejected.
Interactive blas_threads settings affect the R session, so set them deliberately
before other threaded work.

DSS can generate SLURM arrays with backend="slurm". Set slurm_partition,
slurm_module_load, rscript_path, slurm_mem and slurm_time for your cluster.
The DSS shell templates assume simple paths and identifiers without spaces or
shell metacharacters. Generated workers use the selected R executable and cluster
settings.

The inst/examples/record_environment.R script writes the R session, Bioconductor
version, and installed package versions.

## Optional installation library

Use R's default library location:

```sh
Rscript inst/examples/install_dependencies.R
```

Or select a writable library directory. New package installations go there;
downloaded package archives are retained in its downloads subdirectory:

```sh
Rscript inst/examples/install_dependencies.R --lib /your/shared/Rlibs
export R_LIBS_USER="/your/shared/Rlibs"
R CMD INSTALL --library="$R_LIBS_USER" .
```

The installer creates the directories if necessary. It prepends the custom
library for its own R process; export R_LIBS_USER for subsequent terminal commands
and batch jobs. Set it in the job environment as well as interactive terminals.
Compatible dependencies already available through other R library paths may be
reused; this option does not create a fully isolated dependency environment.
With no --lib argument, installation and archive locations use R's normal defaults.
Temporary compilation files continue to use R's temporary directory.

In a fresh RStudio session, add the custom directory before loading the package:

```r
.libPaths(c("/your/shared/Rlibs", .libPaths()))
library(methylomeTools)
```

Optional packages installed separately must also be given lib and, to retain
archives, destdir. For example:

```r
lib <- "/your/shared/Rlibs"
.libPaths(c(lib, .libPaths()))
install.packages("patchwork", lib=lib, destdir=file.path(lib, "downloads"))
```
