library(methylomeTools)
# beta_matrix and metadata are supplied by the user.
metadata <- metadata[match(colnames(beta_matrix), metadata$sample_id), , drop=FALSE]
stopifnot(identical(metadata$sample_id, colnames(beta_matrix)))
rownames(metadata) <- metadata$sample_id
mod <- model.matrix(~ exposure + age + sex + batch, metadata)
mod0 <- model.matrix(~ age + sex + batch, metadata)
stopifnot(nrow(mod) == ncol(beta_matrix), qr(mod)$rank == ncol(mod),
          qr(mod0)$rank == ncol(mod0))
job <- submit_methyl_svs_slurm(beta_matrix, mod, mod0, outdir="outputs/sva",
  cpus_per_task=8, time="12:00:00", partition=NULL, r_module=NULL,
  value_type="beta", feature_strategy="hvf_only", n_hvf=200000,
  override_hvf=TRUE, num_sv_method="be", sva_method="irw",
  submit=FALSE)  # Write batch files without submitting.
# Set r_module to your cluster's module name if needed. The package must be
# installed in an R library visible to the compute node's Rscript process.
# job$result_rds is read only after the job has completed successfully.
