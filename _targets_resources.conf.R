# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Set Resources for Clustermq or Future -----
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


#SBATCH --output=<%= resources$log_path %>   #./log/%x_%A-%a.out

template_path <-  "/home/clz4002/batchtools.slurm.tmpl"

log_path <- here::here("log", "%x_%A-%a.out")


plan <- tweak(
  batchtools_slurm,
  template = template_path,
  resources =
    list(
      log_path = log_path,
      ncpus = 1,
      memory = "16G",
      walltime = 600,
      partition = "scu-cpu",
      ntasks = 1
    )
)


plan_alignment <- tweak(
  batchtools_slurm,
  template = template_path,
  resources =
    list(
      log_path = log_path,
      ncpus = 51,
      memory = "3G",
      walltime = "48:00:00",
      partition = "scu-cpu",
      ntasks = 1
    )
)



resources_all <- tar_resources(
  future = tar_resources_future(plan = plan)
)


resources_alignment <- tar_resources(
  future = tar_resources_future(plan = plan_alignment)
)
