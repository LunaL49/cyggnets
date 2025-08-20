
load("chem_names_ordered.Rdata")
source("compound2info.R")
source("info2options.R")

cp_data = "cp_trimmed_predicted_RNAseq.gctx"
kd_data = "shRNA_predicted_RNAseq_profiles.gctx"
oe_data = "oe_predicted_RNAseq_profiles.gctx"
ko_data = "xpr_predicted_RNAseq_profiles.gctx"
ct_data = "ctl_predicted_RNAseq_profiles.gctx"

library(doParallel)
library(foreach)
cl = makeCluster(7)
registerDoParallel(cl)

conditions = foreach(chem = chem_names, .combine = c, .packages = c("httr", "jsonlite", "rhdf5")) %dopar% {
  out = list()
  info = compound_info(chem)
  targets = info$target_genes # target_genes should never be empty
  if(length(targets) > 3) {
    targets = sample(targets, 3)
  }
  for(target in targets){
    pert = pert_options(target)
    if(is.null(pert)){
      next
    } else {
      for(pert_option in pert){
        cell_lines = cell_line_options(chem, target, pert_option)
        if(length(cell_lines) == 0) {
          next
        } else {
          if(length(cell_lines) > 2) {
            cell_lines = sample(cell_lines, 2)
          }
          for(cell_line in cell_lines){
            doses = dose_options(chem, cell_line) # doses should never be empty 
            if("10 uM" %in% doses) {
              condition_string = paste(chem, target, pert_option, cell_line, "10 uM", sep = "_")
              out[[length(out) + 1]] = condition_string
            }
            if("1.11 uM" %in% doses) {
              condition_string = paste(chem, target, pert_option, cell_line, "1.11 uM", sep = "_")
              out[[length(out) + 1]] = condition_string
            }
          }
        }
      }
    }
  }
  return(out)
}

conditions = unlist(conditions)

save(conditions, file="all_condition_names.RData")

load("all_condition_names.RData")
chems = sample(chem_names, 50)
chems = chems[order(chems)]

pattern = paste0("^(", paste(chems, collapse = "|"), ")")
filtered_conditions = conditions[grepl(pattern, conditions)]
print(length(filtered_conditions))

save(filtered_conditions, file="conditions_short.RData")








