
source("options2idx.R")
source("DE_GSEA.R")

library(fgsea)
library(dplyr)  
library(tibble)
library(limma)
library(rhdf5)

cp_data = "cp_trimmed_predicted_RNAseq.gctx"
kd_data = "shRNA_predicted_RNAseq_profiles.gctx"
oe_data = "oe_predicted_RNAseq_profiles.gctx"
ko_data = "xpr_predicted_RNAseq_profiles.gctx"
ct_data = "ctl_predicted_RNAseq_profiles.gctx"

for(i in 1:length(filtered_conditions)) {
  print(paste0("Currently on index number ", i, " out of ", length(filtered_conditions), "."))
  condition = filtered_conditions[i]
  variables = strsplit(condition, "_")[[1]]
  chem = variables[1]
  target = variables[2]
  pert = variables[3]
  cell = variables[4]
  dose = variables[5]
  
  data_idx = return_data_idx(chem, target, pert, cell, dose)
  if(length(data_idx$cp_treated) < 2) {
    next
  } else if (length(data_idx$cp_control) < 2) {
    next
  } else if (length(data_idx$bp_treated) < 2) {
    next
  } else if (length(data_idx$bp_control) < 2) {
    next
  }
  results = GSEA_analysis_full(data_idx, pert)
  
  saveRDS(results, file=paste0("./results/", condition, ".rds"))
}

# load("remaining_conditions.RData")

# conditions = setdiff(conditions, filtered_conditions)
# save(conditions, file="remaining_conditions.RData")

# chem_names = unique(sub("_.*", "", conditions))

# chems = sample(chem_names, 20)
# chems = chems[order(chems)]
# 
# pattern = paste0("^(", paste(chems, collapse = "|"), ")")
# filtered_conditions = conditions[grepl(pattern, conditions)]
# save(filtered_conditions, file="current_conditions.RData")



