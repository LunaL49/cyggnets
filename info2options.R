
source("compound2info.R")
#load("chem_names.Rdata")
library(rhdf5)

# cp_data = "cp_trimmed_predicted_RNAseq.gctx"
# kd_data = "shRNA_predicted_RNAseq_profiles.gctx"
# oe_data = "oe_predicted_RNAseq_profiles.gctx"
# ko_data = "xpr_predicted_RNAseq_profiles.gctx"
# ct_data = "ctl_predicted_RNAseq_profiles.gctx"

# choose chemical of interest
# chem_name = chem_names[sample(1:length(chem_names), 1)] # choose a random chemical
# chem_info = compound_info(chem_name)
# target_gene = chem_info$target_genes[sample(1:length(chem_info$target_genes), 1)] # choose a random target gene

pert_options = function(target_gene) {
  options = c()
  # shRNA knockdown
  kd_genes = h5read(kd_data, "0/META/COL/pertname")
  if(target_gene %in% kd_genes){
    options = c(options, "shRNA knockdown")
  }
  # CRISPR knockout
  ko_genes = h5read(ko_data, "0/META/COL/pertname")
  if(target_gene %in% ko_genes){
    options = c(options, "CRISPR knockout")
  }
  # overexpression 
  oe_genes = h5read(oe_data, "0/META/COL/pertname")
  if(target_gene %in% oe_genes){
    options = c(options, "overexpression")
  }
  
  # note: it is possible that the target gene is not in any of the three biological 
  # perturbation datasets, the user will have to select a different drug/target gene 
  # in that case.
  return(options)
}

cell_line_options = function(chem_name, target_gene, pert_method) {
  
  # cell lines available in the cp dataset
  cp_pertname = h5read(cp_data, "/0/META/COL/pertname")
  cp_timepoint = h5read(cp_data, "0/META/COL/timepoint")
  
  cp_chem_idx = which(cp_pertname == chem_name)
  cp_timepoint_idx = which(cp_timepoint == "24 h") # 24h timepoint data only
  cp_idx = intersect(cp_chem_idx, cp_timepoint_idx)
  cp_cell_lines = h5read(cp_data, "0/META/COL/cell")[cp_idx]
  cp_cell_lines = unique(cp_cell_lines)
  
  # which biological perturbation dataset to use
  if(pert_method == "overexpression") {
    data = oe_data
  } else if (pert_method == "shRNA knockdown") {
    data = kd_data
  } else if (pert_method == "CRISPR knockout") {
    data = ko_data
  }
  
  # cell lines available in the biological perturbation dataset
  pert_genes = h5read(data, "0/META/COL/pertname")
  pert_gene_idx = which(pert_genes == target_gene)
  
  pert_timepoint = h5read(data, "0/META/COL/timepoint")
  pert_timepoint_idx = which(pert_timepoint == "96 h") # 96h timepoint data only
  pert_idx = intersect(pert_gene_idx, pert_timepoint_idx)
  
  pert_cell_lines = h5read(data, "0/META/COL/cell")[pert_idx]
  pert_cell_lines = unique(pert_cell_lines)
  
  cell_lines = intersect(cp_cell_lines, pert_cell_lines)
  
  return(cell_lines)
}

dose_options = function(chem_name, cell_line) {
  
  cp_pertname = h5read(cp_data, "/0/META/COL/pertname")
  cp_timepoint = h5read(cp_data, "0/META/COL/timepoint")
  
  cp_chem_idx = which(cp_pertname == chem_name)
  cp_timepoint_idx = which(cp_timepoint == "24 h") # 24h timepoint data only
  cp_idx = intersect(cp_chem_idx, cp_timepoint_idx)
  
  cp_cell_lines = h5read(cp_data, "0/META/COL/cell")[cp_idx]
  idx = cp_idx[which(cp_cell_lines == cell_line)]
  doses = h5read(cp_data, "/0/META/COL/dose")[idx]
  doses = unique(doses)
  
  return(doses)
}


