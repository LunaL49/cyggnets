
library(rhdf5)
# chem_name = "vorinostat"
# target_gene = "HDAC3"
# pert_method = "shRNA knockdown"
# cell_line = "A375"
# drug_dose = "10 uM"

# cp_data = "cp_trimmed_predicted_RNAseq.gctx"
# kd_data = "shRNA_predicted_RNAseq_profiles.gctx"
# oe_data = "oe_predicted_RNAseq_profiles.gctx"
# ko_data = "xpr_predicted_RNAseq_profiles.gctx"
# ct_data = "ctl_predicted_RNAseq_profiles.gctx"

# pick out the relevant data columns for control, chemical perturbation, and biological perturbations 
# note: use 24h DMSO control for chemical perturbation
# note: use 96h EMPTY_VECTOR control for all biological perturbation modalities
return_data_idx = function(chem_name, target_gene, pert_method, cell_line, drug_dose) {
  # Treatment data for chemical perturbation
  cp_chem_names = h5read(cp_data, "0/META/COL/pertname")
  cp_chem_names_idx = which(cp_chem_names == chem_name)
  cp_cell_lines = h5read(cp_data, "0/META/COL/cell")
  cp_cell_lines_idx = which(cp_cell_lines == cell_line)
  cp_timepoints = h5read(cp_data, "0/META/COL/timepoint")
  cp_timepoints_idx = which(cp_timepoints == "24 h")
  cp_doses = h5read(cp_data, "0/META/COL/dose")
  cp_doses_idx = which(cp_doses == drug_dose)
  
  cp_treated_idx = Reduce(intersect, list(cp_chem_names_idx, cp_cell_lines_idx, 
                                          cp_timepoints_idx, cp_doses_idx))
  
  # Negative control for chemical perturbation
  # note: there are a LOT of negative control samples for cell lines PC3, A375, A549, MCF7, HT29, & HA1E
  ct_pertnames = h5read(ct_data, "0/META/COL/pertname")
  ct_pertnames_idx = which(ct_pertnames == "DMSO")
  ct_cell_lines = h5read(ct_data, "0/META/COL/cell")
  ct_cell_lines_idx = which(ct_cell_lines == cell_line)
  ct_timepoints = h5read(ct_data, "0/META/COL/timepoint")
  ct_timepoints_idx = which(ct_timepoints == "24 h")
  
  cp_control_idx = Reduce(intersect, list(ct_pertnames_idx, ct_cell_lines_idx, 
                                          ct_timepoints_idx))
  
  # Treatment data for chosen biological perturbation (bp)
  if(pert_method == "overexpression") {
    data = oe_data
  } else if (pert_method == "shRNA knockdown") {
    data = kd_data
  } else if (pert_method == "CRISPR knockout") {
    data = ko_data
  }
  bp_gene_names = h5read(data, "0/META/COL/pertname")
  bp_gene_names_idx = which(bp_gene_names == target_gene)
  bp_cell_lines = h5read(data, "0/META/COL/cell")
  bp_cell_lines_idx = which(bp_cell_lines == cell_line)
  bp_timepoints = h5read(data, "0/META/COL/timepoint")
  bp_timepoints_idx = which(bp_timepoints == "96 h")
  
  bp_treated_idx = Reduce(intersect, list(bp_gene_names_idx, bp_cell_lines_idx, 
                                          bp_timepoints_idx))
  
  # Negative control for biological perturbation (empty vector - "ev" - control)
  ev_pertnames_idx = which(ct_pertnames == "EMPTY_VECTOR")
  ev_cell_lines_idx = ct_cell_lines_idx
  ev_timepoints_idx = which(ct_timepoints == "96 h")
  
  bp_control_idx = Reduce(intersect, list(ev_pertnames_idx, ev_cell_lines_idx, 
                                          ev_timepoints_idx))
  
  return(list(cp_treated = cp_treated_idx, cp_control = cp_control_idx,
              bp_treated = bp_treated_idx, bp_control = bp_control_idx))
}




