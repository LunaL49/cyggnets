
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("rhdf5", quietly = TRUE)) {
  BiocManager::install("rhdf5")
}
library(rhdf5) # needed for reading .gctx files

#### Trim chemical perturbation file by selecting certain cell lines/chemicals ####

# find the list of all relevant cell lines
cp_data = "cp_predicted_RNAseq_profiles.gctx"
kd_data = "shRNA_predicted_RNAseq_profiles.gctx"
oe_data = "oe_predicted_RNAseq_profiles.gctx"
ko_data = "xpr_predicted_RNAseq_profiles.gctx"
ct_data = "ctl_predicted_RNAseq_profiles.gctx"

cp_cell_line = h5read(cp_data, "/0/META/COL/cell")
ct_cell_line = h5read(ct_data, "/0/META/COL/cell")
kd_cell_line = h5read(kd_data, "/0/META/COL/cell")
oe_cell_line = h5read(oe_data, "/0/META/COL/cell")
ko_cell_line = h5read(ko_data, "/0/META/COL/cell")

pert_cell_lines = unique(c(kd_cell_line, oe_cell_line, ko_cell_line))
cp_ct_cell_lines = intersect(cp_cell_line, ct_cell_line)
relevant_cell_lines = intersect(pert_cell_lines, cp_ct_cell_lines)

cp_idx_from_cell_line = which(cp_cell_line %in% relevant_cell_lines) # indices to save from cp_data based on cell line info

# find the list of all relevant compounds
cp_chem_name = h5read(cp_data, "/0/META/COL/pertname")
all_chem_name = unique(cp_chem_name)
filter_chem_name = grep("[0-9]", all_chem_name, invert = TRUE, value = TRUE) # get rid of any investigational drugs with numbers in their names
relevant_chem = filter_chem_name[!grepl("-", filter_chem_name) & nchar(filter_chem_name) >= 5] # get rid of anything with a hyphen in the name or is too short

cp_idx_from_chem_name = which(cp_chem_name %in% relevant_chem)

cp_idx_to_save = intersect(cp_idx_from_cell_line, cp_idx_from_chem_name) # final indices to keep

# write to new .gctx file
input_file = "cp_predicted_RNAseq_profiles.gctx"
output_file = "cp_trimmed_predicted_RNAseq.gctx"
data_path = "/0/DATA/0/matrix"
n_rows = 23614
cols_to_keep = cp_idx_to_save

h5createFile(output_file)
h5createGroup(output_file, "0")
h5createGroup(output_file, "0/DATA")
h5createGroup(output_file, "0/DATA/0")

new_dims = c(n_rows, length(cols_to_keep))
chunk_size = c(n_rows, 10)

h5createDataset(output_file, data_path, 
                dims = new_dims, 
                chunk = chunk_size,
                storage.mode = "double") 

block_size = 10
n_output_cols = length(cols_to_keep)
col_starts = seq(1, n_output_cols, by = block_size)

# this part takes a long time (~4 hours)
for (start_idx in col_starts) {
  end_idx = min(start_idx + block_size - 1, n_output_cols)
  col_block = cols_to_keep[start_idx:end_idx]
  
  # Read block from input
  data_block = h5read(input_file, data_path,
                       index = list(1:n_rows, col_block))
  
  # Write to output file at the correct position
  output_cols = start_idx:end_idx
  h5write(data_block, output_file, data_path,
          index = list(1:n_rows, output_cols))
}

# copy over metadata
row_metadata = h5read(cp_data, "/0/META/ROW/id")
h5createGroup(output_file, "0/META")
h5createGroup(output_file, "0/META/ROW")
h5write(row_metadata, output_file, "/0/META/ROW/id")

h5createGroup(output_file, "0/META/COL")
col_meta_cell = h5read(cp_data, "/0/META/COL/cell")[cols_to_keep]
h5write(col_meta_cell, output_file, "/0/META/COL/cell")
col_meta_dose = h5read(cp_data, "/0/META/COL/dose")[cols_to_keep]
h5write(col_meta_dose, output_file, "/0/META/COL/dose")
col_meta_pertname = h5read(cp_data, "/0/META/COL/pertname")[cols_to_keep]
h5write(col_meta_pertname, output_file, "/0/META/COL/pertname")
col_meta_smiles = h5read(cp_data, "/0/META/COL/smiles")[cols_to_keep]
h5write(col_meta_smiles, output_file, "/0/META/COL/smiles")
col_meta_timepoint = h5read(cp_data, "/0/META/COL/timepoint")[cols_to_keep]
h5write(col_meta_timepoint, output_file, "/0/META/COL/timepoint")



