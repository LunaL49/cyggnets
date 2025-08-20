
#BiocManager::install("limma")
#BiocManager::install("fgsea")
#install.packages("msigdbr")
library(fgsea)
library(dplyr)  
library(tibble)
library(limma)
library(rhdf5)
# chem_name = "vorinostat"
# target_gene = "HDAC3"
# pert_method = "shRNA knockdown"
# cell_line = "A375"
# drug_dose = "10 uM"
# load("one_set_idx.Rdata")

load("valid_genes_mask.Rdata")
load("hallmark_gene_sets.Rdata")
load("C2_CP_gene_sets.Rdata")

# cp_data = "cp_trimmed_predicted_RNAseq.gctx"
# kd_data = "shRNA_predicted_RNAseq_profiles.gctx"
# oe_data = "oe_predicted_RNAseq_profiles.gctx"
# ko_data = "xpr_predicted_RNAseq_profiles.gctx"
# ct_data = "ctl_predicted_RNAseq_profiles.gctx"

# for the data matrix, across each row is the same gene, down each column is one repeat
GSEA_analysis_full = function(data_idx, pert_method) {
  data_path = "/0/DATA/0/matrix"
  n_rows = 23614
  
  cp_treated_data = h5read(cp_data, data_path, 
                           index = list(1:n_rows, data_idx$cp_treated))
  cp_control_data = h5read(ct_data, data_path,
                           index = list(1:n_rows, data_idx$cp_control))
  if(pert_method == "overexpression") {
    data = oe_data
  } else if (pert_method == "shRNA knockdown") {
    data = kd_data
  } else if (pert_method == "CRISPR knockout") {
    data = ko_data
  }
  bp_treated_data = h5read(data, data_path,
                           index = list(1:n_rows, data_idx$bp_treated))
  bp_control_data = h5read(ct_data, data_path,
                           index = list(1:n_rows, data_idx$bp_control))
  
  cp_results = GSEA_analysis(cp_control_data, cp_treated_data)
  bp_results = GSEA_analysis(bp_control_data, bp_treated_data)
  
  return(list(cp = cp_results, bp = bp_results))
}


GSEA_analysis = function(control_data, treated_data) {
  
  #format counts and metadata for limma
  counts = cbind(control_data, treated_data)
  sample_names = paste0("S", seq_len(ncol(counts)))
  colnames(counts) = sample_names
  gene_names = h5read(ct_data, "0/META/ROW/id")
  rownames(counts) = gene_names
  sample_info = data.frame(
    sample = colnames(counts),
    condition = c(rep("control", ncol(control_data)), rep("treated", ncol(treated_data)))
  )
  
  # filter to keep protein-coding genes using mask
  counts = counts[genes_to_keep,]
  
  # pass filtered counts matrix and metadata to limma
  design = model.matrix(~0 + condition, data = sample_info)
  colnames(design) = c("control", "treated")
  
  fit = limma::lmFit(counts, design)
  
  contrast_matrix = makeContrasts(
    Treated_vs_Control = treated - control, # This specifies the comparison
    levels = design
  )
  
  fit2 = contrasts.fit(fit, contrast_matrix)
  fit2 = eBayes(fit2)
  
  de_results = topTable(fit2, coef = "Treated_vs_Control", number = Inf,
                        adjust.method = "fdr")
  
  ranked_genes = de_results$t
  names(ranked_genes) = rownames(de_results)
  ranked_genes = sort(ranked_genes, decreasing = TRUE)
  
  # Hallmark GSEA
  gsea_H = fgsea(pathways = hallmark_gene_sets, stats = ranked_genes,
                       nPermSimple = 1000)
  gsea_H = gsea_H[order(gsea_H$padj), ]
  
  # C2 canonical pathways GSEA
  gsea_C2CP = fgsea(pathways = C2_CP_gene_sets, stats = ranked_genes,
                       nPermSimple = 1000)
  gsea_C2CP = gsea_C2CP[order(gsea_C2CP$padj), ][1:50, ]
  
  # return both GSEA results
  return(list(GSEA_H = gsea_H, GSEA_C2CP = gsea_C2CP))
}

# results = GSEA_analysis_full(one_set_idx, pert_method)
# save(results, file="GSEA_results.Rdata")





