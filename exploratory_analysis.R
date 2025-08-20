
library(rhdf5)
chem_name = "vorinostat"
target_gene = "HDAC3"
pert_method = "shRNA knockdown"
cell_line = "A375"
drug_dose = "10 uM"
load("one_set_idx.Rdata")

cp_data = "cp_trimmed_predicted_RNAseq.gctx"
kd_data = "shRNA_predicted_RNAseq_profiles.gctx"
oe_data = "oe_predicted_RNAseq_profiles.gctx"
ko_data = "xpr_predicted_RNAseq_profiles.gctx"
ct_data = "ctl_predicted_RNAseq_profiles.gctx"

data_path = "/0/DATA/0/matrix"
n_rows = 23614
# for the data matrix, across each row is the same gene, down each column is one repeat
cp_treated_data = h5read(cp_data, data_path, 
                         index = list(1:n_rows, one_set_idx$cp_treated))
cp_control_data = h5read(ct_data, data_path,
                         index = list(1:n_rows, one_set_idx$cp_control))
if(pert_method == "overexpression") {
  data = oe_data
} else if (pert_method == "shRNA knockdown") {
  data = kd_data
} else if (pert_method == "CRISPR knockout") {
  data = ko_data
}
bp_treated_data = h5read(data, data_path,
                         index = list(1:n_rows, one_set_idx$bp_treated))
bp_control_data = h5read(ct_data, data_path,
                         index = list(1:n_rows, one_set_idx$bp_control))

# check total gene counts per repeats and across conditions
cp_treated_gene_count_sum = colSums(cp_treated_data) 
hist(cp_treated_gene_count_sum, breaks=20) # closely clustered around 62000-63000
cp_control_gene_count_sum = colSums(cp_control_data)
hist(cp_control_gene_count_sum, breaks=20) # closely clustered around 63000
bp_treated_gene_count_sum = colSums(bp_treated_data) 
hist(bp_treated_gene_count_sum, breaks=20) # only 8 replicates, but still mostly around 62000-63000
bp_control_gene_count_sum = colSums(bp_control_data) 
hist(bp_control_gene_count_sum, breaks=20) # clustered around 62000
# gene count totals are sufficiently similar across conditions and repeats

# check number of highly expressed genes across repeats and conditions
cp_treated_highly_expressed = colSums(cp_treated_data > 5)
hist(cp_treated_highly_expressed, breaks=20) # clustered around 4800, between 4000-5500
cp_control_highly_expressed = colSums(cp_control_data > 5)
hist(cp_control_highly_expressed, breaks=20) # clustered around 5000, between 4500-6000
bp_treated_highly_expressed = colSums(bp_treated_data > 5)
hist(bp_treated_highly_expressed, breaks=20) # only 8 replicates, around 5200
bp_control_highly_expressed = colSums(bp_control_data > 5)
hist(bp_control_highly_expressed, breaks=20) # clustered around 5000, between 4500-5500
# cp_control seems to have slightly more highly expressed genes, but not hugely unusual

# find the gene with highest mean expression, look at its variance 
most_expressed = rowMeans(bp_control_data)
idx = which.max(most_expressed)
most_expressed[idx] # cp_treated: 12.77 cp_control: 13.74 bp_treated: 13.00 bp_control: 13.55
var(bp_control_data[idx,]) # cp_treated: 1.13 cp_control: 0.52 bp_treated: 1.56 bp_control: 0.87
most_expressed[6] # cp_treated: 2.34 cp_control: 1.99 bp_treated: 2.78 bp_control: 2.27
var(bp_control_data[6,]) # cp_treated: 0.49 cp_control: 0.15 bp_treated: 0.07 bp_contorl: 0.25
# generally true that the most expressed gene has higher variance than a random gene

# log-transformed counts histogram
x = sample(1:ncol(bp_control_data), 1, replace=T)
hist(log(bp_control_data[,x]), breaks=100) 
# fairly consistent across cp_treated, cp_control, bp_treated, and bp_control

# subsample the larger group between test and control, log2(x + 1) transform, then PCA
cp_n = min(c(ncol(cp_control_data), ncol(cp_treated_data)))
bp_n = min(c(ncol(bp_control_data), ncol(bp_treated_data)))

cp_control_subsample = cp_control_data[ , sample(1:ncol(cp_control_data), cp_n, replace=F)]
cp_treated_subsample = cp_treated_data[ , sample(1:ncol(cp_treated_data), cp_n, replace=F)]
bp_control_subsample = bp_control_data[ , sample(1:ncol(bp_control_data), bp_n, replace=F)]
bp_treated_subsample = bp_treated_data[ , sample(1:ncol(bp_treated_data), bp_n, replace=F)]

counts = cbind(bp_control_subsample, bp_treated_subsample)
sample_names = paste0("S", seq_len(2*bp_n))
colnames(counts) = sample_names
sample_info = data.frame(
  sample = colnames(counts),
  condition = c(rep("control", bp_n), rep("treated", bp_n))
)

keep_genes = rowMeans(counts) > 1
filtered_counts = counts[keep_genes, ]

log_counts = log2(filtered_counts + 1)
log_counts_t = t(log_counts) # transpose for PCA

pca = prcomp(log_counts_t, scale. = TRUE)

plot(pca$x[,1:2],
     col = as.factor(sample_info$condition),
     pch = 16,
     xlab = "PC1", ylab = "PC2")

legend("topright",
       legend = levels(as.factor(sample_info$condition)),
       col = 1:2, pch = 16)
# cp: not bad, controls definitely cluster with each other, treated is more spread
# out but still distinct from control. 121 replicates in each arm. 
# bp: slightly problematic, 5-6 of 8 repeats clustered with the controls. 










