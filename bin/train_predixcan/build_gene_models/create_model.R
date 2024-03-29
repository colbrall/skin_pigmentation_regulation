argv <- commandArgs(trailingOnly = TRUE)
# source("GTEx_Tissue_Wide_CV_elasticNet.R")
# source("GTEx_Tissue_Wide_CV_elasticNet_debug.R")
source("Zhang_melanocyte_CV_elasticNet_debug.R")

study <- argv[1]
expression_RDS <- argv[2]
geno_file <- argv[3]
gene_annot_RDS <- argv[4]
snp_annot_RDS <- argv[5]
n_k_folds <- as.numeric(argv[6])
alpha <- as.numeric(argv[7])
out_dir <- argv[8]
chrom <- argv[9]
snpset <- argv[10]
window <- as.numeric(argv[11])

TW_CV_model(expression_RDS, geno_file, gene_annot_RDS, snp_annot_RDS,
    n_k_folds, alpha, out_dir, study, chrom, snpset, window)
