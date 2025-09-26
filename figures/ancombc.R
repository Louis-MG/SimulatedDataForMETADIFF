########################################################
#
#       PACKAGES
#
########################################################

suppressWarnings(suppressPackageStartupMessages({
  library("dplyr")
  library("ANCOMBC")
  library("optparse")
  library("tidyverse")
  library("phyloseq")
  source("/mnt/projects_tn03/KMER_analysis/gitlab/DEG_pval_check.R")
}))
  
########################################################
#
#       ARGS
#
########################################################


option_list = list(
  make_option(c("-m", "--metadata"), type="character", default=NULL,
              help="File with metadata, must contain at least name and group coloumns.", metavar="character"),
  make_option(c("-n", "--noise"), type="character", default=NULL,
              help="Species counts for each noise sample.", metavar="character"),
  make_option(c("-s", "--signal"), type="character", default=NULL,
              help = "Species counts for each fold-change sample.", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output directory. ", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$noise)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (count file for noise).n", call.=FALSE)
}
if (is.null(opt$signal)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (count file for signal).n", call.=FALSE)
}
if (is.null(opt$output)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (output directory file).n", call.=FALSE)
}
if (is.null(opt$metadata)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (metadata).n", call.=FALSE)
}

output_path = opt$output
noise_path = opt$noise
signal_path =opt$signal
metadata_path = opt$metadata

if (!dir.exists(output_path)) {
  # If it does not exist, create the folder
  dir.create(output_path, recursive = TRUE)  # Use recursive = TRUE to create any necessary parent directories
  cat("Folder created:", output_path, "\n")
} else {
  cat("Folder already exists:", output_path, "\n")
}

########################################################
#
#       DATA
#
########################################################

#reads metadata and the two matrices, joins the matrices of species counts and keeps all species 
metadata = read_tsv(metadata_path)

signal = read_tsv(signal_path)
noise = read_tsv(noise_path)
matrix = full_join(signal, noise, by = "S")
#as the joins adds NA for unknown values, it is needed to replcae NAs with 0 counts
matrix[is.na(matrix)] <- 0

########################################################
#
#       TEST
#
########################################################


colnames(metadata) <- c("sample", "group")
matrix <- matrix %>% 
  tibble::column_to_rownames("S")

metadata <- metadata %>% 
  tibble::column_to_rownames("sample") 

#ici data doit etre un objet phyloseq :///
OTU = otu_table(matrix, taxa_are_rows = TRUE)
samples = sample_data(metadata)

phylo_seq_obj <- phyloseq(OTU, samples)

#erreur parce que des taxons ont une variance nulle, mais specifier le group et 
#struct_zero devrait permettre de les retirer
output = ancombc2(data = phylo_seq_obj, 
                  fix_formula = "group", rand_formula = NULL,
                  p_adj_method = "holm", pseudo_sens = TRUE,
                  prv_cut = 0.60, lib_cut = 0, s0_perc = 0.05,
                  struc_zero = FALSE, 
                  alpha = 0.05, n_cl = 1, verbose = TRUE,
                  iter_control = list(tol = 1e-5, max_iter = 20, 
                                      verbose = FALSE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = NULL, 
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100))

pdf_name = paste(output_path, "/ancombc2_hist_pvalues.pdf", sep = "")
pdf(pdf_name)
my_DEG_pval_check$DEG_pval_check(output$res$p_groupnoise)
dev.off()

res <- output$res %>%
  filter(diff_groupnoise == TRUE,
         passed_ss_groupnoise == TRUE)

res$taxon

table_name = paste(output_path, "/differentially_abundant_species.tsv", sep = "")
write.table(res$taxon, table_name)
