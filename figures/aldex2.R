########################################################
#
#       PACKAGES
#
########################################################

suppressWarnings(suppressPackageStartupMessages({
	library("tidyverse")
	library("ALDEx2")
	library("optparse")
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
metadata <- metadata$group

signal = read_tsv(signal_path)
noise = read_tsv(noise_path)
matrix = full_join(signal, noise, by = "S")
#as the joins adds NA for unknown values, it is needed to replcae NAs with 0 counts
matrix[is.na(matrix)] <- 0

S <- matrix$S
matrix$S <- NULL

########################################################
#
#       TEST
#
########################################################

# for testing purposes :D

# metadata = read_tsv("/mnt/projects_tn03/KMER_analysis/DATA/deseq2/rep0_kraken_bracken/FC2_20M/metadata.csv")
# metadata <- metadata$group
# signal = read_tsv("/mnt/projects_tn03/KMER_analysis/DATA/deseq2/rep0_kraken_bracken/FC2_20M/output.S")
# noise = read_tsv("/mnt/projects_tn03/KMER_analysis/DATA/deseq2/rep0_kraken_bracken/noise_20M//output.S")
# matrix = full_join(signal, noise, by = "S")
# matrix[is.na(matrix)] <- 0
# 
# S <- matrix$S
# matrix$S <- NULL

x.aldex <- aldex(matrix, metadata, mc.samples=128, test="kw", effect=TRUE, 
                 include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE, gamma= 1)

pdf_name = paste(output_path, "/aldex2_hist_pval.pdf", sep = "")
pdf(pdf_name)
my_DEG_pval_check$DEG_pval_check(x.aldex$kw.ep) 
dev.off()


species <- S[which(x.aldex$kw.eBH < 0.05)]
table_name = paste(output_path, "/differentially_abundant_species.tsv", sep = "")
write.table(species, table_name)
