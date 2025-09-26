########################################################
#
#       PACKAGES
#
########################################################

suppressWarnings(suppressPackageStartupMessages({
	library("optparse")
	library("tidyverse")
	library("mixOmics")
	source("/mnt/projects_tn03/KMER_analysis/gitlab/DEG_pval_check.R")
}))
# Use of the wilcoxon (unpaired,  so similar to mann-whitney U test) test is because :
# 1. All the observations from both groups are independent of each other,
# 2. The responses are at least ordinal (i.e., one can at least say, of any two observations, which is the greater),
# 3. Under the null hypothesis H0, the distributions of both populations are identical.[3]
# 4. The alternative hypothesis H1 is that the distributions are not identical


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

# Check if the folder exists
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
#removes species with too many 0's
proportion_zeros <- rowSums(matrix == 0) / ncol(matrix)
matrix <- subset(matrix, proportion_zeros <= 0.9)

#saves species info
S <- matrix$S
matrix$S <- NULL

########################################################
#
#       MEAN - VAR plot
#
########################################################

pdf_name = paste(output_path, "/mean_var_plot.pdf", sep="")
pdf(pdf_name)
plot(log(rowMeans(matrix)),log(apply(matrix,1,var)), xlab="Mean species counts", ylab="Var species counts")
abline(0,1, col="red", lwd=2, lty=2)
dev.off()


########################################################
#
#       NORMALISATION
#
########################################################

abun_out.clr = logratio.transfo(as.matrix(t(matrix)), logratio = 'CLR', offset = 1)
abun_out.clr_final = abun_out.clr + min(abun_out.clr) %>% abs()
class(abun_out.clr_final) <- "matrix"
clr_table = as.matrix(abun_out.clr_final) %>% as.data.frame() %>% t() %>% as.data.frame()
table_rounded = round(clr_table *1000, digits = 0)
table_rounded$S <- S


########################################################
#
#       TEST
#
########################################################


wilcoxon <- function(mat) {
  x <- as.numeric(unlist(mat[2:11]))
  y <- as.numeric(unlist(mat[12:21]))
  res <- wilcox.test(x, y)
  return(res$p.value)
}

#saves pvalues
table_rounded$pvalues <- apply(table_rounded, 1, wilcoxon)
#saves histogram of pvalues in a pdf
pdf_name = paste(output_path, "/wilcox_pvalues_histogram.pdf", sep = "")
pdf(pdf_name)
my_DEG_pval_check$DEG_pval_check(table_rounded$pvalues)
dev.off()

significant_species <- table_rounded[which(table_rounded$pvalues < 0.05/dim(matrix)[1]),]$S
write.table(significant_species, paste(output_path, "/differentially_abundant_species.tsv", sep = ""), sep = "\t")

