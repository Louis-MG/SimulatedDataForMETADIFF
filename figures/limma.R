########################################
#
#   PACKAGES
#
########################################

suppressWarnings(suppressPackageStartupMessages({
	library("ggplot2")
	library("ggrepel")
	library("tidyverse")
	library("limma")
	library("edgeR")
	library("optparse")
	source("/mnt/projects_tn03/KMER_analysis/gitlab/DEG_pval_check.R")
	source("/mnt/projects_tn03/KMER_analysis/gitlab/VolcanoPlot.R")
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

# Check if the folder exists
if (!dir.exists(output_path)) {
  # If it does not exist, create the folder
  dir.create(output_path, recursive = TRUE)  # Use recursive = TRUE to create any necessary parent directories
  cat("Folder created:", output_path, "\n")
} else {
  cat("Folder already exists:", output_path, "\n")
}

########################################
#
#   DATA
#
########################################

metadata = read_tsv(metadata_path)
signal = read_tsv(signal_path)
noise = read_tsv(noise_path)
matrix = full_join(signal, noise, by = "S")
matrix[is.na(matrix)] <- 0

S <- matrix$S
matrix$S <- NULL

#proportion_zeros <- rowSums(matrix == 0) / ncol(matrix)
#matrix <- subset(matrix, proportion_zeros <= 0.9)

#transformation en objet DGE
dge <- DGEList(matrix, samples=metadata)
dge <- calcNormFactors(dge) 
# voir ce qui est le defaut pour la normalisation

########################################
#
#   MEAN - VAR PLOT
#
########################################

#voir la tronche de lupus.meta
pdf_name = paste(output_path, "/mean_var_plot.pdf", sep="")
pdf(pdf_name)
plot(log(rowMeans(matrix)),log(apply(matrix,1,var)), xlab="Mean species counts", ylab="Var species counts")
abline(0,1, col="red", lwd=2, lty=2)
dev.off()


########################################
#
#   NORMALISATION
#
########################################

pdf_name = paste(output_path, "/voom_plot.pdf", sep = "")
dge.norm <- voom(dge, normalize.method="quantile", plot=T, save.plot=TRUE)
dev.off()

#voom.plot <- dge.norm %>% data.frame(x=.$voom.xy$x, y=.$voom.xy$y) %>%
#  ggplot(aes(x=x, y=y)) + geom_point() +
#  geom_line(data=data.frame(x=dge.norm$voom.line$x, y=dge.norm$voom.line$y), 
#            aes(x=x, y=y), color="red") +
#  labs(x="log2(count size + 0.5)", 
#       y="Sqrt(standard deviation)", 
#       title="Voom: Mean-variance trend") +
#  theme_bw()
#voom.plot
#dev.off()

########################################
#
#   TEST
#
########################################

design <- model.matrix(~0+dge$sample$group)
colnames(design) <- levels(dge$sample$group)


contrast <- paste0("FC","-","noise")
contrast.matrix <- makeContrasts(contrast, levels=design)

sub_dge.norm <- dge.norm[, dge.norm$targets$group %in% c("FC", "noise")]
sub_dge.norm$targets$group <- as.factor(as.vector(sub_dge.norm$targets$group))
sub_dge.norm$targets$group <- relevel(sub_dge.norm$targets$group, ref="noise")

fit <- lmFit(dge.norm, design)
fit.contrast <- contrasts.fit(fit, contrast.matrix)

res.limma <- eBayes(fit.contrast, trend=T, robust=T)
res.limma$p.value.adj <- p.adjust(res.limma$p.value, method="BH") #Benjamini Hochberg
res.limma <- as.data.frame(res.limma)
colnames(res.limma)[1] <- "log2FC"
res.limma <- res.limma[, c("log2FC", "p.value", "p.value.adj")]
rownames(res.limma) <- rownames(dge.norm)

########################################
#
#   VOLCANO PLOT bof
#
########################################

#volcano <- VolcanoPlot(res.limma, c("noise", "FC"), log2FC_thresholds=c(-0.5, 0.5), p.value="p.value.adj", 
#                       return_DE=F, text_size=6)
#volcano

#df_DE_limma <- VolcanoPlot(res.limma, vals, p.value="p.value.adj", plot=F)

#histogram of pvalues
pdf_name = paste(output_path, "/limma_hist_pvalues.pdf", sep = "")
pdf(pdf_name)
my_DEG_pval_check$DEG_pval_check(res.limma$p.value)
dev.off()

S[which(res.limma$p.value.adj < 0.05)]

table_name = paste(output_path, "/differentially_abundant_species.tsv", sep = "")
write.table(S[which(res.limma$p.value.adj < 0.05)], table_name)
