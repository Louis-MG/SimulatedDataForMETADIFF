########################################################
#
#       PACKAGES
#
########################################################

suppressWarnings(suppressPackageStartupMessages({
	library("mixOmics")
	library("dplyr")
	library("DESeq2")
	library("optparse")
	library("tidyverse")
	library("readr")
	library("tidyr")
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
#       DESEQ2
#
########################################################

metadata = read_tsv("/mnt/projects_tn03/KMER_analysis/DATA/diff_abundances/rep0_kraken_bracken/FC2_20M/metadata.csv")
signal = read_tsv("/mnt/projects_tn03/KMER_analysis/DATA/diff_abundances/rep0_kraken_bracken/FC2_20M/output.S")
noise = read_tsv("/mnt/projects_tn03/KMER_analysis/DATA/diff_abundances/rep0_kraken_bracken/noise_20M//output.S")


#reads metadata and the two matrices, joins the matrices of species counts and keeps all species 
metadata = read_tsv(metadata_path)
 
signal = read_tsv(signal_path)
noise = read_tsv(noise_path)
#matrix = merge(signal, noise, by = "row.names", all.x = TRUE)
matrix = full_join(signal, noise, by = "S")
matrix[is.na(matrix)] <- 0

#print(head(matrix))

tab = matrix %>% column_to_rownames("S") %>%
      dplyr::select(metadata$name)

#hardcoded, change if your groups change 
state1 = "FC"
state2 = "noise"
 
meta3 = metadata %>% dplyr::filter(name %in% colnames(tab)) %>%dplyr::filter(group %in% c(state1, state2))
colData = meta3 %>% column_to_rownames("name")
 
data = tab %>% dplyr::select(all_of(c(rownames(colData))))
#dim(data)
proportion_zeros <- rowSums(data == 0) / ncol(data)
dataf <- subset(data, proportion_zeros <= 0.9)
 
#normalisation
#abun_out.clr = logratio.transfo(as.matrix(t(dataf)), logratio = 'CLR', offset = 1)
#abun_out.clr_final = abun_out.clr + min(abun_out.clr) %>% abs()
#class(abun_out.clr_final) <- "matrix"
#clr_table = as.matrix(abun_out.clr_final) %>% as.data.frame() %>% t() %>% as.data.frame()
#table_rounded = round(clr_table *1000, digits = 0)
 
 
#dataf2 = dplyr::select(as.data.frame(table_rounded), all_of(rownames(colData)))

#countData = table_rounded
dds = DESeqDataSetFromMatrix(countData = dataf,
                             colData = colData,
                             design= ~group ) #+ Animal.ID + Depth.location + Sampling.time the model matrix is not full rank, so the model cannot be fit as specified.
                                                                #One or more variables or interaction terms in the design formula are linear
                                                                #combinations of the others and must be removed.

#saves the mean-var and dispertion estimation plots that must be checked to ensure propoer fit of the model
pdf_name = paste(output_path, "/mean_var_plot.pdf", sep="")
pdf(pdf_name)
plot(log(rowMeans(counts(dds))),log(apply(counts(dds),1,var)), xlab="Mean species counts", ylab="Var species counts")
abline(0,1, col="red", lwd=2, lty=2)
dev.off()


 
#deseq
dds=DESeq(dds, fitType = "parametric")

#dispertion estimation plot
pdf_name = paste(output_path, "/DispEst_plot.pdf", sep = "")
pdf(pdf_name)
par(mfrow=c(2,2))
 
estimateDispersionsFit(dds, fitType = "parametric") %>%
  plotDispEsts(., main= "parametric", finalcol=NA, legend=F)
 
#Choix d'une régression non-paramétrique locale
estimateDispersionsFit(dds, fitType = "local") %>%
  plotDispEsts(., main = "local", finalcol=NA, legend=F)
 
#Choix d'une régression Gamma Poisson
estimateDispersionsFit(dds, fitType = "glmGamPoi") %>%
  plotDispEsts(., main = "glmGamPoi", finalcol=NA, legend=F,)
 
#Choix d'une régression moyenne
estimateDispersionsFit(dds, fitType = "mean") %>%
  plotDispEsts(., main = "mean", finalcol=NA, legend=F)
 
par(mfrow=c(1,1))
dev.off()


#resultsNames(dds)
rld <- rlogTransformation(dds)
 
final_vec_genus = vector()
final_df_genus = data_frame()
#resultsNames(dds)
nomf = str_replace_all(paste0("group_" , state2, "_vs_", state1 ), " ", "_")
nomf = str_replace_all(nomf, "/", "_")
 
res = results(dds, name=nomf, pAdjustMethod = "bonferroni")
res <- res[order(res$padj), ]

#checks the histogram of adjusted pvalues


resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "genomes"

#histogram of adj pvalues
pdf_name = paste(output_path, "/hist_pvalues.pdf", sep = "")
pdf(pdf_name)
my_DEG_pval_check$DEG_pval_check(res$pvalue)
dev.off()

#filters species using adjusted p-value
resdata2 = dplyr::filter(resdata, padj < 0.05)
#final_vec_genus = c(final_vec_genus, resdata2$genomes)
#temp = as.data.frame(counts(dds, normalized=TRUE))
#final_df_genus = rbind(final_df_genus, temp[resdata2$genomes,])
 

#resdata2[is.numeric(resdata2)] <- lapply(resdata2[is.numeric(resdata2)], function(x) format(x, scientific = TRUE))

#print(head(resdata))
#print(head(final_df_genus))
#writes output file
output_file = paste(output_path,"/differentially_abundant_species.tsv", sep = "")
write.table(resdata2$genomes, file = output_file, sep = "\t")
