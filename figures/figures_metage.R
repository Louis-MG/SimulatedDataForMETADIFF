library(ggplot2)
library(cowplot)
library(patchwork)
library(compositions)

df = read.csv("~/Documents/microbiome/figures/benjamini/table_results_metrics.csv", sep = ",", header = TRUE)
df$FC <- as.factor(df$FC)
df$lectures <- as.factor(df$lectures)

###############################################
#
#     GENERALE
#
###############################################

df_counts = read.table("/home/lmgueguen/Documents/microbiome/figures/benjamini/table_kmers_nbrs.tsv", header = TRUE)
df_counts$FC <- as.factor(df_counts$FC)
df_counts$reads <- as.factor(df_counts$reads)

ggplot(df_counts, aes(x = FC, y = counts, color = reads)) + geom_boxplot(fill="white", position = "dodge", outlier.color=NA) + scale_y_continuous(trans="log10") + labs(title = "Number of k-mers detected in cases and controls", x = "Fold-change", y = "Number of k-mers", color = "Sequencing depths") + theme(text=element_text(size=20)) + geom_point(position=position_jitterdodge(dodge.width=0.75)) + facet_grid(cols = vars(conditions))

ggsave(filename = "/home/lmgueguen/Documents/microbiome/figures/benjamini/kmer_counts.png", dpi = 320, width = 10, height = 6)

###############################################
#
#     NMBRE KMER PAR ESPECE et RECALL UNIQUES
#
###############################################

#par especes :

df_1 <- read.table(file = "/home/lmgueguen/Documents/microbiome/figures/benjamini/table_kmers_per_species.csv", header = TRUE, sep =",")
df_1[df_1$kmers == 0,]$kmers <- 1
df_1$FC <- as.factor(df_1$FC)
df_1$Sequencing.depths <- as.factor(df_1$Sequencing.depths)

df_2 = read.csv("~/Documents/microbiome/figures/benjamini/table_recall_uniq.csv", sep = ",", header = TRUE)
df_2$FC <- as.factor(df_2$FC)
df_2$Reads <- as.factor(df_2$Reads)

a <- ggplot(df_1, aes(x = FC, y = kmers, color = Sequencing.depths)) + geom_boxplot(fill="white", position = "dodge", outlier.color=NA) + scale_y_continuous(trans="log10") + labs(title = "Number of enriched k-mers per species", x = "Fold-change", y = "k-mers") + theme(text=element_text(size=10)) + geom_point(position=position_jitterdodge(dodge.width=0.75)) + facet_grid(rows = vars(species)) + guides(color=guide_legend(title="Sequencing depth")) + theme_bw()+ theme(strip.text = element_text(face = "italic"))
b <- ggplot(df_2, aes(x=FC, y = recall, color=Reads)) + geom_boxplot(fill="white", position = "dodge", outlier.color=NA) + labs(y = "Recall", title = "Recall of unique enriched k-mers per species") + scale_fill_brewer(palette = "Set2") + xlab(NULL) + geom_point(position=position_jitterdodge(dodge.width=0.75)) + theme(text=element_text(size=10)) + facet_grid(rows = vars(espece))+ theme_linedraw()+ theme(strip.text = element_text(face = "italic")) + theme_bw()

legend <- get_legend(
  # create some space to the left of the legend
  a + theme(legend.box.margin = margin(0, 0, 0, 10))
)

prow <- plot_grid(a +  theme(legend.position="none"), b +  theme(legend.position="none"),
                  labels = c("A", "B"),
                  hjust = -1,
                  nrow = 1,
                  ncol =2
)

plot_grid(prow, legend, rel_widths = c(3, .7))

ggsave('/home/lmgueguen/Documents/microbiome/figures/benjamini/nbr_kmer_recall_per_sp.png',dpi=400, height = 6, width = 10)
ggsave('/home/lmgueguen/Documents/microbiome/figures/benjamini/nbr_kmer_recall_per_sp.svg', dpi=320, height = 6, width = 10)

###############################################
#
#     PRECISIONS GENERALE
#
###############################################
# A CHANGER
df_recall = read.csv("~/Documents/microbiome/figures/benjamini/table_recall_uniq.csv", sep = ",", header = TRUE)
df_recall$FC <- as.factor(df_recall$FC)
df_recall$Reads <- as.factor(df_recall$Reads)

ggplot(data = df_recall, aes(x=FC, y = recall, color=Reads)) + geom_boxplot(fill="white", position = "dodge", outlier.color=NA) + labs(y = "Recall", title = "Recall of unique enriched k-mers for C.acnes") + scale_fill_brewer(palette = "Set2") + xlab(NULL) + geom_point(position=position_jitterdodge(dodge.width=0.75)) + theme(text=element_text(size=15)) + facet_grid(rows = vars(espece))+ theme_bw()+ theme(strip.text = element_text(face = "italic"))


df_prec = read.csv("~/Documents/microbiome/figures/benjamini/table_precision.csv", sep = ",", header = TRUE)
df_prec$FC <- as.factor(df_prec$FC)
df_prec$Reads <- as.factor(df_prec$Reads)

c <- ggplot(df_prec, aes(x = FC, y = precision, fill=Reads, color = Reads)) + geom_boxplot(fill="white", position = "dodge", outlier.color=NA)  + labs(title = "Precisions of enriched k-mers", x = "", y = "Precision") + theme(text=element_text(size=10)) + geom_point(position=position_jitterdodge(dodge.width=0.75)) + theme_bw()

legend <- get_legend(
  # create some space to the left of the legend
  a + theme(legend.box.margin = margin(0, 0, 0, 10))
)

prow <- plot_grid(a +  theme(legend.position="none"), c +  theme(legend.position="none"),
                  labels = c("A", "B"),
                  hjust = -1,
                  nrow = 1,
                  ncol =2
)

plot_grid(prow, legend, rel_widths = c(3, .7))

ggsave('/home/lmgueguen/Documents/microbiome/figures/nbre_kmer_and_precision_per_sp.png',dpi=400, height = 6, width = 10)
###############################################
#
#     PRECISIONS PAR ESPECE
#
###############################################


df_prec_sp = read.table(file = "~/Documents/microbiome/figures/benjamini/table_precision_per_sp.csv", sep = ",", header = T)
df_prec_sp$species <- as.factor(df_prec_sp$species)
df_prec_sp$Sequencing.depths <- as.factor(df_prec_sp$Sequencing.depths)
df_prec_sp$FC <- as.factor(df_prec_sp$FC)
df_prec_sp$replicate <- as.factor(df_prec_sp$replicate)

d <- ggplot(df_prec_sp, aes(x = FC, y = precision, color = Sequencing.depths)) + geom_boxplot(fill="white", position = "dodge", outlier.color=NA)  + labs(title = "Precisions of unique enriched k-mers per species", x = NULL, y = "Precision") + theme(text=element_text(size=15), ) + geom_point(position=position_jitterdodge(dodge.width=0.75)) + facet_grid(rows = vars(species)) + guides(color=guide_legend(title="Sequencing depth")) + theme_bw()+ theme(strip.text = element_text(face = "italic"))


###############################################
#
#     PATCHWORK PREC. RECALL. (UNIQ)
#
###############################################


(a + b) / (c + d)

legend <- get_legend(
  # create some space to the left of the legend
  a + theme(legend.box.margin = margin(0, 0, 0, 10))
)


prow <- plot_grid(c + theme(legend.position="none"), b + theme(legend.position="none"), a + theme(legend.position="none"), d + theme(legend.position="none"),
                  labels = c("A", "C", "B", "D"),
                  hjust = -1,
                  nrow = 2
)

plot_grid(prow, legend, rel_widths = c(3, .5))

ggsave('/home/lmgueguen/Documents/microbiome/figures/benjamini/plot_4_pannels.svg',dpi=320, height = 12, width = 12)


prow <- plot_grid(b + theme(legend.position="none"), d + theme(legend.position="none"),
                  labels = c("A","B"),
                  hjust = -1,
                  nrow = 1
)

plot_grid(prow, legend, rel_widths = c(3, .5))

ggsave('/home/lmgueguen/Documents/microbiome/figures/plot_2_pannels.png',dpi=320, height = 6, width = 12)

###############################################
#
#     DEMONSTRATION IMPACT PARAMETRES
#
###############################################

df <- read.table(file = "/home/lmgueguen/Documents/microbiome/figures/benjamini/table_kmers_per_species.csv", header = TRUE, sep =",")
df[df$kmers == 0,]$kmers <- 1
df$FC <- as.factor(df$FC)
df$Sequencing.depths <- as.factor(df$Sequencing.depths)

a <- ggplot(data = df[df$species == "C.acnes" & df$FC == 2,], aes(x=Sequencing.depths, y = kmers)) + geom_boxplot(fill="white", position = "dodge", outlier.color=NA) + labs(y = "k-mers", title = expression(paste("Number of k-mers for", italic("C. acnes"), " at FC2", sep = " "))) + scale_fill_brewer(palette = "Set2") + xlab("Sequencing depth") + theme(text=element_text(size=15)) 
b <- ggplot(data = df[df$species == "C.acnes" & df$Sequencing.depths == '20M',], aes(x=FC, y = kmers)) + geom_boxplot(fill="white", position = "dodge", outlier.color=NA) + labs(y = "k-mers", title = expression(paste("Number of k-mers for", italic("C. acnes"), " at 20M reads", sep = " "))) + scale_fill_brewer(palette = "Set2") + xlab("Fold-change") + theme(text=element_text(size=15)) 
c <- ggplot(data = df[df$FC == 2 & df$Sequencing.depths == '20M',], aes(x=species, y = kmers)) + geom_boxplot(fill="white", position = "dodge", outlier.color=NA) + labs(y = "k-mers", title = "Number of k-mers for species at 20M reads and FC2") + scale_fill_brewer(palette = "Set2") + xlab("Species") + theme(text=element_text(size=15), axis.text.x = element_text(face = "italic")) 

prow <- plot_grid(a, b, c,
                  labels = c("A", "B", "C"),
                  hjust = -1,
                  nrow = 3
)

plot_grid(prow, rel_widths = c(3, .75))

ggsave('/home/lmgueguen/Documents/microbiome/figures/benjamini/param_impact_overview.svg',dpi=320, height = 9, width = 9)

############################
#
#   PCA INCONNU ACNE
#
############################

# 
# df = read.table(file = "/home/lmgueguen/Documents/microbiome/glmnet/acne/matrix.tsv", header = T, sep = "\t")
# df$name <- NULL
# 
# Controls <- rep("Healthy", 40)
# Cases <- rep("Acne", 42)
# 
# df <- data.frame(t(df))
# df <- cbind(df, data.frame(c(Controls, Cases)))
# names(df)[names(df) == "c.Controls..Cases."] = "condition"
# 
# pca <- prcomp(df[,-32], scale = FALSE)
# summary(pca)
# 
# three_first_dim = data.frame(pca$x[,c(1,2,3)], df$condition)
# three_first_dim$df.condition <- as.factor(three_first_dim$df.condition)
# a <- ggplot(three_first_dim, aes(x = PC1, y = PC2, color = df$condition),) + geom_point() + theme(text=element_text(size=20)) + labs(x = "", y = "")
# b <- ggplot(three_first_dim, aes(x = PC3, y = PC2, color = df$condition),) + geom_point() + theme(text=element_text(size=20)) + labs(x = "PC3 = 8.7%", y = "PC2 = 14.6%")
# c <- ggplot(three_first_dim, aes(x = PC1, y = PC3, color = df$condition),) + geom_point() + guides(color=guide_legend(title="Condition")) + theme(text=element_text(size=20)) + labs(x = "PC1 = 23.5%", y = "PC3 = 8.7%")
# 
# legend <- get_legend(
#   # create some space to the left of the legend
#   c + theme(legend.box.margin = margin(0, 0, 0, 12))
# )
# 
# prow <- plot_grid(b + theme(legend.position="none"), a + theme(legend.position="none"), legend, c + theme(legend.position="none"),
#                   labels = c("A", "B", "", "C"),
#                   hjust = -1,
#                   nrow = 2
# )
# 
# plot_grid(prow, rel_widths = c(3, .4))
# 
# ggsave('/home/lmgueguen/Documents/microbiome/figures/benjamini/glmnet_pca_acne.svg',dpi=320, height = 12, width = 12)

################################
#
#   PCA INCONNU CRC 1/2 VS CTRL
#
################################

# 
# df = read.table(file = "/home/lmgueguen/Documents/microbiome/glmnet/stage_1_2_vs_ctrl/matrix.tsv", header = T, sep = "\t")
# df$name <- NULL
# 
# Controls <- rep("Healthy", 251)
# Cases <- rep("Stage I/II", 111)
# 
# df <- data.frame(t(df))
# 
# df <- cbind(df, data.frame(c(Controls, Cases)))
# names(df)[names(df) == "c.Controls..Cases."] = "condition"
# 
# pca <- prcomp(df[,-183], scale = FALSE)
# summary(pca)
# 
# three_first_dim = data.frame(pca$x[,c(1,2,3)], df$condition)
# three_first_dim$df.condition <- as.factor(three_first_dim$df.condition)
# a <- ggplot(three_first_dim, aes(x = PC1, y = PC2, color = df$condition),) + geom_point() + theme(text=element_text(size=20)) + labs(x = "", y = "")
# b <- ggplot(three_first_dim, aes(x = PC3, y = PC2, color = df$condition),) + geom_point() + theme(text=element_text(size=20)) + labs(x = "PC3 = 6.9%", y = "PC2 = 8.2%")
# c <- ggplot(three_first_dim, aes(x = PC1, y = PC3, color = df$condition),) + geom_point() + guides(color=guide_legend(title="Condition")) + theme(text=element_text(size=20)) + labs(x = "PC1 = 9.7%", y = "PC3 = 6.9%")
# 
# legend <- get_legend(
#   # create some space to the left of the legend
#   c + theme(legend.box.margin = margin(0, 0, 0, 12))
# )
# 
# prow <- plot_grid(b + theme(legend.position="none"), a + theme(legend.position="none"), legend, c + theme(legend.position="none"),
#                   labels = c("A", "B", "", "C"),
#                   hjust = -1,
#                   nrow = 2
# )
# 
# plot_grid(prow, rel_widths = c(3, .4))
# 
# ggsave('/home/lmgueguen/Documents/microbiome/figures/benjamini/glmnet_pca_stage_1_2_vs_ctrl.svg',dpi=320, height = 12, width = 12)

################################
#
#   PCA INCONNU CRC 1/2 VS 3/4
#
################################

# 
# df = read.table(file = "/home/lmgueguen/Documents/microbiome/glmnet/stage_1_2_vs_stage_3_4/matrix.tsv", header = T, sep = "\t")
# df$name <- NULL
# 
# df <- data.frame(t(df))
# 
# Controls <- rep("Stage I/II", 111)
# Cases <- rep("Stage III/IV", 74)
# 
# df <- cbind(df, data.frame(c(Controls, Cases)))
# names(df)[names(df) == "c.Controls..Cases."] = "condition"
# 
# pca <- prcomp(df[,-170], scale = FALSE)
# summary(pca)
# 
# three_first_dim = data.frame(pca$x[,c(1,2,3)], df$condition)
# three_first_dim$df.condition <- as.factor(three_first_dim$df.condition)
# a <- ggplot(three_first_dim, aes(x = PC1, y = PC2, color = df$condition),) + geom_point() + theme(text=element_text(size=20)) + labs(x = "", y = "")
# b <- ggplot(three_first_dim, aes(x = PC3, y = PC2, color = df$condition),) + geom_point() + theme(text=element_text(size=20)) + labs(x = "PC3 = 6.1%", y = "PC2 = 6.2%")
# c <- ggplot(three_first_dim, aes(x = PC1, y = PC3, color = df$condition),) + geom_point() + guides(color=guide_legend(title="Condition")) + theme(text=element_text(size=20)) + labs(x = "PC1 = 8.1%", y = "PC3 = 6.1%")
# 
# legend <- get_legend(
#   # create some space to the left of the legend
#   c + theme(legend.box.margin = margin(0, 0, 0, 12))
# )
# 
# prow <- plot_grid(b + theme(legend.position="none"), a + theme(legend.position="none"), legend, c + theme(legend.position="none"),
#                   labels = c("A", "B", "", "C"),
#                   hjust = -1,
#                   nrow = 2
# )
# 
# plot_grid(prow, rel_widths = c(3, .4))
# 
# ggsave('/home/lmgueguen/Documents/microbiome/figures/benjamini/glmnet_pca_stage_3_4_vs_1_2.svg',dpi=320, height = 12, width = 12)

###############################################
#
#     SUPP FIGURE SIGNAL DANS CONTROLES
#
###############################################

df = read.table(file = "/home/lmgueguen/Documents/microbiome/figures/benjamini/supp_table_reverse.csv", header = T, sep = ",")
df$replicate <- as.factor(df$replicate)
df$Signal <- as.factor(df$Signal)
df$Condition <- as.factor(df$Condition)

ggplot(data = df, aes(x=Condition, y = kmers, color = Condition)) + geom_boxplot(fill="white", position = "dodge", outlier.color=NA) + labs(y = "Kmers", title = "Number of enriched k-mers for enrichement in case or control") + scale_fill_brewer(palette = "Set2") + theme(text=element_text(size=15)) + facet_grid(cols =  vars(Signal)) + scale_y_continuous(trans="log10") + theme(legend.position="none")

ggsave('/home/lmgueguen/Documents/microbiome/figures/benjamini/supp_reverse_enrichment_comparison.svg',dpi=320, height = 12, width = 12)

# ##############################################
# 
#     SUPP FIGURE SIGNAL DANS CONTROLES ET CAS
# 
# ##############################################

df = read.table(file = "/home/lmgueguen/Documents/microbiome/figures/benjamini/table_depl_enrich_kmers.csv", header = TRUE, sep = ',')
df$condition <- as.factor(df$condition)
df$replicate <- as.factor(df$replicate)

ggplot(data = df, aes(x=condition, y = kmers, color = condition)) + geom_boxplot(fill = "white", position = "dodge", outlier.color=NA) + labs(y = "Kmers", title = "Number of enriched kmers for enrichement in case and control") + scale_fill_brewer(palette = "Set2") + theme(text=element_text(size=15)) + scale_y_continuous(trans="log10") + theme(legend.position="none")

ggsave('/home/lmgueguen/Documents/microbiome/figures/benjamini/supp_depl_enrich_kmers.svg',dpi=320, height = 12, width = 12)

# ##############################################
# 
#     SUPP FIGURE SIGNAL DANS BRUITS
# 
# ##############################################

df = read.table(file = "/home/lmgueguen/Documents/microbiome/figures/benjamini/table_kmers_noise_comparison.tsv", header = T, sep = "\t")
df$comparison <- as.factor(df$comparison)
df$condition <- as.factor(df$condition)

ggplot(data = df, aes(x=condition, y = kmers, color = condition)) + geom_boxplot(fill = "white", position = "dodge", outlier.color=NA) + geom_point(position=position_jitterdodge(dodge.width=0.75)) + labs(y = "Kmers", title = "Number of enriched kmers between controls samples") + scale_fill_brewer(palette = "Set2") + theme(text=element_text(size=15)) + theme(legend.position="none")

ggsave('/home/lmgueguen/Documents/microbiome/figures/benjamini/supp_comparison_noise_samples.svg',dpi=320, height = 12, width = 12)

###############################################
#
#     VERIFICATION APPARTENANCE DES KMERS
#
###############################################

# 
# test <- read.table("~/Documents/microbiome/output_FC15_20M/result_genomes_bifrost_query.tsv", sep = "\t", header = TRUE)
# sum(rowSums(test[,2:155]) == 0)
# test <- read.table("~/Documents/microbiome/output_FC15_40M/result_genomes_bifrost_query.tsv", sep = "\t", header = TRUE)
# sum(rowSums(test[,2:155]) == 0)
# test <- read.table("~/Documents/microbiome/output_FC2_10M/result_genomes_bifrost_query.tsv", sep = "\t", header = TRUE)
# sum(rowSums(test[,2:155]) == 0)
# test <- read.table("~/Documents/microbiome/output_FC2_20M/result_genomes_bifrost_query.tsv", sep = "\t", header = TRUE)
# sum(rowSums(test[,2:155]) == 0)
# test <- read.table("~/Documents/microbiome/output_FC2_40M/result_genomes_bifrost_query.tsv", sep = "\t", header = TRUE)
# sum(rowSums(test[,2:155]) == 0)
# test <- read.table("~/Documents/microbiome/output_FC3_10M/result_genomes_bifrost_query.tsv", sep = "\t", header = TRUE)
# sum(rowSums(test[,2:155]) == 0)
# test <- read.table("~/Documents/microbiome/output_FC3_40M/result_genomes_bifrost_query.tsv", sep = "\t", header = TRUE)
# sum(rowSums(test[,2:155]) == 0)


###############################################
#
#     UNITIGS
#
###############################################
# 
# df <- read.table("/home/lmgueguen/Documents/microbiome/figures/benjamini/table_dist_unitigs.tsv", sep = "\t", header = TRUE)
# df$FC <- as.factor(df$FC)
# df$lectures <- as.factor(df$lectures)
# df$replicate <- as.factor(df$replicate)
# to_remove <- which(df$FC %in% c("1.2", "1.3", "1.4"))
# df <- df[-to_remove,]
# 
# ggplot(df, aes( x = unitig_length, fill = replicate)) + geom_histogram(binwidth=.3, alpha=.85, position = "dodge") + labs(y = "frequency", title = "Distribution of unitig sizes") + facet_grid(lectures~FC) + scale_y_continuous(trans="log10", oob = scales::squish_infinite) + scale_x_continuous(trans="log10",  oob = scales::squish_infinite) + scale_fill_manual(values = c("#FF99FF", "#33CCFF", "#FF6633", "#66CC66", "#FF9900"))

###############################################
#
#     PERCENTAGE OF UNIQUE KMERS FOUND IN FC3 40M
#
###############################################

library(tidyverse)
df <- read.table(file = "/home/lmgueguen/Documents/microbiome/figures/benjamini/table_uniq_kmers_counts_FC3_40M.tsv", sep = "\t", header = TRUE)
uniq_kmers_number <- c(1643682, 5761020, 1497335, 1276446)
percent <- df
percent[,2:6] <- df[,2:6] / uniq_kmers_number
percent <- percent %>% pivot_longer(!espece)
# > mean(unlist(percent[percent$espece == "C.acnes",3]))
# [1] 0.9975928
# > sd(unlist(percent[percent$espece == "C.acnes",3]))
# [1] 1.365562e-05
# > mean(unlist(percent[percent$espece == "M.restricta",3]))
# [1] 0.9973798
# > sd(unlist(percent[percent$espece == "M.restricta",3]))
# [1] 0.0003012527
# > mean(unlist(percent[percent$espece == "S.epidermidis",3]))
# [1] 0.9941994
# > sd(unlist(percent[percent$espece == "S.epidermidis",3]))
# [1] 0.0002301559
# > mean(unlist(percent[percent$espece == "S.pyogenes",3]))
# [1] 0.6221832
# > sd(unlist(percent[percent$espece == "S.pyogenes",3]))
# [1] 0.0129019
percent$espece <-as.factor(percent$espece)
percent$name <- as.factor(percent$name)
ggplot(percent, aes(y = value, x = espece, color = espece)) + geom_boxplot() + theme(text=element_text(size=20)) + labs(title = "Percentage of unique kmers found at FC3 40M for each species", x = "Species", y = "percent")


###############################################
#
#     ranks of the 4 species in simulation align results
#
###############################################
# 
# df <- read.table(file = "/home/lmgueguen/Documents/microbiome/figures/benjamini/ranks.tsv", header = TRUE, sep = "\t")
# df$reads <- as.factor(df$reads)
# #df <- df[which(df$FC %in% c(1.5,2,3)),]
# df$FC <- as.factor(df$FC)
# df$species <- ordered(df$species, levels = c("Malassezia restricta", "Cutibacterium acnes", "Staphylococcus epidermidis", "Streptococcus pyogenes"))
# my_colors <- c("lightgreen", "coral", "lightblue4", "purple")
# names(my_colors) <- levels(df$species)
# 
# a <- ggplot(df[df$FC == 1.5,], aes( y = rank, x = reads, color = species)) + geom_boxplot(position = "dodge", fill = "white", outlier.color=NA) + labs(x = "", y = "", title = "Distribution of species ranks in alignments at FC 1.5") + scale_color_manual(name = "Species", values = my_colors) + theme(text=element_text(size=15)) + geom_point(position=position_jitterdodge(dodge.width=.75)) + coord_trans(y='log10') 
# b <- ggplot(df[df$FC == 2,], aes( y = rank, x = reads, color = species)) + geom_boxplot(position = "dodge", fill = "white", outlier.color=NA) + labs(x = "", y = "", title = "Distribution of species ranks in alignments at FC 2") + scale_color_manual(name = "Species", values = my_colors) + theme(text=element_text(size=15)) + geom_point(position=position_jitterdodge(dodge.width=.75)) + coord_trans(y='log10')
# c <- ggplot(df[df$FC == 3,], aes( y = rank, x = reads, color = species)) + geom_boxplot(position = "dodge", fill = "white", outlier.color=NA) + labs(x = "Sequencing depth", y = "", title = "Distribution of species ranks in alignments at FC 3") + scale_color_manual(name = "Species", values = my_colors) + theme(text=element_text(size=15)) + geom_point(position=position_jitterdodge(dodge.width=.75)) + coord_trans(y='log10')
# 
# legend <- get_legend(
#   # create some space to the left of the legend
#   a + theme(legend.box.margin = margin(0, 0, 0, 12))
# )
# 
# prow <- plot_grid(a + theme(legend.position="none"), b +  theme(legend.position="none"), c + theme(legend.position="none"),
#                   labels = c("A", "B", "C"),
#                   hjust = -1,
#                   nrow = 3
# )
# 
# plot_grid(prow, legend, rel_widths = c(3, .6))
# 
# d <- ggplot(df[df$FC == 12,], aes( y = rank, x = reads, color = species)) + geom_boxplot(position = "dodge", fill = "white", outlier.color=NA) + labs(x = "", title = "Distribution of species ranks in alignments at FC 1.2") + scale_fill_brewer(palette = "Set2", ) + theme(text=element_text(size=15)) + geom_point(position=position_jitterdodge(dodge.width=.75)) + coord_trans(y='log10') 
# e <- ggplot(df[df$FC == 13,], aes( y = rank, x = reads, color = species)) + geom_boxplot(position = "dodge", fill = "white", outlier.color=NA) + labs(x = "", title = "Distribution of species ranks in alignments at FC 1.3") + scale_fill_brewer(palette = "Set2") + theme(text=element_text(size=15)) + geom_point(position=position_jitterdodge(dodge.width=.75)) + coord_trans(y='log10')
# f <- ggplot(df[df$FC == 14,], aes( y = rank, x = reads, color = species)) + geom_boxplot(position = "dodge", fill = "white", outlier.color=NA) + labs(x = "", title = "Distribution of species ranks in alignments at FC 1.4") + scale_fill_brewer(palette = "Set2") + theme(text=element_text(size=15)) + geom_point(position=position_jitterdodge(dodge.width=.75)) + coord_trans(y='log10')
# 
# legend <- get_legend(
#   # create some space to the left of the legend
#   f + theme(legend.box.margin = margin(0, 0, 0, 12))
# )
# 
# prow <- plot_grid(d + theme(legend.position="none"), e +  theme(legend.position="none"), f + theme(legend.position="none"),
#                   labels = c("D", "E", "F"),
#                   hjust = -1,
#                   nrow = 3
# )
# 
# plot_grid(prow, legend, rel_widths = c(3, .6))

#################################
#
#   PRECISSION, RECALL ALT.METH
#
#################################

library("ggplot2")

df = read.table(file = "/mnt/projects_tn03/KMER_analysis/DATA/diff_abundances/all_metrics_table.tsv", header = TRUE, sep = "\t")
df$method <- as.factor(df$method)
df$FC <- as.factor(df$FC)
df$sequencing_depth <- as.factor(df$sequencing_depth)

ggplot(df, mapping = aes(x = method, y = Precision, color = FC)) +
  geom_boxplot() +
  facet_grid(vars(sequencing_depth)) +
  ggtitle("Precision of different methods to detect differential abundances.") +
  theme(text=element_text(size=16)) +
  xlab("Methods") +
  theme_bw()

ggsave(filename = "/mnt/projects_tn03/KMER_analysis/DATA/figures/metadiff_versus_others_precision.pdf", dpi = 300)

ggplot(df, mapping = aes(x = method, y = Recall, color = FC)) +
  geom_boxplot() +
  facet_grid(vars(sequencing_depth)) +
  ggtitle("Recall of different methods to detect differential abundances.") +
  theme(text=element_text(size=16)) +
  xlab("Methods") +
  theme_bw()

ggsave(filename = "/mnt/projects_tn03/KMER_analysis/DATA/figures/metadiff_versus_others_recall.pdf", dpi = 300)

# ################################
# 
#   ALT. DIFF. ABUNDANCES METHODES
# 
# ################################

library("ggplot2")
library("cowplot")

df = read.table(file = "/home/lmgueguen/Documents/microbiome/figures/diff_abundances/all_metrics_table.tsv", header = TRUE, sep = "\t")
# df = read.table(file = "/mnt/projects_tn03/KMER_analysis/DATA/diff_abundances/all_metrics_table.tsv", header = TRUE, sep = "\t")
df$method <- as.factor(df$method)
df$FC <- as.factor(df$FC)
df$sequencing_depth <- as.factor(df$sequencing_depth)

ggplot(df, mapping = aes(x = method, y = Precision, color = FC)) +
  geom_boxplot() +
  facet_grid(vars(sequencing_depth)) +
  ggtitle("Precision of META-DIFF and alternative methods.") +
  theme(text=element_text(size=16)) +
  xlab("Methods") +
  ylab("Precision of enriched species") +
  theme_bw()

ggsave(filename = "/home/lmgueguen/Documents/microbiome/figures/metadiff_versus_others_precision.pdf", dpi = 300, width = 10, height = 6)

ggplot(df, mapping = aes(x = method, y = Recall, color = FC)) +
  geom_boxplot() +
  facet_grid(vars(sequencing_depth)) +
  ggtitle("Recall of META-DIFF and alternative methods.") +
  theme(text=element_text(size=16)) +
  xlab("Methods") +
  ylab("Recall of enriched species") +
  theme_bw()

ggsave(filename = "/home/lmgueguen/Documents/microbiome/figures/metadiff_versus_others_recall.pdf", dpi = 300, height = 6, width = 10)

a <- ggplot(df, mapping = aes(x = method, y = Precision, color = FC)) +
  geom_boxplot() +
  facet_grid(vars(sequencing_depth)) +
  ggtitle("Precision of META-DIFF and alternative methods.") +
  theme(text=element_text(size=16)) +
  xlab("Methods") +
  ylab("Precision of enriched species") +
  theme_bw()

b <- ggplot(df, mapping = aes(x = method, y = Recall, color = FC)) +
  geom_boxplot() +
  facet_grid(vars(sequencing_depth)) +
  ggtitle("Recall of META-DIFF and alternative methods.") +
  theme(text=element_text(size=16)) +
  xlab("Methods") +
  ylab("Recall of enriched species") +
  theme_bw()

legend <- get_legend(
  # create some space to the left of the legend
  a + theme(legend.box.margin = margin(0, 0, 0, 12))
)

prow <- plot_grid(a + theme(legend.position="none"), b +  theme(legend.position="none"),
                  labels = c("E", "F"),
                  hjust = -1,
                  nrow = 1,
                  ncol = 2
)
plot_grid(prow, legend, rel_widths = c(3, .5))
ggsave(filename = "/home/lmgueguen/Documents/microbiome/figures/metadiff_versus_others_precision_recall_and_precision.png", dpi = 300, height = 6, width = 12)


##################################
# 
#   METRICS UNITIGS (FUNCTIONS)
# 
##################################

library("ggplot2")
library("cowplot")


df = read.table("/home/lmgueguen/Documents/microbiome/figures/precision_unitigs.tsv", header = TRUE, sep = "\t")
df$FC <- as.factor(df$FC)
df$sequencing.depth <- as.factor(df$sequencing.depth)
df$replicate <- as.factor(df$replicate)


c <- ggplot(df, aes(x = FC, y = precision, color = sequencing.depth)) +
  geom_boxplot(outlier.color=NA) +
  labs(title = "Precision of unitigs for enriched species", x = "Fold-change", y = "Precision", color = "Sequencing depth") +
  geom_point(position=position_jitterdodge(dodge.width=.75)) +
  theme_bw() + 
  coord_cartesian(ylim = c(0, 1)) +
  theme(text=element_text(size=10))

prow <- plot_grid(c,
                  labels = c("C"),
                  hjust = -1,
                  nrow = 1
)

plot_grid(prow)

ggsave(filename = "/home/lmgueguen/Documents/microbiome/figures/precision_unitigs.png", dpi = 300, height = 6, width = 10)

# tenter de faire les memes metriques mais regardant ce qui est unique a une espece enrichie seulement

df = read.table("/home/lmgueguen/Documents/microbiome/figures/table_recall_functional.tsv", header = TRUE, sep = "\t")
df$Fold.change <- as.factor(df$Fold.change)
df$Sequencing.depth <- as.factor(df$Sequencing.depth)
df$recall <- df$genes.recovered / df$ground_truth

d <- ggplot(df, aes(x = Fold.change, y = recall, color = Sequencing.depth)) +
      geom_boxplot(outlier.color=NA) +
      labs(title = "Recall of enriched functions per species", x = "Fold-change", y = "Recall", color = "Sequencing depth") +
      geom_point(position=position_jitterdodge(dodge.width=.75)) +
      theme_bw() + 
      coord_cartesian(ylim = c(0, 1)) +
      theme(text=element_text(size=10)) +
      facet_grid(rows = vars(Species)) +
      theme(strip.text = element_text(face = "italic"))

d
ggsave(filename = "/home/lmgueguen/Documents/microbiome/figures/recall_functional_unitigs.png", dpi = 300, height = 6, width = 10)

legend <- get_legend(
  # create some space to the left of the legend
  d + theme(legend.box.margin = margin(0, 0, 0, 12))
)

prow <- plot_grid(c + theme(legend.position = "none"), d +  theme(legend.position="none"),
                  labels = c("C","D"),
                  hjust = -1,
                  nrow = 1,
                  ncol = 2
)

plot_grid(prow, legend, rel_widths = c(3, .5))

ggsave(filename = "/home/lmgueguen/Documents/microbiome/figures/2_panels_unitigs.png", dpi = 300, height = 6, width = 12)

#####################################
#
#   big figure all metrics
#     
#####################################


library("ggplot2")
library("cowplot")

# precision unique

df_prec_sp = read.table(file = "~/Documents/microbiome/figures/benjamini/table_precision_per_sp.csv", sep = ",", header = T)
df_prec_sp$species <- as.factor(df_prec_sp$species)
df_prec_sp$Sequencing.depths <- as.factor(df_prec_sp$Sequencing.depths)
df_prec_sp$FC <- as.factor(df_prec_sp$FC)
df_prec_sp$replicate <- as.factor(df_prec_sp$replicate)


a <- ggplot(df_prec_sp, aes(x = FC, y = precision, color = Sequencing.depths)) + 
  geom_boxplot(fill="white", position = "dodge", outlier.color=NA) + 
  labs(title = "Precision of unique enriched k-mers per species", x = NULL, y = "Precision") + 
  theme(text=element_text(size=15), ) + 
  geom_point(position=position_jitterdodge(dodge.width=0.75)) + 
  facet_grid(rows = vars(species)) + 
  guides(color=guide_legend(title="Sequencing depth")) + 
  theme_bw() + 
  theme(strip.text = element_text(face = "italic"))

# recall unique
df_2 = read.csv("~/Documents/microbiome/figures/benjamini/table_recall_uniq.csv", sep = ",", header = TRUE)
df_2$FC <- as.factor(df_2$FC)
df_2$Reads <- as.factor(df_2$Reads)


b <- ggplot(df_2, aes(x=FC, y = recall, color=Reads)) + 
  geom_boxplot(fill="white", position = "dodge", outlier.color=NA) + 
  labs(y = "Recall", title = "Recall of unique enriched k-mers per species") + 
  scale_fill_brewer(palette = "Set2") + 
  xlab(NULL) + 
  geom_point(position=position_jitterdodge(dodge.width=0.75)) + 
  theme(text=element_text(size=15)) + 
  facet_grid(rows = vars(espece)) + 
  theme_linedraw() + 
  theme(strip.text = element_text(face = "italic")) + theme_bw()

prow1 <- plot_grid(a +  theme(legend.position="none"), b +  theme(legend.position="none"),
                  labels = c("A", "B"),
                  hjust = -1,
                  nrow = 1,
                  ncol =2
)

legend1 <- get_legend(
  # create some space to the left of the legend
  a + theme(legend.box.margin = margin(0, 0, 0, 10))
)


# precision unitigs
df = read.table("/home/lmgueguen/Documents/microbiome/figures/precision_unitigs.tsv", header = TRUE, sep = "\t")
df$FC <- as.factor(df$FC)
df$sequencing.depth <- as.factor(df$sequencing.depth)
df$replicate <- as.factor(df$replicate)


c <- ggplot(df, aes(x = FC, y = precision, color = sequencing.depth)) +
  geom_boxplot(outlier.color=NA) +
  labs(title = "Precision of unitigs for enriched species", x = "Fold-change", y = "Precision", color = "Sequencing depth") +
  geom_point(position=position_jitterdodge(dodge.width=.75)) +
  theme_bw() + 
  coord_cartesian(ylim = c(0, 1)) +
  theme(text=element_text(size=11))


df = read.table("/home/lmgueguen/Documents/microbiome/figures/table_recall_functional.tsv", header = TRUE, sep = "\t")
df$Fold.change <- as.factor(df$Fold.change)
df$Sequencing.depth <- as.factor(df$Sequencing.depth)
df$recall <- df$genes.recovered / df$ground_truth

# recall unitigs
d <- ggplot(df, aes(x = Fold.change, y = recall, color = Sequencing.depth)) +
  geom_boxplot(outlier.color=NA) +
  labs(title = "Recall of enriched functions per species", x = "Fold-change", y = "Recall", color = "Sequencing depth") +
  geom_point(position=position_jitterdodge(dodge.width=.75)) +
  theme_bw() + 
  coord_cartesian(ylim = c(0, 1)) +
  theme(text=element_text(size=11)) +
  facet_grid(rows = vars(Species)) +
  theme(strip.text = element_text(face = "italic"))

prow2 <- plot_grid(c +  theme(legend.position="none"), d +  theme(legend.position="none"),
                  labels = c("C", "D"),
                  hjust = -1,
                  nrow = 1,
                  ncol =2
)

legend2 <- get_legend(
  # create some space to the left of the legend
  d + theme(legend.box.margin = margin(0, 0, 0, 10))
)

# precision META-DIFF versus others

df = read.table(file = "/home/lmgueguen/Documents/microbiome/figures/diff_abundances/all_metrics_table.tsv", header = TRUE, sep = "\t")
# df = read.table(file = "/mnt/projects_tn03/KMER_analysis/DATA/diff_abundances/all_metrics_table.tsv", header = TRUE, sep = "\t")
df$method <- as.factor(df$method)
df$FC <- as.factor(df$FC)
df$sequencing_depth <- as.factor(df$sequencing_depth)


e <- ggplot(df, mapping = aes(x = method, y = Precision, color = FC)) +
  geom_boxplot() +
  facet_grid(vars(sequencing_depth)) +
  ggtitle("Precision of META-DIFF and alternative methods.") +
  theme(text=element_text(size=15)) +
  xlab("Methods") +
  ylab("Precision of enriched species") +
  theme_bw()

# recall META-DIFF versus others

f <- ggplot(df, mapping = aes(x = method, y = Recall, color = FC)) +
  geom_boxplot() +
  facet_grid(vars(sequencing_depth)) +
  ggtitle("Recall of META-DIFF and alternative methods.") +
  theme(text=element_text(size=15)) +
  xlab("Methods") +
  ylab("Recall of enriched species") +
  theme_bw()

prow3 <- plot_grid(e +  theme(legend.position="none"), f +  theme(legend.position="none"),
                   labels = c("E", "F"),
                   hjust = -1,
                   nrow = 1,
                   ncol =2
)

legend3 <- get_legend(
  # create some space to the left of the legend
  e + theme(legend.box.margin = margin(0, 0, 0, 10))
)

plot_grid(prow1, legend1, prow2, legend2, prow3, legend3,
          nrow = 3,
          rel_widths = c(1, 0.2, 1, 0.2, 1, 0.2),
          labels = c("", "", "", "", "", ""),
          label_size = 20)

ggsave(filename = "/home/lmgueguen/Documents/microbiome/figures/all_benchmark_metrics.png", dpi = 300, height = 12, width = 12)
