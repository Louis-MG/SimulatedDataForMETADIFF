library(compositions)
library(tidyverse)
library(moments)
library(ggplot2)

#############################################
#
#   FUNCTIONS
#
#############################################


DE <- function(df, fold_change = c(), target = c()){ 
  #
  #' df: data frame. First column must be genomes, second abundances. Columns will be renamed anyway
  #' fold_change: vector. >1 for up, <1 for down. Give target_up or target_down accordingly
  #' target
  #' One of the two above must be provided. In the absence of one or another, the compensation will be distributed evenly across the other values
  #
  if (is.data.frame(df) == FALSE) {
    stop("ERROR: df must be a dataframe.")
  }
  if (is.vector(fold_change) == FALSE) {
    stop("ERROR: fold change must be a vector")
  }
  colnames(df) <- c("genomes", "abundances")
  if ( is.null(target) ) {
    stop("ERROR: provide at least one target to up or down regulate")
  } else { #if both will compensate only on down
    index = target
  }
  if (length(fold_change) != length(target)) {
    stop("ERROR: provide a fold change for each target")
  }
  #applies fold change
  DE_target_abundance <- df$abundances[df$genomes %in% index]*fold_change
  df$abundances[df$genomes %in% index] <- DE_target_abundance
  df$abundances <- df$abundances/sum(df$abundances)
  qplot(x = genomes, y = abundances, data = df, geom = "point")
  return(df)
}

DE_compensate <- function(df, fold_change = c(), target = c(), compensate = c()){ 
  #' df: data frame. First column must be genomes, second abundances. Columns will be renamed anyway
  #' fold_change: vector. >1 for up, <1 for down. Give target_up or target_down accordingly
  #' target: vetcor. Must be the same plength as the gold_change vector
  #' compensate: vector. Targets for abundance compensation.
  #
  if (is.data.frame(df) == FALSE) {
    stop("ERROR: df must be a dataframe.")
  }
  if (is.vector(fold_change) == FALSE) {
    stop("ERROR: fold change must be a vector")
  }
  colnames(df) <- c("genomes", "abundances")
  if ( is.null(target) || is.null(compensate) ) {
    stop("ERROR: provide at least one target to up or down regulate")
  }
  if (length(fold_change) != length(target)) {
    stop("ERROR: provide a fold change for each target")
  }
  #applies fold change
  DE_target_abundance <- df$abundances[df$genomes %in% target]*fold_change
  diff_abundances <- DE_target_abundance - df$abundances[df$genomes %in% target]
  df$abundances[df$genomes %in% target] <- DE_target_abundance
  df$abundances[df$genomes %in% compensate] <- df$abundances[df$genomes %in% compensate] - sum(diff_abundances)/length(compensate)
  print(sum(df$abundances))
  df$abundances <- df$abundances/sum(df$abundances)
  ggplot(data = df, mapping = aes(genomes, abundances))+geom_point()
  return(df)
}

Noise <- function(df, mu = 1, sigma = 0.05) {
  #'
  #' This function creates a biological replicate of the given dataframe. It uses
  #' a set of normaly distributed factors to change change the original values by
  #' multiplying them by the factor. 
  #'
  #' df: dataframe. First column must be genome names/indexes.
  #' mu: mean of the normal distribution.
  #' sigma: standard deviation of the normal distribution
  #'
  if (is.data.frame(df) == FALSE) {
    stop("ERROR: df must be a dataframe.")
  }
  colnames(df) <- c("genomes", "abundances")
  if ( is.numeric(df$abundances) == FALSE ) {
    stop("ERROR: second column of the dataframe must be an integer type.")
  }
  if (is.numeric(c(mu, sigma)) == FALSE) {
    stop("ERROR: mu and sigma must be numeric variables.")
  }
  number_values <- dim(df)[1]
  factor_variation <- rnorm(number_values, mu, sigma)
  df$abundances <- df$abundances * factor_variation
  df$abundances <- df$abundances/sum(df$abundances)
  return(df)
}

#############################################
#
#   DATA
#
#############################################


df = read.table("./SimulatedDataForMETADIFF/taxonomy/abondances/acne_S_2022.tsv", sep = "\t", header = TRUE)
# nettoyage
df <- df[-c(which(df$S == "Homo sapiens")),]
genomes <-  df$S
df$S <- NULL

#############################################
#
#     NORMALISATION CLR
#
#############################################

#cnormalisation
# CLR: Center log ratio transformation

df_mat <- as.matrix(df)
df_mat <- t(df_mat)
df.clr <- clr(df_mat)
apply(df.clr, 1, sum)

#adds absolute value of the minimum value to have all counts >= 0
add_min <- function(x) {
  minimum <- min(x)
  x <- x + abs(minimum)
  return(x)
}

df.clr <- apply(df.clr, 1, function(x) add_min(x))
df.clr <- as.data.frame(df.clr)
df.clr$genomes <- genomes

df_plot <- pivot_longer(df.clr, -c("genomes"))

ggplot(data = df_plot, mapping = aes(x = name, y = value)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 75, hjust =1))


####################
#
#   AUTRE
#
####################


#visualisation og 100 first abundances (sample means) to check 
df_mat <- as.matrix(df.clr[,1:82])
rownames(df_mat) <- genomes
abundances <- apply(df_mat, 1, FUN = mean)
df_plot <- as.data.frame(genomes)
df_plot$abundances <- abundances

# sort values
df_plot <- df_plot[rev(order(df_plot$abundances)), ]
# selects 150 first
df_plot <- df_plot[0:150,]

ggplot(data = df_plot, mapping = aes(x=reorder(genomes, abundances, function(x) sort(x, decreasing = FALSE)), abundances)) + geom_point() + labs(title = "CLR-normalised abundances of genomes", y = "abundances", x = "genomes") +
      theme(axis.text.x = element_text(angle = 75, hjust =1))

##################################
#
#     EXPORT 150 MOST ABUNDANTS
#
##################################

#df_plot$abundances <- df_plot$abundances/sum(df_plot$abundances)
#write.table(df_plot, file ="/home/lmgueguen/Documents/microbiome/acne_150_2022.tsv", sep = "\t", row.names = FALSE, col.names = FALSE)


##################################
#
#     NOISE AND TEST
#
##################################
# add Noise() from simulation_abundance.R

par(mfrow = c(1,2))
# les 2 dataframes
df_1 <- df_plot
df_1$abundances <- df_1$abundances/sum(df_1$abundances)
df_2 <- Noise(df_1, 1, 0.05)

head(df_1)
head(df_2)

#plot des rapport
test1 <- df_1$abundances/df_2$abundances
hist(test1) # pas plus que 1.3 de raport, pas de signal
qqnorm(test1)
qqline(test1)
skewness(test1)
kurtosis(test1)
shapiro.test(test1)

test2 <- df_2$abundances - df_1$abundances
hist(test2)
qqnorm(test2)
qqline(test2)
shapiro.test(test2) #test de normalité des résidus


# dataframe for visualisation

cbind(df_1, as.data.frame(c(rep(1, dim(df_1)[1]))))
df_1_2 <- rbind(cbind(df_1, as.data.frame(list(rep(1, dim(df_1)[1])), col.names = "dataset" )), cbind(df_2, as.data.frame(list(rep(2, dim(df_2)[1])), col.names = "dataset")))
df_1_2$dataset <- as.factor(df_1_2$dataset)

# visulatisation of relative abondances to compare after Noise addition
ggplot(data = df_1_2, mapping = aes(x=genomes, abundances, color = dataset)) + geom_point(alpha = 0.5) + labs(title = "Relative abundance", y = "abundances", x = "genomes") + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1))

df_1_2 <- df_1_2[-c(which(df_1_2$genomes == "Cutibacterium acnes")),]


ggplot(data = df_1_2, mapping = aes(x=genomes, abundances, color = dataset)) + geom_point(alpha = 0.5) + labs(title = "Relative abundance, without C.acnes", y = "abundances", x = "genomes") + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1))

##########################################################
#
#     REPLICATS
#
##########################################################

df <- read.table(file = "./SimulatedDataForMETADIFF/taxonomy/abondances/acne_150_2022.tsv", sep = "\t", header = FALSE)

#ajout du FC
for (j in c(1:5)) {
    for (i in c(1:10)) {
      df_2 <- Noise(df)
      write.table(df_2, file = paste("./rep",j,"/noise/acne_150_rep_", i, ".tsv", sep = ""), sep = "\t", row.names = FALSE, col.names = FALSE)
    }
    
      #CF3
    for (i in c(1:10)) {
      table <- read.table(file = paste("./rep", j,"/noise/acne_150_rep_", i, ".tsv", sep = ""), sep = "\t", header = FALSE)
      table2 <- DE(table, fold_change = c(3,3,3,3), target = c("Cutibacterium_acnes", "Staphylococcus_epidermidis", "Malassezia_restricta", "Streptococcus_pyogenes"))
      write.table(table2, file = paste("./rep",j,"/FC3/acne_150_rep_", i, "_FC3.tsv", sep = ""), sep ="\t", row.names = FALSE, col.names = FALSE)
    }
    
    #FC2
    for (i in c(1:10)) {
      table <- read.table(file = paste("./rep", j,"/noise/acne_150_rep_", i, ".tsv", sep = ""), sep = "\t", header = FALSE)
      table2 <- DE(table, fold_change = c(2,2,2,2), target =  c("Cutibacterium_acnes", "Staphylococcus_epidermidis", "Malassezia_restricta", "Streptococcus_pyogenes"))
      write.table(table2, file = paste("./rep",j,"/FC2/acne_150_rep_", i, "_FC2.tsv", sep = ""), sep ="\t", row.names = FALSE, col.names = FALSE)
    }
    
    #FC1.5
    for (i in c(1:10)) {
      table <- read.table(file = paste("./rep", j,"/noise/acne_150_rep_", i, ".tsv", sep = ""), sep = "\t", header = FALSE)
      table2 <- DE(table, fold_change = c(1.5,1.5,1.5,1.5), target =  c("Cutibacterium_acnes", "Staphylococcus_epidermidis", "Malassezia_restricta", "Streptococcus_pyogenes"))
      write.table(table2, file = paste("./rep", j, "/FC15/acne_150_rep_", i, "_FC15.tsv", sep = ""), sep ="\t", row.names = FALSE, col.names = FALSE)
    }
    
    #FC1.4
    for (i in c(1:10)) {
      table <- read.table(file = paste("./rep", j,"/noise/acne_150_rep_", i, ".tsv", sep = ""), sep = "\t", header = FALSE)
      table2 <- DE(table, fold_change = c(1.4,1.4,1.4,1.4), target =  c("Cutibacterium_acnes", "Staphylococcus_epidermidis", "Malassezia_restricta", "Streptococcus_pyogenes"))
      write.table(table2, file = paste("./rep", j, "/FC14_40/acne_150_rep_", i, "_FC14.tsv", sep = ""), sep ="\t", row.names = FALSE, col.names = FALSE)
    }
    
    #FC1.3
    for (i in c(1:10)) {
      table <- read.table(file = paste("./rep", j, "/noise/acne_150_rep_", i, ".tsv", sep = ""), sep = "\t", header = FALSE)
      table2 <- DE(table, fold_change = c(1.3,1.3,1.3,1.3), target =  c("Cutibacterium_acnes", "Staphylococcus_epidermidis", "Malassezia_restricta", "Streptococcus_pyogenes"))
      write.table(table2, file = paste("./rep", j, "/FC13_40/acne_150_rep_", i, "_FC13.tsv", sep = ""), sep ="\t", row.names = FALSE, col.names = FALSE)
    }
    
    #FC1.2
    for (i in c(1:10)) {
      table <- read.table(file = paste("./rep", j, "/noise/acne_150_rep_", i, ".tsv", sep = ""), sep = "\t", header = FALSE)
      table2 <- DE(table, fold_change = c(1.2,1.2,1.2,1.2), target =  c("Cutibacterium_acnes", "Staphylococcus_epidermidis", "Malassezia_restricta", "Streptococcus_pyogenes"))
      write.table(table2, file = paste("./rep", j, "/FC12_40/acne_150_rep_", i, "_FC12.tsv", sep = ""), sep ="\t", row.names = FALSE, col.names = FALSE)
    }

}

#################################################
#
#     DEPLETION + ENRICHISSEMENT
#
#################################################

for ( j in c(0,1,2,3,4)) {
  for (i in c(1:10)) {
    table <- read.table(file = paste("./rep", j,"/noise/acne_150_rep_", i, ".tsv", sep = ""), sep = "\t", header = FALSE)
    table <- DE(table, fold_change = c(3,1/3), target = c("Cutibacterium_acnes", "Malassezia_restricta"))
    write.table(table, file = paste("./depletion_enrichment/rep",j,"/acne_150_rep_", i, "_deplete_enriched.tsv", sep = ""), sep ="\t", row.names = FALSE, col.names = FALSE)
  }
}
