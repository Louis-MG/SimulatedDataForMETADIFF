#!/usr/bin/env Rscript

# args : outputdir
# usage : Rscript metrics.R ~/output_bifrost

#########################################
#
#     INITIALISATION
#
#########################################

library(caret)

#defaultW <- getOption("warn") 
#options(warn = -1)

args = commandArgs()

setwd(args[6])
results_file <- "./result_genomes_bifrost_query.tsv"
df = read.table(file = results_file, sep = "\t", header = TRUE)

genomes <- c("CutibacteriumAcnes", "Malassezia_restricta", "Staphylococcus_epidermidis", "Streptococcus_pyogenes")
uniques <- c("unique_kmers_CutibacteriumAcnes", "unique_kmers_Malassezia_restricta", "unique_kmers_Staphylococcus_epidermidis", "unique_kmers_Streptococcus_pyogenes")

#########################################
#
#     FUNCTIONS
#
#########################################

verif_results <- function(dataframe, vector_genomes, vector_uniques) {
  #' dataframe: a dataframe
  #' vector_genomes: a vetcor with genome of interest's names, as present in the result file
  #' vector_uniques: a vetcor with unique_genome names, as present in the result file
  #' This function checks that the file is correctly read as a dataframe, that the kmers found to be unique
  #' to a genome really are.  
  if (is.data.frame(dataframe) == FALSE) {
    stop("df must be a dataframe.")
  }  
  if (sum(rowSums(dataframe[, 2:155]) == 0) != 0) {
      not_found <- which(rowSums(dataframe[, 2:155]) == 0)
      print(dataframe[not_found,"query_name"])
      warning("some kmers are not found in any genomes.")
    }
  if (length(which(rowSums(dataframe[, vector_uniques]) > 1)) != 0 ) {
    not_uniques <- which(rowSums(dataframe[, vector_uniques]) > 1 )
    print(not_uniques)
    stop("some unique kmers are not unique ! See index given above.")
  }
  for (i in c(1:4)) {
    index <- which(dataframe[,vector_uniques[i]] == 1)
    if (length(which(dataframe[index, vector_genomes[i]] != 1 )) != 0 ) { #checks if all kmers found in unique_genome_X are also found in genome_X      
      index <- which(dataframe[index, vector_genomes[i]] != 1 )
      print(index)
      stop(paste("some kmers found in unique are not found in the genome of reference ", i, ". See index given above.", sep = ""))
    }
  }
}

graphs <- function(dataframe, vector_genomes, vector_uniques) {
  #' dataframe: a dataframe
  #' vector_genomes: a vetcor with genome of interest's names, as present in the result file
  #' vector_uniques: a vetcor with unique_genome names, as present in the result file
  #' This function produces a graph with 4 panels: specificity, sensibility, recall/precision, and auc/aucpr.
  #' All are saved as png.
  for (i in c(1:4)) {
    tryCatch(expr = {#png(filename = paste("graphs_",vector_genomes[i],".png", sep = ""), width = 1396, height = 791)
                     #par(mfrow=c(2,2))
                     #pred <- prediction(dataframe[,vector_genomes[i]], dataframe[,vector_uniques[i]])
                     #perf <- performance(pred, "tpr", "fpr")
                     #plot(perf, main = vector_genomes[i], colorize = TRUE)
                     #perf <- performance(pred, "prec", "rec")
		     #perf <- performance(pred, "prec")
		     #prec <- round(perf@y.values[[1]], digits = 3)
		     #print(perf@y.values)
                     #plot(perf, main = vector_genomes[i], colorize = FALSE)
	             #perf <- performance(pred, "sens", "spec")
                     #plot(perf, main = vector_genomes[i], colorize = TRUE)
                     #auc <- performance(pred, measure = "auc")
                     #auc <- round(auc@y.values[[1]], digits = 3)
                     #aucpr <- performance(pred, measure = "aucpr")
                     #aucpr <- round(aucpr@y.values[[1]], digits = 3)
                     #plot.new()
                     #text(c(0.45, 0.52, 0.45, 0.52), labels = c("prec = ", prec ))
                     #dev.off()
		     pred <- factor(dataframe[,vector_genomes[i]], levels = c('1','0'))
                     labels <- factor(dataframe[,vector_uniques[i]], levels = c('1','0'))
                     result <- confusionMatrix(pred, labels, mode="prec_recall")
                     precision <- result$byClass["Precision"]
                     print(paste("Results for", vector_genomes[i], ":", sep = " "))
		     print(result)
                     },
              error = function(e){          # Specifying error message
                message(message(paste("Calculation", i, "finished with error:", e$message, sep = " ")))
              },
              warning = function(w){        # Specifying warning message
                message(paste("Calculation", i, "finished with message:", w$message, sep = " "))
              }
            )
  }
}

#########################################
#
#     USE
#
#########################################


verif_results(df, genomes, uniques)
graphs(df, genomes, uniques)