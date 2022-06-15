##### Libraries #####
library('PhosR')
library('limma')
library('statmod')
library('dplyr')
library('lme4')
library('tidyverse')
library('naniar')
library('ggpubr')
library('Biobase')
library("AnnotationDbi")
library("org.Hs.eg.db")
library('ggnewscale')
library('ggridges')
library('enrichplot')
library('ggrepel')
library('dplyr')
library("ggfortify")
library("devtools")
library('readxl')
library('otuSummary')
library('stringi')
library('clusterProfiler')
library('fgsea')
library('reshape2')
library('ggvenn')
library('ggpubr')
library('Rtsne')
library('ggraph')
library("stringi")
library('ReactomePA')
library('rmcorr')
library('tidygraph')
library('eulerr')
library('corrplot')
library('dendextend')
library('foreach')
library('parallel')
library('doParallel')
setwd("/Users/Gerard/Desktop/")
#### Functions ####
#' Z normalization of the data by column
#' 
#' @param input_mat this is a matrix of numeric data and its columns will be z-normalized
#' @return matrix z-normalized by columns
Znorm <- function(input_mat){
  Output_mat = data.frame()
  for(i in 1:ncol(input_mat)){
    col = input_mat[,i]
    mean_col = mean(col,na.rm = T)
    sd_col = sd(col, na.rm = T)
    Output_mat <- rbind(Output_mat,(col-mean_col)/sd_col)
  }
  return(t(Output_mat))
} 

#' Xiao correction
#' 
#' @param matrix_limma this is the resulting table from the topTable limma function with 1 coefficient
#' @result a vector with the Xiao corrected p-values
Xiao_correction <- function(matrix_limma){
  return(matrix_limma[,4] ** abs(matrix_limma[,1]))
}

#' 1D annotation enrichment analysis parallelized
#' 
#' @param matrix_values a matrix with the following columns 1st Annotation, 2nd Protein Name(protein has to be repeated as many times as annotations it has)
#' @param Log_vec a matrix with the first colomn having the protein and the second the values to be enriched (LogFC, pValue, etc.) 
#' @param pval_cutoff a number representing the p-value cutoff to be used
#' @return a matrix with the first column having the enrichment scores for the first group, the second column has the annotation name
#' @example 
#' enrichment_1D_parallel(matrix_annotations_GOBP, r_to_1D)
enrichment_1D_parallel <- function(matrix_values, Log_vec){
  annotations = unique(matrix_values[,1])
  n.cores <- parallel::detectCores() - 1
  my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
  print(my.cluster)
  doParallel::registerDoParallel(cl = my.cluster)
  foreach::getDoParRegistered()
  foreach::getDoParWorkers()

  sigWil <- foreach(anot = annotations) %dopar% {
    make_groups <- matrix_values[matrix_values[,1] == anot,2]
    index <- which(Log_vec[,1] %in% make_groups)
    wilcox <- wilcox.test(as.numeric(Log_vec[index,2]), as.numeric(Log_vec[-index,2]))
    return(c(anot, wilcox$p.value))
  }

  sigWil <- as.data.frame(t(rbind.data.frame(sigWil)))
  colnames(sigWil) <- c("Annotation", "pVal")
  sigWil$p.adj <- p.adjust(sigWil$pVal, method = "BH")
  sigWil <- sigWil[sigWil$p.adj <= 0.1,]

  s <- c()
  for(anot in sigWil[,1]){
    make_groups <- matrix_values[matrix_values[,1] == anot,2]
    index <- which(Log_vec[,1] %in% make_groups)
    ranked_data <- rank(as.numeric(Log_vec[,2]))
    s.calc <- 2*(mean(ranked_data[index]) - mean(ranked_data[-index]))/length(ranked_data)
  s <- append(s,s.calc)
  }
  return(cbind(sigWil, s))
}

#' 2D annotation enrichment analysis parallelized
#' 
#' @param matrix1 a matrix with the following columns 1st Annotation, 2nd Protein Name(protein has to be repeated as many times as annotations it has)
#' @param Log_vec_1 a matrix with the first colomn having the protein and the second the values  (LogFC, pValue, etc.) 
#' @param Log_vec_2 a matrix with the first colomn having the protein and the second the values  (LogFC, pValue, etc.) 
#' @param pval_cutoff a number representing the p-value cutoff to be used
#' @return a matrix with the first column having the enrichment scores for the first group, the second column has the enrichment scores for the second group, the third column has the annotation name
#' @example 
#' Enrichment_2D_parallel(matrix_annotations_GOCC, LogLvsT, LogOvsT, 0.1)
Enrichment_2D_parallel <- function(matrix1,Log_vec_1,Log_vec_2,pval_cutoff){
  annotations <- unique(matrix1[,1])
  
  n.cores <- parallel::detectCores() - 1
  my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
  print(my.cluster)
  doParallel::registerDoParallel(cl = my.cluster)
  foreach::getDoParRegistered()
  foreach::getDoParWorkers()
  
  sigManova <- foreach(anot = annotations) %dopar% {
    make_groups <- matrix1[matrix1[,1] == anot,2]
    Group_p <- numeric(nrow(Log_vec_1))
    Group_p[which(Log_vec_1[,1] %in% make_groups)] <- 1
    Group_p <- as.factor(Group_p)
    Data <- cbind(rank(as.numeric(Log_vec_1[,2])),rank(as.numeric(Log_vec_2[,2])),Group_p)
    res.manova <- manova(cbind(V1,V2) ~ Group_p, as.data.frame(Data))
    summary.man <- summary.aov(res.manova)
    c(anot, summary.man[[1]][1,5], summary.man[[2]][1,5])
  }
  sigManova <- as.data.frame(t(rbind.data.frame(sigManova)))
  colnames(sigManova) <- c("Annotation", "pVal_1", "pVal_2")
  sigManova$p.adj_1 <- p.adjust(as.numeric(sigManova$pVal_1), method = "BH")
  sigManova$p.adj_2 <- p.adjust(as.numeric(sigManova$pVal_2), method = "BH")
  sigManova_short <- sigManova %>% filter(p.adj_1 <= pval_cutoff | p.adj_2 <= pval_cutoff)
  
  sx <- c()
  sy <- c()
  for(anot in sigManova_short[,1]){
    make_groups <- matrix1[matrix1[,1] == anot,2]
    index = which(Log_vec_1[,1] %in% make_groups)
    ranked_data <- cbind(rank(as.numeric(Log_vec_1[,2])), rank(as.numeric(Log_vec_2[,2])))
    s.calc <- 2*(mean(as.numeric(ranked_data[index,1])) - mean(as.numeric(ranked_data[-index,1])))/nrow(Log_vec_1)
    sx <- append(sx,s.calc)
    s.calc <- 2*(mean(as.numeric(ranked_data[index,2])) - mean(as.numeric(ranked_data[-index,2])))/nrow(Log_vec_2)
    sy <- append(sy,s.calc)
  }
  Out_matrix <- as.data.frame(cbind(x = as.numeric(sx), y = as.numeric(sy), annotation = sigManova_short[,1]))
  Out_matrix$x <- as.numeric(Out_matrix$x)
  Out_matrix$y <- as.numeric(Out_matrix$y)
  return(Out_matrix)
}

#### Preprocessing ####
#Setting Color variable
color_plots <- c("#440154FF", "#73D055FF","#404788FF")
### Importing data
#Clinical data
clinical_data <- as.data.frame(read_excel("Gerard_HIIT_data/Gerard_Clinical_data.xlsx"))
clinical_data_noNA <- na.omit(clinical_data)

#Proteomics data
Adipose_proteome <- read.delim("Gerard_HIIT_data/HIIT_adipose_DIA_NN_Library_based.pg_matrix.tsv", 
                               header = TRUE, sep = "\t")
#Read grouping and sample ID
Name_adipose <- read.delim("Gerard_HIIT_data/KH_study_names_for_R.txt", header = F)
Name_adipose_2 <- t(Name_adipose$V2)

colnames(Adipose_proteome)[6:ncol(Adipose_proteome)] <- Name_adipose_2
Exprs_adipose <- Adipose_proteome[,6:ncol(Adipose_proteome)]

### Removing outlier - Only 1.200 proteins quantified in that sample.
Exprs_adipose <- Exprs_adipose[,-grep("T2D_Pre_training_12", colnames(Exprs_adipose))]

### Importing batch-cluster information. Batch effects were found with hierarchical clustering with euclidean distances 
Input_cluster_correction <- read.delim("Gerard_HIIT_data/INPUT_CLUSTER_CORRECTION.txt", header = TRUE) 
reorder_idx <- match(colnames(Exprs_adipose),Input_cluster_correction$Sample_ID)
Input_cluster_correction <- Input_cluster_correction[reorder_idx,]

#Log2 transform expression data
Exprs_adipose = as.matrix(log2(Exprs_adipose))
rownames(Exprs_adipose) <- Adipose_proteome$Protein.Ids

# Defining variables
sample_name = strsplit(gsub("^_", "", colnames(Exprs_adipose)), "_")
df = S4Vectors::DataFrame(
  Group = sapply(sample_name, "[[", 1),
  Condition = sapply(sample_name, "[[", 2),
  replicate = sapply(sample_name, "[[", 4))
rownames(df) = colnames(Exprs_adipose)

subj <- as.factor(df$replicate)
grps <- paste0(df$Group)
Treatment <- paste0(df$Condition)
combined <- paste0(df$Group, df$Condition)
cluster <- as.factor(Input_cluster_correction$Cluster)

#Sample quality
Exprs_adipose_no_filter <- Exprs_adipose
vec <- c()
name <- colnames(Exprs_adipose_no_filter)
for(i in seq_along(as.data.frame(Exprs_adipose_no_filter))){
  len <- length(na.omit(Exprs_adipose_no_filter[,i]))
  vec <-append(vec,len)
}
ggplot(mapping = aes(1:91,y = vec, fill = grps)) + 
  geom_col(alpha = 0.8) +
  theme_minimal() + 
  theme(axis.text = element_text(size= 18),
        axis.line = element_line(color="black", size = 1),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.title = element_text(size = 20), #change legend title font size
        legend.text = element_text(size=15),
        panel.grid = element_blank(),
        axis.title = element_text(size = 30)) +  #change legend text font size
  scale_fill_manual("Groups",values = color_plots) +
  scale_x_continuous(breaks = seq(from = 0, to = 90, by = 5)
                     ,limits = c(0,92), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0,3000), expand = c(0,0)) +
  labs(x = "Sample", y = "Quantified proteins") 

#FILTERING
Exprs_adipose_ranking <- selectGrps(Exprs_adipose, combined, 0.02, n=1)
Exprs_adipose <- selectGrps(Exprs_adipose, combined, 0.5, n=6)
Exprs_adipose_to_PCA <- selectGrps(Exprs_adipose, combined, 1, n=6) # for PCA plot

dim(Exprs_adipose)


#Median normalize
data_median <- apply(Exprs_adipose, 2, median, na.rm=TRUE)
#Exprs_adipose_notNotmalized <- Exprs_adipose
#Exprs_adipose_notImputed <- Exprs_adipose[] - data_median[col(Exprs_adipose)[]]
Exprs_adipose_notImputed <- medianScaling(Exprs_adipose)
Exprs_adipose_scaled <- Exprs_adipose_notImputed 
Exprs_adipose_ranking_scaled <- medianScaling(Exprs_adipose_ranking)

### Removing cluster-batch effect. Only for PCA. 
Group <- factor(grps, levels=c("Lean","Obese", "T2D"))
Training <- factor(Treatment, levels=c("Pre","Post"))
design <- model.matrix(~ 0 + Group*Training)
Exprs_adipose_noBatch_notImp <- removeBatchEffect(Exprs_adipose_notImputed, batch = cluster, design = design)
Exprs_adipose_scaled <- removeBatchEffect(Exprs_adipose_scaled, batch = cluster, design = design)
Exprs_adipose_ranking_scaled_nobatch <- removeBatchEffect(Exprs_adipose_ranking_scaled, batch = cluster, design = design)

data_median <- apply(Exprs_adipose_to_PCA, 2, median, na.rm=TRUE)
Exprs_adipose_to_PCA <- Exprs_adipose_to_PCA[] - data_median[col(Exprs_adipose_to_PCA)[]]
Exprs_adipose_to_PCA <- removeBatchEffect(Exprs_adipose_to_PCA, batch = cluster, design = design)
  
pca_1 <- prcomp(t(Exprs_adipose_to_PCA), scale. = T)
ggplot(mapping = aes(pca_1$x[,1],pca_1$x[,2], color = grps)) + 
  geom_point(size = 2) +
  stat_ellipse(lwd = 0.5, show.legend = F) +
  theme_minimal() +
  theme(axis.line = element_line(),
        axis.text = element_text(size= 18),
        axis.title = element_text(size = 30),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.title = element_text(size=20), #change legend title font size
        legend.text = element_text(size=15), 
        panel.grid = element_blank() ) + 
  labs(x = paste("PC1", round(pca_1$sdev[1]/sum(pca_1$sdev)*100,2),"%", sep = " "), 
       y = paste("PC2", round(pca_1$sdev[2]/sum(pca_1$sdev)*100,2),"%", sep = " ")) + 
  scale_color_manual("Groups",values = color_plots)

geneSymbols <- mapIds(org.Hs.eg.db, keys=rownames(Exprs_adipose_ranking_scaled_nobatch), column="SYMBOL", keytype="ACCNUM", multiVals="first")

protein_median <- apply(Exprs_adipose_ranking_scaled_nobatch,1,FUN = function(x) median(x,na.rm = T))
protein_median <- as.data.frame(cbind(genes = geneSymbols,medians = protein_median))
protein_median <- protein_median[order(as.numeric(protein_median[,2]),decreasing = T),]
protein_median <- as.data.frame(cbind(protein_median, rank =1:3282))

ggplot() +
  geom_point(data = protein_median, 
             mapping = aes(as.numeric(rank), as.numeric(medians)),
             size = 2) +
  geom_point(data = protein_median %>% 
               filter(genes %in% c("HBB", "HBA1", "ALB")), 
             mapping = aes(as.numeric(rank), as.numeric(medians)), 
             color = "red",
             size = 4) +
  geom_point(data = protein_median %>% 
               filter(genes %in% c("FABP4", "PLIN1", "VCL","RAB10","SLC2A4","ADIPOQ","LEP")), 
             mapping = aes(as.numeric(rank), as.numeric(medians)), 
             color = "green",
             size = 4) +
  theme_minimal() + 
  theme(axis.line = element_line(),
        axis.text = element_text(size= 18),
        axis.title = element_text(size = 30), 
        panel.grid = element_blank())+
  labs(x = "Protein abundances ranked", y = "Log2 Intensities") +
  scale_x_continuous(limits = c(-10,3295), 
                     expand = c(0,0)) +
  geom_text_repel(data = protein_median, 
                  aes(x = as.numeric(rank), 
                      y = as.numeric(medians),
                      label=ifelse(genes %in% c("HBB","HBA1","ALB","FABP4","PLIN1", "ADIPOQ","VCL","RAB10","SLC2A4","LEP"),
                                   as.character(protein_median[,1]),'')), 
                  max.overlaps = 3000,
                  size = 8) 
  
#imputation
set.seed(123)
Exprs_adipose_imputed_allLow <- tImpute(Exprs_adipose_notImputed, m=1.6, s=0.6)
Exprs_adipose <- Exprs_adipose_notImputed
#### Plot Imputation ####

Sample6_select = c(60,62,72,27,37,1)
Na_Exprs <- is.na(Exprs_adipose_notImputed[,Sample6_select])
Exprs_adipose_6_I = Exprs_adipose_imputed_allLow[,Sample6_select]

barplot1 <- ggplot() + geom_histogram(aes(Exprs_adipose_6_I[,1], fill = Na_Exprs[,1], alpha = 1000),position = "stack", color = "black") + 
  theme_minimal() + 
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.line = element_line()) + 
  scale_fill_manual(values = c("blue","red")) +
  labs(x = "Intensities", y = "Count")

barplot2 <- ggplot() + geom_histogram(aes(Exprs_adipose_6_I[,2], fill = Na_Exprs[,2], alpha = 1000),position = "stack", color = "black") + 
  theme_minimal() + 
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.line = element_line()) + 
  scale_fill_manual(values = c("blue","red"))  +
  labs(x = "Intensities", y = "")

barplot3 <- ggplot() + geom_histogram(aes(Exprs_adipose_6_I[,3], fill = Na_Exprs[,3], alpha = 1000),position = "stack", color = "black") + 
  theme_minimal() + 
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.line = element_line()) + 
  scale_fill_manual(values = c("blue","red")) +
  labs(x = "Intensities", y = "")

barplot4 <- ggplot() + geom_histogram(aes(Exprs_adipose_6_I[,4], fill = Na_Exprs[,4], alpha = 1000),position = "stack", color = "black") + 
  theme_minimal() + 
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.line = element_line()) + 
  scale_fill_manual(values = c("blue","red")) +
  labs(x = "Intensities", y = "Count")

barplot5 <- ggplot() + geom_histogram(aes(Exprs_adipose_6_I[,5], fill = Na_Exprs[,5], alpha = 1000),position = "stack", color = "black") + 
  theme_minimal() + 
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.line = element_line()) + 
  scale_fill_manual(values = c("blue","red")) +
  labs(x = "Intensities", y = "")

barplot6 <- ggplot() + geom_histogram(aes(Exprs_adipose_6_I[,6], fill = Na_Exprs[,6], alpha = 1000),position = "stack", color = "black") + 
  theme_minimal() + 
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.line = element_line()) + 
  scale_fill_manual(values = c("blue","red")) +
  labs(x = "Intensities", y = "")

ggarrange(barplot1, barplot2, barplot3, barplot4, barplot5, barplot6)
#### Plot Not imputed ####
barplot1 <- ggplot() + geom_histogram(aes(Exprs_adipose_notImputed[,Sample6_select][,1], alpha = 1000),
                                      position = "stack", color = "black", fill = "blue") + 
  theme_minimal() + 
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.line = element_line()) + 
  scale_y_continuous(limits = c(0,250), expand = c(0,0)) +
  labs(y = "Count", x = "Intensities")

barplot2 <- ggplot() + geom_histogram(aes(Exprs_adipose_notImputed[,Sample6_select][,2], alpha = 1000),position = "stack", color = "black", fill = "blue") + 
  theme_minimal() + 
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.line = element_line()) + 
  scale_y_continuous(limits = c(0,270), expand = c(0,0)) +
  labs(y = "", x = "Intensities")

barplot3 <- ggplot() + geom_histogram(aes(Exprs_adipose_notImputed[,Sample6_select][,3], alpha = 1000),position = "stack", color = "black", fill = "blue") + 
  theme_minimal() + 
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.line = element_line()) + 
  scale_y_continuous(limits = c(0,250), expand = c(0,0)) +
  labs(y = "", x = "Intensities")

barplot4 <- ggplot() + geom_histogram(aes(Exprs_adipose_notImputed[,Sample6_select][,4], alpha = 1000),position = "stack", color = "black", fill = "blue") + 
  theme_minimal() + 
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.line = element_line()) + 
  scale_y_continuous(limits = c(0,250), expand = c(0,0)) +
  labs(y = "Count", x = "Intensities")

barplot5 <- ggplot() + geom_histogram(aes(Exprs_adipose_notImputed[,Sample6_select][,5], alpha = 1000),position = "stack", color = "black", fill = "blue") + 
  theme_minimal() + 
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.line = element_line()) + 
  scale_y_continuous(limits = c(0,250), expand = c(0,0)) +
  labs(y = "", x = "Intensities")

barplot6 <- ggplot() + geom_histogram(aes(Exprs_adipose_notImputed[,Sample6_select][,6], alpha = 1000),position = "stack", color = "black", fill = "blue") + 
  theme_minimal() + 
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.line = element_line()) + 
  scale_y_continuous(limits = c(0,250), expand = c(0,0)) +
  labs(y = "", x = "Intensities")

ggarrange(barplot1, barplot2, barplot3, barplot4, barplot5, barplot6)


qqplot1 <- ggplot(mapping = aes(sample = Exprs_adipose_notImputed[,Sample6_select][,1])) + 
  stat_qq_line(size = 1) +
  stat_qq(aes(alpha = 1000),color = "red") + 
  theme_minimal() + 
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.line = element_line()) +
  labs(x = "Theoretical", y = "Sample")

qqplot2 <- ggplot(mapping = aes(sample = Exprs_adipose_notImputed[,Sample6_select][,2])) + 
  stat_qq_line(size = 1) +
  stat_qq(aes(alpha = 1000),color = "red") + 
  theme_minimal() + 
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.line = element_line()) +
  labs(x = "Theoretical", y = "") 

qqplot3 <- ggplot(mapping = aes(sample = Exprs_adipose_notImputed[,Sample6_select][,3])) + 
  stat_qq_line(size = 1) +
  stat_qq(aes(alpha = 1000),color = "red") + 
  theme_minimal() + 
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.line = element_line()) +
  labs(x = "Theoretical", y = "") 

qqplot4 <- ggplot(mapping = aes(sample = Exprs_adipose_notImputed[,Sample6_select][,4])) + stat_qq_line(size = 1) +
  stat_qq(aes(alpha = 1000),color = "red") + 
  theme_minimal() + 
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.line = element_line()) +
  labs(x = "Theoretical", y = "Sample")

qqplot5 <- ggplot(mapping = aes(sample = Exprs_adipose_notImputed[,Sample6_select][,5])) + 
  stat_qq_line(size = 1) +
  stat_qq(aes(alpha = 1000),color = "red") + 
  theme_minimal() + 
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.line = element_line()) +
  labs(x = "Theoretical", y = "")

qqplot6 <- ggplot(mapping = aes(sample = Exprs_adipose_notImputed[,Sample6_select][,6])) + 
  stat_qq_line(size = 1) +
  stat_qq(aes(alpha = 1000),color = "red") + 
  theme_minimal() + 
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.line = element_line()) +
  labs(x = "Theoretical", y = "")

ggarrange(qqplot1,qqplot2,qqplot3,qqplot4,qqplot5,qqplot6)

#### Reproducibility ####
try <- Exprs_adipose_imputed_allLow[,order(grps)]
colnames(try) <- paste0(combined[order(grps)],"_", df$replicate[order(grps)])

Correlation_matrix <- cor(try, method = "pearson")
corrplot(Correlation_matrix, type = "upper",
         tl.col = "black", tl.srt = 45, col.lim = c(0.5,1), is.corr = F)

corr_LeanPre <- cor(Exprs_adipose_notImputed[,combined == "LeanPre"], method = "pearson", use = "complete.obs")
corr_LeanPre <- corr_LeanPre[upper.tri(corr_LeanPre)]
corr_LeanPost <- cor(Exprs_adipose_notImputed[,combined == "LeanPost"], method = "pearson", use = "complete.obs")
corr_LeanPost <- corr_LeanPost[upper.tri(corr_LeanPost)]
corr_ObesePre <- cor(Exprs_adipose_notImputed[,combined == "ObesePre"], method = "pearson", use = "complete.obs")
corr_ObesePre <- corr_ObesePre[upper.tri(corr_ObesePre)]
corr_ObesePost <- cor(Exprs_adipose_notImputed[,combined == "ObesePost"], method = "pearson", use = "complete.obs")
corr_ObesePost <- corr_ObesePost[upper.tri(corr_ObesePost)]
corr_T2DPre <- cor(Exprs_adipose_notImputed[,combined == "T2DPre"], method = "pearson", use = "complete.obs")
corr_T2DPre <- corr_T2DPre[upper.tri(corr_T2DPre)]
corr_T2DPost <- cor(Exprs_adipose_notImputed[,combined == "T2DPost"], method = "pearson", use = "complete.obs")
corr_T2DPost <- corr_T2DPost[upper.tri(corr_T2DPost)]

corr_mat <- rbind(merge("LeanPre",corr_LeanPre), 
                  merge("LeanPost",corr_LeanPost),
                  merge("ObesePre",corr_ObesePre), 
                  merge("ObesePost",corr_ObesePost),
                  merge("T2DPre",corr_T2DPre), 
                  merge("T2DPost",corr_T2DPost))

ggplot(corr_mat,aes(factor(x, levels = unique(combined)),y)) +
  geom_boxplot() +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(), 
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  labs(x = "Group", y = "Correlation Value (r)") 
  
median_vec <- c()
median_vec <- append(median_vec, median(corr_LeanPre,na.rm = T))
median_vec <- append(median_vec, median(corr_LeanPost,na.rm = T))
median_vec <- append(median_vec, median(corr_ObesePre,na.rm = T))
median_vec <- append(median_vec, median(corr_ObesePost,na.rm = T))
median_vec <- append(median_vec, median(corr_T2DPre,na.rm = T))
median_vec <- append(median_vec, median(corr_T2DPost,na.rm = T))

names(median_vec) <- unique(combined)
median_vec <- append(median_vec, median(c(corr_LeanPre,corr_LeanPost,
                                          corr_ObesePre,corr_ObesePost,
                                          corr_T2DPre,corr_T2DPost),na.rm = T))
names(median_vec)[7] <- "Total"
#### Data Transpose ####
#Transpose the Expression Adipose dataframes:
Exprs_adipose_noBatch_notImp <- t(Exprs_adipose_noBatch_notImp)

notinexps <- setdiff(row.names(Exprs_adipose_noBatch_notImp),clinical_data_noNA$New_ID)
for(name in notinexps){
  Exprs_adipose_noBatch_notImp <- Exprs_adipose_noBatch_notImp[-which(row.names(Exprs_adipose_noBatch_notImp) == name),]
}

#### rmCorr ####
Clinical_Z <- Znorm(clinical_data_noNA[5:58])
Clinical_Z <- as.data.frame(cbind(clinical_data_noNA[1:4], Clinical_Z))
colnames(Clinical_Z) = colnames(clinical_data_noNA)

Exprs_Z_noBatch_NI = Znorm(Exprs_adipose_noBatch_notImp)
colnames(Exprs_Z_noBatch_NI) = colnames(Exprs_adipose_noBatch_notImp)
rownames(Exprs_Z_noBatch_NI) = rownames(Exprs_adipose_noBatch_notImp)
Exprs_Z_noBatch_NI_ordered <- Exprs_Z_noBatch_NI[order(match(rownames(Exprs_Z_noBatch_NI),Clinical_Z$New_ID)),]

result <- data.frame()
for(clinical in c("GIR1","VO2max1", "BMI1", "HbA1c1", "FFM1", "FM1")){
  for(prot in colnames(Exprs_Z_noBatch_NI_ordered)){
    if(length(Exprs_Z_noBatch_NI_ordered[,prot][complete.cases(Exprs_Z_noBatch_NI_ordered[,prot])]) >= 10){
      dataset_new <- as.data.frame(cbind(substring(Clinical_Z$ID,4,5), Clinical_Z[,clinical], as.data.frame(Exprs_Z_noBatch_NI_ordered)[,prot]))
      colnames(dataset_new) <- c("Subject", "Clinical", "Protein")
      dataset_new$Clinical <- as.numeric(dataset_new$Clinical)
      dataset_new$Protein <- as.numeric(dataset_new$Protein)
      test <- rmcorr(participant = Subject, measure1 = Clinical, measure2 = Protein, dataset = dataset_new)
      result <- rbind(result, c(clinical, prot, test$r, test$p,test$CI[1],test$CI[2]))
    }
  }
} 

colnames(result) <- c("clinical", "Protein", "correlation", "pVal", "CI Low", "CI high")
result$correlation <- as.numeric(result$correlation)
result$pVal <- as.numeric(result$pVal)
result$correlation <- as.numeric(result$correlation)

result_FFM <- result %>% filter(result$clinical == "FFM1")
result_FFM$p.adj <- p.adjust(result_FFM$pVal, method = "BH")
prot_FFM <- result_FFM %>% filter(pVal < 0.05) %>% filter(abs(correlation) > 0.4)

result_BMI <- result %>% filter(result$clinical == "BMI1")
result_BMI$p.adj <- p.adjust(result_BMI$pVal, method = "BH")
prot_BMI <- result_BMI %>% filter(pVal < 0.05) %>% filter(abs(correlation) > 0.4)

result_GIR <- result %>% filter(result$clinical == "GIR1")
result_GIR$p.adj <- p.adjust(result_GIR$pVal, method = "BH")
prot_GIR <- result_GIR %>% filter(pVal < 0.05) %>% filter(abs(correlation) > 0.4)

result_FM <- result %>% filter(result$clinical == "FM1")
result_FM$p.adj <- p.adjust(result_FM$pVal, method = "BH")
prot_FM <- result_FM %>% filter(pVal < 0.05) %>% filter(abs(correlation) > 0.4)

result_HbA1c <- result %>% filter(result$clinical == "HbA1c1")
result_HbA1c$p.adj <- p.adjust(result_HbA1c$pVal, method = "BH")
prot_HbA1c <- result_HbA1c %>% filter(pVal < 0.05) %>% filter(abs(correlation) > 0.4)

result_VO2max <- result %>% filter(result$clinical == "VO2max1")
result_VO2max$p.adj <- p.adjust(result_VO2max$pVal, method = "BH")
prot_VO2max <- result_VO2max %>% filter(pVal < 0.05) %>% filter(abs(correlation) > 0.4)


list_prot_sig <- rbind(prot_FFM, prot_BMI, prot_GIR, prot_FM, prot_HbA1c, prot_VO2max)

list_prot_sig$sign <- paste(list_prot_sig$clinical,"-")
list_prot_sig$sign[list_prot_sig$correlation > 0] <- paste(list_prot_sig[list_prot_sig$correlation > 0,1],"+")
table(list_prot_sig$sign)

list1 <- list()
Pair_no_leaf <- c()
for(i in 1:6){
  for(j in i:6){
    if(i != j){
      clii <- unique(list_prot_sig$clinical)[i]
      clij <- unique(list_prot_sig$clinical)[j]
      cli1 <- list_prot_sig %>% filter(clinical == clii)
      cli2 <- list_prot_sig %>% filter(clinical == clij)
      inter <- intersect(cli1$Protein,cli2$Protein)
      if(length(inter) != 0){
        list1[stri_join(clii,clij, sep = " ")] <- paste(inter, collapse = ' ')
        Pair_no_leaf <- append(Pair_no_leaf, inter)
      }
    }
  }
}
Pair_no_leaf <- unique(Pair_no_leaf)

geneSymbols <- mapIds(org.Hs.eg.db, keys=rownames(Exprs_adipose), column="SYMBOL", keytype="ACCNUM", multiVals="first")
geneSymbols[1575] <- "PALM2"#Q8IXS6
geneSymbols[1966] <- "AKAP2"#Q9Y2D5 
result_GIR$names <- NA
result_GIR$names[result_GIR$p.adj < 0.1] <- geneSymbols[result_GIR$p.adj < 0.1]
result_GIR$sig <- "NO"
result_GIR$sig[result_GIR$p.adj < 0.1] <- "0.1"
result_GIR$sig[result_GIR$p.adj < 0.05] <- "0.05"

ggplot(result_GIR, aes(as.numeric(correlation), -log10(pVal), color = sig, label = names)) + 
  geom_point(aes(alpha = 0.4), size = 2) + 
  theme_minimal() + 
  labs(y = "-Log10(P.Value)") + 
  geom_hline(yintercept = -log10(1.5e-4), 
             color = "red") + 
  geom_hline(yintercept = -log10(9e-04), 
             color = "purple") +
  scale_x_continuous(limits = c(-0.6,0.6), 
                     seq(-0.6, 0.6, by = 0.3) , 
                     name = "Correlation Coefficient (r)") +
  scale_color_manual(values = c("#39568CFF", "#73D055FF","gray" ), 
                     breaks = c("0.05","0.1")) + 
  guides(alpha = "none") +
  theme(axis.line = element_line(),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 30),
        axis.ticks = element_line(),
        axis.ticks.length = unit(0.5,"cm"))+
  geom_text_repel(show.legend = F, 
                  size = 5) 

#### Correlation plot Sig ####
resultall <- rbind(result_BMI,result_FFM,result_FM,
                   result_GIR[,1:7],result_HbA1c,result_VO2max)
sigprots <- resultall %>% filter(p.adj <= 0.05) %>% filter(Protein != "P02792") %>% filter(Protein != "P02794")


Group = as.factor(clinical_data_noNA$Group)

notinexps <- setdiff(colnames(Exprs_adipose_scaled),clinical_data_noNA$New_ID)
Exprs_adipose_scaled <- Exprs_adipose_scaled[,-which(colnames(Exprs_adipose_scaled) %in% notinexps)]
Exprs_adipose_scaled <- Exprs_adipose_scaled[,order(match(colnames(Exprs_adipose_scaled),clinical_data_noNA$New_ID))]

  
plot2 <- ggplot(mapping = aes(clinical_data_noNA$GIR1,Exprs_adipose_scaled["O95197",])) +
  geom_point(aes(alpha = 0.9,color = Group, shape = clinical_data_noNA$Condition),size = 3) +  
  theme_minimal() +
  geom_smooth(method = "lm") +
  annotate("text", label = paste("r = -0.572"), x = 600, y = 24, size = 10) + 
  annotate("text", label = paste("adj PVal = 0.0215"), x = 550, y = 23.5, size = 10) +   
  labs(x = "GIR", y = "Log2 Intnsities", title = "RTN3") +
  scale_shape_manual("Timepoint",breaks = c("Pre","Post"), values = c(16,17))+
  scale_color_manual("Groups" ,values = color_plots, labels = c("Lean", "Obese", "T2D")) +
  guides(alpha = "none") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 30),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        title = element_text(size = 30))
plot2

plot3 <- ggplot(mapping = aes(clinical_data_noNA$GIR1,Exprs_adipose_scaled["P04075",])) +
  geom_point(aes(color = Group, alpha = 0.9, shape = clinical_data_noNA$Condition), size = 3) +  
  theme_minimal() +
  geom_smooth(method = "lm") +
  annotate("text", label = paste("r = -0.592"), x = 600, y = 26.75,size = 10) + 
  annotate("text", label = paste("adj PVal = 0.0214"), x = 550, y = 26.5,size = 10) +   
  labs(title = "ALDOA", x = "GIR", y = "Log2 Intnsities") +
  scale_shape_manual("Timepoint",breaks = c("Pre","Post"), values = c(16,17))+
  scale_color_manual("Groups" ,values = color_plots, labels = c("Lean", "Obese", "T2D")) +
  guides(alpha = "none") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 30),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        title = element_text(size = 30)) 
plot3

plot4 <- ggplot(mapping = aes(clinical_data_noNA$GIR1,Exprs_adipose_scaled["P10620",])) +
  geom_point(aes(color = Group, alpha = 0.9, shape = clinical_data_noNA$Condition), size = 3) +  
  theme_minimal() +
  geom_smooth(method = "lm") +
  annotate("text", label = paste("r = -0.590"), x = 600, y = 25.75,size = 10) + 
  annotate("text", label = paste("adj PVal = 0.0214"), x = 550, y = 25.5, size = 10) +   
  labs(title = "MGST1", x = "GIR", y = "Log2 Intnsities") +
  scale_shape_manual("Timepoint",breaks = c("Pre","Post"), values = c(16,17))+
  scale_color_manual("Groups" ,values = color_plots, labels = c("Lean", "Obese", "T2D")) +
  guides(alpha = "none") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 30),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        title = element_text(size = 30))
plot4

plot5 <- ggplot(mapping = aes(clinical_data_noNA$GIR1,Exprs_adipose_scaled["P11310",])) +
  geom_point(aes(color = Group,alpha = 0.9, shape = clinical_data_noNA$Condition), size = 3) + 
  theme_minimal() +
  geom_smooth(method = "lm") +
  annotate("text", label = paste("r = -0.582"), x = 600, y = 25,size = 10) + 
  annotate("text", label = paste("adj PVal = 0.0214"), x = 550, y = 24, size = 10) +   
  labs(title = "ACADM", x = "GIR", y = "Log2 Intnsities") +
  scale_shape_manual("Timepoint",breaks = c("Pre","Post"), values = c(16,17))+
  scale_color_manual("Groups" ,values = color_plots, labels = c("Lean", "Obese", "T2D")) +
  guides(alpha = "none") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 30),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        title = element_text(size = 30))
plot5

plot6 <- ggplot(mapping = aes(clinical_data_noNA$GIR1,Exprs_adipose_scaled["Q99685",])) +
  geom_point(aes(color = Group,alpha = 0.9, shape = clinical_data_noNA$Condition), size = 3) + 
  theme_minimal() +
  geom_smooth(method = "lm") +
  annotate("text", label = paste("r = -0.581"), x = 600, y = 25, size = 10) + 
  annotate("text", label = paste("adj PVal = 0.0214"), x = 550, y = 24.5, size = 10) +   
  labs(title = "MGLL", x = "GIR", y = "Log2 Intnsities") +
  scale_shape_manual("Timepoint",breaks = c("Pre","Post"), values = c(16,17))+
  scale_color_manual("Groups" ,values = color_plots, labels = c("Lean", "Obese", "T2D")) +
  guides(alpha = "none") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 30),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        title = element_text(size = 30))
plot6

plot7 <- ggplot(mapping = aes(clinical_data_noNA$GIR1,Exprs_adipose_scaled["Q9NQC3",])) +
  geom_point(aes(color = Group,alpha = 0.9, shape = clinical_data_noNA$Condition), size = 3) + 
  theme_minimal() +
  geom_smooth(method = "lm") +
  annotate("text", label = paste("r = -0.572"), x = 600, y = 25.5, size = 10) + 
  annotate("text", label = paste("adj PVal = 0.0215"), x = 550, y = 25, size = 10) +   
  labs(title = "RTN4", x = "GIR", y = "Log2 Intnsities") +
  scale_shape_manual("Timepoint",breaks = c("Pre","Post"), values = c(16,17))+
  scale_color_manual("Groups" ,values = color_plots, labels = c("Lean", "Obese", "T2D")) +
  guides(alpha = "none") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 30),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        title = element_text(size = 30))
plot7

ggarrange(plot2, plot3, plot4, plot5, plot6, plot7)

#### Graph ####
#plot using pairwise correlation
nodes <- c(unique(list_prot_sig$clinical),unique(list_prot_sig$Protein))

edges <- data.frame()
for(i in 1:nrow(list_prot_sig)){
  cl <- which(nodes == list_prot_sig$clinical[i])
  pr <- which(nodes == list_prot_sig$Protein[i])
  col <- list_prot_sig$correlation[i]
  edges <- rbind(edges,c(cl,pr,col))
}
nodes <- as.data.frame(nodes)
colnames(nodes) <- c("name")
colnames(edges) <- c("CLINICAL", "PROTEIN", "CORRELATION")
replace <- bitr(nodes[7:241,1], fromType = "ACCNUM", toType = "ALIAS",OrgDb =org.Hs.eg.db)
replace_short <- data.frame()
for(prot in unique(replace$ACCNUM)){
  replace_short <- rbind(replace_short, replace$ALIAS[which(replace$ACCNUM == prot)[1]])
}
nodes[7:241,1] <- NA#replace_short
graph <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)

plot_pair <- ggraph(graph = graph) +
  geom_edge_link(aes(color = edges$CORRELATION < 0, edge_width = abs(edges$CORRELATION))) +
  geom_node_point(size = 1) +                                         
  scale_edge_width(range = c(0.25, 1.5)) +
  geom_node_text(aes(label = name, fontface = "bold"),size = 10,nudge_y = 0, nudge_x = 0,)+
  theme_void() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5) ) +
  scale_edge_color_manual(values = c("#73D055FF", "#1F968BFF"))
plot_pair
#### ORA ####

#Preparing universe
universe_vector <- bitr(colnames(Exprs_adipose_noBatch_notImp), fromType = "ACCNUM", toType = "ENTREZID",OrgDb =org.Hs.eg.db )
universe_short <- data.frame()
for(prot in unique(universe_vector$ACCNUM)){
  universe_short <- rbind(universe_short, universe_vector$ENTREZID[which(universe_vector$ACCNUM == prot)[1]])
}

##PAIRWISE
#Separing the different clusters for the later enrichment
community_list_pairwise <- list()
for(i in seq_along(list_prot_sig$clinical)){
  row = list_prot_sig[i,]
  if(!(row$Protein %in% Pair_no_leaf)){
    if(length(community_list_pairwise[[row$clinical]]) == 0){
      community_list_pairwise[[row$clinical]]= c(row$Protein)
    } else {
      community_list_pairwise[[row$clinical]] = append(community_list_pairwise[[row$clinical]],row$Protein)
    }
  }
}
#Change to entrez ID
for(i in seq_along(community_list_pairwise)){
  community_list_pairwise[[i]] = bitr(community_list_pairwise[[i]], fromType = "ACCNUM", toType = "ENTREZID",OrgDb =org.Hs.eg.db )$ENTREZID
}
#enrichment analysis
enrichment_result_pairwise <- list()
for(i in seq_along(community_list_pairwise)){
  tryCatch(
    expr = {
      enrichment_result_pairwise[[names(community_list_pairwise)[i]]] <- enrichPathway(gene = community_list_pairwise[[i]],
                                                                                       universe = universe_short[,1],
                                                                                       pvalueCutoff = 0.05)@result
    }, error = function(e){
      print("Error")
    }
  )
}

enrichment_result_pairwise_GOCC <- list()
for(i in seq_along(community_list_pairwise)){
  tryCatch(
    expr = {
      enrichment_result_pairwise_GOCC[[names(community_list_pairwise)[i]]] <- enrichGO(gene = community_list_pairwise[[i]],
                                                                                  universe = universe_short[,1],
                                                                                  pvalueCutoff = 0.05,
                                                                                  OrgDb = org.Hs.eg.db,
                                                                                  ont = "CC")@result
    }, error = function(e){
      print("Error")
    }
  )
}

enrichment_result_pairwise_GOMF <- list()
for(i in seq_along(community_list_pairwise)){
  tryCatch(
    expr = {
      enrichment_result_pairwise_GOMF[[names(community_list_pairwise)[i]]] <- enrichGO(gene = community_list_pairwise[[i]],
                                                                                       universe = universe_short[,1],
                                                                                       pvalueCutoff = 0.05,
                                                                                       OrgDb = org.Hs.eg.db,
                                                                                       ont = "MF")@result
    }, error = function(e){
      print("Error")
    }
  )
}

enrichment_result_pairwise_GOBP <- list()
for(i in seq_along(community_list_pairwise)){
  tryCatch(
    expr = {
      enrichment_result_pairwise_GOBP[[names(community_list_pairwise)[i]]] <- enrichGO(gene = community_list_pairwise[[i]],
                                                                                       universe = universe_short[,1],
                                                                                       pvalueCutoff = 0.05,
                                                                                       OrgDb = org.Hs.eg.db,
                                                                                       ont = "BP")@result
    }, error = function(e){
      print("Error")
    }
  )
}
#### Plot ORA ####
pair_df_plot <- data.frame()
for(i in seq_along(enrichment_result_pairwise)){
  clinical <- names(enrichment_result_pairwise)[i]
  pair_df_plot <- rbind(pair_df_plot, merge(enrichment_result_pairwise[[i]], clinical))
}

pair_df_plot <- pair_df_plot %>% filter(p.adjust <= 0.1)
pair_df_plot$Description = factor(pair_df_plot$Description, levels=pair_df_plot[order(pair_df_plot$p.adjust), "Description"])


ggplot(pair_df_plot, aes(1, Description)) + 
  geom_point(aes(color = p.adjust, size = as.factor(Count))) + 
  scale_color_viridis_c() + 
  facet_grid(rows = as.factor(pair_df_plot$y), scales = "free_y", space = "free")+
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + 
  labs(y = "Pathways",x = "", title = "ORA pathway") +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 30),
        axis.text.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        panel.background = element_blank(),
        strip.text.y = element_text(size = 18)) +
  scale_size_discrete("Size", 
                    range = c(5,10), 
                    breaks = c(3,6,8))  
  

#GOMF
pair_df_plot_GOMF <- data.frame()
for(i in seq_along(enrichment_result_pairwise_GOMF)){
  clinical <- names(enrichment_result_pairwise_GOMF)[i]
  pair_df_plot_GOMF <- rbind(pair_df_plot_GOMF, merge(enrichment_result_pairwise_GOMF[[i]], clinical))
}

pair_df_plot_GOMF <- pair_df_plot_GOMF %>% filter(p.adjust <= 0.1)

pair_df_plot_GOMF$Description = factor(pair_df_plot_GOMF$Description, levels=pair_df_plot_GOMF[order(pair_df_plot_GOMF$p.adjust), "Description"])


ggplot(pair_df_plot_GOMF, aes(y, Description)) + 
  geom_point(aes(color = p.adjust, size = as.factor(Count))) + 
  scale_color_viridis_c() + 
  facet_grid(rows = as.factor(pair_df_plot_GOMF$y), scales = "free_x",space = "free") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + 
  labs(y = "GOMF",x = "", title = "ORA GOMF") +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 30),
        axis.text.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        panel.background = element_blank(),
        strip.text.y = element_text(size = 18)) +
  scale_size_discrete("Size", 
                      range = c(5,10), 
                      breaks = c(5,6,8,10))  
#GOBP
pair_df_plot_GOBP <- data.frame()
for(i in seq_along(enrichment_result_pairwise_GOBP)){
  clinical <- names(enrichment_result_pairwise_GOBP)[i]
  pair_df_plot_GOBP <- rbind(pair_df_plot_GOBP, merge(enrichment_result_pairwise_GOBP[[i]], clinical))
}

pair_df_plot_GOBP <- pair_df_plot_GOBP %>% filter(p.adjust <= 0.1)

pair_df_plot_GOBP$Description = factor(pair_df_plot_GOBP$Description, levels=pair_df_plot_GOBP[order(pair_df_plot_GOBP$p.adjust), "Description"])

ggplot(pair_df_plot_GOBP, aes(y, Description)) + 
  geom_point(aes(color = p.adjust, size = as.factor(Count))) + 
  scale_color_viridis_c() + 
  facet_grid(rows = as.factor(pair_df_plot_GOBP$y), scales = "free_x") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + 
  labs(y = "GOBP",x = "", title = "ORA GOBP") +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 30),
        axis.text.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        panel.background = element_blank(),
        strip.text.y = element_text(size = 18)) +
  scale_size_discrete("Size", 
                      range = c(4,10), 
                      breaks = c(4,11,17,28))  

#GOCC
pair_df_plot_GOCC <- data.frame()
for(i in seq_along(enrichment_result_pairwise_GOCC)){
  clinical <- names(enrichment_result_pairwise_GOCC)[i]
  pair_df_plot_GOCC <- rbind(pair_df_plot_GOCC, merge(enrichment_result_pairwise_GOCC[[i]], clinical))
}

pair_df_plot_GOCC <- pair_df_plot_GOCC %>% filter(p.adjust <= 0.1) %>% arrange(p.adjust) 

pair_df_plot_GOCC$Description = factor(pair_df_plot_GOCC$Description, levels=pair_df_plot_GOCC[order(pair_df_plot_GOCC$p.adjust), "Description"])
  
ggplot(pair_df_plot_GOCC, aes(1, as.factor(Description))) + 
  geom_point(aes(color = p.adjust, size = as.factor(Count))) + 
  scale_color_viridis_c() + 
  facet_grid(row = as.factor(pair_df_plot_GOCC$y), scales = "free_y", space = "free") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + 
  labs(y = "GOCC",x = "", title = "ORA GOCC") +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 30),
        axis.text.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        panel.background = element_blank(),
        strip.text.y = element_text(size = 18)) +
  scale_size_discrete("Size", 
                      range = c(4,10), 
                      breaks = c(3,19,28))


#### Annotation ####
matrix_annotations <- read.table(gzfile("/Users/Gerard/Desktop/Annotations/mainAnnot.homo_sapiens.txt.gz"), header = T, sep = "\t", fill = T)
matrix_annotations_GOBP = data.frame()
matrix_annotations_GOMF = data.frame()
matrix_annotations_GOCC = data.frame()

for(Protein.IDs in rownames(Exprs_adipose)){
  id = strsplit(x = Protein.IDs, split = ";")[[1]][1]
  row_uni = grep(id, matrix_annotations$UniProt)
  if(length(row_uni) != 0){
    row_to_annotate = matrix_annotations[row_uni,]
    
    GOBP <- strsplit(x = row_to_annotate$GOBP.name, split = ";")[[1]]
    tojoin <- merge(GOBP,id)
    matrix_annotations_GOBP = rbind(matrix_annotations_GOBP,tojoin)
    
    
    GOCC <- strsplit(x = row_to_annotate$GOCC.name, split = ";")[[1]]
    tojoin <- merge(GOCC,id)
    matrix_annotations_GOCC = rbind(matrix_annotations_GOCC,tojoin)
    
    GOMF <- strsplit(x = row_to_annotate$GOMF.name, split = ";")[[1]]
    tojoin <- merge(GOMF,id)
    matrix_annotations_GOMF = rbind(matrix_annotations_GOMF,tojoin)
  }
}

#### LIMMA####
#Get the logFC using the limma package
Group <- factor(combined, levels=c("LeanPre","LeanPost", "ObesePre", "ObesePost", "T2DPre", "T2DPost"))
Training <- factor(Treatment, levels=c("Pre","Post"))

#The following design matrix allows for initial subtraction of main effect of training.
design <- model.matrix(~ 0 + Group + cluster) # Adding cluster-batch effect as covariate
colnames(design)[7:8] <- c("cluser_one","cluster_two")

corfit <- duplicateCorrelation(Exprs_adipose, design, block=df$replicate)
corfit$consensus

fit <- eBayes(lmFit(Exprs_adipose,design,block=df$replicate,correlation=corfit$consensus))

cm <- makeContrasts(GroupLeanPost - GroupLeanPre, 
                    GroupObesePost - GroupObesePre, 
                    GroupT2DPost - GroupT2DPre, 
                    (GroupLeanPost - GroupLeanPre + GroupObesePost - GroupObesePre + GroupT2DPost - GroupT2DPre)/6,
                    (GroupT2DPost - GroupT2DPre) - (GroupObesePost - GroupObesePre), 
                    (GroupObesePost - GroupObesePre) - (GroupLeanPost - GroupLeanPre), 
                    (GroupT2DPost - GroupT2DPre) - (GroupLeanPost - GroupLeanPre), 
                    ((GroupT2DPost - GroupT2DPre)- (GroupObesePost - GroupObesePre) - (GroupLeanPost - GroupLeanPre))/6,
                    GroupObesePre - GroupLeanPre, 
                    GroupT2DPre - GroupLeanPre, 
                    GroupT2DPre - GroupObesePre,
                    levels=design)
fit2 <- eBayes(contrasts.fit(fit, cm))

Effecttrain_Lean <- topTable(fit2, coef = 1, number = Inf, sort.by = "none")
Effecttrain_Lean$Xiao <- Xiao_correction(Effecttrain_Lean) 

Effecttrain_Obese <- topTable(fit2, coef = 2, number = Inf, sort.by = "none")
Effecttrain_Obese$Xiao <- Xiao_correction(Effecttrain_Obese)

Effecttrain_T2D <- topTable(fit2, coef = 3, number = Inf, sort.by = "none")
Effecttrain_T2D$Xiao <- Xiao_correction(Effecttrain_T2D)

Effecttrain_main <- topTable(fit2, coef = 4, number = Inf, sort.by = "none")
Effecttrain_main$Xiao <- Xiao_correction(Effecttrain_main)

Interaction_OvsT <- topTable(fit2, coef = 5, number = Inf, sort.by = "none")
Interaction_OvsT$Xiao <- Xiao_correction(Interaction_OvsT)

Interaction_LvsO <- topTable(fit2, coef = 6, number = Inf, sort.by = "none")
Interaction_LvsO$Xiao <- Xiao_correction(Interaction_LvsO)

Interaction_LvsT <- topTable(fit2, coef = 7, number = Inf, sort.by = "none")
Interaction_LvsT$Xiao <- Xiao_correction(Interaction_LvsT)

Interaction_main <- topTable(fit2, coef = 8, number = Inf, sort.by = "none")
Interaction_main$Xiao <- Xiao_correction(Interaction_main)

OBESE_vs_LEAN <- topTable(fit2, coef = 9, number = Inf, sort.by = "none")
OBESE_vs_LEAN$Xiao <- Xiao_correction(OBESE_vs_LEAN)
logFC_OvsL <- OBESE_vs_LEAN$logFC

T2D_vs_LEAN <- topTable(fit2, coef = 10, number = Inf, sort.by = "none")
T2D_vs_LEAN$Xiao <- Xiao_correction(T2D_vs_LEAN)
logFC_TvsL <- T2D_vs_LEAN$logFC

T2D_vs_OBESE <- topTable(fit2, coef = 11, number = Inf, sort.by = "none")
T2D_vs_OBESE$Xiao <- Xiao_correction(T2D_vs_OBESE)
logFC_TvsO <- T2D_vs_OBESE$logFC

mainEffect_GROUP <- topTable(fit2, coef = 9:11, number = Inf,sort.by = "none")

#### Dataset to send prepare ####
new.data.send <- Exprs_adipose
new.data.send <- cbind(new.data.send,
                       Effecttrain_Lean$logFC,Effecttrain_Lean$P.Value,Effecttrain_Lean$Xiao,
                       Effecttrain_Obese$logFC,Effecttrain_Obese$P.Value,Effecttrain_Obese$Xiao,
                       Effecttrain_T2D$logFC,Effecttrain_T2D$P.Value,Effecttrain_T2D$Xiao,
                       Effecttrain_main$P.Value, Effecttrain_main$adj.P.Val)
colnames(new.data.send)[92:102] <- c("Effect.Training.Lean.LogFC","Effect.Training.Lean.Pval","Effect.Training.Lean.Xiao correction",
                                     "Effect.Training.Obese.LogFC","Effect.Training.Obese.Pval","Effect.Training.Obese.Xiao correction",
                                     "Effect.Training.T2D.LogFC","Effect.Training.T2D.Pval","Effect.Training.T2D.Xiao correction",
                                     "Main.effect.training.Pval", "Main.effect.training.adj Pval")

new.data.send <- cbind(new.data.send, 
                       OBESE_vs_LEAN$logFC, OBESE_vs_LEAN$P.Value, OBESE_vs_LEAN$Xiao,
                       T2D_vs_LEAN$logFC, T2D_vs_LEAN$P.Value, T2D_vs_LEAN$Xiao,
                       T2D_vs_OBESE$logFC, T2D_vs_OBESE$P.Value, T2D_vs_OBESE$Xiao,
                       mainEffect_GROUP$P.Value,mainEffect_GROUP$adj.P.Val)
colnames(new.data.send)[103:113] <- c("Effect.Groups.Obese.vs.Lean.LogFC","Effect.Groups.Obese.vs.Lean.Pval","Effect.Groups.Obese.vs.Lean.Xiao correction",
                                      "Effect.Groups.T2D.vs.Lean.LogFC","Effect.Groups.T2D.vs.Lean.Pval","Effect.Groups.T2D.vs.Lean.Xiao correction",
                                      "Effect.Groups.T2D.vs.Obese.LogFC", "Effect.Groups.T2D.vs.Obese.Pval", "Effect.Groups.T2D.vs.Obese.Xiao correction",
                                      "Main.Effect.Groups.Pval", "Main.Effect.Groups.adj Pval")

new.data.send <- cbind(new.data.send,
                       Interaction_LvsO$logFC,Interaction_LvsO$P.Value,Interaction_LvsO$Xiao,
                       Interaction_LvsT$logFC,Interaction_LvsT$P.Value,Interaction_LvsT$Xiao,
                       Interaction_OvsT$logFC,Interaction_OvsT$P.Value,Interaction_OvsT$Xiao,
                       Interaction_main$P.Value,Interaction_main$adj.P.Val)
colnames(new.data.send)[114:124] <- c("Interaction.Obese.vs.Lean.LogFC","Interaction.Obese.vs.Lean.Pval","Interaction.Obese.vs.Lean.Xiao correction",
                                       "Interaction.T2D.vs.Lean.LogFC","Interaction.T2D.vs.Lean.Pval","Interaction.T2D.vs.Lean.Xiao correction",
                                       "Interaction.T2D.vs.Obese.LogFC","Interaction.T2D.vs.Obese.Pval","Interaction.T2D.vs.Obese.Xiao correction",
                                       "Interaction.main.Pval","Interaction.main.adj Pval")

new.data.send <- cbind(new.data.send, row.names(Exprs_adipose_imputed_allLow), geneSymbols)
colnames(new.data.send)[125:126] <- c("Protein Names", "Gene Names") 

namesprot <- cbind(clinical_data[,2], paste(clinical_data$ID,clinical_data$Condition, sep = "_"))
namesprot <- namesprot[order(match(namesprot[,1], colnames(new.data.send)[1:91])),]
colnames(new.data.send)[1:91] <- namesprot[1:91,2]

write.table(new.data.send,"HIIT_adipose_Allinfo.tsv",append = F, sep = "\t",dec = ".", row.names = F, quote = F)

#### Plots LIMMA ####

#Changing Protein accession ID's to Gene symbols:
geneSymbols <- mapIds(org.Hs.eg.db, keys=rownames(Exprs_adipose), column="SYMBOL", keytype="ACCNUM", multiVals="first")

# Protein accession ID Q8IXS6 and Q9Y2D5 are mapped by madIds function to "PALM2AKAP2"
# However could as well be PALM2 and AKAP2 respectively. To simplify further analyses these
# will be changed.
geneSymbols[1575] <- "PALM2"#Q8IXS6
geneSymbols[1966] <- "AKAP2"#Q9Y2D5 

gene_matrix <- cbind(names(geneSymbols), geneSymbols)
colnames(gene_matrix) = c("y","x")

Effecttrain_main$geneName = geneSymbols
Interaction_main$geneName = geneSymbols
mainEffect_GROUP$geneName = geneSymbols

#Volcano plot pre post
  #Lean
Effecttrain_Lean$sig <- "NO"
Effecttrain_Lean$sig[Effecttrain_Lean$Xiao <= 0.05 & Effecttrain_Lean$logFC >= 0] <- "+" 
Effecttrain_Lean$sig[Effecttrain_Lean$Xiao <= 0.05 & Effecttrain_Lean$logFC <= 0] <- "-"
Effecttrain_Lean$newID <- NA
Effecttrain_Lean$newID[Effecttrain_Lean$Xiao <= 0.05] <- geneSymbols[Effecttrain_Lean$Xiao <= 0.05]  

LeanPlot <- ggplot(Effecttrain_Lean, aes(logFC, -log10(P.Value), color = sig, label = newID)) + 
  geom_point(aes(alpha = 0.99, size = 2)) + 
  geom_text_repel(show.legend = F) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 22),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 25),
        legend.position = "none") +
  labs(x = "Log2 FC (Lean Post - Lean Pre)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual("Significant",
                     values = c("#440154FF", "#55C667FF", "gray"), 
                     breaks = c("-","+"),
                     labels = c("Down", "Up")) + 
  xlim(c(-3,3)) 
LeanPlot
  #Obese
Effecttrain_Obese$sig <- "NO"
Effecttrain_Obese$sig[Effecttrain_Obese$Xiao <= 0.05 & Effecttrain_Obese$logFC >= 0] <- "+" 
Effecttrain_Obese$sig[Effecttrain_Obese$Xiao <= 0.05 & Effecttrain_Obese$logFC <= 0] <- "-" 
Effecttrain_Obese$newID <- NA
Effecttrain_Obese$newID[Effecttrain_Obese$Xiao <= 0.05] <- geneSymbols[Effecttrain_Obese$Xiao <= 0.05]  

ObesePlot <- ggplot(Effecttrain_Obese, aes(logFC, -log10(P.Value), color = sig, label = newID)) + 
  geom_point(aes(alpha = 0.99, size = 2)) + 
  geom_text_repel(show.legend = F) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 22),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 25),
        legend.position = "none")+
  labs(x = "Log2 FC (Obese Post - Obese Pre)") + 
  scale_color_manual("Significant",
                     values = c("#440154FF", "#55C667FF", "gray"), 
                     breaks = c("-","+"),
                     labels = c("Down", "Up")) + 
  xlim(c(-3,3))
ObesePlot
  #T2D
Effecttrain_T2D$sig <- "NO"
Effecttrain_T2D$sig[Effecttrain_T2D$Xiao <= 0.05 & Effecttrain_T2D$logFC >= 0] <- "+" 
Effecttrain_T2D$sig[Effecttrain_T2D$Xiao <= 0.05 & Effecttrain_T2D$logFC <= 0] <- "-" 
Effecttrain_T2D$newID <- NA
Effecttrain_T2D$newID[Effecttrain_T2D$Xiao <= 0.05] <- geneSymbols[Effecttrain_T2D$Xiao <= 0.05]  

T2DPlot <- ggplot(Effecttrain_T2D, aes(logFC, -log10(P.Value), color = sig, label = newID)) + 
  geom_point(aes(alpha = 0.99, size = 2)) + 
  geom_text_repel(show.legend = F) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 22),
        axis.title = element_text(size = 25),
        legend.position = "none") +
  labs(y = "-log10(P.Value)", 
       x = "Log2 FC (T2D Post - T2D Pre)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual("Significant",
                     values = c("#440154FF", "#55C667FF", "gray"), 
                     breaks = c("-","+"),
                     labels = c("Down", "Up")) + 
  xlim(c(-3,3))
T2DPlot
ggarrange(T2DPlot, ObesePlot, LeanPlot, nrow = 1)
#Interaction
  #Lean VS Obese
Interaction_LvsO$sig <- "NO"
Interaction_LvsO$sig[Interaction_LvsO$Xiao <= 0.05 & Interaction_LvsO$logFC >= 0] <- "+" 
Interaction_LvsO$sig[Interaction_LvsO$Xiao <= 0.05 & Interaction_LvsO$logFC <= 0] <- "-"
Interaction_LvsO$newID <- NA
Interaction_LvsO$newID[Interaction_LvsO$Xiao <= 0.05] <- geneSymbols[Interaction_LvsO$Xiao <= 0.05]  

LvsOPlot_I <- ggplot(Interaction_LvsO, aes(logFC, -log10(P.Value), color = sig, label = newID)) + 
  geom_point(aes(alpha = 0.99, size = 2)) + 
  geom_text_repel(show.legend = F) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 22),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 25),
        legend.position = "none") +
  labs(x = "Log2 FC (Obese - Lean)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual("Significant",
                     values = c("#440154FF", "#55C667FF", "gray"), 
                     breaks = c("-","+"),
                     labels = c("Down", "Up")) + 
  xlim(c(-3,3))
LvsOPlot_I

  #Obese vs T2D
Interaction_OvsT$sig <- "NO"
Interaction_OvsT$sig[Interaction_OvsT$Xiao <= 0.05 & Interaction_OvsT$logFC >= 0] <- "+" 
Interaction_OvsT$sig[Interaction_OvsT$Xiao <= 0.05 & Interaction_OvsT$logFC <= 0] <- "-" 
Interaction_OvsT$newID <- NA
Interaction_OvsT$newID[Interaction_OvsT$Xiao <= 0.05] <- geneSymbols[Interaction_OvsT$Xiao <= 0.05]  

OvsTPlot_I <- ggplot(Interaction_OvsT, aes(logFC, -log10(P.Value), color = sig, label = newID)) + 
  geom_point(aes(alpha = 0.99, size = 2)) + 
  geom_text_repel(show.legend = F) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 22),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 25),
        legend.position = "none") +
  labs(x = "Log2 FC (T2D - Obese)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual("Significant",
                     values = c("#440154FF", "#55C667FF", "gray"), 
                     breaks = c("-","+"),
                     labels = c("Down", "Up")) + 
  xlim(c(-4,4)) 
OvsTPlot_I

  #Lean VS T2D
Interaction_LvsT$sig <- "NO"
Interaction_LvsT$sig[Interaction_LvsT$Xiao <= 0.05 & Interaction_LvsT$logFC >= 0] <- "+" 
Interaction_LvsT$sig[Interaction_LvsT$Xiao <= 0.05 & Interaction_LvsT$logFC <= 0] <- "-" 
Interaction_LvsT$newID <- NA
Interaction_LvsT$newID[Interaction_LvsT$Xiao <= 0.05] <- geneSymbols[Interaction_LvsT$Xiao <= 0.05]  

LvsTPlot_I <- ggplot(Interaction_LvsT, aes(logFC, -log10(P.Value), color = sig, label = newID)) + 
  geom_point(aes(alpha = 0.99, size = 2)) + 
  geom_text_repel(show.legend = F) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 22),
        axis.title = element_text(size = 25),
        legend.position = "none") +
  labs(y = "-log10(P.Value)", 
       x = "Log2 FC (T2D - Lean)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual("Significant",
                     values = c("#440154FF", "#55C667FF", "gray"), 
                     breaks = c("-","+"),
                     labels = c("Down", "Up")) + 
  xlim(c(-3,3)) 
LvsTPlot_I  

#Group
#Lean VS Obese
OBESE_vs_LEAN$sig <- "NO"
OBESE_vs_LEAN$sig[OBESE_vs_LEAN$Xiao <= 0.05 & OBESE_vs_LEAN$logFC >= 0] <- "+" 
OBESE_vs_LEAN$sig[OBESE_vs_LEAN$Xiao <= 0.05 & OBESE_vs_LEAN$logFC <= 0] <- "-"
OBESE_vs_LEAN$newID <- NA
OBESE_vs_LEAN$newID[OBESE_vs_LEAN$Xiao <= 0.05] <- geneSymbols[OBESE_vs_LEAN$Xiao <= 0.05]  

OvsLPlot_G <- ggplot(OBESE_vs_LEAN, aes(logFC, -log10(P.Value), color = sig, label = newID)) + 
  geom_point(aes(alpha = 0.99, size = 2)) + 
  geom_text_repel(show.legend = F) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 22),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 25),
        legend.position = "none") +
  labs(x = "Log2 FC (Obese - Lean)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual("Significant",
                     values = c("#440154FF", "#55C667FF", "gray"), 
                     breaks = c("-","+"),
                     labels = c("Down", "Up")) + 
  xlim(c(-4,4)) 

OvsLPlot_G

#Obese vs T2D
T2D_vs_OBESE$sig <- "NO"
T2D_vs_OBESE$sig[T2D_vs_OBESE$Xiao <= 0.05 & T2D_vs_OBESE$logFC >= 0] <- "+" 
T2D_vs_OBESE$sig[T2D_vs_OBESE$Xiao <= 0.05 & T2D_vs_OBESE$logFC <= 0] <- "-" 
T2D_vs_OBESE$newID <- NA
T2D_vs_OBESE$newID[T2D_vs_OBESE$Xiao <= 0.05] <- geneSymbols[T2D_vs_OBESE$Xiao <= 0.05]  

TvsOPlot_G <- ggplot(T2D_vs_OBESE, aes(logFC, -log10(P.Value), color = sig, label = newID)) + 
  geom_point(aes(alpha = 0.99, size = 2)) + 
  geom_text_repel(show.legend = F) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 22),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 25),
        legend.position = "none") +
  labs(x = "Log2 FC (T2D - Obese)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual("Significant",
                     values = c("#440154FF", "#55C667FF", "gray"), 
                     breaks = c("-","+"),
                     labels = c("Down", "Up")) + 
  xlim(c(-3,3)) 
TvsOPlot_G

#Lean VS T2D
T2D_vs_LEAN$sig <- "NO"
T2D_vs_LEAN$sig[T2D_vs_LEAN$Xiao <= 0.05 & T2D_vs_LEAN$logFC >= 0] <- "+" 
T2D_vs_LEAN$sig[T2D_vs_LEAN$Xiao <= 0.05 & T2D_vs_LEAN$logFC <= 0] <- "-" 
T2D_vs_LEAN$newID <- NA
T2D_vs_LEAN$newID[T2D_vs_LEAN$Xiao <= 0.05] <- geneSymbols[T2D_vs_LEAN$Xiao <= 0.05]  

TvsLPlot_G <- ggplot(T2D_vs_LEAN, aes(logFC, -log10(P.Value), color = sig, label = newID)) + 
  geom_point(aes(alpha = 0.99, size = 2)) + 
  geom_text_repel(show.legend = F) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 22),
        axis.title = element_text(size = 25),
        axis.ticks = element_line(),
        axis.ticks.length = unit(0.25,'cm'),
        legend.position = "none") +
  labs(y = "-Log10(P.Value)", 
       x = "Log2 FC (T2D - Lean)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual("Significant",
                     values = c("#440154FF", "#55C667FF", "gray"), 
                     breaks = c("-","+"),
                     labels = c("Down", "Up")) + 
  xlim(c(-3,3)) 
TvsLPlot_G
ggarrange(TvsLPlot_G, OvsLPlot_G,TvsOPlot_G, nrow = 1)
#Venn sig Train
eulerr_options(pointsize=30)

Effecttrain_Lean$newID <- geneSymbols

Effecttrain_Obese$newID<- geneSymbols

Effecttrain_T2D$newID <- geneSymbols

VennIndiv_Training <- ggvenn::ggvenn(list(Lean = Effecttrain_Lean$newID[Effecttrain_Lean$Xiao <= 0.05], Obese = Effecttrain_Obese$newID[Effecttrain_Obese$Xiao <= 0.05], T2D = Effecttrain_T2D$newID[Effecttrain_T2D$Xiao <= 0.05]))

#Proportional Venn
VennT <- euler(c(Lean = 11, Obese = 11, T2D = 12,
                 "Lean&T2D" = 2, "Lean&Obese" = 0, "T2D&Obese" = 1))
VennIndiv_Training_prop <- plot(VennT, quantities = T, 
                                fills = color_plots, 
                                alpha = 0.7)

  
  #Venn sig Inter

Interaction_LvsO$newID <- geneSymbols 

Interaction_LvsT$newID <- geneSymbols

Interaction_OvsT$newID <- geneSymbols  

VennIndiv_Interaction <- ggvenn::ggvenn(list(Lean_vs_Obese = Interaction_LvsO$newID[Interaction_LvsO$Xiao <= 0.05], Lean_vs_T2D = Interaction_LvsT$newID[Interaction_LvsT$Xiao <= 0.05], Obese_vs_T2D = Interaction_OvsT$newID[Interaction_OvsT$Xiao <= 0.05]))
  
#Proportional Venn
VennI <- euler(c("Lean vs Obese" = 11, "Obese vs T2D" = 25, "Lean vs T2D" = 34,
                 "Lean vs Obese&Obese vs T2D" = 4, "Lean vs Obese&Lean vs T2D" = 3, "Lean vs T2D&Obese vs T2D" = 1,
                 "Lean vs Obese&Obese vs T2D&Lean vs T2D" = 1))
VennIndiv_Interaction_prop <- plot(VennI, quantities = T, 
                                   fills = color_plots, 
                                   alpha = 0.7)


  #venn Diagram Groups
OBESE_vs_LEAN$newID <- geneSymbols

T2D_vs_LEAN$newID <- geneSymbols

T2D_vs_OBESE$newID <- geneSymbols

VennIndiv_Groups <- ggvenn::ggvenn(list(OBESE_vs_LEAN = OBESE_vs_LEAN$newID[OBESE_vs_LEAN$Xiao <= 0.05], T2D_vs_LEAN = T2D_vs_LEAN$newID[T2D_vs_LEAN$Xiao <= 0.05], T2D_vs_OBESE = T2D_vs_OBESE$newID[T2D_vs_OBESE$Xiao <= 0.05]))

#Proportional Venn
VennG <- euler(c("Obese vs Lean" = 10, "T2D vs Obese" = 13, "T2D vs Lean" = 38,
                 "Obese vs Lean&T2D vs Obese" = 2, "Obese vs Lean&T2D vs Lean" = 4, "T2D vs Lean&T2D vs Obese" = 11))
VennIndiv_Group_prop <- plot(VennG, quantities = T, 
                             fills = color_plots, 
                             alpha = 0.7)

#Add all gene_names main
Effecttrain_main$newID <- geneSymbols
Interaction_main$newID <- geneSymbols
mainEffect_GROUP$newID <- geneSymbols

#Arrange
ggarrange(VennIndiv_Training, VennIndiv_Interaction, VennIndiv_Groups)
ggarrange(plotlist = list(VennIndiv_Training_prop, ggarrange(VennIndiv_Interaction_prop,VennIndiv_Group_prop)), nrow = 2)

ggarrange(LeanPlot, ObesePlot, T2DPlot, VennIndiv_Training_prop)
ggarrange(LvsOPlot_I, OvsTPlot_I, LvsTPlot_I,VennIndiv_Interaction_prop)
ggarrange(OvsLPlot_G, TvsOPlot_G, TvsLPlot_G, VennIndiv_Group_prop)

## Bar plot baseline difference ##

T2D_vs_LEAN_sig <- T2D_vs_LEAN %>% filter(Xiao <= 0.05) %>% select(logFC,newID)
T2D_vs_OBESE_sig <- T2D_vs_OBESE %>% filter(Xiao <= 0.05) %>% select(logFC,newID)
OBESE_vs_LEAN_sig <- OBESE_vs_LEAN %>% filter(Xiao <= 0.05) %>% select(logFC,newID)

T2D_common <- merge(T2D_vs_LEAN_sig,T2D_vs_OBESE_sig, by = "newID")
Obese_common <- merge(OBESE_vs_LEAN_sig,T2D_vs_OBESE_sig, by = "newID")
Lean_common <- merge(OBESE_vs_LEAN_sig,T2D_vs_LEAN_sig, by = "newID")
common <- rbind(T2D_common,Obese_common,Lean_common)
common$newID <- as.factor(common$newID)

ggplot() +
  geom_col(aes(x = factor(rep(common$newID,2), levels = common$newID), 
               y = c(common$logFC.x,common$logFC.y), 
               fill = as.factor(c(rep(1,11),rep(3,2),rep(3,4),rep(2,11),rep(2,2),rep(1,4)))), 
           position = "dodge") +

  theme_minimal() + 
  scale_fill_manual("Group", 
                    values = color_plots,
                    labels = c("T2D vs Lean",
                               "T2D vs Obese",
                               "Obese vs Lean")) + 
  labs(x = "Gene name", y = "Log2FC") +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        axis.line = element_line(),
        axis.ticks = element_line(),
        axis.ticks.length = unit(0.4,"cm")) 
  #geom_rect(aes(xmin = c("","BMPR2","NDUFV2"), xmax = c("UBXN6","IFT 172", "z"),ymin = -3.2,ymax = 2.5), fill = color_plots, alpha = 0.1)

#### Finding stable protein across all samples ####

#Lowest CV across all samples
Exprs_adipose_medianscaled <- medianScaling(Exprs_adipose_notNotmalized)
Exprs_adipose_noBatch_medianScaled <- removeBatchEffect(Exprs_adipose_medianscaled, batch = cluster, design = design)
proteins_Top100 <- rownames(Exprs_adipose_medianscaled)[order(apply(Exprs_adipose_medianscaled, 1, median, na.rm=TRUE), decreasing = T)][1:100]

Exprs_Top100 <- as.data.frame(cbind(Exprs_adipose_noBatch_medianScaled[proteins_Top100,], mainEffect_GROUP[proteins_Top100,]$P.Value, Effecttrain_main[proteins_Top100,]$P.Value))
Exprs_Top100_filter <- Exprs_Top100 %>% filter(V92 > 0.4) %>% filter(V93 > 0.4)

Exprs_adipose_noLog <- 2**Exprs_Top100_filter[,-c(92,93)]
CV <- apply(Exprs_adipose_noLog,1,FUN = function(x) (sd(x)/mean(x))*100)
CV <- CV[order(CV, decreasing = F)]
head(CV)

#### 1D enrichment ####
# 1D on correlation:
r_to_1D <- cbind(result_GIR$Protein,result_GIR$correlation)
OneD_Enrichment_GOBP_GIR <- enrichment_1D_parallel(matrix_annotations_GOBP, r_to_1D)
rownames(OneD_Enrichment_GOBP_GIR ) <- 1:nrow(OneD_Enrichment_GOBP_GIR)
OneD_Enrichment_GOBP_GIR$ID = "GOBP"
OneD_Enrichment_GOCC_GIR  <- enrichment_1D_parallel(matrix_annotations_GOCC, r_to_1D)
rownames(OneD_Enrichment_GOCC_GIR ) <- 1:nrow(OneD_Enrichment_GOCC_GIR)
OneD_Enrichment_GOCC_GIR$ID = "GOCC"
OneD_Enrichment_GOMF_GIR  <- enrichment_1D_parallel(matrix_annotations_GOMF, r_to_1D)
rownames(OneD_Enrichment_GOMF_GIR ) <- 1:nrow(OneD_Enrichment_GOMF_GIR)
OneD_Enrichment_GOMF_GIR$ID = "GOMF"

OneD_Enrichment_GIR <- rbind(OneD_Enrichment_GOMF_GIR,OneD_Enrichment_GOCC_GIR,OneD_Enrichment_GOBP_GIR)

ggplot(OneD_Enrichment_GIR,aes(s,-log10(as.numeric(pVal)), color = ID, label = Annotation)) +
  geom_vline(xintercept = 0, color = "black") +
  geom_hline(yintercept = 3, color = "black") +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal() +
  scale_color_manual("GO Term", values = color_plots) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18)) +
  scale_x_continuous(breaks = seq(from = -0.9, to = 0.9, by = 0.3), limits = c(-0.9,0.9), labels = c(-0.9,-0.6,-0.3,0,0.3,0.6,0.9)) + 
  labs(x = "Enrichement score", y = "-Log10 pVal") +
  ylim(c(3,8)) + 
  geom_text_repel(show.legend = F)

#### 2D enrichment ####
#Interaction
LogLvsT <- cbind(rownames(Interaction_LvsT), as.numeric(Interaction_LvsT[,1]))
LogOvsT <- cbind(rownames(Interaction_OvsT), as.numeric(Interaction_OvsT[,1]))

Entrichment_2D_Inter_GOCC <- Enrichment_2D_parallel(matrix_annotations_GOCC, LogLvsT, LogOvsT, 0.1)
Entrichment_2D_Inter_GOCC$GO <- "GOCC"
Entrichment_2D_Inter_GOBP <- Enrichment_2D_parallel(matrix_annotations_GOBP, LogLvsT, LogOvsT, 0.1)
Entrichment_2D_Inter_GOBP$GO <- "GOBP"
Entrichment_2D_Inter_GOMF <- Enrichment_2D_parallel(matrix_annotations_GOMF, LogLvsT, LogOvsT, 0.1)
Entrichment_2D_Inter_GOMF$GO <- "GOMF"

Enrichment_2d_Inter <- rbind(Entrichment_2D_Inter_GOCC,Entrichment_2D_Inter_GOBP,Entrichment_2D_Inter_GOMF)
ggplot(Enrichment_2d_Inter, aes(x,y,label = annotation, color = GO)) + 
  geom_vline(xintercept = 0, color = "red") +
  geom_hline(yintercept = 0, color = "red") +
  geom_point(size = 2) +
  geom_text_repel(aes(size = 4), max.overlaps = 9, show.legend = F) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(), 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 30),
        legend.title = element_text(size = 20),
        legend.text =  element_text(size = 18)) +
  scale_color_manual(values =c("#440154FF", "#55C667FF") , "GO term") +
  xlim(c(-0.6,0.6)) +
  ylim(c(-0.6,0.6)) +
  labs(x = "Interaction T2D vs Lean", y = "Interaction T2D vs Obese")

a <- as.data.frame(LogLvsT[which(LogLvsT[,1] %in% matrix_annotations_GOCC$Protein[which(matrix_annotations_GOCC$Annotation == "extracellular region")]),])
b <- as.data.frame(LogOvsT[which(LogOvsT[,1] %in% matrix_annotations_GOCC$Protein[which(matrix_annotations_GOCC$Annotation == "extracellular region")]),])
c <- inner_join(a,b,by = "Protein") 
c$LogFC.x <- as.numeric(c$LogFC.x)
c$LogFC.y <- as.numeric(c$LogFC.y)
str(c)

#Interaction Obesity
LogLvsT <- cbind(rownames(Interaction_LvsT), as.numeric(Interaction_LvsT[,1]))
LogLvsO <- cbind(rownames(Interaction_LvsO), as.numeric(Interaction_LvsO[,1]))

Entrichment_2D_Inter_GOCC_Obesity <- Enrichment_2D_parallel(matrix_annotations_GOCC, LogLvsO, LogLvsT, 0.1)
Entrichment_2D_Inter_GOCC_Obesity$GO <- "GOCC"
Entrichment_2D_Inter_GOBP_Obesity <- Enrichment_2D_parallel(matrix_annotations_GOBP, LogLvsO, LogLvsT, 0.1)
Entrichment_2D_Inter_GOBP_Obesity$GO <- "GOBP"
Entrichment_2D_Inter_GOMF_Obesity <- Enrichment_2D_parallel(matrix_annotations_GOMF, LogLvsO, LogLvsT, 0.1)
Entrichment_2D_Inter_GOMF_Obesity$GO <- "GOMF"

Enrichment_2d_Inter_Obesity <- rbind(Entrichment_2D_Inter_GOCC_Obesity,Entrichment_2D_Inter_GOBP_Obesity,Entrichment_2D_Inter_GOMF_Obesity)
ggplot(Enrichment_2d_Inter, aes(x,y,label = annotation, color = GO)) + 
  geom_vline(xintercept = 0, color = "red") +
  geom_hline(yintercept = 0, color = "red") +
  geom_point(size = 2) +
  geom_text_repel(aes(size = 4), max.overlaps = 9, show.legend = F) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(), 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 30),
        legend.title = element_text(size = 20),
        legend.text =  element_text(size = 18)) +
  scale_color_manual(values =c("#440154FF", "#55C667FF") , "GO term") +
  xlim(c(-0.6,0.6)) +
  ylim(c(-0.6,0.6)) +
  labs(x = "Interaction Obese vs Lean", y = "Interaction T2D vs Lean")

#Group
LogTvsL <- cbind(rownames(T2D_vs_LEAN), T2D_vs_LEAN$logFC)
LogTvsO <- cbind(rownames(T2D_vs_OBESE), T2D_vs_OBESE$logFC)

Entrichment_2D_Group_GOCC <- Enrichment_2D_parallel(matrix_annotations_GOCC, LogTvsL, LogTvsO, 0.1)
Entrichment_2D_Group_GOCC$GO <- "GOCC"
Entrichment_2D_Group_GOBP <- Enrichment_2D_parallel(matrix_annotations_GOBP, LogTvsL, LogTvsO, 0.1)
Entrichment_2D_Group_GOBP$GO <- "GOBP"
Entrichment_2D_Group_GOMF <- Enrichment_2D_parallel(matrix_annotations_GOMF, LogTvsL, LogTvsO, 0.1)
Entrichment_2D_Group_GOMF$GO <- "GOMF"

Enrichment_2d_Group <- rbind(Entrichment_2D_Group_GOCC,Entrichment_2D_Group_GOBP,Entrichment_2D_Group_GOMF)
ggplot(Enrichment_2d_Group, aes(x,y,label = annotation, color = GO)) +
  geom_vline(xintercept = 0, color = "red") +
  geom_hline(yintercept = 0, color = "red") +
  geom_point(size = 2) +
  geom_text_repel(aes(size = 4), max.overlaps = 9, show.legend = F) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(), 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 30),
        legend.title = element_text(size = 20),
        legend.text =  element_text(size = 18)) +
  scale_color_manual(values =c("#440154FF", "#55C667FF") , "GO term") +
  xlim(c(-0.8,0.8)) +
  ylim(c(-0.8,0.8)) +
  labs(x = "T2D vs Lean", y = "T2D vs Obese")

#Effect obesity
LogTvsL <- cbind(rownames(T2D_vs_LEAN), T2D_vs_LEAN$logFC)
LogOvsL <- cbind(rownames(OBESE_vs_LEAN), OBESE_vs_LEAN$logFC)

Entrichment_2D_Group_GOCC_Obese <- Enrichment_2D_parallel(matrix_annotations_GOCC, LogTvsL, LogOvsL, 0.1)
Entrichment_2D_Group_GOCC_Obese$GO <- "GOCC"
Entrichment_2D_Group_GOBP_Obese <- Enrichment_2D_parallel(matrix_annotations_GOBP, LogTvsL, LogOvsL, 0.1)
Entrichment_2D_Group_GOBP_Obese$GO <- "GOBP"
Entrichment_2D_Group_GOMF_Obese <- Enrichment_2D_parallel(matrix_annotations_GOMF, LogTvsL, LogOvsL, 0.1)
Entrichment_2D_Group_GOMF_Obese$GO <- "GOMF"

Enrichment_2d_Group_Obese <- rbind(Entrichment_2D_Group_GOCC_Obese,Entrichment_2D_Group_GOBP_Obese,Entrichment_2D_Group_GOMF_Obese)
ggplot(Enrichment_2d_Group_Obese, aes(x,y,label = annotation, color = GO)) +
  geom_vline(xintercept = 0, color = "red") +
  geom_hline(yintercept = 0, color = "red") +
  geom_point(size = 2) +
  geom_text_repel(aes(size = 4), max.overlaps = 10, show.legend = F) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(), 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 30),
        legend.title = element_text(size = 20),
        legend.text =  element_text(size = 18)) +
  scale_color_manual(values =c("#440154FF", "#55C667FF") , "GO term") +
  xlim(c(-0.8,0.8)) +
  ylim(c(-0.8,0.8)) +
  labs(x = "T2D vs Lean", y = "Obese vs Lean")
