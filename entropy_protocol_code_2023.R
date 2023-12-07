#install packages 
install.packages("ggplot2")
install.packages("reshape2")
install.packages("Matrix")
install.packages("grid")
install.packages("stringr")
install.packages("dplyr")
install.packages("devtools")
devtools::install_github("pcahan1/singleCellNet")

#load packages 
library(ggplot2)
library(reshape2)
library(Matrix)
library(grid)
library(stringr)
library(dplyr)
library(singleCellNet)

#set working directory 
setwd("/Users/Desktop/Entropy/MyWorkSpace")

#load R workspace file 
load("clean_nodatasets_060720.RData")
##optional
##load R workspace file with in vivo reference data to reproduce figures in the protocol paper 
#load("clean_060720.RData")

#load functions 
source("entropy_functions.R")

#prepare counts table and confirm format, cell barcodes as colnames and gene symbols as rownames
counts_table = kannan_ref_data
counts_table[1:4,1:4]

##optional
##convert ENSEMBL ID to gene symbols
##for mouse cells 
#rename_genes(counts_table, species = "mouse")
##for human cells 
#rename_genes(counts_table, species = "human")


#prepare pheno table 
pheno_table = combined_datasets[combined_datasets$data == "kannan_ref_data",]
pheno_table$timepoint[1:10]


#run data_qc() to calculate qc metrics and entropy scores for each cell 

qc_output = data_qc(dataset = "counts_table", #REQUIRED, dataset name must be in quotes 
                    study = " Brilliant Grad Student et al.", #REQUIRED
                    timepoint_list = pheno_table$timepoint, #REQUIRED
                    scn_calc = TRUE, #OPTIONAL, defaulted to TRUE, run SingleCellNet to classify cell types 
                    species = "mouse", #OPTIONAL, defaulted to "mouse", change to "human" for human datasets 
                    sample_type = "in vivo", #OPTIONAL
                    isolation = "LP-FACS", #OPTIONAL
                    sequencing = "mcSCRB-seq", #OPTIONAL
                    mapping = "zUMIs", #OPTIONAL
                    datatype = "UMIs", #OPTIONAL
                    doi = "doi:12345", #OPTIONAL
                    other_meta = NA) #OPTIONAL, other metadata of interest, will be compiled to qc_output

##Optional
##calculate entropy score only, no qc steps
#master_entropy(counts_table)


#Figure 1, entropy score over different timepoints violin plot 
#pick out cardiomyocytes that pass QC, good_cell == "TRUE"
qc_output_goodcell = qc_output[qc_output$good_cell == "TRUE",]
ggplot(qc_output_goodcell , aes(x = timepoint, y = entropy, fill = timepoint)) + geom_jitter(size = 0.5) + geom_violin(scale = "width") + geom_boxplot(position = position_dodge(width = 1)) + xlab("Timepoint") + ylab(expression('Shannon Entropy'~italic(S))) + theme_linedraw() + theme(legend.position = "none")

#Figure 2, gene expression trends over the negative value of entropy score, relative gene expression levels are normalize by counts per million
counts_table_goodcell = counts_table[,rownames(qc_output_goodcell)]
head(counts_table_goodcell)[1:5,1:5]
#calculate depth for each cell
lib_size = colSums(counts_table_goodcell)
lib_size[1:5]
#normalize by counts per million 
counts_table_goodcell_norm = sweep(counts_table_goodcell, 2, lib_size, FUN = '/')
counts_table_goodcell_norm * 10e+06
head(counts_table_goodcell_norm)[1:5,1:5]
dim(counts_table_goodcell_norm)
counts_table_goodcell_norm = t(counts_table_goodcell_norm)
dim(counts_table_goodcell_norm)
head(counts_table_goodcell_norm)[1:5,1:5]
#take the begative value of entropy score, so the more mature state is aligned with the positive direction of the x axis
plot_entropy = -1*qc_output_goodcell$entropy
plot_data = as.data.frame(as.matrix(counts_table_goodcell_norm))
plot_data$plot_entropy = plot_entropy
plot_data$timepoint = qc_output_goodcell$timepoint
#Myh6 as an example
ggplot(plot_data, aes(x = plot_entropy, y = Myh6, color = timepoint)) + geom_jitter(size = 1) + geom_smooth(method = lm, formula = y ~ ns(x, df=3), se = FALSE, color = "black", size = 0.5) + theme_linedraw() + ylim(0,NA) + theme(legend.position="none") + labs(title = "Myh6") + theme(plot.title = element_text(hjust = 0.5)) + xlab("-1*Entropy") +ylab("Relative Expression")

