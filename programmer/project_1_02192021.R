#1 Installed bioconductor and packages

#2 Load packages
biocinstallRepos()
source("https://bioconductor.org/biocLite.R")
biocLite(c("affy", "affyPLM", "sva", "AnnotationDbi", "hgu133plus2.db"))
library("BiocManager")
library("affy")
library("sva")
library("affyPLM")
library("AnnotationDbi")
library("hgu133plus2.db")
library(tidyverse)

#3 Read and normalize CEL files 
ReadAffyData <- ReadAffy(celfile.path="/projectnb/bf528/users/wheeler/project_1/sample/files")                        
RMA_Data <- rma(ReadAffyData)  
#write.table(RMA_Data, file="/projectnb/bf528/users/wheeler/project_1/analysis/analyst/RMA_Data.csv")

#4 Calculate Relative Log Expression (RLE) and Normalized Un-scaled Standard Error (NUSE) scores for microarray samples
RLE_NUSE <- fitPLM(ReadAffyData, normalize = TRUE, background = TRUE)

RLE_Median <- data.frame(median=RLE(RLE_NUSE, type="stats")[1,]) # will compute the summary statistics of RLE (including median) by array
#RLE(RLE_NUSE, type=values) #returns all RLE expression values
ggplot(data = RLE_Median) + geom_histogram(mapping = aes(x=median))

NUSE_Median <- data.frame(median=NUSE(RLE_NUSE, type="stats")[1,]) #will compute the summary stats (including median) of NUSE by array
#NUSE(RLE_NUSE, type=values) #will return all NUSE values
ggplot(data=NUSE_Median) + geom_histogram(mapping=aes(x=median))

#5 Correct for batch effects 
#reading the file with batch correction variable and features variable:
BatchData <- read_csv("/project/bf528/project_1/doc/proj_metadata.csv")
#normalizationcombatbatch <- as.character(BatchData[, "normalizationcombatbatch"])  #column 36 in BatchData #batch effects
#normalizationcombatmod <- BatchData[, "normalizationcombatmod"]  #column 37 in BatchData #features of interest
batch=BatchData$normalizationcombatbatch

mod <- model.matrix(~as.factor(normalizationcombatmod), data=BatchData)
#batch <- model.matrix(~as.factor(normalizationcombatbatch), data=BatchData)
#performing batch corrections while maintaining features of interest-
expression <- exprs(RMA_Data)
batch_corr_data <-ComBat(dat=expression, batch=batch, mod=mod)

#write normalized expression data to a csv file-
write.table(batch_corr_data, file="/projectnb/bf528/users/wheeler/project_1/analysis/analyst/batch_corr_data_final.csv")

#6 Perform PCA on normalized data
transposed_data <-t(batch_corr_data) #transposing the data so the gene is the column and the sample is the row
data_for_PCA <- t(scale(transposed_data, center = TRUE, scale = TRUE)) #need scale the data by gene and then re-transpose prior to PCA
prcomp(data_for_PCA, center=FALSE, scale=FALSE) #already scaled and centered the genes (why had to transpose-otherwise this fxn would scale by sample)
#can view values for ea of the PC's using $rotation attribute in the prcomp object

#7 Plot PC1 vs PC2 and examine for outliers
#examine % variability explained by ea PC by using the $importance attribute
