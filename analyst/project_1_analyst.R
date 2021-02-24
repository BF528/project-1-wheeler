
##################################################################################
#
#  PROGRAMMER TASKS ####
#
##################################################################################

# part 3.1 - packages installation

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()
# 
# 
# BiocManager::install("affy")
# BiocManager::install("affyPLM")
# BiocManager::install("sva")
# BiocManager::install("AnnotationDbi")
# BiocManager::install("hgu133plus2.db")

# part 3.2 - load the R packages

library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)
library(tidyverse)
library(ggplot2)
library(parallel)

# part 3.3 - normalized rma

# path to the cell files
cel_path = "/projectnb/bf528/users/wheeler/project_1/sample/files"

# read in the cell files
cel_data <- ReadAffy(celfile.path=cel_path, verbose=TRUE)

# normalized the cell files
rma_data <- rma(cel_data)

# extract the expression matrix only
rma_normalized <- exprs(rma_data)

# output a csv file of the rma normalized
write.table(rma_normalized, file="/projectnb/bf528/users/wheeler/project_1/analysis/analyst/rma_normalized.csv")

# part 3.4 - compute RLE and NUSE scores
fit <- fitPLM(cel_data, normalize=TRUE, background=TRUE)

# get the RLE scores summary
RLE_score <- RLE(fit, type="stats")

# get the median of RLE scores 
RLE_median <- data.frame(median=RLE_score["median",])

# write out a copy of the median of RLE scores
write.table(RLE_median, file="/projectnb/bf528/users/wheeler/project_1/analysis/analyst/RLE_median.csv")

# histogram of RLE median
RLE_median %>% 
  ggplot(aes(median)) +
  geom_histogram(bins=30, fill='white', color="black") +
  ggtitle("Median RLE Scores")

# get the NUSE scores summary
NUSE_score <- NUSE(fit, type="stats")

# get the median of NUSE scores 
NUSE_median <- data.frame(median=NUSE_score["median",])

# write out a copy of the median of NUSE scores
write.table(NUSE_median, file="/projectnb/bf528/users/wheeler/project_1/analysis/analyst/NUSE_median.csv")

# histogram of NUSE median
NUSE_median %>% 
  ggplot(aes(median)) +
  geom_histogram(bins=30, fill='white', color="black") +
    ggtitle("Median NUSE Scores")

## part 3.5 - correct for batch effects ###

## read in the annotation file ###
annotation <- read.csv("/project/bf528/project_1/doc/proj_metadata.csv", header=T)

## compute combat to correct for batch effects ###
edata = rma_normalized
batch = annotation$normalizationcombatbatch
mod = model.matrix(~as.factor(normalizationcombatmod), data=annotation)
combat_edata = ComBat(dat=edata, batch=batch, mod=mod)

# write out a copy of the combat expression set ###
write.table(combat_edata, "/projectnb/bf528/users/wheeler/project_1/analysis/analyst/combat_edata.csv")

## part 3.6 -- pca ##
rma_normalized <- read.table("/projectnb/bf528/users/wheeler/project_1/analysis/analyst/rma_normalized.csv", quote="\"", comment.char="")

# scale the normalized expression set ##
scale_rma_normalized <- scale(t(rma_normalized), center=TRUE, scale=TRUE) %>% t()

# perform pca on the scale data ###
pca <- prcomp(scale_rma_normalized, center=F, scale=F) %>% summary()

# obtain the values for each of the principal components ###
pca_values <- pca$rotation

# write out a copy of the principal components ###
write.table(pca_values, file="/projectnb/bf528/users/wheeler/project_1/analysis/analyst/pca_rotations.csv")

# the percent variability explained by each principal component 
pca_variance_explained <- pca$importance

# write out a copy of the variability ###
write.table(pca$importance, file="/projectnb/bf528/users/wheeler/project_1/analysis/analyst/pca_importance.csv")

# a plot of PC1 vs PC2 with the percent variability attributed to these principal components shown on the x and y axes labels
pca_values %>% as.data.frame() %>% 
  ggplot() + ggtitle("PC1 vs PC2") +
  geom_point(aes(x=PC1, y=PC2))

##################################################################################
#
#  ANALYSST TASKS ####
#
##################################################################################

## Read in the combat expression dataset
combat_edata <- read.table("/projectnb/bf528/users/wheeler/project_1/analysis/analyst/combat_edata.csv", quote="\"", comment.char="")

# part 4.1 - a function to filter differential expressed genes##
filter_expression <- function(expression_set){
  
  filter_genes <- apply(expression_set, 1, function(x){
    x <- as.numeric(x)
    n <- length(x) # number of expression values
    l <- length(which(x > log2(15))) # number of expression values >= log2(15)
    percent <- l/n  # percent of expression values >= log2(15)
    
    if(percent >= 0.2){
      return(TRUE)
    }else{
      return(FALSE)
    }
  })
  
  # only retain the the gene that has at least 20% of the gene-expression values > log2(15)
  expression_set <- expression_set[filter_genes %in% TRUE,]
  
  return(expression_set)

}

## part 4.2 - a function to calculate the chsq test ##
library(DescTools)

chisq_test <- function(expression_set){
  
  median_variance <- median(apply(expression_set, 1, var))

  chisq_results <- apply(expression_set, 1, function(x){
    x <- as.numeric(x)
    test <- VarTest(x, sigma.squared=median_variance)
    statistic=test$statistic; df=length(x)-1; alpha=0.01;
    chisq_l <- qchisq(p=(alpha/2), df=df)
    chisq_u <- qchisq(p=1-(alpha/2), df=df)
    
    if(statistic < chisq_l | statistic > chisq_u){
      return(TRUE)
    }else{
      return(FALSE)
    }
  })
  
  expression_set <- expression_set[chisq_results %in% TRUE,]
   
  return(expression_set)
  
}

## part 3.3 - a function to filter coefficient of variation  > 0.186 ##
coef_of_variation <- function(expression_set){
  
  CV <- apply(expression_set, 1, function(x){
    x <- as.numeric(x)
    mean_x <- mean(x)
    sd_x <- sd(x)
    cv <- sd_x/mean_x

    if(cv > 0.186){
      return(TRUE)
    }else{
      return(FALSE)
    }
  })
  
  expression_set <- expression_set[CV %in% TRUE,]
  
  return(expression_set)  
  
}

## part 4.4 - Write the filter expression matrix with all the three filter functions ## ##
combat_probesets_all_filtered <- combat_edata %>% 
  filter_expression() %>% 
  chisq_test() %>% 
  coef_of_variation()

write.table(combat_probesets_all_filtered, file="/projectnb/bf528/users/wheeler/project_1/analysis/analyst/combat_probesets_all_filtered.csv")

## part 4.5 - Write the filter expression matrix with 4.2 filter only for biologist ##
combat_probesets_bio_filtered <- combat_edata %>% chisq_test()

write.table(combat_probesets_bio_filtered, file="/projectnb/bf528/users/wheeler/project_1/analysis/analyst/combat_probesets_bio_filtered.csv")

## TASK 5 ##
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms

annotation <- read.csv("/project/bf528/project_1/doc/proj_metadata.csv", header=T)
combat_probesets_all_filtered <- read.table("/projectnb/bf528/users/wheeler/project_1/analysis/analyst/combat_probesets_all_filtered.csv", quote="\"", comment.char="")
combat_probesets_bio_filtered <- read.table("/projectnb/bf528/users/wheeler/project_1/analysis/analyst/combat_probesets_bio_filtered.csv", quote="\"", comment.char="")

## part 5.1 ---

# compute the euclidean distance btw the combat samples ##
combat_distance <- dist(t(combat_probesets_all_filtered), method="euclidean")

# Hierarchical clustering using ward distance
combat_hc <- hclust(combat_distance, method="ward.D")

## part 5.2 - create dendrogram and divide the samples into 2 clusters ###

# Cut tree into 2 groups
combat_sub_grp <- cutree(combat_hc, k=2) 

write.table(combat_sub_grp, file="/projectnb/bf528/users/wheeler/project_1/analysis/analyst/combat_sub_grp.csv")

# Number of members in each cluster
table(combat_sub_grp)

# part 5.3 - heatmap ###
colsidecolor <- lapply(annotation$cit.coloncancermolecularsubtype, function(x){
    if(x == "C3"){
      return("red")
    }else{
      return("blue")
    }
  }) %>% unlist()

## cluster rows (genes) and columns (samples)
hc.row <- hclust(as.dist(1-cor(t(combat_probesets_all_filtered))), method="ward.D2") # genes by correlation
hc.col <- hclust(dist(t(combat_probesets_all_filtered)), method="ward.D2")           # samples by euclidean distance (default)

## color gradient for the expression levels (blue=down-regulated; white=neutral; red=up-regulated)
colGradient <- function(cols, length, cmax=255){
  ramp <- colorRamp(cols)
  rgb( ramp(seq(0,1,length=length)), max=cmax )
}

bwrPalette <- colGradient(c("blue","white","red"), length=13)

# heatmap for combat expression set
heatmap(as.matrix(combat_probesets_all_filtered), Rowv=as.dendrogram(hc.row), Colv=as.dendrogram(hc.col), ColSideColors=colsidecolor, col=bwrPalette, labCol=NA,labRow=NA)

## part 5.4 - function to identify genes differentially expressed between the two clusters using a Welch t-test ###
welch_test <- function(expression_set, sub_grp){
  
  cluster1 <- expression_set[,which(colnames(expression_set) %in% names(sub_grp[sub_grp == 1]))]
  cluster2 <- expression_set[,which(colnames(expression_set) %in% names(sub_grp[sub_grp == 2]))]
  
  results <- data.frame(t=rep(NA, nrow(expression_set)), p=rep(NA, nrow(expression_set)))
  rownames(results) <- rownames(expression_set)
  
  for(r in 1:nrow(results)){
    #r=1;
    x <- as.numeric(cluster1[r,])
    y <- as.numeric(cluster2[r,])
    test <- t.test(x,y)
    results$t[r] = test$statistic
    results$p[r] = test$p.value
  }
  
  return(results)
  
}

## function to calculate log fold change btw the two clusters
log2foldchange <- function(expression_set, sub_grp){
  
  cluster1 <- expression_set[,which(colnames(expression_set) %in% names(sub_grp[sub_grp == 1]))]
  cluster2 <- expression_set[,which(colnames(expression_set) %in% names(sub_grp[sub_grp == 2]))]
  
  results <- data.frame(mean1=rep(NA, nrow(expression_set)), mean2=rep(NA, nrow(expression_set)), logfoldchange=rep(NA, nrow(expression_set)))
  rownames(results) <- rownames(expression_set)
  
  for(r in 1:nrow(results)){
    #r=1;
    mean1 <- mean(as.numeric(cluster1[r,]), na.rm=T)
    mean2 <- mean(as.numeric(cluster2[r,]), na.rm=T)
    results$mean1[r] = mean1
    results$mean2[r] = mean2
    results$logfoldchange[r] = log2(mean1/mean2)
  }
  
  return(results)
  
}

## Compute welch test and the adjusted p-value ###
combat_welch_results <- welch_test(expression_set=combat_probesets_all_filtered, sub_grp=combat_sub_grp) %>% 
  rownames_to_column() %>% 
  mutate(padj=p.adjust(p, method="fdr")) %>% 
  column_to_rownames()

# Export the t-test results computed on the expression matrix  ###
write.table(combat_welch_results, file="/projectnb/bf528/users/wheeler/project_1/analysis/analyst/combat_welch_results.csv")

## number of genes are differentially expressed at adjusted p < 0.05 between the clusters for both lists ##
combat_welch_genes_filtered <- combat_welch_res %>% filter(padj < 0.05)

# Export the t-test results computed on the expression matrix  ###
write.table(combat_welch_genes_filtered, file="/projectnb/bf528/users/wheeler/project_1/analysis/analyst/combat_welch_genes_filtered.csv")

## Compute log2foldchange ###
log2foldchange_results <- log2foldchange(expression_set=combat_probesets_all_filtered, sub_grp=combat_sub_grp)

## filter the adjusted p-value < 0.05 ###
c3_cluster_genes <-  combat_welch_genes_filter %>% mutate(PROBEID = rownames(.)) %>% 
  left_join(log2foldchange_results %>% mutate(PROBEID = rownames(.))) %>% 
  filter(padj < 0.05 & logfoldchange > 0) 

write.table(c3_cluster_genes, file="/projectnb/bf528/users/wheeler/project_1/analysis/analyst/c3_cluster_genes.csv")

c4_cluster_genes <-  combat_welch_genes_filter %>% mutate(PROBEID = rownames(.)) %>% 
  left_join(log2foldchange_results %>% mutate(PROBEID = rownames(.))) %>% 
  filter(padj < 0.05 & logfoldchange < 0) 

write.table(c4_cluster_genes, file="/projectnb/bf528/users/wheeler/project_1/analysis/analyst/c4_cluster_genes.csv")

## part 5.6 ###

## Compute welch test and the adjusted p-value ###
bio_combat_welch_results <- welch_test(expression_set=combat_probesets_bio_filtered, sub_grp=combat_sub_grp) %>% 
  rownames_to_column() %>% 
  mutate(padj=p.adjust(p, method="fdr")) %>% 
  column_to_rownames()

# Export the t-test results computed on the expression matrix on 4.5 for biologist ###
write.table(bio_combat_welch_results, file="/projectnb/bf528/users/wheeler/project_1/analysis/analyst/bio_combat_welch_results.csv")
write.table(bio_combat_welch_results, file="/projectnb/bf528/users/wheeler/project_1/analysis/biologist/differential_expression_results.csv")

## Compute log2foldchange ###
bio_log2foldchange_results <- log2foldchange(expression_set=combat_probesets_bio_filtered, sub_grp=combat_sub_grp)

# Export the t-test results computed on the expression matrix  ###
write.table(bio_log2foldchange_results, file="/projectnb/bf528/users/wheeler/project_1/analysis/analyst/bio_log2foldchange_results.csv")

## number of genes are differentially expressed at adjusted p < 0.05 between the clusters for both lists ##
bio_welch_genes_filtered <- bio_combat_welch_results %>% 
  filter(padj < 0.05)

write.table(bio_welch_genes_filtered, file="/projectnb/bf528/users/wheeler/project_1/analysis/analyst/bio_welch_genes_filtered.csv")

##################################################################################
#
#  BIOLOGIST TASKS ####
#
##################################################################################

## part 6.1 ###
bio_combat_welch_results <- read.table("/projectnb/bf528/users/wheeler/project_1/analysis/analyst/bio_welch_genes_filtered.csv", quote="\"", comment.char="")
bio_log2foldchange_results <- read.table(file="/projectnb/bf528/users/wheeler/project_1/analysis/analyst/bio_log2foldchange_results.csv", quote="\"", comment.char="")

# map the probeset IDs to gene symbols 
gene_symbols <- AnnotationDbi::select(hgu133plus2.db, keys=rownames(bio_combat_welch_results), columns=c("SYMBOL"))

bio_combat_welch_results_2 <- bio_combat_welch_results %>% mutate(PROBEID = rownames(.)) %>% 
  left_join(bio_log2foldchange_results %>% mutate(PROBEID = rownames(.)))

gene_list_no_duplicates <- gene_symbols[!duplicated(gene_symbols$SYMBOL),] %>% 
  left_join(bio_combat_welch_results_2)

write.table(gene_list_no_duplicates, file="/projectnb/bf528/users/wheeler/project_1/analysis/analyst/gene_list_no_duplicates.csv")

## part 6.2 ###

#BiocManager::install("GSEABase")
library(GSEABase)

## part 6.3 - the top 10 up-regulated genes###
bio_combat_welch_up_regulated <- gene_list_no_duplicates %>% 
  filter(logfoldchange > 0) %>% 
  arrange(desc(logfoldchange), padj) %>% 
  slice(1:1000) %>% 
  mutate(genes=SYMBOL, logfoldchange=round(logfoldchange,4), t=round(t,4), p=round(p,4), padj=round(padj,4)) %>% 
  select(genes, logfoldchange, t, p, padj)

write.csv(bio_combat_welch_up_regulated, file="/projectnb/bf528/users/wheeler/project_1/analysis/analyst/top_1000_up_regulated.csv", row.names=F)

## the top 10 down-regulated genes###
bio_combat_welch_down_regulated <- gene_list_no_duplicates %>%  
  filter(logfoldchange < 0) %>% 
  arrange(logfoldchange, padj) %>% 
  slice(1:1000) %>% 
  mutate(genes=SYMBOL, logfoldchange=round(logfoldchange,4), t=round(t,4), p=round(p,4), padj=round(padj,4)) %>% 
  select(genes, logfoldchange, t, p, padj)

write.csv(bio_combat_welch_down_regulated, file="/projectnb/bf528/users/wheeler/project_1/analysis/analyst/top_1000_down_regulated.csv", row.names=F)

## part 6.4 - load gsva ###
library(GSEABase)

##read in the gene sets collection
KEGG <- getGmt(paste0("/projectnb/bf528/users/wheeler/project_1/sample/c2.cp.kegg.v7.2.symbols.gmt")) # 186
GO <- getGmt(paste0("/projectnb/bf528/users/wheeler/project_1/sample//c5.go.v7.2.symbols.gmt")) # 10271
Hallmark <- getGmt(paste0("/projectnb/bf528/users/wheeler/project_1/sample/h.all.v7.2.symbols.gmt")) # 50

## part 6.5 -  a function to do the fisher test ##
Enrichment_Analysis <- function(GeneSetCollection, genelist){
  
  collection_geneset <- names(GeneSetCollection)
  differential_expressed_genes <- genelist %>% filter(padj < 0.05) %>% select(SYMBOL) %>% unlist() %>% as.character()
  not_differential_expressed_genes <- genelist %>% filter(padj >= 0.05) %>% select(SYMBOL) %>% unlist() %>% as.character()
  
  ftest_stat <- seq_along(collection_geneset) %>%
    map_dfr(
      function(g){
        #g=1
        geneset <- GeneSetCollection[[collection_geneset[g]]]@geneIds 
        
        IGS_differential_expressed_genes <- length(which(differential_expressed_genes %in% geneset))
        NIGS_differential_expressed_genes <- length(differential_expressed_genes) - IGS_differential_expressed_genes
        
        IGS_not_differential_expressed_genes <- length(which(not_differential_expressed_genes %in% geneset))
        NIGS_not_differential_expressed_genes <- length(not_differential_expressed_genes) - IGS_not_differential_expressed_genes
        
        table <- matrix(c(IGS_differential_expressed_genes, NIGS_differential_expressed_genes, IGS_not_differential_expressed_genes, NIGS_not_differential_expressed_genes), nrow=2)
        f.test <- fisher.test(table)
        
        results <- data.frame(geneset=collection_geneset[g], statistic=f.test$estimate, p=f.test$p.value)
        
        return(results)
      }
    )
  
  return(ftest_stat)
  
}  

# run the fisher's test and calculate FDR adjusted p-values
KEGG_enrichment <- Enrichment_Analysis(GeneSetCollection=KEGG, genelist=gene_list_no_duplicates) %>% 
  mutate(padj=p.adjust(p, method="fdr")) %>% 
  filter(padj < 0.05) %>% 
  arrange(padj)
GO_enrichment  <- Enrichment_Analysis(GeneSetCollection=GO, genelist=gene_list_no_duplicates) %>% 
  mutate(padj=p.adjust(p, method="fdr")) %>% 
  filter(padj < 0.05) %>% 
  arrange(padj)
Hallmark_enrichment  <- Enrichment_Analysis(GeneSetCollection=Hallmark, genelist=gene_list_no_duplicates) %>% 
  mutate(padj=p.adjust(p, method="fdr")) %>% 
  filter(padj < 0.05) %>% 
  arrange(padj)

# write out a copy to csv
write.csv(KEGG_enrichment, file="/projectnb/bf528/users/wheeler/project_1/analysis/analyst/KEGG_enrichment.csv", row.names=F)
write.csv(GO_enrichment, file="/projectnb/bf528/users/wheeler/project_1/analysis/analyst/GO_enrichment.csv", row.names=F)
write.csv(Hallmark_enrichment, file="/projectnb/bf528/users/wheeler/project_1/analysis/analyst/Hallmark_enrichment.csv", row.names=F)

## read in the enrichment results
KEGG_enrichment <- read.csv(file="/projectnb/bf528/users/wheeler/project_1/analysis/analyst/KEGG_enrichment.csv", header=T)
GO_enrichment <- read.csv(file="/projectnb/bf528/users/wheeler/project_1/analysis/analyst/GO_enrichment.csv", header=T)
Hallmark_enrichment <- read.csv(file="/projectnb/bf528/users/wheeler/project_1/analysis/analyst/Hallmark_enrichment.csv", header=T)

# filter out the top three enriched gene sets
top_3_KEGG <- KEGG_enrichment %>% arrange(padj) %>% slice(1:3)
top_3_GO <- GO_enrichment %>% arrange(padj) %>% slice(1:3)
top_3_Hallmark <- Hallmark_enrichment %>% arrange(padj) %>% slice(1:3)

# write out a copy to csv
write.csv(top_3_KEGG, file="/projectnb/bf528/users/wheeler/project_1/analysis/analyst/top_3_KEGG.csv", row.names=F)
write.csv(top_3_GO, file="/projectnb/bf528/users/wheeler/project_1/analysis/analyst/top_3_GO.csv", row.names=F)
write.csv(top_3_Hallmark, file="/projectnb/bf528/users/wheeler/project_1/analysis/analyst/top_3_Hallmark.csv", row.names=F)




