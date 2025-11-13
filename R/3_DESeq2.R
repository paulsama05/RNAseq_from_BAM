###Goal: use BAM files from TopHat to generate counts and perform differential gene expression
###Section 3: use the count data from the previous script to perform differential gene expression analysis in DESeq2

##Change current location (to match this script)
here::i_am("R/3_DESeq2.R")


##Install DESeq2 (if necessary)
if (!requireNamespace("DESeq2", quietly = TRUE)){renv::install("bioc::DESeq2")}
library(DESeq2)


##Filter out low-count genes (those that don't have a count of at least 10 across all samples)
  #This threshold is arbitrary, but could be determined statistically by looking across the count data
  #For example, by examining the mean expression of each gene in each experimental group and exploring
    #the distribution of that data

fc_filtered <- fc$counts[rowSums(fc$counts)>=10,]


##Create a data frame to define the grouping of samples
  #NOTE: in this analysis the first four samples are "Treated", and the next 4 "Untreated"

grouping <- data.frame(row.names = colnames(fc_filtered), treatment = as.factor(rep(c("Treated","Untreated"), each = 4)))
#Ensure that the "Untreated" samples are prioritized (and will be recognized as baseline)
grouping$treatment <- relevel(grouping$treatment,"Untreated")


##Run DESeq2

#Create the DESeq object from the count data (fc$counts) and the grouping data frame
DESeq_summary <- DESeqDataSetFromMatrix(countData = fc_filtered, colData = grouping, design = ~treatment)
#Run analysis
DESeq_analysis <- DESeq(DESeq_summary)
#Return the results
DESeq_results <- results(DESeq_analysis, contrast = c("treatment", "Treated","Untreated"))
#Peek at the results and output the full data
head(DESeq_results)
write.csv(DESeq_results,"DESeq_results.csv")

##Reorder by adjusted p-value (low to high)
  #The results object is of DESeq2 class and cannot be modified using dplyr calls
  #You can arrange and subset in the "classic" fashion", e.g.  DESeq_ordered <- DESeq_results[order(DESeq_results$padj),]
  #Or, you can convert to a data frame and modify with tidyverse commands

library(dplyr)
library(tidyr)
DESeq_df <- as.data.frame(DESeq_results)
DESeq_df_trimmed <- DESeq_df %>% arrange(padj) %>% drop_na(padj) %>% mutate(log10padj = -log10(padj))
  #The "trimmed" version drops genes without p-value data, organizes by ascending p-value, and  appends a -log10(padj) column
DESeq_df_sig <- DESeq_df_trimmed %>% filter(padj <0.05)
  #The "sig" version is the same, but only significantly differently expressed genes


##Find genes that are exclusively expressed in one group
  #Note: I arbitrarily define this as 0 expression in one group and a mean of > 1 in the other
  #Because the data set was already filtered for total counts >=10, if one group is zero then the other should have an average of ~2.5

#Find mean expression of each gene by group
K_mean <- rowMeans(fc_filtered[,1:4])
head(K_mean)
V_mean <- rowMeans(fc_filtered[,5:8])
#Subset the data set for genes that are only expressed in one group
only_in_treat <- fc_filtered[K_mean>=1 & V_mean ==0,]
  #In my analysis returns 9 genes unique to the treated group
only_in_unt <- fc_filtered[K_mean==0 & V_mean >=1,]
  #In my analysis returns 12 genes unique to the untreated group
#Check to see if genes from the "unique" sets are present in the DESeq output
sum(rownames(rbind(only_in_treat,only_in_unt)) %in% rownames(DESeq_df_sig))
  #Returns '0', which is to be expected because fold change cannot be calculated when one group has no expression


##The above data sets (differential expression from DESeq and "unique" expression) can be used for further analysis
