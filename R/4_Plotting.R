##To do here:
  #Volcano plot of log2fc vs log10 padj using ggplot
  #GO analysis
    #Run variations using gene lists from DESeq data and/or from "unique" expression sets
  #KEGG analysis


##Volcano plot of DESeq output (log2fc vs -log10padj)
 
#Create new category for color data based on fold change and padj
DESeq_df_trimmed$category <- rep("Insignficant/Low",length(nrow(DESeq_df_trimmed)))
#If padj <0.05 and fold change >= 2 (log2fc >=1), set to "Upregulated"
DESeq_df_trimmed$category[DESeq_df_trimmed$log2FoldChange >=1 & DESeq_df_trimmed$padj < 0.05] <- "Upregulated"
#If padj <0.05 and fold change <= 1/2 (log2fc <=-1), set to "Downregulated"
DESeq_df_trimmed$category[DESeq_df_trimmed$log2FoldChange <= -1 & DESeq_df_trimmed$padj < 0.05] <- "Downregulated"
#Set categories as factors
DESeq_df_trimmed$category <- as.factor(DESeq_df_trimmed$category)
#Output the modified file
write.csv(DESeq_df_trimmed, "DESeq_results_categorized.csv")
#Plot
library(ggplot2)
ggplot(DESeq_df_trimmed, aes(log2FoldChange,log10padj)) + geom_point(aes(col=category)) +labs(title = "Differential Gene Expression", subtitle = "Treated vs. Untreated Tumors", x = "log2(FoldChange)", y = "-log10(padj)") + geom_hline(yintercept = -log10(0.05), linetype = "dotted") + geom_vline(xintercept = 1, linetype = "dotted")  + geom_vline(xintercept = -1, linetype = "dotted")
#save 
ggsave("DGE_volcano.pdf")


##Extract genes for GO and KEGG analysis

#Filter for genes that are "Upregulated" or "Downregulated" as defined in "category" and keep only log2fc data
DESeq_df_final <- DESeq_df_trimmed %>% filter(padj < 0.05 & (log2FoldChange >=1 | log2FoldChange <=-1)) %>% select(log2FoldChange)
  #Alternatively, can use: DESeq_df_reduced <- DESeq_df_trimmed %>% filter(category == "Upregulated" | category == "Downregulated")...
#Also, create a vector combining the up/down genes from DESeq and the "uniquely" expressed genes
Genes_combined <- c(rownames(DESeq_df_final), rownames(only_in_treat), rownames(only_in_unt))  


##Download GO terms from ENSEMBL to use for enrichment analysis

#Install/load biomaRt
if (!requireNamespace("biomaRt", quietly = TRUE)){renv::install("bioc::bioMart")}
library(biomaRt)
#Identify the rattus norvegicus data set in ENSEMBL
listDatasets(useMart(biomart = "ENSEMBL_MART_ENSEMBL"))
  #"rnorvegicus_gene_ensembl" is the target data set
#Identify attributes to extract from the data set
listAttributes(useMart("rnorvegicus_gene_ensembl",biomart="ENSEMBL_MART_ENSEMBL"))
  #We will want "ensembl_gene_id" (Gene Stable ID), "go_id" (GO term accession), "name_1006" (GO term name)
#Pull the desired attributes into a data frame
go_ids <- getBM(attributes = c("ensembl_gene_id","go_id","name_1006"),mart = useMart("rnorvegicus_gene_ensembl",biomart="ENSEMBL_MART_ENSEMBL"))
#Check the table head
head(go_ids)


##Run GO enrichment

#Set up term2gene and term2name data frames and the list of all possible unique gene ids
term2gene <- data.frame(term = go_ids$go_id, gene = go_ids$ensembl_gene_id)
term2name <- data.frame(term = go_ids$go_id, name = go_ids$name_1006 )
all_genes <- unique(go_ids$ensembl_gene_id)
#Install/load clusterProfiler
if (!requireNamespace("clusterProfiler", quietly = TRUE)){renv::install("bioc::clusterProfiler")}
library(clusterProfiler)
#Run GO enrichment
enriched <- enricher(rownames(DESeq_df_final), universe = all_genes, TERM2GENE = term2gene, TERM2NAME = term2name)
#Plot
dotplot(enriched)
  #The resulting plot has a "blank" GO category, which means the "name" value was blank
#Now run GO enrichment with the "combined" gene list
enriched_combo <- enricher(Genes_combined, universe = all_genes, TERM2GENE = term2gene, TERM2NAME = term2name)
#Plot
dotplot(enriched_combo)
  #As above, the resulting plot has a "blank" GO category, which means the "name" value was blank

#Note that many rows of the input data have "" (blank) values for the GO terms and names
dim(go_ids[go_ids$go_id=="" | go_ids$name_1006=="",])
#Use dplyr to filter out rows that have blanks in either category (must NOT be blank in both)
go_ids_filtered <- go_ids %>% dplyr::filter(go_id!="" & name_1006!="")
#Use the filtered data frame to recreate the mini dfs for GO enrichment
term2gene_filtered <- data.frame(term = go_ids_filtered$go_id, gene = go_ids_filtered$ensembl_gene_id)
term2name_filtered <- data.frame(term = go_ids_filtered$go_id, name = go_ids_filtered$name_1006 )
all_genes_filtered <- unique(go_ids_filtered$ensembl_gene_id)
#Run GO enrichment
enriched_filtered <- enricher(rownames(DESeq_df_final), universe = all_genes_filtered, TERM2GENE = term2gene_filtered, TERM2NAME = term2name_filtered)
#Plot
pdf("DGE_sig_GO_enrich.pdf")
dotplot(enriched_filtered)
dev.off()
  #Note that in the new graph the "" category (previously the second most enriched) is now gone and a new category has made the top 10
#Regraph the GO enrichment of the combined gene list as well
enriched_combo_filtered <- enricher(Genes_combined, universe = all_genes_filtered, TERM2GENE = term2gene_filtered, TERM2NAME = term2name_filtered)
pdf("Genes_combo_GO_enrich.pdf")
dotplot(enriched_combo_filtered)
dev.off()

##Run KEGG enrichment

#Look up KEGG code for rattus norvegicus
search_kegg_organism("norv", by = "scientific_name")
  #The appropriate code is 'rno'
#Run the enrichment using the "combined" gene set
kegg_result <- enrichKEGG(gene = Genes_combined, organism = "rno", pvalueCutoff = 0.05)
  #The enrichment fails because the gene names are not provided in the correct format
#Remove the "header" from each gene name, which is "ENSRNOG" plus a series of zeroes
Genes_combined_reformatted<-sub("ENSRNOG0+","",Genes_combined)
#Retry KEGG enrichment
kegg_result <- enrichKEGG(gene = Genes_combined_reformatted, organism = "rno", pvalueCutoff = 0.05)
  #The resulting plot is blank because no terms met the p-value cutoff!
#Rerun with "loose" p-value cutoff
kegg_result <- enrichKEGG(gene = Genes_combined_reformatted, organism = "rno", pvalueCutoff = 0.5)
#Look at enriched pathways
kegg_result@result
#Plot
pdf("Genes_combo_KEGG_enrich.pdf")
dotplot(kegg_result)
dev.off()
#And now do the enrichment and graph using just the DESeq genes that met the significance and fold change criteria
DESeq_df_final_reformatted<-sub("ENSRNOG0+","",rownames(DESeq_df_final))
kegg_result_2 <- enrichKEGG(gene = DESeq_df_final_reformatted, organism = "rno", pvalueCutoff = 0.5)
kegg_result_2@result
pdf("DGE_sig_KEGG_enrich.pdf")
dotplot(kegg_result_2)
dev.off()