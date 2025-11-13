
# RNAseq_from_BAM

<!-- badges: start -->
<!-- badges: end -->

This repository provides the scripts necessary to recreate a simple RNAseq workflow in R using the `Rsubread`, `DESeq2`, and `clusterProfiler` packages.

The `.bam` files used in the original analysis are too large to provide in GitHub. However, all output files are provided here.

## Experiment basics and essential data

The core experiment was not published so details are withheld (but are irrelevant to understanding the workflow herein); only a high level overview is provided.

In brief, mice with INS-1 (rat PanNET cell line) xenografts were treated with vehicle ('V' samples) or inhibitor ('K' samples). RNA was purified in bulk and shipped out for sample processing (oligo dT enrichment, rRNA depletion),  paired end sequencing, and initial data processing (trimming and genomic alignment using TopHat).

The resulting alignments were provided as `.bam` files with corresponding sample names.

## General worklow

- Acquire .bam files and the corresponding gene annotation (.gtf) file.

- Use `featureCounts()` to change alignment data into gene-level count data.

- Use `DESeq2` functions to run differential gene expression (DGE) using the gene-level count data and user-defined sample grouping.

- Plot the DGE data using `ggplot2`.

- Use a sub-list of the most significantly upregulated and downregulated genes to perform GO and KEGG enrichments.

## Example results

See the "/outputs/" folder for example files. Where/how these files are generated is detailed below.

## Workflow highlights

Each R script provided in the "/R/" directory encompasses a major step of the workflow:

### 1_env_setup.R
- The basic R environment is setup using the `renv` and `here` packages.

- Local "/R/" and "/data/" directories are created for your R project.

- The relevant, gene-annotated build of the rattus norvegicus genome is acquired as a `.gtf` file (Release 88, which was used for `.bam` generation in TopHat).

- A path is created to the `.bam` files (while ignoring `.bam.bai` files in the same directory).

### 2_featurecounts.R
- `featureCounts()` from the `Rsubread` package uses the `.bam` files (alignments) and the `.gtf` file (gene annotations) to determine gene counts for each sample.

- A quick quality check is performed to make sure a reasonable number of alignments were assigned to genes.

- The resulting count table is output as the file "gene_counts.csv".


### 3_DESeq2.R
- The count table is filtered using an arbitrary threshold.

- A data frame is created to define the grouping of samples (Untreated/'V' vs. Treated/'K') in the count table.

- Differential expression analysis is run by passing the count table and grouping data frame through `DESeq2`.

- The raw output data frame is saved as "DESeq_results.csv".

- The data frame is then filtered to remove genes without associated p-values (NA) before calculating -log10(p-values) for future analysis.

- Separately, genes are identified from the raw data frame that are ONLY expressed in the vehicle-treated ("V") OR inhibitor-treated ("K") group. These samples effectively show differential expression, but do not return a fold change value or p-value in the above differential gene expression data frame.

### 4_Plotting.R
- The genes in the filtered data frame are assigned as "Downregulated", "Upregulated", or "Insignificant/Low" depending on whether they exceed arbitrary p-value and fold change thresholds

- The "categorized", filtered data frame is saved ("DESeq_results_categorized.csv") and passed to `ggplot2` to generate a volcano plot ("DGE_volcano.pdf").

- Genes that meet the "Downregulated" and "Upregulated" criteria are excised into a new sub-list.

- GO terms correspinding to rattus norvegicus are imported by using `biomaRt` to access the relevant ENSEMBL database.

- This sub-list is used to run GO and KEGG enrichment, resulting in the plots "DGE_sig_GO_enrich.pdf" and "DGE_sig_KEGG_enrich.pdf"

- An "extended" version of the sub-list is created by adding in the names of genes that were only expressed in one sample group.

- GO and KEGG enrichment are rerun using this amended sub-list, resulting in the plots "Genes_combo_GO_enrich.pdf" and "Genes_combo_KEGG_enrich.pdf"
