###Goal: use BAM files from TopHat to generate counts and perform differential gene expression
###Section 2: use featurecounts() on designated BAM and GTF files to return transcript counts

##Change current location (to match this script)
here::i_am("R/2_featurecounts.R")

##Get required packages: BiocManager, Rsubread

#Install BiocManager (Bioconductor) if required
if (!requireNamespace("BiocManager", quietly = TRUE)){renv::install("BiocManager")}
#Then install and load Rsubread
renv::install("bioc::Rsubread")
library(Rsubread)


##Run featurecounts() on BAM and GTF files using previously defined addresses
fc <- featureCounts(files = bam_files,
                    annot.ext = annotation,
                    isGTFAnnotationFile = TRUE,
                    GTF.featureType = "exon",
                    GTF.attrType = "gene_id",
                    useMetaFeatures = TRUE,
                    isPairedEnd = TRUE,
                    nthreads = 4)
  #Note that "useMetaFeatures" is set to TRUE to ensure that counts are summarized at the gene level
  #isPairedEnd must also be set to TRUE because the data uses paired end reads. If not, an error is returned


##Check the count data and then write it to a .csv file

head(fc$counts)
  #Note: this file is written to the current working directory, which is the 'data' folder
#Also check summary of assigned vs. unassigned reads by looking at the nested '$stat' data frame
fc$stat
  #In this analysis, ~30% of all reads were unassigned (meets general target of 60-80% assigned reads)
  #These reads are primarily binned to "Unassigned_NoFeatures" (~27%) and "Unassigned_Ambiguity" (~2.7%)
    #The '_NoFeatures" outcome is common, especially in older assemblies with incomplete annotation
    #"_Ambiguity" is also common and is due to overlapping genes or paralogs
write.csv(fc$counts, "gene_counts.csv")