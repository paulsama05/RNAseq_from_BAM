###Goal: use BAM files from TopHat to generate counts and perform differential gene expression
###Section 1: setting up the R environment and file paths

##Get required packages: here

#Install and load "here" if required
if (!requireNamespace("here", quietly = TRUE)){renv::install("here")}
library(here)


##Create new sub-directories, save this script within "/R/", and assign path to this script

sapply(c("R", "data"), dir.create)
#Use 'here" package to assign path to this script
here::i_am("R/1_env_setup.R")
#and move the working directory to the "data" folder to ensure downloaded/created files are deposited there
setwd("./data")


##Download the appropriate GTF (gene annotation file). In this case, Rattus Norvegicus from ENSEMBL release 88
  #Note: this release must be use because it was used to generate the BAM files

#Directly access the desired build/file from ENSEMBL, which is a 'gtf.gz' file
download.file(url = "ftp://ftp.ensembl.org/pub/release-88/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.88.gtf.gz",
  destfile = "Rattus_norvegicus.Rnor_6.0.88.gtf.gz")

#And unzip the file to a '.gtf'
if (!requireNamespace("R.utils", quietly = TRUE)){renv::install("R.utils")}
library(R.utils)
R.utils::gunzip("Rattus_norvegicus.Rnor_6.0.88.gtf.gz", overwrite = TRUE)
  #NOTE: if any step above fails for whatever reason, the final .gtf file is provided in the ./data folder


##Define location of key files

#BAM files (inputs)
bam_files <- list.files(path = here::here("data/bam/"), pattern = "*.bam$", full.names = TRUE)
#GTF (annotation) file for rattus norvegicus whole genome
annotation <- here::here("data/Rattus_norvegicus.Rnor_6.0.88.gtf")
