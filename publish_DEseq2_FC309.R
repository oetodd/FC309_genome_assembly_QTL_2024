#######################################

#__Date:__ December 6th, 2022
#__Author:__ OET 
#__Script:__ DEseq2_FC309
#__Project:__ To analyze RNA-seq data that compares 2 sugar beet lines with differential resistance to Fusarium oxysporum strain F19 over 3 timepoints.
#__Corresponding paper:__ A fully phased, chromosome-scale genome of sugar beet line FC309 enables the discovery of Fusarium yellows resistance QTL. 

# This script is to be used for analysis of each haplome for FC309. It doesn't go through every single comparison, but was used and changed for each necessary comparison.
# 
######################################
######### After use, comment this section ##############

# Install required packages/libraries:

#DESeq2: https://bioconductor.org/packages/release/bioc/html/DESeq2.html
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


################################


###########  LOAD PACKAGE  #####################

# Do this every time.
# Load packages:
library(DESeq2)

#if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
#################################################

###########  READ IN THE DATA  #####################

# Read in the counts data 
#setwd(inputdir) # Set working directory
getwd() # Check the working directory

setwd("~/Desktop/FC309/rnaseq/haplome_analysis")

## Make sure to cut fields 2-6 from the counts data beforehand using terminal.

countsData <- read.table(file = "counts_v1.1.0.txt", header = FALSE, row.names = 1, skip = 2) # import the data

# Check the countsData object

head(countsData)

# Read in the metadata
metadata <- read.table(file = "metadata2.txt", header = FALSE)


# Organize the metadata file to use headers that correspond to the metadata file
colnames(metadata) <- c("fasta1", "fasta2", "names1", "names2", "rep", "time")
metadata

# Organize the countsData file.
# Notice that the countsData file doesn't have any column headers:
head(countsData)

# Let's give countsData some columns names.
as.vector(metadata$names2)

# Name countsData columns headers: (May need to delete fields 2-6 in the counts file)
colnames(countsData) <- c(as.vector(metadata$names2))

# :!: EXERCISE: Now look at the top of countsData again using head():

head(countsData)

################### COUNT MATRIX INPUT ###################

# In this section we will prepare our input data for analysis.
# In the instruction for DESeq2, it states: "We read in a count matrix, which we will name cts, and the sample information table, which we will name coldata."

cts<- as.matrix(countsData)

head(cts)
dim(cts)


# Next we need to make an information called coltable. We can make this out of the metadata table.

# Reorganize the metadata table so the names2 column are now headers
metadata
rownames(metadata)<- metadata$names2
metadata

coldata <- metadata[,c("rep", "time")]

# Now we have coldata!

## One thing we need to explicitly check. The rownames of coldata need to exactly match the colnames of cts.
#Check that the names match  --> Should be TRUE
all(rownames(coldata) == colnames(cts))



# Next we will create an ddsHTSeq object out of cts and coldata:
# This will set a base design for the experiment:
# Load all the _counts.txt files and to attach them to the metadata.

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ time) 

################# PRE-FILTERING -- FILTER FOR PRESENT GENES: ################# 
# Not necessary, but helps keep things fast. 
# Exclude all samples that have 0 reads:
keep <- rowSums(counts(dds)) >= 1
dds <- dds[keep,]

# Exercise: How many did we exclude?
dim(dds)


# Organize the categories based on what makes sense:
coldata
dds$time <- factor(dds$time, levels = c("R_untreated","R_6dpi","R_24h", "S_24h","S_untreated", "S_6dpi"))

#dds$time <- factor(dds$time, levels = c("S_24h", "R_24h","S_6dpi", "R_6dpi"))
any(is.na(dds$time))
table(dds$time)

head(dds$time)

  # PERFORM DESEQ2 analysis:
####### This is where the magic happens! ########### 
# This will transform dds into a specialized object with many more fields filled in.

dds <- DESeq(dds)
resultsNames(dds)

#Access the normalized counts using counts(x, normalized = TRUE)
#Access the raw count info using counts()
head(counts(dds, normalized = TRUE))
head(counts(dds, normalized = FALSE))


############## DIFFERENTIAL EXPRESSION ANALYSIS #####################


# calculate the statistically significant differences between treated and untreated tissue. These results are some of those that were used in the paper,
# some variable names need to be reorganized in the section above to allow for some comparisons. See vignette for details.


resultsNames(dds)
###
res_R24h <- results(dds,
                   lfc = 0.5,
                   alpha = 0.05,
                   contrast=c("time", "R_24h", "R_untreated")) #read as or ""condition", "numerator", "denominator"
summary(res_R24h)
###
res_R6dpi <- results(dds,
                    lfc = 0.5,
                    alpha = 0.05,
                    contrast=c("time", "R_6dpi", "R_untreated"))
summary(res_R6dpi)
###
res_unt <- results(dds,
                     lfc = 0.5,
                     alpha = 0.05,
                     contrast=c("time", "S_untreated", "R_untreated"))
summary(res_unt)

## Perform log fold change
# An input requirement of the lfcShrink function is a coef term. This i pulled from the resultsNames of dds:
resultsNames(dds)
#fold change
res_R24h_LFC <- lfcShrink(dds, coef="time_R_24h_vs_R_untreated", res = res_R24h)
res_R6dpi_LFC <- lfcShrink(dds, coef="time_R_6dpi_vs_R_untreated", res = res_R6dpi)
res_unt_LFC <- lfcShrink(dds, coef="time_S_untreated_vs_R_untreated", res = res_unt)

summary(res_unt_LFC)
###############

# re-organize the data to get the susceptible comparisons
dds$time <- factor(dds$time, levels = c("R_24","S_24h","S_6dpi","R_untreated","R_6dpi","S_untreated"))

#re-run DESeq
dds <- DESeq(dds)
resultsNames(dds)

###
res_S24h <- results(dds,
                    lfc = 0.5,
                    alpha = 0.05,
                    contrast=c("time", "S_24h", "S_untreated")) #read as or ""condition", "numerator", "denominator"
summary(res_S24h)
###
res_S6dpi <- results(dds,
                     lfc = 0.5,
                     alpha = 0.05,
                     contrast=c("time", "S_6dpi", "S_untreated"))
summary(res_S6dpi)

## Perform log fold change
# An input requirement of the lfcShrink function is a coef term. This was pulled from the resultsNames of dds:

resultsNames(dds)
#fold change
res_S24h_LFC <- lfcShrink(dds, coef="time_S_24h_vs_S_untreated", res = res_S24h)
res_S6dpi_LFC <- lfcShrink(dds, coef="time_S_6dpi_vs_S_untreated", res = res_S6dpi)

#####
# re-organize the data to get the last comparisons
dds$time <- factor(dds$time, levels = c("R_24h","S_untreated", "S_24h","S_6dpi","R_untreated","R_6dpi"))

#re-run DESeq
dds <- DESeq(dds)
resultsNames(dds)

##
res_24h <- results(dds,
                   lfc = 0.5,
                   alpha = 0.05,
                   contrast=c("time", "S_24h", "R_24h"))
summary(res_24h)
#fold change
res_24h_LFC <- lfcShrink(dds, coef="time_S_24h_vs_R_24h", res = res_24h)
summary(res_24h_LFC)
########

# re-organize the data to get the suscpetible comparisons
dds$time <- factor(dds$time, levels = c("R_6dpi", "R_24h","S_untreated", "S_24h","S_6dpi","R_untreated"))

#re-run DESeq
dds <- DESeq(dds)
resultsNames(dds)

res_6dpi <- results(dds,
                    lfc = 0.5,
                    alpha = 0.05,
                    contrast=c("time", "S_6dpi", "R_6dpi"))
summary(res_6dpi)
#fold change
res_6dpi_LFC <- lfcShrink(dds, coef="time_S_6dpi_vs_R_6dpi", res = res_6dpi)
summary(res_6dpi_LFC)

############## PERFORM LOG FOLD CHANGE SHRINKAGE FOR VISUALIZATION #####################


head(resLFC)

subset(resLFC, resLFC$log2FoldChange > 1.1)
subset(resLFC, resLFC$padj < 0.001)

##################  Exploring and exporting results ##################
############## GETTING LISTS OF SIGNIFICANTLY CHANGING GENES #####################

# Check the results table:
summary(res_R24h_LFC)
head(res_R24h_LFC)
############################# saving comparisons ############################# 

# Select the significant subset of genes that are up-regulated  ____
Upregulated <- subset(res_R24h_LFC, padj < 0.05 & log2FoldChange > 0.5)
Upregulated <- Upregulated[order(Upregulated$padj),] #order them
head(Upregulated)
write.csv(Upregulated, "R_24h_vs_R_Untreated_upregulated_v1.2.0.csv")

#unfiltered genes:
# Select all genes that are differentially expressed 
Upregulated_allh2 <- res_24h_LFC
write.csv(Upregulated_allh2, "S24vsR24_all_v1.2.0.csv")

# Select all genes that are up-regulated 
sixdpi_all <- res_6dpi_LFC
write.csv(sixdpi_all, "S6dpivsR6dpi_all_v1.2.0.csv")

oneday_all <- res_R24h_LFC
write.csv(oneday_all, "R24hvsRunt_all_v1.2.0.csv")

# Select the significant subset of genes that are down-regulated _____
Downregulated <- subset(res_R24h_LFC, padj < 0.05 & log2FoldChange < -0.5)
Downregulated <- Downregulated[order(Downregulated$padj),]
head(Downregulated)
write.csv(Downregulated, "R_24h_vs_R_Untreated_downregulated_v1.2.0.csv")
#
#
# Select the significant subset of genes that are up-regulated  ____
Upregulated <- subset(res_R6dpi_LFC, padj < 0.05 & log2FoldChange > 0.5)
Upregulated <- Upregulated[order(Upregulated$padj),] #order them
head(Upregulated)
write.csv(Upregulated, "R_6dpi_vs_R_Untreated_upregulated_v1.2.0.csv")

# Select the significant subset of genes that are down-regulated _____
Downregulated <- subset(res_R6dpi_LFC, padj < 0.05 & log2FoldChange < -0.5)
Downregulated <- Downregulated[order(Downregulated$padj),]
head(Downregulated)
write.csv(Downregulated, "R_6dpi_vs_R_Untreated_downregulated_v1.2.0.csv")
#
#
# Select the significant subset of genes that are up-regulated  ____
Upregulated <- subset(res_unt_LFC, padj < 0.05 & log2FoldChange > 0.5)
Upregulated <- Upregulated[order(Upregulated$padj),] #order them
head(Upregulated)
write.csv(Upregulated, "S_Untreated_vs_R_Untreated_upregulated_v1.2.0.csv")

# Select the significant subset of genes that are down-regulated _____
Downregulated <- subset(res_unt_LFC, padj < 0.05 & log2FoldChange < -0.5)
Downregulated <- Downregulated[order(Downregulated$padj),]
head(Downregulated)
write.csv(Downregulated, "S_Untreated_vs_R_Untreated_downregulated_v1.2.0.csv")
#
#
# Select the significant subset of genes that are up-regulated  ____
Upregulated <- subset(res_S24h_LFC, padj < 0.05 & log2FoldChange > 0.5)
Upregulated <- Upregulated[order(Upregulated$padj),] #order them
head(Upregulated)
write.csv(Upregulated, "S_24h_vs_S_Untreated_upregulated_v1.2.0.csv")

# Select the significant subset of genes that are down-regulated _____
Downregulated <- subset(res_S24h_LFC, padj < 0.05 & log2FoldChange < -0.5)
Downregulated <- Downregulated[order(Downregulated$padj),]
head(Downregulated)
write.csv(Downregulated, "S_24h_vs_S_Untreated_downregulated_v1.2.0.csv")
#
#
# Select the significant subset of genes that are up-regulated  ____
Upregulated <- subset(res_S6dpi_LFC, padj < 0.05 & log2FoldChange > 0.5)
Upregulated <- Upregulated[order(Upregulated$padj),] #order them
head(Upregulated)
write.csv(Upregulated, "S_6dpi_vs_S_Untreated_upregulated_v1.2.0.csv")

# Select the significant subset of genes that are down-regulated _____
Downregulated <- subset(res_S6dpi_LFC, padj < 0.05 & log2FoldChange < -0.5)
Downregulated <- Downregulated[order(Downregulated$padj),]
head(Downregulated)
write.csv(Downregulated, "S_6dpi_vs_S_Untreated_downregulated_v1.2.0.csv")
#
#
# Select the significant subset of genes that are up-regulated  ____
Upregulated <- subset(res_24h_LFC, padj < 0.05 & log2FoldChange > 0.5)
Upregulated <- Upregulated[order(Upregulated$padj),] #order them
head(Upregulated)
write.csv(Upregulated, "S_24h_vs_R_24h_upregulated_v1.2.0.csv")

# Select the significant subset of genes that are down-regulated _____
Downregulated <- subset(res_24h_LFC, padj < 0.05 & log2FoldChange < -0.5)
Downregulated <- Downregulated[order(Downregulated$padj),]
head(Downregulated)
write.csv(Downregulated, "S_24h_vs_R_24h_downregulated_v1.2.0.csv")
#
#
# Select the significant subset of genes that are up-regulated  ____
Upregulated <- subset(res_6dpi_LFC, padj < 0.05 & log2FoldChange > 0.5)
Upregulated <- Upregulated[order(Upregulated$padj),] #order them
head(Upregulated)
write.csv(Upregulated, "S_6dpi_vs_R_6dpi_upregulated_v1.2.0.csv")

# Select the significant subset of genes that are down-regulated _____
Downregulated <- subset(res_6dpi_LFC, padj < 0.05 & log2FoldChange < -0.5)
Downregulated <- Downregulated[order(Downregulated$padj),]
head(Downregulated)
write.csv(Downregulated, "S_6dpi_vs_R_6dpi_downregulated_v1.2.0.csv")

# Get versions
sessionInfo()
