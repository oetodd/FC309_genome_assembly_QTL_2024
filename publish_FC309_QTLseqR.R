#######################################

#__Date:__ December 6th, 2022
#__Author:__ KMD 
#__Script:__ QTLseqr_FC309
#__Project:__ To QTL for a sugar beet line with resistance to Fusarium oxysporum strain F19.
#__Corresponding paper:__ A fully phased, chromosome-scale genome of sugar beet line FC309 enables the discovery of Fusarium yellows resistance QTL. 
######################################


#load the package
#install.packages("devtools")
#devtools::install_github("bmansfeld/QTLseqr")

library("QTLseqr")

setwd("~/Desktop")

#Set sample and file names
HighBulk <- "F2_309_R_CKDN230001421-1A_HW3M2DSX5_L2_sorted.bam"
LowBulk <- "F2_309_S_CKDN230001422-1A_HW3M2DSX5_L2_sorted.bam"
file <- "variant_calls.table"

#Choose which chromosomes will be included in the analysis (i.e. exclude smaller contigs)
Chroms <- c("Chr1_Scaffold_1_RagTag", "Chr2_Scaffold_5_RagTag", "Chr3_Scaffold_2_RagTag", "Chr4_Scaffold_3_RagTag", "Chr5_Scaffold_9_RagTag", "Chr6_Scaffold_7_RagTag", "Chr7_Scaffold_4_RagTag", "Chr8_Scaffold_8_RagTag", "Chr9_Scaffold_6_RagTag")

#Import SNP data from file
df <-
  importFromGATK(
    file = file,
    highBulk = HighBulk,
    lowBulk = LowBulk,
    chromList = Chroms
  )

depth = ggplot(data = df) +
  geom_histogram(aes(x = DP.HIGH + DP.LOW)) +
  xlim(0,1000)

ggplot(data = df) +
  geom_histogram(aes(x = REF_FRQ))

high=ggplot(data = df) +
  geom_histogram(aes(x = SNPindex.HIGH))

low=ggplot(data = df) +
  geom_histogram(aes(x = SNPindex.LOW))


#Filter SNPs based on some criteria
df_filt <-
  filterSNPs(
    SNPset = df,
    refAlleleFreq = 0.0,
    minTotalDepth = 100,
    maxTotalDepth = 400,
    minSampleDepth = 40,
    minGQ = 99
  )


#Run G' analysis
df_filt <- runGprimeAnalysis(
  SNPset = df_filt,
  windowSize = 1e6,
  outlierFilter = "deltaSNP")

#Run QTLseq analysis
df_filt <- runQTLseqAnalysis(
  SNPset = df_filt,
  windowSize = 1e6,
  popStruc = "F2",
  bulkSize = c(50, 50),
  replications = 10000,
  intervals = c(95, 99)
)

#Plot
g=plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.001,scaleChroms=TRUE)
delta=plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE,scaleChroms=TRUE)
p=plotQTLStats(SNPset = df_filt, var = "negLog10Pval", plotIntervals = TRUE,scaleChroms=TRUE,plotThreshold=TRUE,q=0.001) + ylim(0,100)

ggsave("g.pdf",g, width=12, height=3, units="in", scale=2)
ggsave("delta.pdf",delta, width=12, height=3, units="in", scale=3)
ggsave("p.pdf",p, width=12, height=5, units="in", scale=2)


#export summary CSV
markers=getQTLTable(SNPset = df_filt, alpha = 0.01, export = TRUE, fileName = "my_BSA_QTL.csv")

