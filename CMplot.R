# R script to create a Circular Manhattan plot of 2 traits 

# Load necessary libraries
library(CMplot)
library(vroom)

# Set working directory
setwd("/data/home/julia.amorim/scripts/data/")

# Load snps
load("mtag_snps.RData")
load("gwas_snps.RData")

# Load data (original and mtag) first iteration: depression
data_orig <- vroom("depression_sumstats.txt", delim = "\t", col_names = TRUE, comment = "#")
data_mtag <- vroom("../mtag/all_results.1NS_trait_2.txt", delim = "\t", col_names = TRUE, comment = "#")

colnames(data_orig) <- c("CHR", "SNP", "BP", "A2", "A1", "EAF", "SE", "P", "CIU", "OR", "CIL", "Z", "N")
colnames(data_mtag) <- c("SNP", "CHR", "BP", "A1", "A2", "Z", "N", "EAF", "BETA", "SE", "Z", "P")
 
# Create a new dataframe with p-values from both traits in CMplot format
munged_data <- merge (data_orig[,c(1:3,9)], data_mtag[,c(1:3,12)], by = c("SNP", "CHR", "BP"), all = T)

# Depression GWAS SNPs
orig_snps <- gwas_snps[3]
orig_snps <- unlist(orig_snps, use.names = F)

# Create a vector with MTAG lead SNPs
lead_snps <- mtag_snps[3]
lead_snps <- unlist(lead_snps, use.names = F)

# CMplot for MTAG
CMplot(munged_data, col=c("#72CC9E","#A09174"), type="p",plot.type="c",r=0.2,cir.axis=FALSE, threshold.col = "#167d57",
      outward=TRUE,cir.axis.col="black", band = 0, axis.cex = 1, threshold = 5e-8, highlight = lead_snps,
      amplify = FALSE, cir.axis.grid = FALSE,cir.chr.h=0.3,chr.den.col=NULL,file="jpg", legend.ncol = 2,
      file.name="Mtag-depression",dpi=700,file.output=TRUE,verbose=TRUE,width=10,height=10, highlight.col = "#9c2858")

# CMplot for GWAS
CMplot(munged_data, col=c("#72CC9E","#A09174"), type="p",plot.type="c",r=0.2,cir.axis=FALSE, threshold.col = "#167d57",
       outward=TRUE,cir.axis.col="black", band = 0, axis.cex = 1, threshold = 5e-8, highlight = orig_snps,
       amplify = FALSE, cir.axis.grid = FALSE,cir.chr.h=0.3,chr.den.col=NULL,file="jpg", legend.ncol = 2,
       file.name="Gwas-depression",dpi=700,file.output=TRUE,verbose=TRUE,width=10,height=10, highlight.col = "#32a89b")


# CMplot qq-plot
#CMplot(munged_data, col=c("#72CC9E","#A09174"), type="p", multraits = TRUE, plot.type="q",threshold = 1e-6,
#       cir.axis.col="black",axis.cex = 1,file="jpg", conf.int = FALSE,
#       file.name="qqplot",dpi=400,file.output=TRUE,verbose=TRUE,width=10,height=10)

# Calculate lambda
#chisq <- qchisq(1-munged_data$Mtag, 1)
#lambda = median(chisq)/qchisq(0.5,1)

