# Some plots of enrichment analysis using different tools

#-----------------------------------------------------
#Ploting Gprofiler enrichment baed on MAGMA full output

library(gprofiler2)

magma_genes <- read.table("magma.genes.out", header = T)

sign_magma_genes <- magma_genes[magma_genes$P<5e-8,]

gostres2 <- gost(sign_magma_genes$GENE, organism = "hsapiens", correction_method = "bonferroni", user_threshold = 0.05)

gostplot(gostres2, capped = TRUE, interactive = TRUE)

#----------------------------------------------------------
# Plotting MAGMA enrichment using first 350 genes from FUMA analysis


magma_enrichment <- read.delim("GS.txt", header = T, sep = "\t")

b <-  magma_enrichment[unique(magma_enrichment$GeneSet),]

library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(forcats)
library(gghighlight)
#---------------------------------------------------------------
# Gwas Catalog


# png("GWAS_catalog.png", width = 1700, height = 1500, res = 250)

c <- magma_enrichment[magma_enrichment$Category == "GWAScatalog",]

ggplot(c,aes(x =-log10(adjP), y = GeneSet, size = N_overlap/N_genes)) +
  geom_point(color = "cadetblue") +
  #gghighlight(GeneSet %in% c("Neurociticism", "Social communication problems", "Schizophrenia", "Depression", "Anorexia nervosa"), use_direct_label = F) +
  #scale_color_gradientn(colours = rev(brewer.pal(11, name = "Spectral"))) +
  xlab("P-value (-log10)") +
  ylab("GWAS Catalog Trait") +
  labs(size = "Gene Ratio", color = "Adjusted\nP-value") +
  theme_classic()

dev.off()

#----------------------------------------------------------
# Plotting MAGMA enrichment
# Gene Ontology

png("Gene_ontology.png", width = 2000, height = 1500, res = 250)

b %>%
  filter(Category %in% c("GO_cc", "GO_mf", "GO_bp")) %>%
  mutate(GeneSet = fct_reorder(GeneSet, N_overlap/N_genes)) %>%
  ggplot(aes(x =N_overlap/N_genes, y = GeneSet, size = N_overlap, color = -log10(adjP))) +
  geom_point() + #colours changed below to a custom pallete w/o white
  scale_color_gradientn(colours = c("#1e3c87","#403784","#57317f","#692977",
                                        "#79216e","#851763","#8f0c56","#960449",
                                        "#9a073b","#9b122d","#9a1e1e"))+
  xlab("Gene Ratio") +
  theme_classic()
dev.off()

#----------------------------------------------------------
# Plotting MAGMA enrichment
# Reactome

png("Reactome.png", width = 3000, height = 4000, res = 250)

magma_enrichment %>%
  filter(Category %in% c("Reactome")) %>%
  mutate(GeneSet = fct_reorder(GeneSet, N_overlap/N_genes)) %>%
  ggplot(aes(x =N_overlap/N_genes, y = GeneSet, size = N_overlap, color = -log10(adjP))) +
  geom_point() +
  scale_color_gradientn(colours = rev(brewer.pal(11, name = "Spectral"))) +
  xlab("Gene Ratio") +
  theme_classic()

dev.off()

#----------------------------------------------------------
# Plotting MAGMA enrichment
# Other genesets

png("Other_genesets.png", width = 2000, height = 3000, res = 250)

magma_enrichment %>%
  filter(Category %in% c("Cell_type_signature", "Chemical_and_Genetic_pertubation",
                         "KEGG", "TF_targets", "Wiki_pathways")) %>%
  mutate(GeneSet = fct_reorder(GeneSet, N_overlap/N_genes)) %>%
  ggplot(aes(x =N_overlap/N_genes, y = GeneSet, size = N_overlap, color = -log10(adjP))) +
  geom_point() +
  scale_color_gradientn(colours = rev(brewer.pal(11, name = "Spectral"))) +
  xlab("Gene Ratio") +
  theme_classic()

dev.off()
#----------------------------------------------------------
# Plotting MAGMA enrichment
# Curated genesets

reactome <- magma_enrichment[magma_enrichment$Category == "Reactome",]
curated_gene_sets <- magma_enrichment[magma_enrichment$Category == "Curated_gene_sets",]

not_reactome_gene_set <- setdiff(curated_gene_sets$GeneSet,reactome$GeneSet)

png("Curated_genesets.png", width = 1500, height = 1500, res = 250)

magma_enrichment %>%
  filter(GeneSet %in% not_reactome_gene_set) %>%
  mutate(GeneSet = fct_reorder(GeneSet, N_overlap/N_genes)) %>%
  ggplot(aes(x =N_overlap/N_genes, y = GeneSet, size = N_overlap, color = -log10(adjP))) +
  geom_point() +
  scale_color_gradientn(colours = rev(brewer.pal(11, name = "Spectral"))) +
  xlab("Gene Ratio") +
  theme_classic()

dev.off()
