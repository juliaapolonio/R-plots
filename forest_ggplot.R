# Create a forest plot using ggplot

library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)

eqtlgen_all=vroom("../GSMR/nextflow/data/results/sdep_eqtlgen/results/results_eqtlgen_sdep.txt")
eqtlgen_bdep=vroom("/data/home/julia.amorim/scripts/GSMR/nextflow/data/results/eqtlgen/results/results_bdep_eqtlgen.txt")

paste0("The number of tested genes before filter is ",nrow(eqtlgen_all))
paste0("The number of tested genes before filter is ",nrow(eqtlgen_bdep))
eqtlgen_result=eqtlgen_all %>% mutate("p"=as.numeric(p)) %>% filter(!is.na(p))
eqtlgen_result=eqtlgen_bdep %>% mutate("p"=as.numeric(p)) %>% filter(!is.na(p))

eqtlgen_result=eqtlgen_result %>% filter(nsnp>3)
paste0("The number of tested genes before filter is ",nrow(eqtlgen_result))
#FDR adjusted p values for remaining probes
eqtlgen_result=eqtlgen_result %>% mutate(Q=p.adjust(p,method="BH"))
eqtlgen_result=eqtlgen_result %>% filter(Q<0.05)
write.table(eqtlgen_result,"eqtlgen_result_IBD_gsmr_BH.tsv",quote=F,sep="\t",row.names=F)


df <- read.table("eqtlgen_result_IBD_gsmr_BH.tsv", T, sep = "\t")
colnames(df) <- c("p1", "p2", "estimate", "stderr", "p", "NSNPs")

df$p1 <- stringr::str_remove(df$p1, "\\_GSMR")

df$upper <- df$estimate+df$stderr
df$lower <- df$estimate-df$stderr

ci.lb = df$lower
ci.ub = df$upper

#taken from forest.default source code
level = 95
alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
vi <- ((ci.ub - ci.lb)/(2 * qnorm(alpha/2, lower.tail = FALSE)))^2
wi <- 1/sqrt(vi)
psize <- wi/sum(wi, na.rm = TRUE)
psize <- (psize - min(psize, na.rm = TRUE))/(max(psize, 
                                                 na.rm = TRUE) - min(psize, na.rm = TRUE))
df$psize <- (psize * 1) + 0.5

ggplot(data=df, aes(y=p1, x=estimate, xmin=estimate-stderr, xmax=estimate+stderr)) +
  geom_point(aes(size=psize)) +
  scale_size_continuous(range = c(0,2)) +
  geom_errorbarh(height=0.5) +
  labs(title='eQTLGen on Strict Depression', x='Effect Size', y = 'Gene') +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  guides(size=FALSE) +
  theme_classic()

ggsave("forest_plot_sdep_ggplot.png",units="in",width=5, height=15, dpi=300)

