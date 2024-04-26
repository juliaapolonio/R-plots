# Create a volcano plot for eQTL MR analysis results

mr_result <- vroom("../GSMR/results_eqtlgen.txt")

library(ggplot2)
library(ggrepel)

mr_result <- mr_result %>%
  mutate(Exposure=gsub(x=Exposure,pattern="_GSMR", replacement=""))

png("eqtlgen_volcano.png", width = 1500, height = 1000)

  ggplot(mr_result, mapping=aes(x=bxy, y=-log10(p), label=Exposure))+
    geom_point(aes(size=nsnp))+
    geom_hline(yintercept=-log10(0.05/16000), linetype="dashed", color="red")+
    geom_text_repel(data=subset(mr_result, (p<0.000005 | abs(bxy)>0.2)), size=4, max.overlaps = 1000)+
    theme_classic()

dev.off()
