# Plot z-scores of sumstats in a cumulative BP x-axis (like a Manhattan plot) 

# Get clumping results
lead_bd_gwas <- read.table("gwas_BD.clumped")
lead_sd_gwas <- read.table("gwas_SD.clumped")
lead_bd_mtag <- read.table("mtag_bd.clumped")
lead_sd_mtag <- read.table("mtag_sd.clumped")


# Tidy clumping data
colnames(lead_bd_gwas) <- lead_bd_gwas[1,]
colnames(lead_sd_gwas) <- lead_sd_gwas[1,]
colnames(lead_bd_mtag) <- lead_bd_mtag[1,]
colnames(lead_sd_mtag) <- lead_sd_mtag[1,]

lead_bd_gwas <- lead_bd_gwas[-1,]
lead_sd_gwas <- lead_sd_gwas[-1,]
lead_bd_mtag <- lead_bd_mtag[-1,]
lead_sd_mtag <- lead_sd_mtag[-1,]


# Get sumstats z-score
z_bdep_gwas <- vroom("depression_sumstats.txt", col_select = c("variant_id", "z_score"))
z_sdep_gwas <- vroom("pgc_depression_sumstats_mtag.txt", col_select = c("variant_id", "z_score"))
z_bdep_mtag <- vroom("outputs/mtag/dep_new_trait_1.txt", col_select = c("SNP", "mtag_z"))
z_sdep_mtag <- vroom("outputs/mtag/sdep_new_trait_1.txt", col_select = c("SNP", "mtag_z"))

# Add z-score to clumped data
lead_bd_gwas <- lead_bd_gwas %>%
  inner_join(z_bdep_gwas, by = c("SNP"="variant_id"))
lead_sd_gwas <- lead_sd_gwas %>%
  inner_join(z_sdep_gwas, by = c("SNP"="variant_id"))
lead_bd_mtag <- lead_bd_mtag %>%
  inner_join(z_bdep_mtag, by = "SNP") 
lead_sd_mtag <- lead_sd_mtag %>%
  inner_join(z_sdep_mtag, by = "SNP")

##################### Create BPcum for all 4 data ###################################################


data <- rbind(lead_bd_gwas, lead_bd_mtag %>% rename(z_score = "mtag_z"), lead_sd_gwas, lead_sd_mtag %>% rename(z_score = "mtag_z"))

### GWAS BDEP
# Create BPcum to find BP position in X axis
data$BP <- as.numeric(data$BP)
data$CHR <- as.numeric(data$CHR)
data <- data %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(data, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)


#################### CREATE AXIS FOR PLOT ##########################################
# Create variables for x axis of Manhattan plot

axisdf = data %>%
  group_by(CHR) %>%
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
             
             
ggplot(data, aes(x=BPcum, y=z_score)) +
  
  # Show all points
  geom_point() +

  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  ylim(-10, 10) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )         
 
