## Zooplankton horizontal trap figures - Summer 2020
#created 15 Jan 2021

#read in packages
pacman::p_load(grid,ggplot2,tidyverse)

#create function to count characters starting at the end of the string
substrEnd <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#read in zoop data 
zoop<- read.csv('SummaryStats/FCR_ZooplanktonSummary2020.csv',header = TRUE)

#select cols to keep
keep <- c("sample_ID","collect_date","DepthOfTow_m","ZoopDensity_No.pL","BiomassConcentration_ugpL",
          "Daphnia_density_NopL","Daphnia_BiomassConcentration_ugpL","Ceriodaphnia_density_NopL","Ceriodaphnia_BiomassConcentration_ugpL",
          "Calanoida_density_NopL","Calanoida_BiomassConcentration_ugpL","Cyclopoida_density_NopL","Cyclopoida_BiomassConcentration_ugpL",
          "Rotifera_density_NopL","Rotifera_BiomassConcentration_ugpL")

zoop_traps <- zoop[is.na(zoop$mesh_size_μm),keep]

#create new column with trap ID
zoop_traps$trap <- substrEnd(zoop_traps$sample_ID,2)

#top vs bottom column
zoop_traps$position <- ifelse(substr(zoop_traps$sample_ID,18,20)=="bot" | substr(zoop_traps$sample_ID,18,20)=="_bo",
                              zoop_traps$position <- "Bottom",zoop_traps$position <- "Top")

#make trap position a factor and switch order
zoop_traps$position <- factor(zoop_traps$position, levels = c("Top","Bottom"))
#------------------------------------------------------------------------------#
# FIGURES

#labels
dat_text <- data.frame(
  label = c("Top","","","Bottom"),
  position = c("Top","Bottom"),
  trap = c("L1","L2","P1","P2"),
  x     = c(4,4,4,4),
  y     = c(400, 25,405,155))

#jpeg("Figures/Horizontal_traps_totaldens.jpg", width = 6, height = 4, units = "in",res = 300)
#removing the first sunset data on 04Aug because P vs. L traps look pretty similar 
ggplot(data=subset(zoop_traps,collect_date!="2020-08-04"),aes(x=trap, y=ZoopDensity_No.pL, fill=trap)) + 
  geom_rect(data=subset(zoop_traps,collect_date == '2020-08-12'),
            aes(fill=collect_date),xmin=-Inf ,xmax = Inf, ymin = -Inf, ymax = Inf, fill = 'black', alpha = 0.053,inherit.aes = FALSE) +
  geom_bar(stat="identity",position=position_dodge()) +
  facet_wrap(~position+collect_date, scales= 'free_y',ncol=2) +
  scale_fill_manual(values=c("#3CBB75FF","#3CBB75FF","#33638DFF","#33638DFF")) +
  theme(strip.text.x = element_blank(),strip.background = element_blank(),
        panel.background = element_blank(), legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
        plot.margin = margin(2.8,0,2,0,unit = "cm"), axis.title.x=element_blank()) +
  geom_text(data = dat_text,  mapping = aes(x = Inf, y = Inf, label = label),
            hjust   = 1.2, vjust   = 1.4, )
#dev.off()

#jpeg("Figures/Horizontal_traps_rotiferdens.jpg", width = 6, height = 4, units = "in",res = 300)
ggplot(data=subset(zoop_traps,collect_date!="2020-08-04"),aes(x=trap, y=Rotifera_density_NopL, fill=trap)) + 
  geom_rect(data=subset(zoop_traps,collect_date == '2020-08-12'), aes(fill=collect_date),xmin=-Inf ,
            xmax = Inf, ymin = -Inf, ymax = Inf, fill = 'black', alpha = 0.053,inherit.aes = FALSE) +
  geom_bar(stat="identity",position=position_dodge()) + ylim(0,105) +
  facet_wrap(~position+collect_date, scales= 'free_y',ncol=2) +
  scale_fill_manual(values=c("#3CBB75FF","#3CBB75FF","#33638DFF","#33638DFF")) +
  theme(strip.text.x = element_blank(),strip.background = element_blank(),
        panel.background = element_blank(), legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
        plot.margin = margin(2.8,0,2,0,unit = "cm"), axis.title.x=element_blank()) +
  geom_text(data = dat_text,  mapping = aes(x = Inf, y = Inf, label = label),
            hjust   = 1.2, vjust   = 1.4, )
#dev.off()



