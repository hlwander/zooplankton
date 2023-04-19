## Zooplankton horizontal trap figures - Summer 2020
#created 15 Jan 2021

#read in packages
pacman::p_load(grid,ggplot2,tidyverse)

#create function to count characters starting at the end of the string
substrEnd <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#Calculates the standard error####
stderr <- function(x) {
  sd(x,na.rm=TRUE)/sqrt(length(na.omit(x)))
}

#read in zoop data 
traps_20<- read.csv('Summer2020-DataAnalysis/SummaryStats/FCR_ZooplanktonSummary2020.csv',header = TRUE) %>% filter(site_no %in% "BVR_trap")
traps_21<- read.csv('Summer2021-DataAnalysis/SummaryStats/FCR_ZooplanktonSummary2021.csv',header = TRUE) %>% filter(site_no %in% "BVR_trap")

#combine data
zoop_traps <- rbind(traps_20,traps_21)

#select cols to keep
keep <- c("sample_ID","collect_date","DepthOfTow_m","ZoopDensity_No.pL","BiomassConcentration_ugpL",
          "Cladocera_density_NopL", "Cladocera_BiomassConcentration_ugpL", "Cladocera_PercentOfTotal",
          "Copepoda_density_NopL", "Copepoda_BiomassConcentration_ugpL", "Copepoda_PercentOfTotal",
          "Rotifera_density_NopL", "Rotifera_BiomassConcentration_ugpL", "Rotifera_PercentOfTotal")

zoop_traps <- zoop_traps[,keep]

#create new rep column 
zoop_traps$rep <- substrEnd(zoop_traps$sample_ID,1)

#create new column with trap ID
zoop_traps$trap <- substr(substrEnd(zoop_traps$sample_ID,2),1,1)

#top vs bottom column
zoop_traps$position <- ifelse(substr(zoop_traps$sample_ID,18,20)=="bot" | substr(zoop_traps$sample_ID,18,20)=="_bo",
                              zoop_traps$position <- "Bottom",zoop_traps$position <- "Top")

#make trap position a factor and switch order
zoop_traps$position <- factor(zoop_traps$position, levels = c("Top","Bottom"))

#add time column
zoop_traps$time <- ifelse(grepl("sunrise",zoop_traps$sample_ID),"sunrise","sunset")

##### Create new df to combine reps over 24 hours
trap_repmeans <- zoop_traps %>% 
  group_by(collect_date , trap, position, time) %>%
  summarise_at(vars(ZoopDensity_No.pL:Rotifera_PercentOfTotal), funs(rep.mean=mean, rep.SE=stderr))

#wide to long
temp_trap <-  trap_repmeans %>%
  gather(metric,value,ZoopDensity_No.pL_rep.mean:Rotifera_PercentOfTotal_rep.mean) 
temp_trap_SE <-  trap_repmeans %>%
  gather(metric,value,ZoopDensity_No.pL_rep.SE:Rotifera_PercentOfTotal_rep.SE) 

trap_repmeans_long <- temp_trap[,c(1:4,16,17)]
trap_repmeans_long$value_SE <- temp_trap_SE$value

#rename metric column so it is shorter
trap_repmeans_long$metric <- substr(trap_repmeans_long$metric,1,nchar(trap_repmeans_long$metric)-14)

#add a column to combine trap and time
trap_repmeans_long$lab <- paste0(trap_repmeans_long$trap,"_",substr(trap_repmeans_long$position,1,3))

#------------------------------------------------------------------------------#
# FIGURES

#rename facet labels
taxa <- c("Total","Total","Cladocera","Cladocera","Cladocera","Copepoda","Copepoda","Copepoda","Rotifera","Rotifera","Rotifera")
names(taxa) <- unique(trap_repmeans_long$metric)

#reorder labels
trap_repmeans_long$lab <- factor(trap_repmeans_long$lab,levels=c("P_Top","L_Top","P_Bot","L_Bot"))

ggplot(subset(trap_repmeans_long, grepl("density",metric,ignore.case = TRUE) & collect_date!="2020-08-04"), 
       aes(x=lab, y=value, fill=collect_date, alpha=position)) +
  geom_rect(data=subset(trap_repmeans_long,time == 'sunset' &grepl("density",metric,ignore.case = TRUE) & collect_date!="2020-08-04"),
            aes(fill=time),xmin=-Inf ,xmax = Inf, ymin = -Inf, ymax = Inf, fill = 'black', alpha = 0.01,inherit.aes = FALSE) +
  geom_bar(stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() + guides(fill=guide_legend(ncol=6))+
  geom_errorbar(aes(ymin=value-value_SE, ymax=value+value_SE), width=.2,position=position_dodge(.9)) + scale_alpha_manual(values = c(0.4, 1)) +
  facet_wrap(time+position~metric, scales= 'free',ncol=4, strip.position = "right", labeller=labeller(metric=taxa)) + ylab("Density (#/L)") +
  scale_fill_manual(values=c("#33CCCC","#339999","#CCCCFF","#9966CC","#FF9966","#FF6633")) +
  theme(strip.text.y = element_text(size = 11 ,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank()) + 
 theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
       axis.text.x=element_text(size=6,family="Times"), plot.title = element_text(hjust = 0.5),
       legend.position = "bottom", legend.margin = margin(0, 1, 2, 1),plot.margin = margin(0,0,0,0.3,unit = "cm"), 
       axis.text.y = element_text(size=13, family="Times"), legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),
       legend.text = element_text(size=10),legend.title = element_blank()) + xlab("") 
#ggsave(file.path(getwd(),"Summer2021-DataAnalysis/Figures/Horizontal_traps_density.jpg"))

ggplot(subset(trap_repmeans_long, grepl("biomass",metric,ignore.case = TRUE) & collect_date!="2020-08-04"), 
       aes(x=lab, y=value, fill=collect_date, alpha=position)) +
  geom_rect(data=subset(trap_repmeans_long,time == 'sunset' &grepl("biomass",metric,ignore.case = TRUE) & collect_date!="2020-08-04"),
            aes(fill=time),xmin=-Inf ,xmax = Inf, ymin = -Inf, ymax = Inf, fill = 'black', alpha = 0.01,inherit.aes = FALSE) +
  geom_bar(stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() + guides(fill=guide_legend(ncol=6))+
  geom_errorbar(aes(ymin=value-value_SE, ymax=value+value_SE), width=.2,position=position_dodge(.9)) + scale_alpha_manual(values = c(0.4, 1)) +
  facet_wrap(time+position~metric, scales= 'free',ncol=4, strip.position = "right", labeller=labeller(metric=taxa)) + ylab("Biomass (ug/L)") +
  scale_fill_manual(values=c("#33CCCC","#339999","#CCCCFF","#9966CC","#FF9966","#FF6633")) +
  theme(strip.text.y = element_text(size = 11 ,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank()) + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text.x=element_text(size=6,family="Times"), plot.title = element_text(hjust = 0.5),
        legend.position = "bottom", legend.margin = margin(0, 1, 2, 1),plot.margin = margin(0,0,0,0.3,unit = "cm"), 
        axis.text.y = element_text(size=13, family="Times"), legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),
        legend.text = element_text(size=10),legend.title = element_blank()) + xlab("") 
#ggsave(file.path(getwd(),"Summer2021-DataAnalysis/Figures/Horizontal_traps_biomass.jpg"))

ggplot(subset(trap_repmeans_long, grepl("percent",metric,ignore.case = TRUE) & collect_date!="2020-08-04"), 
       aes(x=lab, y=value, fill=collect_date, alpha=position)) +
  geom_rect(data=subset(trap_repmeans_long,time == 'sunset' &grepl("percent",metric,ignore.case = TRUE) & collect_date!="2020-08-04"),
            aes(fill=time),xmin=-Inf ,xmax = Inf, ymin = -Inf, ymax = Inf, fill = 'black', alpha = 0.01,inherit.aes = FALSE) +
  geom_bar(stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() + guides(fill=guide_legend(ncol=6))+
  geom_errorbar(aes(ymin=value-value_SE, ymax=value+value_SE), width=.2,position=position_dodge(.9)) + scale_alpha_manual(values = c(0.4, 1)) +
  facet_wrap(time+position~metric, scales= 'free',ncol=3, strip.position = "right", labeller=labeller(metric=taxa)) + ylab("Percent Density (%)") +
  scale_fill_manual(values=c("#33CCCC","#339999","#CCCCFF","#9966CC","#FF9966","#FF6633")) +
  theme(strip.text.y = element_text(size = 11 ,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank()) + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text.x=element_text(size=6,family="Times"), plot.title = element_text(hjust = 0.5),
        legend.position = "bottom", legend.margin = margin(0, 1, 2, 1),plot.margin = margin(0,0,0,0.3,unit = "cm"), 
        axis.text.y = element_text(size=13, family="Times"), legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),
        legend.text = element_text(size=10),legend.title = element_blank()) + xlab("") 
#ggsave(file.path(getwd(),"Summer2021-DataAnalysis/Figures/Horizontal_traps_percent_dens.jpg"))

#------------------------------------------------------------------------------#
#remaking figures for each date pair to better see magnitudes

ggplot(subset(trap_repmeans_long, grepl("percent",metric,ignore.case = TRUE) & 
                collect_date %in% c("2020-08-12", "2020-08-13")), 
       aes(x=lab, y=value, fill=collect_date, alpha=position)) +
  geom_rect(data=subset(trap_repmeans_long,time == 'sunset' &grepl("percent",metric,ignore.case = TRUE) & 
                          collect_date %in% c("2020-08-12", "2020-08-13")),
            aes(fill=time),xmin=-Inf ,xmax = Inf, ymin = -Inf, ymax = Inf, fill = 'black', alpha = 0.01,inherit.aes = FALSE) +
  geom_bar(stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() + guides(fill=guide_legend(ncol=6))+
  geom_errorbar(aes(ymin=value-value_SE, ymax=value+value_SE), width=.2,position=position_dodge(.9)) + scale_alpha_manual(values = c(0.4, 1)) +
  facet_wrap(time+position~metric, scales= 'free',ncol=3, strip.position = "right", labeller=labeller(metric=taxa)) + #change ncol to 4 for dens and biom
  scale_fill_manual(values=c("#33CCCC","#339999")) + ylab("Percent density (%)") +
  theme(strip.text.y = element_text(size = 11 ,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank()) + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text.x=element_text(size=6,family="Times"), plot.title = element_text(hjust = 0.5),
        legend.position = "bottom", legend.margin = margin(0, 1, 2, 1),plot.margin = margin(0,0,0,0.3,unit = "cm"), 
        axis.text.y = element_text(size=13, family="Times"), legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),
        legend.text = element_text(size=10),legend.title = element_blank()) + xlab("") 
#ggsave(file.path(getwd(),"Summer2021-DataAnalysis/Figures/Horizontal_traps_d1_percent_dens.jpg"))

ggplot(subset(trap_repmeans_long, grepl("percent",metric,ignore.case = TRUE) & 
                collect_date %in% c("2021-06-15", "2021-06-16")), 
       aes(x=lab, y=value, fill=collect_date, alpha=position)) +
  geom_rect(data=subset(trap_repmeans_long,time == 'sunset' &grepl("percent",metric,ignore.case = TRUE) & 
                          collect_date %in% c("2021-06-15", "2021-06-16")),
            aes(fill=time),xmin=-Inf ,xmax = Inf, ymin = -Inf, ymax = Inf, fill = 'black', alpha = 0.01,inherit.aes = FALSE) +
  geom_bar(stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() + guides(fill=guide_legend(ncol=6))+
  geom_errorbar(aes(ymin=value-value_SE, ymax=value+value_SE), width=.2,position=position_dodge(.9)) + scale_alpha_manual(values = c(0.4, 1)) +
  facet_wrap(time+position~metric, scales= 'free',ncol=3, strip.position = "right", labeller=labeller(metric=taxa)) + 
  scale_fill_manual(values=c("#CCCCFF","#9966CC")) + ylab("Percent density (%)") +
  theme(strip.text.y = element_text(size = 11 ,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank()) + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text.x=element_text(size=6,family="Times"), plot.title = element_text(hjust = 0.5),
        legend.position = "bottom", legend.margin = margin(0, 1, 2, 1),plot.margin = margin(0,0,0,0.3,unit = "cm"), 
        axis.text.y = element_text(size=13, family="Times"), legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),
        legend.text = element_text(size=10),legend.title = element_blank()) + xlab("") 
#ggsave(file.path(getwd(),"Summer2021-DataAnalysis/Figures/Horizontal_traps_d2_percent_dens.jpg"))

ggplot(subset(trap_repmeans_long, grepl("percent",metric,ignore.case = TRUE) & 
                collect_date %in% c("2021-07-07", "2021-07-08")), 
       aes(x=lab, y=value, fill=collect_date, alpha=position)) +
  geom_rect(data=subset(trap_repmeans_long,time == 'sunset' &grepl("percent",metric,ignore.case = TRUE) & 
                          collect_date %in% c("2021-07-07", "2021-07-08")),
            aes(fill=time),xmin=-Inf ,xmax = Inf, ymin = -Inf, ymax = Inf, fill = 'black', alpha = 0.01,inherit.aes = FALSE) +
  geom_bar(stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() + guides(fill=guide_legend(ncol=6))+
  geom_errorbar(aes(ymin=value-value_SE, ymax=value+value_SE), width=.2,position=position_dodge(.9)) + scale_alpha_manual(values = c(0.4, 1)) +
  facet_wrap(time+position~metric, scales= 'free',ncol=3, strip.position = "right", labeller=labeller(metric=taxa)) + 
  scale_fill_manual(values=c("#FF9966","#FF6633")) + ylab("Percent density (%)") +
  theme(strip.text.y = element_text(size = 11 ,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank()) + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text.x=element_text(size=6,family="Times"), plot.title = element_text(hjust = 0.5),
        legend.position = "bottom", legend.margin = margin(0, 1, 2, 1),plot.margin = margin(0,0,0,0.3,unit = "cm"), 
        axis.text.y = element_text(size=13, family="Times"), legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),
        legend.text = element_text(size=10),legend.title = element_blank()) + xlab("") 
#ggsave(file.path(getwd(),"Summer2021-DataAnalysis/Figures/Horizontal_traps_d3_percent_dens.jpg"))


