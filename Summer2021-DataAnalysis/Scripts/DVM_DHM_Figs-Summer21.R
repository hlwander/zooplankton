##Visualizing all of the data
#Created 17Nov2020
#Updated for 2021 on 20Oct21

### Description of data --> full water column tows and epi tows from FCR + BVR summer 2020 (outside of 12-13Aug 20 MSN)
    #includes samples collected at macrophytes (BVR_l) and pelagic site (BVR_50_p for epi tows collected during MSN ONLY; BVR_50 for full water column tows and tows outside of 24-hour campaigns)
    #samples collected from noon (x1), midnight (x1), sunset (x4), and sunrise (x4)


#TO DO:
#figure out rotifer biomass calcs (simple formula is sig different than other equation with l,w,h params, at least for keratella)

#read in libraries
pacman::p_load(plyr,plotrix,lubridate,dplyr,ggplot2,scales,tidyr,viridis)

#read in zoop summary csv
zoop<- read.csv('./Summer2021-DataAnalysis/SummaryStats/FCR_ZooplanktonSummary2021.csv',header = TRUE)

#create function to count characters starting at the end of the string
substrEnd <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#Calculates the standard error####
stderr <- function(x) {
  sd(x,na.rm=TRUE)/sqrt(length(na.omit(x)))
}

#make sure sample_ID is class character
zoop$sample_ID<- as.character(zoop$sample_ID)

#merge collect_date and hour in a new column
zoop$date<- paste(zoop$collect_date,zoop$Hour,sep=" ")
#get times into date format (character here)
zoop$date<- format(as.POSIXct(zoop$date,format="%Y-%m-%d %H:%M"), format="%Y-%m-%d %H:%M:%S")
#convert to posixct date format
zoop$date<- as.POSIXct(zoop$date, format="%Y-%m-%d %H:%M")

#order by site, then hour
zoop<- arrange(zoop,zoop$site_no,zoop$date)

#pull rep # off as new column
zoop$rep <- ifelse(substrEnd(zoop$sample_ID,4)=="rep1" |substrEnd(zoop$sample_ID,4)=="rep2" | substrEnd(zoop$sample_ID,4)=="rep3"| substrEnd(zoop$sample_ID,4)=="rep4",
   substrEnd(zoop$sample_ID,1),NA)

#drop rep# from sample ID
zoop$sample_ID <- ifelse(substrEnd(zoop$sample_ID,4)=="rep1" |substrEnd(zoop$sample_ID,4)=="rep2" | substrEnd(zoop$sample_ID,4)=="rep3" | substrEnd(zoop$sample_ID,4)=="rep4",
                  substr(zoop$sample_ID,1,nchar(zoop$sample_ID)-5),zoop$sample_ID)

#get hour into character format for grouping
zoop$Hour <- format(round(strptime(paste0(zoop$collect_date, zoop$Hour), format="%Y-%m-%d %H:%M"),units="hours"),format="%H:%M")
#manually change hour of some samples (rounding problems)
zoop$Hour[grepl("midnight",zoop$sample_ID,ignore.case = TRUE)] <- "00:00"
zoop$Hour[grepl("noon",zoop$sample_ID,ignore.case = TRUE)] <- "12:00"

#drop schindler samples and only select BVR_50
zoop <- zoop[substrEnd(zoop$site_no,6)!="schind" & zoop$site_no!="FCR_50",]

#Create new df with more specific taxonomic groupings (Daphnia, Ceriodaphnia, Bosminiidae, calanoida, cyclopoida, nauplius, keratella, kellicottia, Collothecidae ,Conochilidae, Synchaetidae, Trichocercidae)
zoop.repmeans.bytaxa <- zoop %>% select(sample_ID,site_no,collect_date,Hour, Daphnia_density_NopL, Daphnia_BiomassConcentration_ugpL, Daphnia_totalbiomass_ug,
                                        Ceriodaphnia_density_NopL, Ceriodaphnia_BiomassConcentration_ugpL, Ceriodaphnia_totalbiomass_ug,
                                        Bosminidae_density_NopL, Bosminidae_BiomassConcentration_ugpL, Bosminidae_totalbiomass_ug,
                                        Calanoida_density_NopL, Calanoida_BiomassConcentration_ugpL, Calanoida_totalbiomass_ug,
                                        Cyclopoida_density_NopL, Cyclopoida_BiomassConcentration_ugpL,Cyclopoida_totalbiomass_ug,
                                        nauplius_density_NopL, nauplius_BiomassConcentration_ugpL, nauplius_totalbiomass_ug,
                                        Keratella_density_NopL, Keratella_BiomassConcentration_ugpL, Keratella_totalbiomass_ug,
                                        Kellicottia_density_NopL, Kellicottia_BiomassConcentration_ugpL, Kellicottia_totalbiomass_ug,
                                        Collothecidae_density_NopL, Collothecidae_BiomassConcentration_ugpL, Collothecidae_totalbiomass_ug,
                                        Conochilidae_density_NopL, Conochilidae_BiomassConcentration_ugpL, Conochilidae_totalbiomass_ug,
                                        Synchaetidae_density_NopL, Synchaetidae_BiomassConcentration_ugpL, Synchaetidae_totalbiomass_ug,
                                        Trichocercidae_density_NopL, Trichocercidae_BiomassConcentration_ugpL,Trichocercidae_totalbiomass_ug) %>%
                                  group_by(sample_ID, site_no, Hour, collect_date) %>%
                                  summarise_at(vars(Daphnia_density_NopL:Trichocercidae_totalbiomass_ug), funs(rep.mean=mean, rep.SE=stderr))


##### Create new df to combine reps over 24 hours
zoop.repmeans <- zoop %>% select(sample_ID,site_no,collect_date,Hour, ZoopDensity_No.pL, BiomassConcentration_ugpL,
                  TotalBiomass_ug,Cladocera_density_NopL, Cladocera_BiomassConcentration_ugpL, Cladocera_totalbiomass_ug, Cladocera_PercentOfTotal,
                  Cyclopoida_density_NopL, Cyclopoida_BiomassConcentration_ugpL, Cyclopoida_totalbiomass_ug, Cyclopoida_PercentOfTotal,
                  Rotifera_density_NopL, Rotifera_totalbiomass_ug, Rotifera_BiomassConcentration_ugpL,Rotifera_PercentOfTotal,
                  Calanoida_PercentOfTotal, Calanoida_density_NopL, Calanoida_BiomassConcentration_ugpL, Calanoida_totalbiomass_ug,
                  Copepoda_PercentOfTotal, Copepoda_density_NopL, Copepoda_BiomassConcentration_ugpL, Copepoda_totalbiomass_ug,
                  nauplius_PercentOfTotal, nauplius_density_NopL, nauplius_BiomassConcentration_ugpL, nauplius_totalbiomass_ug) %>%
                  group_by(sample_ID, site_no, Hour, collect_date) %>%
                  summarise_at(vars(ZoopDensity_No.pL:nauplius_totalbiomass_ug,), funs(rep.mean=mean, rep.SE=stderr))

#get hour into posixct for graphing
zoop.repmeans$Hour <- strptime(paste0(as.character(zoop.repmeans$collect_date), zoop.repmeans$Hour),format="%Y-%m-%d %H:%M")
zoop.repmeans$Hour <- as.POSIXct(zoop.repmeans$Hour)

zoop.repmeans.bytaxa$Hour <- strptime(paste0(as.character(zoop.repmeans.bytaxa$collect_date), zoop.repmeans.bytaxa$Hour),format="%Y-%m-%d %H:%M")
zoop.repmeans.bytaxa$Hour <- as.POSIXct(zoop.repmeans.bytaxa$Hour)

#make sure zoop.repmeans is a dataframe
zoop.repmeans <- data.frame(zoop.repmeans)
zoop.repmeans.bytaxa <- data.frame(zoop.repmeans.bytaxa)

#order by hour for plotting
zoop.repmeans <- zoop.repmeans[order(zoop.repmeans$Hour),]
zoop.repmeans.bytaxa <- zoop.repmeans.bytaxa[order(zoop.repmeans.bytaxa$Hour),]

## Use round_date to round to nearest hour for plotting
r <- round_date(zoop.repmeans$Hour, "hours")

#subsetting repmeans into new df for DHM analyses/figs
zoop_DHM<- zoop.repmeans[ grepl("epi",zoop.repmeans$sample_ID) |grepl("sunrise",zoop.repmeans$sample_ID) | grepl("sunset",zoop.repmeans$sample_ID) | 
                            zoop.repmeans$site_no=="BVR_l",substrEnd(colnames(zoop.repmeans),14)!="Total_rep.mean" & substrEnd(colnames(zoop.repmeans),12)!="Total_rep.SE"]

zoop_DHM_bytaxa<- zoop.repmeans.bytaxa[ grepl("epi",zoop.repmeans.bytaxa$sample_ID) |grepl("sunrise",zoop.repmeans.bytaxa$sample_ID) | grepl("sunset",zoop.repmeans.bytaxa$sample_ID) | 
                                          zoop.repmeans.bytaxa$site_no=="BVR_l",substrEnd(colnames(zoop.repmeans.bytaxa),14)!="Total_rep.mean" & substrEnd(colnames(zoop.repmeans.bytaxa),12)!="Total_rep.SE"]

#convert new dfs fron tibble to dataframe 
zoop_DHM <- data.frame(zoop_DHM)
zoop_DHM_bytaxa <- data.frame(zoop_DHM_bytaxa)

#-------------------#
#  12-13 Aug 2020   #
#-------------------#

##Fig - time vs. total density and biomass for MSN #1 (15-16 Jun 2022)

variables <- c("ZoopDensity_No.pL_rep.mean","BiomassConcentration_ugpL_rep.mean","Cladocera_density_NopL_rep.mean","Cladocera_BiomassConcentration_ugpL_rep.mean",
               "Cyclopoida_density_NopL_rep.mean","Cyclopoida_BiomassConcentration_ugpL_rep.mean","Rotifera_density_NopL_rep.mean","Rotifera_BiomassConcentration_ugpL_rep.mean",
               "Calanoida_density_NopL_rep.mean","Calanoida_BiomassConcentration_ugpL_rep.mean","Copepoda_density_NopL_rep.mean","Copepoda_BiomassConcentration_ugpL_rep.mean")
SE <- c("ZoopDensity_No.pL_rep.SE","BiomassConcentration_ugpL_rep.SE","Cladocera_density_NopL_rep.SE","Cladocera_BiomassConcentration_ugpL_rep.SE",
        "Cyclopoida_density_NopL_rep.SE","Cyclopoida_BiomassConcentration_ugpL_rep.SE","Rotifera_density_NopL_rep.SE","Rotifera_BiomassConcentration_ugpL_rep.SE",
        "Calanoida_density_NopL_rep.SE","Calanoida_BiomassConcentration_ugpL_rep.SE","Copepoda_density_NopL_rep.SE","Copepoda_BiomassConcentration_ugpL_rep.SE")

#only select density and biomass #/ug / L cols
zoop_DHM <- zoop_DHM[,-c(which(grepl("ug_rep",colnames(zoop_DHM))))]

#convert df from wide to long
zoop_DHM_long <- zoop_DHM %>% gather(metric,value,all_of(variables))
zoop_DHM_long <- zoop_DHM_long %>% gather(metric.SE,value.SE, all_of(SE))

#drop _rep.mean from all metric names
zoop_DHM_long$metric <- substr(zoop_DHM_long$metric,1,nchar(zoop_DHM_long$metric)-9)

sites <- c("Pelagic","Littoral")
names(sites) <- c("BVR_50","BVR_l")

#total density and biomass Jun
ggplot(subset(zoop_DHM_long, metric %in% c("ZoopDensity_No.pL","BiomassConcentration_ugpL") & site_no %in% c("BVR_50","BVR_l") & (collect_date=="2021-06-15" | collect_date=="2021-06-16")), aes(Hour,value,color=site_no)) +
  geom_rect(aes(xmin=as.POSIXct("2021-06-15 11:30:00"),xmax=as.POSIXct("2021-06-15 20:41:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2021-06-15 20:42:00"),xmax=as.POSIXct("2021-06-16 05:59:00"), ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2021-06-16 06:00:00"),xmax=as.POSIXct("2021-06-16 12:30:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + facet_grid(site_no~metric,scales="free_y",labeller = labeller(site_no=sites)) + xlab("")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.background = element_blank(), legend.key = element_blank(), legend.key.height=unit(0.5,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_rect(fill = "transparent"), legend.position = c(0.93,0.95), legend.spacing = unit(0.05, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ scale_x_datetime(expand = c(0,0)) +
  scale_color_manual(values=c("#0099CC","#009966"), guide="none") + geom_line()+ #ylab("Density (Individuals/L)")+
  scale_fill_manual("",values=c("#CCCCCC","white"), guide = guide_legend(keyheight = .7,keywidth = 1.7, order=2, override.aes = list(alpha = 1,color="black")))+
  geom_errorbar(aes(ymin=value-value.SE, ymax=value+value.SE), width=.2,position=position_dodge(.9))
ggsave(file.path(getwd(),"Figures/BVR_2021_jun_zoopdensity_biomass_LittoralvsPelagic.jpg")) 


#taxa density Jun
ggplot(subset(zoop_DHM_long, metric %in% c("Cladocera_density_NopL","Cyclopoida_density_NopL","Rotifera_density_NopL","Calanoida_density_NopL") &
                site_no %in% c("BVR_50","BVR_l") & (collect_date=="2021-06-15" | collect_date=="2021-06-16")), aes(Hour,value,color=site_no)) +
  geom_rect(aes(xmin=as.POSIXct("2021-06-15 11:30:00"),xmax=as.POSIXct("2021-06-15 20:41:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2021-06-15 20:42:00"),xmax=as.POSIXct("2021-06-16 05:59:00"), ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2021-06-16 06:00:00"),xmax=as.POSIXct("2021-06-16 12:30:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + facet_grid(site_no~metric,scales="free_y",labeller = labeller(site_no=sites)) + xlab("")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.background = element_blank(), legend.key = element_blank(), legend.key.height=unit(0.5,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_rect(fill = "transparent"), legend.position = c(0.89,0.15), legend.spacing = unit(0.05, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ scale_x_datetime(expand = c(0,0)) +
  scale_color_manual(values=c("#0099CC","#009966"), guide="none") + geom_line()+ ylab("Density (Individuals/L)")+
  scale_fill_manual("",values=c("#CCCCCC","white"), guide = guide_legend(keyheight = .7,keywidth = 1.7, order=2, override.aes = list(alpha = 1,color="black")))+
  geom_errorbar(aes(ymin=value-value.SE, ymax=value+value.SE), width=.2,position=position_dodge(.9))
ggsave(file.path(getwd(),"Figures/BVR_2021_jun_taxa_density_LittoralvsPelagic.jpg")) 

#taxa biomass Jun
ggplot(subset(zoop_DHM_long, metric %in% c("Cladocera_BiomassConcentration_ugpL","Cyclopoida_BiomassConcentration_ugpL","Rotifera_BiomassConcentration_ugpL","Calanoida_BiomassConcentration_ugpL") &
                site_no %in% c("BVR_50","BVR_l") & (collect_date=="2021-06-15" | collect_date=="2021-06-16")), aes(Hour,value,color=site_no)) +
  geom_rect(aes(xmin=as.POSIXct("2021-06-15 11:30:00"),xmax=as.POSIXct("2021-06-15 20:41:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2021-06-15 20:42:00"),xmax=as.POSIXct("2021-06-16 05:59:00"), ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2021-06-16 06:00:00"),xmax=as.POSIXct("2021-06-16 12:30:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + facet_grid(site_no~metric,scales="free_y",labeller = labeller(site_no=sites)) + xlab("")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.background = element_blank(), legend.key = element_blank(), legend.key.height=unit(0.5,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_rect(fill = "transparent"), legend.position = c(0.89,0.9), legend.spacing = unit(0.05, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ scale_x_datetime(expand = c(0,0)) +
  scale_color_manual(values=c("#0099CC","#009966"), guide="none") + geom_line()+ ylab(expression(paste("Biomass (",mu,"g/L)")))+
  scale_fill_manual("",values=c("#CCCCCC","white"), guide = guide_legend(keyheight = .7,keywidth = 1.7, order=2, override.aes = list(alpha = 1,color="black")))+
  geom_errorbar(aes(ymin=value-value.SE, ymax=value+value.SE), width=.2,position=position_dodge(.9))
ggsave(file.path(getwd(),"Figures/BVR_2021_jun_taxa_biomass_LittoralvsPelagic.jpg"))


#total density and biomass Jul
ggplot(subset(zoop_DHM_long, metric %in% c("ZoopDensity_No.pL","BiomassConcentration_ugpL") & site_no %in% c("BVR_50","BVR_l") & (collect_date=="2021-07-07" | collect_date=="2021-07-08")), aes(Hour,value,color=site_no)) +
  geom_rect(aes(xmin=as.POSIXct("2021-07-07 11:30:00"),xmax=as.POSIXct("2021-07-07 20:42:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2021-07-07 20:43:00"),xmax=as.POSIXct("2021-07-08 06:06:00"), ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2021-07-08 06:07:00"),xmax=as.POSIXct("2021-07-08 12:30:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + facet_grid(site_no~metric,scales="free_y",labeller = labeller(site_no=sites)) + xlab("")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.background = element_blank(), legend.key = element_blank(), legend.key.height=unit(0.5,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_rect(fill = "transparent"), legend.position = c(0.12,0.39), legend.spacing = unit(0.05, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ scale_x_datetime(expand = c(0,0)) +
  scale_color_manual(values=c("#0099CC","#009966"), guide="none") + geom_line()+ #ylab("Density (Individuals/L)")+
  scale_fill_manual("",values=c("#CCCCCC","white"), guide = guide_legend(keyheight = .7,keywidth = 1.7, order=2, override.aes = list(alpha = 1,color="black")))+
  geom_errorbar(aes(ymin=value-value.SE, ymax=value+value.SE), width=.2,position=position_dodge(.9))
ggsave(file.path(getwd(),"Figures/BVR_2021_jul_zoopdensity_biomass_LittoralvsPelagic.jpg")) 


#taxa density Jul
ggplot(subset(zoop_DHM_long, metric %in% c("Cladocera_density_NopL","Cyclopoida_density_NopL","Rotifera_density_NopL","Calanoida_density_NopL") &
                site_no %in% c("BVR_50","BVR_l") & (collect_date=="2021-07-07" | collect_date=="2021-07-08")), aes(Hour,value,color=site_no)) +
  geom_rect(aes(xmin=as.POSIXct("2021-07-07 11:30:00"),xmax=as.POSIXct("2021-07-07 20:42:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2021-07-07 20:43:00"),xmax=as.POSIXct("2021-07-08 06:06:00"), ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2021-07-08 06:07:00"),xmax=as.POSIXct("2021-07-08 12:30:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + facet_grid(site_no~metric,scales="free_y",labeller = labeller(site_no=sites)) + xlab("")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.background = element_blank(), legend.key = element_blank(), legend.key.height=unit(0.5,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_rect(fill = "transparent"), legend.position = c(0.12,0.39), legend.spacing = unit(0.05, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ scale_x_datetime(expand = c(0,0)) +
  scale_color_manual(values=c("#0099CC","#009966"), guide="none") + geom_line()+ ylab("Density (Individuals/L)")+
  scale_fill_manual("",values=c("#CCCCCC","white"), guide = guide_legend(keyheight = .7,keywidth = 1.7, order=2, override.aes = list(alpha = 1,color="black")))+
  geom_errorbar(aes(ymin=value-value.SE, ymax=value+value.SE), width=.2,position=position_dodge(.9))
ggsave(file.path(getwd(),"Figures/BVR_2021_jul_taxa_density_LittoralvsPelagic.jpg")) 

#taxa biomass Jul
ggplot(subset(zoop_DHM_long, metric %in% c("Cladocera_BiomassConcentration_ugpL","Cyclopoida_BiomassConcentration_ugpL","Rotifera_BiomassConcentration_ugpL","Calanoida_BiomassConcentration_ugpL") &
                site_no %in% c("BVR_50","BVR_l") & (collect_date=="2021-07-07" | collect_date=="2021-07-08")), aes(Hour,value,color=site_no)) +
  geom_rect(aes(xmin=as.POSIXct("2021-07-07 11:30:00"),xmax=as.POSIXct("2021-07-07 20:42:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2021-07-07 20:43:00"),xmax=as.POSIXct("2021-07-08 06:06:00"), ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2021-07-08 06:07:00"),xmax=as.POSIXct("2021-07-08 12:30:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + facet_grid(site_no~metric,scales="free_y",labeller = labeller(site_no=sites)) + xlab("")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.background = element_blank(), legend.key = element_blank(), legend.key.height=unit(0.5,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_rect(fill = "transparent"), legend.position = c(0.89,0.9), legend.spacing = unit(0.05, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ scale_x_datetime(expand = c(0,0)) +
  scale_color_manual(values=c("#0099CC","#009966"), guide="none") + geom_line()+ ylab(expression(paste("Biomass (",mu,"g/L)")))+
  scale_fill_manual("",values=c("#CCCCCC","white"), guide = guide_legend(keyheight = .7,keywidth = 1.7, order=2, override.aes = list(alpha = 1,color="black")))+
  geom_errorbar(aes(ymin=value-value.SE, ymax=value+value.SE), width=.2,position=position_dodge(.9))
ggsave(file.path(getwd(),"Figures/BVR_2021_jul_taxa_biomass_LittoralvsPelagic.jpg"))


#-----------------------------------------------------------------------------------------------------#
##Fig - simplified density plot over time to show littoral DHM 

#taxa density Jun
ggplot(subset(zoop_DHM_long, metric %in% c("Cladocera_density_NopL","Cyclopoida_density_NopL","Rotifera_density_NopL","Calanoida_density_NopL") &
                site_no %in% c("BVR_l") & (collect_date=="2021-06-15" | collect_date=="2021-06-16")), aes(Hour,value,color=site_no)) +
  geom_rect(aes(xmin=as.POSIXct("2021-06-15 11:30:00"),xmax=as.POSIXct("2021-06-15 20:41:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2021-06-15 20:42:00"),xmax=as.POSIXct("2021-06-16 05:59:00"), ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2021-06-16 06:00:00"),xmax=as.POSIXct("2021-06-16 12:30:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + facet_wrap(~metric,scales="free_y",labeller = labeller(site_no=sites), ncol=2) + xlab("")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.background = element_blank(), legend.key = element_blank(), legend.key.height=unit(0.5,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_rect(fill = "transparent"), legend.position = c(0.1,0.98), legend.spacing = unit(0.05, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ scale_x_datetime(expand = c(0,0)) +
  scale_color_manual(values=c("#009966"), guide="none") + geom_line()+ ylab("Density (Individuals/L)")+ guides(color="none") +
  scale_fill_manual("",values=c("#CCCCCC","white"), guide = guide_legend(keyheight = .7,keywidth = 1.7, order=2, override.aes = list(alpha = 1,color="black")))+
  geom_errorbar(aes(ymin=value-value.SE, ymax=value+value.SE), width=.2,position=position_dodge(.9))
ggsave(file.path(getwd(),"Figures/BVR_2021_jun_taxa_density_Littoral.jpg")) 


#taxa density Jul
ggplot(subset(zoop_DHM_long, metric %in% c("Cladocera_density_NopL","Cyclopoida_density_NopL","Rotifera_density_NopL","Calanoida_density_NopL") &
                site_no %in% c("BVR_l") & (collect_date=="2021-07-07" | collect_date=="2021-07-08")), aes(Hour,value,color=site_no)) +
  geom_rect(aes(xmin=as.POSIXct("2021-07-07 11:30:00"),xmax=as.POSIXct("2021-07-07 20:42:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2021-07-07 20:43:00"),xmax=as.POSIXct("2021-07-08 06:06:00"), ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2021-07-08 06:07:00"),xmax=as.POSIXct("2021-07-08 12:30:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + facet_wrap(~metric,scales="free_y",labeller = labeller(site_no=sites), ncol = 2) + xlab("")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.background = element_blank(), legend.key = element_blank(), legend.key.height=unit(0.5,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_rect(fill = "transparent"), legend.position = c(0.1,0.98), legend.spacing = unit(0.05, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ scale_x_datetime(expand = c(0,0)) +
  scale_color_manual(values=c("#009966"), guide="none") + geom_line()+ ylab("Density (Individuals/L)")+ guides(color="none") +
  scale_fill_manual("",values=c("#CCCCCC","white"), guide = guide_legend(keyheight = .7,keywidth = 1.7, order=2, override.aes = list(alpha = 1,color="black")))+
  geom_errorbar(aes(ymin=value-value.SE, ymax=value+value.SE), width=.2,position=position_dodge(.9))
ggsave(file.path(getwd(),"Figures/BVR_2021_jul_taxa_density_Littoral.jpg"))


#taxa biomass Jun
ggplot(subset(zoop_DHM_long, metric %in% c("Cladocera_BiomassConcentration_ugpL","Cyclopoida_BiomassConcentration_ugpL","Rotifera_BiomassConcentration_ugpL","Calanoida_BiomassConcentration_ugpL") &
                site_no %in% c("BVR_l") & (collect_date=="2021-06-15" | collect_date=="2021-06-16")), aes(Hour,value,color=site_no)) +
  geom_rect(aes(xmin=as.POSIXct("2021-06-15 11:30:00"),xmax=as.POSIXct("2021-06-15 20:41:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2021-06-15 20:42:00"),xmax=as.POSIXct("2021-06-16 05:59:00"), ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2021-06-16 06:00:00"),xmax=as.POSIXct("2021-06-16 12:30:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + facet_wrap(~metric,scales="free_y",labeller = labeller(site_no=sites), ncol=2) + xlab("")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.background = element_blank(), legend.key = element_blank(), legend.key.height=unit(0.5,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_rect(fill = "transparent"), legend.position = c(0.1,0.98), legend.spacing = unit(0.05, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ scale_x_datetime(expand = c(0,0)) +
  scale_color_manual(values=c("#009966"), guide="none") + geom_line()+ ylab(expression(paste("Biomass (",mu,"g/L)")))+ guides(color="none") +
  scale_fill_manual("",values=c("#CCCCCC","white"), guide = guide_legend(keyheight = .7,keywidth = 1.7, order=2, override.aes = list(alpha = 1,color="black")))+
  geom_errorbar(aes(ymin=value-value.SE, ymax=value+value.SE), width=.2,position=position_dodge(.9))
ggsave(file.path(getwd(),"Figures/BVR_2021_jun_taxa_biomass_Littoral.jpg")) 


#taxa biomass Jul
ggplot(subset(zoop_DHM_long, metric %in% c("Cladocera_BiomassConcentration_ugpL","Cyclopoida_BiomassConcentration_ugpL","Rotifera_BiomassConcentration_ugpL","Calanoida_BiomassConcentration_ugpL") &
                site_no %in% c("BVR_l") & (collect_date=="2021-07-07" | collect_date=="2021-07-08")), aes(Hour,value,color=site_no)) +
  geom_rect(aes(xmin=as.POSIXct("2021-07-07 11:30:00"),xmax=as.POSIXct("2021-07-07 20:42:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2021-07-07 20:43:00"),xmax=as.POSIXct("2021-07-08 06:06:00"), ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2021-07-08 06:07:00"),xmax=as.POSIXct("2021-07-08 12:30:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + facet_wrap(~metric,scales="free_y",labeller = labeller(site_no=sites), ncol = 2) + xlab("")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.background = element_blank(), legend.key = element_blank(), legend.key.height=unit(0.5,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_rect(fill = "transparent"), legend.position = c(0.1,0.98), legend.spacing = unit(0.05, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ scale_x_datetime(expand = c(0,0)) +
  scale_color_manual(values=c("#009966"), guide="none") + geom_line()+ ylab(expression(paste("Biomass (",mu,"g/L)")))+ guides(color="none") +
  scale_fill_manual("",values=c("#CCCCCC","white"), guide = guide_legend(keyheight = .7,keywidth = 1.7, order=2, override.aes = list(alpha = 1,color="black")))+
  geom_errorbar(aes(ymin=value-value.SE, ymax=value+value.SE), width=.2,position=position_dodge(.9))
ggsave(file.path(getwd(),"Figures/BVR_2021_jul_taxa_biomass_Littoral.jpg"))
  
#-------------------------------------------------------------------------------------------------------#
##Fig - more taxa specific density over 24 hrs

variables <- c("Daphnia_density_NopL_rep.mean","Daphnia_BiomassConcentration_ugpL_rep.mean","Ceriodaphnia_density_NopL_rep.mean",
               "Ceriodaphnia_BiomassConcentration_ugpL_rep.mean","Bosminidae_density_NopL_rep.mean","Bosminidae_BiomassConcentration_ugpL_rep.mean",
               "Calanoida_density_NopL_rep.mean","Calanoida_BiomassConcentration_ugpL_rep.mean","Cyclopoida_density_NopL_rep.mean",
               "Cyclopoida_BiomassConcentration_ugpL_rep.mean","nauplius_density_NopL_rep.mean","nauplius_BiomassConcentration_ugpL_rep.mean",
               "Keratella_density_NopL_rep.mean","Keratella_BiomassConcentration_ugpL_rep.mean","Kellicottia_density_NopL_rep.mean",
               "Kellicottia_BiomassConcentration_ugpL_rep.mean","Collothecidae_density_NopL_rep.mean", "Collothecidae_BiomassConcentration_ugpL_rep.mean",
               "Conochilidae_density_NopL_rep.mean","Conochilidae_BiomassConcentration_ugpL_rep.mean","Synchaetidae_density_NopL_rep.mean",
               "Synchaetidae_BiomassConcentration_ugpL_rep.mean","Trichocercidae_density_NopL_rep.mean","Trichocercidae_BiomassConcentration_ugpL_rep.mean")

SE <- c("Daphnia_density_NopL_rep.SE","Daphnia_BiomassConcentration_ugpL_rep.SE","Ceriodaphnia_density_NopL_rep.SE",
        "Ceriodaphnia_BiomassConcentration_ugpL_rep.SE","Bosminidae_density_NopL_rep.SE","Bosminidae_BiomassConcentration_ugpL_rep.SE",
        "Calanoida_density_NopL_rep.SE","Calanoida_BiomassConcentration_ugpL_rep.SE","Cyclopoida_density_NopL_rep.SE",
        "Cyclopoida_BiomassConcentration_ugpL_rep.SE","nauplius_density_NopL_rep.SE","nauplius_BiomassConcentration_ugpL_rep.SE",
        "Keratella_density_NopL_rep.SE","Keratella_BiomassConcentration_ugpL_rep.SE","Kellicottia_density_NopL_rep.SE",
        "Kellicottia_BiomassConcentration_ugpL_rep.SE","Collothecidae_density_NopL_rep.SE", "Collothecidae_BiomassConcentration_ugpL_rep.SE",
        "Conochilidae_density_NopL_rep.SE","Conochilidae_BiomassConcentration_ugpL_rep.SE","Synchaetidae_density_NopL_rep.SE",
        "Synchaetidae_BiomassConcentration_ugpL_rep.SE","Trichocercidae_density_NopL_rep.SE","Trichocercidae_BiomassConcentration_ugpL_rep.SE")

#only select density and biomass #/ug / L cols
zoop_DHM_bytaxa <- zoop_DHM_bytaxa[,-c(which(grepl("ug_rep",colnames(zoop_DHM_bytaxa))))]

#convert df from wide to long (kinda hacky way bc having problems doing this)
df1 <- zoop_DHM_bytaxa %>% gather(metric,value,all_of(variables))
df2 <- zoop_DHM_bytaxa %>% gather(metric.SE,value.SE, all_of(SE))

#cut and paste to merge df
zoop_DHM_bytaxa_long <- df1[,c(1:4,29:30)]
zoop_DHM_bytaxa_long$metric.SE <- df2$metric.SE
zoop_DHM_bytaxa_long$value.SE <- df2$value.SE

#drop _rep.mean from all metric names
zoop_DHM_bytaxa_long$metric <- substr(zoop_DHM_bytaxa_long$metric,1,nchar(zoop_DHM_bytaxa_long$metric)-9)

#change BVR_50_p to BVR_50
zoop_DHM_bytaxa_long$site_no[zoop_DHM_bytaxa_long$site_no=="BVR_50_p"] <- "BVR_50"

metric_taxa <-c("Bosminidae","Bosminidae","Calanoida","Calanoida","Ceriodaphnia","Ceriodaphnia",
                "Collothecidae","Collothecidae","Conochilidae","Conochilidae",          
                "Cyclopoida","Cyclopoida","Daphnia","Daphnia","Kellicottia","Kellicottia",
                "Keratella","Keratella","nauplius","nauplius","Synchaetidae","Synchaetidae",
                "Trichocercidae","Trichocercidae")
names(metric_taxa) <- c(sort(unique(zoop_DHM_bytaxa_long$metric)))

#taxa density Jun 15-16
ggplot(subset(zoop_DHM_bytaxa_long, grepl("density",metric,ignore.case = TRUE) &
                site_no %in% c("BVR_50","BVR_l") & (collect_date=="2021-06-15" | collect_date=="2021-06-16")), aes(Hour,value,color=site_no)) +
  geom_rect(aes(xmin=as.POSIXct("2021-06-15 11:30:00"),xmax=as.POSIXct("2021-06-15 20:41:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2021-06-15 20:42:00"),xmax=as.POSIXct("2021-06-16 05:59:00"), ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2021-06-16 06:00:00"),xmax=as.POSIXct("2021-06-16 12:30:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + facet_grid(site_no~metric,scales="free_y",labeller = labeller(site_no=sites,metric=metric_taxa)) + xlab("")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=6), axis.text = element_text(size=6, color="black"), legend.background = element_blank(), legend.key = element_blank(), legend.key.height=unit(0.5,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_rect(fill = "transparent"), legend.position = c(0.96,0.96), legend.spacing = unit(0.05, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ scale_x_datetime(expand = c(0,0)) +
  scale_color_manual(values=c("#0099CC","#009966"), guide = 'none') + geom_line()+ ylab("Density (Individuals/L)")+
  scale_fill_manual("",values=c("#CCCCCC","white"), guide = guide_legend(keyheight = .5,keywidth = .5, order=2, override.aes = list(alpha = 1,color="black")))+
  geom_errorbar(aes(ymin=value-value.SE, ymax=value+value.SE), width=.2,position=position_dodge(.9))
ggsave(file.path(getwd(),"Figures/BVR_2021_jun_taxa2_density_LittoralvsPelagic.jpg")) 

#taxa biomass Jun 15-16
ggplot(subset(zoop_DHM_bytaxa_long, grepl("biomass",metric,ignore.case = TRUE) &
                site_no %in% c("BVR_50","BVR_l") & (collect_date=="2021-06-15" | collect_date=="2021-06-16")), aes(Hour,value,color=site_no)) +
  geom_rect(aes(xmin=as.POSIXct("2021-06-15 11:30:00"),xmax=as.POSIXct("2021-06-15 20:41:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2021-06-15 20:42:00"),xmax=as.POSIXct("2021-06-16 05:59:00"), ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2021-06-16 06:00:00"),xmax=as.POSIXct("2021-06-16 12:30:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + facet_grid(site_no~metric,scales="free_y",labeller = labeller(site_no=sites,metric=metric_taxa)) + xlab("")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=6), axis.text = element_text(size=6, color="black"), legend.background = element_blank(), legend.key = element_blank(), legend.key.height=unit(0.5,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_rect(fill = "transparent"), legend.position = c(0.96,0.96), legend.spacing = unit(0.05, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ scale_x_datetime(expand = c(0,0)) +
  scale_color_manual(values=c("#0099CC","#009966"), guide = 'none') + geom_line()+ ylab(expression(paste("Biomass (",mu,"g/L)")))+
  scale_fill_manual("",values=c("#CCCCCC","white"), guide = guide_legend(keyheight = .5,keywidth = .5, order=2, override.aes = list(alpha = 1,color="black")))+
  geom_errorbar(aes(ymin=value-value.SE, ymax=value+value.SE), width=.2,position=position_dodge(.9))
ggsave(file.path(getwd(),"Figures/BVR_2021_jun_taxa2_biomass_LittoralvsPelagic.jpg"))

#taxa density Jul 7-8
ggplot(subset(zoop_DHM_bytaxa_long, grepl("density",metric,ignore.case = TRUE) &
                site_no %in% c("BVR_50","BVR_l") & (collect_date=="2021-07-07" | collect_date=="2021-07-08")), aes(Hour,value,color=site_no)) +
  geom_rect(aes(xmin=as.POSIXct("2021-07-07 11:30:00"),xmax=as.POSIXct("2021-07-07 20:42:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2021-07-07 20:43:00"),xmax=as.POSIXct("2021-07-08 06:06:00"), ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2021-07-08 06:07:00"),xmax=as.POSIXct("2021-07-08 12:30:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + facet_grid(site_no~metric,scales="free_y",labeller = labeller(site_no=sites,metric=metric_taxa)) + xlab("")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=6), axis.text = element_text(size=6, color="black"), legend.background = element_blank(), legend.key = element_blank(), legend.key.height=unit(0.5,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_rect(fill = "transparent"), legend.position = c(0.96,0.96), legend.spacing = unit(0.05, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ scale_x_datetime(expand = c(0,0)) +
  scale_color_manual(values=c("#0099CC","#009966"), guide = 'none') + geom_line()+ ylab("Density (Individuals/L)")+
  scale_fill_manual("",values=c("#CCCCCC","white"), guide = guide_legend(keyheight = .5,keywidth = .5, order=2, override.aes = list(alpha = 1,color="black")))+
  geom_errorbar(aes(ymin=value-value.SE, ymax=value+value.SE), width=.2,position=position_dodge(.9))
ggsave(file.path(getwd(),"Figures/BVR_2021_jul_taxa2_density_LittoralvsPelagic.jpg")) 

#taxa biomass Jul 7-8
ggplot(subset(zoop_DHM_bytaxa_long, grepl("biomass",metric,ignore.case = TRUE) &
                site_no %in% c("BVR_50","BVR_l") & (collect_date=="2021-07-07" | collect_date=="2021-07-08")), aes(Hour,value,color=site_no)) +
  geom_rect(aes(xmin=as.POSIXct("2021-07-07 11:30:00"),xmax=as.POSIXct("2021-07-07 20:42:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2021-07-07 20:43:00"),xmax=as.POSIXct("2021-07-08 06:06:00"), ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2021-07-08 06:07:00"),xmax=as.POSIXct("2021-07-08 12:30:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + facet_grid(site_no~metric,scales="free_y",labeller = labeller(site_no=sites,metric=metric_taxa)) + xlab("")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=6), axis.text = element_text(size=6, color="black"), legend.background = element_blank(), legend.key = element_blank(), legend.key.height=unit(0.5,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_rect(fill = "transparent"), legend.position = c(0.96,0.96), legend.spacing = unit(0.05, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ scale_x_datetime(expand = c(0,0)) +
  scale_color_manual(values=c("#0099CC","#009966"), guide = 'none') + geom_line()+ ylab(expression(paste("Biomass (",mu,"g/L)")))+
  scale_fill_manual("",values=c("#CCCCCC","white"), guide = guide_legend(keyheight = .5,keywidth = .5, order=2, override.aes = list(alpha = 1,color="black")))+
  geom_errorbar(aes(ymin=value-value.SE, ymax=value+value.SE), width=.2,position=position_dodge(.9))
ggsave(file.path(getwd(),"Figures/BVR_2021_jul_taxa2_biomass_LittoralvsPelagic.jpg"))


#--------------------------------------#
#       BVR epi vs hypo figures!!      #
#Note: oxycline = epi sample depth (4m)#
#--------------------------------------#

#sum up counts by sample/site/day for DVM analyses + figs
BVR_counts <- zoop %>% select(sample_ID,site_no,collect_date,Hour, OverallCount_n,
  CladoceraCount_n, CyclopoidaCount_n, RotiferaCount_n, CalanoidaCount_n, CopepodaCount_n, naupliusCount_n) %>%
  group_by(sample_ID, site_no, Hour, collect_date) %>%
  summarise_at(vars(OverallCount_n:naupliusCount_n), list(rep.mean=mean, rep.SE=stderr))

#add unadjusted volume
BVR_counts$Volume_unadj<- (zoop %>% select(sample_ID,site_no,collect_date,Hour, Volume_unadj) %>%
  group_by(sample_ID, site_no, Hour, collect_date) %>%
  summarise_at(vars(Volume_unadj),list(rep.mean=mean)))$rep.mean

#add proportional volume for numerator hypo calc
BVR_counts$proportional_vol<- (zoop %>% select(sample_ID,site_no,collect_date,Hour, proportional_vol) %>%
  group_by(sample_ID, site_no, Hour, collect_date) %>%
  summarise_at(vars(proportional_vol),list(rep.mean=mean)))$rep.mean

#get BVR_counts df in same order as zoop.repmeans df
BVR_counts<- BVR_counts[order(match(paste0(BVR_counts$sample_ID,BVR_counts$site_no,BVR_counts$collect_date), 
                                    paste0(zoop.repmeans$sample_ID,zoop.repmeans$site_no,zoop.repmeans$collect_date))),]

#add counts and vol to zoop.repmeans
zoop.repmeans[,paste0(colnames(BVR_counts[5:20]))]<- BVR_counts[5:20]

#new dfs for DVM data (BVR_pelagic_DVM_raw is just raw #/ug; BVR_pelagic_DVM_vol_calculated is #/L and ug/L)
BVR_pelagic_DVM<- zoop.repmeans[(zoop.repmeans$site_no=="BVR_50" |zoop.repmeans$site_no=="BVR_50_p") &
                                 (substrEnd(zoop.repmeans$sample_ID,5)=="night" | substrEnd(zoop.repmeans$sample_ID,4)=="noon" |
                                 substrEnd(zoop.repmeans$sample_ID,9)=="night_epi" | substrEnd(zoop.repmeans$sample_ID,8)=="noon_epi" |
                                 substrEnd(zoop.repmeans$sample_ID,3)=="oxy"),] 

#remove oxy samples 
BVR_pelagic_DVM <- BVR_pelagic_DVM[!c(substrEnd(BVR_pelagic_DVM$sample_ID,3)=="oxy"),]

#only select volume, count, and ug cols (plus SE cols)
BVR_pelagic_DVM_raw<- BVR_pelagic_DVM[,c(1:4,73,74,which(substrEnd(colnames(BVR_pelagic_DVM),10)=="n_rep.mean"),
                                     which(substrEnd(colnames(BVR_pelagic_DVM),11)=="ug_rep.mean"),
                                     which(substrEnd(colnames(BVR_pelagic_DVM),8)=="n_rep.SE"),
                                     which(substrEnd(colnames(BVR_pelagic_DVM),9)=="ug_rep.SE"))]

BVR_pelagic_DVM_vol_calculated <- BVR_pelagic_DVM[,c(1:4,which(substrEnd(colnames(BVR_pelagic_DVM),16)=="y_No.pL_rep.mean"),
                                     which(substrEnd(colnames(BVR_pelagic_DVM),15)=="y_NopL_rep.mean"),
                                     which(substrEnd(colnames(BVR_pelagic_DVM),15)=="n_ugpL_rep.mean"),
                                     which(substrEnd(colnames(BVR_pelagic_DVM),13)=="n_ugpL_rep.SE"),
                                     which(substrEnd(colnames(BVR_pelagic_DVM),14)=="y_No.pL_rep.SE"),
                                     which(substrEnd(colnames(BVR_pelagic_DVM),13)=="y_NopL_rep.SE"))]

#another one for percent calcs
BVR_pelagic_DVM_percent<- BVR_pelagic_DVM[,c(1:4,73,74,which(substrEnd(colnames(BVR_pelagic_DVM),14)=="Total_rep.mean"),
                                     which(substrEnd(colnames(BVR_pelagic_DVM),12)=="Total_rep.SE"))]
  
#initialize df
BVR.DVM.calcs<- data.frame("Hour"=unique(BVR_pelagic_DVM_raw$Hour))
  
#for loop to fill out epi vs hypo calcs 
#hypo density and biomass calculated by subtracting epi raw zoop # from full zoop # and then dividing by the (full volume - epi volume) 
#NOTE: using epi density/L and biomass/L but calculating hypo using raw # and ug values. 
column.names<- colnames(BVR_pelagic_DVM_vol_calculated[,c(5:18)])
variables<- colnames(BVR_pelagic_DVM_raw[,c(7:20)])
percent<- colnames(BVR_pelagic_DVM_percent[,c(7:12)])
for(i in 1:length(variables)){
  BVR.DVM.calcs[,paste0(column.names,"_epi")[i]]<- BVR_pelagic_DVM_vol_calculated[substrEnd(BVR_pelagic_DVM_vol_calculated$sample_ID,3)=="epi",paste0(column.names)[i]]
  BVR.DVM.calcs[,paste0(column.names,"_hypo")[i]] <- (((1/BVR_pelagic_DVM_raw[substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)!="epi" ,"proportional_vol"]) * BVR_pelagic_DVM_raw[substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)!="epi" ,paste0(variables)[i]] * (1/0.051)) - 
                                                        ((1/BVR_pelagic_DVM_raw[substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)=="epi", "proportional_vol"]) *  BVR_pelagic_DVM_raw[substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)=="epi",paste0(variables)[i]]* (1/0.31)))/
                                                     (BVR_pelagic_DVM_raw[substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)!="epi" ,"Volume_unadj"] - BVR_pelagic_DVM_raw[substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)=="epi", "Volume_unadj"])  
  
}
density.percent<- colnames(BVR_pelagic_DVM_vol_calculated[,c(6:11)])
for(i in 1:length(density.percent)){
  for(j in 1:length(unique(BVR.DVM.calcs$Hour))){
    BVR.DVM.calcs[j,paste0(density.percent,"_epi_percent_density")[i]]<- (BVR.DVM.calcs[j,paste0(density.percent,"_epi")][i]/ sum(BVR.DVM.calcs[j,paste0(density.percent,"_epi")[i]],BVR.DVM.calcs[j,paste0(density.percent,"_hypo")[i]])) *100
    BVR.DVM.calcs[j,paste0(density.percent,"_hypo_percent_density")[i]]<- (BVR.DVM.calcs[j,paste0(density.percent,"_hypo")][i]/ sum(BVR.DVM.calcs[j,paste0(density.percent,"_epi")[i]],BVR.DVM.calcs[j,paste0(density.percent,"_hypo")[i]])) * 100
  }       
}



#initialize df
BVR.DVM.calcs.SE<- data.frame("Hour"=unique(BVR_pelagic_DVM_raw$Hour))

#not sure if this is right, but calculating SE of difference between epi mean and hypo mean
SE.diffMean<- function(x,y){
  sqrt((sd(x,na.rm=TRUE)^2/length(na.omit(x))) + 
         (sd(y,na.rm=TRUE)^2/length(na.omit(y))))
}

#pull only noon/midnight samples
DVM_samples_raw <- zoop[(substrEnd(zoop$sample_ID,4)=="noon" | substrEnd(zoop$sample_ID,5)=="night" | substrEnd(zoop$sample_ID,8)=="noon_epi" | substrEnd(zoop$sample_ID,9)=="night_epi") & zoop$site_no!="BVR_l",]
matchingcols <- match(substr(colnames(BVR_pelagic_DVM_raw[1:20]),1,14),substr(colnames(DVM_samples_raw),1,14))
DVM_samples_raw<- DVM_samples_raw[,unique(matchingcols)]

DVM_samples_dens <- zoop[(substrEnd(zoop$sample_ID,4)=="noon" | substrEnd(zoop$sample_ID,5)=="night" | substrEnd(zoop$sample_ID,8)=="noon_epi" | substrEnd(zoop$sample_ID,9)=="night_epi") & zoop$site_no!="BVR_l",]
matchingcols <- match(substr(colnames(BVR_pelagic_DVM_vol_calculated [,c(1:4,6:11)]),1,14),substr(colnames(DVM_samples_dens),1,14))
DVM_samples_dens<- DVM_samples_dens[,unique(matchingcols)] 

#separate full vs epi samples
FullSamples <- DVM_samples_raw$sample_ID[c(3,5,9,15,18,22)]
EpiSamples<- DVM_samples_raw$sample_ID[c(1,11,7,13,17,20)]

#calculate hypo SE 
SEonly<- colnames(DVM_samples_raw)[7:20]
Percentdens <- colnames(DVM_samples_dens)[5:10]

for(i in 1:length(SEonly)){
  BVR.DVM.calcs.SE[,paste0(column.names,"_epi_SE")[i]] <- BVR_pelagic_DVM_vol_calculated[substrEnd(BVR_pelagic_DVM_vol_calculated$sample_ID,3)=="epi",substrEnd(colnames(BVR_pelagic_DVM_vol_calculated),2)=="SE"][i]
  BVR.DVM.calcs.SE[1,paste0(column.names,"_hypo_SE")[i]]<- SE.diffMean(DVM_samples_raw[DVM_samples_raw$sample_ID==FullSamples[1],paste0(SEonly)[i]],
                                                                       DVM_samples_raw[DVM_samples_raw$sample_ID==EpiSamples[1],paste0(SEonly)[i]])
  BVR.DVM.calcs.SE[2,paste0(column.names,"_hypo_SE")[i]]<- SE.diffMean(DVM_samples_raw[DVM_samples_raw$sample_ID==FullSamples[2],paste0(SEonly)[i]],
                                                                       DVM_samples_raw[DVM_samples_raw$sample_ID==EpiSamples[2],paste0(SEonly)[i]])
  BVR.DVM.calcs.SE[3,paste0(column.names,"_hypo_SE")[i]]<- SE.diffMean(DVM_samples_raw[DVM_samples_raw$sample_ID==FullSamples[3],paste0(SEonly)[i]],
                                                                       DVM_samples_raw[DVM_samples_raw$sample_ID==EpiSamples[3],paste0(SEonly)[i]])
} #Note: just repeating hypo SE calc for each date because I can't figure out the for loop across rows and columns...


for (j in 1:length(Percentdens)){
  BVR.DVM.calcs.SE[,paste0(Percentdens,"_epi_percent_density_SE")[j]]<- BVR_pelagic_DVM_vol_calculated[substrEnd(BVR_pelagic_DVM_vol_calculated$sample_ID,3)=="epi",substrEnd(colnames(BVR_pelagic_DVM_vol_calculated),11)=="NopL_rep.SE"][j]
  BVR.DVM.calcs.SE[1,paste0(Percentdens,"_hypo_percent_density_SE")[j]] <- SE.diffMean(DVM_samples_dens[DVM_samples_dens$sample_ID==FullSamples[1],paste0(Percentdens)[j]],
                                                                                       DVM_samples_dens[DVM_samples_dens$sample_ID==EpiSamples[1],paste0(Percentdens)[j]])
  BVR.DVM.calcs.SE[2,paste0(Percentdens,"_hypo_percent_density_SE")[j]] <- SE.diffMean(DVM_samples_dens[DVM_samples_dens$sample_ID==FullSamples[2],paste0(Percentdens)[j]],
                                                                                       DVM_samples_dens[DVM_samples_dens$sample_ID==EpiSamples[2],paste0(Percentdens)[j]])
  BVR.DVM.calcs.SE[3,paste0(Percentdens,"_hypo_percent_density_SE")[j]] <- SE.diffMean(DVM_samples_dens[DVM_samples_dens$sample_ID==FullSamples[3],paste0(Percentdens)[j]],
                                                                                       DVM_samples_dens[DVM_samples_dens$sample_ID==EpiSamples[3],paste0(Percentdens)[j]])
}


#wide to long for both dfs separately
BVR.DVM.calcs.long <-  BVR.DVM.calcs %>%
  gather(metric,value, ZoopDensity_No.pL_rep.mean_epi:nauplius_density_NopL_rep.mean_hypo_percent_density) %>%
  mutate(DateTime = strftime(Hour, "%m-%d-%Y %H:%M"))
BVR.DVM.calcs.SE.long <-  BVR.DVM.calcs.SE %>%
  gather(taxa.metric,SE,ZoopDensity_No.pL_rep.mean_epi_SE:nauplius_density_NopL_hypo_percent_density_SE)


#add the SE column from df2 to df1 for combined df
BVR.DVM.calcs.long$SE <- BVR.DVM.calcs.SE.long$SE

#add watercolumn and hour columns
BVR.DVM.calcs.long$WaterColumn <- ifelse(grepl("epi",BVR.DVM.calcs.long$metric),"epilimnion", ifelse(grepl("eta",BVR.DVM.calcs.long$metric),"metalimnion","hypolimnion"))
BVR.DVM.calcs.long$Hour <- ifelse(substrEnd(BVR.DVM.calcs.long$DateTime,5)=="12:00","noon","midnight")
BVR.DVM.calcs.long$Taxa <- substr(BVR.DVM.calcs.long$metric,1,9)

#shorten date time to just date
BVR.DVM.calcs.long$DateTime <- substr(BVR.DVM.calcs.long$DateTime,1,nchar(BVR.DVM.calcs.long$DateTime)-6)

#replace NAN with 0
BVR.DVM.calcs.long$value[is.nan(BVR.DVM.calcs.long$value)] <- 0

#export 2021 dvm stats
write.csv(BVR.DVM.calcs.long, "./Summer2021-DataAnalysis/SummaryStats/DVM_2021_zoops.csv")


#------------------------------------------------------------------------------#
#reorder water column level for bar order in figure
BVR.DVM.calcs.long$WaterColumn<-factor(BVR.DVM.calcs.long$WaterColumn,levels=c("epilimnion","metalimnion","hypolimnion"))
BVR.DVM.calcs.long<- BVR.DVM.calcs.long[!(BVR.DVM.calcs.long$metric=="ZoopDensity_No.pL_rep.mean_epi" | BVR.DVM.calcs.long$metric=="ZoopDensity_No.pL_rep.mean_meta"| BVR.DVM.calcs.long$metric=="ZoopDensity_No.pL_rep.mean_hypo" | BVR.DVM.calcs.long$metric=="BiomassConcentration_ugpL_rep.mean_epi"| BVR.DVM.calcs.long$metric=="BiomassConcentration_ugpL_rep.mean_meta" |BVR.DVM.calcs.long$metric=="BiomassConcentration_ugpL_rep.mean_hypo" | substrEnd(BVR.DVM.calcs.long$metric,15)=="percent_density"),] 
BVR.DVM.calcs.long$Taxa <-factor(BVR.DVM.calcs.long$Taxa, levels = c("Cladocera","Rotifera_","Cyclopoid","Calanoida"))


#jpeg("Figures/BVR_epimetahypo_density_16Jun2021.jpg", width = 6, height = 4, units = "in",res = 300)
ggplot(subset(BVR.DVM.calcs.long, grepl("density",metric,ignore.case = TRUE)), aes(x=WaterColumn, y=value)) +
  geom_bar(aes(fill=Taxa, alpha=WaterColumn), stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() +
  facet_wrap(~Taxa, scales= 'free',ncol=5, strip.position = "right", labeller=as_labeller(function(variable,value){c("","","","","Noon")})) +
  theme(strip.text.y = element_text(size = 10,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank(), plot.margin = margin(2.8,0,2,0,unit = "cm")) + 
  scale_x_discrete(name="", labels=rep(c("epi","meta","hypo"),10)) + labs(title="16 Jun 2021 Noon") + 
  scale_fill_manual(values=viridis(6),labels=c("Cladocera", "Rotifera", "Cyclopoida","Calanoida")) +
  geom_hline(yintercept=0) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
                                  axis.text.x=element_text(size=7,family="Times"), plot.title = element_text(hjust = 0.5)) + ylab("Density (individual/L)") +
  theme(legend.position = "bottom", legend.margin = margin(0, 1, 2, 1), legend.title = element_text(size=10),legend.key.size = unit(0.3,"cm"),
        legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),legend.text = element_text(size=8),axis.text.y = element_text(size=10, family="Times"))
#dev.off()


density.percent<- colnames(BVR_pelagic_DVM_vol_calculated[,c(6:9)])
for(i in 1:length(density.percent)){
for(j in 1:length(unique(BVR.DVM.calcs$Hour))){
    BVR.DVM.calcs[j,paste0(density.percent,"_epi_percent_density")[i]]<- (BVR.DVM.calcs[j,paste0(density.percent,"_epi")][i]/ sum(BVR.DVM.calcs[j,paste0(density.percent,"_epi")[i]],BVR.DVM.calcs[j,paste0(density.percent,"_hypo")[i]])) *100
    BVR.DVM.calcs[j,paste0(density.percent,"_hypo_percent_density")[i]]<- (BVR.DVM.calcs[j,paste0(density.percent,"_hypo")][i]/ sum(BVR.DVM.calcs[j,paste0(density.percent,"_epi")[i]],BVR.DVM.calcs[j,paste0(density.percent,"_hypo")[i]])) * 100

}       
}

#initialize df
BVR.DVM.calcs.SE<- data.frame("Hour"=unique(BVR_pelagic_DVM_raw$Hour))

#not sure if this is right, but calculating SE of difference between epi mean and hypo mean
SE.diffMean<- function(x,y){
  sqrt((sd(x,na.rm=TRUE)^2/length(na.omit(x))) + 
         (sd(y,na.rm=TRUE)^2/length(na.omit(y))))
}

#pull only noon/midnight samples
DVM_samples_raw <- zoop[(substrEnd(zoop$sample_ID,4)=="noon" | substrEnd(zoop$sample_ID,5)=="night" | substrEnd(zoop$sample_ID,8)=="noon_epi" | substrEnd(zoop$sample_ID,9)=="night_epi") & zoop$site_no!="BVR_l",]
matchingcols <- match(substr(colnames(BVR_pelagic_DVM_raw[1:16]),1,14),substr(colnames(DVM_samples_raw),1,14))
DVM_samples_raw<- DVM_samples_raw[,unique(matchingcols)]
                        
DVM_samples_dens <- zoop[(substrEnd(zoop$sample_ID,4)=="noon" | substrEnd(zoop$sample_ID,5)=="night" | substrEnd(zoop$sample_ID,8)=="noon_epi" | substrEnd(zoop$sample_ID,9)=="night_epi") & zoop$site_no!="BVR_l",]
matchingcols <- match(substr(colnames(BVR_pelagic_DVM_vol_calculated [,c(1:4,6:9)]),1,14),substr(colnames(DVM_samples_dens),1,14))
DVM_samples_dens<- DVM_samples_dens[,unique(matchingcols)] 

#separate full vs epi samples
FullSamples <- unique(DVM_samples_raw$sample_ID)[1:3]
EpiSamples<- unique(DVM_samples_raw$sample_ID)[4:6]

#calculate hypo SE 
SEonly<- colnames(DVM_samples_raw)[7:16]
Percentdens <- colnames(DVM_samples_dens)[5:8]

for(i in 1:length(SEonly)){
  BVR.DVM.calcs.SE[,paste0(column.names,"_epi_SE")[i]] <- BVR_pelagic_DVM_vol_calculated[substrEnd(BVR_pelagic_DVM_vol_calculated$sample_ID,3)=="epi",substrEnd(colnames(BVR_pelagic_DVM_vol_calculated),2)=="SE"][i]
  BVR.DVM.calcs.SE[1,paste0(column.names,"_hypo_SE")[i]]<- SE.diffMean(DVM_samples_raw[DVM_samples_raw$sample_ID==FullSamples[1],paste0(SEonly)[i]],
                                                                      DVM_samples_raw[DVM_samples_raw$sample_ID==EpiSamples[1],paste0(SEonly)[i]])
  BVR.DVM.calcs.SE[2,paste0(column.names,"_hypo_SE")[i]]<- SE.diffMean(DVM_samples_raw[DVM_samples_raw$sample_ID==FullSamples[2],paste0(SEonly)[i]],
                                                                       DVM_samples_raw[DVM_samples_raw$sample_ID==EpiSamples[2],paste0(SEonly)[i]])
  BVR.DVM.calcs.SE[3,paste0(column.names,"_hypo_SE")[i]]<- SE.diffMean(DVM_samples_raw[DVM_samples_raw$sample_ID==FullSamples[3],paste0(SEonly)[i]],
                                                                       DVM_samples_raw[DVM_samples_raw$sample_ID==EpiSamples[3],paste0(SEonly)[i]])
  } #Note: just repeating hypo SE calc for each date because I can't figure out the for loop across rows and columns...


  for (j in 1:length(Percentdens)){
  BVR.DVM.calcs.SE[,paste0(Percentdens,"_epi_percent_density_SE")[j]]<- BVR_pelagic_DVM_vol_calculated[substrEnd(BVR_pelagic_DVM_vol_calculated$sample_ID,3)=="epi",substrEnd(colnames(BVR_pelagic_DVM_vol_calculated),11)=="NopL_rep.SE"][j]
  BVR.DVM.calcs.SE[1,paste0(Percentdens,"_hypo_percent_density_SE")[j]] <- SE.diffMean(DVM_samples_dens[DVM_samples_dens$sample_ID==FullSamples[1],paste0(Percentdens)[j]],
                                                                                      DVM_samples_dens[DVM_samples_dens$sample_ID==EpiSamples[1],paste0(Percentdens)[j]])
  BVR.DVM.calcs.SE[2,paste0(Percentdens,"_hypo_percent_density_SE")[j]] <- SE.diffMean(DVM_samples_dens[DVM_samples_dens$sample_ID==FullSamples[2],paste0(Percentdens)[j]],
                                                                                       DVM_samples_dens[DVM_samples_dens$sample_ID==EpiSamples[2],paste0(Percentdens)[j]])
  BVR.DVM.calcs.SE[3,paste0(Percentdens,"_hypo_percent_density_SE")[j]] <- SE.diffMean(DVM_samples_dens[DVM_samples_dens$sample_ID==FullSamples[3],paste0(Percentdens)[j]],
                                                                                          DVM_samples_dens[DVM_samples_dens$sample_ID==EpiSamples[3],paste0(Percentdens)[j]])
  }

#wide to long for both dfs separately
BVR.DVM.calcs.long <-  BVR.DVM.calcs %>%
    gather(metric,value, ZoopDensity_No.pL_rep.mean_epi:Calanoida_density_NopL_rep.mean_hypo_percent_density) %>%
    mutate(DateTime = strftime(Hour, "%m-%d-%Y %H:%M"))
BVR.DVM.calcs.SE.long <-  BVR.DVM.calcs.SE %>%
gather(taxa.metric,SE,ZoopDensity_No.pL_rep.mean_epi_SE:Calanoida_density_NopL_hypo_percent_density_SE)

#add the SE column from df2 to df1 for combined df
BVR.DVM.calcs.long$SE <- BVR.DVM.calcs.SE.long$SE

#add watercolumn and hour columns
BVR.DVM.calcs.long$WaterColumn <- ifelse((substrEnd(BVR.DVM.calcs.long$metric,3)=="epi" |substrEnd(BVR.DVM.calcs.long$metric,19)=="epi_percent_density") ,"epilimnion","hypolimnion")
BVR.DVM.calcs.long$Hour <- ifelse(substrEnd(BVR.DVM.calcs.long$DateTime,5)=="12:00","noon","midnight")
BVR.DVM.calcs.long$Taxa <- substr(BVR.DVM.calcs.long$metric,1,9)

#shorten date time to just date
BVR.DVM.calcs.long$DateTime <- substr(BVR.DVM.calcs.long$DateTime,1,nchar(BVR.DVM.calcs.long$DateTime)-6)

#ridiculous way to work with these stupid facet labels
facet_labeller_bot <- function(variable, value) {
  rep("",8)
}

facet_labeller_top <- function(variable, value) {
  c("","","","Midnight","","","","Noon")
}

#-------------------------------------------------------#
#        Create new df with more diverse taxa           # 
# subsetting repmeans into new df for DVM analyses/figs #
#-------------------------------------------------------#

#sum up counts by sample/site/day for DVM analyses + figs
BVR_counts_moretaxa <- zoop %>% select(sample_ID,site_no,collect_date,Hour, DaphniaCount_n,
    CeriodaphniaCount_n,BosminidaeCount_n,CalanoidaCount_n,CyclopoidaCount_n,
    naupliusCount_n,KeratellaCount_n,KellicottiaCount_n,CollothecidaeCount_n,
    ConochilidaeCount_n,SynchaetidaeCount_n,TrichocercidaeCount_n) %>%
    group_by(sample_ID, site_no, Hour, collect_date) %>%
    summarise_at(vars(DaphniaCount_n:TrichocercidaeCount_n), funs(rep.mean=mean, rep.SE=stderr))

#add unadjusted volume
BVR_counts_moretaxa$Volume_unadj<- (zoop %>% select(sample_ID,site_no,collect_date,Hour, Volume_unadj) %>%
   group_by(sample_ID, site_no, Hour, collect_date) %>%
   summarise_at(vars(Volume_unadj),funs(rep.mean=mean)))$rep.mean

#add proportional volume (aliquot vol/ total diluted sample volume)
BVR_counts_moretaxa$proportional_vol <- (zoop %>% select(sample_ID,site_no,collect_date,Hour, proportional_vol) %>%
  group_by(sample_ID, site_no, Hour, collect_date) %>%
  summarise_at(vars(proportional_vol),funs(rep.mean=mean)))$rep.mean

#get BVR_counts df in same order as zoop.repmeans df
BVR_counts_moretaxa<- BVR_counts_moretaxa[order(match(paste0(BVR_counts_moretaxa$sample_ID,BVR_counts_moretaxa$site_no,BVR_counts_moretaxa$collect_date), 
   paste0(zoop.repmeans.bytaxa$sample_ID,zoop.repmeans.bytaxa$site_no,zoop.repmeans.bytaxa$collect_date))),]

#add counts and vol and SE to zoop.repmeans
zoop.repmeans.bytaxa[,paste0(colnames(BVR_counts_moretaxa[5:30]))]<- BVR_counts_moretaxa[5:30]

#new dfs for DVM data (BVR_pelagic_DVM_raw is just raw #/ug; BVR_pelagic_DVM_vol_calculated is #/L and ug/L)
BVR_pelagic_taxa_DVM<- zoop.repmeans.bytaxa[(zoop.repmeans.bytaxa$site_no=="BVR_50" |zoop.repmeans.bytaxa$site_no=="BVR_50_p") &
     (substrEnd(zoop.repmeans.bytaxa$sample_ID,5)=="night" | substrEnd(zoop.repmeans.bytaxa$sample_ID,4)=="noon" |
      substrEnd(zoop.repmeans.bytaxa$sample_ID,9)=="night_epi" | substrEnd(zoop.repmeans.bytaxa$sample_ID,8)=="noon_epi" |
      substrEnd(zoop.repmeans.bytaxa$sample_ID,5)=="cline"),] 

#only select volume, count, and ug cols (plus SE cols)
BVR_pelagic_taxa_DVM_raw<- BVR_pelagic_taxa_DVM[,c(1:4,101,102,which(substrEnd(colnames(BVR_pelagic_taxa_DVM),10)=="n_rep.mean"),
     which(substrEnd(colnames(BVR_pelagic_taxa_DVM),11)=="ug_rep.mean"),
     which(substrEnd(colnames(BVR_pelagic_taxa_DVM),8)=="n_rep.SE"),
     which(substrEnd(colnames(BVR_pelagic_taxa_DVM),9)=="ug_rep.SE"))]

BVR_pelagic_taxa_DVM_vol_calculated <- BVR_pelagic_taxa_DVM[,c(1:4,which(substrEnd(colnames(BVR_pelagic_taxa_DVM),15)=="y_NopL_rep.mean"),
     which(substrEnd(colnames(BVR_pelagic_taxa_DVM),15)=="n_ugpL_rep.mean"),
     which(substrEnd(colnames(BVR_pelagic_taxa_DVM),13)=="n_ugpL_rep.SE"),
     which(substrEnd(colnames(BVR_pelagic_taxa_DVM),13)=="y_NopL_rep.SE"))]

#initialize df
BVR.DVM.taxa.calcs<- data.frame("Hour"=unique(BVR_pelagic_taxa_DVM_raw$Hour))

#for loop to fill out epi vs hypo calcs (using 4m as oxycline for MSN 1, but this should be changed to oxycline sample once available!)
column.names<- colnames(BVR_pelagic_taxa_DVM_vol_calculated[,c(5:28)])
variables<- colnames(BVR_pelagic_taxa_DVM_raw[,c(7:30)])
for(i in 1:length(variables)){
  BVR.DVM.taxa.calcs[,paste0(column.names,"_epi")[i]]<- BVR_pelagic_taxa_DVM_vol_calculated[substrEnd(BVR_pelagic_taxa_DVM_vol_calculated$sample_ID,3)=="epi",paste0(column.names)[i]]
  BVR.DVM.taxa.calcs[,paste0(column.names,"_hypo")[i]] <- (((1/BVR_pelagic_taxa_DVM_raw[substrEnd(BVR_pelagic_taxa_DVM_raw$sample_ID,3)!="epi" , "proportional_vol"]) * BVR_pelagic_taxa_DVM_raw[substrEnd(BVR_pelagic_taxa_DVM_raw$sample_ID,3)!="epi" ,paste0(variables)[i]]) - 
                                                             ((1/BVR_pelagic_taxa_DVM_raw[substrEnd(BVR_pelagic_taxa_DVM_raw$sample_ID,3)=="epi", "proportional_vol"]) * BVR_pelagic_taxa_DVM_raw[substrEnd(BVR_pelagic_taxa_DVM_raw$sample_ID,3)=="epi",paste0(variables)[i]])) /
    (BVR_pelagic_taxa_DVM_raw[substrEnd(BVR_pelagic_taxa_DVM_raw$sample_ID,3)!="epi" ,"Volume_unadj"] - BVR_pelagic_taxa_DVM_raw[substrEnd(BVR_pelagic_taxa_DVM_raw$sample_ID,3)=="epi", "Volume_unadj"])  
}

#initialize df
BVR.DVM.taxa.calcs.SE<- data.frame("Hour"=unique(BVR_pelagic_taxa_DVM_raw$Hour))

#pull only noon/midnight samples
DVM_taxa_samples_raw <- zoop[(substrEnd(zoop$sample_ID,4)=="noon" | substrEnd(zoop$sample_ID,5)=="night" | substrEnd(zoop$sample_ID,8)=="noon_epi" | substrEnd(zoop$sample_ID,9)=="night_epi") & zoop$site_no!="BVR_l",]
matchingcols <- match(substr(colnames(BVR_pelagic_taxa_DVM_raw[1:30]),1,15),substr(colnames(DVM_taxa_samples_raw),1,15))
#manually replace NA with 117 (DaphniaCount_n is too short) and manually choose Trichocercidae_totalbiomass_ug (wants to pull percent because string too long; #175)
matchingcols[7]<-117
matchingcols[30]<- 175
DVM_taxa_samples_raw<- DVM_taxa_samples_raw[,unique(matchingcols)]

#calculate epi and hypo SE 
SEonly<- colnames(DVM_taxa_samples_raw)[7:30]
for(i in 1:length(SEonly)){
  BVR.DVM.taxa.calcs.SE[,paste0(column.names,"_epi_SE")[i]] <- BVR_pelagic_taxa_DVM_vol_calculated[substrEnd(BVR_pelagic_taxa_DVM_vol_calculated$sample_ID,3)=="epi",substrEnd(colnames(BVR_pelagic_taxa_DVM_vol_calculated),2)=="SE"][i]
  BVR.DVM.taxa.calcs.SE[1,paste0(column.names,"_hypo_SE")[i]]<- SE.diffMean(DVM_taxa_samples_raw[DVM_taxa_samples_raw$sample_ID==FullSamples[1],paste0(SEonly)[i]],
                                                                       DVM_taxa_samples_raw[DVM_taxa_samples_raw$sample_ID==EpiSamples[1],paste0(SEonly)[i]])
  BVR.DVM.taxa.calcs.SE[2,paste0(column.names,"_hypo_SE")[i]]<- SE.diffMean(DVM_taxa_samples_raw[DVM_taxa_samples_raw$sample_ID==FullSamples[2],paste0(SEonly)[i]],
                                                                       DVM_taxa_samples_raw[DVM_taxa_samples_raw$sample_ID==EpiSamples[2],paste0(SEonly)[i]])
  BVR.DVM.taxa.calcs.SE[3,paste0(column.names,"_hypo_SE")[i]]<- SE.diffMean(DVM_taxa_samples_raw[DVM_taxa_samples_raw$sample_ID==FullSamples[3],paste0(SEonly)[i]],
                                                                       DVM_taxa_samples_raw[DVM_taxa_samples_raw$sample_ID==EpiSamples[3],paste0(SEonly)[i]])
} #Note: just repeating hypo SE calc for each date because I can't figure out the for loop across rows and columns...

#wide to long for both dfs separately
BVR.DVM.taxa.calcs.long <-  BVR.DVM.taxa.calcs %>%
  gather(metric,value, Daphnia_density_NopL_rep.mean_epi:Trichocercidae_BiomassConcentration_ugpL_rep.mean_hypo) %>%
  mutate(DateTime = strftime(Hour, "%m-%d-%Y %H:%M"))
BVR.DVM.taxa.calcs.SE.long <-  BVR.DVM.taxa.calcs.SE %>%
  gather(taxa.metric,SE,Daphnia_density_NopL_rep.mean_epi_SE:Trichocercidae_BiomassConcentration_ugpL_rep.mean_hypo_SE) 

#add the SE column from df2 to df1 for combined df
BVR.DVM.taxa.calcs.long$SE <- BVR.DVM.taxa.calcs.SE.long$SE

#add watercolumn and hour columns
BVR.DVM.taxa.calcs.long$WaterColumn <- ifelse(substrEnd(BVR.DVM.taxa.calcs.long$metric,3)=="epi","epilimnion","hypolimnion")
BVR.DVM.taxa.calcs.long$Hour <- ifelse(substrEnd(BVR.DVM.taxa.calcs.long$DateTime,5)=="12:00","noon","midnight")
BVR.DVM.taxa.calcs.long$Taxa <- substr(BVR.DVM.taxa.calcs.long$metric,1,9)

#shorten date time to just date
BVR.DVM.taxa.calcs.long$DateTime <- substr(BVR.DVM.taxa.calcs.long$DateTime,1,nchar(BVR.DVM.taxa.calcs.long$DateTime)-6)

#-----------------------#
#  Epi vs Hypo figures  #
#-----------------------#

#df for percent density calcs
BVR.DVM.percent.dens.calcs<- BVR.DVM.calcs.long[substrEnd(BVR.DVM.calcs.long$metric,15)=="percent_density",]

#drop total density and biomass rows
BVR.DVM.calcs.long<- BVR.DVM.calcs.long[!(BVR.DVM.calcs.long$metric=="ZoopDensity_No.pL_rep.mean_epi" | BVR.DVM.calcs.long$metric=="ZoopDensity_No.pL_rep.mean_hypo" | BVR.DVM.calcs.long$metric=="BiomassConcentration_ugpL_rep.mean_epi" |BVR.DVM.calcs.long$metric=="BiomassConcentration_ugpL_rep.mean_hypo" | substrEnd(BVR.DVM.calcs.long$metric,15)=="percent_density"),] 
#reorder 
BVR.DVM.calcs.long$Taxa <- factor(BVR.DVM.calcs.long$Taxa, levels = c("Cladocera","Rotifera_","Cyclopoid","Calanoida"))

#Epi vs hypo taxa density for 12-13 Aug 2020
#jpeg("Figures/BVR_epivshypo_density_12Aug2020.jpg", width = 6, height = 4, units = "in",res = 300)
ggplot(subset(BVR.DVM.calcs.long, grepl("density",metric,ignore.case = TRUE) & DateTime =="08-12-2020"), aes(x=WaterColumn, y=value)) +
    geom_bar(aes(fill=Taxa, alpha=WaterColumn), stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() +
    geom_errorbar(aes(ymin=value-SE, ymax=value+SE), width=.2,position=position_dodge(.9)) + scale_alpha_manual(values = c(0.4, 1)) +
    facet_wrap(~Taxa, scales= 'free_y',ncol=4, strip.position = "right", labeller=as_labeller(function(variable,value){c("","","","Noon")})) +
    theme(strip.text.y = element_text(size = 10,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank(), plot.margin = margin(2.8,0,2,0,unit = "cm")) + 
    scale_x_discrete(name="", labels=rep(c("epi","hypo"),10)) + labs(title="12 Aug 2020 Noon") + 
    scale_fill_manual(values=viridis(6),labels=c("Cladocera","Rotifera","Cyclopoida","Calanoida"),drop=TRUE) +
    geom_hline(yintercept=0) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
       axis.text.x=element_text(size=9,family="Times"), plot.title = element_text(hjust = 0.5)) + ylab("Density (individual/L)") +
    theme(legend.position = "bottom", legend.margin = margin(0, 1, 2, 1), legend.title = element_text(size=10),legend.key.size = unit(0.3,"cm"),
        legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),legend.text = element_text(size=8),axis.text.y = element_text(size=10, family="Times"))
#dev.off()

#jpeg("Figures/BVR_epivshypo_biomass_12Aug2020.jpg", width = 6, height = 4, units = "in",res = 300)
ggplot(subset(BVR.DVM.calcs.long, grepl("biomass",metric,ignore.case = TRUE) & DateTime =="08-12-2020"), aes(x=WaterColumn, y=value)) +
  geom_bar(aes(fill=Taxa, alpha=WaterColumn), stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() +
  geom_errorbar(aes(ymin=value-SE, ymax=value+SE), width=.2,position=position_dodge(.9)) + scale_alpha_manual(values = c(0.4, 1)) +
  facet_wrap(~Taxa, scales= 'free',ncol=4, strip.position = "right", labeller=as_labeller(function(variable,value){c("","","","Noon")})) +
  theme(strip.text.y = element_text(size = 10,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank(), plot.margin = margin(2.8,0,2,0,unit = "cm")) + 
  scale_x_discrete(name="", labels=rep(c("epi","hypo"),10)) + labs(title="12 Aug 2020 Noon") + 
  scale_fill_manual(values=viridis(6),labels=c("Cladocera","Rotifera","Cyclopoida","Calanoida")) +
  geom_hline(yintercept=0) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
     axis.text.x=element_text(size=9,family="Times"), plot.title = element_text(hjust = 0.5)) + ylab(expression(paste("Biomass (",mu,"g/L)"))) +
  theme(legend.position = "bottom", legend.margin = margin(0, 1, 2, 1), legend.title = element_text(size=10),legend.key.size = unit(0.3,"cm"),
        legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),legend.text = element_text(size=8),axis.text.y = element_text(size=10, family="Times"))
#dev.off()

  #jpeg("Figures/BVR_epivshypo_density_13Aug2020.jpg", width = 6, height = 4, units = "in",res = 300)
ggplot(subset(BVR.DVM.calcs.long, grepl("density",metric,ignore.case = TRUE) & DateTime=="08-13-2020"), aes(x=WaterColumn, y=value)) +
    geom_rect(data=subset(BVR.DVM.calcs.long,Hour == 'midnight' &grepl("density",metric,ignore.case = TRUE)),
         aes(fill=Hour),xmin=-Inf ,xmax = Inf, ymin = -Inf, ymax = Inf, fill = 'black', alpha = 0.053,inherit.aes = FALSE) +
    geom_bar(aes(fill=Taxa, alpha=WaterColumn), stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() +
    geom_errorbar(aes(ymin=value-SE, ymax=value+SE), width=.2,position=position_dodge(.9)) + scale_alpha_manual(values = c(0.4, 1)) +
    facet_wrap(Hour~Taxa, scales= 'free',ncol=4, strip.position = "right", labeller=labeller(Hour=as_labeller(facet_labeller_bot),Taxa=as_labeller(facet_labeller_top))) +
    theme(strip.text.y = element_text(size = 11 ,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank()) + 
    scale_x_discrete(name="", labels=rep(c("epi","hypo"),10)) + labs(title="13 Aug 2020") + guides(alpha=FALSE) +
    scale_fill_manual("Taxa",values=viridis(6),labels=c("Cladocera","Rotifera","Cyclopoida","Calanoida")) +
    geom_hline(yintercept=0) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
         axis.text.x=element_text(size=6,family="Times"), plot.title = element_text(hjust = 0.5)) + ylab("Density (individual/L)") +
    theme(legend.position = "bottom", legend.margin = margin(0, 1, 2, 1),plot.margin = margin(0,0,0,0.3,unit = "cm"), axis.text.y = element_text(size=13, family="Times"),
          legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),legend.text = element_text(size=10),legend.title = element_text(size=12))
   #dev.off()

 #Epi vs hypo taxa biomass 
#jpeg("Figures/BVR_epivshypo_biomass_13Aug2020.jpg", width = 6, height = 4, units = "in",res = 300)
ggplot(subset(BVR.DVM.calcs.long, grepl("biomass",metric,ignore.case = TRUE) & DateTime=="08-13-2020"), aes(x=WaterColumn, y=value)) +
   geom_rect(data=subset(BVR.DVM.calcs.long,Hour == 'midnight' &grepl("biomass",metric,ignore.case = TRUE)),
           aes(fill=Hour),xmin=-Inf ,xmax = Inf, ymin = -Inf, ymax = Inf, fill = 'black', alpha = 0.053,inherit.aes = FALSE) +
   geom_bar(aes(fill=Taxa, alpha=WaterColumn), stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() +
   geom_errorbar(aes(ymin=value-SE, ymax=value+SE), width=.2,position=position_dodge(.9)) + scale_alpha_manual(values = c(0.4, 1)) +
   facet_wrap(Hour~Taxa, scales= 'free',ncol=4, strip.position = "right", labeller=labeller(Hour=as_labeller(facet_labeller_bot),Taxa=as_labeller(facet_labeller_top))) +
   theme(strip.text.y = element_text(size = 11 ,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank()) + 
   scale_x_discrete(name="", labels=rep(c("epi","hypo"),10)) + labs(title="13 Aug 2020") +  guides(alpha=FALSE) +
   scale_fill_manual("Taxa",values=viridis(6),labels=c("Cladocera","Rotifera","Cyclopoida","Calanoida")) +
   geom_hline(yintercept=0) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
            axis.text.x=element_text(size=6,family="Times"), plot.title = element_text(hjust = 0.5)) + ylab(expression(paste("Biomass (",mu,"g/L)"))) +
      theme(legend.position = "bottom", legend.margin = margin(0, 1, 2, 1),plot.margin = margin(0,0,0,0.3,unit = "cm"), axis.text.y = element_text(size=13, family="Times"),
            legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),legend.text = element_text(size=10),legend.title = element_text(size=12))
  #dev.off()
  
 
#reorder taxa
BVR.DVM.percent.dens.calcs$Taxa <- factor(BVR.DVM.percent.dens.calcs$Taxa, levels = c("Cladocera","Rotifera_","Cyclopoid","Calanoida"))

####Epi vs hypo taxa percent density figs ####
#jpeg("Figures/BVR_epivshypo_percent_density_12Aug2020.jpg", width = 6, height = 4, units = "in",res = 300)
ggplot(subset(BVR.DVM.percent.dens.calcs,  DateTime=="08-12-2020"), aes(x=WaterColumn, y=value)) +
  geom_bar(aes(fill=Taxa, alpha=WaterColumn), stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() +
  geom_errorbar(aes(ymin=value-SE, ymax=value+SE), width=.2,position=position_dodge(.9)) + scale_alpha_manual(values = c(0.4, 1)) +
  facet_wrap(~Taxa, scales= 'free',ncol=4, strip.position = "right", labeller=as_labeller(function(variable,value){c("","","","noon")})) +
  theme(strip.text.y = element_text(size = 13 ,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank()) + 
  scale_x_discrete(name="", labels=rep(c("epi","hypo"),10)) + labs(title="12 Aug 2020") + guides(alpha=FALSE) + 
  scale_fill_manual("Taxa",values=viridis(6),labels=c("Cladocera","Rotifera","Cyclopoida","Calanoida")) +
  geom_hline(yintercept=0) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
                                  axis.text.x=element_text(size=6,family="Times"), plot.title = element_text(hjust = 0.5)) + ylab("Percent Density") +
  theme(legend.position = "bottom", legend.margin = margin(1, 1, 2, 1),plot.margin = margin(0,0,0,0.3,unit = "cm"), axis.text.y = element_text(size=13, family="Times"), axis.title.y = element_text(size=14, family="Times"), 
        legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),legend.text = element_text(size=10),legend.title = element_text(size=12))
#dev.off()

  #jpeg("Figures/BVR_epivshypo_percent_density_13Aug2020.jpg", width = 6, height = 4, units = "in",res = 300)
ggplot(subset(BVR.DVM.percent.dens.calcs,  DateTime=="08-13-2020"), aes(x=WaterColumn, y=value)) +
    geom_rect(data=subset(BVR.DVM.percent.dens.calcs,Hour == 'midnight' &grepl("density",metric,ignore.case = TRUE)),
         aes(fill=Hour),xmin=-Inf ,xmax = Inf, ymin = -Inf, ymax = Inf, fill = 'black', alpha = 0.053,inherit.aes = FALSE) +
    geom_bar(aes(fill=Taxa, alpha=WaterColumn), stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() +
    geom_errorbar(aes(ymin=value-SE, ymax=value+SE), width=.2,position=position_dodge(.9)) + scale_alpha_manual(values = c(0.4, 1)) +
    facet_wrap(Hour~Taxa, scales= 'free',ncol=4, strip.position = "right", labeller=labeller(Hour=as_labeller(facet_labeller_bot),Taxa=as_labeller(facet_labeller_top))) +
    theme(strip.text.y = element_text(size = 13 ,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank()) + 
    scale_x_discrete(name="", labels=rep(c("epi","hypo"),10)) + labs(title="13 Aug 2020") + guides(alpha=FALSE) + 
    scale_fill_manual("Taxa",values=viridis(6),labels=c("Cladocera","Rotifera","Cyclopoida","Calanoida")) +
    geom_hline(yintercept=0) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
         axis.text.x=element_text(size=6,family="Times"), plot.title = element_text(hjust = 0.5)) + ylab("Percent Density") +
    theme(legend.position = "bottom", legend.margin = margin(1, 1, 2, 1),plot.margin = margin(0,0,0,0.3,unit = "cm"), axis.text.y = element_text(size=13, family="Times"), axis.title.y = element_text(size=14, family="Times"), 
         legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),legend.text = element_text(size=10),legend.title = element_text(size=12))
  #dev.off()


###################################################
# more diverse taxa figs

#order taxa to get bars in correct order
BVR.DVM.taxa.calcs.long$Taxa <- factor(BVR.DVM.taxa.calcs.long$Taxa, levels= c("Daphnia_d","Daphnia_B","Ceriodaph","Bosminida",
             "Cyclopoid","Calanoida","nauplius_","Keratella","Kellicott","Collothec","Conochili","Synchaeti","Trichocer"))

#jpeg("Figures/BVR_moretaxa_epivshypo_density_12Aug2020.jpg", width = 6, height = 4, units = "in",res = 300)
ggplot(subset(BVR.DVM.taxa.calcs.long, grepl("density",metric,ignore.case = TRUE) & DateTime =="08-12-2020"), aes(x=WaterColumn, y=value)) +
  geom_bar(aes(fill=Taxa, alpha=WaterColumn), stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() +
  geom_errorbar(aes(ymin=value-SE, ymax=value+SE), width=.2,position=position_dodge(.9)) + scale_alpha_manual(values = c(0.4, 1)) +
  facet_wrap(~Taxa, scales= 'free',ncol=6, strip.position = "right", labeller=as_labeller(function(variable,value){c("","","","","","")}),drop=TRUE) +
  theme(strip.text.y = element_text(size = 10,margin = margin(0, -0.03, 0,0.1, "cm")),strip.background = element_blank(), plot.margin = margin(1,0,1,0,unit = "cm")) + 
  scale_x_discrete(name="", labels=rep(c("epi","hypo"),10)) + labs(title="12 Aug 2020 Noon") + 
  scale_fill_manual(values=rainbow(12),labels=taxa2) +
  geom_hline(yintercept=0) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text.x=element_text(size=8,family="Times"), plot.title = element_text(hjust = 0.5)) + ylab("Density (individual/L)") +
  theme(legend.position = "bottom", legend.margin = margin(0, 1, 2, 1), legend.title = element_text(size=10),legend.key.size = unit(0.3,"cm"),
        legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),legend.text = element_text(size=8),axis.text.y = element_text(size=10, family="Times"))
#dev.off()

#jpeg("Figures/BVR_moretaxa_epivshypo_biomass_12Aug2020.jpg", width = 6, height = 4, units = "in",res = 300)
ggplot(subset(BVR.DVM.taxa.calcs.long, grepl("biomass",metric,ignore.case = TRUE) & DateTime =="08-12-2020"), aes(x=WaterColumn, y=value)) +
  geom_bar(aes(fill=Taxa, alpha=WaterColumn), stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() +
  geom_errorbar(aes(ymin=value-SE, ymax=value+SE), width=.2,position=position_dodge(.9)) + scale_alpha_manual(values = c(0.4, 1)) +
  facet_wrap(~Taxa, scales= 'free',ncol=6, strip.position = "right", labeller=as_labeller(function(variable,value){c("","","","","","")})) +
  theme(strip.text.y = element_text(size = 10,margin = margin(0, -0.03, 0,0.1, "cm")),strip.background = element_blank(), plot.margin = margin(1,0,1,0,unit = "cm")) + 
  scale_x_discrete(name="", labels=rep(c("epi","hypo"),10)) + labs(title="12 Aug 2020 Noon") + 
  scale_fill_manual(values=rainbow(12),labels=taxa2) +
  geom_hline(yintercept=0) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text.x=element_text(size=7,family="Times"), plot.title = element_text(hjust = 0.5)) + ylab(expression(paste("Biomass (",mu,"g/L)"))) +
  theme(legend.position = "bottom", legend.margin = margin(0, 1, 2, 1), legend.title = element_text(size=10),legend.key.size = unit(0.3,"cm"),
        legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),legend.text = element_text(size=8),axis.text.y = element_text(size=10, family="Times"))
#dev.off()

#ridiculous way to work with these stupid facet labels
facet_labeller_bot <- function(variable, value) {
  rep("",24)
}

facet_labeller_top <- function(variable, value) {
  c("","","","","","Midnight","","","","","","Midnight","","","","","","Noon","","","","","","Noon")
}

#Epi vs hypo taxa density and biomass 13 Aug2020
  #jpeg("Figures/BVR_epivshypo_density_13Aug2020.jpg", width = 6, height = 4, units = "in",res = 300)
ggplot(subset(BVR.DVM.taxa.calcs.long, grepl("density",metric,ignore.case = TRUE) & DateTime=="08-13-2020"), aes(x=WaterColumn, y=value)) +
    geom_rect(data=subset(BVR.DVM.taxa.calcs.long,Hour == 'midnight' &grepl("density",metric,ignore.case = TRUE)),
           aes(fill=Hour),xmin=-Inf ,xmax = Inf, ymin = -Inf, ymax = Inf, fill = 'black', alpha = 0.053,inherit.aes = FALSE) +
    geom_bar(aes(fill=Taxa, alpha=WaterColumn), stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() +
    geom_errorbar(aes(ymin=value-SE, ymax=value+SE), width=.2,position=position_dodge(.9)) + scale_alpha_manual(values = c(0.4, 1)) +
    facet_wrap(Hour~Taxa, scales= 'free',ncol=6, strip.position = "right", labeller=labeller(Hour=as_labeller(facet_labeller_bot),Taxa=as_labeller(facet_labeller_top))) +
    theme(strip.text.y = element_text(size = 0 ,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank()) + 
    scale_x_discrete(name="", labels=rep(c("epi","hypo"),10)) + labs(title="13 Aug 2020") + 
    scale_fill_manual("Taxa",values=rainbow(12),labels=taxa2) +
    geom_hline(yintercept=0) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
          axis.text.x=element_text(size=6,family="Times"), plot.title = element_text(hjust = 0.5)) + ylab("Density (individual/L)") +
    theme(legend.position = "bottom", legend.margin = margin(0, 1, 2, 1),plot.margin = margin(0,0.6,0,0.3,unit = "cm"), axis.text.y = element_text(size=9, family="Times"),
          legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),legend.text = element_text(size=8),legend.title = element_text(size=9),legend.key.size = unit(0.3,"cm"))
  #dev.off()

  #jpeg("Figures/BVR_moretaxa_epivshypo_biomass_13Aug2020.jpg", width = 6, height = 4, units = "in",res = 300)
ggplot(subset(BVR.DVM.taxa.calcs.long, grepl("biomass",metric,ignore.case = TRUE) & DateTime=="08-13-2020"), aes(x=WaterColumn, y=value)) +
    geom_rect(data=subset(BVR.DVM.taxa.calcs.long,Hour == 'midnight' &grepl("biomass",metric,ignore.case = TRUE)),
          aes(fill=Hour),xmin=-Inf ,xmax = Inf, ymin = -Inf, ymax = Inf, fill = 'black', alpha = 0.053,inherit.aes = FALSE) +
    geom_bar(aes(fill=Taxa, alpha=WaterColumn), stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() +
    geom_errorbar(aes(ymin=value-SE, ymax=value+SE), width=.2,position=position_dodge(.9)) + scale_alpha_manual(values = c(0.4, 1)) +
    facet_wrap(Hour~Taxa, scales= 'free',ncol=6, strip.position = "right", labeller=labeller(Hour=as_labeller(facet_labeller_bot),Taxa=as_labeller(facet_labeller_top))) +
    theme(strip.text.y = element_text(size = 0 ,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank()) + 
    scale_x_discrete(name="", labels=rep(c("epi","hypo"),10)) + labs(title="13 Aug 2020") + 
    scale_fill_manual("Taxa",values=rainbow(12),labels=taxa2) +
    geom_hline(yintercept=0) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
          axis.text.x=element_text(size=6,family="Times"), plot.title = element_text(hjust = 0.5)) + ylab(expression(paste("Biomass (",mu,"g/L)"))) +
    theme(legend.position = "bottom", legend.margin = margin(0, 1, 2, 1),plot.margin = margin(0,0.6,0,0.3,unit = "cm"), axis.text.y = element_text(size=9, family="Times"),
          legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),legend.text = element_text(size=8),legend.title = element_text(size=9),legend.key.size = unit(0.3,"cm"))
  #dev.off()
