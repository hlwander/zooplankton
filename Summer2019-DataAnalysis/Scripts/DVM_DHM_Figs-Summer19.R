##Figures for GLEON poster, ESA poster, and more...
##Created 10Oct19

### Description of data --> full water column tows and epi tows from BVR summer 2019 (Jun 10-11)
    #includes samples collected at macrophyes (BVR_l), dam (BVR_d) and pelagic site (BVR_50_p for epi tows collected during MSN ONLY; BVR_50 for full water column tows and tows outside of 24-hour campaigns)
    #samples collected from at noon (x1), midnight (x1), sunset (x4), and sunrise (x4)


#TO DO:
#figure out rotifer biomass calcs (simple formula is sig different than other equation with l,w,h params, at least for keratella)

#libarries
pacman::p_load(plyr,plotrix,lubridate,dplyr,ggplot2,scales,tidyr,viridis)

#read in zoop summary csv
zoop<- read.csv('SummaryStats/FCR_ZooplanktonSummary2019.csv',header = TRUE)

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
                  Copepoda_PercentOfTotal, Copepoda_density_NopL, Copepoda_BiomassConcentration_ugpL, Copepoda_totalbiomass_ug,) %>%
                  group_by(sample_ID, site_no, Hour, collect_date) %>%
                  summarise_at(vars(ZoopDensity_No.pL:Copepoda_totalbiomass_ug,), funs(rep.mean=mean, rep.SE=stderr))

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
zoop_DHM<- zoop.repmeans[zoop.repmeans$sample_ID!="B_pel_10Jul19_noon_oxycline" & zoop.repmeans$site_no=="BVR_50_p" | zoop.repmeans$site_no=="BVR_l",substrEnd(colnames(zoop.repmeans),14)!="Total_rep.mean" & substrEnd(colnames(zoop.repmeans),12)!="Total_rep.SE"]
zoop_DHM_bytaxa<- zoop.repmeans.bytaxa[zoop.repmeans.bytaxa$sample_ID!="B_pel_10Jul19_noon_oxycline" & zoop.repmeans.bytaxa$site_no=="BVR_50_p" | zoop.repmeans.bytaxa$site_no=="BVR_l",]

#convert new dfs fron tibble to dataframe 
zoop_DHM <- data.frame(zoop_DHM)
zoop_DHM_bytaxa <- data.frame(zoop_DHM_bytaxa)

#----------------------------------------------------------------------------------#
##Fig - time vs. total density and biomass for MSN #1 and 2 (10-11 Jul, 24-25 Jul 2019)

variables <- c("ZoopDensity_No.pL_rep.mean","BiomassConcentration_ugpL_rep.mean","Cladocera_density_NopL_rep.mean","Cladocera_BiomassConcentration_ugpL_rep.mean",
               "Cyclopoida_density_NopL_rep.mean","Cyclopoida_BiomassConcentration_ugpL_rep.mean","Rotifera_density_NopL_rep.mean","Rotifera_BiomassConcentration_ugpL_rep.mean",
               "Calanoida_density_NopL_rep.mean","Calanoida_BiomassConcentration_ugpL_rep.mean", "Copepoda_density_NopL_rep.mean","Copepoda_BiomassConcentration_ugpL_rep.mean")
SE <- c("ZoopDensity_No.pL_rep.SE","BiomassConcentration_ugpL_rep.SE","Cladocera_density_NopL_rep.SE","Cladocera_BiomassConcentration_ugpL_rep.SE",
        "Cyclopoida_density_NopL_rep.SE","Cyclopoida_BiomassConcentration_ugpL_rep.SE","Rotifera_density_NopL_rep.SE","Rotifera_BiomassConcentration_ugpL_rep.SE",
        "Calanoida_density_NopL_rep.SE","Calanoida_BiomassConcentration_ugpL_rep.SE", "Copepoda_density_NopL_rep.SE","Copepoda_BiomassConcentration_ugpL_rep.SE")

#only select density and biomass #/ug / L cols
zoop_DHM <- zoop_DHM[,-c(which(grepl("ug_rep",colnames(zoop_DHM))))]

#convert df from wide to long (kinda hacky way bc having problems doing this)
df1 <- zoop_DHM %>% gather(metric,value,all_of(variables))
df2 <- zoop_DHM %>% gather(metric.SE,value.SE, all_of(SE))

#cut and paste to merge df
zoop_DHM_long <- df1[,c(1:4,17:18)]
#zoop_DHM_long$metric.SE <- df2$metric.SE
zoop_DHM_long$value.SE <- df2$value.SE

#drop _rep.mean from all metric names
zoop_DHM_long$metric <- substr(zoop_DHM_long$metric,1,nchar(zoop_DHM_long$metric)-9)

#change BVR_50_p to BVR_50
zoop_DHM_long$site_no[zoop_DHM_long$site_no=="BVR_50_p"] <- "BVR_50"

#drop oxy samples
zoop_DHM_long <- zoop_DHM_long[!c(grepl("oxy",zoop_DHM_long$sample_ID)),]

sites <- c("Pelagic","Littoral")
names(sites) <- c("BVR_50","BVR_l")

metric <- c("ZoopDensity_No.pL","BiomassConcentration_ugpL","Cladocera","Cladocera","Cyclopoida","Cyclopoida",
            "Rotifera","Rotifera","Calanoida","Calanoida","Copepoda","Copepoda")
names(metric) <- c("ZoopDensity_No.pL","BiomassConcentration_ugpL","Cladocera_density_NopL",
                   "Cladocera_BiomassConcentration_ugpL","Cyclopoida_density_NopL","Cyclopoida_BiomassConcentration_ugpL",
                   "Rotifera_density_NopL","Rotifera_BiomassConcentration_ugpL","Calanoida_density_NopL",
                   "Calanoida_BiomassConcentration_ugpL","Copepoda_density_NopL","Copepoda_BiomassConcentration_ugpL")

#total density and biomass Jul 10-11
ggplot(subset(zoop_DHM_long, metric %in% c("ZoopDensity_No.pL","BiomassConcentration_ugpL") & site_no %in% c("BVR_50","BVR_l") & (collect_date=="2019-07-10" | collect_date=="2019-07-11")), aes(Hour,value,color=site_no)) +
  geom_rect(aes(xmin=as.POSIXct("2019-07-10 11:30:00"),xmax=as.POSIXct("2019-07-10 20:41:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2019-07-10 20:42:00"),xmax=as.POSIXct("2019-07-11 06:10:00"), ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2019-07-11 06:11:00"),xmax=as.POSIXct("2019-07-11 12:30:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + facet_grid(site_no~metric,scales="free_y",labeller = labeller(site_no=sites, metric=metric)) + xlab("")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.background = element_blank(), legend.key = element_blank(), legend.key.height=unit(0.5,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_rect(fill = "transparent"), legend.position = c(0.9,0.9), legend.spacing = unit(0.05, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ scale_x_datetime(expand = c(0,0)) +
  scale_color_manual(values=c("#0099CC","#009966"), guide = 'none') + geom_line()+ #ylab("Density (Individuals/L)")+
  scale_fill_manual("",values=c("#CCCCCC","white"), guide = guide_legend(keyheight = .7,keywidth = 1.7, order=2, override.aes = list(alpha = 1,color="black")))+
  geom_errorbar(aes(ymin=value-value.SE, ymax=value+value.SE), width=.2,position=position_dodge(.9))
ggsave(file.path(getwd(),"Figures/BVR_2019_jul1_zoopdensity_biomass_LittoralvsPelagic.jpg")) 


#taxa density Jul 10-11
ggplot(subset(zoop_DHM_long, metric %in% c("Cladocera_density_NopL","Cyclopoida_density_NopL","Rotifera_density_NopL","Calanoida_density_NopL") &
                site_no %in% c("BVR_50","BVR_l") & (collect_date=="2019-07-10" | collect_date=="2019-07-11")), aes(Hour,value,color=site_no)) +
  geom_rect(aes(xmin=as.POSIXct("2019-07-10 11:30:00"),xmax=as.POSIXct("2019-07-10 20:41:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2019-07-10 20:42:00"),xmax=as.POSIXct("2019-07-11 06:10:00"), ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2019-07-11 06:11:00"),xmax=as.POSIXct("2019-07-11 12:30:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + facet_grid(site_no~metric,scales="free_y",labeller = labeller(site_no=sites,metric=metric)) + xlab("")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.background = element_blank(), legend.key = element_blank(), legend.key.height=unit(0.5,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_rect(fill = "transparent"), legend.position = c(0.93,0.96), legend.spacing = unit(0.05, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ scale_x_datetime(expand = c(0,0)) +
  scale_color_manual(values=c("#0099CC","#009966"), guide = 'none') + geom_line()+ ylab("Density (Individuals/L)")+
  scale_fill_manual("",values=c("#CCCCCC","white"), guide = guide_legend(keyheight = .7,keywidth = 1.7, order=2, override.aes = list(alpha = 1,color="black")))+
  geom_errorbar(aes(ymin=value-value.SE, ymax=value+value.SE), width=.2,position=position_dodge(.9))
ggsave(file.path(getwd(),"Figures/BVR_2019_jul1_taxa_density_LittoralvsPelagic.jpg")) 

#taxa biomass Jul 10-11
ggplot(subset(zoop_DHM_long, metric %in% c("Cladocera_BiomassConcentration_ugpL","Cyclopoida_BiomassConcentration_ugpL","Rotifera_BiomassConcentration_ugpL","Calanoida_BiomassConcentration_ugpL") &
                site_no %in% c("BVR_50","BVR_l") & (collect_date=="2019-07-10" | collect_date=="2019-07-11")), aes(Hour,value,color=site_no)) +
  geom_rect(aes(xmin=as.POSIXct("2019-07-10 11:30:00"),xmax=as.POSIXct("2019-07-10 20:41:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2019-07-10 20:42:00"),xmax=as.POSIXct("2019-07-11 06:10:00"), ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2019-07-11 06:11:00"),xmax=as.POSIXct("2019-07-11 12:30:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + facet_grid(site_no~metric,scales="free_y",labeller = labeller(site_no=sites,metric=metric)) + xlab("")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.background = element_blank(), legend.key = element_blank(), legend.key.height=unit(0.5,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_rect(fill = "transparent"), legend.position = c(0.93,0.96), legend.spacing = unit(0.05, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ scale_x_datetime(expand = c(0,0)) +
  scale_color_manual(values=c("#0099CC","#009966"), guide = 'none') + geom_line()+ ylab(expression(paste("Biomass (",mu,"g/L)")))+
  scale_fill_manual("",values=c("#CCCCCC","white"), guide = guide_legend(keyheight = .7,keywidth = 1.7, order=2, override.aes = list(alpha = 1,color="black")))+
  geom_errorbar(aes(ymin=value-value.SE, ymax=value+value.SE), width=.2,position=position_dodge(.9))
ggsave(file.path(getwd(),"Figures/BVR_2019_jul1_taxa_biomass_LittoralvsPelagic.jpg"))

#total density and biomass Jul 24-25
ggplot(subset(zoop_DHM_long, metric %in% c("ZoopDensity_No.pL","BiomassConcentration_ugpL") & site_no %in% c("BVR_50","BVR_l") & (collect_date=="2019-07-24" | collect_date=="2019-07-25")), aes(Hour,value,color=site_no)) +
  geom_rect(aes(xmin=as.POSIXct("2019-07-24 11:30:00"),xmax=as.POSIXct("2019-07-24 20:34:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2019-07-24 20:35:00"),xmax=as.POSIXct("2019-07-25 06:19:00"), ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2019-07-25 06:20:00"),xmax=as.POSIXct("2019-07-25 12:30:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + facet_grid(site_no~metric,scales="free_y",labeller = labeller(site_no=sites, metric=metric)) + xlab("")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.background = element_blank(), legend.key = element_blank(), legend.key.height=unit(0.5,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_rect(fill = "transparent"), legend.position = c(0.9,0.9), legend.spacing = unit(0.05, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ scale_x_datetime(expand = c(0,0)) +
  scale_color_manual(values=c("#0099CC","#009966"), guide = 'none') + geom_line()+ #ylab("Density (Individuals/L)")+
  scale_fill_manual("",values=c("#CCCCCC","white"), guide = guide_legend(keyheight = .7,keywidth = 1.7, order=2, override.aes = list(alpha = 1,color="black")))+
  geom_errorbar(aes(ymin=value-value.SE, ymax=value+value.SE), width=.2,position=position_dodge(.9))
ggsave(file.path(getwd(),"Figures/BVR_2019_jul2_zoopdensity_biomass_LittoralvsPelagic.jpg")) 


#taxa density Jul 24-25
ggplot(subset(zoop_DHM_long, metric %in% c("Cladocera_density_NopL","Cyclopoida_density_NopL","Rotifera_density_NopL","Calanoida_density_NopL") &
                site_no %in% c("BVR_50","BVR_l") & (collect_date=="2019-07-24" | collect_date=="2019-07-25")), aes(Hour,value,color=site_no)) +
  geom_rect(aes(xmin=as.POSIXct("2019-07-24 11:30:00"),xmax=as.POSIXct("2019-07-24 20:34:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2019-07-24 20:35:00"),xmax=as.POSIXct("2019-07-25 06:19:00"), ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2019-07-25 06:20:00"),xmax=as.POSIXct("2019-07-25 12:30:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + facet_grid(site_no~metric,scales="free_y",labeller = labeller(site_no=sites,metric=metric)) + xlab("")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.background = element_blank(), legend.key = element_blank(), legend.key.height=unit(0.5,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_rect(fill = "transparent"), legend.position = c(0.93,0.96), legend.spacing = unit(0.05, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ scale_x_datetime(expand = c(0,0)) +
  scale_color_manual(values=c("#0099CC","#009966"), guide = 'none') + geom_line()+ ylab("Density (Individuals/L)")+
  scale_fill_manual("",values=c("#CCCCCC","white"), guide = guide_legend(keyheight = .7,keywidth = 1.7, order=2, override.aes = list(alpha = 1,color="black")))+
  geom_errorbar(aes(ymin=value-value.SE, ymax=value+value.SE), width=.2,position=position_dodge(.9))
ggsave(file.path(getwd(),"Figures/BVR_2019_jul1_taxa_density_LittoralvsPelagic.jpg")) 

#taxa biomass Jul 24-25
ggplot(subset(zoop_DHM_long, metric %in% c("Cladocera_BiomassConcentration_ugpL","Cyclopoida_BiomassConcentration_ugpL","Rotifera_BiomassConcentration_ugpL","Calanoida_BiomassConcentration_ugpL") &
                site_no %in% c("BVR_50","BVR_l") & (collect_date=="2019-07-24" | collect_date=="2019-07-25")), aes(Hour,value,color=site_no)) +
  geom_rect(aes(xmin=as.POSIXct("2019-07-24 11:30:00"),xmax=as.POSIXct("2019-07-24 20:34:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2019-07-24 20:35:00"),xmax=as.POSIXct("2019-07-25 06:19:00"), ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2019-07-25 06:20:00"),xmax=as.POSIXct("2019-07-25 12:30:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + facet_grid(site_no~metric,scales="free_y",labeller = labeller(site_no=sites,metric=metric)) + xlab("")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.background = element_blank(), legend.key = element_blank(), legend.key.height=unit(0.5,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_rect(fill = "transparent"), legend.position = c(0.93,0.96), legend.spacing = unit(0.05, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ scale_x_datetime(expand = c(0,0)) +
  scale_color_manual(values=c("#0099CC","#009966"), guide = 'none') + geom_line()+ ylab(expression(paste("Biomass (",mu,"g/L)")))+
  scale_fill_manual("",values=c("#CCCCCC","white"), guide = guide_legend(keyheight = .7,keywidth = 1.7, order=2, override.aes = list(alpha = 1,color="black")))+
  geom_errorbar(aes(ymin=value-value.SE, ymax=value+value.SE), width=.2,position=position_dodge(.9))
ggsave(file.path(getwd(),"Figures/BVR_2019_jul2_taxa_biomass_LittoralvsPelagic.jpg"))

#------------------------------------------------------------------------------------#
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

#taxa density Jul 10-11
ggplot(subset(zoop_DHM_bytaxa_long, grepl("density",metric,ignore.case = TRUE) &
                site_no %in% c("BVR_50","BVR_l") & (collect_date=="2019-07-10" | collect_date=="2019-07-11")), aes(Hour,value,color=site_no)) +
  geom_rect(aes(xmin=as.POSIXct("2019-07-10 11:30:00"),xmax=as.POSIXct("2019-07-10 20:41:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2019-07-10 20:42:00"),xmax=as.POSIXct("2019-07-11 06:10:00"), ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2019-07-11 06:11:00"),xmax=as.POSIXct("2019-07-11 12:30:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + facet_grid(site_no~metric,scales="free_y",labeller = labeller(site_no=sites,metric=metric_taxa)) + xlab("")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=6), axis.text = element_text(size=6, color="black"), legend.background = element_blank(), legend.key = element_blank(), legend.key.height=unit(0.5,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_rect(fill = "transparent"), legend.position = c(0.96,0.96), legend.spacing = unit(0.05, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ scale_x_datetime(expand = c(0,0)) +
  scale_color_manual(values=c("#0099CC","#009966"), guide = 'none') + geom_line()+ ylab("Density (Individuals/L)")+
  scale_fill_manual("",values=c("#CCCCCC","white"), guide = guide_legend(keyheight = .5,keywidth = .5, order=2, override.aes = list(alpha = 1,color="black")))+
  geom_errorbar(aes(ymin=value-value.SE, ymax=value+value.SE), width=.2,position=position_dodge(.9))
ggsave(file.path(getwd(),"Figures/BVR_2019_jul1_taxa2_density_LittoralvsPelagic.jpg")) 

#taxa biomass Jul 10-11
ggplot(subset(zoop_DHM_bytaxa_long, grepl("biomass",metric,ignore.case = TRUE) &
                site_no %in% c("BVR_50","BVR_l") & (collect_date=="2019-07-10" | collect_date=="2019-07-11")), aes(Hour,value,color=site_no)) +
  geom_rect(aes(xmin=as.POSIXct("2019-07-10 11:30:00"),xmax=as.POSIXct("2019-07-10 20:41:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2019-07-10 20:42:00"),xmax=as.POSIXct("2019-07-11 06:10:00"), ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2019-07-11 06:11:00"),xmax=as.POSIXct("2019-07-11 12:30:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + facet_grid(site_no~metric,scales="free_y",labeller = labeller(site_no=sites,metric=metric_taxa)) + xlab("")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=6), axis.text = element_text(size=6, color="black"), legend.background = element_blank(), legend.key = element_blank(), legend.key.height=unit(0.5,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_rect(fill = "transparent"), legend.position = c(0.96,0.96), legend.spacing = unit(0.05, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ scale_x_datetime(expand = c(0,0)) +
  scale_color_manual(values=c("#0099CC","#009966"), guide = 'none') + geom_line()+ ylab(expression(paste("Biomass (",mu,"g/L)")))+
  scale_fill_manual("",values=c("#CCCCCC","white"), guide = guide_legend(keyheight = .5,keywidth = .5, order=2, override.aes = list(alpha = 1,color="black")))+
  geom_errorbar(aes(ymin=value-value.SE, ymax=value+value.SE), width=.2,position=position_dodge(.9))
ggsave(file.path(getwd(),"Figures/BVR_2019_jul1_taxa2_biomass_LittoralvsPelagic.jpg"))

#taxa density Jul 24-25
ggplot(subset(zoop_DHM_bytaxa_long, grepl("density",metric,ignore.case = TRUE) &
                site_no %in% c("BVR_50","BVR_l") & (collect_date=="2019-07-24" | collect_date=="2019-07-25")), aes(Hour,value,color=site_no)) +
  geom_rect(aes(xmin=as.POSIXct("2019-07-24 11:30:00"),xmax=as.POSIXct("2019-07-24 20:34:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2019-07-24 20:35:00"),xmax=as.POSIXct("2019-07-25 06:19:00"), ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2019-07-25 06:20:00"),xmax=as.POSIXct("2019-07-25 12:30:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + facet_grid(site_no~metric,scales="free_y",labeller = labeller(site_no=sites,metric=metric_taxa)) + xlab("")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=6), axis.text = element_text(size=6, color="black"), legend.background = element_blank(), legend.key = element_blank(), legend.key.height=unit(0.5,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_rect(fill = "transparent"), legend.position = c(0.96,0.96), legend.spacing = unit(0.05, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ scale_x_datetime(expand = c(0,0)) +
  scale_color_manual(values=c("#0099CC","#009966"), guide = 'none') + geom_line()+ ylab("Density (Individuals/L)")+
  scale_fill_manual("",values=c("#CCCCCC","white"), guide = guide_legend(keyheight = .5,keywidth = .5, order=2, override.aes = list(alpha = 1,color="black")))+
  geom_errorbar(aes(ymin=value-value.SE, ymax=value+value.SE), width=.2,position=position_dodge(.9))
ggsave(file.path(getwd(),"Figures/BVR_2019_jul2_taxa2_density_LittoralvsPelagic.jpg")) 

#taxa biomass Jul 24-25
ggplot(subset(zoop_DHM_bytaxa_long, grepl("biomass",metric,ignore.case = TRUE) &
                site_no %in% c("BVR_50","BVR_l") & (collect_date=="2019-07-24" | collect_date=="2019-07-25")), aes(Hour,value,color=site_no)) +
  geom_rect(aes(xmin=as.POSIXct("2019-07-24 11:30:00"),xmax=as.POSIXct("2019-07-24 20:34:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2019-07-24 20:35:00"),xmax=as.POSIXct("2019-07-25 06:19:00"), ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2019-07-25 06:20:00"),xmax=as.POSIXct("2019-07-25 12:30:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + facet_grid(site_no~metric,scales="free_y",labeller = labeller(site_no=sites,metric=metric_taxa)) + xlab("")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=6), axis.text = element_text(size=6, color="black"), legend.background = element_blank(), legend.key = element_blank(), legend.key.height=unit(0.5,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_rect(fill = "transparent"), legend.position = c(0.96,0.96), legend.spacing = unit(0.05, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ scale_x_datetime(expand = c(0,0)) +
  scale_color_manual(values=c("#0099CC","#009966"), guide = 'none') + geom_line()+ ylab(expression(paste("Biomass (",mu,"g/L)")))+
  scale_fill_manual("",values=c("#CCCCCC","white"), guide = guide_legend(keyheight = .5,keywidth = .5, order=2, override.aes = list(alpha = 1,color="black")))+
  geom_errorbar(aes(ymin=value-value.SE, ymax=value+value.SE), width=.2,position=position_dodge(.9))
ggsave(file.path(getwd(),"Figures/BVR_2019_jul2_taxa2_biomass_LittoralvsPelagic.jpg"))


###############################################################################################
  ###############################################################################################
  #### BVR epi vs hypo figures!! 
  #there should be 6 samples per MSN (epi and full tows from 2 noons and 1 midnight)
  #MSN1 oxycline is ~6m; MSN2 oxycline is 5.2(ish)m

#sum up counts by sample/site/day for DVM analyses + figs
BVR_counts <- zoop %>% select(sample_ID,site_no,collect_date,Hour, OverallCount_n,
  CladoceraCount_n, CyclopoidaCount_n, RotiferaCount_n, CalanoidaCount_n, CopepodaCount_n) %>%
  group_by(sample_ID, site_no, Hour, collect_date) %>%
  summarise_at(vars(OverallCount_n:CopepodaCount_n), funs(rep.mean=mean, rep.SE=stderr))

#add unadjusted volume
BVR_counts$Volume_unadj<- (zoop %>% select(sample_ID,site_no,collect_date,Hour, Volume_unadj) %>%
  group_by(sample_ID, site_no, Hour, collect_date) %>%
  summarise_at(vars(Volume_unadj),funs(rep.mean=mean)))$rep.mean

#add proportional volume for numerator hypo calc
BVR_counts$proportional_vol<- (zoop %>% select(sample_ID,site_no,collect_date,Hour, proportional_vol) %>%
  group_by(sample_ID, site_no, Hour, collect_date) %>%
  summarise_at(vars(proportional_vol),funs(rep.mean=mean)))$rep.mean

#get BVR_counts df in same order as zoop.repmeans df
BVR_counts<- BVR_counts[order(match(paste0(BVR_counts$sample_ID,BVR_counts$site_no,BVR_counts$collect_date), 
                                    paste0(zoop.repmeans$sample_ID,zoop.repmeans$site_no,zoop.repmeans$collect_date))),]

#add counts and vol to zoop.repmeans
zoop.repmeans[,paste0(colnames(BVR_counts[5:18]))]<- BVR_counts[5:18]

#new dfs for DVM data (BVR_pelagic_DVM_raw is just raw #/ug; BVR_pelagic_DVM_vol_calculated is #/L and ug/L)
BVR_pelagic_DVM<- zoop.repmeans[(zoop.repmeans$site_no=="BVR_50" |zoop.repmeans$site_no=="BVR_50_p") &
                                 (substrEnd(zoop.repmeans$sample_ID,5)=="night" | substrEnd(zoop.repmeans$sample_ID,4)=="noon" |
                                 substrEnd(zoop.repmeans$sample_ID,9)=="night_epi" | substrEnd(zoop.repmeans$sample_ID,8)=="noon_epi" |
                                 substrEnd(zoop.repmeans$sample_ID,3)=="oxy"),] 

#only select volume, count, and ug cols (plus SE cols)
BVR_pelagic_DVM_raw<- BVR_pelagic_DVM[,c(1:4,63,64,which(substrEnd(colnames(BVR_pelagic_DVM),10)=="n_rep.mean"),
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
BVR_pelagic_DVM_percent<- BVR_pelagic_DVM[,c(1:4,63,64,which(substrEnd(colnames(BVR_pelagic_DVM),14)=="Total_rep.mean"),
                                     which(substrEnd(colnames(BVR_pelagic_DVM),12)=="Total_rep.SE"))]

#drop 24Jul noon sample because missing full water column tow and first noon epi (using oxycline tow instead of 4m for 10Jul)
BVR_pelagic_DVM_raw<- BVR_pelagic_DVM_raw[!(BVR_pelagic_DVM_raw$sample_ID) %in% c("B_pel_24Jul19_noon_epi","B_pel_10Jul19_noon_epi"),]
BVR_pelagic_DVM_vol_calculated<- BVR_pelagic_DVM_vol_calculated[!(BVR_pelagic_DVM_vol_calculated$sample_ID) %in% c("B_pel_24Jul19_noon_epi","B_pel_10Jul19_noon_epi"),]
BVR_pelagic_DVM_percent <- BVR_pelagic_DVM_percent[!(BVR_pelagic_DVM_percent$sample_ID) %in% c("B_pel_24Jul19_noon_epi","B_pel_10Jul19_noon_epi"),]

#rename oxycline tow so ends in epi
BVR_pelagic_DVM_raw$sample_ID[BVR_pelagic_DVM_raw$sample_ID=="B_pel_10Jul19_noon_oxy"]<- "B_pel_10Jul19_noon_oxy_epi"
BVR_pelagic_DVM_vol_calculated$sample_ID[BVR_pelagic_DVM_vol_calculated$sample_ID=="B_pel_10Jul19_noon_oxy"]<- "B_pel_10Jul19_noon_oxy_epi"
BVR_pelagic_DVM_percent$sample_ID[BVR_pelagic_DVM_percent$sample_ID=="B_pel_10Jul19_noon_oxy"]<- "B_pel_10Jul19_noon_oxy_epi"

#calculate hypo SE for 11Jul19 midnight 
SE.hypo.calcs.raw<- zoop[zoop$sample_ID=="B_pel_11Jul19_midnight" |zoop$sample_ID=="B_pel_11Jul19_midnight_epi",]
matchingcols<- match(substr(colnames(BVR_pelagic_DVM_raw[1:18]),1,14),substr(colnames(SE.hypo.calcs.raw),1,14))
SE.hypo.calcs.raw<- SE.hypo.calcs.raw[,unique(matchingcols)]

SE.hypo.calcs.dens <- zoop[zoop$sample_ID=="B_pel_11Jul19_midnight" |zoop$sample_ID=="B_pel_11Jul19_midnight_epi",c(1:5,11,12,which(substrEnd(colnames(zoop),4)=="NopL"))]
matchingcols.dens<- match(substr(colnames(BVR_pelagic_DVM_vol_calculated[,c(1:4,6:10)]),1,14),substr(colnames(SE.hypo.calcs.dens),1,14))
SE.hypo.calcs.dens<- SE.hypo.calcs.dens[,matchingcols.dens]

#not sure if this is right, but calculating SE of difference between epi mean and hypo mean
SE.diffMean<- function(x,y){
  sqrt((sd(x,na.rm=TRUE)^2/length(na.omit(x))) + 
         (sd(y,na.rm=TRUE)^2/length(na.omit(y))))
}
  
#initialize df
BVR.DVM.calcs<- data.frame("Hour"=unique(BVR_pelagic_DVM_raw$Hour))
  
#for loop to fill out epi vs hypo calcs 
#hypo density and biomass calculated by subtracting epi raw zoop # from full zoop # and then dividing by the (full volume - epi volume) 
#NOTE: using epi density/L and biomass/L but calculating hypo using raw # and ug values. 
column.names<- colnames(BVR_pelagic_DVM_vol_calculated[,c(5:16)])
variables<- colnames(BVR_pelagic_DVM_raw[,c(7:18)])
percent<- colnames(BVR_pelagic_DVM_percent[,c(7:11)])
for(i in 1:length(variables)){
  BVR.DVM.calcs[,paste0(column.names,"_epi")[i]]<- BVR_pelagic_DVM_vol_calculated[substrEnd(BVR_pelagic_DVM_vol_calculated$sample_ID,3)=="epi",paste0(column.names)[i]]
  BVR.DVM.calcs[,paste0(column.names,"_hypo")[i]] <- (((1/BVR_pelagic_DVM_raw[substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)!="epi" ,"proportional_vol"]) * BVR_pelagic_DVM_raw[substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)!="epi" ,paste0(variables)[i]] * (1/NetEfficiency2016)) - 
                                                     ((1/BVR_pelagic_DVM_raw[substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)=="epi", "proportional_vol"]) *  BVR_pelagic_DVM_raw[substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)=="epi",paste0(variables)[i]] * (1/NetEfficiency2016)))/
                                                     (BVR_pelagic_DVM_raw[substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)!="epi" ,"Volume_unadj"] - BVR_pelagic_DVM_raw[substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)=="epi", "Volume_unadj"])  

}
density.percent<- colnames(BVR_pelagic_DVM_vol_calculated[,c(6:10)])
for(i in 1:length(density.percent)){
for(j in 1:length(unique(BVR.DVM.calcs$Hour))){
    BVR.DVM.calcs[j,paste0(density.percent,"_epi_percent_density")[i]]<- (BVR.DVM.calcs[j,paste0(density.percent,"_epi")][i]/ sum(BVR.DVM.calcs[j,paste0(density.percent,"_epi")[i]],BVR.DVM.calcs[j,paste0(density.percent,"_hypo")[i]])) *100
    BVR.DVM.calcs[j,paste0(density.percent,"_hypo_percent_density")[i]]<- (BVR.DVM.calcs[j,paste0(density.percent,"_hypo")][i]/ sum(BVR.DVM.calcs[j,paste0(density.percent,"_epi")[i]],BVR.DVM.calcs[j,paste0(density.percent,"_hypo")[i]])) * 100

}       
}


#initialize df
BVR.DVM.calcs.SE<- data.frame("Hour"=unique(BVR_pelagic_DVM_raw$Hour))

#can only calculate hypo SE for 11JUl19 midnight because only time that had reps of epi and full samples
SEonly<- colnames(SE.hypo.calcs.raw)[7:18]
Percentdens <- colnames(SE.hypo.calcs.dens)[5:9]
for(i in 1:length(SEonly)){
  BVR.DVM.calcs.SE[,paste0(column.names,"_epi_SE")[i]] <- BVR_pelagic_DVM_vol_calculated[substrEnd(BVR_pelagic_DVM_vol_calculated$sample_ID,3)=="epi",substrEnd(colnames(BVR_pelagic_DVM_vol_calculated),2)=="SE"][i]
  BVR.DVM.calcs.SE[2,paste0(column.names,"_hypo_SE")[i]]<- SE.diffMean(SE.hypo.calcs.raw[SE.hypo.calcs.raw$sample_ID=="B_pel_11Jul19_midnight",paste0(SEonly)[i]],
                                                           SE.hypo.calcs.raw[SE.hypo.calcs.raw$sample_ID=="B_pel_11Jul19_midnight_epi",paste0(SEonly)[i]])
}
  for (j in 1:length(Percentdens)){
  BVR.DVM.calcs.SE[,paste0(Percentdens,"_epi_percent_density_SE")[j]]<- BVR_pelagic_DVM_vol_calculated[substrEnd(BVR_pelagic_DVM_vol_calculated$sample_ID,3)=="epi",substrEnd(colnames(BVR_pelagic_DVM_vol_calculated),11)=="NopL_rep.SE"][j]
  BVR.DVM.calcs.SE[2,paste0(Percentdens,"_hypo_percent_density_SE")[j]] <- SE.diffMean(SE.hypo.calcs.dens[SE.hypo.calcs.dens$sample_ID=="B_pel_11Jul19_midnight",paste0(Percentdens)[j]],
                                                                                       SE.hypo.calcs.dens[SE.hypo.calcs.dens$sample_ID=="B_pel_11Jul19_midnight_epi",paste0(Percentdens)[j]])
  }
    
#reorder BVR.DVM.calcs so order matches BVR.DVM.calcs.SE
#BVR.DVM.calcs<- BVR.DVM.calcs[,c(1:22,26,23,27,24,28,25,29)]

#wide to long for both dfs separately
BVR.DVM.calcs.long <-  BVR.DVM.calcs %>%
    gather(metric,value, ZoopDensity_No.pL_rep.mean_epi:Copepoda_density_NopL_rep.mean_hypo_percent_density) %>%
    mutate(DateTime = strftime(Hour, "%m-%d-%Y %H:%M"))
BVR.DVM.calcs.SE.long <-  BVR.DVM.calcs.SE %>%
gather(taxa.metric,SE,ZoopDensity_No.pL_rep.mean_epi_SE:Copepoda_density_NopL_hypo_percent_density_SE)

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

#export 2019 dvm stats
write.csv(BVR.DVM.calcs.long,"./SummaryStats/DVM_2019_zoops.csv",row.names = FALSE)

#######################################
#######################################
## Create new df with more diverse taxa
#subsetting repmeans into new df for DVM analyses/figs

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

#drop 24Jul noon sample because missing full water column tow and first noon epi (to replace with oxycline tow)
BVR_pelagic_taxa_DVM_raw<- BVR_pelagic_taxa_DVM_raw[!(BVR_pelagic_taxa_DVM_raw$sample_ID) %in% c("B_pel_24Jul19_noon_epi","B_pel_10Jul19_noon_epi"),]
BVR_pelagic_taxa_DVM_vol_calculated<- BVR_pelagic_taxa_DVM_vol_calculated[!(BVR_pelagic_taxa_DVM_vol_calculated$sample_ID) %in% c("B_pel_24Jul19_noon_epi","B_pel_10Jul19_noon_epi"),]

#rename oxycline tow so ends in epi
BVR_pelagic_taxa_DVM_raw$sample_ID[BVR_pelagic_taxa_DVM_raw$sample_ID=="B_pel_10Jul19_noon_oxycline"]<- 	"B_pel_10Jul19_noon_oxycline_epi"
BVR_pelagic_taxa_DVM_vol_calculated$sample_ID[BVR_pelagic_taxa_DVM_vol_calculated$sample_ID=="B_pel_10Jul19_noon_oxycline"]<- 	"B_pel_10Jul19_noon_oxycline_epi"

#calculate hypo SE for 11Jul19 midnight 
SE.taxa.hypo.calcs.raw<- zoop[zoop$sample_ID=="B_pel_11Jul19_midnight" |zoop$sample_ID=="B_pel_11Jul19_midnight_epi",]
matchingcols<- match(substr(colnames(BVR_pelagic_taxa_DVM_raw[1:30]),1,15),substr(colnames(SE.taxa.hypo.calcs.raw),1,15))
#manually replace NA with 117 (DaphniaCount_n is too short) and manually choose Trichocercidae_totalbiomass_ug (wants to pull percent because string too long; #175)
matchingcols[7]<-117
matchingcols[30]<- 175
SE.taxa.hypo.calcs.raw<- SE.taxa.hypo.calcs.raw[,unique(matchingcols)]

#not sure if this is right, but calculating SE of difference between epi mean and hypo mean
SE.diffMean<- function(x,y){
  sqrt((sd(x,na.rm=TRUE)^2/length(na.omit(x))) + 
         (sd(y,na.rm=TRUE)^2/length(na.omit(y))))
}

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

#can only calculate hypo SE for 11JUl19 midnight because only time that had reps of epi and full samples
SEonly<- colnames(SE.taxa.hypo.calcs.raw[,c(7:30)])
for(i in 1:length(SEonly)){
  BVR.DVM.taxa.calcs.SE[,paste0(column.names,"_epi_SE")[i]] <- BVR_pelagic_taxa_DVM_vol_calculated[substrEnd(BVR_pelagic_taxa_DVM_vol_calculated$sample_ID,3)=="epi",substrEnd(colnames(BVR_pelagic_taxa_DVM_vol_calculated),2)=="SE"][i]
  BVR.DVM.taxa.calcs.SE[2,paste0(column.names,"_hypo_SE")[i]]<- SE.diffMean(SE.taxa.hypo.calcs.raw[SE.taxa.hypo.calcs.raw$sample_ID=="B_pel_11Jul19_midnight",paste0(SEonly)[i]],
                                                                       SE.taxa.hypo.calcs.raw[SE.taxa.hypo.calcs.raw$sample_ID=="B_pel_11Jul19_midnight_epi",paste0(SEonly)[i]])
}

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

############################
############################
#additional df to play around with epi vs meta vs hypo data
#10Jul19 noon is only time we collected tow from 4m, 5.5m, and 10m (all other times were just 4 and 8 :(... )
  #Note: can't calculate hypo SE becasue don't have noon oxycline reps yet

BVR.DVM.10Jul19.tows<- zoop.repmeans[zoop.repmeans$sample_ID=="B_pel_10Jul19_noon" | zoop.repmeans$sample_ID=="B_pel_10Jul19_noon_epi" |
                                       zoop.repmeans$sample_ID=="B_pel_10Jul19_noon_oxycline",]

#initialize df
BVR.DVM.10Jul19.tows.calcs<- data.frame("Hour"=unique(BVR.DVM.10Jul19.tows$Hour))

#replace NAs with

##calculate epi, meta, and hypo (meta and hypo using raw data)
raw_variables<- colnames(BVR.DVM.10Jul19.tows)[c(26,7,29,10,33,14,36,17,42,23)]
variables<- colnames(BVR.DVM.10Jul19.tows)[c(5,6,8,9,12,13,16,18,21,22)] 
for(i in 1:length(variables)){
  BVR.DVM.10Jul19.tows.calcs[,paste0(variables,"_epi")[i]]<- BVR.DVM.10Jul19.tows[substrEnd(BVR.DVM.10Jul19.tows$sample_ID,3)=="epi",paste0(variables)[i]]
  #calculating meta as oxy - epi / oxy vol - epi vol)
  BVR.DVM.10Jul19.tows.calcs[,paste0(variables,"_meta")[i]] <- (((1/BVR.DVM.10Jul19.tows[substrEnd(BVR.DVM.10Jul19.tows$sample_ID,5)=="cline" ,"proportional_vol"]) *  BVR.DVM.10Jul19.tows[substrEnd(BVR.DVM.10Jul19.tows$sample_ID,5)=="cline" ,paste0(variables)[i]]* (1/0.31)) - 
                                                                  ((1/BVR.DVM.10Jul19.tows[substrEnd(BVR.DVM.10Jul19.tows$sample_ID,3)=="epi","proportional_vol"]) *  BVR.DVM.10Jul19.tows[substrEnd(BVR.DVM.10Jul19.tows$sample_ID,3)=="epi",paste0(variables)[i]]* (1/0.31))) /
                                                                  (BVR.DVM.10Jul19.tows[substrEnd(BVR.DVM.10Jul19.tows$sample_ID,5)=="cline" ,"Volume_unadj"] - BVR.DVM.10Jul19.tows[substrEnd(BVR.DVM.10Jul19.tows$sample_ID,3)=="epi", "Volume_unadj"])
  BVR.DVM.10Jul19.tows.calcs[,paste0(variables,"_hypo")[i]] <- (((1/BVR.DVM.10Jul19.tows[substrEnd(BVR.DVM.10Jul19.tows$sample_ID,4)=="noon" ,"proportional_vol"]) * BVR.DVM.10Jul19.tows[substrEnd(BVR.DVM.10Jul19.tows$sample_ID,4)=="noon" ,paste0(variables)[i]]* (1/0.053)) - 
                                                                  ((1/BVR.DVM.10Jul19.tows[substrEnd(BVR.DVM.10Jul19.tows$sample_ID,5)=="cline","proportional_vol"]) * BVR.DVM.10Jul19.tows[substrEnd(BVR.DVM.10Jul19.tows$sample_ID,5)=="cline",paste0(variables)[i]]* (1/0.31))) /
                                                                  (BVR.DVM.10Jul19.tows[substrEnd(BVR.DVM.10Jul19.tows$sample_ID,4)=="noon" ,"Volume_unadj"] - BVR.DVM.10Jul19.tows[substrEnd(BVR.DVM.10Jul19.tows$sample_ID,5)=="cline", "Volume_unadj"])
}

#wide to long df
BVR.DVM.10Jul19.tows.calcs.long <-  BVR.DVM.10Jul19.tows.calcs %>%
  gather(metric,value, ZoopDensity_No.pL_rep.mean_epi:Calanoida_BiomassConcentration_ugpL_rep.mean_hypo) %>%
  mutate(DateTime = strftime(Hour, "%m-%d-%Y %H:%M"))

#add watercolumn and hour columns
BVR.DVM.10Jul19.tows.calcs.long$WaterColumn <- ifelse(substrEnd(BVR.DVM.10Jul19.tows.calcs.long$metric,3)=="epi","epilimnion",ifelse(substrEnd(BVR.DVM.10Jul19.tows.calcs.long$metric,4)=="meta","metalimnion","hypolimnion"))
BVR.DVM.10Jul19.tows.calcs.long$Hour <- ifelse(substrEnd(BVR.DVM.10Jul19.tows.calcs.long$DateTime,5)=="12:00","noon","midnight")
BVR.DVM.10Jul19.tows.calcs.long$Taxa <- substr(BVR.DVM.10Jul19.tows.calcs.long$metric,1,9)


###############################
#### Epi vs Hypo figures #####
###############################
#df for percent density calcs
BVR.DVM.percent.dens.calcs<- BVR.DVM.calcs.long[substrEnd(BVR.DVM.calcs.long$metric,15)=="percent_density",]

#drop total density and biomass rows
BVR.DVM.calcs.long<- BVR.DVM.calcs.long[!(BVR.DVM.calcs.long$metric=="ZoopDensity_No.pL_rep.mean_epi" | BVR.DVM.calcs.long$metric=="ZoopDensity_No.pL_rep.mean_hypo" | BVR.DVM.calcs.long$metric=="BiomassConcentration_ugpL_rep.mean_epi" |BVR.DVM.calcs.long$metric=="BiomassConcentration_ugpL_rep.mean_hypo" | substrEnd(BVR.DVM.calcs.long$metric,15)=="percent_density"),] 
#reorder 
BVR.DVM.calcs.long$Taxa <- factor(BVR.DVM.calcs.long$Taxa, levels = c("Cladocera","Rotifera_","Cyclopoid","Calanoida"))

#jpeg("BVR_epivshypo_density_10Jul2019.jpg", width = 6, height = 4, units = "in",res = 300)
ggplot(subset(BVR.DVM.calcs.long, grepl("density",metric,ignore.case = TRUE) & DateTime =="07-10-2019"), aes(x=WaterColumn, y=value)) +
    geom_bar(aes(fill=Taxa, alpha=WaterColumn), stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() +
    geom_errorbar(aes(ymin=value-SE, ymax=value+SE), width=.2,position=position_dodge(.9)) + scale_alpha_manual(values = c(0.4, 1)) +
    facet_wrap(~Taxa, scales= 'free_y',ncol=4, strip.position = "right", labeller=as_labeller(function(variable,value){c("","","","Noon")})) +
    theme(strip.text.y = element_text(size = 10,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank(), plot.margin = margin(2.8,0,2,0,unit = "cm")) + 
    scale_x_discrete(name="", labels=rep(c("epi","hypo"),10)) + labs(title="10 Jul 2019 Noon") + 
    scale_fill_manual(values=viridis(6),labels=c("Cladocera","Rotifera","Cyclopoida","Calanoida"),drop=TRUE) +
    geom_hline(yintercept=0) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
       axis.text.x=element_text(size=9,family="Times"), plot.title = element_text(hjust = 0.5)) + ylab("Density (individual/L)") +
    theme(legend.position = "bottom", legend.margin = margin(0, 1, 2, 1), legend.title = element_text(size=10),legend.key.size = unit(0.3,"cm"),
        legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),legend.text = element_text(size=8),axis.text.y = element_text(size=10, family="Times"))
#dev.off()

#jpeg("BVR_epivshypo_biomass_10Jul2019.jpg", width = 6, height = 4, units = "in",res = 300)
ggplot(subset(BVR.DVM.calcs.long, grepl("biomass",metric,ignore.case = TRUE) & DateTime =="07-10-2019"), aes(x=WaterColumn, y=value)) +
  geom_bar(aes(fill=Taxa, alpha=WaterColumn), stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() +
  geom_errorbar(aes(ymin=value-SE, ymax=value+SE), width=.2,position=position_dodge(.9)) + scale_alpha_manual(values = c(0.4, 1)) +
  facet_wrap(~Taxa, scales= 'free',ncol=4, strip.position = "right", labeller=as_labeller(function(variable,value){c("","","","Noon")})) +
  theme(strip.text.y = element_text(size = 10,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank(), plot.margin = margin(2.8,0,2,0,unit = "cm")) + 
  scale_x_discrete(name="", labels=rep(c("epi","hypo"),10)) + labs(title="10 Jul 2019 Noon") + 
  scale_fill_manual(values=viridis(6),labels=c("Cladocera","Rotifera","Cyclopoida","Calanoida")) +
  geom_hline(yintercept=0) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
     axis.text.x=element_text(size=9,family="Times"), plot.title = element_text(hjust = 0.5)) + ylab(expression(paste("Biomass (",mu,"g/L)"))) +
  theme(legend.position = "bottom", legend.margin = margin(0, 1, 2, 1), legend.title = element_text(size=10),legend.key.size = unit(0.3,"cm"),
        legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),legend.text = element_text(size=8),axis.text.y = element_text(size=10, family="Times"))
#dev.off()


#loop for last 2 dates
dates<- c("11 Jul 2019","25 Jul 2019")

#Epi vs hypo taxa density
for(d in 1:length(dates)){
  #jpeg(paste0("BVR_epivshypo_density_",dates[d],".jpg"), width = 6, height = 4, units = "in",res = 300)
  plot<- ggplot(subset(BVR.DVM.calcs.long, grepl("density",metric,ignore.case = TRUE) & DateTime %in%
         unique(DateTime)[2:3][d]), aes(x=WaterColumn, y=value)) +
    geom_rect(data=subset(BVR.DVM.calcs.long,Hour == 'midnight' &grepl("density",metric,ignore.case = TRUE)),
         aes(fill=Hour),xmin=-Inf ,xmax = Inf, ymin = -Inf, ymax = Inf, fill = 'black', alpha = 0.053,inherit.aes = FALSE) +
    geom_bar(aes(fill=Taxa, alpha=WaterColumn), stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() +
    geom_errorbar(aes(ymin=value-SE, ymax=value+SE), width=.2,position=position_dodge(.9)) + scale_alpha_manual(values = c(0.4, 1)) +
    facet_wrap(Hour~Taxa, scales= 'free',ncol=4, strip.position = "right", labeller=labeller(Hour=as_labeller(facet_labeller_bot),Taxa=as_labeller(facet_labeller_top))) +
    theme(strip.text.y = element_text(size = 11 ,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank()) + 
    scale_x_discrete(name="", labels=rep(c("epi","hypo"),10)) + labs(title=paste0(dates[d])) + 
    scale_fill_manual("Taxa",values=viridis(6),labels=c("Cladocera","Rotifera","Cyclopoida","Calanoida")) +
    geom_hline(yintercept=0) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
         axis.text.x=element_text(size=6,family="Times"), plot.title = element_text(hjust = 0.5)) + ylab("Density (individual/L)") +
    theme(legend.position = "bottom", legend.margin = margin(0, 1, 2, 1),plot.margin = margin(0,0,0,0.3,unit = "cm"), axis.text.y = element_text(size=13, family="Times"),
          legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),legend.text = element_text(size=10),legend.title = element_text(size=12))
 print(plot)
   #dev.off()
}

#create new column to convert negative values to 0
BVR.DVM.calcs.long$value_no_neg <- ifelse(BVR.DVM.calcs.long$value < 0, BVR.DVM.calcs.long$value_no_neg <- 0,
                                          BVR.DVM.calcs.long$value_no_neg <- BVR.DVM.calcs.long$value)

 #Epi vs hypo taxa biomass 
  for(d in 1:length(dates)){
#jpeg(paste0("BVR_epivshypo_biomass_",dates[d],".jpg"), width = 6, height = 4, units = "in",res = 300)
plot<- ggplot(subset(BVR.DVM.calcs.long, grepl("biomass",metric,ignore.case = TRUE) & DateTime %in%
           unique(DateTime)[2:3][d]), aes(x=WaterColumn, y=value_no_neg)) +
   geom_rect(data=subset(BVR.DVM.calcs.long,Hour == 'midnight' &grepl("biomass",metric,ignore.case = TRUE)),
           aes(fill=Hour),xmin=-Inf ,xmax = Inf, ymin = -Inf, ymax = Inf, fill = 'black', alpha = 0.053,inherit.aes = FALSE) +
   geom_bar(aes(fill=Taxa, alpha=WaterColumn), stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() +
   geom_errorbar(aes(ymin=value_no_neg-SE, ymax=value_no_neg+SE), width=.2,position=position_dodge(.9)) + scale_alpha_manual(values = c(0.4, 1)) +
   facet_wrap(Hour~Taxa, scales= 'free',ncol=4, strip.position = "right", labeller=labeller(Hour=as_labeller(facet_labeller_bot),Taxa=as_labeller(facet_labeller_top))) +
   theme(strip.text.y = element_text(size = 12, face="bold" ,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank()) + 
   scale_x_discrete(name="", labels=rep(c("epi","hypo"),10)) + labs(title=paste0(dates[d])) + guides(alpha=FALSE) +
   scale_fill_manual("Taxa",values=viridis(6),labels=c("Cladocera","Rotifera","Cyclopoida","Calanoida")) +
   geom_hline(yintercept=0) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
            axis.text.x=element_text(size=7,family="Times"), plot.title = element_text(hjust = 0.5)) + ylab(expression(paste("Biomass (",mu,"g/L)"))) +
      theme(legend.position = "bottom", legend.margin = margin(0, 1, 2, 1),plot.margin = margin(0,0,0,0.3,unit = "cm"), axis.text.y = element_text(size=13, family="Times"),
            legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),legend.text = element_text(size=10),legend.title = element_text(size=12), axis.title.y = element_text(size = 13))
    print(plot)
    #dev.off()
  }
 
#reorder taxa
BVR.DVM.percent.dens.calcs$Taxa <- factor(BVR.DVM.percent.dens.calcs$Taxa, levels = c("Cladocera","Rotifera_","Cyclopoid","Calanoida"))

####Epi vs hypo taxa percent density fig ####
for(d in 1:length(dates)){
  #jpeg(paste0("BVR_epivshypo_percent_density_",dates[d],".jpg"), width = 6, height = 4, units = "in",res = 300)
  plot<- ggplot(subset(BVR.DVM.percent.dens.calcs,  DateTime %in% unique(DateTime)[2:3][d]), aes(x=WaterColumn, y=value)) +
    geom_rect(data=subset(BVR.DVM.percent.dens.calcs,Hour == 'midnight' &grepl("density",metric,ignore.case = TRUE)),
         aes(fill=Hour),xmin=-Inf ,xmax = Inf, ymin = -Inf, ymax = Inf, fill = 'black', alpha = 0.053,inherit.aes = FALSE) +
    geom_bar(aes(fill=Taxa, alpha=WaterColumn), stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() +
    geom_errorbar(aes(ymin=value-SE, ymax=value+SE), width=.2,position=position_dodge(.9)) + scale_alpha_manual(values = c(0.4, 1)) +
    facet_wrap(Hour~Taxa, scales= 'free',ncol=4, strip.position = "right", labeller=labeller(Hour=as_labeller(facet_labeller_bot),Taxa=as_labeller(facet_labeller_top))) +
    theme(strip.text.y = element_text(size = 13 ,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank()) + 
    scale_x_discrete(name="", labels=rep(c("epi","hypo"),10)) + labs(title=paste0(dates[d])) + guides(alpha=FALSE) + 
    scale_fill_manual("Taxa",values=viridis(6),labels=c("Cladocera","Rotifera","Cyclopoida","Calanoida")) +
    geom_hline(yintercept=0) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
         axis.text.x=element_text(size=6,family="Times"), plot.title = element_text(hjust = 0.5)) + ylab("Percent Density") +
    theme(legend.position = "bottom", legend.margin = margin(1, 1, 2, 1),plot.margin = margin(0,0,0,0.3,unit = "cm"), axis.text.y = element_text(size=13, family="Times"), axis.title.y = element_text(size=14, family="Times"), 
         legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),legend.text = element_text(size=10),legend.title = element_text(size=12))
  print(plot)
  #dev.off()
}

###################################################
# more diverse taxa figs

#order taxa to get bars in correct order
BVR.DVM.taxa.calcs.long$Taxa <- factor(BVR.DVM.taxa.calcs.long$Taxa, levels= c("Daphnia_d","Daphnia_B","Ceriodaph","Bosminida",
             "Cyclopoid","Calanoida","nauplius_","Keratella","Kellicott","Collothec","Conochili","Synchaeti","Trichocer"))

#jpeg("BVR_moretaxa_epivshypo_density_10Jul2019.jpg", width = 6, height = 4, units = "in",res = 300)
ggplot(subset(BVR.DVM.taxa.calcs.long, grepl("density",metric,ignore.case = TRUE) & DateTime =="07-10-2019"), aes(x=WaterColumn, y=value)) +
  geom_bar(aes(fill=Taxa, alpha=WaterColumn), stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() +
  geom_errorbar(aes(ymin=value-SE, ymax=value+SE), width=.2,position=position_dodge(.9)) + scale_alpha_manual(values = c(0.4, 1)) +
  facet_wrap(~Taxa, scales= 'free',ncol=6, strip.position = "right", labeller=as_labeller(function(variable,value){c("","","","","","")}),drop=TRUE) +
  theme(strip.text.y = element_text(size = 10,margin = margin(0, -0.03, 0,0.1, "cm")),strip.background = element_blank(), plot.margin = margin(1,0,1,0,unit = "cm")) + 
  scale_x_discrete(name="", labels=rep(c("epi","hypo"),10)) + labs(title="10 Jul 2019 Noon") + 
  scale_fill_manual(values=rainbow(12),labels=taxa2) +
  geom_hline(yintercept=0) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text.x=element_text(size=8,family="Times"), plot.title = element_text(hjust = 0.5)) + ylab("Density (individual/L)") +
  theme(legend.position = "bottom", legend.margin = margin(0, 1, 2, 1), legend.title = element_text(size=10),legend.key.size = unit(0.3,"cm"),
        legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),legend.text = element_text(size=8),axis.text.y = element_text(size=10, family="Times"))
#dev.off()

#jpeg("BVR_moretaxa_epivshypo_biomass_10Jul2019.jpg", width = 6, height = 4, units = "in",res = 300)
ggplot(subset(BVR.DVM.taxa.calcs.long, grepl("biomass",metric,ignore.case = TRUE) & DateTime =="07-10-2019"), aes(x=WaterColumn, y=value)) +
  geom_bar(aes(fill=Taxa, alpha=WaterColumn), stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() +
  geom_errorbar(aes(ymin=value-SE, ymax=value+SE), width=.2,position=position_dodge(.9)) + scale_alpha_manual(values = c(0.4, 1)) +
  facet_wrap(~Taxa, scales= 'free',ncol=6, strip.position = "right", labeller=as_labeller(function(variable,value){c("","","","","","")})) +
  theme(strip.text.y = element_text(size = 10,margin = margin(0, -0.03, 0,0.1, "cm")),strip.background = element_blank(), plot.margin = margin(1,0,1,0,unit = "cm")) + 
  scale_x_discrete(name="", labels=rep(c("epi","hypo"),10)) + labs(title="10 Jul 2019 Noon") + 
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

#loop for last 2 dates
dates<- c("11 Jul 2019","25 Jul 2019")

#Epi vs hypo taxa density
for(d in 1:length(dates)){
  #jpeg(paste0("BVR_epivshypo_density_",dates[d],".jpg"), width = 6, height = 4, units = "in",res = 300)
  plot<- ggplot(subset(BVR.DVM.taxa.calcs.long, grepl("density",metric,ignore.case = TRUE) & DateTime %in%
           unique(DateTime)[2:3][d]), aes(x=WaterColumn, y=value)) +
    geom_rect(data=subset(BVR.DVM.taxa.calcs.long,Hour == 'midnight' &grepl("density",metric,ignore.case = TRUE)),
           aes(fill=Hour),xmin=-Inf ,xmax = Inf, ymin = -Inf, ymax = Inf, fill = 'black', alpha = 0.053,inherit.aes = FALSE) +
    geom_bar(aes(fill=Taxa, alpha=WaterColumn), stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() +
    geom_errorbar(aes(ymin=value-SE, ymax=value+SE), width=.2,position=position_dodge(.9)) + scale_alpha_manual(values = c(0.4, 1)) +
    facet_wrap(Hour~Taxa, scales= 'free',ncol=6, strip.position = "right", labeller=labeller(Hour=as_labeller(facet_labeller_bot),Taxa=as_labeller(facet_labeller_top))) +
    theme(strip.text.y = element_text(size = 0 ,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank()) + 
    scale_x_discrete(name="", labels=rep(c("epi","hypo"),10)) + labs(title=paste0(dates[d])) + 
    scale_fill_manual("Taxa",values=rainbow(12),labels=taxa2) +
    geom_hline(yintercept=0) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
          axis.text.x=element_text(size=6,family="Times"), plot.title = element_text(hjust = 0.5)) + ylab("Density (individual/L)") +
    theme(legend.position = "bottom", legend.margin = margin(0, 1, 2, 1),plot.margin = margin(0,0.6,0,0.3,unit = "cm"), axis.text.y = element_text(size=9, family="Times"),
          legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),legend.text = element_text(size=8),legend.title = element_text(size=9),legend.key.size = unit(0.3,"cm"))
  print(plot)
  #dev.off()
}

#reorder 
#BVR.DVM.taxa.calcs.long$Taxa<- factor(BVR.DVM.calcs.long$Taxa, levels = taxa2)

#Epi vs hypo taxa biomass 
for(d in 1:length(dates)){
  #jpeg(paste0("BVR_moretaxa_epivshypo_biomass_",dates[d],".jpg"), width = 6, height = 4, units = "in",res = 300)
  plot<- ggplot(subset(BVR.DVM.taxa.calcs.long, grepl("biomass",metric,ignore.case = TRUE) & DateTime %in%
         unique(DateTime)[2:3][d]), aes(x=WaterColumn, y=value)) +
    geom_rect(data=subset(BVR.DVM.taxa.calcs.long,Hour == 'midnight' &grepl("biomass",metric,ignore.case = TRUE)),
          aes(fill=Hour),xmin=-Inf ,xmax = Inf, ymin = -Inf, ymax = Inf, fill = 'black', alpha = 0.053,inherit.aes = FALSE) +
    geom_bar(aes(fill=Taxa, alpha=WaterColumn), stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() +
    geom_errorbar(aes(ymin=value-SE, ymax=value+SE), width=.2,position=position_dodge(.9)) + scale_alpha_manual(values = c(0.4, 1)) +
    facet_wrap(Hour~Taxa, scales= 'free',ncol=6, strip.position = "right", labeller=labeller(Hour=as_labeller(facet_labeller_bot),Taxa=as_labeller(facet_labeller_top))) +
    theme(strip.text.y = element_text(size = 0 ,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank()) + 
    scale_x_discrete(name="", labels=rep(c("epi","hypo"),10)) + labs(title=paste0(dates[d])) + 
    scale_fill_manual("Taxa",values=rainbow(12),labels=taxa2) +
    geom_hline(yintercept=0) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
          axis.text.x=element_text(size=6,family="Times"), plot.title = element_text(hjust = 0.5)) + ylab(expression(paste("Biomass (",mu,"g/L)"))) +
    theme(legend.position = "bottom", legend.margin = margin(0, 1, 2, 1),plot.margin = margin(0,0.6,0,0.3,unit = "cm"), axis.text.y = element_text(size=9, family="Times"),
          legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),legend.text = element_text(size=8),legend.title = element_text(size=9),legend.key.size = unit(0.3,"cm"))
  print(plot)
  #dev.off()
}

###################################################
#reorder water column level for bar order in figure
BVR.DVM.10Jul19.tows.calcs.long$WaterColumn<-factor(BVR.DVM.10Jul19.tows.calcs.long$WaterColumn,levels=c("epilimnion","metalimnion","hypolimnion"))
BVR.DVM.10Jul19.tows.calcs.long<- BVR.DVM.10Jul19.tows.calcs.long[!(BVR.DVM.10Jul19.tows.calcs.long$metric=="ZoopDensity_No.pL_rep.mean_epi" | BVR.DVM.10Jul19.tows.calcs.long$metric=="ZoopDensity_No.pL_rep.mean_meta"| BVR.DVM.10Jul19.tows.calcs.long$metric=="ZoopDensity_No.pL_rep.mean_hypo" | BVR.DVM.10Jul19.tows.calcs.long$metric=="BiomassConcentration_ugpL_rep.mean_epi"| BVR.DVM.10Jul19.tows.calcs.long$metric=="BiomassConcentration_ugpL_rep.mean_meta" |BVR.DVM.10Jul19.tows.calcs.long$metric=="BiomassConcentration_ugpL_rep.mean_hypo" | substrEnd(BVR.DVM.10Jul19.tows.calcs.long$metric,15)=="percent_density"),] 
BVR.DVM.10Jul19.tows.calcs.long$Taxa <-factor(BVR.DVM.10Jul19.tows.calcs.long$Taxa, levels = c("Cladocera","Rotifera_","Cyclopoid","Calanoida"))

#BVR 10Jul19 epi meta hypo
#jpeg("BVR_epimetahypo_density_10Jul2019.jpg", width = 6, height = 4, units = "in",res = 300)
ggplot(subset(BVR.DVM.10Jul19.tows.calcs.long, grepl("density",metric,ignore.case = TRUE)), aes(x=WaterColumn, y=value)) +
  geom_bar(aes(fill=Taxa, alpha=WaterColumn), stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() +
  facet_wrap(~Taxa, scales= 'free',ncol=5, strip.position = "right", labeller=as_labeller(function(variable,value){c("","","","","Noon")})) +
  theme(strip.text.y = element_text(size = 10,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank(), plot.margin = margin(2.8,0,2,0,unit = "cm")) + 
  scale_x_discrete(name="", labels=rep(c("epi","meta","hypo"),10)) + labs(title="10 Jul 2019 Noon") + 
  scale_fill_manual(values=viridis(6),labels=c("Cladocera","Cyclopoida","Calanoida","Rotifera")) +
  geom_hline(yintercept=0) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
      axis.text.x=element_text(size=7,family="Times"), plot.title = element_text(hjust = 0.5)) + ylab("Density (individual/L)") +
  theme(legend.position = "bottom", legend.margin = margin(0, 1, 2, 1), legend.title = element_text(size=10),legend.key.size = unit(0.3,"cm"),
        legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),legend.text = element_text(size=8),axis.text.y = element_text(size=10, family="Times"))
#dev.off()

#jpeg("BVR_epimetahypo_biomass_10Jul2019.jpg", width = 6, height = 4, units = "in",res = 300)
ggplot(subset(BVR.DVM.10Jul19.tows.calcs.long, grepl("biomass",metric,ignore.case = TRUE)), aes(x=WaterColumn, y=value)) +
  geom_bar(aes(fill=Taxa, alpha=WaterColumn), stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() +
  facet_wrap(~Taxa, scales= 'free',ncol=5, strip.position = "right", labeller=as_labeller(function(variable,value){c("","","","","Noon")})) +
  theme(strip.text.y = element_text(size = 10,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank(), plot.margin = margin(2.8,0,2,0,unit = "cm")) + 
  scale_x_discrete(name="", labels=rep(c("epi","meta","hypo"),10)) + labs(title="10 Jul 2019 Noon") + 
  scale_fill_manual(values=viridis(6),labels=c("Cladocera","Cyclopoida","Calanoida","Rotifera","Total biomass")) +
  geom_hline(yintercept=0) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
      axis.text.x=element_text(size=7,family="Times"), plot.title = element_text(hjust = 0.5)) + ylab(expression(paste("Biomass (",mu,"g/L)"))) +
  theme(legend.position = "bottom", legend.margin = margin(0, 1, 2, 1), legend.title = element_text(size=10),legend.key.size = unit(0.3,"cm"),
      legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),legend.text = element_text(size=8),axis.text.y = element_text(size=10, family="Times"))
#dev.off()
  
  ###############################################################################################
  ###############################################################################################
  #### FCR figures 
      #simple boxplot comparing epi vs hypo in the presence of MOM vs anoxia

  #df for all FCR samples
  FCR_pelagic_DVM<- zoop.repmeans[substr(zoop.repmeans$sample_ID,1,1)=="F",]
  #remove 2019-06-13 because missing full water column data
  FCR_pelagic_DVM <- FCR_pelagic_DVM[!FCR_pelagic_DVM$collect_date== "2019-06-13",]

  #initialize df
  FCR.DVM.calcs<- data.frame("Hour"=unique(FCR_pelagic_DVM$Hour))
  
  #for loop to fill out epi vs hypo calcs 
  variables<- colnames(FCR_pelagic_DVM)[c(5:18)] #change range
  for(i in 1:length(variables)){
    FCR.DVM.calcs[,paste0(variables,"_epi")[i]]<- FCR_pelagic_DVM[substrEnd(FCR_pelagic_DVM$sample_ID,4)!="9.0m" & substrEnd(FCR_pelagic_DVM$sample_ID,4)!="8.0m",paste0(variables)[i]]
    FCR.DVM.calcs[,paste0(variables,"_hypo")[i]] <- FCR_pelagic_DVM[substrEnd(FCR_pelagic_DVM$sample_ID,4)=="9.0m" | substrEnd(FCR_pelagic_DVM$sample_ID,4)=="8.0m",paste0(variables)[i]] - FCR.DVM.calcs[,paste0(variables,"_epi")[i]]
  }
  
  #wide to long df
  FCR.DVM.calcs.long <-  FCR.DVM.calcs %>%
    gather(metric,value, ZoopDensity_No.pL_rep.mean_epi:Calanoida_BiomassConcentration_ugpL_rep.mean_hypo) %>%
    mutate(DateTime = strftime(Hour, "%m-%d-%Y %H:%M"))
  
  #add watercolumn and hour columns
  FCR.DVM.calcs.long$WaterColumn <- ifelse(substrEnd(FCR.DVM.calcs.long$metric,3)=="epi","epilimnion","hypolimnion")
  FCR.DVM.calcs.long$Hour <- ifelse(substrEnd(FCR.DVM.calcs.long$DateTime,5)=="12:00","noon","midnight")
  FCR.DVM.calcs.long$Taxa <- substr(FCR.DVM.calcs.long$metric,1,9)

#unique dates for for loop
dates<- unique(FCR.DVM.calcs.long$DateTime)
  
for(d in 1:length(dates)){
#jpeg("FCR_epivshypo_density.jpg", width = 6, height = 5, units = "in",res = 300)
  plot<- ggplot(subset(FCR.DVM.calcs.long, grepl("density",metric,ignore.case = TRUE) & DateTime %in%
           unique(DateTime)[d]), aes(x=WaterColumn, y=value)) +
    geom_bar(aes(fill=Taxa, alpha=WaterColumn), stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() +
    facet_wrap(Hour~Taxa, scales= 'free',ncol=5, strip.position = "right", labeller=as_labeller(function(variable,value){c("","","","","Noon")})) +
    theme(strip.text.y = element_text(size = 11 ,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank()) + 
    scale_x_discrete(name="", labels=rep(c("epi","hypo"),10)) + labs(title=paste0(dates[d])) + 
    scale_fill_manual("Taxa",values=viridis(6),labels=c("Cladocera","Cyclopoida","Calanoida","Rotifera","Total biomass")) +
    geom_hline(yintercept=0) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
          axis.text.x=element_text(size=9,family="Times"), plot.title = element_text(hjust = 0.5)) + ylab("Density (indiviual/L)") +
    theme(legend.position = "bottom", legend.margin = margin(0, 1, 2, 1),plot.margin = margin(0,0,0,0.3,unit = "cm"), axis.text.y = element_text(size=13, family="Times"),
          legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),legend.text = element_text(size=10),legend.title = element_text(size=12))
  print(plot)
  #dev.off()
}

#reorder taxa levels for figure
FCR.DVM.calcs.long$Taxa<-factor(FCR.DVM.calcs.long$Taxa,levels=c("Cladocera","Cyclopoid","Calanoida","Rotifera_","BiomassCo","ZoopDensi"))

for(d in 1:length(dates)){
  #jpeg("FCR_epivshypo_biomass.jpg", width = 6, height = 5, units = "in",res = 300)
  plot<- ggplot(subset(FCR.DVM.calcs.long, grepl("biomass",metric,ignore.case = TRUE) & DateTime %in%
            unique(DateTime)[d]), aes(x=WaterColumn, y=value)) +
    geom_bar(aes(fill=Taxa, alpha=WaterColumn), stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() +
    facet_wrap(Hour~Taxa, scales= 'free',ncol=5, strip.position = "right", labeller=as_labeller(function(variable,value){c("","","","","Noon")})) +
    theme(strip.text.y = element_text(size = 11 ,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank()) + 
    scale_x_discrete(name="", labels=rep(c("epi","hypo"),10)) + labs(title=paste0(dates[d])) + 
    scale_fill_manual("Taxa",values=viridis(6),labels=c("Cladocera","Cyclopoida","Calanoida","Rotifera","Total biomass")) +
    geom_hline(yintercept=0) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
           axis.text.x=element_text(size=9,family="Times"), plot.title = element_text(hjust = 0.5)) + ylab(expression(paste("Biomass (",mu,"g/L)"))) +
    theme(legend.position = "bottom", legend.margin = margin(0, 1, 2, 1),plot.margin = margin(0,0,0,0.3,unit = "cm"), axis.text.y = element_text(size=13, family="Times"),
          legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),legend.text = element_text(size=10),legend.title = element_text(size=12))
  print(plot)
  #dev.off()
}

###############################################################################################
###############################################################################################
######################### Super simple plot of temp and DO during MSNs

#read in CTD data
CTD1<- read.csv('BVR_CTD_10Jul19.csv', header=TRUE)
CTD2<- read.csv('BVR_CTD_24Jul19.csv', header=TRUE)
  
#order df for plotting purposes
CTD1 <- CTD1 %>% group_by(Date) %>% arrange(Depth_m, by_group=TRUE)
CTD2 <- CTD2 %>% group_by(Date) %>% arrange(Depth_m, by_group=TRUE)

#jpeg("2019campaigns_Temp_O2_profile.jpg", width = 6, height = 5, units = "in",res = 300)
par(mfrow=c(1,2))
par(mar = c(4,0,0,0))
par(oma = c(2,4,4,4))

plot(CTD1$Depth_m[CTD1$time=="10 Jul 11:33"]~CTD1$DO_mgL[CTD1$time=="10 Jul 11:33"], ylim=c(11,0), pch=16, type='o',col="blue",xlim=c(0,7.5), ylab="", xlab="",cex.axis=1.5, cex.lab=1.5)
par(new=TRUE)
plot(CTD1$Depth_m[CTD1$time=="10 Jul 11:33"]~CTD1$Temp_C[CTD1$time=="10 Jul 11:33"],pch=16,type='o',col="red", yaxt='n',xaxt='n',xlab=" ",ylab=" ",ylim=c(11,0),xlim=c(8,28))
axis(3,at=seq(round(min(CTD1$Temp_C),-1),round(max(CTD1$Temp_C),-1),length.out=6), cex.axis=1.5, cex.lab=1.5)
text(15,0,"10 Jul 2019", cex=1.3)

plot(CTD2$Depth_m[CTD2$time=="24 Jul 11:30"]~CTD2$DO_mgL[CTD2$time=="24 Jul 11:30"], ylim=c(11,0), pch=16, type='o',col="blue",xlim=c(0,7.5),ylab="", xlab="",cex.axis=1.5, cex.lab=1.5, yaxt='n')
par(new=TRUE)
plot(CTD2$Depth_m[CTD2$time=="24 Jul 11:30"]~CTD2$Temp_C[CTD2$time=="24 Jul 11:30"],pch=16,type='o',col="red", yaxt='n',xaxt='n',xlab=" ",ylab=" ",ylim=c(11,0),xlim=c(8,28))
axis(3,at=seq(round(min(CTD2$Temp_C),-1),round(max(CTD2$Temp_C),-1),length.out=6),cex.axis=1.5, cex.lab=1.5)
text(15,0,"24 Jul 2019", cex=1.3)
axis(4, at=pretty(CTD2$Depth_m[CTD2$time=="24 Jul 11:30"]),cex.axis=1.5, cex.lab=1.5)
mtext(expression(paste("Temperature (",degree,"C)")),side=3, line=2.4, cex=1.5, outer = TRUE)
mtext("Dissolved Oxygen (mg/L)",side=1, line=-1.8, cex=1.5, outer = TRUE)
mtext("          Depth (m)",side=2, line=2.4, cex=1.5, outer = TRUE)
mtext("          Depth (m)",side=4, line=2.4, cex=1.5, outer = TRUE)
legend("bottomright", legend=c("Dissolved Oxygen","Temperature"), col=c("blue", "red"), 
       cex=1.1, pch=16, box.lty=0,bg="transparent",xjust=1)
#dev.off()
