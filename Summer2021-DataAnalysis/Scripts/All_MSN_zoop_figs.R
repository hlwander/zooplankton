#Zoop density change over 24hrs for all 5 MSNs
#created 16 Oct 2022

#read in libraries
pacman::p_load(dplyr, vegan, labdsv, goeveg, rLakeAnalyzer, ggplot2,tidyr,lubridate, scales, colorblindcheck)

#create function to count characters starting at the end of the string
substrEnd <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#Calculates the standard error####
stderr <- function(x) {
  sd(x,na.rm=TRUE)/sqrt(length(na.omit(x)))
}

#read in zoop data from all 3 years
zoops2019<- read.csv(file.path(getwd(),'Summer2019-DataAnalysis/SummaryStats/FCR_ZooplanktonSummary2019.csv'),header = TRUE)
zoops2020<- read.csv(file.path(getwd(),'Summer2020-DataAnalysis/SummaryStats/FCR_ZooplanktonSummary2020.csv'),header = TRUE)
zoops2021<- read.csv(file.path(getwd(),'Summer2021-DataAnalysis/SummaryStats/FCR_ZooplanktonSummary2021.csv'),header = TRUE)

#drop holopedium becuase only in 2019
zoops2019 <- zoops2019[,!c(grepl("Holopedium",colnames(zoops2019)))]

#combine all zoop datasets
zoop <- rbind(zoops2019,zoops2020,zoops2021)

#make sure sample_ID is class character
zoop$sample_ID<- as.character(zoop$sample_ID)

#change BVR_50_p --> BVR_50
zoop$site_no[zoop$site_no=="BVR_50_p"] <- "BVR_50"

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

#drop 20um samples
zoop <- zoop[!c(substrEnd(zoop$sample_ID,4)=="filt" | substrEnd(zoop$sample_ID,2)=="20"),]

#get hour into character format for grouping
zoop$Hour <- format(round(strptime(paste0(zoop$collect_date, zoop$Hour), format="%Y-%m-%d %H:%M"),units="hours"),format="%H:%M")
#manually change hour of some samples (rounding problems)
zoop$Hour[grepl("midnight",zoop$sample_ID,ignore.case = TRUE)] <- "00:00"
zoop$Hour[grepl("noon",zoop$sample_ID,ignore.case = TRUE)] <- "12:00"
zoop$Hour[grepl("sunrise_h1",zoop$sample_ID,ignore.case = TRUE) |
            grepl("sunrise_epi_h1",zoop$sample_ID,ignore.case = TRUE)] <- "04:00"
zoop$Hour[grepl("sunrise_h2",zoop$sample_ID,ignore.case = TRUE) |
            grepl("sunrise_epi_h2",zoop$sample_ID,ignore.case = TRUE)] <- "05:00"
zoop$Hour[grepl("sunrise_h3",zoop$sample_ID,ignore.case = TRUE) |
            grepl("sunrise_epi_h3",zoop$sample_ID,ignore.case = TRUE)] <- "06:00"
zoop$Hour[grepl("sunrise_h4",zoop$sample_ID,ignore.case = TRUE) |
            grepl("sunrise_epi_h4",zoop$sample_ID,ignore.case = TRUE)] <- "07:00"
zoop$Hour[grepl("sunset_h1",zoop$sample_ID,ignore.case = TRUE) | 
            grepl("sunset_epi_h1",zoop$sample_ID,ignore.case = TRUE)] <- "18:00"
zoop$Hour[grepl("sunset_h2",zoop$sample_ID,ignore.case = TRUE) | 
            grepl("sunset_epi_h2",zoop$sample_ID,ignore.case = TRUE)] <- "19:00"
zoop$Hour[grepl("sunset_h3",zoop$sample_ID,ignore.case = TRUE) |
            grepl("sunset_epi_h3",zoop$sample_ID,ignore.case = TRUE)] <- "20:00"
zoop$Hour[grepl("sunset_h4",zoop$sample_ID,ignore.case = TRUE) |
            grepl("sunset_epi_h4",zoop$sample_ID,ignore.case = TRUE)] <- "21:00"

#drop the 06-29 samples (these were test samples)
zoop <- zoop[substrEnd(zoop$sample_ID,4)!="filt",]

#drop schindler samples and only select BVR_50 epi samples
zoop <- zoop[substrEnd(zoop$site_no,6)!="schind" & zoop$site_no!="BVR_d" & zoop$site_no!="BVR_dam" & zoop$site_no!="BVR_trap" & zoop$site_no!="FCR_50",]

##### Create new df to combine reps over 24 hours
zoop.repmeans <- zoop %>% select(sample_ID,site_no,collect_date,Hour, ZoopDensity_No.pL, BiomassConcentration_ugpL,
                                 Cladocera_density_NopL, Cladocera_BiomassConcentration_ugpL, Cladocera_PercentOfTotal,
                                 Cyclopoida_density_NopL, Cyclopoida_BiomassConcentration_ugpL, Cyclopoida_PercentOfTotal,
                                 Rotifera_density_NopL, Rotifera_BiomassConcentration_ugpL,Rotifera_PercentOfTotal,
                                 Calanoida_PercentOfTotal, Calanoida_density_NopL, Calanoida_BiomassConcentration_ugpL,
                                 Copepoda_PercentOfTotal, Copepoda_density_NopL, Copepoda_BiomassConcentration_ugpL) %>%
  group_by(sample_ID, site_no, Hour, collect_date) %>%
  summarise_at(vars(ZoopDensity_No.pL:Copepoda_BiomassConcentration_ugpL), funs(rep.mean=mean, rep.SE=stderr))

#merge collect_date and hour in a new column
zoop.repmeans$datetime<- paste(zoop.repmeans$collect_date,zoop.repmeans$Hour,sep=" ")
#get times into date format (character here)
zoop.repmeans$datetime<- format(as.POSIXct(zoop.repmeans$datetime,format="%Y-%m-%d %H:%M"), format="%Y-%m-%d %H:%M:%S")
#convert to posixct date format
zoop.repmeans$datetime<- as.POSIXct(zoop.repmeans$datetime, format="%Y-%m-%d %H:%M")

zoop.repmeans$dates <- ifelse(zoop.repmeans$collect_date=="2019-07-10" | zoop.repmeans$collect_date=="2019-07-24"| 
                                zoop.repmeans$collect_date=="2020-08-12" | zoop.repmeans$collect_date=="2021-06-15" |
                                zoop.repmeans$collect_date=="2021-07-07","2022-10-15","2022-10-16")

#only select hour and then add arbitrary dates for plotting
zoop.repmeans$Hour <- format(zoop.repmeans$Hour, format='%H:%M')

#combine hour and date
zoop.repmeans$Hour <- strptime(paste0(as.character(zoop.repmeans$dates), zoop.repmeans$Hour),format="%Y-%m-%d %H:%M")
zoop.repmeans$Hour <- as.POSIXct(zoop.repmeans$Hour)

#make sure zoop.repmeans is a dataframe
zoop.repmeans <- data.frame(zoop.repmeans)

#order by hour for plotting
zoop.repmeans <- zoop.repmeans[order(zoop.repmeans$Hour),]

#Export all zoop data df
write.csv(zoop.repmeans,"./Summer2021-DataAnalysis/SummaryStats/All_MSN_zoops.csv",row.names = FALSE)

#only select epi samples for DHM plots
zoop_epi <- zoop.repmeans[grepl("epi",zoop.repmeans$sample_ID) |grepl("sunrise",zoop.repmeans$sample_ID) | grepl("sunset",zoop.repmeans$sample_ID) | zoop.repmeans$site_no=="BVR_l",]

#convert new dfs from tibble to dataframe 
zoop_DHM <- data.frame(zoop_epi)

#convert df from wide to long (kinda hacky way bc having problems doing this)
df1 <- zoop_DHM %>% gather(metric,value,ZoopDensity_No.pL_rep.mean:Copepoda_BiomassConcentration_ugpL_rep.mean)
df2 <- zoop_DHM %>% gather(metric.SE,value.SE, ZoopDensity_No.pL_rep.SE:Copepoda_BiomassConcentration_ugpL_rep.SE)

##cut and paste to merge df
zoop_DHM_long <- df1[,c(1:4,24,25)]
#zoop_DHM_long$metric.SE <- df2$metric.SE #use this as a check to make sure rows match up
zoop_DHM_long$value.SE <- df2$value.SE

#Export DHM csv
write.csv(zoop_DHM_long,"./Summer2021-DataAnalysis/SummaryStats/All_MSN_DHM.csv",row.names = FALSE)

#reset DHM df so can only look at density
zoop_DHM_long <- NA

variables <- c("ZoopDensity_No.pL_rep.mean","Cladocera_density_NopL_rep.mean", "Cladocera_PercentOfTotal_rep.mean",
               "Cyclopoida_density_NopL_rep.mean","Cyclopoida_PercentOfTotal_rep.mean","Rotifera_density_NopL_rep.mean",
               "Rotifera_PercentOfTotal_rep.mean","Calanoida_density_NopL_rep.mean","Calanoida_PercentOfTotal_rep.mean",
               "Copepoda_density_NopL_rep.mean","Copepoda_PercentOfTotal_rep.mean")
SE <- c("ZoopDensity_No.pL_rep.SE","Cladocera_density_NopL_rep.SE", "Cladocera_PercentOfTotal_rep.SE",
        "Cyclopoida_density_NopL_rep.SE","Cyclopoida_PercentOfTotal_rep.SE","Rotifera_density_NopL_rep.SE",
        "Rotifera_PercentOfTotal_rep.SE","Calanoida_density_NopL_rep.SE","Calanoida_PercentOfTotal_rep.SE",
        "Copepoda_density_NopL_rep.SE","Copepoda_PercentOfTotal_rep.SE")

#select density columns
zoop_DHM <- zoop_DHM[,-c(which(grepl("ug_rep",colnames(zoop_DHM))))]

#remove biomass too
zoop_DHM <- zoop_DHM[,-c(which(grepl("ugpL",colnames(zoop_DHM))))]

#convert df from wide to long (kinda hacky way bc having problems doing this)
df1 <- zoop_DHM %>% gather(metric,value,all_of(variables))
df2 <- zoop_DHM %>% gather(metric.SE,value.SE, all_of(SE))

#cut and paste to merge df
zoop_DHM_long <- df1[,c(1:4,16,18:19)]
#zoop_DHM_long$metric.SE <- df2$metric.SE #use this as a check to make sure rows match up
zoop_DHM_long$value.SE <- df2$value.SE

#drop _rep.mean from all metric names
zoop_DHM_long$metric <- substr(zoop_DHM_long$metric,1,nchar(zoop_DHM_long$metric)-9)

#add column for MSN #
zoop_DHM_long$MSN <- ifelse(zoop_DHM_long$collect_date=="2019-07-10" | zoop_DHM_long$collect_date=="2019-07-11",1,
                            ifelse(zoop_DHM_long$collect_date=="2019-07-24" | zoop_DHM_long$collect_date=="2019-07-25",2,
                                   ifelse(zoop_DHM_long$collect_date=="2020-08-12" | zoop_DHM_long$collect_date=="2020-08-13",3,
                                          ifelse(zoop_DHM_long$collect_date=="2021-06-15" | zoop_DHM_long$collect_date=="2021-06-16",4,5))))

metric_taxa <-c("ZoopDensity","Cladocera","Cladocera","Cyclopoida",
                "Cyclopoida", "Rotifera", "Rotifera", "Calanoida",
                "Calanoida","Copepoda","Copepoda")
names(metric_taxa) <- c(unique(zoop_DHM_long$metric))

sites <- c("Pelagic","Littoral")
names(sites) <- c("BVR_50","BVR_l")

#cb_friendly_2 <- c("#8C510A", "#BF812D","#DFC27D", "#C7EAE5", "#35978F")

#Figure for zoop density for each MSN 24-hours 
ggplot(subset(zoop_DHM_long, metric %in% c("Cladocera_density_NopL","Cyclopoida_density_NopL","Rotifera_density_NopL")),
                aes(Hour,value, color=as.factor(MSN))) + 
  geom_rect(aes(xmin=as.POSIXct("2022-10-15 11:30:00"),xmax=as.POSIXct("2022-10-15 20:41:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2022-10-15 20:42:00"),xmax=as.POSIXct("2022-10-16 06:10:00"), ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2022-10-16 06:11:00"),xmax=as.POSIXct("2022-10-16 12:30:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + facet_grid(site_no~metric,scales="free_y",labeller = labeller(metric=metric_taxa, site_no=sites)) + xlab("")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=8), axis.text = element_text(size=7, color="black"), legend.background = element_blank(), legend.key = element_blank(), legend.key.height=unit(0.3,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_rect(fill = "transparent"), legend.position = c(0.1,0.92), legend.spacing = unit(-0.5, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), legend.key.width =unit(0.7,"line"))+ scale_x_datetime(expand = c(0,0),labels = date_format("%H-%M",tz="EST5EDT")) +
  #scale_colour_viridis_d("",option="viridis", labels=c("10-11 Jul 2019","24-25 Jul 2019","12-13 Aug 2020","15-16 Jun 2021","7-8 Jul 2021")) + 
  scale_color_manual("",values=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"), labels=c("10-11 Jul 2019","24-25 Jul 2019","12-13 Aug 2020","15-16 Jun 2021","7-8 Jul 2021"), guide=guide_legend(order=1)) + 
  geom_line()+ ylab("Density (Individuals/L)") + scale_fill_manual("",values=c("#CCCCCC","white"), guide = "none")+
  geom_errorbar(aes(ymin=value-value.SE, ymax=value+value.SE), width=.2,position=position_dodge(.9))
ggsave(file.path(getwd(),"Summer2021-DataAnalysis/Figures/BVR_MSNs_taxa_density.jpg"), width=5, height=3) 

ggplot(subset(zoop_DHM_long, metric %in% c("Cladocera_PercentOfTotal","Cyclopoida_PercentOfTotal","Rotifera_PercentOfTotal")), #removing calanoids because they are super low %
                aes(Hour,value, color=as.factor(MSN))) + 
  geom_rect(aes(xmin=as.POSIXct("2022-10-15 11:30:00"),xmax=as.POSIXct("2022-10-15 20:41:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2022-10-15 20:42:00"),xmax=as.POSIXct("2022-10-16 06:10:00"), ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2022-10-16 06:11:00"),xmax=as.POSIXct("2022-10-16 12:30:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + facet_grid(site_no~metric,scales="free_y",labeller = labeller(metric=metric_taxa, site_no=sites)) + xlab("")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=10), axis.text = element_text(size=9, color="black"), legend.background = element_blank(), legend.key = element_blank(), legend.key.height=unit(0.3,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_rect(fill = "transparent"), legend.position = c(0.12,0.35), legend.spacing = unit(-0.5, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), legend.key.width =unit(0.7,"line"),legend.spacing.y=unit(-1,'cm'), legend.spacing.x=unit(0,'cm')) + 
  scale_x_datetime(expand = c(0,0),labels = date_format("%H-%M",tz="EST5EDT"))+
  #scale_colour_viridis_d("",option="viridis", labels=c("10-11 Jul 2019","24-25 Jul 2019","12-13 Aug 2020","15-16 Jun 2021","7-8 Jul 2021")) + 
  #scale_color_manual("",values=hcl.colors(5,"Geyser"), labels=c("10-11 Jul 2019","24-25 Jul 2019","12-13 Aug 2020","15-16 Jun 2021","7-8 Jul 2021"), guide=guide_legend(order=1)) + 
  scale_color_manual("",values=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"), labels=c("10-11 Jul 2019","24-25 Jul 2019","12-13 Aug 2020","15-16 Jun 2021","7-8 Jul 2021"), guide=guide_legend(order=1)) + 
  geom_line()+ ylab("% Density") + scale_fill_manual("",values=c("#CCCCCC","white"), guide = "none")+ 
  geom_errorbar(aes(ymin=value-value.SE, ymax=value+value.SE), width=.2,position=position_dodge(.9))
ggsave(file.path(getwd(),"Summer2021-DataAnalysis/Figures/BVR_MSNs_taxa_percent_density.jpg"), width=5, height=4) 

#-------------------------------------------------------------------------------------#
# DVM calcs for epi and hypo density across all 5 campaigns

#read in all 3 DVM tables from 2019-2021
DVM_2019 <- read.csv(file.path(getwd(),'Summer2019-DataAnalysis/SummaryStats/DVM_2019_zoops.csv'),header = TRUE)
DVM_2020 <- read.csv(file.path(getwd(),'Summer2020-DataAnalysis/SummaryStats/DVM_2020_zoops.csv'),header = TRUE) %>% select(-X)
DVM_2021 <- read.csv(file.path(getwd(),'Summer2021-DataAnalysis/SummaryStats/DVM_2021_zoops.csv'),header = TRUE) %>% select(-X)

#merge all zoop files
DVM_all <- rbind(DVM_2019,DVM_2020,DVM_2021)


