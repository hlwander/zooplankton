#multivariate stats for 2019 - 2021 zoop data
#created 25Nov2021

#read in libraries
pacman::p_load(dplyr, vegan, labdsv, goeveg, rLakeAnalyzer, ggplot2,tidyr,viridis)

#function to count characters starting at the end of the string
substrEnd <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#function to find the closest value in a df
closest<-function(x,y){
  x[which(abs(x-y)==min(abs(x-y)))] }

#read in zoop data from all 3 years
zoops2019<- read.csv(file.path(getwd(),'Summer2019-DataAnalysis/SummaryStats/FCR_ZooplanktonSummary2019.csv'),header = TRUE)
zoops2020<- read.csv(file.path(getwd(),'Summer2020-DataAnalysis/SummaryStats/FCR_ZooplanktonSummary2020.csv'),header = TRUE)
zoops2021<- read.csv(file.path(getwd(),'Summer2021-DataAnalysis/SummaryStats/FCR_ZooplanktonSummary2021.csv'),header = TRUE)

#select density cols to keep
zoops2019 <- zoops2019 %>% select("sample_ID","site_no","collect_date","DepthOfTow_m","Hour","mesh_size_μm","ZoopDensity_No.pL","Calanoida_density_NopL",
                                   "Cyclopoida_density_NopL", "Keratella_density_NopL","Kellicottia_density_NopL", "Bosmina_density_NopL",
                                   "Daphnia_density_NopL", "Ceriodaphnia_density_NopL","nauplius_density_NopL", "Collothecidae_density_NopL",
                                   "Synchaetidae_density_NopL", "Conochilidae_density_NopL")

zoops2020 <- zoops2020 %>% select("sample_ID","site_no","collect_date","DepthOfTow_m","Hour","mesh_size_μm","ZoopDensity_No.pL","Calanoida_density_NopL",
                                  "Cyclopoida_density_NopL", "Keratella_density_NopL","Kellicottia_density_NopL", "Bosmina_density_NopL",
                                  "Daphnia_density_NopL", "Ceriodaphnia_density_NopL","nauplius_density_NopL", "Collothecidae_density_NopL",
                                  "Synchaetidae_density_NopL", "Conochilidae_density_NopL")

zoops2021 <- zoops2021 %>% select("sample_ID","site_no","collect_date","DepthOfTow_m","Hour","mesh_size_μm","ZoopDensity_No.pL","Calanoida_density_NopL",
                                  "Cyclopoida_density_NopL", "Keratella_density_NopL","Kellicottia_density_NopL", "Bosmina_density_NopL",
                                  "Daphnia_density_NopL", "Ceriodaphnia_density_NopL","nauplius_density_NopL", "Collothecidae_density_NopL",
                                  "Synchaetidae_density_NopL", "Conochilidae_density_NopL")

#combine all zoop datasets
zoops <- rbind(zoops2019,zoops2020,zoops2021)

#ignore 20 um and horizontal trap samples
zoops <- zoops %>% filter(mesh_size_μm >20, na.rm=TRUE)

#create df for temporal epi tows
zoop_epi_tows <- zoops[zoops$site_no!="FCR_50"& zoops$site_no!="BVR_d" & zoops$site_no!="BVR_dam" & (grepl("epi",zoops$sample_ID) |grepl("sunrise",zoops$sample_ID) | 
                        grepl("sunset",zoops$sample_ID) | zoops$site_no=="BVR_l"), ] %>%
  mutate(Hour=substr(Hour,1,2)) %>% 
  select(!c(DepthOfTow_m, sample_ID)) %>% group_by(site_no,collect_date,Hour) %>%
  summarise(across(everything(),list(mean)))

zoop_epi_tows$time <-ifelse(zoop_epi_tows$Hour=="12" | zoop_epi_tows$Hour=="11", "noon", ifelse(
   zoop_epi_tows$Hour =="0:" | zoop_epi_tows$Hour =="23", "midnight",ifelse(zoop_epi_tows$Hour=="18"|
   zoop_epi_tows$Hour=="19" | zoop_epi_tows$Hour=="20" | zoop_epi_tows$Hour=="21", "sunset", "sunrise")))

zoop_epi_tows$site <- ifelse(substrEnd(zoop_epi_tows$site_no,1)=="l","lit","pel")

#also add columns to group times by sampling event and time and sites together
zoop_epi_tows$groups <- ifelse(zoop_epi_tows$collect_date=="2019-07-10" | zoop_epi_tows$collect_date=="2019-07-11","1",
                        ifelse(zoop_epi_tows$collect_date=="2019-07-24" | zoop_epi_tows$collect_date=="2019-07-25","2",
                        ifelse(zoop_epi_tows$collect_date=="2020-08-12" | zoop_epi_tows$collect_date=="2020-08-13","3",
                        ifelse(zoop_epi_tows$collect_date=="2021-06-15" | zoop_epi_tows$collect_date=="2021-06-16","4","5"))))

zoop_epi_tows$timesite <- paste0(zoop_epi_tows$time,zoop_epi_tows$site)

zoop_epi_tows$timegroup <- paste0(zoop_epi_tows$time,zoop_epi_tows$groups)

#------------------------------------------------------------------------------#
# set up data for NMDS
#relies on rank orders for ordination, no assumptions of linear relationship
zoop_temporal_dens <- zoop_epi_tows[,c(grepl("density_NopL",colnames(zoop_epi_tows)))]  

#transforming data - hellinger transformation because gives low weight to low/zero values
#converts species abundances from absolute to relative - use w/ bray curtis (jaccard might be another one to try later on to look at absolute)
zoop_temporal_dens_trans <- hellinger(zoop_temporal_dens)

#-------------------------------------------------------------------------------#
#                  Bray-curtis dissimilarity --> NMDS figs                      #
#-------------------------------------------------------------------------------#
#scree plot to choose dimension #
dimcheckMDS(zoop_temporal_dens_trans, distance = "bray", k = 6, trymax = 20, autotransform = TRUE)

#choosing 4 dimensions (0.07) so stress is between 0.05 and 0.1
NMDS_temporal_bray <- metaMDS(zoop_temporal_dens_trans, distance='bray', k=4, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_temporal_bray$stress

#  Points divided into time groups
#jpeg("Figures/2019-2021_NMDS_1v2_bray_days.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_temporal_bray,display = c('sites','species'),choices = c(1,2),type = "n") 
points(NMDS_temporal_bray$points[zoop_epi_tows$groups=="1",1], NMDS_temporal_bray$points[zoop_epi_tows$groups=="1",2], pch=21,bg="#008585")
points(NMDS_temporal_bray$points[zoop_epi_tows$groups=="2",1], NMDS_temporal_bray$points[zoop_epi_tows$groups=="2",2], pch=21,bg="#9BBAA0")
points(NMDS_temporal_bray$points[zoop_epi_tows$groups=="3",1], NMDS_temporal_bray$points[zoop_epi_tows$groups=="3",2], pch=21,bg="#F2E2B0")
points(NMDS_temporal_bray$points[zoop_epi_tows$groups=="4",1], NMDS_temporal_bray$points[zoop_epi_tows$groups=="4",2], pch=21,bg="#DEA868")
points(NMDS_temporal_bray$points[zoop_epi_tows$groups=="5",1], NMDS_temporal_bray$points[zoop_epi_tows$groups=="5",2], pch=21,bg="#C7522B")
legend("bottomright", legend=c('Day1','Day2','Day3','Day4','Day5'), pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"), bty = "n") 
ordihull(ord, zoop_epi_tows$groups, display = "sites", draw = c("polygon"),
         col = c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"), alpha = 75,cex = 2)
text(NMDS_temporal_bray$species[,2],NMDS_temporal_bray$species[,4], labels = c("Calanoida","Cyclopoida","Keratella","Kellicottia",
     "Bosmina","Daphnia","Ceriodaphnia","nauplius","Collothecidae","Synchaetidae","Conochilidae"), cex=0.9)
#dev.off()

#jpeg("Figures/2019-2021_NMDS_1v3_bray_days.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_temporal_bray,display = c('sites','species'),choices = c(1,3),type = "n") 
points(NMDS_temporal_bray$points[zoop_epi_tows$groups=="1",1], NMDS_temporal_bray$points[zoop_epi_tows$groups=="1",3], pch=21,bg="#008585")
points(NMDS_temporal_bray$points[zoop_epi_tows$groups=="2",1], NMDS_temporal_bray$points[zoop_epi_tows$groups=="2",3], pch=21,bg="#9BBAA0")
points(NMDS_temporal_bray$points[zoop_epi_tows$groups=="3",1], NMDS_temporal_bray$points[zoop_epi_tows$groups=="3",3], pch=21,bg="#F2E2B0")
points(NMDS_temporal_bray$points[zoop_epi_tows$groups=="4",1], NMDS_temporal_bray$points[zoop_epi_tows$groups=="4",3], pch=21,bg="#DEA868")
points(NMDS_temporal_bray$points[zoop_epi_tows$groups=="5",1], NMDS_temporal_bray$points[zoop_epi_tows$groups=="5",3], pch=21,bg="#C7522B")
legend("bottomright", legend=c('Day1','Day2','Day3','Day4','Day5'), pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"), bty = "n") 
ordihull(ord, zoop_epi_tows$groups, display = "sites", draw = c("polygon"),
         col = c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"), alpha = 75,cex = 2)
text(NMDS_temporal_bray$species[,1],NMDS_temporal_bray$species[,3], labels = c("Calanoida","Cyclopoida","Keratella","Kellicottia",
                                                                               "Bosmina","Daphnia","Ceriodaphnia","nauplius","Collothecidae","Synchaetidae","Conochilidae"), cex=0.9)
#dev.off()

#jpeg(file.path(getwd(),"Summer2021-DataAnalysis/Figures/2019-2021_NMDS_1v2_bray_sites.jpg"), width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_temporal_bray,display = c('sites','species'),choices = c(1,2),type = "n") 
points(NMDS_temporal_bray$points[zoop_epi_tows$site=="pel",1], NMDS_temporal_bray$points[zoop_epi_tows$site=="pel",2], pch=21,bg="slateblue4")
points(NMDS_temporal_bray$points[zoop_epi_tows$site=="lit",1], NMDS_temporal_bray$points[zoop_epi_tows$site=="lit",2], pch=21,bg="palegreen3")
legend("bottomright", legend=c("Pelagic","Littoral"), pch=21, pt.bg=c("slateblue4","palegreen3"),bty = "n") 
ordihull(ord, zoop_epi_tows$site, display = "sites", draw = c("polygon"),
         col = c("palegreen3","slateblue4"), alpha = 75,cex = 2)
text(NMDS_temporal_bray$species[,1],NMDS_temporal_bray$species[,2], labels = c("Calanoida","Cyclopoida","Keratella","Kellicottia",
     "Bosmina","Daphnia","Ceriodaphnia","nauplius","Collothecidae","Synchaetidae","Conochilidae"), cex=0.9)
#dev.off()

#jpeg("Figures/2019-2021_NMDS_1v3_bray_sites.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_temporal_bray,display = c('sites','species'),choices = c(1,3),type = "n") 
points(NMDS_temporal_bray$points[zoop_epi_tows$site=="pel",1], NMDS_temporal_bray$points[zoop_epi_tows$site=="pel",3], pch=21,bg="slateblue4")
points(NMDS_temporal_bray$points[zoop_epi_tows$site=="lit",1], NMDS_temporal_bray$points[zoop_epi_tows$site=="lit",3], pch=21,bg="palegreen3")
legend("bottomright", legend=c("Pelagic","Littoral"), pch=21, pt.bg=c("slateblue4","palegreen3"),bty = "n") 
ordihull(ord, zoop_epi_tows$site, display = "sites", draw = c("polygon"),
         col = c("palegreen3","slateblue4"), alpha = 75,cex = 2)
text(NMDS_temporal_bray$species[,1],NMDS_temporal_bray$species[,3], labels = c("Calanoida","Cyclopoida","Keratella","Kellicottia",
     "Bosmina","Daphnia","Ceriodaphnia","nauplius","Collothecidae","Synchaetidae","Conochilidae"), cex=0.9)
#dev.off()


#jpeg("Figures/2019-2021_NMDS_1v4_bray_time.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_temporal_bray,display = c('sites','species'),choices = c(1,4),type = "n") 
points(NMDS_temporal_bray$points[zoop_epi_tows$time=='sunrise',1], NMDS_temporal_bray$points[zoop_epi_tows$time=='sunrise',4], pch=21,bg="#CC79A7")
points(NMDS_temporal_bray$points[zoop_epi_tows$time=='noon',1], NMDS_temporal_bray$points[zoop_epi_tows$time=='noon',4], pch=21,bg="#F0E442")
points(NMDS_temporal_bray$points[zoop_epi_tows$time=='sunset',1], NMDS_temporal_bray$points[zoop_epi_tows$time=='sunset',4], pch=21,bg="#009E73")
points(NMDS_temporal_bray$points[zoop_epi_tows$time=='midnight',1], NMDS_temporal_bray$points[zoop_epi_tows$time=='midnight',4], pch=21,bg="#0072B2")
legend("bottomright", legend=c('sunrise','noon', 'sunset', 'midnight'), pch=21, pt.bg=c('#CC79A7', '#F0E442','#009E73','#0072B2'),bty = "n") 
ordihull(ord, zoop_epi_tows$time, display = "sites", draw = c("polygon"),
         col = c("#0072B2", "#F0E442","#CC79A7","#009E73"), alpha = 75,cex = 2)
text(NMDS_temporal_bray$species[,1],NMDS_temporal_bray$species[,4], labels = c("Calanoida","Cyclopoida","Keratella","Kellicottia",
     "Bosmina","Daphnia","Ceriodaphnia","nauplius","Collothecidae","Synchaetidae","Conochilidae"), cex=0.9)
#dev.off()

#jpeg("Figures/2019-2021_NMDS_1v3_bray_time.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_temporal_bray,display = c('sites','species'),choices = c(1,3),type = "n") 
points(NMDS_temporal_bray$points[zoop_epi_tows$time=='sunrise',1], NMDS_temporal_bray$points[zoop_epi_tows$time=='sunrise',3], pch=21,bg="#CC79A7")
points(NMDS_temporal_bray$points[zoop_epi_tows$time=='noon',1], NMDS_temporal_bray$points[zoop_epi_tows$time=='noon',3], pch=21,bg="#F0E442")
points(NMDS_temporal_bray$points[zoop_epi_tows$time=='sunset',1], NMDS_temporal_bray$points[zoop_epi_tows$time=='sunset',3], pch=21,bg="#009E73")
points(NMDS_temporal_bray$points[zoop_epi_tows$time=='midnight',1], NMDS_temporal_bray$points[zoop_epi_tows$time=='midnight',3], pch=21,bg="#0072B2")
legend("bottomright", legend=c('sunrise','noon', 'sunset', 'midnight'), pch=21, pt.bg=c('#CC79A7', '#F0E442','#009E73','#0072B2'),bty = "n") 
ordihull(ord, zoop_epi_tows$time, display = "sites", draw = c("polygon"),
         col = c("#0072B2", "#F0E442","#CC79A7","#009E73"), alpha = 75,cex = 2)
text(NMDS_temporal_bray$species[,1],NMDS_temporal_bray$species[,3], labels = c("Calanoida","Cyclopoida","Keratella","Kellicottia",
     "Bosmina","Daphnia","Ceriodaphnia","nauplius","Collothecidae","Synchaetidae","Conochilidae"), cex=0.9)
#dev.off()

#jpeg("Figures/2019-2021_NMDS_1v2_bray_timegroups.jpg", width = 6, height = 5, units = "in",res = 300)
#ord <- ordiplot(NMDS_temporal_bray,display = c('sites','species'),choices = c(1,2),type = "n") 
#points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunrise1',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunrise1',2], pch=21,bg="cornsilk")
#points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunrise2',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunrise2',2], pch=21,bg="yellow")
#points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunrise3',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunrise3',2], pch=21,bg="yellow3")
#points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='noon1',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='noon1',2], pch=21,bg="lightpink")
#points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='noon2',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='noon2',2], pch=21,bg="indianred")
#points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='noon3',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='noon3',2], pch=21,bg="red")
#points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunset1',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunset1',2], pch=21,bg="palegreen")
#points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunset2',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunset2',2], pch=21,bg="mediumseagreen")
#points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunset3',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunset3',2], pch=21,bg="darkgreen")
#points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='midnight1',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='midnight1',2], pch=21,bg="lightsteelblue")
#points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='midnight2',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='midnight2',2], pch=21,bg="skyblue")
#points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='midnight3',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='midnight3',2], pch=21,bg="steelblue")
#ordiellipse(ord, zoop_epi_tows$timegroup, display = "sites", kind="sd",draw="lines", col = c("lightsteelblue", "skyblue","steelblue","lightpink","indianred","red","cornsilk","yellow","yellow3","palegreen","mediumseagreen","darkgreen"), alpha = 75,cex = 2)
#legend("bottomright", legend=c('sunrise1', 'sunrise2', 'sunrise3', 'noon1', 'noon2', 'noon3', 'sunset1', 'sunset2', 'sunset3', 'midnight1', 'midnight2', 'midnight3'), pch=21,
#       pt.bg=c('cornsilk', 'yellow','yellow3','lightpink','indianred','red','palegreen','mediumseagreen','darkgreen','lightsteelblue','skyblue','steelblue'),bty = "n") 
#dev.off()


#-------------------------------------------------------------------------------#
#           Averaging zoops by time and campaign/day for new NMDS               #
#-------------------------------------------------------------------------------#
#first average times for each 24-hour campaign so there are 11 points per day (basically just averaging noon and midnight)
zoop_epi_tows$order <- ifelse(zoop_epi_tows$Hour=="11" | zoop_epi_tows$Hour=="12",1, 
                              ifelse(zoop_epi_tows$Hour=="18",2, ifelse(zoop_epi_tows$Hour=="19",3,
                              ifelse(zoop_epi_tows$Hour=="20",4, ifelse(zoop_epi_tows$Hour=="21",5,
                              ifelse(zoop_epi_tows$Hour=="0:" | zoop_epi_tows$Hour=="23",6,
                              ifelse(zoop_epi_tows$Hour=="4:" | zoop_epi_tows$Hour=="3:",7,
                              ifelse(zoop_epi_tows$Hour=="5:",8,
                              ifelse(zoop_epi_tows$Hour=="6:",9,10)))))))))

#add order 11 for noon2
zoop_epi_tows$order[zoop_epi_tows$order==1 & (zoop_epi_tows$collect_date=="2019-07-10" | zoop_epi_tows$collect_date=="2019-07-24" |
                                                zoop_epi_tows$collect_date=="2020-08-12" | zoop_epi_tows$collect_date=="2021-06-15" |
                                                zoop_epi_tows$collect_date=="2021-07-07")] <- 11

#now specify whether it is noon1 or noon2
zoop_epi_tows$time[zoop_epi_tows$order==1] <- "noon1"
zoop_epi_tows$time[zoop_epi_tows$order==11] <- "noon2"

#average by MSN, site, then hour
zoop_avg <- zoop_epi_tows %>% group_by(groups,site,order) %>%
  summarise_at(vars(Calanoida_density_NopL_1:Conochilidae_density_NopL_1), list(mean = mean))

#time df
zoop_avg_time <- zoop_epi_tows %>% group_by(order) %>%
  summarise_at(vars(Calanoida_density_NopL_1:Conochilidae_density_NopL_1), list(mean = mean))

#df with only noon and midnight
zoop_avg_noon_midnight <- zoop_epi_tows %>% group_by(groups,site,order) %>% filter(order %in% c(1,6,11)) %>%
  summarise_at(vars(Calanoida_density_NopL_1:Conochilidae_density_NopL_1), list(mean = mean))

#pelagic vs littoral dfs
zoop_pel <- zoop_avg[zoop_avg$site=="pel",]
zoop_lit <- zoop_avg[zoop_avg$site=="lit",]


#only select data cols
zoop_temporal_avg_dens <- zoop_avg[,c(grepl("mean",colnames(zoop_avg)))] 
zoop_avg_time_dens <- zoop_avg_time[,c(grepl("mean",colnames(zoop_avg_time)))] 
zoop_avg_noon_midnight_dens <- zoop_avg_noon_midnight[,c(grepl("mean",colnames(zoop_avg_noon_midnight)))] 
zoop_pel_dens <- zoop_pel[,c(grepl("mean",colnames(zoop_pel)))]
zoop_lit_dens <- zoop_lit[,c(grepl("mean",colnames(zoop_lit)))]

#transforming data - hellinger transformation because gives low weight to low/zero values
#converts species abundances from absolute to relative - use w/ bray curtis
zoop_temporal_dens_avg_trans <- hellinger(zoop_temporal_avg_dens)
zoop_avg_time_dens_trans <- hellinger(zoop_avg_time_dens)
zoop_avg_noon_midnight_dens_trans <- hellinger(zoop_avg_noon_midnight_dens)
zoop_pel_dens_trans <- hellinger(zoop_pel_dens)
zoop_lit_dens_trans <- hellinger(zoop_lit_dens)


#scree plot to choose dimension #
dimcheckMDS(zoop_temporal_dens_avg_trans, distance = "bray", k = 6, trymax = 20, autotransform = TRUE)
dimcheckMDS(zoop_avg_time_dens_trans, distance = "bray", k = 6, trymax = 20, autotransform = TRUE)

#now do NMDS using averages w/ 4 dimensions for consistency
NMDS_temporal_avg_bray <- metaMDS(zoop_temporal_dens_avg_trans, distance='bray', k=4, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_temporal_avg_bray$stress

#only doing 2 dimensions here because not enough data (and stress is already low)
NMDS_avg_time_bray <- metaMDS(zoop_avg_time_dens_trans, distance='bray', k=2, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_avg_time_bray$stress

NMDS_noon_midnight_avg_bray <- metaMDS(zoop_avg_noon_midnight_dens_trans, distance='bray', k=4, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_noon_midnight_avg_bray$stress

NMDS_pel_bray <- metaMDS(zoop_pel_dens_trans, distance='bray', k=4, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_pel_bray$stress

NMDS_lit_bray <- metaMDS(zoop_lit_dens_trans, distance='bray', k=4, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_lit_bray$stress

#-------------------------------------------------------------------------------#
#                                 Tracking figs                                 #
#-------------------------------------------------------------------------------#
#jpeg(file.path(getwd(),"Summer2021-DataAnalysis/Figures/2019-2021_NMDS_1v2_bray_temporal_avg_days.jpg"), width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_temporal_avg_bray,display = c('sites','species'),choices = c(1,2),type = "n")
ordihull(ord, zoop_avg$groups, display = "sites", draw = c("polygon"),
         col = c("#008585","#89B199","#EFECBF","#DB9B5A","#C7522B"), alpha = 75,cex = 2)
lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="1" & zoop_avg$site=="pel",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="1" & zoop_avg$site=="pel",2], col="#008585")
lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="1" & zoop_avg$site=="lit",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="1" & zoop_avg$site=="lit",2], col="#369187")
points(NMDS_temporal_avg_bray$points[zoop_avg$groups=="1",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="1",2], pch=rep(c(0,1,2,3,4,5,6,7,8,9,10),2),col=c(rep("#369187",10),rep("#008585",10)))

lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="2" & zoop_avg$site=="pel",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="2" & zoop_avg$site=="pel",2], col="#89B199")
lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="2" & zoop_avg$site=="lit",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="2" & zoop_avg$site=="lit",2], col="#C4D5B2")
points(NMDS_temporal_avg_bray$points[zoop_avg$groups=="2",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="2",2], pch=rep(c(0,1,2,3,4,5,6,7,8,9,10),2),col=c(rep("#C4D5B2",10),rep("#89B199",10)))

lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="3" & zoop_avg$site=="pel",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="3" & zoop_avg$site=="pel",2], col="#EFECBF")
lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="3" & zoop_avg$site=="lit",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="3" & zoop_avg$site=="lit",2], col="#EFDBA7")
points(NMDS_temporal_avg_bray$points[zoop_avg$groups=="3",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="3",2], pch=rep(c(0,1,2,3,4,5,6,7,8,9,10),2),col=c(rep("#EFDBA7",10),rep("#EFECBF",10)))

lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="4" & zoop_avg$site=="pel",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="4" & zoop_avg$site=="pel",2], col="#DB9B5A")
lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="4" & zoop_avg$site=="lit",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="4" & zoop_avg$site=="lit",2], col="#E4BC80")
points(NMDS_temporal_avg_bray$points[zoop_avg$groups=="4",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="4",2], pch=rep(c(0,1,2,3,4,5,6,7,8,9,10),2),col=c(rep("#E4BC80",10),rep("#DB9B5A",10)))

lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="5" & zoop_avg$site=="pel",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="5" & zoop_avg$site=="pel",2], col="#C7522B")
lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="5" & zoop_avg$site=="lit",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="5" & zoop_avg$site=="lit",2], col="#CA602E")
points(NMDS_temporal_avg_bray$points[zoop_avg$groups=="5",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="5",2], pch=rep(c(0,1,2,3,4,5,6,7,8,9,10),2),col=c(rep("#CA602E",10),rep("#C7522B",10)))

legend("bottomright", legend=c('pelagic day1','littoral day1','pelagic day2','littoral day2','pelagic day3','littoral day3', "pelagic day4", "littoral day4", "pelagic day5", "littoral day5"), pch=21, 
       pt.bg=c("#008585","#369187","#89B199","#C4D5B2","#EFECBF","#EFDBA7","#DB9B5A","#E4BC80", "#C7522B","#CA602E"),bty = "n",cex=0.8) 
legend("bottomleft", legend=c('noon1','sunset1','sunset2','sunset3','sunset4','midnight','sunrise1','sunrise2','sunrise3','sunrise4','noon2'), pch=c(0,1,2,3,4,5,6,7,8,9,10) ,bty = "n",cex=0.8) 
#dev.off()

#hours
ord <- ordiplot(NMDS_temporal_avg_bray,display = c('sites','species'),choices = c(1,2),type = "n") 
ordihull(ord, zoop_avg$order, display = "sites", draw = c("polygon"),
         col =  hcl.colors(11, palette = "Earth"), alpha = 75,cex = 2)

#sites
ord <- ordiplot(NMDS_temporal_avg_bray,display = c('sites','species'),choices = c(1,2),type = "n") 
ordihull(ord, zoop_avg$site, display = "sites", draw = c("polygon"),
         col =  hcl.colors(2, palette = "Earth"), alpha = 75,cex = 2)
#-----------------------------------------------------------------------------------------#
#time tracking noon vs midnight across all sites and days
#jpeg(file.path(getwd(),"Summer2021-DataAnalysis/Figures/2019-2021_NMDS_1v2_bray_density_noon_midnight.jpg"), width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_noon_midnight_avg_bray,display = c('sites','species'),choices = c(1,2),type = "n") 
ordihull(ord, zoop_avg_noon_midnight$order, display = "sites", draw = c("polygon"),
         col = c("#FF99CC","#440154FF","#FDE725FF"), alpha = 75,cex = 4)

points(NMDS_noon_midnight_avg_bray$points[zoop_avg_noon_midnight$order=="1",1], NMDS_noon_midnight_avg_bray$points[zoop_avg_noon_midnight$order=="1",2], pch=19,col="#FF99CC")
points(NMDS_noon_midnight_avg_bray$points[zoop_avg_noon_midnight$order=="6",1], NMDS_noon_midnight_avg_bray$points[zoop_avg_noon_midnight$order=="6",2], pch=19,col="#440154FF")
points(NMDS_noon_midnight_avg_bray$points[zoop_avg_noon_midnight$order=="11",1], NMDS_noon_midnight_avg_bray$points[zoop_avg_noon_midnight$order=="11",2], pch=19,col="#FDE725FF")

legend("bottomleft", legend=c('noon1','midnight','noon2'), pch=19,
       col = c("#FF99CC","#440154FF","#FDE725FF"),bty = "n",cex=1.2 )
#dev.off() 


#pelagic only tracking density through time
#jpeg(file.path(getwd(),"Summer2021-DataAnalysis/Figures/2019-2021_NMDS_1v2_bray_pelagic_tracking_density_time.jpg"), width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_pel_bray,display = c('sites','species'),choices = c(1,2),type = "n") 
ordihull(ord, zoop_pel$groups, display = "sites", draw = c("polygon"),
         col = c("#008585","#89B199","#EFECBF","#DB9B5A","#C7522B"), alpha = 75,cex = 4)
lines(NMDS_pel_bray$points[zoop_pel$groups=="1",1], NMDS_pel_bray$points[zoop_pel$groups=="1",2], col="#008585")
points(NMDS_pel_bray$points[zoop_pel$groups=="1",1], NMDS_pel_bray$points[zoop_pel$groups=="1",2], pch=c(0,1,2,3,4,5,6,7,8,9,10),col="#008585")

lines(NMDS_pel_bray$points[zoop_pel$groups=="2",1], NMDS_pel_bray$points[zoop_pel$groups=="2",2], col="#89B199")
points(NMDS_pel_bray$points[zoop_pel$groups=="2",1], NMDS_pel_bray$points[zoop_pel$groups=="2",2], pch=c(0,1,2,3,4,5,6,7,8,9,10),col="#89B199")

lines(NMDS_pel_bray$points[zoop_pel$groups=="3",1], NMDS_pel_bray$points[zoop_pel$groups=="3",2], col="#EFECBF")
points(NMDS_pel_bray$points[zoop_pel$groups=="3",1], NMDS_pel_bray$points[zoop_pel$groups=="3",2], pch=c(0,1,2,3,4,5,6,7,8,9,10),col="#EFECBF")

lines(NMDS_pel_bray$points[zoop_pel$groups=="4",1], NMDS_pel_bray$points[zoop_pel$groups=="4",2], col="#DB9B5A")
points(NMDS_pel_bray$points[zoop_pel$groups=="4",1], NMDS_pel_bray$points[zoop_pel$groups=="4",2], pch=c(0,1,2,3,4,5,6,7,8,9,10),col="#DB9B5A")

lines(NMDS_pel_bray$points[zoop_pel$groups=="5",1], NMDS_pel_bray$points[zoop_pel$groups=="5",2], col="#C7522B")
points(NMDS_pel_bray$points[zoop_pel$groups=="5",1], NMDS_pel_bray$points[zoop_pel$groups=="5",2], pch=c(0,1,2,3,4,5,6,7,8,9,10),col="#C7522B")

legend("bottomright", legend=c('10-11 Jul 2019', '24-25 Jul 2019','12-13 Aug 2020','15-16 Jun 2021', '7-8 Jul 2021'), pch=21, 
       pt.bg=c("#008585","#89B199","#EFECBF","#DB9B5A", "#C7522B"),bty = "n",cex=1.2) 
legend("bottomleft", legend=c('12pm','6pm','7pm','8pm','9pm','12am','4am','5am','6am','7am','12pm'), pch=c(0,1,2,3,4,5,6,7,8,9,10) ,bty = "n",cex=1.2) 
#dev.off() 

#littoral only tracking density through time
#jpeg(file.path(getwd(),"Summer2021-DataAnalysis/Figures/2019-2021_NMDS_1v2_bray_littoral_tracking_density_time.jpg"), width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_lit_bray,display = c('sites','species'),choices = c(1,2),type = "n") 
ordihull(ord, zoop_lit$groups, display = "sites", draw = c("polygon"),
         col = c("#008585","#89B199","#EFECBF","#DB9B5A","#C7522B"), alpha = 75,cex = 2)
lines(NMDS_lit_bray$points[zoop_lit$groups=="1",1], NMDS_lit_bray$points[zoop_lit$groups=="1",2], col="#369187")
points(NMDS_lit_bray$points[zoop_lit$groups=="1",1], NMDS_lit_bray$points[zoop_lit$groups=="1",2], pch=c(0,1,2,3,4,5,6,7,8,9,10),col="#369187")

lines(NMDS_lit_bray$points[zoop_lit$groups=="2",1], NMDS_lit_bray$points[zoop_lit$groups=="2",2], col="#C4D5B2")
points(NMDS_lit_bray$points[zoop_lit$groups=="2",1], NMDS_lit_bray$points[zoop_lit$groups=="2",2], pch=c(0,1,2,3,4,5,6,7,8,9,10),col="#C4D5B2")

lines(NMDS_lit_bray$points[zoop_lit$groups=="3",1], NMDS_lit_bray$points[zoop_lit$groups=="3",2], col="#EFDBA7")
points(NMDS_lit_bray$points[zoop_lit$groups=="3",1], NMDS_lit_bray$points[zoop_lit$groups=="3",2], pch=c(0,1,2,3,4,5,6,7,8,9,10),col=c(rep("#EFDBA7",10),rep("#EFECBF",10)))

lines(NMDS_lit_bray$points[zoop_lit$groups=="4",1], NMDS_lit_bray$points[zoop_lit$groups=="4",2], col="#E4BC80")
points(NMDS_lit_bray$points[zoop_lit$groups=="4",1], NMDS_lit_bray$points[zoop_lit$groups=="4",2], pch=c(0,1,2,3,4,5,6,7,8,9,10),col="#E4BC80")

lines(NMDS_lit_bray$points[zoop_lit$groups=="5",1], NMDS_lit_bray$points[zoop_lit$groups=="5",2], col="#CA602E")
points(NMDS_lit_bray$points[zoop_lit$groups=="5",1], NMDS_lit_bray$points[zoop_lit$groups=="5",2], pch=c(0,1,2,3,4,5,6,7,8,9,10),col="#CA602E")

legend("topright", legend=c('10-11 Jul 2019', '24-25 Jul 2019','12-13 Aug 2020','15-16 Jun 2021', '7-8 Jul 2021'), pch=21, 
       pt.bg=c("#369187","#C4D5B2","#EFDBA7","#E4BC80","#CA602E"),bty = "n",cex=0.8) 
legend("bottomleft", legend=c('12pm','6pm','7pm','8pm','9pm','12am','4am','5am','6am','7am','12pm'), pch=c(0,1,2,3,4,5,6,7,8,9,10) ,bty = "n",cex=0.8) 
#dev.off()


#-------------------------------------------------------------------------------#
#                     Calculating euclidean distance                            #
#-------------------------------------------------------------------------------#
#so I'm technically calculating the euclidean distances from bray-curtis distance matrices

#step 1: take NMDS output for each site using NMDS coordinates
zoop_pel_euc <- as.matrix(vegdist(NMDS_pel_bray$points, method='euclidean')) 
zoop_lit_euc <- as.matrix(vegdist(NMDS_lit_bray$points, method='euclidean'))

zoop_sites_euc <- as.matrix(vegdist(NMDS_temporal_bray$points, method='euclidean'))

#step 2: select and sum the 9 distances between connecting points for each of the 5 days
pel_day1 <- sum(zoop_pel_euc[1,2],zoop_pel_euc[2,3],zoop_pel_euc[3,4],zoop_pel_euc[4,5],zoop_pel_euc[5,6],
                zoop_pel_euc[6,7],zoop_pel_euc[7,8],zoop_pel_euc[8,9],zoop_pel_euc[9,10],zoop_pel_euc[10,11])

pel_day2 <- sum(zoop_pel_euc[12,13], zoop_pel_euc[13,14],zoop_pel_euc[14,15],zoop_pel_euc[15,16],zoop_pel_euc[16,17],
                zoop_pel_euc[17,18],zoop_pel_euc[18,19],zoop_pel_euc[19,20],zoop_pel_euc[20,21],zoop_pel_euc[21,22])

pel_day3 <- sum(zoop_pel_euc[23,24], zoop_pel_euc[24,25],zoop_pel_euc[25,26],zoop_pel_euc[26,27],zoop_pel_euc[27,28],
                zoop_pel_euc[28,29],zoop_pel_euc[29,30],zoop_pel_euc[30,31],zoop_pel_euc[31,32],zoop_pel_euc[32,33])

pel_day4 <- sum(zoop_pel_euc[34,35], zoop_pel_euc[35,36],zoop_pel_euc[36,37],zoop_pel_euc[37,38],zoop_pel_euc[38,39],
                zoop_pel_euc[39,40],zoop_pel_euc[40,41],zoop_pel_euc[41,42],zoop_pel_euc[42,43],zoop_pel_euc[43,44])

pel_day5 <- sum(zoop_pel_euc[45,46], zoop_pel_euc[46,47],zoop_pel_euc[47,48],zoop_pel_euc[48,49],zoop_pel_euc[49,50],
                zoop_pel_euc[50,51],zoop_pel_euc[51,52],zoop_pel_euc[52,53],zoop_pel_euc[53,54],zoop_pel_euc[54,55])

mean(pel_day1, pel_day2, pel_day3, pel_day4, pel_day5)

lit_day1 <- sum(zoop_lit_euc[1,2],zoop_lit_euc[2,3],zoop_lit_euc[3,4],zoop_lit_euc[4,5],zoop_lit_euc[5,6],
                zoop_lit_euc[6,7],zoop_lit_euc[7,8],zoop_lit_euc[8,9],zoop_lit_euc[9,10],zoop_lit_euc[10,11])

lit_day2 <- sum(zoop_lit_euc[12,13], zoop_lit_euc[13,14],zoop_lit_euc[14,15],zoop_lit_euc[15,16],zoop_lit_euc[16,17],
                zoop_lit_euc[17,18],zoop_lit_euc[18,19],zoop_lit_euc[19,20],zoop_lit_euc[20,21],zoop_lit_euc[21,22])

lit_day3 <- sum(zoop_lit_euc[23,24], zoop_lit_euc[24,25],zoop_lit_euc[25,26],zoop_lit_euc[26,27],zoop_lit_euc[27,28],
                zoop_lit_euc[28,29],zoop_lit_euc[29,30],zoop_lit_euc[30,31],zoop_lit_euc[31,32],zoop_lit_euc[32,33])

lit_day4 <- sum(zoop_lit_euc[34,35], zoop_lit_euc[35,36],zoop_lit_euc[36,37],zoop_lit_euc[37,38],zoop_lit_euc[38,39],
                zoop_lit_euc[39,40],zoop_lit_euc[40,41],zoop_lit_euc[41,42],zoop_lit_euc[42,43],zoop_lit_euc[43,44])

lit_day5 <- sum(zoop_lit_euc[45,46], zoop_lit_euc[46,47],zoop_lit_euc[47,48],zoop_lit_euc[48,49],zoop_lit_euc[49,50],
                zoop_lit_euc[50,51],zoop_lit_euc[51,52],zoop_lit_euc[52,53],zoop_lit_euc[53,54],zoop_lit_euc[54,55])

mean(lit_day1, lit_day2, lit_day3, lit_day4, lit_day5)

#convert ED matrix back into distance structure for next steps
zoop_pel_euc <- vegdist(NMDS_pel_bray$points, method='euclidean')
zoop_lit_euc <- vegdist(NMDS_lit_bray$points, method='euclidean')

zoop_sites_euc <- vegdist(NMDS_temporal_bray$points, method='euclidean')

#Now calculate the centroids of each polygon and calculate the Euclidean distance between each centroid (n=10)
centroids_pel <- betadisper(zoop_pel_euc, group = as.factor(zoop_pel$groups), type="centroid")
centroids_lit <- betadisper(zoop_lit_euc, group = as.factor(zoop_lit$groups), type="centroid")

centroids_sites <- betadisper(zoop_sites_euc, group = as.factor(zoop_epi_tows$site), type="centroid")

mean(dist(centroids_pel$centroids))
mean(dist(centroids_lit$centroids))

mean(dist(centroids_sites$centroids))
#but this is wrong because it doesn't account for the negative eigenvalues????
#do I need to square the negative ones?



#-------------------------------------------------------------------------------
#step 3: make a dataset of data
euc_distances_df <- data.frame(pelagic=c(pel_day1,pel_day2,pel_day3,pel_day4,pel_day5),
                               littoral=c(lit_day1,lit_day2,lit_day3,lit_day4,lit_day5))
write.csv(euc_distances_df, file.path(getwd(),"/Summer2021-DataAnalysis/SummaryStats/Daily_Euclidean_distances_pelvslit.csv"))

#plot littoral vs pelagic euclidean distances
#jpeg(file.path(getwd(),"Summer2021-DataAnalysis/Figures/2019-2020_pelagic_vs_littoral_euclidean_dist_daily_sums.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euc_distances_df$littoral,euc_distances_df$pelagic, xlab="littoral", ylab="pelagic", 
     main="Daily euclidean distance sums", cex=2.8, pch=21, cex.lab = 1.5)
points(euc_distances_df$littoral[1],euc_distances_df$pelagic[1], bg="#008585",pch=21,cex=3)
points(euc_distances_df$littoral[2],euc_distances_df$pelagic[2], bg="#9BBAA0" ,pch=21,cex=3)
points(euc_distances_df$littoral[3],euc_distances_df$pelagic[3], bg="#F2E2B0",pch=21,cex=3)
points(euc_distances_df$littoral[4],euc_distances_df$pelagic[4], bg="#DEA868",pch=21,cex=3)
points(euc_distances_df$littoral[5],euc_distances_df$pelagic[5], bg="#C7522B",pch=21,cex=3)
legend("bottomleft",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

