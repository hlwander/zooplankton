#multivariate stats for 2019 - 2021 zoop data
#created 25Nov2021

#read in libraries
pacman::p_load(dplyr, vegan, labdsv)

#function to count characters starting at the end of the string
substrEnd <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#read in zoop data from all 3 years
zoops2019<- read.csv('/Users/heatherwander/Documents/VirginiaTech/research/zooplankton/Summer2019-DataAnalysis/SummaryStats/FCR_ZooplanktonSummary2019.csv',header = TRUE)
zoops2020<- read.csv('/Users/heatherwander/Documents/VirginiaTech/research/zooplankton/Summer2020-DataAnalysis/SummaryStats/FCR_ZooplanktonSummary2020.csv',header = TRUE)
zoops2021<- read.csv('SummaryStats/FCR_ZooplanktonSummary2021.csv',header = TRUE)

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
                                  "Daphnia_density_NopL", "Ceriodaphnia_density_NopL","nauplius_density_NopL",
                                  "Synchaetidae_density_NopL", "Conochilidae_density_NopL")

#add collotheca column for 2021 (none found so column was not created
zoops2021$Collothecidae_density_NopL <- 0

#combine all zoop datasets
zoops <- rbind(zoops2019,zoops2020,zoops2021)

#ignore 20 um and horizontal trap samples
zoops <- zoops %>% filter(mesh_size_μm >20, na.rm=TRUE)

#create df for temporal epi tows
zoop_epi_tows <- zoops[zoops$site_no=="BVR_l" | zoops$site_no=="BVR_50_p", ] %>%
  mutate(Hour=substr(Hour,1,2)) %>% 
  select(!c(DepthOfTow_m, sample_ID)) %>% group_by(site_no,collect_date,Hour) %>%
  summarise(across(everything(),list(mean)))

zoop_epi_tows$time <-ifelse(zoop_epi_tows$Hour=="12" | zoop_epi_tows$Hour=="11", "noon", ifelse(
   zoop_epi_tows$Hour =="0:" | zoop_epi_tows$Hour =="23", "midnight",ifelse(zoop_epi_tows$Hour=="18"|
   zoop_epi_tows$Hour=="19" | zoop_epi_tows$Hour=="20" | zoop_epi_tows$Hour=="21", "sunset", "sunrise")))

zoop_epi_tows$site <- ifelse(substrEnd(zoop_epi_tows$site_no,1)=="p","pel","lit")

#also add columns to group times by sampling event and time and siites tigether
zoop_epi_tows$groups <- ifelse(zoop_epi_tows$collect_date=="2019-07-10" | zoop_epi_tows$collect_date=="2019-07-11","1",
                        ifelse(zoop_epi_tows$collect_date=="2019-07-24" | zoop_epi_tows$collect_date=="2019-07-25","2","3"))

zoop_epi_tows$timesite <- paste0(zoop_epi_tows$time,zoop_epi_tows$site)

zoop_epi_tows$timegroup <- paste0(zoop_epi_tows$time,zoop_epi_tows$groups)

#------------------------------------------------------------------------------#
# set up data for NMDS
#relies on rank orders for ordination, no assumptions of linear relationship
zoop_temporal_dens <- zoop_epi_tows[,c(6:16)]

#transforming data - hellinger transformation because gives low weight to low/zero values
#converts species abundances from absolute to relative - use w/ bray curtis
zoop_temporal_dens_trans <- hellinger(zoop_temporal_dens)

#-------------------------------------------------------------------------------#
#                  Bray-curtis dissimilarity --> NMDS figs                      #
#-------------------------------------------------------------------------------#
par(ask=TRUE)
NMDS_temporal_bray <- metaMDS(zoop_temporal_dens_trans, distance='bray', k=3, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_temporal_bray$stress

#  Points divided into time groups
#jpeg("Figures/2019-2020_NMDS_1v2_bray_days.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_temporal_bray,display = c('sites','species'),choices = c(1,2),type = "n") 
points(NMDS_temporal_bray$points[zoop_epi_tows$groups=="1",1], NMDS_temporal_bray$points[zoop_epi_tows$groups=="1",2], pch=21,bg="tomato3")
points(NMDS_temporal_bray$points[zoop_epi_tows$groups=="2",1], NMDS_temporal_bray$points[zoop_epi_tows$groups=="2",2], pch=21,bg="tan1")
points(NMDS_temporal_bray$points[zoop_epi_tows$groups=="3",1], NMDS_temporal_bray$points[zoop_epi_tows$groups=="3",2], pch=21,bg="steelblue4")
legend("bottomright", legend=c('Day1','Day2','Day3'), pch=21, pt.bg=c("tomato3","tan1","steelblue4"),bty = "n") 
ordihull(ord, zoop_epi_tows$groups, display = "sites", draw = c("polygon"),
         col = c("tomato3","tan1","steelblue4"), alpha = 75,cex = 2)
text(NMDS_temporal_bray$species[,1],NMDS_temporal_bray$species[,2], labels = c("Calanoida","Cyclopoida","Keratella","Kellicottia",
     "Bosmina","Daphnia","Ceriodaphnia","nauplius","Collothecidae","Synchaetidae","Conochilidae"), cex=0.9)
#dev.off()

#jpeg("Figures/2019-2020_NMDS_1v3_bray_days.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_temporal_bray,display = c('sites','species'),choices = c(1,3),type = "n") 
points(NMDS_temporal_bray$points[zoop_epi_tows$groups=="1",1], NMDS_temporal_bray$points[zoop_epi_tows$groups=="1",3], pch=21,bg="tomato3")
points(NMDS_temporal_bray$points[zoop_epi_tows$groups=="2",1], NMDS_temporal_bray$points[zoop_epi_tows$groups=="2",3], pch=21,bg="tan1")
points(NMDS_temporal_bray$points[zoop_epi_tows$groups=="3",1], NMDS_temporal_bray$points[zoop_epi_tows$groups=="3",3], pch=21,bg="steelblue4")
legend("bottomright", legend=c('Day1','Day2','Day3'), pch=21, pt.bg=c("tomato3","tan1","steelblue4"),bty = "n") 
ordihull(ord, zoop_epi_tows$groups, display = "sites", draw = c("polygon"),
         col = c("tomato3","tan1","steelblue4"), alpha = 75,cex = 2)
text(NMDS_temporal_bray$species[,1],NMDS_temporal_bray$species[,2], labels = c("Calanoida","Cyclopoida","Keratella","Kellicottia",
     "Bosmina","Daphnia","Ceriodaphnia","nauplius","Collothecidae","Synchaetidae","Conochilidae"), cex=0.9)
#dev.off()

#jpeg("Figures/2019-2020_NMDS_1v2_bray_sites.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_temporal_bray,display = c('sites','species'),choices = c(1,2),type = "n") 
points(NMDS_temporal_bray$points[zoop_epi_tows$site=="pel",1], NMDS_temporal_bray$points[zoop_epi_tows$site=="pel",2], pch=21,bg="slateblue4")
points(NMDS_temporal_bray$points[zoop_epi_tows$site=="lit",1], NMDS_temporal_bray$points[zoop_epi_tows$site=="lit",2], pch=21,bg="palegreen3")
legend("bottomright", legend=c("Pelagic","Littoral"), pch=21, pt.bg=c("slateblue4","palegreen3"),bty = "n") 
ordihull(ord, zoop_epi_tows$site, display = "sites", draw = c("polygon"),
         col = c("palegreen3","slateblue4"), alpha = 75,cex = 2)
text(NMDS_temporal_bray$species[,1],NMDS_temporal_bray$species[,2], labels = c("Calanoida","Cyclopoida","Keratella","Kellicottia",
     "Bosmina","Daphnia","Ceriodaphnia","nauplius","Collothecidae","Synchaetidae","Conochilidae"), cex=0.9)
#dev.off()

#jpeg("Figures/2019-2020_NMDS_1v3_bray_sites.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_temporal_bray,display = c('sites','species'),choices = c(1,3),type = "n") 
points(NMDS_temporal_bray$points[zoop_epi_tows$site=="pel",1], NMDS_temporal_bray$points[zoop_epi_tows$site=="pel",3], pch=21,bg="slateblue4")
points(NMDS_temporal_bray$points[zoop_epi_tows$site=="lit",1], NMDS_temporal_bray$points[zoop_epi_tows$site=="lit",3], pch=21,bg="palegreen3")
legend("bottomright", legend=c("Pelagic","Littoral"), pch=21, pt.bg=c("slateblue4","palegreen3"),bty = "n") 
ordihull(ord, zoop_epi_tows$site, display = "sites", draw = c("polygon"),
         col = c("palegreen3","slateblue4"), alpha = 75,cex = 2)
text(NMDS_temporal_bray$species[,1],NMDS_temporal_bray$species[,2], labels = c("Calanoida","Cyclopoida","Keratella","Kellicottia",
     "Bosmina","Daphnia","Ceriodaphnia","nauplius","Collothecidae","Synchaetidae","Conochilidae"), cex=0.9)
#dev.off()


#jpeg("Figures/2019-2020_NMDS_1v2_bray_time.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_temporal_bray,display = c('sites','species'),choices = c(1,2),type = "n") 
points(NMDS_temporal_bray$points[zoop_epi_tows$time=='sunrise',1], NMDS_temporal_bray$points[zoop_epi_tows$time=='sunrise',2], pch=21,bg="#CC79A7")
points(NMDS_temporal_bray$points[zoop_epi_tows$time=='noon',1], NMDS_temporal_bray$points[zoop_epi_tows$time=='noon',2], pch=21,bg="#F0E442")
points(NMDS_temporal_bray$points[zoop_epi_tows$time=='sunset',1], NMDS_temporal_bray$points[zoop_epi_tows$time=='sunset',2], pch=21,bg="#009E73")
points(NMDS_temporal_bray$points[zoop_epi_tows$time=='midnight',1], NMDS_temporal_bray$points[zoop_epi_tows$time=='midnight',2], pch=21,bg="#0072B2")
legend("bottomright", legend=c('sunrise','noon', 'sunset', 'midnight'), pch=21, pt.bg=c('#CC79A7', '#F0E442','#009E73','#0072B2'),bty = "n") 
ordihull(ord, zoop_epi_tows$time, display = "sites", draw = c("polygon"),
         col = c("#0072B2", "#F0E442","#CC79A7","#009E73"), alpha = 75,cex = 2)
text(NMDS_temporal_bray$species[,1],NMDS_temporal_bray$species[,2], labels = c("Calanoida","Cyclopoida","Keratella","Kellicottia",
     "Bosmina","Daphnia","Ceriodaphnia","nauplius","Collothecidae","Synchaetidae","Conochilidae"), cex=0.9)
#dev.off()

#jpeg("Figures/2019-2020_NMDS_1v3_bray_time.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_temporal_bray,display = c('sites','species'),choices = c(1,3),type = "n") 
points(NMDS_temporal_bray$points[zoop_epi_tows$time=='sunrise',1], NMDS_temporal_bray$points[zoop_epi_tows$time=='sunrise',3], pch=21,bg="#CC79A7")
points(NMDS_temporal_bray$points[zoop_epi_tows$time=='noon',1], NMDS_temporal_bray$points[zoop_epi_tows$time=='noon',3], pch=21,bg="#F0E442")
points(NMDS_temporal_bray$points[zoop_epi_tows$time=='sunset',1], NMDS_temporal_bray$points[zoop_epi_tows$time=='sunset',3], pch=21,bg="#009E73")
points(NMDS_temporal_bray$points[zoop_epi_tows$time=='midnight',1], NMDS_temporal_bray$points[zoop_epi_tows$time=='midnight',3], pch=21,bg="#0072B2")
legend("bottomright", legend=c('sunrise','noon', 'sunset', 'midnight'), pch=21, pt.bg=c('#CC79A7', '#F0E442','#009E73','#0072B2'),bty = "n") 
ordihull(ord, zoop_epi_tows$time, display = "sites", draw = c("polygon"),
         col = c("#0072B2", "#F0E442","#CC79A7","#009E73"), alpha = 75,cex = 2)
text(NMDS_temporal_bray$species[,1],NMDS_temporal_bray$species[,2], labels = c("Calanoida","Cyclopoida","Keratella","Kellicottia",
     "Bosmina","Daphnia","Ceriodaphnia","nauplius","Collothecidae","Synchaetidae","Conochilidae"), cex=0.9)
#dev.off()


#-------------------------------------------------------------------------------#
#              want to make some kind of tracking time NMDS fig                 #
#-------------------------------------------------------------------------------#

#jpeg("Figures/2019-2020_NMDS_1v2_bray_timegroups.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_temporal_bray,display = c('sites','species'),choices = c(1,2),type = "n") 
points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunrise1',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunrise1',2], pch=21,bg="cornsilk")
points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunrise2',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunrise2',2], pch=21,bg="yellow")
points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunrise3',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunrise3',2], pch=21,bg="yellow3")
points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='noon1',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='noon1',2], pch=21,bg="lightpink")
points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='noon2',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='noon2',2], pch=21,bg="indianred")
points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='noon3',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='noon3',2], pch=21,bg="red")
points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunset1',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunset1',2], pch=21,bg="palegreen")
points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunset2',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunset2',2], pch=21,bg="mediumseagreen")
points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunset3',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunset3',2], pch=21,bg="darkgreen")
points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='midnight1',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='midnight1',2], pch=21,bg="lightsteelblue")
points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='midnight2',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='midnight2',2], pch=21,bg="skyblue")
points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='midnight3',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='midnight3',2], pch=21,bg="steelblue")
ordiellipse(ord, zoop_epi_tows$timegroup, display = "sites", kind="sd",draw="lines", col = c("lightsteelblue", "skyblue","steelblue","lightpink","indianred","red","cornsilk","yellow","yellow3","palegreen","mediumseagreen","darkgreen"), alpha = 75,cex = 2)
legend("bottomright", legend=c('sunrise1', 'sunrise2', 'sunrise3', 'noon1', 'noon2', 'noon3', 'sunset1', 'sunset2', 'sunset3', 'midnight1', 'midnight2', 'midnight3'), pch=21,
       pt.bg=c('cornsilk', 'yellow','yellow3','lightpink','indianred','red','palegreen','mediumseagreen','darkgreen','lightsteelblue','skyblue','steelblue'),bty = "n") 
#dev.off()

#jpeg("Figures/2019-2020_NMDS_1v3_bray_timegroups.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_temporal_bray,display = c('sites','species'),choices = c(1,3),type = "n") 
points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunrise1',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunrise1',3], pch=21,bg="cornsilk")
points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunrise2',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunrise2',3], pch=21,bg="yellow")
points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunrise3',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunrise3',3], pch=21,bg="yellow3")
points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='noon1',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='noon1',3], pch=21,bg="lightpink")
points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='noon2',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='noon2',3], pch=21,bg="indianred")
points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='noon3',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='noon3',3], pch=21,bg="red")
points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunset1',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunset1',3], pch=21,bg="palegreen")
points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunset2',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunset2',3], pch=21,bg="mediumseagreen")
points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunset3',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='sunset3',3], pch=21,bg="darkgreen")
points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='midnight1',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='midnight1',3], pch=21,bg="lightsteelblue")
points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='midnight2',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='midnight2',3], pch=21,bg="skyblue")
points(NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='midnight3',1], NMDS_temporal_bray$points[zoop_epi_tows$timegroup=='midnight3',3], pch=21,bg="steelblue")
ordiellipse(ord, zoop_epi_tows$timegroup, display = "sites", kind="sd",draw="lines", col = c("lightsteelblue", "skyblue","steelblue","lightpink","indianred","red","cornsilk","yellow","yellow3","palegreen","mediumseagreen","darkgreen"), alpha = 75,cex = 2)
legend("bottomright", legend=c('sunrise1', 'sunrise2', 'sunrise3', 'noon1', 'noon2', 'noon3', 'sunset1', 'sunset2', 'sunset3', 'midnight1', 'midnight2', 'midnight3'), pch=21,
       pt.bg=c('cornsilk', 'yellow','yellow3','lightpink','indianred','red','palegreen','mediumseagreen','darkgreen','lightsteelblue','skyblue','steelblue'),bty = "n") 
#dev.off()
