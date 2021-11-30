#multivariate stats for 2019 - 2021 zoop data
#created 25Nov2021

#read in libraries
pacman::p_load(dplyr, vegan, labdsv, goeveg)

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

#add collotheca column for 2021 (none found so column was not created)
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


#-------------------------------------------------------------------------------#
#           Averaging zoops by time and campaign/day for new NMDS               #
#-------------------------------------------------------------------------------#

#first average times for each 24-hour campaign so there are 10 points per day (basically just averaging noon and midnight)
zoop_epi_tows$order <- ifelse(zoop_epi_tows$Hour=="4:",1,ifelse(zoop_epi_tows$Hour=="5:",2,ifelse(
  zoop_epi_tows$Hour=="6:",3,ifelse(zoop_epi_tows$Hour=="7:",4,ifelse(zoop_epi_tows$time=="noon",5,ifelse(
    zoop_epi_tows$Hour=="18",6,ifelse(zoop_epi_tows$Hour=="19",7,ifelse(zoop_epi_tows$Hour=="20",8,ifelse(zoop_epi_tows$Hour=="21",9,10)))))))))

#average by groups, site, then order
zoop_avg <- zoop_epi_tows %>% group_by(groups,site,order) %>%
  summarise_at(vars(Calanoida_density_NopL_1:Conochilidae_density_NopL_1), list(mean = mean))

#pelagic df
zoop_pel <- zoop_avg[zoop_avg$site=="pel",]
zoop_lit <- zoop_avg[zoop_avg$site=="lit",]

#only select data cols
zoop_temporal_avg_dens <- zoop_avg[,c(4:14)]
zoop_pel_dens <- zoop_pel[,c(4:14)]
zoop_lit_dens <- zoop_lit[,c(4:14)]

#transforming data - hellinger transformation because gives low weight to low/zero values
#converts species abundances from absolute to relative - use w/ bray curtis
zoop_temporal_dens_avg_trans <- hellinger(zoop_temporal_avg_dens)
zoop_pel_dens_trans <- hellinger(zoop_pel_dens)
zoop_lit_dens_trans <- hellinger(zoop_lit_dens)

#now do NMDS using averages
NMDS_temporal_avg_bray <- metaMDS(zoop_temporal_dens_avg_trans, distance='bray', k=3, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_temporal_avg_bray$stress

NMDS_pel_bray <- metaMDS(zoop_pel_dens_trans, distance='bray', k=3, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_pel_bray$stress

NMDS_lit_bray <- metaMDS(zoop_lit_dens_trans, distance='bray', k=3, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_lit_bray$stress

#-------------------------------------------------------------------------------#
#                                 Tracking figs                                 #
#-------------------------------------------------------------------------------#

#jpeg("Figures/2019-2020_NMDS_1v2_bray_temporal_avg_days.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_temporal_avg_bray,display = c('sites','species'),choices = c(1,2),type = "n") 
ordihull(ord, zoop_avg$groups, display = "sites", draw = c("polygon"),
         col = c("tomato3","tan1","steelblue4"), alpha = 75,cex = 2)
lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="1" & zoop_avg$site=="pel",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="1" & zoop_avg$site=="pel",2], col="tomato3")
lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="1" & zoop_avg$site=="lit",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="1" & zoop_avg$site=="lit",2], col="tomato")
points(NMDS_temporal_avg_bray$points[zoop_avg$groups=="1",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="1",2], pch=rep(c(0,1,2,3,4,5,6,7,8,9),2),col=c(rep("tomato",10),rep("tomato3",10)))

lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="2" & zoop_avg$site=="pel",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="2" & zoop_avg$site=="pel",2], col="tan4")
lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="2" & zoop_avg$site=="lit",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="2" & zoop_avg$site=="lit",2], col="tan1")
points(NMDS_temporal_avg_bray$points[zoop_avg$groups=="2",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="2",2], pch=rep(c(0,1,2,3,4,5,6,7,8,9),2),col=c(rep("tan1",10),rep("tan4",10)))

lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="3" & zoop_avg$site=="pel",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="3" & zoop_avg$site=="pel",2], col="steelblue4")
lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="3" & zoop_avg$site=="lit",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="3" & zoop_avg$site=="lit",2], col="steelblue1")
points(NMDS_temporal_avg_bray$points[zoop_avg$groups=="3",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="3",2], pch=rep(c(0,1,2,3,4,5,6,7,8,9),2),col=c(rep("steelblue1",10),rep("steelblue4",10)))

legend("bottomright", legend=c('pelagic day1','littoral day1','pelagic day2','littoral day2','pelagic day3','littoral day3'), pch=21, pt.bg=c("tomato3","tomato","tan4","tan1", "steelblue4","steelblue1"),bty = "n",cex=0.8) 
legend("bottomleft", legend=c('sunrise1','sunrise2','sunrise3','sunrise4','noon','sunset1','sunset2','sunset3','sunset4','midnight'), pch=c(0,1,2,3,4,5,6,7,8,9) ,bty = "n",cex=0.8) 
#dev.off()

#jpeg("Figures/2019-2020_NMDS_1v3_bray_temporal_avg_days.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_temporal_avg_bray,display = c('sites','species'),choices = c(1,3),type = "n") 
ordihull(ord, zoop_avg$groups, display = "sites", draw = c("polygon"),
         col = c("tomato3","tan1","steelblue4"), alpha = 75,cex = 2)
lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="1" & zoop_avg$site=="pel",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="1" & zoop_avg$site=="pel",3], col="tomato3")
lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="1" & zoop_avg$site=="lit",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="1" & zoop_avg$site=="lit",3], col="tomato")
points(NMDS_temporal_avg_bray$points[zoop_avg$groups=="1",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="1",3], pch=rep(c(0,1,2,3,4,5,6,7,8,9),2),col=c(rep("tomato",10),rep("tomato3",10)))

lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="2" & zoop_avg$site=="pel",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="2" & zoop_avg$site=="pel",3], col="tan4")
lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="2" & zoop_avg$site=="lit",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="2" & zoop_avg$site=="lit",3], col="tan1")
points(NMDS_temporal_avg_bray$points[zoop_avg$groups=="2",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="2",3], pch=rep(c(0,1,2,3,4,5,6,7,8,9),2),col=c(rep("tan1",10),rep("tan4",10)))

lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="3" & zoop_avg$site=="pel",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="3" & zoop_avg$site=="pel",3], col="steelblue4")
lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="3" & zoop_avg$site=="lit",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="3" & zoop_avg$site=="lit",3], col="steelblue1")
points(NMDS_temporal_avg_bray$points[zoop_avg$groups=="3",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="3",3], pch=rep(c(0,1,2,3,4,5,6,7,8,9),2),col=c(rep("steelblue1",10),rep("steelblue4",10)))

legend("bottomright", legend=c('pelagic day1','littoral day1','pelagic day2','littoral day2','pelagic day3','littoral day3'), pch=21, pt.bg=c("tomato3","tomato","tan4","tan1", "steelblue4","steelblue1"),bty = "n",cex=0.8) 
legend("bottomleft", legend=c('sunrise1','sunrise2','sunrise3','sunrise4','noon','sunset1','sunset2','sunset3','sunset4','midnight'), pch=c(0,1,2,3,4,5,6,7,8,9) ,bty = "n",cex=0.8) 
#dev.off()

#pelagic only tracking density through time
#jpeg("Figures/2019-2020_NMDS_1v2_bray_pelagic_tracking_density_time.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_pel_bray,display = c('sites','species'),choices = c(1,2),type = "n") 
ordihull(ord, zoop_pel$groups, display = "sites", draw = c("polygon"),
         col = c("tomato3","tan1","steelblue4"), alpha = 75,cex = 2)
lines(NMDS_pel_bray$points[zoop_pel$groups=="1",1], NMDS_pel_bray$points[zoop_pel$groups=="1",2], col="tomato3")
points(NMDS_pel_bray$points[zoop_pel$groups=="1",1], NMDS_pel_bray$points[zoop_pel$groups=="1",2], pch=c(0,1,2,3,4,5,6,7,8,9),col="tomato3")

lines(NMDS_pel_bray$points[zoop_pel$groups=="2",1], NMDS_pel_bray$points[zoop_pel$groups=="2",2], col="tan4")
points(NMDS_pel_bray$points[zoop_pel$groups=="2",1], NMDS_pel_bray$points[zoop_pel$groups=="2",2], pch=c(0,1,2,3,4,5,6,7,8,9),col="tan4")

lines(NMDS_pel_bray$points[zoop_pel$groups=="3",1], NMDS_pel_bray$points[zoop_pel$groups=="3",2], col="steelblue4")
points(NMDS_pel_bray$points[zoop_pel$groups=="3",1], NMDS_pel_bray$points[zoop_pel$groups=="3",2], pch=c(0,1,2,3,4,5,6,7,8,9),col="steelblue4")

legend("bottomright", legend=c('day1','day2','day3'), pch=21, pt.bg=c("tomato3","tan4","steelblue4"),bty = "n",cex=0.8) 
legend("bottomleft", legend=c('sunrise1','sunrise2','sunrise3','sunrise4','noon','sunset1','sunset2','sunset3','sunset4','midnight'), pch=c(0,1,2,3,4,5,6,7,8,9) ,bty = "n",cex=0.8) 
#dev.off()

#jpeg("Figures/2019-2020_NMDS_1v3_bray_pelagic_tracking_density_time.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_pel_bray,display = c('sites','species'),choices = c(1,3),type = "n") 
ordihull(ord, zoop_pel$groups, display = "sites", draw = c("polygon"),
         col = c("tomato3","tan1","steelblue4"), alpha = 75,cex = 2)
lines(NMDS_pel_bray$points[zoop_pel$groups=="1",1], NMDS_pel_bray$points[zoop_pel$groups=="1",3], col="tomato3")
points(NMDS_pel_bray$points[zoop_pel$groups=="1",1], NMDS_pel_bray$points[zoop_pel$groups=="1",3], pch=c(0,1,2,3,4,5,6,7,8,9),col="tomato3")

lines(NMDS_pel_bray$points[zoop_pel$groups=="2",1], NMDS_pel_bray$points[zoop_pel$groups=="2",3], col="tan4")
points(NMDS_pel_bray$points[zoop_pel$groups=="2",1], NMDS_pel_bray$points[zoop_pel$groups=="2",3], pch=c(0,1,2,3,4,5,6,7,8,9),col="tan4")

lines(NMDS_pel_bray$points[zoop_pel$groups=="3",1], NMDS_pel_bray$points[zoop_pel$groups=="3",3], col="steelblue4")
points(NMDS_pel_bray$points[zoop_pel$groups=="3",1], NMDS_pel_bray$points[zoop_pel$groups=="3",3], pch=c(0,1,2,3,4,5,6,7,8,9),col="steelblue4")

legend("bottomright", legend=c('day1','day2','day3'), pch=21, pt.bg=c("tomato3","tan4","steelblue4"),bty = "n",cex=0.8) 
legend("bottomleft", legend=c('sunrise1','sunrise2','sunrise3','sunrise4','noon','sunset1','sunset2','sunset3','sunset4','midnight'), pch=c(0,1,2,3,4,5,6,7,8,9) ,bty = "n",cex=0.8) 
#dev.off()

#jpeg("Figures/2019-2020_NMDS_1v2_bray_littoral_tracking_density_time.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_lit_bray,display = c('sites','species'),choices = c(1,2),type = "n") 
ordihull(ord, zoop_lit$groups, display = "sites", draw = c("polygon"),
         col = c("tomato3","tan1","steelblue4"), alpha = 75,cex = 2)
lines(NMDS_lit_bray$points[zoop_lit$groups=="1",1], NMDS_lit_bray$points[zoop_lit$groups=="1",2], col="tomato")
points(NMDS_lit_bray$points[zoop_lit$groups=="1",1], NMDS_lit_bray$points[zoop_lit$groups=="1",2], pch=c(0,1,2,3,4,5,6,7,8,9),col="tomato")

lines(NMDS_lit_bray$points[zoop_lit$groups=="2",1], NMDS_lit_bray$points[zoop_lit$groups=="2",2], col="tan1")
points(NMDS_lit_bray$points[zoop_lit$groups=="2",1], NMDS_lit_bray$points[zoop_lit$groups=="2",2], pch=c(0,1,2,3,4,5,6,7,8,9),col="tan1")

lines(NMDS_lit_bray$points[zoop_lit$groups=="3",1], NMDS_lit_bray$points[zoop_lit$groups=="3",2], col="steelblue1")
points(NMDS_lit_bray$points[zoop_lit$groups=="3",1], NMDS_lit_bray$points[zoop_lit$groups=="3",2], pch=c(0,1,2,3,4,5,6,7,8,9),col="steelblue1")

legend("bottomright", legend=c('day1','day2','day3'), pch=21, pt.bg=c("tomato","tan1","steelblue1"),bty = "n",cex=0.8) 
legend("bottomleft", legend=c('sunrise1','sunrise2','sunrise3','sunrise4','noon','sunset1','sunset2','sunset3','sunset4','midnight'), pch=c(0,1,2,3,4,5,6,7,8,9) ,bty = "n",cex=0.8) 
#dev.off()

#jpeg("Figures/2019-2020_NMDS_1v3_bray_littoral_tracking_density_time.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_lit_bray,display = c('sites','species'),choices = c(1,3),type = "n") 
ordihull(ord, zoop_lit$groups, display = "sites", draw = c("polygon"),
         col = c("tomato3","tan1","steelblue4"), alpha = 75,cex = 2)
lines(NMDS_lit_bray$points[zoop_lit$groups=="1",1], NMDS_lit_bray$points[zoop_lit$groups=="1",3], col="tomato")
points(NMDS_lit_bray$points[zoop_lit$groups=="1",1], NMDS_lit_bray$points[zoop_lit$groups=="1",3], pch=c(0,1,2,3,4,5,6,7,8,9),col="tomato")

lines(NMDS_lit_bray$points[zoop_lit$groups=="2",1], NMDS_lit_bray$points[zoop_lit$groups=="2",3], col="tan1")
points(NMDS_lit_bray$points[zoop_lit$groups=="2",1], NMDS_lit_bray$points[zoop_lit$groups=="2",3], pch=c(0,1,2,3,4,5,6,7,8,9),col="tan1")

lines(NMDS_lit_bray$points[zoop_lit$groups=="3",1], NMDS_lit_bray$points[zoop_lit$groups=="3",3], col="steelblue1")
points(NMDS_lit_bray$points[zoop_lit$groups=="3",1], NMDS_lit_bray$points[zoop_lit$groups=="3",3], pch=c(0,1,2,3,4,5,6,7,8,9),col="steelblue1")

legend("bottomright", legend=c('day1','day2','day3'), pch=21, pt.bg=c("tomato","tan1","steelblue1"),bty = "n",cex=0.8) 
legend("bottomleft", legend=c('sunrise1','sunrise2','sunrise3','sunrise4','noon','sunset1','sunset2','sunset3','sunset4','midnight'), pch=c(0,1,2,3,4,5,6,7,8,9) ,bty = "n",cex=0.8) 
#dev.off()


#-------------------------------------------------------------------------------#
#                     Calculating euclidean distance                            #
#-------------------------------------------------------------------------------#

#step 1: take NMDS output for each site using NMDS coordinates
zoop_pel_euc <- as.matrix(vegdist(NMDS_pel_bray$points, method='euclidean'))
zoop_lit_euc <- as.matrix(vegdist(NMDS_lit_bray$points, method='euclidean'))

#step 2: select and sum the 10 distances between connecting points for each of the 3 days
pel_day1 <- sum(zoop_pel_euc[1,2],zoop_pel_euc[2,3],zoop_pel_euc[3,4],zoop_pel_euc[4,5],zoop_pel_euc[5,6],
                zoop_pel_euc[6,7],zoop_pel_euc[7,8],zoop_pel_euc[8,9],zoop_pel_euc[9,10])

pel_day2 <- sum(zoop_pel_euc[11,12], zoop_pel_euc[12,13],zoop_pel_euc[13,14],zoop_pel_euc[14,15],zoop_pel_euc[15,16],
                zoop_pel_euc[16,17],zoop_pel_euc[17,18],zoop_pel_euc[18,19],zoop_pel_euc[19,20])

pel_day3 <- sum(zoop_pel_euc[21,22], zoop_pel_euc[22,23],zoop_pel_euc[23,24],zoop_pel_euc[24,25],zoop_pel_euc[25,26],
                zoop_pel_euc[26,27],zoop_pel_euc[27,28],zoop_pel_euc[28,29],zoop_pel_euc[29,30])

lit_day1 <- sum(zoop_lit_euc[1,2],zoop_lit_euc[2,3],zoop_lit_euc[3,4],zoop_lit_euc[4,5],zoop_lit_euc[5,6],
                zoop_lit_euc[6,7],zoop_lit_euc[7,8],zoop_lit_euc[8,9],zoop_lit_euc[9,10])

lit_day2 <- sum(zoop_lit_euc[11,12],zoop_lit_euc[12,13],zoop_lit_euc[13,14],zoop_lit_euc[14,15],zoop_lit_euc[15,16],
                zoop_lit_euc[16,17],zoop_lit_euc[17,18],zoop_lit_euc[18,19],zoop_lit_euc[19,20])

lit_day3 <- sum(zoop_lit_euc[21,22],zoop_lit_euc[22,23],zoop_lit_euc[23,24],zoop_lit_euc[24,25],zoop_lit_euc[25,26],
                zoop_lit_euc[26,27],zoop_lit_euc[27,28],zoop_lit_euc[28,29],zoop_lit_euc[29,30])

#step 3: make a dataset of data
euc_distances_df <- data.frame(pelagic=c(pel_day1,pel_day2,pel_day3),littoral=c(lit_day1,lit_day2,lit_day3))

#plot littoral vs pelagic euclidean distances
#jpeg("Figures/2019-2020_pelagic_vs_littoral_euclidean_dist_daily_sums.jpg", width = 6, height = 5, units = "in",res = 300)
plot(euc_distances_df$littoral,euc_distances_df$pelagic, xlab="littoral", ylab="pelagic", 
     main="Daily euclidean distance sums", cex=2.8, pch=21, cex.lab = 1.5)
points(euc_distances_df$littoral[1],euc_distances_df$pelagic[1], bg="orangered3",pch=21,cex=3)
points(euc_distances_df$littoral[2],euc_distances_df$pelagic[2], bg="cadetblue",pch=21,cex=3)
points(euc_distances_df$littoral[3],euc_distances_df$pelagic[3], bg="mediumpurple4",pch=21,cex=3)
legend("bottomleft",legend=c("day1","day2","day3"),pch=21, pt.bg=c("orangered3","cadetblue","mediumpurple4"),bty = "n",cex=1.4)
#dev.off()

#driver data df
euclidean_drivers_df <- data.frame(site = c("pel_day1","pel_day2","pel_day3","lit_day1","lit_day2","lit_day3"),
                                   euc_dist = c(euc_distances_df$pelagic,euc_distances_df$littoral))

#read in DO and temp from sampling days
a <- read.csv("RawData/071019_bvr50_a.csv")
b <- read.csv("RawData/072419_bvr50_a.csv")
c <- read.csv("RawData/081220_bvr50.csv")

#pull 0.1 and 9m 
a <- a[a$Depth_m >=0.1 & a$Depth_m<=0.13 | a$Depth_m >=9 & a$Depth_m<=9.03,]
b <- b[b$Depth_m >=0.1 & b$Depth_m<0.112 | b$Depth_m >=9 & b$Depth_m<=9.03,][c(1,3),]
c <- c[c$Depth_m ==0.1 | c$Depth_m >=9 & c$Depth_m<=9.01,]

ctd <- rbind(a,b,c)

#add ctd cols to df
euclidean_drivers_df$DO_0.1m <- c(ctd$DO_mgL[c(1,3,5)],NA,NA,NA)
euclidean_drivers_df$DO_9m <- c(ctd$DO_mgL[c(2,4,6)],NA,NA,NA)
euclidean_drivers_df$temp_0.1m <- c(ctd$Temp_C[c(1,3,5)],NA,NA,NA)
euclidean_drivers_df$temp_9m <- c(ctd$Temp_C[c(2,4,6)],NA,NA,NA)
euclidean_drivers_df$chl_0.1m <- c(ctd$Chla_ugL[c(1,3,5)],NA,NA,NA)
euclidean_drivers_df$chl_9m <- c(ctd$Chla_ugL[c(2,4,6)],NA,NA,NA)
euclidean_drivers_df$sp_cond_0.1m <- c(ctd$Spec_Cond_uScm[c(1,3,5)],NA,NA,NA)
euclidean_drivers_df$sp_cond_9m <- c(ctd$Spec_Cond_uScm[c(2,4,6)],NA,NA,NA)
euclidean_drivers_df$percent_DO_0.1m <- c(ctd$DO_pSat[c(1,3,5)],NA,NA,NA)
euclidean_drivers_df$percent_DO_9m <- c(ctd$DO_pSat[c(2,4,6)],NA,NA,NA)
  
#calculate avg density col for each day at both sites
euclidean_drivers_df$avg_dens <- c(mean(zoop_epi_tows$ZoopDensity_No.pL_1[zoop_epi_tows$site=="pel" & zoop_epi_tows$groups==1]),
                                   mean(zoop_epi_tows$ZoopDensity_No.pL_1[zoop_epi_tows$site=="pel" & zoop_epi_tows$groups==2]),
                                   mean(zoop_epi_tows$ZoopDensity_No.pL_1[zoop_epi_tows$site=="pel" & zoop_epi_tows$groups==3]),
                                   mean(zoop_epi_tows$ZoopDensity_No.pL_1[zoop_epi_tows$site=="lit" & zoop_epi_tows$groups==1]),
                                   mean(zoop_epi_tows$ZoopDensity_No.pL_1[zoop_epi_tows$site=="lit" & zoop_epi_tows$groups==2]),
                                   mean(zoop_epi_tows$ZoopDensity_No.pL_1[zoop_epi_tows$site=="lit" & zoop_epi_tows$groups==3]))

#plot response variable (euclidean distances) against environmental/biological data
#jpeg("Figures/2019-2020_euclidean_dist_daily_sums_vs_epiDO.jpg", width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$DO_0.1m,euclidean_drivers_df$euc_dist, xlab="Epi DO", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$DO_0.1m[1],euclidean_drivers_df$euc_dist[1], bg="orangered3",pch=21,cex=3)
points(euclidean_drivers_df$DO_0.1m[2],euclidean_drivers_df$euc_dist[2], bg="cadetblue",pch=21,cex=3)
points(euclidean_drivers_df$DO_0.1m[3],euclidean_drivers_df$euc_dist[3], bg="mediumpurple4",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3"),pch=21, pt.bg=c("orangered3","cadetblue","mediumpurple4"),bty = "n",cex=1.4)
#dev.off()

#jpeg("Figures/2019-2020_euclidean_dist_daily_sums_vs_hypoDO.jpg", width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$DO_9m,euclidean_drivers_df$euc_dist, xlab="Hypo DO", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$DO_9m[1],euclidean_drivers_df$euc_dist[1], bg="orangered3",pch=21,cex=3)
points(euclidean_drivers_df$DO_9m[2],euclidean_drivers_df$euc_dist[2], bg="cadetblue",pch=21,cex=3)
points(euclidean_drivers_df$DO_9m[3],euclidean_drivers_df$euc_dist[3], bg="mediumpurple4",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3"),pch=21, pt.bg=c("orangered3","cadetblue","mediumpurple4"),bty = "n",cex=1.4)
#dev.off()

#jpeg("Figures/2019-2020_euclidean_dist_daily_sums_vs_epi_temp.jpg", width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$temp_0.1m,euclidean_drivers_df$euc_dist, xlab="Epi temp", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$temp_0.1m[1],euclidean_drivers_df$euc_dist[1], bg="orangered3",pch=21,cex=3)
points(euclidean_drivers_df$temp_0.1m[2],euclidean_drivers_df$euc_dist[2], bg="cadetblue",pch=21,cex=3)
points(euclidean_drivers_df$temp_0.1m[3],euclidean_drivers_df$euc_dist[3], bg="mediumpurple4",pch=21,cex=3)
legend("bottomleft",legend=c("day1","day2","day3"),pch=21, pt.bg=c("orangered3","cadetblue","mediumpurple4"),bty = "n",cex=1.4)
#dev.off()

#jpeg("Figures/2019-2020_euclidean_dist_daily_sums_vs_hypo_temp.jpg", width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$temp_9m,euclidean_drivers_df$euc_dist, xlab="Hypo temp", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$temp_9m[1],euclidean_drivers_df$euc_dist[1], bg="orangered3",pch=21,cex=3)
points(euclidean_drivers_df$temp_9m[2],euclidean_drivers_df$euc_dist[2], bg="cadetblue",pch=21,cex=3)
points(euclidean_drivers_df$temp_9m[3],euclidean_drivers_df$euc_dist[3], bg="mediumpurple4",pch=21,cex=3)
legend("topright",legend=c("day1","day2","day3"),pch=21, pt.bg=c("orangered3","cadetblue","mediumpurple4"),bty = "n",cex=1.4)
#dev.off()

#jpeg("Figures/2019-2020_euclidean_dist_daily_sums_vs_epi_sp_cond.jpg", width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$sp_cond_0.1m,euclidean_drivers_df$euc_dist, xlab="Epi sp. cond", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$sp_cond_0.1m[1],euclidean_drivers_df$euc_dist[1], bg="orangered3",pch=21,cex=3)
points(euclidean_drivers_df$sp_cond_0.1m[2],euclidean_drivers_df$euc_dist[2], bg="cadetblue",pch=21,cex=3)
points(euclidean_drivers_df$sp_cond_0.1m[3],euclidean_drivers_df$euc_dist[3], bg="mediumpurple4",pch=21,cex=3)
legend("bottomleft",legend=c("day1","day2","day3"),pch=21, pt.bg=c("orangered3","cadetblue","mediumpurple4"),bty = "n",cex=1.4)
#dev.off()

#jpeg("Figures/2019-2020_euclidean_dist_daily_sums_vs_hypo_sp_cond.jpg", width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$sp_cond_9m,euclidean_drivers_df$euc_dist, xlab="Hypo sp. cond", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$sp_cond_9m[1],euclidean_drivers_df$euc_dist[1], bg="orangered3",pch=21,cex=3)
points(euclidean_drivers_df$sp_cond_9m[2],euclidean_drivers_df$euc_dist[2], bg="cadetblue",pch=21,cex=3)
points(euclidean_drivers_df$sp_cond_9m[3],euclidean_drivers_df$euc_dist[3], bg="mediumpurple4",pch=21,cex=3)
legend("bottomleft",legend=c("day1","day2","day3"),pch=21, pt.bg=c("orangered3","cadetblue","mediumpurple4"),bty = "n",cex=1.4)
#dev.off()

#jpeg("Figures/2019-2020_euclidean_dist_daily_sums_vs_epi_chla.jpg", width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$chl_0.1m,euclidean_drivers_df$euc_dist, xlab="Epi chla", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$chl_0.1m[1],euclidean_drivers_df$euc_dist[1], bg="orangered3",pch=21,cex=3)
points(euclidean_drivers_df$chl_0.1m[2],euclidean_drivers_df$euc_dist[2], bg="cadetblue",pch=21,cex=3)
points(euclidean_drivers_df$chl_0.1m[3],euclidean_drivers_df$euc_dist[3], bg="mediumpurple4",pch=21,cex=3)
legend("bottomleft",legend=c("day1","day2","day3"),pch=21, pt.bg=c("orangered3","cadetblue","mediumpurple4"),bty = "n",cex=1.4)
#dev.off()

#jpeg("Figures/2019-2020_euclidean_dist_daily_sums_vs_hypo_chla.jpg", width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$chl_9m,euclidean_drivers_df$euc_dist, xlab="Hypo chla", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$chl_9m[1],euclidean_drivers_df$euc_dist[1], bg="orangered3",pch=21,cex=3)
points(euclidean_drivers_df$chl_9m[2],euclidean_drivers_df$euc_dist[2], bg="cadetblue",pch=21,cex=3)
points(euclidean_drivers_df$chl_9m[3],euclidean_drivers_df$euc_dist[3], bg="mediumpurple4",pch=21,cex=3)
legend("bottomleft",legend=c("day1","day2","day3"),pch=21, pt.bg=c("orangered3","cadetblue","mediumpurple4"),bty = "n",cex=1.4)
#dev.off()

#jpeg("Figures/2019-2020_euclidean_dist_daily_sums_vs_pel_zoop_dens.jpg", width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$avg_dens[1:3],euclidean_drivers_df$euc_dist[1:3], xlab="Avg pelagic zoop dens", ylab="distance",cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$avg_dens[1],euclidean_drivers_df$euc_dist[1], bg="orangered3",pch=21,cex=3)
points(euclidean_drivers_df$avg_dens[2],euclidean_drivers_df$euc_dist[2], bg="cadetblue",pch=21,cex=3)
points(euclidean_drivers_df$avg_dens[3],euclidean_drivers_df$euc_dist[3], bg="mediumpurple4",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3"),pch=21, pt.bg=c("orangered3","cadetblue","mediumpurple4"),bty = "n",cex=1.4)
#dev.off()

#jpeg("Figures/2019-2020_euclidean_dist_daily_sums_vs_lit_zoop_dens.jpg", width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$avg_dens[4:6],euclidean_drivers_df$euc_dist[4:6], xlab="Avg littoral zoop dens", ylab="distance",cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$avg_dens[4],euclidean_drivers_df$euc_dist[4], bg="orangered3",pch=21,cex=3)
points(euclidean_drivers_df$avg_dens[5],euclidean_drivers_df$euc_dist[5], bg="cadetblue",pch=21,cex=3)
points(euclidean_drivers_df$avg_dens[6],euclidean_drivers_df$euc_dist[6], bg="mediumpurple4",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3"),pch=21, pt.bg=c("orangered3","cadetblue","mediumpurple4"),bty = "n",cex=1.4)
#dev.off()

#should have done this at the very beginning, but creating a scree plot to make sure 3 is the right number of axes to be using 
dimcheckMDS(zoop_temporal_dens_trans, distance = "bray", k = 6, trymax = 20, autotransform = TRUE)

  