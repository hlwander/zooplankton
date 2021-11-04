#NMDS with 2019 data
#created 25Nov2021

#read in libraries
pacman::p_load(dplyr, vegan, labdsv)

#function to count characters starting at the end of the string
substrEnd <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#read in 2019 data
zoop<- read.csv('SummaryStats/FCR_ZooplanktonSummary2019.csv',header = TRUE)

#ignore 20 um and horizontal trap samples
zoop <- zoop %>% filter(mesh_size_μm >20, na.rm=TRUE)

#select zoop density cols
zoop <- zoop %>% select("sample_ID","site_no","collect_date","DepthOfTow_m","Hour","ZoopDensity_No.pL","Calanoida_density_NopL",
                          "Cyclopoida_density_NopL", "Keratella_density_NopL","Kellicottia_density_NopL", "Bosmina_density_NopL",
                          "Daphnia_density_NopL", "Ceriodaphnia_density_NopL","nauplius_density_NopL", "Collothecidae_density_NopL",
                          "Synchaetidae_density_NopL", "Conochilidae_density_NopL")

#create df for temporal epi tows
zoop_epi_tows <- zoop[zoop$site_no=="BVR_l" | zoop$site_no=="BVR_50_p", ] %>%
  mutate(Hour=substr(Hour,1,2)) %>% 
  select(!c(DepthOfTow_m, sample_ID)) %>% group_by(site_no,collect_date,Hour) %>%
  summarise(across(everything(),list(mean)))

zoop_epi_tows$time <-ifelse(zoop_epi_tows$Hour=="12" | zoop_epi_tows$Hour=="11", "noon", ifelse(
  zoop_epi_tows$Hour =="0:" | zoop_epi_tows$Hour =="23", "midnight",ifelse(zoop_epi_tows$Hour=="18"|
  zoop_epi_tows$Hour=="19" | zoop_epi_tows$Hour=="20" | zoop_epi_tows$Hour=="21", "sunset", "sunrise")))

zoop_epi_tows$site <- ifelse(substrEnd(zoop_epi_tows$site_no,1)=="p","pel","lit")

zoop_epi_tows$day <- ifelse(zoop_epi_tows$collect_date=="2019-07-10" | 
                           zoop_epi_tows$collect_date=="2019-07-11","1", "2")

#subset zoop_epi_tows by day
zoop_epi_tows1 <- zoop_epi_tows[zoop_epi_tows$day=="1",]
zoop_epi_tows2 <- zoop_epi_tows[zoop_epi_tows$day=="2",]

#------------------------------------------------------------------------------#
# set up data for NMDS
#relies on rank orders for ordination, no assumptions of linear relationship
zoop_temporal_dens1 <- zoop_epi_tows[zoop_epi_tows$collect_date=="2019-07-10" | 
                                     zoop_epi_tows$collect_date=="2019-07-11" ,c(5:15)]

zoop_temporal_dens2 <- zoop_epi_tows[zoop_epi_tows$collect_date=="2019-07-24" | 
                                     zoop_epi_tows$collect_date=="2019-07-25" ,c(5:15)]

#transforming data - hellinger transformation because gives low weight to low/zero values
#converts species abundances from absolute to relative - use w/ bray curtis
zoop_temporal_dens_trans1 <- hellinger(zoop_temporal_dens1)
zoop_temporal_dens_trans2 <- hellinger(zoop_temporal_dens2)

#-------------------------------------------------------------------------------#
#                  Bray-curtis dissimilarity --> NMDS figs                      #
#-------------------------------------------------------------------------------#
par(ask=TRUE)
NMDS_temporal_bray1 <- metaMDS(zoop_temporal_dens_trans1, distance='bray', k=3, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_temporal_bray1$stress

NMDS_temporal_bray2 <- metaMDS(zoop_temporal_dens_trans2, distance='bray', k=3, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_temporal_bray2$stress


#  Points divided into time groups
#jpeg("Figures/2019_NMDS_1v2_bray_temporal_day1.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_temporal_bray1,display = c('sites','species'),choices = c(1,2),type = "n",xlim=c(-0.6,1)) 
points(NMDS_temporal_bray1$points[zoop_epi_tows1$time=='sunrise',1],NMDS_temporal_bray1$points[zoop_epi_tows1$time=='sunrise',2], pch=21,bg="#CC79A7")
points(NMDS_temporal_bray1$points[zoop_epi_tows1$time=='noon',1], NMDS_temporal_bray1$points[zoop_epi_tows1$time=='noon',2], pch=21,bg="#F0E442")
points(NMDS_temporal_bray1$points[zoop_epi_tows1$time=='sunset',1], NMDS_temporal_bray1$points[zoop_epi_tows1$time=='sunset',2], pch=21,bg="#009E73")
points(NMDS_temporal_bray1$points[zoop_epi_tows1$time=='midnight',1], NMDS_temporal_bray1$points[zoop_epi_tows1$time=='midnight',2], pch=21,bg="#0072B2")
legend("bottomright", legend=c('sunrise','noon', 'sunset', 'midnight'), pch=21, pt.bg=c('#CC79A7', '#F0E442','#009E73','#0072B2'),bty = "n") 
ordihull(ord, zoop_epi_tows1$time, display = "sites", draw = c("polygon"),
         col = c("#0072B2", "#F0E442","#CC79A7","#009E73"), alpha = 75,cex = 2)
text(NMDS_temporal_bray1$species[,1],NMDS_temporal_bray1$species[,2], labels = c("Calanoida","Cyclopoida","Keratella","Kellicottia",
     "Bosmina","Daphnia","Ceriodaphnia","nauplius","Collothecidae","Synchaetidae","Conochilidae"), cex=0.9)
#dev.off()

#jpeg("Figures/2019_NMDS_1v3_bray_temporal_day1.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_temporal_bray1,display = c('sites','species'),choices = c(1,3),type = "n",xlim=c(-0.6,1)) 
points(NMDS_temporal_bray1$points[zoop_epi_tows1$time=='sunrise',1],NMDS_temporal_bray1$points[zoop_epi_tows1$time=='sunrise',3], pch=21,bg="#CC79A7")
points(NMDS_temporal_bray1$points[zoop_epi_tows1$time=='noon',1], NMDS_temporal_bray1$points[zoop_epi_tows1$time=='noon',3], pch=21,bg="#F0E442")
points(NMDS_temporal_bray1$points[zoop_epi_tows1$time=='sunset',1], NMDS_temporal_bray1$points[zoop_epi_tows1$time=='sunset',3], pch=21,bg="#009E73")
points(NMDS_temporal_bray1$points[zoop_epi_tows1$time=='midnight',1], NMDS_temporal_bray1$points[zoop_epi_tows1$time=='midnight',3], pch=21,bg="#0072B2")
legend("bottomright", legend=c('sunrise','noon', 'sunset', 'midnight'), pch=21, pt.bg=c('#CC79A7', '#F0E442','#009E73','#0072B2'),bty = "n") 
ordihull(ord, zoop_epi_tows1$time, display = "sites", draw = c("polygon"),
         col = c("#0072B2", "#F0E442","#CC79A7","#009E73"), alpha = 75,cex = 2)
text(NMDS_temporal_bray1$species[,1],NMDS_temporal_bray1$species[,3], labels = c("Calanoida","Cyclopoida","Keratella","Kellicottia",
     "Bosmina","Daphnia","Ceriodaphnia","nauplius","Collothecidae","Synchaetidae","Conochilidae"), cex=0.9)
#dev.off()

#jpeg("Figures/2019_NMDS_1v2_bray_temporal_day2.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_temporal_bray2,display = c('sites','species'),choices = c(1,2),type = "n",xlim=c(-0.6,1)) 
points(NMDS_temporal_bray2$points[zoop_epi_tows2$time=='sunrise',1],NMDS_temporal_bray2$points[zoop_epi_tows2$time=='sunrise',2], pch=21,bg="#CC79A7")
points(NMDS_temporal_bray2$points[zoop_epi_tows2$time=='noon',1], NMDS_temporal_bray2$points[zoop_epi_tows2$time=='noon',2], pch=21,bg="#F0E442")
points(NMDS_temporal_bray2$points[zoop_epi_tows2$time=='sunset',1], NMDS_temporal_bray2$points[zoop_epi_tows2$time=='sunset',2], pch=21,bg="#009E73")
points(NMDS_temporal_bray2$points[zoop_epi_tows2$time=='midnight',1], NMDS_temporal_bray2$points[zoop_epi_tows2$time=='midnight',2], pch=21,bg="#0072B2")
legend("bottomright", legend=c('sunrise','noon', 'sunset', 'midnight'), pch=21, pt.bg=c('#CC79A7', '#F0E442','#009E73','#0072B2'),bty = "n") 
ordihull(ord, zoop_epi_tows2$time, display = "sites", draw = c("polygon"),
         col = c("#0072B2", "#F0E442","#CC79A7","#009E73"), alpha = 75,cex = 2)
text(NMDS_temporal_bray2$species[,1],NMDS_temporal_bray2$species[,2], labels = c("Calanoida","Cyclopoida","Keratella","Kellicottia",
     "Bosmina","Daphnia","Ceriodaphnia","nauplius","Collothecidae","Synchaetidae","Conochilidae"), cex=0.9)
#dev.off()

#jpeg("Figures/2019_NMDS_1v3_bray_temporal_day2.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_temporal_bray2,display = c('sites','species'),choices = c(1,3),type = "n",xlim=c(-0.65,1)) 
points(NMDS_temporal_bray2$points[zoop_epi_tows2$time=='sunrise',1],NMDS_temporal_bray2$points[zoop_epi_tows2$time=='sunrise',3], pch=21,bg="#CC79A7")
points(NMDS_temporal_bray2$points[zoop_epi_tows2$time=='noon',1], NMDS_temporal_bray2$points[zoop_epi_tows2$time=='noon',3], pch=21,bg="#F0E442")
points(NMDS_temporal_bray2$points[zoop_epi_tows2$time=='sunset',1], NMDS_temporal_bray2$points[zoop_epi_tows2$time=='sunset',3], pch=21,bg="#009E73")
points(NMDS_temporal_bray2$points[zoop_epi_tows2$time=='midnight',1], NMDS_temporal_bray2$points[zoop_epi_tows2$time=='midnight',3], pch=21,bg="#0072B2")
legend("bottomright", legend=c('sunrise','noon', 'sunset', 'midnight'), pch=21, pt.bg=c('#CC79A7', '#F0E442','#009E73','#0072B2'),bty = "n") 
ordihull(ord, zoop_epi_tows2$time, display = "sites", draw = c("polygon"),
         col = c("#0072B2", "#F0E442","#CC79A7","#009E73"), alpha = 75,cex = 2)
text(NMDS_temporal_bray2$species[,1],NMDS_temporal_bray2$species[,2], labels = c("Calanoida","Cyclopoida","Keratella","Kellicottia",
     "Bosmina","Daphnia","Ceriodaphnia","nauplius","Collothecidae","Synchaetidae","Conochilidae"), cex=0.9)
#dev.off()
                           
                           