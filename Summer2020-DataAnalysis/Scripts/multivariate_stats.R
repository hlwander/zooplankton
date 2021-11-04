#playing around with multivariate analyses using BVR 2020 zoop data

pacman::p_load(dplyr, vegan, labdsv)

#function to count characters starting at the end of the string
substrEnd <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#read in zoop summary file and get into format for analyses
zoop<- read.csv('SummaryStats/FCR_ZooplanktonSummary2020.csv',header = TRUE)

#ignore 20 um and horizontal trap samples
zoop <- zoop %>% filter(mesh_size_μm >20, na.rm=TRUE)

#select density and biomass cols to keep
zoop <- zoop %>% select("sample_ID","site_no","collect_date","DepthOfTow_m","Hour","ZoopDensity_No.pL","BiomassConcentration_ugpL","Calanoida_density_NopL","Calanoida_BiomassConcentration_ugpL",
  "Cyclopoida_density_NopL","Cyclopoida_BiomassConcentration_ugpL", "Keratella_density_NopL","Keratella_BiomassConcentration_ugpL",
  "Kellicottia_density_NopL","Kellicottia_BiomassConcentration_ugpL", "Bosmina_density_NopL","Bosmina_BiomassConcentration_ugpL",
  "Daphnia_density_NopL","Daphnia_BiomassConcentration_ugpL", "Ceriodaphnia_density_NopL","Ceriodaphnia_BiomassConcentration_ugpL",
  "nauplius_density_NopL","nauplius_BiomassConcentration_ugpL", "Collothecidae_density_NopL","Collothecidae_BiomassConcentration_ugpL",
  "Synchaetidae_density_NopL","Synchaetidae_BiomassConcentration_ugpL", "Conochilidae_density_NopL","Conochilidae_BiomassConcentration_ugpL")

#create dfs by depth (schindler) and time (tows)
zoop_schind <- zoop[zoop$site_no=="BVR_schind", ] %>% select(!c(site_no,Hour, collect_date, sample_ID)) %>%
  group_by(DepthOfTow_m) %>% summarise(across(everything(),list(mean))) %>% rename(Depth_m = DepthOfTow_m)

zoop_epi_tows <- zoop[zoop$site_no=="BVR_l" | zoop$site_no=="BVR_50_p", ] %>%
  mutate(Hour=substr(Hour,1,2)) %>% 
  select(!c(DepthOfTow_m, sample_ID)) %>% group_by(site_no,collect_date,Hour) %>%
  summarise(across(everything(),list(mean)))

zoop_epi_tows$time <-ifelse(zoop_epi_tows$Hour=="12", "noon", ifelse(
  zoop_epi_tows$Hour =="0:", "midnight",ifelse(zoop_epi_tows$Hour=="18"|zoop_epi_tows$Hour=="19" |
  zoop_epi_tows$Hour=="20" | zoop_epi_tows$Hour=="21", "sunset", "sunrise")))

zoop_epi_tows$site <- ifelse(substrEnd(zoop_epi_tows$site_no,1)=="p","pel","lit")

#------------------------------------------------------------------------------#
#pull in driver data

#now read in the ctd and fp data from 12-13 Aug
fp <- read.csv('RawData/Fluora.csv') %>% select(!c(X, Date, DOY)) %>% filter(DateTime=="2020-08-12 16:00:00") 
ctd<- read.csv('./RawData/081220_bvr50.csv', header=TRUE) %>% mutate(Date = as.Date(as.character(Date))) %>% 
  rename(DateTime = Date) %>% select(!Flag)

depths = c(0.1, seq(1, 10, by = 1))
fp.final<-data.frame()
ctd.final<-data.frame() 

for (i in 1:length(depths)){
  fp_layer <- fp %>% group_by(DateTime) %>% slice(which.min(abs(as.numeric(Depth_m) - depths[i])))
  ctd_layer <- ctd %>% group_by(DateTime) %>% slice(which.min(abs(as.numeric(Depth_m) - depths[i])))
  # Bind each of the data layers together.
  fp.final = bind_rows(fp.final, fp_layer)
  ctd.final = bind_rows(ctd.final, ctd_layer)
}

#manually set depths so can merge datasets
fp.final$Depth_m <- c(0.1, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0)
ctd.final$Depth_m <- c(0.1, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0)

#drop ph because it definitely was NOT 14 ever
ctd.final <- ctd.final[-c(10)]

#add fp and ctd data to zoop_schind
depth_df <- merge(zoop_schind, fp.final, by="Depth_m")
depth_df <- merge(depth_df, ctd.final, by="Depth_m")

#add in epi vs hypo and crust vs rot for NMDS
depth_df$zone <- ifelse(depth_df$Depth_m > 4, "hypo","epi")

#------------------------------------------------------------------------------#
#Nonmetric multidimensional scaling (NMDS) 
#relies on rank orders for ordination, no assumptions of linear relationship
zoop_dens <- depth_df[,c(4,6,8,10,12,14,16,18,20,22,24)]
zoop_biom <- depth_df[,c(5,7,9,11,13,15,17,19,21,23,25)]
zoop_temporal_dens <- zoop_epi_tows[,c(6,8,10,12,14,16,18,20,22,24,26)]
zoop_temporal_biom <- zoop_epi_tows[,c(7,9,11,13,15,17,19,21,23,25,27)]

#transforming data - hellinger transformation because gives low weight to low/zero values
#converts species abundances from absolute to relative - use w/ bray curtis
zoop_dens_trans <- hellinger(zoop_dens)
zoop_biom_trans <- hellinger(zoop_biom)
zoop_temporal_dens_trans <- hellinger(zoop_temporal_dens)
zoop_temporal_biom_trans <- hellinger(zoop_temporal_biom)

#-------------------------------------------------------------------------------#
#              Bray-curtis dissimilarity --> depth NMDS figs                    #
#-------------------------------------------------------------------------------#
par(ask=TRUE)
NMDS_depth_bray <- metaMDS(zoop_dens_trans, distance='bray', k=3, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_depth_bray$stress

#  Points divided into epi vs. hypo
#jpeg("Figures/2020_NMDS_1v2_bray_epivshypo.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- plot(NMDS_depth_bray, display=c('sites',"species"),choices=c(1,2), type='n')
points(NMDS_depth_bray$points[depth_df$zone=='epi',1], NMDS_depth_bray$points[depth_df$zone=='epi',2], pch=21,bg='#B47846')
points(NMDS_depth_bray$points[depth_df$zone=='hypo',1], NMDS_depth_bray$points[depth_df$zone=='hypo',2], pch=21,bg='#4682B4')
text(NMDS_depth_bray$species[,1],NMDS_depth_bray$species[,2], labels = c("Calanoida","Cyclopoida","Keratella","Kellicottia",
    "Bosmina","Daphnia","Ceriodaphnia","nauplius","Synchaetidae","Conochilidae"), cex=0.9)
legend("bottomright", legend=c('epi','hypo'), pch=21, pt.bg=c('#B47846', '#4682B4'),bty='n') 

ordihull(ord, depth_df$zone, display = "sites", draw = c("polygon"),
         col = c("#B47846","#4682B4"), alpha = 75,cex = 2)
#dev.off()

#jpeg("Figures/2020_NMDS_1v3_bray_epivshypo.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- plot(NMDS_depth_bray, display=c('sites',"species"),choices=c(1,3), type='n', ylim=c(-0.4,0.9))
points(NMDS_depth_bray$points[depth_df$zone=='epi',1], NMDS_depth_bray$points[depth_df$zone=='epi',3], pch=21,bg='#B47846')
points(NMDS_depth_bray$points[depth_df$zone=='hypo',1], NMDS_depth_bray$points[depth_df$zone=='hypo',3], pch=21,bg='#4682B4')
text(NMDS_depth_bray$species[,1],NMDS_depth_bray$species[,3], labels = c("Calanoida","Cyclopoida","Keratella","Kellicottia",
     "Bosmina","Daphnia","Ceriodaphnia","nauplius","Synchaetidae","Conochilidae"), cex=0.9)
legend("bottomright", legend=c('epi','hypo'), pch=21, pt.bg=c('#B47846', '#4682B4'),bty='n') 

ordihull(ord, depth_df$zone, display = "sites", draw = c("polygon"),
         col = c("#B47846","#4682B4"), alpha = 75,cex = 2)
#dev.off()

#jpeg("Figures/2020_NMDS_2v3_bray_epivshypo.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- plot(NMDS_depth_bray, display=c('sites',"species"),choices=c(2,3), type='n')
points(NMDS_depth_bray$points[depth_df$zone=='epi',2], NMDS_depth_bray$points[depth_df$zone=='epi',3], pch=21,bg='#B47846')
points(NMDS_depth_bray$points[depth_df$zone=='hypo',2], NMDS_depth_bray$points[depth_df$zone=='hypo',3], pch=21,bg='#4682B4')
text(NMDS_depth_bray$species[,2],NMDS_depth_bray$species[,3], labels = c("Calanoida","Cyclopoida","Keratella","Kellicottia",
     "Bosmina","Daphnia","Ceriodaphnia","nauplius","Synchaetidae","Conochilidae"), cex=0.9)
legend("bottomright", legend=c('epi','hypo'), pch=21, pt.bg=c('#B47846', '#4682B4'),bty='n') 

ordihull(ord, depth_df$zone, display = "sites", draw = c("polygon"),
         col = c("#B47846","#4682B4"), alpha = 75,cex = 2)
#dev.off()

#looking at how species comp changes across sites
#ord <- ordiplot(NMDS_depth_bray,display = c('sites'),choices = c(1,2),type = "points")
#points(NMDS_depth_bray$points[depth_df$zone=='epi',1], NMDS_depth_bray$points[depth_df$zone=='epi',3], pch=21,bg='black')
#points(NMDS_depth_bray$points[depth_df$zone=='hypo',1], NMDS_depth_bray$points[depth_df$zone=='hypo',3], pch=21,bg='grey75')
#legend("topleft", legend=c('epi','hypo'), pch=21, pt.bg=c('black', 'grey75'),bty = "n") 
#legend("bottomright",legend = c("k = 3"),cex = 1.4,bty = 'n')

#look at environmental data 
fit <- envfit(NMDS_depth_bray, depth_df[,c(27:32,35:44)], permu=999)

#jpeg("Figures/2020_NMDS_1v2_bray_epivshypo_envcorr.jpg", width = 6, height = 5, units = "in",res = 300)
plot(NMDS_depth_bray, display="sites",choices=c(1,2))
points(NMDS_depth_bray$points[depth_df$zone=='epi',1], NMDS_depth_bray$points[depth_df$zone=='epi',2], pch=21,bg="#B47846")
points(NMDS_depth_bray$points[depth_df$zone=='hypo',1], NMDS_depth_bray$points[depth_df$zone=='hypo',2], pch=21,bg="#4682B4")
plot(fit, p.max=0.05) # only display variables that are significant
legend("topleft", legend=c('epi','hypo'), pch=21, pt.bg=c("#B47846", "#4682B4"),bty ='n') 
#dev.off()

#jpeg("Figures/2020_NMDS_2v3_bray_epivshypo_envcorr.jpg", width = 6, height = 5, units = "in",res = 300)
plot(NMDS_depth_bray, display="sites",choices=c(2,3))
points(NMDS_depth_bray$points[depth_df$zone=='epi',2], NMDS_depth_bray$points[depth_df$zone=='epi',3], pch=21,bg="#B47846")
points(NMDS_depth_bray$points[depth_df$zone=='hypo',2], NMDS_depth_bray$points[depth_df$zone=='hypo',3], pch=21,bg="#4682B4")
plot(fit, p.max=0.05) # only display variables that are significant
legend("topleft", legend=c('epi','hypo'), pch=21, pt.bg=c("#B47846", "#4682B4"),bty ='n') 
#dev.off()

#jpeg("Figures/2020_NMDS_1v3_bray_epivshypo_envcorr.jpg", width = 6, height = 5, units = "in",res = 300)
plot(NMDS_depth_bray, display="sites",choices=c(1,3))
points(NMDS_depth_bray$points[depth_df$zone=='epi',1], NMDS_depth_bray$points[depth_df$zone=='epi',3], pch=21,bg="#B47846")
points(NMDS_depth_bray$points[depth_df$zone=='hypo',1], NMDS_depth_bray$points[depth_df$zone=='hypo',3], pch=21,bg="#4682B4")
plot(fit, p.max=0.05) # only display variables that are significant
legend("topleft", legend=c('epi','hypo'), pch=21, pt.bg=c("#B47846", "#4682B4"),bty ='n') 
#dev.off()


#-------------------------------------------------------------------------------#
#             Morisita horn dissimilarity --> depth NMDS figs                   #
#-------------------------------------------------------------------------------#
NMDS_depth_horn <- metaMDS(zoop_dens_trans, distance='horn', k=3, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_depth_horn$stress

#jpeg("Figures/2020_NMDS_1v2_mhorn_epivshypo.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- plot(NMDS_depth_horn, display=c('sites',"species"),choices=c(1,2), type='n')
points(NMDS_depth_horn$points[depth_df$zone=='epi',1], NMDS_depth_horn$points[depth_df$zone=='epi',2], pch=21,bg='#B47846')
points(NMDS_depth_horn$points[depth_df$zone=='hypo',1], NMDS_depth_horn$points[depth_df$zone=='hypo',2], pch=21,bg='#4682B4')
text(NMDS_depth_horn$species[,1],NMDS_depth_horn$species[,2], labels = c("Calanoida","Cyclopoida","Keratella","Kellicottia",
      "Bosmina","Daphnia","Ceriodaphnia","nauplius","Synchaetidae","Conochilidae"), cex=0.9)
legend("bottomright", legend=c('epi','hypo'), pch=21, pt.bg=c('#B47846', '#4682B4'),bty='n') 

ordihull(ord, depth_df$zone, display = "sites", draw = c("polygon"),
         col = c("#B47846","#4682B4"), alpha = 75,cex = 2)
#dev.off()

#jpeg("Figures/2020_NMDS_1v3_mhorn_epivshypo.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- plot(NMDS_depth_horn, display=c('sites',"species"),choices=c(1,3), type='n', ylim=c(-0.4,0.9))
points(NMDS_depth_horn$points[depth_df$zone=='epi',1], NMDS_depth_horn$points[depth_df$zone=='epi',3], pch=21,bg='#B47846')
points(NMDS_depth_horn$points[depth_df$zone=='hypo',1], NMDS_depth_horn$points[depth_df$zone=='hypo',3], pch=21,bg='#4682B4')
text(NMDS_depth_horn$species[,1],NMDS_depth_horn$species[,3], labels = c("Calanoida","Cyclopoida","Keratella","Kellicottia",
     "Bosmina","Daphnia","Ceriodaphnia","nauplius","Synchaetidae","Conochilidae"), cex=0.9)
legend("bottomright", legend=c('epi','hypo'), pch=21, pt.bg=c('#B47846', '#4682B4'),bty='n') 

ordihull(ord, depth_df$zone, display = "sites", draw = c("polygon"),
         col = c("#B47846","#4682B4"), alpha = 75,cex = 2)
#dev.off()

#jpeg("Figures/2020_NMDS_2v3_mhorn_epivshypo.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- plot(NMDS_depth_horn, display=c('sites',"species"),choices=c(2,3), type='n')
points(NMDS_depth_horn$points[depth_df$zone=='epi',2], NMDS_depth_horn$points[depth_df$zone=='epi',3], pch=21,bg='#B47846')
points(NMDS_depth_horn$points[depth_df$zone=='hypo',2], NMDS_depth_horn$points[depth_df$zone=='hypo',3], pch=21,bg='#4682B4')
text(NMDS_depth_horn$species[,2],NMDS_depth_horn$species[,3], labels = c("Calanoida","Cyclopoida","Keratella","Kellicottia",
      "Bosmina","Daphnia","Ceriodaphnia","nauplius","Synchaetidae","Conochilidae"), cex=0.9)
legend("bottomright", legend=c('epi','hypo'), pch=21, pt.bg=c('#B47846', '#4682B4'),bty='n') 

ordihull(ord, depth_df$zone, display = "sites", draw = c("polygon"),
         col = c("#B47846","#4682B4"), alpha = 75,cex = 2)
#dev.off()

#look at environmental data 
fit_horn <- envfit(NMDS_depth_horn, depth_df[,c(27:32,35:44)], permu=999)

#jpeg("Figures/2020_NMDS_1v2_mhorn_epivshypo_envcorr.jpg", width = 6, height = 5, units = "in",res = 300)
plot(NMDS_depth_horn, display="sites",choices=c(1,2))
points(NMDS_depth_horn$points[depth_df$zone=='epi',1], NMDS_depth_horn$points[depth_df$zone=='epi',2], pch=21,bg="#B47846")
points(NMDS_depth_horn$points[depth_df$zone=='hypo',1], NMDS_depth_horn$points[depth_df$zone=='hypo',2], pch=21,bg="#4682B4")
plot(fit_horn, p.max=0.05) # only display variables that are significant
legend("topleft", legend=c('epi','hypo'), pch=21, pt.bg=c("#B47846", "#4682B4"),bty ='n') 
#dev.off()

#jpeg("Figures/2020_NMDS_2v3_mhorn_epivshypo_envcorr.jpg", width = 6, height = 5, units = "in",res = 300)
plot(NMDS_depth_horn, display="sites",choices=c(2,3))
points(NMDS_depth_horn$points[depth_df$zone=='epi',2], NMDS_depth_horn$points[depth_df$zone=='epi',3], pch=21,bg="#B47846")
points(NMDS_depth_horn$points[depth_df$zone=='hypo',2], NMDS_depth_horn$points[depth_df$zone=='hypo',3], pch=21,bg="#4682B4")
plot(fit_horn, p.max=0.05) # only display variables that are significant
legend("topleft", legend=c('epi','hypo'), pch=21, pt.bg=c("#B47846", "#4682B4"),bty ='n') 
#dev.off()

#jpeg("Figures/2020_NMDS_1v3_mhorn_epivshypo_envcorr.jpg", width = 6, height = 5, units = "in",res = 300)
plot(NMDS_depth_horn, display="sites",choices=c(1,3))
points(NMDS_depth_horn$points[depth_df$zone=='epi',1], NMDS_depth_horn$points[depth_df$zone=='epi',3], pch=21,bg="#B47846")
points(NMDS_depth_horn$points[depth_df$zone=='hypo',1], NMDS_depth_horn$points[depth_df$zone=='hypo',3], pch=21,bg="#4682B4")
plot(fit_horn, p.max=0.05) # only display variables that are significant
legend("topleft", legend=c('epi','hypo'), pch=21, pt.bg=c("#B47846", "#4682B4"),bty ='n') 
#dev.off()


#-------------------------------------------------------------------------------#
#             Bray-curtis dissimilarity --> temporal NMDS figs                  #
#-------------------------------------------------------------------------------#
par(ask=TRUE)
NMDS_temporal_bray <- metaMDS(zoop_temporal_dens_trans, distance='bray', k=3, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_temporal_bray$stress

#  Points divided into time groups
#jpeg("Figures/2020_NMDS_1v2_bray_temporal.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_temporal_bray,display = c('sites','species'),choices = c(1,2),type = "n",xlim=c(-0.6,1)) 
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

#jpeg("Figures/2020_NMDS_1v3_bray_temporal.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_temporal_bray,display = c('sites','species'),choices = c(1,3),type = "n",xlim=c(-0.6,1)) 
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

#jpeg("Figures/2020_NMDS_2v3_bray_temporal.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_temporal_bray,display = c('sites','species'),choices = c(2,3),type = "n") 
points(NMDS_temporal_bray$points[zoop_epi_tows$time=='sunrise',2], NMDS_temporal_bray$points[zoop_epi_tows$time=='sunrise',3], pch=21,bg="#CC79A7")
points(NMDS_temporal_bray$points[zoop_epi_tows$time=='noon',2], NMDS_temporal_bray$points[zoop_epi_tows$time=='noon',3], pch=21,bg="#F0E442")
points(NMDS_temporal_bray$points[zoop_epi_tows$time=='sunset',2], NMDS_temporal_bray$points[zoop_epi_tows$time=='sunset',3], pch=21,bg="#009E73")
points(NMDS_temporal_bray$points[zoop_epi_tows$time=='midnight',2], NMDS_temporal_bray$points[zoop_epi_tows$time=='midnight',3], pch=21,bg="#0072B2")
text(NMDS_temporal_bray$species[,2],NMDS_temporal_bray$species[,3], labels = c("Calanoida","Cyclopoida","Keratella","Kellicottia",
      "Bosmina","Daphnia","Ceriodaphnia","nauplius","Collothecidae","Synchaetidae","Conochilidae"), cex=0.9)
legend("bottomright", legend=c('sunrise','noon', 'sunset', 'midnight'), pch=21, pt.bg=c('#CC79A7', '#F0E442','#009E73','#0072B2'),bty = "n") 

ordihull(ord, zoop_epi_tows$time, display = "sites", draw = c("polygon"),
         col = c("#0072B2", "#F0E442","#CC79A7","#009E73"), alpha = 75,cex = 2)
#dev.off()

#-------------------------------------------------------------------------------#
#            Morisita horn dissimilarity --> temporal NMDS figs                 #
#-------------------------------------------------------------------------------#
NMDS_temporal_horn <- metaMDS(zoop_temporal_dens, distance='horn', k=3, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_temporal_horn$stress

#jpeg("Figures/2020_NMDS_1v2_mhorn_temporal.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_temporal_horn,display = c('sites','species'),choices = c(1,2),type = "n",xlim=c(-0.6,1)) 
points(NMDS_temporal_horn$points[zoop_epi_tows$time=='sunrise',1], NMDS_temporal_horn$points[zoop_epi_tows$time=='sunrise',2], pch=21,bg="#CC79A7")
points(NMDS_temporal_horn$points[zoop_epi_tows$time=='noon',1], NMDS_temporal_horn$points[zoop_epi_tows$time=='noon',2], pch=21,bg="#F0E442")
points(NMDS_temporal_horn$points[zoop_epi_tows$time=='sunset',1], NMDS_temporal_horn$points[zoop_epi_tows$time=='sunset',2], pch=21,bg="#009E73")
points(NMDS_temporal_horn$points[zoop_epi_tows$time=='midnight',1], NMDS_temporal_horn$points[zoop_epi_tows$time=='midnight',2], pch=21,bg="#0072B2")
text(NMDS_temporal_horn$species[,1],NMDS_temporal_horn$species[,2], labels = c("Calanoida","Cyclopoida","Keratella","Kellicottia",
     "Bosmina","Daphnia","Ceriodaphnia","nauplius","Collothecidae","Synchaetidae","Conochilidae"), cex=0.9)
legend("bottomright", legend=c('sunrise','noon', 'sunset', 'midnight'), pch=21, pt.bg=c('#CC79A7', '#F0E442','#009E73','#0072B2'),bty = "n") 

ordihull(ord, zoop_epi_tows$time, display = "sites", draw = c("polygon"),
         col = c("#0072B2", "#F0E442","#CC79A7","#009E73"), alpha = 75,cex = 2)
#dev.off()

#jpeg("Figures/2020_NMDS_1v3_mhorn_temporal.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_temporal_horn,display = c('sites','species'),choices = c(1,3),type = "n",xlim=c(-0.6,1)) 
points(NMDS_temporal_horn$points[zoop_epi_tows$time=='sunrise',1], NMDS_temporal_horn$points[zoop_epi_tows$time=='sunrise',3], pch=21,bg="#CC79A7")
points(NMDS_temporal_horn$points[zoop_epi_tows$time=='noon',1], NMDS_temporal_horn$points[zoop_epi_tows$time=='noon',3], pch=21,bg="#F0E442")
points(NMDS_temporal_horn$points[zoop_epi_tows$time=='sunset',1], NMDS_temporal_horn$points[zoop_epi_tows$time=='sunset',3], pch=21,bg="#009E73")
points(NMDS_temporal_horn$points[zoop_epi_tows$time=='midnight',1], NMDS_temporal_horn$points[zoop_epi_tows$time=='midnight',3], pch=21,bg="#0072B2")
text(NMDS_temporal_horn$species[,1],NMDS_temporal_horn$species[,3], labels = c("Calanoida","Cyclopoida","Keratella","Kellicottia",
     "Bosmina","Daphnia","Ceriodaphnia","nauplius","Collothecidae","Synchaetidae","Conochilidae"), cex=0.9)
legend("bottomright", legend=c('sunrise','noon', 'sunset', 'midnight'), pch=21, pt.bg=c('#CC79A7', '#F0E442','#009E73','#0072B2'),bty = "n") 

ordihull(ord, zoop_epi_tows$time, display = "sites", draw = c("polygon"),
         col = c("#0072B2", "#F0E442","#CC79A7","#009E73"), alpha = 75,cex = 2)
#dev.off()

#jpeg("Figures/2020_NMDS_2v3_mhorn_temporal.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_temporal_horn,display = c('sites','species'),choices = c(2,3),type = "n") 
points(NMDS_temporal_horn$points[zoop_epi_tows$time=='sunrise',2], NMDS_temporal_horn$points[zoop_epi_tows$time=='sunrise',3], pch=21,bg="#CC79A7")
points(NMDS_temporal_horn$points[zoop_epi_tows$time=='noon',2], NMDS_temporal_horn$points[zoop_epi_tows$time=='noon',3], pch=21,bg="#F0E442")
points(NMDS_temporal_horn$points[zoop_epi_tows$time=='sunset',2], NMDS_temporal_horn$points[zoop_epi_tows$time=='sunset',3], pch=21,bg="#009E73")
points(NMDS_temporal_horn$points[zoop_epi_tows$time=='midnight',2], NMDS_temporal_horn$points[zoop_epi_tows$time=='midnight',3], pch=21,bg="#0072B2")
text(NMDS_temporal_horn$species[,2],NMDS_temporal_horn$species[,3], labels = c("Calanoida","Cyclopoida","Keratella","Kellicottia",
     "Bosmina","Daphnia","Ceriodaphnia","nauplius","Collothecidae","Synchaetidae","Conochilidae"), cex=0.9)
legend("bottomright", legend=c('sunrise','noon', 'sunset', 'midnight'), pch=21, pt.bg=c('#CC79A7', '#F0E442','#009E73','#0072B2'),bty = "n") 

ordihull(ord, zoop_epi_tows$time, display = "sites", draw = c("polygon"),
         col = c("#0072B2", "#F0E442","#CC79A7","#009E73"), alpha = 75,cex = 2)
#dev.off()


