#multivariate stats for 2019 - 2021 zoop data
#created 25Nov2021

#read in libraries
pacman::p_load(dplyr, vegan, labdsv, goeveg, rLakeAnalyzer, ggplot2,tidyr)

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
#converts species abundances from absolute to relative - use w/ bray curtis
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
#jpeg("Figures/2019-2020_NMDS_2v4_bray_days.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_temporal_bray,display = c('sites','species'),choices = c(2,4),type = "n") 
points(NMDS_temporal_bray$points[zoop_epi_tows$groups=="1",2], NMDS_temporal_bray$points[zoop_epi_tows$groups=="1",4], pch=21,bg="#008585")
points(NMDS_temporal_bray$points[zoop_epi_tows$groups=="2",2], NMDS_temporal_bray$points[zoop_epi_tows$groups=="2",4], pch=21,bg="#9BBAA0")
points(NMDS_temporal_bray$points[zoop_epi_tows$groups=="3",2], NMDS_temporal_bray$points[zoop_epi_tows$groups=="3",4], pch=21,bg="#FBF2C4")
points(NMDS_temporal_bray$points[zoop_epi_tows$groups=="4",2], NMDS_temporal_bray$points[zoop_epi_tows$groups=="4",4], pch=21,bg="#DEA868")
points(NMDS_temporal_bray$points[zoop_epi_tows$groups=="5",2], NMDS_temporal_bray$points[zoop_epi_tows$groups=="5",4], pch=21,bg="#C7522B")
legend("bottomright", legend=c('Day1','Day2','Day3','Day4','Day5'), pch=21, pt.bg=hcl.colors(5,"Geyser"), bty = "n") 
ordihull(ord, zoop_epi_tows$groups, display = "sites", draw = c("polygon"),
         col = hcl.colors(5,"Geyser"), alpha = 75,cex = 2)
text(NMDS_temporal_bray$species[,2],NMDS_temporal_bray$species[,4], labels = c("Calanoida","Cyclopoida","Keratella","Kellicottia",
     "Bosmina","Daphnia","Ceriodaphnia","nauplius","Collothecidae","Synchaetidae","Conochilidae"), cex=0.9)
#dev.off()

#jpeg("Figures/2019-2020_NMDS_1v3_bray_days.jpg", width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_temporal_bray,display = c('sites','species'),choices = c(1,4),type = "n") 
points(NMDS_temporal_bray$points[zoop_epi_tows$groups=="1",1], NMDS_temporal_bray$points[zoop_epi_tows$groups=="1",3], pch=21,bg="#008585")
points(NMDS_temporal_bray$points[zoop_epi_tows$groups=="2",1], NMDS_temporal_bray$points[zoop_epi_tows$groups=="2",3], pch=21,bg="#9BBAA0")
points(NMDS_temporal_bray$points[zoop_epi_tows$groups=="3",1], NMDS_temporal_bray$points[zoop_epi_tows$groups=="3",3], pch=21,bg="#FBF2C4")
points(NMDS_temporal_bray$points[zoop_epi_tows$groups=="4",1], NMDS_temporal_bray$points[zoop_epi_tows$groups=="4",3], pch=21,bg="#DEA868")
points(NMDS_temporal_bray$points[zoop_epi_tows$groups=="5",1], NMDS_temporal_bray$points[zoop_epi_tows$groups=="5",3], pch=21,bg="#C7522B")
legend("bottomright", legend=c('Day1','Day2','Day3','Day4','Day5'), pch=21, pt.bg=hcl.colors(5,"Geyser"), bty = "n") 
ordihull(ord, zoop_epi_tows$groups, display = "sites", draw = c("polygon"),
         col = hcl.colors(5,"Geyser"), alpha = 75,cex = 2)
text(NMDS_temporal_bray$species[,1],NMDS_temporal_bray$species[,3], labels = c("Calanoida","Cyclopoida","Keratella","Kellicottia",
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
text(NMDS_temporal_bray$species[,1],NMDS_temporal_bray$species[,3], labels = c("Calanoida","Cyclopoida","Keratella","Kellicottia",
     "Bosmina","Daphnia","Ceriodaphnia","nauplius","Collothecidae","Synchaetidae","Conochilidae"), cex=0.9)
#dev.off()


#jpeg("Figures/2019-2020_NMDS_1v4_bray_time.jpg", width = 6, height = 5, units = "in",res = 300)
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

#jpeg("Figures/2019-2020_NMDS_1v3_bray_time.jpg", width = 6, height = 5, units = "in",res = 300)
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

#jpeg("Figures/2019-2020_NMDS_1v2_bray_timegroups.jpg", width = 6, height = 5, units = "in",res = 300)
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
#first average times for each 24-hour campaign so there are 10 points per day (basically just averaging noon and midnight)
zoop_epi_tows$order <- ifelse(zoop_epi_tows$Hour=="4:" | zoop_epi_tows$Hour=="3:",1,ifelse(zoop_epi_tows$Hour=="5:",2,ifelse(
  zoop_epi_tows$Hour=="6:",3,ifelse(zoop_epi_tows$Hour=="7:",4,ifelse(zoop_epi_tows$time=="noon",5,ifelse(
    zoop_epi_tows$Hour=="18",6,ifelse(zoop_epi_tows$Hour=="19",7,ifelse(zoop_epi_tows$Hour=="20",8,ifelse(zoop_epi_tows$Hour=="21",9,10)))))))))

#average by groups, site, then order
zoop_avg <- zoop_epi_tows %>% group_by(groups,site,order) %>%
  summarise_at(vars(Calanoida_density_NopL_1:Conochilidae_density_NopL_1), list(mean = mean))

#pelagic vs littoral dfs
zoop_pel <- zoop_avg[zoop_avg$site=="pel",]
zoop_lit <- zoop_avg[zoop_avg$site=="lit",]

#only select data cols
zoop_temporal_avg_dens <- zoop_avg[,c(grepl("mean",colnames(zoop_avg)))] 
zoop_pel_dens <- zoop_pel[,c(grepl("mean",colnames(zoop_pel)))]
zoop_lit_dens <- zoop_lit[,c(grepl("mean",colnames(zoop_lit)))]

#transforming data - hellinger transformation because gives low weight to low/zero values
#converts species abundances from absolute to relative - use w/ bray curtis
zoop_temporal_dens_avg_trans <- hellinger(zoop_temporal_avg_dens)
zoop_pel_dens_trans <- hellinger(zoop_pel_dens)
zoop_lit_dens_trans <- hellinger(zoop_lit_dens)


#scree plot to choose dimension #
dimcheckMDS(zoop_temporal_dens_avg_trans, distance = "bray", k = 6, trymax = 20, autotransform = TRUE)

#now do NMDS using averages w/ 4 dimensions for consistency
NMDS_temporal_avg_bray <- metaMDS(zoop_temporal_dens_avg_trans, distance='bray', k=4, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_temporal_avg_bray$stress

NMDS_pel_bray <- metaMDS(zoop_pel_dens_trans, distance='bray', k=4, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_pel_bray$stress

NMDS_lit_bray <- metaMDS(zoop_lit_dens_trans, distance='bray', k=4, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_lit_bray$stress

#-------------------------------------------------------------------------------#
#                                 Tracking figs                                 #
#-------------------------------------------------------------------------------#
#jpeg(file.path(getwd(),"Summer2021-DataAnalysis/Figures/2019-2020_NMDS_1v2_bray_temporal_avg_days.jpg"), width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_temporal_avg_bray,display = c('sites','species'),choices = c(1,2),type = "n") 
ordihull(ord, zoop_avg$groups, display = "sites", draw = c("polygon"),
         col = c("#008585","#89B199","#EFECBF","#DB9B5A","#C7522B"), alpha = 75,cex = 2)
lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="1" & zoop_avg$site=="pel",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="1" & zoop_avg$site=="pel",2], col="#008585")
lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="1" & zoop_avg$site=="lit",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="1" & zoop_avg$site=="lit",2], col="#369187")
points(NMDS_temporal_avg_bray$points[zoop_avg$groups=="1",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="1",2], pch=rep(c(0,1,2,3,4,5,6,7,8,9),2),col=c(rep("#369187",10),rep("#008585",10)))

lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="2" & zoop_avg$site=="pel",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="2" & zoop_avg$site=="pel",2], col="#89B199")
lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="2" & zoop_avg$site=="lit",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="2" & zoop_avg$site=="lit",2], col="#C4D5B2")
points(NMDS_temporal_avg_bray$points[zoop_avg$groups=="2",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="2",2], pch=rep(c(0,1,2,3,4,5,6,7,8,9),2),col=c(rep("#C4D5B2",10),rep("#89B199",10)))

lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="3" & zoop_avg$site=="pel",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="3" & zoop_avg$site=="pel",2], col="#EFECBF")
lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="3" & zoop_avg$site=="lit",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="3" & zoop_avg$site=="lit",2], col="#EFDBA7")
points(NMDS_temporal_avg_bray$points[zoop_avg$groups=="3",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="3",2], pch=rep(c(0,1,2,3,4,5,6,7,8,9),2),col=c(rep("#EFDBA7",10),rep("#EFECBF",10)))

lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="4" & zoop_avg$site=="pel",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="4" & zoop_avg$site=="pel",2], col="#DB9B5A")
lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="4" & zoop_avg$site=="lit",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="4" & zoop_avg$site=="lit",2], col="#E4BC80")
points(NMDS_temporal_avg_bray$points[zoop_avg$groups=="4",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="4",2], pch=rep(c(0,1,2,3,4,5,6,7,8,9),2),col=c(rep("#E4BC80",10),rep("#DB9B5A",10)))

lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="5" & zoop_avg$site=="pel",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="5" & zoop_avg$site=="pel",2], col="#C7522B")
lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="5" & zoop_avg$site=="lit",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="5" & zoop_avg$site=="lit",2], col="#CA602E")
points(NMDS_temporal_avg_bray$points[zoop_avg$groups=="5",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="5",2], pch=rep(c(0,1,2,3,4,5,6,7,8,9),2),col=c(rep("#CA602E",10),rep("#C7522B",10)))

legend("bottomright", legend=c('pelagic day1','littoral day1','pelagic day2','littoral day2','pelagic day3','littoral day3', "pelagic day4", "littoral day4", "pelagic day5", "littoral day5"), pch=21, 
       pt.bg=c("#008585","#369187","#89B199","#C4D5B2","#EFECBF","#EFDBA7","#DB9B5A","#E4BC80", "#C7522B","#CA602E"),bty = "n",cex=0.8) 
legend("bottomleft", legend=c('sunrise1','sunrise2','sunrise3','sunrise4','noon','sunset1','sunset2','sunset3','sunset4','midnight'), pch=c(0,1,2,3,4,5,6,7,8,9) ,bty = "n",cex=0.8) 
#dev.off()

#jpeg(file.path(getwd(),"Summer2021-DataAnalysis/Figures/2019-2020_NMDS_1v3_bray_temporal_avg_days.jpg"), width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_temporal_avg_bray,display = c('sites','species'),choices = c(1,3),type = "n") 
ordihull(ord, zoop_avg$groups, display = "sites", draw = c("polygon"),
         col = c("#008585","#89B199","#EFECBF","#DB9B5A","#C7522B"), alpha = 75,cex = 2)
lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="1" & zoop_avg$site=="pel",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="1" & zoop_avg$site=="pel",3], col="#008585")
lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="1" & zoop_avg$site=="lit",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="1" & zoop_avg$site=="lit",3], col="#369187")
points(NMDS_temporal_avg_bray$points[zoop_avg$groups=="1",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="1",3], pch=rep(c(0,1,2,3,4,5,6,7,8,9),2),col=c(rep("#369187",10),rep("#008585",10)))

lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="2" & zoop_avg$site=="pel",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="2" & zoop_avg$site=="pel",3], col="#89B199")
lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="2" & zoop_avg$site=="lit",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="2" & zoop_avg$site=="lit",3], col="#C4D5B2")
points(NMDS_temporal_avg_bray$points[zoop_avg$groups=="2",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="2",3], pch=rep(c(0,1,2,3,4,5,6,7,8,9),2),col=c(rep("#C4D5B2",10),rep("#89B199",10)))

lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="3" & zoop_avg$site=="pel",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="3" & zoop_avg$site=="pel",3], col="#EFECBF")
lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="3" & zoop_avg$site=="lit",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="3" & zoop_avg$site=="lit",3], col="#EFDBA7")
points(NMDS_temporal_avg_bray$points[zoop_avg$groups=="3",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="3",3], pch=rep(c(0,1,2,3,4,5,6,7,8,9),2),col=c(rep("#EFDBA7",10),rep("#EFECBF",10)))

lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="4" & zoop_avg$site=="pel",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="4" & zoop_avg$site=="pel",3], col="#DB9B5A")
lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="4" & zoop_avg$site=="lit",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="4" & zoop_avg$site=="lit",3], col="#E4BC80")
points(NMDS_temporal_avg_bray$points[zoop_avg$groups=="4",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="4",3], pch=rep(c(0,1,2,3,4,5,6,7,8,9),2),col=c(rep("#E4BC80",10),rep("#DB9B5A",10)))

lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="5" & zoop_avg$site=="pel",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="5" & zoop_avg$site=="pel",3], col="#C7522B")
lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="5" & zoop_avg$site=="lit",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="5" & zoop_avg$site=="lit",3], col="#CA602E")
points(NMDS_temporal_avg_bray$points[zoop_avg$groups=="5",1], NMDS_temporal_avg_bray$points[zoop_avg$groups=="5",3], pch=rep(c(0,1,2,3,4,5,6,7,8,9),2),col=c(rep("#CA602E",10),rep("#C7522B",10)))

legend("bottomright", legend=c('pelagic day1','littoral day1','pelagic day2','littoral day2','pelagic day3','littoral day3', "pelagic day4", "littoral day4", "pelagic day5", "littoral day5"), pch=21, 
       pt.bg=c("#008585","#369187","#89B199","#C4D5B2","#EFECBF","#EFDBA7","#DB9B5A","#E4BC80", "#C7522B","#CA602E"),bty = "n",cex=0.8) 
legend("bottomleft", legend=c('sunrise1','sunrise2','sunrise3','sunrise4','noon','sunset1','sunset2','sunset3','sunset4','midnight'), pch=c(0,1,2,3,4,5,6,7,8,9) ,bty = "n",cex=0.8) 
#dev.off()


#pelagic only tracking density through time
#jpeg(file.path(getwd(),"Summer2021-DataAnalysis/Figures/2019-2020_NMDS_1v2_bray_pelagic_tracking_density_time.jpg"), width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_pel_bray,display = c('sites','species'),choices = c(1,2),type = "n") 
ordihull(ord, zoop_pel$groups, display = "sites", draw = c("polygon"),
         col = c("#008585","#89B199","#EFECBF","#DB9B5A","#C7522B"), alpha = 75,cex = 2)
lines(NMDS_pel_bray$points[zoop_pel$groups=="1",1], NMDS_pel_bray$points[zoop_pel$groups=="1",2], col="#008585")
points(NMDS_pel_bray$points[zoop_pel$groups=="1",1], NMDS_pel_bray$points[zoop_pel$groups=="1",2], pch=c(0,1,2,3,4,5,6,7,8,9),col="#008585")

lines(NMDS_pel_bray$points[zoop_pel$groups=="2",1], NMDS_pel_bray$points[zoop_pel$groups=="2",2], col="#89B199")
points(NMDS_pel_bray$points[zoop_pel$groups=="2",1], NMDS_pel_bray$points[zoop_pel$groups=="2",2], pch=c(0,1,2,3,4,5,6,7,8,9),col="#89B199")

lines(NMDS_pel_bray$points[zoop_pel$groups=="3",1], NMDS_pel_bray$points[zoop_pel$groups=="3",2], col="#EFECBF")
points(NMDS_pel_bray$points[zoop_pel$groups=="3",1], NMDS_pel_bray$points[zoop_pel$groups=="3",2], pch=c(0,1,2,3,4,5,6,7,8,9),col="#EFECBF")

lines(NMDS_pel_bray$points[zoop_pel$groups=="4",1], NMDS_pel_bray$points[zoop_pel$groups=="4",2], col="#DB9B5A")
points(NMDS_pel_bray$points[zoop_pel$groups=="4",1], NMDS_pel_bray$points[zoop_pel$groups=="4",2], pch=c(0,1,2,3,4,5,6,7,8,9),col="#DB9B5A")

lines(NMDS_pel_bray$points[zoop_pel$groups=="5",1], NMDS_pel_bray$points[zoop_pel$groups=="5",2], col="#C7522B")
points(NMDS_pel_bray$points[zoop_pel$groups=="5",1], NMDS_pel_bray$points[zoop_pel$groups=="5",2], pch=c(0,1,2,3,4,5,6,7,8,9),col="#C7522B")

legend("bottomright", legend=c('day1', 'day2','day3','day4', 'day5'), pch=21, 
       pt.bg=c("#008585","#89B199","#EFECBF","#DB9B5A", "#C7522B"),bty = "n",cex=0.8) 
legend("bottomleft", legend=c('sunrise1','sunrise2','sunrise3','sunrise4','noon','sunset1','sunset2','sunset3','sunset4','midnight'), pch=c(0,1,2,3,4,5,6,7,8,9) ,bty = "n",cex=0.8) 
#dev.off()

#jpeg(file.path(getwd(),"Summer2021-DataAnalysis/Figures/2019-2020_NMDS_1v3_bray_pelagic_tracking_density_time.jpg"), width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_pel_bray,display = c('sites','species'),choices = c(1,3),type = "n") 
ordihull(ord, zoop_pel$groups, display = "sites", draw = c("polygon"),
         col = c("#008585","#89B199","#EFECBF","#DB9B5A","#C7522B"), alpha = 75,cex = 2)
lines(NMDS_pel_bray$points[zoop_pel$groups=="1",1], NMDS_pel_bray$points[zoop_pel$groups=="1",3], col="#008585")
points(NMDS_pel_bray$points[zoop_pel$groups=="1",1], NMDS_pel_bray$points[zoop_pel$groups=="1",3], pch=c(0,1,2,3,4,5,6,7,8,9),col="#008585")

lines(NMDS_pel_bray$points[zoop_pel$groups=="2",1], NMDS_pel_bray$points[zoop_pel$groups=="2",3], col="#89B199")
points(NMDS_pel_bray$points[zoop_pel$groups=="2",1], NMDS_pel_bray$points[zoop_pel$groups=="2",3], pch=c(0,1,2,3,4,5,6,7,8,9),col="#89B199")

lines(NMDS_pel_bray$points[zoop_pel$groups=="3",1], NMDS_pel_bray$points[zoop_pel$groups=="3",3], col="#EFECBF")
points(NMDS_pel_bray$points[zoop_pel$groups=="3",1], NMDS_pel_bray$points[zoop_pel$groups=="3",3], pch=c(0,1,2,3,4,5,6,7,8,9),col="#EFECBF")

lines(NMDS_pel_bray$points[zoop_pel$groups=="4",1], NMDS_pel_bray$points[zoop_pel$groups=="4",3], col="#DB9B5A")
points(NMDS_pel_bray$points[zoop_pel$groups=="4",1], NMDS_pel_bray$points[zoop_pel$groups=="4",3], pch=c(0,1,2,3,4,5,6,7,8,9),col="#DB9B5A")

lines(NMDS_pel_bray$points[zoop_pel$groups=="5",1], NMDS_pel_bray$points[zoop_pel$groups=="5",3], col="#C7522B")
points(NMDS_pel_bray$points[zoop_pel$groups=="5",1], NMDS_pel_bray$points[zoop_pel$groups=="5",3], pch=c(0,1,2,3,4,5,6,7,8,9),col="#C7522B")

legend("bottomright", legend=c('day1', 'day2','day3','day4', 'day5'), pch=21, 
       pt.bg=c("#008585","#89B199","#EFECBF","#DB9B5A", "#C7522B"),bty = "n",cex=0.8) 
legend("bottomleft", legend=c('sunrise1','sunrise2','sunrise3','sunrise4','noon','sunset1','sunset2','sunset3','sunset4','midnight'), pch=c(0,1,2,3,4,5,6,7,8,9) ,bty = "n",cex=0.8) 
#dev.off()

#littoral only tracking density through time
#jpeg(file.path(getwd(),"Summer2021-DataAnalysis/Figures/2019-2020_NMDS_1v2_bray_littoral_tracking_density_time.jpg"), width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_lit_bray,display = c('sites','species'),choices = c(1,2),type = "n") 
ordihull(ord, zoop_lit$groups, display = "sites", draw = c("polygon"),
         col = c("#008585","#89B199","#EFECBF","#DB9B5A","#C7522B"), alpha = 75,cex = 2)
lines(NMDS_lit_bray$points[zoop_lit$groups=="1",1], NMDS_lit_bray$points[zoop_lit$groups=="1",2], col="#369187")
points(NMDS_lit_bray$points[zoop_lit$groups=="1",1], NMDS_lit_bray$points[zoop_lit$groups=="1",2], pch=c(0,1,2,3,4,5,6,7,8,9),col="#369187")

lines(NMDS_lit_bray$points[zoop_lit$groups=="2",1], NMDS_lit_bray$points[zoop_lit$groups=="2",2], col="#C4D5B2")
points(NMDS_lit_bray$points[zoop_lit$groups=="2",1], NMDS_lit_bray$points[zoop_lit$groups=="2",2], pch=c(0,1,2,3,4,5,6,7,8,9),col="#C4D5B2")

lines(NMDS_lit_bray$points[zoop_lit$groups=="3",1], NMDS_lit_bray$points[zoop_lit$groups=="3",2], col="#EFDBA7")
points(NMDS_lit_bray$points[zoop_lit$groups=="3",1], NMDS_lit_bray$points[zoop_lit$groups=="3",2], pch=c(0,1,2,3,4,5,6,7,8,9),col=c(rep("#EFDBA7",10),rep("#EFECBF",10)))

lines(NMDS_lit_bray$points[zoop_lit$groups=="4",1], NMDS_lit_bray$points[zoop_lit$groups=="4",2], col="#E4BC80")
points(NMDS_lit_bray$points[zoop_lit$groups=="4",1], NMDS_lit_bray$points[zoop_lit$groups=="4",2], pch=c(0,1,2,3,4,5,6,7,8,9),col="#E4BC80")

lines(NMDS_lit_bray$points[zoop_lit$groups=="5",1], NMDS_lit_bray$points[zoop_lit$groups=="5",2], col="#CA602E")
points(NMDS_lit_bray$points[zoop_lit$groups=="5",1], NMDS_lit_bray$points[zoop_lit$groups=="5",2], pch=c(0,1,2,3,4,5,6,7,8,9),col="#CA602E")

legend("bottomright", legend=c('day1', 'day2','day3','day4', 'day5'), pch=21, 
       pt.bg=c("#369187","#C4D5B2","#EFDBA7","#E4BC80","#CA602E"),bty = "n",cex=0.8) 
legend("bottomleft", legend=c('sunrise1','sunrise2','sunrise3','sunrise4','noon','sunset1','sunset2','sunset3','sunset4','midnight'), pch=c(0,1,2,3,4,5,6,7,8,9) ,bty = "n",cex=0.8) 
#dev.off()

#jpeg(file.path(getwd(),"Summer2021-DataAnalysis/Figures/2019-2020_NMDS_1v3_bray_littoral_tracking_density_time.jpg"), width = 6, height = 5, units = "in",res = 300)
ord <- ordiplot(NMDS_lit_bray,display = c('sites','species'),choices = c(1,3),type = "n") 
ordihull(ord, zoop_lit$groups, display = "sites", draw = c("polygon"),
         col = c("#008585","#89B199","#EFECBF","#DB9B5A","#C7522B"), alpha = 75,cex = 2)
lines(NMDS_lit_bray$points[zoop_lit$groups=="1",1], NMDS_lit_bray$points[zoop_lit$groups=="1",3], col="#369187")
points(NMDS_lit_bray$points[zoop_lit$groups=="1",1], NMDS_lit_bray$points[zoop_lit$groups=="1",3], pch=c(0,1,2,3,4,5,6,7,8,9),col="#369187")

lines(NMDS_lit_bray$points[zoop_lit$groups=="2",1], NMDS_lit_bray$points[zoop_lit$groups=="2",3], col="#C4D5B2")
points(NMDS_lit_bray$points[zoop_lit$groups=="2",1], NMDS_lit_bray$points[zoop_lit$groups=="2",3], pch=c(0,1,2,3,4,5,6,7,8,9),col="#C4D5B2")

lines(NMDS_lit_bray$points[zoop_lit$groups=="3",1], NMDS_lit_bray$points[zoop_lit$groups=="3",3], col="#EFDBA7")
points(NMDS_lit_bray$points[zoop_lit$groups=="3",1], NMDS_lit_bray$points[zoop_lit$groups=="3",3], pch=c(0,1,2,3,4,5,6,7,8,9),col="#EFDBA7")

lines(NMDS_lit_bray$points[zoop_lit$groups=="4",1], NMDS_lit_bray$points[zoop_lit$groups=="4",3], col="#E4BC80")
points(NMDS_lit_bray$points[zoop_lit$groups=="4",1], NMDS_lit_bray$points[zoop_lit$groups=="4",3], pch=c(0,1,2,3,4,5,6,7,8,9),col="#E4BC80")

lines(NMDS_lit_bray$points[zoop_lit$groups=="5",1], NMDS_lit_bray$points[zoop_lit$groups=="5",3], col="#CA602E")
points(NMDS_lit_bray$points[zoop_lit$groups=="5",1], NMDS_lit_bray$points[zoop_lit$groups=="5",3], pch=c(0,1,2,3,4,5,6,7,8,9),col="#CA602E")

legend("bottomright", legend=c('day1', 'day2','day3','day4', 'day5'), pch=21, 
       pt.bg=c("#369187","#C4D5B2","#EFDBA7","#E4BC80","#CA602E"),bty = "n",cex=0.8) 
legend("bottomleft", legend=c('sunrise1','sunrise2','sunrise3','sunrise4','noon','sunset1','sunset2','sunset3','sunset4','midnight'), pch=c(0,1,2,3,4,5,6,7,8,9) ,bty = "n",cex=0.8) 
#dev.off()

#-------------------------------------------------------------------------------#
#                     Calculating euclidean distance                            #
#-------------------------------------------------------------------------------#
#so I'm technically calculating the euclidean distances from bray-curtis distance matrices

#step 1: take NMDS output for each site using NMDS coordinates
zoop_pel_euc <- as.matrix(vegdist(NMDS_pel_bray$points, method='euclidean')) 
zoop_lit_euc <- as.matrix(vegdist(NMDS_lit_bray$points, method='euclidean'))

#step 2: select and sum the 9 distances between connecting points for each of the 5 days
pel_day1 <- sum(zoop_pel_euc[1,2],zoop_pel_euc[2,3],zoop_pel_euc[3,4],zoop_pel_euc[4,5],zoop_pel_euc[5,6],
                zoop_pel_euc[6,7],zoop_pel_euc[7,8],zoop_pel_euc[8,9],zoop_pel_euc[9,10])

pel_day2 <- sum(zoop_pel_euc[11,12], zoop_pel_euc[12,13],zoop_pel_euc[13,14],zoop_pel_euc[14,15],zoop_pel_euc[15,16],
                zoop_pel_euc[16,17],zoop_pel_euc[17,18],zoop_pel_euc[18,19],zoop_pel_euc[19,20])

pel_day3 <- sum(zoop_pel_euc[21,22], zoop_pel_euc[22,23],zoop_pel_euc[23,24],zoop_pel_euc[24,25],zoop_pel_euc[25,26],
                zoop_pel_euc[26,27],zoop_pel_euc[27,28],zoop_pel_euc[28,29],zoop_pel_euc[29,30])

pel_day4 <- sum(zoop_pel_euc[31,32], zoop_pel_euc[32,33],zoop_pel_euc[33,34],zoop_pel_euc[34,35],zoop_pel_euc[35,36],
                zoop_pel_euc[36,37],zoop_pel_euc[37,38],zoop_pel_euc[38,39],zoop_pel_euc[39,40])

pel_day5 <- sum(zoop_pel_euc[41,42], zoop_pel_euc[42,43],zoop_pel_euc[43,44],zoop_pel_euc[44,45],zoop_pel_euc[45,46],
                zoop_pel_euc[46,47],zoop_pel_euc[47,48],zoop_pel_euc[48,49],zoop_pel_euc[49,50])

lit_day1 <- sum(zoop_lit_euc[1,2],zoop_lit_euc[2,3],zoop_lit_euc[3,4],zoop_lit_euc[4,5],zoop_lit_euc[5,6],
                zoop_lit_euc[6,7],zoop_lit_euc[7,8],zoop_lit_euc[8,9],zoop_lit_euc[9,10])

lit_day2 <- sum(zoop_lit_euc[11,12],zoop_lit_euc[12,13],zoop_lit_euc[13,14],zoop_lit_euc[14,15],zoop_lit_euc[15,16],
                zoop_lit_euc[16,17],zoop_lit_euc[17,18],zoop_lit_euc[18,19],zoop_lit_euc[19,20])

lit_day3 <- sum(zoop_lit_euc[21,22],zoop_lit_euc[22,23],zoop_lit_euc[23,24],zoop_lit_euc[24,25],zoop_lit_euc[25,26],
                zoop_lit_euc[26,27],zoop_lit_euc[27,28],zoop_lit_euc[28,29],zoop_lit_euc[29,30])

lit_day4 <- sum(zoop_lit_euc[31,32],zoop_lit_euc[32,33],zoop_lit_euc[33,34],zoop_lit_euc[34,35],zoop_lit_euc[35,36],
                zoop_lit_euc[36,37],zoop_lit_euc[37,38],zoop_lit_euc[38,39],zoop_lit_euc[39,40])

lit_day5 <- sum(zoop_lit_euc[41,42],zoop_lit_euc[42,43],zoop_lit_euc[43,44],zoop_lit_euc[44,45],zoop_lit_euc[45,46],
                zoop_lit_euc[46,47],zoop_lit_euc[47,48],zoop_lit_euc[48,49],zoop_lit_euc[49,50])

#step 3: make a dataset of data
euc_distances_df <- data.frame(pelagic=c(pel_day1,pel_day2,pel_day3,pel_day4,pel_day5),
                               littoral=c(lit_day1,lit_day2,lit_day3,lit_day4,lit_day5))

#plot littoral vs pelagic euclidean distances
#jpeg(file.path(getwd(),"Summer2021-DataAnalysis/Figures/2019-2020_pelagic_vs_littoral_euclidean_dist_daily_sums.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euc_distances_df$littoral,euc_distances_df$pelagic, xlab="littoral", ylab="pelagic", 
     main="Daily euclidean distance sums", cex=2.8, pch=21, cex.lab = 1.5)
points(euc_distances_df$littoral[1],euc_distances_df$pelagic[1], bg="#008585",pch=21,cex=3)
points(euc_distances_df$littoral[2],euc_distances_df$pelagic[2], bg="#9BBAA0" ,pch=21,cex=3)
points(euc_distances_df$littoral[3],euc_distances_df$pelagic[3], bg="#FBF2C4",pch=21,cex=3)
points(euc_distances_df$littoral[4],euc_distances_df$pelagic[4], bg="#DEA868",pch=21,cex=3)
points(euc_distances_df$littoral[5],euc_distances_df$pelagic[5], bg="#C7522B",pch=21,cex=3)
legend("bottomleft",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#FBF2C4","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#-----------------------------------------------------------------------------------------#
# read in all driver data (YSI, CTD, flora, chem?)
#day1 (10-11 Jul 2019), day2 (24-25 Jul 2019), day3 (12-13 Aug 2020), day4 (15-16 Jun 2021), day5 (7-8 Jul 2021)

#df to compile euclidean distances for each day/site
euclidean_drivers_df <- data.frame(site = c("day1","day2","day3","day4","day5"),
                                   pel_euc_dist = euc_distances_df$pelagic,
                                   lit_euc_dist = euc_distances_df$littoral)


#EDI YSI data
inUrl2  <- "https://pasta.lternet.edu/package/data/eml/edi/198/10/b3bd353312f9e37ca392e2a5315cc9da" 
infile2 <- tempfile()
try(download.file(inUrl2,infile2,method="curl"))
if (is.na(file.size(infile2))) download.file(inUrl2,infile2,method="auto")

ysi <-read.csv(infile2,header=F 
               ,skip=1
               ,sep=","  
               , col.names=c(
                 "Reservoir",     
                 "Site",     
                 "DateTime",     
                 "Depth_m",     
                 "Temp_C",     
                 "DO_mgL",     
                 "DOSat",     
                 "Cond_uScm",     
                 "Sp_cond_uScm",     
                 "PAR_umolm2s",     
                 "ORP_mV",     
                 "pH",     
                 "Flag_DateTime",     
                 "Flag_Temp",     
                 "Flag_DO",     
                 "Flag_DOSat",     
                 "Flag_Cond",     
                 "Flag_Sp_Cond",     
                 "Flag_PAR",     
                 "Flag_ORP",     
                 "Flag_pH"    ), check.names=TRUE)

unlink(infile2)

#some of the sites are wrong for 7-8 Jul 2021
ysi$Site[ysi$DateTime=="2021-07-07 11:37:00"] <- 50
ysi$Site[ysi$DateTime=="2021-07-08 11:45:00"] <- 50

#change date format
ysi$DateTime <- as.Date(ysi$DateTime)

#separate ysi by day
#msn1_ysi <-  ysi[(ysi$DateTime=="2019-07-10" | ysi$DateTime=="2019-07-11")  & ysi$Reservoir=="BVR",] #no ysi data
#msn2_ysi <-  ysi[(ysi$DateTime=="2019-07-24" | ysi$DateTime=="2019-07-25")  & ysi$Reservoir=="BVR",] #no ysi data
msn3_ysi <-  ysi[(ysi$DateTime=="2020-08-12" | ysi$DateTime=="2020-08-13")  & ysi$Reservoir=="BVR",]
msn4_ysi <-  ysi[(ysi$DateTime=="2021-06-15" | ysi$DateTime=="2020-06-16")  & ysi$Reservoir=="BVR",]
msn5_ysi <-  ysi[(ysi$DateTime=="2021-07-07" | ysi$DateTime=="2021-07-08")  & ysi$Reservoir=="BVR" & ysi$Site==50,]

msn5_ysi_lit <- ysi[(ysi$DateTime=="2021-07-07" | ysi$DateTime=="2021-07-08")  & ysi$Reservoir=="BVR" & ysi$Site==51,]

#read in CTD EDI file
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/200/12/0a62d1946e8d9a511bc1404e69e59b8c" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")

ctd <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               , col.names=c(
                 "Reservoir",     
                 "Site",     
                 "Date",     
                 "Depth_m",     
                 "Temp_C",     
                 "DO_mgL",     
                 "DO_pSat",     
                 "Cond_uScm",     
                 "Spec_Cond_uScm",     
                 "Chla_ugL",     
                 "Turb_NTU",     
                 "pH",     
                 "ORP_mV",     
                 "PAR_umolm2s",     
                 "Desc_rate",     
                 "Flag_Temp",     
                 "Flag_DO",     
                 "Flag_Cond",     
                 "Flag_SpecCond",     
                 "Flag_Chla",     
                 "Flag_Turb",     
                 "Flag_pH",     
                 "Flag_ORP",     
                 "Flag_PAR",     
                 "Flag_DescRate"    ), check.names=TRUE)

unlink(infile1)

#change date format
ctd$Date <- as.Date(ctd$Date)

#read in DO and temp from sampling days, select 0.1 and 9m, and then average both day DO
msn1_ctd <- ctd %>% filter((Date=="2019-07-10" | Date=="2019-07-11") & Reservoir=="BVR") %>%
  group_by(Date) %>% filter(abs(Depth_m - 0.1) == min(abs(Depth_m - 0.1)) | (abs(Depth_m - 9) == min(abs(Depth_m - 9)))) %>%
  distinct(Depth_m, .keep_all = TRUE)

msn2_ctd <- ctd %>% filter((Date=="2019-07-24" | Date=="2019-07-25") & Reservoir=="BVR") %>%
  group_by(Date) %>% filter(abs(Depth_m - 0.1) == min(abs(Depth_m - 0.1)) | (abs(Depth_m - 9) == min(abs(Depth_m - 9)))) %>%
  distinct(Depth_m, .keep_all = TRUE)

msn3_ctd <- ctd %>% filter(Date=="2020-08-12" | Date=="2020-08-13" & Reservoir=="BVR") %>%
  group_by(Date) %>% filter(abs(Depth_m - 0.1) == min(abs(Depth_m - 0.1)) | (abs(Depth_m - 9) == min(abs(Depth_m - 9)))) %>%
  distinct(Depth_m, .keep_all = TRUE)

msn4_ctd <- ctd %>% filter(Date=="2021-06-15" | Date=="2021-06-16" & Reservoir=="BVR") %>%
  group_by(Date) %>% filter(abs(Depth_m - 0.1) == min(abs(Depth_m - 0.1)) | (abs(Depth_m - 9) == min(abs(Depth_m - 9)))) %>%
  distinct(Depth_m, .keep_all = TRUE)

before_msn5_ctd <- ctd[(ctd$Date=="2021-06-28")  & ctd$Reservoir=="BVR",] #no ctd on day 5
before_msn5_ctd <-  before_msn5_ctd[before_msn5_ctd$Depth_m==closest(before_msn5_ctd$Depth_m, 0.1) | before_msn5_ctd$Depth_m==closest(before_msn5_ctd$Depth_m, 9),] 

after_msn5_ctd <- ctd[(ctd$Date=="2021-07-12")  & ctd$Reservoir=="BVR",] #no ctd on day 5
after_msn5_ctd <-  after_msn5_ctd[after_msn5_ctd$Depth_m==closest(after_msn5_ctd$Depth_m, 0.1) | after_msn5_ctd$Depth_m==closest(after_msn5_ctd$Depth_m, 9),] 

#so for now, I'm just going to take the do from 4 days after because I'm assuming that DO doesn't change that much in 4 days...
msn5_0.1_ctd <- after_msn5_ctd$DO_mgL[after_msn5_ctd$Depth_m<2]
msn5_9_ctd <- after_msn5_ctd$DO_mgL[after_msn5_ctd$Depth_m>2]

#select every 0.5m from casts
ctd_final <- ctd %>%
  dplyr::mutate(rdepth = plyr::round_any(Depth_m, 0.5)) %>% 
  dplyr::group_by(Date, rdepth, Reservoir, Site) %>%
  dplyr::summarise(value = mean(Temp_C)) %>% 
  dplyr::rename(depth = rdepth) 

#calcualte thermocline depth now
msn1_therm_depth <- ctd_final %>% filter((Date=="2019-07-10" | Date=="2019-07-11") & Reservoir=="BVR") %>%
  group_by(Date) %>% filter(depth > 0) %>% mutate(therm_depth = thermo.depth(value,depth))

msn2_therm_depth <- ctd_final %>% filter((Date=="2019-07-24"  | Date=="2019-07-25") & Reservoir=="BVR") %>%
  group_by(Date) %>% filter(depth > 0) %>% mutate(therm_depth = thermo.depth(value,depth))

msn3_therm_depth <- ctd_final %>% filter((Date=="2020-08-12"  | Date=="2020-08-13") & Reservoir=="BVR") %>%
  group_by(Date) %>% filter(depth > 0) %>% mutate(therm_depth = thermo.depth(value,depth))

msn4_therm_depth <- ctd_final %>% filter((Date=="2021-06-15"  | Date=="2021-06-16") & Reservoir=="BVR") %>%
  group_by(Date) %>% filter(depth > 0) %>% mutate(therm_depth = thermo.depth(value,depth))

msn5_therm_depth <- ctd_final %>% filter((Date=="2021-07-12"  | Date=="2021-07-12") & Reservoir=="BVR") %>%
  group_by(Date) %>% filter(depth > 0) %>% mutate(therm_depth = thermo.depth(value,depth)) #4 days after msn


#read in fp data from edi
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/272/6/6b3151c0fdd913e02641363c2b00ae57" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")


fp <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               , col.names=c(
                 "Reservoir",     
                 "Site",     
                 "DateTime",     
                 "Depth_m",     
                 "GreenAlgae_ugL",     
                 "Bluegreens_ugL",     
                 "BrownAlgae_ugL",     
                 "MixedAlgae_ugL",     
                 "YellowSubstances_ugL",     
                 "TotalConc_ugL",     
                 "Transmission",     
                 "Temp_degC",     
                 "RFU_525nm",     
                 "RFU_570nm",     
                 "RFU_610nm",     
                 "RFU_370nm",     
                 "RFU_590nm",     
                 "RFU_470nm",     
                 "Flag_GreenAlgae",     
                 "Flag_BluegreenAlgae",     
                 "Flag_BrownAlgae",     
                 "Flag_MixedAlgae",     
                 "Flag_TotalConc",     
                 "Flag_Temp",     
                 "Flag_Transmission",     
                 "Flag_525nm",     
                 "Flag_570nm",     
                 "Flag_610nm",     
                 "Flag_370nm",     
                 "Flag_590nm",     
                 "Flag_470nm"    ), check.names=TRUE)

unlink(infile1)

#change date format for fp data
fp$DateTime <- as.Date(fp$DateTime)

#fp for all msns
msn1_fp <- fp %>% filter((DateTime=="2019-07-10" | DateTime=="2019-07-11") & Reservoir=="BVR") %>%
  group_by(DateTime) %>% filter(abs(Depth_m - 0.1) == min(abs(Depth_m - 0.1)) | (abs(Depth_m - 9) == min(abs(Depth_m - 9)))) %>%
  distinct(Depth_m, .keep_all = TRUE)

msn2_fp <- fp %>% filter((DateTime=="2019-07-24" | DateTime=="2019-07-25") & Reservoir=="BVR") %>%
  group_by(DateTime) %>% filter(abs(Depth_m - 0.1) == min(abs(Depth_m - 0.1)) | (abs(Depth_m - 9) == min(abs(Depth_m - 9)))) %>%
  distinct(Depth_m, .keep_all = TRUE) #nope

msn3_fp <- fp %>% filter((DateTime=="2020-08-12" | DateTime=="2020-08-13") & Reservoir=="BVR") %>%
  group_by(DateTime) %>% filter(abs(Depth_m - 0.1) == min(abs(Depth_m - 0.1)) | (abs(Depth_m - 9) == min(abs(Depth_m - 9)))) %>%
  distinct(Depth_m, .keep_all = TRUE)

msn4_fp <- fp %>% filter((DateTime=="2021-06-15" | DateTime=="2021-06-16") & Reservoir=="BVR") %>%
  group_by(DateTime) %>% filter(abs(Depth_m - 0.1) == min(abs(Depth_m - 0.1)) | (abs(Depth_m - 9) == min(abs(Depth_m - 9)))) %>%
  distinct(Depth_m, .keep_all = TRUE)

msn5_fp <- fp %>% filter((DateTime=="2021-07-07" | DateTime=="2021-07-08") & Reservoir=="BVR") %>%
  group_by(DateTime) %>% filter(abs(Depth_m - 0.1) == min(abs(Depth_m - 0.1)) | (abs(Depth_m - 9) == min(abs(Depth_m - 9)))) %>%
  distinct(Depth_m, .keep_all = TRUE) #nope

#read in nutreint data from edi
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/199/10/aa2ccc23688fc908f9d61cb217210a3d" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")


chem <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               , col.names=c(
                 "Reservoir",     
                 "Site",     
                 "DateTime",     
                 "Depth_m",     
                 "Rep",     
                 "TN_ugL",     
                 "TP_ugL",     
                 "NH4_ugL",     
                 "NO3NO2_ugL",     
                 "SRP_ugL",     
                 "DOC_mgL",     
                 "DIC_mgL",     
                 "DC_mgL",     
                 "DN_mgL",     
                 "Flag_DateTime",     
                 "Flag_TN",     
                 "Flag_TP",     
                 "Flag_NH4",     
                 "Flag_NO3NO2",     
                 "Flag_SRP",     
                 "Flag_DOC",     
                 "Flag_DIC",     
                 "Flag_DC",     
                 "Flag_DN"    ), check.names=TRUE)

unlink(infile1)

#change date format for fp data
chem$DateTime <- as.Date(chem$DateTime)

#chem for all msns
msn1_chem <- chem %>% filter((DateTime=="2019-07-10" | DateTime=="2019-07-11") & Reservoir=="BVR") %>%
  group_by(DateTime) %>% filter(Depth_m== 0.1 | Depth_m== 9)

msn2_chem <- chem %>% filter((DateTime=="2019-07-24" | DateTime=="2019-07-25") & Reservoir=="BVR") %>%
  group_by(DateTime) %>% filter(Depth_m== 0.1 | Depth_m== 9)

msn3_chem <- chem %>% filter((DateTime=="2020-08-12" | DateTime=="2020-08-13") & Reservoir=="BVR") %>%
  group_by(DateTime) %>% filter(Depth_m== 0.1 | Depth_m== 9)

msn4_chem <- chem %>% filter((DateTime=="2021-06-15" | DateTime=="2021-06-16") & Reservoir=="BVR" & Site==50) %>%
  group_by(DateTime) %>% filter(Depth_m== 0.1 | Depth_m== 9)

msn5_chem <- chem %>% filter((DateTime=="2021-07-11" | DateTime=="2021-07-12") & Reservoir=="BVR" & Site==50) %>%
  group_by(DateTime) %>% filter(Depth_m== 0.1 | Depth_m== 9)

#add environmental drivers to euclidean distance df
euclidean_drivers_df$DO_0.1m <- c(mean(msn1_ctd$DO_mgL[msn1_ctd$Depth_m<2]), mean(msn2_ctd$DO_mgL[msn2_ctd$Depth_m<2]),
                                  mean(msn3_ctd$DO_mgL[msn3_ctd$Depth_m<2]),mean(msn4_ctd$DO_mgL[msn4_ctd$Depth_m<2]),msn5_0.1_ctd)

euclidean_drivers_df$DO_9.0m <- c(mean(msn1_ctd$DO_mgL[msn1_ctd$Depth_m>2]), mean(msn2_ctd$DO_mgL[msn2_ctd$Depth_m>2]),
                                  mean(msn3_ctd$DO_mgL[msn3_ctd$Depth_m>2]),mean(msn4_ctd$DO_mgL[msn4_ctd$Depth_m>2]),msn5_9_ctd)

euclidean_drivers_df$DOsat_0.1m <- c(mean(msn1_ctd$DO_pSat[msn1_ctd$Depth_m<2]), mean(msn2_ctd$DO_pSat[msn2_ctd$Depth_m<2]),
                                     mean(msn3_ctd$DO_pSat[msn3_ctd$Depth_m<2]),mean(msn4_ctd$DO_pSat[msn4_ctd$Depth_m<2]),
                                     after_msn5_ctd$DO_pSat[after_msn5_ctd$Depth_m<2])

euclidean_drivers_df$DOsat_9.0m <- c(mean(msn1_ctd$DO_pSat[msn1_ctd$Depth_m>2]), mean(msn2_ctd$DO_pSat[msn2_ctd$Depth_m>2]),
                                     mean(msn3_ctd$DO_pSat[msn3_ctd$Depth_m>2]),mean(msn4_ctd$DO_pSat[msn4_ctd$Depth_m>2]),
                                     after_msn5_ctd$DO_pSat[after_msn5_ctd$Depth_m>2])

euclidean_drivers_df$Spcond_0.1m <- c(mean(msn1_ctd$Spec_Cond_uScm[msn1_ctd$Depth_m<2]), mean(msn2_ctd$Spec_Cond_uScm[msn2_ctd$Depth_m<2]),
                                      mean(msn3_ctd$Spec_Cond_uScm[msn3_ctd$Depth_m<2]),mean(msn4_ctd$Spec_Cond_uScm[msn4_ctd$Depth_m<2]),
                                      after_msn5_ctd$Spec_Cond_uScm[after_msn5_ctd$Depth_m<2])

euclidean_drivers_df$Spcond_9.0m <- c(mean(msn1_ctd$Spec_Cond_uScm[msn1_ctd$Depth_m>2]), mean(msn2_ctd$Spec_Cond_uScm[msn2_ctd$Depth_m>2]),
                                      mean(msn3_ctd$Spec_Cond_uScm[msn3_ctd$Depth_m>2]),mean(msn4_ctd$Spec_Cond_uScm[msn4_ctd$Depth_m>2]),
                                      after_msn5_ctd$Spec_Cond_uScm[after_msn5_ctd$Depth_m>2])

euclidean_drivers_df$chla_0.1m <- c(mean(msn1_ctd$Chla_ugL[msn1_ctd$Depth_m<2]), mean(msn2_ctd$Chla_ugL[msn2_ctd$Depth_m<2]),
                                    mean(msn3_ctd$Chla_ugL[msn3_ctd$Depth_m<2]),mean(msn4_ctd$Chla_ugL[msn4_ctd$Depth_m<2]),
                                    after_msn5_ctd$Chla_ugL[after_msn5_ctd$Depth_m<2])

euclidean_drivers_df$chla_9.0m <- c(mean(msn1_ctd$Chla_ugL[msn1_ctd$Depth_m>2]), mean(msn2_ctd$Chla_ugL[msn2_ctd$Depth_m>2]),
                                    mean(msn3_ctd$Chla_ugL[msn3_ctd$Depth_m>2]),mean(msn4_ctd$Chla_ugL[msn4_ctd$Depth_m>2]),
                                    after_msn5_ctd$Chla_ugL[after_msn5_ctd$Depth_m>2])

euclidean_drivers_df$turb_0.1m <- c(mean(msn1_ctd$Turb_NTU[msn1_ctd$Depth_m<2]), mean(msn2_ctd$Turb_NTU[msn2_ctd$Depth_m<2]),
                                    mean(msn3_ctd$Turb_NTU[msn3_ctd$Depth_m<2]),mean(msn4_ctd$Turb_NTU[msn4_ctd$Depth_m<2]),
                                    after_msn5_ctd$Turb_NTU[after_msn5_ctd$Depth_m<2])

euclidean_drivers_df$turb_9.0m <- c(mean(msn1_ctd$Turb_NTU[msn1_ctd$Depth_m>2]), mean(msn2_ctd$Turb_NTU[msn2_ctd$Depth_m>2]),
                                    mean(msn3_ctd$Turb_NTU[msn3_ctd$Depth_m>2]),mean(msn4_ctd$Turb_NTU[msn4_ctd$Depth_m>2]),
                                    after_msn5_ctd$Turb_NTU[after_msn5_ctd$Depth_m>2])

euclidean_drivers_df$temp_0.1m <- c(mean(msn1_ctd$Temp_C[msn1_ctd$Depth_m<2]), mean(msn2_ctd$Temp_C[msn2_ctd$Depth_m<2]),
                                    mean(msn3_ctd$Temp_C[msn3_ctd$Depth_m<2]),mean(msn4_ctd$Temp_C[msn4_ctd$Depth_m<2]),
                                    after_msn5_ctd$Temp_C[after_msn5_ctd$Depth_m<2])

euclidean_drivers_df$temp_9.0m <- c(mean(msn1_ctd$Temp_C[msn1_ctd$Depth_m>2]), mean(msn2_ctd$Temp_C[msn2_ctd$Depth_m>2]),
                                    mean(msn3_ctd$Temp_C[msn3_ctd$Depth_m>2]),mean(msn4_ctd$Temp_C[msn4_ctd$Depth_m>2]),
                                    after_msn5_ctd$Temp_C[after_msn5_ctd$Depth_m>2])

euclidean_drivers_df$par_0.1m <- c(mean(msn1_ctd$PAR_umolm2s[msn1_ctd$Depth_m<2]), mean(msn2_ctd$PAR_umolm2s[msn2_ctd$Depth_m<2]),
                                    mean(msn3_ctd$PAR_umolm2s[msn3_ctd$Depth_m<2]),mean(msn4_ctd$PAR_umolm2s[msn4_ctd$Depth_m<2]),
                                    after_msn5_ctd$PAR_umolm2s[after_msn5_ctd$Depth_m<2])

euclidean_drivers_df$par_9.0m <- c(mean(msn1_ctd$PAR_umolm2s[msn1_ctd$Depth_m>2]), mean(msn2_ctd$PAR_umolm2s[msn2_ctd$Depth_m>2]),
                                   mean(msn3_ctd$PAR_umolm2s[msn3_ctd$Depth_m>2]),mean(msn4_ctd$PAR_umolm2s[msn4_ctd$Depth_m>2]),
                                   after_msn5_ctd$PAR_umolm2s[after_msn5_ctd$Depth_m>2])


euclidean_drivers_df$thermo_depth <- c(mean(c(first(msn1_therm_depth$therm_depth),last(msn1_therm_depth$therm_depth))),
                                       mean(c(first(msn2_therm_depth$therm_depth),last(msn2_therm_depth$therm_depth))),
                                       mean(c(first(msn3_therm_depth$therm_depth),last(msn3_therm_depth$therm_depth))),
                                       mean(c(first(msn4_therm_depth$therm_depth),last(msn4_therm_depth$therm_depth))),
                                       mean(c(first(msn5_therm_depth$therm_depth),last(msn5_therm_depth$therm_depth))))

euclidean_drivers_df$TN_surf <- c(msn1_chem$TN_ugL[msn1_chem$Depth_m==0.1], 
                                  msn2_chem$TN_ugL[msn2_chem$Depth_m==0.1],
                                  msn3_chem$TN_ugL[msn3_chem$Depth_m==0.1],
                                  msn4_chem$TN_ugL[msn4_chem$Depth_m==0.1],
                                  msn5_chem$TN_ugL[msn5_chem$Depth_m==0.1]) 

euclidean_drivers_df$TP_surf <- c(msn1_chem$TP_ugL[msn1_chem$Depth_m==0.1], 
                                  msn2_chem$TP_ugL[msn2_chem$Depth_m==0.1],
                                  msn3_chem$TP_ugL[msn3_chem$Depth_m==0.1],
                                  msn4_chem$TP_ugL[msn4_chem$Depth_m==0.1],
                                  msn5_chem$TP_ugL[msn5_chem$Depth_m==0.1]) 

euclidean_drivers_df$NH4_surf <- c(msn1_chem$NH4_ugL[msn1_chem$Depth_m==0.1], 
                                  msn2_chem$NH4_ugL[msn2_chem$Depth_m==0.1],
                                  msn3_chem$NH4_ugL[msn3_chem$Depth_m==0.1],
                                  msn4_chem$NH4_ugL[msn4_chem$Depth_m==0.1],
                                  msn5_chem$NH4_ugL[msn5_chem$Depth_m==0.1]) 

euclidean_drivers_df$NO3NO2_surf <- c(msn1_chem$NO3NO2_ugL[msn1_chem$Depth_m==0.1], 
                                  msn2_chem$NO3NO2_ugL[msn2_chem$Depth_m==0.1],
                                  msn3_chem$NO3NO2_ugL[msn3_chem$Depth_m==0.1],
                                  msn4_chem$NO3NO2_ugL[msn4_chem$Depth_m==0.1],
                                  msn5_chem$NO3NO2_ugL[msn5_chem$Depth_m==0.1]) 

euclidean_drivers_df$SRP_surf <- c(msn1_chem$SRP_ugL[msn1_chem$Depth_m==0.1], 
                                  msn2_chem$SRP_ugL[msn2_chem$Depth_m==0.1],
                                  msn3_chem$SRP_ugL[msn3_chem$Depth_m==0.1],
                                  msn4_chem$SRP_ugL[msn4_chem$Depth_m==0.1],
                                  msn5_chem$SRP_ugL[msn5_chem$Depth_m==0.1]) 

euclidean_drivers_df$DOC_surf <- c(msn1_chem$DOC_mgL[msn1_chem$Depth_m==0.1], 
                                  msn2_chem$DOC_mgL[msn2_chem$Depth_m==0.1],
                                  msn3_chem$DOC_mgL[msn3_chem$Depth_m==0.1],
                                  msn4_chem$DOC_mgL[msn4_chem$Depth_m==0.1],
                                  msn5_chem$DOC_mgL[msn5_chem$Depth_m==0.1]) 


#----------------------------------------------------------------------------------------#
#calculate avg density col for each day at both sites
euclidean_drivers_df$pel_avg_dens <- c(mean(zoop_epi_tows$ZoopDensity_No.pL_1[zoop_epi_tows$site=="pel" & zoop_epi_tows$groups==1]),
                                   mean(zoop_epi_tows$ZoopDensity_No.pL_1[zoop_epi_tows$site=="pel" & zoop_epi_tows$groups==2]),
                                   mean(zoop_epi_tows$ZoopDensity_No.pL_1[zoop_epi_tows$site=="pel" & zoop_epi_tows$groups==3]),
                                   mean(zoop_epi_tows$ZoopDensity_No.pL_1[zoop_epi_tows$site=="pel" & zoop_epi_tows$groups==4]),
                                   mean(zoop_epi_tows$ZoopDensity_No.pL_1[zoop_epi_tows$site=="pel" & zoop_epi_tows$groups==5]))


euclidean_drivers_df$lit_avg_dens <- c(mean(zoop_epi_tows$ZoopDensity_No.pL_1[zoop_epi_tows$site=="lit" & zoop_epi_tows$groups==1]),
                                  mean(zoop_epi_tows$ZoopDensity_No.pL_1[zoop_epi_tows$site=="lit" & zoop_epi_tows$groups==2]),
                                  mean(zoop_epi_tows$ZoopDensity_No.pL_1[zoop_epi_tows$site=="lit" & zoop_epi_tows$groups==3]),
                                  mean(zoop_epi_tows$ZoopDensity_No.pL_1[zoop_epi_tows$site=="lit" & zoop_epi_tows$groups==4]),
                                  mean(zoop_epi_tows$ZoopDensity_No.pL_1[zoop_epi_tows$site=="lit" & zoop_epi_tows$groups==5]))

euclidean_drivers_df$oxycline_depth <- c(6.5,5.5,3.5,6.5,5) #see calcs below

#plot response variable (euclidean distances) against environmental/biological data
#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_epiDO.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$DO_0.1m,euclidean_drivers_df$pel_euc_dist, xlab="Epi DO", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5, ylim=c(0,3))
points(euclidean_drivers_df$DO_0.1m[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$DO_0.1m[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$DO_0.1m[3],euclidean_drivers_df$pel_euc_dist[3], bg="#FBF2C4",pch=21,cex=3)
points(euclidean_drivers_df$DO_0.1m[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$DO_0.1m[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=hcl.colors(5,"Geyser"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_hypoDO.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$DO_9.0m,euclidean_drivers_df$pel_euc_dist, xlab="Hypo DO", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$DO_9.0m[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$DO_9.0m[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$DO_9.0m[3],euclidean_drivers_df$pel_euc_dist[3], bg="#FBF2C4",pch=21,cex=3)
points(euclidean_drivers_df$DO_9.0m[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$DO_9.0m[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=hcl.colors(5,"Geyser"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_epitemp.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$temp_0.1m,euclidean_drivers_df$pel_euc_dist, xlab="Epi Temp", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$temp_0.1m[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$temp_0.1m[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$temp_0.1m[3],euclidean_drivers_df$pel_euc_dist[3], bg="#FBF2C4",pch=21,cex=3)
points(euclidean_drivers_df$temp_0.1m[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$temp_0.1m[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=hcl.colors(5,"Geyser"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_hypotemp.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$temp_9.0m,euclidean_drivers_df$pel_euc_dist, xlab="Hypo Temp", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$temp_9.0m[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$temp_9.0m[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$temp_9.0m[3],euclidean_drivers_df$pel_euc_dist[3], bg="#FBF2C4",pch=21,cex=3)
points(euclidean_drivers_df$temp_9.0m[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$temp_9.0m[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=hcl.colors(5,"Geyser"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_epispcond.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$Spcond_0.1m,euclidean_drivers_df$pel_euc_dist, xlab="Epi Sp Cond", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$Spcond_0.1m[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$Spcond_0.1m[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$Spcond_0.1m[3],euclidean_drivers_df$pel_euc_dist[3], bg="#FBF2C4",pch=21,cex=3)
points(euclidean_drivers_df$Spcond_0.1m[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$Spcond_0.1m[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=hcl.colors(5,"Geyser"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_hypospcond.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$Spcond_9.0m,euclidean_drivers_df$pel_euc_dist, xlab="Hypo Sp Cond", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$Spcond_9.0m[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$Spcond_9.0m[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$Spcond_9.0m[3],euclidean_drivers_df$pel_euc_dist[3], bg="#FBF2C4",pch=21,cex=3)
points(euclidean_drivers_df$Spcond_9.0m[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$Spcond_9.0m[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=hcl.colors(5,"Geyser"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_epichl.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$chla_0.1m,euclidean_drivers_df$pel_euc_dist, xlab="Epi Chl a", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$chla_0.1m[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$chla_0.1m[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$chla_0.1m[3],euclidean_drivers_df$pel_euc_dist[3], bg="#FBF2C4",pch=21,cex=3)
points(euclidean_drivers_df$chla_0.1m[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$chla_0.1m[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=hcl.colors(5,"Geyser"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_hypochl.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$chla_9.0m,euclidean_drivers_df$pel_euc_dist, xlab="Hypo Chl a", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$chla_9.0m[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$chla_9.0m[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$chla_9.0m[3],euclidean_drivers_df$pel_euc_dist[3], bg="#FBF2C4",pch=21,cex=3)
points(euclidean_drivers_df$chla_9.0m[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$chla_9.0m[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=hcl.colors(5,"Geyser"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_epiturb.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$turb_0.1m,euclidean_drivers_df$pel_euc_dist, xlab="Epi Turb", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$turb_0.1m[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$turb_0.1m[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$turb_0.1m[3],euclidean_drivers_df$pel_euc_dist[3], bg="#FBF2C4",pch=21,cex=3)
points(euclidean_drivers_df$turb_0.1m[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$turb_0.1m[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=hcl.colors(5,"Geyser"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_hypoturb.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$turb_9.0m,euclidean_drivers_df$pel_euc_dist, xlab="Hypo Turb", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$turb_9.0m[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$turb_9.0m[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$turb_9.0m[3],euclidean_drivers_df$pel_euc_dist[3], bg="#FBF2C4",pch=21,cex=3)
points(euclidean_drivers_df$turb_9.0m[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$turb_9.0m[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=hcl.colors(5,"Geyser"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_epipar.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$par_0.1m,euclidean_drivers_df$pel_euc_dist, xlab="Epi PAR", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$par_0.1m[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$par_0.1m[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$par_0.1m[3],euclidean_drivers_df$pel_euc_dist[3], bg="#FBF2C4",pch=21,cex=3)
points(euclidean_drivers_df$par_0.1m[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$par_0.1m[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=hcl.colors(5,"Geyser"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_hypopar.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$par_9.0m,euclidean_drivers_df$pel_euc_dist, xlab="Hypo PAR", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$par_9.0m[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$par_9.0m[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$par_9.0m[3],euclidean_drivers_df$pel_euc_dist[3], bg="#FBF2C4",pch=21,cex=3)
points(euclidean_drivers_df$par_9.0m[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$par_9.0m[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=hcl.colors(5,"Geyser"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_pelzoopdens.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$pel_avg_dens,euclidean_drivers_df$pel_euc_dist, xlab="Pel Zoop Dens", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$pel_avg_dens[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$pel_avg_dens[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$pel_avg_dens[3],euclidean_drivers_df$pel_euc_dist[3], bg="#FBF2C4",pch=21,cex=3)
points(euclidean_drivers_df$pel_avg_dens[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$pel_avg_dens[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=hcl.colors(5,"Geyser"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_litzoopdens.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$lit_avg_dens,euclidean_drivers_df$pel_euc_dist, xlab="Lit Zoop Dens", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$lit_avg_dens[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$lit_avg_dens[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$lit_avg_dens[3],euclidean_drivers_df$pel_euc_dist[3], bg="#FBF2C4",pch=21,cex=3)
points(euclidean_drivers_df$lit_avg_dens[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$lit_avg_dens[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=hcl.colors(5,"Geyser"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_thermdepth.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$thermo_depth,euclidean_drivers_df$pel_euc_dist, xlab="Thermocline Depth", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$thermo_depth[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$thermo_depth[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$thermo_depth[3],euclidean_drivers_df$pel_euc_dist[3], bg="#FBF2C4",pch=21,cex=3)
points(euclidean_drivers_df$thermo_depth[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$thermo_depth[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=hcl.colors(5,"Geyser"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_TN.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$TN_surf,euclidean_drivers_df$pel_euc_dist, xlab="TN (ug/L)", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$TN_surf[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$TN_surf[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$TN_surf[3],euclidean_drivers_df$pel_euc_dist[3], bg="#FBF2C4",pch=21,cex=3)
points(euclidean_drivers_df$TN_surf[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$TN_surf[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=hcl.colors(5,"Geyser"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_TP.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$TP_surf,euclidean_drivers_df$pel_euc_dist, xlab="TP (ug/L)", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$TP_surf[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$TP_surf[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$TP_surf[3],euclidean_drivers_df$pel_euc_dist[3], bg="#FBF2C4",pch=21,cex=3)
points(euclidean_drivers_df$TP_surf[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$TP_surf[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=hcl.colors(5,"Geyser"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_NH4.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$NH4_surf,euclidean_drivers_df$pel_euc_dist, xlab="NH4 (ug/L)", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$NH4_surf[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$NH4_surf[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$NH4_surf[3],euclidean_drivers_df$pel_euc_dist[3], bg="#FBF2C4",pch=21,cex=3)
points(euclidean_drivers_df$NH4_surf[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$NH4_surf[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=hcl.colors(5,"Geyser"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_NO3NO2.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$NO3NO2_surf,euclidean_drivers_df$pel_euc_dist, xlab="NO3NO2 (ug/L)", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$NO3NO2_surf[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$NO3NO2_surf[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$NO3NO2_surf[3],euclidean_drivers_df$pel_euc_dist[3], bg="#FBF2C4",pch=21,cex=3)
points(euclidean_drivers_df$NO3NO2_surf[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$NO3NO2_surf[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=hcl.colors(5,"Geyser"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_SRP.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$SRP_surf,euclidean_drivers_df$pel_euc_dist, xlab="SRP (ug/L)", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$SRP_surf[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$SRP_surf[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$SRP_surf[3],euclidean_drivers_df$pel_euc_dist[3], bg="#FBF2C4",pch=21,cex=3)
points(euclidean_drivers_df$SRP_surf[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$SRP_surf[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=hcl.colors(5,"Geyser"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_DOC.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$DOC_surf,euclidean_drivers_df$pel_euc_dist, xlab="DOC (mg/L)", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$DOC_surf[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$DOC_surf[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$DOC_surf[3],euclidean_drivers_df$pel_euc_dist[3], bg="#FBF2C4",pch=21,cex=3)
points(euclidean_drivers_df$DOC_surf[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$DOC_surf[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=hcl.colors(5,"Geyser"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_oxydepth.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$oxycline_depth,euclidean_drivers_df$pel_euc_dist, xlab="Oxycline Depth", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$oxycline_depth[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$oxycline_depth[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$oxycline_depth[3],euclidean_drivers_df$pel_euc_dist[3], bg="#FBF2C4",pch=21,cex=3)
points(euclidean_drivers_df$oxycline_depth[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$oxycline_depth[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=hcl.colors(5,"Geyser"),bty = "n",cex=1.4)
#dev.off()


#ctd 
ctd_clean <- ctd %>%
  dplyr::mutate(rdepth = plyr::round_any(Depth_m, 0.5)) %>% 
  dplyr::filter(!is.na(DO_mgL)) %>%
  dplyr::group_by(Reservoir,Date,Site, rdepth) %>%
  dplyr::summarise(value = mean(DO_mgL)) %>% 
  dplyr::rename(depth = rdepth)  %>%
  dplyr::group_by(Reservoir,Site, Date) %>% 
  dplyr::mutate(oxy = min(depth[value<=2], na.rm=TRUE))

depths <- ctd_clean$oxy[(ctd_clean$Date=="2019-07-10" | ctd_clean$Date=="2019-07-24" |
                          ctd_clean$Date=="2020-08-12" | ctd_clean$Date=="2021-06-16" | ctd_clean$Date=="2021-07-12") & ctd_clean$Reservoir=="BVR" & ctd_clean$Site==50 & ctd_clean$depth>0]


#oxygen plots for all 5 MSNs
ggplot(subset(ctd_clean, depth > 0 & Reservoir=="BVR" & Site==50 & Date %in% c(as.Date("2019-07-10"), as.Date("2019-07-24"),
                                             as.Date("2020-08-12"), as.Date("2021-06-16"), as.Date("2021-07-12"))), aes(value,depth,color=as.factor(Date))) + 
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin= depths, ymax=Inf), fill="red",alpha=0.03) +
  scale_color_manual("",values=hcl.colors(5,"Geyser"), guide="none") +  xlab("DO (mg/L)") + ylab("Depth (m)") +
  ylim(10,0) + geom_point() + geom_path() + theme_bw() + facet_wrap(~Date, ncol=5) +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = c(0.76,0.02),
        legend.background = element_blank(),legend.direction = "horizontal", panel.grid.minor = element_blank(), legend.key=element_rect(fill=NA),
        plot.margin = unit(c(0,0.05,0,0), "cm"),legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),
        legend.title = element_text(size = 4.5),legend.text  = element_text(size = 4.5), panel.spacing=unit(0, "cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6), axis.text.y = element_text(size=6)) 
#ggsave(file.path(getwd(),"Summer2021-DataAnalysis/Figures/DO_profiles.jpg"), width=4, height=3) 

       