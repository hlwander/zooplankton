#multivariate stats for 2019 - 2021 zoop data
#created 25Nov2021

#read in libraries
pacman::p_load(dplyr, vegan, labdsv, goeveg, rLakeAnalyzer, ggplot2,tidyr,viridis, egg, ggordiplots, splancs, ggpubr, FSA, rcompanion)

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

#manually change a noon and midnight hour so  summarizing below actually works
zoops$Hour[zoops$sample_ID=="B_pel_10Jul19_noon_epi_rep1"] <- "12:00"
zoops$Hour[zoops$sample_ID=="B_pel_08Jul21_midnight_epi_rep1"] <- "0:00"

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
ord <- ordiplot(NMDS_temporal_bray,display = c('sites'),choices = c(1,2),type = "n") 
points(NMDS_temporal_bray$points[zoop_epi_tows$groups=="1",1], NMDS_temporal_bray$points[zoop_epi_tows$groups=="1",2], pch=21,bg="#008585")
points(NMDS_temporal_bray$points[zoop_epi_tows$groups=="2",1], NMDS_temporal_bray$points[zoop_epi_tows$groups=="2",2], pch=21,bg="#9BBAA0")
points(NMDS_temporal_bray$points[zoop_epi_tows$groups=="3",1], NMDS_temporal_bray$points[zoop_epi_tows$groups=="3",2], pch=21,bg="#F2E2B0")
points(NMDS_temporal_bray$points[zoop_epi_tows$groups=="4",1], NMDS_temporal_bray$points[zoop_epi_tows$groups=="4",2], pch=21,bg="#DEA868")
points(NMDS_temporal_bray$points[zoop_epi_tows$groups=="5",1], NMDS_temporal_bray$points[zoop_epi_tows$groups=="5",2], pch=21,bg="#C7522B")
legend("bottomright", legend=c('Day1','Day2','Day3','Day4','Day5'), pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"), bty = "n") 
ordihull(ord, zoop_epi_tows$groups, display = "sites", draw = c("polygon"),
         col = c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"), alpha = 75,cex = 2)
#text(NMDS_temporal_bray$species[,2],NMDS_temporal_bray$species[,4], labels = c("Calanoida","Cyclopoida","Keratella","Kellicottia",
#     "Bosmina","Daphnia","Ceriodaphnia","nauplius","Collothecidae","Synchaetidae","Conochilidae"), cex=0.9)
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
#                                 NMDS ms figs                                  #
#-------------------------------------------------------------------------------#

ord <- ordiplot(NMDS_temporal_avg_bray,display = c('sites','species'),choices = c(1,2),type = "n")
sites <- gg_ordiplot(ord, zoop_avg$site, kind = "ehull", 
                    ellipse=FALSE, hull = TRUE, plot = FALSE, pt.size=0.9) 

NMDS_site <- sites$plot + geom_point() + theme_bw() + 
                geom_polygon(data = sites$df_hull, aes(x = x, y = y, fill = Group), alpha=0.2) +
                xlim(c(-0.53, 0.55)) + ylim(c(-0.7,0.64)) +
                theme(text = element_text(size=7), axis.text = element_text(size=7, color="black"), 
                      legend.background = element_blank(), 
                      legend.key.height=unit(0.3,"line"),
                      legend.key = element_blank(),
                      axis.text.x = element_text(vjust = 0.5), 
                      strip.background = element_rect(fill = "transparent"), 
                      legend.position = c(0.14,0.1), legend.spacing = unit(-0.5, 'cm'),
                      plot.margin = unit(c(0,-0.1,0,0), 'lines'),
                      panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
                      legend.key.width =unit(0.1,"line")) + guides(fill="none") +
                scale_fill_manual("",values=c("#882255","#3399CC"))+
                scale_color_manual("",values=c("#882255","#3399CC"),
                                   label=c('littoral','pelagic'))


days <- gg_ordiplot(ord, zoop_avg$groups, kind = "ehull", 
                         ellipse=FALSE, hull = TRUE, plot = FALSE, pt.size=0.9) 

NMDS_day <- days$plot + geom_point() + theme_bw() + geom_path() + ylab(NULL) +
              geom_polygon(data = days$df_hull, aes(x = x, y = y, fill = Group), alpha=0.2) +
              xlim(c(-0.53, 0.55)) + ylim(c(-0.7,0.64)) +
              theme(text = element_text(size=7), axis.text = element_text(size=7, color="black"), 
                    legend.background = element_blank(), 
                    legend.key.height=unit(0.3,"line"), 
                    legend.key = element_blank(),
                    axis.text.x = element_text(vjust = 0.5), 
                    axis.text.y=element_blank(),
                    axis.ticks.y = element_blank(),
                    strip.background = element_rect(fill = "transparent"), 
                    legend.position = c(0.25,0.16), legend.spacing = unit(-0.5, 'cm'),
                    plot.margin = unit(c(0,-0.1,0,-0.1), 'lines'),
                    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
                    legend.key.width =unit(0.1,"line")) + guides(fill="none") +
              scale_fill_manual("",values=c("#008585","#89B199","#EFECBF","#DB9B5A","#C7522B"))+
              scale_color_manual("",values=c("#008585","#89B199","#EFECBF","#DB9B5A","#C7522B"),
                                 label=c('10-11 Jul 2019', '24-25 Jul 2019','12-13 Aug 2020',
                                        '15-16 Jun 2021', '7-8 Jul 2021'))
  
  
hours <- gg_ordiplot(ord, zoop_avg$order, kind = "ehull", 
                    ellipse=FALSE, hull = TRUE, plot = FALSE, pt.size=0.9) 

NMDS_hour <- hours$plot + geom_point() + theme_bw() + geom_path() + ylab(NULL) +
              geom_polygon(data = hours$df_hull, aes(x = x, y = y, fill = Group), alpha=0.2) +
              xlim(c(-0.53, 0.55)) + ylim(c(-0.7,0.64)) +
              theme(text = element_text(size=7), axis.text = element_text(size=7, color="black"), 
                    legend.background = element_blank(), 
                    legend.key.height=unit(0.3,"line"), 
                    legend.key = element_blank(),
                    axis.text.x = element_text(vjust = 0.5), 
                    axis.text.y=element_blank(),
                    axis.ticks.y = element_blank(),
                    strip.background = element_rect(fill = "transparent"), 
                    legend.position = c(0.12,0.28), legend.spacing = unit(-0.5, 'cm'),
                    plot.margin = unit(c(0,0,0,-0.1), 'lines'),
                    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
                    legend.key.width =unit(0.1,"line")) + guides(fill="none") +
              scale_fill_manual("",values=hcl.colors(11,"sunset"))+
              scale_color_manual("",values=hcl.colors(11,"sunset"),
                                 label=c('12pm','6pm','7pm','8pm','9pm','12am',
                                         '4am','5am','6am','7am','12pm'))

fig5 <- egg::ggarrange(NMDS_site, NMDS_day, NMDS_hour, nrow=1)
ggsave(file.path(getwd(),"Summer2021-DataAnalysis/Figures/NMDS_multipanel_2v1.jpg"),
       fig5, width=5, height=2) 

#ggsave(file.path(getwd(),"Summer2021-DataAnalysis/Figures/NMDS_hours_2v1.jpg"),
#       NMDS_hour, width=5, height=2) 

#-----------------------------------------------------------------------------------------#
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
#technically calculating the euclidean distances from bray-curtis distance matrices

#step 1: take NMDS output for each site using NMDS coordinates
zoop_euc <- as.matrix(vegdist(NMDS_temporal_avg_bray$points, method='euclidean'))

#get collect_date into correct format
zoop_epi_tows$collect_date <- as.Date(zoop_epi_tows$collect_date)

#order zoop epi tows by hour, MSN, and site
zoop_epi_tows <- zoop_epi_tows %>% dplyr::arrange(site, groups, order)

#step 2: select and sum the 9 distances between connecting points for each of the 5 days

lit_day1 <- sum(zoop_euc[1,2],zoop_euc[2,3],zoop_euc[3,4],zoop_euc[4,5],zoop_euc[5,6],
                zoop_euc[6,7],zoop_euc[7,8],zoop_euc[8,9],zoop_euc[9,10],zoop_euc[10,11])

lit_day2 <- sum(zoop_euc[12,13], zoop_euc[13,14],zoop_euc[14,15],zoop_euc[15,16],zoop_euc[16,17],
                zoop_euc[17,18],zoop_euc[18,19],zoop_euc[19,20],zoop_euc[20,21],zoop_euc[21,22])

lit_day3 <- sum(zoop_euc[23,24], zoop_euc[24,25],zoop_euc[25,26],zoop_euc[26,27],zoop_euc[27,28],
                zoop_euc[28,29],zoop_euc[29,30],zoop_euc[30,31],zoop_euc[31,32],zoop_euc[32,33])

lit_day4 <- sum(zoop_euc[34,35], zoop_euc[35,36],zoop_euc[36,37],zoop_euc[37,38],zoop_euc[38,39],
                zoop_euc[39,40],zoop_euc[40,41],zoop_euc[41,42],zoop_euc[42,43],zoop_euc[43,44])

lit_day5 <- sum(zoop_euc[45,46], zoop_euc[46,47],zoop_euc[47,48],zoop_euc[48,49],zoop_euc[49,50],
                zoop_euc[50,51],zoop_euc[51,52],zoop_euc[52,53],zoop_euc[53,54],zoop_euc[54,55])

mean(lit_day1,lit_day2,lit_day3,lit_day4,lit_day5)

pel_day1 <- sum(zoop_euc[56,57],zoop_euc[57,58],zoop_euc[58,59],zoop_euc[59,60],zoop_euc[60,61],
                zoop_euc[61,62],zoop_euc[62,63],zoop_euc[63,64],zoop_euc[64,65],zoop_euc[65,66])

pel_day2 <- sum(zoop_euc[67,68], zoop_euc[68,69],zoop_euc[69,70],zoop_euc[70,71],zoop_euc[71,72],
                zoop_euc[72,73],zoop_euc[73,74],zoop_euc[74,75],zoop_euc[75,76],zoop_euc[76,77])

pel_day3 <- sum(zoop_euc[78,79], zoop_euc[79,80],zoop_euc[80,81],zoop_euc[81,82],zoop_euc[82,83],
                zoop_euc[83,84],zoop_euc[84,85],zoop_euc[85,86],zoop_euc[86,87],zoop_euc[87,88])

pel_day4 <- sum(zoop_euc[89,90], zoop_euc[90,91],zoop_euc[91,92],zoop_euc[92,93],zoop_euc[93,94],
                zoop_euc[94,95],zoop_euc[95,96],zoop_euc[96,97],zoop_euc[97,98],zoop_euc[98,99])

pel_day5 <- sum(zoop_euc[100,101], zoop_euc[101,102],zoop_euc[102,103],zoop_euc[103,104],zoop_euc[104,105],
                zoop_euc[105,106],zoop_euc[106,107],zoop_euc[107,108],zoop_euc[108,109],zoop_euc[109,110])

mean(pel_day1,pel_day2,pel_day3,pel_day4,pel_day5) #community structure is more variable at pelagic site than littoral

#convert ED matrix back into distance structure for next steps
zoop_euc <- vegdist(NMDS_temporal_avg_bray$points, method='euclidean')

#Now calculate the centroids of each polygon AND the avg distance of each point to its polygon centroid
centroids_sites <- betadisper(zoop_euc, group = as.factor(zoop_epi_tows$site), type="centroid")
centroids_hours <- betadisper(zoop_euc, group = as.factor(zoop_epi_tows$order), type="centroid")
centroids_days <-  betadisper(zoop_euc, group = as.factor(zoop_epi_tows$groups), type="centroid")

#-------------------------------------------------------------------------------#
#METHOD 1:average distance of each point to polygon centroid (dispersion approach)

#site variability 
disp_site <- mean(centroids_sites$group.distances)
disp_site_sd <- sd(centroids_sites$group.distances)

#hourly variability- MOST VARIABLE!
disp_hour <- mean(centroids_hours$group.distances)
disp_hour_sd <- sd(centroids_hours$group.distances)

#daily variability - LEAST VARIABLE!
disp_day <- mean(centroids_days$group.distances)
disp_day_sd <- sd(centroids_days$group.distances)

#-------------------------------------------------------------------------------#
#METHOD 2: average distance between all combinations of centroids (pairwise approach)

#site variability - MOST VARIABLE
pair_site <- mean(dist(centroids_sites$centroids))
pair_site_sd <- sd(dist(centroids_sites$centroids)) #NA

#hourly variability
pair_hour <- mean(dist(centroids_hours$centroids))
pair_hour_sd <- sd(dist(centroids_hours$centroids))

#annual variability - LEAST VARIABLE!
pair_day <- mean(dist(centroids_days$group.distances))
pair_day_sd <- mean(dist(centroids_days$group.distances))

#-------------------------------------------------------------------------------#
#METHOD 3: average areas of polygons (avg_area) - use NMDS df for this

#site areas             
area_lit <- areapl(cbind(sites$df_hull$x[sites$df_hull$Group=="lit"],
                             sites$df_hull$y[sites$df_hull$Group=="lit"]))

area_pel <- areapl(cbind(sites$df_hull$x[sites$df_hull$Group=="pel"],
                             sites$df_hull$y[sites$df_hull$Group=="pel"]))

avg_area_site <- mean(c(area_lit, area_pel))
avg_sd_site <- sd(c(area_lit, area_pel))
            
#day areas
area_day1 <- areapl(cbind(days$df_hull$x[days$df_hull$Group==1],
                               days$df_hull$y[days$df_hull$Group==1]))

area_day2 <- areapl(cbind(days$df_hull$x[days$df_hull$Group==2],
                               days$df_hull$y[days$df_hull$Group==2]))

area_day3 <- areapl(cbind(days$df_hull$x[days$df_hull$Group==3],
                               days$df_hull$y[days$df_hull$Group==3]))

area_day4 <- areapl(cbind(days$df_hull$x[days$df_hull$Group==4],
                               days$df_hull$y[days$df_hull$Group==4]))

area_day5 <- areapl(cbind(days$df_hull$x[days$df_hull$Group==5],
                               days$df_hull$y[days$df_hull$Group==5]))

avg_area_day <- mean(c(area_day1, area_day2, area_day3, area_day4, area_day5))
avg_sd_day <- sd(c(area_day1, area_day2, area_day3, area_day4, area_day5))

#hour areas
area_hour1 <- areapl(cbind(hours$df_hull$x[hours$df_hull$Group==1],
                           hours$df_hull$y[hours$df_hull$Group==1]))

area_hour2 <- areapl(cbind(hours$df_hull$x[hours$df_hull$Group==2],
                           hours$df_hull$y[hours$df_hull$Group==2]))

area_hour3 <- areapl(cbind(hours$df_hull$x[hours$df_hull$Group==3],
                           hours$df_hull$y[hours$df_hull$Group==3]))

area_hour4 <- areapl(cbind(hours$df_hull$x[hours$df_hull$Group==4],
                           hours$df_hull$y[hours$df_hull$Group==4]))

area_hour5 <- areapl(cbind(hours$df_hull$x[hours$df_hull$Group==5],
                           hours$df_hull$y[hours$df_hull$Group==5]))

area_hour6 <- areapl(cbind(hours$df_hull$x[hours$df_hull$Group==6],
                           hours$df_hull$y[hours$df_hull$Group==6]))

area_hour7 <- areapl(cbind(hours$df_hull$x[hours$df_hull$Group==7],
                           hours$df_hull$y[hours$df_hull$Group==7]))

area_hour8 <- areapl(cbind(hours$df_hull$x[hours$df_hull$Group==8],
                           hours$df_hull$y[hours$df_hull$Group==8]))

area_hour9 <- areapl(cbind(hours$df_hull$x[hours$df_hull$Group==9],
                           hours$df_hull$y[hours$df_hull$Group==9]))

area_hour10 <- areapl(cbind(hours$df_hull$x[hours$df_hull$Group==10],
                            hours$df_hull$y[hours$df_hull$Group==10]))

area_hour11 <- areapl(cbind(hours$df_hull$x[hours$df_hull$Group==11],
                            hours$df_hull$y[hours$df_hull$Group==11]))

avg_area_hour <- mean(c(area_hour1, area_hour2, area_hour3, area_hour4, area_hour5, area_hour6,
                        area_hour7, area_hour8, area_hour9, area_hour10, area_hour11))
avg_sd_hour <- sd(c(area_hour1, area_hour2, area_hour3, area_hour4, area_hour5, area_hour6,
                    area_hour7, area_hour8, area_hour9, area_hour10, area_hour11))

#-------------------------------------------------------------------------------#
#METHOD 4: max area - min area of polygons (area_diff) 
#note - don't think I can calculate sd of difference because idk what the sample mean is

area_diff_site <- area_lit - area_pel

area_diff_day <- area_day2 - area_day4

area_diff_hour <- area_hour7 - area_hour11

#-------------------------------------------------------------------------------
#put variability values into a dataset
euc_distances_df <- data.frame("Method" = c("Dispersion","Pairwise", "Average Area"),
                               "Site" = c(paste0(round(disp_site,2)," ± ",round(disp_site_sd,2)),
                                          paste0(round(pair_site,2)),
                                          paste0(round(avg_area_site,2), " ± ", round(avg_sd_site,2))),
                               "Day" = c(paste0(round(disp_day,2)," ± ",round(disp_day_sd,2)),
                                         paste0(round(pair_day,2)," ± ",round(pair_day_sd,2)),
                                         paste0(round(avg_area_day,2), " ± ", round(avg_sd_day,2))),
                               "Hour" = c(paste0(round(disp_hour,2)," ± ",round(disp_hour_sd,2)),
                                          paste0(round(pair_hour,2)," ± ",round(pair_hour_sd,2)),
                                          paste0(round(avg_area_hour,2), " ± ", round(avg_sd_hour,2))))
  
#write.csv(euc_distances_df, file.path(getwd(),"/Summer2021-DataAnalysis/SummaryStats/Euclidean_distances.csv"))


#-------------------------------------------------------------------------------#
#Kruskal-wallis test to determine if group means are significant (can only do with dispersion and area)
disp_site_df <- data.frame("Group" = c("site","site"),
                           "dist" = c(centroids_sites$group.distances),
                           "avg_area" = c(area_lit, area_pel)) 

disp_hours_df <- data.frame("Group"=c(rep("hour",11)), 
                            "dist" = c(centroids_hours$group.distances),
                            "avg_area" = c(area_hour1, area_hour2, area_hour3, area_hour4, area_hour5, area_hour6,
                                           area_hour7, area_hour8, area_hour9, area_hour10, area_hour11))

disp_days_df <- data.frame("Group"=c(rep("day",5)), 
                           "dist" = c(centroids_days$group.distances),
                            "avg_area" = c(area_day1, area_day2, area_day3,
                                           area_day4, area_day5))

disp_df <- rbind(disp_site_df,disp_hours_df, disp_days_df)


#KW test to see if groups are significant
kruskal.test(dist ~ Group, data = disp_df) #not significant
kruskal.test(avg_area ~ Group, data = disp_df) #significant sometimes??

#now dunn test to determine which areas are different from each other
dunnTest(avg_area ~ as.factor(Group),
         data=disp_df,
         method="bonferroni")
#site avg areas are sig larger from day avg area

ggboxplot(disp_df, x = "Group", y = "avg_area", 
          color = "Group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("site", "hour", "day"),
          ylab = "Average Area", xlab = "Group")


#-------------------------------------------------------------------------------#
#Calculate within site/day/hour variability as the dispersion and area of each polygon from above (n=2 methods)

#polygon dispersion
lit_disp <- centroids_sites$group.distances[1]
pel_disp <- centroids_sites$group.distances[2]

day1_disp <- centroids_days$group.distances[1]
day2_disp <- centroids_days$group.distances[2]
day3_disp <- centroids_days$group.distances[3]
day4_disp <- centroids_days$group.distances[4]
day5_disp <- centroids_days$group.distances[5]

hour1_disp <- centroids_hours$group.distances[1]
hour2_disp <- centroids_hours$group.distances[2]
hour3_disp <- centroids_hours$group.distances[3]
hour4_disp <- centroids_hours$group.distances[4]
hour5_disp <- centroids_hours$group.distances[5]
hour6_disp <- centroids_hours$group.distances[6]
hour7_disp <- centroids_hours$group.distances[7]
hour8_disp <- centroids_hours$group.distances[8]
hour9_disp <- centroids_hours$group.distances[9]
hour10_disp <- centroids_hours$group.distances[10]
hour11_disp <- centroids_hours$group.distances[11]

#-------------------------------------------------------------------------------
#dataframe with area and dispersion for each polygon
polygons_df <- data.frame("Polygon" = c("Littoral","Pelagic",
                                        "10-11 Jul 2019","24-25 Jul 2019","12-13 Aug 2020",
                                        "15-16 Jun 2021","7-8 Jul 2021",
                                        "12pm","6pm","7pm","8pm","9pm","12am",
                                        "4am","5am","6am","7am","12pm"),
                          "Area" = c(area_lit, area_pel, area_day1, area_day2, area_day3,
                                     area_day4, area_day5, area_hour1, area_hour2, area_hour3,
                                     area_hour4, area_hour5, area_hour6, area_hour7, area_hour8, 
                                     area_hour9, area_hour10, area_hour11),
                          "Dispersion" = c(lit_disp, pel_disp, day1_disp, day2_disp, day3_disp,
                                           day4_disp, day5_disp, hour1_disp, hour2_disp, hour3_disp, 
                                           hour4_disp, hour5_disp, hour6_disp, hour7_disp, hour8_disp,
                                           hour9_disp, hour10_disp, hour11_disp))

#write.csv(polygons_df, file.path(getwd(),"/Summer2021-DataAnalysis/SummaryStats/polygon_area_dispersion.csv"))

#calculate the range of all areas and dispersion values
range(polygons_df$Area[1:2])[2] - range(polygons_df$Area[1:2])[1] #smallest
range(polygons_df$Area[3:7])[2] - range(polygons_df$Area[3:7])[1] #largest
range(polygons_df$Area[8:18])[2] - range(polygons_df$Area[8:18])[1]

range(polygons_df$Dispersion[1:2])[2] - range(polygons_df$Dispersion[1:2])[1] #largest
range(polygons_df$Dispersion[3:7])[2] - range(polygons_df$Dispersion[3:7])[1] #smallest
range(polygons_df$Dispersion[8:18])[2] - range(polygons_df$Dispersion[8:18])[1]


#-------------------------------------------------------------------------------#
#dfs to calculate significance within sites, days, and hours
within_site_dist <- data.frame("group" = c(rep("lit",55),rep("pel",55)),
                                "dist" = c(centroids_sites$distances[zoop_epi_tows$site=="lit"],
                                          centroids_sites$distances[zoop_epi_tows$site=="pel"]))

within_day_dist <- data.frame("group" = c(rep("day1",22),rep("day2",22),rep("day3",22),
                                          rep("day4",22),rep("day5",22)),
                              "dist" = c(centroids_days$distances[zoop_epi_tows$groups==1],
                                         centroids_days$distances[zoop_epi_tows$groups==2],
                                         centroids_days$distances[zoop_epi_tows$groups==3],
                                         centroids_days$distances[zoop_epi_tows$groups==4],
                                         centroids_days$distances[zoop_epi_tows$groups==5]))

within_hour_dist <- data.frame("group" = c(rep("hour1",10),rep("hour2",10),rep("hour3",10),
                                           rep("hour4",10), rep("hour5",10), rep("hour6",10),
                                           rep("hour7",10), rep("hour8",10), rep("hour9",10),
                                           rep("hour10",10), rep("hour11",10)),
                               "dist" = c(centroids_hours$distances[zoop_epi_tows$order==1],
                                          centroids_hours$distances[zoop_epi_tows$order==2],
                                          centroids_hours$distances[zoop_epi_tows$order==3],
                                          centroids_hours$distances[zoop_epi_tows$order==4],
                                          centroids_hours$distances[zoop_epi_tows$order==5],
                                          centroids_hours$distances[zoop_epi_tows$order==6],
                                          centroids_hours$distances[zoop_epi_tows$order==7],
                                          centroids_hours$distances[zoop_epi_tows$order==8],
                                          centroids_hours$distances[zoop_epi_tows$order==9],
                                          centroids_hours$distances[zoop_epi_tows$order==10],
                                          centroids_hours$distances[zoop_epi_tows$order==11]))


#now kw tests for significance 
kruskal.test(dist ~ group, data = within_site_dist) #sig! - so pelagic and littoral are different
kruskal.test(dist ~ group, data = within_day_dist) #sig
kruskal.test(dist ~ group, data = within_hour_dist) #nope

#dunn test for days
dunn_day_disp <- dunnTest(dist ~ as.factor(group),
                  data = within_day_dist,
                  method="bonferroni")
#after adj, p-values are all not significant, but before, day1+3, day2+3, day3+4, day3+5 are all different

#letters
cldList(P.unadj ~ Comparison, data=dunn_day_disp$res, threshold = 0.05) #honestly don't really believe this - n=22 isn't a tremendous amount of points...

ggboxplot(within_day_dist, x = "group", y = "dist", 
          color = "group", palette = hcl.colors(5,"earth"),
          order = c("day1", "day2", "day3", "day4", "day5"),
          ylab = "Distance to centroid", xlab = "")



#-------------------------------------------------------------------------------#
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

