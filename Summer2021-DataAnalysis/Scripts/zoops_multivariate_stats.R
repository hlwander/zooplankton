#multivariate stats for 2019 - 2021 zoop data
#created 25Nov2021

#read in libraries
pacman::p_load(dplyr, vegan, labdsv, goeveg, rLakeAnalyzer, ggplot2,tidyr,
               viridis, egg, ggordiplots, splancs, ggpubr, FSA, rcompanion)

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

#turn transformed community data into Euclidean distance matrix
zoop_euc <- as.matrix(vegdist(zoop_temporal_dens_avg_trans, method='euclidean'))
#zoop_bray <- as.matrix(vegdist(zoop_temporal_dens_avg_trans, method='bray'))

#scree plot to choose dimension #
dimcheckMDS(zoop_euc, distance = "bray", k = 6, trymax = 20, autotransform = TRUE)

#now do NMDS using averages w/ 4 dimensions for consistency
NMDS_temporal_avg_bray <- metaMDS(zoop_euc, distance='bray', k=4, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_temporal_avg_bray$stress

#-------------------------------------------------------------------------------#
#                                 NMDS ms figs                                  #
#-------------------------------------------------------------------------------#

ord <- ordiplot(NMDS_temporal_avg_bray,display = c('sites','species'),choices = c(1,2),type = "n")
sites <- gg_ordiplot(ord, zoop_avg$site, kind = "ehull", 
                     ellipse=FALSE, hull = TRUE, plot = FALSE, pt.size=0.9) 

NMDS_site <- sites$plot + geom_point() + theme_bw() + 
  geom_polygon(data = sites$df_hull, aes(x = x, y = y, fill = Group), alpha=0.2) +
  geom_point(data=sites$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2, fill=c("#882255","#3399CC")) +
  # xlim(c(-0.63, 0.6)) + ylim(c(-0.7,0.64)) +
  theme(text = element_text(size=7), axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), 
        legend.key.height=unit(0.3,"line"),
        legend.key = element_blank(),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-0,-0,-0,-0),
        legend.direction = "horizontal",
        axis.text.x = element_text(vjust = 0.5), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "top", legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(0,-0.1,0,0), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) + guides(fill="none") +
        annotate("text", x=-0.42, y=0.6, label= "bolditalic(a)", parse=T, size = 3) +
        scale_fill_manual("",values=c("#882255","#3399CC"))+
        scale_color_manual("",values=c("#882255","#3399CC"),
                     label=c('littoral','pelagic')) 

days <- gg_ordiplot(ord, zoop_avg$groups, kind = "ehull", 
                    ellipse=FALSE, hull = TRUE, plot = FALSE, pt.size=0.9) 

NMDS_day <- days$plot + geom_point() + theme_bw() + geom_path() + ylab(NULL) +
  geom_polygon(data = days$df_hull, aes(x = x, y = y, fill = Group), alpha=0.2) +
  geom_point(data=days$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2, 
             fill=c("#008585","#89B199","#EFECBF","#DB9B5A","#C7522B")) +
  # xlim(c(-0.63, 0.6)) + ylim(c(-0.7,0.64)) +
  theme(text = element_text(size=7), axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), 
        legend.key.height=unit(0.3,"line"), 
        legend.key = element_blank(),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-0,-0,-0,-0),
        legend.direction = "horizontal",
        axis.text.x = element_text(vjust = 0.5), 
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "top", legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(0,-0.1,0,-0.1), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) + 
        annotate("text", x=-0.42, y=0.6, label= "bolditalic(b)", parse=T, size = 3) +
        guides(fill="none", color = guide_legend(ncol=2)) +
        scale_fill_manual("",values=c("#008585","#89B199","#EFECBF","#DB9B5A","#C7522B"))+
        scale_color_manual("",values=c("#008585","#89B199","#EFECBF","#DB9B5A","#C7522B"),
                     label=c('10-11 Jul 2019', '24-25 Jul 2019','12-13 Aug 2020',
                             '15-16 Jun 2021', '7-8 Jul 2021'))


hours <- gg_ordiplot(ord, zoop_avg$order, kind = "ehull", 
                     ellipse=FALSE, hull = TRUE, plot = FALSE, pt.size=0.9) 
#order hours properly
hours$df_hull$Group <- factor(hours$df_hull$Group, levels = c(unique(hours$df_hull$Group)))

NMDS_hour <- hours$plot + geom_point() + theme_bw() + geom_path() + ylab(NULL) +
  geom_polygon(data = hours$df_hull, aes(x = x, y = y, fill = Group), alpha=0.2) +
  geom_point(data=hours$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2, fill=hcl.colors(11,"sunset")) +
  theme(text = element_text(size=7), axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), 
        legend.key.height=unit(0.3,"line"), 
        legend.key = element_blank(),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-0,-0,-0,-0),
        legend.direction = "horizontal",
        axis.text.x = element_text(vjust = 0.5), 
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "top", legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(0,0,0,-0.1), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) + guides(fill="none") +
        annotate("text", x=-0.42, y=0.6, label= "bolditalic(c)", parse=T, size=3) +
        scale_fill_manual("",values=hcl.colors(11,"sunset"))+
        scale_color_manual("",values=hcl.colors(11,"sunset"),
                     label=c('12pm','6pm','7pm','8pm','9pm','12am',
                             '4am','5am','6am','7am','12pm'))

fig5 <- egg::ggarrange(NMDS_site, NMDS_day, NMDS_hour, nrow=1)
ggsave(file.path(getwd(),"Summer2021-DataAnalysis/Figures/NMDS_multipanel_2v1.jpg"),
       fig5, width=5, height=3) 

#ggsave(file.path(getwd(),"Summer2021-DataAnalysis/Figures/NMDS_hours_2v1.jpg"),
#       NMDS_hour, width=5, height=2) 

#-------------------------------------------------------------------------------#
#track hours for each campaign
ord <- ordiplot(NMDS_temporal_avg_bray,display = c('sites','species'),choices = c(1,2),type = "n")
ordihull(ord, zoop_avg$groups, display = "sites", draw = c("polygon"),
         col = c("#008585","#89B199","#EFECBF","#DB9B5A","#C7522B"), alpha = 75,cex = 2)
lines(NMDS_temporal_avg_bray$points[zoop_avg$groups=="3"  & zoop_avg$site=="lit" ,1], 
      NMDS_temporal_avg_bray$points[(zoop_avg$groups=="3" & zoop_avg$site=="lit"),2], col="#EFECBF")

points(NMDS_temporal_avg_bray$points[(zoop_avg$groups=="3" & zoop_avg$site=="lit"),1], 
       NMDS_temporal_avg_bray$points[(zoop_avg$groups=="3" & zoop_avg$site=="lit"),2], 
       pch=as.character(zoop_avg$order),col="#EFECBF")

#legend("bottomright", legend=c('10-11 Jul 2019', '24-25 Jul 2019','12-13 Aug 2020','15-16 Jun 2021', '7-8 Jul 2021'), pch=21, 
#       pt.bg=c("#008585","#89B199","#EFECBF","#DB9B5A", "#C7522B"),bty = "n",cex=1.2) 
#legend("bottomleft", legend=c('12pm','6pm','7pm','8pm','9pm','12am','4am','5am','6am','7am','12pm'), pch=c(0,1,2,3,4,5,6,7,8,9,10) ,bty = "n",cex=1.2) 

#-------------------------------------------------------------------------------#
#                     Calculating euclidean distance                            #
#-------------------------------------------------------------------------------#
#step 1: use Euclidean distance matrix from transformed community data (zoop_euc)

#get collect_date into correct format
zoop_epi_tows$collect_date <- as.Date(zoop_epi_tows$collect_date)

#order zoop epi tows by hour, MSN, and site
zoop_epi_tows <- zoop_epi_tows %>% dplyr::arrange(site, groups, order)

#convert ED matrix back into distance structure for next steps
zoop_euc <- vegdist(zoop_temporal_dens_avg_trans, method='euclidean', upper = TRUE)

#Now calculate the centroids of each polygon AND the avg distance of each point to its polygon centroid
centroids_sites <- betadisper(zoop_euc, group = as.factor(zoop_epi_tows$site), type="centroid")
centroids_hours <- betadisper(zoop_euc, group = as.factor(zoop_epi_tows$order), type="centroid")
centroids_days <-  betadisper(zoop_euc, group = as.factor(zoop_epi_tows$groups), type="centroid")

#-------------------------------------------------------------------------------#
#METHOD 1:average distance of each point to polygon centroid (dispersion approach)
#METHOD 2: average distance between all combinations of centroids (pairwise approach)

#"bootstrapping" or randomly select 10 points per group and 2 groups and then calculating values using methods above

#create a df to store all the different variability values
var_results <- data.frame("site_disp"=rep(NA,500))

for (i in 1:500){ 
  #randomly select 10 points in each group
  ord_sub <- sample(unique(zoop_epi_tows$order), 2)
  groups_sub <-  sample(unique(zoop_epi_tows$groups), 2)
  
  zoop_sub <-  zoop_epi_tows %>% group_by(order) %>% 
               filter(order %in% c(ord_sub), groups %in% c(groups_sub)) %>%
               slice_sample(n=10)
  
  #only select data cols
  zoop_dens_sub <- zoop_sub[,c(grepl("density",colnames(zoop_sub)))] 
  
  #hellinger transformation
  zoop_dens_sub_trans <- hellinger(zoop_dens_sub)
  
  #convert ED matrix back into distance structure for next steps
  zoop_euc_sub <- vegdist(zoop_dens_sub_trans, method='euclidean', upper = TRUE)
  
  #Now calculate the centroids of each polygon AND the avg distance of each point to its polygon centroid
  centroids_sites_sub <- betadisper(zoop_euc_sub, group = as.factor(zoop_sub$site), type="centroid")
  centroids_days_sub <-  betadisper(zoop_euc_sub, group = as.factor(zoop_sub$groups), type="centroid")
  centroids_hours_sub <- betadisper(zoop_euc_sub, group = as.factor(zoop_sub$order), type="centroid")
  
  #distance calcs
  var_results$site_disp[i] <- mean(centroids_sites_sub$group.distances)
  var_results$day_disp[i] <- mean(centroids_days_sub$group.distances)
  var_results$hour_disp[i] <- mean(centroids_hours_sub$group.distances)
  
  var_results$site_pair[i] <- mean(dist(centroids_sites_sub$centroids))
  var_results$day_pair[i] <- mean(dist(centroids_days_sub$centroids))
  var_results$hour_pair[i] <- mean(dist(centroids_hours_sub$centroids))
  
  }
  
#-------------------------------------------------------------------------------#
#Kruskal-wallis test to determine if group means are significant

#first convert wide to long
disp_df <- var_results[,grepl("disp",colnames(var_results))] %>%  pivot_longer(everything(), names_to="group")
pair_df <- var_results[,grepl("pair",colnames(var_results))] %>%  pivot_longer(everything(), names_to="group")

#now kw test
kw_disp <- kruskal.test(value ~ group, data = disp_df) #significant
kw_pair <- kruskal.test(value ~ group, data = pair_df) #significant

#now dunn test to determine which groups are different from each other
dunnTest(value ~ as.factor(group),
         data=disp_df,
         method="bonferroni")

dunnTest(value ~ as.factor(group),
         data=pair_df,
         method="bonferroni")

#sites, days, and hours are all different from each other!

ggboxplot(disp_df, x = "group", y = "value", 
          color = "group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("site_disp", "day_disp", "hour_disp"),
          ylab = "Dispersion", xlab = "Group")

ggboxplot(pair_df, x = "group", y = "value", 
          color = "group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("site_pair", "day_pair", "hour_pair"),
          ylab = "Pairwise", xlab = "Group")


#calculate the range of all dispersion and pairwise values
range(disp_df$value[disp_df$group=="site_disp"])[2] - range(disp_df$value[disp_df$group=="site_disp"])[1]
range(disp_df$value[disp_df$group=="day_disp"])[2] - range(disp_df$value[disp_df$group=="day_disp"])[1]
range(disp_df$value[disp_df$group=="hour_disp"])[2] - range(disp_df$value[disp_df$group=="hour_disp"])[1]

range(pair_df$value[pair_df$group=="site_pair"])[2] - range(pair_df$value[pair_df$group=="site_pair"])[1]
range(pair_df$value[pair_df$group=="day_pair"])[2] - range(pair_df$value[pair_df$group=="day_pair"])[1]
range(pair_df$value[pair_df$group=="hour_pair"])[2] - range(pair_df$value[pair_df$group=="hour_pair"])[1]

#create table for kw test results
kw_results <- data.frame("Variability" = c(" ", "Dispersion", " ", " ", "Pairwise", " "),
                         "Scale" = rep(c("Sites", "Sampling campaigns", "Hours of the day"),2),
                         "n" = c(rep(500,6)),
                         "mean" = c(round(mean(var_results$site_disp),2), round(mean(var_results$day_disp),2),
                                    round(mean(var_results$hour_disp),2), round(mean(var_results$site_pair),2),
                                    round(mean(var_results$day_pair),2), round(mean(var_results$hour_pair),2)),
                         "sd" = c(round(sd(var_results$site_disp),2), round(sd(var_results$day_disp),2),
                                  round(sd(var_results$hour_disp),2), round(sd(var_results$site_pair),2),
                                  round(sd(var_results$day_pair),2), round(sd(var_results$hour_pair),2)),
                         "df" = c(" ",kw_disp$parameter, " ", " ", kw_pair$parameter, " "),
                         "χ2" = c(" ",round(kw_disp$statistic,3), " ", " ", round(kw_pair$statistic,3), " "),
                         "p-value" = c(" ",kw_disp$p.value, " ", " ", kw_pair$p.value, " "))
#write.csv(kw_results, file.path(getwd(),"/Summer2021-DataAnalysis/SummaryStats/Euclidean_distances_bootstrapped_kw_results.csv"),row.names = FALSE)

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
kw_sites <- kruskal.test(dist ~ group, data = within_site_dist) #sig! - so pelagic and littoral are different
kw_days <- kruskal.test(dist ~ group, data = within_day_dist) #sig
kw_hours <- kruskal.test(dist ~ group, data = within_hour_dist) #nope

#dunn post-hoc test
dunn_within_days <- dunnTest(dist ~ as.factor(group),
                        data=within_day_dist,
                        method="bonferroni")

cldList(P.adj ~ Comparison, data=dunn_within_days$res, threshold = 0.05)

ggboxplot(within_day_dist, x = "group", y = "dist", 
          color = "group", palette = c("#008585", "#9BBAA0", "#F2E2B0","#DEA868","#C7522B"),
          order = c("day1", "day2", "day3","day4","day5"),
          ylab = "average dispersion", xlab = "Sampling Camapigns")
#day 3 has greatest dispersion!

#create table for within scale kw test results
kw_results_disp <- data.frame("Group" = c("Littoral", "Pelagic", "10-11 Jul 2019", "24-25 Jul 2019", 
                                     "12-13 Aug 2020", "15-16 Jun 2021", "7-8 Jul 2021", "12pm",
                                     "6pm", "7pm", "8pm", "9pm", "12am", "4am", "5am", "6am", "7am", "12pm"),
                         "n" = c(55, 55, 22, 22, 22, 22, 22, rep(10,11)),
                         "mean" = c(round(mean(within_site_dist$dist[within_site_dist$group=="lit"]),2), 
                                    round(mean(within_site_dist$dist[within_site_dist$group=="pel"]),2),
                                    round(mean(within_day_dist$dist[within_day_dist$group=="day1"]),2),
                                    round(mean(within_day_dist$dist[within_day_dist$group=="day2"]),2),
                                    round(mean(within_day_dist$dist[within_day_dist$group=="day3"]),2),
                                    round(mean(within_day_dist$dist[within_day_dist$group=="day4"]),2),
                                    round(mean(within_day_dist$dist[within_day_dist$group=="day5"]),2),
                                    round(mean(within_hour_dist$dist[within_hour_dist$group=="hour1"]),2),
                                    round(mean(within_hour_dist$dist[within_hour_dist$group=="hour2"]),2),
                                    round(mean(within_hour_dist$dist[within_hour_dist$group=="hour3"]),2),
                                    round(mean(within_hour_dist$dist[within_hour_dist$group=="hour4"]),2),
                                    round(mean(within_hour_dist$dist[within_hour_dist$group=="hour5"]),2),
                                    round(mean(within_hour_dist$dist[within_hour_dist$group=="hour6"]),2),
                                    round(mean(within_hour_dist$dist[within_hour_dist$group=="hour7"]),2),
                                    round(mean(within_hour_dist$dist[within_hour_dist$group=="hour8"]),2),
                                    round(mean(within_hour_dist$dist[within_hour_dist$group=="hour9"]),2),
                                    round(mean(within_hour_dist$dist[within_hour_dist$group=="hour10"]),2),
                                    round(mean(within_hour_dist$dist[within_hour_dist$group=="hour11"]),2)),
                         "sd" = c(round(sd(within_site_dist$dist[within_site_dist$group=="lit"]),2), 
                                  round(sd(within_site_dist$dist[within_site_dist$group=="pel"]),2),
                                  round(sd(within_day_dist$dist[within_day_dist$group=="day1"]),2),
                                  round(sd(within_day_dist$dist[within_day_dist$group=="day2"]),2),
                                  round(sd(within_day_dist$dist[within_day_dist$group=="day3"]),2),
                                  round(sd(within_day_dist$dist[within_day_dist$group=="day4"]),2),
                                  round(sd(within_day_dist$dist[within_day_dist$group=="day5"]),2),
                                  round(sd(within_hour_dist$dist[within_hour_dist$group=="hour1"]),2),
                                  round(sd(within_hour_dist$dist[within_hour_dist$group=="hour2"]),2),
                                  round(sd(within_hour_dist$dist[within_hour_dist$group=="hour3"]),2),
                                  round(sd(within_hour_dist$dist[within_hour_dist$group=="hour4"]),2),
                                  round(sd(within_hour_dist$dist[within_hour_dist$group=="hour5"]),2),
                                  round(sd(within_hour_dist$dist[within_hour_dist$group=="hour6"]),2),
                                  round(sd(within_hour_dist$dist[within_hour_dist$group=="hour7"]),2),
                                  round(sd(within_hour_dist$dist[within_hour_dist$group=="hour8"]),2),
                                  round(sd(within_hour_dist$dist[within_hour_dist$group=="hour9"]),2),
                                  round(sd(within_hour_dist$dist[within_hour_dist$group=="hour10"]),2),
                                  round(sd(within_hour_dist$dist[within_hour_dist$group=="hour11"]),2)),
                         "df" = c(kw_sites$parameter, " ", rep(" ",2), kw_days$parameter, rep(" ",2),
                                  rep(" ", 5), kw_hours$parameter, rep(" ", 5)),
                         "χ2" = c(round(kw_sites$statistic,3), " ", rep(" ",2), round(kw_days$statistic,3), rep(" ",2),
                                  rep(" ", 5), round(kw_hours$statistic,3), rep(" ", 5)),
                         "p-value" = c(kw_sites$p.value, " ", rep(" ",2), kw_days$p.value, rep(" ",2),
                                       rep(" ", 5), kw_hours$p.value, rep(" ", 5)))
#write.csv(kw_results_disp, file.path(getwd(),"/Summer2021-DataAnalysis/SummaryStats/within_group_dispersion_kw_results.csv"),row.names = FALSE)


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

#------------------------------------------------------------------------------#
#pull in driver data --> can't use ysi or fp because only have 2 or 3 casts out of 5

chem <- read.csv('~/Documents/VirginiaTech/research/Reservoirs/Data/DataAlreadyUploadedToEDI/EDIProductionFiles/MakeEMLChemistry/2022/chemistry_2013_2022.csv') %>% 
  mutate(DateTime = as.Date(DateTime)) %>%
  filter(Reservoir =="BVR" & Site==50 & 
           DateTime %in% c(as.Date("2019-07-10"), as.Date("2019-07-11"), as.Date("2019-07-24"),
                           as.Date("2019-07-25"), as.Date("2020-08-12"), as.Date("2020-08-13"),
                           as.Date("2021-06-15"), as.Date("2021-06-16"), as.Date("2021-07-07"), 
                           as.Date("2021-07-08"))) #can only do TN/TP for all 5 MSNs

#add msn #
chem$msn <- ifelse(chem$DateTime=="2019-07-10" | chem$DateTime=="2019-07-11","1",
                        ifelse(chem$DateTime=="2019-07-24" | chem$DateTime=="2019-07-25","2",
                               ifelse(chem$DateTime=="2020-08-12" | chem$DateTime=="2020-08-13","3",
                                      ifelse(chem$DateTime=="2021-06-15" | chem$DateTime=="2021-06-16","4","5"))))

#average across msn and depth
chem_avg <- chem %>% group_by(Reservoir, Depth_m, msn) %>% summarise(across(everything(), list(mean)))


secchi <- read.csv('~/Documents/VirginiaTech/research/Reservoirs/Data/DataAlreadyUploadedToEDI/EDIProductionFiles/MakeEMLYSI_PAR_secchi/2022/Data/Secchi_depth_2013-2022.csv') %>% 
  mutate(DateTime = as.Date(DateTime)) %>%
  filter(Reservoir =="BVR" & Site==50 & 
           DateTime %in% c(as.Date("2019-07-10"), as.Date("2019-07-11"), as.Date("2019-07-24"),
                           as.Date("2019-07-25"), as.Date("2020-08-12"), as.Date("2020-08-13"),
                           as.Date("2021-06-15"), as.Date("2021-06-16"), as.Date("2021-07-07"), 
                           as.Date("2021-07-08"))) #can only do TN/TP for all 5 MSNs

#add secchi for 5th MSN - 1.8m
secchi[5,] <- c("BVR",50,"2021-07-08",1.8,0,0)

ctd<- read.csv('~/Documents/VirginiaTech/research/BVR_GLM/bvr_glm/field_data/CTD_final_2013_2022.csv') %>% 
  mutate(DateTime = as.Date(DateTime)) %>%
  filter(Reservoir =="BVR" & 
           DateTime %in% c(as.Date("2019-07-10"), as.Date("2019-07-11"), as.Date("2019-07-24"),
                           as.Date("2019-07-25"), as.Date("2020-08-12"), as.Date("2020-08-13"),
                           as.Date("2021-06-15"), as.Date("2021-06-16"), as.Date("2021-07-07"), 
                           as.Date("2021-07-08"), as.Date("2021-07-12"))) #ctd cast from 4 days after MSN 5

#select every 0.5m from casts
ctd_final <- ctd %>%
  dplyr::mutate(rdepth = plyr::round_any(Depth_m, 0.5)) %>% 
  dplyr::group_by(DateTime, rdepth, Reservoir, Site) %>%
  dplyr::summarise(value = mean(Temp_C)) %>% 
  dplyr::rename(depth = rdepth) 


#calculate thermocline depth by date
ctd_thermo_depth <- ctd_final %>% group_by(DateTime) %>% filter(depth > 0) %>% 
  mutate(therm_depth = thermo.depth(value,depth))

#add msn # 
ctd_thermo_depth$msn <- ifelse(ctd_thermo_depth$DateTime=="2019-07-10" | ctd_thermo_depth$DateTime=="2019-07-11","1",
                        ifelse(ctd_thermo_depth$DateTime=="2019-07-24" | ctd_thermo_depth$DateTime=="2019-07-25","2",
                               ifelse(ctd_thermo_depth$DateTime=="2020-08-12" | ctd_thermo_depth$DateTime=="2020-08-13","3",
                                      ifelse(ctd_thermo_depth$DateTime=="2021-06-15" | ctd_thermo_depth$DateTime=="2021-06-16","4","5"))))



depths = c(0.1, seq(1, 10, by = 1))
ctd.final<-data.frame() 

for (i in 1:length(depths)){
  ctd_layer <- ctd %>% group_by(DateTime) %>% slice(which.min(abs(as.numeric(Depth_m) - depths[i])))
  # Bind each of the data layers together.
  ctd.final = bind_rows(ctd.final, ctd_layer)
}

#round to one decimal place
ctd.final$Depth_m <- round(ctd.final$Depth_m,1)

#add msn # 
ctd.final$msn <- ifelse(ctd.final$DateTime=="2019-07-10" | ctd.final$DateTime=="2019-07-11","1",
                        ifelse(ctd.final$DateTime=="2019-07-24" | ctd.final$DateTime=="2019-07-25","2",
                               ifelse(ctd.final$DateTime=="2020-08-12" | ctd.final$DateTime=="2020-08-13","3",
                                      ifelse(ctd.final$DateTime=="2021-06-15" | ctd.final$DateTime=="2021-06-16","4","5"))))

#average across msn and depth
ctd.final_avg <- ctd.final %>% group_by(Reservoir, Depth_m, msn) %>% summarise(across(everything(), list(mean)))

#make an environmental driver df
msn_drivers <- data.frame("groups" = as.character(1:5), #epi is 0.1m, hypo is 10m - avg for both noons of msn when available
                          "epi_temp" = c(ctd.final_avg$Temp_C_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==1],
                                         ctd.final_avg$Temp_C_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==2],
                                         ctd.final_avg$Temp_C_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==3],
                                         ctd.final_avg$Temp_C_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==4],
                                         ctd.final_avg$Temp_C_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==5]), 
                          "hypo_temp" = c(ctd.final_avg$Temp_C_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==1],
                                          ctd.final_avg$Temp_C_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==2],
                                          ctd.final_avg$Temp_C_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==3],
                                          ctd.final_avg$Temp_C_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==4],
                                          ctd.final_avg$Temp_C_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==5]),
                          "epi_DO" = c(ctd.final_avg$DO_mgL_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==1],
                                       ctd.final_avg$DO_mgL_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==2],
                                       ctd.final_avg$DO_mgL_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==3],
                                       ctd.final_avg$DO_mgL_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==4],
                                       ctd.final_avg$DO_mgL_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==5]), 
                          "hypo_DO" = c(ctd.final_avg$DO_mgL_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==1],
                                        ctd.final_avg$DO_mgL_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==2],
                                        ctd.final_avg$DO_mgL_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==3],
                                        ctd.final_avg$DO_mgL_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==4],
                                        ctd.final_avg$DO_mgL_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==5]),
                          "epi_DOsat" = c(ctd.final_avg$DOsat_percent_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==1],
                                          ctd.final_avg$DOsat_percent_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==2],
                                          ctd.final_avg$DOsat_percent_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==3],
                                          ctd.final_avg$DOsat_percent_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==4],
                                          ctd.final_avg$DOsat_percent_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==5]),
                          "hypo_DOsat" = c(ctd.final_avg$DOsat_percent_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==1],
                                           ctd.final_avg$DOsat_percent_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==2],
                                           ctd.final_avg$DOsat_percent_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==3],
                                           ctd.final_avg$DOsat_percent_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==4],
                                           ctd.final_avg$DOsat_percent_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==5]),
                          "epi_spcond" = c(ctd.final_avg$SpCond_uScm_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==1],
                                           ctd.final_avg$SpCond_uScm_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==2],
                                           ctd.final_avg$SpCond_uScm_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==3],
                                           ctd.final_avg$SpCond_uScm_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==4],
                                           ctd.final_avg$SpCond_uScm_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==5]),
                          "hypo_spcond" = c(ctd.final_avg$SpCond_uScm_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==1],
                                            ctd.final_avg$SpCond_uScm_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==2],
                                            ctd.final_avg$SpCond_uScm_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==3],
                                            ctd.final_avg$SpCond_uScm_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==4],
                                            ctd.final_avg$SpCond_uScm_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==5]),
                          "epi_chl" = c(ctd.final_avg$Chla_ugL_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==1],
                                        ctd.final_avg$Chla_ugL_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==2],
                                        ctd.final_avg$Chla_ugL_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==3],
                                        ctd.final_avg$Chla_ugL_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==4],
                                        ctd.final_avg$Chla_ugL_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==5]),
                          "hypo_chl" = c(ctd.final_avg$Chla_ugL_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==1],
                                         ctd.final_avg$Chla_ugL_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==2],
                                         ctd.final_avg$Chla_ugL_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==3],
                                         ctd.final_avg$Chla_ugL_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==4],
                                         ctd.final_avg$Chla_ugL_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==5]),
                          "epi_par" = c(ctd.final_avg$PAR_umolm2s_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==1],
                                        ctd.final_avg$PAR_umolm2s_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==2],
                                        ctd.final_avg$PAR_umolm2s_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==3],
                                        ctd.final_avg$PAR_umolm2s_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==4],
                                        ctd.final_avg$PAR_umolm2s_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==5]),
                          "epi_tn" = c(chem_avg$TN_ugL_1[chem_avg$Depth_m==0.1 & chem_avg$msn==1],
                                       chem_avg$TN_ugL_1[chem_avg$Depth_m==0.1 & chem_avg$msn==2],
                                       chem_avg$TN_ugL_1[chem_avg$Depth_m==0.1 & chem_avg$msn==3],
                                       chem_avg$TN_ugL_1[chem_avg$Depth_m==0.1 & chem_avg$msn==4],
                                       chem_avg$TN_ugL_1[chem_avg$Depth_m==0.1 & chem_avg$msn==5]), 
                          "epi_tp" = c(chem_avg$TP_ugL_1[chem_avg$Depth_m==0.1 & chem_avg$msn==1],
                                       chem_avg$TP_ugL_1[chem_avg$Depth_m==0.1 & chem_avg$msn==2],
                                       chem_avg$TP_ugL_1[chem_avg$Depth_m==0.1 & chem_avg$msn==3],
                                       chem_avg$TP_ugL_1[chem_avg$Depth_m==0.1 & chem_avg$msn==4],
                                       chem_avg$TP_ugL_1[chem_avg$Depth_m==0.1 & chem_avg$msn==5]), 
                          "secchi" = secchi$Secchi_m,
                          "thermo_depth" = c(mean(ctd_thermo_depth$therm_depth[ctd_thermo_depth$msn==1]),
                                            mean(ctd_thermo_depth$therm_depth[ctd_thermo_depth$msn==2]),
                                            mean(ctd_thermo_depth$therm_depth[ctd_thermo_depth$msn==3]),
                                            mean(ctd_thermo_depth$therm_depth[ctd_thermo_depth$msn==4]),
                                            mean(ctd_thermo_depth$therm_depth[ctd_thermo_depth$msn==5])))

#join driver and env dfs
zoops_plus_drivers <- left_join(zoop_avg, msn_drivers)

#fit environmental drivers onto ordination
fit_env <- envfit(ord$sites, zoops_plus_drivers[,c(15:29)])

#pull out vectors
scores <- as.data.frame(scores(fit_env, display = "vectors"))
scores <- cbind(scores, env = rownames(scores))

#plot drivers w/ NMDS
NMDS_day_env <- days$plot + geom_point() + theme_bw() + geom_path() +
  geom_polygon(data = days$df_hull, aes(x = x, y = y, fill = Group), alpha=0.2) +
  geom_point(data=days$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2, 
             fill=c("#008585","#89B199","#EFECBF","#DB9B5A","#C7522B")) +
  theme(text = element_text(size=7), axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), 
        legend.key.height=unit(0.3,"line"), 
        legend.key = element_blank(),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-0,-0,-0,-0),
        legend.direction = "horizontal",
        axis.text.x = element_text(vjust = 0.5), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "right", legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(0,-0.1,0,-0.1), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) + 
  guides(fill="none", color = guide_legend(ncol=1)) +
  scale_fill_manual("",values=c("#008585","#89B199","#EFECBF","#DB9B5A","#C7522B"))+
  scale_color_manual("",values=c("#008585","#89B199","#EFECBF","#DB9B5A","#C7522B"),
                     label=c('10-11 Jul 2019', '24-25 Jul 2019','12-13 Aug 2020',
                             '15-16 Jun 2021', '7-8 Jul 2021')) +
  geom_segment(data = scores,
              aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), linewidth= 0.3,
              arrow = arrow(length = unit(0.1, "cm")), colour = "black") +
  geom_text(data = scores, aes(x = NMDS1, y = NMDS2, label = env),
            size = 2)
ggsave(file.path(getwd(),"Summer2021-DataAnalysis/Figures/NMDS_2v1_envfit.jpg"),
              NMDS_day_env, width=3.5, height=4) 

#------------------------------------------------------------------------------#
#plot env variable vs NMDS2 centroids

plot(msn_drivers$epi_DO ~ days$df_mean.ord$y, pch=19, cex=2,
     col=c("#008585","#89B199","#EFECBF","#DB9B5A","#C7522B"), 
     xlab="NMDS2", ylab="Epilimnetic DO (mgL)")
