#script to look at space vs time euclidean distances 

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
#first average times for each 24-hour campaign so there are 10 points per day (basically just averaging noon and midnight)
zoop_epi_tows$order <- ifelse(zoop_epi_tows$Hour=="11" | zoop_epi_tows$Hour=="12",1, ifelse(zoop_epi_tows$Hour=="18",2,ifelse(zoop_epi_tows$Hour=="19",3,
                      ifelse(zoop_epi_tows$Hour=="20",4,ifelse(zoop_epi_tows$Hour=="21",5,
                      ifelse(zoop_epi_tows$Hour=="0:" | zoop_epi_tows$Hour=="23",6,
                      ifelse(zoop_epi_tows$Hour=="4:" | zoop_epi_tows$Hour=="3:",7,ifelse(zoop_epi_tows$Hour=="5:",8,
                      ifelse(zoop_epi_tows$Hour=="6:",9,10)))))))))

#add order 11 for noon2
zoop_epi_tows$order[zoop_epi_tows$order==1 & (zoop_epi_tows$collect_date=="2019-07-10" | zoop_epi_tows$collect_date=="2019-07-24" |
                      zoop_epi_tows$collect_date=="2020-08-12" | zoop_epi_tows$collect_date=="2021-06-15" |
                      zoop_epi_tows$collect_date=="2021-07-07")] <- 11

#now specify whether it is noon1 or noon2
zoop_epi_tows$time[zoop_epi_tows$order==1] <- "noon"
zoop_epi_tows$time[zoop_epi_tows$order==11] <- "noon2"

#average by groups, site, then order
zoop_avg <- zoop_epi_tows %>% group_by(groups,site,order) %>%
  summarise_at(vars(Calanoida_density_NopL_1:Conochilidae_density_NopL_1), list(mean = mean))

#pelagic vs littoral dfs
zoop_pel <- zoop_avg[zoop_avg$site=="pel",]
zoop_lit <- zoop_avg[zoop_avg$site=="lit",]

#10 dfs for time (sunrise 1-4, noon, sunset 1-4, midnight)
zoop_noon1 <- zoop_avg[zoop_avg$order==1,]
zoop_sunset1 <- zoop_avg[zoop_avg$order==2,]
zoop_sunset2 <- zoop_avg[zoop_avg$order==3,]
zoop_sunset3 <- zoop_avg[zoop_avg$order==4,]
zoop_sunset4 <- zoop_avg[zoop_avg$order==5,]
zoop_midnight <- zoop_avg[zoop_avg$order==6,]
zoop_sunrise1 <- zoop_avg[zoop_avg$order==7,]
zoop_sunrise2 <- zoop_avg[zoop_avg$order==8,]
zoop_sunrise3 <- zoop_avg[zoop_avg$order==9,]
zoop_sunrise4 <- zoop_avg[zoop_avg$order==10,]
zoop_noon2 <- zoop_avg[zoop_avg$order==11,]

#5 dfs for MSN #
zoop_msn1 <- zoop_avg[zoop_avg$groups==1,]
zoop_msn2 <- zoop_avg[zoop_avg$groups==2,]
zoop_msn3 <- zoop_avg[zoop_avg$groups==3,]
zoop_msn4 <- zoop_avg[zoop_avg$groups==4,]
zoop_msn5 <- zoop_avg[zoop_avg$groups==5,]

#only select data cols
zoop_temporal_avg_dens <- zoop_avg[,c(grepl("mean",colnames(zoop_avg)))] 
zoop_pel_dens <- zoop_pel[,c(grepl("mean",colnames(zoop_pel)))]
zoop_lit_dens <- zoop_lit[,c(grepl("mean",colnames(zoop_lit)))]
zoop_sunrise1_dens <- zoop_sunrise1[,c(grepl("mean",colnames(zoop_sunrise1)))]
zoop_sunrise2_dens <- zoop_sunrise2[,c(grepl("mean",colnames(zoop_sunrise2)))]
zoop_sunrise3_dens <- zoop_sunrise3[,c(grepl("mean",colnames(zoop_sunrise3)))]
zoop_sunrise4_dens <- zoop_sunrise4[,c(grepl("mean",colnames(zoop_sunrise4)))]
zoop_noon1_dens <- zoop_noon1[,c(grepl("mean",colnames(zoop_noon1)))]
zoop_noon2_dens <- zoop_noon2[,c(grepl("mean",colnames(zoop_noon2)))]
zoop_sunset1_dens <- zoop_sunset1[,c(grepl("mean",colnames(zoop_sunset1)))]
zoop_sunset2_dens <- zoop_sunset2[,c(grepl("mean",colnames(zoop_sunset2)))]
zoop_sunset3_dens <- zoop_sunset3[,c(grepl("mean",colnames(zoop_sunset3)))]
zoop_sunset4_dens <- zoop_sunset4[,c(grepl("mean",colnames(zoop_sunset4)))]
zoop_midnight_dens <- zoop_midnight[,c(grepl("mean",colnames(zoop_midnight)))]
zoop_msn1_dens <- zoop_msn1[,c(grepl("mean",colnames(zoop_msn1)))]
zoop_msn2_dens <- zoop_msn2[,c(grepl("mean",colnames(zoop_msn2)))]
zoop_msn3_dens <- zoop_msn3[,c(grepl("mean",colnames(zoop_msn3)))]
zoop_msn4_dens <- zoop_msn4[,c(grepl("mean",colnames(zoop_msn4)))]
zoop_msn5_dens <- zoop_msn5[,c(grepl("mean",colnames(zoop_msn5)))]

#transforming data - hellinger transformation because gives low weight to low/zero values
#converts species abundances from absolute to relative - use w/ bray curtis
zoop_temporal_dens_avg_trans <- hellinger(zoop_temporal_avg_dens)
zoop_pel_dens_trans <- hellinger(zoop_pel_dens)
zoop_lit_dens_trans <- hellinger(zoop_lit_dens)
zoop_sunrise1_dens_trans <- hellinger(zoop_sunrise1_dens)
zoop_sunrise2_dens_trans <- hellinger(zoop_sunrise2_dens)
zoop_sunrise3_dens_trans <- hellinger(zoop_sunrise3_dens)
zoop_sunrise4_dens_trans <- hellinger(zoop_sunrise4_dens)
zoop_noon1_dens_trans <- hellinger(zoop_noon1_dens)
zoop_noon2_dens_trans <- hellinger(zoop_noon2_dens)
zoop_sunset1_dens_trans <- hellinger(zoop_sunset1_dens)
zoop_sunset2_dens_trans <- hellinger(zoop_sunset2_dens)
zoop_sunset3_dens_trans <- hellinger(zoop_sunset3_dens)
zoop_sunset4_dens_trans <- hellinger(zoop_sunset4_dens)
zoop_midnight_dens_trans <- hellinger(zoop_midnight_dens)
zoop_msn1_dens_trans <- hellinger(zoop_msn1_dens)
zoop_msn2_dens_trans <- hellinger(zoop_msn2_dens)
zoop_msn3_dens_trans <- hellinger(zoop_msn3_dens)
zoop_msn4_dens_trans <- hellinger(zoop_msn4_dens)
zoop_msn5_dens_trans <- hellinger(zoop_msn5_dens)

#-------------------------------------------------------------------------------#
#              NMDS on each datafram (time, day, and site)                      #
#-------------------------------------------------------------------------------#
#scree plot to choose dimension #
dimcheckMDS(zoop_temporal_dens_avg_trans, distance = "bray", k = 6, trymax = 20, autotransform = TRUE)

#sites NMDS
NMDS_pel_bray <- metaMDS(zoop_pel_dens_trans, distance='bray', k=4, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_pel_bray$stress

NMDS_lit_bray <- metaMDS(zoop_lit_dens_trans, distance='bray', k=4, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_lit_bray$stress

#time NMDS (doing k=3 for these because fewer data points and stress is similar as pelagic NMDS)
NMDS_sunrise1_bray <- metaMDS(zoop_sunrise1_dens_trans, distance='bray', k=3, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_sunrise1_bray$stress

NMDS_sunrise2_bray <- metaMDS(zoop_sunrise2_dens_trans, distance='bray', k=3, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_sunrise2_bray$stress

NMDS_sunrise3_bray <- metaMDS(zoop_sunrise3_dens_trans, distance='bray', k=3, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_sunrise3_bray$stress

NMDS_sunrise4_bray <- metaMDS(zoop_sunrise4_dens_trans, distance='bray', k=3, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_sunrise4_bray$stress

NMDS_noon1_bray <- metaMDS(zoop_noon1_dens_trans, distance='bray', k=3, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_noon1_bray$stress

NMDS_noon2_bray <- metaMDS(zoop_noon2_dens_trans, distance='bray', k=3, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_noon2_bray$stress

NMDS_sunset1_bray <- metaMDS(zoop_sunset1_dens_trans, distance='bray', k=3, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_sunset1_bray$stress

NMDS_sunset2_bray <- metaMDS(zoop_sunset2_dens_trans, distance='bray', k=3, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_sunset2_bray$stress

NMDS_sunset3_bray <- metaMDS(zoop_sunset3_dens_trans, distance='bray', k=3, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_sunset3_bray$stress

NMDS_sunset4_bray <- metaMDS(zoop_sunset4_dens_trans, distance='bray', k=3, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_sunset4_bray$stress

NMDS_midnight_bray <- metaMDS(zoop_midnight_dens_trans, distance='bray', k=3, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_midnight_bray$stress

#day NMDS (back to k=4)
NMDS_msn1_bray <- metaMDS(zoop_msn1_dens_trans, distance='bray', k=4, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_msn1_bray$stress

NMDS_msn2_bray <- metaMDS(zoop_msn2_dens_trans, distance='bray', k=4, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_msn2_bray$stress

NMDS_msn3_bray <- metaMDS(zoop_msn3_dens_trans, distance='bray', k=4, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_msn3_bray$stress

NMDS_msn4_bray <- metaMDS(zoop_msn4_dens_trans, distance='bray', k=4, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_msn4_bray$stress

NMDS_msn5_bray <- metaMDS(zoop_msn5_dens_trans, distance='bray', k=4, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_msn5_bray$stress


#-------------------------------------------------------------------------------#
#                     Calculating euclidean distance                            #
#-------------------------------------------------------------------------------#
#so I'm technically calculating the euclidean distances from bray-curtis distance matrices

#step 1: take NMDS output for each df using NMDS coordinates
zoop_pel_euc <- as.matrix(vegdist(NMDS_pel_bray$points, method='euclidean')) 
zoop_lit_euc <- as.matrix(vegdist(NMDS_lit_bray$points, method='euclidean'))

zoop_noon1_euc <- as.matrix(vegdist(NMDS_noon1_bray$points, method='euclidean'))
zoop_sunset1_euc <- as.matrix(vegdist(NMDS_sunset1_bray$points, method='euclidean'))
zoop_sunset2_euc <- as.matrix(vegdist(NMDS_sunset2_bray$points, method='euclidean'))
zoop_sunset3_euc <- as.matrix(vegdist(NMDS_sunset3_bray$points, method='euclidean'))
zoop_sunset4_euc <- as.matrix(vegdist(NMDS_sunset4_bray$points, method='euclidean'))
zoop_midnight_euc <- as.matrix(vegdist(NMDS_midnight_bray$points, method='euclidean'))
zoop_sunrise1_euc <- as.matrix(vegdist(NMDS_sunrise1_bray$points, method='euclidean'))
zoop_sunrise2_euc <- as.matrix(vegdist(NMDS_sunrise2_bray$points, method='euclidean'))
zoop_sunrise3_euc <- as.matrix(vegdist(NMDS_sunrise3_bray$points, method='euclidean'))
zoop_sunrise4_euc <- as.matrix(vegdist(NMDS_sunrise4_bray$points, method='euclidean'))
zoop_noon2_euc <- as.matrix(vegdist(NMDS_noon2_bray$points, method='euclidean'))

zoop_msn1_euc <- as.matrix(vegdist(NMDS_msn1_bray$points, method='euclidean'))
zoop_msn2_euc <- as.matrix(vegdist(NMDS_msn2_bray$points, method='euclidean'))
zoop_msn3_euc <- as.matrix(vegdist(NMDS_msn3_bray$points, method='euclidean'))
zoop_msn4_euc <- as.matrix(vegdist(NMDS_msn4_bray$points, method='euclidean'))
zoop_msn5_euc <- as.matrix(vegdist(NMDS_msn5_bray$points, method='euclidean'))


#step 2: site, time of day, year distances
pel_dist <- sum(zoop_pel_euc[1,2],zoop_pel_euc[2,3],zoop_pel_euc[3,4],zoop_pel_euc[4,5],zoop_pel_euc[5,6],
                zoop_pel_euc[6,7],zoop_pel_euc[7,8],zoop_pel_euc[8,9],zoop_pel_euc[9,10],zoop_pel_euc[10,11],
                zoop_pel_euc[11,12],zoop_pel_euc[12,13],zoop_pel_euc[13,14],zoop_pel_euc[14,15],zoop_pel_euc[15,16],
                zoop_pel_euc[16,17],zoop_pel_euc[17,18],zoop_pel_euc[18,19],zoop_pel_euc[19,20],zoop_pel_euc[20,21],
                zoop_pel_euc[21,22],zoop_pel_euc[22,23],zoop_pel_euc[23,24],zoop_pel_euc[24,25],zoop_pel_euc[25,26],
                zoop_pel_euc[26,27],zoop_pel_euc[27,28],zoop_pel_euc[28,29],zoop_pel_euc[29,30],zoop_pel_euc[30,31],
                zoop_pel_euc[31,32],zoop_pel_euc[32,33],zoop_pel_euc[33,34],zoop_pel_euc[34,35],zoop_pel_euc[35,36],
                zoop_pel_euc[36,37],zoop_pel_euc[37,38],zoop_pel_euc[38,39],zoop_pel_euc[39,40],zoop_pel_euc[40,41],
                zoop_pel_euc[41,42],zoop_pel_euc[42,43],zoop_pel_euc[43,44],zoop_pel_euc[44,45],zoop_pel_euc[45,46],
                zoop_pel_euc[46,47],zoop_pel_euc[47,48],zoop_pel_euc[48,49],zoop_pel_euc[49,50],zoop_pel_euc[50,51],
                zoop_pel_euc[51,52],zoop_pel_euc[52,53],zoop_pel_euc[53,54],zoop_pel_euc[54,55])

lit_dist <- sum(zoop_lit_euc[1,2],zoop_lit_euc[2,3],zoop_lit_euc[3,4],zoop_lit_euc[4,5],zoop_lit_euc[5,6],
                zoop_lit_euc[6,7],zoop_lit_euc[7,8],zoop_lit_euc[8,9],zoop_lit_euc[9,10],zoop_lit_euc[10,11],
                zoop_lit_euc[11,12],zoop_lit_euc[12,13],zoop_lit_euc[13,14],zoop_lit_euc[14,15],zoop_lit_euc[15,16],
                zoop_lit_euc[16,17],zoop_lit_euc[17,18],zoop_lit_euc[18,19],zoop_lit_euc[19,20],zoop_lit_euc[20,21],
                zoop_lit_euc[21,22],zoop_lit_euc[22,23],zoop_lit_euc[23,24],zoop_lit_euc[24,25],zoop_lit_euc[25,26],
                zoop_lit_euc[26,27],zoop_lit_euc[27,28],zoop_lit_euc[28,29],zoop_lit_euc[29,30],zoop_lit_euc[30,31],
                zoop_lit_euc[31,32],zoop_lit_euc[32,33],zoop_lit_euc[33,34],zoop_lit_euc[34,35],zoop_lit_euc[35,36],
                zoop_lit_euc[36,37],zoop_lit_euc[37,38],zoop_lit_euc[38,39],zoop_lit_euc[39,40],zoop_lit_euc[40,41],
                zoop_lit_euc[41,42],zoop_lit_euc[42,43],zoop_lit_euc[43,44],zoop_lit_euc[44,45],zoop_lit_euc[45,46],
                zoop_lit_euc[46,47],zoop_lit_euc[47,48],zoop_lit_euc[48,49],zoop_lit_euc[49,50],zoop_lit_euc[50,51],
                zoop_lit_euc[51,52],zoop_lit_euc[52,53],zoop_lit_euc[53,54],zoop_lit_euc[54,55]) 

sunrise1_dist <- sum(zoop_sunrise1_euc[1,2],zoop_sunrise1_euc[2,3],zoop_sunrise1_euc[3,4],zoop_sunrise1_euc[4,5],
                     zoop_sunrise1_euc[5,6],zoop_sunrise1_euc[6,7],zoop_sunrise1_euc[7,8],zoop_sunrise1_euc[8,9],zoop_sunrise1_euc[9,10])

sunrise2_dist <- sum(zoop_sunrise2_euc[1,2],zoop_sunrise2_euc[2,3],zoop_sunrise2_euc[3,4],zoop_sunrise2_euc[4,5],
                     zoop_sunrise2_euc[5,6],zoop_sunrise2_euc[6,7],zoop_sunrise2_euc[7,8],zoop_sunrise2_euc[8,9],zoop_sunrise2_euc[9,10])

sunrise3_dist <- sum(zoop_sunrise3_euc[1,2],zoop_sunrise3_euc[2,3],zoop_sunrise3_euc[3,4],zoop_sunrise3_euc[4,5],
                     zoop_sunrise3_euc[5,6],zoop_sunrise3_euc[6,7],zoop_sunrise3_euc[7,8],zoop_sunrise3_euc[8,9],zoop_sunrise3_euc[9,10])

sunrise4_dist <- sum(zoop_sunrise4_euc[1,2],zoop_sunrise4_euc[2,3],zoop_sunrise4_euc[3,4],zoop_sunrise4_euc[4,5],
                     zoop_sunrise4_euc[5,6],zoop_sunrise4_euc[6,7],zoop_sunrise4_euc[7,8],zoop_sunrise4_euc[8,9],zoop_sunrise4_euc[9,10])

noon1_dist <- sum(zoop_noon1_euc[1,2],zoop_noon1_euc[2,3],zoop_noon1_euc[3,4],zoop_noon1_euc[4,5],
                  zoop_noon1_euc[5,6],zoop_noon1_euc[6,7],zoop_noon1_euc[7,8],zoop_noon1_euc[8,9],zoop_noon1_euc[9,10])

noon2_dist <- sum(zoop_noon2_euc[1,2],zoop_noon2_euc[2,3],zoop_noon2_euc[3,4],zoop_noon2_euc[4,5],
                  zoop_noon2_euc[5,6],zoop_noon2_euc[6,7],zoop_noon2_euc[7,8],zoop_noon2_euc[8,9],zoop_noon2_euc[9,10])

sunset1_dist <- sum(zoop_sunset1_euc[1,2],zoop_sunset1_euc[2,3],zoop_sunset1_euc[3,4],zoop_sunset1_euc[4,5],
                    zoop_sunset1_euc[5,6],zoop_sunset1_euc[6,7],zoop_sunset1_euc[7,8],zoop_sunset1_euc[8,9],zoop_sunset1_euc[9,10])

sunset2_dist <- sum(zoop_sunset2_euc[1,2],zoop_sunset2_euc[2,3],zoop_sunset2_euc[3,4],zoop_sunset2_euc[4,5],
                    zoop_sunset2_euc[5,6],zoop_sunset2_euc[6,7],zoop_sunset2_euc[7,8],zoop_sunset2_euc[8,9],zoop_sunset2_euc[9,10])

sunset3_dist <- sum(zoop_sunset3_euc[1,2],zoop_sunset3_euc[2,3],zoop_sunset3_euc[3,4],zoop_sunset3_euc[4,5],
                    zoop_sunset3_euc[5,6],zoop_sunset3_euc[6,7],zoop_sunset3_euc[7,8],zoop_sunset3_euc[8,9],zoop_sunset3_euc[9,10])

sunset4_dist <- sum(zoop_sunset4_euc[1,2],zoop_sunset4_euc[2,3],zoop_sunset4_euc[3,4],zoop_sunset4_euc[4,5],
                    zoop_sunset4_euc[5,6],zoop_sunset4_euc[6,7],zoop_sunset4_euc[7,8],zoop_sunset4_euc[8,9],zoop_sunset4_euc[9,10])

midnight_dist <- sum(zoop_midnight_euc[1,2],zoop_midnight_euc[2,3],zoop_midnight_euc[3,4],zoop_midnight_euc[4,5],
                     zoop_midnight_euc[5,6],zoop_midnight_euc[6,7],zoop_midnight_euc[7,8],zoop_midnight_euc[8,9],zoop_midnight_euc[9,10])

msn1_dist <- sum(zoop_msn1_euc[1,2],zoop_msn1_euc[2,3],zoop_msn1_euc[3,4],zoop_msn1_euc[4,5],
                 zoop_msn1_euc[5,6],zoop_msn1_euc[6,7],zoop_msn1_euc[7,8],zoop_msn1_euc[8,9],
                 zoop_msn1_euc[9,10],zoop_msn1_euc[10,11],zoop_msn1_euc[11,12],zoop_msn1_euc[12,13],
                 zoop_msn1_euc[13,14],zoop_msn1_euc[14,15],zoop_msn1_euc[15,16],zoop_msn1_euc[16,17],
                 zoop_msn1_euc[17,18],zoop_msn1_euc[18,19],zoop_msn1_euc[19,20],zoop_msn1_euc[20,21],
                 zoop_msn1_euc[21,22])

msn2_dist <- sum(zoop_msn2_euc[1,2],zoop_msn2_euc[2,3],zoop_msn2_euc[3,4],zoop_msn2_euc[4,5],
                 zoop_msn2_euc[5,6],zoop_msn2_euc[6,7],zoop_msn2_euc[7,8],zoop_msn2_euc[8,9],
                 zoop_msn2_euc[9,10],zoop_msn2_euc[10,11],zoop_msn2_euc[11,12],zoop_msn2_euc[12,13],
                 zoop_msn2_euc[13,14],zoop_msn2_euc[14,15],zoop_msn2_euc[15,16],zoop_msn2_euc[16,17],
                 zoop_msn2_euc[17,18],zoop_msn2_euc[18,19],zoop_msn2_euc[19,20],zoop_msn2_euc[20,21],
                 zoop_msn2_euc[21,22])

msn3_dist <- sum(zoop_msn3_euc[1,2],zoop_msn3_euc[2,3],zoop_msn3_euc[3,4],zoop_msn3_euc[4,5],
                 zoop_msn3_euc[5,6],zoop_msn3_euc[6,7],zoop_msn3_euc[7,8],zoop_msn3_euc[8,9],
                 zoop_msn3_euc[9,10],zoop_msn3_euc[10,11],zoop_msn3_euc[11,12],zoop_msn3_euc[12,13],
                 zoop_msn3_euc[13,14],zoop_msn3_euc[14,15],zoop_msn3_euc[15,16],zoop_msn3_euc[16,17],
                 zoop_msn3_euc[17,18],zoop_msn3_euc[18,19],zoop_msn3_euc[19,20],zoop_msn3_euc[20,21],
                 zoop_msn3_euc[21,22])

msn4_dist <- sum(zoop_msn4_euc[1,2],zoop_msn4_euc[2,3],zoop_msn4_euc[3,4],zoop_msn4_euc[4,5],
                 zoop_msn4_euc[5,6],zoop_msn4_euc[6,7],zoop_msn4_euc[7,8],zoop_msn4_euc[8,9],
                 zoop_msn4_euc[9,10],zoop_msn4_euc[10,11],zoop_msn4_euc[11,12],zoop_msn4_euc[12,13],
                 zoop_msn4_euc[13,14],zoop_msn4_euc[14,15],zoop_msn4_euc[15,16],zoop_msn4_euc[16,17],
                 zoop_msn4_euc[17,18],zoop_msn4_euc[18,19],zoop_msn4_euc[19,20],zoop_msn4_euc[20,21],
                 zoop_msn4_euc[21,22])

msn5_dist <- sum(zoop_msn5_euc[1,2],zoop_msn5_euc[2,3],zoop_msn5_euc[3,4],zoop_msn5_euc[4,5],
                 zoop_msn5_euc[5,6],zoop_msn5_euc[6,7],zoop_msn5_euc[7,8],zoop_msn5_euc[8,9],
                 zoop_msn5_euc[9,10],zoop_msn5_euc[10,11],zoop_msn5_euc[11,12],zoop_msn5_euc[12,13],
                 zoop_msn5_euc[13,14],zoop_msn5_euc[14,15],zoop_msn5_euc[15,16],zoop_msn5_euc[16,17],
                 zoop_msn5_euc[17,18],zoop_msn5_euc[18,19],zoop_msn5_euc[19,20],zoop_msn5_euc[20,21],
                 zoop_msn5_euc[21,22])


#step 3: make a dataset of data
euc_distances_spacevtime_df <- data.frame(site=sum(pel_dist,lit_dist), 
                               time=sum(sunrise1_dist,sunrise2_dist,sunrise3_dist,sunrise4_dist,
                                        noon1_dist,sunset1_dist,sunset2_dist,sunset3_dist,sunset4_dist,
                                        midnight_dist,noon2_dist),
                               day=sum(msn1_dist,msn2_dist,msn3_dist,msn4_dist,msn5_dist))
#Take-home: hourly variability is greater than year to year and site variability

#SO, create a df for just time euclidean distances
time_distances <- data.frame(time= c("sunrise1","sunrise2","sunrise3","sunrise4","noon","noon2",
                                     "sunset1","sunset2","sunset3","sunset4","midnight"),
                               ED = c(sunrise1_dist, sunrise2_dist, sunrise3_dist, sunrise4_dist,
                                    noon1_dist, noon2_dist, sunset1_dist, sunset2_dist,
                                    sunset3_dist, sunset4_dist, midnight_dist))

#read in zoop hourly proportion in pel vs lit
total_prop<- read.csv(paste0(getwd(),"/Summer2021-DataAnalysis/SummaryStats/Hourly_proportions_pelvslit.csv"), header=T)

#get hourly just for total dens and biom  
hourly_prop <- total_prop %>%  
  filter(metric %in% c("ZoopDensity_No.pL","BiomassConcentration_ugpL")) %>%
  group_by(Hour) %>%
  summarise(prop_pel = mean(proportion_pel),
            prop_lit = mean(proportion_lit))

#get both dfs in same order
hourly_prop <- hourly_prop[match(time_distances$time,hourly_prop$Hour),]

#add in pelagic proportion of density and biomass to df
time_distances$pel_dens_prop <- hourly_prop$prop_pel
time_distances$lit_dens_prop <- hourly_prop$prop_lit

 #plot hourly euclidean distances vs pelagic total zoop proportion
#jpeg(file.path(getwd(),"Summer2021-DataAnalysis/Figures/2019-2021_HourlyEDvspel_proportion.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(time_distances$pel_dens_prop,time_distances$ED, xlab="Pel zoop proportion", 
     ylab="Hourly ED", cex=2.8, pch=19, cex.lab = 1.5,col="darkblue")
#points(time_distances$pel_dens_prop[1],time_distances$ED[1],bg=hcl.colors(11,"sunset")[1],pch=21,cex=3)
#points(time_distances$pel_dens_prop[2],time_distances$ED[2],bg=hcl.colors(11,"sunset")[2],pch=21,cex=3)
#points(time_distances$pel_dens_prop[3],time_distances$ED[3],bg=hcl.colors(11,"sunset")[3],pch=21,cex=3)
#points(time_distances$pel_dens_prop[4],time_distances$ED[4],bg=hcl.colors(11,"sunset")[4],pch=21,cex=3)
#points(time_distances$pel_dens_prop[5],time_distances$ED[5],bg=hcl.colors(11,"sunset")[5],pch=21,cex=3)
#points(time_distances$pel_dens_prop[6],time_distances$ED[6],bg=hcl.colors(11,"sunset")[6],pch=21,cex=3)
#points(time_distances$pel_dens_prop[7],time_distances$ED[7],bg=hcl.colors(11,"sunset")[7],pch=21,cex=3)
#points(time_distances$pel_dens_prop[8],time_distances$ED[8],bg=hcl.colors(11,"sunset")[8],pch=21,cex=3)
#points(time_distances$pel_dens_prop[9],time_distances$ED[9],bg=hcl.colors(11,"sunset")[9],pch=21,cex=3)
#points(time_distances$pel_dens_prop[10],time_distances$ED[10],bg=hcl.colors(11,"sunset")[10],pch=21,cex=3)
#points(time_distances$pel_dens_prop[11],time_distances$ED[11],bg=hcl.colors(11,"sunset")[11],pch=21,cex=3)
#legend("bottomright",legend=c("sunrise1","sunrise2","sunrise3","sunrise4",
#                              "noon1","noon2","sunset1","sunset2","sunset3","sunset4"),
#                        pch=21, pt.bg=hcl.colors(11,"sunset"),bty = "n",cex=1.4,ncol=2)
#dev.off()

#jpeg(file.path(getwd(),"Summer2021-DataAnalysis/Figures/2019-2021_HourlyEDvslit_proportion.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(time_distances$lit_dens_prop,time_distances$ED, xlab="Lit zoop proportion", 
     ylab="Hourly ED", cex=2.8, pch=19, cex.lab = 1.5, col="darkgreen")
#points(time_distances$lit_dens_prop[1],time_distances$ED[1],bg=hcl.colors(11,"sunset")[1],pch=21,cex=3)
#points(time_distances$lit_dens_prop[2],time_distances$ED[2],bg=hcl.colors(11,"sunset")[2],pch=21,cex=3)
#points(time_distances$lit_dens_prop[3],time_distances$ED[3],bg=hcl.colors(11,"sunset")[3],pch=21,cex=3)
#points(time_distances$lit_dens_prop[4],time_distances$ED[4],bg=hcl.colors(11,"sunset")[4],pch=21,cex=3)
#points(time_distances$lit_dens_prop[5],time_distances$ED[5],bg=hcl.colors(11,"sunset")[5],pch=21,cex=3)
#points(time_distances$lit_dens_prop[6],time_distances$ED[6],bg=hcl.colors(11,"sunset")[6],pch=21,cex=3)
#points(time_distances$lit_dens_prop[7],time_distances$ED[7],bg=hcl.colors(11,"sunset")[7],pch=21,cex=3)
#points(time_distances$lit_dens_prop[8],time_distances$ED[8],bg=hcl.colors(11,"sunset")[8],pch=21,cex=3)
#points(time_distances$lit_dens_prop[9],time_distances$ED[9],bg=hcl.colors(11,"sunset")[9],pch=21,cex=3)
#points(time_distances$lit_dens_prop[10],time_distances$ED[10],bg=hcl.colors(11,"sunset")[10],pch=21,cex=3)
#points(time_distances$lit_dens_prop[11],time_distances$ED[11],bg=hcl.colors(11,"sunset")[11],pch=21,cex=3)
#legend("bottomleft",legend=c("sunrise1","sunrise2","sunrise3","sunrise4",
#                              "noon1","noon2","sunset1","sunset2","sunset3","sunset4"),
#                       pch=21, pt.bg=hcl.colors(11,"sunset"),bty = "n",cex=1.4,ncol=2)
#dev.off()

