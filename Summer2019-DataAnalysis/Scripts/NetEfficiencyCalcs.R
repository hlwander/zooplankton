#script for calculating zooplankton net efficiency in Beaverdam Reservoir
#Created: 22July2020

#load in packages
pacman::p_load(dplyr,ggplot2)

#read in most recent zooplankton tow vs. schindler trap data
zoop<- read.csv("2014-2016_zooplankton.csv", header=TRUE)

#add in Date column
zoop$Date <- format(as.POSIXct(zoop$DateTime, format="%Y-%m-%d %H:%M:%S"), format="%Y-%m-%d")

#only select BVR data
zoop_methods_comparison<- zoop[zoop$Reservoir=="BVR",]

#subset data to days when schindler and tows were collected
zoop_methods_comparison <- zoop_methods_comparison[zoop_methods_comparison$Date=="2016-08-03" | zoop_methods_comparison$Date=="2016-08-04",]

#calculate volume and count
zoop_methods_comparison_calcs <- zoop_methods_comparison %>% group_by(DateTime,StartDepth_m,CollectionMethod)  %>%
  mutate(Volume_L=MeanWeight_ug/Biomass_ugL)  %>% mutate(Count_n=Density_IndPerL * Volume_L)

#drop chaoborus
zoop_methods_comparison_calcs <- zoop_methods_comparison_calcs[!zoop_methods_comparison_calcs$Taxon=="Chaoborus",]

#sum count across datetime, StartDepth, and collectionMethod
zoop_methods_comparison_calcs <- zoop_methods_comparison_calcs %>% group_by(DateTime,StartDepth_m,CollectionMethod)  %>%
  mutate(TotalCount= sum(Count_n, na.rm=TRUE))

#then average totalcount in new df used for calcs below --> averaging takes into account the NAs
zoop_totalcount <- zoop_methods_comparison_calcs %>% select(Reservoir,DateTime,StartDepth_m,CollectionMethod,TotalCount) %>%
  group_by(DateTime,StartDepth_m,CollectionMethod)  %>% summarise(TotalCount_n=mean(TotalCount))

#add column for number of pooled depths and unadjusted volume (pi*r^2*depth * 1000L / 1m^3)
netArea <- pi * (0.1524)^2 *1000 
zoop_totalcount <- zoop_totalcount %>% group_by(DateTime,CollectionMethod)  %>% 
  mutate(Num_pooled_Depths = length(StartDepth_m))  %>% mutate(Vol_unadj = netArea * StartDepth_m)

#replace vol_unadj column with NA if collectionMethod==Schindler
zoop_totalcount$Vol_unadj[zoop_totalcount$CollectionMethod=="Schindler"] <- NA

#separate tow data from schindler data
Schindler_totalCount <- zoop_totalcount[zoop_totalcount$CollectionMethod=="Schindler",]
Tow_totalCount_final <- zoop_totalcount[zoop_totalcount$CollectionMethod=="Tow",colnames(zoop_totalcount)!="Num_pooled_Depths" & colnames(zoop_totalcount)!="CollectionMethod"]

#sum total count for DateTime 
Schindler_totalCount_final <- Schindler_totalCount %>% select(DateTime,StartDepth_m,TotalCount_n) %>%
  group_by(DateTime)  %>% summarise(StartDepth_m=last(StartDepth_m),TotalCount_n=sum(TotalCount_n))
  
#get rid of 8,9,10 m to sum all 7 m data
Schindler_totalCount_7 <- Schindler_totalCount[Schindler_totalCount$StartDepth_m<=7,]

Schindler_totalCount_7 <- Schindler_totalCount_7 %>% select(DateTime,StartDepth_m,TotalCount_n) %>%
  group_by(DateTime)  %>% summarise(StartDepth_m=last(StartDepth_m),TotalCount_n=sum(TotalCount_n))

#add column for totalCount_7m
Schindler_totalCount_final$TotalCount_7m <- Schindler_totalCount_7$TotalCount_n

#initialize df
Density.neteff <- data.frame("DateTime"=unique(Tow_totalCount_final$DateTime))

#### Calculate APPARENT density from vertical tows ####
#count / volume_unadj 
for(i in 1:length(Density.neteff$DateTime)){
  Density.neteff$Apparent_dens_7m[i] <- Tow_totalCount_final$TotalCount_n[Tow_totalCount_final$StartDepth_m==7][i] / Tow_totalCount_final$Vol_unadj[Tow_totalCount_final$StartDepth_m==7][i]
  Density.neteff$Apparent_dens_10m[i] <-  Tow_totalCount_final$TotalCount_n[Tow_totalCount_final$StartDepth_m==10][i] / Tow_totalCount_final$Vol_unadj[Tow_totalCount_final$StartDepth_m==10][i]
}

#### Calculate ACTUAL density from schindler traps ####
#count / 30L (vol of schindler trap) * # of samples (pooled from all depths)
for(i in 1:length(Density.neteff$DateTime)){
  Density.neteff$Actual_dens_7m[i] <- Schindler_totalCount_final$TotalCount_7m[i] / (30*8)
  Density.neteff$Actual_dens_10m[i] <-  Schindler_totalCount_final$TotalCount_n[i] / (30*11)
}

#### Net Efficiency = APPARENT density / ACTUAL density ####
for(i in 1:length(Density.neteff$DateTime)){
  Density.neteff$NetEff_7m[i] <- Density.neteff$Apparent_dens_7m[i]/Density.neteff$Actual_dens_7m[i]
  Density.neteff$NetEff_10m[i] <- Density.neteff$Apparent_dens_10m[i]/Density.neteff$Actual_dens_10m[i]
}

#Going to take the average 10m net efficiency across noon samples only
NetEfficiency2016 <- mean(Density.neteff$NetEff_10m[Density.neteff$DateTime=="2016-08-03 12:00:00" | Density.neteff$DateTime=="2016-08-04 12:00:00"])

#---------------------------#
#visualize density by depth #
#---------------------------#
#order depth by decreasing number
Schindler_totalCount <- Schindler_totalCount[with(Schindler_totalCount,order(StartDepth_m)),]

#jpeg("Schindler_density_vs_depth.jpg", width = 6, height = 5, units = "in",res = 300)
ggplot(data=Schindler_totalCount,aes(x=TotalCount_n/30, y=StartDepth_m,color=DateTime)) + geom_point() +
  scale_y_reverse() + geom_path()
#dev.off()

##########################################
###########  CCR as a check  #############
##########################################

#only select CCR data
zoop_methods_comparison_CCR<- zoop[zoop$Reservoir=="CCR",]

#subset data to days when schindler and tows were collected
zoop_methods_comparison_CCR <- zoop_methods_comparison_CCR[zoop_methods_comparison_CCR$Date=="2016-05-27" | zoop_methods_comparison_CCR$Date=="2016-05-28" | 
                                                         zoop_methods_comparison_CCR$Date=="2016-08-12" | zoop_methods_comparison_CCR$Date=="2016-08-13",]

#if collection method = tow and end depth = 7 then ignore
zoop_methods_comparison_CCR <- zoop_methods_comparison_CCR[!(zoop_methods_comparison_CCR$CollectionMethod=="Tow" & zoop_methods_comparison_CCR$EndDepth_m==7),]

#get rid of 8m schindlers and 21m tows 
zoop_methods_comparison_CCR <- zoop_methods_comparison_CCR[!(zoop_methods_comparison_CCR$CollectionMethod=="Schindler" & zoop_methods_comparison_CCR$StartDepth_m==8),]
zoop_methods_comparison_CCR <- zoop_methods_comparison_CCR[!(zoop_methods_comparison_CCR$CollectionMethod=="Tow" & zoop_methods_comparison_CCR$StartDepth_m==21),]

#calculate volume and count
zoop_methods_comparison_calcs_CCR <- zoop_methods_comparison_CCR %>% group_by(DateTime,StartDepth_m,CollectionMethod)  %>%
  mutate(Volume_L=MeanWeight_ug/Biomass_ugL)  %>% mutate(Count_n=Density_IndPerL * Volume_L)

#drop chaoborus
zoop_methods_comparison_calcs_CCR <- zoop_methods_comparison_calcs_CCR[!zoop_methods_comparison_calcs_CCR$Taxon=="Chaoborus",]

#sum count across datetime, StartDepth, and collectionMethod
zoop_methods_comparison_calcs_CCR <- zoop_methods_comparison_calcs_CCR %>% group_by(DateTime,StartDepth_m,CollectionMethod)  %>%
  mutate(TotalCount= sum(Count_n, na.rm=TRUE))

#then average totalcount in new df used for calcs below
zoop_totalcount_CCR <- zoop_methods_comparison_calcs_CCR %>% select(Reservoir,DateTime,StartDepth_m,CollectionMethod,TotalCount) %>%
  group_by(DateTime,StartDepth_m,CollectionMethod)  %>% summarise(TotalCount_n=mean(TotalCount))

#add column for number of pooled depths and unadjusted volume (pi*r^2*depth * 1000L / 1m^3)
netArea <- pi * (0.1524)^2 *1000 
zoop_totalcount_CCR <- zoop_totalcount_CCR %>% group_by(DateTime,CollectionMethod)  %>% 
  mutate(Num_pooled_Depths = length(StartDepth_m))  %>% mutate(Vol_unadj = netArea * StartDepth_m)

#replace vol_unadj column with NA if collectionMethod==Schindler
zoop_totalcount_CCR$Vol_unadj[zoop_totalcount_CCR$CollectionMethod=="Schindler"] <- NA

#separate tow data from schindler data
Schindler_totalCount_CCR <- zoop_totalcount_CCR[zoop_totalcount_CCR$CollectionMethod=="Schindler",]
Tow_totalCount_final_CCR <- zoop_totalcount_CCR[zoop_totalcount_CCR$CollectionMethod=="Tow",colnames(zoop_totalcount_CCR)!="Num_pooled_Depths" & colnames(zoop_totalcount_CCR)!="CollectionMethod"]

#sum total count for DateTime 
Schindler_totalCount_final_CCR <- Schindler_totalCount_CCR %>% select(DateTime,StartDepth_m,TotalCount_n) %>%
  group_by(DateTime)  %>% summarise(StartDepth_m=last(StartDepth_m),TotalCount_n=sum(TotalCount_n))

#initialize df
Density.neteff.CCR <- data.frame("DateTime"=unique(Tow_totalCount_final_CCR$DateTime))

#### Calculate APPARENT density from vertical tows ####
#count / volume_unadj 
for(i in 1:length(Density.neteff.CCR$DateTime)){
  Density.neteff.CCR$Apparent_dens_7m[i] <- Tow_totalCount_final_CCR$TotalCount_n[i] / Tow_totalCount_final_CCR$Vol_unadj[i]
}

#### Calculate ACTUAL density from schindler traps ####
#count / 30L (vol of schindler trap) * # of samples (pooled from all depths)
for(i in 1:length(Density.neteff.CCR$DateTime)){
  Density.neteff.CCR$Actual_dens_7m[i] <- Schindler_totalCount_final_CCR$TotalCount_n[i] / (30*8)
}

#### Net Efficiency = APPARENT density / ACTUAL density ####
for(i in 1:length(Density.neteff.CCR$DateTime)){
  Density.neteff.CCR$NetEff_7m[i] <- Density.neteff.CCR$Apparent_dens_7m[i]/Density.neteff.CCR$Actual_dens_7m[i]
}

mean(Density.neteff.CCR$NetEff_7m)
sd(Density.neteff.CCR$NetEff_7m)

