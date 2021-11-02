#script for calculating zooplankton net efficiency in Beaverdam Reservoir
#Created: 20 November 2020
#for some reason, did not collect full water column tow at the same time as schindlers (??!!)
#options are to either average noon and midnight full densities or just use one or the other

#load in packages
pacman::p_load(dplyr,ggplot2)

#function to count characters starting at the end of the string
substrEnd <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#read in zoop data 
zoop<- read.csv('SummaryStats/FCR_ZooplanktonSummary2020.csv',header = TRUE)

#select cols to keep
keep <- c("sample_ID","site_no","collect_date","Hour","DepthOfTow_m","Zooplankton_No.","Volume_L","Volume_unadj","mesh_size_μm")

zoop_totalcount <- zoop[(!is.na(zoop$mesh_size_μm) & zoop$mesh_size_μm==61 | zoop$sample_ID=="B_pel_12Aug20_noon_rep1" | zoop$sample_ID=="B_pel_12Aug20_noon_rep2" |
                           zoop$sample_ID=="B_pel_13Aug20_midnight_rep1" | zoop$sample_ID=="B_pel_13Aug20_midnight_rep2" | zoop$sample_ID=="F11Sep20_9.0m" |
                           zoop$sample_ID=="F10Sep20_9.0m" | zoop$sample_ID=="F10Sep20_2.5m" | zoop$sample_ID=="F11Sep20_2.2m"),keep]

#round to nearest hour
zoop_totalcount$Hour <- format(strptime(zoop_totalcount$Hour,format="%H"),"%H")
zoop_totalcount$Hour <- ifelse(zoop_totalcount$Hour==20, 19, ifelse(zoop_totalcount$Hour==11, 12, zoop_totalcount$Hour))

#replace vol_unadj column with NA for schindler traps
zoop_totalcount$Volume_unadj[zoop_totalcount$mesh_size_μm ==61] <- NA

#------------------------------------------------------------------------------#
#### BVR and FCR 2020 net efficiency ####

#separate tow data from schindler data
Schindler_totalCount <- zoop_totalcount[zoop_totalcount$mesh_size_μm==61,]
Tow_totalCount_final <- zoop_totalcount[zoop_totalcount$mesh_size_μm==80,]

#add a column for rep
Schindler_totalCount$Rep <- substrEnd(Schindler_totalCount$sample_ID,1)

#sum total count for rep1 vs rep2 for each reservoir
Schindler_totalCount_final <- Schindler_totalCount %>% select(collect_date, site_no,Hour,DepthOfTow_m,Zooplankton_No.,Rep) %>%
  group_by(site_no,Rep)  %>% summarise(DepthOfTow_m=max(DepthOfTow_m),TotalCount_n=sum(Zooplankton_No.))
  
#initialize df
Density.neteff <- data.frame("sample_ID"=unique(Tow_totalCount_final$sample_ID))

#### Calculate APPARENT density from vertical tows ####
#count / volume_unadj (because we want to whole total volume not volume counted)
for(i in 1:length(Density.neteff$sample_ID)){
  Density.neteff$Apparent_dens[i] <-  Tow_totalCount_final$Zooplankton_No.[i] / Tow_totalCount_final$Volume_unadj[i]
}

#### Calculate ACTUAL density from schindler traps ####
#count / 30L (vol of schindler trap) * # of samples (pooled from all depths)
#repeating these for noon and midnight samples to see if there is a difference
for(i in 1:length(Density.neteff$sample_ID)){
  Density.neteff$Actual_dens[i] <-  Schindler_totalCount_final$TotalCount_n[i] / (30*11)
  Density.neteff$Actual_dens[3] <- Density.neteff$Actual_dens[1]
  Density.neteff$Actual_dens[4] <- Density.neteff$Actual_dens[2]
  Density.neteff$Actual_dens[5:8] <- Schindler_totalCount_final$TotalCount_n[3] / (30*11)
  
}

#### Net Efficiency = APPARENT density / ACTUAL density ####
for(i in 1:length(Density.neteff$sample_ID)){
  Density.neteff$NetEff[i] <- Density.neteff$Apparent_dens[i]/Density.neteff$Actual_dens[i]
}

#Going to take the average net efficiency across both reps because values are super close to each other
NetEfficiency2020 <- c(mean(Density.neteff$NetEff[1:4]),mean(Density.neteff$NetEff[6],Density.neteff$NetEff[8]), mean(Density.neteff$NetEff[5],Density.neteff$NetEff[7]))
  
#---------------------------#
#visualize density by depth #
#---------------------------#
#order depth by decreasing number
Schindler_totalCount <- Schindler_totalCount[with(Schindler_totalCount,order(DepthOfTow_m)),]

#jpeg("Figures/Schindler_density_vs_depth.jpg", width = 6, height = 5, units = "in",res = 300)
ggplot(data=Schindler_totalCount,aes(x=Zooplankton_No./30, y=DepthOfTow_m,color=Hour)) + geom_point() +
  scale_y_reverse() + geom_path() + facet_grid(~site_no)
#dev.off()

#summarize schindler_totalCount so one # per depth
Schindler_avgCount <- Schindler_totalCount %>% group_by(site_no, DepthOfTow_m) %>% summarise(mean_num = mean(Zooplankton_No.))

#jpeg("Figures/AvgSchindler_density_vs_depth.jpg", width = 6, height = 5, units = "in",res = 300)
ggplot(data=Schindler_avgCount,aes(x=mean_num, y=DepthOfTow_m)) + geom_point() +
  scale_y_reverse() + geom_path() + xlab("Zooplankton (#)") + ylab("Depth (m)")  + facet_grid(~site_no)
#dev.off()

