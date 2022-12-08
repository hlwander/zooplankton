#FCR MOM figs - noon and midnight tows 

#load libraries
pacman::p_load(ggplot2, tidyverse,viridis)

#----------------------------------------------------------------------#
#functions

#count characters starting at the end of the string
substrEnd <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#calculate the standard error
stderr <- function(x) {
  sd(x,na.rm=TRUE)/sqrt(length(na.omit(x)))
}

#not sure if this is right, but calculating SE of difference between epi mean and hypo mean
SE.diffMean<- function(x,y){
  sqrt((sd(x,na.rm=TRUE)^2/length(na.omit(x))) + 
         (sd(y,na.rm=TRUE)^2/length(na.omit(y))))
}

#----------------------------------------------------------------------#
#read in zoop DVM stats
fcr_zoops_2020<- read.csv('Summer2020-DataAnalysis/SummaryStats/FCR_ZooplanktonSummary2020.csv',header = TRUE) %>% filter(site_no %in% c("FCR_50", "FCR_schind"))
fcr_zoops_2021<- read.csv('Summer2021-DataAnalysis/SummaryStats/FCR_ZooplanktonSummary2021.csv',header = TRUE) %>% filter(site_no %in% c("FCR_50", "FCR_schind"))

#combine 2020 and 2021
fcr_zoops_mom <- rbind(fcr_zoops_2020, fcr_zoops_2021)

#drop 20 um samples
fcr_zoops_mom <- fcr_zoops_mom[!c(fcr_zoops_mom$mesh_size_μm=="20"),]

#manually set hour to midnight or noon
fcr_zoops_mom$Hour <- ifelse(grepl("midnight",fcr_zoops_mom$sample_ID),00, 12)

#pull rep # off as new column
fcr_zoops_mom$rep <- ifelse(substrEnd(fcr_zoops_mom$sample_ID,4)=="rep1" |substrEnd(fcr_zoops_mom$sample_ID,4)=="rep2",
                            substrEnd(fcr_zoops_mom$sample_ID,1),NA)

#drop reps from sampleIDs
fcr_zoops_mom$sample_ID <- ifelse(substrEnd(fcr_zoops_mom$sample_ID,4)=="rep1" |substrEnd(fcr_zoops_mom$sample_ID,4)=="rep2",
                                  substr(fcr_zoops_mom$sample_ID,1,nchar(fcr_zoops_mom$sample_ID)-5),fcr_zoops_mom$sample_ID)


#remove sample where sizeid was not digitized (n=1, F10Jun21_midnight_epi_rep1 sizeid lost in digitization process)
fcr_zoops_mom <- fcr_zoops_mom[!is.na(fcr_zoops_mom$OverallMeanSize_mm),]

#Replace NAs with 0
fcr_zoops_mom[20:204][is.na(fcr_zoops_mom[20:204])] <- 0

#drop oxy sample
fcr_zoops_mom <- fcr_zoops_mom[!c(substrEnd(fcr_zoops_mom$sample_ID,3)=="oxy"),]

##### Create new df to combine reps over 24 hours
zoop.repmeans <- fcr_zoops_mom %>% select(sample_ID,site_no,collect_date,Hour, Volume_L, Volume_unadj, proportional_vol, ZoopDensity_No.pL, OverallCount_n, TotalBiomass_ug,
                                          BiomassConcentration_ugpL,Cladocera_density_NopL, Cladocera_BiomassConcentration_ugpL, CladoceraCount_n, Cladocera_totalbiomass_ug, Cladocera_PercentOfTotal,
                                          Cyclopoida_density_NopL, Cyclopoida_BiomassConcentration_ugpL, CyclopoidaCount_n, Cyclopoida_totalbiomass_ug, Cyclopoida_PercentOfTotal,
                                          Rotifera_density_NopL,Rotifera_BiomassConcentration_ugpL, RotiferaCount_n, Rotifera_totalbiomass_ug, Rotifera_PercentOfTotal, 
                                          Calanoida_density_NopL, Calanoida_BiomassConcentration_ugpL, CalanoidaCount_n, Calanoida_PercentOfTotal, Calanoida_totalbiomass_ug,
                                          Copepoda_density_NopL, Copepoda_BiomassConcentration_ugpL, CopepodaCount_n, Copepoda_PercentOfTotal, Copepoda_totalbiomass_ug) %>%
  group_by(sample_ID, site_no, collect_date ,Hour) %>%
  summarise_at(vars(Volume_L:Copepoda_totalbiomass_ug,), funs(rep.mean=mean, rep.SE=stderr))

#separate tow vs schindler data
zoop_repmeans_tows <- zoop.repmeans[zoop.repmeans$site_no=="FCR_50",]
zoop_repmeans_schind <- zoop.repmeans[zoop.repmeans$site_no=="FCR_schind",]


#two dfs for raw (no net efficiency calcs needed) and volume calculated values
FCR_taxa_DVM_raw<- zoop_repmeans_tows[,c(1:4,6,7,which(substrEnd(colnames(zoop_repmeans_tows),10)=="n_rep.mean"),
                                         which(substrEnd(colnames(zoop_repmeans_tows),11)=="ug_rep.mean"),
                                         which(substrEnd(colnames(zoop_repmeans_tows),8)=="n_rep.SE"),
                                         which(substrEnd(colnames(zoop_repmeans_tows),9)=="ug_rep.SE"))]


FCR_DVM_vol_calculated <- zoop_repmeans_tows[,c(1:4,which(substrEnd(colnames(zoop_repmeans_tows),16)=="y_No.pL_rep.mean"),
                                                        which(substrEnd(colnames(zoop_repmeans_tows),15)=="y_NopL_rep.mean"),
                                                        which(substrEnd(colnames(zoop_repmeans_tows),15)=="n_ugpL_rep.mean"),
                                                        which(substrEnd(colnames(zoop_repmeans_tows),13)=="n_ugpL_rep.SE"),
                                                        which(substrEnd(colnames(zoop_repmeans_tows),14)=="y_No.pL_rep.SE"),
                                                        which(substrEnd(colnames(zoop_repmeans_tows),13)=="y_NopL_rep.SE"))]


#another one for percent calcs
FCR_DVM_percent<- zoop_repmeans_tows[,c(1:4,6,7,which(substrEnd(colnames(zoop_repmeans_tows),14)=="Total_rep.mean"),
                                        which(substrEnd(colnames(zoop_repmeans_tows),12)=="Total_rep.SE"))]

#----------------------------------------------------------------------#
#        DFs to calculate hypo and epi density and biomass             #
#----------------------------------------------------------------------#
fcr_zoops_mom_tows <- fcr_zoops_mom[fcr_zoops_mom$site_no=="FCR_50",]
fcr_zoops_mom_schind <- fcr_zoops_mom[fcr_zoops_mom$site_no=="FCR_schind",]

#df for calculating epi vs hypo density/biomass and then reorder
FCR_DVM_calcs_raw<- data.frame("SampleID"=unique(substr(fcr_zoops_mom_tows$sample_ID,1,13))) %>% arrange(SampleID,unique(substr(FCR_taxa_DVM_raw$sample_ID,1,13)))
FCR_DVM_calcs_vol_calc<- data.frame("SampleID"=unique(substr(fcr_zoops_mom_tows$sample_ID,1,13))) %>% arrange(SampleID,unique(substr(FCR_taxa_DVM_raw$sample_ID,1,13)))
FCR_DVM_calcs_percent <- data.frame("SampleID"=unique(substr(fcr_zoops_mom_tows$sample_ID,1,13))) %>% arrange(SampleID,unique(substr(FCR_taxa_DVM_raw$sample_ID,1,13)))

#for loop to fill out epi vs hypo calcs RAW counts and ug
variables<- colnames(FCR_taxa_DVM_raw)[c(7:18)]
for(i in 1:length(variables)){
  FCR_DVM_calcs_raw[,paste0(variables,"_epi")[i]]<- FCR_taxa_DVM_raw[substrEnd(FCR_taxa_DVM_raw$sample_ID,3)=="epi",paste0(variables)[i]]
  FCR_DVM_calcs_raw[,paste0(variables,"_hypo")[i]] <- FCR_taxa_DVM_raw[substrEnd(FCR_taxa_DVM_raw$sample_ID,3)!="epi",paste0(variables)[i]] - FCR_DVM_calcs_raw[,paste0(variables,"_epi")[i]]
}

#Net efficiencies for fcr 2020 and 2021; note that I'm assuming that neteff is 100% for epi tows (>95% in 2021 and >100% in 2020)
netefficiency <- c(0.05279169, 0.06146684)

temp <- FCR_DVM_calcs_vol_calc

#hypo density and biomass calculated by subtracting epi raw zoop # from full zoop # and then dividing by the (full volume - epi volume) 
#NOTE: using epi density/L and biomass/L but calculating hypo using raw # and ug values. AND including the net efficiency to correct for the tows when using raw counts
column.names<- colnames(FCR_DVM_vol_calculated[,c(5:16)])
variables<- colnames(FCR_taxa_DVM_raw)[c(7:18)]
percent<- colnames(FCR_DVM_percent[,c(7:11)])
for(i in 1:length(variables)){
  FCR_DVM_calcs_vol_calc[,paste0(column.names,"_epi")[i]]<- FCR_DVM_vol_calculated[substrEnd(FCR_DVM_vol_calculated$sample_ID,3)=="epi",paste0(column.names)[i]]
  FCR_DVM_calcs_vol_calc[,paste0(column.names,"_hypo")[i]] <- (((FCR_taxa_DVM_raw[substrEnd(FCR_taxa_DVM_raw$sample_ID,3)!="epi" ,"proportional_vol_rep.mean"]) * FCR_taxa_DVM_raw[substrEnd(FCR_taxa_DVM_raw$sample_ID,3)!="epi" ,paste0(variables)[i]]) * (1/netefficiency[1])) - 
                                                               ((FCR_taxa_DVM_raw[substrEnd(FCR_taxa_DVM_raw$sample_ID,3)=="epi", "proportional_vol_rep.mean"]) *  FCR_taxa_DVM_raw[substrEnd(FCR_taxa_DVM_raw$sample_ID,3)=="epi",paste0(variables)[i]])/ 
                                                               (FCR_taxa_DVM_raw[substrEnd(FCR_taxa_DVM_raw$sample_ID,3)!="epi" ,"Volume_unadj_rep.mean"] - FCR_taxa_DVM_raw[substrEnd(FCR_taxa_DVM_raw$sample_ID,3)=="epi", "Volume_unadj_rep.mean"])                                                                     
  temp[,paste0(column.names,"_hypo")[i]] <- (((FCR_taxa_DVM_raw[substrEnd(FCR_taxa_DVM_raw$sample_ID,3)!="epi" ,"proportional_vol_rep.mean"]) * FCR_taxa_DVM_raw[substrEnd(FCR_taxa_DVM_raw$sample_ID,3)!="epi" ,paste0(variables)[i]]) * (1/netefficiency[2])) - 
                                                               ((FCR_taxa_DVM_raw[substrEnd(FCR_taxa_DVM_raw$sample_ID,3)=="epi", "proportional_vol_rep.mean"]) *  FCR_taxa_DVM_raw[substrEnd(FCR_taxa_DVM_raw$sample_ID,3)=="epi",paste0(variables)[i]])/ 
                                                               (FCR_taxa_DVM_raw[substrEnd(FCR_taxa_DVM_raw$sample_ID,3)!="epi" ,"Volume_unadj_rep.mean"] - FCR_taxa_DVM_raw[substrEnd(FCR_taxa_DVM_raw$sample_ID,3)=="epi", "Volume_unadj_rep.mean"])                                                                                  
} 
#now manually adjusting rows 2+3 because these are the 21 samples and have a different net efficiency
FCR_DVM_calcs_vol_calc[2,paste0(column.names,"_hypo")] <- temp[2,c(2:13)]
FCR_DVM_calcs_vol_calc[3,paste0(column.names,"_hypo")] <- temp[3,c(2:13)]

#percent density
density.percent<- colnames(FCR_DVM_vol_calculated[,c(6:10)]) 
for(i in 1:length(density.percent)){
  for(j in 1:length(unique(FCR_DVM_calcs_vol_calc$SampleID))){
    FCR_DVM_calcs_vol_calc[j,paste0(density.percent,"_epi_percent_density")[i]]<- (FCR_DVM_calcs_vol_calc[j,paste0(density.percent,"_epi")][i]/ sum(FCR_DVM_calcs_vol_calc[j,paste0(density.percent,"_epi")[i]],FCR_DVM_calcs_vol_calc[j,paste0(density.percent,"_hypo")[i]])) *100
    FCR_DVM_calcs_vol_calc[j,paste0(density.percent,"_hypo_percent_density")[i]]<- (FCR_DVM_calcs_vol_calc[j,paste0(density.percent,"_hypo")][i]/ sum(FCR_DVM_calcs_vol_calc[j,paste0(density.percent,"_epi")][i],FCR_DVM_calcs_vol_calc[j,paste0(density.percent,"_hypo")][i])) * 100
    
  }       
}

#################################
#           SE dfs              #
#################################

#need reps from zoop df to calculate SE for raw and vol_calculated
matchingcols <- match(substr(colnames(FCR_taxa_DVM_raw[,c(1:4,8:18)]),1,14),substr(colnames(fcr_zoops_mom_tows),1,14))
DVM_samples_raw<- fcr_zoops_mom_tows[,unique(matchingcols)] 

matchingcols <- match(substr(colnames(FCR_DVM_vol_calculated[,c(1:16)]),1,14),substr(colnames(fcr_zoops_mom_tows),1,14))
DVM_samples_dens<- fcr_zoops_mom_tows[,unique(matchingcols)] 

#separate full vs epi samples
FullSamples <- data.frame(unique(DVM_samples_raw$sample_ID[substrEnd(DVM_samples_raw$sample_ID,3)!="epi"])) %>% arrange(FullSamples,unique(substr(FCR_taxa_DVM_raw$sample_ID,1,13)))
EpiSamples<- data.frame(unique(DVM_samples_raw$sample_ID[substrEnd(DVM_samples_raw$sample_ID,3)=="epi"])) %>% arrange(EpiSamples,unique(substr(FCR_taxa_DVM_raw$sample_ID,1,13)))
#check bvr dvm calcs because I think the order pof dfs might be slightly off

#initialize df
FCR_DVM_calcs_SE<- data.frame("SampleID"=unique(substr(fcr_zoops_mom_tows$sample_ID,1,13))) %>% arrange(SampleID,unique(substr(FCR_taxa_DVM_raw$sample_ID,1,13)))

#calculate hypo SE 
SEonly<- colnames(FCR_taxa_DVM_raw)[7:18]
SEonly <- substr(SEonly,1,nchar(SEonly)-9)

Percentdens <- colnames(FCR_DVM_vol_calculated)[6:10]
Percentdens <- substr(Percentdens,1,nchar(Percentdens)-9)

for(i in 1:length(SEonly)){ #hypo only works if there is a noon and epi rep (which there is not...)
  FCR_DVM_calcs_SE[,paste0(column.names,"_epi_SE")[i]] <- FCR_DVM_vol_calculated[substrEnd(FCR_DVM_vol_calculated$sample_ID,3)=="epi",substrEnd(colnames(FCR_DVM_vol_calculated),2)=="SE"][i]
  FCR_DVM_calcs_SE[1,paste0(column.names,"_hypo_SE")[i]]<- SE.diffMean(DVM_samples_raw[DVM_samples_raw$sample_ID==FullSamples[[1]][1],paste0(SEonly)[i]],
                                                                       DVM_samples_raw[DVM_samples_raw$sample_ID==EpiSamples[[1]][1],paste0(SEonly)[i]])
  FCR_DVM_calcs_SE[2,paste0(column.names,"_hypo_SE")[i]]<- SE.diffMean(DVM_samples_raw[DVM_samples_raw$sample_ID==FullSamples[[1]][2],paste0(SEonly)[i]],
                                                                       DVM_samples_raw[DVM_samples_raw$sample_ID==EpiSamples[[1]][2],paste0(SEonly)[i]])
  FCR_DVM_calcs_SE[3,paste0(column.names,"_hypo_SE")[i]]<- SE.diffMean(DVM_samples_raw[DVM_samples_raw$sample_ID==FullSamples[[1]][3],paste0(SEonly)[i]],
                                                                       DVM_samples_raw[DVM_samples_raw$sample_ID==EpiSamples[[1]][3],paste0(SEonly)[i]])
  FCR_DVM_calcs_SE[4,paste0(column.names,"_hypo_SE")[i]]<- SE.diffMean(DVM_samples_raw[DVM_samples_raw$sample_ID==FullSamples[[1]][4],paste0(SEonly)[i]],
                                                                       DVM_samples_raw[DVM_samples_raw$sample_ID==EpiSamples[[1]][4],paste0(SEonly)[i]])
  FCR_DVM_calcs_SE[5,paste0(column.names,"_hypo_SE")[i]]<- SE.diffMean(DVM_samples_raw[DVM_samples_raw$sample_ID==FullSamples[[1]][5],paste0(SEonly)[i]],
                                                                       DVM_samples_raw[DVM_samples_raw$sample_ID==EpiSamples[[1]][5],paste0(SEonly)[i]])
  FCR_DVM_calcs_SE[6,paste0(column.names,"_hypo_SE")[i]]<- SE.diffMean(DVM_samples_raw[DVM_samples_raw$sample_ID==FullSamples[[1]][6],paste0(SEonly)[i]],
                                                                       DVM_samples_raw[DVM_samples_raw$sample_ID==EpiSamples[[1]][6],paste0(SEonly)[i]])
  FCR_DVM_calcs_SE[7,paste0(column.names,"_hypo_SE")[i]]<- SE.diffMean(DVM_samples_raw[DVM_samples_raw$sample_ID==FullSamples[[1]][7],paste0(SEonly)[i]],
                                                                       DVM_samples_raw[DVM_samples_raw$sample_ID==EpiSamples[[1]][7],paste0(SEonly)[i]])
  FCR_DVM_calcs_SE[8,paste0(column.names,"_hypo_SE")[i]]<- SE.diffMean(DVM_samples_raw[DVM_samples_raw$sample_ID==FullSamples[[1]][8],paste0(SEonly)[i]],
                                                                       DVM_samples_raw[DVM_samples_raw$sample_ID==EpiSamples[[1]][8],paste0(SEonly)[i]])
  FCR_DVM_calcs_SE[9,paste0(column.names,"_hypo_SE")[i]]<- SE.diffMean(DVM_samples_raw[DVM_samples_raw$sample_ID==FullSamples[[1]][9],paste0(SEonly)[i]],
                                                                       DVM_samples_raw[DVM_samples_raw$sample_ID==EpiSamples[[1]][9],paste0(SEonly)[i]])
  FCR_DVM_calcs_SE[10,paste0(column.names,"_hypo_SE")[i]]<- SE.diffMean(DVM_samples_raw[DVM_samples_raw$sample_ID==FullSamples[[1]][10],paste0(SEonly)[i]],
                                                                       DVM_samples_raw[DVM_samples_raw$sample_ID==EpiSamples[[1]][10],paste0(SEonly)[i]])
  FCR_DVM_calcs_SE[11,paste0(column.names,"_hypo_SE")[i]]<- SE.diffMean(DVM_samples_raw[DVM_samples_raw$sample_ID==FullSamples[[1]][11],paste0(SEonly)[i]],
                                                                       DVM_samples_raw[DVM_samples_raw$sample_ID==EpiSamples[[1]][11],paste0(SEonly)[i]])
  FCR_DVM_calcs_SE[12,paste0(column.names,"_hypo_SE")[i]]<- SE.diffMean(DVM_samples_raw[DVM_samples_raw$sample_ID==FullSamples[[1]][12],paste0(SEonly)[i]],
                                                                       DVM_samples_raw[DVM_samples_raw$sample_ID==EpiSamples[[1]][12],paste0(SEonly)[i]])
  
  
} #Note: just repeating hypo SE calc for each date because I can't figure out the for loop across rows and columns...


for (j in 1:length(Percentdens)){
  FCR_DVM_calcs_SE[,paste0(Percentdens,"_epi_percent_density_SE")[j]]<- FCR_DVM_vol_calculated[substrEnd(FCR_DVM_vol_calculated$sample_ID,3)=="epi",substrEnd(colnames(FCR_DVM_vol_calculated),11)=="NopL_rep.SE"][j]
  FCR_DVM_calcs_SE[1,paste0(Percentdens,"_hypo_percent_density_SE")[j]] <- SE.diffMean(DVM_samples_dens[DVM_samples_dens$sample_ID==FullSamples[[1]][1],paste0(Percentdens)[j]],
                                                                                       DVM_samples_dens[DVM_samples_dens$sample_ID==EpiSamples[[1]][1],paste0(Percentdens)[j]])
  FCR_DVM_calcs_SE[2,paste0(Percentdens,"_hypo_percent_density_SE")[j]] <- SE.diffMean(DVM_samples_dens[DVM_samples_dens$sample_ID==FullSamples[[1]][2],paste0(Percentdens)[j]],
                                                                                       DVM_samples_dens[DVM_samples_dens$sample_ID==EpiSamples[[1]][2],paste0(Percentdens)[j]])
  FCR_DVM_calcs_SE[3,paste0(Percentdens,"_hypo_percent_density_SE")[j]] <- SE.diffMean(DVM_samples_dens[DVM_samples_dens$sample_ID==FullSamples[[1]][3],paste0(Percentdens)[j]],
                                                                                       DVM_samples_dens[DVM_samples_dens$sample_ID==EpiSamples[[1]][3],paste0(Percentdens)[j]])
  FCR_DVM_calcs_SE[4,paste0(Percentdens,"_hypo_percent_density_SE")[j]] <- SE.diffMean(DVM_samples_dens[DVM_samples_dens$sample_ID==FullSamples[[1]][4],paste0(Percentdens)[j]],
                                                                                       DVM_samples_dens[DVM_samples_dens$sample_ID==EpiSamples[[1]][4],paste0(Percentdens)[j]])
  FCR_DVM_calcs_SE[5,paste0(Percentdens,"_hypo_percent_density_SE")[j]] <- SE.diffMean(DVM_samples_dens[DVM_samples_dens$sample_ID==FullSamples[[1]][5],paste0(Percentdens)[j]],
                                                                                       DVM_samples_dens[DVM_samples_dens$sample_ID==EpiSamples[[1]][5],paste0(Percentdens)[j]])
  FCR_DVM_calcs_SE[6,paste0(Percentdens,"_hypo_percent_density_SE")[j]] <- SE.diffMean(DVM_samples_dens[DVM_samples_dens$sample_ID==FullSamples[[1]][6],paste0(Percentdens)[j]],
                                                                                       DVM_samples_dens[DVM_samples_dens$sample_ID==EpiSamples[[1]][6],paste0(Percentdens)[j]])
  FCR_DVM_calcs_SE[7,paste0(Percentdens,"_hypo_percent_density_SE")[j]] <- SE.diffMean(DVM_samples_dens[DVM_samples_dens$sample_ID==FullSamples[[1]][7],paste0(Percentdens)[j]],
                                                                                       DVM_samples_dens[DVM_samples_dens$sample_ID==EpiSamples[[1]][7],paste0(Percentdens)[j]])
  FCR_DVM_calcs_SE[8,paste0(Percentdens,"_hypo_percent_density_SE")[j]] <- SE.diffMean(DVM_samples_dens[DVM_samples_dens$sample_ID==FullSamples[[1]][8],paste0(Percentdens)[j]],
                                                                                       DVM_samples_dens[DVM_samples_dens$sample_ID==EpiSamples[[1]][8],paste0(Percentdens)[j]])
  FCR_DVM_calcs_SE[9,paste0(Percentdens,"_hypo_percent_density_SE")[j]] <- SE.diffMean(DVM_samples_dens[DVM_samples_dens$sample_ID==FullSamples[[1]][9],paste0(Percentdens)[j]],
                                                                                       DVM_samples_dens[DVM_samples_dens$sample_ID==EpiSamples[[1]][9],paste0(Percentdens)[j]])
  FCR_DVM_calcs_SE[10,paste0(Percentdens,"_hypo_percent_density_SE")[j]] <- SE.diffMean(DVM_samples_dens[DVM_samples_dens$sample_ID==FullSamples[[1]][10],paste0(Percentdens)[j]],
                                                                                       DVM_samples_dens[DVM_samples_dens$sample_ID==EpiSamples[[1]][10],paste0(Percentdens)[j]])
  FCR_DVM_calcs_SE[11,paste0(Percentdens,"_hypo_percent_density_SE")[j]] <- SE.diffMean(DVM_samples_dens[DVM_samples_dens$sample_ID==FullSamples[[1]][11],paste0(Percentdens)[j]],
                                                                                       DVM_samples_dens[DVM_samples_dens$sample_ID==EpiSamples[[1]][11],paste0(Percentdens)[j]])
  FCR_DVM_calcs_SE[12,paste0(Percentdens,"_hypo_percent_density_SE")[j]] <- SE.diffMean(DVM_samples_dens[DVM_samples_dens$sample_ID==FullSamples[[1]][12],paste0(Percentdens)[j]],
                                                                                       DVM_samples_dens[DVM_samples_dens$sample_ID==EpiSamples[[1]][12],paste0(Percentdens)[j]])
}

#------------------------------------------------------------------------------#
#wide to long df
FCR_DVM_calcs_long <-   FCR_DVM_calcs_vol_calc %>%
  gather(metric,value, ZoopDensity_No.pL_rep.mean_epi:Copepoda_density_NopL_rep.mean_hypo_percent_density) 

FCR_DVM_calcs_SE_long <-  FCR_DVM_calcs_SE %>%
  gather(metric,SE,ZoopDensity_No.pL_rep.mean_epi_SE:Copepoda_density_NopL_hypo_percent_density_SE) 

#add the SE column from df2 to df1 for combined df
FCR_DVM_calcs_long$SE <- FCR_DVM_calcs_SE_long$SE

#add watercolumn and hour columns
FCR_DVM_calcs_long$WaterColumn <- ifelse(substrEnd(FCR_DVM_calcs_long$metric,3)=="epi" |substrEnd(FCR_DVM_calcs_long$metric,19)=="epi_percent_density" ,"epilimnion","hypolimnion")
FCR_DVM_calcs_long$Hour <- ifelse(substr(FCR_DVM_calcs_long$SampleID,10,13)=="noon","noon","midnight")
FCR_DVM_calcs_long$Taxa <- substr(FCR_DVM_calcs_long$metric,1,9)

#Export FCR MOM DVM stats
write.csv(FCR_DVM_calcs_long, "Summer2021-DataAnalysis/SummaryStats/FCR_DVM_2020-2021_zoops.csv")

#----------------------------------------------------------------------#
#                           FCR figures                                #
#simple boxplot comparing epi vs hypo in the presence of MOM vs anoxia #
#----------------------------------------------------------------------#
#ridiculous way to work with these stupid facet labels
facet_labeller_bot <- function(variable, value) {
  rep("",8)
}

facet_labeller_top <- function(variable, value) {
  c("","","","Midnight","","","","Noon")
}

#drop total zoop density/biomass
FCR_DVM_calcs_long <- FCR_DVM_calcs_long[substr(FCR_DVM_calcs_long$metric,1,17)!="ZoopDensity_No.pL" & substr(FCR_DVM_calcs_long$metric,1,17)!="BiomassConcentrat",]

#reorder taxa
FCR_DVM_calcs_long$Taxa <- factor(FCR_DVM_calcs_long$Taxa, levels = c("Cladocera","Rotifera_","Cyclopoid","Calanoida")) #"ZoopDensi", "BiomassCo"))

#taxa list
taxa <- c("Cladocera", "Rotifera_", "Cyclopoid", "Calanoida")

#jpeg("Figures/FCR_epivshypo_density_10-11Jun2021.jpg", width = 6, height = 4, units = "in",res = 300)
ggplot(subset(FCR_DVM_calcs_long, Taxa %in%taxa & grepl("density",metric,ignore.case = TRUE) & substrEnd(metric,7)!="density" &
                SampleID %in% c("F14Sep20_midn", "F15Sep20_noon")), aes(x=WaterColumn, y=value)) +
  geom_rect(data=subset(FCR_DVM_calcs_long, Taxa %in%taxa & Hour == 'midnight' &grepl("density",metric,ignore.case = TRUE)),
            aes(fill=Hour),xmin=-Inf ,xmax = Inf, ymin = -Inf, ymax = Inf, fill = 'black', alpha = 0.053,inherit.aes = FALSE) +
  geom_bar(aes(fill=Taxa, alpha=WaterColumn), stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() +
  geom_errorbar(aes(ymin=value-SE, ymax=value+SE), width=.2,position=position_dodge(.9)) + scale_alpha_manual(values = c(0.4, 1)) +
  facet_wrap(Hour~Taxa, scales= 'free',ncol=4, strip.position = "right", labeller=labeller(Hour=as_labeller(facet_labeller_bot),Taxa=as_labeller(facet_labeller_top))) +
  theme(strip.text.y = element_text(size = 11 ,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank()) + 
  scale_x_discrete(name="", labels=rep(c("epi","hypo"),10)) + guides(alpha=FALSE) +
  scale_fill_manual("Taxa",values=viridis(6),labels=c("Cladocera","Rotifera","Cyclopoida","Calanoida")) +
  geom_hline(yintercept=0) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
                                  axis.text.x=element_text(size=8,family="Times"), plot.title = element_text(hjust = 0.5)) + ylab("Density (individual/L)") +
  theme(legend.position = "bottom", legend.margin = margin(0, 1, 2, 1),plot.margin = margin(0,0,0,0.3,unit = "cm"), axis.text.y = element_text(size=13, family="Times"),
        legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),legend.text = element_text(size=10),legend.title = element_text(size=12))
#dev.off()

#---------------------------------------------------------------------------------------#
#also do same thing but with schindelr trap data

#add anoxic depth
zoop_repmeans_schind$anoxic_top <- ifelse(zoop_repmeans_schind$collect_date=="2020-09-11",2.5,
                                          ifelse(zoop_repmeans_schind$collect_date=="2020-09-15",2.7,4.3))

zoop_repmeans_schind$anoxic_bot <- ifelse(zoop_repmeans_schind$collect_date=="2020-09-11",5,
                                          ifelse(zoop_repmeans_schind$collect_date=="2020-09-15",5,10))


#new df to calculate hypo and epi zoop density 
schindler_mom <- data.frame("SampleID"=unique(substr(zoop_repmeans_schind$sample_ID,1,18)))

#select columns that I want to average
dens <- c(colnames(zoop_repmeans_schind)[c(8,12,17,22,27,32)])

#calculate epi as mean >= anoxic top #percent?
for (i in 1:length(dens)){
  schindler_mom[,paste0(dens,"epi_dens")[i]] <- mean(zoop_repmeans_schind[which(zoop_repmeans_schind$anoxic_top >= as.numeric(substrEnd(zoop_repmeans_schind$sample_ID,3))[i]), 
                                                                          paste0(dens)[i]])
}











#----------------------------------------------#
# Super simple plot of temp and DO during MSNs #
#----------------------------------------------#

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


ysi <- ysi[,c(1:7)]
ysi$DateTime <- as.Date(ysi$DateTime)
ysi <- ysi[ysi$DateTime=="2022-06-30" | ysi$DateTime=="2022-07-01",]


#jpeg("Figures/2021_Temp_O2_profile.jpg", width = 6, height = 5, units = "in",res = 300)
par(mfrow=c(1,1))
par(mar = c(4,0,0,0))
par(oma = c(2,4,4,4))

plot(ysi$Depth_m~ysi$DO_mgL, ylim=c(11,0), pch=16, type='o',col="blue",xlim=c(0,7.5), ylab="", xlab="",cex.axis=1.5, cex.lab=1.5)
par(new=TRUE)
plot(ysi$Depth_m~ysi$Temp_C,pch=16,type='o',col="red", yaxt='n',xaxt='n',xlab=" ",ylab=" ",ylim=c(11,0),xlim=c(8,28))
#axis(3,at=seq(round(min(CTD$Temp_C),-1),round(max(CTD$Temp_C),-1),length.out=6), cex.axis=1.5, cex.lab=1.5)
#text(11.5,0.5,"12 Aug 2020", cex=1.6)

mtext(expression(paste("Temperature (",degree,"C)")),side=3, line=2.4, cex=1.5, outer = TRUE)
mtext("Dissolved Oxygen (mg/L)",side=1, line=-1.8, cex=1.5, outer = TRUE)
mtext("          Depth (m)",side=2, line=2.4, cex=1.5, outer = TRUE)
legend("bottomright", legend=c("Dissolved Oxygen","Temperature"), col=c("blue", "red"), 
       cex=1.1, pch=16, box.lty=0,bg="transparent",xjust=1)
#dev.off()

#------------------------------------------------------------------------------#
#Look at schindler data because tows aren't super interesting

schindlers <- zoop[zoop$mesh_size_μm==61 & !is.na(zoop$mesh_size_μm) & zoop$site_no=="FCR_schind",]

#order depth by decreasing number
schindlers <- schindlers[with(schindlers,order(DepthOfTow_m)),]

#select the specific taxa
schindlers <- schindlers[,c(1:6,8,38,46,54,62)]

#wide to long
schindlers_long <- schindlers %>% gather(metric,value,Zooplankton_No.:RotiferaCount_n)

#jpeg("Figures/Schindler_density_vs_depth.jpg", width = 6, height = 5, units = "in",res = 300)
ggplot(data=schindlers_long,aes(x=value/30, y=DepthOfTow_m,color=collect_date)) + geom_point() +
  scale_y_reverse() + geom_path() + facet_grid(metric~site_no+collect_date)
#dev.off()

#summarize schindler_totalCount so one # per depth
Schindler_avgCount <- schindlers_long %>% group_by(site_no, DepthOfTow_m, collect_date, metric) %>% summarise(mean_num = mean(value))

#jpeg("Figures/AvgSchindler_density_vs_depth.jpg", width = 6, height = 5, units = "in",res = 300)
ggplot(data=Schindler_avgCount,aes(x=mean_num/30, y=DepthOfTow_m,color=collect_date)) + geom_point() +
  scale_y_reverse() + geom_path() + xlab("Zooplankton (#/L)") + ylab("Depth (m)")  + facet_wrap(metric~site_no, scales="free")
#dev.off()

