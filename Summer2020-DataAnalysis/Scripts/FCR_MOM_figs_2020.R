#FCR MOM figs 2020 - noon and midnight tows 

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
#read in zoop summary csv
zoop<- read.csv('SummaryStats/FCR_ZooplanktonSummary2020.csv',header = TRUE)

#only select FCR tows
zoop <- zoop[zoop$site_no=="FCR_50"& zoop$mesh_size_μm==80,]

#manually set hour to midnight or noon
zoop$Hour <- ifelse(substr(zoop$sample_ID,10,17)=="midnight", 00,12)

#pull rep # off as new column
zoop$rep <- ifelse(substrEnd(zoop$sample_ID,4)=="rep1" |substrEnd(zoop$sample_ID,4)=="rep2",
                   substrEnd(zoop$sample_ID,1),NA)

#drop reps from sampleIDs
zoop$sample_ID <- ifelse(substrEnd(zoop$sample_ID,4)=="rep1" |substrEnd(zoop$sample_ID,4)=="rep2",
                         substr(zoop$sample_ID,1,nchar(zoop$sample_ID)-5),zoop$sample_ID)

#Replace NAs with 0
zoop[20:204][is.na(zoop[20:204])] <- 0

#select noon/midnight pairs
zoop <- zoop[zoop$collect_date== "2020-06-28"|zoop$collect_date== "2020-06-29"|
               zoop$collect_date== "2020-09-10"|zoop$collect_date== "2020-09-11"|
               zoop$collect_date== "2020-09-14"|zoop$collect_date== "2020-09-15",]

#drop oxy samples(s)
zoop <- zoop[substrEnd(zoop$sample_ID,3)!="oxy",]

##### Create new df to combine reps over 24 hours
zoop.repmeans <- zoop %>% select(sample_ID,site_no,collect_date,Hour, Volume_L, Volume_unadj, proportional_vol, ZoopDensity_No.pL, OverallCount_n, TotalBiomass_ug,
                                 BiomassConcentration_ugpL,Cladocera_density_NopL, Cladocera_BiomassConcentration_ugpL, CladoceraCount_n, Cladocera_totalbiomass_ug, Cladocera_PercentOfTotal,
                                 Cyclopoida_density_NopL, Cyclopoida_BiomassConcentration_ugpL, CyclopoidaCount_n, Cyclopoida_totalbiomass_ug, Cyclopoida_PercentOfTotal,
                                 Rotifera_density_NopL,Rotifera_BiomassConcentration_ugpL, RotiferaCount_n, Rotifera_totalbiomass_ug, Rotifera_PercentOfTotal, 
                                 Calanoida_density_NopL, Calanoida_BiomassConcentration_ugpL, CalanoidaCount_n, Calanoida_PercentOfTotal, Calanoida_totalbiomass_ug) %>%
  group_by(sample_ID, site_no, collect_date, Hour) %>%
  summarise_at(vars(Volume_L:Calanoida_totalbiomass_ug,), funs(rep.mean=mean, rep.SE=stderr))



#----------------------------------------------------------------------#
#                       FCR noon + midnight                            #
#       creating new dfs for density and raw count/biomass data        #
#----------------------------------------------------------------------#

#two dfs for raw (no net efficiency calcs needed) and volume calculated values
FCR_taxa_DVM_raw<- zoop.repmeans[,c(1:4,6,7,which(substrEnd(colnames(zoop.repmeans),10)=="n_rep.mean"),
                                    which(substrEnd(colnames(zoop.repmeans),11)=="ug_rep.mean"),
                                    which(substrEnd(colnames(zoop.repmeans),8)=="n_rep.SE"),
                                    which(substrEnd(colnames(zoop.repmeans),9)=="ug_rep.SE"))]


#df for vol_calculated
FCR_taxa_DVM_vol_calculated <- zoop.repmeans[,c(1:4,8,35,which(substrEnd(colnames(zoop.repmeans),15)=="y_NopL_rep.mean"),
                                                which(substrEnd(colnames(zoop.repmeans),15)=="n_ugpL_rep.mean"),
                                                which(substrEnd(colnames(zoop.repmeans),13)=="n_ugpL_rep.SE"),
                                                which(substrEnd(colnames(zoop.repmeans),13)=="y_NopL_rep.SE"))]

#another one for percent calcs
FCR_DVM_percent<- zoop.repmeans[,c(1:4,which(substrEnd(colnames(zoop.repmeans),14)=="Total_rep.mean"),
                                   which(substrEnd(colnames(zoop.repmeans),12)=="Total_rep.SE"))]


#df for calculating epi vs hypo density/biomass
FCR_DVM_calcs_raw<- data.frame("SampleID"=unique(substr(zoop$sample_ID,1,13)))
FCR_DVM_calcs_vol_calc<- data.frame("SampleID"=unique(substr(zoop$sample_ID,1,13)))

#for loop to fill out epi vs hypo calcs RAW counts and ug
variables<- colnames(FCR_taxa_DVM_raw)[c(7:16)] #change range
for(i in 1:length(variables)){
  FCR_DVM_calcs_raw[,paste0(variables,"_epi")[i]]<- FCR_taxa_DVM_raw[substrEnd(FCR_taxa_DVM_raw$sample_ID,3)=="epi",paste0(variables)[i]]
  FCR_DVM_calcs_raw[,paste0(variables,"_hypo")[i]] <- FCR_taxa_DVM_raw[substrEnd(FCR_taxa_DVM_raw$sample_ID,3)!="epi",paste0(variables)[i]] - FCR_DVM_calcs_raw[,paste0(variables,"_epi")[i]]
}

#Net efficiencies from NetEfficiencyCalcs script
NetEfficiency2020 <- c(0.03491837, 0.05757148)

#hypo density and biomass calculated by subtracting epi raw zoop # from full zoop # and then dividing by the (full volume - epi volume) 
#NOTE: using epi density/L and biomass/L but calculating hypo using raw # and ug values. Note including the net efficiency to correct for the tows when using raw counts
column.names<- colnames(FCR_taxa_DVM_vol_calculated[,c(5,7:15)])
variables<- colnames(FCR_taxa_DVM_raw)[c(7:16)]
percent<- colnames(FCR_DVM_percent[,c(5:8)])
for(i in 1:length(variables)){
  FCR_DVM_calcs_vol_calc[,paste0(column.names,"_epi")[i]]<- FCR_taxa_DVM_vol_calculated[substrEnd(FCR_taxa_DVM_vol_calculated$sample_ID,3)=="epi",paste0(column.names)[i]]
  FCR_DVM_calcs_vol_calc[,paste0(column.names,"_hypo")[i]] <- (((1/FCR_taxa_DVM_raw[substrEnd(FCR_taxa_DVM_raw$sample_ID,3)!="epi" ,"proportional_vol_rep.mean"]) * FCR_taxa_DVM_raw[substrEnd(FCR_taxa_DVM_raw$sample_ID,3)!="epi" ,paste0(variables)[i]]) * (1/NetEfficiency2020[2])) - 
    ((1/FCR_taxa_DVM_raw[substrEnd(FCR_taxa_DVM_raw$sample_ID,3)=="epi", "proportional_vol_rep.mean"]) *  FCR_taxa_DVM_raw[substrEnd(FCR_taxa_DVM_raw$sample_ID,3)=="epi",paste0(variables)[i]]) / #both these net efficiencies are calculated from 2020 data (since these tows are <4m, net efficiency ~100%!)
    (FCR_taxa_DVM_raw[substrEnd(FCR_taxa_DVM_raw$sample_ID,3)!="epi" ,"Volume_unadj_rep.mean"] - FCR_taxa_DVM_raw[substrEnd(FCR_taxa_DVM_raw$sample_ID,3)=="epi", "Volume_unadj_rep.mean"])  
} #Note - not really sure why, but total density w/ net efficiency taken into account results in a million zoops in the hypo (and some of the other taxa seem unreasonably high)
#second note is that I was initially doing 1/proportional volume, which is right - see notebook from 250ct21


#percent density
density.percent<- colnames(FCR_taxa_DVM_vol_calculated[,c(7:10)]) 
for(i in 1:length(density.percent)){
  for(j in 1:length(unique(FCR_DVM_calcs_vol_calc$SampleID))){
    FCR_DVM_calcs_vol_calc[j,paste0(density.percent,"_epi_percent_density")[i]]<- (FCR_DVM_calcs_vol_calc[j,paste0(density.percent,"_epi")][i]/ sum(FCR_DVM_calcs_vol_calc[j,paste0(density.percent,"_epi")[i]],FCR_DVM_calcs_vol_calc[j,paste0(density.percent,"_hypo")[i]])) *100
    FCR_DVM_calcs_vol_calc[j,paste0(density.percent,"_hypo_percent_density")[i]]<- (FCR_DVM_calcs_vol_calc[j,paste0(density.percent,"_hypo")][i]/ sum(FCR_DVM_calcs_vol_calc[j,paste0(density.percent,"_epi")][i],FCR_DVM_calcs_vol_calc[j,paste0(density.percent,"_hypo")][i])) * 100
    
  }       
}


# now calculate SE dfs

#need reps from zoop df to calculate SE for raw and vol_calculated (but don't have reps here...)
matchingcols <- match(substr(colnames(FCR_taxa_DVM_raw[,c(1:4,6:16)]),1,14),substr(colnames(zoop.repmeans),1,14))
DVM_samples_raw<- zoop[,unique(matchingcols)]

matchingcols <- match(substr(colnames(FCR_taxa_DVM_vol_calculated[,c(1:4,7:10)]),1,14),substr(colnames(zoop.repmeans),1,14))
DVM_samples_dens<- zoop[,unique(matchingcols)] 

#separate full vs epi samples
FullSamples <- unique(DVM_samples_raw$sample_ID[substrEnd(DVM_samples_raw$sample_ID,3)!="epi"])
EpiSamples<- unique(DVM_samples_raw$sample_ID[substrEnd(DVM_samples_raw$sample_ID,3)=="epi"])

#initialize df
FCR_DVM_calcs_SE<- data.frame("SampleID"=unique(substr(zoop.repmeans$sample_ID,1,13)))

#calculate hypo SE 
SEonly<- colnames(FCR_taxa_DVM_raw)[7:16]
SEonly <- substr(SEonly,1,nchar(SEonly)-9)

Percentdens <- colnames(FCR_taxa_DVM_vol_calculated)[7:10]
Percentdens <- substr(Percentdens,1,nchar(Percentdens)-9)

for(i in 1:length(SEonly)){
  FCR_DVM_calcs_SE[,paste0(column.names,"_epi_SE")[i]] <- FCR_taxa_DVM_vol_calculated[substrEnd(FCR_taxa_DVM_vol_calculated$sample_ID,3)=="epi",substrEnd(colnames(FCR_taxa_DVM_vol_calculated),2)=="SE"][i]
  FCR_DVM_calcs_SE[1,paste0(column.names,"_hypo_SE")[i]]<- SE.diffMean(DVM_samples_raw[DVM_samples_raw$sample_ID==FullSamples[1],paste0(SEonly)[i]],
                                                                       DVM_samples_raw[DVM_samples_raw$sample_ID==EpiSamples[1],paste0(SEonly)[i]])
  FCR_DVM_calcs_SE[2,paste0(column.names,"_hypo_SE")[i]]<- SE.diffMean(DVM_samples_raw[DVM_samples_raw$sample_ID==FullSamples[2],paste0(SEonly)[i]],
                                                                       DVM_samples_raw[DVM_samples_raw$sample_ID==EpiSamples[2],paste0(SEonly)[i]])
} #Note: just repeating hypo SE calc for each date because I can't figure out the for loop across rows and columns...


for (j in 1:length(Percentdens)){
  FCR_DVM_calcs_SE[,paste0(Percentdens,"_epi_percent_density_SE")[j]]<- FCR_taxa_DVM_vol_calculated[substrEnd(FCR_taxa_DVM_vol_calculated$sample_ID,3)=="epi",substrEnd(colnames(FCR_taxa_DVM_vol_calculated),11)=="NopL_rep.SE"][j]
  FCR_DVM_calcs_SE[1,paste0(Percentdens,"_hypo_percent_density_SE")[j]] <- SE.diffMean(DVM_samples_dens[DVM_samples_dens$sample_ID==FullSamples[1],paste0(Percentdens)[j]],
                                                                                       DVM_samples_dens[DVM_samples_dens$sample_ID==EpiSamples[1],paste0(Percentdens)[j]])
  FCR_DVM_calcs_SE[2,paste0(Percentdens,"_hypo_percent_density_SE")[j]] <- SE.diffMean(DVM_samples_dens[DVM_samples_dens$sample_ID==FullSamples[2],paste0(Percentdens)[j]],
                                                                                       DVM_samples_dens[DVM_samples_dens$sample_ID==EpiSamples[2],paste0(Percentdens)[j]])
}

#------------------------------------------------------------------------------#
#wide to long df
FCR_DVM_calcs_long <-   FCR_DVM_calcs_vol_calc %>%
  gather(metric,value, ZoopDensity_No.pL_rep.mean_epi:Calanoida_density_NopL_rep.mean_hypo_percent_density) 

FCR_DVM_calcs_SE_long <-  FCR_DVM_calcs_SE %>%
  gather(metric,SE,ZoopDensity_No.pL_rep.mean_epi_SE:Calanoida_density_NopL_hypo_percent_density_SE) 

#add the SE column from df2 to df1 for combined df
FCR_DVM_calcs_long$SE <- FCR_DVM_calcs_SE_long$SE

#add watercolumn and hour columns
FCR_DVM_calcs_long$WaterColumn <- ifelse(substrEnd(FCR_DVM_calcs_long$metric,3)=="epi" |substrEnd(FCR_DVM_calcs_long$metric,19)=="epi_percent_density" ,"epilimnion","hypolimnion")
FCR_DVM_calcs_long$Hour <- ifelse(substrEnd(FCR_DVM_calcs_long$SampleID,4)=="midn","midnight","noon")
FCR_DVM_calcs_long$Taxa <- substr(FCR_DVM_calcs_long$metric,1,9)

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
FCR_DVM_calcs_long$Taxa <- factor(FCR_DVM_calcs_long$Taxa, levels = c("Cladocera","Rotifera_","Cyclopoid","Calanoida"))

#sample ids
sep1<- c("F10Sep20_midn", "F11Sep20_noon")
sep2 <- c("F14Sep20_midn", "F15Sep20_noon")
noonjun<- c("F28Jun20_noon", "F29Jun20_noon")

#jpeg("Figures/FCR_epivshypo_density_10-11Sep2020.jpg", width = 6, height = 4, units = "in",res = 300)
ggplot(subset(FCR_DVM_calcs_long, grepl("density",metric,ignore.case = TRUE) & substrEnd(metric,7)!="density"& SampleID %in% sep1), aes(x=WaterColumn, y=value)) +
  geom_rect(data=subset(FCR_DVM_calcs_long,Hour == 'midnight' &grepl("density",metric,ignore.case = TRUE)),
            aes(fill=Hour),xmin=-Inf ,xmax = Inf, ymin = -Inf, ymax = Inf, fill = 'black', alpha = 0.053,inherit.aes = FALSE) +
  geom_bar(aes(fill=Taxa, alpha=WaterColumn), stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() +
  geom_errorbar(aes(ymin=value-SE, ymax=value+SE), width=.2,position=position_dodge(.9)) + scale_alpha_manual(values = c(0.4, 1)) +
  facet_wrap(Hour~Taxa, scales= 'free',ncol=4, strip.position = "right", labeller=labeller(Hour=as_labeller(facet_labeller_bot),Taxa=as_labeller(facet_labeller_top))) +
  theme(strip.text.y = element_text(size = 11 ,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank()) + 
  scale_x_discrete(name="", labels=rep(c("epi","hypo"),10)) + labs(title="10-11 Sep 2020") + guides(alpha=FALSE) +
  scale_fill_manual("Taxa",values=viridis(6),labels=c("Cladocera","Rotifera","Cyclopoida","Calanoida")) +
  geom_hline(yintercept=0) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
                                  axis.text.x=element_text(size=8,family="Times"), plot.title = element_text(hjust = 0.5)) + ylab("Density (individual/L)") +
  theme(legend.position = "bottom", legend.margin = margin(0, 1, 2, 1),plot.margin = margin(0,0,0,0.3,unit = "cm"), axis.text.y = element_text(size=13, family="Times"),
        legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),legend.text = element_text(size=10),legend.title = element_text(size=12))
#dev.off()

#jpeg("Figures/FCR_epivshypo_biomass_10-11Sep2020.jpg", width = 6, height = 4, units = "in",res = 300)
ggplot(subset(FCR_DVM_calcs_long, grepl("biomass",metric,ignore.case = TRUE) & substrEnd(metric,7)!="density" & SampleID %in% sep1), aes(x=WaterColumn, y=value)) +
  geom_rect(data=subset(FCR_DVM_calcs_long,Hour == 'midnight' &grepl("biomass",metric,ignore.case = TRUE)),
            aes(fill=Hour),xmin=-Inf ,xmax = Inf, ymin = -Inf, ymax = Inf, fill = 'black', alpha = 0.053,inherit.aes = FALSE) +
  geom_bar(aes(fill=Taxa, alpha=WaterColumn), stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() +
  geom_errorbar(aes(ymin=value-SE, ymax=value+SE), width=.2,position=position_dodge(.9)) + scale_alpha_manual(values = c(0.4, 1)) +
  facet_wrap(Hour~Taxa, scales= 'free',ncol=4, strip.position = "right", labeller=labeller(Hour=as_labeller(facet_labeller_bot),Taxa=as_labeller(facet_labeller_top))) +
  theme(strip.text.y = element_text(size = 11 ,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank()) + 
  scale_x_discrete(name="", labels=rep(c("epi","hypo"),10)) + labs(title="10-11 Sep 2020") + guides(alpha=FALSE) +
  scale_fill_manual("Taxa",values=viridis(6),labels=c("Cladocera","Rotifera","Cyclopoida","Calanoida")) +
  geom_hline(yintercept=0) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
                                  axis.text.x=element_text(size=8,family="Times"), plot.title = element_text(hjust = 0.5)) + ylab(expression(paste("Biomass (",mu,"g/L)"))) +
  theme(legend.position = "bottom", legend.margin = margin(0, 1, 2, 1),plot.margin = margin(0,0,0,0.3,unit = "cm"), axis.text.y = element_text(size=13, family="Times"),
        legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),legend.text = element_text(size=10),legend.title = element_text(size=12))
#dev.off()

#jpeg("Figures/FCR_epivshypo_percent_dens_10-11Sep2020.jpg", width = 6, height = 4, units = "in",res = 300)
ggplot(subset(FCR_DVM_calcs_long, grepl("density",metric,ignore.case = TRUE) & substrEnd(metric,7)=="density"& SampleID %in% sep1), aes(x=WaterColumn, y=value)) +
  geom_rect(data=subset(FCR_DVM_calcs_long,Hour == 'midnight' &grepl("density",metric,ignore.case = TRUE)),
            aes(fill=Hour),xmin=-Inf ,xmax = Inf, ymin = -Inf, ymax = Inf, fill = 'black', alpha = 0.053,inherit.aes = FALSE) +
  geom_bar(aes(fill=Taxa, alpha=WaterColumn), stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() +
  geom_errorbar(aes(ymin=value-SE, ymax=value+SE), width=.2,position=position_dodge(.9)) + scale_alpha_manual(values = c(0.4, 1)) + ylim(-20,150) +
  facet_wrap(Hour~Taxa, scales= 'free',ncol=4, strip.position = "right", labeller=labeller(Hour=as_labeller(facet_labeller_bot),Taxa=as_labeller(facet_labeller_top))) +
  theme(strip.text.y = element_text(size = 11 ,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank()) + 
  scale_x_discrete(name="", labels=rep(c("epi","hypo"),10)) + labs(title="10-11 Sep 2020") + guides(alpha=FALSE) +
  scale_fill_manual("Taxa",values=viridis(6),labels=c("Cladocera","Rotifera","Cyclopoida","Calanoida")) +
  geom_hline(yintercept=0) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
                                  axis.text.x=element_text(size=8,family="Times"), plot.title = element_text(hjust = 0.5)) + ylab("Percent Density") +
  theme(legend.position = "bottom", legend.margin = margin(0, 1, 2, 1),plot.margin = margin(0,0,0,0.3,unit = "cm"), axis.text.y = element_text(size=13, family="Times"),
        legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),legend.text = element_text(size=10),legend.title = element_text(size=12))
#dev.off()


#------------------------------------------------------------------------------#
#Look at schindler data because tows aren't super interesting

schindlers <- zoop[zoop$mesh_size_μm==61 & !is.na(zoop$mesh_size_μm) & zoop$site_no=="FCR_schind",]

#order depth by decreasing number
schindlers <- schindlers[with(schindlers,order(DepthOfTow_m)),]

#jpeg("Figures/Schindler_density_vs_depth.jpg", width = 6, height = 5, units = "in",res = 300)
ggplot(data=schindlers,aes(x=Zooplankton_No./30, y=DepthOfTow_m,color=collect_date)) + geom_point() +
  scale_y_reverse() + geom_path() + facet_grid(~site_no+collect_date)
#dev.off()













#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#extra code that needs to be cleaned up or deleted


#initialize df
FCR.DVM.calcs<- data.frame("Hour"=unique(FCR_pelagic_DVM$Hour))

#for loop to fill out epi vs hypo calcs 
variables<- colnames(FCR_pelagic_DVM)[c(5:18)] #change range
for(i in 1:length(variables)){
  FCR.DVM.calcs[,paste0(variables,"_epi")[i]]<- FCR_pelagic_DVM[substrEnd(FCR_pelagic_DVM$sample_ID,4)!="9.0m" & substrEnd(FCR_pelagic_DVM$sample_ID,4)!="8.0m",paste0(variables)[i]]
  FCR.DVM.calcs[,paste0(variables,"_hypo")[i]] <- FCR_pelagic_DVM[substrEnd(FCR_pelagic_DVM$sample_ID,4)=="9.0m" | substrEnd(FCR_pelagic_DVM$sample_ID,4)=="8.0m",paste0(variables)[i]] - FCR.DVM.calcs[,paste0(variables,"_epi")[i]]
}

#wide to long df
FCR.DVM.calcs.long <-  FCR.DVM.calcs %>%
  gather(metric,value, ZoopDensity_No.pL_rep.mean_epi:Calanoida_BiomassConcentration_ugpL_rep.mean_hypo) %>%
  mutate(DateTime = strftime(Hour, "%m-%d-%Y %H:%M"))

#add watercolumn and hour columns
FCR.DVM.calcs.long$WaterColumn <- ifelse(substrEnd(FCR.DVM.calcs.long$metric,3)=="epi","epilimnion","hypolimnion")
FCR.DVM.calcs.long$Hour <- ifelse(substrEnd(FCR.DVM.calcs.long$DateTime,5)=="12:00","noon","midnight")
FCR.DVM.calcs.long$Taxa <- substr(FCR.DVM.calcs.long$metric,1,9)

#unique dates for for loop
dates<- unique(FCR.DVM.calcs.long$DateTime)

for(d in 1:length(dates)){
  #jpeg("Figures/FCR_epivshypo_density.jpg", width = 6, height = 5, units = "in",res = 300)
  plot<- ggplot(subset(FCR.DVM.calcs.long, grepl("density",metric,ignore.case = TRUE) & DateTime %in%
                         unique(DateTime)[d]), aes(x=WaterColumn, y=value)) +
    geom_bar(aes(fill=Taxa, alpha=WaterColumn), stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() +
    facet_wrap(Hour~Taxa, scales= 'free',ncol=5, strip.position = "right", labeller=as_labeller(function(variable,value){c("","","","","Noon")})) +
    theme(strip.text.y = element_text(size = 11 ,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank()) + 
    scale_x_discrete(name="", labels=rep(c("epi","hypo"),10)) + labs(title=paste0(dates[d])) + 
    scale_fill_manual("Taxa",values=viridis(6),labels=c("Cladocera","Cyclopoida","Calanoida","Rotifera","Total biomass")) +
    geom_hline(yintercept=0) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
                                    axis.text.x=element_text(size=9,family="Times"), plot.title = element_text(hjust = 0.5)) + ylab("Density (indiviual/L)") +
    theme(legend.position = "bottom", legend.margin = margin(0, 1, 2, 1),plot.margin = margin(0,0,0,0.3,unit = "cm"), axis.text.y = element_text(size=13, family="Times"),
          legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),legend.text = element_text(size=10),legend.title = element_text(size=12))
  print(plot)
  #dev.off()
}

#reorder taxa levels for figure
FCR.DVM.calcs.long$Taxa<-factor(FCR.DVM.calcs.long$Taxa,levels=c("Cladocera","Cyclopoid","Calanoida","Rotifera_","BiomassCo","ZoopDensi"))

for(d in 1:length(dates)){
  #jpeg("Figures/FCR_epivshypo_biomass.jpg", width = 6, height = 5, units = "in",res = 300)
  plot<- ggplot(subset(FCR.DVM.calcs.long, grepl("biomass",metric,ignore.case = TRUE) & DateTime %in%
                         unique(DateTime)[d]), aes(x=WaterColumn, y=value)) +
    geom_bar(aes(fill=Taxa, alpha=WaterColumn), stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() +
    facet_wrap(Hour~Taxa, scales= 'free',ncol=5, strip.position = "right", labeller=as_labeller(function(variable,value){c("","","","","Noon")})) +
    theme(strip.text.y = element_text(size = 11 ,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank()) + 
    scale_x_discrete(name="", labels=rep(c("epi","hypo"),10)) + labs(title=paste0(dates[d])) + 
    scale_fill_manual("Taxa",values=viridis(6),labels=c("Cladocera","Cyclopoida","Calanoida","Rotifera","Total biomass")) +
    geom_hline(yintercept=0) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
                                    axis.text.x=element_text(size=9,family="Times"), plot.title = element_text(hjust = 0.5)) + ylab(expression(paste("Biomass (",mu,"g/L)"))) +
    theme(legend.position = "bottom", legend.margin = margin(0, 1, 2, 1),plot.margin = margin(0,0,0,0.3,unit = "cm"), axis.text.y = element_text(size=13, family="Times"),
          legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),legend.text = element_text(size=10),legend.title = element_text(size=12))
  print(plot)
  #dev.off()
}




