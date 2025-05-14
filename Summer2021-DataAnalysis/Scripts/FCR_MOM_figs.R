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
  summarise_at(vars(Volume_L:Copepoda_totalbiomass_ug,), list(rep.mean=mean, rep.SE=stderr))

#separate tow vs schindler data
zoop_repmeans_tows <- zoop.repmeans[zoop.repmeans$site_no=="FCR_50",]
zoop_repmeans_schind <- zoop.repmeans[zoop.repmeans$site_no=="FCR_schind",]


#two dfs for raw (no net efficiency calcs needed) and volume calculated values (vol_unadj and prop_vol)
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
} #this isn't corrected for the net efficiency so hypo is almost always negative (I think this df is not useful, so flagging for now to potentially delete later)

#mean net efficiency for fcr 2020, 2021, and 2022; note that I'm assuming that neteff is 100% for epi tows (>95% in 2021 and >100% in 2020)
netefficiency <- mean(c(0.05279169, 0.05513654))

#hypo density and biomass calculated by subtracting epi raw zoop # from full zoop # and then dividing by the (full volume - epi volume) 
#NOTE: using epi density/L and biomass/L but calculating hypo using raw # and ug values. AND including the net efficiency to correct for the tows when using raw counts
column.names<- colnames(FCR_DVM_vol_calculated[,c(5:16)])
variables<- colnames(FCR_taxa_DVM_raw)[c(7:18)]
percent<- colnames(FCR_DVM_percent[,c(7:11)])
for(i in 1:length(variables)){
  FCR_DVM_calcs_vol_calc[,paste0(column.names,"_epi")[i]]<- FCR_DVM_vol_calculated[substrEnd(FCR_DVM_vol_calculated$sample_ID,3)=="epi",paste0(column.names)[i]]
  FCR_DVM_calcs_vol_calc[,paste0(column.names,"_hypo")[i]] <- (((FCR_taxa_DVM_raw[substrEnd(FCR_taxa_DVM_raw$sample_ID,3)!="epi" ,"proportional_vol_rep.mean"]) * FCR_taxa_DVM_raw[substrEnd(FCR_taxa_DVM_raw$sample_ID,3)!="epi" ,paste0(variables)[i]]) * (1/netefficiency)) - 
                                                               ((FCR_taxa_DVM_raw[substrEnd(FCR_taxa_DVM_raw$sample_ID,3)=="epi", "proportional_vol_rep.mean"]) *  FCR_taxa_DVM_raw[substrEnd(FCR_taxa_DVM_raw$sample_ID,3)=="epi",paste0(variables)[i]])/ 
                                                               (FCR_taxa_DVM_raw[substrEnd(FCR_taxa_DVM_raw$sample_ID,3)!="epi" ,"Volume_unadj_rep.mean"] - FCR_taxa_DVM_raw[substrEnd(FCR_taxa_DVM_raw$sample_ID,3)=="epi", "Volume_unadj_rep.mean"])                                                                     
} 


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
FullSamples <-data.frame("SampleID"=unique(DVM_samples_raw$sample_ID[substrEnd(DVM_samples_raw$sample_ID,3)!="epi"])) %>% arrange(SampleID,unique(substr(FCR_taxa_DVM_raw$sample_ID,1,13)))
EpiSamples<- data.frame("SampleID"=unique(DVM_samples_raw$sample_ID[substrEnd(DVM_samples_raw$sample_ID,3)=="epi"])) %>% arrange(SampleID,unique(substr(FCR_taxa_DVM_raw$sample_ID,1,13)))

#initialize df
FCR_DVM_calcs_SE<- data.frame("SampleID"=unique(substr(fcr_zoops_mom_tows$sample_ID,1,13))) %>% arrange(SampleID,unique(substr(FCR_taxa_DVM_raw$sample_ID,1,13)))

#calculate hypo SE 
SEonly<- colnames(FCR_taxa_DVM_raw)[7:18]
SEonly <- substr(SEonly,1,nchar(SEonly)-9)

Percentdens <- colnames(FCR_DVM_vol_calculated)[6:10]
Percentdens <- substr(Percentdens,1,nchar(Percentdens)-9)

for(i in 1:length(SEonly)){ 
  for(j in 1:nrow(EpiSamples)){ #hypo only works if there are reps for both full and epi tows (n=3)
  FCR_DVM_calcs_SE[,paste0(column.names,"_epi_SE")[i]] <- FCR_DVM_vol_calculated[substrEnd(FCR_DVM_vol_calculated$sample_ID,3)=="epi",substrEnd(colnames(FCR_DVM_vol_calculated),2)=="SE"][i]
  FCR_DVM_calcs_SE[j,paste0(column.names,"_hypo_SE")[i]]<- SE.diffMean(DVM_samples_raw[DVM_samples_raw$sample_ID==FullSamples[[1]][j],paste0(SEonly)[i]],
                                                                       DVM_samples_raw[DVM_samples_raw$sample_ID==EpiSamples[[1]][j],paste0(SEonly)[i]])
  }
} #double check but I think this works



for (i in 1:length(Percentdens)){
  for(j in 1:nrow(EpiSamples)){
  FCR_DVM_calcs_SE[,paste0(Percentdens,"_epi_percent_density_SE")[i]]<- FCR_DVM_vol_calculated[substrEnd(FCR_DVM_vol_calculated$sample_ID,3)=="epi",substrEnd(colnames(FCR_DVM_vol_calculated),11)=="NopL_rep.SE"][i]
  FCR_DVM_calcs_SE[j,paste0(Percentdens,"_hypo_percent_density_SE")[i]] <- SE.diffMean(DVM_samples_dens[DVM_samples_dens$sample_ID==FullSamples[[1]][j],paste0(Percentdens)[i]],
                                                                                       DVM_samples_dens[DVM_samples_dens$sample_ID==EpiSamples[[1]][j],paste0(Percentdens)[i]])
  }
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

#drop days that don't have a noon/midnight pair
FCR_DVM_calcs_long <- FCR_DVM_calcs_long %>% 
  filter(!SampleID %in% c("F08Jun20_noon","F20Jul20_noon","F27Jul20_noon","F30Sep20_noon"))

#pull out date
FCR_DVM_calcs_long$DateTime <- as.Date(substr(FCR_DVM_calcs_long$SampleID,2,8), format="%d%b%y")

#replace NAN with 0
FCR_DVM_calcs_long$value[is.nan(FCR_DVM_calcs_long$value)] <- 0

#Export FCR MOM DVM stats
#write.csv(FCR_DVM_calcs_long, "Summer2021-DataAnalysis/SummaryStats/FCR_DVM_2020-2022_zoops.csv")


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

#add a column for campaign #
FCR_DVM_calcs_long$campaign_num <- ifelse(FCR_DVM_calcs_long$SampleID %in% c("F28Jun20_midn", "F29Jun20_noon"), 1,
                                          ifelse(FCR_DVM_calcs_long$SampleID %in% c("F10Sep20_midn", "F11Sep20_noon"), 2,
                                                 ifelse(FCR_DVM_calcs_long$SampleID %in% c("F14Sep20_midn", "F15Sep20_noon"), 3,
                                                        ifelse(FCR_DVM_calcs_long$SampleID %in% c("F10Jun21_midn", "F10Jun21_noon"), 4, 5))))


date_list <- c("28-29Jun2020", "10-11Sep2020", "14-15Sep2020", "10-11Jun2021", "30Jun-01Jul2022")

for(i in 1:length(date_list)){ 
plot <- ggplot(subset(FCR_DVM_calcs_long, Taxa %in%taxa & grepl("density",metric,ignore.case = TRUE) & 
                        substrEnd(metric,7)!="density" &
                      campaign_num== c(1:5)[i]), aes(x=WaterColumn, y=value)) +
              geom_rect(data=subset(FCR_DVM_calcs_long, Taxa %in%taxa & 
                                      Hour == 'midnight' &grepl("density",metric,ignore.case = TRUE) & 
                                      substrEnd(metric,7)!="density"),
                        aes(fill=Hour),xmin=-Inf ,xmax = Inf, ymin = -Inf, ymax = Inf, fill = 'gray', 
                        alpha = 0.053,inherit.aes = FALSE) +
              geom_bar(aes(fill=Taxa, alpha=WaterColumn), stat="identity", 
                       position=position_dodge(),show.legend = TRUE) + theme_bw() +
              geom_errorbar(aes(ymin=value-SE, ymax=value+SE), width=.2,
                            position=position_dodge(.9)) + scale_alpha_manual(values = c(0.4, 1)) +
              facet_wrap(Hour~Taxa, ncol=4, strip.position = "right", 
                         labeller=labeller(Hour=as_labeller(facet_labeller_bot),
                                           Taxa=as_labeller(facet_labeller_top))) +
              scale_x_discrete(name="", labels=rep(c("epi","hypo"),10)) + 
              guides(alpha=FALSE) + ggtitle(date_list[i]) +
              scale_fill_manual("Taxa",values=viridis(6),
                                labels=c("Cladocera","Rotifera","Cyclopoida","Calanoida")) +
              geom_hline(yintercept=0) +
              ylab("Density (individuals/L)") +
              theme(legend.position = "bottom", legend.margin = margin(0, 1, 2, 1),
                    plot.margin = margin(0,0,0,0.3,unit = "cm"), 
                    axis.text.y = element_text(size=13, family="Times"),
                    legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),
                    legend.text = element_text(size=6),legend.title = element_text(size=6),
                    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
                    axis.text.x=element_text(size=8,family="Times"), 
                    plot.title = element_text(hjust = 0.5, size = 8),
                    strip.text.y = element_text(size = 11 ,margin = margin(0, -0.01, 0,0.1, "cm")),
                    strip.background = element_blank(),
                    legend.key.size = unit(0.4, 'cm'))
print(plot)
#ggsave(paste0(getwd(),"/Summer2021-DataAnalysis/Figures/fcr_mom/FCR_dvm_dens_", date_list[i],".jpg"), width=4, height=3)

}


#plot taxa on x and color by epi vs hypo
for(i in 1:length(date_list)){ 
  plot <- ggplot(subset(FCR_DVM_calcs_long, Taxa %in%taxa & grepl("density",metric,ignore.case = TRUE) & 
                          substrEnd(metric,7)!="density" &
                          campaign_num== c(1:5)[i]), aes(x=Taxa, y=value, fill=WaterColumn)) +
    geom_rect(data=subset(FCR_DVM_calcs_long, Taxa %in%taxa & 
                            Hour == 'midnight' &grepl("density",metric,ignore.case = TRUE) & 
                            substrEnd(metric,7)!="density"),
              aes(fill=Hour),xmin=-Inf ,xmax = Inf, ymin = -Inf, ymax = Inf, fill = 'gray', 
              alpha = 0.053,inherit.aes = FALSE) +
    geom_bar(aes(fill=WaterColumn), stat="identity", 
             position=position_dodge(),show.legend = TRUE) + theme_bw() +
    geom_errorbar(aes(ymin=value-SE, ymax=value+SE), width=.2,
                  position=position_dodge(.9)) + scale_alpha_manual(values = c(0.4, 1)) +
    facet_grid(~Hour) + xlab("") +
    ggtitle(date_list[i]) +
    scale_fill_manual("",values=viridis(2)) +
    geom_hline(yintercept=0) +
    ylab("Density (individuals/L)") +
    theme(legend.position = "bottom", legend.margin = margin(0, 1, 2, 1),
          plot.margin = margin(0,0,0,0.3,unit = "cm"), 
          axis.text.y = element_text(size=13, family="Times"),
          legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),
          legend.text = element_text(size=6),legend.title = element_text(size=6),
          panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
          axis.text.x=element_text(vjust=0.5,angle=40,size=8,family="Times"), 
          plot.title = element_text(hjust = 0.5, size = 8),
          strip.text.y = element_text(size = 11 ,margin = margin(0, -0.01, 0,0.1, "cm")),
          strip.background = element_blank(),
          legend.key.size = unit(0.4, 'cm'))
  print(plot)
 # ggsave(paste0(getwd(),"/Summer2021-DataAnalysis/Figures/fcr_mom/FCR_epi_hypo_dens_", date_list[i],".jpg"), width=4, height=3)
}

#-------------------------------------------------------------------------------
####                         DVM METRICS                                    ####
#-------------------------------------------------------------------------------
#work with raw density and biomass to calculate metrics
FCR_DVM_calcs_long <- FCR_DVM_calcs_long[!grepl("percent",FCR_DVM_calcs_long$metric),]

#shorten metric name
FCR_DVM_calcs_long$metric <- ifelse(FCR_DVM_calcs_long$WaterColumn=="epilimnion", substr(FCR_DVM_calcs_long$metric,1,nchar(FCR_DVM_calcs_long$metric)-13),substr(FCR_DVM_calcs_long$metric,1,nchar(FCR_DVM_calcs_long$metric)-14))

DVM_proportion <-  plyr::ddply(FCR_DVM_calcs_long, c("metric", "campaign_num", "Hour","DateTime"), function(x) {
  data.frame(
    proportion_epi = x$value[x$WaterColumn=="epilimnion"] / sum(x$value)
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 

#DVM metrics for a single day --> DVM = (Depi / Depi + Dhypo)Night - (Depi / Depi + Dhypo)Day
DVM_metric_df <-  plyr::ddply(DVM_proportion, c("metric", "campaign_num"), function(x) {
  data.frame(
    DVM_metric = x$proportion_epi[x$Hour=="midnight"] - x$proportion_epi[x$Hour=="noon"]
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 

#replace NAN with 0 bc none of that taxa were found
DVM_metric_df$DVM_metric[is.nan(DVM_metric_df$DVM_metric)] <- 0

#---------------------------------------------------------------------------------------#
dens_list <- c("Cladocera_density_NopL","Copepoda_density_NopL","Rotifera_density_NopL")

metric_taxa <-c("Calanoida","Calanoida","Cladocera","Cladocera",
                "Copepoda","Copepoda","Cyclopoida", "Cyclopoida",
                "Rotifera", "Rotifera")
names(metric_taxa) <- c(unique(DVM_metric_df$metric))

#plot migration metrics
ggplot(subset(DVM_metric_df, grepl("density",metric, ignore.case=T) & 
                metric %in% dens_list), 
       aes(x=as.factor(campaign_num), y=DVM_metric, color=as.factor(campaign_num))) + 
  geom_point(position=position_dodge(.9)) + theme_bw() + geom_hline(yintercept = 0, linetype="dotted")+
  scale_x_discrete(breaks=c("1","2","3","4","5"),
                   labels=date_list) +  xlab("") +
  theme(text = element_text(size=8), axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), legend.key = element_blank(), 
        legend.key.height=unit(0.3,"line"), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        strip.background = element_rect(fill = "transparent"), 
        plot.margin = unit(c(0,3,0,0), 'lines'),
        legend.position = c(0.92,0.94), legend.spacing = unit(-0.5, 'cm'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.7,"line")) + guides(color="none") +
  scale_color_manual("",values=c("#008585","#89B199","#EFECBF","#DB9B5A","#C7522B"))+
  facet_wrap(~metric, labeller = labeller(metric=metric_taxa)) + ylab("Density DVM metric") +
  geom_text(aes(x=5.9, y=c(rep(0,13),0.02,-0.02), label=c(rep(NA,13),"Normal \nMigration", "Reverse \nMigration")), 
            hjust = 0, size = 3, color="black") + coord_cartesian(xlim = c(1, 5), clip = 'off')
#ggsave(file.path(getwd(),"Summer2021-DataAnalysis/Figures/fcr_mom/FCR_MSNs_DVM_metric_dens_3taxa.jpg"), width=5, height=4) 
#I think this counts as NO MIGRATION bc magnitude never exceeds 0.1 (bvr ranges from ~ 0.6 to -0.4)

#---------------------------------------------------------------------------------------#
#also do same thing but with Schindler trap data

#add anoxic depth
zoop_repmeans_schind$anoxic_top <- ifelse(zoop_repmeans_schind$collect_date=="2020-09-11",2.5,
                                          ifelse(zoop_repmeans_schind$collect_date=="2020-09-15",2.7,
                                                 ifelse(zoop_repmeans_schind$collect_date=="2022-07-01",2.3,4.2)))

zoop_repmeans_schind$anoxic_bot <- ifelse(zoop_repmeans_schind$collect_date=="2020-09-11",5,
                                          ifelse(zoop_repmeans_schind$collect_date=="2020-09-15",5,10)) #only have a mom for 2 days

#select columns that I want to average
dens <- c(colnames(zoop_repmeans_schind)[c(8,12,17,22,27,32)])

#new df to calculate hypo and epi zoop density 
schindler_mom <- data.frame("SampleID"=unique(substr(zoop_repmeans_schind$sample_ID,1,18)))

#list of anoxic depths for each date
anoxic_depths <- c(2.3,2.3,4.2,4.2,2.5,2.7) #match up with schindler_mom dates

#change 2022 midnight date so that all dates are unique for for loop below
zoop_repmeans_schind$collect_date[zoop_repmeans_schind$collect_date=="2022-07-01" & 
                                    zoop_repmeans_schind$Hour==0] <- "2022-06-30"

#calculate epi as mean >= anoxic top 
for (i in 1:length(dens)){
  for(j in 1:length(schindler_mom$SampleID)){
  schindler_mom[j,paste0(dens,"epi_dens")[i]] <- lapply(zoop_repmeans_schind[which(anoxic_depths[j]  >= as.numeric(substrEnd(zoop_repmeans_schind$sample_ID,3)) &
                                                        zoop_repmeans_schind$collect_date==unique(zoop_repmeans_schind$collect_date)[j]),paste0(dens[i])],mean, na.rm=T)   
  schindler_mom[j,paste0(dens,"hypo_dens")[i]] <- lapply(zoop_repmeans_schind[which(anoxic_depths[j]  < as.numeric(substrEnd(zoop_repmeans_schind$sample_ID,3)) &
                                                                                     zoop_repmeans_schind$collect_date==unique(zoop_repmeans_schind$collect_date)[j]),paste0(dens[i])],mean, na.rm=T)   
  schindler_mom[j,paste0(dens,"per_epid")[i]] <- (schindler_mom[j,paste0(dens,"epi_dens")[i]]/ sum(schindler_mom[j,paste0(dens,"epi_dens")[i]],schindler_mom[j,paste0(dens,"hypo_dens")[i]])) *100
  schindler_mom[j,paste0(dens,"per_hypod")[i]] <- (schindler_mom[j,paste0(dens,"hypo_dens")[i]]/ sum(schindler_mom[j,paste0(dens,"epi_dens")[i]],schindler_mom[j,paste0(dens,"hypo_dens")[i]])) *100
  }
} 

dens_npl <- colnames(schindler_mom[,which(substrEnd(colnames(schindler_mom),8)=="epi_dens"| substrEnd(colnames(schindler_mom),9)=="hypo_dens")])
per_dens <- colnames(schindler_mom[,which(substrEnd(colnames(schindler_mom),4)=="epid"| substrEnd(colnames(schindler_mom),5)=="hypod")])
    
#wide to long
temp_dens <-   schindler_mom %>% gather(metric.raw,value.raw, all_of(dens_npl))
temp_percentdens <- schindler_mom %>% gather(metric.per,value.per, all_of(per_dens))

#cut and paste to merge df
schindler_mom_long <- temp_dens[,c(1,14,15)]
#schindler_mom_long$metric.per <- temp_percentdens$metric.per
schindler_mom_long$value.per <- temp_percentdens$value.per

#add date and watercolumn cols
schindler_mom_long$WaterColumn <- ifelse(substrEnd(schindler_mom_long$metric.raw,8)=="epi_dens" | substrEnd(schindler_mom_long$metric.raw,8)=="per_epid","epilimnion","hypolimnion")
schindler_mom_long$CollectDate <- as.Date(substr(schindler_mom_long$SampleID,7,13), format="%d%b%y")

#change collect_date for 2022 midnight samples
schindler_mom_long$CollectDate[schindler_mom_long$SampleID=="F_pel_01Jul22_midn"] <- "2022-06-30"

#make metric names shorter
schindler_mom_long$metric.raw <- ifelse(schindler_mom_long$WaterColumn=="epilimnion",
                                    substr(schindler_mom_long$metric.raw,1,nchar(schindler_mom_long$metric.raw)-30),
                                    substr(schindler_mom_long$metric.raw,1,nchar(schindler_mom_long$metric.raw)-31))

#now make plots to look at epi vs hypo boxplots
ggplot(schindler_mom_long, aes(metric.raw, value.raw,fill=WaterColumn)) + ylab("Density (#/L)") + xlab("")+
  geom_bar(stat="identity", position=position_dodge()) + facet_wrap(~CollectDate, ncol=3)+
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = c(0.72,0.94),
        legend.background = element_blank(),legend.direction = "horizontal", panel.grid.minor = element_blank(), legend.key=element_rect(fill=NA),
        plot.margin = unit(c(0,0.05,0,0), "cm"),legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),
        legend.title = element_blank(),legend.text  = element_text(size = 4.5), panel.spacing=unit(0, "cm"),
        axis.text.x = element_text(angle=45, vjust = 0.7,size=6, hjust = 0.6), axis.text.y = element_text(size=6)) 
#ggsave(file.path(getwd(),"Summer2021-DataAnalysis/Figures/fcr_mom/FCR_schind_taxa_dens.jpg"), width=4, height=3) 

ggplot(schindler_mom_long, aes(metric.raw, value.per,fill=WaterColumn)) + ylab("Percent Density (%)") + xlab("")+
  geom_bar(stat="identity", position=position_dodge()) + facet_wrap(~CollectDate, ncol=3, scales="free_y")+
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = c(0.85,0.96),
        legend.background = element_blank(),legend.direction = "horizontal", panel.grid.minor = element_blank(), legend.key=element_rect(fill=NA),
        plot.margin = unit(c(0,0.05,0,0), "cm"),legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),
        legend.title = element_blank(),legend.text  = element_text(size = 4.5), panel.spacing=unit(0, "cm"),
        axis.text.x = element_text(angle=45, vjust = 0.7,size=6, hjust = 0.6), axis.text.y = element_text(size=6)) 
#ggsave(file.path(getwd(),"Summer2021-DataAnalysis/Figures/fcr_mom/FCR_schind_taxa_percent_dens.jpg"), width=4, height=3) 


#----------------------------------------------#
# Super simple plot of temp and DO during MSNs #
#----------------------------------------------#

inUrl2  <-  "https://pasta.lternet.edu/package/data/eml/edi/198/11/6e5a0344231de7fcebbe6dc2bed0a1c3" 
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

dates <- c("2020-06-28", "2020-06-29", "2020-09-10", "2020-09-11", "2020-09-14", "2020-09-15", 
           "2021-06-10", "2021-06-11", "2022-06-30", "2022-07-01")

ysi <- ysi[,c(1:7)]
ysi$DateTime <- as.Date(ysi$DateTime)
ysi <- ysi[ysi$DateTime %in% as.Date(dates) & ysi$Reservoir=="FCR",]

#order depths 
ysi <- ysi[order(ysi$Depth_m),]

for(i in 1:length(unique(ysi$DateTime))){
plot <- ggplot(data = subset(ysi, DateTime==unique(ysi$DateTime)[i]), aes(DO_mgL, Depth_m, col="DO (mg/L)")) + 
    geom_point() + geom_path() + theme_bw() + ylim(rev(range(ysi$Depth_m))) +
    geom_vline(xintercept=2, lty=2)+
    geom_point(data = subset(ysi, DateTime==unique(ysi$DateTime)[i]), aes(Temp_C, Depth_m, col="Temp (C)")) +
    geom_path(data = subset(ysi, DateTime==unique(ysi$DateTime)[i]), aes(Temp_C, Depth_m, col="Temp (C)")) +
    scale_x_continuous(sec.axis = sec_axis(~ . * 0.9, name = "Temperature")) +
    ggtitle(unique(ysi$DateTime)[i]) +
    scale_color_manual(values = c("black", "red")) +
    theme(text = element_text(size=10), axis.text = element_text(size=8, color="black"), 
          legend.position = c(0.72,0.94),
          legend.background = element_blank(),legend.direction = "horizontal", 
          panel.grid.minor = element_blank(), legend.key=element_rect(fill=NA),
          plot.margin = unit(c(0,0.05,0,0), "cm"),legend.key.size = unit(0.5, "lines"), 
          panel.grid.major = element_blank(), legend.title = element_blank(),
          legend.text  = element_text(size = 6), panel.spacing=unit(0, "cm"),
          axis.text.x = element_text(vjust = 0.5,size=8), axis.text.y = element_text(size=6)) 
print(plot)
#ggsave(paste0(getwd(),"/Figures/fcr_mom/2021_Temp_O2_",unique(ysi$DateTime)[i],"_ysi.jpg"))
}

#check the ctd data to confirm anoxia in 2021 and 2022
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/200/13/27ceda6bc7fdec2e7d79a6e4fe16ffdf" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")

ctd<- read.csv(infile1) |> 
  mutate(DateTime = as.Date(DateTime)) |> 
  filter(Reservoir =="FCR" & Depth_m > 0 &
           DateTime %in% c(as.Date("2020-09-10"), as.Date("2020-09-11"), as.Date("2020-09-14"),
                           as.Date("2020-09-15"), as.Date("2021-06-10"), as.Date("2021-06-11"),
                           as.Date("2022-06-30"), as.Date("2022-07-01")))


depths <- seq(0,10, by = 1)
newDepths <- depths
df.final.raw<- ctd %>%
  group_by(DateTime) %>%
  slice(which.min(abs(as.numeric(Depth_m) - depths[1]))) #Create a new dataframe
df.final.raw$Depth_m <- newDepths[1]
#loop through all depths and add the closest values to the final dataframe
for (i in 2:length(depths)){
  ctd_atThisDepth <- ctd %>%
    group_by(DateTime) %>%
    slice(which.min(abs(as.numeric(Depth_m) - depths[i])))
  ctd_atThisDepth <- ctd_atThisDepth %>%
    #only include if the measured depth is within 0.1 of the depth label
    filter(abs(Depth_m-newDepths[i])<0.1)
  ctd_atThisDepth$Depth_m <- newDepths[i]
  df.final.raw <- rbind(df.final.raw,ctd_atThisDepth)
}

for(i in 1:length(unique(df.final.raw$DateTime))){
  plot <- ggplot(data = subset(df.final.raw, DateTime==unique(df.final.raw$DateTime)[i]), aes(DO_mgL, Depth_m, col="DO (mg/L)")) + 
    geom_point() + geom_path() + theme_bw() + ylim(rev(range(df.final.raw$Depth_m))) +
    geom_vline(xintercept=2, lty=2)+
    geom_point(data = subset(df.final.raw, DateTime==unique(df.final.raw$DateTime)[i]), aes(Temp_C, Depth_m, col="Temp (C)")) +
    geom_path(data = subset(df.final.raw, DateTime==unique(df.final.raw$DateTime)[i]), aes(Temp_C, Depth_m, col="Temp (C)")) +
    scale_x_continuous(sec.axis = sec_axis(~ . * 0.9, name = "Temperature")) +
    ggtitle(unique(df.final.raw$DateTime)[i]) +
    scale_color_manual(values = c("black", "red")) +
    theme(text = element_text(size=10), axis.text = element_text(size=8, color="black"), 
          legend.position = c(0.72,0.94),
          legend.background = element_blank(),legend.direction = "horizontal", 
          panel.grid.minor = element_blank(), legend.key=element_rect(fill=NA),
          plot.margin = unit(c(0,0.05,0,0), "cm"),legend.key.size = unit(0.5, "lines"), 
          panel.grid.major = element_blank(), legend.title = element_blank(),
          legend.text  = element_text(size = 6), panel.spacing=unit(0, "cm"),
          axis.text.x = element_text(vjust = 0.5,size=8), axis.text.y = element_text(size=6)) 
  print(plot)
  #ggsave(paste0(getwd(),"/Figures/fcr_mom/2021_Temp_O2_",unique(df.final.raw$DateTime)[i],"_ctd.jpg"))
}

#date list with both ctd and ysi casts
dates_ctd_ysi <- c("2020-09-11","2020-09-15","2021-06-11")

#ysi + ctd DO on same plot
for(i in 1:length(dates_ctd_ysi)){
  plot <- ggplot(data = subset(df.final.raw, DateTime==dates_ctd_ysi[i]), aes(DO_mgL, Depth_m, col="CTD DO (mg/L)")) + 
    geom_point() + geom_path() + theme_bw() + ylim(rev(range(df.final.raw$Depth_m))) +
    geom_point(data = subset(ysi, DateTime==dates_ctd_ysi[i]), aes(DO_mgL, Depth_m, col="YSI DO (mg/L)")) +
    geom_path(data = subset(ysi, DateTime==dates_ctd_ysi[i]), aes(DO_mgL, Depth_m, col="YSI DO (mg/L)")) +
    geom_vline(xintercept=2, lty=2)+
    ggtitle(dates_ctd_ysi[i]) +
    scale_color_manual(values = c("black", "red")) +
    theme(text = element_text(size=10), axis.text = element_text(size=8, color="black"), 
          legend.position = c(0.72,0.94),
          legend.background = element_blank(),legend.direction = "horizontal", 
          panel.grid.minor = element_blank(), legend.key=element_rect(fill=NA),
          plot.margin = unit(c(0,0.05,0,0), "cm"),legend.key.size = unit(0.5, "lines"), 
          panel.grid.major = element_blank(), legend.title = element_blank(),
          legend.text  = element_text(size = 6), panel.spacing=unit(0, "cm"),
          axis.text.x = element_text(vjust = 0.5,size=8), axis.text.y = element_text(size=6)) 
  print(plot)
  #ggsave(paste0(getwd(),"/Figures/fcr_mom/2021_CTD_YSI_DO_comp_",dates_ctd_ysi[i],".jpg"))
}


#------------------------------------------------------------------------------#
#Look at schindler data because tows aren't super interesting

#select the specific taxa
schindlers <- zoop_repmeans_schind[,c(1:7, which(grepl("_n_",colnames(zoop_repmeans_schind))))]

#add depth column
schindlers$depth_m <- as.numeric(substrEnd(schindlers$sample_ID,3))

#wide to long
df1 <- schindlers %>% gather(metric,value,OverallCount_n_rep.mean:CopepodaCount_n_rep.mean)
df2 <- schindlers %>% gather(metric,value,OverallCount_n_rep.SE:CopepodaCount_n_rep.SE)

#cut and paste to merge df
schindlers_long <- df1[,c(1:7,14,15,16)]
#schindlers_long$metric.SE <- df2$metric
schindlers_long$value.SE <- df2$value

#make metric names shorter
schindlers_long$metric <-  substr(schindlers_long$metric,1,nchar(schindlers_long$metric)-16)

ggplot(data=schindlers_long,aes(x=value/30, y=depth_m,color=collect_date)) + geom_point() +
  scale_y_reverse() + geom_path() + facet_grid(metric~collect_date) + guides(color = FALSE) +
  theme(text = element_text(size=8), 
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5,size=6), 
        axis.text.y = element_text(size=6)) 
#ggsave(paste0(getwd(),"/Summer2021-DataAnalysis/Figures/fcr_mom/Schindler_density_vs_depth.jpg"))
#looks like rotifers are migrating!

#summarize schindler_totalCount so one # per depth
Schindler_avgCount <- schindlers_long %>% group_by(depth_m, collect_date, metric) %>% summarise(mean_num = mean(value))

ggplot(data=subset(Schindler_avgCount, !collect_date %in% c("2020-09-11","2020-09-15")),
       aes(x=mean_num/30, y=depth_m,color=collect_date)) + 
  geom_point() + scale_y_reverse() + geom_path() + xlab("Zooplankton (#/L)") + 
  ylab("Depth (m)")  + facet_wrap(~metric, scales="free") +
  theme(text = element_text(size=8), 
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5,size=6), 
        axis.text.y = element_text(size=6)) 
#ggsave(paste0(getwd(),"/Summer2021-DataAnalysis/Figures/fcr_mom/AvgSchindler_density_vs_depth_2021-2022.jpg"))
#wow this could actually be really cool! look at copepods - some are migrating and others aren't??

#-------------------------------------------------------------------------------#
#new dfs for multivariate stats
zoops_summary <- fcr_zoops_mom %>% select(sample_ID,site_no,collect_date,Hour, 
                                          Bosmina_density_NopL, Bosmina_BiomassConcentration_ugpL,
                                          Daphnia_density_NopL, Daphnia_BiomassConcentration_ugpL,
                                          Ceriodaphnia_density_NopL, Ceriodaphnia_BiomassConcentration_ugpL,
                                          Cyclopoida_density_NopL, Cyclopoida_BiomassConcentration_ugpL, 
                                          Keratella_density_NopL,Keratella_BiomassConcentration_ugpL, 
                                          Kellicottia_density_NopL,Kellicottia_BiomassConcentration_ugpL, 
                                          Ploima_density_NopL, Ploima_BiomassConcentration_ugpL,
                                          Gastropidae_density_NopL, Gastropidae_BiomassConcentration_ugpL,
                                          Collothecidae_density_NopL, Collothecidae_BiomassConcentration_ugpL,
                                          Conochilidae_density_NopL, Conochilidae_BiomassConcentration_ugpL,
                                          Synchaetidae_density_NopL, Synchaetidae_BiomassConcentration_ugpL,
                                          Trichocercidae_density_NopL, Trichocercidae_BiomassConcentration_ugpL,
                                          Lepadella_density_NopL, Lepadella_BiomassConcentration_ugpL,
                                          Monostyla_density_NopL, Monostyla_BiomassConcentration_ugpL,
                                          Lecane_density_NopL, Lecane_BiomassConcentration_ugpL,
                                          Calanoida_density_NopL, Calanoida_BiomassConcentration_ugpL, 
                                          nauplius_density_NopL, nauplius_BiomassConcentration_ugpL) |>
  mutate(depth = substrEnd(fcr_zoops_mom$sample_ID,3)) |>
  filter(site_no %in% c("FCR_schind"), collect_date %in% c(dates)) |>
  group_by(sample_ID, collect_date ,Hour, depth) |>
  summarise_at(vars(Bosmina_density_NopL:nauplius_BiomassConcentration_ugpL,), funs(rep.mean=mean))
  

zoops_summary_3groups <- fcr_zoops_mom %>% select(sample_ID,site_no,collect_date,Hour, 
                                          Cladocera_density_NopL, Cladocera_BiomassConcentration_ugpL,
                                          Copepoda_density_NopL, Copepoda_BiomassConcentration_ugpL,
                                          Rotifera_density_NopL, Rotifera_BiomassConcentration_ugpL) |>
  mutate(depth = substrEnd(fcr_zoops_mom$sample_ID,3)) |>
  filter(site_no %in% c("FCR_schind"), collect_date %in% c(dates)) |>
  group_by(sample_ID, collect_date ,Hour, depth) |>
  summarise_at(vars(Cladocera_density_NopL:Rotifera_BiomassConcentration_ugpL,), funs(rep.mean=mean))


#change 2022 midnight date so that all dates are unique for for loop below
zoops_summary$collect_date[zoops_summary$collect_date=="2022-07-01" & 
                             zoops_summary$Hour==0] <- "2022-06-30"

zoops_summary_3groups$collect_date[zoops_summary_3groups$collect_date=="2022-07-01" & 
                                     zoops_summary_3groups$Hour==0] <- "2022-06-30"

#add a column for campaign #
zoops_summary$campaign_num <- ifelse(zoops_summary$collect_date %in% c("2020-09-11"), 1,
                                          ifelse(zoops_summary$collect_date %in% c("2020-09-15"), 2,
                                                 ifelse(zoops_summary$collect_date %in% c("2021-06-10","2021-06-11"), 3, 4)))

zoops_summary_3groups$campaign_num <- ifelse(zoops_summary_3groups$collect_date %in% c("2020-09-11"), 1,
                                     ifelse(zoops_summary_3groups$collect_date %in% c("2020-09-15"), 2,
                                            ifelse(zoops_summary_3groups$collect_date %in% c("2021-06-10","2021-06-11"), 3, 4)))


#create dfs for multivariate zoop script (depth specific)
schind_dens <- zoops_summary[ ,c(1:4,19, which(grepl("density",colnames(zoops_summary))))]
schind_dens_3groups <- zoops_summary_3groups[ ,c(1:4,11, which(grepl("density",colnames(zoops_summary_3groups))))]

schind_biom <- zoops_summary[ ,c(1:4,19, which(grepl("Biomass",colnames(zoops_summary))))]
schind_biom_3groups <- zoops_summary_3groups[ ,c(1:4,11, which(grepl("Biomass",colnames(zoops_summary_3groups))))]

write.csv(schind_dens, "Summer2021-DataAnalysis/SummaryStats/FCR_MOM_schind_2020-2022_zoopdens.csv", row.names = FALSE)
write.csv(schind_biom, "Summer2021-DataAnalysis/SummaryStats/FCR_MOM_schind_2020-2022_zoopbiom.csv", row.names = FALSE)

write.csv(schind_dens_3groups, "Summer2021-DataAnalysis/SummaryStats/FCR_MOM_schind_2020-2022_zoopdens_3groups.csv", row.names = FALSE)
write.csv(schind_biom_3groups, "Summer2021-DataAnalysis/SummaryStats/FCR_MOM_schind_2020-2022_zoopbiom_3groups.csv", row.names = FALSE)


#look at schindlers for all taxa

schind_dens_long <- schind_dens |> pivot_longer(Bosmina_density_NopL_rep.mean:nauplius_density_NopL_rep.mean)
  
ggplot(data=subset(schind_dens_long, !collect_date %in% c("2020-09-11","2020-09-15")),
       aes(x=value, y=as.numeric(depth),color=collect_date)) + 
  geom_point() + scale_y_reverse() + geom_path() + xlab("Density (#/L)") + 
  ylab("Depth (m)")  + facet_wrap(~name, scales="free") +
  theme(text = element_text(size=8), 
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5,size=6), 
        axis.text.y = element_text(size=6)) 
#ggsave("Summer2021-DataAnalysis/Figures/fcr_mom/Schindler_density_vs_depth_2021-2022_all_taxa.jpg")

