#Script to calculate DVM and DHM metrics
#Created 5Dec2022

#read in libraries
pacman::p_load(dplyr, vegan, labdsv, goeveg, rLakeAnalyzer, ggplot2,tidyr,lubridate, scales, colorblindcheck)

#Calculates the standard error####
stderr <- function(x) {
  sd(x,na.rm=TRUE)/sqrt(length(na.omit(x)))
}

#read in DHM csv
all_DHM <- read.csv(paste0(getwd(),"/Summer2021-DataAnalysis/SummaryStats/All_MSN_DHM.csv"))

#convert DHM hour to character
all_DHM$Hour <- ifelse(grepl("noon",all_DHM$sample_ID),"noon",
                       ifelse(grepl("midnight",all_DHM$sample_ID),"midnight",
                              ifelse(grepl("sunrise",all_DHM$sample_ID),"sunrise","sunset")))

#add taxa column
all_DHM$Taxa <- substr(all_DHM$metric,1,9)

#rename collect_date to DateTime
names(all_DHM)[names(all_DHM)=="collect_date"] <- "DateTime"

#drop sample_id
all_DHM <- all_DHM %>% select(!sample_ID)

#work with raw density and biomass to calculate metrics
all_DHM <- all_DHM[!grepl("percent",all_DHM$metric, ignore.case = TRUE),]

#add MSN# column
all_DHM$MSN <- ifelse(all_DHM$DateTime=="2019-07-10" | all_DHM$DateTime=="2019-07-11",1,
                      ifelse(all_DHM$DateTime=="2019-07-24" | all_DHM$DateTime=="2019-07-25",2,
                             ifelse(all_DHM$DateTime=="2020-08-12" | all_DHM$DateTime=="2020-08-13",3,
                                    ifelse(all_DHM$DateTime=="2021-06-15" | all_DHM$DateTime=="2021-06-16",4,5))))

#change second noon for each MSN to noon2
all_DHM$Hour <- ifelse(all_DHM$Hour=="noon" & (all_DHM$DateTime=="2019-07-11" | 
                                                 all_DHM$DateTime=="2019-07-25" |
                                                 all_DHM$DateTime=="2020-08-13" |
                                                 all_DHM$DateTime=="2021-06-16" |
                                                 all_DHM$DateTime=="2021-07-08"), "noon2", all_DHM$Hour)

#read in DVM annual csvs
DVM_2021 <- read.csv(paste0(getwd(),"/Summer2021-DataAnalysis/SummaryStats/DVM_2021_zoops.csv")) %>% select(-X)
DVM_2020 <- read.csv(paste0(getwd(),"/Summer2020-DataAnalysis/SummaryStats/DVM_2020_zoops.csv")) %>% select(-X)
DVM_2019 <- read.csv(paste0(getwd(),"/Summer2019-DataAnalysis/SummaryStats/DVM_2019_zoops.csv")) 

#combine annual DVM dfs
all_DVM <- rbind(DVM_2019,DVM_2020,DVM_2021)

#change date format
all_DVM$DateTime <- format(as.Date(all_DVM$DateTime, "%m-%d-%Y"), "%Y-%m-%d")

#work with raw density and biomass to calculate metrics
all_DVM <- all_DVM[!grepl("percent",all_DVM$metric),]

#add MSN# column
all_DVM$MSN <- ifelse(all_DVM$DateTime=="2019-07-10" | all_DVM$DateTime=="2019-07-11",1,
                            ifelse(all_DVM$DateTime=="2019-07-24" | all_DVM$DateTime=="2019-07-25",2,
                                   ifelse(all_DVM$DateTime=="2020-08-12" | all_DVM$DateTime=="2020-08-13",3,
                                          ifelse(all_DVM$DateTime=="2021-06-15" | all_DVM$DateTime=="2021-06-16",4,5))))

#change second noon for each MSN to noon2
all_DVM$Hour <- ifelse(all_DVM$Hour=="noon" & (all_DVM$DateTime=="2019-07-11" | 
                                                 all_DVM$DateTime=="2019-07-25" |
                                                 all_DVM$DateTime=="2020-08-13" |
                                                 all_DVM$DateTime=="2021-06-16" |
                                                 all_DVM$DateTime=="2021-07-08"), "noon2", all_DVM$Hour)

#not sure if this is the right call, but going to set all negative values to 0 bc they are not real and due to additive error across all steps
all_DVM$value[all_DVM$value < 0] <- 0 #n=20

#shorten metric name
all_DVM$metric <- ifelse(all_DVM$WaterColumn=="epilimnion", substr(all_DVM$metric,1,nchar(all_DVM$metric)-13),substr(all_DVM$metric,1,nchar(all_DVM$metric)-14))

#--------------------------------------------------------------------------------------------------#
#calculate DVM and DHM metrics
#note that I think we are missing  noon tow from 24jul19 --> will look into this

#Calculate proportion of zoops in epi at noon and midnight
DVM_proportion <-  plyr::ddply(all_DVM, c("metric", "MSN", "Hour","DateTime"), function(x) {
  data.frame(
    proportion_epi = x$value[x$WaterColumn=="epilimnion"] / sum(x$value)
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 

#initialize df for DVM metrics for each day
#migration_df <- data.frame("MSN"= unique(DVM_proportion$MSN))

#DVM metrics for a single day --> DVM = (Depi / Depi + Dhypo)Night - (Depi / Depi + Dhypo)Day
DVM_noon1 <-  plyr::ddply(DVM_proportion, c("metric", "MSN"), function(x) {
  data.frame(
    DVM1 = x$proportion_epi[x$Hour=="midnight"] - x$proportion_epi[x$Hour=="noon"]
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 


DVM_noon2 <-  plyr::ddply(DVM_proportion, c("metric", "MSN"), function(x) {
  data.frame(
DVM2 = x$proportion_epi[x$Hour=="midnight"] - x$proportion_epi[x$Hour=="noon2"]
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 

#merge both dfs 
DVM_metrics <- right_join(DVM_noon1,DVM_noon2, by=c("metric","MSN"))


#Calculate proportion of zoops in epi at noon and midnight
DHM_proportion <-  plyr::ddply(all_DHM, c("metric", "MSN", "Hour","DateTime"), function(x) {
  data.frame(
    proportion_epi = x$value[x$site_no=="BVR_50"] / sum(x$value)
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 


#DHM metrics for a single day --> DHM = (Dpelepi / Dpelepi + Dlit)Night - (Dpelepi / Dpelepi + Dlit)Day
DHM_metrics <-  plyr::ddply(DHM_proportion, c("metric", "MSN"), function(x) {
  data.frame(
    DHM1 = x$proportion_epi[x$Hour=="midnight"] - x$proportion_epi[x$Hour=="noon"],
    DHM2 = x$proportion_epi[x$Hour=="midnight"] - x$proportion_epi[x$Hour=="noon2"],
    #DHM3 = x$proportion_epi[x$Hour=="sunset"] - x$proportion_epi[x$Hour=="sunrise"]
    )
}, .progress = plyr::progress_text(), .parallel = FALSE) 

#initialize final migration df
migration_df <- data.frame(metric = DHM_metrics$metric, MSN = DHM_metrics$MSN)

#get DVM df into same order as DHM df
DVM_metrics <- DVM_metrics[order(DVM_metrics$metric, DVM_metrics$MSN),]

#add average and SE of DVM and DHM calcs
migration_df$DVM_avg <- rowMeans(DVM_metrics[,c(3,4)], na.rm=TRUE)
migration_df$DHM_avg <- rowMeans(DHM_metrics[,c(3,4)], na.rm=TRUE)

#calculate stderr across for migration metrics
migration_df$DVM_SE <- apply(DVM_metrics[,c(3,4)], 1, stderr)
migration_df$DHM_SE <- apply(DHM_metrics[,c(3,4)], 1, stderr)

#replace NAN with 0 bc none of that taxa were found
migration_df$DHM_avg[is.nan(migration_df$DHM_avg)] <- 0

#now convert from wide to long
df1 <- migration_df %>% gather(migration, value, DVM_avg:DHM_avg)
df2 <- migration_df %>% gather(migration, value, DVM_SE:DHM_SE)

migration_long <- df1[,c(1,2,5,6)]
migration_long$SE <- df2$value

#export migration metrics
write.csv(migration_long,"./Summer2021-DataAnalysis/SummaryStats/migration_metrics.csv",row.names = FALSE)

#-------------------------------------------------------------------------------#

metric_taxa <-c("Total","Calanoida","Calanoida","Cladocera","Cladocera",
                "Copepoda","Copepoda","Cyclopoida", "Cyclopoida",
                 "Rotifera", "Rotifera", "Total")
names(metric_taxa) <- c(unique(migration_long$metric))

#plot migration metrics
ggplot(subset(migration_long, grepl("density",metric)), 
              aes(x=MSN, y=value, color=metric, shape=migration)) + 
  geom_point() + theme_bw() + geom_hline(yintercept = 0, linetype="dotted")+
  theme(text = element_text(size=8), axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), legend.key = element_blank(), 
        legend.key.height=unit(0.3,"line"), 
        axis.text.x = element_text(vjust = 0.5, hjust=1), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = c(0.92,0.94), legend.spacing = unit(-0.5, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.7,"line")) + guides(color=FALSE) +
  scale_color_manual("",values=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"))+
  facet_wrap(~metric, labeller = labeller(metric=metric_taxa)) + ylab("Density migration metric")
ggsave(file.path(getwd(),"Summer2021-DataAnalysis/Figures/BVR_MSNs_migration_metrics.jpg"), width=5, height=4) 

