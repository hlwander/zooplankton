#read in libraries
pacman::p_load(tidyverse)

#playing around with schindler data
inUrl1  <- "https://pasta-s.lternet.edu/package/data/eml/edi/1090/14/c7a04035b0a99adc489f5b6daec1cd52" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")

zoops <-read.csv(infile1) |> filter(CollectionMethod=="Schindler" &
                                      Reservoir %in% c("FCR","BVR")) |> 
  select(-c(Site,EndDepth_m,CollectionMethod))

#new df for 2016-2018 zoop data
zoops_2016_2018 <- zoops[as.Date(zoops$DateTime)<"2019-01-01",]
zoops_2019_2021 <- zoops[as.Date(zoops$DateTime)>="2019-01-01",]

#just look at cladocerans, copepods, and rotifers here
zoops_3groups_pre <- zoops_2016_2018 |> group_by(Reservoir, DateTime, StartDepth_m) |> 
  summarise(Cladocera_Density_IndPerL = sum(Density_IndPerL[Taxon %in% c("Bosmina","D. catawba",
                                                      "Chydorus","D. ambigua",
                                                      "Diaphanosoma","Ceriodaphnia")]),
            Copepoda_Density_IndPerL = sum(Density_IndPerL[Taxon %in% c("Diaptomus","Nauplii",
                                                         "Cyclopoids")]),
            Rotifera_Density_IndPerL = sum(Density_IndPerL[Taxon %in% c("Total Rotifers")]),
            Cladocera_MeanLength_mm = mean(MeanLength_mm[Taxon %in% c("Bosmina","D. catawba",
                                                        "Chydorus","D. ambigua",
                                                        "Diaphanosoma","Ceriodaphnia")], 
                             na.rm=T),
            Copepoda_MeanLength_mm = mean(MeanLength_mm[Taxon %in% c("Diaptomus","Nauplii",
                                                          "Cyclopoids")],
                             na.rm=T),
            Rotifera_MeanLength_mm = mean(MeanLength_mm[Taxon %in% c("Total Rotifers")],
                            na.rm=T))

#wide to long
zoops_3groups_pre_dens_long <- zoops_3groups_pre |> 
  group_by(Reservoir, DateTime, StartDepth_m) |> 
  pivot_longer(cols=Cladocera_Density_IndPerL:Rotifera_Density_IndPerL,
               names_to = c("Taxon"),
               values_to = "Density_IndPerL")

zoops_3groups_pre_size_long <- zoops_3groups_pre |> 
  group_by(Reservoir, DateTime, StartDepth_m) |> 
  pivot_longer(cols=Cladocera_MeanLength_mm:Rotifera_MeanLength_mm,
               names_to = c("Taxon"),
               values_to = "MeanLength_mm")

#create long df with both dens and length
zoops_final_pre <- zoops_3groups_pre_size_long |> 
  select(-c("Cladocera_Density_IndPerL","Copepoda_Density_IndPerL",
            "Rotifera_Density_IndPerL")) |> 
  mutate(Taxon = str_extract(Taxon, "[^_]+")) 

#kinda hacky way to add density data
zoops_final_pre$Density_IndPerL <- zoops_3groups_pre_dens_long$Density_IndPerL

#drop midnight samples
zoops_final_pre <- zoops_final_pre |> 
  mutate(DateTime = as.POSIXct(DateTime, format="%Y-%m-%d %H:%M:%S", tz="UTC")) |> 
  filter(hour(DateTime) %in% c(9,10,11,12,13,14)) |>  #drop nighttime samples
  mutate(DateTime = as.Date(DateTime)) 

#average rep 1 and 2 when appropriate
zoops_final_post <- zoops_2019_2021 |> 
  mutate(DateTime = as.POSIXct(DateTime, format="%Y-%m-%d %H:%M:%S", tz="UTC")) |> 
  filter(hour(DateTime) %in% c(9,10,11,12,13,14)) |> #drop nighttime samples
  filter(Taxon %in% c("Cladocera","Copepoda","Rotifera")) |> 
  mutate(DateTime = as.Date(DateTime)) |> 
  group_by(Reservoir, StartDepth_m, DateTime, Taxon) |> 
  summarise(Density_IndPerL = mean(Density_IndPerL),
            MeanLength_mm = mean(MeanLength_mm))

#combine 2016-2018 + 2019-2021 data
all_zoops_final <- bind_rows(zoops_final_pre, zoops_final_post) |> 
  mutate_all(~replace(., is.nan(.), NA)) #replace NAN with NA

#plot schindler density for each day/reservoir
dates <-c("2016-06-08", "2016-08-03", "2021-06-15", "2021-07-07",
          "2020-09-11", "2020-09-15", "2021-06-10", "2022-07-01")

for(i in 1:length(dates)){
plot <- ggplot(data=subset(all_zoops_final, DateTime==dates[i]), 
       aes(Density_IndPerL, StartDepth_m)) + geom_point() +
  geom_path() + theme_bw() + ylim(rev(range(all_zoops_final$StartDepth_m))) +
  facet_wrap(~Taxon, scales="free_x") +
    ggtitle(paste0(all_zoops_final$Reservoir[all_zoops_final$DateTime==dates[i]],
                   " ",dates[i])) 
print(plot)
ggsave(paste0(getwd(),"/Figures/ch4/schind_zoops_",dates[i],".jpg"))

  }

#------------------------------------------------------------------------------#
#bring in ctd DO data

inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/200/13/27ceda6bc7fdec2e7d79a6e4fe16ffdf" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")

ctd_dates <- c("2016-06-08","2016-07-28","2021-06-16","2021-07-12", 
               "2020-09-11","2020-09-15","2021-06-10","2022-07-01")
#note that ctd casts were either taken on the day of zoop sampling or within a day
#EXCEPT for bvr #4 taken 5 days after zoop sampling (2021-07-07) and
#bvr #2 taken 6 days before zoop sampling (2016-08-03)

ctd<- read.csv(infile1) |>  
  select(Reservoir, DateTime, Site, Depth_m, DO_mgL) |> 
  mutate(DateTime = as.Date(DateTime)) |> 
  filter( Depth_m > 0 &
           DateTime %in% ctd_dates)           

#select 1m depth intervals and drop upcast
depths <- seq(0,10, by = 1)
ctd.final.raw<- ctd %>%
  group_by(DateTime, Reservoir) %>%
  slice(which.min(abs(as.numeric(Depth_m) - depths[1]))) #Create a new dataframe
ctd.final.raw$Depth_m <- depths[1]
#loop through all depths and add the closest values to the final dataframe
for (i in 2:length(depths)){
  ctd_atThisDepth <- ctd %>%
    group_by(DateTime, Reservoir) %>%
    slice(which.min(abs(as.numeric(Depth_m) - depths[i])))
  ctd_atThisDepth <- ctd_atThisDepth %>%
    #only include if the measured depth is within 0.1 of the depth label
    filter(abs(Depth_m-depths[i])<0.1)
  ctd_atThisDepth$Depth_m <- depths[i]
  ctd.final.raw <- rbind(ctd.final.raw,ctd_atThisDepth)
}
        
#reservoir for DO plots
res_temp <- c("BVR","BVR","BVR","BVR", "FCR","FCR","FCR","FCR")

#plot ctd and/or ysi data when available
for(i in 1:length(dates)){
plot <- ggplot(data = subset(ctd.final.raw, DateTime==ctd_dates[i] &
                               Reservoir==res_temp[i]), 
         aes(DO_mgL, Depth_m, col="DO (mg/L)")) + 
         geom_point() + geom_path() + theme_bw() + 
         ylim(rev(range(ctd.final.raw$Depth_m))) +
        geom_vline(xintercept=2, lty=2)+
        ggtitle(paste0(res_temp[i]," ",ctd_dates[i])) +
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
#ggsave(paste0(getwd(),"/Figures/ch4/CTD_DO_",ctd_dates[i],".jpg"))
  }




inUrl2  <-  "https://pasta.lternet.edu/package/data/eml/edi/198/11/6e5a0344231de7fcebbe6dc2bed0a1c3" 
infile2 <- tempfile()
try(download.file(inUrl2,infile2,method="curl"))
if (is.na(file.size(infile2))) download.file(inUrl2,infile2,method="auto")

ysi<- read.csv(infile2) |>  
  select(Reservoir, DateTime, Site, Depth_m, DO_mgL) |> 
  mutate(DateTime = as.Date(DateTime)) |> 
  filter(DateTime %in% c("2016-08-04", "2021-07-07") & #"2016-07-28"
           Reservoir == "BVR" )

#change site numbers for a few that were labeled wrong
ysi$Site[c(13:25)] <- 50

ggplot(data=subset(ysi, Site==50), aes(DO_mgL, Depth_m, col="DO (mg/L)")) + 
  geom_point() + geom_path() + theme_bw() + 
  ylim(rev(range(ysi$Depth_m))) +
  geom_vline(xintercept=2, lty=2)+
  facet_wrap(~DateTime) +
  scale_color_manual(values = c("black", "red")) +
  theme(text = element_text(size=10), axis.text = element_text(size=8, color="black"), 
        legend.position = c(0.72,0.94),
        legend.background = element_blank(),legend.direction = "horizontal", 
        panel.grid.minor = element_blank(), legend.key=element_rect(fill=NA),
        plot.margin = unit(c(0,0.05,0,0), "cm"),legend.key.size = unit(0.5, "lines"), 
        panel.grid.major = element_blank(), legend.title = element_blank(),
        legend.text  = element_text(size = 6), panel.spacing=unit(0, "cm"),
        axis.text.x = element_text(vjust = 0.5,size=8), axis.text.y = element_text(size=6)) 


#compare ysi and ctd to see if 2016-08-03 could have actually been neg heterograde
ctd_ysi_comp <- "2016-07-28"

ggplot(data = subset(ysi, DateTime==ctd_ysi_comp), 
       aes(DO_mgL, Depth_m, col="YSI DO (mg/L)")) + 
  geom_point() + geom_path() + theme_bw() + 
  ylim(rev(range(ysi$Depth_m))) +
  geom_vline(xintercept=2, lty=2)+
  geom_point(data = subset(ctd.final.raw, DateTime==ctd_ysi_comp & Reservoir=="BVR"), 
             aes(DO_mgL, Depth_m, col="CTD DO (mg/L))")) +
  geom_path(data = subset(ctd.final.raw, DateTime==ctd_ysi_comp& Reservoir=="BVR"), 
            aes(DO_mgL, Depth_m, col="CTD DO (mg/L))")) +
  ggtitle(paste0("BVR ", ctd_ysi_comp)) +
  scale_color_manual(values = c("black", "red")) +
  theme(text = element_text(size=10), axis.text = element_text(size=8, color="black"), 
        legend.position = c(0.72,0.94),
        legend.background = element_blank(),legend.direction = "horizontal", 
        panel.grid.minor = element_blank(), legend.key=element_rect(fill=NA),
        plot.margin = unit(c(0,0.05,0,0), "cm"),legend.key.size = unit(0.5, "lines"), 
        panel.grid.major = element_blank(), legend.title = element_blank(),
        legend.text  = element_text(size = 6), panel.spacing=unit(0, "cm"),
        axis.text.x = element_text(vjust = 0.5,size=8), axis.text.y = element_text(size=6)) 
  
#------------------------------------------------------------------------------#
#size and density histogram for fcr vs. bvr

zoops_tow <-read.csv(infile1) |> filter(CollectionMethod=="Tow" &
                                      Reservoir %in% c("FCR","BVR")) |> 
  select(-c(Site,EndDepth_m,CollectionMethod))

ggplot(data=subset(zoops_tow, Taxon %in% c("Cladocera","Copepoda","Rotifera")), 
       aes(MeanLength_mm, fill=Taxon))+ 
  geom_density(alpha=0.2) + theme_bw() +
  facet_wrap(~Reservoir)
#note that this is prob just 2019-2021 bc JPD used different taxa groupings     

ggplot(data=subset(zoops_tow, Taxon %in% c("Cladocera","Copepoda","Rotifera")), 
       aes(Density_IndPerL, fill=Taxon)) + xlim(0,100) + 
  geom_density(alpha=0.2) + theme_bw() +
  facet_wrap(~Reservoir)
