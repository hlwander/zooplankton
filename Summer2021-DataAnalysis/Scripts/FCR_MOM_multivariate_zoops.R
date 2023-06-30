#script for multivariate stats using FCR MOM + environmental correlate data
#created 26 May 2023

#read in libraries
pacman::p_load(tidyverse, labdsv, goeveg, vegan, ggordiplots, viridis, RColorBrewer,
               rLakeAnalyzer, ggrepel)

#read in FCR MOM schindler data
zoop_dens <- read.csv("Summer2021-DataAnalysis/SummaryStats/FCR_MOM_schind_2020-2022_zoopdens.csv", header=TRUE)
zoop_biom <- read.csv("Summer2021-DataAnalysis/SummaryStats/FCR_MOM_schind_2020-2022_zoopbiom.csv", header=TRUE)

#add anoxic vs mom col
zoop_dens$DO <- ifelse(zoop_dens$collect_date %in% c("2020-09-11", "2020-09-15"), "MOM", "anoxic")
zoop_biom$DO <- ifelse(zoop_biom$collect_date %in% c("2020-09-11", "2020-09-15"), "MOM", "anoxic")

#add msn #
zoop_dens$msn <- ifelse(zoop_dens$collect_date=="2020-09-11", "1",
                        ifelse(zoop_dens$collect_date=="2020-09-15", "2",
                               ifelse(zoop_dens$collect_date=="2021-06-10" | 
                                        zoop_dens$collect_date=="2021-06-11" , "3", "4")))

zoop_biom$msn <- ifelse(zoop_biom$collect_date=="2020-09-11", "1",
                        ifelse(zoop_biom$collect_date=="2020-09-15", "2",
                               ifelse(zoop_biom$collect_date=="2021-06-10" | 
                                        zoop_biom$collect_date=="2021-06-11" , "3", "4")))

#split data up to look and noons alone and then day/night pairs
#zoop_dens_noons <- zoop_dens

#only select data cols for NMDS
density_only <- zoop_dens[, c(grepl("density", colnames(zoop_dens)))]
biomass_only <- zoop_biom[, c(grepl("Biomass", colnames(zoop_biom)))]

#hellinger transformation to convert abundances from absolute to relative (and give low weight to low/0 values)
zoop_dens_trans <- hellinger(density_only)
zoop_biom_trans <- hellinger(biomass_only)

#scree plot to choose # of dimensions --> using transformed data instead of distance matrix for now
dimcheckMDS(zoop_dens_trans, distance = "bray", k = 6, trymax = 20, autotransform = TRUE)

#now do NMDS w/ 4 dimensions
set.seed(1)
NMDS_dens_bray <- metaMDS(zoop_dens_trans, distance='bray', k=4, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_dens_bray$stress

date_list <- c("28-29Jun2020", "10-11Sep2020", "14-15Sep2020", "10-11Jun2021", "30Jun-01Jul2022")

#plot NMDS
ord <- ordiplot(NMDS_dens_bray,display = c('sites','species'),choices = c(1,2),type = "n")
depth <- gg_ordiplot(ord, zoop_dens$depth, kind = "ehull", 
                     ellipse=FALSE, hull = TRUE, plot = FALSE, pt.size=0.9) 

NMDS_depth <- depth$plot + geom_point() + theme_bw() + 
  geom_polygon(data = depth$df_hull, aes(x = x, y = y, fill = Group), alpha=0.2) +
  geom_point(data=depth$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2,
             fill = brewer.pal(10, "BrBG")) +
  theme(text = element_text(size=7), axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), 
        legend.key.height=unit(0.5,"line"),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(0,-0,-0,-0),
        legend.direction = "horizontal",
        legend.text = element_text(size=9),
        axis.text.x = element_text(vjust = 0.5), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "top",
        plot.margin = unit(c(1,-0.1,0,0), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.5,"line")) +
  guides(fill="none", color = guide_legend(override.aes = list(size = 2))) +
     #   annotate("text", x=-0.42, y=0.6, label= "bolditalic(a)", parse=T, size = 3) +
        scale_fill_manual("",values=brewer.pal(10, "BrBG"))+
        scale_color_manual("",values=brewer.pal(10, "BrBG"),
                     label=c(unique(zoop_dens$depth))) 
#ggsave("Summer2021-DataAnalysis/Figures/fcr_mom/NMDS_dens_depths_fcr_2v1.jpg")

time <- gg_ordiplot(ord, zoop_dens$Hour, kind = "ehull", 
                     ellipse=FALSE, hull = TRUE, plot = FALSE, pt.size=0.9) 

NMDS_time <- time$plot + geom_point() + theme_bw() + 
  geom_polygon(data = time$df_hull, aes(x = x, y = y, fill = Group), alpha=0.2) +
  geom_point(data=time$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2,
             fill = c("#2b4261","#45d3f7")) +
  theme(text = element_text(size=7), axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), 
        legend.key.height=unit(0.3,"line"),
        legend.key = element_blank(),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-0,-0,-0,-0),
        legend.direction = "horizontal",
        legend.text = element_text(size=9),
        axis.text.x = element_text(vjust = 0.5), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "top", legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(1,-0.1,0,0), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) + 
  guides(fill="none", color = guide_legend(override.aes = list(size = 2))) +
  #annotate("text", x=-0.42, y=0.6, label= "bolditalic(a)", parse=T, size = 3) +
  scale_fill_manual("",values=c("#2b4261","#45d3f7"))+
  scale_color_manual("",values=c("#2b4261","#45d3f7"),
                     label=c('Midnight','Noon')) 
#ggsave("Summer2021-DataAnalysis/Figures/fcr_mom/NMDS_dens_time_fcr_2v1.jpg")

day <- gg_ordiplot(ord, zoop_dens$collect_date, kind = "ehull", 
                     ellipse=FALSE, hull = TRUE, plot = FALSE, pt.size=0.9) 
                     
NMDS_day <- day$plot + geom_point() + theme_bw() + 
  geom_polygon(data = day$df_hull, aes(x = x, y = y, fill = Group), alpha=0.2) +
  geom_point(data=day$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2,
             fill = c("#111258","#4C1559","#B1235E","#F0855F","#F8C48A","#F8F4A9")) +
  theme(text = element_text(size=7), axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), 
        legend.key.height=unit(0.3,"line"),
        legend.key = element_blank(),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-0,-0,-0,-0),
        legend.direction = "horizontal",
        legend.text = element_text(size=9),
        axis.text.x = element_text(vjust = 0.5), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "top", legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(1,-0.1,0,0), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) + 
  guides(fill="none", color = guide_legend(override.aes = list(size = 2))) +
  #annotate("text", x=-0.42, y=0.6, label= "bolditalic(a)", parse=T, size = 3) +
  scale_fill_manual("",values=c("#111258","#4C1559","#B1235E","#F0855F","#F8C48A","#F8F4A9"))+
  scale_color_manual("",values=c("#111258","#4C1559","#B1235E","#F0855F","#F8C48A","#F8F4A9"),
                     label=c("2020-09-11", "2020-09-15", "2021-06-10", "2021-06-11",
                             "2022-06-30", "2022-07-01")) 
#ggsave("Summer2021-DataAnalysis/Figures/fcr_mom/NMDS_dens_day_fcr_2v1.jpg")

do <- gg_ordiplot(ord, zoop_dens$DO, kind = "ehull", 
                   ellipse=FALSE, hull = TRUE, plot = FALSE, pt.size=0.9) 

NMDS_do <- do$plot + geom_point() + theme_bw() + 
  geom_polygon(data = do$df_hull, aes(x = x, y = y, fill = Group), alpha=0.2) +
  geom_point(data=do$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2,
             fill = c("#d33131","#f1d8a7")) +
  theme(text = element_text(size=7), axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), 
        legend.key.height=unit(0.3,"line"),
        legend.key = element_blank(),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-0,-0,-0,-0),
        legend.direction = "horizontal",
        legend.text = element_text(size=9),
        axis.text.x = element_text(vjust = 0.5), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "top", legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(1,-0.1,0,0), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) + 
  guides(fill="none", color = guide_legend(override.aes = list(size = 2))) +
  #annotate("text", x=-0.42, y=0.6, label= "bolditalic(a)", parse=T, size = 3) +
  scale_fill_manual("",values=c("#d33131","#f1d8a7"))+
  scale_color_manual("",values=c("#d33131","#f1d8a7"),
                     label=c("anoxic","MOM")) 
#ggsave("Summer2021-DataAnalysis/Figures/fcr_mom/NMDS_dens_do_fcr_2v1.jpg")

#------------------------------------------------------------------------------#
#same for biomass

#now do NMDS w/ 4 dimensions
set.seed(1)
NMDS_biom_bray <- metaMDS(zoop_biom_trans, distance='bray', k=4, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_biom_bray$stress

#plot NMDS
ord <- ordiplot(NMDS_biom_bray,display = c('sites','species'),choices = c(1,2),type = "n")
depth <- gg_ordiplot(ord, zoop_biom$depth, kind = "ehull", 
                     ellipse=FALSE, hull = TRUE, plot = FALSE, pt.size=0.9) 

NMDS_depth <- depth$plot + geom_point() + theme_bw() + 
  geom_polygon(data = depth$df_hull, aes(x = x, y = y, fill = Group), alpha=0.2) +
  geom_point(data=depth$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2,
             fill = brewer.pal(10, "BrBG")) +
  theme(text = element_text(size=7), axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), 
        legend.key.height=unit(0.5,"line"),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(0,-0,-0,-0),
        legend.direction = "horizontal",
        legend.text = element_text(size=9),
        axis.text.x = element_text(vjust = 0.5), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "top",
        plot.margin = unit(c(1,-0,0,0), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.5,"line")) +
  guides(fill="none", color = guide_legend(override.aes = list(size = 2))) +
  #   annotate("text", x=-0.42, y=0.6, label= "bolditalic(a)", parse=T, size = 3) +
  scale_fill_manual("",values=brewer.pal(10, "BrBG"))+
  scale_color_manual("",values=brewer.pal(10, "BrBG"),
                     label=c(unique(zoop_biom$depth))) 
#ggsave("Summer2021-DataAnalysis/Figures/fcr_mom/NMDS_biom_depths_fcr_2v1.jpg")

time <- gg_ordiplot(ord, zoop_biom$Hour, kind = "ehull", 
                    ellipse=FALSE, hull = TRUE, plot = FALSE, pt.size=0.9) 

NMDS_time <- time$plot + geom_point() + theme_bw() + 
  geom_polygon(data = time$df_hull, aes(x = x, y = y, fill = Group), alpha=0.2) +
  geom_point(data=time$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2,
             fill = c("#2b4261","#45d3f7")) +
  theme(text = element_text(size=7), axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), 
        legend.key.height=unit(0.3,"line"),
        legend.key = element_blank(),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-0,-0,-0,-0),
        legend.direction = "horizontal",
        legend.text = element_text(size=9),
        axis.text.x = element_text(vjust = 0.5), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "top", legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(1,0,0,0), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) + 
  guides(fill="none", color = guide_legend(override.aes = list(size = 2))) +
  #annotate("text", x=-0.42, y=0.6, label= "bolditalic(a)", parse=T, size = 3) +
  scale_fill_manual("",values=c("#2b4261","#45d3f7"))+
  scale_color_manual("",values=c("#2b4261","#45d3f7"),
                     label=c('Midnight','Noon')) 
#ggsave("Summer2021-DataAnalysis/Figures/fcr_mom/NMDS_biom_time_fcr_2v1.jpg")

day <- gg_ordiplot(ord, zoop_biom$collect_date, kind = "ehull", 
                   ellipse=FALSE, hull = TRUE, plot = FALSE, pt.size=0.9) 

NMDS_day <- day$plot + geom_point() + theme_bw() + 
  geom_polygon(data = day$df_hull, aes(x = x, y = y, fill = Group), alpha=0.2) +
  geom_point(data=day$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2,
             fill = c("#111258","#4C1559","#B1235E","#F0855F","#F8C48A","#F8F4A9")) +
  theme(text = element_text(size=7), axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), 
        legend.key.height=unit(0.3,"line"),
        legend.key = element_blank(),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-0,-0,-0,-0),
        legend.direction = "horizontal",
        legend.text = element_text(size=9),
        axis.text.x = element_text(vjust = 0.5), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "top", legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(1,0,0,0), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) + 
  guides(fill="none", color = guide_legend(override.aes = list(size = 2))) +
  #annotate("text", x=-0.42, y=0.6, label= "bolditalic(a)", parse=T, size = 3) +
  scale_fill_manual("",values=c("#111258","#4C1559","#B1235E","#F0855F","#F8C48A","#F8F4A9"))+
  scale_color_manual("",values=c("#111258","#4C1559","#B1235E","#F0855F","#F8C48A","#F8F4A9"),
                     label=c("2020-09-11", "2020-09-15", "2021-06-10", "2021-06-11",
                             "2022-06-30", "2022-07-01"))  
#ggsave("Summer2021-DataAnalysis/Figures/fcr_mom/NMDS_biom_day_fcr_2v1.jpg")

do <- gg_ordiplot(ord, zoop_biom$DO, kind = "ehull", 
                  ellipse=FALSE, hull = TRUE, plot = FALSE, pt.size=0.9) 

NMDS_do <- do$plot + geom_point() + theme_bw() + 
  geom_polygon(data = do$df_hull, aes(x = x, y = y, fill = Group), alpha=0.2) +
  geom_point(data=do$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2,
             fill = c("#d33131","#f1d8a7")) +
  theme(text = element_text(size=7), axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), 
        legend.key.height=unit(0.3,"line"),
        legend.key = element_blank(),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-0,-0,-0,-0),
        legend.direction = "horizontal",
        legend.text = element_text(size=9),
        axis.text.x = element_text(vjust = 0.5), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "top", legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(1,0,0,0), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) + 
  guides(fill="none", color = guide_legend(override.aes = list(size = 2))) +
  #annotate("text", x=-0.42, y=0.6, label= "bolditalic(a)", parse=T, size = 3) +
  scale_fill_manual("",values=c("#d33131","#f1d8a7"))+
  scale_color_manual("",values=c("#d33131","#f1d8a7"),
                     label=c("anoxic","MOM")) 
#ggsave("Summer2021-DataAnalysis/Figures/fcr_mom/NMDS_biom_do_fcr_2v1.jpg")


#------------------------------------------------------------------------------#
#pull in env data for driver analysis

ctd<- read.csv('~/Documents/VirginiaTech/research/BVR_GLM/bvr_glm/field_data/CTD_final_2013_2022.csv') %>% 
  mutate(DateTime = as.Date(DateTime)) %>%
  filter(Reservoir =="FCR" & Depth_m > 0 &
           DateTime %in% c(as.Date("2020-09-10"), as.Date("2020-09-11"), as.Date("2020-09-14"),
                           as.Date("2020-09-15"), as.Date("2021-06-10"), as.Date("2021-06-11"),
                           as.Date("2022-06-30"), as.Date("2022-07-01")))

#select every 0.5m from casts
ctd_final_temp <- ctd %>%
  dplyr::mutate(rdepth = plyr::round_any(Depth_m, 0.5)) %>% 
  dplyr::group_by(DateTime, rdepth, Reservoir, Site) %>%
  dplyr::summarise(value = mean(Temp_C)) %>% 
  dplyr::rename(depth = rdepth) 

#do the same for DO values now
ctd_final_DO <- ctd %>%
  dplyr::mutate(rdepth = plyr::round_any(Depth_m, 0.5)) %>% 
  dplyr::group_by(DateTime, rdepth, Reservoir, Site) %>%
  dplyr::summarise(value = mean(DO_mgL)) %>% 
  dplyr::rename(depth = rdepth) 


#calculate thermocline depth and oxycline depth by date
ctd_thermo_depth <- ctd_final_temp %>% group_by(DateTime) %>% filter(depth > 0) %>% 
  mutate(therm_depth = thermo.depth(value,depth))

#add msn #
ctd_thermo_depth$msn  <- ifelse(ctd_thermo_depth$DateTime=="2020-09-11", "1",
                            ifelse(ctd_thermo_depth$DateTime=="2020-09-15", "2",
                                   ifelse(ctd_thermo_depth$DateTime=="2021-06-10" |
                                            ctd_thermo_depth$DateTime=="2021-06-11" , "3", "4")))

ctd_oxy_depth <- ctd_final_DO %>% group_by(DateTime) %>% filter(depth > 0) %>% 
  mutate(oxy_depth = thermo.depth(value,depth))

#add msn #
ctd_oxy_depth$msn  <- ifelse(ctd_oxy_depth$DateTime=="2020-09-11", "1",
                         ifelse(ctd_oxy_depth$DateTime=="2020-09-15", "2",
                                ifelse(ctd_oxy_depth$DateTime=="2021-06-10" | 
                                         ctd_oxy_depth$DateTime=="2021-06-11" , "3", "4")))

depths = c(0.1, seq(1, 9, by = 1))
ctd.final<-data.frame() 

for (i in 1:length(depths)){
  ctd_layer <- ctd %>% group_by(DateTime) %>% slice(which.min(abs(as.numeric(Depth_m) - depths[i])))
  # Bind each of the data layers together.
  ctd.final = bind_rows(ctd.final, ctd_layer)
}

#round to one decimal place
ctd.final$Depth_m <- round(ctd.final$Depth_m,1)

#add msn # 
ctd.final$msn <-  ifelse(ctd.final$DateTime=="2020-09-11", "1",
                         ifelse(ctd.final$DateTime=="2020-09-15", "2",
                                ifelse(ctd.final$DateTime=="2021-06-10" | 
                                         ctd.final$DateTime=="2021-06-11", "3", "4")))

#average across msn and depth
ctd.final_avg <- ctd.final %>% group_by(Reservoir, Depth_m, msn) %>% summarise(across(everything(), list(mean), na.rm=T))


chem <- read.csv('~/Documents/VirginiaTech/research/Reservoirs/Data/DataAlreadyUploadedToEDI/EDIProductionFiles/MakeEMLChemistry/2022/chemistry_2013_2022.csv') %>% 
  mutate(DateTime = as.Date(DateTime)) %>%
  filter(Reservoir =="FCR" & Site==50 & 
           DateTime %in% c(as.Date("2020-09-10"), as.Date("2020-09-11"), as.Date("2020-09-14"),
                           as.Date("2020-09-15"), as.Date("2021-06-10"), as.Date("2021-06-11"),
                           as.Date("2022-06-30"), as.Date("2022-07-01"), as.Date("2022-06-27"), 
                           as.Date("2021-06-07"))) # 2021 and 2022 chem samples were taken 3 days earlier than MSN

#add msn #
chem$msn <- ifelse(chem$DateTime=="2020-09-11", "1",
                   ifelse(chem$DateTime=="2020-09-15", "2",
                          ifelse(chem$DateTime=="2021-06-07", "3", "4")))

#average across msn and depth
chem_avg <- chem %>% group_by(Reservoir, Depth_m, msn) %>% summarise(across(everything(), list(mean)))



secchi <- read.csv('~/Documents/VirginiaTech/research/Reservoirs/Data/DataAlreadyUploadedToEDI/EDIProductionFiles/MakeEMLYSI_PAR_secchi/2022/Data/Secchi_depth_2013-2022.csv') %>% 
  mutate(DateTime = as.Date(DateTime)) %>%
  filter(Reservoir =="FCR" & Site==50 & 
           DateTime %in% c(as.Date("2020-09-10"), as.Date("2020-09-11"), as.Date("2020-09-14"),
                           as.Date("2020-09-15"), as.Date("2021-06-10"), as.Date("2021-06-11"),
                           as.Date("2022-06-30"), as.Date("2022-07-01"), as.Date("2022-06-27"), 
                           as.Date("2021-06-07")))


ysi <- read.csv('~/Documents/VirginiaTech/research/Reservoirs/Data/DataAlreadyUploadedToEDI/EDIProductionFiles/MakeEMLYSI_PAR_secchi/2022/Data/YSI_PAR_profiles_2013-2022.csv') %>% 
  mutate(DateTime = as.Date(DateTime)) %>%
  filter(Reservoir =="FCR" & Site==50 & 
           DateTime %in% c(as.Date("2020-09-10"), as.Date("2020-09-11"), as.Date("2020-09-14"),
                           as.Date("2020-09-15"), as.Date("2021-06-10"), as.Date("2021-06-11"),
                           as.Date("2022-06-30"), as.Date("2022-07-01"), as.Date("2022-06-27"))) # 2022 ysi profile collected 3 days earlier than MSN


fp <- read.csv('~/Documents/VirginiaTech/research/Reservoirs/Data/DataAlreadyUploadedToEDI/EDIProductionFiles/MakeEMLFluoroProbe/2022/FluoroProbe_2014_2022_temp.csv') %>% 
  mutate(DateTime = as.Date(DateTime)) %>%
  filter(Reservoir =="FCR" & Site==50 & Depth_m < 9.1 &
           DateTime %in% c(as.Date("2020-09-10"), as.Date("2020-09-11"),
                           as.Date("2020-09-15"), as.Date("2021-06-10"), as.Date("2021-06-11"),
                           as.Date("2022-06-30"), as.Date("2022-07-01"), as.Date("2021-06-07")))  # 2021 fp profile collected 3 days earlier than MSN

#select every 1m from casts
fp_final <- fp %>%
  dplyr::mutate(rdepth = plyr::round_any(Depth_m, 1)) %>% 
  dplyr::group_by(DateTime, rdepth, Reservoir, Site) %>%
  dplyr::summarise(across(c(everything()),mean, na.rm=T)) %>% 
  dplyr::rename(depth = rdepth) 

#only keep first 11 cols in fp df
fp_final <- fp_final[,c(1:11)]

#add msn col
fp_final$msn <- ifelse(fp_final$DateTime=="2020-09-11", "1",
                       ifelse(fp_final$DateTime=="2020-09-15", "2",
                              ifelse(fp_final$DateTime=="2021-06-07", "3", "4")))


#------------------------------------------------------------------------------
#env fit with daily drivers

daily_drivers  <- data.frame("msn" = as.character(1:4), #mean full water column vals
                          "temp" = c(mean(ctd.final_avg$Temp_C_1[ctd.final_avg$msn==1]),
                                     mean(ctd.final_avg$Temp_C_1[ctd.final_avg$msn==2]),
                                     mean(ctd.final_avg$Temp_C_1[ctd.final_avg$msn==3]),
                                     mean(ctd.final_avg$Temp_C_1[ctd.final_avg$msn==4])), 
                          "DO_conc" = c(mean(ctd.final_avg$DO_mgL_1[ctd.final_avg$msn==1]),
                                        mean(ctd.final_avg$DO_mgL_1[ctd.final_avg$msn==2]),
                                        mean(ctd.final_avg$DO_mgL_1[ctd.final_avg$msn==3]),
                                        mean(ctd.final_avg$DO_mgL_1[ctd.final_avg$msn==4])),
                          "DOsat" = c(mean(ctd.final_avg$DOsat_percent_1[ctd.final_avg$msn==1]),
                                      mean(ctd.final_avg$DOsat_percent_1[ctd.final_avg$msn==2]),
                                      mean(ctd.final_avg$DOsat_percent_1[ctd.final_avg$msn==3]),
                                      mean(ctd.final_avg$DOsat_percent_1[ctd.final_avg$msn==4])),
                          "Spcond" = c(mean(ctd.final_avg$SpCond_uScm_1[ctd.final_avg$msn==1]),
                                       mean(ctd.final_avg$SpCond_uScm_1[ctd.final_avg$msn==2]),
                                       mean(ctd.final_avg$SpCond_uScm_1[ctd.final_avg$msn==3]),
                                       mean(ctd.final_avg$SpCond_uScm_1[ctd.final_avg$msn==4])),
                        "Green_algae" = c(mean(fp_final$GreenAlgae_ugL[fp_final$msn==1]),
                                     mean(fp_final$GreenAlgae_ugL[fp_final$msn==2]),
                                     mean(fp_final$GreenAlgae_ugL[fp_final$msn==3]),
                                     mean(fp_final$GreenAlgae_ugL[fp_final$msn==4])),
                        "Bluegreen_algae" = c(mean(fp_final$Bluegreens_ugL[fp_final$msn==1]),
                                          mean(fp_final$Bluegreens_ugL[fp_final$msn==2]),
                                          mean(fp_final$Bluegreens_ugL[fp_final$msn==3]),
                                          mean(fp_final$Bluegreens_ugL[fp_final$msn==4])),
                        "Brown_algae" = c(mean(fp_final$BrownAlgae_ugL[fp_final$msn==1]),
                                              mean(fp_final$BrownAlgae_ugL[fp_final$msn==2]),
                                              mean(fp_final$BrownAlgae_ugL[fp_final$msn==3]),
                                              mean(fp_final$BrownAlgae_ugL[fp_final$msn==4])),
                        "Mixed_algae" = c(mean(fp_final$MixedAlgae_ugL[fp_final$msn==1]),
                                          mean(fp_final$MixedAlgae_ugL[fp_final$msn==2]),
                                          mean(fp_final$MixedAlgae_ugL[fp_final$msn==3]),
                                          mean(fp_final$MixedAlgae_ugL[fp_final$msn==4])),
                        "Yellow_substances" = c(mean(fp_final$YellowSubstances_ugL[fp_final$msn==1]),
                                          mean(fp_final$YellowSubstances_ugL[fp_final$msn==2]),
                                          mean(fp_final$YellowSubstances_ugL[fp_final$msn==3]),
                                          mean(fp_final$YellowSubstances_ugL[fp_final$msn==4])),
                        "Total_algae" = c(mean(fp_final$TotalConc_ugL[fp_final$msn==1]),
                                                mean(fp_final$TotalConc_ugL[fp_final$msn==2]),
                                                mean(fp_final$TotalConc_ugL[fp_final$msn==3]),
                                                mean(fp_final$TotalConc_ugL[fp_final$msn==4])),
                          "PAR" = c(mean(ctd.final_avg$PAR_umolm2s_1[ctd.final_avg$msn==1]),
                                        mean(ctd.final_avg$PAR_umolm2s_1[ctd.final_avg$msn==2]),
                                        mean(ctd.final_avg$PAR_umolm2s_1[ctd.final_avg$msn==3]),
                                        mean(ctd.final_avg$PAR_umolm2s_1[ctd.final_avg$msn==4])),
                          "TN" = c(mean(chem_avg$TN_ugL_1[chem_avg$msn==1]),
                                       mean(chem_avg$TN_ugL_1[chem_avg$msn==2]),
                                       mean(chem_avg$TN_ugL_1[chem_avg$msn==3]),
                                       mean(chem_avg$TN_ugL_1[chem_avg$msn==4])), 
                          "TP" = c(mean(chem_avg$TP_ugL_1[chem_avg$msn==1]),
                                   mean(chem_avg$TP_ugL_1[chem_avg$msn==2]),
                                   mean(chem_avg$TP_ugL_1[chem_avg$msn==3]),
                                   mean(chem_avg$TP_ugL_1[chem_avg$msn==4])), 
                        "NO3" = c(mean(chem_avg$NO3NO2_ugL_1[chem_avg$msn==1]),
                                 mean(chem_avg$NO3NO2_ugL_1[chem_avg$msn==2]),
                                 mean(chem_avg$NO3NO2_ugL_1[chem_avg$msn==3]),
                                 mean(chem_avg$NO3NO2_ugL_1[chem_avg$msn==4])), 
                        "NH4" = c(mean(chem_avg$NH4_ugL_1[chem_avg$msn==1]),
                                  mean(chem_avg$NH4_ugL_1[chem_avg$msn==2]),
                                  mean(chem_avg$NH4_ugL_1[chem_avg$msn==3]),
                                  mean(chem_avg$NH4_ugL_1[chem_avg$msn==4])), 
                        "PO4" = c(mean(chem_avg$SRP_ugL_1[chem_avg$msn==1]),
                                  mean(chem_avg$SRP_ugL_1[chem_avg$msn==2]),
                                  mean(chem_avg$SRP_ugL_1[chem_avg$msn==3]),
                                  mean(chem_avg$SRP_ugL_1[chem_avg$msn==4])), 
                        "DOC" = c(mean(chem_avg$DOC_mgL_1[chem_avg$msn==1]),
                                  mean(chem_avg$DOC_mgL_1[chem_avg$msn==2]),
                                  mean(chem_avg$DOC_mgL_1[chem_avg$msn==3]),
                                  mean(chem_avg$DOC_mgL_1[chem_avg$msn==4])), 
                          "secchi" = as.numeric(secchi$Secchi_m),
                          "thermo_depth" = c(mean(ctd_thermo_depth$therm_depth[ctd_thermo_depth$msn==1]),
                                             mean(ctd_thermo_depth$therm_depth[ctd_thermo_depth$msn==2]),
                                             mean(ctd_thermo_depth$therm_depth[ctd_thermo_depth$msn==3]),
                                             mean(ctd_thermo_depth$therm_depth[ctd_thermo_depth$msn==4])))


#join driver and env dfs
zoops_plus_daily_drivers <- left_join(zoop_dens, daily_drivers)

#fit environmental drivers onto ordination
fit_env <- envfit(ord$sites, zoops_plus_daily_drivers[,c(25:43)])

#pull out vectors
scores <- as.data.frame(scores(fit_env, display = "vectors"))
scores <- cbind(scores, env = rownames(scores))

#plot drivers w/ NMDS
ord <- ordiplot(NMDS_dens_bray,display = c('sites','species'),choices = c(1,2),type = "n")

day <- gg_ordiplot(ord, zoop_dens$collect_date, kind = "ehull", 
                   ellipse=FALSE, hull = TRUE, plot = FALSE, pt.size=0.9) 

NMDS_day_env <- day$plot + geom_point() + theme_bw() + 
  geom_polygon(data = day$df_hull, aes(x = x, y = y, fill = Group), alpha=0.2) +
  geom_point(data=day$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2,
             fill = c("#111258","#4C1559","#B1235E","#F0855F","#F8C48A","#F8F4A9")) +
  theme(text = element_text(size=7), axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), 
        legend.key.height=unit(0.3,"line"),
        legend.key = element_blank(),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-0,-0,-0,-0),
        legend.direction = "horizontal",
        legend.text = element_text(size=9),
        axis.text.x = element_text(vjust = 0.5), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "top", legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(1,0,0,0), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) + 
  guides(fill="none", color = guide_legend(override.aes = list(size = 2))) +
  #annotate("text", x=-0.42, y=0.6, label= "bolditalic(a)", parse=T, size = 3) +
  scale_fill_manual("",values=c("#111258","#4C1559","#B1235E","#F0855F","#F8C48A","#F8F4A9"))+
  scale_color_manual("",values=c("#111258","#4C1559","#B1235E","#F0855F","#F8C48A","#F8F4A9"),
                     label=c("2020-09-11", "2020-09-15", "2021-06-10", "2021-06-11",
                             "2022-06-30", "2022-07-01"))  +
  geom_segment(data = scores,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), linewidth= 0.3,
               arrow = arrow(length = unit(0.1, "cm")), colour = "black") +
  geom_text_repel(data=scores, aes(x = NMDS1, y = NMDS2, label = env), size=2)
#ggsave("Summer2021-DataAnalysis/Figures/fcr_mom/NMDS_fcr_dens_2v1_envfit_days.jpg", NMDS_day_env, width=3, height=4) 


#------------------------------------------------------------------------------
#env fit with depth-specific drivers

depth_drivers  <- data.frame("msn" = c(rep(as.character(1),10), rep(as.character(2),10),
                                         rep(as.character(3),10), rep(as.character(4),10)), #mean full water column vals
                             "depth" = c(0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                             "temp" = c(ctd.final_avg$Temp_C_1[ctd.final_avg$msn==1],
                                        ctd.final_avg$Temp_C_1[ctd.final_avg$msn==2],
                                        ctd.final_avg$Temp_C_1[ctd.final_avg$msn==3],
                                        ctd.final_avg$Temp_C_1[ctd.final_avg$msn==4]), 
                             "DO_conc" = c(ctd.final_avg$DO_mgL_1[ctd.final_avg$msn==1],
                                           ctd.final_avg$DO_mgL_1[ctd.final_avg$msn==2],
                                           ctd.final_avg$DO_mgL_1[ctd.final_avg$msn==3],
                                           ctd.final_avg$DO_mgL_1[ctd.final_avg$msn==4]),
                             "DOsat" = c(ctd.final_avg$DOsat_percent_1[ctd.final_avg$msn==1],
                                         ctd.final_avg$DOsat_percent_1[ctd.final_avg$msn==2],
                                         ctd.final_avg$DOsat_percent_1[ctd.final_avg$msn==3],
                                         ctd.final_avg$DOsat_percent_1[ctd.final_avg$msn==4]),
                             "Spcond" = c(ctd.final_avg$SpCond_uScm_1[ctd.final_avg$msn==1],
                                          ctd.final_avg$SpCond_uScm_1[ctd.final_avg$msn==2],
                                          ctd.final_avg$SpCond_uScm_1[ctd.final_avg$msn==3],
                                          ctd.final_avg$SpCond_uScm_1[ctd.final_avg$msn==4]),
                             "Green_algae" = c(c(NA,NA,fp_final$GreenAlgae_ugL[fp_final$msn==1]),
                                               c(NA,fp_final$GreenAlgae_ugL[fp_final$msn==2]),
                                               fp_final$GreenAlgae_ugL[fp_final$msn==3],
                                               fp_final$GreenAlgae_ugL[fp_final$msn==4]),
                             "Bluegreen_algae" = c(c(NA,NA,fp_final$Bluegreens_ugL[fp_final$msn==1]),
                                                   c(NA,fp_final$Bluegreens_ugL[fp_final$msn==2]),
                                                   fp_final$Bluegreens_ugL[fp_final$msn==3],
                                                   fp_final$Bluegreens_ugL[fp_final$msn==4]),
                             "Brown_algae" = c(c(NA,NA,fp_final$BrownAlgae_ugL[fp_final$msn==1]),
                                               c(NA,fp_final$BrownAlgae_ugL[fp_final$msn==2]),
                                               fp_final$BrownAlgae_ugL[fp_final$msn==3],
                                               fp_final$BrownAlgae_ugL[fp_final$msn==4]),
                             "Mixed_algae" = c(c(NA,NA,fp_final$MixedAlgae_ugL[fp_final$msn==1]),
                                               c(NA,fp_final$MixedAlgae_ugL[fp_final$msn==2]),
                                               fp_final$MixedAlgae_ugL[fp_final$msn==3],
                                               fp_final$MixedAlgae_ugL[fp_final$msn==4]),
                             "Yellow_substances" = c(c(NA,NA,fp_final$YellowSubstances_ugL[fp_final$msn==1]),
                                                     c(NA,fp_final$YellowSubstances_ugL[fp_final$msn==2]),
                                                     fp_final$YellowSubstances_ugL[fp_final$msn==3],
                                                     fp_final$YellowSubstances_ugL[fp_final$msn==4]),
                             "Total_algae" = c(c(NA,NA,fp_final$TotalConc_ugL[fp_final$msn==1]),
                                               c(NA,fp_final$TotalConc_ugL[fp_final$msn==2]),
                                               fp_final$TotalConc_ugL[fp_final$msn==3],
                                               fp_final$TotalConc_ugL[fp_final$msn==4]),
                             "PAR" = c(ctd.final_avg$PAR_umolm2s_1[ctd.final_avg$msn==1],
                                       ctd.final_avg$PAR_umolm2s_1[ctd.final_avg$msn==2],
                                       ctd.final_avg$PAR_umolm2s_1[ctd.final_avg$msn==3],
                                       ctd.final_avg$PAR_umolm2s_1[ctd.final_avg$msn==4]))


#join driver and env dfs
zoops_plus_depth_drivers <- left_join(zoop_dens, depth_drivers)

#fit environmental drivers onto ordination
fit_env <- envfit(ord$sites, zoops_plus_depth_drivers[,c(25:35)], na.rm=T) #some missing values so make sure can remove them like this

#pull out vectors
scores <- as.data.frame(scores(fit_env, display = "vectors"))
scores <- cbind(scores, env = rownames(scores))

#plot drivers w/ NMDS
ord <- ordiplot(NMDS_dens_bray,display = c('sites','species'),choices = c(1,2),type = "n")
depth <- gg_ordiplot(ord, zoop_dens$depth, kind = "ehull", 
                     ellipse=FALSE, hull = TRUE, plot = FALSE, pt.size=0.9) 

NMDS_depth_env <- depth$plot + geom_point() + theme_bw() + 
  geom_polygon(data = depth$df_hull, aes(x = x, y = y, fill = Group), alpha=0.2) +
  geom_point(data=depth$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2,
             fill = brewer.pal(10, "BrBG")) +
  theme(text = element_text(size=7), axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), 
        legend.key.height=unit(0.5,"line"),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(0,-0,-0,-0),
        legend.direction = "horizontal",
        legend.text = element_text(size=9),
        axis.text.x = element_text(vjust = 0.5), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "top",
        plot.margin = unit(c(1,-0,0,0), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.5,"line")) +
  guides(fill="none", color = guide_legend(override.aes = list(size = 2))) +
  #   annotate("text", x=-0.42, y=0.6, label= "bolditalic(a)", parse=T, size = 3) +
  scale_fill_manual("",values=brewer.pal(10, "BrBG"))+
  scale_color_manual("",values=brewer.pal(10, "BrBG"),
                     label=c(unique(zoop_biom$depth))) +
  geom_segment(data = scores,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), linewidth= 0.3,
               arrow = arrow(length = unit(0.1, "cm")), colour = "black") +
  geom_text_repel(data=scores, aes(x = NMDS1, y = NMDS2, label = env), size=2)
#ggsave("Summer2021-DataAnalysis/Figures/fcr_mom/NMDS_fcr_dens_2v1_envfit_depths.jpg", NMDS_depth_env, width=3, height=4) 



