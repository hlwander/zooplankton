#Code for visualizing fluoroprobe casts from 2016, 2019, and 2020 MSNs
#Created 24 Dec 2020; modified from MEL and RPM code

#load packages
pacman::p_load(tidyverse, lubridate, akima, reshape2, 
               gridExtra, grid, colorRamps,RColorBrewer, rLakeAnalyzer, cowplot)

#read in fp casts
col_names <- names(read_tsv("./RawData/FluoroprobeData/20160803_BVR_50_a.txt",n_max=0))

raw_fp <- dir(path = "./RawData/FluoroprobeData", pattern = "_BVR_50") %>% 
  map_df(~ read_tsv(file.path(path = "./RawData/FluoroprobeData", .), col_types = cols(.default = "c"), col_names = col_names, skip = 2))

fp <- raw_fp %>%
  mutate(DateTime = `Date/Time`, GreenAlgae_ugL = as.numeric(`Green Algae`), Bluegreens_ugL = as.numeric(`Bluegreen`),
         Browns_ugL = as.numeric(`Diatoms`), Mixed_ugL = as.numeric(`Cryptophyta`), YellowSubstances_ugL = as.numeric(`Yellow substances`),
         TotalConc_ugL = as.numeric(`Total conc.`), Transmission_perc = as.numeric(`Transmission`), Depth_m = as.numeric(`Depth`)) %>%
  select(DateTime, GreenAlgae_ugL, Bluegreens_ugL, Browns_ugL, Mixed_ugL, YellowSubstances_ugL,
         TotalConc_ugL, Transmission_perc, Depth_m) %>%
  mutate(DateTime = as.POSIXct(as_datetime(DateTime, tz = "", format = "%m/%d/%Y %I:%M:%S %p"))) %>%
  mutate(Date = date(DateTime), DOY = yday(DateTime))

#not sure why 8-17-20 is in the 8-12-20 file, but drop this date
fp <- fp[fp$Date!="2020-08-17",]

#round to nearest hour
fp$DateTime <- ceiling_date(fp$DateTime,"hour")

# filter out depths in the fp cast that are closest to specified values.
  depths = seq(0.1, 10.3, by = 0.3)
  df.final<-data.frame()
  
  for (i in 1:length(depths)){
    fp_layer <- fp %>% group_by(DateTime) %>% slice(which.min(abs(as.numeric(Depth_m) - depths[i])))
    # Bind each of the data layers together.
    df.final = bind_rows(df.final, fp_layer)
  }

  # Re-arrange the data frames by date
  fp_new <- arrange(df.final, DateTime)

  # Round each extracted depth to the nearest 10th. 
  fp_new$Depth_m <- round(as.numeric(fp_new$Depth_m), digits = 0.5)
  
  #save fp data
  write.csv(fp_new, paste0("./RawData/Fluora.csv"))
#-----------------------------------------------------------------------------#  
#look at casts
for (i in 1:length(unique(fp_new$DateTime))){
  profile = subset(fp_new, DateTime == unique(fp_new$DateTime)[i])
  castname = profile$DateTime[1]
  fp_casts = profile %>%
    select(Depth_m, GreenAlgae_ugL, Bluegreens_ugL, Browns_ugL, Mixed_ugL, TotalConc_ugL)%>%
    gather(GreenAlgae_ugL:TotalConc_ugL, key = spectral_group, value = ugL)
  profile_plot <- ggplot(data = fp_casts, aes(x = ugL, y = Depth_m, group = spectral_group, color = spectral_group))+
    geom_path(size = 1)+  scale_y_reverse()+ ggtitle(castname)+ theme_bw()
  filename = paste0(getwd(),"/RawData/FluoroprobeData/Figures/",castname,".png")
  ggsave(filename = filename, plot = profile_plot, device = "png")
  
}

  #casts look similar with cmax between 6.5 and 7.5m except for 2019-07-11 third cast and 2020-08-12
  #these cases actually have brown and green algae peaks (respectively) at weird depths instead of bluegreen peaks at ~7m
  #figure out what times these are - assuming the 2020 one is noon ish
  
  
  