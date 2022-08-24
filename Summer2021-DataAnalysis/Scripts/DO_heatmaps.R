library(tidyverse)
library(akima)
library(colorRamps)
library(lubridate)
library(ggpubr)

on <- as.Date(c("2019-06-03", "2019-07-08", "2019-08-05","2019-09-02","2020-06-29","2020-07-13","2020-07-23","2020-09-25", "2021-06-11","2021-07-12","2021-07-26","2021-09-14"))
off <- as.Date(c("2019-06-17", "2019-07-18", "2019-08-19",NA,"2020-07-12","2020-07-22","2020-09-11","2020-12-10", "2021-06-26","2021-07-14","2021-09-06",NA))
on_off = data.frame(on,off,Reservoir="FCR")%>%
  mutate(Year = year(on))

#Load data
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/200/12/0a62d1946e8d9a511bc1404e69e59b8c" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")


dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               , col.names=c(
                 "Reservoir",     
                 "Site",     
                 "Date",     
                 "Depth_m",     
                 "Temp_C",     
                 "DO_mgL",     
                 "DO_pSat",     
                 "Cond_uScm",     
                 "Spec_Cond_uScm",     
                 "Chla_ugL",     
                 "Turb_NTU",     
                 "pH",     
                 "ORP_mV",     
                 "PAR_umolm2s",     
                 "Desc_rate",     
                 "Flag_Temp",     
                 "Flag_DO",     
                 "Flag_Cond",     
                 "Flag_SpecCond",     
                 "Flag_Chla",     
                 "Flag_Turb",     
                 "Flag_pH",     
                 "Flag_ORP",     
                 "Flag_PAR",     
                 "Flag_DescRate"    ), check.names=TRUE)

unlink(infile1)

do = dt1%>%
  mutate(Date = as.Date(Date))%>%
  filter(!is.na(DO_mgL),
         Site == 50,
         year(Date)!=2021|
           (Date>=as.Date("2021-03-01")&Date<as.Date("2021-07-01")),
         Date != as.Date("2019-07-18")|Reservoir=="FCR")
#Trim dataset
depths <- seq(0,11, by = .1)
newDepths <- depths
df.final<- do %>% group_by(Date, Reservoir) %>% slice(which.min(abs(as.numeric(Depth_m) - depths[1]))) #Create a new dataframe
df.final$Depth_m <- newDepths[1]
for (i in 2:length(depths)){ #loop through all depths and add the closest values to the final dataframe
  ctd_atThisDepth <- do %>% group_by(Date, Reservoir) %>% slice(which.min(abs(as.numeric(Depth_m) - depths[i])))
  ctd_atThisDepth$Depth_m <- newDepths[i]
  df.final <- rbind(df.final,ctd_atThisDepth)
}


inUrl2  <- "https://pasta.lternet.edu/package/data/eml/edi/198/9/b3bd353312f9e37ca392e2a5315cc9da" 
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

unlink(infile2)

df.final2 = df.final%>%
  full_join(ysi%>%select("Reservoir",     
                         "Site",     
                         "DateTime",     
                         "Depth_m",  
                         "DO_mgL")%>%
              filter(!is.na(DO_mgL),
                     Site == 50)%>%
              mutate(Date = as.Date(DateTime))%>%
              filter((Date>=as.Date("2021-07-01")&Date<as.Date("2021-12-01"))|
                       (Date>=as.Date("2019-03-01")&Date<as.Date("2019-12-01"))|
                       (Date>=as.Date("2017-05-01")&Date<as.Date("2017-11-01"))|
                       (Date>=as.Date("2020-03-01")&Date<as.Date("2020-06-18"))
              ))

ysi%>%
  filter(Site==50)%>%
  ggplot(aes(x = as.Date(DateTime), y = DO_mgL))+
  geom_point()

interp_do = data.frame(x = NA, y = NA, z = NA, Reservoir = NA)
years = c(2020:2021) #change this range for DO before 2020
for(year in years){
  do_fcr<- df.final2%>%
    mutate(Year = year(Date))%>%
    filter(Year == year,
           Reservoir == "FCR",
           Depth_m <= 9)
  do_bvr<- df.final2%>%
    mutate(Year = year(Date))%>%
    filter(Year == year,
           Reservoir == "BVR",
           Depth_m <= 11)
  do_fcr_interp <-
    interp2xyz(interp(do_fcr$Date,do_fcr$Depth_m,do_fcr$DO_mgL,xo = seq(min(do_fcr$Date), max(do_fcr$Date),1 ),yo = seq(min(do_fcr$Depth_m), max(do_fcr$Depth_m), by = .01), duplicate = "mean"),data.frame = T)%>%
    mutate(Reservoir = "FCR")
  do_bvr_interp <- interp2xyz(interp(do_bvr$Date,do_bvr$Depth_m,do_bvr$DO_mgL,xo = seq(min(do_bvr$Date), max(do_bvr$Date), 1),yo = seq(min(do_bvr$Depth_m), max(do_bvr$Depth_m), by = .01), duplicate = "mean"),data.frame = T)%>%
    mutate(Reservoir = "BVR")
  interp_do = interp_do%>%
    full_join(do_fcr_interp)%>%
    full_join(do_bvr_interp)
}
interp_do = interp_do%>%
  filter(!is.na(Reservoir))

#subset df to only include FCR data
interp_do <- interp_do[interp_do$Reservoir=="FCR" & !is.na(interp_do$z),]

#noon/midnight zoop sampling days
fcr_zoops <- c("2020-09-10","2020-09-14","2021-06-10")

#add manual groups for heatmap 
interp_do$groups <- cut(interp_do$z,    
                       breaks = c(0, 2, 4, 6, 8, 10, 12, 14, 16))


#DO heatmap
#jpeg("./Figures/FCR_DO_discrete_heatmap_2020_2021.jpeg",width = 8, height = 4, units = "in", res = 300)
interp_do%>%
  mutate(x = as.Date(x, origin = "1970-01-01"),
         Year = year(x))%>%
  ggplot(aes(x=x, y=y,fill=groups))+
  geom_tile()+ #geom_raster(aes(fill=z)) + #
  scale_y_reverse(expand = c(0,0))+
  labs(x = "", y = "Depth (m)", title = "Dissolved Oxygen ",fill=expression('(mg/L)'))+
  scale_fill_manual(breaks=levels(interp_do$groups),values=c("#CC3300","#FF9900","#CCFF00","#66FF66","#66FF99","#33FFFF","#3399FF","#0066CC"))+ #scale_fill_gradientn(colors=rev(blue2green2red(100)))+
  geom_vline(aes(xintercept = on), data= on_off)+
  geom_vline(aes(xintercept = off),lty=2, data = on_off)+
  geom_vline(xintercept = as.Date(fcr_zoops),lty=1,cex=0.6, col="dark blue")+
  scale_x_date(expand = c(0,0))+
  theme(panel.border = element_rect(fill = NA))+
  facet_grid(#rows= vars(Reservoir),
             cols = vars(Year),
             scales = "free",
             space = "free_x"
  )
#dev.off()
