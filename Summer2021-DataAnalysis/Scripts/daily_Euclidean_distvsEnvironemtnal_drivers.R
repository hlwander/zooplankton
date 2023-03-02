#plotting Euclidean distances vs. environemtnal drivers

#read in Euclidean distance df
read.csv(file.path(getwd(),"/Summer2021-DataAnalysis/SummaryStats/Daily_Euclidean_distances_pelvslit.csv"))

# read in all driver data (YSI, CTD, flora, chem?)
#day1 (10-11 Jul 2019), day2 (24-25 Jul 2019), day3 (12-13 Aug 2020), day4 (15-16 Jun 2021), day5 (7-8 Jul 2021)

#df to compile euclidean distances for each day/site
euclidean_drivers_df <- data.frame(site = c("day1","day2","day3","day4","day5"),
                                   pel_euc_dist = euc_distances_df$pelagic,
                                   lit_euc_dist = euc_distances_df$littoral)


#EDI YSI data
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

unlink(infile2)

#some of the sites are wrong for 7-8 Jul 2021
ysi$Site[ysi$DateTime=="2021-07-07 11:37:00"] <- 50
ysi$Site[ysi$DateTime=="2021-07-08 11:45:00"] <- 50

#change date format
ysi$DateTime <- as.Date(ysi$DateTime)

#separate ysi by day
#msn1_ysi <-  ysi[(ysi$DateTime=="2019-07-10" | ysi$DateTime=="2019-07-11")  & ysi$Reservoir=="BVR",] #no ysi data
#msn2_ysi <-  ysi[(ysi$DateTime=="2019-07-24" | ysi$DateTime=="2019-07-25")  & ysi$Reservoir=="BVR",] #no ysi data
msn3_ysi <-  ysi[(ysi$DateTime=="2020-08-12" | ysi$DateTime=="2020-08-13")  & ysi$Reservoir=="BVR",]
msn4_ysi <-  ysi[(ysi$DateTime=="2021-06-15" | ysi$DateTime=="2020-06-16")  & ysi$Reservoir=="BVR",]
msn5_ysi <-  ysi[(ysi$DateTime=="2021-07-07" | ysi$DateTime=="2021-07-08")  & ysi$Reservoir=="BVR" & ysi$Site==50,]

msn5_ysi_lit <- ysi[(ysi$DateTime=="2021-07-07" | ysi$DateTime=="2021-07-08")  & ysi$Reservoir=="BVR" & ysi$Site==51,]

#read in CTD EDI file
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/200/12/0a62d1946e8d9a511bc1404e69e59b8c" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")

ctd <-read.csv(infile1,header=F 
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

#change date format
ctd$Date <- as.Date(ctd$Date)

#read in DO and temp from sampling days, select 0.1 and 9m, and then average both day DO
msn1_ctd <- ctd %>% filter((Date=="2019-07-10" | Date=="2019-07-11") & Reservoir=="BVR") %>%
  group_by(Date) %>% filter(abs(Depth_m - 0.1) == min(abs(Depth_m - 0.1)) | (abs(Depth_m - 9) == min(abs(Depth_m - 9)))) %>%
  distinct(Depth_m, .keep_all = TRUE)

msn2_ctd <- ctd %>% filter((Date=="2019-07-24" | Date=="2019-07-25") & Reservoir=="BVR") %>%
  group_by(Date) %>% filter(abs(Depth_m - 0.1) == min(abs(Depth_m - 0.1)) | (abs(Depth_m - 9) == min(abs(Depth_m - 9)))) %>%
  distinct(Depth_m, .keep_all = TRUE)

msn3_ctd <- ctd %>% filter(Date=="2020-08-12" | Date=="2020-08-13" & Reservoir=="BVR") %>%
  group_by(Date) %>% filter(abs(Depth_m - 0.1) == min(abs(Depth_m - 0.1)) | (abs(Depth_m - 9) == min(abs(Depth_m - 9)))) %>%
  distinct(Depth_m, .keep_all = TRUE)

msn4_ctd <- ctd %>% filter(Date=="2021-06-15" | Date=="2021-06-16" & Reservoir=="BVR") %>%
  group_by(Date) %>% filter(abs(Depth_m - 0.1) == min(abs(Depth_m - 0.1)) | (abs(Depth_m - 9) == min(abs(Depth_m - 9)))) %>%
  distinct(Depth_m, .keep_all = TRUE)

before_msn5_ctd <- ctd[(ctd$Date=="2021-06-28")  & ctd$Reservoir=="BVR",] #no ctd on day 5
before_msn5_ctd <-  before_msn5_ctd[before_msn5_ctd$Depth_m==closest(before_msn5_ctd$Depth_m, 0.1) | before_msn5_ctd$Depth_m==closest(before_msn5_ctd$Depth_m, 9),] 

after_msn5_ctd <- ctd[(ctd$Date=="2021-07-12")  & ctd$Reservoir=="BVR",] #no ctd on day 5
after_msn5_ctd <-  after_msn5_ctd[after_msn5_ctd$Depth_m==closest(after_msn5_ctd$Depth_m, 0.1) | after_msn5_ctd$Depth_m==closest(after_msn5_ctd$Depth_m, 9),] 

#so for now, I'm just going to take the do from 4 days after because I'm assuming that DO doesn't change that much in 4 days...
msn5_0.1_ctd <- after_msn5_ctd$DO_mgL[after_msn5_ctd$Depth_m<2]
msn5_9_ctd <- after_msn5_ctd$DO_mgL[after_msn5_ctd$Depth_m>2]

#select every 0.5m from casts
ctd_final <- ctd %>%
  dplyr::mutate(rdepth = plyr::round_any(Depth_m, 0.5)) %>% 
  dplyr::group_by(Date, rdepth, Reservoir, Site) %>%
  dplyr::summarise(value = mean(Temp_C)) %>% 
  dplyr::rename(depth = rdepth) 

#calcualte thermocline depth now
msn1_therm_depth <- ctd_final %>% filter((Date=="2019-07-10" | Date=="2019-07-11") & Reservoir=="BVR") %>%
  group_by(Date) %>% filter(depth > 0) %>% mutate(therm_depth = thermo.depth(value,depth))

msn2_therm_depth <- ctd_final %>% filter((Date=="2019-07-24"  | Date=="2019-07-25") & Reservoir=="BVR") %>%
  group_by(Date) %>% filter(depth > 0) %>% mutate(therm_depth = thermo.depth(value,depth))

msn3_therm_depth <- ctd_final %>% filter((Date=="2020-08-12"  | Date=="2020-08-13") & Reservoir=="BVR") %>%
  group_by(Date) %>% filter(depth > 0) %>% mutate(therm_depth = thermo.depth(value,depth))

msn4_therm_depth <- ctd_final %>% filter((Date=="2021-06-15"  | Date=="2021-06-16") & Reservoir=="BVR") %>%
  group_by(Date) %>% filter(depth > 0) %>% mutate(therm_depth = thermo.depth(value,depth))

msn5_therm_depth <- ctd_final %>% filter((Date=="2021-07-12"  | Date=="2021-07-12") & Reservoir=="BVR") %>%
  group_by(Date) %>% filter(depth > 0) %>% mutate(therm_depth = thermo.depth(value,depth)) #4 days after msn


#read in fp data from edi
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/272/6/6b3151c0fdd913e02641363c2b00ae57" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")


fp <-read.csv(infile1,header=F 
              ,skip=1
              ,sep=","  
              , col.names=c(
                "Reservoir",     
                "Site",     
                "DateTime",     
                "Depth_m",     
                "GreenAlgae_ugL",     
                "Bluegreens_ugL",     
                "BrownAlgae_ugL",     
                "MixedAlgae_ugL",     
                "YellowSubstances_ugL",     
                "TotalConc_ugL",     
                "Transmission",     
                "Temp_degC",     
                "RFU_525nm",     
                "RFU_570nm",     
                "RFU_610nm",     
                "RFU_370nm",     
                "RFU_590nm",     
                "RFU_470nm",     
                "Flag_GreenAlgae",     
                "Flag_BluegreenAlgae",     
                "Flag_BrownAlgae",     
                "Flag_MixedAlgae",     
                "Flag_TotalConc",     
                "Flag_Temp",     
                "Flag_Transmission",     
                "Flag_525nm",     
                "Flag_570nm",     
                "Flag_610nm",     
                "Flag_370nm",     
                "Flag_590nm",     
                "Flag_470nm"    ), check.names=TRUE)

unlink(infile1)

#change date format for fp data
fp$DateTime <- as.Date(fp$DateTime)

#fp for all msns
msn1_fp <- fp %>% filter((DateTime=="2019-07-10" | DateTime=="2019-07-11") & Reservoir=="BVR") %>%
  group_by(DateTime) %>% filter(abs(Depth_m - 0.1) == min(abs(Depth_m - 0.1)) | (abs(Depth_m - 9) == min(abs(Depth_m - 9)))) %>%
  distinct(Depth_m, .keep_all = TRUE)

msn2_fp <- fp %>% filter((DateTime=="2019-07-24" | DateTime=="2019-07-25") & Reservoir=="BVR") %>%
  group_by(DateTime) %>% filter(abs(Depth_m - 0.1) == min(abs(Depth_m - 0.1)) | (abs(Depth_m - 9) == min(abs(Depth_m - 9)))) %>%
  distinct(Depth_m, .keep_all = TRUE) #nope

msn3_fp <- fp %>% filter((DateTime=="2020-08-12" | DateTime=="2020-08-13") & Reservoir=="BVR") %>%
  group_by(DateTime) %>% filter(abs(Depth_m - 0.1) == min(abs(Depth_m - 0.1)) | (abs(Depth_m - 9) == min(abs(Depth_m - 9)))) %>%
  distinct(Depth_m, .keep_all = TRUE)

msn4_fp <- fp %>% filter((DateTime=="2021-06-15" | DateTime=="2021-06-16") & Reservoir=="BVR") %>%
  group_by(DateTime) %>% filter(abs(Depth_m - 0.1) == min(abs(Depth_m - 0.1)) | (abs(Depth_m - 9) == min(abs(Depth_m - 9)))) %>%
  distinct(Depth_m, .keep_all = TRUE)

msn5_fp <- fp %>% filter((DateTime=="2021-07-07" | DateTime=="2021-07-08") & Reservoir=="BVR") %>%
  group_by(DateTime) %>% filter(abs(Depth_m - 0.1) == min(abs(Depth_m - 0.1)) | (abs(Depth_m - 9) == min(abs(Depth_m - 9)))) %>%
  distinct(Depth_m, .keep_all = TRUE) #nope

#read in nutreint data from edi
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/199/10/aa2ccc23688fc908f9d61cb217210a3d" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")


chem <-read.csv(infile1,header=F 
                ,skip=1
                ,sep=","  
                , col.names=c(
                  "Reservoir",     
                  "Site",     
                  "DateTime",     
                  "Depth_m",     
                  "Rep",     
                  "TN_ugL",     
                  "TP_ugL",     
                  "NH4_ugL",     
                  "NO3NO2_ugL",     
                  "SRP_ugL",     
                  "DOC_mgL",     
                  "DIC_mgL",     
                  "DC_mgL",     
                  "DN_mgL",     
                  "Flag_DateTime",     
                  "Flag_TN",     
                  "Flag_TP",     
                  "Flag_NH4",     
                  "Flag_NO3NO2",     
                  "Flag_SRP",     
                  "Flag_DOC",     
                  "Flag_DIC",     
                  "Flag_DC",     
                  "Flag_DN"    ), check.names=TRUE)

unlink(infile1)

#change date format for fp data
chem$DateTime <- as.Date(chem$DateTime)

#chem for all msns
msn1_chem <- chem %>% filter((DateTime=="2019-07-10" | DateTime=="2019-07-11") & Reservoir=="BVR") %>%
  group_by(DateTime) %>% filter(Depth_m== 0.1 | Depth_m== 9)

msn2_chem <- chem %>% filter((DateTime=="2019-07-24" | DateTime=="2019-07-25") & Reservoir=="BVR") %>%
  group_by(DateTime) %>% filter(Depth_m== 0.1 | Depth_m== 9)

msn3_chem <- chem %>% filter((DateTime=="2020-08-12" | DateTime=="2020-08-13") & Reservoir=="BVR") %>%
  group_by(DateTime) %>% filter(Depth_m== 0.1 | Depth_m== 9)

msn4_chem <- chem %>% filter((DateTime=="2021-06-15" | DateTime=="2021-06-16") & Reservoir=="BVR" & Site==50) %>%
  group_by(DateTime) %>% filter(Depth_m== 0.1 | Depth_m== 9)

msn5_chem <- chem %>% filter((DateTime=="2021-07-11" | DateTime=="2021-07-12") & Reservoir=="BVR" & Site==50) %>%
  group_by(DateTime) %>% filter(Depth_m== 0.1 | Depth_m== 9)

#add environmental drivers to euclidean distance df
euclidean_drivers_df$DO_0.1m <- c(mean(msn1_ctd$DO_mgL[msn1_ctd$Depth_m<2]), mean(msn2_ctd$DO_mgL[msn2_ctd$Depth_m<2]),
                                  mean(msn3_ctd$DO_mgL[msn3_ctd$Depth_m<2]),mean(msn4_ctd$DO_mgL[msn4_ctd$Depth_m<2]),msn5_0.1_ctd)

euclidean_drivers_df$DO_9.0m <- c(mean(msn1_ctd$DO_mgL[msn1_ctd$Depth_m>2]), mean(msn2_ctd$DO_mgL[msn2_ctd$Depth_m>2]),
                                  mean(msn3_ctd$DO_mgL[msn3_ctd$Depth_m>2]),mean(msn4_ctd$DO_mgL[msn4_ctd$Depth_m>2]),msn5_9_ctd)

euclidean_drivers_df$DOsat_0.1m <- c(mean(msn1_ctd$DO_pSat[msn1_ctd$Depth_m<2]), mean(msn2_ctd$DO_pSat[msn2_ctd$Depth_m<2]),
                                     mean(msn3_ctd$DO_pSat[msn3_ctd$Depth_m<2]),mean(msn4_ctd$DO_pSat[msn4_ctd$Depth_m<2]),
                                     after_msn5_ctd$DO_pSat[after_msn5_ctd$Depth_m<2])

euclidean_drivers_df$DOsat_9.0m <- c(mean(msn1_ctd$DO_pSat[msn1_ctd$Depth_m>2]), mean(msn2_ctd$DO_pSat[msn2_ctd$Depth_m>2]),
                                     mean(msn3_ctd$DO_pSat[msn3_ctd$Depth_m>2]),mean(msn4_ctd$DO_pSat[msn4_ctd$Depth_m>2]),
                                     after_msn5_ctd$DO_pSat[after_msn5_ctd$Depth_m>2])

euclidean_drivers_df$Spcond_0.1m <- c(mean(msn1_ctd$Spec_Cond_uScm[msn1_ctd$Depth_m<2]), mean(msn2_ctd$Spec_Cond_uScm[msn2_ctd$Depth_m<2]),
                                      mean(msn3_ctd$Spec_Cond_uScm[msn3_ctd$Depth_m<2]),mean(msn4_ctd$Spec_Cond_uScm[msn4_ctd$Depth_m<2]),
                                      after_msn5_ctd$Spec_Cond_uScm[after_msn5_ctd$Depth_m<2])

euclidean_drivers_df$Spcond_9.0m <- c(mean(msn1_ctd$Spec_Cond_uScm[msn1_ctd$Depth_m>2]), mean(msn2_ctd$Spec_Cond_uScm[msn2_ctd$Depth_m>2]),
                                      mean(msn3_ctd$Spec_Cond_uScm[msn3_ctd$Depth_m>2]),mean(msn4_ctd$Spec_Cond_uScm[msn4_ctd$Depth_m>2]),
                                      after_msn5_ctd$Spec_Cond_uScm[after_msn5_ctd$Depth_m>2])

euclidean_drivers_df$chla_0.1m <- c(mean(msn1_ctd$Chla_ugL[msn1_ctd$Depth_m<2]), mean(msn2_ctd$Chla_ugL[msn2_ctd$Depth_m<2]),
                                    mean(msn3_ctd$Chla_ugL[msn3_ctd$Depth_m<2]),mean(msn4_ctd$Chla_ugL[msn4_ctd$Depth_m<2]),
                                    after_msn5_ctd$Chla_ugL[after_msn5_ctd$Depth_m<2])

euclidean_drivers_df$chla_9.0m <- c(mean(msn1_ctd$Chla_ugL[msn1_ctd$Depth_m>2]), mean(msn2_ctd$Chla_ugL[msn2_ctd$Depth_m>2]),
                                    mean(msn3_ctd$Chla_ugL[msn3_ctd$Depth_m>2]),mean(msn4_ctd$Chla_ugL[msn4_ctd$Depth_m>2]),
                                    after_msn5_ctd$Chla_ugL[after_msn5_ctd$Depth_m>2])

euclidean_drivers_df$turb_0.1m <- c(mean(msn1_ctd$Turb_NTU[msn1_ctd$Depth_m<2]), mean(msn2_ctd$Turb_NTU[msn2_ctd$Depth_m<2]),
                                    mean(msn3_ctd$Turb_NTU[msn3_ctd$Depth_m<2]),mean(msn4_ctd$Turb_NTU[msn4_ctd$Depth_m<2]),
                                    after_msn5_ctd$Turb_NTU[after_msn5_ctd$Depth_m<2])

euclidean_drivers_df$turb_9.0m <- c(mean(msn1_ctd$Turb_NTU[msn1_ctd$Depth_m>2]), mean(msn2_ctd$Turb_NTU[msn2_ctd$Depth_m>2]),
                                    mean(msn3_ctd$Turb_NTU[msn3_ctd$Depth_m>2]),mean(msn4_ctd$Turb_NTU[msn4_ctd$Depth_m>2]),
                                    after_msn5_ctd$Turb_NTU[after_msn5_ctd$Depth_m>2])

euclidean_drivers_df$temp_0.1m <- c(mean(msn1_ctd$Temp_C[msn1_ctd$Depth_m<2]), mean(msn2_ctd$Temp_C[msn2_ctd$Depth_m<2]),
                                    mean(msn3_ctd$Temp_C[msn3_ctd$Depth_m<2]),mean(msn4_ctd$Temp_C[msn4_ctd$Depth_m<2]),
                                    after_msn5_ctd$Temp_C[after_msn5_ctd$Depth_m<2])

euclidean_drivers_df$temp_9.0m <- c(mean(msn1_ctd$Temp_C[msn1_ctd$Depth_m>2]), mean(msn2_ctd$Temp_C[msn2_ctd$Depth_m>2]),
                                    mean(msn3_ctd$Temp_C[msn3_ctd$Depth_m>2]),mean(msn4_ctd$Temp_C[msn4_ctd$Depth_m>2]),
                                    after_msn5_ctd$Temp_C[after_msn5_ctd$Depth_m>2])

euclidean_drivers_df$par_0.1m <- c(mean(msn1_ctd$PAR_umolm2s[msn1_ctd$Depth_m<2]), mean(msn2_ctd$PAR_umolm2s[msn2_ctd$Depth_m<2]),
                                   mean(msn3_ctd$PAR_umolm2s[msn3_ctd$Depth_m<2]),mean(msn4_ctd$PAR_umolm2s[msn4_ctd$Depth_m<2]),
                                   after_msn5_ctd$PAR_umolm2s[after_msn5_ctd$Depth_m<2])

euclidean_drivers_df$par_9.0m <- c(mean(msn1_ctd$PAR_umolm2s[msn1_ctd$Depth_m>2]), mean(msn2_ctd$PAR_umolm2s[msn2_ctd$Depth_m>2]),
                                   mean(msn3_ctd$PAR_umolm2s[msn3_ctd$Depth_m>2]),mean(msn4_ctd$PAR_umolm2s[msn4_ctd$Depth_m>2]),
                                   after_msn5_ctd$PAR_umolm2s[after_msn5_ctd$Depth_m>2])


euclidean_drivers_df$thermo_depth <- c(mean(c(first(msn1_therm_depth$therm_depth),last(msn1_therm_depth$therm_depth))),
                                       mean(c(first(msn2_therm_depth$therm_depth),last(msn2_therm_depth$therm_depth))),
                                       mean(c(first(msn3_therm_depth$therm_depth),last(msn3_therm_depth$therm_depth))),
                                       mean(c(first(msn4_therm_depth$therm_depth),last(msn4_therm_depth$therm_depth))),
                                       mean(c(first(msn5_therm_depth$therm_depth),last(msn5_therm_depth$therm_depth))))

euclidean_drivers_df$TN_surf <- c(msn1_chem$TN_ugL[msn1_chem$Depth_m==0.1], 
                                  msn2_chem$TN_ugL[msn2_chem$Depth_m==0.1],
                                  msn3_chem$TN_ugL[msn3_chem$Depth_m==0.1],
                                  msn4_chem$TN_ugL[msn4_chem$Depth_m==0.1],
                                  msn5_chem$TN_ugL[msn5_chem$Depth_m==0.1]) 

euclidean_drivers_df$TP_surf <- c(msn1_chem$TP_ugL[msn1_chem$Depth_m==0.1], 
                                  msn2_chem$TP_ugL[msn2_chem$Depth_m==0.1],
                                  msn3_chem$TP_ugL[msn3_chem$Depth_m==0.1],
                                  msn4_chem$TP_ugL[msn4_chem$Depth_m==0.1],
                                  msn5_chem$TP_ugL[msn5_chem$Depth_m==0.1]) 

euclidean_drivers_df$NH4_surf <- c(msn1_chem$NH4_ugL[msn1_chem$Depth_m==0.1], 
                                   msn2_chem$NH4_ugL[msn2_chem$Depth_m==0.1],
                                   msn3_chem$NH4_ugL[msn3_chem$Depth_m==0.1],
                                   msn4_chem$NH4_ugL[msn4_chem$Depth_m==0.1],
                                   msn5_chem$NH4_ugL[msn5_chem$Depth_m==0.1]) 

euclidean_drivers_df$NO3NO2_surf <- c(msn1_chem$NO3NO2_ugL[msn1_chem$Depth_m==0.1], 
                                      msn2_chem$NO3NO2_ugL[msn2_chem$Depth_m==0.1],
                                      msn3_chem$NO3NO2_ugL[msn3_chem$Depth_m==0.1],
                                      msn4_chem$NO3NO2_ugL[msn4_chem$Depth_m==0.1],
                                      msn5_chem$NO3NO2_ugL[msn5_chem$Depth_m==0.1]) 

euclidean_drivers_df$SRP_surf <- c(msn1_chem$SRP_ugL[msn1_chem$Depth_m==0.1], 
                                   msn2_chem$SRP_ugL[msn2_chem$Depth_m==0.1],
                                   msn3_chem$SRP_ugL[msn3_chem$Depth_m==0.1],
                                   msn4_chem$SRP_ugL[msn4_chem$Depth_m==0.1],
                                   msn5_chem$SRP_ugL[msn5_chem$Depth_m==0.1]) 

euclidean_drivers_df$DOC_surf <- c(msn1_chem$DOC_mgL[msn1_chem$Depth_m==0.1], 
                                   msn2_chem$DOC_mgL[msn2_chem$Depth_m==0.1],
                                   msn3_chem$DOC_mgL[msn3_chem$Depth_m==0.1],
                                   msn4_chem$DOC_mgL[msn4_chem$Depth_m==0.1],
                                   msn5_chem$DOC_mgL[msn5_chem$Depth_m==0.1]) 


#calculate avg density col for each day at both sites
euclidean_drivers_df$pel_avg_dens <- c(mean(zoop_epi_tows$ZoopDensity_No.pL_1[zoop_epi_tows$site=="pel" & zoop_epi_tows$groups==1]),
                                       mean(zoop_epi_tows$ZoopDensity_No.pL_1[zoop_epi_tows$site=="pel" & zoop_epi_tows$groups==2]),
                                       mean(zoop_epi_tows$ZoopDensity_No.pL_1[zoop_epi_tows$site=="pel" & zoop_epi_tows$groups==3]),
                                       mean(zoop_epi_tows$ZoopDensity_No.pL_1[zoop_epi_tows$site=="pel" & zoop_epi_tows$groups==4]),
                                       mean(zoop_epi_tows$ZoopDensity_No.pL_1[zoop_epi_tows$site=="pel" & zoop_epi_tows$groups==5]))


euclidean_drivers_df$lit_avg_dens <- c(mean(zoop_epi_tows$ZoopDensity_No.pL_1[zoop_epi_tows$site=="lit" & zoop_epi_tows$groups==1]),
                                       mean(zoop_epi_tows$ZoopDensity_No.pL_1[zoop_epi_tows$site=="lit" & zoop_epi_tows$groups==2]),
                                       mean(zoop_epi_tows$ZoopDensity_No.pL_1[zoop_epi_tows$site=="lit" & zoop_epi_tows$groups==3]),
                                       mean(zoop_epi_tows$ZoopDensity_No.pL_1[zoop_epi_tows$site=="lit" & zoop_epi_tows$groups==4]),
                                       mean(zoop_epi_tows$ZoopDensity_No.pL_1[zoop_epi_tows$site=="lit" & zoop_epi_tows$groups==5]))

euclidean_drivers_df$oxycline_depth <- c(6.5,5.5,3.5,6.5,5) #see calcs below

#read in DVM and DHM metrics
migration_metrics <- read.csv(file.path(getwd(),"Summer2021-DataAnalysis/SummaryStats/migration_metrics.csv"), header = TRUE)

#long to wide
migration_metrics_wide <-  migration_metrics %>% pivot_wider(
  id_cols = "metric", id_expand = TRUE,
  names_from = c("migration","MSN"), values_from = value)

euclidean_drivers_df$Calanoida_dens_DVM_metric <-c(migration_metrics_wide$DVM_avg_1[migration_metrics_wide$metric=="Calanoida_density_NopL"],
                                                   migration_metrics_wide$DVM_avg_2[migration_metrics_wide$metric=="Calanoida_density_NopL"],
                                                   migration_metrics_wide$DVM_avg_3[migration_metrics_wide$metric=="Calanoida_density_NopL"],
                                                   migration_metrics_wide$DVM_avg_4[migration_metrics_wide$metric=="Calanoida_density_NopL"],
                                                   migration_metrics_wide$DVM_avg_5[migration_metrics_wide$metric=="Calanoida_density_NopL"])

euclidean_drivers_df$Calanoida_biom_DVM_metric <-c(migration_metrics_wide$DVM_avg_1[migration_metrics_wide$metric=="Calanoida_BiomassConcentration_ugpL"],
                                                   migration_metrics_wide$DVM_avg_2[migration_metrics_wide$metric=="Calanoida_BiomassConcentration_ugpL"],
                                                   migration_metrics_wide$DVM_avg_3[migration_metrics_wide$metric=="Calanoida_BiomassConcentration_ugpL"],
                                                   migration_metrics_wide$DVM_avg_4[migration_metrics_wide$metric=="Calanoida_BiomassConcentration_ugpL"],
                                                   migration_metrics_wide$DVM_avg_5[migration_metrics_wide$metric=="Calanoida_BiomassConcentration_ugpL"])

euclidean_drivers_df$Cladocera_dens_DVM_metric <-c(migration_metrics_wide$DVM_avg_1[migration_metrics_wide$metric=="Cladocera_density_NopL"],
                                                   migration_metrics_wide$DVM_avg_2[migration_metrics_wide$metric=="Cladocera_density_NopL"],
                                                   migration_metrics_wide$DVM_avg_3[migration_metrics_wide$metric=="Cladocera_density_NopL"],
                                                   migration_metrics_wide$DVM_avg_4[migration_metrics_wide$metric=="Cladocera_density_NopL"],
                                                   migration_metrics_wide$DVM_avg_5[migration_metrics_wide$metric=="Cladocera_density_NopL"])

euclidean_drivers_df$Cladocera_biom_DVM_metric <-c(migration_metrics_wide$DVM_avg_1[migration_metrics_wide$metric=="Cladocera_BiomassConcentration_ugpL"],
                                                   migration_metrics_wide$DVM_avg_2[migration_metrics_wide$metric=="Cladocera_BiomassConcentration_ugpL"],
                                                   migration_metrics_wide$DVM_avg_3[migration_metrics_wide$metric=="Cladocera_BiomassConcentration_ugpL"],
                                                   migration_metrics_wide$DVM_avg_4[migration_metrics_wide$metric=="Cladocera_BiomassConcentration_ugpL"],
                                                   migration_metrics_wide$DVM_avg_5[migration_metrics_wide$metric=="Cladocera_BiomassConcentration_ugpL"])

euclidean_drivers_df$Copepoda_dens_DVM_metric <-c(migration_metrics_wide$DVM_avg_1[migration_metrics_wide$metric=="Copepoda_density_NopL"],
                                                  migration_metrics_wide$DVM_avg_2[migration_metrics_wide$metric=="Copepoda_density_NopL"],
                                                  migration_metrics_wide$DVM_avg_3[migration_metrics_wide$metric=="Copepoda_density_NopL"],
                                                  migration_metrics_wide$DVM_avg_4[migration_metrics_wide$metric=="Copepoda_density_NopL"],
                                                  migration_metrics_wide$DVM_avg_5[migration_metrics_wide$metric=="Copepoda_density_NopL"])

euclidean_drivers_df$Copepoda_biom_DVM_metric <-c(migration_metrics_wide$DVM_avg_1[migration_metrics_wide$metric=="Copepoda_BiomassConcentration_ugpL"],
                                                  migration_metrics_wide$DVM_avg_2[migration_metrics_wide$metric=="Copepoda_BiomassConcentration_ugpL"],
                                                  migration_metrics_wide$DVM_avg_3[migration_metrics_wide$metric=="Copepoda_BiomassConcentration_ugpL"],
                                                  migration_metrics_wide$DVM_avg_4[migration_metrics_wide$metric=="Copepoda_BiomassConcentration_ugpL"],
                                                  migration_metrics_wide$DVM_avg_5[migration_metrics_wide$metric=="Copepoda_BiomassConcentration_ugpL"])

euclidean_drivers_df$Cyclopoida_dens_DVM_metric <-c(migration_metrics_wide$DVM_avg_1[migration_metrics_wide$metric=="Cyclopoida_density_NopL"],
                                                    migration_metrics_wide$DVM_avg_2[migration_metrics_wide$metric=="Cyclopoida_density_NopL"],
                                                    migration_metrics_wide$DVM_avg_3[migration_metrics_wide$metric=="Cyclopoida_density_NopL"],
                                                    migration_metrics_wide$DVM_avg_4[migration_metrics_wide$metric=="Cyclopoida_density_NopL"],
                                                    migration_metrics_wide$DVM_avg_5[migration_metrics_wide$metric=="Cyclopoida_density_NopL"])

euclidean_drivers_df$Cyclopoida_biom_DVM_metric <-c(migration_metrics_wide$DVM_avg_1[migration_metrics_wide$metric=="Cyclopoida_BiomassConcentration_ugpL"],
                                                    migration_metrics_wide$DVM_avg_2[migration_metrics_wide$metric=="Cyclopoida_BiomassConcentration_ugpL"],
                                                    migration_metrics_wide$DVM_avg_3[migration_metrics_wide$metric=="Cyclopoida_BiomassConcentration_ugpL"],
                                                    migration_metrics_wide$DVM_avg_4[migration_metrics_wide$metric=="Cyclopoida_BiomassConcentration_ugpL"],
                                                    migration_metrics_wide$DVM_avg_5[migration_metrics_wide$metric=="Cyclopoida_BiomassConcentration_ugpL"])

euclidean_drivers_df$Rotifera_dens_DVM_metric <-c(migration_metrics_wide$DVM_avg_1[migration_metrics_wide$metric=="Rotifera_density_NopL"],
                                                  migration_metrics_wide$DVM_avg_2[migration_metrics_wide$metric=="Rotifera_density_NopL"],
                                                  migration_metrics_wide$DVM_avg_3[migration_metrics_wide$metric=="Rotifera_density_NopL"],
                                                  migration_metrics_wide$DVM_avg_4[migration_metrics_wide$metric=="Rotifera_density_NopL"],
                                                  migration_metrics_wide$DVM_avg_5[migration_metrics_wide$metric=="Rotifera_density_NopL"])

euclidean_drivers_df$Rotifera_biom_DVM_metric <-c(migration_metrics_wide$DVM_avg_1[migration_metrics_wide$metric=="Rotifera_BiomassConcentration_ugpL"],
                                                  migration_metrics_wide$DVM_avg_2[migration_metrics_wide$metric=="Rotifera_BiomassConcentration_ugpL"],
                                                  migration_metrics_wide$DVM_avg_3[migration_metrics_wide$metric=="Rotifera_BiomassConcentration_ugpL"],
                                                  migration_metrics_wide$DVM_avg_4[migration_metrics_wide$metric=="Rotifera_BiomassConcentration_ugpL"],
                                                  migration_metrics_wide$DVM_avg_5[migration_metrics_wide$metric=="Rotifera_BiomassConcentration_ugpL"])

euclidean_drivers_df$Calanoida_dens_DHM_metric <-c(migration_metrics_wide$DHM_avg_1[migration_metrics_wide$metric=="Calanoida_density_NopL"],
                                                   migration_metrics_wide$DHM_avg_2[migration_metrics_wide$metric=="Calanoida_density_NopL"],
                                                   migration_metrics_wide$DHM_avg_3[migration_metrics_wide$metric=="Calanoida_density_NopL"],
                                                   migration_metrics_wide$DHM_avg_4[migration_metrics_wide$metric=="Calanoida_density_NopL"],
                                                   migration_metrics_wide$DHM_avg_5[migration_metrics_wide$metric=="Calanoida_density_NopL"])

euclidean_drivers_df$Calanoida_biom_DHM_metric <-c(migration_metrics_wide$DHM_avg_1[migration_metrics_wide$metric=="Calanoida_BiomassConcentration_ugpL"],
                                                   migration_metrics_wide$DHM_avg_2[migration_metrics_wide$metric=="Calanoida_BiomassConcentration_ugpL"],
                                                   migration_metrics_wide$DHM_avg_3[migration_metrics_wide$metric=="Calanoida_BiomassConcentration_ugpL"],
                                                   migration_metrics_wide$DHM_avg_4[migration_metrics_wide$metric=="Calanoida_BiomassConcentration_ugpL"],
                                                   migration_metrics_wide$DHM_avg_5[migration_metrics_wide$metric=="Calanoida_BiomassConcentration_ugpL"])

euclidean_drivers_df$Cladocera_dens_DHM_metric <-c(migration_metrics_wide$DHM_avg_1[migration_metrics_wide$metric=="Cladocera_density_NopL"],
                                                   migration_metrics_wide$DHM_avg_2[migration_metrics_wide$metric=="Cladocera_density_NopL"],
                                                   migration_metrics_wide$DHM_avg_3[migration_metrics_wide$metric=="Cladocera_density_NopL"],
                                                   migration_metrics_wide$DHM_avg_4[migration_metrics_wide$metric=="Cladocera_density_NopL"],
                                                   migration_metrics_wide$DHM_avg_5[migration_metrics_wide$metric=="Cladocera_density_NopL"])

euclidean_drivers_df$Cladocera_biom_DHM_metric <-c(migration_metrics_wide$DHM_avg_1[migration_metrics_wide$metric=="Cladocera_BiomassConcentration_ugpL"],
                                                   migration_metrics_wide$DHM_avg_2[migration_metrics_wide$metric=="Cladocera_BiomassConcentration_ugpL"],
                                                   migration_metrics_wide$DHM_avg_3[migration_metrics_wide$metric=="Cladocera_BiomassConcentration_ugpL"],
                                                   migration_metrics_wide$DHM_avg_4[migration_metrics_wide$metric=="Cladocera_BiomassConcentration_ugpL"],
                                                   migration_metrics_wide$DHM_avg_5[migration_metrics_wide$metric=="Cladocera_BiomassConcentration_ugpL"])

euclidean_drivers_df$Copepoda_dens_DHM_metric <-c(migration_metrics_wide$DHM_avg_1[migration_metrics_wide$metric=="Copepoda_density_NopL"],
                                                  migration_metrics_wide$DHM_avg_2[migration_metrics_wide$metric=="Copepoda_density_NopL"],
                                                  migration_metrics_wide$DHM_avg_3[migration_metrics_wide$metric=="Copepoda_density_NopL"],
                                                  migration_metrics_wide$DHM_avg_4[migration_metrics_wide$metric=="Copepoda_density_NopL"],
                                                  migration_metrics_wide$DHM_avg_5[migration_metrics_wide$metric=="Copepoda_density_NopL"])

euclidean_drivers_df$Copepoda_biom_DHM_metric <-c(migration_metrics_wide$DHM_avg_1[migration_metrics_wide$metric=="Copepoda_BiomassConcentration_ugpL"],
                                                  migration_metrics_wide$DHM_avg_2[migration_metrics_wide$metric=="Copepoda_BiomassConcentration_ugpL"],
                                                  migration_metrics_wide$DHM_avg_3[migration_metrics_wide$metric=="Copepoda_BiomassConcentration_ugpL"],
                                                  migration_metrics_wide$DHM_avg_4[migration_metrics_wide$metric=="Copepoda_BiomassConcentration_ugpL"],
                                                  migration_metrics_wide$DHM_avg_5[migration_metrics_wide$metric=="Copepoda_BiomassConcentration_ugpL"])

euclidean_drivers_df$Cyclopoida_dens_DHM_metric <-c(migration_metrics_wide$DHM_avg_1[migration_metrics_wide$metric=="Cyclopoida_density_NopL"],
                                                    migration_metrics_wide$DHM_avg_2[migration_metrics_wide$metric=="Cyclopoida_density_NopL"],
                                                    migration_metrics_wide$DHM_avg_3[migration_metrics_wide$metric=="Cyclopoida_density_NopL"],
                                                    migration_metrics_wide$DHM_avg_4[migration_metrics_wide$metric=="Cyclopoida_density_NopL"],
                                                    migration_metrics_wide$DHM_avg_5[migration_metrics_wide$metric=="Cyclopoida_density_NopL"])

euclidean_drivers_df$Cyclopoida_biom_DHM_metric <-c(migration_metrics_wide$DHM_avg_1[migration_metrics_wide$metric=="Cyclopoida_BiomassConcentration_ugpL"],
                                                    migration_metrics_wide$DHM_avg_2[migration_metrics_wide$metric=="Cyclopoida_BiomassConcentration_ugpL"],
                                                    migration_metrics_wide$DHM_avg_3[migration_metrics_wide$metric=="Cyclopoida_BiomassConcentration_ugpL"],
                                                    migration_metrics_wide$DHM_avg_4[migration_metrics_wide$metric=="Cyclopoida_BiomassConcentration_ugpL"],
                                                    migration_metrics_wide$DHM_avg_5[migration_metrics_wide$metric=="Cyclopoida_BiomassConcentration_ugpL"])

euclidean_drivers_df$Rotifera_dens_DHM_metric <-c(migration_metrics_wide$DHM_avg_1[migration_metrics_wide$metric=="Rotifera_density_NopL"],
                                                  migration_metrics_wide$DHM_avg_2[migration_metrics_wide$metric=="Rotifera_density_NopL"],
                                                  migration_metrics_wide$DHM_avg_3[migration_metrics_wide$metric=="Rotifera_density_NopL"],
                                                  migration_metrics_wide$DHM_avg_4[migration_metrics_wide$metric=="Rotifera_density_NopL"],
                                                  migration_metrics_wide$DHM_avg_5[migration_metrics_wide$metric=="Rotifera_density_NopL"])

euclidean_drivers_df$Rotifera_biom_DHM_metric <-c(migration_metrics_wide$DHM_avg_1[migration_metrics_wide$metric=="Rotifera_BiomassConcentration_ugpL"],
                                                  migration_metrics_wide$DHM_avg_2[migration_metrics_wide$metric=="Rotifera_BiomassConcentration_ugpL"],
                                                  migration_metrics_wide$DHM_avg_3[migration_metrics_wide$metric=="Rotifera_BiomassConcentration_ugpL"],
                                                  migration_metrics_wide$DHM_avg_4[migration_metrics_wide$metric=="Rotifera_BiomassConcentration_ugpL"],
                                                  migration_metrics_wide$DHM_avg_5[migration_metrics_wide$metric=="Rotifera_BiomassConcentration_ugpL"])


#----------------------------------------------------------------------------------------#
#Now calculate a correlation matrix for euclidean_drivers_df (and just select row 1 for the pel distance vs driver )
distVsdriver_cor <- cor(euclidean_drivers_df[2:27], method="spearman")[c(1,2),]

#plot response variable (euclidean distances) against environmental/biological data
#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_epiDO.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$DO_0.1m,euclidean_drivers_df$pel_euc_dist, xlab="Epi DO", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5, ylim=c(0,3))
points(euclidean_drivers_df$DO_0.1m[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$DO_0.1m[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$DO_0.1m[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$DO_0.1m[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$DO_0.1m[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_hypoDO.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$DO_9.0m,euclidean_drivers_df$pel_euc_dist, xlab="Hypo DO", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$DO_9.0m[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$DO_9.0m[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$DO_9.0m[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$DO_9.0m[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$DO_9.0m[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_epitemp.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$temp_0.1m,euclidean_drivers_df$pel_euc_dist, xlab="Epi Temp", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$temp_0.1m[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$temp_0.1m[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$temp_0.1m[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$temp_0.1m[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$temp_0.1m[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_hypotemp.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$temp_9.0m,euclidean_drivers_df$pel_euc_dist, xlab="Hypo Temp", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$temp_9.0m[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$temp_9.0m[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$temp_9.0m[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$temp_9.0m[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$temp_9.0m[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_epispcond.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$Spcond_0.1m,euclidean_drivers_df$pel_euc_dist, xlab="Epi Sp Cond", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$Spcond_0.1m[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$Spcond_0.1m[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$Spcond_0.1m[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$Spcond_0.1m[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$Spcond_0.1m[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_hypospcond.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$Spcond_9.0m,euclidean_drivers_df$pel_euc_dist, xlab="Hypo Sp Cond", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$Spcond_9.0m[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$Spcond_9.0m[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$Spcond_9.0m[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$Spcond_9.0m[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$Spcond_9.0m[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_epichl.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$chla_0.1m,euclidean_drivers_df$pel_euc_dist, xlab="Epi Chl a", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$chla_0.1m[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$chla_0.1m[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$chla_0.1m[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$chla_0.1m[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$chla_0.1m[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_hypochl.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$chla_9.0m,euclidean_drivers_df$pel_euc_dist, xlab="Hypo Chl a", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$chla_9.0m[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$chla_9.0m[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$chla_9.0m[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$chla_9.0m[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$chla_9.0m[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_epiturb.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$turb_0.1m,euclidean_drivers_df$pel_euc_dist, xlab="Epi Turb", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$turb_0.1m[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$turb_0.1m[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$turb_0.1m[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$turb_0.1m[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$turb_0.1m[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_hypoturb.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$turb_9.0m,euclidean_drivers_df$pel_euc_dist, xlab="Hypo Turb", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$turb_9.0m[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=4)
points(euclidean_drivers_df$turb_9.0m[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=4)
points(euclidean_drivers_df$turb_9.0m[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=4)
points(euclidean_drivers_df$turb_9.0m[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=4)
points(euclidean_drivers_df$turb_9.0m[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=4)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_epipar.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$par_0.1m,euclidean_drivers_df$pel_euc_dist, xlab="Epi PAR", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$par_0.1m[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$par_0.1m[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$par_0.1m[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$par_0.1m[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$par_0.1m[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_hypopar.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$par_9.0m,euclidean_drivers_df$pel_euc_dist, xlab="Hypo PAR", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$par_9.0m[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=4)
points(euclidean_drivers_df$par_9.0m[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=4)
points(euclidean_drivers_df$par_9.0m[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=4)
points(euclidean_drivers_df$par_9.0m[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=4)
points(euclidean_drivers_df$par_9.0m[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=4)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_pelzoopdens.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$pel_avg_dens,euclidean_drivers_df$pel_euc_dist, xlab="Pel Zoop Dens", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$pel_avg_dens[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$pel_avg_dens[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$pel_avg_dens[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$pel_avg_dens[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$pel_avg_dens[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_litzoopdens.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$lit_avg_dens,euclidean_drivers_df$pel_euc_dist, xlab="Lit Zoop Dens", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$lit_avg_dens[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$lit_avg_dens[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$lit_avg_dens[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$lit_avg_dens[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$lit_avg_dens[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_thermdepth.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$thermo_depth,euclidean_drivers_df$pel_euc_dist, xlab="Thermocline Depth", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$thermo_depth[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$thermo_depth[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$thermo_depth[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$thermo_depth[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$thermo_depth[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_TN.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$TN_surf,euclidean_drivers_df$pel_euc_dist, xlab="TN (ug/L)", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$TN_surf[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$TN_surf[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$TN_surf[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$TN_surf[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$TN_surf[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_TP.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$TP_surf,euclidean_drivers_df$pel_euc_dist, xlab="TP (ug/L)", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$TP_surf[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$TP_surf[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$TP_surf[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$TP_surf[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$TP_surf[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_NH4.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$NH4_surf,euclidean_drivers_df$pel_euc_dist, xlab="NH4 (ug/L)", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$NH4_surf[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$NH4_surf[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$NH4_surf[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$NH4_surf[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$NH4_surf[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_NO3NO2.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$NO3NO2_surf,euclidean_drivers_df$pel_euc_dist, xlab="NO3NO2 (ug/L)", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$NO3NO2_surf[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$NO3NO2_surf[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$NO3NO2_surf[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$NO3NO2_surf[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$NO3NO2_surf[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_SRP.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$SRP_surf,euclidean_drivers_df$pel_euc_dist, xlab="SRP (ug/L)", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$SRP_surf[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$SRP_surf[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$SRP_surf[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$SRP_surf[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$SRP_surf[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_DOC.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$DOC_surf,euclidean_drivers_df$pel_euc_dist, xlab="DOC (mg/L)", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$DOC_surf[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$DOC_surf[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$DOC_surf[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$DOC_surf[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$DOC_surf[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_oxydepth.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$oxycline_depth,euclidean_drivers_df$pel_euc_dist, xlab="Oxycline Depth", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$oxycline_depth[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$oxycline_depth[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$oxycline_depth[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$oxycline_depth[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$oxycline_depth[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#-------------------------------------------------------------------------------
# Euclidean distance vs migration metrics

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_calanodiaDVMdens.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$Calanoida_dens_DVM_metric,euclidean_drivers_df$pel_euc_dist, xlab="Calanoida DVM density", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$Calanoida_dens_DVM_metric[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$Calanoida_dens_DVM_metric[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$Calanoida_dens_DVM_metric[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$Calanoida_dens_DVM_metric[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$Calanoida_dens_DVM_metric[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_calanodiaDVMbiom.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$Calanoida_biom_DVM_metric,euclidean_drivers_df$pel_euc_dist, xlab="Calanoida DVM biomass", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$Calanoida_biom_DVM_metric[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$Calanoida_biom_DVM_metric[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$Calanoida_biom_DVM_metric[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$Calanoida_biom_DVM_metric[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$Calanoida_biom_DVM_metric[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_cladoceraDVMdens.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$Cladocera_dens_DVM_metric,euclidean_drivers_df$pel_euc_dist, xlab="Cladocera DVM density", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$Cladocera_dens_DVM_metric[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$Cladocera_dens_DVM_metric[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$Cladocera_dens_DVM_metric[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$Cladocera_dens_DVM_metric[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$Cladocera_dens_DVM_metric[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_cladoceraDVMbiom.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$Cladocera_biom_DVM_metric,euclidean_drivers_df$pel_euc_dist, xlab="Cladocera DVM biomass", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$Cladocera_biom_DVM_metric[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$Cladocera_biom_DVM_metric[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$Cladocera_biom_DVM_metric[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$Cladocera_biom_DVM_metric[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$Cladocera_biom_DVM_metric[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_cyclopoidaDVMdens.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$Cyclopoida_dens_DVM_metric,euclidean_drivers_df$pel_euc_dist, xlab="Cyclopoida DVM density", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$Cyclopoida_dens_DVM_metric[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$Cyclopoida_dens_DVM_metric[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$Cyclopoida_dens_DVM_metric[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$Cyclopoida_dens_DVM_metric[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$Cyclopoida_dens_DVM_metric[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_cyclopoidaDVMbiom.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$Cyclopoida_biom_DVM_metric,euclidean_drivers_df$pel_euc_dist, xlab="Cyclopoida DVM biomass", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$Cyclopoida_biom_DVM_metric[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$Cyclopoida_biom_DVM_metric[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$Cyclopoida_biom_DVM_metric[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$Cyclopoida_biom_DVM_metric[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$Cyclopoida_biom_DVM_metric[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_rotiferaDVMdens.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$Rotifera_dens_DVM_metric,euclidean_drivers_df$pel_euc_dist, xlab="Rotifera DVM density", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$Rotifera_dens_DVM_metric[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$Rotifera_dens_DVM_metric[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$Rotifera_dens_DVM_metric[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$Rotifera_dens_DVM_metric[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$Rotifera_dens_DVM_metric[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_rotiferaDVMbiom.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$Rotifera_biom_DVM_metric,euclidean_drivers_df$pel_euc_dist, xlab="Rotifera DVM biomass", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$Rotifera_biom_DVM_metric[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$Rotifera_biom_DVM_metric[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$Rotifera_biom_DVM_metric[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$Rotifera_biom_DVM_metric[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$Rotifera_biom_DVM_metric[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_copepodaDVMdens.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$Copepoda_dens_DVM_metric,euclidean_drivers_df$pel_euc_dist, xlab="Copepoda DVM density", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$Copepoda_dens_DVM_metric[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$Copepoda_dens_DVM_metric[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$Copepoda_dens_DVM_metric[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$Copepoda_dens_DVM_metric[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$Copepoda_dens_DVM_metric[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_copepodaDVMbiom.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$Copepoda_biom_DVM_metric,euclidean_drivers_df$pel_euc_dist, xlab="Copepoda DVM biomass", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$Copepoda_biom_DVM_metric[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$Copepoda_biom_DVM_metric[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$Copepoda_biom_DVM_metric[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$Copepoda_biom_DVM_metric[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$Copepoda_biom_DVM_metric[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_calanodiaDHMdens.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$Calanoida_dens_DHM_metric,euclidean_drivers_df$pel_euc_dist, xlab="Calanoida DHM density", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$Calanoida_dens_DHM_metric[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$Calanoida_dens_DHM_metric[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$Calanoida_dens_DHM_metric[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$Calanoida_dens_DHM_metric[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$Calanoida_dens_DHM_metric[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_calanodiaDHMbiom.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$Calanoida_biom_DHM_metric,euclidean_drivers_df$pel_euc_dist, xlab="Calanoida DHM biomass", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$Calanoida_biom_DHM_metric[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$Calanoida_biom_DHM_metric[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$Calanoida_biom_DHM_metric[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$Calanoida_biom_DHM_metric[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$Calanoida_biom_DHM_metric[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_cladoceraDHMdens.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$Cladocera_dens_DHM_metric,euclidean_drivers_df$pel_euc_dist, xlab="Cladocera DHM density", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$Cladocera_dens_DHM_metric[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$Cladocera_dens_DHM_metric[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$Cladocera_dens_DHM_metric[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$Cladocera_dens_DHM_metric[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$Cladocera_dens_DHM_metric[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_cladoceraDHMbiom.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$Cladocera_biom_DHM_metric,euclidean_drivers_df$pel_euc_dist, xlab="Cladocera DHM biomass", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$Cladocera_biom_DHM_metric[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$Cladocera_biom_DHM_metric[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$Cladocera_biom_DHM_metric[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$Cladocera_biom_DHM_metric[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$Cladocera_biom_DHM_metric[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_cyclopoidaDHMdens.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$Cyclopoida_dens_DHM_metric,euclidean_drivers_df$pel_euc_dist, xlab="Cyclopoida DHM density", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$Cyclopoida_dens_DHM_metric[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$Cyclopoida_dens_DHM_metric[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$Cyclopoida_dens_DHM_metric[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$Cyclopoida_dens_DHM_metric[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$Cyclopoida_dens_DHM_metric[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_cyclopoidaDHMbiom.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$Cyclopoida_biom_DHM_metric,euclidean_drivers_df$pel_euc_dist, xlab="Cyclopoida DHM biomass", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$Cyclopoida_biom_DHM_metric[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$Cyclopoida_biom_DHM_metric[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$Cyclopoida_biom_DHM_metric[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$Cyclopoida_biom_DHM_metric[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$Cyclopoida_biom_DHM_metric[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_rotiferaDHMdens.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$Rotifera_dens_DHM_metric,euclidean_drivers_df$pel_euc_dist, xlab="Rotifera DHM density", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$Rotifera_dens_DHM_metric[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$Rotifera_dens_DHM_metric[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$Rotifera_dens_DHM_metric[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$Rotifera_dens_DHM_metric[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$Rotifera_dens_DHM_metric[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_rotiferaDHMbiom.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$Rotifera_biom_DHM_metric,euclidean_drivers_df$pel_euc_dist, xlab="Rotifera DHM biomass", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$Rotifera_biom_DHM_metric[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$Rotifera_biom_DHM_metric[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$Rotifera_biom_DHM_metric[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$Rotifera_biom_DHM_metric[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$Rotifera_biom_DHM_metric[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_copepodaDHMdens.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$Copepoda_dens_DHM_metric,euclidean_drivers_df$pel_euc_dist, xlab="Copepoda DHM density", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$Copepoda_dens_DHM_metric[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$Copepoda_dens_DHM_metric[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$Copepoda_dens_DHM_metric[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$Copepoda_dens_DHM_metric[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$Copepoda_dens_DHM_metric[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()

#jpeg(file.path(getwd(), "Summer2021-DataAnalysis/Figures/2019-2021_pel_euclidean_dist_daily_sums_vs_copepodaDHMbiom.jpg"), width = 6, height = 5, units = "in",res = 300)
plot(euclidean_drivers_df$Copepoda_biom_DHM_metric,euclidean_drivers_df$pel_euc_dist, xlab="Copepoda DHM biomass", ylab="distance", cex=2.8, pch=21, cex.lab = 1.5)
points(euclidean_drivers_df$Copepoda_biom_DHM_metric[1],euclidean_drivers_df$pel_euc_dist[1], bg="#008585",pch=21,cex=3)
points(euclidean_drivers_df$Copepoda_biom_DHM_metric[2],euclidean_drivers_df$pel_euc_dist[2], bg="#9BBAA0",pch=21,cex=3)
points(euclidean_drivers_df$Copepoda_biom_DHM_metric[3],euclidean_drivers_df$pel_euc_dist[3], bg="#F2E2B0",pch=21,cex=3)
points(euclidean_drivers_df$Copepoda_biom_DHM_metric[4],euclidean_drivers_df$pel_euc_dist[4], bg="#DEA868",pch=21,cex=3)
points(euclidean_drivers_df$Copepoda_biom_DHM_metric[5],euclidean_drivers_df$pel_euc_dist[5], bg="#C7522B",pch=21,cex=3)
legend("bottomright",legend=c("day1","day2","day3","day4","day5"),pch=21, pt.bg=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"),bty = "n",cex=1.4)
#dev.off()


#ctd 
ctd_clean <- ctd %>%
  dplyr::mutate(rdepth = plyr::round_any(Depth_m, 0.5)) %>% 
  dplyr::filter(!is.na(DO_mgL)) %>%
  dplyr::group_by(Reservoir,Date,Site, rdepth) %>%
  dplyr::summarise(value = mean(DO_mgL)) %>% 
  dplyr::rename(depth = rdepth)  %>%
  dplyr::group_by(Reservoir,Site, Date) %>% 
  dplyr::mutate(oxy = min(depth[value<=2], na.rm=TRUE))

depths <- ctd_clean$oxy[(ctd_clean$Date=="2019-07-10" | ctd_clean$Date=="2019-07-24" |
                           ctd_clean$Date=="2020-08-12" | ctd_clean$Date=="2021-06-16" | ctd_clean$Date=="2021-07-12") & ctd_clean$Reservoir=="BVR" & ctd_clean$Site==50 & ctd_clean$depth>0]


#oxygen plots for all 5 MSNs
ggplot(subset(ctd_clean, depth > 0 & Reservoir=="BVR" & Site==50 & Date %in% c(as.Date("2019-07-10"), as.Date("2019-07-24"),
                                                                               as.Date("2020-08-12"), as.Date("2021-06-16"), as.Date("2021-07-12"))), aes(value,depth,color=as.factor(Date))) + 
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin= depths, ymax=Inf), fill="red",alpha=0.03) +
  scale_color_manual("",values=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"), guide="none") +  xlab("DO (mg/L)") + ylab("Depth (m)") +
  ylim(10,0) + geom_point() + geom_path() + theme_bw() + facet_wrap(~Date, ncol=5) +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = c(0.76,0.02),
        legend.background = element_blank(),legend.direction = "horizontal", panel.grid.minor = element_blank(), legend.key=element_rect(fill=NA),
        plot.margin = unit(c(0,0.05,0,0), "cm"),legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),
        legend.title = element_text(size = 4.5),legend.text  = element_text(size = 4.5), panel.spacing=unit(0, "cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6), axis.text.y = element_text(size=6)) 
#ggsave(file.path(getwd(),"Summer2021-DataAnalysis/Figures/DO_profiles.jpg"), width=4, height=3) 

