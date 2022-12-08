#Script to calculate DVM and DHM metrics
#Created 5Dec2022

#read in libraries
pacman::p_load(dplyr, vegan, labdsv, goeveg, rLakeAnalyzer, ggplot2,tidyr,lubridate, scales, colorblindcheck)

#cyclopoid, rotifer,calanoid, copepod density (#/L)
boxplot(all_zoops[,c(12,16,21,25)])

#read in the DVM and DHM csvs (note: these are just density and % of total)
all_DHM <- read.csv(paste0(getwd(),"/Summer2021-DataAnalysis/SummaryStats/All_MSN_DHM.csv"))
all_DVM <- read.csv()

#DVM metrics for a single day --> DVM = (Depi / Depi + Dhypo)Night - (Depi / Depi + Dhypo)Day


#DHM metrics for a single day --> DHM = (Dpelepi / Dpelepi + Dlit)Night - (Dpelepi / Dpelepi + Dlit)Day



