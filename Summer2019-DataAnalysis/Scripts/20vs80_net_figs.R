# Figs comparing 20 vs 80 um zoop net data
# 22Jun2023

#read in 2019 data
zoops2019<- read.csv(file.path(getwd(),'Summer2019-DataAnalysis/SummaryStats/FCR_ZooplanktonSummary2019.csv'),header = TRUE)

#pull data for days where we collected both 20um and 80um net samples
