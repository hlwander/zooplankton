#script for getting all summary files into EDI format

#read in libraries
pacman::p_load(tidyr, dplyr)

#function to count characters starting at the end of the string
substrEnd <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#read in zoop data from all 3 years
zoops2019<- read.csv(file.path(getwd(),'Summer2019-DataAnalysis/SummaryStats/FCR_ZooplanktonSummary2019.csv'),header = TRUE)
zoops2020<- read.csv(file.path(getwd(),'Summer2020-DataAnalysis/SummaryStats/FCR_ZooplanktonSummary2020.csv'),header = TRUE)
zoops2021<- read.csv(file.path(getwd(),'Summer2021-DataAnalysis/SummaryStats/FCR_ZooplanktonSummary2021.csv'),header = TRUE)

#drop holopedium from 2019 so can merge files below
zoops2019 <- zoops2019[,!grepl("Holopedium", colnames(zoops2019))]

#merge all files
all_zoops <- rbind(zoops2019, zoops2020, zoops2021)

#remove horizontal trap data
all_zoops <- all_zoops[!c(grepl("top", all_zoops$sample_ID) | 
                           grepl("bot", all_zoops$sample_ID)),]

#remove 20um samples
all_zoops <- all_zoops[all_zoops$mesh_size_μm!=20,]

#drop daphniidae and bsominidae bc I already have the genera
all_zoops <- all_zoops[,!c(grepl("Daphniidae", colnames(all_zoops)) |
                      grepl("Bosminidae", colnames(all_zoops)))]

#drop percent of total cols, count cols, and se cols
all_zoops <- all_zoops[,!c(grepl("Percent", colnames(all_zoops)) |
                             grepl("Count", colnames(all_zoops)) |
                             grepl("SE", colnames(all_zoops)))]

#drop the total zoop dens/biomass cols
all_zoops <- all_zoops[ ,-c(14:17)]

#now convert from wide to long for dens
all_zoops_dens <- all_zoops |> 
  pivot_longer(cols = colnames(all_zoops[grepl("density", colnames(all_zoops))]), names_to = "Taxon",
               values_to = "Density_IndPerL") 

#convert from wide to long for size
all_zoops_size <- all_zoops |> 
  pivot_longer(cols = colnames(all_zoops[grepl("MeanSize", colnames(all_zoops))]), names_to = "Taxon",
               values_to = "MeanLength_mm") 

#convert from wide to long for weight
all_zoops_weight <- all_zoops |> 
  pivot_longer(cols = colnames(all_zoops[grepl("totalbiomass_ug", colnames(all_zoops))]), names_to = "Taxon",
               values_to = "MeanWeight_ug") 

#convert from wide to long for biomass
all_zoops_biom <- all_zoops |> 
  pivot_longer(cols = colnames(all_zoops[grepl("BiomassConcentration", colnames(all_zoops))]), names_to = "Taxon",
               values_to = "Biomass_ugL") 

#combine dfs
all_zoops_final <- all_zoops_dens[,c(1:6, 77, 78)]
all_zoops_final$MeanLength_mm <- all_zoops_size$MeanLength_mm
all_zoops_final$MeanWeight_ug <- all_zoops_weight$MeanWeight_ug
all_zoops_final$Biomass_ugL <- all_zoops_biom$Biomass_ugL

all_zoops_final$Reservoir <- ifelse(substr(all_zoops_final$sample_ID,1,1) == "B", "BVR", "FCR")
all_zoops_final$Site <- ifelse(grepl("dam", all_zoops_final$sample_ID), 49,
                               ifelse(grepl("mac", all_zoops_final$sample_ID), 51, 50))
all_zoops_final$DateTime <- as.POSIXct(paste(all_zoops_final$collect_date, all_zoops_final$Hour), 
                                       format = "%Y-%m-%d %H:%M")

#drop density_nopl from end of taxon
all_zoops_final$Taxon <- substr(all_zoops_final$Taxon,1,nchar(all_zoops_final$Taxon)-13)

#add collection method col
all_zoops_final$CollectionMethod <- ifelse(grepl("schind", all_zoops_final$sample_ID), "Schindler", "Tow")

#rename depth of tow to StartDepth_m
all_zoops_final <- all_zoops_final |> rename(StartDepth_m = DepthOfTow_m)

#end depth is 0 for tows and start depth for schindlers 
all_zoops_final$EndDepth_m <- ifelse(all_zoops_final$CollectionMethod=="Tow", 0,
                                     all_zoops_final$StartDepth_m)

#drop unnecessary cols
all_zoops_final <- all_zoops_final |> select(-c(sample_ID, Project, site_no, collect_date, Hour))

#change order of cols
all_zoops_final <- all_zoops_final |> select(Reservoir, Site, DateTime, StartDepth_m, EndDepth_m,
                                             CollectionMethod, Taxon, Density_IndPerL, MeanLength_mm,
                                             MeanWeight_ug, Biomass_ugL)

#if density is 0, set mean length to 0
all_zoops_final$MeanLength_mm <- ifelse(all_zoops_final$Density_IndPerL==0, 0, all_zoops_final$MeanLength_mm)

#if density is not 0, but weight/biomass is 0, set to NA
all_zoops_final$MeanWeight_ug <- ifelse(all_zoops_final$Density_IndPerL!=0 & all_zoops_final$MeanWeight_ug==0, 
                                        NA, all_zoops_final$MeanWeight_ug)

all_zoops_final$Biomass_ugL <- ifelse(all_zoops_final$Density_IndPerL!=0 & all_zoops_final$Biomass_ugL==0, 
                                        NA, all_zoops_final$Biomass_ugL)

#add flag cols
all_zoops_final$Flag_Length <- ifelse(is.na(all_zoops_final$MeanLength_mm), 1, 0)
all_zoops_final$Flag_Weight <- ifelse(is.na(all_zoops_final$MeanWeight_ug), 1, 0)
all_zoops_final$Flag_Biomass <- ifelse(is.na(all_zoops_final$Biomass_ugL), 1, 0)

#export file 
write.csv(all_zoops_final, file.path(getwd(),'Summer2021-DataAnalysis/SummaryStats/EDI_zoop_taxa_2019-2022.csv'), row.names = FALSE)

#-------------------------------------------------------------------------------#
#raw density data sheet for EDI

#read in density csv from all 3 years
zoopdens19 <- read.csv(file.path(getwd(),'Summer2019-DataAnalysis/RawData/FCR2019_ZooplanktonCounting_Density_DataEntry.csv'),header = TRUE)
zoopdens20<- read.csv(file.path(getwd(),'Summer2020-DataAnalysis/RawData/FCR2020_ZooplanktonCounting_Density_DataEntry.csv'),header = TRUE)
zoopdens21<- read.csv(file.path(getwd(),'Summer2021-DataAnalysis/RawData/FCR2021_ZooplanktonCounting_Density_DataEntry.csv'),header = TRUE)

#combine files
zoopdens <- rbind(zoopdens19, zoopdens20, zoopdens21)

#drop 20um samples and horizontal trap samples
zoopdens <- zoopdens[!c(is.na(zoopdens$mesh_size_μm) | zoopdens$mesh_size_μm==20),]

#drop other test sample that doesn't have an hour
zoopdens <- zoopdens[zoopdens$Hour!="",]

#add reservoir
zoopdens$Reservoir <- ifelse(substr(zoopdens$sample_ID,1,1) == "B", "BVR", "FCR")

#add site
zoopdens$Site <- zoopdens$Site <- ifelse(grepl("dam", zoopdens$sample_ID), 49,
                                                ifelse(grepl("mac", zoopdens$sample_ID), 51, 50))

#add datetime
zoopdens$DateTime <- as.POSIXct(paste(zoopdens$collect_date, zoopdens$Hour), 
                                       format = "%d-%b-%y %H:%M")

#add collection method
zoopdens$CollectionMethod <- ifelse(grepl("schind", zoopdens$sample_ID), "Schindler", "Tow")

#start depth
zoopdens <- zoopdens |> rename(StartDepth_m = DepthOfTow_m)

#end depth is 0 for tows and start depth for schindlers 
zoopdens$EndDepth_m <- ifelse(zoopdens$CollectionMethod=="Tow", 0,
                              zoopdens$StartDepth_m)

#add column for rep # - need this bc not all reps have different collection times recorded...
zoopdens$Rep <- ifelse(grepl("rep",zoopdens$sample_ID), substrEnd(zoopdens$sample_ID,1), 1)

#drop unnecessary cols
zoopdens <- zoopdens |> select(-c(sample_ID, Project, site_no, collect_date, 
                                  Hour, date_processed, mesh_size_μm, INT, Notes,
                                  PhantomMidgePresence_YesNo))

#change order of cols
zoopdens <- zoopdens |> select(Reservoir, Site, DateTime, StartDepth_m, EndDepth_m, Rep,
                                             CollectionMethod, InitialSampleVolume_mL, Subsample,
                                             SubsampleVolume_mL, Zooplankton_No.)

#export density file
write.csv(zoopdens, file.path(getwd(),'Summer2021-DataAnalysis/SummaryStats/EDI_zoop_raw_dens_2019-2022.csv'), row.names = FALSE)

#-------------------------------------------------------------------------------#
#raw biomass data sheet for EDI

#read in density csv from all 3 years
zoopbiom19 <- read.csv(file.path(getwd(),'Summer2019-DataAnalysis/RawData/FCR2019_ZooplanktonCounting_SizeID_DataEntry.csv'),header = TRUE)
zoopbiom20<- read.csv(file.path(getwd(),'Summer2020-DataAnalysis/RawData/FCR2020_ZooplanktonCounting_SizeID_DataEntry.csv'),header = TRUE)
zoopbiom21<- read.csv(file.path(getwd(),'Summer2021-DataAnalysis/RawData/FCR2021_ZooplanktonCounting_SizeID_DataEntry.csv'),header = TRUE)

#drop 17jun19 bc was a practice sample
zoopbiom19 <- zoopbiom19[zoopbiom19$collect_date!="17-Jun-19",]

#combine files
zoopbiom <- rbind(zoopbiom19, zoopbiom20, zoopbiom21)

#add column for rep #
zoopbiom$Rep <- ifelse(grepl("rep",zoopbiom$sample_ID), substrEnd(zoopbiom$sample_ID,1), 1)

#drop 20um samples and horizontal trap samples
zoopbiom <- zoopbiom[!c(grepl("top", zoopbiom$sample_ID) |
                          grepl("bot", zoopbiom$sample_ID) |
                          grepl("_20_", zoopbiom$sample_ID)),]

#add reservoir
zoopbiom$Reservoir <- ifelse(substr(zoopbiom$sample_ID,1,1) == "B", "BVR", "FCR")

#add site
zoopbiom$Site <- zoopbiom$Site <- ifelse(grepl("dam", zoopbiom$sample_ID), 49,
                                         ifelse(grepl("mac", zoopbiom$sample_ID), 51, 50))

#add datetime
zoopbiom$DateTime <- as.POSIXct(paste(zoopbiom$collect_date, zoopbiom$Hour), 
                                format = "%d-%b-%y %H:%M")

#add collection method
zoopbiom$CollectionMethod <- ifelse(grepl("schind", zoopbiom$sample_ID), "Schindler", "Tow")

#drop unnecessary cols
zoopbiom <- zoopbiom |> select(-c(sample_ID, Project, site_no, collect_date, 
                                  Hour, MarksInOcularMicrometer_Width_No., 
                                  MarksInOcularMicrometer_Height_No., Initials, Notes))
                                  

#change order of cols
zoopbiom <- zoopbiom |> select(Reservoir, Site, DateTime, Rep,
                               CollectionMethod, Subsample, LowestTaxonomicLevelOfID,
                               TaxaID, Nauplius, ObjectiveMagnification, 
                               MarksInOcularMicrometer_No.)

#export density file
write.csv(zoopbiom, file.path(getwd(),'Summer2021-DataAnalysis/SummaryStats/EDI_zoop_raw_biom_2019-2022.csv'), row.names = FALSE)

