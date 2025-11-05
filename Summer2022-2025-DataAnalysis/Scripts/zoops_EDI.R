#script for getting all summary files into EDI format 2022-2025
# 7 August 2025

#read in libraries
pacman::p_load(tidyr, dplyr, stringr)

#set wd
setwd("/Users/heatherwander/Documents/VirginiaTech/research/zooplankton/Summer2022-2025-DataAnalysis")

#function to count characters starting at the end of the string
substrEnd <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#read in zoop data from most recent year
zoops_new <- read.csv(file.path(getwd(),"SummaryStats/BVR_ZooplanktonSummary2022-2025.csv"))

#drop daphniidae and bsominidae bc I already have the genera
zoops_new <- zoops_new[,!c(grepl("Daphniidae", colnames(zoops_new)) |
                      grepl("Bosminidae", colnames(zoops_new)))]

#drop percent of total cols, count cols, and se cols
zoops_new <- zoops_new[,!c(grepl("Percent", colnames(zoops_new)) |
                             grepl("Count", colnames(zoops_new)) |
                             grepl("SE", colnames(zoops_new)))]

#drop the total zoop dens/biomass cols
zoops_new <- zoops_new[ ,-c(14:17)]

#now convert from wide to long for dens
zoops_dens <- zoops_new |> 
  pivot_longer(cols = colnames(zoops_new[grepl("density", colnames(zoops_new))]), names_to = "Taxon",
               values_to = "Density_IndPerL") 

#convert from wide to long for size
zoops_size <- zoops_new |> 
  pivot_longer(cols = colnames(zoops_new[grepl("MeanSize", colnames(zoops_new))]), names_to = "Taxon",
               values_to = "MeanLength_mm") 

#convert from wide to long for weight
zoops_weight <- zoops_new |> 
  pivot_longer(cols = colnames(zoops_new[grepl("totalbiomass_ug", colnames(zoops_new))]), names_to = "Taxon",
               values_to = "MeanWeight_ug") 

#convert from wide to long for biomass
zoops_biom <- zoops_new |> 
  pivot_longer(cols = colnames(zoops_new[grepl("BiomassConcentration", colnames(zoops_new))]), names_to = "Taxon",
               values_to = "Biomass_ugL") 

#combine dfs
zoops_final <- zoops_dens[,c(1:6, 80, 81)]
zoops_final$MeanLength_mm <- zoops_size$MeanLength_mm
zoops_final$MeanWeight_ug <- zoops_weight$MeanWeight_ug
zoops_final$Biomass_ugL <- zoops_biom$Biomass_ugL

zoops_final$Reservoir <- ifelse(substr(zoops_final$sample_ID,1,1) == "B", "BVR", "FCR")
zoops_final$Site <- 50
zoops_final$DateTime <- as.POSIXct(paste(zoops_final$collect_date, zoops_final$Hour), 
                                       format = "%Y-%m-%d %H:%M")

#drop density_nopl from end of taxon
zoops_final$Taxon <- substr(zoops_final$Taxon,1,nchar(zoops_final$Taxon)-13)

#add rep col (all these tows only have 1 rep but consider changing this to an ifelse statement for future sampling campaigns)
zoops_final$Rep <- 1

#add collection method col
zoops_final$CollectionMethod <- "Tow" #all tows for these years

#rename depth of tow to StartDepth_m
zoops_final <- zoops_final |> rename(StartDepth_m = DepthOfTow_m)

#end depth is 0 for tows and start depth for schindlers 
zoops_final$EndDepth_m <- 0

#drop unnecessary cols
zoops_final <- zoops_final |> select(-c(sample_ID, Project, site_no, collect_date, Hour))

#change order of cols
zoops_final <- zoops_final |> select(Reservoir, Site, DateTime, StartDepth_m, EndDepth_m, Rep,
                                             CollectionMethod, Taxon, Density_IndPerL, MeanLength_mm,
                                             MeanWeight_ug, Biomass_ugL)

#if density is 0, set mean length to 0
zoops_final$MeanLength_mm <- ifelse(zoops_final$Density_IndPerL==0, 0, zoops_final$MeanLength_mm)

#if density is not 0, but weight/biomass is 0, set to NA
zoops_final$MeanWeight_ug <- ifelse(zoops_final$Density_IndPerL!=0 & zoops_final$MeanWeight_ug==0, 
                                        NA, zoops_final$MeanWeight_ug)

zoops_final$Biomass_ugL <- ifelse(zoops_final$Density_IndPerL!=0 & zoops_final$Biomass_ugL==0, 
                                        NA, zoops_final$Biomass_ugL)

#add flag cols
zoops_final$Flag_MeanLength_mm <- ifelse(is.na(zoops_final$MeanLength_mm), 1, 0)
zoops_final$Flag_MeanWeight_ug <- ifelse(is.na(zoops_final$MeanWeight_ug), 1, 0)
zoops_final$Flag_Biomass_ugL <- ifelse(is.na(zoops_final$Biomass_ugL), 1, 0)

#export file 
write.csv(zoops_final, file.path(getwd(),'SummaryStats/EDI_zoop_summary_2022-2025.csv'), row.names = FALSE)

#-------------------------------------------------------------------------------#
#raw density data sheet for EDI

#read in density csv for most recent year
zoops_dens_new <- read.csv(file.path(getwd(),"RawData/BVR2022-2025_ZooplanktonCounting_Density_DataEntry.csv"),header = TRUE)

#add reservoir
zoops_dens_new$Reservoir <- ifelse(substr(zoops_dens_new$sample_ID,1,1) == "B", "BVR", "FCR")

#add site
zoops_dens_new$Site <- 50

#add datetime
zoops_dens_new$DateTime <- as.POSIXct(paste(zoops_dens_new$collect_date, zoops_dens_new$Hour), 
                                       format = "%d-%b-%y %H:%M")

#add collection method
zoops_dens_new$CollectionMethod <- "Tow"

#start depth
zoops_dens_new <- zoops_dens_new |> rename(StartDepth_m = DepthOfTow_m)

#end depth is 0 for tows and start depth for schindlers 
zoops_dens_new$EndDepth_m <- 0

#add column for rep # - need this bc not all reps have different collection times recorded...
zoops_dens_new$Rep <- 1

#drop unnecessary cols
zoops_dens_new <- zoops_dens_new |> select(-c(sample_ID, Project, site_no, collect_date, 
                                  Hour, date_processed, mesh_size_μm, INT, Notes,
                                  PhantomMidgePresence_YesNo))

#change order of cols
zoops_dens_new <- zoops_dens_new |> select(Reservoir, Site, DateTime, StartDepth_m, EndDepth_m, Rep,
                                             CollectionMethod, InitialSampleVolume_mL, Subsample,
                                             SubsampleVolume_mL, Zooplankton_No.)

#export density file
write.csv(zoops_dens_new, file.path(getwd(),'SummaryStats/EDI_zoop_raw_dens_2022-2025.csv'), row.names = FALSE)

#-------------------------------------------------------------------------------#
#raw biomass data sheet for EDI

#read in biom csv for most recent year
zoops_biom_new <- read.csv(file.path(getwd(),"RawData/BVR2022-2025_ZooplanktonCounting_SizeID_DataEntry.csv"),header = TRUE)

#add column for rep #
zoops_biom_new$Rep <- 1

#add reservoir
zoops_biom_new$Reservoir <- ifelse(substr(zoops_biom_new$sample_ID,1,1) == "B", "BVR", "FCR")

#add site
zoops_biom_new$Site <- 50

#add datetime
zoops_biom_new$DateTime <- as.POSIXct(paste(zoops_biom_new$collect_date, zoops_biom_new$Hour), 
                                format = "%d-%b-%y %H:%M")

#start depth
zoops_biom_new <- zoops_biom_new |>
  mutate(StartDepth_m = substrEnd(sample_ID,4) |>
           str_remove(".*_") |> # remove everything before underscore
           str_sub(1, -2))      # drop the m at the end of each depth
  
#end depth is 0 for tows and start depth for schindlers 
zoops_biom_new$EndDepth_m <- 0.1

#add collection method
zoops_biom_new$CollectionMethod <- "Tow"

zoops_biom_new$OcularMagnification <- "10x"

#drop unnecessary cols
zoops_biom_new <- zoops_biom_new |> select(-c(sample_ID, Project, site_no, collect_date, 
                                  Hour, MarksInOcularMicrometer_Width_No., 
                                  MarksInOcularMicrometer_Height_No., Initials, Notes))
                                  

#change order of cols
zoops_biom_new <- zoops_biom_new |> select(Reservoir, Site, DateTime, StartDepth_m,
                                           EndDepth_m, Rep, CollectionMethod, 
                                           Subsample, LowestTaxonomicLevelOfID,
                                            TaxaID, Nauplius, ObjectiveMagnification, 
                                           OcularMagnification, MarksInOcularMicrometer_No.)

#export density file
write.csv(zoops_biom_new, file.path(getwd(),'SummaryStats/EDI_zoop_raw_biom_2022-2025.csv'), row.names = FALSE)

