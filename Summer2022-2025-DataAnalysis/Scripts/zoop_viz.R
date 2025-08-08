# Zoop visualization script 
# Created 07 Aug 2025 - HLW

pacman::p_load(tidyverse, ggplot2, lubridate)

new_zoops <- read.csv("SummaryStats/FCR_ZooplanktonSummary2022-2025.csv") |>
  select(-c(site_no,sample_ID,Project,InitialSampleVolume_mL,Zooplankton_No.,INT,Volume_L,
            Volume_unadj,proportional_vol,mesh_size_μm))

#visualize total zoop metrics

#biomass
ggplot(new_zoops, aes(x = as.Date(collect_date), y = BiomassConcentration_ugpL)) +
  geom_line() +
  geom_point() +
  labs(title = "Zooplankton Biomass Concentration Over Time",
       x = "Date",
       y = "Biomass (µg/L)") +
  theme_minimal()

#density
ggplot(new_zoops, aes(x = as.Date(collect_date), y = ZoopDensity_No.pL)) +
  geom_line() +
  geom_point() +
  labs(title = "Zooplankton Density Over Time",
       x = "Date",
       y = "Density (#/L)") +
  theme_minimal()

#size
ggplot(new_zoops, aes(x = as.Date(collect_date), y = OverallMeanSize_mm)) +
  geom_point() +
  geom_errorbar(aes(ymin = OverallMeanSize_mm - OverallSESize_mm,
                    ymax = OverallMeanSize_mm + OverallSESize_mm),
                width = 0.2) +
  labs(title = "Zooplankton Size Over Time",
       x = "Date",
       y = "Mean Size (mm)") +
  theme_minimal()

#changes in zoop biomass w/ tow depth
ggplot(new_zoops, aes(x = as.Date(collect_date), y = DepthOfTow_m, size = BiomassConcentration_ugpL)) +
  geom_point(alpha = 0.7) +
  scale_y_reverse() +
  labs(title = "Vertical Distribution of Biomass",
       x = "Date",
       y = "Tow Depth (m)",
       size = "Biomass (µg/L)") +
  theme_minimal()

#changes in zoop density w/ tow depth
ggplot(new_zoops, aes(x = as.Date(collect_date), y = DepthOfTow_m, size = ZoopDensity_No.pL)) +
  geom_point(alpha = 0.7) +
  scale_y_reverse() +
  labs(title = "Vertical Distribution of Density",
       x = "Date",
       y = "Tow Depth (m)",
       size = "Density (#/L)") +
  theme_minimal()

#another way of visualizing how density changes with depth of tow
ggplot(new_zoops, aes(x = ZoopDensity_No.pL, y = DepthOfTow_m)) +
  geom_line() +
  geom_point() +
  scale_y_reverse() +
  labs(title = "Zooplankton Density Profiles",
       x = "Density (No./pL)",
       y = "Tow Depth (m)") +
  theme_minimal()

#biomass vs. density - not really a 1:1 relationship
ggplot(new_zoops, aes(x = ZoopDensity_No.pL, y = BiomassConcentration_ugpL)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Biomass vs. Density",
       x = "Density (No./pL)",
       y = "Biomass (µg/L)") +
  theme_minimal()

#-----------------------------------------------------------------------------#
#list of taxa 
taxa<-c("Daphniidae","Copepoda","Calanoida","Cladocera","Cyclopoida","Rotifera",
        "Keratella","Kellicottia","Crustacea","Bosminidae","nauplius",
        "Ceriodaphnia","Daphnia","Bosmina","Ploima","Gastropidae",
        "Conochilidae","Synchaetidae","Trichocercidae","Monostyla","Lecane")

taxa_groupings <- data.frame(
  Taxon = c("Crustacea","Daphniidae", "Cladocera", "Bosminidae", "Ceriodaphnia", "Daphnia", "Bosmina", 
            "Copepoda", "Calanoida", "Cyclopoida", "nauplius",
            "Rotifera", "Keratella", "Kellicottia", "Ploima", "Gastropidae", 
            "Conochilidae", "Synchaetidae", "Trichocercidae", "Monostyla", "Lecane"),
  Group = c(rep("NA", 1),
            rep("Cladocera", 6),
            rep("Copepod", 4),
            rep("Rotifer", 10)))

# Reshape size
size_long <- new_zoops %>%
  select(collect_date, DepthOfTow_m, Hour, contains("MeanSize_mm")) %>%
  pivot_longer(
    cols = contains("MeanSize_mm"),
    names_to = "Taxon",
    values_to = "MeanSize_mm"
  ) %>%
  mutate(Taxon = gsub("MeanSize_mm", "", Taxon)) |>
  filter(!Taxon %in% "Overall",
         !is.na(MeanSize_mm)) |>
  left_join(taxa_groupings, by = "Taxon")

# Reshape biomass
biomass_long <- new_zoops %>%
  select(collect_date, DepthOfTow_m, Hour, contains("_totalbiomass_ug")) %>%
  pivot_longer(
    cols = contains("_totalbiomass_ug"),
    names_to = "Taxon",
    values_to = "TotalBiomass_ug"
  ) %>%
  mutate(Taxon = gsub("_totalbiomass_ug", "", Taxon)) |>
  left_join(taxa_groupings, by = "Taxon")

# Reshape density
density_long <- new_zoops %>%
  select(collect_date, DepthOfTow_m, Hour, contains("_density_NopL")) %>%
  pivot_longer(
    cols = contains("_density_NopL"),
    names_to = "Taxon",
    values_to = "Density_NopL"
  ) %>%
  mutate(Taxon = gsub("_density_NopL", "", Taxon)) |>
  left_join(taxa_groupings, by = "Taxon")

#now visualize density
ggplot(density_long, aes(x = as.Date(collect_date), y = Density_NopL, color = Group)) +
  geom_point(alpha = 0.6) +
  geom_line() +
  facet_wrap(~ Taxon, scales = "free_y") +
  labs(title = "Zooplankton Density Over Time",
       x = "Date",
       y = "Density (#/L)") +
  theme_minimal()

#biomass
ggplot(biomass_long, aes(x = as.Date(collect_date), y = TotalBiomass_ug, color = Group)) +
  geom_point(alpha = 0.6) +
  geom_line() +
  facet_wrap(~ Taxon, scales = "free_y") +
  labs(title = "Zooplankton Biomass Over Time",
       x = "Date",
       y = "Biomass (µg)") +
  theme_minimal()

#size
ggplot(size_long, aes(x = as.Date(collect_date), y = MeanSize_mm, color = Group)) +
  geom_point(alpha = 0.6) +
  geom_line() +
  facet_wrap(~ Taxon, scales = "free_y") +
  labs(title = "Mean Zooplankton Size Over Time",
       x = "Date",
       y = "Mean Size (mm)") +
  theme_minimal()

#------------------------------------------------------------------------------#
# visualize the full zooplankton dataset (2016-2025)








