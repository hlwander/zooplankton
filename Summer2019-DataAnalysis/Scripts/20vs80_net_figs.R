# Figs comparing 20 vs 80 um zoop net data
# 23Jun2023

#read in 2019 data
zoops2020<- read.csv(file.path(getwd(),'Summer2020-DataAnalysis/SummaryStats/FCR_ZooplanktonSummary2020.csv'),header = TRUE)

#pull data for days where we collected both 20um and 80um net samples (29Jun20)
zoop_comp <- zoops2020[(zoops2020$collect_date=="2020-06-29" | zoops2020$collect_date=="2020-06-28") &
                         zoops2020$DepthOfTow_m==2.6,]

#only keep the midnight comparison tow and filtered 20um reps
zoop_comp <- zoop_comp[grepl("_filt",zoop_comp$sample_ID) | zoop_comp$collect_date=="2020-06-28",]

#average the 3 20um reps
zoops_20v80 <- zoop_comp |>  
  group_by(site_no, Hour, collect_date) |>
  summarise_at(vars(ZoopDensity_No.pL:Synchaetidae_BiomassConcentration_ugpL), list(mean=mean)) 

zoops_20v80 <- zoops_20v80 |> 
  select(site_no, Hour, collect_date, which(grepl("density", colnames(zoops_20v80), ignore.case=TRUE)))

#add mesh size 
zoops_20v80$mesh <- c(20,80)

zoops_20v80_long <- zoops_20v80 |> 
  pivot_longer(cols = ZoopDensity_No.pL_mean:Synchaetidae_density_NopL_mean, names_to = "variable")

#change facet labels
metric_taxa <-c("Total", "Daphniidae", "Copepoda", "Calanoida", "Cladocera", "Cyclopoida", "Rotifera",
                "Keratella", "Kellicottia","Crustacea", "Bosminidae","Nauplius","Ceriodaphnia",
                "Daphnia","Bosmina", "Ploima", "Gastropidae", "Collothecidae", "Conochilidae",
                "Synchaetidae")
names(metric_taxa) <- c(unique(zoops_20v80_long$variable))

#now plot densities
mesh_comp <- ggplot(zoops_20v80_long, aes(x=as.factor(mesh), y=value)) + geom_bar(stat="identity") +
                facet_wrap(~variable, labeller = labeller(variable=metric_taxa)) + 
                theme_bw() + xlab("mesh") + ylab("Density (#/L)") +
                theme(text = element_text(size=6), axis.text = element_text(size=5, color="black"), 
                      legend.background = element_blank(), legend.key = element_blank(), 
                      strip.background = element_rect(fill = "transparent"), 
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank())
ggsave(file.path(getwd(),"Summer2021-DataAnalysis/Figures/20vs80_density_comp.jpg"),
       mesh_comp, width=3, height=3) 
#take-home is that 20um gets clogged way to quickly, even for a 2.6m tow (so full water column must be pretty bad...)


