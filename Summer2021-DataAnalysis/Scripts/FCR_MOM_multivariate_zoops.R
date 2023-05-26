#script for multivariate stats using FCR MOM + environmental correlate data
#created 26 May 2023

#read in libraries
pacman::p_load(tidyverse, labdsv, goeveg, vegan, ggordiplots, viridis)

#read in FCR MOM schindler data
zoop_dens <- read.csv("Summer2021-DataAnalysis/SummaryStats/FCR_MOM_schind_2020-2022_zoopdens.csv", header=TRUE)
zoop_biom <- read.csv("Summer2021-DataAnalysis/SummaryStats/FCR_MOM_schind_2020-2022_zoopbiom.csv", header=TRUE)

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
             color="black", pch=21, size=2) +
  theme(text = element_text(size=7), axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), 
        legend.key.height=unit(0.3,"line"),
        legend.key = element_blank(),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-0,-0,-0,-0),
        legend.direction = "horizontal",
        axis.text.x = element_text(vjust = 0.5), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "top", legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(0,-0.1,0,0), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) + guides(fill="none") +
        annotate("text", x=-0.42, y=0.6, label= "bolditalic(a)", parse=T, size = 3) +
        scale_fill_manual("",values=viridis(9, option="F"))+
        scale_color_manual("",values=viridis(9, option="F"),
                     label=c('littoral','pelagic')) 
