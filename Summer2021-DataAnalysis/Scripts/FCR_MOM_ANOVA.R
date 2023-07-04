# FCR MOM BACI
# 4 Jul 2023

#Note: trying this out, but not necessarily a before/after case because it goes back and forth between anoxia and mom most years depending on sss usage
#using MOM schindler density data bc I trust that more given high net efficiency

#A bit of a spin on the traditional BACI design where before = noon, after = midnight, control = anoxia, impact = MOM

#read in libraries
pacman::p_load(dplyr, ggpubr, ARTool)

#read in FCR MOM schindler data
zoop_dens <- read.csv("Summer2021-DataAnalysis/SummaryStats/FCR_MOM_schind_2020-2022_zoopdens_3groups.csv", header=TRUE)

#add anoxic vs mom col
zoop_dens$DO <- ifelse(zoop_dens$collect_date %in% c("2020-09-11", "2020-09-15"), "MOM", "anoxic")

#histogram of densities 
hist(zoop_dens[, which(colnames(zoop_dens) %in% c("Cladocera_density_NopL_rep.mean"))])
hist(zoop_dens[, which(colnames(zoop_dens) %in% c("Copepoda_density_NopL_rep.mean"))])
hist(zoop_dens[, which(colnames(zoop_dens) %in% c("Rotifera_density_NopL_rep.mean"))])

#wide to long df
zoop_dens_long <- zoop_dens |> 
  pivot_longer(Cladocera_density_NopL_rep.mean:Rotifera_density_NopL_rep.mean)

#get rid of midnight samples so can test DO * depth
zoop_dens_long <- zoop_dens_long |> filter(Hour %in% c(12))

#visualize the data
ggplot(zoop_dens_long, aes(value, depth, color = as.factor(collect_date))) + geom_path() + 
  facet_wrap(~name+DO, nrow = 3) + ylim(9,0)

#check that data are normally distributed
ggpubr::ggqqplot(log(zoop_dens_long$value[zoop_dens_long$name=="Cladocera_density_NopL_rep.mean"]))
ggpubr::ggqqplot(log(zoop_dens_long$value[zoop_dens_long$name=="Copepoda_density_NopL_rep.mean"]))
ggpubr::ggqqplot(log(zoop_dens_long$value[zoop_dens_long$name=="Rotifera_density_NopL_rep.mean"]))

shapiro.test(log(zoop_dens_long$value))
# p < 0.05 so not normal + need nonparametric test

#Aligned Ranks Transformation (ART) ANOVA (depth, time, and DO as categorical variables)

#step1: make all independent variables factors
zoop_dens_long$depth <- as.factor(zoop_dens_long$depth)
zoop_dens_long$Hour <- as.factor(zoop_dens_long$Hour)
zoop_dens_long$DO <- as.factor(zoop_dens_long$DO)

#now used linear mixed model syntax

m <- art(value ~ DO * depth + (1|name), data = zoop_dens_long)

anova(m)
#cool - anoxia vs MOM affects zoop density

mean(zoop_dens_long$value[zoop_dens_long$DO=="anoxic"])
mean(zoop_dens_long$value[zoop_dens_long$DO=="MOM"])
#so greater density during anoxic conditions interesting


