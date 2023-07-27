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

#-------------------------------------------------------------------------------#
#richness, evenness, and diversity during MOM vs anoxic days

#read in 2020-2022 zoop diversity data
zoop_div_2020 <- read.csv("Summer2020-DataAnalysis/SummaryStats/zoop_diversity_2020.csv", header=TRUE)
zoop_div_2021 <- read.csv("Summer2021-DataAnalysis/SummaryStats/zoop_diversity_2021.csv", header=TRUE)

#combine dfs
zoop_div <- rbind(zoop_div_2020, zoop_div_2021) 

#only pull schindler trap data
zoop_div_schind <- zoop_div |> filter(site_no %in% c("FCR_schind"), 
                                      sample_ID %in% c(zoop_div$sample_ID[grepl("noon",zoop_div$sample_ID)]))
                                
#add DO col
zoop_div_schind$DO <- ifelse(zoop_div_schind$collect_date %in% c("2020-09-11", "2020-09-15"), "MOM", "anoxic")

#summarise across DO and depth (note - taking out depth bc no crazy differences w/ depth)
zoop_div_schind_avg <- zoop_div_schind |> group_by(DO) |> #DepthOfTow_m
  summarise_at(vars(GenusRichness:OrderSubOrderEvenness), list(mean=mean))

#convert from wide to long
zoop_div_schind_avg_long <- zoop_div_schind_avg |> 
  pivot_longer(GenusRichness_mean:OrderSubOrderEvenness_mean)

#plot MOM vs anoxic diversity values (unfortunately not super interesting...)
ggplot(data=subset(zoop_div_schind_avg_long, name %in% 
                     c(zoop_div_schind_avg_long$name[grepl("Genus", zoop_div_schind_avg_long$name)])),
       aes(x=name, y=value, fill=name)) + geom_bar(stat="identity") +
  facet_wrap(~DO) + theme_bw() + 
  theme(text = element_text(size=8), axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), legend.key = element_blank(), 
        legend.key.height=unit(0.3,"line"),
        axis.text.x = element_text(angle = 45, hjust=1), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = c(0.1,0.44), legend.spacing = unit(-0.5, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.7,"line"))+ 
  guides(fill="none") + xlab("")

ggplot(data=subset(zoop_div_schind_avg_long, name %in% 
                     c(zoop_div_schind_avg_long$name[grepl("Family", zoop_div_schind_avg_long$name)])),
       aes(x=name, y=value, fill=name)) + geom_bar(stat="identity") +
  facet_wrap(~DO) + theme_bw() + 
  theme(text = element_text(size=8), axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), legend.key = element_blank(), 
        legend.key.height=unit(0.3,"line"),
        axis.text.x = element_text(angle = 45, hjust=1), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = c(0.1,0.44), legend.spacing = unit(-0.5, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.7,"line"))+ 
  guides(fill="none") + xlab("")

ggplot(data=subset(zoop_div_schind_avg_long, name %in% 
                     c(zoop_div_schind_avg_long$name[grepl("Order", zoop_div_schind_avg_long$name)])),
       aes(x=name, y=value, fill=name)) + geom_bar(stat="identity") +
  facet_wrap(~DO) + theme_bw() + 
  theme(text = element_text(size=8), axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), legend.key = element_blank(), 
        legend.key.height=unit(0.3,"line"),
        axis.text.x = element_text(angle = 45, hjust=1), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = c(0.1,0.44), legend.spacing = unit(-0.5, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.7,"line"))+ 
  guides(fill="none") + xlab("")

