#Simulation of random data with different numbers of groups
#step 1 - randomly select 100 values 
#step 2 - randomly group these points using different number of groups
#step 3 - calculate variability within each group 
#step 4 - plot variability vs group #

#read in libraries
pacman::p_load(ggplot2, vegan, ggordiplots, splancs)

vec <- data.frame("vals1" = runif(110, 0.0, 1),
                  "vals2" = runif(110, 0.0, 1),
                  "vals3" = runif(110, 0.0, 1),
                  "vals4" = runif(110, 0.0, 1),
                  "vals5" = runif(110, 0.0, 1),
                  "vals6" = runif(110, 0.0, 1),
                  "vals7" = runif(110, 0.0, 1),
                  "vals8" = runif(110, 0.0, 1),
                  "vals9" = runif(110, 0.0, 1),
                  "vals10" = runif(110, 0.0, 1),
                  "vals11" = runif(110, 0.0, 1))

grps <- data.frame("n"=rep(NA,110))

for(j in 2:100){
    temp_grp <- as.numeric(cut_number(vec$vals1, j))
      grps[,paste0("grp_",j)] <- temp_grp
}

#create new df to save betadisper results for all different groupings
var_results <- data.frame("disp" = rep(NA,100),
                          "pairwise" = rep(NA,100),
                          "area" = rep(NA,100))

for (i in 2:length(grps)){
  message(paste0("Running group ",i))
  temp_dist <- vegdist(vec, method="euclidean")
  var <- betadisper(temp_dist, group = as.factor(grps[,i]), type="centroid")
  var_results$disp[i] <- mean(var$group.distances)
  var_results$pairwise[i] <- mean(dist(var$centroids))
  nmds_temp <- metaMDS(temp_dist, distance='bray', k=4, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
  if(i<50){
    message(paste0("Area ",i))
    ord_temp <- ordiplot(nmds_temp,display = c('sites','species'),choices = c(1,2),type = "n")
    grp_temp <- gg_ordiplot(ord_temp, grps[,i], kind = "ehull", 
                           ellipse=FALSE, hull = TRUE, plot = FALSE, pt.size=0.9)
    areas = c()
    for(j in 1:i){
      area_temp <- areapl(cbind(grp_temp$df_hull$x[grp_temp$df_hull$Group==j],
                                grp_temp$df_hull$y[grp_temp$df_hull$Group==j]))
      areas = c(areas,area_temp)
    }
    var_results$area[i] <- mean(areas)
  }else{
    var_results$area[i] <- 0
  }
}


#test to make sure for loop is working the way I think it is
#temp_dist <- vegdist(vec[,7], method="euclidean")
#var <- betadisper(temp_dist, group = as.factor(vec[,7]), type="centroid")
#test$disp[7] <- mean(var$group.distances)
#test$pairwise[7] <- mean(dist(var$centroids))

#add group num col
var_results$grp_num <- c(1:100)

#plot variability vs number of groups (n = 2-100)
ggplot(var_results, aes(grp_num, disp)) + geom_point() + theme_bw()

ggplot(var_results, aes(grp_num, pairwise)) + geom_point() + theme_bw()

ggplot(var_results, aes(grp_num, area)) + geom_point() + theme_bw()

#sooo definitely need to bootstrap because the more groups in a df, the lower variability is...

#------------------------------------------------------------------------------#
# trying again with my data rather than randomn data
grps <- data.frame("n"=rep(NA,110))

for(j in 2:100){
  temp_grp <- as.numeric(cut_number(zoop_temporal_dens_avg_trans$Cyclopoida_density_NopL_1_mean, j))
  grps[,paste0("grp_",j)] <- temp_grp
}

#create new df to save betadisper results for all different groupings
var_results <- data.frame("disp"=rep(NA,100))

for (i in 2:length(grps)){
  temp_dist <- vegdist(zoop_temporal_dens_avg_trans, method="euclidean")
  var <- betadisper(temp_dist, group = as.factor(grps[,i]), type="centroid")
  var_results$disp[i] <- mean(var$group.distances)
  var_results$pairwise[i] <- mean(dist(var$centroids))
}

#add group num col
var_results$grp_num <- c(1:100)

#plot variability vs number of groups (n = 2-100)
ggplot(var_results, aes(grp_num, disp)) + geom_point() + theme_bw()

ggplot(var_results, aes(grp_num, pairwise)) + geom_point() + theme_bw()

