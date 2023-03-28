
# potential for apparent competition or facilitation
# for every pair of species across every pair of sites

# -------------------------------------------------------------------------

library(tidyverse)

# -------------------------------------------------------------------------
library(foreach)
library(doParallel)
# 
workers <- 4
cl <- makeCluster(workers)
# register the cluster for using foreach
registerDoParallel(cl)

# -------------------------------------------------------------------------

site.info <- read.csv2("data/site_info_id.csv")
d1 <- read.csv2("data/melbourne_all_interactions.csv")
load("data/melbourne_network_list.RData")
sp.info <- read.csv2("data/melbourne_species_info.csv")
site.distances <- read.csv("data/site.distances.csv")
site.distances$neardist.norm <- scales::rescale(site.distances$neardist,to = c(0,1))
all.sites <- sort(site.info$site_id)
interaction.types <- names(net.list[[1]])

all.sp <- sort(unique(sp.info$species))
plant.sp <- subset(sp.info, functional.group == "Plant")
plant.sp <- plant.sp$species

pol.sp <- subset(sp.info, functional.group == "Pollinator")
pol.sp <- pol.sp$species

herb.sp <- subset(sp.info, functional.group == "Herbivore")
herb.sp <- herb.sp$species
# -------------------------------------------------------------------------
# template.df <- expand.grid(sp_1 = plant.sp,site_1 = all.sites,
#                            sp_2 = plant.sp,site_2 = all.sites)

# this part involves checking every combination of plant sp pairs
# at every site, and mediated by every animal species
# therefore, I need to code many nested loops

# 1 - pollination

# i.pol <- i.plant.1 <- i.site.1 <- 1
# i.plant.2 <- i.site.2 <- 2
# plant.sp <- plant.sp[1:32]

pol.list <- foreach(i.plant.1 = 1:length(plant.sp),.packages = "dplyr") %dopar%
  {
    # for(i.plant.1 in 1:length(plant.sp)){
    plant.1.list <- list()

    for(i.site.1 in 1:nrow(site.info)){
      my.net.1 <- net.list[[i.site.1]][["Pollinator"]]

      for(i.plant.2 in 1:length(plant.sp)){
        for(i.site.2 in 1:nrow(site.info)){
          my.net.2 <- net.list[[i.site.2]][["Pollinator"]]

          # the overall potential for apparent interaction is
          # the sum of individual contributions
          sum.app.int <- 0

          # calculate contributions from each animal species
          for(i.pol in 1:length(pol.sp)){

            # does the interaction occur at site 1 with i.plant.1?
            interaction_1_occurring <- which(my.net.1$vegspecies == plant.sp[i.plant.1] &
                                               my.net.1$species == pol.sp[i.pol])
            # site 2 with i.plant.2?
            interaction_2_occurring <- which(my.net.2$vegspecies == plant.sp[i.plant.2] &
                                               my.net.2$species == pol.sp[i.pol])

            if(length(interaction_1_occurring) > 0 &
               length(interaction_2_occurring) > 0){
              # if animal "i.pol" interacts with "i.plant.1" at site.1
              # and with "i.plant.2" at site.2
              # I move on to calculate the apparent interaction

              # first, I need the distance between both sites
              if(i.site.1 == i.site.2){
                my.dist <- 0
              }else{
                my.dist <- site.distances$neardist.norm[which(site.distances$siteA == site.info$site[i.site.1] &
                                                                site.distances$siteB == site.info$site[i.site.2])]
              }

              # second, transform this distance to an inter-layer strength
              # TODO this is just for testing, strength is simply 1-distance
              my.interlayer <- 1 - my.dist

              # third, obtain the other parts of the formula
              # relative effect of i.pol on i.plant.1 from site.1
              den.e1 <- my.net.1$freq[my.net.1$vegspecies == plant.sp[i.plant.1] &
                                        my.net.1$species == pol.sp[i.pol]]
              num.e1 <- sum(my.net.1$freq[my.net.1$vegspecies == plant.sp[i.plant.1]])

              # relative effect of i.plant.2 on i.pol from site.2
              den.e2 <- my.net.2$freq[my.net.2$vegspecies == plant.sp[i.plant.2] &
                                        my.net.2$species == pol.sp[i.pol]]
              num.e2 <- sum(my.net.2$freq[my.net.2$species == pol.sp[i.pol]])

              # add the apparent interaction from this specific pollinator i.pol
              sum.app.int <- sum.app.int + ((num.e1/den.e1) * (num.e2/den.e2) * my.interlayer)
            }# if interaction_1_occurring and interaction_2_occurring
          }# for i.pol

          # return(data.frame(v1 = plant.sp[i.plant.1],
          #        v2 = plant.sp[i.plant.2]))

          if(sum.app.int>0){
            # pol.list[[length(pol.list)+1]] <-
            plant.1.list[[length(plant.1.list)+1]] <- data.frame(sp_1 = plant.sp[i.plant.1],
                                                                 site_1 = site.info$site[i.site.1],
                                                                 sp_2 = plant.sp[i.plant.2],
                                                                 site_2 = site.info$site[i.site.2],
                                                                 interaction_type = "Pollination",
                                                                 apparent_interaction_value = sum.app.int)

          }# if apparent interaction > 0

        }# for i.site.2
      }# for i.plant.2
    }# for i.site

    plant.1.df <- dplyr::bind_rows(plant.1.list)
    return(plant.1.df)
  }# for i.plant.1

pol.df <- bind_rows(pol.list)

# -------------------------------------------------------------------------

# 1 - herbivory
# herb.list2 <- list()

# i.herb <- i.plant.1 <- i.site.1 <- 1
# i.plant.2 <- i.site.2 <- 2

# for(i.plant.1 in 1:length(plant.sp)){
herb.list <- foreach(i.plant.1 = 1:length(plant.sp),.packages = "dplyr") %dopar% 
  {
    
    plant.1.list <- list()
    
    for(i.site.1 in 1:nrow(site.info)){
      my.net.1 <- net.list[[i.site.1]][["Herbivore"]]
      
      for(i.plant.2 in 1:length(plant.sp)){
        for(i.site.2 in 1:nrow(site.info)){
          my.net.2 <- net.list[[i.site.2]][["Herbivore"]]
          
          # the overall potential for apparent interaction is
          # the sum of individual contributions
          sum.app.int <- 0
          
          # calculate contributions from each animal species
          for(i.herb in 1:length(herb.sp)){
            
            # does the interaction occur at site 1 with i.plant.1?
            interaction_1_occurring <- which(my.net.1$vegspecies == plant.sp[i.plant.1] &
                                               my.net.1$species == herb.sp[i.herb])
            # site 2 with i.plant.2?
            interaction_2_occurring <- which(my.net.2$vegspecies == plant.sp[i.plant.2] &
                                               my.net.2$species == herb.sp[i.herb])
            
            if(length(interaction_1_occurring) > 0 & 
               length(interaction_2_occurring) > 0){
              # if animal "i.herb" interacts with "i.plant.1" at site.1
              # and with "i.plant.2" at site.2
              # I move on to calculate the apparent interaction
              
              # first, I need the distance between both sites
              if(i.site.1 == i.site.2){
                my.dist <- 0
              }else{
                my.dist <- site.distances$neardist.norm[which(site.distances$siteA == site.info$site[i.site.1] &
                                                                site.distances$siteB == site.info$site[i.site.2])]
              }
              
              # second, transform this distance to an inter-layer strength
              # TODO this is just for testing, strength as simply 1-distance
              my.interlayer <- 1 - my.dist
              
              # third, obtain the other parts of the formula
              # relative effect of i.herb on i.plant.1 from site.1
              den.e1 <- my.net.1$freq[my.net.1$vegspecies == plant.sp[i.plant.1] &
                                        my.net.1$species == herb.sp[i.herb]]
              num.e1 <- sum(my.net.1$freq[my.net.1$vegspecies == plant.sp[i.plant.1]])
              
              # relative effect of i.plant.2 on i.herb from site.2
              den.e2 <- my.net.2$freq[my.net.2$vegspecies == plant.sp[i.plant.2] &
                                        my.net.2$species == herb.sp[i.herb]]
              num.e2 <- sum(my.net.2$freq[my.net.2$species == herb.sp[i.herb]])
              
              # add the apparent interaction from this specific pollinator i.herb
              sum.app.int <- sum.app.int + ((num.e1/den.e1) * (num.e2/den.e2) * my.interlayer)
            }# if interaction_1_occurring and interaction_2_occurring
          }# for i.herb
          
          if(sum.app.int>0){
            # herb.list[[length(herb.list)+1]] <- 
            plant.1.list[[length(plant.1.list)+1]] <- data.frame(sp_1 = plant.sp[i.plant.1],
                                                                        site_1 = site.info$site[i.site.1],
                                                                        sp_2 = plant.sp[i.plant.2],
                                                                        site_2 = site.info$site[i.site.2],
                                                                        interaction_type = "Herbivory",
                                                                        apparent_interaction_value = sum.app.int)
          }# if apparent interaction > 0
        }# for i.site.2
      }# for i.plant.2
    }# for i.site
    plant.1.df <- dplyr::bind_rows(plant.1.list)
    return(plant.1.df)
  }# for i.plant.1

herb.df <- bind_rows(herb.list)

# -------------------------------------------------------------------------
stopCluster(cl)
# -------------------------------------------------------------------------

write.csv2(pol.df,"results/potential_apparent_facilitation.csv",row.names = F)
write.csv2(herb.df,"results/potential_apparent_competition.csv",row.names = F)









