
library(tidyverse)
library(igraph)
library(bipartite)
source("R/matrix_from_edge_list.R")

# -------------------------------------------------------------------------

d1 <- read.csv2("data/melbourne_all_interactions.csv")
load("data/melbourne_network_list.RData")
sp.info <- read.csv2("data/melbourne_species_info.csv")

sites <- names(net.list)
interaction.types <- names(net.list[[1]])
all.sp <- sort(unique(sp.info$species))

# -------------------------------------------------------------------------
# for each network, obtain a series of metrics
# i.site <- i.type <- 2

metrics.list <- list()
for(i.site in 1:length(sites)){
  metrics.list[[i.site]] <- list()
  for(i.type in 1:length(interaction.types)){
    metrics.list[[i.site]][[i.type]] <- list()
  } # for i.type
  names(metrics.list[[i.site]]) <- interaction.types
} # for i.site
names(metrics.list) <- sites

for(i.site in 1:length(sites)){
  for(i.type in 1:length(interaction.types)){
    
    # my.graph <- graph_from_data_frame(net.list[[i.site]][[i.type]])
    my.net <- as.data.frame(net.list[[i.site]][[i.type]])
    names(my.net) <- c("node_to","node_from","link_value")
    
    if(nrow(my.net)>0){
      
      my.matrix <- matrix_from_edge_list(my.net,weighted = T)
      
      if(ncol(my.matrix) > 1 & nrow(my.matrix) > 1){
        my.metrics <- bipartite::specieslevel(my.matrix,
                                              index = c("normalised degree",
                                                        "species strength",
                                                        "betweenness",
                                                        "closeness",
                                                        "d"),
                                              level = "both")
      }else{
        my.metrics <- bipartite::specieslevel(my.matrix,
                                              index = c("normalised degree",
                                                        "species strength",
                                                        # "betweenness",
                                                        # "closeness",
                                                        "d"),
                                              level = "both")
        my.metrics[[1]]$betweenness <- NA
        my.metrics[[2]]$betweenness <- NA
        my.metrics[[1]]$weighted.betweenness <- NA
        my.metrics[[2]]$weighted.betweenness <- NA
        my.metrics[[1]]$closeness <- NA
        my.metrics[[2]]$closeness <- NA
        my.metrics[[1]]$weighted.closeness <- NA
        my.metrics[[2]]$weighted.closeness <- NA
      }
    }else{
      my.metrics <- NA
    }

    metrics.list[[i.site]][[i.type]] <- my.metrics
    
  }# i.type
}# i.site

# -------------------------------------------------------------------------
# transform nested list to dataframe

metrics.df.list <- list()
metric.names <- names(metrics.list[[1]][[1]][[1]])

# i.site <- i.type <- 2
# i.level <- "higher level"

for(i.site in 1:length(sites)){
  for(i.type in 1:length(interaction.types)){
    for(i.level in c("lower level","higher level")){
      
      my.sp <-unique(net.list[[i.site]][[i.type]]$vegspecies)
      if(i.level == "higher level"){
        my.sp <- unique(net.list[[i.site]][[i.type]]$species)
      }
      
      if(length(my.sp)>0){
        
        net.df <- data.frame(site = sites[i.site],
                             interaction_type = interaction.types[i.type],
                             species = my.sp)
        
        for(i.metric in metric.names){
          net.df[,ncol(net.df)+1] <- metrics.list[[i.site]][[i.type]][[i.level]][[i.metric]]
        }# for i.metric
        names(net.df)[4:ncol(net.df)] <- metric.names
        
        metrics.df.list[[length(metrics.df.list)+1]] <- net.df
      }# if species in the network
      
    }# for i.level
  }# for i.type
}# for i.site


metrics.df <- bind_rows(metrics.df.list)
metrics.df.long <- pivot_longer(metrics.df,cols = 4:ncol(metrics.df),names_to = "metric",values_to = "value")

metrics.df.long$species_functional_group <- sp.info$functional.group[match(metrics.df.long$species,sp.info$species)]

# -------------------------------------------------------------------------
# clean up a bit
metrics.df.long$interaction_type[metrics.df.long$interaction_type == "Herbivore"] <- "Herbivory"
metrics.df.long$interaction_type[metrics.df.long$interaction_type == "Pollinator"] <- "Pollination"

metrics.df.long <- metrics.df.long[,c("site","interaction_type","species","species_functional_group","metric","value")]

# -------------------------------------------------------------------------

write.csv2(metrics.df.long,file = "results/melbourne_species_metrics.csv",row.names = F)



