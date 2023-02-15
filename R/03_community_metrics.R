
library(igraph) # for centralization
library(maxnodf) # for nestedness
library(tidyverse)
library(purrr)
library(MultitrophicFun)
source("R/metric_functions.R")
source("R/matrix_from_edge_list.R")
source("R/reshuffle_matrix.R")

# metrics and, optionally, z-scores for all networks
set.seed(6373)

# -------------------------------------------------------------------------
# library(foreach)
# library(doParallel)
# # 
# workers <- 1
# cl <- makeCluster(workers)
# # register the cluster for using foreach
# registerDoParallel(cl)

# -------------------------------------------------------------------------
# calculate z-scores?
calc.z <- T
null.reps <- 100

# -------------------------------------------------------------------------

d1 <- read.csv2("data/melbourne_all_interactions.csv")
load("data/melbourne_network_list.RData")
sp.info <- read.csv2("data/melbourne_species_info.csv")
site.distances <- read.csv("data/site.distances.csv")
sites <- names(net.list)
interaction.types <- names(net.list[[1]])
all.sp <- sort(unique(sp.info$species))

# -------------------------------------------------------------------------
# for each network, obtain a series of metrics

metrics.list <- list()
for(i.site in 1:length(sites)){
  metrics.list[[i.site]] <- list()
  for(i.type in 1:length(interaction.types)){
    metrics.list[[i.site]][[i.type]] <- list()
  } # for i.type
  names(metrics.list[[i.site]]) <- interaction.types
} # for i.site
names(metrics.list) <- sites

# i.site <- i.type <- 2
for(i.site in 1:length(sites)){
  for(i.type in 1:length(interaction.types)){
    
    my.metrics <- list()
    
    # my.graph <- graph_from_data_frame(net.list[[i.site]][[i.type]])
    my.net <- as.data.frame(net.list[[i.site]][[i.type]])
    names(my.net) <- c("node_to","node_from","link_value")
    
    num.to <- length(unique(my.net$node_to))
    num.from <- length(unique(my.net$node_from))
    
    if(nrow(my.net)>1 & num.to > 1 & num.from > 1){
      
      rich <- richness.edge.list(my.net)
      con <- connectance.edge.list(my.net)
      ld <- link.density.edge.list(my.net)
      
      deg <- NA
      try(deg <- degree.edge.list(my.net))
      
      mod.eigen <- NA
      try(mod.eigen <- modularity.edge.list(my.net,modularity.type = "leading_eigen"))
      
      mod.betweenness <- NA
      try(mod.betweenness <- modularity.edge.list(my.net,modularity.type = "edge_betweenness"))
      
      # mod.infomap <- NA
      # try(mod.infomap <- modularity.edge.list(my.net,modularity.type = "infomap"))
      
      nest <- 0
      if(con < 1){
        nest <- NA
        try(nest <- nestedness.edge.list(my.net))
      }
      
      cent.degree <- NA
      try(cent.degree <- centralization.edge.list(my.net,centrality.type = "degree"))
      
      cent.betweenness <- NA
      try(cent.betweenness <- centralization.edge.list(my.net, centrality.type = "betweenness"))
      
      cent.eigen <- NA
      try(cent.eigen <- centralization.edge.list(my.net, centrality.type = "eigen"))
      
      avg.overlap <- NA
      try(avg.overlap <- interaction.overlap.edge.list(my.net,method = "horn"))
      if(inherits(avg.overlap,"data.frame")){
        avg.overlap <- mean(avg.overlap$overlap)
      }
      
      metrics.list[[i.site]][[i.type]]$richness <- rich
      metrics.list[[i.site]][[i.type]]$connectance <- con
      metrics.list[[i.site]][[i.type]]$link_density <- ld
      metrics.list[[i.site]][[i.type]]$degree_heterogeneity <- deg
      metrics.list[[i.site]][[i.type]]$modularity_eigen <- mod.eigen
      metrics.list[[i.site]][[i.type]]$modularity_betweenness <- mod.betweenness
      metrics.list[[i.site]][[i.type]]$nestedness <- nest
      metrics.list[[i.site]][[i.type]]$centrality_degree <- cent.degree
      metrics.list[[i.site]][[i.type]]$centrality_betweenness <- cent.betweenness
      metrics.list[[i.site]][[i.type]]$centrality_eigen <- cent.eigen
      metrics.list[[i.site]][[i.type]]$interaction_overlap <- avg.overlap 
    }else{
      metrics.list[[i.site]][[i.type]]$richness <- NA
      metrics.list[[i.site]][[i.type]]$connectance <- NA
      metrics.list[[i.site]][[i.type]]$link_density <- NA
      metrics.list[[i.site]][[i.type]]$degree_heterogeneity <- NA
      metrics.list[[i.site]][[i.type]]$modularity_eigen <- NA
      metrics.list[[i.site]][[i.type]]$modularity_betweenness <- NA
      metrics.list[[i.site]][[i.type]]$nestedness <- NA
      metrics.list[[i.site]][[i.type]]$centrality_degree <- NA
      metrics.list[[i.site]][[i.type]]$centrality_betweenness <- NA
      metrics.list[[i.site]][[i.type]]$centrality_eigen <- NA
      metrics.list[[i.site]][[i.type]]$interaction_overlap <- NA 
    }
  }# i.type
}# i.site

# -------------------------------------------------------------------------
# repeat the process for the randomised networks
if(calc.z){
  
  metrics.null.list <- list()
  for(i.site in 1:length(sites)){
    metrics.null.list[[i.site]] <- list()
    for(i.type in 1:length(interaction.types)){
      metrics.null.list[[i.site]][[i.type]] <- list()
      for(i.rep in 1:null.reps){
        metrics.null.list[[i.site]][[i.type]][[i.rep]] <- list()
      }
      names(metrics.null.list[[i.site]][[i.type]]) <- as.character(1:null.reps)
    } # for i.type
    names(metrics.null.list[[i.site]]) <- interaction.types
  } # for i.site
  names(metrics.null.list) <- sites
  
  # i.site <- i.type <- 2
  for(i.site in 1:length(sites)){
    for(i.type in 1:length(interaction.types)){
      for(i.rep in 1:null.reps){
        
        my.metrics <- list()
        
        # my.graph <- graph_from_data_frame(net.list[[i.site]][[i.type]])
        my.net <- as.data.frame(net.list[[i.site]][[i.type]])
        names(my.net) <- c("node_to","node_from","link_value")
        num.to <- length(unique(my.net$node_to))
        num.from <- length(unique(my.net$node_from))
        
        if(nrow(my.net)>1 & num.to > 1 & num.from > 1){
          
          # this randomizes the network
          my.null <- reshuffle.edge.list(my.net)
          
          rich <- richness.edge.list(my.null)
          con <- connectance.edge.list(my.null)
          ld <- link.density.edge.list(my.null)
          
          deg <- NA
          try(deg <- degree.edge.list(my.null))
          
          mod.eigen <- NA
          try(mod.eigen <- modularity.edge.list(my.null,modularity.type = "leading_eigen"))
          
          mod.betweenness <- NA
          try(mod.betweenness <- modularity.edge.list(my.null,modularity.type = "edge_betweenness"))
          
          # mod.infomap <- NA
          # try(mod.infomap <- modularity.edge.list(my.null,modularity.type = "infomap"))
          
          nest <- 0
          if(con < 1){
            nest <- NA
            try(nest <- nestedness.edge.list(my.null))
          }
          
          cent.degree <- NA
          try(cent.degree <- centralization.edge.list(my.null,centrality.type = "degree"))
          
          cent.betweenness <- NA
          try(cent.betweenness <- centralization.edge.list(my.null, centrality.type = "betweenness"))
          
          cent.eigen <- NA
          try(cent.eigen <- centralization.edge.list(my.null, centrality.type = "eigen"))
          
          avg.overlap <- NA
          try(avg.overlap <- interaction.overlap.edge.list(my.null,method = "horn"))
          if(inherits(avg.overlap,"data.frame")){
            avg.overlap <- mean(avg.overlap$overlap)
          }
          
          metrics.null.list[[i.site]][[i.type]][[i.rep]]$richness <- rich
          metrics.null.list[[i.site]][[i.type]][[i.rep]]$connectance <- con
          metrics.null.list[[i.site]][[i.type]][[i.rep]]$link_density <- ld
          metrics.null.list[[i.site]][[i.type]][[i.rep]]$degree_heterogeneity <- deg
          metrics.null.list[[i.site]][[i.type]][[i.rep]]$modularity_eigen <- mod.eigen
          metrics.null.list[[i.site]][[i.type]][[i.rep]]$modularity_betweenness <- mod.betweenness
          metrics.null.list[[i.site]][[i.type]][[i.rep]]$nestedness <- nest
          metrics.null.list[[i.site]][[i.type]][[i.rep]]$centrality_degree <- cent.degree
          metrics.null.list[[i.site]][[i.type]][[i.rep]]$centrality_betweenness <- cent.betweenness
          metrics.null.list[[i.site]][[i.type]][[i.rep]]$centrality_eigen <- cent.eigen
          metrics.null.list[[i.site]][[i.type]][[i.rep]]$interaction_overlap <- avg.overlap 
        }else{
          metrics.null.list[[i.site]][[i.type]][[i.rep]]$richness <- NA
          metrics.null.list[[i.site]][[i.type]][[i.rep]]$connectance <- NA
          metrics.null.list[[i.site]][[i.type]][[i.rep]]$link_density <- NA
          metrics.null.list[[i.site]][[i.type]][[i.rep]]$degree_heterogeneity <- NA
          metrics.null.list[[i.site]][[i.type]][[i.rep]]$modularity_eigen <- NA
          metrics.null.list[[i.site]][[i.type]][[i.rep]]$modularity_betweenness <- NA
          metrics.null.list[[i.site]][[i.type]][[i.rep]]$nestedness <- NA
          metrics.null.list[[i.site]][[i.type]][[i.rep]]$centrality_degree <- NA
          metrics.null.list[[i.site]][[i.type]][[i.rep]]$centrality_betweenness <- NA
          metrics.null.list[[i.site]][[i.type]][[i.rep]]$centrality_eigen <- NA
          metrics.null.list[[i.site]][[i.type]][[i.rep]]$interaction_overlap <- NA 
        }# if-else valid network
      }# i.rep
    }# i.type
  }# i.site
}# if calc.z

# -------------------------------------------------------------------------
# transform nested list to dataframe
# there should be better ways, but this works

# i.site <- i.type <- 2
metrics.df.list <- list()
metric.names <- names(metrics.list[[1]][[1]])

for(i.site in 1:length(sites)){
  for(i.type in 1:length(interaction.types)){

        net.df <- data.frame(site = sites[i.site],
                             interaction_type = interaction.types[i.type])

        for(i.metric in metric.names){
          net.df[,ncol(net.df)+1] <- metrics.list[[i.site]][[i.type]][[i.metric]]
        }# for i.metric
        names(net.df)[3:ncol(net.df)] <- metric.names

        metrics.df.list[[length(metrics.df.list)+1]] <- net.df
      # }# if species in the network
      
  }# for i.type
}# for i.site

metrics.df <- bind_rows(metrics.df.list)
metrics.df.long <- pivot_longer(metrics.df,cols = 3:ncol(metrics.df),names_to = "metric",values_to = "value")

# i.site <- i.type <- i.rep <- 2
metrics.df.null.list <- list()
for(i.site in 1:length(sites)){
  for(i.type in 1:length(interaction.types)){
    for(i.rep in 1:null.reps){
      
    net.df <- data.frame(site = sites[i.site],
                         interaction_type = interaction.types[i.type],
                         rep = i.rep)

    for(i.metric in metric.names){
      net.df[,ncol(net.df)+1] <- metrics.null.list[[i.site]][[i.type]][[i.rep]][[i.metric]]
    }# for i.metric
    names(net.df)[4:ncol(net.df)] <- metric.names

    metrics.df.null.list[[length(metrics.df.null.list)+1]] <- net.df
    # }# if species in the network
    
    }# for i.rep
  }# for i.type
}# for i.site

metrics.null.df <- bind_rows(metrics.df.null.list)
metrics.null.df.long <- pivot_longer(metrics.null.df,cols = 4:ncol(metrics.null.df),names_to = "metric",values_to = "value")

# -------------------------------------------------------------------------

write.csv2(metrics.df.long,file = "results/melbourne_community_metrics.csv",row.names = F)
write.csv2(metrics.null.df.long,file = "results/melbourne_community_null_metrics.csv",row.names = F)


