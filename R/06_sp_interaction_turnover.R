
library(tidyverse)
library(igraph)
library(betalink)

# -------------------------------------------------------------------------

load("data/melbourne_network_list.RData")
site.distances <- read.csv("data/site.distances.csv")
sp.info <- read.csv2("data/melbourne_species_info.csv")

# -------------------------------------------------------------------------
# different lists for herbivore networks and for pollinator networks
# networks must be igraph objects, I think
herb.list <- list()
pol.list <- list()

for(i.site in 1:length(net.list)){
  i.herb <- net.list[[i.site]][["Herbivore"]]
  i.pol <- net.list[[i.site]][["Pollinator"]]
  
  herb.list[[i.site]] <- igraph::graph_from_data_frame(i.herb)
  pol.list[[i.site]] <- igraph::graph_from_data_frame(i.pol)
}
names(herb.list) <- names(net.list)
names(pol.list) <- names(net.list)
# -------------------------------------------------------------------------
# s: dissimilarity in composition
# wn: dissimilarity in interactions (all interactions)
# os: dissimilarity in interactions between common species to both networks
# st: dissimilarity in interactions due to species turnover

herb.dissim <- betalink::network_betadiversity(herb.list)
pol.dissim <- betalink::network_betadiversity(pol.list)

names(herb.dissim)[1:2] <- c("siteA","siteB")
names(pol.dissim)[1:2] <- c("siteA","siteB")

herb.dissim.dist <- left_join(herb.dissim,site.distances[,c("siteA","siteB","neardist")]) %>%
  select(siteA,siteB,neardist,S,WN)
pol.dissim.dist <- left_join(pol.dissim,site.distances[,c("siteA","siteB","neardist")]) %>%
  select(siteA,siteB,neardist,S,WN)
names(herb.dissim.dist)[4:5] <- c("species_turnover","interaction_turnover")
names(pol.dissim.dist)[4:5] <- c("species_turnover","interaction_turnover")

herb.dd.long <- pivot_longer(herb.dissim.dist,species_turnover:interaction_turnover,names_to = "metric")
herb.dd.long$interaction <- "herbivory"
pol.dd.long <- pivot_longer(pol.dissim.dist,species_turnover:interaction_turnover,names_to = "metric")
pol.dd.long$interaction <- "pollination"

all.dd.long <- bind_rows(herb.dd.long,pol.dd.long)

# -------------------------------------------------------------------------

dissim.dist.plot <- ggplot(all.dd.long, aes(x = neardist, y = value, group = metric)) +
  geom_line(aes(linetype = metric)) +
  facet_grid(interaction~.) +
  NULL
# dissim.dist.plot

# -------------------------------------------------------------------------

# ggsave(filename = "results/images/turnover_distance.pdf",plot = dissim.dist.plot,
#        device = cairo_pdf,
#        width = 8,height = 5,dpi = 300)

