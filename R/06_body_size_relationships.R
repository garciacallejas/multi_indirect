
library(tidyverse)

# -------------------------------------------------------------------------

d1 <- read.csv2("data/melbourne_all_interactions.csv")
load("data/melbourne_network_list.RData")
site.distances <- read.csv("data/site.distances.csv")

sp.info <- read.csv2("data/melbourne_species_info.csv")
sp.info$species_functional_group <- sp.info$functional.group
sp.info$functional.group <- NULL
sp.metrics <- read.csv2("results/melbourne_species_metrics.csv")

# -------------------------------------------------------------------------
# spatial turnover/cooccurrences
sp.turnover <- subset(sp.info, species_functional_group %in% c("Herbivore","Pollinator"))

turnover.list <- list()
for(i.pair in 1:nrow(site.distances)){
  
  herb.site.1 <- unique(net.list[[site.distances$siteA[i.pair]]][["Herbivore"]]$species)
  herb.site.2 <- unique(net.list[[site.distances$siteB[i.pair]]][["Herbivore"]]$species)
  
  pol.site.1 <- unique(net.list[[site.distances$siteA[i.pair]]][["Pollinator"]]$species)
  pol.site.2 <- unique(net.list[[site.distances$siteB[i.pair]]][["Pollinator"]]$species)
  
  for(i.sp in 1:nrow(sp.turnover)){
    
    my.sp <- sp.turnover$species[i.sp]
    site.1.present <- F
    site.2.present <- F
    
    if(sp.turnover$species_functional_group[i.sp] == "Herbivore"){
      site.1.present <- my.sp %in% herb.site.1
      site.2.present <- my.sp %in% herb.site.2
    }else{
      site.1.present <- my.sp %in% pol.site.1
      site.2.present <- my.sp %in% pol.site.2
    }
    turnover.list[[length(turnover.list)+1]] <- data.frame(species = my.sp,
                                                           guild = sp.turnover$species_functional_group[i.sp],
                                                           body.size = sp.turnover$size[i.sp],
                                                           site1 = site.distances$siteA[i.pair],
                                                           site2 = site.distances$siteB[i.pair],
                                                           site1.presence = site.1.present,
                                                           site2.presence = site.2.present,
                                                           cooccurrence = site.1.present & site.2.present,
                                                           distance = site.distances$neardist[i.pair])
  }# for i.sp
}# for i.pair
turnover.df <- bind_rows(turnover.list)
turnover.df <- turnover.df %>%
  group_by(guild) %>%
  mutate(body.size.interval = cut_interval(body.size, n = 5))
cooc.distances.df <- turnover.df %>% group_by(guild,distance,body.size.interval) %>%
  summarise(num.cooc = sum(cooccurrence))
# -------------------------------------------------------------------------
cooc.body.size.distance.plot <- ggplot(cooc.distances.df, aes(x = distance, y = num.cooc, group = body.size.interval)) +
  geom_line(aes(color = body.size.interval)) + 
  facet_grid(guild~.,scales = "free_y") +
  NULL
# cooc.body.size.distance.plot

# -------------------------------------------------------------------------
avg.metrics <- sp.metrics %>% group_by(species,species_functional_group,metric) %>%
  summarise(avg.value = mean(value))

metrics.body.size <- left_join(avg.metrics, 
                               sp.info[,c("species","species_functional_group","size")]) %>%
  filter(species_functional_group != "Plant")

# -------------------------------------------------------------------------

plot.metrics.body.size <- ggplot(metrics.body.size,aes(x = size, y = avg.value)) + 
  geom_point(aes(fill = species_functional_group), shape = 21) +
  facet_wrap(metric~species_functional_group, scales = "free",drop = T) +
  theme(legend.position = "none") +
  NULL
# plot.metrics.body.size


