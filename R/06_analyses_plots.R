
# analyze the potential for apparent interactions 
# and other spatial patterns

# -------------------------------------------------------------------------

library(tidyverse)

# -------------------------------------------------------------------------

app.fac <- read.csv2("results/potential_apparent_facilitation.csv")
app.com <- read.csv2("results/potential_apparent_competition.csv")
site.distances <- read.csv("data/site.distances.csv")
sp.info <- read.csv2("data/melbourne_species_info.csv")
load("data/melbourne_network_list.RData")

# -------------------------------------------------------------------------
# spatial turnover
sp.turnover <- subset(sp.info, functional.group %in% c("Herbivore","Pollinator"))

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
    
    if(sp.turnover$functional.group[i.sp] == "Herbivore"){
      site.1.present <- my.sp %in% herb.site.1
      site.2.present <- my.sp %in% herb.site.2
    }else{
      site.1.present <- my.sp %in% pol.site.1
      site.2.present <- my.sp %in% pol.site.2
    }
    turnover.list[[length(turnover.list)+1]] <- data.frame(species = my.sp,
                                                           guild = sp.turnover$functional.group[i.sp],
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
# -------------------------------------------------------------------------

num.cooc.sp <- turnover.df %>%
  group_by(distance,guild) %>%
  summarise(num.cooccurring.sp = sum(cooccurrence))

dist.cooc.plot <- ggplot(num.cooc.sp, aes(x = distance, y = num.cooccurring.sp, group = guild)) + 
  geom_line(aes(colour = guild)) +
  theme_bw() +
  NULL
# dist.cooc.plot

cooc.sp <- subset(turnover.df,cooccurrence == T)
body.size.dist <- ggplot(cooc.sp,aes(x = distance, y = body.size)) + 
  stat_summary(aes(color = guild)) +
  scale_x_continuous(breaks = c(100,1000,2000,3000,4000,5000,6000)) +
  NULL
# body.size.dist

# -------------------------------------------------------------------------
# potential for apparent interactions

site.dist.clean <- site.distances[,c("siteA","siteB","neardist")]
names(site.dist.clean)[1:2] <- c("site_1","site_2")

app.fac.dist <- left_join(app.fac,site.dist.clean) %>% replace_na(list(neardist = 0))
app.com.dist <- left_join(app.com,site.dist.clean) %>% replace_na(list(neardist = 0))

app.all <- bind_rows(app.fac.dist,app.com.dist)

# -------------------------------------------------------------------------

app.dist.plot <- ggplot(app.all,aes(x = neardist, y = apparent_interaction_value)) +
  geom_point() + 
  facet_grid(interaction_type~.,scales = "free_y") +
  theme_bw() +
  NULL
# app.dist.plot

# -------------------------------------------------------------------------

ggsave(filename = "results/images/cooccurrence_distance.pdf",plot = dist.cooc.plot,
       device = cairo_pdf,
       width = 8,height = 5,dpi = 300)
ggsave(filename = "results/images/body_size_distance.pdf",plot = body.size.dist,
       device = cairo_pdf,
       width = 8,height = 5,dpi = 300)
ggsave(filename = "results/images/apparent_interaction_test.pdf",plot = app.dist.plot,
       device = cairo_pdf,
       width = 8,height = 5,dpi = 300)





