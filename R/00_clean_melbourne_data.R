
library(tidyverse)

# -------------------------------------------------------------------------

# tests of melbourne datasets
# 1 - original

# herb <- read.csv("../datasets/Melbourne_urban_parks/herb.final.csv")
# pol <- read.csv("../datasets/Melbourne_urban_parks/pol.final.csv")
# 
# # 2 - mata 2020
# full <- read.csv("../datasets/Melbourne_urban_parks/mata2020.csv")

# 3 - mata revised
full.2 <- read.csv2("../datasets/Melbourne_urban_parks/mata_revised.csv")
full.2$size <- as.numeric(full.2$size)

site.info <- read.csv("data/site.data.csv")
site.info$site_id <- paste("site_",1:nrow(site.info),sep="")

# add butterfly data
butterflies <- read.csv2("data/melbourne_butterflies_clean.csv")

# -------------------------------------------------------------------------

full.3 <- full.2[,c("site","vegspecies","buggroup","family","species","size","pol2","her2")]

# fix hoverflies
full.3$her2[full.3$family == "Syrphidae"] <- 0

# pol.her <- subset(full.2,pol2 == 1 & her2 == 1)
pol.only <- subset(full.3, pol2 == 1 & her2 == 0)
her.only <- subset(full.3, pol2 == 0 & her2 == 1)
none <- subset(full.3, pol2 == 0 & her2 == 0)

# -------------------------------------------------------------------------

full.3$functional.group <- "Other"
full.3$functional.group[full.3$pol2 == 1] <- "Pollinator"
full.3$functional.group[full.3$her2 == 1] <- "Herbivore"

full.3 <- full.3[,c("site","vegspecies","buggroup","family","species","size","functional.group")]

butter.2 <- butterflies[,c("site","vegspecies","butterfly")]
butter.2$functional.group <- "Pollinator"
butter.2$buggroup <- "Butterflies"

butter.2$family <- NA
butter.2$family[butter.2$butterfly == "Brown butterfly group"] <- "Nymphalidae"
butter.2$family[butter.2$butterfly == "Eurema smilax"] <- "Pieridae"
butter.2$family[butter.2$butterfly == "Hesperiidae group"] <- "Hesperiidae"
butter.2$family[butter.2$butterfly == "Junonia villida"] <- "Nymphalidae"
butter.2$family[butter.2$butterfly == "Little blue butterfly group"] <- "Lycaenidae"
butter.2$family[butter.2$butterfly == "Papilio anactus"] <- "Papilionidae"
butter.2$family[butter.2$butterfly == "Pieris rapae"] <- "Pieridae"
butter.2$family[butter.2$butterfly %in% c("Vanessa itea","Vanessa kershawi")] <- "Nymphalidae"

# from kirk et al.'s report
butter.2$wingspan <- NA
butter.2$wingspan[butter.2$butterfly == "Brown butterfly group"] <- 46.5
butter.2$wingspan[butter.2$butterfly == "Eurema smilax"] <- 31.5
butter.2$wingspan[butter.2$butterfly == "Hesperiidae group"] <- 16.5
butter.2$wingspan[butter.2$butterfly == "Junonia villida"] <- 41.5
butter.2$wingspan[butter.2$butterfly == "Little blue butterfly group"] <- 29
butter.2$wingspan[butter.2$butterfly == "Papilio anactus"] <- 69.5
butter.2$wingspan[butter.2$butterfly == "Pieris rapae"] <- 44
butter.2$wingspan[butter.2$butterfly == "Vanessa itea"] <- 50
butter.2$wingspan[butter.2$butterfly == "Vanessa kershawi"] <- 45

# TODO how to adjust size between butterflies and the rest of the insects?
butter.2$size <- butter.2$wingspan
butter.2$wingspan <- NULL

butter.2$species <- butter.2$butterfly
butter.2$butterfly <- NULL

butter.3 <- butter.2[,c("site","vegspecies","buggroup","family","species","size","functional.group")]

# -------------------------------------------------------------------------
full.4 <- bind_rows(full.3,butter.3)

# -------------------------------------------------------------------------
# store networks in a list
sites <- sort(unique(full.4$site))
interaction.types <- c("Pollinator","Herbivore")

# -------------------------------------------------------------------------
# lists of networks (edge lists)
# one list for herbivory networks
# and another for pollination networks

net.list <- list()
for(i.site in 1:length(sites)){
  net.list[[i.site]] <- list()
  for(i.type in 1:length(interaction.types)){
    net.list[[i.site]][[i.type]] <- NA
  }# i.type
  names(net.list[[i.site]]) <- interaction.types
}# i.site
names(net.list) <- sites

# -------------------------------------------------------------------------

# i.site <- i.type <- 2

for(i.site in 1:length(sites)){
  for(i.type in 1:length(interaction.types)){
    
    my.int <- subset(full.4,site == sites[i.site] & 
                       functional.group == interaction.types[i.type])
    
    net.list[[i.site]][[i.type]] <- my.int %>% 
      group_by(vegspecies,species) %>%
      summarise(freq = n())
    
    # my.graph <- graph_from_data_frame(int.freq)
  }# for i.type
}# i.site

# -------------------------------------------------------------------------
# species attributes
plant.sp <- data.frame(species = sort(unique(full.4$vegspecies)),
                       functional.group = "Plant",
                       size = NA)
animal.sp <- unique(full.4[,c("species","functional.group","size")])

all.sp <- bind_rows(plant.sp,animal.sp)
all.sp <- all.sp %>%
  group_by(functional.group) %>%
  mutate(sp_id = paste(functional.group,"_",row_number(),sep = ""))

# -------------------------------------------------------------------------
# pd <- .1
# size.plot <- ggplot(full.3,aes(x = functional.group,y = size)) + 
#   geom_point(position = position_jitter(pd)) +
#   geom_boxplot(alpha = .4) + 
#   NULL
# size.plot

# -------------------------------------------------------------------------

save(net.list,file="data/melbourne_network_list.RData")
write.csv2(all.sp,"data/melbourne_species_info.csv",row.names = F)
write.csv2(full.4,"data/melbourne_all_interactions.csv",row.names = F)
write.csv2(site.info,"data/site_info_id.csv",row.names = F)
