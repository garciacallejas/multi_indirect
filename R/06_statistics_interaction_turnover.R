
# analyze the potential for apparent interactions 
# and other spatial patterns

# -------------------------------------------------------------------------

library(tidyverse)
library(vegan)
library(betareg)
library(DHARMa)
library(scales)
library(performance)

# -------------------------------------------------------------------------
all.interactions <- read.csv2("data/melbourne_all_interactions.csv")
app.fac <- read.csv2("results/potential_apparent_facilitation.csv")
app.com <- read.csv2("results/potential_apparent_competition.csv")
site.distances <- read.csv("data/site.distances.csv")
sp.info <- read.csv2("data/melbourne_species_info.csv")
load("data/melbourne_network_list.RData")
dissimilarity.metrics <- read.csv2("results/dissimilarity_between_sites.csv")
site.info <- read.csv("data/site.data.csv")
site.info$aream2 <- as.numeric(gsub(",","",site.info$aream2))

# -------------------------------------------------------------------------
# metrics across pairs of sites:
# 1 - how similar are two sites in habitat composition?
# 2 - what is their combined area? and the difference in areas?

site.pairs <- site.distances[,c("siteA","siteB","neardist")] 
site.pairs$habitat.similarity <- NA
site.pairs$combined.area <- NA
site.pairs$diff.area <- NA

for(i.pair in 1:nrow(site.pairs)){
  site1 <- site.pairs$siteA[i.pair]
  site2 <- site.pairs$siteB[i.pair]
  
  site.pairs$combined.area[i.pair] <- site.info$aream2[site.info$site == site1] + 
    site.info$aream2[site.info$site == site2]
  site.pairs$diff.area[i.pair] <- abs(site.info$aream2[site.info$site == site1] - 
                     site.info$aream2[site.info$site == site2])
  lawn <- site.info$lawn[site.info$site == site1] * site.info$lawn[site.info$site == site2]
  tree <- site.info$tree[site.info$site == site1] * site.info$tree[site.info$site == site2]
  grassland <- site.info$grassland[site.info$site == site1] * site.info$grassland[site.info$site == site2]
  midstorey <- site.info$midstorey[site.info$site == site1] * site.info$midstorey[site.info$site == site2]
  
  site.pairs$habitat.similarity[i.pair] <- lawn+tree+grassland+midstorey
  
}
# -------------------------------------------------------------------------
# plant composition of each site
all.plants <- all.interactions %>% group_by(site,vegspecies) %>% summarise(presence = n())
all.plants$presence <- 1
all.plants2 <- expand_grid(site = unique(all.plants$site),vegspecies = unique(all.plants$vegspecies))
all.plants2 <- left_join(all.plants2,all.plants[,c("site","vegspecies","presence")])
all.plants2$presence[which(is.na(all.plants2$presence))] <- 0

plant.community <- all.plants2 %>% pivot_wider(names_from = vegspecies,values_from = presence)
plant.comm.matrix <- as.matrix(plant.community[,2:ncol(plant.community)])
rownames(plant.comm.matrix) <- plant.community$site

site.plant.diss.wide <- vegan::betadiver(plant.comm.matrix,method = "w")
site.plant.diss.df <- as.data.frame(as.matrix(site.plant.diss.wide))
site.plant.diss.df$siteA <- rownames(site.plant.diss.df) 
site.plant.diss.long <- site.plant.diss.df %>% pivot_longer(cols = 1:(ncol(site.plant.diss.df)-1),
                                                            names_to = "siteB",
                                                            values_to = "plant_beta_div")

site.pairs <- left_join(site.pairs,site.plant.diss.long)

# -------------------------------------------------------------------------
# add dissimilarity metrics
diss.wide <- dissimilarity.metrics %>% pivot_wider(names_from = metric,values_from = value)
diss.wide.herb <- subset(diss.wide,interaction == "herbivory")
diss.wide.pol <- subset(diss.wide,interaction == "pollination")

site.pairs.herb <- left_join(diss.wide.herb,site.pairs)
site.pairs.pol <- left_join(diss.wide.pol,site.pairs)

# will use a beta regression, for outcomes [0,1]
# to avoid zeros and ones, transform (mentioned in the betareg documentation)
site.pairs.herb$interaction_turnover_transf <- (site.pairs.herb$interaction_turnover * (nrow(site.pairs.herb)-1) + .5)/nrow(site.pairs.herb)
site.pairs.pol$interaction_turnover_transf <- (site.pairs.pol$interaction_turnover * (nrow(site.pairs.pol)-1) + .5)/nrow(site.pairs.pol)
site.pairs.herb$species_turnover_transf <- (site.pairs.herb$species_turnover * (nrow(site.pairs.herb)-1) + .5)/nrow(site.pairs.herb)
site.pairs.pol$species_turnover_transf <- (site.pairs.pol$species_turnover * (nrow(site.pairs.pol)-1) + .5)/nrow(site.pairs.pol)

site.pairs.herb$habitat.similarity <- as.factor(site.pairs.herb$habitat.similarity)
site.pairs.pol$habitat.similarity <- as.factor(site.pairs.pol$habitat.similarity)

# -------------------------------------------------------------------------
# model
# apparently in a beta regression residuals do not need to be normally distributed?
# Espinheira, P. L., Ferrari, S. L., & Cribari-Neto, F. (2008a). Influence diagnostics in beta regression. Computational Statistics & Data Analysis, 52 (9), 4417–4431.
# Espinheira, P. L., Ferrari, S. L., & Cribari-Neto, F. (2008b). On beta regression residuals. Journal of Applied Statistics, 35(4), 407–419.
# Espinheira, P. L., Santos, E. G., & Cribari-Neto, F. (2017). On nonlinear beta regression residuals. Biometrical Journal, 59(3), 445–461.
# Ferrari, S. L., & Cribari-Neto, F. (2004). Beta regression for modelling rates and proportions. Journal of Applied Statistics, 31(7), 799–815.
# Pereira, G. H. (2017). On quantile residuals in beta regression. Communications in Statistics - Simulation and Computation, 1–15.

herb.interaction.turnover <- betareg(interaction_turnover_transf ~ scale(neardist) + habitat.similarity + scale(combined.area) + plant_beta_div,
                                     data = site.pairs.herb)
summary(herb.interaction.turnover)
performance::check_model(herb.interaction.turnover)
performance::model_performance(herb.interaction.turnover)

pol.interaction.turnover <- betareg(interaction_turnover_transf ~ scale(neardist) + habitat.similarity + scale(combined.area) + plant_beta_div,
                                     data = site.pairs.pol)
summary(pol.interaction.turnover)
performance::check_model(pol.interaction.turnover)
performance::model_performance(pol.interaction.turnover)

# -------------------------------------------------------------------------
# this below does not make much sense since the outcome is partly given by turnover in plant species alone
# TODO I should calculate turnover in pollinator/herbivore composition independently
# herb.sp.turnover <- betareg(species_turnover_transf ~ scale(neardist) + habitat.similarity + scale(combined.area) + plant_beta_div,
#                                      data = site.pairs.herb)
# summary(herb.sp.turnover)
# 
# pol.sp.turnover <- betareg(species_turnover_transf ~ scale(neardist) + habitat.similarity + scale(combined.area) + plant_beta_div,
#                                     data = site.pairs.pol)
# summary(pol.sp.turnover)
