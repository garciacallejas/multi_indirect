
library(tidyverse)
library(lmerTest)
library(DHARMa)
library(scales)
library(betareg)
library(colorblindr)
library(GGally)
library(ggfortify)
library(ggConvexHull)
library(patchwork)
# -------------------------------------------------------------------------

d1 <- read.csv2("data/melbourne_all_interactions.csv")
sp.info <- read.csv2("data/melbourne_species_info.csv")
sp.metrics <- read.csv2("results/melbourne_species_metrics.csv")

# -------------------------------------------------------------------------
# 1 - ubiquity-body size

# number of presences
freq.data <- d1 %>% 
  filter(functional.group != "Other") %>%
  group_by(species,functional.group, buggroup,size) %>%
  summarise(num.presences = n())
freq.by.site <- d1 %>% 
  filter(functional.group != "Other") %>%
  group_by(species,site,functional.group,buggroup) %>%
  summarise(num.presences = n(),mean.body.size = mean(size))
# m.freq <- glm(num.presences ~ size*functional.group, family=poisson, data=freq.data)

# per group
# herb.data <- subset(freq.data,functional.group == "Herbivore")
# m.herb.freq <- glm(num.presences ~ size, family=poisson, data=herb.data)
# 
# pol.data <- subset(freq.data,functional.group == "Pollinator")
# m.pol.freq <- glm(num.presences ~ size, family=poisson, data=pol.data)

# -------------------------------------------------------------------------
# body size - species metrics

sp.full.data <- left_join(freq.by.site,sp.metrics)
sp.body.size.data <- subset(sp.full.data, !is.na(mean.body.size)) %>% 
  select(site,species,functional.group,buggroup,metric,value,mean.body.size,num.presences)

sp.bs.wide <- pivot_wider(sp.body.size.data, names_from = "metric",values_from = "value")
sp.bs.wide.clean <- sp.bs.wide[complete.cases(sp.bs.wide),]

# -------------------------------------------------------------------------
# visualization
pd <- .2
size.plot <- ggplot(sp.body.size.data, aes(x = functional.group, y = mean.body.size)) + 
  geom_point(aes(color = functional.group), position = position_jitter(pd)) +
  geom_boxplot(aes(fill = functional.group), alpha = .4) + 
  theme(legend.position="none") +
  NULL
# size.plot

metrics.plot <- ggplot(sp.body.size.data, aes(x = functional.group, y = value)) + 
  geom_point(aes(color = functional.group), position = position_jitter(pd)) +
  geom_boxplot(aes(fill = functional.group), alpha = .4) + 
  facet_wrap(metric~., scales = "free_y") +
  theme(legend.position="none") +
  NULL
# metrics.plot
sp.bsd.2 <- sp.body.size.data
sp.bsd.2$apis <- sp.bsd.2$species == "Apis mellifera"
size.rel.plots <- ggplot(sp.bsd.2, aes(x = value, y = mean.body.size)) + 
  geom_point(aes(color = apis)) + 
  facet_wrap(metric~functional.group, scales = "free") +
  NULL
# size.rel.plots

# ggsave(filename = paste("results/images/body_size_distribution.pdf",sep=""),
#        plot = size.plot,
#        device = cairo_pdf,
#        width = 6, height = 4,dpi = 300)
# 
# ggsave(filename = paste("results/images/species_metrics_distribution.pdf",sep=""),
#        plot = metrics.plot,
#        device = cairo_pdf,
#        width = 10, height = 6,dpi = 300)

# ggsave(filename = paste("results/images/scatterplot_metrics_size.pdf",sep=""),
#        plot = size.rel.plots,
#        device = cairo_pdf,
#        width = 10, height = 8,dpi = 300)

# -------------------------------------------------------------------------
# how do size-metric plots look like?
my.metrics <- c("normalised.degree","species.strength",
                "weighted.betweenness","weighted.closeness","d")
plot.list <- list()
for(i.metric in 1:length(my.metrics)){
  my.column <- sym(my.metrics[i.metric])
  plot.list[[i.metric]] <- ggplot(sp.bs.wide.clean, aes(y = !!my.column,x = mean.body.size)) + 
    geom_point(aes(color = functional.group)) + 
    theme_bw() +
    ggtitle(my.metrics[i.metric]) +
    NULL
}
# patchwork::wrap_plots(plot.list,ncol = 2)


# -------------------------------------------------------------------------
# prepare data

# will use a beta regression, for outcomes [0,1]
# to avoid zeros and ones, transform (mentioned in the betareg documentation)
sp.bs.wide.clean$normalised.degree.transf <- (sp.bs.wide.clean$normalised.degree * (nrow(sp.bs.wide.clean)-1) + .5)/nrow(sp.bs.wide.clean)
sp.bs.wide.clean$buggroup.simplified <- sp.bs.wide.clean$buggroup
sp.bs.wide.clean$buggroup.simplified[which(sp.bs.wide.clean$buggroup.simplified %in% c("Cicadas","Flies","Parasitoid wasps","Sawflies"))] <- "Other"

# -------------------------------------------------------------------------
# relationship between the different metrics
# 1 - pca
pca.metrics <- sp.bs.wide.clean %>% 
  ungroup() %>%
  dplyr::select(species.strength:normalised.degree.transf) %>%
  # dplyr::select(-modularity.infomap) %>%
  # drop_na() %>%
  prcomp(scale = TRUE)

# nicer plot with ggfortify
pca.plot <- autoplot(pca.metrics, data = sp.bs.wide.clean, 
                     fill = 'buggroup.simplified', 
                     colour = "buggroup.simplified",
                     shape = 21,
                     # size = .8,
                     loadings = TRUE, 
                     loadings.colour = 'darkgrey',
                     loadings.label = TRUE, 
                     loadings.label.size = 3.5,
                     loadings.label.colour = "darkred") +
  theme_bw() + 
  # scale_fill_OkabeIto(name = "Interaction type", darken = .2) +
  # scale_colour_OkabeIto(name = "Interaction type", darken = .2) +
  # geom_convexhull(aes(color = buggroup.simplified),
  #                 alpha = 0.2) +
  NULL
# pca.plot

# correlations between pairs of metrics
metrics.pairs <- ggpairs(sp.bs.wide.clean, columns = c(7,8,10,12,13), 
                         ggplot2::aes(colour=functional.group))
# metrics.pairs

# -------------------------------------------------------------------------
pol.data <- subset(sp.bs.wide.clean, functional.group == "Pollinator")
herb.data <- subset(sp.bs.wide.clean, functional.group == "Herbivore")

# -------------------------------------------------------------------------
herb.degree <- betareg(normalised.degree.transf ~ scale(mean.body.size) + scale(num.presences) + buggroup.simplified,
                                     data = herb.data)
summary(herb.degree)

pol.degree <- betareg(normalised.degree.transf ~ scale(mean.body.size) + scale(num.presences) + buggroup.simplified,
                       data = pol.data)
summary(pol.degree)
# -------------------------------------------------------------------------
# these do not fit very well
herb.strength <- glm(species.strength ~ scale(mean.body.size) + scale(num.presences) + buggroup.simplified,
                       data = herb.data, family = Gamma(link = "log"))
summary(herb.strength)
DHARMa::plotResiduals(herb.strength)
DHARMa::testResiduals(herb.strength)

pol.strength <- glmer(species.strength ~ scale(mean.body.size) + scale(num.presences) + buggroup.simplified + (1|site),
                     data = pol.data, family = Gamma(link = "log"))
summary(pol.strength)
DHARMa::plotResiduals(pol.strength)
DHARMa::testResiduals(pol.strength)

# 
# pol.strength.1 <- glmer(size ~ species.strength + (1|site), 
#             data = pol.data, 
#             family = Gamma(link = "log"))
# herb.strength.1 <- glmer(size ~ species.strength + (1|site), 
#             data = herb.data, 
#             family = Gamma(link = "log"))
# 
# testResiduals(pol.strength.1)
# testResiduals(herb.strength.1)
# 
# testDispersion(pol.strength.1)
# testDispersion(herb.strength.1)
# 
# plotResiduals(pol.strength.1)
# plotResiduals(herb.strength.1)
# 
# pol.strength.1 <- glmer(size ~ species.strength + (1|site), 
#                         data = pol.data, 
#                         family = Gamma(link = "log"))
# herb.strength.1 <- glmer(size ~ species.strength + (1|site), 
#                          data = herb.data, 
#                          family = Gamma(link = "log"))
# 
# testResiduals(pol.strength.1)
# testResiduals(herb.strength.1)
# 
# testDispersion(pol.strength.1)
# testDispersion(herb.strength.1)
# 
# plotResiduals(pol.strength.1)
# plotResiduals(herb.strength.1)
# 
# # -------------------------------------------------------------------------
# 
# herb.full.1 <- glmer(size ~ scale(species.strength) + 
#                        scale(d) + 
#                        scale(normalised.degree) +
#                        scale(weighted.closeness) +
#                        scale(weighted.betweenness) + 
#                        (1|site), 
#                      data = herb.data, 
#                      family = Gamma(link = "log"))
# 
# plotResiduals(herb.full.1)
# 
# herb.full.2 <- glm(size ~ scale(species.strength) + 
#                        scale(d) + 
#                        scale(normalised.degree) +
#                        scale(weighted.closeness) +
#                        scale(weighted.betweenness),
#                      data = herb.data, 
#                      family = Gamma(link = "log"))
# 
# plotResiduals(herb.full.2)
# testResiduals(herb.full.2)
# 
# herb.full.3 <- lmer(size ~ scale(species.strength) + 
#                      scale(d) + 
#                      scale(normalised.degree) +
#                      scale(weighted.closeness) +
#                      scale(weighted.betweenness) + 
#                       (1|site), 
#                    data = herb.data)
# 
# plotResiduals(herb.full.3)
# testResiduals(herb.full.3)
