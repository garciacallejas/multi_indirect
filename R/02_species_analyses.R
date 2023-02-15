
library(tidyverse)
library(lmerTest)
library(DHARMa)
library(scales)

# -------------------------------------------------------------------------

d1 <- read.csv2("data/melbourne_all_interactions.csv")
sp.info <- read.csv2("data/melbourne_species_info.csv")
sp.metrics <- read.csv2("results/melbourne_species_metrics.csv")

# -------------------------------------------------------------------------
# 1 - ubiquity-body size

# number of presences
freq.data <- d1 %>% 
  group_by(species,functional.group,size) %>%
  summarise(num.presences = length(unique(site)))

# m.freq <- glm(num.presences ~ size*functional.group, family=poisson, data=freq.data)

# per group
herb.data <- subset(freq.data,functional.group == "Herbivore")
m.herb.freq <- glm(num.presences ~ size, family=poisson, data=herb.data)

pol.data <- subset(freq.data,functional.group == "Pollinator")
m.pol.freq <- glm(num.presences ~ size, family=poisson, data=pol.data)

# -------------------------------------------------------------------------
# body size - species metrics

sp.full.data <- left_join(sp.metrics,sp.info)
sp.body.size.data <- subset(sp.full.data, !is.na(size))

sp.bs.wide <- pivot_wider(sp.body.size.data, names_from = "metric",values_from = "value")
sp.bs.wide.clean <- sp.bs.wide[complete.cases(sp.bs.wide),]

# -------------------------------------------------------------------------
# visualization
pd <- .2
size.plot <- ggplot(sp.body.size.data, aes(x = functional.group, y = size)) + 
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
size.rel.plots <- ggplot(sp.bsd.2, aes(x = value, y = size)) + 
  geom_point(aes(color = apis)) + 
  facet_wrap(metric~species_functional_group, scales = "free") +
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

ggsave(filename = paste("results/images/scatterplot_metrics_size.pdf",sep=""),
       plot = size.rel.plots,
       device = cairo_pdf,
       width = 10, height = 8,dpi = 300)
# -------------------------------------------------------------------------

pol.data <- subset(sp.bs.wide.clean, species_functional_group == "Pollinator")
herb.data <- subset(sp.bs.wide.clean, species_functional_group == "Herbivore")

pol.strength.1 <- glmer(size ~ species.strength + (1|site), 
            data = pol.data, 
            family = Gamma(link = "log"))
herb.strength.1 <- glmer(size ~ species.strength + (1|site), 
            data = herb.data, 
            family = Gamma(link = "log"))

testResiduals(pol.strength.1)
testResiduals(herb.strength.1)

testDispersion(pol.strength.1)
testDispersion(herb.strength.1)

plotResiduals(pol.strength.1)
plotResiduals(herb.strength.1)

pol.strength.1 <- glmer(size ~ species.strength + (1|site), 
                        data = pol.data, 
                        family = Gamma(link = "log"))
herb.strength.1 <- glmer(size ~ species.strength + (1|site), 
                         data = herb.data, 
                         family = Gamma(link = "log"))

testResiduals(pol.strength.1)
testResiduals(herb.strength.1)

testDispersion(pol.strength.1)
testDispersion(herb.strength.1)

plotResiduals(pol.strength.1)
plotResiduals(herb.strength.1)

# -------------------------------------------------------------------------

herb.full.1 <- glmer(size ~ scale(species.strength) + 
                       scale(d) + 
                       scale(normalised.degree) +
                       scale(weighted.closeness) +
                       scale(weighted.betweenness) + 
                       (1|site), 
                     data = herb.data, 
                     family = Gamma(link = "log"))

plotResiduals(herb.full.1)

herb.full.2 <- glm(size ~ scale(species.strength) + 
                       scale(d) + 
                       scale(normalised.degree) +
                       scale(weighted.closeness) +
                       scale(weighted.betweenness),
                     data = herb.data, 
                     family = Gamma(link = "log"))

plotResiduals(herb.full.2)
testResiduals(herb.full.2)

herb.full.3 <- lmer(size ~ scale(species.strength) + 
                     scale(d) + 
                     scale(normalised.degree) +
                     scale(weighted.closeness) +
                     scale(weighted.betweenness) + 
                      (1|site), 
                   data = herb.data)

plotResiduals(herb.full.3)
testResiduals(herb.full.3)
