
library(tidyverse)
# -------------------------------------------------------------------------

d1 <- read.csv("../datasets/Melbourne_urban_parks/mata2020.csv")

d2 <- d1[,c("site","date","habitat","po","vegspecies","buggroup","order","family","genus","species","size","pol","her")]

# -------------------------------------------------------------------------
# subset species that are classified as both pollinator and herbivore
# there may be species/taxa that are neither pollinator nor herbivore

insect.sp <- unique(d2[,c("buggroup","order","family","genus","species","size","pol","her")])
insect.sp$category <- NA

to.classify <- subset(insect.sp, pol == her)

# -------------------------------------------------------------------------
# diptera
# classification from wikipedia and online search
# diptera.families <- sort(unique(to.classify$family[which(to.classify$order == "Diptera")]))

insect.sp$category[which(insect.sp$family == "Asilidae")] <- "Predator"
insect.sp$category[which(insect.sp$family == "Chloropidae")] <- "Other"
insect.sp$category[which(insect.sp$family == "Dolichopodidae")] <- "Predator"

# doubtful
insect.sp$category[which(insect.sp$family == "Drosophilidae")] <- "Other"

insect.sp$category[which(insect.sp$family == "Empididae")] <- "Predator"
insect.sp$category[which(insect.sp$family == "Ephydridae")] <- "Other"
insect.sp$category[which(insect.sp$family == "Heleomyzidae")] <- "Other"
insect.sp$category[which(insect.sp$family == "Lauxaniidae")] <- "Other"
insect.sp$category[which(insect.sp$family == "Lonchopteridae")] <- "Other"
insect.sp$category[which(insect.sp$family == "Phoridae")] <- "Pollinator"
insect.sp$category[which(insect.sp$family == "Pyrgotidae")] <- "Other"
insect.sp$category[which(insect.sp$family == "Scenopinidae")] <- "Pollinator"
insect.sp$category[which(insect.sp$family == "Sepsidae")] <- "Other"
insect.sp$category[which(insect.sp$family == "Syrphidae")] <- "Pollinator"
insect.sp$category[which(insect.sp$family == "Therevidae")] <- "Pollinator"

# -------------------------------------------------------------------------
# coleoptera.families <- sort(unique(to.classify$family[which(to.classify$order == "Coleoptera")]))

insect.sp$category[which(insect.sp$family == "Buprestidae")] <- "Herbivore"
insect.sp$category[which(insect.sp$family == "Cantharidae")] <- "Pollinator"
insect.sp$category[which(insect.sp$family == "Carabidae")] <- "Predator"
insect.sp$category[which(insect.sp$family == "Cleridae")] <- "Predator"
insect.sp$category[which(insect.sp$family == "Coccinellidae")] <- "Predator"
insect.sp$category[which(insect.sp$family == "Dermestidae")] <- "Pollinator"
insect.sp$category[which(insect.sp$family == "Hydrophilidae")] <- "Herbivore"
insect.sp$category[which(insect.sp$family == "Mordellidae")] <- "Pollinator"
insect.sp$category[which(insect.sp$family == "Nitidulidae")] <- "Herbivore"
insect.sp$category[which(insect.sp$family == "Phalacridae")] <- "Pollinator"
insect.sp$category[which(insect.sp$family == "Ptinidae")] <- "Other" # adults do not feed?
insect.sp$category[which(insect.sp$family == "Staphylinidae")] <- "Predator"

# -------------------------------------------------------------------------
# hemiptera.families <- sort(unique(to.classify$family[which(to.classify$order == "Hemiptera")]))

insect.sp$category[which(insect.sp$family == "Geocoridae")] <- "Predator"
insect.sp$category[which(insect.sp$family == "Lygaeidae")] <- "Herbivore"
insect.sp$category[which(insect.sp$family == "Nabidae")] <- "Predator" 
insect.sp$category[which(insect.sp$family == "Pentatomidae")] <- "Herbivore"
insect.sp$category[which(insect.sp$family == "Reduviidae")] <- "Predator"

# -------------------------------------------------------------------------
# hymenoptera.families <- sort(unique(to.classify$family[which(to.classify$order == "Hymenoptera")]))

insect.sp$category[which(insect.sp$family == "Agaonidae")] <- "Pollinator"
insect.sp$category[which(insect.sp$family == "Apidae")] <- "Pollinator"
insect.sp$category[which(insect.sp$family == "Bethylidae")] <- "Predator" 
insect.sp$category[which(insect.sp$family == "Braconidae")] <- "Other" # unknown
insect.sp$category[which(insect.sp$family == "Chalcididae")] <- "Other" # unknown
insect.sp$category[which(insect.sp$family == "Colletidae")] <- "Pollinator"
insect.sp$category[which(insect.sp$family == "Crabronidae")] <- "Pollinator"
insect.sp$category[which(insect.sp$family == "Diapriidae")] <- "Other" # unknown
insect.sp$category[which(insect.sp$family == "Encyrtidae")] <- "Other" # unknown
insect.sp$category[which(insect.sp$family == "Formicidae")] <- "Predator" # stretch?
insect.sp$category[which(insect.sp$family == "Gasteruptiidae")] <- "Pollinator"
# typo
insect.sp$family[which(insect.sp$family == "Halticidae")] <- "Halictidae"
insect.sp$category[which(insect.sp$family == "Halictidae")] <- "Pollinator"
insect.sp$category[which(insect.sp$family == "Ichneumonidae")] <- "Pollinator"
insect.sp$category[which(insect.sp$family == "Mutillidae")] <- "Pollinator"
insect.sp$category[which(insect.sp$family == "Pompilidae")] <- "Pollinator"
insect.sp$category[which(insect.sp$family == "Tiphiidae")] <- "Pollinator"
insect.sp$category[which(insect.sp$family == "Torymidae")] <- "Herbivore" 
insect.sp$category[which(insect.sp$family == "Vespidae")] <- "Predator"


