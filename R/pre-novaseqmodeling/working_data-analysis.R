
## Working script for data analysis 
# EVS 7/2021

# get data
require(tidyverse)
require(phyloseq)
require(microViz)
require(vegan)
require(ggpubr)
load("./data/ps-decontamed-samples.RData")

### ---- Prelim: Is there variability throughout time? ---

# Look at control group to remove influence of metformin
cont <- psdecon %>% 
  ps_filter(Treatment %in% c("Control", "Males")) %>% 
  ps_mutate(Time = factor(Time, ordered = TRUE, levels = c(1, 2, 3)))

ord <- ordinate(cont, method = "PCoA", distance = "bray")
plot_ordination(cont, ord, type = "samples", color = "Time")

### ---- Q1: Is there a difference between metformin treatment and control ----

# subset ps to time 3 with only Control and T4
sub <- psdecon %>% 
  ps_filter(Time == 3) %>% 
  ps_filter(Treatment %in% c("Control", "T4"))

# how many birds per group?
length(sub@sam_data$Bird[sub@sam_data$Treatment == "Control"]) # 9 birds
length(sub@sam_data$Bird[sub@sam_data$Treatment == "T4"]) # 8 birds

## ALPHA DIVERSITY

# get sample df
samp <- data.frame(sub@sam_data) %>% 
  rownames_to_column(var = "ID")

# get sample data and diversity metrics
alpha <- estimate_richness(sub, measures = c("Shannon", "Simpson", "Observed")) %>% 
  rownames_to_column(var = "ID") %>% 
  merge(samp, by = "ID") %>% 
  # make vertical
  pivot_longer(cols = c("Observed", "Shannon", "Simpson"), names_to = "Index", values_to = "value")

# plot diagnostics
ggdensity(alpha, x = "value", fill = "Treatment", facet.by = "Index", scales = "free")
ggqqplot(alpha, x = "value", color = "Treatment", facet.by = "Index", scales = "free",
         add = "qqline")

# fit a linear model
t.test(value ~ Treatment, data = alpha[alpha$Index == "Shannon",])
t.test(value ~ Treatment, data = alpha[alpha$Index == "Simpson",])

# visualize
ggboxplot(data = alpha, x = "Treatment", y = "value", fill = "Treatment", facet.by = "Index", scales = "free")

### PERMANOVA

## pre-process: are there any singletons
any(taxa_sums(sub) == 1)
ps <- prune_taxa(taxa_sums(sub) > 1, sub)

## transform to relative abundance
psT <- transform_sample_counts(ps, function(x)  x/sum(x))

# calculate Bray-Curtis distance
dis <- distance(psT, method = "bray", na.rm = TRUE)

# get sample data
sampdf <- data.frame(sample_data(psT))

# test with adonis
adonis(dis ~ Treatment, data = sampdf) # no significance

#### ORDINATE

# use pre-processed ps object above
ord <- ordinate(ps, method = "PCoA", distance = "bray")
plot_ordination(ps, ord, type = "samples", color = "Treatment")

## get to biom format for Songbird
require(biomformat)
otu <- t(data.frame(otu_table(sub)))

data <- make_biom(otu)
write_biom(data, biom_file = "./biom-Cont-T4-Time3.biom")

# get metadata for SongBird
write.table(samp, file = "./Control-T4-Time3.txt", sep = "\t", row.names = FALSE)
