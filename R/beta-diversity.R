## beta diversity

# Emily Van Syoc 12/2021

require(tidyverse)
require(phyloseq)
require(ggpubr)
require(vegan)
require(microViz)
require(pairwiseAdonis)
require(ggordiplots)

load("./data/ps-decontam-filtered-counts.RData")
ps <- pscount %>% 
  ps_mutate(Time = as.factor(recode(Sample.Date, '8/8/19' = 1, '10/16/19' = 2, '12/20/19' = 3)))


### ---- test 1: difference between control and T4 ----

# get data
t1 <- pscount %>% 
  ps_filter(Treatment %in% c("Control", "T4")) %>% 
  # transform to relative abundance
  tax_transform("compositional", keep_counts = FALSE) %>% 
  dist_calc("bray")

# test beta dispersion
bd <- t1 %>% dist_bdisp(variables = "Treatment") %>% bdisp_get() # significance
# plot
plot(bd$Treatment$model)

# test
mod1 <- t1 %>% 
  dist_permanova(
    seed = 123,
    variables = c("Treatment"),
    n_perms = 9999
  )

# plot
mod1 %>% 
  ord_calc(method = "PCoA") %>% 
  ord_plot(color = "Treatment") +
  stat_ellipse(aes(linetype = Treatment, color = Treatment))

## ---- test 2: difference in metformin treatment over time ----

# get data
t2 <- ps %>% 
  ps_filter(!Treatment %in% c("Control", "Males")) %>% 
  # transform to relative abundance
  tax_transform("compositional", keep_counts = FALSE) %>% 
  dist_calc("bray")

# test beta dispersion
bd <- t2 %>% dist_bdisp(variables = "Sample.Date") %>% bdisp_get() # not significant
# plot
plot(bd$Sample.Date$model)

# test
mod2 <- t2 %>% 
  dist_permanova(
    seed = 123,
    variables = c("Time"), # significant
    n_perms = 9999
  )

# plot
mod2 %>% 
  ord_calc(method = "PCoA") %>% 
  ord_plot(color = "Time") +
  stat_ellipse(aes(linetype = Time, color = Time))

## ---- test 3: test interaction between treatment & time ----

# get data
t1 <- pscount %>% 
  ps_filter(Treatment %in% c("Control", "T4")) %>% 
  # transform to relative abundance
  tax_transform("compositional", keep_counts = FALSE) %>% 
  dist_calc("bray")

# test beta dispersion
bd <- t1 %>% dist_bdisp(variables = "Treatment") %>% bdisp_get() # significance
# plot
plot(bd$Treatment$model)

# test
mod1 <- t1 %>% 
  dist_permanova(
    seed = 123,
    variables = c("Treatment"),
    n_perms = 9999
  )

# plot
mod1 %>% 
  ord_calc(method = "PCoA") %>% 
  ord_plot(color = "Treatment") +
  stat_ellipse(aes(linetype = Treatment, color = Treatment))

## ---- test 2: difference in metformin treatment over time ----

# get data
t3 <- ps %>% 
  # transform to relative abundance
  tax_transform("compositional", keep_counts = FALSE) %>% 
  dist_calc("bray")

# test beta dispersion
bd <- t3 %>% dist_bdisp(variables = c("Treatment", "Time")) %>% bdisp_get() # not significant
# plot
plot(bd$Treatment$model)

# test
mod3 <- t3 %>% 
  dist_permanova(
    seed = 123,
    variables = "Time + Treatment + Time*Treatment", # interaction is significant
    n_perms = 9999
  )

# plot
mod3 %>% 
  ord_calc(method = "PCoA") %>% 
  ord_plot(color = "Treatment", shape = "Time") 
