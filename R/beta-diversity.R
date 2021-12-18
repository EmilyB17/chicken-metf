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
t1 <- ps %>% 
  ps_filter(Treatment %in% c("Control", "T4")) %>% 
  # transform to relative abundance
  tax_transform("compositional", keep_counts = FALSE) %>% 
  dist_calc("bray")

# test beta dispersion
bd <- t1 %>% dist_bdisp(variables = "Treatment") %>% bdisp_get() # significance
# plot
#plot(bd$Treatment$model)

# test
mod1 <- t1 %>% 
  dist_permanova(
    seed = 123,
    variables = c("Treatment"),
    n_perms = 9999
  ) # significant

# plot
mod1 %>% 
  ord_calc(method = "auto") %>% 
  ord_plot(color = "Treatment") +
  stat_ellipse(aes(linetype = Treatment, color = Treatment)) +
  geom_text(aes(x = -0.8, y = -1.2), label = "adonis p value = 0.0026") +
  ggtitle("T4 versus control")

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
#plot(bd$Sample.Date$model)

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
  stat_ellipse(aes(linetype = Time, color = Time)) +
  geom_text(aes(x = -0.8, y = -1.6), label = "adonis p value = 0.0043") +
  ggtitle("All treatments versus time")

## ---- test 3: metformin treatment dose response ----

# get data
t1 <- ps %>% 
  ps_filter(!Treatment %in% c("Control", "Males")) %>% 
  # transform to relative abundance
  tax_transform("compositional", keep_counts = FALSE) %>% 
  dist_calc("bray")

# test beta dispersion
bd <- t1 %>% dist_bdisp(variables = "Treatment") %>% bdisp_get() # significance
# plot
#plot(bd$Treatment$model)

# test
mod4 <- t1 %>% 
  dist_permanova(
    seed = 123,
    variables = c("Treatment"), # significant
    n_perms = 9999
  )

# plot
mod4 %>% 
  ord_calc(method = "PCoA") %>% 
  ord_plot(color = "Treatment") +
  stat_ellipse(aes(linetype = Treatment, color = Treatment)) +
  geom_text(aes(x = -0.8, y = -1.4), label = "adonis p value = 2e-4") +
  ggtitle("Dose response")

