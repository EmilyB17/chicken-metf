## Figure 3: Differential relative abundance

# EVS 2/2022

require(microViz)
require(tidyverse)
require(phyloseq)
require(ggpubr)

load("data/ps-decontam-filtered-counts.RData")

# load the model information
source("R/relative-abundance/log2-linearmod.R")

ps <- pscount %>% 
  ps_mutate(Time = as.factor(recode(Sample.Date, '8/8/19' = '40 weeks', '10/16/19' = '50 weeks', '12/20/19' = '60 weeks')),
            
            Treatment = recode_factor(Treatment,
                                      Control = "0 mg/kg",
                                      T2 = "25 mg/kg",
                                      T3 = "50 mg/kg",
                                      T4 = "75 mg/kg",
                                      .ordered = TRUE)) %>% 
  
  tax_fix()

# get colors
cols <- qualitative_hcl("Dark2", n = 3)

## ---- lollipop plot ----

# combine mod1, mod2, and mod3 (very few significant taxa)
df<- mod1$taxatree_stats %>% 
  # add identifying column
  mutate(model = "75 mg/kg vs 0 mg/kg") %>% 
  rbind(mod2$taxatree_stats %>% mutate(model = "60 weeks vs 40 weeks")) %>%  
  rbind(mod3$taxatree_stats %>% mutate(model = "75 mg/kg vs 25 mg/kg")) %>% 
  # get only sig genera
  filter(p.adj.Bon < 0.05) 

## get full taxonomic names for the 4 significant genera
names <- mod1 %>% ps_get() %>% tax_select(c("Acinetobacter",
                                            "Alloiococcus",
                                            "Cellulosilyticum",
                                            "Lachnospiraceae UCG-010",
                                            "[Ruminococcus] gnavus group")) %>% 
  ps_melt() %>% 
  select(Phylum, Class, Order, Family, Genus) %>% 
  distinct() %>% 
  # make name from Class & Genus
  mutate(taxa = paste(
    str_remove(Class, "C: "),
    str_remove(Genus, "G: "),
    sep = " "
  ))

## join
plot_dat <- df %>% full_join(names, by = c("taxon" = "Genus"))

p <- ggdotchart(plot_dat, x = "taxa", y = "estimate",
           color = "model", shape = "model",
           sorting = "descending",
           add = "segments", add.params = list(size = 2),
           dot.size = 5,
           rotate = TRUE,
           xlab = "",
           ylab = "log2 fold change") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  # remove legend label
  labs(color = "", shape = "")  +
  scale_color_manual(values = cols)

# change font
p1 <- ggpar(p, font.ytickslab = "italic")

# save
ggsave(filename = "R/figures/fig3-rel-abund.png", plot = p1, dpi = 600)

