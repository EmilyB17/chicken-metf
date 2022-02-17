## FIGURE 1: taxa summary

# EVS 2/2022

require(microViz)
require(ggpubr)
require(tidyverse)
require(RColorBrewer)

load("data/ps-decontam-filtered-counts.RData")

# remove males
ps <- pscount %>% 
  ps_filter(Treatment != "Males")

## ---- barplot ----

ps %>% 
  tax_fix() %>% 
  comp_barplot(
    tax_level = "Phylum",
    label = NULL, # don't label each sample (too many!)
    n_taxa = length(get_taxa_unique(ps, "Phylum")), # plot all taxa
    bar_outline_color = NA
    ) +
  coord_flip()

# save
ggsave(filename = "R/figures/fig1-tax-summary.png", plot = last_plot(), dpi = 600)

## ---- top genera ----

# plot top 10 most abundant genera
psg <- ps %>% 
  tax_glom("Genus")

top <- prune_taxa(names(sort(taxa_sums(psg), TRUE))[1:10], psg)

top %>% 
  # recode metformin treatment
  ps_mutate(newTreat = recode(Treatment,
                              Control = "0 mg/kg",
                              T2 = "25 mg/kg",
                              T3 = "50 mg/kg", 
                              T4 = "75 mg/kg")) %>% 
  tax_fix() %>% 
  merge_samples(group = "newTreat") %>% 
  comp_barplot(tax_level = "Genus",
               n_taxa = 10,
               bar_width = 0.95,
               bar_outline_colour = "grey5"
  ) + 
  coord_flip()

# save
ggsave(filename = "R/figures/fig1-genus-taxsum.png", plot = last_plot(), dpi = 600)


## which ranks do these belong to?
tax_table(prune_taxa(names(sort(taxa_sums(psg), TRUE))[1:6], psg))

## ---- heatmap ----

# recode metformin for coloring 
ps <- ps %>% 
  # recode metformin treatment
  ps_mutate(Treatment = recode(Treatment,
                              Control = "0 mg/kg",
                              T2 = "25 mg/kg",
                              T3 = "50 mg/kg", 
                              T4 = "75 mg/kg")) 
# assign colors
###cols <- brewer.pal("Set2", n = 4) ## deuces RColorBrewerr
require(colorspace)
cols <- sequential_hcl(palette = "TealGrn", n = 4)
#names(cols) <- unique(samdat_tbl(ps)$Treatment)
names(cols) <- c("75 mg/kg", "50 mg/kg", "25 mg/kg", "0 mg/kg")


## for heatmap: view color palettes from colorspace: https://colorspace.r-forge.r-project.org/articles/colorspace.html

p <- ps %>% 
  tax_fix() %>% 
  tax_transform("clr", rank = "Phylum") %>% 
  #tax_filter(min_prevalence = 0.5) %>%  # keep only taxa in most samples 
  comp_heatmap(
    colors = heat_palette(palette = "Blue-Yellow", rev = TRUE), name = "CLR",
    tax_anno = taxAnnotation(
      Prevalence = anno_tax_prev(bar_width = 0.3, size = grid::unit(1, "cm"))
    ),
    sample_anno = sampleAnnotation(
      Treatment = anno_sample("Treatment"),
      col = list(Treatment = cols), border = FALSE
                             )
    )

# save - complex heatmap so save as pdf
pdf(file = "R/figures/fig1-heatmap-phyla.pdf")
p
dev.off()
