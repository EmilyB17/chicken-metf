## rel abundance - log2 linear model 

require(microViz)
require(tidyverse)
require(phyloseq)

load("data/ps-decontam-filtered-counts.RData")
ps <- pscount %>% 
  ps_mutate(Time = as.factor(recode(Sample.Date, '8/8/19' = '40 weeks', '10/16/19' = '50 weeks', '12/20/19' = '60 weeks'))) %>% 
  tax_fix()

## ---- test1: control v t4 ----

# subset 
phylo <- ps %>% 
  ps_filter(Treatment %in% c("Control", "T4")) 


# model on Genus
lm_models <- phylo %>% 
  tax_fix() %>% 
  tax_prepend_ranks() %>% 
  tax_transform("compositional", rank = "Genus") %>% 
  tax_transform(
    trans = "log2", 
    # chain allows for multiple transformations (TSS then log2)
    chain = TRUE, 
    # for log, replace 0's with half of min value
    zero_replace = "halfmin"
  ) %>% 
  taxatree_models(
    type = lm, 
    ranks = "Genus",
    variables = c("Treatment")
  )

# get stats
lm_stats <- taxatree_models2stats(lm_models)
lm_stats$taxatree_stats

# adjust for multiple comparisons
lm_stats <- taxatree_stats_p_adjust(
  data = lm_stats, method = "BH", grouping = "rank"
)

# which are significant?
lm_stats$taxatree_stats %>% filter(p.adj.BH.rank < 0.05)

# what are the effect sizes?
hist(lm_stats$taxatree_stats$estimate)

# which taxa are changing?
unique(lm_stats$taxatree_stats$taxon[lm_stats$taxatree_stats$p.adj.BH.rank < 0.05])

# save for comparison
mod1 <- lm_stats

## ---- test2: difference between metformin over time ----

# only pairwise comparisons, so use times 1 and 3
# get data
phylo <- ps %>% 
  ps_filter(!Treatment %in% c("Control", "Males")) %>% 
  ps_filter(Time %in% c("40 weeks", "60 weeks")) %>% 
  tax_fix()

# model on Genus
lm_models <- phylo %>% 
  tax_fix() %>% 
  tax_prepend_ranks() %>% 
  tax_transform("compositional", rank = "Genus") %>% 
  tax_transform(
    trans = "log2", 
    # chain allows for multiple transformations (TSS then log2)
    chain = TRUE, 
    # for log, replace 0's with half of min value
    zero_replace = "halfmin"
  ) %>% 
  taxatree_models(
    type = lm, 
    ranks = "Genus",
    variables = c("Time")
  )

# get stats
lm_stats <- taxatree_models2stats(lm_models)
lm_stats$taxatree_stats

# adjust for multiple comparisons
lm_stats <- taxatree_stats_p_adjust(
  data = lm_stats, method = "BH", grouping = "rank"
)

# which are significant?
lm_stats$taxatree_stats %>% filter(p.adj.BH.rank < 0.05)

# what are the effect sizes?
hist(lm_stats$taxatree_stats$estimate)

# which taxa are changing?
unique(lm_stats$taxatree_stats$taxon[lm_stats$taxatree_stats$p.adj.BH.rank < 0.05])

# save for later
mod2 <- lm_stats

### ---- sanity check: control over time ----

# get data
phylo <- ps %>% 
  ps_filter(Treatment %in% "Control") %>% 
  ps_filter(Time %in% c("40 weeks", "60 weeks")) %>% 
  tax_fix()

# model on Genus
lm_models <- phylo %>% 
  tax_fix() %>% 
  tax_prepend_ranks() %>% 
  tax_transform("compositional", rank = "Genus") %>% 
  tax_transform(
    trans = "log2", 
    # chain allows for multiple transformations (TSS then log2)
    chain = TRUE, 
    # for log, replace 0's with half of min value
    zero_replace = "halfmin"
  ) %>% 
  taxatree_models(
    type = lm, 
    ranks = "Genus",
    variables = c("Time")
  )

# get stats
lm_stats <- taxatree_models2stats(lm_models)
lm_stats$taxatree_stats

# adjust for multiple comparisons
lm_stats <- taxatree_stats_p_adjust(
  data = lm_stats, method = "BH", grouping = "rank"
)

# which are significant?
lm_stats$taxatree_stats %>% filter(p.adj.BH.rank < 0.05)

# what are the effect sizes?
hist(lm_stats$taxatree_stats$estimate)

# which taxa are changing?
unique(lm_stats$taxatree_stats$taxon[lm_stats$taxatree_stats$p.adj.BH.rank < 0.05])

# save for later
mod2.5 <- lm_stats

## ---- test 3: dose response ----

# get data
phylo <- ps %>% 
  ps_filter(!Treatment %in% c("Control", "Males")) %>% 
  tax_fix()

# model on Genus
lm_models <- phylo %>% 
  tax_fix() %>% 
  tax_prepend_ranks() %>% 
  tax_transform("compositional", rank = "Genus") %>% 
  tax_transform(
    trans = "log2", 
    # chain allows for multiple transformations (TSS then log2)
    chain = TRUE, 
    # for log, replace 0's with half of min value
    zero_replace = "halfmin"
  ) %>% 
  taxatree_models(
    type = lm, 
    ranks = "Genus",
    variables = c("Treatment")
  )

# get stats
lm_stats <- taxatree_models2stats(lm_models)
lm_stats$taxatree_stats

# adjust for multiple comparisons
lm_stats <- taxatree_stats_p_adjust(
  data = lm_stats, method = "BH", grouping = "rank"
)

# which are significant?
lm_stats$taxatree_stats %>% filter(p.adj.BH.rank < 0.05)

# what are the effect sizes?
hist(lm_stats$taxatree_stats$estimate)

# which taxa are changing?
unique(lm_stats$taxatree_stats$taxon[lm_stats$taxatree_stats$p.adj.BH.rank < 0.05])

# save for later
mod3 <- lm_stats
