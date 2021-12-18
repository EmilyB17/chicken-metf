## Preliminary data report 

# 1/2022

# EVS

## ---- getData ----

# load packages
require(tidyverse)
require(ggpubr)
require(vegan)
require(phyloseq)
require(microViz)
require(rstatix)
require(DESeq2)

# load RData
load("./data/ps-decontam-filtered-counts.RData")
# add Time
ps <- pscount %>% 
  ps_mutate(Time = fct_recode(Sample.Date, "1" = "8/8/19", "2" = "10/16/19", "3" = "12/20/19")) %>% 
  tax_fix()

### ---- `relAbun` ----

# number of taxa
#table(tax_table(ps)[, "Phylum"], exclude = NULL) 
#length(get_taxa_unique(ps, "Genus"))
#length(get_taxa_unique(ps, "Phylum"))


## plot
ggplot(data = psmelt(ps), mapping = aes_string(x = "Treatment", y = "Abundance")) +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill") +
  labs(x = NULL, y = "Relative Abundance",
       title = "Relative abundance of all phylum by time & treatment") +
  theme_bw() +
  facet_grid(~Time)

# get top 10 most abundance taxa
top <- prune_taxa(names(sort(taxa_sums(ps), TRUE))[1:10], ps)
ggplot(data = psmelt(top), mapping = aes_string(x = "Treatment", y = "Abundance")) +
  geom_bar(aes(fill = Genus), stat = "identity", position = "fill") +
  labs(x = NULL, y = "Relative Abundance",
       title = "Top 10 most relatively abundant genera") +
  theme_bw() +
  facet_grid(~Time)

## ---- AlphaDiversity ----

# get metrics
sampdf <- sample_data(ps) %>% 
  data.frame() %>% 
  rownames_to_column(var = "ID")

# get diversity
alpha <- estimate_richness(ps, measures = c("Shannon", "Simpson", "Observed")) %>% 
  # make ID
  rownames_to_column(var = "ID") %>% 
  merge(sampdf, by = "ID") %>% 
  # make Time an ordered factor
  mutate(Time = factor(Time, ordered = TRUE, levels = c(1, 2, 3)))

# plot Shannon's sig diffs
ggboxplot(filter(alpha, !Treatment %in% c("Males", "Control")),
          x = "Time", y = "Shannon",
          add = "dotplot", add.params = list(color = "Treatment"),
          title = "Shannon's Diversity of all treatments compared to time") +
  stat_compare_means(method = "anova", label.y = 4.8)

## ---- BetaDiversity ----

### test 1: difference between control and T4

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

##  test 2: difference in metformin treatment over time 

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

## test 3: metformin treatment dose response 

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

## ---- RelativeAbundance ----

ps <- ps %>% tax_fix()

## test1: difference between control and T4 

# get data
t1 <- ps %>% 
  ps_filter(Treatment %in% c("Control", "T4")) %>% 
  tax_fix()

# phyloseq-DEseq
dds <- phyloseq_to_deseq2(t1, ~ Treatment)

# perform DESeq
ds <- DESeq(dds, test = "Wald", fitType = "parametric")

# get results
res = results(ds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(t1)[rownames(sigtab), ], "matrix"))


## plot
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=4) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  ggtitle("Metformin T4 versus Control")

## test2: difference between metformin over time 

# only pairwise comparisons, so use times 1 and 3
# get data
t1 <- ps %>% 
  ps_filter(!Treatment %in% c("Control", "Males")) %>% 
  ps_filter(!Time == 2) %>% 
  tax_fix()

# phyloseq-DEseq
dds <- phyloseq_to_deseq2(t1, ~ Time)

# perform DESeq
ds <- DESeq(dds, test = "Wald", fitType = "parametric")

# get results
res = results(ds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(t1)[rownames(sigtab), ], "matrix"))


## plot
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=4) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  ggtitle("Metformin treatments T3 versus T1")

### sanity check: control over time 

# get data
t1 <- ps %>% 
  ps_filter(Treatment %in% "Control") %>% 
  ps_filter(!Time == 2) %>% 
  tax_fix()

# phyloseq-DEseq
dds <- phyloseq_to_deseq2(t1, ~ Time)

# perform DESeq
ds <- DESeq(dds, test = "Wald", fitType = "parametric")

# get results
res = results(ds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(t1)[rownames(sigtab), ], "matrix"))


## plot
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=4) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  ggtitle("Control over time")

## test 3: dose response 

# get data
t1 <- ps %>% 
  ps_filter(!Treatment %in% c("Control", "Males")) %>% 
  ps_filter(!Treatment %in% c("T2")) %>% 
  tax_fix()

# phyloseq-DEseq
dds <- phyloseq_to_deseq2(t1, ~ Treatment)

# perform DESeq
ds <- DESeq(dds, test = "Wald", fitType = "parametric")

# get results
res = results(ds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(t1)[rownames(sigtab), ], "matrix"))


## plot
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=4) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  ggtitle("Dose Response T4 versus T2")
