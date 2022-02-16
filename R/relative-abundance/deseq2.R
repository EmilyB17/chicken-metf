## relative abundance

## EVS 12/2021


require(tidyverse)
require(phyloseq)
require(microViz)
require(rstatix)
require(DESeq2)

load("./data/ps-decontam-filtered-counts.RData")
ps <- pscount %>% 
  ps_mutate(Time = as.factor(recode(Sample.Date, '8/8/19' = '40 weeks', '10/16/19' = '50 weeks', '12/20/19' = '60 weeks'))) %>% 
  tax_fix()

## plot
ggplot(data = psmelt(ps), mapping = aes_string(x = "Treatment", y = "Abundance")) +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill") +
  labs(x = NULL, y = "Relative Abundance",
       title = "Phylum") +
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

# what phyla do the most abundant genera belong to?
t <- psmelt(top)

## ---- test1: difference between control and T4 ----

# get data
t1 <- ps %>% 
  ps_filter(Treatment %in% c("Control", "T4")) %>% 
  # dummy code
  ps_mutate(tdummy = recode(Treatment, Control = 0, T4 = 1)) %>%
  tax_glom("Genus") %>% 
  tax_fix()

# phyloseq-DEseq
dds <- phyloseq_to_deseq2(t1, ~ tdummy)

# perform DESeq
ds <- DESeq(dds, test = "Wald", fitType = "parametric")

# get results
res = results(ds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(t1)[rownames(sigtab), ], "matrix"))

# how many were more abundant in T4 than control
length(which(sigtab$log2FoldChange > 0))
unique(sigtab$Phylum[sigtab$log2FoldChange > 0])

# how many were more abundant in control than T4
length(which(sigtab$log2FoldChange < 0))
unique(sigtab$Phylum[sigtab$log2FoldChange < 0])


## plot
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  ggtitle("75 mg/kg metformin compared to control")

## ---- test2: difference between metformin over time ----

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

# which were higher in 60 weeks than 40 weeks
length(which(sigtab$log2FoldChange > 0))

# get these genera and log2fold change to compare to the control treatment
mets <- sigtab %>% select(Genus, log2FoldChange) %>% mutate(type = "metformin")

## plot
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

### ---- sanity check: control over time ----

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

# which were higher in 60 weeks
length(which(sigtab$log2FoldChange > 0))

# get the genera and log2fold change to compare to metformin birds
con <- sigtab %>% select(Genus, log2FoldChange) %>% mutate(type = "control")


## plot
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

### compare metformin and control time responses
both <- rbind(mets, con) %>% 
  mutate(Genus = paste(Genus, rownames(both), sep = "_")) %>% 
  pivot_wider(names_from = "type", values_from = "log2FoldChange", values_fill = NA)


## ---- test 3: dose response ----

# get data
t1 <- ps %>% 
  ps_filter(!Treatment %in% c("Control", "Males")) %>% 
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

# which are higher in T4 than T2
length(which(sigtab$log2FoldChange > 0))
length(which(sigtab$log2FoldChange < 0))

## plot
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  ggtitle("75 mg/kg metformin compared to 25 mg/kg metformin")
