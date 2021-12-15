## Decontamination and filtering

# Emily Van Syoc 11/21

# load packages
require(tidyverse)
require(phyloseq)
require(dada2)
# install latest version of microViz
#devtools::install_github("david-barnett/microViz@0.7.1")
require(microViz)
#require(BiocManager)
#BiocManager::install("decontam")
require(decontam)

# get data
load("./data/raw-ps.RData")

## ---- basic info ----

ps 

ntaxa(ps) # 24,146 taxa

get_taxa_unique(ps, "Phylum") # 30 different phyla including NA
get_taxa_unique(ps, "Kingdom") # 4 kingdoms - need to filter this first

# ---- filter for Bacteria & remove NA phyla ----

# get only bacteria
psf <- subset_taxa(ps, Kingdom == "Bacteria")

# remove NA phylum
psf1 <- subset_taxa(psf, !is.na(Phylum) & !Phylum %in% c("", "NA"))

# validate
psf2 <- tax_fix(psf1)
psf3 <- phyloseq_validate(psf2, remove_undetected = TRUE)

## ---- filter for relative abundance ----

# transform to relative abundance
pst <- transform_sample_counts(psf3, function(x) x / sum(x) )

# remove taxa with total relative abundance less than 10e-5
psr <- filter_taxa(pst, function(x) mean(x) > 1e-5, TRUE)

## ---- decontam ----

## workflow from: https://benjjneb.github.io/decontam/vignettes/decontam_intro.html

# add sampling variables
psd <- psr %>% 
  ps_mutate(
    Time = recode(Sample.Date, '8/8/19' = 1, '10/16/19' = 2, '12/20/19' = 3),
    SampleBinary = if_else(Bird == "Sequence-Control", true = "Control", false = "Sample")
  ) 

## remove positive controls for this 
## update - keep and remove contaminants like a sample
#psnopos <- psd %>% ps_filter(!Treatment %in% c("PCP2", "PCP3", "PCP1", "PCISO1"))

## inspect library sizes 
sampdf <- data.frame(sample_data(psd))
sampdf$LibrarySize <- sample_sums(psd)
sampdf <- sampdf[order(sampdf$LibrarySize), ]
sampdf$Index <- seq(nrow(sampdf))
ggplot(data = sampdf, aes(x = Index, y = LibrarySize, color = SampleBinary)) + geom_point()
## interesting - two controls are scattered throughout

## use negative controls to identify contaminants
# add distinguisher for positive controls and treat as a sample
psp <- psd %>% ps_mutate(SampleBinary = case_when(
  SampleBinary == "Control" & str_detect(Treatment, "PC") ~ "Sample",
  SampleBinary == "Control" & !str_detect(Treatment, "PC") ~ "Control",
  SampleBinary == "Sample" ~ "Sample"
))
pds <- psp %>% ps_mutate(is.neg = if_else(SampleBinary == "Control", true = TRUE, false = FALSE))
contamdf.prev <- isContaminant(pds, method="prevalence", neg="is.neg", threshold = 0.5) # more aggressive threshold is 0.5
table(contamdf.prev$contaminant) # this identifies 673 contaminants and 27059 'true' taxa

## identify the contaminant taxa and remove them from the phyloseq object
psdecon <- prune_taxa(!contamdf.prev$contaminant, pds)

## remove the negative controls and save the decontaminated phyloseq object
pssave <- ps_filter(psdecon, SampleBinary == "Sample") %>% 
  # remove positive controls
  ps_filter(!str_detect(Treatment, "PC")) %>% 
  ps_select(-c(is.neg, SampleBinary))

# write
save(pssave, file = "./data/ps-decontam-filtered.RData")

## ----- positive control -----

## get positive controls only
pspos <- psdecon %>% ps_filter(Treatment %in% c("PCP2", "PCP3", "PCP1", "PCISO1"))

# there are 981 unique taxa
ntaxa(pspos)

pspos <- tax_fix(pspos)
# visualize in relative abundance barplot
pspos %>%
  comp_barplot(
    tax_level = "Family", n_taxa = 15,  
    merge_other = FALSE 
    # set merge_other = TRUE (the default) to remove outlines from inside "other" category
  ) +
  facet_wrap("Treatment", scales = "free") +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

## ---- get phyloseq of counts ----

## we need the count data with the taxa present in the decontaminated and filtered PS
pscount <- subset_taxa(ps, taxa_names(ps) %in% taxa_names(pssave)) %>% 
  # remove controls
  ps_filter(!Bird == "Sequence-Control")

# save
save(pscount, file = "./data/ps-decontam-filtered-counts.RData")

## ---- Relative abundance barplot ----


## all samples
psdecon %>% 
  comp_barplot(
    tax_level = "Phylum", n_taxa = 10,
    merge_other = TRUE,
  )+
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) 


