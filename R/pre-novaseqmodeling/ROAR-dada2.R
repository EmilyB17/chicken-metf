
## dada2 results from ROAR

# Emily Van Syoc, 6/2021

# load packages
require(dada2)
require(phyloseq)
require(tidyverse)

## load data from ROAR
load("../dada2/asv-table.RData")
load("../dada2/tax-table.RData")


# make the sample names match
otu <- otu_table(seqtab.nochim, taxa_are_rows = FALSE)
rownames(otu) <- str_extract(rownames(otu), "M(\\d){1,3}")

# read in sample data 
dat <- read.table("./data/sampleIDs-with-sequencing-controls.txt", sep = "\t", header = TRUE) %>% 
  column_to_rownames(var = "NovaGeneID")

# make phyloseq object
ps <- phyloseq(otu,
               tax_table(tax),
               sample_data(dat))

# change DNA strings to ASV number
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

# save the object to RData
save(ps, file = "./data/ps-object.RData")

#### ---- Edit phyloseq for downstream use

# remove all taxa that are completely unassigned through the Kingdom level (n = 241)

psEdit <- subset_taxa(ps, !Kingdom == "NA")

# remove all non-bacterial reads (n = 97)
psEdit <- subset_taxa(ps, Kingdom == "Bacteria")

# save the object
save(psEdit, file = "./data/ps-object-noNA.RData")
