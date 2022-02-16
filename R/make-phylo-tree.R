## make phylogenetic tree
# EVS 2/2022

require(phyloseq)
# load phyloseq
load("data/ps-decontam-filtered-counts.RData")
load("data/raw-ps.RData")

psf <- subset_taxa(ps, taxa_names(ps) %in% taxa_names(pscount))

# get DNA strings as taxa names
otu <- otu_table(psf)
taxa_names(otu) <- psf@refseq[,'seq']


## from: https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html#introduction
library("knitr")
BiocManager::install("BiocStyle")
library("BiocStyle")
.cran_packages <- c("ggplot2", "gridExtra")
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  #source("http://bioconductor.org/biocLite.R")
  #biocLite(.bioc_packages[!.inst], ask = F)
  BiocManager::install(.bioc_packages[!.inst], ask = F)
}
# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

# make phylo tree
seqs <- getSequences(otu)


## subset to test
seqs <- seqs[1:100]
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)

phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)