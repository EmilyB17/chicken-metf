## wrangle data for Access

# require package
require(tidyverse)
require(readxl)

# get the key
key <- read_excel("../selected-Cage bird numbers for microbiome collab.xlsx") %>% 
  select(Treatment, Bird)

# get the key for 12/20 samples
key.dec <- read_excel("../SampleIDs.xlsx", sheet = "temp-12-20-key" ) %>% 
  # remove dead birds
  filter(!Sample == "dead") %>% 
  mutate(Bird = as.character(Bird))

# get data
dat <- read_excel("../SampleIDs.xlsx") %>% 
  # get the columns we want
  select(NovaGeneID, SamplingTubeID, ExtraTubeID, LabID)

# remove sequencing controls
datnocon <- dat %>% 
  # remove sequence controls
  filter(!str_detect(SamplingTubeID, "PosControl")) %>% 
  filter(!str_detect(SamplingTubeID, "NegControl"))

## ---- 12/20/19 labels ----

## get 12-20-19 sampling
dec <- datnocon %>% 
  filter(str_detect(SamplingTubeID, "12/20/19")) %>% 
  # get tube ID and the handwritten tube ID
  mutate(labelID = str_extract(SamplingTubeID, "^(\\d){1,3}")) %>% 
  select(-LabID)

# match to the key
matchedid <- dec %>% 
  # the tubeID is the Sample unless it is over 50, then it is the bird
  inner_join(key.dec, by = c("labelID" = "Sample")) %>% 
  rename(ExtraTubeID = labelID) %>% 
  select(NovaGeneID, SamplingTubeID, LabID, Bird, Treatment)

# match by bird
matched2 <- dec %>% 
  inner_join(key.dec, by = c("labelID" = "Bird")) %>% 
  rename(Bird = labelID) %>% 
  select(NovaGeneID, SamplingTubeID, LabID, Bird, Treatment)

# did anything not match?
nomatch <- dec %>% 
  anti_join(key.dec, by = c("labelID" = "Sample")) %>% 
  anti_join(key.dec, by = c("labelID" = "Bird")) %>% 
  # i checked this by hand and it should be bird 85 (bird 84 died)
  mutate(bird.correct = "85") %>% 
  inner_join(key.dec, by = c("bird.correct" = "Bird")) %>% 
  rename(Bird = bird.correct) %>% 
  select(NovaGeneID, SamplingTubeID, LabID, Bird, Treatment)

## join together
all.dec <- rbind(matchedid, matched2, nomatch) %>% 
  mutate(Sample.Date = "12/20/19") %>% 
  select(NovaGeneID, Bird, Treatment, Sample.Date)

## ---- Other times ----

# remove the wonky 12/19
datnodec <- datnocon %>% 
  filter(!str_detect(SamplingTubeID, "12/20/19-50$")) %>% 
  filter(!str_detect(SamplingTubeID, "12/20")) %>% 
  # get bird name
  mutate(Bird = as.numeric(str_extract(SamplingTubeID, "^(\\d){1,3}")),
         # get treatment
         Treatment = case_when(
           str_detect(SamplingTubeID, "T(\\d)") ~ str_extract(SamplingTubeID, "T(\\d)"),
           str_detect(SamplingTubeID, "Males") ~ "Males",
           str_detect(SamplingTubeID, "Control") ~ "Control"
         ),
         # get sampling time
         Sample.Date = str_extract(SamplingTubeID, "(\\d){1,2}/(\\d){1,2}/19")) 


## ---- clean up and export ----

# clean up
clean <- datnodec %>% 
  # get sampling time from SamplingTubeID
  mutate(Sample.Date = str_extract(SamplingTubeID, "(\\d){1,2}/(\\d){1,2}/19$")) %>% 
  select(NovaGeneID, Bird, Treatment, Sample.Date) %>% 
  # add December samples
  rbind(all.dec) %>% 
  mutate(Treatment = str_replace_all(Treatment, "Treatment ", "T"),
         Bird = as.numeric(Bird))
# from Evelyn: Bird 74 is a Control
clean$Treatment[clean$Bird == "74"] <- "Control"
# There is a TYPO: Bird 241 is TREATMENT 2
clean$Treatment[clean$Bird == "241"] <- "T2"

# write to table
#write.table(clean, file = "./data/sampleIDs.txt", sep = "\t", row.names = FALSE)

# add sequencing controls back
cont <- read_excel("../SampleIDs.xlsx") %>% 
  # get the columns we want
  select(NovaGeneID, SamplingTubeID, ExtraTubeID, LabID) %>% 
  # get controls
  filter(str_detect(SamplingTubeID, "Control")) %>% 
  # remove control birds
  filter(!str_detect(SamplingTubeID, "-Control-(\\d)")) %>% 
  # add correct columns
  mutate(Bird = "Sequence-Control",
         Treatment = ExtraTubeID,
         Sample.Date = "NA") %>% 
  select(NovaGeneID, Bird, Treatment, Sample.Date) %>% 
  # add the samples
  rbind(clean)

# write to table
#write.table(cont, file = "./data/sampleIDs-with-sequencing-controls.txt", sep = "\t", row.names = FALSE)

## ---- check against Evelyn's key ----

ekey <- read_excel("../SampleIDs.xlsx", sheet = "ids-from-Evelyn") %>% 
  rename(Cage = `Cage#`) 

# check
check <- ekey %>% 
  anti_join(clean, by = c("Cage" = "Bird")) # 4 don't match, 2 controls and 2 males

matches <- ekey %>% 
  inner_join(clean, by = c("Cage" = "Bird"))

# anti check
check1 <- clean %>% 
  anti_join(ekey, by = c("Bird" = "Cage"))


