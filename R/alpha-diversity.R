## alpha diversity analysis

# EVS 12/2021


require(tidyverse)
require(phyloseq)
require(microViz)
require(rstatix)
require(ggpubr)

load("./data/ps-decontam-filtered-counts.RData")
ps <- pscount %>% 
  ps_mutate(Time = recode(Sample.Date, '8/8/19' = 1, '10/16/19' = 2, '12/20/19' = 3))

## ---- get alpha diversity metrics -----

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

## ----- Observed ----

#exploratory plots
obs <- alpha %>% select(-c(Shannon, Simpson))
hist(obs$Observed)
ggboxplot(obs, x = "Treatment", y = "Observed", facet.by = "Time")

## ---- Shannon ----

#exploratory plots
sh <- alpha %>% select(-c(Observed, Simpson))
hist(sh$Shannon)
ggboxplot(sh, x = "Treatment", y = "Shannon", facet.by = "Time")

# build test model for linear assumption
testm <- lm(Shannon ~ Treatment, data = filter(sh, Treatment %in% c("Control", "T4")))
summary(testm)
plot(testm)
hist(resid(testm))
qqnorm(resid(testm))
qqline(resid(testm))
shapiro.test(resid(testm)) # everything looks good for linear tests

## test 1: is there a difference between control and T4 (control and metformin)?
mod1 <- t.test(Shannon ~ Treatment, data = filter(sh, Treatment %in% c("Control", "T4") & Time == 3))
ggboxplot(filter(sh, Treatment %in% c("Control", "T4") & Time == 3),
          x = "Treatment", y = "Shannon") +
  stat_compare_means(method = 't.test')
# no difference

## test 2: is there a difference in metformin treatment over time


mod2 <- aov(Shannon ~ Time, data = filter(sh, !Treatment %in% c("Males", "Control")))
summary(mod2)
ggboxplot(filter(sh, !Treatment %in% c("Males", "Control")),
          x = "Time", y = "Shannon") +
  stat_compare_means(method = "anova")
# post hoc test: Tukey
TukeyHSD(mod2) # differences between 2 and 3 (marginal between 1 and 3)

## sanity check: is there a difference in control over time
mod3 <- aov(Shannon ~ Time, data = filter(sh, Treatment == "Control"))
summary(mod3) # no difference
ggboxplot(filter(sh, Treatment == "Control"), x = "Time", y = "Shannon")

## test 3: is there a dose response
mod4 <- aov(Shannon ~ Treatment, data = filter(sh, Treatment %in% c("T2", "T3", "T4")))
summary(mod4) # none

## ---- Simpson ----

# exploratory plots
si <- alpha %>% select(-c(Observed, Shannon))
hist(si$Simpson)
ggboxplot(si, x = "Treatment", y = "Simpson", facet.by = "Time")

# build test model for linear assumption
testm <- lm(Simpson ~ Treatment, data = filter(si, Treatment %in% c("Control", "T4")))
summary(testm)
plot(testm)
hist(resid(testm))
qqnorm(resid(testm))
qqline(resid(testm))
shapiro.test(resid(testm)) # a little iffy

## test 1: is there a difference between control and T4 (control and metformin)?
mod1 <- t.test(Simpson ~ Treatment, data = filter(si, Treatment %in% c("Control", "T4") & Time == 3))
ggboxplot(filter(si, Treatment %in% c("Control", "T4") & Time == 3),
          x = "Treatment", y = "Simpson") +
  stat_compare_means(method = 't.test')
# no difference

## test 2: is there a difference in metformin treatment over time
mod2 <- aov(Simpson ~ Time, data = filter(si, !Treatment %in% c("Males", "Control")))
summary(mod2)
ggboxplot(filter(si, !Treatment %in% c("Males", "Control")),
          x = "Time", y = "Simpson") +
  stat_compare_means(method = "anova")
# no difference

## sanity check: is there a difference in control over time
mod3 <- aov(Simpson ~ Time, data = filter(si, Treatment == "Control"))
summary(mod3) # no difference
ggboxplot(filter(si, Treatment == "Control"), x = "Time", y = "Simpson")

## test 3: is there a dose response
mod4 <- aov(Simpson ~ Treatment, data = filter(si, Treatment %in% c("T2", "T3", "T4")))
summary(mod4) # none
