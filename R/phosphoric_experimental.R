library(dplyr)
library(ggplot2)

# Titration of 25 mL of H3PO4 (0.04389 M)
# with NaOH (0.09948 M)
# (Veq1=11.03 mL; Veq2=22.06 mL).

# dubious_scaling_factor <- 2.47
# exp.result.df <- read.csv("data/titrPho.csv") %>%
#   mutate(Na=dubious_scaling_factor*(Volume.mL-min(Volume.mL))/max(Volume.mL)) %>% 
#   select(Na, pH)

NaOH.Conc <- 0.09948
Acid.Vol <- 25

exp.result.df <- read.csv("data/titrPho.csv") %>%
  dplyr::rename(Na.Vol=Volume.mL)  %>% 
  mutate(Na = Na.Vol*NaOH.Conc/(Na.Vol+Acid.Vol)) %>% 
  select(Na, pH)
