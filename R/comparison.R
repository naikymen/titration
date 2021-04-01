library(dplyr)
library(ggplot2)

# Overlay to numerical solution ####
source("R/phosphoric_titration_analytic.R")
pH.seq <- seq(from=1.5,to=13,length.out = 1000)
Na.seq <- Na(H = 10^-pH.seq, P_ca=0.04389)

ann.result.df <- data.frame(pH=pH.seq, Na=Na.seq)

num.result.df <- readRDS("output/result.df.0.04389.RDS") %>%
  mutate(pH = -log10(H)) %>% select(pH, Na)

# Titration of 25 mL of H3PO4 (0.04389 M)
# with NaOH (0.09948 M)
# (Veq1=11.03 mL; Veq2=22.06 mL).
source("R/phosphoric_experimental.R")
# dubious_scaling_factor <- 1.8
# dubious_offset_factor <- 0.4
# exp.result.df <- exp.result.df %>% 
#   mutate(Na = Na*dubious_scaling_factor,
#          pH = pH + dubious_offset_factor)

titration.curves <- bind_rows(
  analytic=ann.result.df,
  numeric=num.result.df,
  experimental=exp.result.df,
  .id = "method")

ggplot(titration.curves, mapping=aes(x=Na, y=pH)) + 
  geom_path(aes(color=method)) +
  theme_minimal() + 
  ggtitle("Titration curves for 0.04389 M Phosphoric acid") +
  xlab("[Na+]")

ggsave("output/comparison.pdf")
