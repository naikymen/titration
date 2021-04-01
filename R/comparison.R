library(dplyr)
library(ggplot2)

# Overlay solutions from diferent methods ####

# Analytical methods
# source("R/phosphoric_titration_analytic.R")
source("R/phosphoric_titration_analytic.vol_correct.R")
pH.seq <- seq(from=1,to=12,length.out = 1000)
Na.seq <- Na(H = 10^-pH.seq, P_ca=0.04389)
Na.seq.adj <- Na.adj(H=10^-pH.seq, P.ca=0.04389)

ann.result.df <- data.frame(pH=pH.seq, Na=Na.seq$Na)
ann.result.df.adj <- data.frame(pH=pH.seq, Na=Na.seq.adj$Na)

# Numerical methods
num.result.df <- readRDS("output/result.df.0.04389.RDS") %>%
  mutate(pH = -log10(H)) %>% select(pH, Na)

# Experimental data
# Titration of 25 mL of H3PO4 (0.04389 M), with NaOH (0.09948 M).
# Veq1=11.03 mL, Veq2=22.06 mL.
source("R/phosphoric_experimental.R")
exp.result.df <- exp.result()


titration.curves <- bind_rows(
  analytic=ann.result.df,
  analytic.adj=ann.result.df.adj,
  numeric=num.result.df,
  experimental=exp.result.df,
  .id = "method")

ggplot(titration.curves, mapping=aes(x=Na, y=pH)) + 
  geom_path(aes(color=method), size=2, alpha=.5) +
  theme_minimal() + 
  ggtitle("Titration curves for 0.04389 M Phosphoric acid",
          "Note: negative [Na+] values were removed") +
  xlab("[Na+]") +
  xlim(c(0,NA))

# ggsave("output/comparison2.png")
