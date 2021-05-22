library(dplyr)
library(ggplot2)

# Overlay solutions from diferent methods ####
P.vol=0.025    # 25 mL
P.ca=0.04389   # 0.04389 M
P.mass <- P.ca*P.vol # moles
Na.ca=0.09948  # 0.09948
pH.seq <- seq(from=1, to=12, length.out = 1000)

# Analytical methods
# source("R/phosphoric_titration_analytic.R")
source("R/phosphoric_titration_analytic.vol_correct.R")
# Not corrected
Na.seq      <- (Na.aprox(H=10^-pH.seq, P.ca=P.ca)$Na * P.vol) / Na.ca
# Numeric correction
Na.seq.adj  <- Na.adj(pH.seq=pH.seq, P.ca = P.ca, P.vol = P.vol, Na.ca = Na.ca, n.its = 200)$Na.vol
# Analytic correction
Na.seq.adj2 <- Na.adj2(pH.seq=pH.seq, P.mass = P.mass, P.vol = P.vol, Na.ca = Na.ca)$Na.vol

ann.result.df     <- data.frame(pH=pH.seq, Na.vol=Na.seq)
ann.result.df.adj <- data.frame(pH=pH.seq, Na.vol=Na.seq.adj)
ann.result.df.adj2 <- data.frame(pH=pH.seq, Na.vol=Na.seq.adj2)

# Numerical methods
num.result.df <- readRDS("output/result.df.0.04389.RDS") %>%
  mutate(pH = -log10(H)) %>% 
  mutate(Na.vol = (Na * P.vol) / Na.ca) %>% 
  select(pH, Na.vol)

# Experimental data
# Titration of 25 mL of H3PO4 (0.04389 M), with NaOH (0.09948 M).
# Veq1=11.03 mL, Veq2=22.06 mL.
source("R/phosphoric_experimental.R")
exp.result.df <- exp.result()


titration.curves <- bind_rows(
  numeric.uncorrected  = num.result.df,
  analytic.uncorrected = ann.result.df,
  analytic.corrected      = ann.result.df.adj2,
  analytic.corrected.iter = ann.result.df.adj,
  experimental = exp.result.df,
  .id = "method")

titration.curves %>% 
  filter(Na.vol >= 0 & Na.vol <= 0.03) %>% 
  ggplot(mapping=aes(x=Na.vol, y=pH)) + 
  geom_path(size=1,
            alpha=.7, data = exp.result.df) +
  geom_path(aes(color=method,size=as.factor(method)), 
            size=1,
            alpha=.7) +
  theme_minimal() + theme(legend.position = "top") +
  ggtitle("Titration curves for 0.04389 M Phosphoric acid",
          "Note: negative Na.vol values were removed") +
  facet_wrap(~method)

titration.curves %>% 
  filter(Na.vol >= 0 & Na.vol <= 0.03) %>% 
  ggplot(mapping=aes(x=Na.vol, y=pH)) + 
  geom_hline(yintercept = -log10(c(Ka.1, Ka.2, Ka.3)), alpha = .5, color = "black") +
  geom_path(aes(color=method,size=as.factor(method)), 
            size=2,
            alpha=.5) +
  theme_minimal() + theme(legend.position = "top") +
  ggtitle("Titration curves for 0.04389 M Phosphoric acid",
          "pKas are drawn as gray horizontal lines")

# ggsave("output/comparison4.w_experimental.png", width = 8)
