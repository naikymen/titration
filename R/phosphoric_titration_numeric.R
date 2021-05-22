library(dplyr)
library(tidyr)
library(ggplot2)

# Numerical solution for the triprotic buffer system
# Using 1M phosphoric acid
# pseudo-Titrating with NaOH (not volume corrected, i.e. assumes solid NaOH)

Ka.1 <- 7.1 * 10^-3
Ka.2 <- 6.3 * 10^-8
Ka.3 <- 4.5 * 10^-13
Kw <- 10^-14

P_ca <- 0.04389
ci.start <- c(H=0.1, H3A=0.9, H2A=0.1, HA=0.1, A=0.1)

Na.seq <- seq(from=0,to=4*P_ca,by=P_ca/1000)

# SETUP ####
balance.jabobian <- function(x, convert.fun, ...){
  n <- length(x)
  Df <- matrix(numeric(n*n),n,n)
  x <- convert.fun(x)
  
  # d(BM) H3A + H2A + HA + A - P_ca
  Df[1,1] <- 0    # /d(H)
  Df[1,2] <- 1    # /d(H3A)
  Df[1,3] <- 1    # /d(H2A)
  Df[1,4] <- 1    # /d(HA)
  Df[1,5] <- 1    # /d(A)
  
  # d(BC) H + Na - Kw/H - H2A - 2*HA - 3*A
  Df[2,1] <- 1+Kw/(x[1]^2)   # /d(H)
  Df[2,2] <- 0               # /d(H3A)
  Df[2,3] <- -1              # /d(H2A)
  Df[2,4] <- -2              # /d(HA)
  Df[2,5] <- -3              # /d(A)
  
  # d(EQ1) H * H2A / Ka.1 - H3A
  Df[3,1] <-  x[3]/Ka.1  # /d(H)
  Df[3,2] <-  -1         # /d(H3A)
  Df[3,3] <-  x[1]/Ka.1  # /d(H2A)
  Df[3,4] <-  0          # /d(HA)
  Df[3,5] <-  0          # /d(A)
  
  # d(EQ2) H * HA  / Ka.2 - H2A
  Df[4,1] <-  x[4]/Ka.2  # /d(H)
  Df[4,2] <-  0          # /d(H3A)
  Df[4,3] <-  -1         # /d(H2A)
  Df[4,4] <-  x[1]/Ka.2  # /d(HA)
  Df[4,5] <-  0          # /d(A)
  
  # d(EQ3) H * A   / Ka.3 - HA
  Df[5,1] <-  x[5]/Ka.3  # /d(H)
  Df[5,2] <-  0          # /d(H3A)
  Df[5,3] <-  0          # /d(H2A)
  Df[5,4] <-  -1         # /d(HA)
  Df[5,5] <-  x[1]/Ka.3  # /d(A)
  
  Df
}

balance <- function(vars, Na_ca, P_ca, convert.fun){
  # Apply positive only constraint
  vars <- convert.fun(vars)
  
  H <- vars[1]
  H3A <- vars[2]
  H2A <- vars[3]
  HA <- vars[4]
  A <- vars[5]
  
  Na <- convert.fun(Na_ca)
  
  eq.system <- c(
    BM = H3A + H2A + HA + A - P_ca,
    # safe.sum(H3A, + H2A, + HA, + A, - P_ca),
    BC = H + Na - Kw/H - H2A - 2*HA - 3*A,
    # safe.sum(H, +Na, -Kw/H, -H2A, -2*HA, -3*A),
    
    A1 = H * H2A/Ka.1 - H3A,
    # A1 = H * H2A - H3A*Ka.1,
    # A1 = (H * H2A) / (Ka.1*H3A) - 1,
    # safe.eq(H3A,H2A,Ka.1,H),
    
    A2 = H * HA/Ka.2 - H2A,
    # A2 = H * HA - H2A*Ka.2,
    # A2 = (H * HA) / (Ka.2*H2A) - 1,
    # safe.eq(H2A,HA,Ka.2,H),
    
    A3 = H * A/Ka.3 - HA
    # A3 = H * A - HA*Ka.3
    # A3 = (H * A) / (Ka.3*HA) - 1
    # safe.eq(HA,A,Ka.3,H)
  )
  return(eq.system)
}

# RUN SOLVER ####

varnames <- c("Na", "H",  "H3A", "H2A", "HA",  "A", "termcd")

result.m <- matrix(ncol = length(varnames), nrow = length(Na.seq))
colnames(result.m) <- varnames

result.m[,1] <- Na.seq

convert.fun <- function(x) abs(x)
convert.fun.i <- function(x) abs(x)

for(i in 1:length(Na.seq)){
  
  Na_ca <- result.m[i,1]
  
  idx.good <- which(result.m[,7] == 1)
  
  if(i == 1){                          # Si es la primera iteración,
    ci <- ci.start                     # usar los valores "start" como C.I.
  } else if(length(idx.good)>0){       # Si hay algun valor bueno en la matriz
    ci <- result.m[max(idx.good),2:6]  # usar el más reciente como C.I.
  } else {                             # Si no,
    ci <- result.m[i-1, 2:6]           # usar los valores de la solución anterior
  }
  
  iter.trace <- 0
  # iter.trace <- 1
  # cat("Iteration ",i,"\n\n")
  result <- nleqslv::nleqslv(x = convert.fun.i(ci), 
                             fn = balance, 
                             jac = balance.jabobian,
                             Na = Na_ca, P = P_ca,
                             convert.fun = convert.fun,
                             method="Newton",
                             # method="Broyden",
                             global="dbldog",
                             # xscalm = "auto",
                             control = list(allowSingular=TRUE,
                                            xtol=1e-10,
                                            maxit=1000,
                                            # scalex=rep(1e-12, 5),
                                            scalex=1/ci,
                                            # delta=1E,
                                            # stepmax=1E5,
                                            trace=iter.trace))
  
  # cat("\n\n ",result$message,"\n\n") 
  
  result$x <- convert.fun(result$x)
  
  result.m[i,2:6] <- result$x
  result.m[i,7] <- result$termcd
  
  stopifnot(all(result$x >= 0))
  
} # END LOOP

# tail(result.m[!is.na(result.m[,2]),], n = 10)
result.df <- as.data.frame(result.m)
# saveRDS(result.df, "output/result.df.RDS")
# saveRDS(result.df, "output/result.df.0.04389.RDS")

# PLOT SOLVER RETURN CODES ####
plot(termcd~Na, result.df,pch=20)
stopifnot(all(result.df$termcd==1))


# PLOT TITRATION CURVE ####
result.df <- as.data.frame(result.m)
result.df$termcd <- c("ok", "xtol", "stall", "maxit", "jacobad")[result.df$termcd]

result.df.long <- result.df %>%
  mutate(pH = -log10(H)) %>% select(-H) %>%
  pivot_longer(c(-Na, -pH, -termcd), names_to="species", values_to="concentration")

p.titration <- ggplot(result.df.long, 
            aes(x=Na, y=pH)
) + 
  geom_path() +
  geom_point(aes(color=termcd)) +
  theme_minimal()

# p.titration

# PLOT MINIMIZATION OF EQUATIONS ####
sols <- 
  apply(result.m, MARGIN=1, FUN = function(sols){
    sols <- sols[1:6]
    
    sols <- 
      with(as.data.frame(t(sols)),{
        c(
          BM = H3A + H2A + HA + A - P_ca,
          # BM=safe.sum(H3A, + H2A, + HA, + A, - P_ca), 
          BC = H + Na - Kw/H - H2A - 2*HA - 3*A,
          # BC=safe.sum(H, +Na, -Kw/H, -H2A, -2*HA, -3*A), 
          A1 = H * H2A / Ka.1 - H3A,
          # A1=safe.eq(H3A,H2A,Ka.1,H),
          A2 = H * HA  / Ka.2 - H2A,
          # A2=safe.eq(H2A,HA,Ka.2,H),
          A3 = H * A   / Ka.3 - HA
          # A3=safe.eq(HA,A,Ka.3,H)
        )
      })
    
    # t(sols)
    sols
  })

sols.df <- as.data.frame(t(sols)) %>% mutate(Na= result.df$Na, termcd=result.df$termcd)

sols.df.long <- sols.df %>%
  pivot_longer(-c(Na, termcd), names_to="eqs", values_to="vals")

p.equations <- ggplot(sols.df.long, 
            aes(x=Na, y=vals, 
                color=eqs, group=eqs)
) + 
  geom_path() +
  # geom_point(aes(shape=termcd), color = "black") +
  theme_minimal()

# p.equations


# PLOT SPECIATION DIAGRAM ####
result.df <- as.data.frame(result.m)
result.df$termcd <- c("ok", "xtol", "stall", "maxit", "jacobad")[result.df$termcd]

result.df.long <- result.df %>%
  mutate(pH = -log10(H)) %>% select(-H, -termcd) %>%
  pivot_longer(c(-pH, -Na), names_to="species", values_to="concentration") #%>% filter(concentration > 10^-4)

p.speciation <- 
  ggplot(result.df.long, 
         aes(x=pH, y=concentration,
             color=species, group=species, 
             text=paste("Na:", Na))) +
  geom_path(alpha=1, size=1) +
  theme_minimal()

# p.speciation


# PLOT SPECIATION VS Na ####
result.df <- as.data.frame(result.m)
# result.df <- readRDS("output/result.df.0.04389.RDS")
result.df$termcd <- c("ok", "xtol", "stall", "maxit", "jacobad")[result.df$termcd]

result.df.long <- result.df %>%
  mutate(pH = -log10(H)) %>% select(-H, -termcd) %>%
  pivot_longer(c(-Na, -pH), names_to="species", values_to="concentration")

p.species_vs_Na <-
  ggplot(result.df.long, 
            aes(x=Na, y=concentration, 
                color=species, group=species,
                text=paste("pH:", pH))) + 
  geom_path() +
  theme_minimal()

# p.species_vs_Na


# PLOT ALL ####
# pdf("output/numerical.solution.pdf")
gridExtra::grid.arrange(p.titration,
                        p.equations,
                        p.speciation,
                        p.species_vs_Na,
                        ncol=2, nrow =2)
# dev.off()


