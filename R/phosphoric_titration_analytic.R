# Sumamente facil ####

Ka.1 <- 7.1 * 10^-3
Ka.2 <- 6.3 * 10^-8
Ka.3 <- 4.5 * 10^-13
Kw <- 10^-14

a1 <- function(H) H^3
a2 <- function(H) H^2 * Ka.1
a3 <- function(H) H * Ka.1 * Ka.2
a4 <- function(H) Ka.1 * Ka.2 * Ka.3
cd <- function(H) a1(H) + a2(H) + a3(H) + a4(H)

H3A <- function(H, P_ca) P_ca * a1(H) / cd(H)
H2A <- function(H, P_ca) P_ca * a2(H) / cd(H)
HA <- function(H, P_ca) P_ca * a3(H) / cd(H)
A <- function(H, P_ca) P_ca * a4(H) / cd(H)

OH <- function(H) Kw/H

Na <- function(H=10^-pH.seq, P_ca=1)
  H2A(H, P_ca) + 2*HA(H, P_ca) + 3*A(H, P_ca) + OH(H) - H

# Uncomment to compute
# pH.seq <- seq(from=0,to=14,length.out = 1000)
# Na.seq <- Na(H = 10^-pH.seq, P_ca=1)
# plot(Na.seq, pH.seq)
