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

Na <- function(H=10^-pH.seq, P_ca=1)list(
  H = H,
  Na = H2A(H, P_ca) + 2*HA(H, P_ca) + 3*A(H, P_ca) + OH(H) - H,
  H3A = H3A(H, P_ca),
  H2A = H2A(H, P_ca),
  HA = HA(H, P_ca),
  A = A(H, P_ca)
)

Na.adj.one <- function(H=10^-3, P.ca=1, P.vol=1, Na.ca=1, n.its=1:100){
  result <- matrix(ncol = 4, nrow = length(n.its))
  colnames(result) <- c("H", "Na", "Vol", "P.ca")
  
  P.mass <- P.ca*P.vol
  
  for(i in seq_along(n.its)){
    
    Na = H2A(H, P.ca) + 2*HA(H, P.ca) + 3*A(H, P.ca) + OH(H) - H
    Vol <- P.vol + Na/Na.ca
    P.ca <- P.mass/Vol
    
    result[i,] <- c(H, Na, Vol, P.ca)
    
    if(isTRUE(all.equal(result[i-1,], result[i,]))){
      break
    } else if(i==max(seq_along(n.its))){
      print("pH: ")
      print(-log10(H))
      stop("Calculation did not converge")
    }
  }
  # d <- as.data.frame(result[!is.na(result[,1]),])
  # plot(d$Na)
  # plot(d$Vol)
  # plot(d$P.ca)
  
  res.adj <- Na(H = H, P_ca=P.mass/Vol)
  
  return(res.adj)
}

Na.adj <- function(H, ...){
  res.matrix <- matrix(nrow = length(pH.seq), ncol = 6)
  colnames(res.matrix) <- c("H", "Na", "H3A", "H2A", "HA", "A")
  
  for(i in seq_along(pH.seq)){
    res.matrix[i,] <- unlist(Na.adj.one(H=H[i], ...))
  }
  
  return(as.data.frame(res.matrix))
}

# Uncomment to compute
# pH.seq <- seq(from=1, to=14, length.out = 1000)
# Na.seq <- Na(H=10^-pH.seq, P_ca=1)
# Na.seq.adj <- Na.adj(H=10^-pH.seq, P.ca=1)
# 
# plot(Na.seq.adj$Na, pH.seq, type="l", xlim = c(0,4), col=2, 
#      xlab = "Na", ylab = "pH", main = "Titration curves")
# lines(Na.seq$Na, pH.seq, col=3)
# legend("right", legend = c("Adjusted", "Original"),
#        box.lty = 0, bg="transparent",
#        lty = rep(1,2),
#        col = c(2,3))
