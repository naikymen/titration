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

H3A <- function(H, P.ca) P.ca * a1(H) / cd(H)
H2A <- function(H, P.ca) P.ca * a2(H) / cd(H)
HA <- function(H, P.ca) P.ca * a3(H) / cd(H)
A <- function(H, P.ca) P.ca * a4(H) / cd(H)

.H3A <- function(H) a1(H) / cd(H)
.H2A <- function(H) a2(H) / cd(H)
.HA <- function(H) a3(H) / cd(H)
.A <- function(H) a4(H) / cd(H)

OH <- function(H) Kw/H

Na.aprox <- function(H=10^-pH.seq, P.ca=1){
  H <- unname(H)
  P.ca <- unname(P.ca)
  list(
    H = H,
    Na = H2A(H, P.ca) + 2*HA(H, P.ca) + 3*A(H, P.ca) + OH(H) - H,
    H3A = H3A(H, P.ca),
    H2A = H2A(H, P.ca),
    HA = HA(H, P.ca),
    A = A(H, P.ca)
  )
}
#pH.seq <- seq(from=1,to=12,length.out = 1000)
#10^-pH.seq
# Na.aprox(H = 10^-7)

# Constantes
# P.vol=1
# P.ca=2
# P.mass <- P.ca*P.vol # moles
# H=10^-2
# Na.ca=1

# # Solución analítica con ajuste de volumen:
# 
# # [Na+] en funcion del titulante agregado (parte 1):
# V.tot = (P.vol + Na.vol)
# Na = Na.mass / V.tot
# ######
# Na = Na.mass / (P.vol + Na.vol)  # ok
# 
# 
# # [Na+] en funcion del titulante agregado (parte 2):
# Na = Na.mass / (P.vol + Na.vol)
# Na.mass = Na.vol * Na.ca
# ######
# Na = Na.vol * Na.ca / (P.vol + Na.vol)  # ok
# 
# 
# # Sacar "P.ca" como factor común de las ecuaciones de equilibrio:
# Na = H2A(H, P.ca) + 2*HA(H, P.ca) + 3*A(H, P.ca) + OH(H) - H
# #####
# Na = P.ca*.H2A(H) + 2*P.ca*.HA(H) + 3*P.ca*.A(H) + OH(H) - H
# #####
# Na = P.ca * (.H2A(H) + 2*.HA(H) + 3*.A(H)) + OH(H) - H
# 
# 
# # Reemplazar [Na+]:
# Na = Na.vol * Na.ca / (P.vol + Na.vol)
# Na = P.ca * (.H2A(H) + 2*.HA(H) + 3*.A(H)) + OH(H) - H
# #####
# Na.vol * Na.ca / (P.vol + Na.vol) = P.ca * (.H2A(H) + 2*.HA(H) + 3*.A(H)) + OH(H) - H
# 
# # "P.ca" en funcion del titulante agregado:
# P.ca = P.mass / V.tot
# V.tot = P.vol + Na.vol
# #####
# P.ca = P.mass / (P.vol + Na.vol)  # ok
# 
# # Reemplazar "P.ca":
# P.ca = P.mass / (P.vol + Na.vol)
# Na.vol * Na.ca / (P.vol + Na.vol) = P.ca * (.H2A(H) + 2*.HA(H) + 3*.A(H)) + OH(H) - H
# #####
# Na.vol * Na.ca / (P.vol + Na.vol) = (P.mass / (P.vol + Na.vol)) * (.H2A(H) + 2*.HA(H) + 3*.A(H)) + OH(H) - H
# ##### Multiplico por (P.vol + Na.vol)
# Na.vol * Na.ca = P.mass * (.H2A(H) + 2*.HA(H) + 3*.A(H)) + (P.vol + Na.vol) * (OH(H) - H)
# ##### Distribuyo (OH(H) - H)
# Na.vol * Na.ca = P.mass * (.H2A(H) + 2*.HA(H) + 3*.A(H)) + P.vol * (OH(H) - H) + Na.vol * (OH(H) - H) 
# #####
# Na.vol * Na.ca - Na.vol * (OH(H) - H) = P.mass * (.H2A(H) + 2*.HA(H) + 3*.A(H)) + P.vol * (OH(H) - H)
# 
# Luego despejamos Na.vol y listo!
# 
# Armamps una solucón analítica:
Na.adj.analitico <- function(H, P.mass, P.vol, Na.ca){
  Na.vol <- (P.mass * (.H2A(H) + 2*.HA(H) + 3*.A(H)) + P.vol * (OH(H) - H)) / (Na.ca - (OH(H) - H))
  return(Na.vol)
}

# Uno numérico: es divertido!
Na.adj.one <- function(H=10^-3, P.ca=1, P.vol=1, Na.ca=1, n.its=100){
  result <- matrix(ncol = 4, nrow = n.its)
  colnames(result) <- c("H", "Na", "Na.vol", "P.ca")
  
  P.mass <- P.ca*P.vol
  
  Na <- Na.aprox(H = H, P.ca = P.ca)[["Na"]]
  Na.mass <- Na * P.vol
  Na.vol <- Na.mass/Na.ca
  Vol <- P.vol + Na.vol
  P.ca.new <- P.mass/Vol
  
  result[1,] <- c(H, Na, Vol, P.ca.new)
  
  for(i in 2:n.its){
    
    Na <- Na.aprox(H = H, P.ca = P.ca.new)[["Na"]]
    Na.mass <- Na * Vol
    Na.vol <- Na.mass/Na.ca
    Vol <- P.vol + Na.vol
    P.ca.new <- P.mass/Vol
    
    result[i,] <- c(H, Na, Na.vol, P.ca.new)
    
    if(isTRUE(all.equal(result[i-1,], result[i,]))){
      break
      # print("got it!")
      # print(i)
    } else if(i==n.its){
      warning(paste("Calculation did not converge at pH: ", -log10(H)))
      return(rep(NA, 7))
    }
  }
  
  return(c(
    unlist(Na.aprox(H=H, P.ca=P.ca.new)),
    Na.vol=Na.vol
  ))
}

# Uno numérico: es divertido!
# P.vol=1
# P.ca=2
# P.mass <- P.ca*P.vol # moles
# Na.ca=1
# H=10^-1
# plot.convergence = T
Na.adj.numeric <- function(H, P.ca, P.vol, Na.ca, n.its=100, plot.convergence=F){
  H <- unname(H)
  P.ca <- unname(P.ca)
  P.vol <- unname(P.vol)
  Na.ca <- unname(Na.ca)
  
  P.mass <- P.ca * P.vol
  
  # Primero imaginamos que agregamos NaOH(s) y no hay incremento de volumen.
  # En ese caso:
  Na <- Na.aprox(H = H, P.ca = P.ca)[["Na"]]
  
  # Sin embargo, en el laboratorio hacemos la titulación con NaOH(l),
  # y, a diferencia de NaOH(s), eso cambia el pH de dos maneras:
  # 1) Por agregado de una base fuerte (el NaOH).
  # 2) Por dilución del ácido! (porque al agregar titulante acuoso se agrega agua, y eso aumenta volumen de la solución).
  
  # Imaginemos entonces que agregamos el NaOH(s) para llegar a pH=7,
  # ¡Pero que nos olvidamos de agregar el agua! ¿Cuánta agua olvidamos agregar?
  # Si teníamos Na.ca=1M, podemos calcular la masa de NaOH(s) 
  Na.mass <- Na * P.vol
  # y el volumen de NaOH(l) correspondiente:
  Na.vol <- Na.mass/Na.ca
  
  # Esta es nuestra apuesta respecto al volumen de titulante que debíamos agregar.
  
  # Guardar data de las iteraciones
  iter.data <- matrix(ncol=2, nrow=n.its+1, dimnames = list(NULL, c("P.ca", "Na.vol")))
  iter.data[1,] <- c(P.ca, Na.vol)
  
  for(i in 1:n.its){
    # Y luego el volumen final de la solución:
    P.vol.new <- P.vol + Na.vol # * 2
    
    # Imaginemos que es tarde y queremos irnos a casa, y tiramos una hipótesis:
    # "Agregar el agua no cambió nada! y me puedo ir a casa en paz".
    
    # Pero JTP nos dice:
    # "Cambiar el volumen diluye al ácido, así que chequealo antes de irte.".
    # "Podés fijarte si da igual agregar primero el agua y después el NaOH(s)".
    
    # Entonces primero diluimos el ácido:
    P.ca.new <- P.mass / P.vol.new
    
    # Ahora nos fijamos la concentración de NaOH requerida para llegar al mismo pH
    # en la solución de ácido diluída:
    Na.new <- Na.aprox(H = H, P.ca = P.ca.new)[["Na"]]
    # Y la masa y volumen correspondiente de titulante:
    Na.mass.new <- Na.new * P.vol.new
    Na.vol.new <- Na.mass.new/Na.ca
    
    # Si agregar el volumen primero no cambiaba nada,
    # entonces nuestra apuesta inicial para la cantidad de titulante
    # debería ser buena! e igual a la que calculamos haciendo el TP al revés.
    # ¿Lo es?
    
    # Guardar data de las iteraciones
    iter.data[i+1,] <- c(P.ca.new, Na.vol.new)
    
    if(!isTRUE(all.equal(Na.vol, Na.vol.new))){
      # Si no son iguales... ¿qué significa?
      # En ese caso la dilución del ácido cambia el volumen de titulante
      # necesario para llevar la solución a cierto pH (o sea, JTP tenía razón).
      
      # El objetivo es encontrar una dilución del ácido, tal que, a cierto pH,
      # el volumen de NaOH 1M necesario coincida con el volumen que diluyó al ácido.
      
      # Si nuestra aproximación inicial no fue suficiente:
      # ¿qué valor de Na.vol probamos a continuación?
      
      # Primero notamos que "Na.vol.new" no puede ser más grande que "Na.vol".
      # Esto es porque si el agua diluye al ácido, en ninguna circunstancia
      # voy a tener que agregar más NaOH para llegar al mismo pH.
      # "Es como que la dilución por sí misma puede titular" (?)
      
      # Ahora sabemos que "Na.vol.new" tiene cota superior.
      # ¿Pero que pasa si adivinamos un valor muy chico para probar?
      # En ese caso el volumen agregado sería pequeño y el volumen total
      # se acercaría al de la solución original, que es el lugar de donde partimos.
      # O sea, si probamos con un valor muy bajo, el ácido está concentrado,
      # y la corrección de volumen va a dar un "Na.vol.new" alto otra vez.
      
      # A partir de esto podemos decir que hay una tendencia hacia el centro
      # la pregunta que queda es ¿converge realmente?
      
      # print(Na.vol.new)
      Na.vol <- Na.vol.new
    } else {

      # Graficar el proceso
      if(plot.convergence){
        
        iter.data.ok <- iter.data[!is.na(iter.data[,1]),]
        
        plot(x=iter.data[,1], y = iter.data[,2], type = "l",
             xlim = range(iter.data[,1], na.rm = T),
             ylim=range(iter.data[,2], na.rm = T),
             xlab="P.ca.new", ylab="Na.vol.new")
        points(x=iter.data[1,1], y = iter.data[1,2], col = "red")
        points(x=iter.data[-1,1], y = iter.data[-1,2], col = "blue")
        title(main = "Apuesta inicial (rojo) e iteraciones (azul)")
        
        # plot(iter.data.ok[,1], type = "l",
        #      ylim=range(iter.data.ok, na.rm = T))
        # lines(iter.data.ok[,2])
        # points(iter.data[,1], col = "red")
        # points(iter.data[,2], col = "blue")
        # title(main = "Apuesta inicial (rojo) e iteraciones (azul)")
      }
      
      # Si son iguales, podemos irnos a casa :)
      Na.aprox.new <- Na.aprox(H = H, P.ca = P.ca.new)
      return(c(
        unlist(Na.aprox.new), 
        Vol=Na.vol.new
        ))
    }
    
    if(i==n.its){
      stop(paste("Calculation did not converge at pH: ", -log10(H)))
    }
  }
  
}

# Ejemplo:
# 
# Na.adj.result <- Na.adj.one()
# 
# Na.adj.result
# 
# unlist(Na.aprox(H = Na.adj.result["H"], P_ca = Na.adj.result["P.ca"]))

# Numeric (old function)
Na.adj <- function(pH.seq, P.ca, P.vol, Na.ca, n.its){
  res.columns <- c("H", "Na", "H3A", "H2A", "HA", "A", "Na.vol")
  res.matrix <- matrix(nrow = length(pH.seq), ncol = length(res.columns))
  colnames(res.matrix) <- res.columns
  
  for(i in seq_along(pH.seq)){
    h <- 10^(-pH.seq[i])
    # res.matrix[i,] <- unlist(Na.adj.numeric(H=h, P.ca, P.vol, Na.ca, n.its = n.its))
    res.matrix[i,] <- unlist(Na.adj.one(H=h, P.ca, P.vol, Na.ca))
  }
  
  return(as.data.frame(res.matrix))
}

# Analytic
Na.adj2 <- function(pH.seq, P.mass, P.vol, Na.ca){
  res.columns <- c("H", "Na", "H3A", "H2A", "HA", "A", "Na.vol")
  res.matrix <- matrix(nrow = length(pH.seq), ncol = length(res.columns))
  colnames(res.matrix) <- res.columns
  
  for(i in seq_along(pH.seq)){
    h <- 10^(-pH.seq[i])
    Na.vol <- Na.adj.analitico(H=h, P.mass, P.vol, Na.ca)
    res.matrix[i,"Na.vol"] <- Na.vol
    
    P.ca.new <- P.mass / (P.vol + Na.vol)
    especies <- unlist(Na.aprox(H=h, P.ca = P.ca.new))
    res.matrix[i,-7] <- especies
  }
  
  return(as.data.frame(res.matrix))
}

# Ejemplos:
if(F){
  # Diagrama de Especiacion (sin ajustar volumen)
  P.vol=1
  P.ca=2
  P.mass <- P.ca*P.vol # moles
  H=10^-2
  Na.ca=1
  pH.seq <- seq(from=0, to=14, length.out = 1000)
  Species.seq <- Na.aprox(H=10^-pH.seq, P.ca=P.ca)
  
  plot(x=-log10(Species.seq$H), y = Species.seq$H3A/P.ca, 
       main = "Diagrama de especiación del ácido fosfórico",
       xlab = "pH", ylab = "alpha",
       xlim = c(0,14), ylim=c(0,1), type = "l", col=1)
  lines(x=-log10(Species.seq$H), y = Species.seq$H2A/P.ca, col=2)
  lines(x=-log10(Species.seq$H), y = Species.seq$HA/P.ca, col=3)
  lines(x=-log10(Species.seq$H), y = Species.seq$A/P.ca, col=4)
  legend(x = "right", box.lty = 0, bg="transparent",
         legend = c("H3A","H2A","HA","A"), col = 1:4, lty = rep(1,4))
  
  # Titulacion (sin ajustar, ajuste analítico y ajuste numérico)
  P.vol=1
  P.ca=2
  P.mass <- P.ca*P.vol # moles
  Na.ca=1
  pH.seq <- seq(from=1, to=13.25, length.out = 1000)
  
  # Experiment conditions:
  # P.vol=0.025    # 25 mL
  # P.ca=0.04389   # 0.04389 M
  # P.mass <- P.ca*P.vol # moles
  # Na.ca=0.09948  # 0.09948
  # pH.seq <- seq(from=1, to=12, length.out = 1000)
  
  # Todo en volumen de [NaOH]
  Na.seq      <- (Na.aprox(H=10^-pH.seq, P.ca=P.ca)$Na * P.vol) / Na.ca
  Na.seq.adj  <- Na.adj(pH.seq=pH.seq, P.ca = P.ca, P.vol = P.vol, Na.ca = Na.ca, n.its = 200)$Na.vol
  Na.seq.adj2 <- Na.adj2(pH.seq=pH.seq, P.mass = P.mass, P.vol = P.vol, Na.ca = Na.ca)$Na.vol
  
  # Plot
  plot(x = Na.seq, y = pH.seq, type="l", col=2, xlab = "Na.vol", ylab = "pH", main = "Titration curves")
  lines(x = Na.seq.adj, y = pH.seq, col=3)
  lines(x = Na.seq.adj2,     y = pH.seq, col=4)
  legend("right", legend = c("Original", "Adjusted (numeric)", "Adjusted (analytic)"),
         box.lty = 0, bg="transparent",
         lty = rep(1,3),
         col = 2:4)
}

