rm(list = ls())

#equations from Munz, P. et al. (2009) ‘Zombies!’, Science.
#parameters (greek letters anglicised)
P <- 0 # birth rate
d <- 0.00001 #natural death rate percent per day
B <- 0.0095 #transmission percent per day
G <- 0.0001 #resurrect percent per day
A <- 0.0001 #destroy percent per day


#install.packages("odeintr")
library("odeintr")

f <- function(y,t){
  #start with values in y[0:2]
Si <- y[1] #initial susceptible amount (for this iteration)
Zi <- y[2] #initial zombie amount
Ri <- y[3] #initial resurrected amount
#equations from Munz, P. et al. (2009) ‘Zombies!’, Science.
#S′ <- Π−βSZ−δS
Srate <- P - B*Si*Zi - d*Si
#Z′ <- βSZ+ζR−αSZ
Zrate <- B*Si*Zi + G*Ri - A*Si*Zi
#R′ <- δS+αSZ−ζR.
Rrate <- d*Si + A*Si*Zi - G*Ri
return(c(Srate, Zrate, Rrate))
}

#intialise
S0 <- 500 #initial susceptibles
Z0 <- 0 #initial zombies
R0 <- 0 #initial death population

y0 <- c(S0,Z0,R0)

soln = integrate_sys(f, y0, 5, 0.01)

plot(soln)

# dxdt = function(x, t) x * (1 - x)
# 
# integrate_sys(dxdt, 0.001, 15, 0.01)
# 
# integrate_sys(dxdt, 0.01, 15, 0.01)
# 
# dxdt(5,2)
