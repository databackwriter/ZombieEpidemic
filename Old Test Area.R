
#####old test area#####
S.params=config.Susceptible
Z.params=config.Zombie
outputfilename=config.plot$outputfilename
# ){

Ni=S.params["n"]
xinit<-S.params["xinit"]
yinit<-S.params["yinit"]
Xt<-rep(xinit,Ni) #set of start points (one for each humanoid)
Yt<-rep(yinit,Ni) #set of start points (one for each humanoid)
Z<-rwalk2dhumanoids(
  x=S.params,
  Ni,
  Xt,
  Yt,
  humanoidtype=prefix.Susceptible,
  startind = 1)

Ni=Z.params["n"]
xinit<-Z.params["xinit"]
yinit<-Z.params["yinit"]
Xt<-rep(xinit,Ni) #set of start points (one for each humanoid)
Yt<-rep(yinit,Ni) #set of start points (one for each humanoid)
startind<-length(Z[,1,1])+1
Z1<-rwalk2dhumanoids(
  x=Z.params,
  Ni,
  Xt,
  Yt,
  humanoidtype=prefix.Zombie,
  startind=startind)

#combine them to get a humanoid array
Z<-arrayunion(Z,Z1)
Z<-plotoutbreakXYTasgif(Z,
                        S.params,
                        Z.params,
                        unleashplot = TRUE)

