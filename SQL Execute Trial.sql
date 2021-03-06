﻿DECLARE @BigRs NVARCHAR(MAX)= '
#####global constants
prefix.Susceptible<-"Susceptible"
prefix.Zombie<-"Zombie"
prefix.Removed<-"Removed"
prefix.Newborn<-"ProtoLiving"
prefix.Resurrected<-"BackFromDead"
prime.SusceptibleWins<-3
prime.ZombieWins<-2
prime.Resurrection<-5
prime.NaturalDeath<-7
#####configuration area####
configuration.SRZModel<-function(P=0,# birth rate
                                 d=0.0001,# natural death percent (per unit time)
                                 B=0.0095,# transmission percent  (per unit time)
                                 G=0.0001,# resurrect percent (per unit time)
                                 a=0.0001,# destroy
                                 rcritical=2,# critical Euclidean distance below which a zombie and human will go into battle
                                 zombiewinratio=0.95# ratio at which Zombie wins
){
  x<-c(P,d,B,G,a,rcritical,zombiewinratio)
  names(x)<-c("P","d","B","G","a","rcritical","zombiewinratio")
  return(x)
}
configuration.Humanoid<-function(xinit,#set a starting co-ordinate for humanoids on the x-dimension
                                 yinit,#set a starting co-ordinate for humanoids on the y-dimension
                                 n,#initial number of humanoids
                                 mu,#mean footstep length for humanoids IMPORTANT: because we have unit time this is average humanoids speed
                                 sigma#standard deviation of footstep lengths for humanoids
){
  x<-c(xinit,
       yinit,
       n,
       mu,
       sigma)
  names(x)<-c("xinit",
              "yinit",
              "n",
              "mu",
              "sigma")
  return(x)
}
config.plot.global<-function(outputfilename="/Users/peterfmoore/Documents/GitHub/ZombieEpidemic/zombieanimation.gif"){
  #define axes limits on the whole data set
  xmin<-min(config.Susceptible["xinit"], config.Zombie["xinit"])
  ymin<-min(config.Susceptible["yinit"], config.Zombie["yinit"])
  xmax<-max(config.Susceptible["xinit"], config.Zombie["xinit"])
  ymax<-max(config.Susceptible["yinit"], config.Zombie["yinit"])
  xlimlocal<-c(xmin-1*floor(sqrt(time.ntimesteps+1)),xmax+ceiling(sqrt(time.ntimesteps+1))) #using square root of n as arough approximatei
  ylimlocal<-c(ymin-1*floor(sqrt(time.ntimesteps+1)),ymax+ceiling(sqrt(time.ntimesteps+1)))
  x<-list(xlimlocal,ylimlocal,outputfilename)
  names(x)<-c("xlimlocal","ylimlocal","outputfilename")
  return(x)
}

#####model area####
model.timestep.assessment<-function(Scount,Zcount){
  breakme<-FALSE
  if (Scount==0){
    breakme<-TRUE
  } else if (Zcount==0){
    breakme<-TRUE
  }
  return(breakme)
}
model.run.params<-function(
  S.params,
  Z.params,
  outputfilename=config.plot$outputfilename,
  mode="matrix",
  iterations=c(1:1)
){
  #runs a model of zombies atacking humans and vice versa, has two modes:
  #1. matrix: the default, produces a matrix of Zombie and Susceptible ubers at each timestep
  #2. gif: maps the same information in a gif animation
  if (mode=="matrix"){
    
    p<-plot(x=c(),y=c(),type="l", col="blue",xlab="time", ylab="Humanoid Count", main="Zombies vs Susceptibles",
            frame.plot = FALSE,xlim=c(0,length(time.vtimesteps)),ylim=c(0,(S.params["n"]+Z.params["n"])))
    p<-legend("bottom", c("Susceptibles", "Zombies"), xpd = TRUE, horiz = TRUE, 
              inset = c(0, 0), pch=15,bty = "n", col = c("blue","red"), cex = 0.5)
  }
  for(mi in iterations){ 
   #####setup####
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
    #####iterate#####
    m.M<- matrix(data=NA,nrow = length(time.vtimesteps),ncol=2,dimnames = list(paste0("t",time.vtimesteps),c("SusceptibleCount", "ZombieCount")))
    m.M["t0",]<-c(S.params["n"],Z.params["n"])
    if (mode=="matrix"){
      for (i in time.vtimesteps){
        j<-i+1
        
        #for each timestep we need to remove any humanoids that get too near, too born, or too tired (i.e. they die)
        Z<-reconfigurepopulation(Z=Z,
                                 S.params = S.params,
                                 Z.params = Z.params,
                                 timesteplocal=j)
        #look at just the suceptibles
        m.Susceptible<-arrayStartWith(Z,prefix.Susceptible)
        #look at just the zombies
        m.Zombie<-arrayStartWith(Z,prefix.Zombie)
        m.Susceptible.count<-length(rownames(m.Susceptible))
        m.Zombie.count<-length(rownames(m.Zombie))
        m.M[j,]<-c(m.Susceptible.count,m.Zombie.count)
        x<-(paste0("Final timestep recorded at t=t",j))
        if(model.timestep.assessment(m.Susceptible.count,m.Zombie.count)){
          break #if there are no zombies or suscpetibles left then leave
        }
      }
      print(x)
      time.axis<-time.vtimesteps[1:length(time.vtimesteps)]
      m.M<-subset(m.M,rowSums(m.M)>=0)
      p<-lines(x=time.axis,y=m.M[,"SusceptibleCount"],type="l", col="blue")
      p<-lines(x=time.axis,y=m.M[,"ZombieCount"],type="l", col="red")
  
      } else if (mode=="gif"){
        Z<-plotoutbreakXYTasgif(Z,
                                S.params,
                                Z.params,
                                unleashplot = TRUE)
      }
  }
  return(m.M)
}

#####random walk area####
rwalk1dhumanoids<-function(Ni, mu, sigma, Xt,humanoidtype="Humanoid",startind=1,humanoidnames=c(),offset=0){
  #builds a one-d random walk along  nominal time axis
  #where mu=0 this is akin to Brownian motion
  #inspired by ideas at http://www.phytools.org/eqg/Exercise_4.1/
  #time.ntimesteps-number of timesteps, Ni-size of initial Humaoid population,
  #mu,sigma-mean and standard deviation for footstep length
  #xinit-initial x-position
  #time.vtimesteps=number of timestepos
  #humanoidtype: typically human or zombie
  #startind-starting ordinal for zombies
  Nlarge<-time.ntimesteps*Ni # number of time steps multiplied by initial number of humanoids
  #Xt<-rep(xinit,Ni) #set of start points (one for each humanoid)
  footsteps<-rnorm(n=Nlarge,mean=mu,sd=sigma) #lets now get a super long set of footsteps
  X<-matrix(data=footsteps,nrow=Ni,ncol=time.ntimesteps) #now get a matrix where we get time.ntimesteps columns for each of our Ni humanoids
  X<-cbind(Xt,X) #now place our initial co-ordinate against x (literally x vs t)
  X<-apply(X=X,MARGIN=1,FUN=cumsum)
  X<-t(X) #SUBTLETY: transpose
  if(length(humanoidnames)==0){
    humanoidnames<-paste(humanoidtype,c(startind:(Ni+startind-1)))
  }
  dimnames(X)=list(humanoidnames,paste0("t", time.vtimesteps+offset))
  return(X)
}
rwalk2dhumanoids<-function(x,Ni, Xt, Yt, humanoidtype="Humanoid",startind=1,humanoidnames=c(),offset=0){
  mu<-x["mu"]
  sigma<-x["sigma"]
  # xinit<-x["xinit"]
  # yinit<-x["yinit"]
  # Xt<-rep(xinit,Ni) #set of start points (one for each humanoid)
  # Yt<-rep(yinit,Ni) #set of start points (one for each humanoid)
  X<-rwalk1dhumanoids( Ni=Ni, mu=mu, sigma=sigma, Xt=Xt,humanoidtype=humanoidtype,startind=startind,humanoidnames=humanoidnames,offset=offset) # get Ni rows of humanoids X co-ordinates
  Y<-rwalk1dhumanoids(Ni=Ni, mu=mu, sigma=sigma, Xt=Yt,humanoidtype=humanoidtype,startind=startind,humanoidnames=humanoidnames,offset=offset) # get Ni rows of humanoids Y co-ordinates
  Z<-cbind(X,Y)
  Zdim<-c(Ni,time.ntimesteps+1,2)#humanoid,time,x-y
  humanoidnames<-rownames(X)#paste(humanoidtype,c(1:Ni))
  W<-array(data= Z,
           dim = Zdim,
           dimnames=list(humanoidnames,paste0("t", time.vtimesteps+offset),c("x-displacement","y-displacement"))
           )
  return(W)
}
#####linear algebra area####
arrayunion<-function(H,J, topology="NTXY"){
  #takes two arrays where expected dimensions are planar co-ordinates in time (for humanoids)
  #these are then deconstructed into the X and Y co-ordinates
  #and reconstructed
  #the result is a union in the n-dimensional vector space

  Hx<-H[,,1]
  Hy<-H[,,2]
  Jx<-J[,,1]
  Jy<-J[,,2]

  xcoords<-rbind(Hx,Jx)
  ycoords<-rbind(Hy,Jy)
  allcoords<-cbind(xcoords,ycoords)

  Nih<-length(H[,1,1])
  Nij<-length(J[,1,1])
  Zdim<-c(Nih+Nij,time.ntimesteps+1,2)#humanoid,time,x-y
  Z<-array(data= allcoords,
           dim = Zdim,
           dimnames=list(rownames(allcoords),paste0("t", time.vtimesteps),c("x-displacement","y-displacement"))
  )
  return(Z)
}
arrayChangeRowNames<-function(Z,oldnames,newnames){
  #takes an array, Z, and replaces a vector of oldnames with a vector of new ones
  myrownames<-rownames(Z)
  for (i in 1:length(oldnames)){ #blimey, there has to be a more efficient way than this
    myrownames<-replace(myrownames, myrownames == oldnames[i], newnames[i])
  }
  rownames(Z)<-myrownames
  return(Z)
}
arrayStartWith<-function(Z,prefix){
  nm<-rownames(Z)
  startwithlocal<-paste("^",prefix,sep="")
  vstartwith<- nm %in% grep(startwithlocal, nm, value = TRUE)
  nm<-subset(nm,vstartwith)
  Z<-Z[nm,1:(time.ntimesteps+1),1:2,drop=FALSE]

  return(Z)
}
euclideandistancebetweentwomatrices<-function(H,J){
  #function that takes in two matrices EXPECTING col 1 and col2 two be pairwise co-ordinates in Euclidean space
  #note also the importance of the naming here, the function itself is humanoid agnostic, this information is contained in
  #the meta data of the matrix (i.e. the row and column names)
  #these are then returned and this is non-trivial
  library(pdist)
  dists<-as.matrix(pdist(H,J))
  distsnames<-list(rownames(H), rownames(J))
  dimnames(dists)<-distsnames
  return(dists)
}
isolaterowminima<-function(M, removeNA = TRUE){
  #takes a matrix and returns the values that match the minima for that row and leaves the rest as NA
  library(matrixStats)
  x<-rowMins(M,na.rm = removeNA)
  L<-ifelse(M==x,M,NA)
  return(L)
}
isolatecolminima<-function(M, removeNA = TRUE){
  #takes a matrix and returns the values that match the minima for that column and leaves the rest as NA
  return(t(isolaterowminima(t(M),removeNA))) #note that we transpose twice
}
isolateminima<-function(M, removeNA = TRUE){
  #takes a matrix and returns the values that match the minima for that row and column and leaves the rest as NA
  return(isolatecolminima(isolaterowminima(M, removeNA), removeNA))
}
matrixnn<-function(K){
  booGood<-FALSE
  if (!is.null(rownames(K))){
    booGood<-TRUE
  }
  return(booGood)
}
#####repopulation area####
deathmatrix<-function(Z,timestep){
  # a matrix of zombie winners (2s), zombie losers (3s), and live to fight another day-ers (NAs)
  # NB it is important that these numbers are prime.SusceptibleWins
  interfacedistance = config.model["rcritical"]
  zombiewinrate = config.model["zombiewinratio"]
  # naturaldeathrate=config.model["d"]

  H<-arrayStartWith(Z,prefix.Susceptible)#Susceptibles
  J<-arrayStartWith(Z,prefix.Zombie)#Zombies
  xyfootprintH<-as.matrix(H[,timestep,])
  xyfootprintJ<-as.matrix(J[,timestep,])

  #nearest neighbour killings
  pdistHJ<-euclideandistancebetweentwomatrices(xyfootprintH, xyfootprintJ)
  L<-isolateminima(pdistHJ) #isolate minima (ASSUMPTION only nearest neighbours do anything)
  x<-ifelse(runif(length(L))>zombiewinrate, prime.SusceptibleWins, prime.ZombieWins) #CONVENTION 2 means zombie wins, 3 means human wins (zombie doesnt), NA means do nothing
  K<-ifelse(L<interfacedistance,x,NA) #create a matrix of winners, losers, and live to fight another day-ers


  # #natural causes killings
  # y<-ifelse(runif(length(K))<naturaldeathrate,prime.NaturalDeath,NA)
  # K<-ifelse(is.na(K),y,K)
  return(K)
}
deatmatrixisolatehumanoids<-function(deathmatrixlocal,deathprime){
  dmatrowsums<-rowSums(deathmatrixlocal, na.rm=TRUE)
  namelist<-(subset(dmatrowsums,dmatrowsums==deathprime))
  namelist<-names(namelist)

  if (length(subset(namelist,namelist=="x-displacement"))+length(subset(namelist,namelist=="y-displacement"))>0){
    namelist<-c()
  }

  return(namelist)
}
reconfigurepopulation.changenames<-function(Z,namelist,prefixfind,prefixreplace){

  deaths<-sort(unique(namelist))
  oldnames<-sort(unique(deaths))
  prefix<-paste("^", prefixfind,sep="")
  newnames<-sort(sub(prefix, prefixreplace, oldnames))
  Z<-arrayChangeRowNames(Z=Z,oldnames=oldnames,newnames=newnames)
  return(Z)
}
reconfigurepopulation.MCSnames<-function(Z,prefix,r,mode="names"){
  #returns the names of humanoids whose status is aboutto change (e.g. Susceptible->Dead)
  #mode="names" means return a set of names, "numbers" means return lenght
  K<-arrayStartWith(Z,prefix=prefix)
  Knames<-rownames(K)
  u<-runif(1:length(Knames))
  u<-ifelse(u<r,1,0)
  names(u)<-Knames
  x<-subset(u,u==1)
  x<-names(x)
  if(mode=="numbers"){
    x<-length(x)
  }
  return(x)
}
reconfigurepopulation.introduce<-function(Z,S.params,Ni,Xt,Yt,startind,humanoidtype,humanoidnames){
  Z1<-rwalk2dhumanoids(
    x=S.params,
    Ni,
    Xt,
    Yt,
    humanoidtype,
    startind,
    humanoidnames)
  unames<-union(rownames(Z),rownames(Z1))
  Z<-arrayunion(Z,Z1)
  rownames(Z)<-unames

  return(Z)
}
reconfigurepopulation.revive<-function(Z,Zparams,Ni,Xt,Yt,startind,humanoidtype,humanoidnames,offset,timesteplocal){
  K<-rwalk2dhumanoids(
    x=Zparams,
    Ni=Ni,
    Xt=Xt,
    Yt=Yt,
    humanoidtype=humanoidtype,
    startind,
    humanoidnames=humanoidnames,
    offset=offset)
  
  rownames(K)<-humanoidnames
  Knames<-rownames(K)
  vt<-paste0("t",(timesteplocal-1):(length(time.vtimesteps)-1))
  
  for (i in (1:2)){
    Z[Knames,vt,i]<-K[Knames,vt,i]
  }
  return(Z)
}

exhumeremoved<-function(Z,timesteplocal,prefix,replacement,mode="stun"){
  #when humanoids are removed something happens to their position
  #the mode handles this:
  #mode="stun"=>dead humanoids stay still so they get stuck with the x-y of their last known whereabouts for eternity
  #mode="start"=>newborns already have their co-ordinates so we dont need to do anything
  #mode="revive"=>resurrected were previously stunned so they need a new set of coords
  K<-arrayStartWith(Z,prefix)
  if(matrixnn(K)){
    Knames<-rownames(K)
    if (mode=="stun"){
      for (i in (1:2)){
            Z[Knames,,i]<-K[,timesteplocal-1,i]#SUBTLETY: timestep-1
        }
    }
    prefix<-paste(prefix,".",sep="")
    pm<-sub(prefix,replacement,Knames)
    Z<-arrayChangeRowNames(Z,Knames,pm)
  }
  return(Z)
}
reconfigurepopulation<-function(Z,
                                S.params,
                                Z.params,
                        timesteplocal){

    #anything with the name RemovedX was removed last time
    #now it must join the official removed class
    #this means that the name is changed to remove and the co-ordinates are locked in time
    Z<-exhumeremoved(Z,timesteplocal,prefix="RemovedX[ZD]",replacement="Removed ")
    Z<-exhumeremoved(Z,timesteplocal,prefix="BackFromDeadXR",replacement="Zombie ")
    Z<-exhumeremoved(Z,timesteplocal,prefix="RemovedX[S]",replacement="Zombie ", mode="revive")
    Z<-exhumeremoved(Z,timesteplocal,prefix=prefix.Newborn,replacement=prefix.Susceptible, mode="start")

    #start with the resurrections
    introducednames.Zombies.Resurrected<-reconfigurepopulation.MCSnames(Z=Z,prefix=paste(prefix.Removed,""),r=config.model["G"])
    if(length(introducednames.Zombies.Resurrected)>0){
      Ni<-length(introducednames.Zombies.Resurrected)
      Xt<-Z[introducednames.Zombies.Resurrected,timesteplocal,1]
      Yt<-Z[introducednames.Zombies.Resurrected,timesteplocal,2]
      startind<-1
      humanoidtype=prefix.Resurrected
      humanoidnames<-introducednames.Zombies.Resurrected #build a new random walk for these guys but update dont insert
      offset<-max(timesteplocal-1,0)
      Z<-reconfigurepopulation.revive(Z=Z,Zparams=Z.params,Ni=Ni,Xt=Xt,Yt=Yt,startind=startind,humanoidtype=humanoidtype,humanoidnames=humanoidnames,offset=offset,timesteplocal=timesteplocal)
      Z<-reconfigurepopulation.changenames(Z=Z,namelist=introducednames.Zombies.Resurrected,prefixfind=prefix.Removed,prefixreplace="BackFromDeadXR")
    }
    
    #death by contact (nearest neighbour deaths)
    deathmatrixlocal<-deathmatrix(Z,timestep=timesteplocal)
    s<-sum(deathmatrixlocal,na.rm =TRUE)
    if(s > 0){
      #rename the removed class with a placeholder name of RemovedX
      deathnames.Susceptibles<-deatmatrixisolatehumanoids(deathmatrixlocal,prime.ZombieWins)
      # print(deathnames.Susceptibles)
      if (length(deathnames.Susceptibles)>0){
        # print(deathnames.Susceptibles)
        # Ni<-length(deathnames.Susceptibles)
        # Xt<-Z[deathnames.Susceptibles,timesteplocal,1]
        # Yt<-Z[deathnames.Susceptibles,timesteplocal,2]
        # startind<-1
        # humanoidtype="Removed"
        # humanoidnames<-deathnames.Susceptibles #build a new random walk for these guys but update dont insert
        # offset<-max(timesteplocal-1,0)
        # Z<-reconfigurepopulation.revive(Z=Z,Zparams=Z.params,Ni=Ni,Xt=Xt,Yt=Yt,startind=startind,humanoidtype=humanoidtype,humanoidnames=deathnames.Susceptibles,offset=offset,timesteplocal=timesteplocal)
        Z<-reconfigurepopulation.changenames(Z=Z,namelist=deathnames.Susceptibles,prefixfind=prefix.Susceptible,prefixreplace="RemovedXS")
      }
      
      deathnames.Zombies<-deatmatrixisolatehumanoids(t(deathmatrixlocal),prime.SusceptibleWins)
      if(length(deathnames.Zombies)>0){
        Z<-reconfigurepopulation.changenames(Z=Z,namelist=deathnames.Zombies,prefixfind=prefix.Zombie,prefixreplace="RemovedXZ")
        }
      }

    #nautral deaths
    deathnames.Susceptibles.Natural<-reconfigurepopulation.MCSnames(Z=Z,prefix=prefix.Susceptible,r=config.model["d"])
    if(length(deathnames.Susceptibles.Natural)>0){
     Z<-reconfigurepopulation.changenames(Z=Z,namelist=deathnames.Susceptibles.Natural,prefixfind=prefix.Susceptible,prefixreplace="RemovedXD")
     }
    
    #natural births
    introducednames.Susceptibles.Natural<-reconfigurepopulation.MCSnames(Z=Z,prefix=prefix.Susceptible,r=config.model["P"],mode="numbers")
    introducednames.Susceptibles.Natural.Names<-reconfigurepopulation.MCSnames(Z=Z,prefix=prefix.Susceptible,r=config.model["P"])
    if(introducednames.Susceptibles.Natural>0){#SUBTLETY: length function omitted
      Ni<-introducednames.Susceptibles.Natural
      xinit<-S.params["xinit"]
      yinit<-S.params["yinit"]
      Xt<-rep(xinit,Ni) #set of start points (one for each humanoid)
      Yt<-rep(yinit,Ni) #set of start points (one for each humanoid)
      startind<-length(Z[,1,1])+1
      humanoidtype<-prefix.Newborn
      humanoidnames<-introducednames.Susceptibles.Natural.Names
      Z<-reconfigurepopulation.introduce(Z,S.params,Ni,Xt,Yt,startind,humanoidtype, humanoidnames)
    }
 


     return(Z)
}
#####plot area####
plotsub<-function(plotmatrix, timestep, singlecolour){
  nameextent<-1:length(rownames(plotmatrix))
  xyfootprint<-(plotmatrix[nameextent,timestep,1:2,drop=FALSE])
 #get ourselves a nice two-D matrix--,1,,drop=FALSE
  xfootprint<-xyfootprint[,,1] #a nice x-displacement vector
  yfootprint<-xyfootprint[,,2] #and a nice y-displacement vector
  nlocal<-length(xfootprint)
  collocal<-rep(singlecolour,nlocal)
  #danger CONVENTION: axes labels are asumed the same because of the way H and J are created (means it wonl;t work for any old matrix)
  xlocal<-colnames(xyfootprint)[1]
  ylocal<-colnames(xyfootprint)[2]
  l1<-list(xyfootprint,xfootprint,yfootprint,nlocal,collocal,xlocal,ylocal)
  names(l1)<-c("xyfootprint","xfootprint","yfootprint","nlocal","collocal","xlocal","ylocal")

  return(l1)
}
plotoutbreakXYall<-function(Z,
                            timestep
                            ){

  #humanoids outbreak for x and y displacements for two matrices of humanoids
  #this is effectively a four dimensional vector space contained in Z
  #this function uses that to plot the X-Y displacement for humans and zombies in x-y x-coordinate for this epoch
  #by calling it repeatedly we can build an animation per plotourbreakXYTasGIf

  #Susceptibles
  m.Susceptible<-arrayStartWith(Z,prefix.Susceptible)
  #Zombies
  m.Zombie<-arrayStartWith(Z,prefix.Zombie)
  #Those killed on this iteration
  m.Death.Susceptibles.Killed<-arrayStartWith(Z,"RemovedXS")
  m.Death.Susceptibles.NaturalDeath<-arrayStartWith(Z,"RemovedXD")
  m.Death.Zombies.Killed<-arrayStartWith(Z,"RemovedXZ")
  #Mark those already dead
  m.Dead<-arrayStartWith(Z,"Removed[^X]")
  #and those just born
  m.Born<-arrayStartWith(Z,prefix.Newborn)
  
  #and those back from the dead
  m.Back<-arrayStartWith(Z,prefix.Resurrected)

  susplot<-plotsub(m.Susceptible,timestep = timestep,singlecolour = "blue")
  zomplot<-plotsub(m.Zombie,timestep = timestep,singlecolour = "red")


  #generic
  mlocal<-paste("x-y coordinates for ", susplot$nlocal, " ", prefix.Susceptible, "s,", sep="")
  mlocal<-paste(mlocal, " and ", zomplot$nlocal, " ", prefix.Zombie, "s,", sep="")
  mlocal<-paste(mlocal, " at t=t", (timestep-1), sep="")
  plot(x=c(), xlab="x-displacement", ylab="y-displacement", main=mlocal, cex.main = 1.0,  frame.plot = FALSE, xlim = config.plot$xlimlocal, ylim = config.plot$ylimlocal)
 
  
  if(matrixnn(m.Susceptible)){
    points(x=susplot$xfootprint,y=susplot$yfootprint,type="p",pch=19, col=susplot$collocal)
  }
  
  if(matrixnn(m.Zombie)){
   points(x=zomplot$xfootprint,y=zomplot$yfootprint,type="p", pch=15, col=zomplot$collocal)
  }
  
  # if(matrixnn(m.Death)){
  #   justplot<-plotsub(m.Death, timestep = timestep,singlecolour = "red")#plotsub(plotmatrix=m.Death,timestep = timestep,singlecolour = "red")
  #   points(x=justplot$xfootprint,y=justplot$yfootprint,type="p", col=justplot$collocal, pch=8, cex=2)
  # }
  
  if(matrixnn(m.Death.Susceptibles.Killed)){
    justplot<-plotsub(m.Death.Susceptibles.Killed, timestep = timestep,singlecolour = "blue")#plotsub(plotmatrix=m.Death,timestep = timestep,singlecolour = "red")
    points(x=justplot$xfootprint,y=justplot$yfootprint,type="p", col=justplot$collocal, pch=8, cex=2)
  }
  
  if(matrixnn(m.Death.Susceptibles.NaturalDeath)){
    justplot<-plotsub(m.Death.Susceptibles.NaturalDeath, timestep = timestep,singlecolour = "black")#plotsub(plotmatrix=m.Death,timestep = timestep,singlecolour = "red")
    points(x=justplot$xfootprint,y=justplot$yfootprint,type="p", col=justplot$collocal, pch=8, cex=2)
  }
  
  if(matrixnn(m.Death.Zombies.Killed)){
    justplot<-plotsub(m.Death.Zombies.Killed, timestep = timestep,singlecolour = "red")#plotsub(plotmatrix=m.Death,timestep = timestep,singlecolour = "red")
    points(x=justplot$xfootprint,y=justplot$yfootprint,type="p", col=justplot$collocal, pch=8, cex=2)
  }
   
  if(matrixnn(m.Dead)){
    deadplot<-plotsub(m.Dead,timestep = timestep,singlecolour = "black")#plotsub(plotmatrix=m.Dead,timestep = timestep,singlecolour = "black")
    points(x=deadplot$xfootprint,y=deadplot$yfootprint,type="p", col=deadplot$collocal, pch=13, cex=1)
  }
  
  if(matrixnn(m.Born)){
    bornplot<-plotsub(m.Born,timestep = timestep,singlecolour = "pink")#plotsub(plotmatrix=m.Born,timestep = timestep,singlecolour = "black")
    points(x=bornplot$xfootprint,y=bornplot$yfootprint,type="p", col=bornplot$collocal, pch=11, cex=1)
  }
  
  if(matrixnn(m.Back)){
    backplot<-plotsub(m.Back,timestep = timestep,singlecolour = "green")#plotsub(plotmatrix=m.Back,timestep = timestep,singlecolour = "black")
    points(x=backplot$xfootprint,y=backplot$yfootprint,type="p", col=backplot$collocal, pch=6, cex=1)
  }
  
  legend("bottom", c("Susceptibles", "Zombies", "Removed", "Birth", "Resurrection"), xpd = TRUE, horiz = TRUE, 
         inset = c(0, 0),bty = "n", col = c("blue","red", "black", "pink", "green"), cex = 0.5,pch=c(19,15,13,8,8))
}

plotoutbreakXYTasgif<-function(Z,
                               S.params,
                               Z.params,
                               outputfilename=config.plot$outputfilename,
                               unleashplot=TRUE){
  #builds a gif of the zombie vs susceptible action and outputs it to the specified filename (Downloads folder)
  if(unleashplot){
    library(animation)

    saveGIF({
      for (i in time.vtimesteps){
        j<-i+1
          
        #for each timestep we need to remove any humanoids that get too near, too born, or too tired (i.e. they die)
        Z<-reconfigurepopulation(Z=Z,
                                 S.params,
                                 Z.params,
                         timesteplocal=j)
        #look at just the suceptibles
        m.Susceptible<-arrayStartWith(Z,prefix.Susceptible)
        #look at just the zombies
        m.Zombie<-arrayStartWith(Z,prefix.Zombie)
        m.Susceptible.count<-length(rownames(m.Susceptible))
        m.Zombie.count<-length(rownames(m.Zombie))
        if(model.timestep.assessment(m.Susceptible.count,m.Zombie.count)){
          print(paste("Game Over at t=t",j))
          break #if there are no zombies or suscpetibles left then leave
        }
        plotoutbreakXYall(Z=Z,
                        timestep=j)
      }
     }, movie.name=outputfilename)
  }

    return(Z)
}
#####have some fun
time.ntimesteps<-99  #time steps
time.vtimesteps<-c(0:time.ntimesteps) #set up a zero-based vector, t, that will represent time, NB unit time in steps of one
# epidemiology parameters (note that these have been anglicised from Munz et al, e.g. Π is represented as P)
config.model<-configuration.SRZModel(P=0,d=0.0001,B=0.0095,G=0.0001,a=0.0001,rcritical=1,zombiewinratio=0.95)
config.Susceptible<-configuration.Humanoid(xinit=10,yinit=10,n=500,mu=0,sigma=1)# Susceptible configuration
config.Zombie<-configuration.Humanoid(xinit= 0,yinit=0,n=2,mu=0,sigma=1)# Zombie configuration
config.plot<-config.plot.global()#axes configuration for gif

#zombie playground :-)#
OutputDataSet<-model.run.params(S.params = config.Susceptible, Z.params = config.Zombie, mode="matrix", iterations = c(1:1))

'



    EXEC sys.sp_execute_external_script @language = N'R'
                                      , @script = N'OutputDataSet<-data.frame(installed.packages())'--@BigRs
                                      --, @input_data_1 = @SQLInput
                                      --, @input_data_1_name = @SQLInputName
                                      --, @params = N'@modelbin varbinary(max) OUTPUT'
                                      --, @modelbin = @model OUTPUT;