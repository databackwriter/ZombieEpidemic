#keep R tidy
rm(list=ls()) #clear variables
cat("\014")  #clear console
#dev.off(dev.list()["RStudioGD"]) # clear phumanoids

# library("plot3D")
# library("odeintr")

#nominal definitions
#x-"x" displacement in a two-dimensional plane where x and y are orthogonal
#y-"y" displacement in a two-dimensional plane where x and y are orthogonal
#t-time coordinate
#together x,y,t make a three-dimensional vector space in R3

#global parameters
prefixS<-"Susceptible"
prefixZ<-"Zombie"
primeS<-3
primeZ<-2

#time parameters
ntimesteps<-3 #time steps
vtimesteps<-c(0:ntimesteps) #set up a zero-based vector, t, that will represent time, NB unit time in steps of one

# epidemiology parameters (note that these have been anglicised from Munz et al, e.g. Π is represented as P)
P<-0      # birth rate
d<-0.0001  # natural death percent (per unit time)
B<-0.0095  # transmission percent  (per unit time)
G<-0.0001  # resurrect percent (per unit time)
a<-0.0001  # destroy percent  (per unit time)

# Zombie configuration
xinitZ<-0 #set a starting co-ordinate for Zombies on the x-dimension
yinitZ<-0 #set a starting co-ordinate for Zombies on the y-dimension
nZ<-5 #initial number of Zombies
muZ<-0 #mean footstep length for Zombies IMPORTANT: because we have unit time this is average zombie speed
sigmaZ<-1 #standard deviation of footstep lengths for Zombies

# Susceptible configuration
xinitS<-10#set a starting co-ordinate for Susceptibles on the x-dimension
yinitS<-10#set a starting co-ordinate for Susceptibles on the y-dimension
nS<-6 #initial number of Susceptibles
muS<-0 #mean footstep length for Susceptibles
sigmaS<-1 #standard deviation of footstep lengths for Susceptibles (observation: for Susceptibles, sigmaS is akin to a drunkenness factor!)

# #Removed configuration
# nRemoved<-0 #initial number of Removeds ##ASSUMPTION: start with zero in all cases
# matRemoved <- matrix((NA),nrow=0,ncol=0) #start with an empty matrix to collect Removals
# 
# nResurrectabled<-0 #initial number of Removeds ##ASSUMPTION: start with zero in all cases
# matResurrectable <- rwalk1dhumanoids(1,1,1,1,1) #start with an empty matrix to collect Resurrectable


#killing parameters
interfacedistance<-2 #distance at which zombies infect susceptibles or susceptiles kill zombies
zombiewinrate<-0.95

#random walk functions
rwalk1dhumanoids<-function(Nt, Ni, mu, sigma, xinit, timestepslocal, humanoidtype="Humanoid"){
  #builds a one-d random walk along  nominal time axis
  #where mu=0 this is akin to Brownian motion
  #inspired by ideas at http://www.phytools.org/eqg/Exercise_4.1/
  #Nt-number of timesteps, Ni-size of initial Humaoid population,
  #mu,sigma-mean and standard deviation for footstep length
  #xinit-initial x-position
  Nlarge<-Nt*Ni # number of time steps multiplied by initial number of humanoids
  Xt<-rep(xinit,Ni) #set of start points (one for each humanoid)
  footsteps<-rnorm(n=Nlarge,mean=mu,sd=sigma) #lets now get a super long set of footsteps
  X<-matrix(data=footsteps,nrow=Ni,ncol=Nt) #now get a matrix where we get Nt columns for each of our Ni humanoids
  X<-cbind(Xt,X) #now place our initial co-ordinate against x (literally x vs t)
  X<-apply(X=X,MARGIN=1,FUN=cumsum) 
  X<-t(X) #SUBTLETY: transpose
  dimnames(X)=list(paste(humanoidtype,c(1:Ni)),paste0("t", timestepslocal))
  return(X)
}
#get your two arrays, break them into x and y co-ords (Ax,Bx,Ay,By), adapt lines 83 to 90
humanx<-rwalk1dhumanoids(Nt=7, Ni=3, mu = 1, sigma= 1, xinit = 10, timestepslocal=c(0:7),humanoidtype=prefixS)
humany<-rwalk1dhumanoids(Nt=7, Ni=3, mu = 1, sigma= 1, xinit = 10, timestepslocal=c(0:7),humanoidtype=prefixS)

zombiex<-rwalk1dhumanoids(Nt=7, Ni=6, mu = 0, sigma= 1, xinit = 1, timestepslocal=c(0:7),humanoidtype=prefixZ)
zombiey<-rwalk1dhumanoids(Nt=7, Ni=6, mu = 0, sigma= 1, xinit = 1, timestepslocal=c(0:7),humanoidtype=prefixZ)

xcoords<-rbind(humanx,zombiex)
ycoords<-rbind(humany,zombiey)
allcoords<-cbind(xcoords,ycoords)
Zdim<-c(3+6,7+1,2)#humanoid,time,x-y
W<-array(data= allcoords,
         dim = Zdim,
         dimnames=list(rownames(allcoords),paste0("t", c(0:7)),c("x-displacement","y-displacement"))
)

rwalk2dhumanoids<-function(Nt, Ni, mu, sigma, xinit, yinit, timestepslocal, humanoidtype="Humanoid"){
  
  #example usage: Hs<-rwalk2dhumanoids(Nt=ntimesteps, Ni=nS, mu=muS, sigma=sigmaS, xinit=xinitS, yinit=yinitS, timestepslocal = vtimesteps, useunif=FALSE, humanoidtype=prefixS)
  X<-rwalk1dhumanoids(Nt=Nt, Ni=Ni, mu=mu, sigma=sigma, xinit=xinit,timestepslocal=timestepslocal,humanoidtype=humanoidtype) # get Ni rows of humanoids X co-ordinates
  Y<-rwalk1dhumanoids(Nt=Nt, Ni=Ni, mu=mu, sigma=sigma, xinit=yinit,timestepslocal=timestepslocal,humanoidtype=humanoidtype) # get Ni rows of humanoids Y co-ordinates
  Z<-cbind(X,Y)
  Zdim<-c(Ni,Nt+1,2)#humanoid,time,x-y
  W<-array(data= Z,
           dim = Zdim,
           dimnames=list(paste(humanoidtype,c(1:Ni)),paste0("t", timestepslocal),c("x-displacement","y-displacement"))
           )
  return(W)
}


#linear algebra functions
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

deleterowzero<-function(M, removeNA = TRUE){
  #takes a matrix and returns the values that match the zero for that row and leaves the rest as NA
  Mrs<-rowSums(M, na.rm=TRUE)
  L<-subset(M,Mrs>0)
  return(L)
}
deletecolzero<-function(M, removeNA = TRUE){
  #takes a matrix and returns the values that match the zero for that column and leaves the rest as NA
  return(t(deleterowzero(t(M),removeNA))) #note that we transpose twice
}
deletezero<-function(M, removeNA = TRUE){
  #takes a matrix and returns the values that match the zero for that row and column and leaves the rest as NA
  return(deletecolzero(deleterowzero(M, removeNA), removeNA))
}

#death functions
deathmatrix<-function(H,J,timestep,interfacedistance, zombiewinrate, sprime, zprime){
  # a matrix of zombie winners (2's), zombie losers (3's), and live to fight another day-ers (NA's)
  # it is important that these numbers are primes
  xyfootprintH<-as.matrix(H[,timestep,]) 
  xyfootprintJ<-as.matrix(J[,timestep,]) 
  
  #important step, get distances between humanoids
  pdistHJ<-euclideandistancebetweentwomatrices(xyfootprintH, xyfootprintJ)
  
  L<-isolateminima(pdistHJ) #isolate minima (ASSUMPTION only nearest neighbours do anything)
  x<-ifelse(runif(length(L))>zombiewinrate, sprime, zprime) #CONVENTION 2 means zombie wins, 3 means human wins (zombie doesn't), NA means do nothing
  K<-ifelse(L<interfacedistance,x,NA) #create a matrix of winners, losers, and live to fight another day-ers
  # K<-deletezero(K)
  return(K)
}
humanoiddeathlist<-function(dm, susceptiblename, zombiename, humanoidtype, sprime, zprime){
  # internal function purpose: take an n x m matrix of humanoids to die
  # return a list of what's about to die
  deathprime<-0
  if (humanoidtype==zombiename){
    deathprime<-sprime #if we are looking at a matrix of zombies then the sprime kills them
  } else if (humanoidtype==susceptiblename){
    deathprime<-zprime # and vice versa
  } else{
    deathprime<-NA
  }
  
  dmatrowsums<-rowSums(dm, na.rm=TRUE) 
  rowdeathlist<-subset(dmatrowsums,dmatrowsums==deathprime) 
  rowdeathlist<-names(rowdeathlist)
  return(rowdeathlist)
}
humanoidtype<-function(hmat,susceptiblename, zombiename){
  humanoidnamesH<-rownames(hmat)
  humanoidtypeH<-sub(" .*","",humanoidnamesH[1])
  if(humanoidtypeH==susceptiblename){ #kill the susceptibles
    Htype<-susceptiblename
  } else if (humanoidtypeH==zombiename){
    Htype<-zombiename
  }
}

#plot functions
plotoutbreakXT<-function(H,verticalaxis=1){
  #plot humanoids outbreak versus time
  #verticalaxiss:
  #1 means plot x-displacement versus time
  #2 means plot y-displacement versus time
  nlocal<-length(H[,1,1]) #time and displacement fixed, so variable length is number of humanoids
  tlen<-length(H[1,,1])-1 #humanoid held constant, axis is constant (time varies)
  tlocal<-c(0:tlen)
  H1<-H[1,,verticalaxis]
  H2n<-H[2:nlocal,,verticalaxis]
  ylocal<-colnames(H[1,,])[verticalaxis] 
  ylimlocal<-c((-1*tlen),tlen)
  humanoidtype<-sub(" .*","",rownames(H[,1,])[verticalaxis])
  mlocal<-paste(ylocal, "versus time for a", humanoidtype)
  plot(x=tlocal,y=H1,type="l", xlab="time", ylab=ylocal, main=mlocal, frame.plot = FALSE, ylim = ylimlocal)
  apply(X=H2n,MARGIN=1,FUN=lines,x=tlocal)
  invisible(H)
}
plotoutbreakXY<-function(H,timestep){
  #phumanoids outbreak for x and y displacements for a specified time

  xyfootprint<-as.matrix(H[,timestep,]) #get ourselves a nice two-D matrix
  
  xfootprint<-xyfootprint[,1] #a nice x-displacement vector
  yfootprint<-xyfootprint[,2] #and a nice y-displacement vector
  
  nlocal<-length(xfootprint) # take length of xfootprint as our num ber of Zombies (could just as easily take y becvause time is fixed and we have an x-y for all)

  xlocal<-colnames(xyfootprint)[1] 
  ylocal<-colnames(xyfootprint)[2]
  
  xlimlocal<-c(floor(min(H[,,1])),ceiling(max(H[,,1]))) #note here that we are ignoring time altogether
  ylimlocal<-c(floor(min(H[,,2])),ceiling(max(H[,,2])))
  
  humanoidnames<-rownames(xyfootprint)
  humanoidtype<-sub(" .*","",humanoidnames[1])
  
  collocal<-rainbow(nlocal,end=1/6)
  
  mlocal<-paste("x-y coordinates for ", nlocal, " ", humanoidtype, "s, at t=t", (timestep-1), sep="")
  p<-plot(x=xfootprint,y=yfootprint,type="p", xlab=xlocal, ylab=ylocal, main=mlocal, frame.plot = FALSE,col=collocal, xlim = xlimlocal, ylim = ylimlocal)
  # text(x=xfootprint,y=yfootprint,labels = humanoidnames, pos = 4)
  # invisible(H)
}
plotoutbreakXYall<-function(H,
                            J,
                            timestep, 
                            deathmatrixlocal,
                            susceptiblename = "Susceptible", 
                            zombiename = "Zombie",
                            sprime=3,
                            zprime=2
                            ){

  #humanoids outbreak for x and y displacements for two matrices of humanoids
  #this is effectively a four dimensional vector space
  #first know (internally) what type of humanoids we have
  humanoidtypeH<-humanoidtype(hmat=H,susceptiblename=susceptiblename, zombiename=zombiename)
  humanoidtypeJ<-humanoidtype(hmat=J,susceptiblename=susceptiblename, zombiename=zombiename)

  #whatever is in the H variable will be in the rows of the death matrix 

  
  #before we start this iteration, we need to see if anything got killed on the last iteration
  deathsH<-humanoiddeathlist(dm = deathmatrixlocal,
        susceptiblename=susceptiblename,
        zombiename=zombiename,
        humanoidtype=humanoidtypeH,
        sprime=sprime,
        zprime=zprime)
  deathsJ<-humanoiddeathlist(dm = t(deathmatrixlocal), #SUBTLETY: transpose the matrix
                             susceptiblename=susceptiblename,
                             zombiename=zombiename,
                             humanoidtype=humanoidtypeJ,
                             sprime=sprime,
                             zprime=zprime)
  

  
  # s<-sum(deathmatrixlocal,na.rm =TRUE)
  # if(s > 0){
  #   removaltimesteps<-c(timestep:ntimesteps+1)
  #   retentiontimesteps<-c(0:timestep-1)
  #   print(timestep)
  #   print(deathmatrixlocal)
  #   if(length(deathsH)>0){
  #     removedH<-(H[deathsH,removaltimesteps,])
  #     print("removed H")
  #     print(removedH)
  #   }
  #   if(length(deathsJ)>0){
  #     removedJ<-(J[deathsJ,removaltimesteps,])
  #     print("removed J")
  #     print(removedJ)
  #   }
  # }

  #now start to visualise
  #get the x-y details for the matrix H
  xyfootprintH<-as.matrix(H[,timestep,]) #get ourselves a nice two-D matrix
  xfootprintH<-xyfootprintH[,1] #a nice x-displacement vector
  yfootprintH<-xyfootprintH[,2] #and a nice y-displacement vector
  nlocalH<-length(xfootprintH) # take length of xfootprint as our num ber of Zombies (could just as easily take y becvause time is fixed and we have an x-y for all)
  # xlocalH<-colnames(xyfootprint)[1] 
  # ylocalH<-colnames(xyfootprint)[2]
  xspreadH<-H[,,1]
  yspreadH<-H[,,2]
  collocalH<-rainbow(nlocalH,start=0,end=1/6)


  
  #get the x-y details for the matrix H
  xyfootprintJ<-as.matrix(J[,timestep,]) #get ourselves a nice two-D matrix
  xfootprintJ<-xyfootprintJ[,1] #a nice x-displacement vector
  yfootprintJ<-xyfootprintJ[,2] #and a nice y-displacement vector
  nlocalJ<-length(xfootprintJ) # take length of xfootprint as our num ber of Zombies (could just as easily take y becvause time is fixed and we have an x-y for all)
  # xlocalJ<-colnames(xyfootprint)[1] 
  # ylocalJ<-colnames(xyfootprint)[2]
  xspreadJ<-J[,,1]
  yspreadJ<-J[,,2]
  collocalJ<-rainbow(nlocalJ,start=5/6, end = 1)#NB this is the only strongly-typed difference and we would make more gerneric were there more types
 
   if(humanoidtypeJ==susceptiblename){ #kill the susceptibles
    print("||")
  }
  
  
  #define axes limits on the whole data set
  xspreadall<-union(xspreadH,xspreadJ)
  yspreadall<-union(yspreadH,yspreadJ)
  xlimlocal<-c(floor(min(xspreadall)),ceiling(max(xspreadall))) #note here that we are ignoring time altogether
  ylimlocal<-c(floor(min(yspreadall)),ceiling(max(yspreadall)))

  #danger CONVENTION: axes labels are asumed the same because of the way H and J are created (means it wonl;t work for any old matrix)
  xlocal<-colnames(xyfootprintH)[1] 
  ylocal<-colnames(xyfootprintH)[2]
  
  #generic
  mlocal<-paste("x-y coordinates for ", nlocalH, " ", humanoidtypeH, "s,", sep="")
  mlocal<-paste(mlocal, " and ",  nlocalJ, " ", humanoidtypeJ, "s,", sep="")
  mlocal<-paste(mlocal, " at t=t", (timestep-1), sep="")
  plot(x=xfootprintH,y=yfootprintH,type="p", xlab=xlocal, ylab=ylocal, main=mlocal, cex.main = 1.0, frame.plot = FALSE,col=collocalH, xlim = xlimlocal, ylim = ylimlocal)
  points(x=xfootprintJ,y=yfootprintJ,type="p", col=collocalJ)
  
  
  # text(x=xfootprint,y=yfootprint,labels = humanoidnames, pos = 4)
  # invisible(H)
}

plotoutbreakXYTasgif<-function(H,
                               J,
                               criticaldistance = 0, 
                               winratez = 0.95, 
                               outputfilename="/Users/petermoore/Documents/GitHub/ZombieEpidemic/zombieanimation.gif",
                               sprime=3,
                               zprime=2){
  #builds a gif of the zombie vs susceptible action and outputs it to the specified filename (Downloads folder)
  library(animation)
  saveGIF({
    for (i in vtimesteps){
      j<-i+1
      dm<-deathmatrix(H=H,
                      J=J,
                      timestep=j,
                      interfacedistance=interfacedistance, 
                      zombiewinrate=zombiewinrate,
                      sprime=sprime,
                      zprime=zprime)
      plotoutbreakXYall(H,J,j,dm)
    }
  }, movie.name=outputfilename)
}

#experimental area
#get a matrix of zombies
Hs<-rwalk2dhumanoids(Nt=ntimesteps, Ni=nS, mu=muS, sigma=sigmaS, xinit=xinitS, yinit=yinitS, timestepslocal = vtimesteps, humanoidtype=prefixS)
Hz<-rwalk2dhumanoids(Nt=ntimesteps, Ni=nZ, mu=muZ, sigma=sigmaZ, xinit=xinitZ, yinit=yinitZ, timestepslocal = vtimesteps, humanoidtype=prefixZ)

plotoutbreakXYTasgif(H=Hs,J=Hz,criticaldistance = interfacedistance, winratez = zombiewinrate)







H<-Hs
J<-Hz
criticaldistance = interfacedistance 
winratez = zombiewinrate
outputfilename="/Users/petermoore/Documents/GitHub/ZombieEpidemic/zombieanimation.gif"
sprime=3
zprime=2


H[1:2,1:2,1:2]


H[1,1:2,1:2] #ignores humanoids (gives only Susceptible 1)
H[1:2,1,1:2] #ignores time (gives only t0)
H[1:2,1:2,1] #ignores y-axis (gives only x-displacement)

myrownames<-rownames(J)
oldnames<-c("Susceptible 2", "Zombie 3")
newnames<-c("Rod", "Jane")
for (i in 1:length(oldnames)){ #blimey, there has to be a more efficient way than this
  myrownames<-replace(myrownames, myrownames == oldnames[i], newnames[i])
}
rownames(J)<-myrownames
J


nm<-rownames(J)
nm
startwithZombie<- nm %in% grep("^Z", nm, value = TRUE)
nm<-subset(nm,startwithZombie)
nm
J[nm,,]





















