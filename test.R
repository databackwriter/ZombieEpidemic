rm(list=ls()) #clear variables
cat("\014")  #clear console
#dev.off(dev.list()["RStudioGD"]) # clear phumanoids
library(memisc)
m <- matrix(1,2,2)
rownames(m) <- letters[1:2]
colnames(m) <- LETTERS[1:2]
m
dimrename(m,1,a="first",b="second")
dimrename(m,1,A="first",B="second")
dimrename(m,2,"A"="first",B="second")
m
ddd<a
rowrename(m,ddd="first",b="second")
colrename(m,"A"="first",B="second")



x <- c(a=1, b=2)
x

rename(x,a="A",b="B")
x
rename(iris,
           Sepal.Length="Sepal_Length",
           Sepal.Width ="Sepal_Width",
           Petal.Length="Petal_Length",
           Petal.Width ="Petal_Width"
)
rename(iris,
           al.="_"
           ,gsub=TRUE)




x <- cbind(x1 = 3, x2 = c(4:1, 2:5))
dimnames(x)[[1]] <- letters[1:8]
apply(x, 2, mean, trim = .2)










Nt=ntimesteps
Ni=nS
mu=muS
# sigma=sigmaS, xinit=xinitS, yinit=yinitS, timestepslocal = vtimesteps, humanoidtype=prefixS
Nt=ntimesteps
Ni=nZ
mu=muZ
# , sigma=sigmaZ, xinit=xinitZ, yinit=yinitZ, timestepslocal = vtimesteps, humanoidtype=prefixZ

H<-Hs
J<-Hz
criticaldistance = interfacedistance 
winratez = zombiewinrate
outputfilename="/Users/petermoore/Documents/GitHub/ZombieEpidemic/zombieanimation.gif"
sprime=3
zprime=2




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








