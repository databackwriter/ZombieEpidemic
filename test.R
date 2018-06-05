rm(list=ls()) #clear variables
cat("\014")  #clear console
#dev.off(dev.list()["RStudioGD"]) # clear phumanoids

library("plot3D")
library("odeintr")

#nominal definitions
#x-"x" displacement in a two-dimensional plane where x and y are orthogonal
#y-"y" displacement in a two-dimensional plane where x and y are orthogonal
#t-time coordinate
#together x,y,t make a three-dimensional vector space in R3
