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
