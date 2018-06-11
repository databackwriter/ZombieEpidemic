
#install.packages("RODBC")

# install.packages("RODBC", type = "source")
# install.packages("RODBC")
library(RODBC)
cn <- odbcConnect("LAYDSQL_Zombie")


tabledata <- sqlFetch(cn, "DimModel")
print(names(tabledata))


 ?odbc


df<-data.frame("hi","bye")
names(df)<-c("hello","goodbye")

de<-data.frame("hola","ciao")
names(de)<-c("hello","goodbye")

newdf <- rbind(df, de)