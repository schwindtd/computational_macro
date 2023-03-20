library(mvQuad)
library(mvtnorm)
setwd("/Users/danielschwindt/Documents/UMD/Fall 2022/ECON630/PS1/")

ghq3 = createNIGrid(dim=2, type="GHe", level=c(3,3))
ghq7 = createNIGrid(dim=2, type="GHe", level=c(7,7))
ghq15 = createNIGrid(dim=2, type="GHe", level=c(15,15))

data3 = as.data.frame(cbind(ghq3$weights, ghq3$nodes))
data7 = as.data.frame(cbind(ghq7$weights, ghq7$nodes))
data15 = as.data.frame(cbind(ghq15$weights, ghq15$nodes))

names(data3) = c("w_i", "x_i", "y_i")
names(data7) = c("w_i", "x_i", "y_i")
names(data15) = c("w_i", "x_i", "y_i")

write.csv(data3, file="ghq_2d_n3.csv", row.names=F)
write.csv(data7, file="ghq_2d_n7.csv", row.names=F)
write.csv(data15, file="ghq_2d_n15.csv", row.names=F)

## Test Code
C=matrix(c(4,-1,-1,1), ncol=2)
m=c(0,5)
myGrid = createNIGrid(dim=2, type="GHe", level=3)
rescale(myGrid, m=m, C=C, dec.type=0)
myFun = function(x){
  x[,1]*x[,2]
}
quadrature(myFun, myGrid)