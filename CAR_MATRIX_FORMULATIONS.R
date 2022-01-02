Neighbor_Mat_Big <- matrix(0, nrow = 11, ncol = 11)
for(i in 1:ncol(Neighbor_Mat_Big)){
  Neighbor_Mat_Big[i,i+1] <- 1
  Neighbor_Mat_Big[i+1,i] <- 1
}

Neighbor_Mat_Big[4,5] <- 0
Neighbor_Mat_Big[5,4] <- 0
Neighbor_Mat_Big[4,11] <- 1
Neighbor_Mat_Big[5,11] <- 1
Neighbor_Mat_Big[11,10] <- 0
Neighbor_Mat_Big[10,11] <- 0
Neighbor_Mat_Big[11,4] <- 1
Neighbor_Mat_Big[11,5] <- 1

D_Big= diag(rowSums(Neighbor_Mat_Big))

Neighbor_Mat_Small <- matrix(0, nrow = 10, ncol = 10)
for(i in 1:ncol(Neighbor_Mat_Small)){
  Neighbor_Mat_Small[i,i+1] <- 1
  Neighbor_Mat_Small[i+1,i] <- 1
}
Neighbor_Mat_Small[4,5] <- .5
Neighbor_Mat_Small[5,4] <- .5

D_Small= diag(rowSums(Neighbor_Mat_Small))

var <- 0.5 # Variance
phi <- 0.2 # Correlation parameter

Var_Cov_Big <- solve(var*(D_Big - phi * Neighbor_Mat_Big))
Var_Cov_Small <- solve(var*(D_Small - phi * Neighbor_Mat_Small))