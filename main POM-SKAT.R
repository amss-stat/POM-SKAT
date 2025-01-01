#### Input
## genotype matrix G, dim: n*m
## phenotype vector Y, dim: n
## covariate matrix X, dim: n*p

#### package
library(MASS)
library(PearsonDS)

#### POM-SKAT

n <- nrow(G)
M <- ncol(G)
K <- if (is.matrix(X)) ncol(X) else 1
Yfac <- as.factor(Y)
J <- nlevels(Yfac)
datafac <- data.frame(Yfac,X,G)
colnames(datafac) <- paste("X", 1:ncol(datafac), sep = "")
datanum <- data.frame(Y,X,G)
colnames(datanum) <- paste("X", 1:ncol(datanum), sep = "")

### null
formula <- as.formula(paste("X1 ~", paste("X", 2:(K+1), sep = "", collapse = " + ")))
remod <- polr(formula,data = datafac,Hess = TRUE)
hats <- summary(remod)$coefficients[,"Value"]
hatb <- hats[1:K]
hatajs <- sapply(1:(J-1), function(j) hats[K + j])

### estimator
mu0 <- rowSums(sapply(hatajs, function(hataj) {
  exp(X %*% hatb - hataj) / (1 + exp(X %*% hatb - hataj))
})) 
v0 <- diag(sqrt(as.vector(
  mu0 - mu0^2 + rowSums(sapply(2:(J-1), function(j) {
    2 * (j-1) * exp(X %*% hatb - hatajs[j]) / (1 + exp(X %*% hatb - hatajs[j]))
  }))
)))

Vb <- diag((dbeta(colMeans(G)/2,1/2,1/2))^2)
GN <- G
for (i in 1:m) {
  formula <- as.formula(paste( paste("X", i+1+K, sep = ""),paste("~"), paste("X", 2:(ncol(X)+1), sep = "", collapse = " + ")))
  mG <- lm(formula, data = datanum)
  GN[,i] <- predict(mG)
}
Gn <- G-GN

### test
U <- t(Y-mu0)%*%Gn%*%Vb%*%t(Gn)%*%(Y-mu0)/2
eige <- v0%*%Gn%*%Vb%*%t(Gn)%*%v0/2  
lambda <- eigen(eige)$values[which(eigen(eige)$values>0.0001)]
  
EU <- sum(lambda)
VU <- 2*sum(lambda^2)
TU <- 8*sum(lambda^3)/(VU^(3/2))
b = (2/TU)^2
a = TU/2
c = -2/TU
t = (U-EU)/sqrt(VU)
p = 1-ppearsonIII(t,b,c,a)
  
