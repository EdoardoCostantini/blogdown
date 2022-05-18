# Project:   blogdown
# Objective: TODO
# Author:    Edoardo Costantini
# Created:   2022-05-18
# Modified:  2022-05-18

# Prepare environment ----------------------------------------------------------

library(psych)
library(PCAmixdata)

library("FactoMineR")
library("factoextra")

data(tea)

head(tea)

lapply(tea[, c("where", "how", "SPC")], nlevels)

sapply(tea, nlevels)

x <- tea[, c("where", "how", "Tea")]

CTD <- tab.disjonctif(x)
pk <- colMeans(CTD)

CTD[1, ] / pk
CTD[4, ] / pk
1/.64
1/.56666667

CTD_t <- t(apply(CTD, 1, function(r) {t(r)/pk} ))
colMeans(CTD_t)


N <- tab.disjonctif(x)
1/nrow(tea)

# Correspondance analysis based on contingency table ---------------------------
x <- tea[, c("SPC", "where")]
n <- nrow(x)
r <- nlevels(x[, 1]) 
C <- nlevels(x[, 2]) 
N <- table(x)

# Contingency table
N

# Indicator matrix
Z <- tab.disjonctif(x)
Z1 <- Z[, 1:r]
Z2 <- Z[, -c(1:r)]
N - t(Z1) %*% Z2

# Correspondance matrix
P <- 1/n * N

# From Greenacre1984
N

# Column and row sums
r_bold <- rowSums(N)
c_bold <- colSums(N)

drop(N %*% rep(1, ncol(N))) - r_bold
drop(t(N) %*% rep(1, nrow(N))) - c_bold

D_r <- diag(r_bold)
D_c <- diag(c_bold)

# Matrices of profiles

R <- solve(D_r) %*% P
C <- solve(D_c) %*% t(P)

# Centroids
r <- t(C) %*% c_bold
c <- t(R) %*% r_bold

# Generalized SVD of P - r_bold t(c_bold)

P - r_bold %*% t(c_bold)
A <- svd(P - r %*% t(c))$u
B <- svd(P - r %*% t(c))$v

t(A) %*% solve(D_r) %*% A
t(B) %*% solve(D_c) %*% B

# From Jolliffe p. 37
r_bold <- rep(1, r)
c_bold <- rep(1, C)
D_r <- diag(r_bold)
D_c <- diag(c_bold)
Omega <- solve(D_r)
Psi <- solve(D_c)
X <- P - r_bold %*% t(c_bold)

V <- svd(X)$u
M <- diag(svd(X)$d)
B <- svd(X)$v

round(t(V) %*% Omega %*% V, 3)
round(t(B) %*% Psi %*% B, 3)

X_til <- sqrt(Omega) %*% X %*% sqrt(Psi)

W <- svd(X_til)$u
K <- diag(svd(X_til)$d)
C <- svd(X_til)$v

solve(sqrt(Omega)) %*% W - V
solve(sqrt(Omega)) %*% C - B

W %*% K

# Row profiles
D_r
