### Object: Performing Multiple Correspondance Analysis in R 
### Source: Practical Guide to Principal Component Methods 
###         (Chapter 5 Code)

# Load Packages
  library("FactoMineR") # for analysis
  library("factoextra") # for visualizarion
  library("httpg")

# Data Prep ---------------------------------------------------------------

  data("poison")
  head(poison[, 1:7], 3) # survey style data
  
# Subset active individuals and variables
  poison.active <- poison[1:55, 5:15] 
  head(poison.active[, 1:6], 3)
  
# Summaries
  str(poison.active)
  for (i in 1:4) {
    plot(poison.active[,i], main=colnames(poison.active)[i],
         ylab = "Count", col="steelblue", las = 2)
  }
  
# The analysis ------------------------------------------------------------
  
  res.mca <- FactoMineR::MCA(X = poison.active,
                             ncp = 5,
                             graph = TRUE)

# Visualization -----------------------------------------------------------

  # Eigenvalues / Variances
  eig.val <- get_eigenvalue(res.mca)
  eig.val # proportion of variances retained by dimensions
  
  # Percentages of Inertia explained by MCA
  fviz_screeplot(res.mca, addlabels = TRUE, ylim = c(0,45))
  
  # Biplot
  fviz_mca_biplot(res.mca,
                  repel = TRUE, # avoid text overlapping
                  ggtheme = theme_minimal())
    # Rows (individuals) are represented by blue points;
    # Columns (variable categories) by red triangles.
  
  # Graph of variables
  var <- get_mca_var(res.mca) 
  var
  # Coordinates 
  head(var$coord)
  # Cos2: quality on the factore map 
  head(var$cos2)
  # Contributions to the principal components
  head(var$contrib)
  
  # Graph of individuals
  ind <- get_mca_ind(res.mca) # extract the results for individuals
  ind
  # Coordinates of column points
  head(ind$coord)
  # Quality of representation
  head(ind$cos2)
  # Contributions
  head(ind$contrib)
  # BIplot for individuals only (no vars)
  fviz_mca_ind(res.mca, col.ind = "cos2",
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE, ggtheme = theme_minimal())
  fviz_mca_ind(res.mca, 
               label = "none",
               habillage = "Vomiting", # color by groups defined by variable
               palette = c("#00AFBB", "#E7B800"),
               addEllipses = TRUE, ellipse.type = "confidence",
               repel = TRUE, ggtheme = theme_minimal())
  # More than 1 grouping variable
  fviz_ellipses(res.mca, c("Vomiting", "Fever"),
                geom = "point")
  # Bar plot for Cos2 of individuals 
  fviz_cos2(res.mca, choice = "ind", axes = 1:2, top = 20)
  # Contribution of individuals to the dimensions 
  fviz_contrib(res.mca, choice = "ind", axes = 1:2, top = 20)

# Supplementary elements --------------------------------------------------

  res.mca <- MCA(poison, 
                 ind.sup = 53:55)
  ind <- get_mca_ind(res.mca)
  ind$coord
  res.mca$ind.sup # Supplementary individuals 
  
# Doing it by hand --------------------------------------------------------
  
  # Data
  data(wine)
  MCA.dt <- wine[, sapply(wine, is.factor)]
  sapply(MCA.dt, nlevels)
  
  # Goal: Find the coordinates of 
  npcs <- 4
  res.mca <- FactoMineR::MCA(X = MCA.dt, ncp = npcs, graph = FALSE)
  res.mca$ind$coord # Coordinates of the individuals on the dimensions
  res.mca$svd$V
  
# Following AudigierEtAl2016_MICAT ----------------------------------------
  
  I <- nrow(MCA.dt)     # numebr of individuals
  K <- ncol(MCA.dt)     # number of categorical predictors
  qk <- sapply(MCA.dt,  # number of levels per categorical variable 
               nlevels) 
  J <- sum(qk)          # total number of categories
  
  # Disjunctive table
  Z <- tab.disjonctif(MCA.dt)
  # Notice the relationship between the disjunctive table and a 
  # contingency table 
  N <- t(Z[, 1:3]) %*% Z[, 4:7]
  N - table(MCA.dt)
  
  # Distance Metric Matrix
  pxkqk <- colSums(Z)/I # props ids taking category value on variable
  D_Sigma <- diag(pxkqk)
  
  # Weight Matrix
  W_mat <- diag(rep(1, I))/I
  
  # M matrix (center matrix)
  M <- matrix((rep(colMeans(Z), nrow(Z))), 
              nrow = nrow(Z), 
              byrow = TRUE)
  
  # SVD of triplet (Z-M, D_Sigma, W_mat)
  SVD.trip <- svd.triplet(X = Z - M,
                          row.w = diag(W_mat),
                          col.w = diag(1/K*solve(D_Sigma)),
                          ncp = npcs)
  
  # Manual SVD triplet
  SVD.man <- svd(sqrt(W_mat) %*% (Z - M) %*% sqrt(1/K*solve(D_Sigma)))
  
  # Convert back to correct scales (according to Chaven 2017 p. 3)
  V.man <- (solve(sqrt(1/K*solve(D_Sigma))) %*% SVD.man$v)[, 1:npcs]
  U.man <- (solve(sqrt(W_mat)) %*% SVD.man$u)[, 1:npcs]
  L.man <- SVD.man$d[1:npcs]

  # Compare SVD triplet and manual SVD of weighted matrix
  round(abs(SVD.trip$V) - abs(V.man), 5)
  round(abs(SVD.trip$U) - abs(U.man), 5)
  round(abs(SVD.trip$vs[1:npcs]) - abs(L.man), 5)
  
  # Reconstruction Formula
  d_hat <- SVD.trip$vs[1:npcs] # matrix of the singular values 
                            # (Squared would be eigenvalues of Z)
  u_hat <- SVD.trip$U # Left singular vectors matrix
  v_hat <- SVD.trip$V # Right singular vectors matrix
  
  z_hat <- u_hat %*% diag(d_hat) %*% t(v_hat) + M
  colSums(z_hat)
  colSums(Z)
  
  # Compare SVD matrices
  # Matrix of singular values
  res.mca$svd$vs[1:npcs] -
    d_hat[1:npcs]

  res.mca$svd$vs[1:npcs] -
    L.man
  
  # Left Singular Vectors Matrix
  round(
    res.mca$svd$U - u_hat,
    3
  )

  round(
    res.mca$svd$U - U.man,
    3
  )

  # Right Singular Vectors Matrix
  round(
    res.mca$svd$V -
      v_hat, 
    3)
  # Correlation between columns
  # And look into the PCAmixdata package
  round(cor(v_hat, res.mca$svd$V), 1)
  
  # Coordinates on Dimensions are recovered
  round(
    res.mca$ind$coord -
      u_hat %*% diag(d_hat),
    3
  )

# Following JosseHusson2016 -----------------------------------------------

  I <- nrow(MCA.dt)
  J <- ncol(MCA.dt)
  X <- tab.disjonctif(MCA.dt)
  rowMarg <- rowSums(X) # = J
  colMarg <- colSums(X) # = number of ids in a category
  D_Sigma <- diag(colMarg)
  D <- 1/I * diag(rep(1, I)) # rowMasses
  SVD.trip <- svd.triplet(X = I %*% X %*% solve(D_Sigma),
                          row.w = diag( D ),
                          col.w = diag( 1/(I*J)*D_Sigma ),
                          ncp=2
                          )
  svd(I * X %*% solve(D_Sigma))
  
  SVD.trip$vs
  round(SVD.trip$vs[2:3] - res.mca$svd$vs[1:2], 3)
  
  SVD.trip$U
  res.mca$svd$U
  round(SVD.trip$U[, 2:3] - res.mca$svd$U, 3)
  res.mca$svd$U
  
  SVD.trip$V
  round(SVD.trip$V[, 2:3] - res.mca$svd$V, 3)
  res.mca$svd$V
  
  

  