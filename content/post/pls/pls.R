# Project:   blogdown
# Objective: Partial Least Square ideas
# Author:    Edoardo Costantini
# Created:   2022-06-13
# Modified:  2022-06-24

# Load packages ----------------------------------------------------------------

    # load pls package
    library(pls)
    library(plsdof)
    library(PCovR)

# Preapre data -----------------------------------------------------------------

        data(alexithymia)

    # 1. Select vairables for model
    dt <- cbind(yarn[[2]], yarn[[1]])
    dim(dt)

    # 2. Check scales
    # Means
    round(colMeans(dt), 3)

    # Variances
    round(apply(dt, 2, var), 3)

    # 3. Objects needed
    p <- ncol(dt) - 1
    n <- nrow(dt)
    dt_sd <- cbind(
        yarn[[2]],
        scale(yarn[[1]]) # center and scale to variance = 1
    )

    # Checks
    round(colMeans(dt_sd), 3)
    round(apply(dt_sd, 2, var), 3)

# Estimation -------------------------------------------------------------------

    # Helper function
    orthogonalize <- function(vec1, vec2) {
        v <- vec1
        u <- vec2

        newv <- v - drop(t(u) %*% v / (t(u) %*% u)) * u

        return(newv)
    }

    # Parms
    M         <- 10 # number of "components"
    X         <- lapply(1:M, matrix, nrow = n, ncol = p)
    X[[1]]    <- dt_sd[, -1]
    y         <- dt_sd[, 1]
    y_hat     <- cbind(
        mean(dt_sd[, 1]),
        matrix(rep(NA, n * (M - 1)), nrow = n)
    )
    z         <- matrix(NA, nrow = n, ncol = M)
    theta_hat <- rep(NA, M)

    # PLS Algorithm following HastieEtAl2017 p 81 (Algorithm 3.3)
    for (m in 2:M) {
        # 2a
        store_2a <- matrix(NA, nrow = n, ncol = p)
        for (j in 1:p) {
            rho_hat_mj <- t(X[[m - 1]][, j]) %*% y
            store_2a[, j] <- rho_hat_mj %*% X[[m - 1]][, j]
        }        
        z[, m] <- rowSums(store_2a)

        # 2b
        theta_hat[m] <- drop(t(z[, m]) %*% y / t(z[, m]) %*% z[, m])

        # 2c
        y_hat[, m] <- y_hat[, m - 1] + theta_hat[m] * z[, m]

        # 2d orthogonalize all columns
        for (j in 1:p) {
            X[[m]][, j] <- orthogonalize(X[[m-1]][, j], z[, m])
        }
    }

    # PLS Algorithm following MevikWehrens2007 (PLS package)
    E <- X
    F <- vector("list", M)
        F[[1]] <- y
    T <- matrix(nrow = n, ncol = M)
    W <- P <- matrix(nrow = p, ncol = M)
    P <- matrix(nrow = p, ncol = M)
    Q <- matrix(nrow = 1, ncol = M)

    for (m in 1:M) {
        # SVD of cross product
        S <- t(E[[m]]) %*% F[[m]]
        svdS <- svd(S)
        w <- drop(svdS$u)
        q <- drop(svdS$v)

        # Scores
        t <- E[[m]] %*% w
        t <- t / drop(sqrt(t(t) %*% t))
        u <- y * q

        # Weights
        p <- t(E[[m]]) %*% t
        q <- t(F[[m]]) %*% t
        
        # Deflation
        E[[m + 1]] <- E[[m]] - t %*% t(p)
        F[[m + 1]] <- F[[m]] - drop(t %*% t(q))
    }

    # Fit PCR model w/ pls package
    pls_fit_pls <- plsr(
        dt[, 1] ~ dt[, -1],
        ncomp = M,
        scale = FALSE,
        center = FALSE,
        method = "oscorespls",
        validation = "none"
    )

    # Fit PCR model w/ plsdof package
    pls_fit_plsdof <- pls.model(X[[1]], y)

    # Copmare y_hats
    m <- 2
    data.frame(
        pls = round(as.data.frame(fitted(pls_fit_pls)), 3)[, m],
        plsdof = round(pls_fit_plsdof$Yhat, 3)[, m],
        man = round(y_hat, 3)[, m]
    )

# Types of DVs -----------------------------------------------------------------

    data(oliveoil)
    sens.pcr <- pcr(sensory ~ chemical, ncomp = 4, scale = TRUE, data = oliveoil)
    sens.pls <- plsr(sensory ~ chemical, ncomp = 4, scale = TRUE, data = oliveoil)

    oliveoil$sensory

    cppls.fito

# PLS degrees of freedom -------------------------------------------------------

    library(plsdof)

    data(Boston)
    X <- as.matrix(Boston[, -14])
    y <- as.vector(Boston[, 14])

    dim()

    my.pls1 <- pls.model(X, y)
    
    my.pls1

    