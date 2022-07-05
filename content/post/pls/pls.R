# Project:   blogdown
# Objective: Partial Least Square ideas
# Author:    Edoardo Costantini
# Created:   2022-06-13
# Modified:  2022-07-05

# Prepare environment ----------------------------------------------------------

    # Load packages
    library(pls)
    library(plsdof)
    library(PCovR)

    # Load personal funcitons
    source("./content/post/pls/dfSDI.R")
    source("./content/post/pls/pls.manual.R")

# Preapre data -----------------------------------------------------------------

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

# Prediction -------------------------------------------------------------------

    n <- 50 # number of observations
    p <- 15 # number of variables
    X <- matrix(rnorm(n * p), ncol = p)
    y <- rnorm(n)
    M <- 10 # number of "components"

    ntest <- 200 #
    Xtest <- matrix(rnorm(ntest * p), ncol = p) # test data
    ytest <- rnorm(ntest) # test data

    # Fit alternative PLSs
    out_pls <- plsr(
        y ~ X,
        ncomp = M,
        scale = TRUE,
        center = TRUE,
        method = "oscorespls",
        validation = "none"
    )
    out_plsdof <- pls.model(X, y, compute.DoF = TRUE, Xtest = Xtest, ytest = NULL)

    # Manual fit
    outpls_man <- pls.manual(ivs = X, dv = y, m = M)

    # Obtain predictions on new data
    round(
        cbind(
            PLS = predict(out_pls, newdata = Xtest)[, , M],
            PLSdof = out_plsdof$prediction[, M]
        ), 5
    )

# PLS degrees of freedom -------------------------------------------------------

    library(plsdof)
    set.seed(1234)

    # Generate data data
    n <- 100 # number of observations
    p <- 15 # number of variables
    m <- 15
    X <- matrix(rnorm(n * p), ncol = p)
    y <- rnorm(n)

    # Fit model with package
    outpls <- pls.model(X, y, compute.DoF = TRUE)
    outpls$DoF
    outpls.internal <- linear.pls.fit(X, y, m, DoF.max = min(n - 1, p + 1))

    # Fit model with person PLS function
    outpls_man <- pls.manual(ivs = X, dv = y, m = m)

    # Y hats
    round(outpls_man$Yhat - outpls$Yhat, 5)

    # T scores
    j <- 1
    cbind(
        PLSTT = outpls.internal$TT[, j],
        manualTs = outpls_man$Ts[, j],
        manualTsn = outpls_man$Tsn[, j]
    )

    # Degrees of freedom
    DoF_manual <- dofPLS(
        X,
        y,
        TT = outpls_man$Tsn,
        Yhat = outpls_man$Yhat[, 2:(m + 1)],
        DoF.max = m + 1
    )

    cbind(
        PLS = outpls$DoF,
        PLS.manual = DoF_manual,
        diff = round(outpls$DoF - DoF_manual, 5)
    )
