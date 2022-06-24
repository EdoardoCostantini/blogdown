# Project:   blogdown
# Objective: Partial Least Square ideas
# Author:    Edoardo Costantini
# Created:   2022-06-13
# Modified:  2022-06-24

# Load packages ----------------------------------------------------------------

    # install pls package (if not already installed)
    install.packages("pls")

    # load pls package
    library(pls)
    library(plsdof)

# Preapre data -----------------------------------------------------------------

    # 1. Select vairables for model
    dt <- cbind(yarn[[2]], yarn[[1]])

    # 2. Scale
    # Means
    round(colMeans(dt), 3)

    # 3. Objects needed
    p <- ncol(dt) - 1
    n <- nrow(dt)

    # Variance
    round(apply(dt, 2, var), 3)

    # center and scale to variance = 1
    dt_sd <- cbind(yarn[[2]], scale(yarn[[1]]))

    # Checks
    round(colMeans(dt_sd), 3)
    round(apply(dt_sd, 2, var), 3)

# Estimation -------------------------------------------------------------------

    Xmj <- dt_sd[, -1]
    y <- dt_sd[, 1]
    M <- 10 # number of "components"
    y_hat <- cbind(mean(dt_sd[, 1]), matrix(rep(NA, n * (M - 1)), nrow = n))

    z <- matrix(NA, nrow = n, ncol = M + 1)
    theta_hat <- rep(NA, M)

    orthogonalize <- function(vec1, vec2) {
        v <- vec1
        u <- vec2

        newv <- v - drop(t(u) %*% v / (t(u) %*% u)) * u

        return(newv)
    }

    for (m in 2:M) {
        # 2a
        IP <- matrix(NA, nrow = n, ncol = p)
        for (j in 1:p) {
            IP[, j] <- t(Xmj[, j]) %*% y %*% Xmj[, j]
        }
        z[, m] <- rowSums(IP)

        # 2b
        theta_hat[m] <- drop(z[, m] %*% y / z[, m] %*% z[, m])

        # 2c
        y_hat[, m] <- y_hat[, m - 1] + theta_hat[m] * z[, m]

        # 2d orthogonalize all columns
        for (j in 1:p) {
            Xmj[, j] <- orthogonalize(Xmj[, j], z[, m])
        }
    }

    # fit PCR model
    model <- plsr(
        dt[, 1] ~ dt[, -1],
        ncomp = M,
        scale = FALSE,
        center = FALSE,
        method = "oscorespls",
        validation = "none"
    )

    round(as.data.frame(fitted(model)), 3)

    round(y_hat, 3)

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

    my.pls1 <- pls.model(X, y, m = 5, compute.DoF = TRUE)
    my.pls1

    my.pls1$DoF