# Pls function manual

pls.manual <- function(ivs, dv, m = 1L){

    # Parms
    # M   <- 10 # number of "components"
    # ivs <- yarn[[1]]
    # dv <- yarn[[2]]

    # Scale data
    M <- m + 1
    p <- ncol(ivs)
    n <- nrow(ivs)
    X      <- lapply(1:M, matrix, nrow = n, ncol = p)
    X[[1]] <- scale(ivs)
    y      <- dv
    y_hat <- cbind(
        mean(y),
        matrix(rep(NA, n * (M - 1)), nrow = n)
    )
    z         <- matrix(NA, nrow = n, ncol = M)
    theta_hat <- rep(NA, M)
    W <- matrix(nrow = ncol(ivs), ncol = M)

    # PLS Algorithm following HastieEtAl2017 p 81 (Algorithm 3.3)
    for (m in 2:M) {
        # 2a
        store_2a <- matrix(NA, nrow = n, ncol = p)
        for (j in 1:p) {
            rho_hat_mj <- t(X[[m - 1]][, j]) %*% y
            store_2a[, j] <- rho_hat_mj %*% X[[m - 1]][, j]
            W[j, m] <- rho_hat_mj
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

    # Normalize the component scores
    Tsn <- apply(z[, -1], 2, function(j) j / sqrt(sum(j^2)))

    # Return
    return(
        list(
            Ts = z[, -1],
            Tsn = Tsn,
            Yhat = y_hat,
            W = W[, -1]
        )
    )

}

# Orthogonalize two vectors
orthogonalize <- function(vec1, vec2) {
    v <- vec1
    u <- vec2

    newv <- v - drop(t(u) %*% v / (t(u) %*% u)) * u

    return(newv)
}