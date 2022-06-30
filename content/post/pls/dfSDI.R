# Degrees of freedom for supervised derived input models

dofPLS <- function(X, y, TT, Yhat, m = ncol(X), DoF.max = ncol(X) + 1){
    # Example inputs
    # X = scale(mtcars[, -1])
    # y = mtcars[, 1]
    # m = ncol(X)
    # DoF.max = m + 1
    # TT <- linear.pls.fit(X, y, m, DoF.max = DoF.max)$TT # normalizezs PC scores
    # Yhat <- linear.pls.fit(X, y, m, DoF.max = DoF.max)$Yhat[, 2:(m + 1)]

    # Body
    n <- nrow(X)

    # Scale data
    mean.X <- apply(X, 2, mean)
    sd.X <- apply(X, 2, sd)
    sd.X[sd.X == 0] <- 1
    X <- X - rep(1, nrow(X)) %*% t(mean.X)
    X <- X / (rep(1, nrow(X)) %*% t(sd.X))
    K <- X %*% t(X)

    # pls.dof
    DoF.max <- DoF.max - 1
    TK <- matrix(, m, m)
    KY <- krylov(K, K %*% y, m)
    lambda <- eigen(K)$values
    tr.K <- vector(length = m)
    for (i in 1:m) {
        tr.K[i] <- sum(lambda^i)
    }
    BB <- t(TT) %*% KY
    BB[row(BB) > col(BB)] <- 0
    b <- t(TT) %*% y
    DoF <- vector(length = m)
    Binv <- backsolve(BB, diag(m))
    tkt <- rep(0, m)
    ykv <- rep(0, m)
    KjT <- array(dim = c(m, n, m))
    dummy <- TT
    for (i in 1:m) {
        dummy <- K %*% dummy
        KjT[i, , ] <- dummy
    }
    trace.term <- rep(0, m)

    for (i in 1:m) {
        Binvi <- Binv[1:i, 1:i, drop = FALSE]
        ci <- Binvi %*% b[1:i]
        Vi <- TT[, 1:i, drop = FALSE] %*% t(Binvi)
        trace.term[i] <- sum(ci * tr.K[1:i])
        ri <- y - Yhat[, i]
        for (j in 1:i) {
            KjTj <- KjT[j, , ]
            tkt[i] <- tkt[i] + ci[j] * tr(t(TT[, 1:i, drop = FALSE]) %*%
                KjTj[, 1:i, drop = FALSE])
            ri <- K %*% ri
            ykv[i] <- ykv[i] + sum(ri * Vi[, j])
        }
    }

    DoF <- trace.term + 1:m - tkt + ykv

    DoF[DoF > DoF.max] <- DoF.max
    DoF <- c(0, DoF) + 1
    DoF

}
