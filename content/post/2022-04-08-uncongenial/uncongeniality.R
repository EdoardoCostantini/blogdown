# Project:   blogdown
# Objective: Description of main uncongenial situations
# Author:    Edoardo Costantini
# Created:   2022-04-08
# Modified:  2022-04-08

# Set up environement ----------------------------------------------------------

library(mice)
library(tidyverse)

# Uncongenial form 1: imputation model is less complex than true model ---------

# Meng 1994 --------------------------------------------------------------------

# Data
set.seed(20220408)
n <- 1e4
x <- sample(c(0, 1), size = n, replace = TRUE)
y <- rnorm(n, sd = 1) + 5 * x

# Variances in the groups are the same
tapply(X = y, INDEX = x, FUN = var)

# Group means are different
tapply(X = y, INDEX = x, FUN = mean)

# Impose MAR missingnes on y
ampute_out <- ampute(cbind(y, x),
                     patterns = c(0, 1),
                     mech = "MAR")
y_NAs <- ampute_out$amp[, "y"]

# Imputer assumes more ---------------------------------------------------------

# Analyst's parameter of interest estiamted on complete data

  Pog <- c(theta     = mean(y[x == 1]),
           var_theta = sum(x == 1)^(-1))

# Define useful objects

  Ri        <- !is.na(y_NAs)   # R = 1 observe, R = 0 missing
  n_exc     <- sum(Ri)         # number of obserserved cases on y
  y_bar_exc <- mean(y_NAs[Ri]) # mean of y using only observed cases

  # Group objects
  n_j0        <- sum(x == 0)              # sample size in x = 0
  n_j0_exc    <- sum(x == 0 & Ri)         # observed y size in x = 0
  n_j0_que    <- n_j0 - n_j0_exc          # missing y size in x = 0
  ybar_j0_exc <- mean(y_NAs[x == 0 & Ri]) # mean observed in x = 0

  n_j1        <- sum(x == 1)              # sample size in x = 1
  n_j1_exc    <- sum(x == 1 & Ri)         # observed y size in x = 1
  n_j1_que    <- n_j1 - n_j1_exc          # missing y size in x = 1
  ybar_j1_exc <- mean(y_NAs[x == 1 & Ri]) # mean observed in x = 1

# Incomplete data procedure (listwise deletion)

  Pobs <- c(theta = ybar_j1_exc,
            var_theta = n_j1_exc^(-1))

# Complete data procedure through (uncongenial) imputation

  # Define the number of imputations
  m <- 10

  # Obtain multiply imputed ys
  y_imp <- lapply(1:m, function (i){
    # Sample imputaions
    y_imps <- rnorm(n    = sum(!Ri),
                    mean = mean(y_NAs[Ri]),
                    sd   = 1)

    # Create copy of original y with missing values
    y_fill <- y_NAs

    # Fill in the imputation
    y_fill[!Ri] <- y_imps

    # Return the imputed variable
    return(y_fill)
  })

  # Perform a complete data procedure (Pcom) on every imputated dataset
  Pcoms <- sapply(y_imp, function (j){
    c(theta     = mean(j[x == 1]),        # mean of y in group x = 1
      var_theta = length(j[x == 1])^(-1)) # standard error of mean of y in group x = 1
  } )

  # Pool the estimates to obtain Pm
  theta_bar_m <- mean(Pcoms[1, ])
  U_bar_m     <- mean(Pcoms[2, ])
  Bm          <- 1/(m-1) * sum((Pcoms[1, ] - theta_bar_m)^2)
  Tm          <- U_bar_m + (1+1/m) * Bm
  Pm_unco     <- c(theta     = theta_bar_m,
                   var_theta = Tm)

# Compare results

  round(cbind(Pog, Pobs, Pm_unco), 5)

  # Super-efficiency
  # set.seed(20220408)
  # y <- rnorm(n, sd = 1) + 0 * x # group means are the same!
  # imputer knows!
  # analyst doesn't

# Imputer assumes less ---------------------------------------------------------

# Analyst's parameter of interest estiamted on complete data

  Pog <- c(theta     = mean(y),
           var_theta = n^(-1))

# > a) MCAR Data - Missing data is missing at completely at random -------------

# Impose missing values

  set.seed(20220411)
  ym <- ampute(cbind(y, x),
               patterns = c(0, 1),
               mech = "MCAR")$amp[, "y"]

# Compute all of the objects of interest

  Ri        <- !is.na(ym)   # R = 1 observe, R = 0 missing
  n_exc     <- sum(Ri)      # number of obserserved cases on y
  y_bar_exc <- mean(ym[Ri]) # mean of y using only observed cases

  n_j0        <- sum(x == 0)              # sample size in x = 0
  n_j0_exc    <- sum(x == 0 & Ri)         # observed y size in x = 0
  n_j0_que    <- n_j0 - n_j0_exc          # missing y size in x = 0
  ybar_j0_exc <- mean(ym[x == 0 & Ri]) # mean observed in x = 0

  n_j1        <- sum(x == 1)              # sample size in x = 1
  n_j1_exc    <- sum(x == 1 & Ri)         # observed y size in x = 1
  n_j1_que    <- n_j1 - n_j1_exc          # missing y size in x = 1
  ybar_j1_exc <- mean(ym[x == 1 & Ri]) # mean observed in x = 1

# Incomplete data procedure

  # Incomplete data procedure 1: unkown x
  Pobs_1 <- c(Qh = y_bar_exc,
              U  = n_exc^(-1))

    # the computation of the mean is equivalent to computing
    (n_j1_exc / n_exc * ybar_j1_exc + n_j0_exc / n_exc * ybar_j0_exc) - Pobs_1["Qh"]

  # Incomplete data procedure 2: known x
  Pobs_2 <- c(Qh = n_j1 / n * ybar_j1_exc + n_j0 / n * ybar_j0_exc,
              U  = n^(-1))

# Complete data procedure with imputations under more general model

  # Impute values
  y_imp <- lapply(1:m, function (i){

    # For x = 1
    theta1 <- rnorm(1, ybar_j1_exc, n_j1_exc^(-1))
    y_imps_j1 <- rnorm(sum(!Ri & x == 1), mean = theta1, sd = 1)

    # For x = 0
    theta0 <- rnorm(1, ybar_j0_exc, n_j0_exc^(-1))
    y_imps_j0 <- rnorm(sum(!Ri & x == 0), mean = theta0, sd = 1)

    # Fill missing vlaues
    y_fill <- ym
    y_fill[x == 0 & !Ri] <- y_imps_j0
    y_fill[x == 1 & !Ri] <- y_imps_j1

    # Returns
    return(y_fill)
  })

  # Perform a complete data procedure (Pcom) on every imputated dataset
  Pcoms <- sapply(y_imp, function (j){
    c(Qh = mean(j),        # mean of y in the population
      U  = length(j)^(-1)) # standard error of population mean of y
  } )

  # Pool the estimates to obtain Pm
  theta_bar_m <- mean(Pcoms[1, ])
  U_bar_m <- mean(Pcoms[2, ])
  Bm <- 1/(m-1) * sum((Pcoms[1, ] - theta_bar_m)^2)
  Tm <- U_bar_m + (1+1/m) * Bm
  Pm_cong <- c(Qh = theta_bar_m,
               U = Tm)

# Compare results

  MCAR_results <- round(cbind(Pog, Pobs_1, Pobs_2, Pm_cong), 5)

# > b) MAR Data - Missing values are missing at random (depend on X) -----------

# Impose missing values

  set.seed(20220411)
  ym <- ampute(cbind(y, x),
               patterns = c(0, 1),
               mech = "MAR")$amp[, "y"] # MAR instead of MCAR

# Compute all of the objects of interest

  Ri        <- !is.na(ym)   # R = 1 observe, R = 0 missing
  n_exc     <- sum(Ri)      # number of obserserved cases on y
  y_bar_exc <- mean(ym[Ri]) # mean of y using only observed cases

  n_j0        <- sum(x == 0)              # sample size in x = 0
  n_j0_exc    <- sum(x == 0 & Ri)         # observed y size in x = 0
  n_j0_que    <- n_j0 - n_j0_exc          # missing y size in x = 0
  ybar_j0_exc <- mean(ym[x == 0 & Ri]) # mean observed in x = 0

  n_j1        <- sum(x == 1)              # sample size in x = 1
  n_j1_exc    <- sum(x == 1 & Ri)         # observed y size in x = 1
  n_j1_que    <- n_j1 - n_j1_exc          # missing y size in x = 1
  ybar_j1_exc <- mean(ym[x == 1 & Ri]) # mean observed in x = 1

# Incomplete data procedure

  # Incomplete data procedure 1: unkown x
  Pobs_1 <- c(Qh = y_bar_exc,
              U  = n_exc^(-1))

    # the computation of the mean is equivalent to computing
    (n_j1_exc / n_exc * ybar_j1_exc + n_j0_exc / n_exc * ybar_j0_exc) - Pobs_1["Qh"]

  # Incomplete data procedure 2: known x
  Pobs_2 <- c(Qh = n_j1 / n * ybar_j1_exc + n_j0 / n * ybar_j0_exc,
              U  = n^(-1))

# Complete data procedure with imputations under more general model

  # Impute values
  y_imp <- lapply(1:m, function (i){

    # For x = 1
    theta1 <- rnorm(1, ybar_j1_exc, n_j1_exc^(-1))
    y_imps_j1 <- rnorm(sum(!Ri & x == 1), mean = theta1, sd = 1)

    # For x = 0
    theta0 <- rnorm(1, ybar_j0_exc, n_j0_exc^(-1))
    y_imps_j0 <- rnorm(sum(!Ri & x == 0), mean = theta0, sd = 1)

    # Fill missing vlaues
    y_fill <- ym
    y_fill[x == 0 & !Ri] <- y_imps_j0
    y_fill[x == 1 & !Ri] <- y_imps_j1

    # Returns
    return(y_fill)
  })

  # Perform a complete data procedure (Pcom) on every imputated dataset
  Pcoms <- sapply(y_imp, function (j){
    c(Qh = mean(j),        # mean of y in the population
      U  = length(j)^(-1)) # standard error of population mean of y
  } )

  # Pool the estimates to obtain Pm
  theta_bar_m <- mean(Pcoms[1, ])
  U_bar_m <- mean(Pcoms[2, ])
  Bm <- 1/(m-1) * sum((Pcoms[1, ] - theta_bar_m)^2)
  Tm <- U_bar_m + (1+1/m) * Bm
  Pm_cong <- c(Qh = theta_bar_m,
               U = Tm)

# Compare results

  MAR_results <- round(cbind(Pog, Pobs_1, Pobs_2, Pm_cong), 5)

# > When is what statistically valid? ------------------------------------------
MCAR_results
MAR_results