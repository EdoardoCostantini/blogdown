# Project:   blogdown
# Objective: Study of superefficiency as presented by Rubin 1996 and Meng 1994
# Author:    Edoardo Costantini
# Created:   2022-04-14
# Modified:  2022-04-14

rm(list = ls())

# Pacakges
library(splitstackshape)

# Define population
N <- 1e4

dat <- data.frame(y = rnorm(N, mean = 0, sd = 1),
                 x = c(rep("a", N/2), rep("b", N/2)))

table(dat$x)


# Sample
n <- 100 # sample stratum size
dat_start <- stratified(indt = dat, group = "x", size = n)
smp_og <- smp <- as.data.frame(dat_start)

# Missing values
n0 <- .4 * n # number of missing cases in each group

smp[smp$x == "a", ][sample(1:n, n0), "y"] <- NA
smp[smp$x == "b", ][sample(1:n, n0), "y"] <- NA

# Estimands
y_bar_a <- mean(smp[smp$x == "a", "y"], na.rm = TRUE)
y_bar_b <- mean(smp[smp$x == "a", "y"], na.rm = TRUE)
y_bar <- (y_bar_a + y_bar_b)/2
D_bar <- y_bar_a - y_bar_b

# Multiply impute the data with the knowledge of groups having same distribution

Ri <- !is.na(smp$y)   # R = 1 observe, R = 0 missing
m <- 50
y_bar_obs <- mean(smp$y, na.rm = TRUE)

# Obtain multiply imputed ys
  y_imp <- lapply(1:m, function (i){
    # Sample imputaions
    y_imps <- rnorm(n    = sum(!Ri),
                    mean = 0,
                    sd   = 1)

    # Create copy of original y with missing values
    y_fill <- smp$y

    # Fill in the imputation
    y_fill[!Ri] <- y_imps

    # Return the imputed variable
    return(y_fill)
  })

# Perform a complete data procedure (Pcom) on every imputated dataset
  Pcoms <- sapply(y_imp, function (j){
    c(theta     = mean(j),        # mean of y in group x = 1
      var_theta = length(j)^(-1)) # standard error of mean of y in group x = 1
  } )

# Pool the estimates to obtain Pm
  theta_bar_m <- mean(Pcoms[1, ])
  U_bar_m     <- mean(Pcoms[2, ])
  Bm          <- 1/(m-1) * sum((Pcoms[1, ] - theta_bar_m)^2)
  Tm          <- U_bar_m + (1+1/m) * Bm
  Pm_unco     <- c(theta     = theta_bar_m,
                   var_theta = Tm)

# Estimates with comeplte data
  Pog <- c(theta     = mean(smp_og$y),
           var_theta = nrow(smp_og)^(-1))

data.frame(Pog = Pog,
           Pm_unco = Pm_unco)