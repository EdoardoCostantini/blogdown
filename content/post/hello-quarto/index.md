---
title: Hello, Quarto
date: "2012-04-06"
slug: quarto
categories: ["Drafts", "Tutorials"]
tags: ["statistics", "regression", "penalty", "high-dimensional data"]
subtitle: ''
summary: ''
authors: ["admin"]
format: 
  hugo:
    toc: true
    toc-depth: 3
    number-sections: true
---



-   <a href="#introduction" id="toc-introduction"><span class="toc-section-number">1</span> <span class="header-section-number">1</span> Introduction</a>
-   <a href="#learn-by-coding" id="toc-learn-by-coding"><span class="toc-section-number">2</span> <span class="header-section-number">2</span> Learn by coding</a>
    -   <a href="#fitting-ridge-regression-manually" id="toc-fitting-ridge-regression-manually"><span class="toc-section-number">2.1</span> <span class="header-section-number">2.1</span> Fitting ridge regression manually</a>
        -   <a href="#an-alternative-way-to-avoid-penalising-the-intercept" id="toc-an-alternative-way-to-avoid-penalising-the-intercept"><span class="toc-section-number">2.1.1</span> <span class="header-section-number">2.1.1</span> An alternative way to avoid penalising the intercept</a>
    -   <a href="#fit-ridge-regression-with-glmnet" id="toc-fit-ridge-regression-with-glmnet"><span class="toc-section-number">2.2</span> <span class="header-section-number">2.2</span> Fit ridge regression with glmnet</a>
        -   <a href="#use-the-biased-estimation-of-variance" id="toc-use-the-biased-estimation-of-variance"><span class="toc-section-number">2.2.1</span> <span class="header-section-number">2.2.1</span> Use the biased estimation of variance</a>
        -   <a href="#return-the-unstandardized-coefficients" id="toc-return-the-unstandardized-coefficients"><span class="toc-section-number">2.2.2</span> <span class="header-section-number">2.2.2</span> Return the unstandardized coefficients</a>
        -   <a href="#adjust-the-parametrization-of-lambda-for-glmnet" id="toc-adjust-the-parametrization-of-lambda-for-glmnet"><span class="toc-section-number">2.2.3</span> <span class="header-section-number">2.2.3</span> Adjust the parametrization of <img style="vertical-align:middle" src="https://latex.codecogs.com/svg.latex?%5Ctextstyle%20%5Clambda" alt="\lambda" title="\lambda" class="math inline" /> for <code>glmnet</code></a>
        -   <a href="#compare-manual-and-glmnet-ridge-regression-output" id="toc-compare-manual-and-glmnet-ridge-regression-output"><span class="toc-section-number">2.2.4</span> <span class="header-section-number">2.2.4</span> Compare manual and <code>glmnet</code> ridge regression output</a>
-   <a href="#tldr-just-give-me-the-code" id="toc-tldr-just-give-me-the-code"><span class="toc-section-number">3</span> <span class="header-section-number">3</span> TL;DR, just give me the code!</a>

# <span class="header-section-number">1</span> Introduction

When there are many correlated predictors in a linear regression model, their regression coefficients can become poorly determined and exhibit high variance.
This problem can be alleviated by imposing a size constraint (or penalty) on the coefficients.
Ridge regression shrinks the regression coefficients by imposing a penalty on their size.
The ridge coefficients values minimize a penalized residual sum of squares:

![\hat{\beta}^{\text{ridge}} = \text{argmin}\_{\beta} \left\\{ \sum\_{i=1}^{N} \left( y_i - \beta_0 - \sum\_{j=1}^{p} x\_{ij}\beta_j \right)^2 + \lambda \sum\_{j=1}^{p}\beta_j^2 \right\\}](https://latex.codecogs.com/svg.latex?%5Chat%7B%5Cbeta%7D%5E%7B%5Ctext%7Bridge%7D%7D%20%3D%20%5Ctext%7Bargmin%7D_%7B%5Cbeta%7D%20%5Cleft%5C%7B%20%5Csum_%7Bi%3D1%7D%5E%7BN%7D%20%5Cleft%28%20y_i%20-%20%5Cbeta_0%20-%20%5Csum_%7Bj%3D1%7D%5E%7Bp%7D%20x_%7Bij%7D%5Cbeta_j%20%5Cright%29%5E2%20%2B%20%5Clambda%20%5Csum_%7Bj%3D1%7D%5E%7Bp%7D%5Cbeta_j%5E2%20%5Cright%5C%7D "\hat{\beta}^{\text{ridge}} = \text{argmin}_{\beta} \left\{ \sum_{i=1}^{N} \left( y_i - \beta_0 - \sum_{j=1}^{p} x_{ij}\beta_j \right)^2 + \lambda \sum_{j=1}^{p}\beta_j^2 \right\}")

The ridge solutions are not equivariant under scaling of the inputs.
Therefore, it is recommended to standardize the inputs before solving the minimization problem.

Notice that the intercept ![\beta_0](https://latex.codecogs.com/svg.latex?%5Cbeta_0 "\beta_0") has been left out of the penalty term.
Penalization of the intercept would make the procedure depend on the origin chosen for ![Y](https://latex.codecogs.com/svg.latex?Y "Y").
Furthermore, by centering the predictors, we can separate the solution to the [minimazion problem](https://www.notion.so/Ridge-regression-8134d8babda5413ab182df645c6196a8) into two parts:

1.  Intercept

![\hat{\beta}\_0 = \bar{y}=\frac{1}{N}\sum\_{i = 1}^{N} y_i](https://latex.codecogs.com/svg.latex?%5Chat%7B%5Cbeta%7D_0%20%3D%20%5Cbar%7By%7D%3D%5Cfrac%7B1%7D%7BN%7D%5Csum_%7Bi%20%3D%201%7D%5E%7BN%7D%20y_i "\hat{\beta}_0 = \bar{y}=\frac{1}{N}\sum_{i = 1}^{N} y_i")

2.  Penalised regression coefficients

![\hat{\beta}^{\text{ridge}}=(\mathbf{X}^T\mathbf{X} + \lambda \mathbf{I})^{-1}\mathbf{X}^Ty](https://latex.codecogs.com/svg.latex?%5Chat%7B%5Cbeta%7D%5E%7B%5Ctext%7Bridge%7D%7D%3D%28%5Cmathbf%7BX%7D%5ET%5Cmathbf%7BX%7D%20%2B%20%5Clambda%20%5Cmathbf%7BI%7D%29%5E%7B-1%7D%5Cmathbf%7BX%7D%5ETy "\hat{\beta}^{\text{ridge}}=(\mathbf{X}^T\mathbf{X} + \lambda \mathbf{I})^{-1}\mathbf{X}^Ty")

which is the regular way of estimating regression coefficients with a penalty term (![\lambda](https://latex.codecogs.com/svg.latex?%5Clambda "\lambda")) added on the diagonal (![\mathbf{I}](https://latex.codecogs.com/svg.latex?%5Cmathbf%7BI%7D "\mathbf{I}")) of the cross-product matrix (![\mathbf{X}^T\mathbf{X}](https://latex.codecogs.com/svg.latex?%5Cmathbf%7BX%7D%5ET%5Cmathbf%7BX%7D "\mathbf{X}^T\mathbf{X}")) to make it invertible (![(\...)^{-1}](https://latex.codecogs.com/svg.latex?%28...%29%5E%7B-1%7D "(...)^{-1}")).

# <span class="header-section-number">2</span> Learn by coding

The `glmnet` package can be used to obtain the ridge regression estimates of the regression coefficients.
In this section, we will first see how to obtain these estimates "manually", that is coding every step on our own, and then we will see how to obtain the same results using the `glmnet` package.

Let's start by setting up the R environment.
In this post, we will work with the `mtcars` data.
If you are not familiar with it, just look up the R help file on it.
We will use the first column of the dataset (variable named `mpg`) as a dependent variable and the remaining ones as predictors in a linear regression.

``` r
# Set up -----------------------------------------------------------------------

# Load packages
library(glmnet)

# Take the mtcars data
y <- mtcars[, "mpg"]
X <- mtcars[, -1]

# Create a few shorthands we will use
n <- nrow(X)
p <- ncol(X)
```

## <span class="header-section-number">2.1</span> Fitting ridge regression manually

First, let's make sure the predictors are centered on the mean and scaled to have a variance of 1.

``` r
# Fitting ridge regression manually --------------------------------------------

# Scale the data (standardize)
X_scale <- scale(X, center = TRUE, scale = TRUE)
```

Then, we want to **fit the ridge regression** manually by separating the intercept and the regression coefficients estimation (two-step approach):

1.  Estimate the intercept (![\hat{\beta}\_0](https://latex.codecogs.com/svg.latex?%5Chat%7B%5Cbeta%7D_0 "\hat{\beta}_0"))

``` r
# Estimate the intercept
b0_hat_r <- mean(y)
```

2.  Estimate the ridge regression coefficients (![\hat{\beta}^{\text{ridge}}](https://latex.codecogs.com/svg.latex?%5Chat%7B%5Cbeta%7D%5E%7B%5Ctext%7Bridge%7D%7D "\hat{\beta}^{\text{ridge}}")).

<!-- -->

1.  Compute the cross-product matrix of the predictors.

    This is the same step we would take if we wanted to compute the OLS estimates.

    ::: {.cell}

    ``` r
    # Compute the cross-product matrix of the data
    XtX <- t(X_scale) %*% X_scale
    ```

    :::

2.  Define a value of ![\lambda](https://latex.codecogs.com/svg.latex?%5Clambda "\lambda").

    This value is usually chosen by cross-validation from a grid of possible values.
    However, here we are only interested in how ![\lambda](https://latex.codecogs.com/svg.latex?%5Clambda "\lambda") is used in the computation, so we can simply give it a fixed value.

    ::: {.cell}

    ``` r
    # Define a lambda value
    lambda <- .1
    ```

    :::

3.  Compute ![\hat{\beta}^{\text{ridge}}](https://latex.codecogs.com/svg.latex?%5Chat%7B%5Cbeta%7D%5E%7B%5Ctext%7Bridge%7D%7D "\hat{\beta}^{\text{ridge}}").

    ::: {.cell}

    ``` r
    # Estimate the regression coefficients with the ridge penalty
    bs_hat_r <- solve(XtX + lambda * diag(p)) %*% t(X_scale) %*% y
    ```

    :::

    where `diag(p)` is the identity matrix ![\mathbf{I}](https://latex.codecogs.com/svg.latex?%5Cmathbf%7BI%7D "\mathbf{I}").

Finally, let's print the results:

``` r
# Print the results
round(
  data.frame(twostep = c(b0 = b0_hat_r,
                         b = bs_hat_r)),
  3
)
```

        twostep
    b0   20.091
    b1   -0.194
    b2    1.366
    b3   -1.373
    b4    0.438
    b5   -3.389
    b6    1.361
    b7    0.162
    b8    1.243
    b9    0.496
    b10  -0.460

It is important to note the effect of centering and scaling.
When fitting ridge regression, many sources recommend centering the data.
This allows to separate the estimation of the intercept from the estimation of the regression coefficients.
As a result, only the regression coefficients are penalised.
To understand the effect of centering, consider what happens in regular OLS estimation when **predictors are centered**:

``` r
# Centering in regular OLS -----------------------------------------------------

# Create a version of X that is centered
X_center <- scale(X, center = TRUE, scale = FALSE)

# Fit an regular linear model
lm_ols <- lm(y ~ X_center)

# Check that b0 is equal to the mean of y
coef(lm_ols)["(Intercept)"] - mean(y)
```

      (Intercept) 
    -3.552714e-15 

Furthermore, let's see what would have happened if we had penalised the intercept as well.

``` r
# Consequence of penalising the intercept --------------------------------------

# Add a vector of 1s to penalise the intercept
X_scale_w1 <- cbind(1, X_scale)

# Compute the cross-product matrix of the data
XtX <- t(X_scale_w1) %*% X_scale_w1

# Estimate the regression coefficients with the ridge penalty
bs_hat_r_w1 <- solve(XtX + lambda * diag(p+1)) %*% t(X_scale_w1) %*% y

# Print the results
round(
  data.frame(twostep = c(b0 = b0_hat_r,
                         b = bs_hat_r),
             onestep = c(b0 = bs_hat_r_w1[1],
                         b = bs_hat_r_w1[-1])),
  3
)
```

        twostep onestep
    b0   20.091  20.028
    b1   -0.194  -0.194
    b2    1.366   1.366
    b3   -1.373  -1.373
    b4    0.438   0.438
    b5   -3.389  -3.389
    b6    1.361   1.361
    b7    0.162   0.162
    b8    1.243   1.243
    b9    0.496   0.496
    b10  -0.460  -0.460

As you see, the intercept would be shrunk toward zero, without any benefit.
As a result, any prediction would also be offset by the same amount.

### <span class="header-section-number">2.1.1</span> An alternative way to avoid penalising the intercept

It can be handy to obtain estimates of the regression coefficients and intercept in one step.
We can use matrix algebra and R to simplify the two-step procedure to a single step.
In particular, we can avoid the penalisation of the intercept by setting to 0 the first element of the "penalty" matrix `lambda * diag(p + 1)`.

``` r
# Alternative to avoid penalization of the intercept ---------------------------

# Compute cross-product matrix
XtX <- crossprod(X_scale_w1)

# Create penalty matrix
pen <- lambda * diag(p + 1)

# replace first element with 0
pen[1, 1] <- 0

# Obtain standardized estimates
bs_hat_r2 <- solve(XtX + pen) %*% t(X_scale_w1) %*% (y)

# Compare
round(
        data.frame(
                twostep = c(b0 = b0_hat_r, b = bs_hat_r),
                onestep = c(b0 = bs_hat_r2[1], b = bs_hat_r2[-1])
        ),
        3
)
```

        twostep onestep
    b0   20.091  20.091
    b1   -0.194  -0.194
    b2    1.366   1.366
    b3   -1.373  -1.373
    b4    0.438   0.438
    b5   -3.389  -3.389
    b6    1.361   1.361
    b7    0.162   0.162
    b8    1.243   1.243
    b9    0.496   0.496
    b10  -0.460  -0.460

## <span class="header-section-number">2.2</span> Fit ridge regression with glmnet

The most popular R package to fit regularised regression is `glmnet`.
Let's see how we can replicate the results we obtained with the manual approach with glmnet.
There are three important differences to consider:

-   `glmnet` uses the [biased sample variance estimate](https://en.wikipedia.org/wiki/Variance#Biased_sample_variance) when scaling the predictors;
-   `glmnet` returns the unstandardized regression coefficients;
-   `glmnet` uses a different parametrization for ![\lambda](https://latex.codecogs.com/svg.latex?%5Clambda "\lambda").

To obtain the same results with the manual approach and `glmnet` we need to account for these.

### <span class="header-section-number">2.2.1</span> Use the biased estimation of variance

First, let's use the biased sample variance estimate in computing ![\hat{\beta}^{\text{ridge}}](https://latex.codecogs.com/svg.latex?%5Chat%7B%5Cbeta%7D%5E%7B%5Ctext%7Bridge%7D%7D "\hat{\beta}^{\text{ridge}}") with the manual approach:

``` r
# Fitting ridge manually with biased variance estimation -----------------------

# Standardize X
X_scale <- sapply(1:p, function (j){
  muj <- mean(X[, j])                  # mean
  sj <- sqrt( var(X[, j]) * (n-1) / n) # (biased) sd
  (X[, j] - muj) / sj                  # center and scale
})

# Craete the desing matrix
X_scale_dm <- cbind(1, X_scale)

# Compute cross-product matrix
XtX <- crossprod(X_scale_dm)

# Create penalty matrix
pen <- lambda * diag(p + 1)
pen[1, 1] <- 0

# Obtain standardized estimates
bs_hat_r3 <- solve(XtX + pen) %*% t(X_scale_dm) %*% (y)

# Print results
round(
      data.frame(
              manual = c(b0 = bs_hat_r3[1], b = bs_hat_r3[-1])
      ),
      3
)
```

        manual
    b0  20.091
    b1  -0.191
    b2   1.353
    b3  -1.354
    b4   0.430
    b5  -3.343
    b6   1.343
    b7   0.159
    b8   1.224
    b9   0.488
    b10 -0.449

### <span class="header-section-number">2.2.2</span> Return the unstandardized coefficients

Next, we need to revert these regression coefficients to their original scale.
Since we are estimating the regression coefficients on the scaled data, they are computed on the standardized scale.

``` r
# Return the  unstandardized coefficients --------------------------------------

# Extract the original mean and standard deviations of all X variables
mean_x <- colMeans(X)
sd_x <- sqrt(apply(X, 2, var) * (n - 1) / n) # biased version

# Revert to original scale
bs_hat_r4 <- c(bs_hat_r3[1] - crossprod(mean_x, bs_hat_r3[-1] / sd_x),
               bs_hat_r3[-1] / sd_x)

# Compare manual standardized and unstandardized results
round(
      data.frame(
              standardized = c(b0 = bs_hat_r3[1], b = bs_hat_r3[-1]),
              unstandardized = c(b0 = bs_hat_r4[1], b = bs_hat_r4[-1])
      ),
      3
)
```

        standardized unstandardized
    b0        20.091         12.908
    b1        -0.191         -0.109
    b2         1.353          0.011
    b3        -1.354         -0.020
    b4         0.430          0.818
    b5        -3.343         -3.471
    b6         1.343          0.764
    b7         0.159          0.320
    b8         1.224          2.491
    b9         0.488          0.672
    b10       -0.449         -0.282

### <span class="header-section-number">2.2.3</span> Adjust the parametrization of ![\lambda](https://latex.codecogs.com/svg.latex?%5Clambda "\lambda") for `glmnet`

Next, we need to understand the relationship between the ![\lambda](https://latex.codecogs.com/svg.latex?%5Clambda "\lambda") parametrization we used and the one used by `glmnet`.
The following code shows that if we want to use a given value of `lambda` in `glmnet` we need to multiply it by the standard deviation of the dependent variable (`sd_y`) and divide it by the sample size (`n`).

``` r
# Adjust the parametrization of lambda -----------------------------------------

# Extract the original mean and standard deviations of y (for lambda parametrization)
mean_y <- mean(y)
sd_y <- sqrt(var(y) * (n - 1) / n)

# Compute the value glmnet wants for your target lambda
lambda_glmnet <- sd_y * lambda / n
```

### <span class="header-section-number">2.2.4</span> Compare manual and `glmnet` ridge regression output

Finally, we can compare the results:

``` r
# Fitting ridge regression with glmnet -----------------------------------------

# Fit glmnet
fit_glmnet_s <- glmnet(x = X,
                       y = y,
                       alpha = 0,
                       lambda = lambda_glmnet, # correction for how penalty is used
                       thresh = 1e-20)

bs_glmnet <- coef(fit_glmnet_s)

# Compare estimated coefficients
round(
      data.frame(
        manual = c(b0 = bs_hat_r4[1], b = bs_hat_r4[-1]),
        glmnet = c(b0 = bs_glmnet[1], b = bs_glmnet[-1])
      ),
      3
)
```

           manual glmnet
    b0     12.908 12.908
    b.cyl  -0.109 -0.109
    b.disp  0.011  0.011
    b.hp   -0.020 -0.020
    b.drat  0.818  0.818
    b.wt   -3.471 -3.471
    b.qsec  0.764  0.764
    b.vs    0.320  0.320
    b.am    2.491  2.491
    b.gear  0.672  0.672
    b.carb -0.282 -0.282

# <span class="header-section-number">3</span> TL;DR, just give me the code!
