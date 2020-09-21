
# allcoeff <- list()
# class(allcoeff) <- "linreg"

linreg <- function(formula, data)
{
  X <- model.matrix(formula, data)
  y1 <- data[all.vars(formula)[1]]
  y <- as.matrix(y1)

  #1.2.1  Computations using ordinary least squares
  #Regression Coefficients:
  beta_a <- solve(t(X) %*% X)
  beta_b <- (t(X)) %*% y
  beta_cap <- (beta_a) %*% (beta_b)

  print(beta_cap)

  #The fitted values:
  y_cap <- X %*% beta_cap

  print(y_cap)

  #The residuals:
  e_cap <- y - y_cap

  print(e_cap)

  #The degrees of freedom:
  n <- nrow(X)
  p <- ncol(X)
  df <- n - p

  print(e_cap)

  #The residual variance
  sigmasq_cap_a <- (t(e_cap) %*% e_cap)
  sigmasq_cap <-  sigmasq_cap_a/df

  print(sigmasq_cap)

  #The variance of the regression coefficients:
  Varofbeta_cap <- (as.numeric(sigmasq_cap)) * beta_a

  print(Varofbeta_cap)

  # The t-values for each coefficient:
  t_beta <- beta_cap/(sqrt(diag(Varofbeta_cap)))

  pt_beta <- 2*pt(abs(t_beta), df)

  print(pt_beta)

  allcoeff <- list(Coefficients=beta_cap, fitted_values=y_cap, residuals=e_cap, degreeoffreedom=df, residvariance=sigmasq_cap, varofregcoeff=Varofbeta_cap, tvalue=pt_beta)
  attr(allcoeff,"class") <- "linreg"
  return(allcoeff)
}

linreg(formula = Petal.Length~Species, data = iris)

#Defining Methods

  print.linreg <- function(x){
  aa <- x$Coefficients
  print(mod_object)

}

plot.linreg <- function(x){

  mod_object <- lm(formula, data)
  print(mod_object)

}

xxx <- linreg(formula = Petal.Length~Species, data = iris)
print(xxx)


