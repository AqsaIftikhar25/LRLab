#' Linear Regression
#' @name Linear Regression
#'
#' @param formula linear expression
#' @param data data source
#'
#' @description calculate coefficients of linear regression
#'
#'
#' @return return the coefficients and coefficient names
#' @usage linreg (formula,data)
#'
#' @examples
#' linreg(Petal.Length~Species, data = iris)
#' linreg(mpg~cyl, data = mtcars)
#' @export


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

  #print(beta_cap)

  #The fitted values:
  y_cap <- X %*% beta_cap

  #print(y_cap)

  #The residuals:
  e_cap <- as.vector(y - y_cap)

  #print(e_cap)

  #The degrees of freedom:
  n <- nrow(X)
  p <- ncol(X)
  df <- n - p

  #print(df)

  #The residual variance
  sigmasq_cap_a <- (t(e_cap) %*% e_cap)
  sigmasq_cap <-  sigmasq_cap_a/df

  #print(sigmasq_cap)

  #The variance of the regression coefficients:
  Varofbeta_cap <- (as.numeric(sigmasq_cap)) * beta_a

  #print(Varofbeta_cap)

  # The t-values for each coefficient:
  t_beta <- beta_cap/(sqrt(diag(Varofbeta_cap)))

  # print(t_beta)
  pt_beta <- 2*pt(abs(t_beta), df)

  #print(pt_beta)

  allcoeff <- list(Coefficients=beta_cap, fitted_values=y_cap, residuals=e_cap, degreeoffreedom=df,
                residvariance=sigmasq_cap, varofregcoeff=Varofbeta_cap,
                t_value=t_beta, p_value=pt_beta, formula1=formula, data1=data)
  # attr(allcoeff,"class") <- "linreg"
  class(allcoeff) <- 'linreg'
  return(allcoeff)
}

xx = linreg(formula = Petal.Length~Species, data = iris)

#Defining Methods

print.linreg <- function(x){

  coeff = as.vector(x$Coefficients)
  names(coeff) = rownames(x$Coefficients)
  cat("Call:\n", paste(x$formula1 , sep = "\n"),"\n\n", sep = "")
  # cat("linreg(formula = Petal.Length ~ Sepal.Width + Sepal.Length, data = iris)", '\n')

  if (length(x$Coefficients)) {
    cat("Coefficients:\n")
    print.default(format(coeff), print.gap = 2L,quote = FALSE)
  }
  else
  {
    cat("No coefficients\n")
  }


}


plot.linreg <- function(x){

  stand_resid = sqrt(abs(x$residuals/sd(x$residuals)))
  df = data.frame(x$fitted_values,x$residuals,stand_resid)
  colnames(df) = c('x', 'y','y1')

  p1 = ggplot(data=df,mapping = aes(x,y)) + geom_point(size = 5, shape = 1)+
    stat_summary(fun = median, color = 'red', geom = 'line', size = 1)+
    labs(y= "Residuals", x = "Fitted values \n lm(Petal.Length ~ Species)", title = "Residuals vs Fitted")+
    theme(plot.title = element_text(hjust = 0.5))

  p2 = ggplot(data=df,mapping = aes(x,y1)) + geom_point(size = 5, shape = 1)+
    stat_summary(fun = mean, color = 'red', geom = 'line', size = 1)+
    labs(y= expression(sqrt('|Standardized Residuals|')), x = "Fitted values \n lm(Petal.Length ~ Species)", title = "Scale-Location")+
    theme(plot.title = element_text(hjust = 0.5))

  print(p1)
  print(p2)
}

resid.linreg <- function(x)
  {
  return(as.vector(x$residuals))
  }

pred <- function(x)
{
  UseMethod("pred")
}

pred.linreg <- function(x)
{
  return(as.vector(x$fitted_values))
}

coef <- function(x)
{
  UseMethod("coef")
}

coef.linreg <- function(x)
{
  cat("Coefficients:",'\n')
  coeff <- as.vector(x$Coefficients)
  names(coeff) = rownames(x$Coefficients)
  return(coeff)
}

summary.linreg <- function(x, formula)
{
  emptyvect <- c("***", "***", "***")
  coef_matrix <- cbind(x$Coefficients, sqrt(diag(x$varofregcoeff)), x$t_value, x$p_value, emptyvect)
  colnames(coef_matrix) <- c("Coefficients", "Standard Error", "t values", "p values", " ")
  print.table(coef_matrix)

  cat('Residual standard error:', sqrt(x$residvariance))
  cat(' on ')
  cat(x$degreeoffreedom,"degrees of freedom\n")


}

library(ggplot2)

print(xx)
plot(xx)
resid(xx)
pred(xx)
coef(xx)
summary(xx)
