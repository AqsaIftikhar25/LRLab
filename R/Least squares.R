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
  e_cap <- y - y_cap

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

  print(t_beta)
  pt_beta <- 2*pt(abs(t_beta), df)

  #print(pt_beta)

  allcoeff <- list(Coefficients=beta_cap, fitted_values=y_cap, residuals=e_cap, degreeoffreedom=df, residvariance=sigmasq_cap, varofregcoeff=Varofbeta_cap, p_value=pt_beta)
  attr(allcoeff,"class") <- "linreg"
  return(allcoeff)
}

xx = linreg(formula = Petal.Length~Species, data = iris)

#Defining Methods

print.linreg <- function(x){
  "print function like lm print"
  vec_coe = as.vector(x$Coefficients)
  names(vec_coe) = rownames(x$Coefficients)
  #cat("\nCall:\n", paste(deparse(.self$call), sep = "\n", collapse = "\n"),
  #  "\n\n", sep = "")
  if (length(x$Coefficients)) {
    cat("Coefficients:\n")
    print.default(format(vec_coe), print.gap = 2L,quote = FALSE)
  }
  else
  {
    cat("No coefficients\n")
  }

}


# plot.linreg <- function(x){
#
#   stand_r = sqrt(abs(x$residuals/sqrt(as.vector(x$residvariance))))
#   df = data.frame(x$fitted_values,x$residuals,stand_r)
#   colnames(df) = c('x', 'y','y1')
#   df
# #   p1 = ggplot(data=df,mapping = aes(x,y)) + geom_point()+
# #     geom_smooth( method = 'lm',color="red", se=FALSE)
# #
# #
# #   p2 = ggplot(data=df,mapping = aes(x,y1)) + geom_point() +theme(plot.title = element_text(hjust = 0.5))+
# #     geom_smooth( method = 'lm', se=FALSE)
# #
#   print(p1)
#   print(p2)
# }
#
# plot(xx)
#
# print(xx)


