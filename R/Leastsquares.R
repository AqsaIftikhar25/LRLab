#' Linear Regression
#' @name Linear Regression
#' @param formula linear expression
#' @param data data source
#' @description calculate coefficients of linear regression
#' @return return the coefficients and coefficient names
#' @usage linreg (formula,data)
#' @import ggplot2
#' @import gridExtra
#' @importFrom stats median model.matrix pt sd
#' @examples
#' l1 <- linreg(Petal.Length~Species, data = iris)
#' l1 <- linreg(mpg~cyl, data = mtcars)
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
                t_value=t_beta, p_value=pt_beta, formula1=formula, data1=substitute(data),
                formula_Call = match.call())
  # attr(allcoeff,"class") <- "linreg"
  class(allcoeff) <- 'linreg'
  return(allcoeff)
}


# xx = linreg(formula = Petal.Length~Species, data = iris)

#Defining Methods

#print.linreg <- function(x){
#
#  coeff = as.vector(x$Coefficients)
#  names(coeff) = rownames(x$Coefficients)
#  cat("linreg(formula = ",format(x$formula1), ", data = ", x$data1, ")\n\n", sep = "")
#
#  if (length(x$Coefficients)) {
#    cat("Coefficients:\n")
#    print.default(format(coeff), print.gap = 2L,quote = FALSE)
#  }
#  else
#  {
#    cat("No coefficients\n")
#  }
#}

#' @rdname print
#' @method print linreg
#' @export

print.linreg <- function (x){

  cat("call:\n")
  base::print(x$formula_Call)
  cat("\ncoefficients:\n")
  base::print(structure(as.vector(t(x$Coefficients)), names = row.names(x$Coefficients)))

}

#' @rdname plot
#' @method plot linreg
#' @export

plot.linreg <- function(obj){

  #Creating a theme using Linkoping university colors
  theme_linkoping <- function() {
    font <- "Times"
    windowsFonts("Times" = windowsFont("Times"))
    theme(

      panel.background = element_rect(colour = '#00b5e4', fill = '#a0dbed'),
      panel.border = element_rect(colour = '#00b5e4', fill = NA),
      axis.text =    element_text(colour = '#3a3b3b', family = font, size = 12),
      axis.title =   element_text(colour = '#3a3b3b', family = font, size = 12),
      plot.title =   element_text(colour = '#3a3b3b', family = font, face = 'bold', size = 16),)
  }

  stand_resid = sqrt(abs(obj$residuals/sd(obj$residuals)))
  df = data.frame(obj$fitted_values,obj$residuals,stand_resid)
  colnames(df) = c('x', 'y','y1')

  p1 = ggplot(data=df,mapping = aes(x,y)) + geom_point(size = 5, shape = 1)+
    stat_summary(fun = median, color = 'red', geom = 'line', size = 1)+
    labs(y= "Residuals", x = "Fitted values \n lm(Petal.Length ~ Species)", title = "Residuals vs Fitted")+
    theme(plot.title = element_text(hjust = 0.5))+ theme_linkoping()

  p2 = ggplot(data=df,mapping = aes(x,y1)) + geom_point(size = 5, shape = 1)+
    stat_summary(fun = mean, color = 'red', geom = 'line', size = 1)+
    labs(y= expression(sqrt('|Standardized Residuals|')), x = "Fitted values \n lm(Petal.Length ~ Species)", title = "Scale-Location")+
    theme(plot.title = element_text(hjust = 0.5))+ theme_linkoping()

  # print(p1)
  # print(p2)

  gridExtra::grid.arrange(p1,p2,nrow = 2)
}

#' @rdname resid
#' @method resid linreg
#' @export

resid <- function(x){
  UseMethod("resid",x)

}

resid.linreg <- function(x)
{
  return(as.vector(x$residuals))
}

#' @rdname pred
#' @method pred linreg
#' @export

pred <- function(x)
{
  UseMethod("pred",x)
}

pred.linreg <- function(x)
{
  return(as.vector(x$fitted_values))
}

#' @rdname coef
#' @method coef linreg
#' @export

coef.linreg <- function(obj)
{
  cat("Coefficients:",'\n')
  coeff <- as.vector(obj$Coefficients)
  names(coeff) = rownames(obj$Coefficients)
  return(coeff)
}

#' @rdname summary
#' @method summary linreg
#' @export

summary.linreg <- function(obj, formula)
{
  emptyvect <- c("***", "***", "***")
  coef_matrix <- cbind(obj$Coefficients, sqrt(diag(obj$varofregcoeff)), obj$t_value, obj$p_value, emptyvect)
  colnames(coef_matrix) <- c("Coefficients", "Standard Error", "t values", "p values", " ")
  print.table(coef_matrix)

  cat('Residual standard error:', sqrt(obj$residvariance))
  cat(' on ')
  cat(obj$degreeoffreedom,"degrees of freedom\n")


}


# print(xx)
# plot(xx)
# resid(xx)
# pred(xx)
# coef(xx)
# summary(xx)

