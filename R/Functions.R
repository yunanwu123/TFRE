#' Estimate the tuning parameter for Rank Lasso
#'
#' @description This function estimates the tuning parameter of the Rank Lasso
#' regression given the covariate matrix X.
#' @param X Input matrix, of dimension n_obs x n_vars; each row is an observation vector.
#' @param alpha0 The level to estimate the tuning parameter. Default value is 0.1.
#' See more details in "Details".
#' @param const_lambda The constant to estimate the tuning parameter, should be
#' greater than 1. Default value is 1.01. See more details in "Details".
#' @param times The size of simulated samples to estimate the tuning parameter.
#' Default value is 500.
#' @details In Rank Lasso regressions, the tuning parameter can be estimated independent
#' of errors. In Wang \emph{et al.} (2020), the following tuning parameter is suggested:
#' \deqn{\lambda^* = const\_lambda * G^{-1}_{||\bm{S}_n||_\infty}(1-alpha0)},
#' where \eqn{\bm{S}_n = -2[n(n-1)]^{-1}\sum_{j=1}^n\bm{x}_j[2r_j-(n+1)]}, \eqn{r_1,\ldots,r_n}
#' follows the uniform distribution on the per-mutations of the integers \eqn{\{1,\ldots,n\}},
#' and \eqn{G^{-1}_{||\bm{S}_n||_\infty}(1-alpha0)} denotes the \eqn{(1-alpha0)}-quantile
#' of the distribution of \eqn{||\bm{S}_n||_\infty}.
#' @return The estimated tuning parameter of the Rank Lasso regression given X.
#' @author Yunan Wu and Lan Wang\cr Maintainer:
#' Yunan Wu <yunan.wu@@utdallas.edu>
#' @seealso \code{\link{RankReg}}
#' @references Wang, L., Peng, B., Bradic, J., Li, R. and Wu, Y. (2020),
#'  \emph{A Tuning-free Robust and Efficient Approach to High-dimensional Regression,
#'  Journal of the American Statistical Association, 115:532, 1700-1714},
#'  \doi{10.1080/01621459.2020.1840989}.
#' @examples
#' n <- 100; p <- 400
#' X <- matrix(rnorm(n*p),n)
#' est_lambda(X)
#'
est_lambda<-function(X, alpha0 = 0.1, const_lambda = 1.01, times = 500){
  if(missing(X)|!is.matrix(X)){
    stop("Please supply the covariate matrix X to estimate the tuning parameter for Rank Lasso")
  }
  n <- nrow(X)
  res <- NULL
  for (i in 1:times){
    epsilon_rank <- sample(1:n,n)
    xi <- 2*epsilon_rank-(n+1)
    S <- (-2/n/(n-1))*crossprod(X,xi)
    res <- c(res,max(abs(S)))
  }
  return(quantile(res,1-alpha0)*const_lambda)
}

p_diff <- function(theta,second_stage, lambda, a = 3.7){
  if(second_stage=="scad"){
    less <- theta<=lambda
    y <- pmax(a*lambda-theta,0)
    res <- lambda*less+y*(1-less)/(a-1)
  }else{
    res <- pmin(lambda-abs(theta)/a,lambda)
  }
  return(pmax(res,1e-3))
}

hbic.rank_second <- function(newx, newy, n, beta_int, second_stage, lambda_list, a,
                          thresh, maxin, maxout, const){
  hbic <- NULL
  Beta <- NULL
  C_n <- log(log(n))/n
  logpn <- log(ncol(newx))

  for(i in 1:length(lambda_list)){
    penalty <- p_diff(beta_int,second_stage,lambda_list[i],a)
    x_update <- t(apply(newx,1,function(x)x/penalty))
    obj <- QICD(x_update,newy,lambda_list = 1,thresh = thresh, maxin = maxin, maxout = maxout)
    beta_2nd <- as.vector(obj$beta_final)/penalty
    Beta <- rbind(Beta, beta_2nd)
    HBIC <- log(sum(abs(newy - newx %*% beta_2nd))) + logpn * sum(abs(beta_2nd) > 1e-06) * C_n/const
    hbic <- c(hbic, HBIC)
  }
  return(list(hbic = hbic, Beta = Beta))
}

#' Fit a Rank regression model with Lasso, SCAD or MCP regularization
#'
#' @description This function fits a Rank Lasso model and/or a Rank SCAD or MCP model.
#' The Rank regression models are fitted via QICD algorithm and \emph{Incomplete
#' U-statistics} resampling technique (optional). The tuning parameter of Rank
#' Lasso regression is estimated by the covariate matrix X. The Rank SCAD / MCP
#' regressions are computed at a grid of values for the tuning parameter eta. High
#' dimensional BIC (HBIC) will be used as the criterion on the Rank SCAD / MCP
#' tuning parameter searching.
#' @param X Input matrix, of dimension n_obs x n_vars; each row is an observation vector.
#' @param y Response variable.
#' @param alpha0 The level to estimate the tuning parameter. Default value is 0.1.
#' See more details in the "Details" section of \code{\link{est_lambda}}.
#' @param const_lambda The constant to estimate the tuning parameter, should be
#' greater than 1. Default value is 1.01. See more details in the "Details" section
#' of \code{\link{est_lambda}}.
#' @param times The size of simulated samples to estimate the tuning parameter.
#' Default value is 500.
#' @param incomplete Logical. If \code{incomplete = TRUE}, the \emph{Incomplete
#' U-statistics} resampling technique would be applied; if \code{incomplete = FALSE},
#' the complete U-statistics would be used in computation. See more details in
#' Clémençon, Colin and Bellet (2016).
#' @param const_incomplete The constant for the \emph{Incomplete U-statistics}
#' resampling technique. If \code{incomplete = TRUE}, \code{const_incomplete} x n_obs
#' samples will be randomly selected in the coefficient estimation. Default value is 10.
#' See more details in Clémençon, Colin and Bellet (2016).
#' @param thresh Convergence threshold for QICD algorithm. Default value is 1e-6.
#' See more details in Peng and Wang (2015).
#' @param maxin Maximum number of inner coordiante descent iterations in QICD
#' algorithm; default is 100. See more details in Peng and Wang (2015).
#' @param maxout Maximum number of outter Majoriaztion Minimization step (MM)
#' iterations in QICD algorithm; default is 20. See more details in Peng and Wang (2015).
#' @param second_stage Penalty function for the second stage model. Character vector,
#' which can be "scad", "mcp" and "none". If \code{second_stage = "scad"}, the Rank
#' SCAD regression would be fitted; if \code{second_stage = "mcp"}, the Rank MCP
#' regression would be fitted; if \code{scad = "none"}, only the Rank Lasso regression
#' outputs would be returned.
#' @param a an unknown parameter in SCAD and MCP penalty functions. The default
#' value is 3.7, suggested by Fan and Li (2001).
#' @param eta_list A numerical vector for the tuning parameters to be used in the
#' Rank SCAD or MCP regression. Cannot be \code{NULL} if \code{second_stage = "scad"}
#' or \code{"mcp"}.
#' @param const_hbic The constant to be used in calculating HBIC in the Rank SCAD
#' regression. Default value is 6. See more details in "Details".
#' @details Wang \emph{et al.} (2020) proposed the Rank Lasso estimator for high-dimensional
#' linear regressions with heavy-tailed errors as below:
#' \deqn{\widehat{\bm{\beta}}(\lambda^*) = \arg\min_{\bm{\beta}}\frac{1}{n(n-1)}{\sum\sum}_{i\neq j}\left|(Y_i-\bm{x}_i^T\bm{\beta})-(Y_j-\bm{x}_j^T\bm{\beta})\right| + \lambda^*\sum_{k=1}^p|\beta_k|,}
#' where \eqn{\lambda^*} is the tuning parameter estimated by \code{\link{est_lambda}}.
#' The Rank Lasso model is fitted by QICD algorithm proposed in Peng and Wang (2015).
#' To overcome the computational barrier arising from the U-statistics structure of
#' the aforementioned loss function, we apply the \emph{Incomplete U-statistics}
#' resampling technique which was first proposed in Clémençon, Colin and Bellet (2016).\cr
#' Wang \emph{et al.} (2020) also proposed a second-stage enhancement by using the
#' Rank Lasso estimator \eqn{\widehat{\bm{\beta}}(\lambda^*)} as an initial estimator.
#' It is defined as:
#' \deqn{\widetilde{\bm{\beta}}^{(1)} = \arg\min_{\bm{\beta}}\frac{1}{n(n-1)}{\sum\sum}_{i\neq j}\left|(Y_i-\bm{x}_i^T\bm{\beta})-(Y_j-\bm{x}_j^T\bm{\beta})\right| + \sum_{k=1}^pp_{\eta}'(|\widehat{\beta}_k(\lambda^*)|)|\beta_k|,}
#' where \eqn{p'_{\eta}(\cdot)} denotes the derivative of some nonconvex penalty
#' function \eqn{p_{\eta}(\cdot)}, \eqn{\eta > 0} is a tuning parameter. This
#' function implements the second-stage enhancement with two popular nonconvex
#' penalty functions: SCAD and MCP. The modified high-dimensional BIC criterion
#' in Wang \emph{et al.} (2020) is employed for selecting \eqn{\eta}. Define:
#' \deqn{HBIC(\eta) = \log\left\{{\sum\sum}_{i\neq j}\left|(Y_i-\bm{x}_i^T\widetilde{\bm{\beta}}_{\eta})-(Y_j-\bm{x}_j^T\widetilde{\bm{\beta}}_{\eta})\right|\right\} + |A_\eta|\frac{\log\log n}{n* const\_hbic}\log p,}
#' where \eqn{\widetilde{\bm{\beta}}_{\eta}} denotes the second-stage estimator with
#' the tuning parameter value \eqn{\eta}, and \eqn{|A_\eta|} denotes the cardinality
#' of the index set of the selected model. This function selects the value of \eqn{\eta}
#'  that minimizes HBIC(\eqn{\eta}).
#' @return An object of class "RankReg", which is a list containing at least the
#' following components:
#' \item{X}{The input matrix used.}
#' \item{y}{The response variable used.}
#' \item{incomplete}{Logical. \code{TRUE} if the \emph{Incomplete U-statistics}
#' resampling technique is applied, and \code{FALSE} if not.}
#' \item{beta_RankLasso}{The estimated coefficient vector of the Rank Lasso regression.}
#' \item{rank_lambda}{The estimated tuning parameter of the Rank Lasso regression.}
#' \item{second_stage}{Character vector, \code{"scad"} if the Rank SCAD regression
#' is fitted,  \code{"mcp"} if the Rank MCP regression is fitted, \code{"none"} if
#' only the Rank Lasso regression is fitted.}
#' If \code{second_stage = "scad"}, then the fitted RankReg object will also contain
#' an object named as "Rankscad", which is a list containing the following components:
#' \item{Beta_Rankscad}{The estimated coefficient matrix of the Rank SCAD regression.
#' The diminsion is n_eta x (p+1) with the first column to be the intercepts,
#' where n_eta is the length of \code{eta_list} vector.}
#' \item{df_Rankscad}{The number of nonzero coefficients (intercept excluded) for
#' each value in \code{eta_list}.}
#' \item{eta_list}{The tuning parameter vector used in the Rank SCAD regressions}
#' \item{hbic}{A numerical vector of HBIC values for the Rank SCAD model corresponding
#' to each value in \code{eta_list}.}
#' \item{eta_min}{The eta value which yields the smallest HBIC value in the Rank
#' SCAD regression.}
#' \item{Beta_Rankscad_min}{The estimated coefficient vector which employs \code{eta_min}
#' as the eta value in the Rank SCAD regression.}
#' If \code{second_stage = "mcp"}, then the fitted RankReg object will also contain
#' an object named as "Rankmcp", which is a list containing the following components:
#' \item{Beta_Rankmcp}{The estimated coefficient matrix of the Rank MCP regression.
#' The diminsion is n_eta x (p+1) with the first column to be the intercepts,
#' where n_eta is the length of \code{eta_list} vector.}
#' \item{df_Rankmcp}{The number of nonzero coefficients (intercept excluded) for
#' each value in \code{eta_list}.}
#' \item{eta_list}{The tuning parameter vector used in the Rank MCP regressions}
#' \item{hbic}{A numerical vector of HBIC values for the Rank MCP model corresponding
#' to each value in \code{eta_list}.}
#' \item{eta_min}{The eta value which yields the smallest HBIC value in the Rank
#' MCP regression.}
#' \item{Beta_Rankmcp_min}{The estimated coefficient vector which employs \code{eta_min}
#' as the eta value in the Rank MCP regression.}
#' @author Yunan Wu and Lan Wang\cr Maintainer:
#' Yunan Wu <yunan.wu@@utdallas.edu>
#' @seealso \code{\link{predict.RankReg}}, \code{\link{coef.RankReg}}, \code{\link{plot.RankReg}}, \code{\link{est_lambda}}
#' @references Wang, L., Peng, B., Bradic, J., Li, R. and Wu, Y. (2020),
#' \emph{A Tuning-free Robust and Efficient Approach to High-dimensional Regression,
#' Journal of the American Statistical Association, 115:532, 1700-1714},
#' \doi{10.1080/01621459.2020.1840989}.\cr
#' Peng, B. and Wang, L. (2015),
#' \emph{An Iterative Coordinate Descent Algorithm for High-Dimensional Nonconvex
#' Penalized Quantile Regression, Journal of Computational and Graphical Statistics,
#' 24:3, 676-694}, \doi{10.1080/10618600.2014.913516}. \cr
#' Clémençon, S., Colin, I., and Bellet, A. (2016),
#' \emph{Scaling-up empirical risk minimization: optimization of incomplete u-statistics.
#' The Journal of Machine Learning Research, 17(1):2682–2717}. \cr
#' Fan, J. and Li, R. (2001),
#' \emph{Variable Selection via Nonconcave Penalized Likelihood and its Oracle
#' Properties, Journal of the American Statistical Association, 96:456, 1348-1360},
#' \doi{10.1198/016214501753382273}.
#' @examples
#' n <- 100; p <- 400
#' beta0 <- c(1.5,-1.25,1,-0.75,0.5,rep(0,p-5))
#' eta_list = 0.1*8:20*sqrt(log(p)/n)
#' X <- matrix(rnorm(n*p),n)
#' y <- X %*% beta0 + rt(n,4)
#'
#' Obj_RankLasso <- RankReg(X, y, second_stage = "none")
#' Obj_RankLasso$beta_RankLasso
#'
#' Obj_RankSCAD <- RankReg(X, y, eta_list = eta_list)
#' Obj_RankSCAD$Rankscad$hbic
#' Obj_RankSCAD$Rankscad$df_Rankscad
#' Obj_RankSCAD$Rankscad$Beta_Rankscad_min
#'
#' Obj_RankMCP <- RankReg(X, y, second_stage = "mcp", eta_list = eta_list)
#' Obj_RankMCP$Rankmcp$hbic
#' Obj_RankMCP$Rankmcp$df_Rankmcp
#' Obj_RankMCP$Rankmcp$Beta_Rankmcp_min
#'

RankReg <- function(X, y, alpha0 = 0.1, const_lambda = 1.01, times = 500,
                    incomplete = TRUE, const_incomplete = 10, thresh = 1e-6,
                    maxin = 100, maxout = 20, second_stage = "scad", a = 3.7,
                    eta_list = NULL, const_hbic = 6){
  if(missing(X)|missing(y)){
    stop("Please supply the data (X, y) for the regression")
  }
  if(nrow(X)!=length(y)){
    stop("The row number of X and the length of y should be consistent")
  }
  if(second_stage!="none" & is.null(eta_list)){
    stop("Please supply the tuning parameter list for the Rank SCAD or MCP regression")
  }

  elm_diff <- function(index,v){
    return(v[index[1],]-v[index[2],])
  }
  n <- length(y)
  index_list <- combn(1:n,2)
  x_diff <- t(apply(index_list,2,elm_diff,v=X))
  y_diff <- apply(index_list,2,elm_diff,v=y)
  n_diff <- length(y_diff)
  xbar <- colMeans(X)
  ybar <- mean(y)
  lam_lasso <- est_lambda(X, alpha0, const_lambda, times)

  if(incomplete){
    id <- sample(1:n_diff,const_incomplete*n)
    newx <- x_diff[id,]
    newy <- y_diff[id]
  } else{
    newx <- x_diff
    newy <- y_diff
  }
  obj_RankLasso <- QICD(newx, newy, lambda_list = lam_lasso, thresh = thresh,
                        maxin = maxin, maxout = maxout)
  beta_RL <- as.vector(obj_RankLasso$beta_final)

  intercpt_RL <- ybar- crossprod(beta_RL,xbar)
  beta_RankLasso <- c(intercpt_RL,beta_RL)

  res <- list(X=X,y=y,incomplete=incomplete,beta_RankLasso=beta_RankLasso,
              rank_lambda=lam_lasso,second_stage=second_stage)

  if(second_stage == "scad"){
    obj_Rankscad <- hbic.rank_second(newx, newy, n, beta_RL, second_stage, eta_list, a = a,
                                  thresh = thresh, maxin = maxin, maxout = maxout,
                                  const = const_hbic)
    intercpt_Rs <- ybar - (obj_Rankscad$Beta%*%xbar)
    Beta_Rankscad <- cbind(intercpt_Rs,obj_Rankscad$Beta)
    min_ind <- which.min(obj_Rankscad$hbic)
    df <- colSums(t(obj_Rankscad$Beta)!=0)
    res_scad <- list(Beta_Rankscad = Beta_Rankscad, df_Rankscad = df,
                     eta_list = eta_list, hbic = obj_Rankscad$hbic,
                     eta_min = eta_list[min_ind], Beta_Rankscad_min = Beta_Rankscad[min_ind,])

    res <- append(res, list(Rankscad = res_scad))
  } else if(second_stage == "mcp"){
    obj_Rankmcp <- hbic.rank_second(newx, newy, n, beta_RL, second_stage, eta_list, a = a,
                                  thresh = thresh, maxin = maxin, maxout = maxout,
                                  const = const_hbic)
    intercpt_Rm <- ybar - (obj_Rankmcp$Beta%*%xbar)
    Beta_Rankmcp <- cbind(intercpt_Rm,obj_Rankmcp$Beta)
    min_ind <- which.min(obj_Rankmcp$hbic)
    df <- colSums(t(obj_Rankmcp$Beta)!=0)
    res_mcp <- list(Beta_Rankmcp = Beta_Rankmcp, df_Rankmcp = df,
                     eta_list = eta_list, hbic = obj_Rankmcp$hbic,
                     eta_min = eta_list[min_ind], Beta_Rankmcp_min = Beta_Rankmcp[min_ind,])

    res <- append(res, list(Rankmcp = res_mcp))
  }else if(second_stage != "none"){
    warning("'second_stage' should be one of 'none', 'scad' and 'mcp'")
  }
  attr(res,"class") <- "RankReg"
  return(res)
}


#' Make predictions from a 'RankReg' object
#'
#' @description This function make predictions for new X values from a fitted
#' Rank Lasso or Rank SCAD model.
#' @param object Fitted "RankReg" model object.
#' @param newX Matrix of new values for X at which predictions are to be made.
#' @param s Regression model to use for prediction. Should be one of "1st" and
#' "2nd". See more details in "Details".
#' @param ... Not used. Other arguments to predict.
#' @details If \code{object$second_stage = "none"}, \code{s} cannot be "2nd". If
#' \code{object$second_stage = "none"} and \code{s = "2nd"}, the function will
#' return the predictions based on the Rank Lasso regression. If \code{object$second_stage = "scad"}
#' or \code{"mcp"}, and \code{s = "2nd"}, the function will return the predictions
#' based on the Rank SCAD or MCP regression with the smallest HBIC.
#' @return A vector of predictions for the new X values given the fitted RankReg model.
#' @author Yunan Wu and Lan Wang\cr Maintainer:
#' Yunan Wu <yunan.wu@@utdallas.edu>
#' @seealso \code{\link{RankReg}}, \code{\link{coef.RankReg}}, \code{\link{plot.RankReg}}
#' @references Wang, L., Peng, B., Bradic, J., Li, R. and Wu, Y. (2020),
#' \emph{A Tuning-free Robust and Efficient Approach to High-dimensional Regression,
#' Journal of the American Statistical Association, 115:532, 1700-1714},
#' \doi{10.1080/01621459.2020.1840989}.
#' @examples
#' n <- 100; p <- 400
#' beta0 <- c(1.5,-1.25,1,-0.75,0.5,rep(0,p-5))
#' eta_list = 0.1*8:20*sqrt(log(p)/n)
#' X <- matrix(rnorm(n*p),n)
#' y <- X %*% beta0 + rt(n,4)
#'
#' Obj_RankLasso <- RankReg(X, y, second_stage = "none")
#' Obj_RankSCAD <- RankReg(X, y, eta_list = eta_list)
#' Obj_RankMCP <- RankReg(X, y, second_stage = "mcp", eta_list = eta_list)
#'
#' newX <- matrix(rnorm(n*p),n)
#' predict(Obj_RankLasso, newX, "1st")
#' predict(Obj_RankLasso, newX, "2nd")
#' predict(Obj_RankSCAD, newX, "1st")
#' predict(Obj_RankSCAD, newX, "2nd")
#' predict(Obj_RankMCP, newX, "1st")
#' predict(Obj_RankMCP, newX, "2nd")
#'
#' @method predict RankReg
#'
predict.RankReg<-function(object, newX, s, ...){
  if (missing(newX)) {
    stop("Please supply a value for 'newX'")
  }
  if(attr(object, "class")!="RankReg"){
    stop("Please supply a valid 'RankReg' object")
  }
  if(s=="1st"){
    pred <- cbind(1,newX)%*%object$beta_RankLasso
  }else if(s=="2nd"){
    if(object$second_stage == "scad"){
      pred <- cbind(1,newX)%*%object$Rankscad$Beta_Rankscad_min
    } else if(object$second_stage == "mcp"){
      pred <- cbind(1,newX)%*%object$Rankmcp$Beta_Rankmcp_min
    }else{
      pred <- cbind(1,newX)%*%object$beta_RankLasso
      warning("The object doesn't include a second stage model. Return the predicted values according to the Rank Lasso regression")
    }
  }else{
    stop("s should be one of '1st' and '2nd")
  }
  return(as.numeric(pred))
}



#' Extract coefficients from a 'RankReg' object
#'
#' @description This function extracts the coefficient vector from a fitted Rank
#' Lasso, SCAD or MCP model.
#' @param object Fitted "RankReg" model object.
#' @param s Regression model to use for coefficient extraction. Should be one of
#' "1st" and "2nd". See more details in "Details".
#' @param ... Not used. Other arguments to extract coefficients.
#' @details If \code{object$second_stage = "none"}, \code{s} cannot be "2nd". If
#' \code{object$second_stage = "none"} and \code{s = "2nd"}, the function will
#' return the coefficient vector from the Rank Lasso regression. If \code{object$second_stage = "scad"}
#' or \code{"mcp"}, and \code{s = "2nd"}, the function will return the coefficient
#' vector from the Rank SCAD or MCP regression with the smallest HBIC.
#' @return The coefficient vector from the fitted RankReg model, with the first
#' element as the intercept.
#' @author Yunan Wu and Lan Wang\cr Maintainer:
#' Yunan Wu <yunan.wu@@utdallas.edu>
#' @seealso \code{\link{RankReg}}, \code{\link{predict.RankReg}}, \code{\link{plot.RankReg}}
#' @references Wang, L., Peng, B., Bradic, J., Li, R. and Wu, Y. (2020),
#' \emph{A Tuning-free Robust and Efficient Approach to High-dimensional Regression,
#' Journal of the American Statistical Association, 115:532, 1700-1714},
#' \doi{10.1080/01621459.2020.1840989}.
#' @examples
#' n <- 100; p <- 400
#' beta0 <- c(1.5,-1.25,1,-0.75,0.5,rep(0,p-5))
#' eta_list = 0.1*8:20*sqrt(log(p)/n)
#' X <- matrix(rnorm(n*p),n)
#' y <- X %*% beta0 + rt(n,4)
#'
#' Obj_RankLasso <- RankReg(X, y, second_stage = "none")
#' Obj_RankSCAD <- RankReg(X, y, eta_list = eta_list)
#' Obj_RankMCP <- RankReg(X, y, second_stage = "mcp", eta_list = eta_list)
#'
#' coef(Obj_RankLasso, "1st")
#' coef(Obj_RankLasso, "2nd")
#' coef(Obj_RankSCAD, "1st")
#' coef(Obj_RankSCAD, "2nd")
#' coef(Obj_RankMCP, "1st")
#' coef(Obj_RankMCP, "2nd")
#'
#' @method coef RankReg
#'
coef.RankReg<-function(object, s, ...){
  if(attr(object, "class")!="RankReg"){
    stop("Please supply a valid 'RankReg' object")
  }
  if(s=="1st"){
    coef <- object$beta_RankLasso
  }else if(s=="2nd"){
    if(object$second_stage == "scad"){
      coef <- object$Rankscad$Beta_Rankscad_min
    }else if(object$second_stage == "mcp"){
      coef <- object$Rankmcp$Beta_Rankmcp_min
    }else{
      coef <- object$beta_RankLasso
      warning("The object doesn't include a second stage model. Return the coefficient vector from the Rank Lasso regression")
    }
  }else{
    stop("s should be one of '1st' and '2nd")
  }
  return(as.numeric(coef))
}


#' Plot the second stage model curve for a 'RankReg' object
#'
#' @description This function plots the HBIC curve and the model size curve as a
#' function of the \code{eta} values used, from a fitted Rank SCAD or MCP model.
#' @param x A fitted "RankReg" model object. It should contain a second stage model.
#' @param ... Not used. Other arguments to be passed through plotting functions.
#' @details In the output plot, the red line represents the HBIC curve as a function
#' of \code{eta} values, the blue line represents the number of nonzero coefficients
#' as a function of  \code{eta} values, and the purple vertical dashed line denotes
#' the model selected with the smallest HBIC.\cr
#' This function cannot plot the object if \code{object$second_stage = "none"}.
#' @author Yunan Wu and Lan Wang\cr Maintainer:
#' Yunan Wu <yunan.wu@@utdallas.edu>
#' @seealso \code{\link{RankReg}}, \code{\link{predict.RankReg}}, \code{\link{coef.RankReg}}
#' @references Wang, L., Peng, B., Bradic, J., Li, R. and Wu, Y. (2020),
#' \emph{A Tuning-free Robust and Efficient Approach to High-dimensional Regression,
#' Journal of the American Statistical Association, 115:532, 1700-1714},
#' \doi{10.1080/01621459.2020.1840989}.
#' @examples
#' n <- 100; p <- 400
#' beta0 <- c(1.5,-1.25,1,-0.75,0.5,rep(0,p-5))
#' eta_list = 0.1*8:20*sqrt(log(p)/n)
#' X <- matrix(rnorm(n*p),n)
#' y <- X %*% beta0 + rt(n,4)
#'
#' Obj_RankLasso <- RankReg(X, y, second_stage = "none")
#' Obj_RankSCAD <- RankReg(X, y, eta_list = eta_list)
#' Obj_RankMCP <- RankReg(X, y, second_stage = "mcp", eta_list = eta_list)
#'
#' \dontrun{plot(Obj_RankLasso)}
#' plot(Obj_RankSCAD)
#' plot(Obj_RankMCP)
#'
#' @method plot RankReg
#'
plot.RankReg<-function(x, ...){
  if(attr(x, "class")!="RankReg"){
    stop("Please supply a valid 'RankReg' object")
  }
  if(x$second_stage == "scad"){
    par(mar = c(5, 4, 4, 4) + 0.3)
    plot(x$Rankscad$eta_list, x$Rankscad$hbic,type = "l", col = "red",
         xlab = "eta value", ylab = "HBIC")
    par(new = TRUE)
    plot(x$Rankscad$eta_list, x$Rankscad$df_Rankscad, type = "l",
         col = "blue",  axes = FALSE, xlab = "", ylab = "")
    abline(v = x$Rankscad$eta_min, lty = 2, col = "purple")
    axis(side = 4, at = pretty(range(x$Rankscad$df_Rankscad)))
    mtext("df", side = 4, line = 3)
  }else if(x$second_stage == "mcp"){
    par(mar = c(5, 4, 4, 4) + 0.3)
    plot(x$Rankmcp$eta_list, x$Rankmcp$hbic, type = "l", col = "red",
         xlab = "eta value", ylab = "HBIC")
    par(new = TRUE)
    plot(x$Rankmcp$eta_list, x$Rankmcp$df_Rankmcp, type="l", col = "blue",
         axes = FALSE, xlab = "", ylab = "")
    abline(v = x$Rankmcp$eta_min, lty = 2, col = "purple")
    axis(side = 4, at = pretty(range(x$Rankmcp$df_Rankmcp)))
    mtext("df", side = 4, line = 3)
  }else{
    stop("Please supply a valid 'RankReg' object with a second stage model")
  }
}
