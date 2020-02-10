#' Create data according to simulation setting
#'
#' @param lista A list with the setup functions of the DGP
#' @return A list of the data, for X, y0, y1, yobs, z, tau
#' @examples
#' dgp(list("N"=1231,"p"=4, "covariate"="linear", "y_mean"="linear", "ite"="linear")
#' @export

dgp <- function(lista){
  if (!is.list(lista) & is.null(lista)){
    stop("dgp requires a list of setup.")
  }
  if (!any(c("N", "covariate", "p", "y_mean", "ite")  %in% names(lista))){
    stop("A setup condition is missing")
  }
  N <- lista[['N']]
  covariate <- lista[['covariate']]
  p <- lista[['p']]
  y_mean <- lista[['y_mean']]
  ite <- lista[['ite']]
  if (!is.numeric(lista[["N"]]) | !is.numeric(lista[["p"]])){
    stop("N/p should be a number")
  }
  if (!any(c("linear", "linear_normal", "zaidi_lower", "zaidi", "lu") %in% lista[["covariate"]])){
    stop("covariate should be one of the following: linear, linear_normal, zaidi_lower, zaidi, lu")
  }
  if (!any(c("linear", "zero", "zaidi_lower", "zaidi", "lu") %in% lista[["y_mean"]])){
    stop("y_mean should be one of the following: linear, zero, zaidi_lower, zaidi, lu")
  }
  if (!any(c("linear", "zero", "square", "athey", "lu", "sin") %in% lista[["ite"]])){
    stop("y_mean should be one of the following: linear, zero, square, athey, lu, sin")
  }
  #### X ####
  if (covariate == "linear"){
    X1 <- matrix(rnorm(N*(p/2), 0, 1), nrow=N)
    X2 <- matrix(rbinom((N*(p/2)), 1, 0.5), nrow=N)
    X  <- cbind(X1, X2)
  } else if (covariate == "linear_normal"){
    X <- mvtnorm::rmvnorm(N, mean=c(rep(0, p)), sigma = matrix(diag(p), ncol=p))
  } else if (covariate == "zaidi_lower") {
    X1 <- matrix(rnorm(N*3),nrow=N)
    X4 <- purrr::rbernoulli(N, 0.25)
    X5 <- rbinom(N, 2, 0.5)
    X <- cbind(X1, X4, X5)
  } else if (covariate=="zaidi"){
    X <- matrix(rnorm(N*15),nrow=N)
    X <- cbind(X, matrix(runif(N*15),nrow=N))
    X <- cbind(X, pnorm(X[, 1:5] - X[, 16:20]))
    for (k in 36:40){
      lambda_temp <- 5 + 0.75 * X[, (k-35)]* (X[, (k-20)] + X[, (k-5)])
      X <- cbind(X, matrix(rpois(N, lambda_temp),nrow=N))
    }
  } else if (covariate == "lu"){
    X1_10 <- matrix(rnorm(N*10, 0, 1), ncol=10)
    X10_20 <- matrix(rbinom(N*10, 1, 0.5), ncol=10)
    X <- cbind(X1_10, X10_20)
  }
  #### Y mean ####
  if (y_mean == "linear"){
    beta <- rnorm(p, 0, 1)
    y0 <- X%*%matrix(beta)
  } else if (y_mean == "zero"){
    y0 <- rep(0, N) + ifelse(2*X[, 1]*2*X[,2] > -2 & 2*X[, 1]*2*X[,2] < 0, rnorm(N, 3), 0)
  } else if (y_mean == "zaidi_lower"){
    y0 <- -6 + ifelse(X[, 5]==0, 2, ifelse(X[,5]==1, -1, -4)) + abs(X[,3]-1)
  } else if (y_mean == "zaidi"){
    fx_temp <- rowSums(X[,16:19] * exp(X[,30:33]))
    fx <- fx_temp / (1+fx_temp)
    y0 <- 0.15*sum(X[, 1:5]) + 1.5*exp(1+1.5*fx) - 5
  } else if (y_mean == "lu"){
    y0 <- 2.455 - (.4*X[, 1] + .154*X[, 2] - .152*X[, 11] - .126*X[, 12])
  }

  #### ITE ####
  if (ite == "zero"){
    tau <- rep(0, N)
  } else if (ite=="square"){
    tau <- (X[, 1])^2 - 2
  } else if (ite == "linear"){
    tau <- 2*X[, 1] + 2*X[,2]  - 3
  } else if (ite == "athey"){
    tau <- zeta(X[, 1])*zeta(X[,2])*ifelse(X[, 1] < 0, -1, 1)
  } else if (ite=="sin"){
    tau <- 2*sin(2*X[, 1]+2*X[,2])
  } else if (ite=="lu"){
    gx <- .254*X[,2]^2 - .152*X[, 11] - .4*X[,11]^2 - .126*X[, 12]
    tau <- (.4*X[, 1] + .154*X[, 2] - .152*X[, 11] - .126*X[, 12]) - ifelse(gx>0, 1,0)
  }
  #### Z ####
  z <- rbinom(N, 1, 0.5)
  #### Return ####
  colnames(X) <- paste0("X", c(1:ncol(X)))
  return(list(
    "X" = X,
    "z" = z,
    "y0" = y0,
    "y1" = y0 + tau,
    "yobs" = ifelse(z==1, y0+tau, y0),
    "tau" = tau
  ))
}

#' Split the data into train and test set
#'
#' @param data data to split to experiment-roll-out
#' @param n1 number of units in experimentation phase
#' @export
split_data <- function(data, n1){
  X <-   data[["X"]]
  y <-   data[["yobs"]]
  z <-   data[["z"]]
  y0 <-  data[["y0"]]
  y1 <-  data[["y1"]]
  tau <- data[["tau"]]
  N <- nrow(data["X"])
  ix <- sample(1:nrow(X), n1)
  X.train <- as.matrix(X[ix,], ncol=ncol(X))
  z.train <- z[ix]
  y.train <- y[ix]
  y0.train <- y0[ix]
  y1.train <- y1[ix]
  tau.train <- tau[ix]
  X.test <- as.matrix(X[-ix,], ncol=ncol(X))
  z.test <- z[-ix]
  y.test <- y[-ix]
  y0.test <- y0[-ix]
  y1.test <- y1[-ix]
  tau.test <- tau[-ix]
  return(list(
    "experimentation" = list("X" = X.train,
                   "yobs" = y.train,
                   "z" = z.train,
                   "y0"=y0.train,
                   "y1"=y1.train,
                   "tau"=tau.train),
    "rollout"  = list("X" = X.test,
                   "yobs" = y.test,
                   "z" = z.test,
                   "y0"=y0.test,
                   "y1"=y1.test,
                   "tau"=tau.test)))
}


#' Create n_1+n_2 dataset
#'
#' @param data data to add together
#' @param indexes indexes to move from rollout to sampling
#' @export

sampling_data <- function(data, indexes){
  data[["sampling"]][["X"]] <- rbind(data[["experimentation"]][["X"]],
                                     data[["rollout"]][["X"]][indexes,])
  data[["sampling"]][["yobs"]] <- c(data[["experimentation"]][["yobs"]],
                                    data[["rollout"]][["yobs"]][indexes])
  data[["sampling"]][["z"]] <- c(data[["experimentation"]][["z"]],
                                 data[["rollout"]][["z"]][indexes])
  data[["sampling"]][["y1"]] <- c(data[["experimentation"]][["y1"]],
                                  data[["rollout"]][["y1"]][indexes])
  data[["sampling"]][["y0"]] <- c(data[["experimentation"]][["y0"]],
                                  data[["rollout"]][["y0"]][indexes])
  data[["sampling"]][["tau"]] <- c(data[["experimentation"]][["tau"]],
                                   data[["rollout"]][["tau"]][indexes])
  data[["rollout"]][["X"]] <- data[["rollout"]][["X"]][-indexes,]
  data[["rollout"]][["yobs"]] <- data[["rollout"]][["yobs"]][-indexes]
  data[["rollout"]][["z"]] <- data[["rollout"]][["z"]][-indexes]
  data[["rollout"]][["tau"]] <- data[["rollout"]][["tau"]][-indexes]
  data[["rollout"]][["y1"]] <- data[["rollout"]][["y1"]][-indexes]
  data[["rollout"]][["y0"]] <- data[["rollout"]][["y0"]][-indexes]
  return(data)
}


.onLoad <- function(libname, pkgname) {
  # use superassignment to update global reference to
  # sgd_py <<- reticulate::import("sklearn.linear_model", delay_load = TRUE)
  # print(system("which python3"))
  # use_python(Sys.which("python3"))
  np <<- reticulate::import("numpy", delay_load = TRUE)
  sk <<- reticulate::import("sklearn", delay_load = TRUE)
}
