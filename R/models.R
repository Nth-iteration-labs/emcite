#' Train BART model for ITE estimation
#'
#' @param data Data to train on, previously splitted into train/test set
#' @return outcome Outcome psBART
#' @return ps propensity score BART model
#' @export
train_model <- function(data){
    x <- data[["X"]]
    z <- data[['z']]
    y <- data[['yobs']]
  # Continuous or binary outcome
  if (length(unique(y)) == 2){
    outcome_linear <- F
  } else {outcome_linear <- T}
    # Train prpensity-score model
    ps <- pbart(x.train = x,
                y.train = z,
                rm.const=F)
    est_ps <- ps$prob.train.mean
    if (length(unique(est_ps))==1){
      est_ps <- est_ps + rnorm(length(est_ps), 0, 0.001)
    }
    # possibility to add external propensity score
    if (!is.null(data$prop_scores)){
      est_ps <- data$prop_scores
    }
    x <- cbind(x, est_ps, z)
    # Train normal model
    if (outcome_linear){
      outcome <- wbart(x, y, rm.const=F)
    } else {
      outcome <- mc.pbart(x.train = x,
                       y.train = y,
      rm.const=F)
    }
  # Return both models
  list("outcome"= outcome,
       "ps"     = ps)
}

#' Predict ITE from trained BART model on given data
#'
#' @param data data to predict on
#' @param models models that were trained previously
#' @return predictions of y1, y0, tau and propensity score
#' @export
predict_ite <- function(data, models){
  x.test <- data.frame(data[["X"]])
  outcome <- models[['outcome']]
  ps      <- models[['ps']]
  # Add predicted PS-score:
  x.test$ps  <- predict(ps, x.test)$prob.test.mean
  x.test.control <- copy(x.test)
  x.test.treat <- copy(x.test)

  x.test.control$z <- 0
  x.test.treat$z   <- 1
  y_hat_treat <- predict(outcome, x.test.treat)
  y_hat_control <- predict(outcome, x.test.control)
  # Check binary or continuous outcome
  if (class(outcome) == "pbart"){
    outcome_lin <- "logit"
  } else {outcome_lin <- "linear"}
  preds_tau <- switch(outcome_lin,
                  "logit" = y_hat_treat$prob.test - y_hat_control$prob.test,
                  "linear" = y_hat_treat - y_hat_control)
  if (outcome_lin == "linear"){
    preds <- list("tau" = preds_tau,
                  "train.ps"  = ps$prob.train.mean,
                  "y1"  = y_hat_treat,
                  "y0"  = y_hat_control)
  } else if (outcome_lin=="logit") {
    preds <- list("tau" = preds_tau,
                  "train.ps"  = ps$prob.train.mean,
                  "y1"  = y_hat_treat$prob.test,
                  "y0"  = y_hat_control$prob.test)
  }
  return(preds)
}

#' Fit stochastic gradient descent
#'
#' @param X the feature matrix
#' @param tau_prediction ITE predictions
#' @param weight_types value for $\alpha$
#'
#' @return vector of weights
#' @export
fit_gradient_descent <- function(X, tau_predictions, weight_types, parallel=F, R=F) {
  # Calculate predicted Type-S error
  prob_of_mistake <- apply(tau_predictions,
                           2,
                           function(x)
                             weight_types * (ifelse(mean(x) > 0, sum(x < 0), sum(x > 0)) / length(x))) + 1
  tau_mean <- apply(tau_predictions, 2, mean)
  if (R){
    m.glm <- cv.glmnet(X, tau_mean, weights = prob_of_mistake, alpha=0)
    thetas <- coef(m.glm, s = "lambda.min")
    lambda.m1 <- m.glm$lambda.min
  }
  if (!R){
    df.tau.model <- cbind(1, X)
  # if (!parallel){
    thetas <- sgd_py(X = df.tau.model, y = tau_mean, w = prob_of_mistake, theta=NULL)
    lambda.m1 <- NA
  }
  # } else {
  #   thetas <- sgd_py2(X = df.tau.model, y = tau_mean, w = prob_of_mistake, theta=NULL)
  # }
  return(list(thetas=thetas, lambda=lambda.m1))
}

#' Run model in python
#'
#' @param X feature matrix
#' @param y outcomes
#' @param w weight vector
#' @param theta initial parameter vector
#' @export
sgd_py <- function(X, y, w, theta=NULL){
  model = sk$linear_model$SGDRegressor()
  if (is.null(theta)){
    model$fit(X, y, sample_weight=w)
  } else {
    original_weight = np$copy(theta)
    model$fit(X, y, sample_weight=w, coef_init=original_weight)
  }
  return(model$coef_)
}

sgd_py_p <- function(lista){
  X <- lista[["X"]]
  y <- lista[["y"]]
  w <- lista[["w"]]
  l1 <- lista[["lambda"]]
  theta = lista[["theta"]]
  # model = sk$linear_model$SGDRegressor()
  if (is.null(theta)){
    model$fit(X, y, sample_weight=w)
  } else {
    # original_weight = np$copy(theta)
    # model$fit(X, y, sample_weight=w, coef_init=original_weight)
    m.glm <- cv.glmnet(X, y, weights = w)
    thetas <- coef(m.glm, s = m.glm$lambda.min)
  }
  return(thetas)
}

#' Acquisition function, selects $n_2$ units from the test set
#'
#' @param data list containing train and test covariates
#' @param tau_train estimated ITE for train set
#' @param tau_test  predicted ITE for test set
#' @param theta initialization of theta
#' @param n2 number of samples to select
#' @param type which method to run: random/variance/type-s/emcite
#' @param weight_types weight of type-s error in EMCITE
#' @param B number of bootstrap predictions in EMCITE
#'
#' @return vector of selected indexes
#' @export
af <- function(data, tau_train, tau_test, theta, n2,
               type="emcite",
               weight_types=5,
               B=2, parallel=F, lambda=NULL, R=F){
  print(paste0("start of ", type))
  types <- apply(tau_test,
                 2,
                 function(x) (ifelse(mean(x) > 0,
                                     sum(x<0), sum(x>0))/length(x)))
  if (type=="random"){
    print(sample(1:ncol(tau_test), n2))
    sel.vec <- sample(1:ncol(tau_test), n2)
  } else if (type=="variance"){
    v_tau <- apply(tau_test, 2, var)
    sel.vec <- order(v_tau, decreasing = T)[1:n2]
  } else if (type=="type-s"){
    # types <- apply(tau_test,
    #       2,
    #       function(x) (ifelse(mean(x) > 0,
    #                                        sum(x<0), sum(x>0))/length(x)))
    sel.vec <- order(types, decreasing = T)[1:n2]
  } else if (type=="emcite") {
    if (!R){
      theta <- theta$thetas
    }
    X.test <- data$rollout$X
    X.train <- data$experimentation$X
    tau_train_mean <- apply(tau_train, 2, mean)
    tau_test_mean <- apply(tau_test, 2, mean)
    prob_of_mistake_train <- apply(tau_train,
                                     2,
                                     function(x) weight_types*(ifelse(mean(x) > 0,
                                                                      sum(x<0), sum(x>0))/length(x))) + 1
    prob_of_mistake <- apply(tau_test,
                               2,
                               function(x) weight_types*(ifelse(mean(x) > 0,
                                                                sum(x<0), sum(x>0))/length(x))) + 1
    #X.model.train <- data.frame(cbind(1, X.train))
    #X.model <- data.frame(cbind(1, X.test))
    #names(X.model.train) <- names(X.model) <- paste0("X", c(1:ncol(X.model.train)))
    sel.vec <- c()
    id <- 1:ncol(tau_test)
    # Create empty matrix for all the draws <- draw B row of simulation to get a bootstrap estimate
    if (B>1){
      draws <- matrix(nrow=B, ncol=ncol(tau_test))
      # Sample randomly 5
      for (dr in 1:ncol(tau_test)){
        draws[, dr] <- sample(tau_test[, dr], B)
      }
    } else {
        draws <- c(rep(NA,ncol(tau_test)))
      for (dr in 1:ncol(tau_test)){
        draws[dr] <- sample(tau_test[, dr], 1)
      }
    }
    # EMCITE, ALgorithm 1
    for (i in 1:n2){
      change <- list()
      # Gradients to update every time <- the one with the highest norm will be selected
      gradients <- list()
      # if (ncol(tau_test) > 1000){
      if (R){
        run.numbers <- sample(id, 100, prob = types)
      } else{run.numbers <- id}
      #} else {
      #  run.numbers <- id
      #}
      # print(id[(length(id)-10):length(id)])
      # print(run.numbers)
      if (!parallel){
        for (ix in seq_along(run.numbers)){
          ix.no <- run.numbers[ix]
          # Initialize vector to store the changes, mean of this will be the measurement
          mean_change <- c()
            for (simulation in c(1:B)){
              sim.value <- ifelse(B>1, draws[simulation, ix.no],  draws[ix.no])
              # In Python, with automated function

              new_theta <- sgd_py(X=as.matrix(cbind(1, rbind(X.train, X.test[ix.no, ]))),
                                       y=c(tau_train_mean, sim.value),
                                       w=c(prob_of_mistake_train, prob_of_mistake[ix.no]),
                                       theta = theta)
              # print(new_theta)
              grad <- theta - new_theta
              mean_change <- c(mean_change, sum(abs(grad)))
              # Calculate mean of gradient to update weights
              if (length(gradients) < ix){
                gradients[[ix]] <- new_theta
              } else {
                gradients[[ix]] <- gradients[[ix]] + (new_theta-gradients[[ix]])/(simulation)
              }
          }
          change[[ix]] <- mean(mean_change)}
          selected <- which(unlist(change)==max(unlist(change)))
          if (length(selected) > 1){
            selected <- sample(selected, 1)
          }
          tau_train_mean <- c(tau_train_mean, ifelse(B>1, mean(draws[, run.numbers[selected]]),
                                                     draws[run.numbers[selected]]))
          theta <- gradients[[selected]]
          sel.vec <- c(sel.vec, run.numbers[selected])
          X.train <- rbind(X.train,X.test[run.numbers[selected], ])
          prob_of_mistake_train <- c(prob_of_mistake_train, prob_of_mistake[run.numbers[selected]])
          # X.test <- X.test[-run.numbers[selected], ]
          # tau_test <- tau_test[, -run.numbers[selected]]
          # if (B>1){
          #   draws <- draws[, -run.numbers[selected]]
          # } else {draws <- draws[-run.numbers[selected]]}
          id <- id[-run.numbers[selected]]
          } else {
        lista_to_p <- list()
        list_res   <- list()
          for (ix in seq_along(run.numbers)){
            ix.no <- run.numbers[ix]
            # print(as.matrix(rbind(X.train,
            #                      X.test[ix.no, ])))
            for (simulation in 1:B){
              lista_to_p[[paste0("ix_", ix.no, "_sim_",
                                 simulation)]] <- list(X=as.matrix(rbind(X.train,
                                                              X.test[ix.no, ])),
                                                      y=c(tau_train_mean, ifelse(B>1,
                                                          draws[simulation, ix.no],
                                                          draws[ix.no])),
                                                      w=c(prob_of_mistake_train,
                                                          prob_of_mistake[ix.no]),
                                                      theta = theta,
                                                      lambda=lambda)
            }
          }
        # new_theta_l <- mclapply(lista_to_p, sgd_py_p, mc.cores=(parallel::detectCores()-1))
        new_theta_l <- lapply(lista_to_p, sgd_py_p)
        # test <- mclapply(new_theta_l, function(x) data.table(values=sum(abs(theta - x))),
        #                  mc.cores = (parallel::detectCores()-1))
        test <- lapply(new_theta_l, function(x) data.table(values=sum(abs(theta - x))))
        # print(test)

        td <- data.table(names=names(test), rbindlist(test))
        td[, c("ix", "sim") := tstrsplit(names, "_")[c(2,4)]]
        selected <- as.integer(td[, .(mv = mean(values)), by=ix][order(-mv)][1, ix])
        print(selected)
        tau_train_mean <- c(tau_train_mean, ifelse(B>1, mean(draws[, selected]),
                                                   draws[selected]))
        # print(grep(paste0("ix_",selected, "_sim_"), names(new_theta_l)))
        # print(new_theta_l[grepl(paste0("ix_",selected, "_sim_"), names(new_theta_l))])
        theta <- colMeans(rbindlist(lapply(new_theta_l[grepl(paste0("ix_",selected, "_sim_"), names(new_theta_l))],
               function(x) data.table(t(as.matrix(x))))))
        sel.vec <- c(sel.vec, selected)
        X.train <- rbind(X.train,X.test[selected, ])
        prob_of_mistake_train <- c(prob_of_mistake_train, prob_of_mistake[selected])
        id.sel <- which(id==selected)
        print(id.sel)
        id <- id[-id.sel]
        types <- types[-id.sel]
        }
      # prob_of_mistake <- prob_of_mistake[-run.numbers[selected]]
    }
  }
  return(sel.vec)
}

#' Retrain and evaluate the model based on PEHE
#'
#' @param data data that contains the sampling phase
#' @export
retrain_and_metrics <- function(data){
  m.new <- train_model(data[["sampling"]])
  tau.preds <- predict_ite(data[["rollout"]], m.new)
  tau.train.preds <- predict_ite(data[["sampling"]], m.new)
  if (any(!is.na(data[["rollout"]][["y1"]]))){
    optimal.y.experimentation <- ifelse(data[["experimentation"]][['tau']] > 0,
                                        data[["experimentation"]][['y1']],
                                        data[["experimentation"]][['y0']])
    optimal.y.train <- ifelse(data[["sampling"]][['tau']] > 0,
                              data[["sampling"]][['y1']],
                              data[["sampling"]][['y0']])
    optimal.y.test  <- ifelse(data[["rollout"]][['tau']] > 0,
                              data[["rollout"]][['y1']],
                              data[["rollout"]][['y0']])
    regret.train <- sum(optimal.y.train - data[["sampling"]][['yobs']]) - sum(optimal.y.experimentation-
                                                                                data[["experimentation"]][["yobs"]])
    regret.test  <- sum(optimal.y.test - ifelse(apply(tau.preds$tau,2,mean)>0,
                                                data[["rollout"]][['y1']],
                                                data[["rollout"]][['y0']]))
    pehe.train   <-  mean((data[["sampling"]][['tau']] - apply(tau.train.preds$tau,2,mean))^2)
    pehe.test    <- mean((data[["rollout"]][['tau']] - apply(tau.preds$tau,2,mean))^2)
    return(data.table(pehe_sampling   = pehe.train,
                      pehe_rollout    = pehe.test,
                      regret_sampling = regret.train,
                      regret_rollout  = regret.test))
  } else {
    predictions <- apply(tau.preds$tau, 2 , mean)
    perf <- performance(apply(tau.preds$y1, 2, mean),
                        apply(tau.preds$y0, 2, mean),
                        data[["rollout"]][["yobs"]],
                        data[["rollout"]][["z"]],
                        direction=1
    )
    # Risk of policy with rejection sampling
   #  if (data[["active"]]){
   #    risk_policy <- sapply(seq(0, 1, .1), function(x) {
   #      1 - (sum(data[["rollout"]][["yobs"]][data[["rollout"]][["z"]]==(
   #        predictions > quantile(predictions, 1-x))])*(1/data[["rollout"]][["prop_scores"]]))/sum(
   #          data[["rollout"]][["z"]]==(predictions > quantile(predictions, 1-x)))
   #      })
   #    } else {
        risk_policy <- sapply(seq(0, 1, .1), function(x) {
          1 - (sum(data[["rollout"]][["yobs"]][data[["rollout"]][["z"]]==(
          predictions > quantile(predictions, 1-x))])*(1/0.5))/sum(data[["rollout"]][["z"]]==(predictions > quantile(predictions, 1-x)))
    })
    return(list("perf"=perf,
                "risk" = risk_policy))
  }
}
#'
#'
#'
