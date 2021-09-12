# library(glmnet)
# library(data.table)


train_glmnet_and_evaluate <- function(
  x_train, y_train, x_test, y_test, lambda, ...)
{
  args <- list(...)
  if ("family" %in% names(args) && args$family == "binomial") {
    y_train <- cbind(1 - y_train, y_train)
  }
  
  model <- glmnet(
    x_train, y_train, lambda = lambda, ...)
  if (is.null(lambda)) {
    lambda <- model$lambda
  }

  xtb_test <- predict(model, x_test, type = "response")
  deviance <- apply(xtb_test, 2, function(v) get_deviance(y_test, v))
  deviance <- c(deviance, rep.int(
    deviance[length(deviance)], length(lambda) - length(deviance)))
  list(model = model,
       evaluation = data.table(lambda = lambda, deviance = deviance))
}


train_cv_glmnet <- function(
  x, y_train, y_eval, foldid, ...)
{
  model <- train_glmnet_and_evaluate(
    x, y_train, x, y_eval, lambda = NULL, ...)$model
  lambda <- model$lambda
  
  nfolds <- max(foldid)
  cv <- rbindlist(lapply(seq_len(nfolds), function(fold) {
    i <- (foldid != fold)
    out <- train_glmnet_and_evaluate(
      x[i, , drop = FALSE], y_train[i],
      x[!i, , drop = FALSE], y_eval[!i],
      lambda = lambda, ...)$evaluation
    out[, fold := fold]
    out
  }))
  
  cv_summary <- cv[, list(
    nfolds = .N,
    deviance_mean = mean(deviance),
    deviance_se = sd(deviance) / sqrt(.N),
    chosen = 0L
    ), keyby = lambda]
  cv_summary[which.min(deviance_mean), chosen := 1L]
  
  lambda_chosen <- cv_summary[chosen == 1L, lambda]
  est <- as.double(coef(model, s = lambda_chosen))
  alpha <- est[1L]
  beta <- est[-1L]
  
  list(model = model, 
       cv = cv, cv_summary = cv_summary,
       estimate = list(alpha = alpha, beta = beta))
}


prior_lasso_core_fit <- function(x, s, y, yp, foldid, family, 
  nlambdas = 50, eta = c(0, 2 ^ (-4 : 4)))
{
  z  <- cbind(s, x)
  fit <- lapply(eta, function(eta1) {
    ytilde <- (y + eta1 * yp) / (1 + eta1)
    fit <- train_cv_glmnet(z, ytilde, y, 
      foldid = foldid, family = family, 
      nlambda = nlambdas, lambda.min.ratio = 0.005,
      thresh = 1e-6)
    fit
  })
  cv <- rbindlist(lapply(seq_along(eta), function(i) 
      data.table(
        index = i, eta = eta[i], fit[[i]]$cv_summary)))
  i <- which.min(cv[, deviance_mean])
  beta <- as.double(coef(
    fit[[cv[i, index]]]$model, s = cv[i, lambda]))

  beta
}


prior_lasso_using_support <- function(x, s, y, foldid, family, alpha.hat)
{
  z <- cbind(s, x)
  j <- c(1, 1 + which(abs(alpha.hat) > 1e-8))
  factor <- rep.int(1.0, ncol(z))
  factor[j] <- 0.0001
  
  lambda <- glmnet(
    z, y, family = family,
    nlambda = 50, lambda.min.ratio = 0.005,
    thresh = 1e-2)$lambda
  lambda <- c(2.0 * max(lambda) / min(factor), lambda)
  fit <- cv.glmnet(
    z, y, foldid = foldid, family = family, 
    lambda = lambda * mean(factor),
    penalty.factor = factor / mean(factor),
    thresh = 1e-6)
  beta <- as.double(coef(fit, s = "lambda.min"))
  yp <- plogis(as.double(beta[1] + z %*% beta[-1]))
  
  prior_lasso_core_fit(x, s, y, yp, foldid, family)
}


prior_lasso_using_value <- function(x, s, y, foldid, family, alpha.hat)
{
  z <- cbind(s, x %*% alpha.hat)
  fit <- glm(y ~ z, family = family)
  simple <- coef(fit)
  yp <- plogis(as.double(simple[1] + z %*% simple[-1]))
  
  prior_lasso_core_fit(x, s, y, yp, foldid, family)
}

