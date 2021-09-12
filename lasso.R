library(glmnet)


balanced_cv_fold <- function(kfolds, 
  n = NULL, y = NULL, seed = NULL)
{
  stopifnot(!is.null(n) || !is.null(y))
  stopifnot(kfolds >= 2)

  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (!is.null(y)) {
    n <- length(y)
  } else {
    n <- as.integer(n)
    y <- runif(n)
  }
  kfolds <- as.integer(kfolds)

  quo <- n %/% kfolds  # quotient
  rem <- n %% kfolds   # remainder
  offset <- rep.int(seq_len(quo), 
    c(rep_len(kfolds + 1L, rem), rep_len(kfolds, quo - rem)))
  u <- runif(n, 0, 1) + offset
  
  i <- sort.list(y)     # get sorting indexes
  i <- i[sort.list(u)]  # permute indexes only within similar y values
  fold <- integer(n)
  fold[i] <- rep_len(seq_len(kfolds), n)
  
  fold
}


get_lambda_sequence <- function(
  x, y, foldid = NULL,
  nlambda = 50L, lambda_max_multiplier = 1.5,
  lambda_min_ratio = 1 / length(y))
{
  if (is.null(foldid)) {
    lambda_max <- max(abs(crossprod(x, y - mean(y)) / length(y)))
  } else {
    lambda_max <- max(sapply(unique(foldid), function(f) {
      i <- which(foldid != f)
      x1 <- x[i, ]
      y1 <- y[i]
      max(abs(crossprod(x1, y1 - mean(y1)) / length(y1)))
    }))
  }
  lambda_max <- lambda_max * lambda_max_multiplier
  lambda_min <- pmax(lambda_max * lambda_min_ratio, 1e-6)
  lambda <- exp(seq(log(lambda_max), log(lambda_min), 
                    length.out = nlambda))
  lambda
}


get_estimate <- function(
  z, y, family, l1_fraction, lambda, factor, exclude, 
  tuning, foldid, df_multiplier)
{
  if (tuning == "cv") {
    model <- cv.glmnet(
      z, y, foldid = foldid,
      family = family, alpha = l1_fraction, 
      standardize = FALSE, intercept = TRUE,
      lambda = lambda * mean(factor),
      penalty.factor = factor / mean(factor), exclude = exclude)
    beta <- as.double(coef(model, s = "lambda.min"))
    tuning <- data.frame(
      lambda = model$lambda, nonzero = model$nzero + 1L,
      cvm = model$cvm, cvsd = model$cvsd)
  } else if (tuning == "ic") {
    model <- glmnet(
      z, y, 
      family = family, alpha = l1_fraction, 
      standardize = FALSE, intercept = TRUE,
      lambda = lambda * mean(factor),
      penalty.factor = factor / mean(factor), exclude = exclude)
    nonzero <- sapply(predict(model, type = "nonzero"), length) + 1L
    ic <- deviance(model) + df_multiplier * nonzero
    beta <- as.double(coef(model, s = model$lambda[which.min(ic)]))
    tuning <- data.frame(
      lambda = model$lambda, nonzero = nonzero,
      ic = ic)
  } else {
    stop("Invalid value for 'tuning'")
  }
  list(beta = beta, tuning = tuning)
}


lasso_fit <- function(
  x, y, family = "gaussian", l1_fraction = 1.0, 
  factor = 1.0, exclude = integer(0L), 
  tuning = c("cv", "ic"), kfolds = 10L, foldid = 123456L, 
  df_multiplier = "log(n)",
  ...)
{
  tuning <- match.arg(tuning)
  if (tuning == "cv") {
    if (length(foldid) == 1L) {
      foldid <- balanced_cv_fold(kfolds, y = y, seed = foldid)
    } else if (length(foldid) != nrow(x)) {
      foldid <- rep_len(foldid, nrow(x))
    }
  } else {
    foldid <- NULL
    df_multiplier <- eval(
      parse(text = df_multiplier), list(n = length(y)))
  }

  if (length(factor) != ncol(x)) {
    factor <- rep_len(factor, ncol(x))
  }
  z <- scale(x)
  z[is.na(z)] <- 0.0
  lambda <- get_lambda_sequence(z, y, foldid)

  result <- get_estimate(
    z, y, family, l1_fraction, lambda, factor, exclude, 
    tuning, foldid, df_multiplier)
  beta <- result$beta

  intercept <- beta[1L]
  slope <- beta[-1L]
  slope <- ifelse(
    attr(z, "scaled:scale") > 1e-8, 
    slope / attr(z, "scaled:scale"), 0.0)
  intercept <- intercept - sum(attr(z, "scaled:center") * slope)
  beta <- c(intercept, slope)

  beta <- structure(beta, foldid = foldid, tuning = result$tuning)
  beta
}


adaptive_lasso_fit <- function(
  x, y, family = "gaussian", l1_fraction = 1.0, power = 1.0, 
  factor = 1.0, exclude = integer(0L),
  tuning = c("cv", "ic"), kfolds = 10L, foldid = 123456L, 
  df_multiplier = "log(n)",
  ...)
{
  tuning <- match.arg(tuning)
  if (tuning == "cv") {
    if (length(foldid) == 1L) {
      foldid <- balanced_cv_fold(kfolds, y = y, seed = foldid)
    } else if (length(foldid) != nrow(x)) {
      foldid <- rep_len(foldid, nrow(x))
    }
  } else {
    foldid <- NULL
    df_multiplier <- eval(
      parse(text = df_multiplier), list(n = length(y)))
  }

  if (length(factor) != ncol(x)) {
    factor <- rep_len(factor, ncol(x))
  }
  z <- scale(x)
  z[is.na(z)] <- 0.0
  lambda <- get_lambda_sequence(z, y, foldid)
  
  result <- get_estimate(
    z, y, family, l1_fraction, lambda, factor, exclude, 
    tuning, foldid, df_multiplier)
  beta <- result$beta

  # drop intercept
  if (sum(abs(beta[-1L]) > 1e-6) == 0L) {
    beta[-1L] <- 0.0
    beta <- structure(beta, foldid = foldid, tuning = result$tuning)
    return(beta)
  }

  factor <- 1.0 / abs(beta[-1L]) ** power
  factor <- factor / min(factor)
  threshold <- min(1e6, max(factor) * 0.99)
  exclude <- which(factor > threshold)
  factor[exclude] <- threshold
  
  result <- get_estimate(
    z, y, family, l1_fraction, lambda, factor, exclude, 
    tuning, foldid, df_multiplier)
  beta <- result$beta

  intercept <- beta[1L]
  slope <- beta[-1L]
  slope <- ifelse(
    attr(z, "scaled:scale") > 1e-8, 
    slope / attr(z, "scaled:scale"), 0.0)
  intercept <- intercept - sum(attr(z, "scaled:center") * slope)
  beta <- c(intercept, slope)
  
  beta <- structure(beta, foldid = foldid, tuning = result$tuning)
  return(beta)
}


linear_model_predict <- function(beta, x, probability = FALSE)
{
  y <- drop(x %*% beta[-1L]) + beta[1L]
  if (probability) {
    y <- plogis(y)
  }
  y
}

