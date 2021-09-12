source("lasso.R")
source("prior_lasso.R")


supervised_lasso <- function(
  x, y, kfolds, fold_seed, 
  ...)
{
  labeled <- which(!is.na(y))
  lasso_fit(
    x[labeled, ], y[labeled], 
    family = "binomial", tuning = "cv",
    kfolds = kfolds, foldid = fold_seed)
}


supervised_adaptive_lasso <- function(
  x, y, kfolds, fold_seed, 
  ...)
{
  labeled <- which(!is.na(y))
  adaptive_lasso_fit(
    x[labeled, ], y[labeled], 
    family = "binomial", tuning = "cv",
    kfolds = kfolds, foldid = fold_seed)
}


surrogate_only <- function(
  x, y, surrogate_index, 
  ...)
{
  surrogate_x <- x[, surrogate_index]
  
  labeled <- which(!is.na(y))
  model <- glm(
    y[labeled] ~ surrogate_x[labeled], 
    family = binomial())
  gamma <- as.double(coef(model))
  
  beta <- double(ncol(x) + 1L)
  beta[c(1L, surrogate_index + 1L)] <- gamma
  beta
}


perfect_surrogate <- function(
  x, y, surrogate_index, kfolds, fold_seed,
  alpha = NULL,
  ...)
{
  surrogate_x <- x[, surrogate_index]
  other_x <- x[, -surrogate_index]
  
  if (is.null(alpha)) {
    alpha <- adaptive_lasso_fit(
      other_x, surrogate_x, 
      family = "gaussian", tuning = "ic",
      # kfolds = kfolds, foldid = fold_seed,
      df_multiplier = "log(n)")
  } else {
    stopifnot(length(alpha) == ncol(x))
  }
  alpha <- alpha[-1L]  # drop intercept
  augmented_x <- as.double(other_x %*% alpha)

  labeled <- which(!is.na(y))
  
  if (sum(augmented_x^2) == 0){
    model <- glm(
      y[labeled] ~ surrogate_x[labeled], 
      family = binomial())
    gamma <- c(as.double(coef(model)), 0)
  }else{
    model <- glm(
      y[labeled] ~ surrogate_x[labeled] + augmented_x[labeled], 
      family = binomial())
    gamma <- as.double(coef(model))
  }
  beta <- double(ncol(x) + 1L)
  beta[c(1L, surrogate_index + 1L)] <- gamma[c(1L, 2L)]
  beta[-c(1L, surrogate_index + 1L)] <- gamma[3L] * alpha
  beta
}


naive_semisupervised <- function(
  x, y, surrogate_index, kfolds, fold_seed, 
  alpha = NULL,
  ...)
{
  surrogate_x <- x[, surrogate_index]
  other_x <- x[, -surrogate_index]
  
  if (is.null(alpha)) {
    alpha <- adaptive_lasso_fit(
      other_x, surrogate_x, 
      family = "gaussian", tuning = "ic",
      # kfolds = kfolds, foldid = fold_seed,
      df_multiplier = "log(n)")
  } else {
    stopifnot(length(alpha) == ncol(x))
  }
  alpha <- alpha[-1L]  # drop intercept
  augmented_x <- as.double(other_x %*% alpha)
  
  labeled <- which(!is.na(y))
  labeled_z <- cbind(surrogate_x[labeled], augmented_x[labeled], 
                     other_x[labeled, ], deparse.level = 0)
  labeled_y <- y[labeled]
  gamma <- lasso_fit(
    labeled_z, labeled_y,
    family = "binomial", tuning = "cv",
    kfolds = kfolds, foldid = fold_seed,
    penalty_factor = c(0, 0, rep_len(1, ncol(other_x))))

  beta <- double(ncol(x) + 1L)
  beta[c(1L, surrogate_index + 1L)] <- gamma[c(1L, 2L)]
  beta[-c(1L, surrogate_index + 1L)] <- (
    gamma[3L] * alpha + gamma[-c(1L, 2L, 3L)])
  beta
}


get_kappa_sequence <- function(
  x, 
  nkappa = 10L, 
  kappa_multiplier_range = c(1.0, 8.0))
{
  kappa_base <- log(max(nrow(x), ncol(x))) / log(nrow(x))
  kappa_range <- kappa_base * kappa_multiplier_range
  kappa <- exp(seq(log(kappa_range[2]), log(kappa_range[1]),
    length.out = nkappa))
  kappa
}


adaptive_semisupervised <- function(
  x, y, surrogate_index, kfolds, fold_seed, 
  alpha = NULL,
  ...)
{
  surrogate_x <- x[, surrogate_index]
  other_x <- x[, -surrogate_index]
  
  if (is.null(alpha)) {
    alpha <- adaptive_lasso_fit(
      other_x, surrogate_x, 
      family = "gaussian", tuning = "ic",
      # kfolds = kfolds, foldid = fold_seed,
      df_multiplier = "log(n)")
  } else {
    stopifnot(length(alpha) == ncol(x))
  }
  alpha <- alpha[-1L]  # drop intercept
  augmented_x <- as.double(other_x %*% alpha)

  z <- scale(cbind(surrogate_x, augmented_x, 
                   other_x, deparse.level = 0))
  z[is.na(z)] <- 0.0
  labeled <- which(!is.na(y))
  labeled_z <- z[labeled, ]
  labeled_y <- y[labeled]
  foldid <- balanced_cv_fold(
    kfolds, y = labeled_y + 0.1 * labeled_z[, 1] / max(labeled_z[, 1]) +
      0.1 * labeled_z[, 2] / max(labeled_z[, 2]), seed = fold_seed)

  abs_alpha <- abs(alpha)
  active <- c(1L, 2L, 2L + which(abs_alpha >= 1e-8))
  inactive <- 2L + which(abs_alpha < 1e-8)
  factor <- c(0.0001, 0.0001, rep.int(1, ncol(x) - 1L))
  kappa <- get_kappa_sequence(labeled_z)
  lambda <- get_lambda_sequence(labeled_z, labeled_y, foldid)
  
  fit <- lapply(kappa, function(kappa1) {
    factor[inactive] <- kappa1
    cv.glmnet(
      labeled_z, labeled_y, foldid = foldid, 
      family = "binomial", alpha = 1.0,
      standardize = FALSE, intercept = TRUE,
      lambda = lambda * mean(factor),
      penalty.factor = factor / mean(factor))
  })
  tuning <- do.call(
    "rbind", lapply(seq_along(kappa), function(i) 
      data.frame(
        index = i, kappa = kappa[i], 
        lambda = fit[[i]]$lambda, nonzero = fit[[i]]$nzero,
        cvm = fit[[i]]$cvm, cvsd = fit[[i]]$cvsd)))
  good <- min(tuning$cvm + tuning$cvsd * 0.01)
  i <- which(tuning$cvm <= good)[1]
  
  gamma <- as.double(
    coef(fit[[tuning$index[i]]], s = tuning$lambda[i]))
  intercept <- gamma[1L]
  slope <- gamma[-1L]
  slope <- ifelse(
    attr(z, "scaled:scale") > 1e-8, 
    slope / attr(z, "scaled:scale"), 0.0)
  intercept <- intercept - sum(attr(z, "scaled:center") * slope)
  gamma <- c(intercept, slope)
  
  beta <- double(ncol(x) + 1L)
  beta[c(1L, surrogate_index + 1L)] <- gamma[c(1L, 2L)]
  beta[-c(1L, surrogate_index + 1L)] <- (
    gamma[3L] * alpha + gamma[-c(1L, 2L, 3L)])
  attr(beta, "tuning") <- tuning
  
  beta
}


prior_lasso_support <- function(
  x, y, surrogate_index, kfolds, fold_seed, 
  alpha = NULL,
  ...)
{
  if (is.null(alpha)) {
    alpha <- adaptive_lasso_fit(
      x[, -surrogate_index], x[, surrogate_index], 
      family = "gaussian", tuning = "ic",
      # kfolds = kfolds, foldid = fold_seed,
      df_multiplier = "log(n)")
  } else {
    stopifnot(length(alpha) == ncol(x))
  }
  alpha <- alpha[-1L]  # drop intercept

  labeled <- which(!is.na(y))
  beta <- prior_lasso_using_support(
    x = x[labeled, -surrogate_index], 
    s = x[labeled, surrogate_index], 
    y = y[labeled], 
    foldid = balanced_cv_fold(
      kfolds, n = length(labeled), seed = fold_seed), 
    family = "binomial", 
    alpha.hat = alpha)
  beta
}


U_lasso_direction <- function(
  x, y, surrogate_index,
  ...)
{
  labeled <- which(!is.na(y))
  surrogate_x <- x[- labeled, surrogate_index]
  other_x <- x[-labeled, -surrogate_index]
  prev_y <- mean(y[labeled])
  
  cutoff1 <- quantile(surrogate_x, 0.98)
  cutoff2 <- quantile(surrogate_x, 0.02)
  fit_set <- which(surrogate_x > cutoff1 | surrogate_x < cutoff2)
  x_train <- other_x[fit_set, ]
  surrogate_train <- surrogate_x[fit_set]
  Y_star <- ifelse(surrogate_train >= cutoff1, 1, 0)
  
  U_lasso_direction <- lasso_fit(x_train, Y_star, 
                                 family = "binomial", tuning = "cv")
  return(as.vector(U_lasso_direction))
}


U_lasso <- function(x, y, surrogate_index, alpha,...)
{
  alpha <- alpha[-1]
  labeled <- which(!is.na(y))
  y <- y[labeled]
  s <- x[labeled, surrogate_index]
  x_direct <- as.vector(x[labeled, -surrogate_index] %*% alpha)
  data_label <- data.frame(y = y, x_direct = x_direct, s = s)
  mod.fit <- glm(y ~ x_direct, data = data_label, family = binomial())
  fit_coef <- coef(mod.fit)
  
  beta <- double(ncol(x) + 1)
  beta[c(1L, surrogate_index + 1L)] <- c(0, 0)
  beta[-c(1L, surrogate_index + 1L)] <- alpha
  
  beta
}

U_lasso_SS <- function(x, y, surrogate_index, alpha,...){
  alpha <- alpha[-1]
  labeled <- which(!is.na(y))
  y <- y[labeled]
  s <- x[labeled, surrogate_index]
  x_direct <- as.vector(x[labeled, -surrogate_index] %*% alpha)
  data_label <- data.frame(y = y, x_direct = x_direct, s = s)
  mod.fit <- glm(y ~ s + x_direct, data = data_label, family = binomial())
  fit_coef <- coef(mod.fit)
  
  beta <- double(ncol(x) + 1)
  
  beta[c(1L, surrogate_index + 1L)] <- c(fit_coef[1], fit_coef[2])
  beta[-c(1L, surrogate_index + 1L)] <- fit_coef[3] * alpha
  
  beta
}




prior_lasso_value <- function(
  x, y, surrogate_index, kfolds, fold_seed, 
  alpha = NULL,
  ...)
{
  if (is.null(alpha)) {
    alpha <- adaptive_lasso_fit(
      x[, -surrogate_index], x[, surrogate_index], 
      family = "gaussian", tuning = "ic",
      # kfolds = kfolds, foldid = fold_seed,
      df_multiplier = "log(n)")
  } else {
    stopifnot(length(alpha) == ncol(x))
  }
  alpha <- alpha[-1L]  # drop intercept

  labeled <- which(!is.na(y))
  beta <- prior_lasso_using_value(
    x = x[labeled, -surrogate_index], 
    s = x[labeled, surrogate_index], 
    y = y[labeled], 
    foldid = balanced_cv_fold(
      kfolds, n = length(labeled), seed = fold_seed), 
    family = "binomial", 
    alpha.hat = alpha)
  beta
}

