library(compiler)
enableJIT(3L)
setCompilerOptions(optimize = 3L, suppressAll = TRUE)
library(data.table)
library(glmnet)
library(Matrix)

source("roc_auc.R")
source("methods.R")


## data generation ----

simulate_data <- function(scenario, nl = 100, nt = 10000, p = 500, seed = 1L)
{
  # nl <- 100    # number of observations with labels
  # nt <- 10000  # number of observations in total
  # p <- 100     # dimension of x
  # nl = 100; nt = 50000; p = 50; seed = 1; scenario = 1

  # set alpha and beta
  v1 <- c(0.5, 1, -0.8, 0.6, 0.2)
  v2 <- c(0.1, -0.2, -0.2, 0.2, 0.7)
  w1 <- c(-0.05, -0.5, 1.4, 0.5, -0.6)
  w2 <- c(0.02, 0.05, 0.02, -0.02, -0.05)
  if (scenario == 1) {
    # exactly multiple
    alpha <- c(v1, v2, double(p - 10))
    beta <- c(v1, v2, double(p - 10)) * 1.5
  } else if (scenario == 2) {
    # some deviation, more sparse delta
    alpha <- c(v1, v2, double(p - 10))
    beta <- c(v1 + w1, v2 + w2, double(p - 10)) * 1.5
  } else if (scenario == 3) {
    # alpha is more dense
    alpha <- c(v1, v2, v2, v2, double(p - 20))
    beta <- c(v1 + w1, v2 + w2, double(p - 10)) * 1.5
  } else if (scenario == 4) {
    # beta is more dense
    alpha <- c(v1, double(p - 5))
    beta <- c(v1 + w1, v2 + w2, double(p - 10)) * 1.5
  } else if (scenario == 5) {
    # same support but different magnitude
    alpha <- c(v1, v2, double(p - 10))
    beta <- c(v2, v1, double(p - 10)) * 1.5
  } else if (scenario == 6) {
    # different support and different magnitude
    alpha <- c(v1, v2, double(p - 10))
    beta <- c(v2, double(5), v1, double(p - 15)) * 1.5
  } else {
    alpha <- c(v1, v2, double(p - 10))
    beta <- alpha * 0
  }

  # generate x, s, and y
  set.seed(seed)
  rho <- 0.5
  sigma <- outer(seq_len(p), seq_len(p), function(a, b) rho ^ abs(a - b))
  sigma_cholesky <- chol(sigma)
  x <- t(crossprod(sigma_cholesky, matrix(rnorm(p * nt), p, nt)))
  x <- log1p(round(exp(2 * x)))
  
  s <- 1 + as.double(x %*% alpha) + 2 * rnorm(nt)
  s <- log1p(round(exp(s)))
  
  eta <- -4 + 0.5 * s + as.double(x %*% beta)
  y <- rbinom(nt, 1, plogis(eta))
  
  # organize x, s and y into list
  data <- list(y = y, x = I(cbind(s, x, deparse.level = 0)))
  if (nt > nl) {
    data$y[seq.int(nl + 1L, nt)] <- NA_real_
  }
  data$surrogate_index <- 1L

  # save parameter information
  k <- beta / alpha
  k <- k[is.finite(k)]
  rho <- k[
    which.min(sapply(k, function(rho) 
      sum(abs(beta - rho * alpha)) +
        1e-4 * sum(abs(beta - rho * alpha) > 1e-8)))]
  info <- list(
    alpha_true = c(1, alpha),
    beta_true = c(-4, 0.5, beta),
    rho_true = rho,
    relative_l1_distance = sum(abs(beta - rho * alpha)),
    relative_sparsity = sum(abs(beta - rho * alpha) > 1e-8),
    auc_eta = get_auc(y, eta), 
    auc_s = get_auc(y, s))
  
  list(data = data, info = info)
}


simulate_data_misspecified  <- function(scenario, nl = 100, nt = 10000,
                                        p = 500, seed = 1L, beta.fit = T)
{
  beta_x <- c(c(0.8, 1, -1, 0.8, 0.4), rep(0, p - 5))
  delta_1 <- c(c(c(0.02, 0.04, 0.04, -0.02, -0.02), rep(0, 15)), rep(0, p - 20))
  delta_2 <- c(c(10 * c(0.02, 0.04, 0.04, -0.02, -0.02), rep(0, 15)), rep(0, p - 20))
  delta_3 <- c(10 * rep(c(0.02, 0.04, 0.04, -0.02, -0.02), 4), rep(0, p - 20))
  
  # Generate binary x
  set.seed(seed)
  rho <- 0.5
  sigma <- outer(seq_len(20), seq_len(20), function(a, b) rho ^ abs(a - b))
  sigma_20 <- chol(sigma)
  sigma <- outer(seq_len(p - 20), seq_len(p - 20), function(a, b) rho ^ abs(a - b))
  sigma_rem <- chol(sigma)
  
  sigma_full <- bdiag(sigma_20,  sigma_rem) 
  x_latent <- t(crossprod(sigma_full, matrix(rnorm(p * nt), p, nt)))
  x <- matrix(pnorm(as.vector(x_latent), 0, 1), nrow = nt, ncol = p)  
  x <- 2 * x - 1
  eta <- as.double(x %*% beta_x)
  y <- I(eta + rnorm(nt, 0, 1) > 0)
  
  if (scenario == 1){
    alpha <- rep(0, p)
    gamma <- rep(0, p)
    mu <- 1
  }
  
  if (scenario == 2){
    alpha <- c(c(0.6, -0.4, 0.4, 0.5, -0.5), rep(0, p - 5))
    gamma <- c(c(0.3, 0.4, 0.6, -0.5, - 0.5), rep(0, p - 5))
    mu <- 1.5
  }
  
  if (scenario == 3){
    alpha <- c(rep(c(0.6, -0.4, 0.4, 0.5, -0.5), 3), rep(0, p - 15))
    gamma <- c(rep(c(0.3, 0.4, 0.6, -0.5, - 0.5), 3), rep(0, p - 15))
    mu <- 2
  }

  mean_s <- x %*% alpha + x %*% gamma * y + mu * y
  s <- mean_s + rnorm(nt, 0, 1)
  
  data <- list(y = y, x = I(cbind(s, x, deparse.level = 0)))
  if (nt > nl) {
    data$y[seq.int(nl + 1L, nt)] <- NA_real_
  }
  data$surrogate_index <- 1L
  
  # save parameter information
  
  info <- NULL
  
  if (beta.fit == T){
    # Obtain the true beta
    data.all <- list(y = y, x = I(cbind(s, x[,1:20], deparse.level = 0)))
    
    oracle.fit <- glm(y~., data = data.all, family = binomial())
    beta <- as.vector(coef(oracle.fit))
    oracle.fit.alpha <- lm(s ~ x[,1:20])
    alpha <- as.vector(coef(oracle.fit.alpha))
    
    beta.true <- as.vector(c(beta, rep(0, p - 20)))
    alpha.true <- as.vector(c(alpha, rep(0, p - 20)))
    
    k <- exp(seq(log(0.01), log(100), length.out = 10000))
    rho <- k[
      which.min(sapply(k, function(rho) 
        sum(abs(beta[-c(1,2)] - rho * alpha[-c(1)])) +
          1e-4 * sum(abs(beta[-c(1,2)] - rho * alpha[-c(1)]) > 1e-8)))]
    
    info <- list(
      alpha_true = alpha.true, beta_true = beta.true, rho_true = rho,
      relative_l1_distance = sum(abs(beta[-c(1,2)] - rho * alpha[-c(1)])),
      relative_sparsity = sum(abs(beta[-c(1,2)] - rho * alpha[-c(1)]) > 1e-8),
      cor.coef = cor(beta[-c(1,2)], alpha[-c(1)]))
    
  }
  list(data = data, info = info)
  
}




## evaluations ----

assess_quality <- function(data, beta_true, beta_hat)
{
  link_true <- as.double(data$x %*% beta_true[-1L]) + beta_true[1L]
  response_true <- plogis(link_true)

  link_hat <- as.double(data$x %*% beta_hat[-1L]) + beta_hat[1L]
  response_hat <- plogis(link_hat)
  
  c(auc = get_auc(data$y, response_hat),
    excess_risk = 0.5 * (get_deviance(data$y, response_hat) - 
                    get_deviance(data$y, response_true)),
    mse_link = mean((link_hat - link_true) ** 2),
    mse_response = mean((response_hat - response_true) ** 2),
    bias_response = mean(response_hat - response_true),
    var_response = var(response_hat - response_true),
    var_ratio_response = var(response_hat) / var(response_true),
    l1_error = sum(abs(beta_hat - beta_true)),
    l2_error = sqrt(sum(beta_hat - beta_true) ** 2),
    l0_norm = sum(abs(beta_hat) > 1e-8),
    l1_norm = sum(abs(beta_hat)))
}


## simulations ----


run_simulation <- function(
  method = c("supervised_lasso", "supervised_adaptive_lasso", 
             "surrogate_only", "perfect_surrogate", 
             "naive_semisupervised", "adaptive_semisupervised",
             "prior_lasso_support", "prior_lasso_value"),
  seeds = seq_len(10), scenario = 1, 
  nl = 100, nt = 10000, p = 100)
{
  # method = "supervised_lasso"; seeds = 1 : 5; i = 1
  # scenario = 2; nl = 100; nt = 10000; p = 100
  method <- match.arg(method)
  compute_estimate <- get(method)
  
  quality <- vector("list", length(seeds))
  for (i in seq_along(seeds)) {
    data <- simulate_data(scenario = scenario, 
                          nl = nl, nt = nt, p = p, 
                          seed = seeds[i])$data
    
    beta_hat <- compute_estimate(
      x = data$x, y = data$y, 
      surrogate_index = data$surrogate_index,
      kfolds = 10L, fold_seed = seeds[i] + 10000L)
    
    test <- simulate_data(scenario = scenario,
                          nl = 10000L, nt = 10000L, p = p, 
                          seed = seeds[i] + 100000L)
    quality[[i]] <- assess_quality(
      test$data, test$info$beta_true, beta_hat)
    
    print(i)
  }
  
  quality <- do.call("rbind", quality)
  quality <- data.table(
    scenario = scenario, method = method, 
    nl = nl, nt = nt, p = p,
    quality)
  quality
}





run_simulation_mis <- function(
  method = c("supervised_lasso", "supervised_adaptive_lasso", 
             "surrogate_only", "perfect_surrogate", 
             "naive_semisupervised", "adaptive_semisupervised",
             "prior_lasso_support", "prior_lasso_value"),
  seeds = seq_len(10), scenario = 1, 
  nl = 100, nt = 10000, p = 100)
{
  # method = "supervised_lasso"; seeds = 1 : 5; i = 1
  # scenario = 2; nl = 100; nt = 10000; p = 100
  method <- match.arg(method)
  compute_estimate <- get(method)
  
  quality <- vector("list", length(seeds))
  for (i in seq_along(seeds)) {
    data.gen <- simulate_data_misspecified(scenario = scenario, 
                                           nl = nl, nt = nt, p = p, 
                                           seed = seeds[i], beta.fit = F)
    data <- data.gen$data
    beta_hat <- compute_estimate(
      x = data$x, y = data$y, 
      surrogate_index = data$surrogate_index,
      kfolds = 10L, fold_seed = seeds[i] + 10000L)
    
    test <- simulate_data_misspecified(scenario = scenario,
                                       nl = 10000L, nt = 10000L, p = p, 
                                       seed = seeds[i] + 100000L, beta.fit = T)
    quality[[i]] <- assess_quality(test$data, test$info$beta_true, beta_hat)
    
    print(i)
  }
  
  quality <- do.call("rbind", quality)
  quality <- data.table(
    scenario = scenario, method = method, 
    nl = nl, nt = nt, p = p,
    quality)
  quality
}


