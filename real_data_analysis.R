
source("simulations.R")


assess_quality_real <- function(data, beta_true, beta_hat)
{
  link_true <- as.double(data$x %*% beta_true[-1L]) + beta_true[1L]
  response_true <- plogis(link_true)
  
  link_hat <- as.double(data$x %*% beta_hat[-1L]) + beta_hat[1L]
  response_hat <- plogis(link_hat)
  
  c(auc = get_auc(data$y, response_hat),
    excess_risk = 0.5 * (get_deviance(data$y, response_hat) - 
                           get_deviance(data$y, response_true)),
    mse_link = mean((link_hat - link_true) ** 2),
    mse_response = mean((data$y - response_hat) ** 2),
    bias_response = mean(response_hat - response_true),
    var_response = var(response_hat - response_true),
    var_ratio_response = var(response_hat) / var(response_true),
    l1_error = sum(abs(beta_hat - beta_true)),
    l2_error = sqrt(sum(beta_hat - beta_true) ** 2),
    l0_norm = sum(abs(beta_hat) > 1e-8),
    l1_norm = sum(abs(beta_hat)))
}


adjust_note_count <- function(u, note_count)
{
  u <- u - mean(u)
  note_count <- note_count - mean(note_count)
  gamma <- sum(u * note_count) / sum(note_count ** 2)
  u - gamma * note_count
}


clean_data <- function(dt, count_normalization)
{
  dt <- dt[!is.na(note_count), ]
  dt[, surrogate := icd_main + nlp_main]
  dt[, icd_main := NULL]
  dt[, nlp_main := NULL]
  
  xn <- grep("^(surrogate|note_count|C[0-9]{7})", names(dt), 
             ignore.case = TRUE, value = TRUE)
  for (j in xn) {
    set(dt, j = j, value = log1p(dt[[j]]))
  }

  nonzero <- dt[!is.na(label), sapply(.SD, function(z) sum(z > 0.5))]
  nonzero <- names(nonzero)[nonzero >= 1L]

  xn <- intersect(xn, nonzero)
  if (count_normalization) {
    for (j in setdiff(xn, "note_count")) {
      set(dt, j = j, value = adjust_note_count(dt[[j]], dt[["note_count"]]))
    }
  }
  
  y <- dt[, label]
  x <- as.matrix(
    dt[, c("surrogate", setdiff(xn, "surrogate")), with = FALSE])
  
  label_table <- table(y, exclude = NULL)
  surrogate_auc <- get_auc(y[!is.na(y)], x[!is.na(y), 1L])

  list(y = y, x = x,
       label_table = label_table, surrogate_auc = surrogate_auc)
}


fit_full_labeled_data <- function(
  x, y, 
  method = c("supervised_lasso", "supervised_adaptive_lasso", 
             "surrogate_only", "perfect_surrogate", 
             "naive_semisupervised", "adaptive_semisupervised",
             "prior_lasso_support", "prior_lasso_value"))
{
  method <- match.arg(method)
  compute_estimate <- get(method)

  beta <- compute_estimate(
    x, y, surrogate_index = 1L, kfolds = 10L, fold_seed = 1L)
  beta <- as.double(beta)
  
  beta
}




train_and_evaluate <- function(
  x, y, 
  method = c('U_lasso', "supervised_lasso", "supervised_adaptive_lasso", 
             "surrogate_only", "perfect_surrogate", 
             "naive_semisupervised", "adaptive_semisupervised",
             "prior_lasso_support", "prior_lasso_value", 'U_lasso_SS'),
  beta_good = NULL, num_replications = 5, num_bootstrap = 20, train_size = 50)
{
  method <- match.arg(method)
  compute_estimate <- get(method)
  
  labeled <- which(!is.na(y))
  unlabeled <- which(is.na(y))
  alpha <- double(ncol(x))
  if (method %in% c("perfect_surrogate", 
                    "naive_semisupervised", "adaptive_semisupervised",
                    "prior_lasso_support", "prior_lasso_value")) {
    alpha <- adaptive_lasso_fit(
      x[, -1L], x[, 1L],
      family = "gaussian", tuning = "ic",
      # kfolds = kfolds, foldid = fold_seed,
      df_multiplier = "log(n)")
  }
  if (method %in% c('U_lasso', 'U_lasso_SS')){
    alpha <- U_lasso_direction(x, y, 1L)
  }
  
  names(alpha) <- c("intercept", colnames(x)[-1L])
  
  auc_x_alpha <- get_auc(
    y[labeled],
    drop(x[labeled, -1L] %*% alpha[-1L]))
  auc_s <- get_auc(
    y[labeled],
    x[labeled, 1L])
  auc_list <- c(auc_x_alpha = auc_x_alpha, auc_s = auc_s)
  
  n <- length(labeled)
  
  if (is.null(beta_good)) {
    beta_good <- double(ncol(x) + 1L)
  }
  
  result <- lapply(seq_len(num_replications), function(k) {
    #print(k)
    set.seed(1234L + k)
    sample_id <- sample(labeled, length(labeled), replace = F)
    fold_size <- as.integer(0.25 * length(labeled))
    
    eva_mat <- c()
    for (fold in 1:4) {
      i_test <- sample_id[(1 + (fold - 1) * fold_size): (fold * fold_size)]
      i_train_all <- setdiff(sample_id, i_test)
      
      eva_mat_fold <- c()
      for (j in 1:num_bootstrap) {
        set.seed(66666 + j + 1000 * k)
        i_train <- sample(i_train_all, train_size, replace = F)
        i_train <- c(i_train, unlabeled)
        beta <- compute_estimate(
          x[i_train, ], y[i_train], surrogate_index = 1L,
          kfolds = 10L, fold_seed = 5678L + k, alpha = alpha)
        
        test_data <- list(x = x[i_test, ], y = y[i_test])
        evaluate <- list(quality = assess_quality_real(test_data, beta_good, beta), 
                         beta = beta)
        eva_mat_fold <- rbind(eva_mat_fold, evaluate$quality)
      }
      eva_mat_fold <- as.data.frame(eva_mat_fold)
      eva_mat <- rbind(eva_mat,
                       c(mean(eva_mat_fold$auc), sd(eva_mat_fold$auc),
                         mean(eva_mat_fold$mse_response), sd(eva_mat_fold$mse_response)))
    }
    colMeans(eva_mat)
  })
  
  metric <- rbindlist(lapply(result, function(r) as.data.table(rbind((r)))))
  metric <- colMeans(metric)
  print(method)
  return(metric)
}




