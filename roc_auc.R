get_roc <- function(y_true, y_score) 
{
  i <- order(y_score, decreasing = TRUE)
  y_score <- y_score[i]
  y_true <- y_true[i]
  
  n <- length(y_true)
  t1 <- sum(y_true)
  t0 <- n - t1
  p1 <- seq_len(n)
  p0 <- n - p1

  true_positive <- cumsum(y_true)
  false_positive <- p1 - true_positive
  true_negative <- t0 - false_positive

  cut <- y_score
  pct <- p1 / n
  acc <- (true_positive + true_negative) / n
  
  tpr <- true_positive / t1
  fpr <- false_positive / t0
  tnr <- true_negative / t0
  
  ppv <- true_positive / p1
  fdr <- false_positive / p1
  npv <- true_negative / p0

  sen <- tpr
  rec <- tpr
  spec <- tnr
  prec <- ppv
  f1 <- prec * rec / (prec + rec)

  i <- which(diff(cut) < -1e-8)
  data.frame(
    cut = 0.5 * (cut[i] + cut[i + 1]), pct = pct[i], acc = acc[i],
    tpr = tpr[i], fpr = fpr[i], tnr = tnr[i],
    ppv = ppv[i], fdr = fdr[i], npv = npv[i],
    sen = sen[i], spec = spec[i], prec = prec[i], rec = rec[i],
    f1 = f1[i])
}


get_auc <- function(y_true, y_score) 
{
  t1 <- sum(y_true)
  t0 <- length(y_true) - t1
  r <- rank(y_score, ties.method = "average")
  auc <- (sum(r[y_true > 0.5]) - t1 * (t1 + 1) / 2) / t1 / t0
  auc
}


get_confusion_metrics <- function(y_true, y_pred)
{
  pct <- mean(y_pred > 0.5)
  acc <- mean((y_true > 0.5 & y_pred > 0.5) | 
                (y_true < 0.5 & y_pred < 0.5))
  
  tpr <- sum(y_true > 0.5 & y_pred > 0.5) / sum(y_true > 0.5)
  sen <- tpr
  rec <- sen
  tnr <- sum(y_true < 0.5 & y_pred < 0.5) / sum(y_true < 0.5)
  spec <- tnr
  fpr <- 1.0 - spec
  
  ppv <- sum(y_true > 0.5 & y_pred > 0.5) / sum(y_pred > 0.5)
  prec <- ppv
  fdr <- 1.0 - ppv
  npv <- sum(y_true < 0.5 & y_pred < 0.5) / sum(y_pred < 0.5)
  
  f1 <- 1.0 / (1.0 / rec + 1.0 / prec)

  c(pct = pct, acc = acc, 
    tpr = tpr, fpr = fpr, tnr = tnr, 
    ppv = ppv, fdr = fdr, npv = npv, 
    sen = sen, spec = spec, prec = prec, rec = rec,
    f1 = f1)
}


get_roc_at_threshold <- function(y_true, y_score, cut)
{
  result <- sapply(cut, function(c) 
    get_confusion_metrics(y_true, as.double(y_score > c)))
  result <- data.frame(cut = cut, t(result))
  result
}


get_deviance <- function(y_true, y_score)
{
  y_score <- pmin(pmax(y_score, 1e-8), 1 - 1e-8)
  mean(-2 * (y_true * log(y_score) + (1 - y_true) * log(1 - y_score)))
}

