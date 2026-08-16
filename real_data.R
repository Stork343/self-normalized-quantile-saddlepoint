## ================================================================
##  VaR 实证（完整精简版）：HS / FHS / Delta / EVT / MCt / QR / GARCH
##  - SNQESA 离散精确区间（HS/FHS）
##  - 其它方法的 CI 用 bootstrap + pbmcmapply 并行
##  - 统一回测：POF / IND / CC / Traffic Light
##  - 输出统一，含 time_sec
##  仅依赖：stats （可选 pbmcapply、tibble）
## ================================================================
library(pbmcapply)
library(parallel)
library(progress)
library(quantileCI)
library(quantreg)
library(rugarch)
library(ggplot2)
library(gridExtra)
library(scales)
library(patchwork)
library(dplyr)


## ---------- 0) 并行与极简小工具 ----------
# 基于你统一的 nwin、原始 r 与 idx 的定义方式：
make_L_next <- function(r, nwin){
  r <- as.numeric(r)
  idx <- seq_len(length(r) - nwin)
  L   <- -r
  # 与各 fit_* 里 out 的行数一致（长度 = length(idx)）
  L_next <- L[nwin + idx]
  return(L_next)
}

# 你的 L_next 是基于 nwin 后的“下一日损失”，其日期=原始 dates[nwin + idx]
make_aligned_dates <- function(dates, nwin){
  idx <- seq_len(length(dates) - nwin)
  as.Date(dates[nwin + idx])
}

#L_next <- make_L_next(spx$r, nwin)  # <<< 全局一次性构造

# 输入：aligned_dates（与 L_next 同长度）、命名的危机期列表
# 输出：一个命名 list，每个元素是与 aligned_dates 等长的 TRUE/FALSE 向量
build_crisis_masks <- function(aligned_dates, crisis_periods) {
  masks <- list()
  for (nm in names(crisis_periods)) {
    pr <- crisis_periods[[nm]]
    start <- as.Date(pr[1]); end <- as.Date(pr[2])
    masks[[nm]] <- (aligned_dates >= start) & (aligned_dates <= end)
  }
  masks
}


.has_pbmc <- function(){
  requireNamespace("pbmcapply", quietly = TRUE) && .Platform$OS.type != "windows"
}

boot_ci_parallel <- function(stat_sampler, B = 2000, alpha = 0.05, cores = 2, seed = NULL){
  stopifnot(is.function(stat_sampler), B >= 100, alpha > 0, alpha < 1)
  if (!is.null(seed)) set.seed(seed)
  if (.has_pbmc() && cores > 1L){
    vals <- as.numeric(parallel::mcmapply(function(i) stat_sampler(), 1:B,
                                             SIMPLIFY = TRUE, mc.cores = cores))
  } else {
    vals <- sapply(1:B, function(i) stat_sampler())
  }
  qs <- stats::quantile(vals, c(alpha/2, 1 - alpha/2), type = 8, na.rm = TRUE, names = FALSE)
  list(lower = qs[1], upper = qs[2], draws = vals)
}

snqesa_ci_discrete <- function(x, tau = 0.99, alpha = 0.05, split = c("equal","minlength")){
  split <- match.arg(split)
  x <- sort(as.numeric(x)); n <- length(x)
  pick <- function(aL){
    aL <- max(0, min(alpha, aL)); aR <- alpha - aL
    L <- stats::qbinom(aL, n, tau); U <- stats::qbinom(1 - aR, n, tau) + 1
    L <- min(n, max(1L, as.integer(L))); U <- min(n, max(1L, as.integer(U)))
    c(lower = x[L], upper = x[U])
  }
  if (split == "equal") pick(alpha/2) else {
    grid <- seq(0, alpha, length.out = max(31, round(10 + 4/alpha)))
    out  <- t(sapply(grid, pick)); out[which.min(out[,"upper"] - out[,"lower"]), ]
  }
}

ewma_sigma <- function(r, lambda = 0.94){
  r <- as.numeric(r); n <- length(r); s2 <- numeric(n)
  s2[1] <- stats::var(r, na.rm = TRUE)
  for (t in 2:n) s2[t] <- lambda*s2[t-1] + (1-lambda)*r[t-1]^2
  list(sigma = sqrt(s2), sigma_next = sqrt(lambda*s2 + (1-lambda)*r^2))
}

backtest_pof <- function(exceed, alpha_level){
  y <- as.integer(exceed); y <- y[is.finite(y)]
  n <- length(y); x <- sum(y); if (n <= 0) return(c(stat=NA, p.value=NA, pi=NA))
  phat <- x / n
  logL <- function(p) x*log(p) + (n-x)*log(1-p)
  LR   <- 2 * (logL(phat) - logL(alpha_level))
  pval <- 1 - stats::pchisq(LR, df = 1)
  c(stat = unname(LR), p.value = unname(pval), pi = unname(phat))
}

backtest_ind_cc <- function(exceed, alpha_level){
  y <- as.integer(exceed); y <- y[is.finite(y)]
  n <- length(y); if (n <= 1) return(c(LR_ind=NA, p_ind=NA, LR_cc=NA, p_cc=NA))
  ylag <- c(NA, head(y, -1)); keep <- is.finite(ylag); y <- y[keep]; ylag <- ylag[keep]
  n00 <- sum(ylag==0 & y==0); n01 <- sum(ylag==0 & y==1)
  n10 <- sum(ylag==1 & y==0); n11 <- sum(ylag==1 & y==1)
  n0 <- n00 + n01; n1 <- n10 + n11
  p01 <- if (n0>0) n01/n0 else 0; p11 <- if (n1>0) n11/n1 else 0
  phat <- (n01 + n11) / (n0 + n1 + 1e-15)
  L1 <- (1 - p01)^n00 * p01^n01 * (1 - p11)^n10 * p11^n11
  L0 <- (1 - phat)^(n00 + n10) * phat^(n01 + n11)
  LR_ind <- -2 * log((L0 + 1e-32)/(L1 + 1e-32))
  p_ind  <- 1 - stats::pchisq(LR_ind, df = 1)
  LR_pof <- backtest_pof(y, alpha_level = alpha_level)["stat"]
  LR_cc  <- as.numeric(LR_pof) + LR_ind
  p_cc   <- 1 - stats::pchisq(LR_cc, df = 2)
  c(LR_ind = LR_ind, p_ind = p_ind, LR_cc = LR_cc, p_cc = p_cc)
}

traffic_light <- function(exceed, window = 250, tau = 0.99){
  # Basel 交通灯：对 NA 做稳健处理（当作 0），避免 if 条件出现 NA
  thr <- if (abs(tau - 0.99) < 1e-8) {
    c(green = 4, yellow = 9)
  } else {
    c(green = ceiling((1 - tau) * window * 1.25),
      yellow = ceiling((1 - tau) * window * 2.25))
  }
  
  y <- as.integer(exceed)
  # 关键：把 NA 当作 0（不计入超损），避免 sum 得到 NA
  y[!is.finite(y)] <- 0L
  
  m <- length(y)
  if (m < window) {
    return(list(share = c(green = NA_real_, yellow = NA_real_, red = NA_real_)))
  }
  
  col <- character(m - window + 1)
  for (i in seq_len(m - window + 1)) {
    # 关键：na.rm = TRUE（保险）
    k <- sum(y[i:(i + window - 1)], na.rm = TRUE)
    col[i] <- if (k <= thr["green"]) "green" else if (k <= thr["yellow"]) "yellow" else "red"
  }
  sh <- prop.table(table(factor(col, levels = c("green","yellow","red"))))
  # 确保返回三色
  sh <- sh[c("green","yellow","red")]
  sh[is.na(sh)] <- 0
  list(share = sh)
}


## ---------- 1) 7 个方法（独立函数体） ----------
# 1) HS: 历史模拟（SNQESA 离散区间）
fit_HS <- function(r, dates, tau, nwin, alpha){
  t0 <- proc.time()[3]
  L <- -as.numeric(r); idx <- seq_len(length(L) - nwin)
  out <- data.frame(date = as.Date(dates[nwin + idx]),
                    VaR = NA_real_, ES = NA_real_, CI_low = NA_real_, CI_up = NA_real_)
  for (k in seq_along(idx)){
    win <- L[idx[k]:(idx[k]+nwin-1)]; win <- win[is.finite(win)]
    if (length(win) < 5) next
    q  <- as.numeric(stats::quantile(win, tau, type = 8, na.rm = TRUE))
    ci <- snqesa_ci_discrete(win, tau = tau, alpha = alpha, split = "minlength")
    es <- mean(win[win > q], na.rm = TRUE)
    out$VaR[k] <- q; out$ES[k] <- es; out$CI_low[k] <- ci["lower"]; out$CI_up[k] <- ci["upper"]
  }
  L_next <- L[nwin + idx]; ex <- as.integer(L_next > out$VaR)
  t1 <- proc.time()[3]
  
  list(
    name   = "HS", series = out, time_sec = unname(t1 - t0),
    backtest = list(POF = backtest_pof(ex, 1 - tau),
                    IND = backtest_ind_cc(ex, 1 - tau),
                    TL  = traffic_light(ex, 250, tau))
  )
}

res.HS = fit_HS(spx$r, spx$date, 0.99, 250, 0.05)

# format_result(res.HS)
# report_var_ci(res.HS)
# report_var_ci_timeseries(res.HS)
visualize_var_ci(res.HS, "HS")

# 2) FHS: EWMA 过滤历史模拟（SNQESA 离散区间）
fit_FHS <- function(r, dates, tau, nwin, alpha, lambda){
  t0 <- proc.time()[3]
  L <- -as.numeric(r); idx <- seq_len(length(L) - nwin)
  out <- data.frame(date = as.Date(dates[nwin + idx]),
                    VaR = NA_real_, ES = NA_real_, CI_low = NA_real_, CI_up = NA_real_)
  for (k in seq_along(idx)){
    i <- idx[k]
    winL <- L[i:(i+nwin-1)]
    ew   <- ewma_sigma(r[i:(i+nwin-1)], lambda = lambda)
    z    <- winL / ew$sigma; z <- z[is.finite(z)]
    if (length(z) < 5) next
    qz  <- as.numeric(stats::quantile(z, tau, type = 8, na.rm = TRUE))
    ciZ <- snqesa_ci_discrete(z, tau = tau, alpha = alpha, split = "minlength")
    sig_next <- ew$sigma_next[nwin]
    out$VaR[k]    <- qz * sig_next
    out$ES[k]     <- mean(z[z > qz], na.rm = TRUE) * sig_next
    out$CI_low[k] <- ciZ["lower"] * sig_next
    out$CI_up[k]  <- ciZ["upper"] * sig_next
  }
  L_next <- L[nwin + idx]; ex <- as.integer(L_next > out$VaR)
  t1 <- proc.time()[3]
  list(
    name   = "FHS", series = out, time_sec = unname(t1 - t0),
    backtest = list(POF = backtest_pof(ex, 1 - tau),
                    IND = backtest_ind_cc(ex, 1 - tau),
                    TL  = traffic_light(ex, 250, tau))
  )
}

res.FHS = fit_FHS(spx$r, spx$date, 0.99, 250, 0.05, 0.971)

# format_result(res.FHS)
# report_var_ci(res.FHS)
# report_var_ci_timeseries(res.FHS)
visualize_var_ci(res.FHS, "FHS")


# 3) Delta-Normal：正态参数法（bootstrap CI）
fit_Delta <- function(r, dates, tau, nwin, alpha, B_ci = 2000, cores_ci = 10){
  t0 <- proc.time()[3]
  L <- -as.numeric(r); idx <- seq_len(length(L) - nwin)
  z <- stats::qnorm(tau)
  out <- data.frame(date = as.Date(dates[nwin + idx]),
                    VaR = NA_real_, ES = NA_real_, CI_low = NA_real_, CI_up = NA_real_)
  cat(sprintf("开始处理 %d 个滚动窗口...\n", length(idx)))
  start_time <- Sys.time()
  
  # 使用progress包，它会自动计算和显示剩余时间
  pb <- progress_bar$new(
    format = "处理中 [:bar] :percent | 已用: :elapsed | 剩余: :eta",
    total = length(idx),
    clear = FALSE,
    width = 80,
    show_after = 0  # 立即显示
  )
  
  for (k in seq_along(idx)){
    win <- L[idx[k]:(idx[k]+nwin-1)]; win <- win[is.finite(win)]
    if (length(win) < 5) next
    mu <- mean(win); sdv <- stats::sd(win)
    out$VaR[k] <- mu + z*sdv
    out$ES[k]  <- mu + sdv * stats::dnorm(z) / (1 - tau)
    n  <- length(win)
    stat_sampler <- function(){
      xb <- mu + sdv * stats::rnorm(n)
      mean(xb) + z*stats::sd(xb)
    }
    ci <- boot_ci_parallel(stat_sampler, B = B_ci, alpha = alpha, cores = cores_ci)
    out$CI_low[k] <- ci$lower; out$CI_up[k] <- ci$upper
    
    pb$tick()
  }
  
  L_next <- L[nwin + idx]; ex <- as.integer(L_next > out$VaR)
  total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  cat(sprintf("\n✅ 完成! 总耗时: %.1f 秒\n", total_time))
  
  t1 <- proc.time()[3]
  list(
    name   = "Delta", series = out, time_sec = unname(t1 - t0),
    backtest = list(POF = backtest_pof(ex, 1 - tau),
                    IND = backtest_ind_cc(ex, 1 - tau),
                    TL  = traffic_light(ex, 250, tau))
  )
}

res.Delta = fit_Delta(spx$r, spx$date, 0.99, 250, 0.05, 1000, 10)

format_result(res.Delta)
report_var_ci(res.Delta)
report_var_ci_timeseries(res.Delta)
visualize_var_ci(res.Delta, "Delta")


# 4) MCt：t 分布蒙特卡洛（bootstrap CI）
fit_MCt <- function(r, dates, tau, nwin, df = 5, alpha = 0.05, B_ci = 2000, cores_ci = 4){
  if (df <= 2) stop("自由度df必须大于2")
  if (tau <= 0 || tau >= 1) stop("tau必须在0和1之间")
  
  
  t0 <- proc.time()[3]
  L <- -as.numeric(r); idx <- seq_len(length(L) - nwin)
  tcut <- stats::qt(tau, df = df)
  out <- data.frame(date = as.Date(dates[nwin + idx]),
                    VaR = NA_real_, ES = NA_real_, CI_low = NA_real_, CI_up = NA_real_)
  
  # 检查是否安装了progress包
  if (requireNamespace("progress", quietly = TRUE)) {
    library(progress)
    
    cat(sprintf("开始处理 %d 个滚动窗口 (MCt方法, df=%d)...\n", length(idx), df))
    
    # 使用progress包创建进度条
    pb <- progress_bar$new(
      format = "处理中 [:bar] :percent | 已用: :elapsed | 剩余: :eta",
      total = length(idx),
      clear = FALSE,
      width = 80,
      show_after = 0  # 立即显示
    )
    
    for (k in seq_along(idx)){
      win <- L[idx[k]:(idx[k]+nwin-1)]; win <- win[is.finite(win)]
      if (length(win) < 5) next
      mu <- mean(win); sdv <- stats::sd(win)
      out$VaR[k] <- mu + sdv * tcut
      out$ES[k]  <- mu + sdv * ((df + tcut^2)/((df - 1)*(1 - tau))) * stats::dt(tcut, df = df)
      n  <- length(win)
      stat_sampler <- function(){
        xb <- mu + sdv * stats::rt(n, df = df)
        mean(xb) + tcut*stats::sd(xb)
      }
      ci <- boot_ci_parallel(stat_sampler, B = B_ci, alpha = alpha, cores = cores_ci)
      out$CI_low[k] <- ci$lower; out$CI_up[k] <- ci$upper
      
      pb$tick()  # 更新进度条
    }
    
  } else {
    # 回退到基础版本（无progress包）
    cat(sprintf("开始处理 %d 个滚动窗口 (MCt方法, df=%d)...\n", length(idx), df))
    
    # 创建基础进度条
    pb <- txtProgressBar(min = 0, max = length(idx), style = 3, width = 50)
    
    for (k in seq_along(idx)){
      win <- L[idx[k]:(idx[k]+nwin-1)]; win <- win[is.finite(win)]
      if (length(win) < 5) next
      mu <- mean(win); sdv <- stats::sd(win)
      out$VaR[k] <- mu + sdv * tcut
      out$ES[k]  <- mu + sdv * ((df + tcut^2)/((df - 1)*(1 - tau))) * stats::dt(tcut, df = df)
      n  <- length(win)
      stat_sampler <- function(){
        xb <- mu + sdv * stats::rt(n, df = df)
        mean(xb) + tcut*stats::sd(xb)
      }
      ci <- boot_ci_parallel(stat_sampler, B = B_ci, alpha = alpha, cores = cores_ci)
      out$CI_low[k] <- ci$lower; out$CI_up[k] <- ci$upper
      
      setTxtProgressBar(pb, k)  # 更新基础进度条
    }
    
    close(pb)  # 关闭基础进度条
  }
  
  L_next <- L[nwin + idx]; ex <- as.integer(L_next > out$VaR)
  t1 <- proc.time()[3]
  
  # 输出完成信息
  cat(sprintf("✅ MCt方法完成! 总耗时: %.1f 秒\n", unname(t1 - t0)))
  
  list(
    name   = "MCt", 
    series = out, 
    time_sec = unname(t1 - t0),
    parameters = list(tau = tau, nwin = nwin, df = df, alpha = alpha),
    backtest = list(
      POF = backtest_pof(ex, 1 - tau),
      IND = backtest_ind_cc(ex, 1 - tau),
      TL  = traffic_light(ex, 250, tau)
    )
  )
}

res.MCt = fit_MCt(spx$r, spx$date, 0.99, 250, B_ci = 1000, cores_ci = 10)

format_result(res.MCt)
report_var_ci(res.MCt)
report_var_ci_timeseries(res.MCt)
visualize_var_ci(res.MCt, "MCt")


# 5) QR：简化分位回归（两步法：LS + 残差分位数；bootstrap CI）
#     L_t ≈ b0 + b1 * sigma_EWMA + e；预测 VaR = b0 + b1*sigma_next + q_tau(e)
fit_QR <- function(r, dates, tau, nwin, alpha, lambda = 0.95, B_ci = 2000, cores_ci = 10){
  if (lambda <= 0 || lambda >= 1) stop("lambda必须在0和1之间")
  if (tau <= 0 || tau >= 1) stop("tau必须在0和1之间")
  t0 <- proc.time()[3]
  L <- -as.numeric(r); idx <- seq_len(length(L) - nwin)
  out <- data.frame(date = as.Date(dates[nwin + idx]),
                    VaR = NA_real_, ES = NA_real_, CI_low = NA_real_, CI_up = NA_real_)
  
  # 检查是否安装了progress包
  if (requireNamespace("progress", quietly = TRUE)) {
    library(progress)
    
    cat(sprintf("开始处理 %d 个滚动窗口 (QR方法, lambda=%.2f)...\n", length(idx), lambda))
    
    # 使用progress包创建进度条
    pb <- progress_bar$new(
      format = "处理中 [:bar] :percent | 已用: :elapsed | 剩余: :eta",
      total = length(idx),
      clear = FALSE,
      width = 80,
      show_after = 0  # 立即显示
    )
    
    for (k in seq_along(idx)){
      i <- idx[k]
      winL <- L[i:(i+nwin-1)]
      ew   <- ewma_sigma(r[i:(i+nwin-1)], lambda = lambda)
      s    <- ew$sigma
      keep <- is.finite(winL) & is.finite(s)
      winL <- winL[keep]; s <- s[keep]
      if (length(winL) < 5) next
      fit <- stats::lm(winL ~ s)  # LS slope
      e   <- fit$residuals
      q_e <- as.numeric(stats::quantile(e, tau, type = 8, na.rm = TRUE))
      sig_next <- ew$sigma_next[nwin]
      pred_mu  <- coef(fit)[1] + coef(fit)[2] * sig_next
      out$VaR[k] <- pred_mu + q_e
      out$ES[k]  <- pred_mu + mean(e[e > q_e], na.rm = TRUE)
      # pairs bootstrap (L,s) -> LS -> residual τ-quantile -> VaR_next
      n  <- length(winL)
      stat_sampler <- function(){
        id <- sample.int(n, n, replace = TRUE)
        fb <- stats::lm(winL[id] ~ s[id])
        eb <- winL[id] - (coef(fb)[1] + coef(fb)[2]*s[id])
        qb <- as.numeric(stats::quantile(eb, tau, type = 8, na.rm = TRUE))
        (coef(fb)[1] + coef(fb)[2]*sig_next) + qb
      }
      ci <- boot_ci_parallel(stat_sampler, B = B_ci, alpha = alpha, cores = cores_ci)
      out$CI_low[k] <- ci$lower; out$CI_up[k] <- ci$upper
      
      pb$tick()  # 更新进度条
    }
    
  } else {
    # 回退到基础版本（无progress包）
    cat(sprintf("开始处理 %d 个滚动窗口 (QR方法, lambda=%.2f)...\n", length(idx), lambda))
    
    # 创建基础进度条
    pb <- txtProgressBar(min = 0, max = length(idx), style = 3, width = 50)
    
    for (k in seq_along(idx)){
      i <- idx[k]
      winL <- L[i:(i+nwin-1)]
      ew   <- ewma_sigma(r[i:(i+nwin-1)], lambda = lambda)
      s    <- ew$sigma
      keep <- is.finite(winL) & is.finite(s)
      winL <- winL[keep]; s <- s[keep]
      if (length(winL) < 5) next
      fit <- stats::lm(winL ~ s)  # LS slope
      e   <- fit$residuals
      q_e <- as.numeric(stats::quantile(e, tau, type = 8, na.rm = TRUE))
      sig_next <- ew$sigma_next[nwin]
      pred_mu  <- coef(fit)[1] + coef(fit)[2] * sig_next
      out$VaR[k] <- pred_mu + q_e
      out$ES[k]  <- pred_mu + mean(e[e > q_e], na.rm = TRUE)
      # pairs bootstrap (L,s) -> LS -> residual τ-quantile -> VaR_next
      n  <- length(winL)
      stat_sampler <- function(){
        id <- sample.int(n, n, replace = TRUE)
        fb <- stats::lm(winL[id] ~ s[id])
        eb <- winL[id] - (coef(fb)[1] + coef(fb)[2]*s[id])
        qb <- as.numeric(stats::quantile(eb, tau, type = 8, na.rm = TRUE))
        (coef(fb)[1] + coef(fb)[2]*sig_next) + qb
      }
      ci <- boot_ci_parallel(stat_sampler, B = B_ci, alpha = alpha, cores = cores_ci)
      out$CI_low[k] <- ci$lower; out$CI_up[k] <- ci$upper
      
      setTxtProgressBar(pb, k)  # 更新基础进度条
    }
    
    close(pb)  # 关闭基础进度条
  }
  
  L_next <- L[nwin + idx]; ex <- as.integer(L_next > out$VaR)
  t1 <- proc.time()[3]
  
  # 输出完成信息
  cat(sprintf("✅ QR方法完成! 总耗时: %.1f 秒\n", unname(t1 - t0)))
  
  list(
    name   = "QR", 
    series = out, 
    time_sec = unname(t1 - t0),
    parameters = list(tau = tau, nwin = nwin, lambda = lambda, alpha = alpha),
    backtest = list(
      POF = backtest_pof(ex, 1 - tau),
      IND = backtest_ind_cc(ex, 1 - tau),
      TL  = traffic_light(ex, 250, tau)
    )
  )
}

res.QR = fit_QR(spx$r, spx$date, 0.99, 250, 0.05, 0.95, 1000, 10)


format_result(res.QR)
report_var_ci(res.QR)



# 6) EVT：POT + GPD（MLE），高阈值 u=Q_{evt_uq}；bootstrap CI（重拟合）
#     VaR_tau = u + beta/xi * ( (( (1 - tau)/(1 - F(u)) )^(-xi) - 1) ),  |xi|>0
#     ES_tau  = ( VaR_tau - xi*u + beta ) / (1 - xi),  xi < 1
fit_EVT <- function(r, dates, tau, nwin, alpha, evt_uq = 0.9, B_ci = 400, cores_ci = 10){
  t0 <- proc.time()[3]
  L <- -as.numeric(r); idx <- seq_len(length(L) - nwin)
  out <- data.frame(date = as.Date(dates[nwin + idx]),
                    VaR = NA_real_, ES = NA_real_, CI_low = NA_real_, CI_up = NA_real_)
  
  gpd_mle <- function(y){ # y: exceedances > 0
    if (length(y) < 30) return(c(xi=NA, beta=NA))
    # 参数 (xi, logbeta) 以避免 beta<=0
    nll <- function(par){
      xi <- par[1]; lb <- par[2]; beta <- exp(lb)
      if (beta <= 0) return(1e12)
      t1 <- 1 + xi*y/beta
      if (any(t1 <= 0)) return(1e12)
      sum(lb + (1/xi + 1)*log(t1))
    }
    par0 <- c(0.1, log(stats::sd(y)))
    fit  <- tryCatch(stats::optim(par0, nll, method="Nelder-Mead"),
                     error = function(e) NULL)
    if (is.null(fit) || fit$convergence != 0) return(c(xi=NA, beta=NA))
    xi <- fit$par[1]; beta <- exp(fit$par[2])
    c(xi=xi, beta=beta)
  }
  
  q_evt <- function(u, xi, beta, p_tail_ratio){
    # p_tail_ratio = (1 - tau)/(1 - F(u))
    if (!is.finite(xi) || !is.finite(beta) || beta <= 0) return(NA_real_)
    if (abs(xi) < 1e-6){
      u + beta*log(1/p_tail_ratio)
    } else {
      u + (beta/xi)*(p_tail_ratio^(-xi) - 1)
    }
  }
  
  es_evt <- function(VaR, u, xi, beta){
    if (!is.finite(VaR) || !is.finite(xi) || !is.finite(beta)) return(NA_real_)
    if (xi >= 1) return(NA_real_)
    (VaR - xi*u + beta)/(1 - xi)
  }
  
  # 检查是否安装了progress包
  if (requireNamespace("progress", quietly = TRUE)) {
    library(progress)
    
    cat(sprintf("开始处理 %d 个滚动窗口 (EVT方法, u_quantile=%.2f)...\n", length(idx), evt_uq))
    
    # 使用progress包创建进度条
    pb <- progress_bar$new(
      format = "处理中 [:bar] :percent | 已用: :elapsed | 剩余: :eta",
      total = length(idx),
      clear = FALSE,
      width = 80,
      show_after = 0  # 立即显示
    )
    
    for (k in seq_along(idx)){
      win <- L[idx[k]:(idx[k]+nwin-1)]
      win <- win[is.finite(win)]; n <- length(win)
      if (n < 60) {
        pb$tick()  # 即使跳过也要更新进度条
        next
      }
      u <- as.numeric(stats::quantile(win, evt_uq, type = 8, na.rm = TRUE))
      y <- win[win > u] - u; m <- length(y)
      if (m < 30) {
        pb$tick()
        next
      }
      par <- gpd_mle(y); xi <- par["xi"]; beta <- par["beta"]
      if (!is.finite(xi) || !is.finite(beta)) {
        pb$tick()
        next
      }
      p_tail_ratio <- (1 - tau) / (1 - evt_uq)
      VaR <- q_evt(u, xi, beta, p_tail_ratio)
      ES  <- es_evt(VaR, u, xi, beta)
      out$VaR[k] <- VaR; out$ES[k] <- ES
      # parametric bootstrap: 重采样 y* ~ GPD(xi,beta)，重拟合 -> VaR*
      rGPD <- function(n, xi, beta){
        u <- stats::runif(n)
        if (abs(xi) < 1e-6){
          -beta*log(1 - u)
        } else {
          beta/xi * ((1 - u)^(-xi) - 1)
        }
      }
      stat_sampler <- function(){
        yb <- rGPD(m, xi, beta)
        par_b <- gpd_mle(yb)
        xi_b <- par_b["xi"]; beta_b <- par_b["beta"]
        q_evt(u, xi_b, beta_b, p_tail_ratio)
      }
      ci <- boot_ci_parallel(stat_sampler, B = B_ci, alpha = alpha, cores = cores_ci)
      out$CI_low[k] <- ci$lower; out$CI_up[k] <- ci$upper
      
      pb$tick()  # 更新进度条
    }
    
  } else {
    # 回退到基础版本（无progress包）
    cat(sprintf("开始处理 %d 个滚动窗口 (EVT方法, u_quantile=%.2f)...\n", length(idx), evt_uq))
    
    # 创建基础进度条
    pb <- txtProgressBar(min = 0, max = length(idx), style = 3, width = 50)
    
    for (k in seq_along(idx)){
      win <- L[idx[k]:(idx[k]+nwin-1)]
      win <- win[is.finite(win)]; n <- length(win)
      if (n < 60) {
        setTxtProgressBar(pb, k)
        next
      }
      u <- as.numeric(stats::quantile(win, evt_uq, type = 8, na.rm = TRUE))
      y <- win[win > u] - u; m <- length(y)
      if (m < 30) {
        setTxtProgressBar(pb, k)
        next
      }
      par <- gpd_mle(y); xi <- par["xi"]; beta <- par["beta"]
      if (!is.finite(xi) || !is.finite(beta)) {
        setTxtProgressBar(pb, k)
        next
      }
      p_tail_ratio <- (1 - tau) / (1 - evt_uq)
      VaR <- q_evt(u, xi, beta, p_tail_ratio)
      ES  <- es_evt(VaR, u, xi, beta)
      out$VaR[k] <- VaR; out$ES[k] <- ES
      # parametric bootstrap: 重采样 y* ~ GPD(xi,beta)，重拟合 -> VaR*
      rGPD <- function(n, xi, beta){
        u <- stats::runif(n)
        if (abs(xi) < 1e-6){
          -beta*log(1 - u)
        } else {
          beta/xi * ((1 - u)^(-xi) - 1)
        }
      }
      stat_sampler <- function(){
        yb <- rGPD(m, xi, beta)
        par_b <- gpd_mle(yb)
        xi_b <- par_b["xi"]; beta_b <- par_b["beta"]
        q_evt(u, xi_b, beta_b, p_tail_ratio)
      }
      ci <- boot_ci_parallel(stat_sampler, B = B_ci, alpha = alpha, cores = cores_ci)
      out$CI_low[k] <- ci$lower; out$CI_up[k] <- ci$upper
      
      setTxtProgressBar(pb, k)  # 更新基础进度条
    }
    
    close(pb)  # 关闭基础进度条
  }
  
  L_next <- L[nwin + idx]; ex <- as.integer(L_next > out$VaR)
  t1 <- proc.time()[3]
  
  # 输出完成信息
  cat(sprintf("✅ EVT方法完成! 总耗时: %.1f 秒\n", unname(t1 - t0)))
  
  list(
    name   = "EVT", 
    series = out, 
    time_sec = unname(t1 - t0),
    backtest = list(
      POF = backtest_pof(ex, 1 - tau),
      IND = backtest_ind_cc(ex, 1 - tau),
      TL  = traffic_light(ex, 250, tau)
    )
  )
}


res.EVT = fit_EVT(spx$r, spx$date, 0.99, 250, evt_uq = 0.7, 0.05)

format_result(res.EVT)
report_var_ci(res.EVT)


# 7) GARCH(1,1)（Quasi-ML, Normal），bootstrap CI（参数重估）
fit_GARCH <- function(r, dates, tau, nwin, alpha, B_ci = 400, cores_ci = 10){
  t0 <- proc.time()[3]
  idx <- seq_len(length(r) - nwin)
  out <- data.frame(date = as.Date(dates[nwin + idx]),
                    VaR = NA_real_, ES = NA_real_, CI_low = NA_real_, CI_up = NA_real_)
  qz <- stats::qnorm(tau)
  
  garch_fit <- function(x){ # 简化 QML 正态
    x <- as.numeric(x); n <- length(x); if (n < 30) return(NULL)
    nll <- function(par){
      w  <- exp(par[1]); a <- stats::plogis(par[2]) * 0.999
      b  <- stats::plogis(par[3]) * (0.999 - a)
      h  <- numeric(n); h[1] <- stats::var(x)
      for (t in 2:n) h[t] <- w + a*x[t-1]^2 + b*h[t-1]
      if (any(h <= 0) || !is.finite(sum(h))) return(1e12)
      0.5*sum(log(2*pi) + log(h) + x^2/h)
    }
    par0 <- c(log(var(x)*0.01), 0, 0)
    fit  <- tryCatch(stats::optim(par0, nll, method="Nelder-Mead"),
                     error=function(e) NULL)
    if (is.null(fit) || fit$convergence != 0) return(NULL)
    w <- exp(fit$par[1]); a <- stats::plogis(fit$par[2]) * 0.999
    b <- stats::plogis(fit$par[3]) * (0.999 - a)
    list(omega = w, alpha = a, beta = b)
  }
  
  garch_sigma_next <- function(x, par){
    x <- as.numeric(x); n <- length(x); h <- numeric(n)
    h[1] <- stats::var(x)
    for (t in 2:n) h[t] <- par$omega + par$alpha*x[t-1]^2 + par$beta*h[t-1]
    sqrt(par$omega + par$alpha*x[n]^2 + par$beta*h[n])
  }
  
  # 检查是否安装了progress包
  if (requireNamespace("progress", quietly = TRUE)) {
    library(progress)
    
    cat(sprintf("开始处理 %d 个滚动窗口 (GARCH方法)...\n", length(idx)))
    
    # 使用progress包创建进度条
    pb <- progress_bar$new(
      format = "处理中 [:bar] :percent | 已用: :elapsed | 剩余: :eta",
      total = length(idx),
      clear = FALSE,
      width = 80,
      show_after = 0  # 立即显示
    )
    
    for (k in seq_along(idx)){
      x <- as.numeric(r[idx[k]:(idx[k]+nwin-1)]); 
      if (sum(is.finite(x)) < 30) {
        pb$tick()
        next
      }
      par <- garch_fit(x); 
      if (is.null(par)) {
        pb$tick()
        next
      }
      sig_next <- garch_sigma_next(x, par)
      VaR <- qz * sig_next * 1  # mean≈0，VaR 在损失上 = -return 的右尾
      ES  <- sig_next * stats::dnorm(qz) / (1 - tau)
      out$VaR[k] <- VaR; out$ES[k] <- ES
      # parametric bootstrap: 重抽样 N(0,1)，重估参数，再算 VaR_next
      n  <- length(x)
      stat_sampler <- function(){
        xb <- stats::rnorm(n) * sig_next  # 简化：用当前 sig_next 近似波动（更快）
        par_b <- garch_fit(xb)
        if (is.null(par_b)) return(NA_real_)
        sig_b <- garch_sigma_next(xb, par_b)
        qz * sig_b
      }
      ci <- boot_ci_parallel(stat_sampler, B = B_ci, alpha = alpha, cores = cores_ci)
      out$CI_low[k] <- ci$lower; out$CI_up[k] <- ci$upper
      
      pb$tick()  # 更新进度条
    }
    
  } else {
    # 回退到基础版本（无progress包）
    cat(sprintf("开始处理 %d 个滚动窗口 (GARCH方法)...\n", length(idx)))
    
    # 创建基础进度条
    pb <- txtProgressBar(min = 0, max = length(idx), style = 3, width = 50)
    
    for (k in seq_along(idx)){
      x <- as.numeric(r[idx[k]:(idx[k]+nwin-1)]); 
      if (sum(is.finite(x)) < 30) {
        setTxtProgressBar(pb, k)
        next
      }
      par <- garch_fit(x); 
      if (is.null(par)) {
        setTxtProgressBar(pb, k)
        next
      }
      sig_next <- garch_sigma_next(x, par)
      VaR <- qz * sig_next * 1  # mean≈0，VaR 在损失上 = -return 的右尾
      ES  <- sig_next * stats::dnorm(qz) / (1 - tau)
      out$VaR[k] <- VaR; out$ES[k] <- ES
      # parametric bootstrap: 重抽样 N(0,1)，重估参数，再算 VaR_next
      n  <- length(x)
      stat_sampler <- function(){
        xb <- stats::rnorm(n) * sig_next  # 简化：用当前 sig_next 近似波动（更快）
        par_b <- garch_fit(xb)
        if (is.null(par_b)) return(NA_real_)
        sig_b <- garch_sigma_next(xb, par_b)
        qz * sig_b
      }
      ci <- boot_ci_parallel(stat_sampler, B = B_ci, alpha = alpha, cores = cores_ci)
      out$CI_low[k] <- ci$lower; out$CI_up[k] <- ci$upper
      
      setTxtProgressBar(pb, k)  # 更新基础进度条
    }
    
    close(pb)  # 关闭基础进度条
  }
  
  L <- -as.numeric(r); L_next <- L[nwin + idx]; ex <- as.integer(L_next > out$VaR)
  t1 <- proc.time()[3]
  
  # 输出完成信息
  cat(sprintf("✅ GARCH方法完成! 总耗时: %.1f 秒\n", unname(t1 - t0)))
  
  list(
    name   = "GARCH", 
    series = out, 
    time_sec = unname(t1 - t0),
    backtest = list(
      POF = backtest_pof(ex, 1 - tau),
      IND = backtest_ind_cc(ex, 1 - tau),
      TL  = traffic_light(ex, 250, tau)
    )
  )
}

res.GARCH = fit_GARCH(spx$r, spx$date, 0.99, 250, 0.05)


format_result(res.GARCH)
report_var_ci(res.GARCH)


# 8) EWHS：指数加权历史模拟（weighted quantile；bootstrap CI）
fit_EWHS <- function(r, dates, tau, nwin, alpha,
                     lambda = 0.97, B_ci = 1000, cores_ci = 10){
  stopifnot(tau > 0, tau < 1, lambda > 0, lambda < 1)
  t0 <- proc.time()[3]
  L   <- -as.numeric(r)
  idx <- seq_len(length(L) - nwin)
  out <- data.frame(date = as.Date(dates[nwin + idx]),
                    VaR = NA_real_, ES = NA_real_,
                    CI_low = NA_real_, CI_up = NA_real_)
  
  # 加权分位数（线性插值，严格单调的累计权重）
  wq <- function(x, w, p){
    stopifnot(length(x) == length(w))
    o   <- order(x); x <- x[o]; w <- w[o]
    w   <- w / sum(w)
    cw  <- cumsum(w)
    if (p <= cw[1]) return(x[1])
    if (p >= cw[length(cw)]) return(x[length(x)])
    j   <- findInterval(p, cw)             # cw[j] <= p < cw[j+1]
    tlo <- cw[j]; thi <- cw[j+1]
    xlo <- x[j];  xhi <- x[j+1]
    xlo + (p - tlo) / (thi - tlo) * (xhi - xlo)
  }
  
  # 指数权重：越近权重越大（最后一条=1）
  make_weights <- function(n){
    w <- lambda^rev(seq_len(n) - 1L)
    w / sum(w)
  }
  
  # 进度条
  use_prog <- requireNamespace("progress", quietly = TRUE)
  if (use_prog) {
    pb <- progress::progress_bar$new(
      format = "EWHS [:bar] :percent | :elapsed eta: :eta",
      total = length(idx), clear = FALSE, width = 80, show_after = 0
    )
  }
  
  for (k in seq_along(idx)){
    i   <- idx[k]
    win <- L[i:(i+nwin-1)]
    win <- win[is.finite(win)]
    n   <- length(win)
    if (n < 10){
      if (use_prog) pb$tick()
      next
    }
    w <- make_weights(n)
    
    # VaR & ES（加权）
    q  <- wq(win, w, tau)
    m  <- win > q
    es <- if (any(m)) sum(w[m] * win[m]) / sum(w[m]) else NA_real_
    
    # CI：加权自助法（按 w 概率有放回抽样，重算加权分位数）
    stat_sampler <- function(){
      id <- sample.int(n, n, replace = TRUE, prob = w)
      xb <- win[id]
      # 抽样后用“等权”分位数 ≡ 原加权经验分布的自助重采样
      # 也可用相同权重 w[id] 再做一次 wq，二者在大样本下等价
      wq(xb, rep(1/n, n), tau)
    }
    ci <- boot_ci_parallel(stat_sampler, B = B_ci, alpha = alpha, cores = cores_ci)
    
    out$VaR[k]    <- q
    out$ES[k]     <- es
    out$CI_low[k] <- ci$lower
    out$CI_up[k]  <- ci$upper
    
    if (use_prog) pb$tick()
  }
  
  L_next <- L[nwin + idx]
  ex     <- as.integer(L_next > out$VaR)
  t1     <- proc.time()[3]
  
  list(
    name   = "EWHS",
    series = out,
    time_sec = unname(t1 - t0),
    parameters = list(tau = tau, nwin = nwin, lambda = lambda, alpha = alpha),
    backtest = list(
      POF = backtest_pof(ex, 1 - tau),
      IND = backtest_ind_cc(ex, 1 - tau),
      TL  = traffic_light(ex, 250, tau)
    )
  )
}


res.EWHS = fit_EWHS(spx$r, spx$date, 0.99, 250, 0.05, lambda = 0.96)


format_result(res.EWHS)
report_var_ci(res.EWHS)


## ---------- 2) 统一调度与汇报 ----------
format_result <- function(results) {
  cat("=== 风险模型回测结果 ===\n\n")
  
  # 计算时间
  cat(sprintf("计算时间: %.3f 秒\n\n", results$time_sec))
  
  # 失败比例检验
  cat("1. 失败比例检验 (POF):\n")
  cat(sprintf("   统计量(POF.stat): %.4f\n", results$backtest$POF[[1]]))
  cat(sprintf("   P-value: %.4f\n", results$backtest$POF[[2]]))
  cat(sprintf("   实际失败率(pi): %.2f%%\n\n", results$backtest$POF[[3]] * 100))
  
  # 独立性检验
  cat("2. 独立性检验 (IND):\n")
  cat(sprintf("   独立性统计量(IND.stat): %.4f\n", results$backtest$IND[[1]]))
  cat(sprintf("   P-value: %.4f\n", results$backtest$IND[[2]]))
  cat(sprintf("   条件覆盖统计量(LR_cc): %.4f\n", results$backtest$IND[[3]]))
  cat(sprintf("   P-value: %.4f\n\n", results$backtest$IND[[4]]))
  
  # 交通灯检验
  cat("3. 交通灯检验 (TL):\n")
  cat(sprintf("   绿色区域: %.2f%%\n", results$backtest$TL$share[[1]] * 100))
  cat(sprintf("   黄色区域: %.2f%%\n", results$backtest$TL$share[[2]] * 100))
  cat(sprintf("   红色区域: %.2f%%\n", results$backtest$TL$share[[3]] * 100))
}

report_var_ci <- function(result, method_name = NULL) {
  if (is.null(method_name)) {
    method_name <- result$name
  }
  
  series <- result$series
  series_clean <- series[complete.cases(series$VaR, series$CI_low, series$CI_up), ]
  
  cat("=== ", method_name, "方法 VaR 置信区间报告 ===\n\n")
  
  # 基本统计
  cat("1. 置信区间基本统计:\n")
  ci_width <- series_clean$CI_up - series_clean$CI_low
  cat(sprintf("   平均置信区间宽度: %.4f\n", mean(ci_width, na.rm = TRUE)))
  cat(sprintf("   置信区间宽度标准差: %.4f\n", sd(ci_width, na.rm = TRUE)))
  cat(sprintf("   最小置信区间宽度: %.4f\n", min(ci_width, na.rm = TRUE)))
  cat(sprintf("   最大置信区间宽度: %.4f\n", max(ci_width, na.rm = TRUE)))
  cat(sprintf("   覆盖观测点数: %d\n", nrow(series_clean)))
  
  # 最近时期的置信区间
  cat("\n2. 最近时期 VaR 置信区间:\n")
  last_point <- tail(series_clean, 1)
  cat(sprintf("   日期: %s\n", last_point$date))
  cat(sprintf("   VaR 估计: %.4f\n", last_point$VaR))
  cat(sprintf("   95%% 置信区间: [%.4f, %.4f]\n", 
              last_point$CI_low, last_point$CI_up))
  cat(sprintf("   区间宽度: %.4f\n", last_point$CI_up - last_point$CI_low))
  
  return(invisible(list(
    summary = list(
      mean_width = mean(ci_width, na.rm = TRUE),
      sd_width = sd(ci_width, na.rm = TRUE),
      n_obs = nrow(series_clean)
    ),
    recent = last_point
  )))
}

report_var_ci_timeseries <- function(result) {
  series <- result$series
  series_clean <- series[complete.cases(series$VaR, series$CI_low, series$CI_up), ]
  
  # 计算置信区间宽度
  series_clean$CI_width <- series_clean$CI_up - series_clean$CI_low
  
  cat("=== VaR 置信区间时间序列特征 ===\n\n")
  
  # 分位数统计
  cat("置信区间宽度分位数:\n")
  quantiles <- quantile(series_clean$CI_width, 
                        probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), 
                        na.rm = TRUE)
  for (i in seq_along(quantiles)) {
    cat(sprintf("   %s 分位数: %.4f\n", 
                names(quantiles)[i], 
                quantiles[i]))
  }

  # 高波动时期识别
  high_vol_threshold <- quantile(series_clean$CI_width, 0.9, na.rm = TRUE)
  high_vol_periods <- series_clean[series_clean$CI_width > high_vol_threshold, ]
  
  cat(sprintf("\n高不确定性时期 (前10%%): %d 个观测点\n", nrow(high_vol_periods)))
  cat("高不确定性时期示例:\n")
  if (nrow(high_vol_periods) > 0) {
    for (i in 1:min(3, nrow(high_vol_periods))) {
      cat(sprintf("   %s: VaR=%.4f, CI=[%.4f, %.4f], 宽度=%.4f\n",
                  high_vol_periods$date[i], high_vol_periods$VaR[i],
                  high_vol_periods$CI_low[i], high_vol_periods$CI_up[i],
                  high_vol_periods$CI_width[i]))
    }
  }
  
  return(invisible(list(
    quantiles = quantiles,
    high_vol_periods = high_vol_periods
  )))
}

compare_var_ci_methods <- function(results_list) {
  cat("=== 多方法 VaR 置信区间比较 ===\n\n")
  
  comparison <- data.frame()
  
  for (i in seq_along(results_list)) {
    result <- results_list[[i]]
    method_name <- result$name
    
    series_clean <- result$series[complete.cases(
      result$series$VaR, result$series$CI_low, result$series$CI_up), ]
    
    ci_width <- series_clean$CI_up - series_clean$CI_low
    
    method_stats <- data.frame(
      Method = method_name,
      Mean_CI_Width = mean(ci_width, na.rm = TRUE),
      SD_CI_Width = sd(ci_width, na.rm = TRUE),
      Median_CI_Width = median(ci_width, na.rm = TRUE),
      Min_CI_Width = min(ci_width, na.rm = TRUE),
      Max_CI_Width = max(ci_width, na.rm = TRUE),
      N_Observations = nrow(series_clean),
      Coverage_Rate = mean(!is.na(ci_width)),
      stringsAsFactors = FALSE
    )
    
    comparison <- rbind(comparison, method_stats)
  }
  
  # 打印比较表格
  print(comparison, row.names = FALSE)
  
  cat("\n关键发现:\n")
  # 找出最精确的方法（最小平均宽度）
  most_precise <- comparison[which.min(comparison$Mean_CI_Width), "Method"]
  cat(sprintf("• 最精确的方法: %s (平均置信区间宽度最小)\n", most_precise))
  
  # 找出最稳定的方法（最小标准差）
  most_stable <- comparison[which.min(comparison$SD_CI_Width), "Method"]
  cat(sprintf("• 最稳定的方法: %s (置信区间宽度波动最小)\n", most_stable))
  
  return(invisible(comparison))
}

visualize_var_ci <- function(result, title_suffix = "",
                             save = NULL, width = 10, height = 7, dpi = 300,
                             base_family = "Helvetica") {
  series <- result$series
  series_clean <- series[complete.cases(series$VaR, series$CI_low, series$CI_up), ]
  series_clean$date <- as.Date(series_clean$date)
  series_clean$CI_width <- with(series_clean, CI_up - CI_low)
  
  # 配色与风格
  col_var <- "#4C72B0"; col_hist <- "#55A868"; col_scatter <- "#8172B2"; col_mean <- "red"
  
  nature_theme <- theme_minimal(base_size = 13, base_family = base_family) +
    theme(
      panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold"),
      plot.title  = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 6)),
      plot.margin = margin(8, 12, 8, 12),
      legend.position = "none"   # 🚫 彻底去掉所有图例
    )
  
  # 1) VaR + CI
  p1 <- ggplot(series_clean, aes(date)) +
    geom_ribbon(aes(ymin = CI_low, ymax = CI_up), fill = col_var, alpha = 0.20) +
    geom_line(aes(y = VaR), color = col_var, linewidth = 0.9) +
    labs(title = paste0("VaR and 95% CI ", title_suffix), x = "Date", y = "VaR") +
    nature_theme
  
  # 2) CI 宽度分布
  mean_w <- mean(series_clean$CI_width, na.rm = TRUE)
  p2 <- ggplot(series_clean, aes(CI_width)) +
    geom_histogram(bins = 30, fill = col_hist, color = "white", alpha = 0.85) +
    geom_vline(xintercept = mean_w, color = col_mean, linetype = "dashed", linewidth = 0.8) +
    annotate("label", x = mean_w, y = Inf, label = paste0("Mean = ", round(mean_w, 4)),
             vjust = 1.6, color = col_mean, fill = scales::alpha("white", 0.9), label.size = 0) +
    labs(title = "Distribution of CI width", x = "CI width", y = "Count") +
    nature_theme
  
  # 3) CI 宽度时间变化
  p3 <- ggplot(series_clean, aes(date, CI_width)) +
    geom_line(linewidth = 0.8, color = col_hist) +
    geom_hline(yintercept = mean_w, color = col_mean, linetype = "dashed") +
    labs(title = "Temporal change of CI width", x = "Date", y = "CI width") +
    nature_theme
  
  # 4) VaR vs CI width
  corr <- suppressWarnings(cor(series_clean$VaR, series_clean$CI_width, use = "complete.obs"))
  p4 <- ggplot(series_clean, aes(VaR, CI_width)) +
    geom_point(alpha = 0.7, size = 2, color = col_scatter) +
    labs(title = sprintf("VaR vs CI width (corr = %.3f)", corr),
         x = "VaR", y = "CI width") +
    nature_theme
  
  # 组合四图（整齐对齐、无外部图例）
  final_plot <- (p1 + p2) / (p3 + p4)
  
  if (!is.null(save)) {
    ggsave(save, final_plot, width = width, height = height, dpi = dpi, limitsize = FALSE)
  }
  
  print(final_plot)
  
  invisible(list(mean_ci_width = mean_w, correlation = corr))
}


## ---------- 2.1) 统一评价 ----------
evaluate_backtest_performance <- function(results_list) {
  cat("=== 回测性能综合评价 ===\n\n")
  
  performance_df <- data.frame()
  
  for (result in results_list) {
    method <- result$name
    backtest <- result$backtest
    
    # 失败比例检验
    pof_pvalue <- backtest$POF[['p.value']]
    actual_failure_rate <- backtest$POF[['pi']]
    theoretical_failure_rate <- 0.01  # 对于tau=0.99
    
    # 条件覆盖检验
    cc_pvalue <- backtest$IND[['p_cc']]
    
    # 交通灯检验
    tl_green <- backtest$TL$share[["green"]]
    tl_red <- backtest$TL$share[["red"]]
    
    method_stats <- data.frame(
      Method = method,
      POF_Pvalue = pof_pvalue,
      Actual_Failure_Rate = actual_failure_rate,
      Theoretical_Failure_Rate = theoretical_failure_rate,
      CC_Pvalue = cc_pvalue,
      Green_Zone = tl_green,
      Red_Zone = tl_red,
      POF_Pass = pof_pvalue > 0.05,  # 是否通过失败比例检验
      CC_Pass = cc_pvalue > 0.05,     # 是否通过条件覆盖检验
      TL_Pass = tl_red < 0.1,         # 红色区域是否小于10%
      stringsAsFactors = FALSE
    )
    
    performance_df <- rbind(performance_df, method_stats)
  }
  
  # 计算综合得分
  performance_df$Composite_Score <- 
    (performance_df$POF_Pass * 0.4) + 
    (performance_df$CC_Pass * 0.3) + 
    (performance_df$TL_Pass * 0.3) +
    (1 - abs(performance_df$Actual_Failure_Rate - performance_df$Theoretical_Failure_Rate)) * 0.5
  
  print(performance_df, row.names = FALSE)
  return(invisible(performance_df))
}

evaluate_stability <- function(results_list) {
  cat("=== 风险预测稳定性评估 ===\n\n")
  
  stability_df <- data.frame()
  
  for (result in results_list) {
    series <- result$series
    method <- result$name
    
    # 移除缺失值
    series_clean <- series[complete.cases(series$VaR), ]
    
    # VaR序列的波动性
    var_volatility <- sd(series_clean$VaR, na.rm = TRUE)
    
    # VaR变化的波动性 (衡量预测的平滑度)
    var_changes <- diff(series_clean$VaR)
    change_volatility <- sd(var_changes, na.rm = TRUE)
    
    # 最大回撤 (最大连续恶化)
    max_drawdown <- max(cummax(series_clean$VaR) - series_clean$VaR, na.rm = TRUE)
    
    # 转折点比率 (衡量预测的稳定性)
    turning_points <- sum(diff(sign(diff(series_clean$VaR))) != 0, na.rm = TRUE)
    turning_ratio <- turning_points / (nrow(series_clean) - 2)
    
    method_stats <- data.frame(
      Method = method,
      VaR_Volatility = var_volatility,
      Change_Volatility = change_volatility,
      Max_Drawdown = max_drawdown,
      Turning_Ratio = turning_ratio,
      Stability_Score = 1 / (change_volatility + 0.1 * turning_ratio),  # 稳定性得分
      stringsAsFactors = FALSE
    )
    
    stability_df <- rbind(stability_df, method_stats)
  }
  
  print(stability_df, row.names = FALSE)
  return(invisible(stability_df))
}

comprehensive_model_evaluation <- function(results_list, L_next) {
  cat("=== 风险模型综合评价 ===\n\n")
  
  # 获取各项评估结果
  backtest_perf <- evaluate_backtest_performance(results_list)
  stability <- evaluate_stability(results_list)
  extreme_perf <- evaluate_extreme_event_performance(results_list, L_next)
  ci_comparison <- compare_var_ci_methods(results_list)  # 您已有的函数
  
  # 合并所有评估指标
  comprehensive_df <- backtest_perf %>%
    left_join(stability, by = "Method") %>%
    left_join(extreme_perf, by = "Method") %>%
    left_join(ci_comparison, by = "Method")
  
  # 计算最终综合得分 (可根据需求调整权重)
  comprehensive_df$Final_Score <- 
    comprehensive_df$Composite_Score * 0.3 +      # 回测性能 30%
    comprehensive_df$Stability_Score * 0.2 +      # 稳定性 20%
    comprehensive_df$Extreme_Score * 0.25 +       # 极端事件表现 25%
    (1 / comprehensive_df$Mean_CI_Width) * 0.25   # 精度 25%
  
  # 标准化最终得分
  comprehensive_df$Final_Score_Normalized <- 
    comprehensive_df$Final_Score / max(comprehensive_df$Final_Score)
  
  # 按最终得分排序
  comprehensive_df <- comprehensive_df[order(-comprehensive_df$Final_Score_Normalized), ]
  
  cat("\n最终排名:\n")
  print(comprehensive_df[, c("Method", "Final_Score_Normalized", 
                             "POF_Pass", "CC_Pass", "TL_Pass",
                             "Mean_CI_Width", "Stability_Score", "Extreme_Score")], 
        row.names = FALSE)
  
  return(invisible(comprehensive_df))
}

# results_list: list(res.HS, res.FHS, ...)
# L_next: 用 make_L_next(...) 得到的长度与各方法相同的损失序列
# k_top: 取前多少个极端损失
# scoring: "ratio" 或 "penalty"（见下）
evaluate_extreme_event_performance <- function(results_list,
                                               L_next,
                                               crisis_periods = NULL,
                                               k_top = 10,
                                               scoring = c("ratio","penalty")) {
  cat("=== 极端事件表现评估 ===\n\n")
  scoring <- match.arg(scoring)
  extreme_df <- data.frame(stringsAsFactors = FALSE)
  
  # 可选：若传入危机期掩码（长度与 L_next 相同的 TRUE/FALSE 向量），只在危机期内评估
  if (!is.null(crisis_periods)) {
    if (is.logical(crisis_periods) && length(crisis_periods) == length(L_next)) {
      L_mask <- crisis_periods
    } else {
      warning("crisis_periods 忽略：需为与 L_next 等长的逻辑向量。")
      L_mask <- rep(TRUE, length(L_next))
    }
  } else {
    L_mask <- rep(TRUE, length(L_next))
  }
  
  L_next <- as.numeric(L_next)
  
  for (result in results_list) {
    method <- result$name
    series <- result$series
    
    if (!("VaR" %in% names(series))) {
      warning(sprintf("Method '%s' 缺少 VaR 列，跳过。", method))
      next
    }
    
    VaR <- as.numeric(series$VaR)
    
    # 1) 长度对齐
    n <- min(length(L_next), length(VaR))
    Ln <- L_next[seq_len(n)]
    Vr <- VaR   [seq_len(n)]
    Mk <- L_mask[seq_len(n)]
    
    # 2) 有效索引：Ln/Vr 同时有限，且在掩码范围内
    valid_idx <- which(is.finite(Ln) & is.finite(Vr) & Mk)
    if (length(valid_idx) == 0) {
      extreme_df <- rbind(extreme_df, data.frame(
        Method = method,
        Extreme_Failure_Rate = NA_real_,
        Conservative_Gap     = NA_real_,
        Conservativeness_Ratio = NA_real_,
        Extreme_Score        = NA_real_,
        K_Top                = k_top,
        stringsAsFactors = FALSE
      ))
      next
    }
    
    # 3) 在有效索引里选前 k_top 个最大损失（索引排序，避免 match/浮点问题）
    ord  <- order(Ln[valid_idx], decreasing = TRUE)
    take <- min(k_top, length(ord))
    top_idx <- valid_idx[ord[seq_len(take)]]
    
    # 4) 指标计算
    ex_top <- as.integer(Ln[top_idx] > Vr[top_idx])                # 是否突破
    extreme_failure_rate <- mean(ex_top, na.rm = TRUE)             # 失败率
    gaps <- Ln[top_idx] - Vr[top_idx]                              # 损失- VaR（负值=过保守）
    conservative_gap <- mean(gaps, na.rm = TRUE)
    
    # 5) 评分（截断到 [0,1]）
    if (scoring == "ratio") {
      # A: 简单稳健 —— VaR/Loss，>1 截断成 1
      ratio <- Vr[top_idx] / pmax(Ln[top_idx], 1e-12)
      comp_part <- mean(pmin(ratio, 1), na.rm = TRUE)              # ∈ [0,1]
    } else {
      # B: 不同惩罚 —— 低估重惩罚、过保守轻惩罚（权重可调）
      over  <- pmax(Vr[top_idx] - Ln[top_idx], 0) / pmax(Ln[top_idx], 1e-12) # 过保守比例
      under <- pmax(Ln[top_idx] - Vr[top_idx], 0) / pmax(Ln[top_idx], 1e-12) # 低估比例
      comp_part <- 1 - mean(under, na.rm = TRUE) - 0.25 * mean(over, na.rm = TRUE)
      comp_part <- max(0, min(comp_part, 1))
    }
    
    # 终分：越不失败越高 + 越接近不过度/不低估越高
    extreme_score <- (1 - extreme_failure_rate) * 0.7 + comp_part * 0.3
    extreme_score <- max(0, min(extreme_score, 1))
    
    # 记录一个可解释的“保守度比例”输出（便于横向比较）
    conservativeness_ratio <- mean(pmin(Vr[top_idx] / pmax(Ln[top_idx], 1e-12), 1), na.rm = TRUE)
    
    extreme_df <- rbind(extreme_df, data.frame(
      Method = method,
      Extreme_Failure_Rate = extreme_failure_rate,
      Conservative_Gap     = conservative_gap,
      Conservativeness_Ratio = conservativeness_ratio,
      Extreme_Score        = extreme_score,
      K_Top                = take,
      stringsAsFactors = FALSE
    ))
  }
  
  print(extreme_df, row.names = FALSE)
  invisible(extreme_df)
}


# 逐危机期循环：在每个掩码内仅选取该期内的极端损失进行评估
evaluate_extremes_by_period <- function(results_list,
                                        L_next,
                                        aligned_dates,
                                        crisis_periods = default_crisis_periods,
                                        k_top = 10,
                                        scoring = c("ratio","penalty")) {
  scoring <- match.arg(scoring)
  masks <- build_crisis_masks(aligned_dates, crisis_periods)
  out_all <- data.frame(stringsAsFactors = FALSE)
  
  for (nm in names(masks)) {
    cat(sprintf("\n=== 危机期：%s ===\n", nm))
    msk <- masks[[nm]]
    # 若该期内有效样本过少，跳过
    if (sum(msk, na.rm = TRUE) < 5) {
      cat("  （样本过少，跳过）\n")
      next
    }
    tab <- evaluate_extreme_event_performance(
      results_list,
      L_next = L_next,
      crisis_periods = msk,
      k_top = k_top,
      scoring = scoring
    )
    tab$Period <- nm
    out_all <- rbind(out_all, tab[, c("Period", names(tab))])
  }
  # 调整列顺序便于阅读
  cols <- c("Period","Method","Extreme_Failure_Rate","Conservative_Gap",
            "Conservativeness_Ratio","Extreme_Score","K_Top")
  cols <- cols[cols %in% names(out_all)]
  out_all <- out_all[, cols, drop = FALSE]
  rownames(out_all) <- NULL
  out_all
}


## ---------- 3) 数据获取（FRED S&P500，简版） ----------
fetch_spx_fred <- function(){
  url <- "https://fred.stlouisfed.org/graph/fredgraph.csv?id=SP500"
  dat <- tryCatch(utils::read.csv(url, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(dat)) stop("FRED 下载失败")
  names(dat) <- c("date","SPX")
  dat$date <- as.Date(dat$date); dat <- dat[order(dat$date), ]
  dat$SPX <- as.numeric(dat$SPX); dat <- dat[is.finite(dat$SPX), ]
  r <- diff(log(dat$SPX))
  list(date = dat$date[-1], r = as.numeric(r))
}

## ---------- 4) 示例运行 ----------
spx   <- fetch_spx_fred()
tau   <- 0.99
nwin  <- 250
alpha <- 0.05
lambda<- 0.95
res   <- run_methods(spx$r, spx$date, tau=tau, nwin=nwin, alpha=alpha,
                     lambda=lambda, df=5, evt_uq=0.9,
                     methods = c("HS","FHS","Delta","EVT","MCt","QR","GARCH"),
                     B_ci = 1000, cores_ci = 6)

options(scipen = 10)

results_list = list(res.HS, res.FHS, res.MCt, res.QR, res.EVT, res.GARCH, res.EWHS)

compare_var_ci_methods(results_list)

evaluate_backtest_performance(results_list)

evaluate_stability(results_list)

# # （1）先全局构造 L_next（只做一次）
# L_next <- make_L_next(spx$r, nwin)
# 
# 
# # （3）直接评估 —— 推荐 ratio 评分，取前 10 个极端损失
# extreme_tab <- evaluate_extreme_event_performance(
#   results_list,
#   L_next = L_next,
#   k_top = 10,
#   scoring = "ratio"      # 或 "penalty"
# )

# 典型危机/高波动期（S&P 500 视角）
default_crisis_periods <- list(
  `2015 China Deval & Black Monday` = c("2015-08-17", "2015-09-04"),
  `2016 Early Selloff`              = c("2016-01-04", "2016-02-11"),
  `2016 Brexit`                     = c("2016-06-20", "2016-06-30"),
  `2018 Volmageddon`                = c("2018-02-01", "2018-02-15"),
  `2018 Q4 Drawdown`                = c("2018-10-01", "2018-12-24"),
  `2020 COVID Crash`                = c("2020-02-19", "2020-06-23"),
  `2020 US Election Vol`            = c("2020-10-26", "2020-11-06"),
  `2022 Fed Tightening Bear`        = c("2022-01-03", "2022-10-14"),
  `2023 US Regional Banks`          = c("2023-03-08", "2023-03-31")
)


## 1) 对齐构造（只需一次）
L_next        <- make_L_next(spx$r, nwin)
aligned_dates <- make_aligned_dates(spx$date, nwin)

## 2) 模型结果列表（你已有）
results_list  <- list(res.HS, res.FHS, res.MCt, res.QR, res.EVT, res.GARCH, res.EWHS)

## 3) 使用默认危机期，评分用 "ratio"，每期取前 10 个极端损失
by_period <- evaluate_extremes_by_period(
  results_list,
  L_next        = L_next,
  aligned_dates = aligned_dates,
  crisis_periods = default_crisis_periods,  # 也可传你自定义的列表
  k_top = 10,
  scoring = "ratio"                          # 或 "penalty"
)

## 4) 查看结果（按 Period, Extreme_Score 由高到低排序）
by_period[order(by_period$Period, -by_period$Extreme_Score), ]


# 每个危机期取 Extreme_Score 最大的模型
best_by_period <- by_period %>%
  group_by(Period) %>%
  slice_max(Extreme_Score, n = 1, with_ties = FALSE) %>%
  arrange(Period)

best_by_period



ggplot(by_period, aes(x = Period, y = Method, fill = Extreme_Score)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "#3182bd") +
  labs(title = "Extreme Event Performance by Crisis Period",
       x = "Crisis Period", y = "Method", fill = "Extreme Score") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(by_period, aes(x = Period, y = Method, fill = Extreme_Score)) +
  geom_tile(color = "grey80") +
  scale_fill_gradient(low = "white", high = "#08519c", limits = c(0,1)) +
  labs(title = "Extreme Event Performance by Crisis Period",
       x = NULL, y = NULL, fill = "Extreme Score") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




ggplot(results_list$HS$series, aes(date)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_up), fill = col_var, alpha = 0.20) +
  geom_line(aes(y = VaR), color = col_var, linewidth = 0.9) +
  
  geom_ribbon(data = results_list$FHS$series, aes(ymin = CI_low, ymax = CI_up), fill = col_var, alpha = 0.20) +
  geom_line(data = results_list$FHS$series, aes(x = date, y = VaR), color = 2, linewidth = 0.9) +
  labs(title = paste0("VaR and 95% CI ", 'HS'), x = "Date", y = "VaR") + 
    
    
  geom_ribbon(data = results_list$Delta$series, aes(ymin = CI_low, ymax = CI_up), fill = col_var, alpha = 0.20) +
  geom_line(data = results_list$Delta$series, aes(x = date, y = VaR), color = 3, linewidth = 0.9) +
  labs(title = paste0("VaR and 95% CI ", 'HS'), x = "Date", y = "VaR") +
    
  geom_ribbon(data = results_list$MCt$series, aes(ymin = CI_low, ymax = CI_up), fill = col_var, alpha = 0.20) +
  geom_line(data = results_list$MCt$series, aes(x = date, y = VaR), color = 4, linewidth = 0.9) +
  labs(title = paste0("VaR and 95% CI ", 'HS'), x = "Date", y = "VaR") + 
  
  geom_ribbon(data = results_list$QR$series, aes(ymin = CI_low, ymax = CI_up), fill = col_var, alpha = 0.20) +
  geom_line(data = results_list$QR$series, aes(x = date, y = VaR), color = 5, linewidth = 0.9) +
  labs(title = paste0("VaR and 95% CI ", 'HS'), x = "Date", y = "VaR") + 
  
  nature_theme


