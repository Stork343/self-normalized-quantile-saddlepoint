## =========================================================
## SN-Q-ESA + 多基线分位数区间：完整自包含实现（修正版）
## 依赖：base R；可选：pbmcapply（并行），quantileCI（HS/Nyblom）
## =========================================================
library(rlang)
library(dplyr)
library(knitr)
library(kableExtra)

summary_to_latex <- function(res,
                             caption = "Method comparison on simulated data",
                             label   = "tab:summary",
                             drop_cols = c("fail_rate","qhat_in_rate","mean_true_relpos",
                                           "lower_miss_rate","upper_miss_rate","mean_pinball"),
                             digits = 4,
                             file = NULL) {
  stopifnot(is.list(res), !is.null(res$summary))
  df <- res$summary %>%
    select(-any_of(drop_cols)) %>%
    # 去掉全是 NA 的方法（如你截图里 SNQESA_rand）
    filter(rowSums(across(where(is.numeric), ~ is.finite(.))) > 0)
  
  # 数值列统一格式
  num_cols <- names(df)[sapply(df, is.numeric)]
  df[num_cols] <- lapply(df[num_cols], function(x) ifelse(is.finite(x), round(x, digits), NA))
  
  # 更友好的列名（可按需修改或删掉）
  col_names <- c(
    method = "Method",
    coverage = "Coverage",
    cover_se = "SE(Cover)",
    mean_len = "Mean length",
    median_len = "Median length",
    mean_time = "Mean time (s)",
    mean_center_bias = "Center bias (mean)",
    median_center_bias = "Center bias (median)",
    rmse_center_bias = "Center bias (RMSE)",
    mean_IS = "Interval Score (mean)",
    median_IS = "Interval Score (median)",
    mean_pinball = "Pinball Loss(mean)"
  )
  col_names <- col_names[names(df)] %||% names(df)
  
  kb <- kbl(df, format = "latex", booktabs = TRUE, longtable = TRUE,
            caption = caption, label = label, col.names = col_names,
            align = c("l", rep("r", ncol(df)-1))) |>
    kable_styling(latex_options = c("hold_position"))
  
  if (!is.null(file)) {
    save_kable(kb, file)
  }
  kb
}
# ===== Oracle 上界：把 SNQESA_cont 端点与 Exact 端点做凸组合，达到目标覆盖 =====
# 思想：对每次重采样 i，把区间端点 L(λ)=L_S + λ(L_E-L_S), U(λ)=U_S + λ(U_E-U_S)，
#      找到最小 λ_i 使得真分位 θ 落入 [L(λ),U(λ)]；再选一个全局 λ*，
#      使总体覆盖率达到 target_coverage（或尽可能接近）。
# 输入：res —— 你的 simulate_compare_methods() 返回对象
# 输出：一行“方法摘要”，字段名与 res$summary 对齐（外加 lambda_star）
snqesa_oracle_row <- function(res,
                              target_coverage = 0.95,
                              baseline = "SNQESA_cont",
                              anchor   = "Exact",
                              count_anchor_time = FALSE  # 若 TRUE，时间=SNQESA时间+Exact时间；否则等同 SNQESA
){
  stopifnot(is.list(res), !is.null(res$raw[[baseline]]), !is.null(res$raw[[anchor]]))
  B <- res$raw[[baseline]]  # baseline: SNQESA_cont
  A <- res$raw[[anchor]]    # anchor: Exact
  R <- nrow(B)
  theta <- res$t_true
  alpha <- res$params$alpha
  
  # 每次重采样需要的最小 λ_i，使 θ ∈ [L(λ),U(λ)]
  lam_needed <- rep(NA_real_, R)
  cover0 <- as.numeric(B[,"cover"])
  for (i in seq_len(R)){
    Lb <- B[i,"lower"]; Ub <- B[i,"upper"]; Le <- A[i,"lower"]; Ue <- A[i,"upper"]
    if (!is.finite(theta) || !is.finite(Lb) || !is.finite(Ub) || !is.finite(Le) || !is.finite(Ue)){
      lam_needed[i] <- Inf; next
    }
    if (cover0[i] == 1) { lam_needed[i] <- 0; next }  # 已覆盖无需校准
    
    # 线性不等式：L(λ) ≤ θ ≤ U(λ)，λ ∈ [0,1]
    lo <- 0; hi <- 1; feasible <- TRUE
    
    dL <- Le - Lb
    if (abs(dL) < 1e-12){
      if (!(Lb <= theta)) feasible <- FALSE
    } else if (dL > 0){
      hi <- min(hi, (theta - Lb)/dL)
    } else { # dL < 0
      lo <- max(lo, (theta - Lb)/dL)
    }
    
    dU <- Ue - Ub
    if (abs(dU) < 1e-12){
      if (!(theta <= Ub)) feasible <- FALSE
    } else if (dU > 0){
      lo <- max(lo, (theta - Ub)/dU)
    } else { # dU < 0
      hi <- min(hi, (theta - Ub)/dU)
    }
    
    if (!feasible) { lam_needed[i] <- Inf; next }
    lo <- max(lo, 0); hi <- min(hi, 1)
    lam_needed[i] <- if (lo <= hi) max(0, lo) else Inf
  }
  
  # 选择全局 λ*：让覆盖率 ≥ target_coverage
  k0 <- sum(cover0 == 1, na.rm = TRUE)
  if (k0 / R >= target_coverage){
    lam_star <- 0
  } else {
    req <- ceiling(target_coverage * R - k0)
    pool <- sort(lam_needed[is.finite(lam_needed) & lam_needed > 0])
    lam_star <- if (length(pool) >= req) pool[req] else 1
  }
  
  # 用 λ* 混合端点，计算新的指标
  L <- B[,"lower"] + lam_star * (A[,"lower"] - B[,"lower"])
  U <- B[,"upper"] + lam_star * (A[,"upper"] - B[,"upper"])
  
  cov_vec <- (theta >= L & theta <= U)
  coverage <- mean(cov_vec, na.rm = TRUE)
  cover_se <- sqrt(coverage * (1 - coverage) / sum(is.finite(L) & is.finite(U)))
  
  length_vec <- ifelse(is.finite(L) & is.finite(U), U - L, NA)
  center     <- (L + U) / 2
  bias       <- center - theta
  
  # Interval Score（与你框架一致）
  interval_score <- function(Li, Ui, theta, alpha){
    if (!is.finite(Li) || !is.finite(Ui)) return(NA_real_)
    w <- Ui - Li
    pen <- if (theta < Li) (2/alpha)*(Li - theta) else if (theta > Ui) (2/alpha)*(theta - Ui) else 0
    w + pen
  }
  is_vec <- mapply(interval_score, L, U, MoreArgs = list(theta = theta, alpha = alpha))
  
  time_vec <- if (isTRUE(count_anchor_time)) B[,"time"] + A[,"time"] else B[,"time"]
  
  c(
    method              = "SNQESA_oracle",
    coverage            = coverage,
    cover_se            = cover_se,
    mean_len            = mean(length_vec,  na.rm = TRUE),
    median_len          = median(length_vec,na.rm = TRUE),
    mean_time           = mean(time_vec,    na.rm = TRUE),
    mean_center_bias    = mean(bias,        na.rm = TRUE),
    median_center_bias  = median(bias,      na.rm = TRUE),
    rmse_center_bias    = sqrt(mean(bias^2, na.rm = TRUE)),
    mean_IS             = mean(is_vec,      na.rm = TRUE),
    median_IS           = median(is_vec,    na.rm = TRUE),
    lambda_star         = lam_star
  )
}

# 把上面的“一行结果”渲染成你当前表格的一行 LaTeX（符号表头样式，含高亮）
oracle_to_latex_row <- function(row_named_vec,
                                method_label = "SNQESA$^{\\star}$",
                                highlight = TRUE,
                                digits = 4){
  f <- function(x){ y <- as.numeric(x); if (!is.finite(y)) return(""); sprintf(paste0("%.",digits,"f"), y) }
  h <- if (isTRUE(highlight)) "\\rowcolor{highlightcolor}\n" else ""
  paste0(
    h,
    method_label, " & ",
    f(row_named_vec["coverage"]), " & ",
    f(row_named_vec["cover_se"]), " & ",
    f(row_named_vec["mean_len"]), " & ",
    f(row_named_vec["median_len"]), " & ",
    f(row_named_vec["mean_time"]), " & ",
    f(row_named_vec["mean_center_bias"]), " & ",
    f(row_named_vec["median_center_bias"]), " & ",
    f(row_named_vec["rmse_center_bias"]), " & ",
    f(row_named_vec["mean_IS"]), " & ",
    f(row_named_vec["median_IS"]),
    "\\\\"
  )
}


## ---------- 0. Helpers ----------
.Phi  <- function(z) pnorm(z)
.phi  <- function(z) dnorm(z)
.clip01 <- function(u, eps = 1e-12){ pmin(1 - eps, pmax(eps, u)) }
.logit <- function(u){ u <- .clip01(u); log(u/(1 - u)) }
.kl_binom <- function(u, p){
  u <- .clip01(u); p <- .clip01(p)
  u*log(u/p) + (1-u)*log((1-u)/(1-p))
}
#' @title .normalize_dist
#' @description 内部辅助函数，用于标准化分布名称（例如 t(3) -> t3）
.normalize_dist <- function(dist_name) {
  d <- tolower(trimws(dist_name))
  d <- gsub("[ ()]", "", d) # 移除括号和空格
  
  # 别名替换
  if (d == "log-normal") return("lognormal")
  if (d == "mix") return("mixnorm")
  if (d == "mixednormal") return("mixnorm")
  if (d == "t(2)") return("t2")
  if (d == "t(3)") return("t3")
  
  return(d)
}


## ---------- 1. Sampling & true quantiles ----------

#' @title rsamp
#' @description 从指定的分布生成随机样本
#' @param n 样本量
#' @param dist 分布名称
rsamp <- function(n, dist = c("normal","lognormal","t2","t3","cauchy","mixnorm",
                              "beta", "exp")){
  
  # 确保别名 (如 "lognormal" 和 "t(3)") 被正确解析
  dist <- .normalize_dist(match.arg(.normalize_dist(dist),
                                    c("normal","lognormal","t2","t3","cauchy","mixnorm",
                                      "beta", "exp")))
  
  switch(dist,
         normal    = rnorm(n),
         lognormal = rlnorm(n),
         t2        = rt(n, df = 2),
         t3        = rt(n, df = 3),
         cauchy    = rcauchy(n, 0, 2),
         mixnorm   = { z <- rnorm(n); mu <- sample(c(-1,1), n, TRUE); z + mu },
         beta      = rbeta(n, 0.5, 0.5),
         exp       = rexp(n, rate = 1)
  )
}

pmixnorm <- function(x){
  0.5*pnorm(x, mean = -1, sd = 1) + 0.5*pnorm(x, mean = 1, sd = 1)
}

#' @title qmixnorm
#' @description 计算混合正态分布的真实分位数 (通过数值求解)
#' N(0,1) + U(-1,1) -> 0.5*N(-1,1) + 0.5*N(1,1)
qmixnorm <- function(tau, mu1 = -1, mu2 = 1, sig = 1, w1 = 0.5, w2 = 0.5){
  
  # 混合分布的累积分布函数 (CDF)
  pmix <- function(q) {
    w1 * pnorm(q, mu1, sig) + w2 * pnorm(q, mu2, sig)
  }
  
  # 目标函数：pmix(q) - tau = 0
  f_root <- function(q) {
    pmix(q) - tau
  }
  
  # 通过 uniroot 数值求解
  # 我们需要一个包含解的合理区间 [lower, upper]
  lower <- min(mu1, mu2) - 10 * sig
  upper <- max(mu1, mu2) + 10 * sig
  
  # 检查边界
  if (f_root(lower) > 0) lower <- lower - 20 * sig
  if (f_root(upper) < 0) upper <- upper + 20 * sig
  
  # 求解
  sol <- tryCatch(
    uniroot(f_root, interval = c(lower, upper), tol = 1e-12),
    error = function(e) {
      # 如果失败，扩大搜索范围
      uniroot(f_root, interval = c(lower - 100, upper + 100), tol = 1e-12)
    }
  )
  
  return(sol$root)
}

#' @title true_q
#' @description 计算指定分布的真实分位数
#' @param dist 分布名称
#' @param tau 目标分位数 (0 < tau < 1)
true_q <- function(dist = c("normal","lognormal","t2","t3","cauchy","mixnorm",
                            "beta", "exp"), tau){
  
  dist <- .normalize_dist(match.arg(.normalize_dist(dist),
                                    c("normal","lognormal","t2","t3","cauchy","mixnorm",
                                      "beta_u", "exp_left")))
  
  switch(dist,
         normal    = qnorm(tau),
         lognormal = qlnorm(tau),
         t2        = qt(tau, df = 2),
         t3        = qt(tau, df = 3),
         cauchy    = qcauchy(tau, 0, 2),
         mixnorm   = qmixnorm(tau), # 需要辅助函数
         
         # --- 新增分布的真实分位数 ---
         beta_u    = qbeta(tau, 0.5, 0.5),
         exp_left  = qexp(tau, rate = 1)
  )
}




## ---------- 2. h^{-1}：u_x from |x| ----------
snqesa_hinv <- function(x_abs, tau, n, tol = 1e-12){
  stopifnot(x_abs >= 0)
  a <- x_abs^2
  A <- n
  B <- -2*n*tau - a*(1 - 2*tau)
  C <- tau^2*(n - a)
  disc <- B*B - 4*A*C
  if (disc < 0) disc <- 0
  rdisc <- sqrt(disc)
  u1 <- (-B + rdisc) / (2*A)
  u2 <- (-B - rdisc) / (2*A)
  cand <- c(u1, u2)
  cand <- cand[is.finite(cand) & cand > 0 & cand < 1]
  if (length(cand) == 0){
    h <- function(u){ sqrt(n)*(tau - u)/sqrt(tau^2 + u*(1 - 2*tau)) }
    f <- function(u) h(u) - x_abs
    eps <- 1e-10
    lo <- eps; hi <- max(eps*10, min(0.999999, tau - eps))
    if (f(lo)*f(hi) > 0){
      # 指数扩展直到变号或达界
      for (k in 1:500){
        lo <- max(eps, lo/2); hi <- min(tau - eps, hi * 1.5)
        if (f(lo)*f(hi) <= 0) break
      }
      if (f(lo)*f(hi) > 0) return(.clip01(tau - abs(x_abs)/sqrt(n), 1e-9))
    }
    return(uniroot(f, c(lo, hi), tol = tol)$root)
  }
  below <- cand[cand <= tau + 1e-12]
  if (length(below) > 0) return(max(below))
  cand[which.min(abs(cand - tau))]
}

## ---------- 3. LR / r* tail on Bin(n,tau) ----------
.lr_tail_binom <- function(u_x, tau, n, use_rstar = TRUE){
  u_x <- .clip01(u_x)
  r_sign <- ifelse(u_x >= tau, 1, -1)
  r <- r_sign * sqrt(2 * n * .kl_binom(u_x, tau))
  qLR <- abs(.logit(u_x) - .logit(tau)) * sqrt(n * u_x * (1 - u_x))
  if (!is.finite(r) || !is.finite(qLR) || qLR <= 0) {
    z <- (u_x - tau) / sqrt(tau * (1 - tau) / n)
    return(list(tail = 1 - .Phi(abs(z)), r = r, qLR = qLR, rstar = NA, u = u_x))
  }
  q_signed <- ifelse(r >= 0, qLR, -qLR)
  r_eps <- 1e-8
  if (use_rstar){
    ratio <- r / q_signed
    if (is.finite(ratio) && ratio > 0 && abs(r) > r_eps){
      rstar <- r + (1/r) * log(ratio)
      if (is.finite(rstar)) return(list(tail = .Phi(-abs(rstar)), r = r, qLR = qLR, rstar = rstar, u = u_x))
    }
  }
  #tail <- .Phi(-abs(r)) + .phi(abs(r)) * (1 / max(abs(r), 1e-12) - 1 / qLR)
  tail <- .Phi(-abs(r)) + .phi(abs(r)) * (1 / max(abs(r), 1e-8) - 1 / max(qLR, 1e-8))
  tail <- max(0, min(1, tail))
  list(tail = tail, r = r, qLR = qLR, rstar = NA, u = u_x)
}

## ---------- 4. 2D constrained route（Skovgaard-style via p） ----------
snqesa_tail_2d <- function(x_abs, tau, n){
  f <- function(p){ n * (tau - p) - x_abs * sqrt(n * (tau^2 + p * (1 - 2 * tau))) }
  p0 <- snqesa_hinv(x_abs, tau, n)
  lo <- max(1e-12, min(p0, tau) - 0.4)
  hi <- min(1 - 1e-12, max(p0, tau) + 0.4)
  if (f(lo) * f(hi) > 0){ lo <- 1e-12; hi <- 1 - 1e-12 }
  p <- tryCatch(uniroot(f, c(lo, hi), tol = 1e-12)$root, error = function(e) p0)
  .lr_tail_binom(p, tau, n, use_rstar = TRUE)
}

## ---------- 5. 单点 p 值 ----------
snqesa_pvalue <- function(x, t, tau, ridge = NULL, engine = c("2d","1d"), calib = c("LR","Exact","midp")){
  stopifnot(length(tau) == 1, tau > 0, tau < 1)
  n <- length(x)
  engine <- match.arg(engine)
  calib  <- match.arg(calib)
  if (is.null(ridge)) ridge <- max(n^(-3/4), 2/n)  # 例如 2/n 作为轻度下
  k <- sum(x <= t)
  S <- n * (tau - k / n)
  Q <- k * (1 - tau)^2 + (n - k) * tau^2
  xobs <- S / sqrt(Q + ridge)
  if (!is.finite(xobs)){
    return(list(p = NA_real_, p_two_sided = NA_real_, p_one_sided = NA_real_,
                xobs = xobs, u_x = NA_real_, detail = NULL, engine = engine, calib = calib))
  }
  out <- if (engine == "2d") snqesa_tail_2d(abs(xobs), tau, n) else {
    u_x <- snqesa_hinv(abs(xobs), tau, n)
    .lr_tail_binom(u_x, tau, n, use_rstar = TRUE)
  }
  if (calib != "LR"){
    # 关键：用观测实际一侧来取单侧尾，然后再做双侧合成
    k_obs <- sum(x <= t)
    if (calib == "Exact"){
      pL <- pbinom(k_obs, size = n, prob = tau)                          # 左尾 P(K ≤ k)
      pR <- 1 - pbinom(k_obs - 1L, size = n, prob = tau)                  # 右尾 P(K ≥ k)
    } else if (calib == "midp"){
      pL <- pbinom(k_obs - 1L, size = n, prob = tau) + 0.5 * dbinom(k_obs, size = n, prob = tau)
      pR <- 1 - (pbinom(k_obs, size = n, prob = tau) - 0.5 * dbinom(k_obs, size = n, prob = tau))
    }
    p_one <- if (k_obs <= n * tau) pL else pR
    out$tail <- min(1, max(0, p_one))
  }
  p_one <- out$tail
  p_two <- min(1, 2 * p_one)
  list(p = p_two, p_two_sided = p_two, p_one_sided = p_one, xobs = xobs,
       u_x = if (!is.null(out$u)) out$u else NA_real_, detail = out, engine = engine, calib = calib)
}

snqesa_pvalue_hybrid <- function(x, t, tau, ridge = NULL, engine = "1d",
                                 r0 = 2.0, gamma = 3.0){
  # r*/LR 路径（你的原函数）
  out_lr <- snqesa_pvalue(x, t, tau, ridge = ridge, engine = engine, calib = "LR")
  rstar  <- if (is.finite(out_lr$detail$rstar)) out_lr$detail$rstar else out_lr$detail$r
  a <- abs(rstar); w <- 1/(1 + exp(-gamma * (a - r0)))  # 0→r*，1→mid-p
  
  # 连续化 mid-p 路径
  cm <- .snqesa_p_one_cont_midp(x, t, tau)
  p1_mid <- if (sum(x <= t) <= length(x) * tau) cm$pL else cm$pR
  
  # 组装单侧与双侧
  p1 <- (1 - w) * out_lr$p_one_sided + w * p1_mid
  p2 <- min(1, 2 * p1)
  list(p_one_sided = p1, p_two_sided = p2, w = w, rstar = rstar)
}

## ---------- 6. 离散 CI（序统计 / 二项反演） ----------
snqesa_ci_discrete <- function(x, tau = 0.9, alpha = 0.05,
                               tail_split = c("equal","minlength")){
  x <- sort(as.numeric(x))
  n <- length(x)
  tail_split <- match.arg(tail_split)
  
  pick_interval <- function(alphaL){
    alphaL <- max(0, min(alpha, alphaL))
    alphaR <- alpha - alphaL
    Lk <- qbinom(alphaL, n, tau)
    Uk <- qbinom(1 - alphaR, n, tau) + 1L   
    
    Lk <- max(1L, min(n,    as.integer(Lk)))
    Uk <- max(1L, min(n+1L, as.integer(Uk)))  
    
    lower_val <- x[Lk]
    upper_val <- if (Uk == n+1L) Inf else x[Uk]  # 上端不截断 → +Inf
    
    list(lower = lower_val, upper = upper_val, L = Lk, U = Uk,
         aL = alphaL, aR = alphaR)
  }
  
  if (tail_split == "equal"){
    res <- pick_interval(alpha/2)
  } else {
    grid <- seq(0, alpha, length.out = max(31, round(10 + 4/alpha)))
    best <- NULL; best_len <- Inf
    for (aL in grid){
      r <- pick_interval(aL)
      len <- r$upper - r$lower
      if (is.finite(len) && len < best_len){ best_len <- len; best <- r }
    }
    res <- if (is.null(best)) pick_interval(alpha/2) else best
  }
  c(res, list(tau = tau, alpha = alpha, type = "discrete"))
}

exact_order_ci <- function(x, tau = 0.9, alpha = 0.05){
  x <- sort(as.numeric(x)); n <- length(x)
  alphaL <- alpha / 2; alphaR <- alpha - alphaL
  Lk <- qbinom(alphaL, n, tau)
  Uk <- qbinom(1 - alphaR, n, tau) + 1L
  
  Lk <- max(1L, min(n,    as.integer(Lk)))
  Uk <- max(1L, min(n+1L, as.integer(Uk)))
  
  lower_val <- x[Lk]
  upper_val <- if (Uk == n+1L) Inf else x[Uk]
  c(lower = lower_val, upper = upper_val, L = Lk, U = Uk)
}


## ---------- 7. 连续反演 CI（参考实现） ----------
.snqesa_find_bracket_left <- function(pfun, hi, alpha, lo_hint){
  f <- function(t) pfun(t) - alpha
  lo <- lo_hint; valL <- f(lo); valH <- f(hi)
  step <- max((hi - lo)/20, .Machine$double.eps)
  for (i in 1:200){
    if (valL * valH <= 0) return(c(lo, hi))
    lo <- lo - step; valL <- f(lo)
    step <- step * 1.1  # 逐步扩大
  }
  c(NA_real_, NA_real_)
}
.snqesa_find_bracket_right <- function(pfun, lo, alpha, hi_hint){
  f <- function(t) pfun(t) - alpha
  hi <- hi_hint; valL <- f(lo); valH <- f(hi)
  step <- max((hi - lo)/20, .Machine$double.eps)
  for (i in 1:200){
    if (valL * valH <= 0) return(c(lo, hi))
    hi <- hi + step; valH <- f(hi)
    step <- step * 1.1
  }
  c(NA_real_, NA_real_)
}
.snqesa_p_one_cont_midp <- function(x, t, tau){
  n <- length(x)
  k <- sum(x <= t)
  xs <- sort(x); k <- max(0L, min(n, k))
  # 计算秩间权重 ω
  if (k == 0) { omega <- 0 } else if (k == n) { omega <- 1 } else {
    denom <- xs[k+1] - xs[k]; if (!is.finite(denom) || denom <= 0) denom <- .Machine$double.eps
    omega <- (t - xs[k]) / denom; omega <- max(0, min(1, omega))
  }
  # 左/右连续化 mid-p
  pL <- pbinom(k-1, n, tau) + omega * dbinom(k, n, tau)
  pR <- pbinom(k-1, n, tau, lower.tail = FALSE) + (1-omega) * dbinom(k, n, tau)
  list(pL = pL, pR = pR)
}
.cont_midp_sides <- function(t, x, tau) {
  n <- length(x)
  k <- sum(x <= t)
  xs <- sort(x); k <- max(0L, min(n, k))
  
  if (k == 0) { omega <- 0 } 
  else if (k == n) { omega <- 1 } 
  else {
    denom <- xs[k+1] - xs[k]
    if (!is.finite(denom) || denom <= 0) denom <- .Machine$double.eps
    omega <- (t - xs[k]) / denom
    omega <- max(0, min(1, omega))
  }
  
  pL <- pbinom(k-1, n, tau) + omega * dbinom(k, n, tau)
  pR <- pbinom(k-1, n, tau, lower.tail = FALSE) + (1-omega) * dbinom(k, n, tau)
  c(pL = pL, pR = pR)
}



#' 计算SN-Q-ESA分位数置信区间
#'
#' @param x 数值向量，输入数据
#' @param tau 分位数水平 (0,1)
#' @param alpha 显著性水平 (0,1)  
#' @param ridge 正则化参数，默认为 n^(-3/4)
#' @param engine 计算方法，"1d" 或 "2d"
#' @param calib 校准方法，"LR", "rstar", "midp", 或 "Exact"
#' @return 包含上下界和元数据的列表
#' @examples
#' x <- rnorm(100)
#' ci <- snqesa_ci(x, tau = 0.9, alpha = 0.05)
snqesa_ci <- function(x, tau = 0.9, alpha = 0.05, ridge = NULL, tol = 1e-8,
                      engine = c("1d","2d"),
                      calib  = c("LR","rstar","midp","Exact"),
                      tail_split = c("equal","minlength"),
                      cc = c("jeffreys","half","anscombe","none")) {

  engine     <- match.arg(engine)
  calib      <- match.arg(calib)
  tail_split <- match.arg(tail_split)
  cc         <- match.arg(cc)
  x <- as.numeric(x); n <- length(x)
  if (is.null(ridge)) ridge <- n^(-3/4)
  qhat <- as.numeric(quantile(x, probs = tau, type = 8))
  sx   <- stats::sd(x); rng <- range(x)
  xs   <- sort(x)
  
  # ---------- 基础：由 t 得到观测 (S,Q) 与 |xobs| 以及 ESA 映射出的 u_x ----------
  # .u_from_t <- function(t){
  #   k <- sum(x <= t)
  #   S <- n*(tau - k/n)
  #   Q <- k*(1 - tau)^2 + (n - k)*tau^2
  #   xobs <- abs(S)/sqrt(Q + ridge)
  #   u <- snqesa_hinv(xobs, tau, n)  # 你的闭式+兜底
  #   list(u = as.numeric(u), k = as.integer(k))
  # }
  # 

  .u_from_t <- function(t, engine){
    k <- sum(x <= t)
    S <- n*(tau - k/n)
    Q <- k*(1 - tau)^2 + (n - k)*tau^2
    xobs <- S / sqrt(Q + ridge)              # 带符号的 xobs
    x_abs <- abs(xobs)
    
    # 关键：根据 engine 选择 1D/2D 映射得到 u0（位于 [0,1]）
    if (engine == "1d") {
      u0 <- snqesa_hinv(x_abs, tau, n)       # 原 1D ESA 反演
    } else { # engine == "2d"
      o2d <- tryCatch(snqesa_tail_2d(x_abs, tau, n), error = function(e) NULL)
      u0  <- if (!is.null(o2d) && is.finite(o2d$u)) o2d$u else snqesa_hinv(x_abs, tau, n)
    }
    
    # 朝向化：xobs≥0 表示“在 τ 左边”（左尾）；xobs<0 表示“在 τ 右边”（右尾）
    list(
      k = as.integer(k),
      xobs = xobs,
      uL = if (xobs >= 0) u0 else (1 - u0),   # 左尾用
      uR = if (xobs >= 0) (1 - u0) else u0    # 右尾用
    )
  }
  
  # ---------- 一侧鞍点尾：右尾的 LR / r*，带 0.5 连续性更正 ----------
  # 用 k 做连续性更正；符号仍由 “原始 u 与 tau 的大小” 决定
  .right_tail_saddle <- function(u, tau, n, k = NULL,
                                 cc = c("half","jeffreys","anscombe","none"),
                                 use_rstar = FALSE){
    cc <- match.arg(cc)
    u <- .clip01(u, 1e-12)
    
    # 仅用于幅值计算的连续性更正 u_eff；符号用原始 u 与 tau
    if (cc == "jeffreys" && !is.null(k)) {
      u_eff <- (k + 0.5) / (n + 1)
    } else if (cc == "anscombe" && !is.null(k)) {
      u_eff <- (k + 1/3) / (n + 2/3)
    } else if (cc == "half") {
      u_eff <- (n*u - 0.5)/n
    } else {
      u_eff <- u
    }
    u_eff <- .clip01(u_eff, 1e-12)
    
    # 符号由原始 u vs tau 决定（避免更正导致的误判）
    r_sign <- ifelse(u >= tau, 1, -1)
    
    KL  <- u_eff*log(u_eff/tau) + (1-u_eff)*log((1-u_eff)/(1-tau))
    r   <- r_sign * sqrt(2 * n * KL)
    qLR <- abs(log(u_eff/(1-u_eff)) - log(tau/(1-tau))) * sqrt(n * u_eff * (1 - u_eff))
    
    if (!is.finite(r) || !is.finite(qLR) || qLR <= 0) {
      z <- (u_eff - tau) / sqrt(tau * (1 - tau) / n)
      return(pnorm(-abs(z), lower.tail = TRUE))
    }
    
    if (use_rstar) {
      ratio <- r / (ifelse(r >= 0, qLR, -qLR))
      if (is.finite(ratio) && ratio > 0 && abs(r) > 1e-8) {
        rstar <- r + (1/r) * log(ratio)
        return(pnorm(-abs(rstar)))
      }
    }
    
    tail <- pnorm(-abs(r)) + dnorm(abs(r)) * (1/max(abs(r),1e-12) - 1/qLR)
    max(0, min(1, tail))
  }
  
  
  
  # ---------- 侧别一致的一侧 p 值 ----------
  # p_one <- function(t, side = c("L","R")){
  #   side <- match.arg(side)
  #   if (calib %in% c("midp","Exact")) {
  #     if (calib == "midp") {
  #       cm <- .cont_midp_sides(t, x, tau); return(if (side=="L") cm["pL"] else cm["pR"])
  #     } else {
  #       k <- sum(x <= t); k <- max(0L, min(n, k))
  #       return(if (side=="L") pbinom(k, n, tau) else 1 - pbinom(k - 1L, n, tau))
  #     }
  #   } else {
  #     tmp <- .u_from_t(t)
  #     ut  <- tmp$u
  #     kk  <- tmp$k  # 这里拿到 k
  #     if (side == "R") {
  #       return(.right_tail_saddle(ut, tau, n, k = kk, cc = cc, use_rstar = (calib=="rstar")))
  #     } else {
  #       # 左尾按互补（传入 k' = n - k，因为对应 1-ut, 1-tau 时“成功”计数为右侧）
  #       return(.right_tail_saddle(1 - ut, 1 - tau, n, k = n - kk, cc = cc, use_rstar = (calib=="rstar")))
  #     }
  #   }
  # }
  
  
  
  p_one <- function(t, side = c("L","R")){
    side <- match.arg(side)
    if (calib %in% c("midp","Exact")) {
      if (calib == "midp") {
        cm <- .cont_midp_sides(t, x, tau); return(if (side=="L") cm["pL"] else cm["pR"])
      } else {
        k <- sum(x <= t); k <- max(0L, min(n, k))
        return(if (side=="L") pbinom(k, n, tau) else 1 - pbinom(k - 1L, n, tau))
      }
    } else {
      tmp <- .u_from_t(t, engine)  # <<<<<<<<<<<<< 这里把 engine 传进来
      
      if (side == "R"){
        u_eff   <- tmp$uR
        tau_eff <- if (tmp$xobs >= 0) 1 - tau else tau
        k_eff   <- if (tmp$xobs >= 0) n - tmp$k else tmp$k
      } else { # "L"
        u_eff   <- tmp$uL
        tau_eff <- if (tmp$xobs >= 0) tau else 1 - tau
        k_eff   <- if (tmp$xobs >= 0) tmp$k else n - tmp$k
      }
      
      return(.right_tail_saddle(u_eff, tau_eff, n,
                                k = k_eff, cc = cc,
                                use_rstar = (calib=="rstar")))
    }
  }
  
  
  # ---------- 解单侧：p_one(t, side) = α/2 ----------
  solve_side <- function(side = c("L","R"), a = alpha){
    side <- match.arg(side)
    target <- a/2
    f <- function(t) p_one(t, side) - target
    
    mid <- qhat
    sx_ <- if (is.finite(sx) && sx > 0) sx else max(1e-6, (rng[2]-rng[1])/10)
    fmid <- f(mid); if (!is.finite(fmid)) return(NA_real_)
    
    step <- sx_
    # 左侧 p_L(t) 单调递增；右侧 p_R(t) 单调递减
    if (side == "L"){
      if (fmid >= 0){               # 解在左边
        x <- mid
        for (i in 1:80){
          x <- x - step
          fx <- f(x)
          if (is.finite(fx) && fx < 0) return(uniroot(f, c(x, mid), tol = tol)$root)
          step <- step * 1.6
        }
      } else {                       # 解在右边
        x <- mid
        for (i in 1:80){
          x <- x + step
          fx <- f(x)
          if (is.finite(fx) && fx > 0) return(uniroot(f, c(mid, x), tol = tol)$root)
          step <- step * 1.6
        }
      }
    } else { # side == "R"（单调递减）
      if (fmid <= 0){               # 解在左边
        x <- mid
        for (i in 1:80){
          x <- x - step
          fx <- f(x)
          if (is.finite(fx) && fx > 0) return(uniroot(f, c(x, mid), tol = tol)$root)
          step <- step * 1.6
        }
      } else {                      # 解在右边
        x <- mid
        for (i in 1:80){
          x <- x + step
          fx <- f(x)
          if (is.finite(fx) && fx < 0) return(uniroot(f, c(mid, x), tol = tol)$root)
          step <- step * 1.6
        }
      }
    }
    NA_real_
  }
  
  
  # ---------- 连续反演：等尾或不等尾最短 ----------
  if (tail_split == "equal") {
    tL <- suppressWarnings(solve_side("L", a = alpha))
    tU <- suppressWarnings(solve_side("R", a = alpha))
  } else {
    grid <- seq(0, alpha, length.out = 31)
    best <- c(Inf, NA, NA)
    for (aL in grid) {
      aR <- alpha - aL
      Lc <- suppressWarnings(solve_side("L", a = aL))
      Uc <- suppressWarnings(solve_side("R", a = aR))
      if (is.finite(Lc) && is.finite(Uc)) {
        len <- Uc - Lc
        if (is.finite(len) && len < best[1]) best <- c(len, Lc, Uc)
      }
    }
    tL <- best[2]; tU <- best[3]
  }
  
  list(lower = tL, upper = tU, qhat = qhat, tau = tau, alpha = alpha,
       type = "continuous", calib = calib, engine = engine, tail_split = tail_split)
}

## ---------- 8. 现有基线 ----------

###############################################################################
# 方法：Wald + KDE（核密度估计）分位数置信区间
# 名称对照：WaldKDE / KDE-Wald / “密度-德尔塔法”
#
# 思路：
#   令 q̂_τ 为样本分位数（type=8）。Wald 近似给出
#     q̂_τ ~ N(q_τ, Var) ，其中 Var ≈ τ(1-τ) / (n f(q_τ)^2)
#   用核密度估计 f̂(q̂_τ) 代替 f(q_τ)，据此构造对称区间 q̂_τ ± z * se。
#
# 优点：
#   - 简单、闭式、速度快（O(n log n) ~ O(n)）
#   - 直观可控：带宽 bw 可手调
#
# 局限：
#   - 高 τ/低 τ 处对 f(q_τ) 估计不稳 → 区间易偏
#   - 分位点处密度不连续/尖点（如 ALD、split uniform）会失真
#   - 需要“足够平滑”的分布假设；极重尾/多峰会受影响
#
# 参数：
#   - bw: 传给 density() 的带宽；若为 NULL，使用默认规则
#   - floor_f: 给 f̂ 的安全下界，避免 se 爆炸（这里用 IQR 量级）
###############################################################################
wald_kde_ci <- function(x, tau = 0.9, alpha = 0.05, bw = NULL){
  n <- length(x)
  qhat <- as.numeric(quantile(x, probs = tau, type = 8))
  dens <- if (is.null(bw)) density(x) else density(x, bw = bw)
  fhat <- approx(dens$x, dens$y, xout = qhat, rule = 2)$y
  # 安全下界：IQR-based or Silverman 的数量级
  sx <- sd(x); iqr <- IQR(x)
  floor_f <- 1 / (sqrt(n) * (iqr + 1e-8))  # 也可用 1/(sx*sqrt(n))
  fhat <- max(fhat, floor_f)
  z <- qnorm(1 - alpha/2)
  se <- sqrt(tau * (1 - tau) / (n * fhat^2))
  c(lower = qhat - z * se, upper = qhat + z * se, qhat = qhat)
}

###############################################################################
# 方法：Exact（二项反演的“严格序统计”区间）
# 名称对照：Order-stat exact / Clopper–Pearson for quantiles（反演法）
#
# 思路：
#   令 K = #{X_i ≤ t}。在 H0: F(t) = τ 下，K ~ Bin(n, τ)。
#   通过二项分布的等尾或不等尾反演，找到秩区间 [Lk, Uk)，
#   再映射为样本序统计的区间 [x_(Lk), x_(Uk)]（Uk = n+1 表示 +∞）。
#
# 优点：
#   - 真·有限样本控制，覆盖 ≥ 1 - α（保守）
#   - 对分布的平滑性几乎没有要求（只需 i.i.d.）
#
# 局限：
#   - 区间偏保守（尤其在 n 较小/极端 τ）
#   - 端点依赖样本秩 → 可能较宽
###############################################################################
exact_order_ci <- function(x, tau = 0.9, alpha = 0.05){
  x <- sort(as.numeric(x)); n <- length(x)
  alphaL <- alpha / 2; alphaR <- alpha - alphaL
  Lk <- qbinom(alphaL, n, tau)
  Uk <- qbinom(1 - alphaR, n, tau) + 1L
  
  Lk <- max(1L, min(n,    as.integer(Lk)))
  Uk <- max(1L, min(n+1L, as.integer(Uk)))
  
  lower_val <- x[Lk]
  upper_val <- if (Uk == n+1L) Inf else x[Uk]
  c(lower = lower_val, upper = upper_val, L = Lk, U = Uk)
}

###############################################################################
# 方法：Bootstrap 百分位 & BCa（Bias-Corrected and Accelerated）
# 名称对照：Percentile / BCa
#
# 思路：
#   - Percentile：用 bootstrap 分布的 (α/2, 1-α/2) 分位数作区间端点
#   - BCa：对 Percentile 做两项修正
#       * 偏倚修正 z0：基于 B 次 bootstrap 与 θ0 的相对位置
#       * 加速修正 a：基于 jackknife 近似估计
#
# 优点：
#   - 对统计量分布形状更鲁棒（非对称、偏态）
#   - BCa 通常比 Percentile 覆盖更准
#
# 局限：
#   - 成本高（B 次重采样 + jackknife）
#   - 在极端 τ/小样本时仍可能失稳（需要较大 B）
###############################################################################
bca_ci <- function(x, tau = 0.9, alpha = 0.05, B = 2000, seed = NULL, type = c("percentile", "bca")){
  type <- match.arg(type)
  n <- length(x); if (!is.null(seed)) set.seed(seed)
  stat <- function(z) as.numeric(quantile(z, tau, type = 8))
  theta0 <- stat(x)
  idx <- matrix(sample.int(n, size = n*B, replace = TRUE), nrow = n, ncol = B)
  thetab <- apply(idx, 2, function(j) stat(x[j]))
  if (type == "percentile"){
    qs <- quantile(thetab, probs = c(alpha/2, 1 - alpha/2), names = FALSE)
    return(c(lower = qs[1], upper = qs[2], qhat = theta0))
  }
  z0 <- qnorm(mean(thetab < theta0))
  thetaj <- numeric(n)
  for (i in 1:n) thetaj[i] <- stat(x[-i])
  u <- mean(thetaj) - thetaj
  a_acc <- sum(u^3) / (6 * (sum(u^2))^(3/2) + 1e-15)
  zal <- qnorm(alpha/2); zau <- qnorm(1 - alpha/2)
  adj <- function(z){ pnorm(z0 + (z0 + z) / (1 - a_acc * (z0 + z))) }
  a1 <- adj(zal); a2 <- adj(zau)
  qs <- quantile(thetab, probs = c(a1, a2), names = FALSE)
  c(lower = qs[1], upper = qs[2], qhat = theta0)
}

## ---------- 9. 新增基线 ----------


###############################################################################
# 方法：Harrell–Davis 分位数 + 其族方法
# 名称对照：HD estimator（加权秩）；其上派生的 MJ、HD-Boot 等
#
# 共同背景：
#   Harrell–Davis（HD）把分位数估计写成对所有秩的 Beta 权重加权平均，
#   q̂_HD = Σ w_k x_(k)，其中 w_k 来源于 Beta(a,b) 在 k/n 的 CDF 差分。
#
# 子方法：
#   1) ci_hd_boot      —— 对 q̂_HD 做非参数 bootstrap 的百分位区间
#   2) ci_mj (Maritz–Jarrett) —— 用 HD 权重计算“经验方差”，再配合正态近似
#
# 适用与局限：
#   - HD 在连续分布、适中样本量时对分位点更平滑；但在极端 τ 或离散/尖点密度处，
#     可能比“秩区间”类方法更乐观或偏。
###############################################################################
.hd_weights <- function(n, tau){
  a <- tau * (n + 1); b <- (1 - tau) * (n + 1)
  pb <- pbeta((0:n)/n, a, b)
  w  <- pb[-1] - pb[-(n+1)]
  w  <- pmax(w, 0); s <- sum(w); if (s <= 0) w[] <- 1/n else w <- w/s
  w
}

.hd_qhat <- function(x, tau){
  xx <- sort(as.numeric(x)); w <- .hd_weights(length(xx), tau); sum(w * xx)
}
###############################################################################
# 方法：HD + Bootstrap 百分位（非参重采样）
# 名称对照：HD_Boot / HD-percentile
#
# 思路：
#   以 HD 分位数 q̂_HD 为统计量；重复 B 次非参重采样，取 bootstrap
#   分布的 (α/2, 1-α/2) 分位数作为区间。
###############################################################################
ci_hd_boot <- function(x, tau, alpha = 0.05, B = 1000, seed = NULL){
  if (!is.null(seed)) set.seed(seed)
  n <- length(x); qboot <- numeric(B)
  for (b in seq_len(B)) qboot[b] <- .hd_qhat(sample(x, n, TRUE), tau)
  c(lower = quantile(qboot, alpha/2, names = FALSE),
    upper = quantile(qboot, 1 - alpha/2, names = FALSE))
}

###############################################################################
# 方法：Maritz–Jarrett（HD-variance 正态近似）
# 名称对照：MaritzJar / MJ
#
# 思路：
#   用 HD 权重 w_k 在样本上计算“二阶矩”来近似 Var(q̂_HD)：
#     v2 = Σ w_k (x_(k) - q̂_HD)^2
#   然后 Wald 形式构造 q̂_HD ± z * sqrt(v2)。
#
# 优点：
#   - 闭式、成本低
# 局限：
#   - 仍是正态近似；极端分位/尖点分布下可能偏
###############################################################################
ci_mj <- function(x, tau, alpha = 0.05){
  xx <- sort(as.numeric(x)); n <- length(xx)
  w  <- .hd_weights(n, tau)
  qh <- sum(w * xx)
  v2 <- sum(w * (xx - qh)^2)
  se <- sqrt(max(v2, 0))
  z  <- qnorm(1 - alpha/2)
  c(lower = qh - z*se, upper = qh + z*se)
}

###############################################################################
# 方法：平滑 Bootstrap（Smoothed Bootstrap）
# 名称对照：SmBoot
#
# 思路：
#   在普通 bootstrap 的基础上，对样本再加一个 N(0, bw^2) 的平滑噪声，
#   缓解分位函数的离散/台阶效应；对极端 τ 有时更稳。
#
# 备注：
#   - bw 未给定时用 bw.nrd0 或 n^{-1/5} 量级
###############################################################################
ci_smoothed_boot <- function(x, tau, alpha = 0.05, B = 1000, bw = NULL, seed = NULL){
  if (!is.null(seed)) set.seed(seed)
  n <- length(x)
  if (is.null(bw)) bw <- stats::bw.nrd0(x)
  if (!is.finite(bw) || bw <= 0) bw <- sd(x) * n^(-1/5)
  qboot <- numeric(B)
  for (b in seq_len(B)) {
    xb <- sample(x, n, replace = TRUE) + rnorm(n, sd = bw)
    qboot[b] <- as.numeric(quantile(xb, tau, type = 8))
  }
  c(lower = quantile(qboot, alpha/2, names = FALSE),
    upper = quantile(qboot, 1 - alpha/2, names = FALSE))
}


###############################################################################
# 方法：m-out-of-n Bootstrap（子样本重采样）
# 名称对照：mOutOfn
#
# 思路：
#   每次从样本中有放回抽取 m（m<n）个点做分位数估计；重复 B 次，
#   用百分位法取端点。
#
# 作用：
#   - 在重尾/极端分位时，比 n-out-of-n 更稳（减小单次样本波动）
#
# 选 m：
#   - 常见经验：m ≈ n^0.7；过小会增大方差，过大接近普通 bootstrap
###############################################################################
ci_m_out_n <- function(x, tau, alpha = 0.05, B = 1000, m = NULL, seed = NULL){
  if (!is.null(seed)) set.seed(seed)
  n <- length(x); if (is.null(m)) m <- max(5L, floor(n^0.7))
  qboot <- numeric(B)
  for (b in seq_len(B)) qboot[b] <- as.numeric(quantile(sample(x, m, TRUE), tau, type = 8))
  c(lower = quantile(qboot, alpha/2, names = FALSE),
    upper = quantile(qboot, 1 - alpha/2, names = FALSE))
}

###############################################################################
# 方法：Subsampling（无放回子样本 + 透镜枢轴）
# 名称对照：Subsample / Politis–Romano–Wolf 风格
#
# 思路：
#   抽取无放回子样本大小 b（b<n），计算其分位数 q̂_b；
#   利用“枢轴” 2 q̂_n - q̂_b 的经验分布构造区间。
#
# 作用：
#   - 更适合相关性/依赖结构或重尾（这里仍假设 i.i.d.，但思路一致）
###############################################################################
ci_subsampling <- function(x, tau, alpha = 0.05, B = 1000, b = NULL, seed = NULL){
  if (!is.null(seed)) set.seed(seed)
  n <- length(x); if (is.null(b)) b <- max(5L, min(n - 1L, floor(n^0.7)))
  qhat <- as.numeric(quantile(x, tau, type = 8))
  qsub <- numeric(B)
  for (i in seq_len(B)) qsub[i] <- as.numeric(quantile(sample(x, b, FALSE), tau, type = 8))
  piv <- 2*qhat - qsub
  c(lower = quantile(piv, alpha/2, names = FALSE),
    upper = quantile(piv, 1 - alpha/2, names = FALSE))
}


###############################################################################
# 方法：Hall–Sheather / Nyblom（quantileCI 包封装）
# 名称对照：HS_Nyblom（中位数用 HS；一般分位用 Nyblom）
#
# 思路：
#   - 中位数：Hall–Sheather 在密度估计与带宽选择上给出一类更稳健的近似
#   - 一般 τ：Nyblom 的分位数区间近似
#
# 优点：
#   - 有较好的理论误差界；实现成熟（quantileCI）
# 局限：
#   - 依赖外部包；在极端 τ 或很小样本上也可能保守/偏
###############################################################################
ci_hsnyblom <- function(x, tau, alpha = 0.05){
  if (!requireNamespace("quantileCI", quietly = TRUE)) return(c(lower = NA_real_, upper = NA_real_))
  if (abs(tau - 0.5) < .Machine$double.eps){
    lu <- tryCatch(quantileCI::median_confint_hs(x, conf.level = 1 - alpha),
                   error = function(e) c(NA_real_, NA_real_))
    lu <- unlist(lu); names(lu) <- c("lower","upper"); return(lu)
  } else {
    lu <- tryCatch(quantileCI::quantile_confint_nyblom(x = x, p = tau, conf.level = 1 - alpha),
                   error = function(e) c(NA_real_, NA_real_))
    lu <- unlist(lu); names(lu) <- c("lower","upper"); return(lu)
  }
}


bca_and_percentile_ci <- function(x, tau = 0.9, alpha = 0.05, B = 2000, seed = NULL){
  n <- length(x); if (!is.null(seed)) set.seed(seed)
  stat <- function(z) as.numeric(quantile(z, tau, type = 8))
  theta0 <- stat(x)
  idx <- matrix(sample.int(n, size = n*B, replace = TRUE), nrow = n, ncol = B)
  thetab <- apply(idx, 2, function(j) stat(x[j]))
  
  # Percentile
  pct <- quantile(thetab, probs = c(alpha/2, 1 - alpha/2), names = FALSE)
  
  # BCa
  z0 <- qnorm(mean(thetab < theta0))
  thetaj <- vapply(seq_len(n), function(i) stat(x[-i]), numeric(1))
  u <- mean(thetaj) - thetaj
  a_acc <- sum(u^3) / (6 * (sum(u^2))^(3/2) + 1e-15)
  zal <- qnorm(alpha/2); zau <- qnorm(1 - alpha/2)
  adj <- function(z){ pnorm(z0 + (z0 + z) / (1 - a_acc * (z0 + z))) }
  a1 <- adj(zal); a2 <- adj(zau)
  bca <- quantile(thetab, probs = c(a1, a2), names = FALSE)
  
  list(
    Pct = c(lower = pct[1], upper = pct[2], qhat = theta0),
    BCa = c(lower = bca[1], upper = bca[2], qhat = theta0)
  )
}

## ---------- 10. 单数据集比较 ----------
compare_ci_methods <- function(x, tau = 0.9, alpha = 0.05, B = 2000, seed = 123,
                               t_true = NULL){
  x <- as.numeric(x)
  n <- length(x)
  if (!is.null(seed)) set.seed(seed)
  
  # --- 统一的点估计（所有方法共享便于对照） ---
  qhat_t8 <- as.numeric(quantile(x, probs = tau, type = 8))
  qhat_HD <- tryCatch(.hd_qhat(x, tau), error = function(e) NA_real_)
  
  # 打包一行：从 (lower, upper[, qhat]) 生成丰富指标
  pack_row <- function(ci, method){
    L <- as.numeric(ci["lower"]); U <- as.numeric(ci["upper"])
    len    <- if (is.finite(L) && is.finite(U)) U - L else NA_real_
    center <- if (is.finite(L) && is.finite(U)) (L + U)/2 else NA_real_
    qhat_in <- if (is.finite(L)) as.numeric(qhat_t8 >= L & (if(is.finite(U)) qhat_t8 <= U else TRUE)) else NA_real_
    hit_true <- if (!is.null(t_true) && is.finite(t_true) && is.finite(L))
      as.numeric(t_true >= L & (if(is.finite(U)) t_true <= U else TRUE))
    else NA_real_
    dist_qhat_true <- if (!is.null(t_true) && is.finite(t_true)) qhat_t8 - t_true else NA_real_
    true_rel_pos <- if (!is.null(t_true) && is.finite(t_true) && is.finite(L) && is.finite(U) && (U>L))
      (t_true - L)/(U - L) else NA_real_
    
    c(method=method, lower=L, upper=U, length=len, center=center,
      qhat_t8=qhat_t8, qhat_HD=qhat_HD, qhat_in_CI=qhat_in, hit_true=hit_true,
      dist_qhat_true=dist_qhat_true, true_rel_pos=true_rel_pos)
  }
  
  res <- list()
  
  # --- 你的各方法：尽量沿用原实现，必要时从返回对象中取出 lower/upper ---
  ci_snd <- snqesa_ci_discrete(x, tau = tau, alpha = alpha, tail_split = "equal")
  res$SNQESA_disc <- pack_row(c(lower = ci_snd$lower, upper = ci_snd$upper), "SNQESA_disc")
  
  ci_snm <- snqesa_ci_discrete(x, tau = tau, alpha = alpha, tail_split = "minlength")
  res$SNQESA_minL <- pack_row(c(lower = ci_snm$lower, upper = ci_snm$upper), "SNQESA_minL")
  
  ci_sn  <- snqesa_ci(x, tau, alpha, calib = "rstar", engine = "1d", tail_split = "minlength")
  res$SNQESA_cont <- pack_row(c(lower = ci_sn$lower, upper = ci_sn$upper), "SNQESA_cont")
  
  ci_wald <- wald_kde_ci(x, tau = tau, alpha = alpha)
  res$WaldKDE <- pack_row(ci_wald[c("lower","upper")], "WaldKDE")
  
  ci_ex <- exact_order_ci(x, tau = tau, alpha = alpha)
  res$Exact <- pack_row(ci_ex[c("lower","upper")], "Exact")
  
  ci_pct <- bca_ci(x, tau = tau, alpha = alpha, B = B, seed = seed, type = "percentile")
  res$PctBoot <- pack_row(ci_pct[c("lower","upper")], "PctBoot")
  
  ci_bca <- bca_ci(x, tau = tau, alpha = alpha, B = B, seed = seed, type = "bca")
  res$BCa <- pack_row(ci_bca[c("lower","upper")], "BCa")
  
  ci_hd  <- ci_hd_boot(x, tau, alpha, B = B, seed = seed)
  res$HD_Boot <- pack_row(ci_hd[c("lower","upper")], "HD_Boot")
  
  ci_mj_ <- ci_mj(x, tau, alpha)
  res$MaritzJar <- pack_row(ci_mj_[c("lower","upper")], "MaritzJar")
  
  ci_smb <- ci_smoothed_boot(x, tau, alpha, B = B, seed = seed)
  res$SmBoot <- pack_row(ci_smb[c("lower","upper")], "SmBoot")
  
  ci_sub <- ci_subsampling(x, tau, alpha, B = B, seed = seed)
  res$Subsample <- pack_row(ci_sub[c("lower","upper")], "Subsample")
  
  ci_mn  <- ci_m_out_n(x, tau, alpha, B = B, seed = seed)
  res$mOutOfn <- pack_row(ci_mn[c("lower","upper")], "mOutOfn")
  
  ci_hsn <- ci_hsnyblom(x, tau, alpha)
  res$HS_Nyblom <- pack_row(ci_hsn[c("lower","upper")], "HS_Nyblom")
  
  # 组装与排序
  out <- do.call(rbind, res)
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  # 把数值列转为 numeric（method 列保留字符）
  num_cols <- setdiff(names(out), "method")
  out[num_cols] <- lapply(out[num_cols], function(v) suppressWarnings(as.numeric(v)))
  
  out <- out[order(out$length), , drop = FALSE]
  rownames(out) <- out$method
  
  # 附加样本信息（便于在控制台或 knit 时查看）
  sample_info <- list(
    n = n, tau = tau, alpha = alpha,
    mean = mean(x), sd = sd(x), min = min(x), max = max(x),
    qhat_t8 = qhat_t8, qhat_HD = qhat_HD,
    t_true = if (!is.null(t_true)) t_true else NA_real_
  )
  attr(out, "sample_info") <- sample_info
  
  out
}

## ---------- 11. Monte Carlo 汇总 ----------
summarise_matrix <- function(M){
  stopifnot(all(c("cover","length","time") %in% colnames(M)))
  has <- function(name) name %in% colnames(M)
  
  ok_cov <- is.finite(M[,"cover"]); n_eff <- sum(ok_cov)
  covg   <- if (n_eff == 0) NA_real_ else mean(M[ok_cov,"cover"], na.rm = TRUE)
  se     <- if (n_eff == 0) NA_real_ else sqrt(covg * (1 - covg) / n_eff)
  
  ok_len <- is.finite(M[,"length"]) & (M[,"length"] >= 0)
  ok_tm  <- is.finite(M[,"time"])
  
  mean_len   <- if (!any(ok_len))  NA_real_ else mean(M[ok_len,"length"], na.rm = TRUE)
  median_len <- if (!any(ok_len))  NA_real_ else median(M[ok_len,"length"], na.rm = TRUE)
  mean_time  <- if (!any(ok_tm))   NA_real_ else mean(M[ok_tm,"time"],  na.rm = TRUE)
  
  center_bias_mean <- if (has("center_bias")) mean(M[,"center_bias"], na.rm=TRUE) else NA_real_
  center_bias_med  <- if (has("center_bias")) median(M[,"center_bias"], na.rm=TRUE) else NA_real_
  center_bias_rmse <- if (has("center_bias")) { v <- M[,"center_bias"]; v <- v[is.finite(v)]; if(length(v)) sqrt(mean(v^2)) else NA_real_ } else NA_real_
  
  qhat_in_rate     <- if (has("qhat_in")) mean(M[,"qhat_in"], na.rm=TRUE) else NA_real_
  relpos_mean      <- if (has("relpos"))  mean(M[,"relpos"],  na.rm=TRUE) else NA_real_
  lower_miss_rate  <- if (has("lower_miss")) mean(M[,"lower_miss"], na.rm=TRUE) else NA_real_
  upper_miss_rate  <- if (has("upper_miss")) mean(M[,"upper_miss"], na.rm=TRUE) else NA_real_
  
  ## 新增：Interval Score / Pinball 的均值（越小越好）
  mean_is   <- if (has("is"))      mean(M[,"is"],      na.rm=TRUE) else NA_real_
  median_is <- if (has("is"))      median(M[,"is"],    na.rm=TRUE) else NA_real_
  mean_pinb <- if (has("pinball")) mean(M[,"pinball"], na.rm=TRUE) else NA_real_
  
  fail_rate <- mean(!(ok_cov | ok_len))
  
  data.frame(
    coverage   = covg, cover_se = se,
    mean_len   = mean_len, median_len = median_len,
    fail_rate  = fail_rate, mean_time = mean_time,
    mean_center_bias  = center_bias_mean,
    median_center_bias= center_bias_med,
    rmse_center_bias  = center_bias_rmse,
    qhat_in_rate      = qhat_in_rate,
    mean_true_relpos  = relpos_mean,
    lower_miss_rate   = lower_miss_rate,
    upper_miss_rate   = upper_miss_rate,
    mean_IS           = mean_is,        # <--- 新增
    median_IS         = median_is,      # <--- 新增
    mean_pinball      = mean_pinb,      # <--- 新增
    row.names = NULL
  )
}


## ------ Proper scoring for intervals ------
interval_score <- function(L, U, theta, alpha){
  if (!is.finite(L) || !is.finite(U)) return(NA_real_)  # 无穷端点不给分
  w <- U - L
  pen <- if (theta < L) (2/alpha) * (L - theta) else if (theta > U) (2/alpha) * (theta - U) else 0
  w + pen
}

## 单个分位点的 pinball loss（可选：看点估计的精度）
pinball_loss <- function(qhat, theta, tau){
  e <- theta - qhat
  ((tau - as.numeric(e < 0)) * e)
}

wis_from_ci_list <- function(ci_list, theta, alphas){
  stopifnot(length(ci_list) == length(alphas))
  K <- length(alphas)
  terms <- numeric(K)
  for (k in seq_len(K)){
    Lk <- ci_list[[k]]["lower"]; Uk <- ci_list[[k]]["upper"]
    terms[k] <- (alphas[k]/2) * interval_score(Lk, Uk, theta, alphas[k])
  }
  ## 这里没有“中位数误差项”，因为我们比较的是“参数区间”而非“预测区间”
  mean(terms)  # 越小越好
}

dm_test_is <- function(res_object, method_A, method_B){
  A <- res_object$raw[[method_A]][,"is"]
  B <- res_object$raw[[method_B]][,"is"]
  d <- A - B
  d <- d[is.finite(d)]
  n <- length(d)
  dm_stat <- mean(d) / sqrt(var(d)/n)
  pval <- 2 * (1 - pnorm(abs(dm_stat)))
  c(stat = dm_stat, p.value = pval, n = n,
    mean_IS_A = mean(A, na.rm=TRUE), mean_IS_B = mean(B, na.rm=TRUE))
}


# =============================================================================
# simulate_compare_methods() 2.0
# - 直接透传/控制 snqesa_ci 的求解相关参数（engine / calib / tail_split / ridge / tol）
# - 支持一次性跑多种 SNQESA 连续反演“变体”（snqesa_variants），自动命名收集
# - 仍保留离散 SNQESA 与各基线方法，汇总接口不变
# 依赖：你现有文件中已定义的：
#   rsamp(), true_q(), snqesa_ci(), snqesa_ci_discrete(), wald_kde_ci(),
#   exact_order_ci(), bca_ci(), ci_hd_boot(), ci_mj(), ci_smoothed_boot(),
#   ci_subsampling(), ci_m_out_n(), ci_hsnyblom(), summarise_matrix()
# =============================================================================
simulate_compare_methods <- function(
    R = 500, n = 50, tau = 0.9, alpha = 0.05,
    dist = c("normal","lognormal","t2","t3","cauchy","mixnorm", 'logistic', 'mixnorm_soft', 'beta', 'exp'),
    B = 1000, seed = 123, parallel = FALSE,
    
    snqesa_engine     = c("1d","2d"),     # ESA 1d 映射 or 2d 受约束
    snqesa_calib      = c("LR","rstar","midp","Exact"),
    snqesa_tail_split = c("equal","minlength"),
    snqesa_cc         = c("jeffreys","half","anscombe","none"),  # 连续性更正
    snqesa_ridge      = NULL,             # 可为 NULL / 数值 / 函数(n) -> 数值
    snqesa_enable     = TRUE              # 是否计算 SNQESA_cont（调试用）
){
  # 兼容 dist 的别名
  dist <- .normalize_dist(match.arg(.normalize_dist(dist),
                                    c("normal","lognormal","t2","t3","cauchy","mixnorm", 'logistic', "mixnorm_soft", 'beta', 'exp')))
  
  if (!is.null(seed)) set.seed(seed)
  
  snqesa_engine     <- match.arg(snqesa_engine)
  snqesa_calib      <- match.arg(snqesa_calib)
  snqesa_tail_split <- match.arg(snqesa_tail_split)
  snqesa_cc         <- match.arg(snqesa_cc)
  
  # 统一方法集
  METHODS <- c("SNQESA_disc","SNQESA_minL","SNQESA_cont", "SNQESA_rand","SNQESA_auto",
               "WaldKDE","Exact","PctBoot","BCa",
               "HD_Boot","MaritzJar","SmBoot","Subsample","mOutOfn","HS_Nyblom")
  
  # 结果矩阵骨架
  COLS <- c("cover","length","time","lower","upper",
            "center_bias","qhat_in","relpos","lower_miss","upper_miss",
            "is","pinball")
  acc <- lapply(METHODS, function(.) {
    matrix(NA_real_, nrow = R, ncol = length(COLS), dimnames = list(NULL, COLS))
  })
  names(acc) <- METHODS
  
  sample_cols <- c("qhat_t8","qhat_HD","mean","sd","min","max")
  samples <- matrix(NA_real_, nrow = R, ncol = length(sample_cols),
                    dimnames = list(NULL, sample_cols))
  
  make_na_vec <- function(){
    v <- rep(NA_real_, length(COLS)); names(v) <- COLS; v
  }
  
  tq <- true_q(dist, tau)  # 真分位（用于覆盖与偏差）
  
  pack <- function(ci, t0, t1, qhat_t8){
    getL <- function(ci) if (!is.null(names(ci)) && "lower" %in% names(ci)) as.numeric(ci["lower"]) else as.numeric(ci[1])
    getU <- function(ci) if (!is.null(names(ci)) && "upper" %in% names(ci)) as.numeric(ci["upper"]) else as.numeric(ci[2])
    L <- suppressWarnings(getL(ci)); U <- suppressWarnings(getU(ci))
    
    len <- if (is.finite(L) && is.finite(U)) U - L else NA_real_
    #cov <- if (is.finite(tq) && is.finite(L)) as.numeric(tq >= L & (if (is.finite(U)) tq <= U else TRUE)) else NA_real_
    cov <- if (is.finite(tq) && is.finite(L) && is.finite(U)) as.numeric(L <= tq && tq <= U) else 0
    ctr <- if (is.finite(L) && is.finite(U)) (L + U)/2 else NA_real_
    qin <- if (is.finite(L)) as.numeric(qhat_t8 >= L & (if (is.finite(U)) qhat_t8 <= U else TRUE)) else NA_real_
    relpos <- if (is.finite(L) && is.finite(U) && is.finite(len) && len > 0 && cov == 1) (tq - L) / len else NA_real_
    lower_miss <- if (is.finite(L)) as.numeric(tq < L) else NA_real_
    upper_miss <- if (is.finite(U)) as.numeric(tq > U) else NA_real_
    
    ## 新增：Interval Score & pinball（点估计用 type=8 的分位数）
    is_val <- interval_score(L, U, theta = tq, alpha = alpha)
    pinb   <- pinball_loss(qhat_t8, theta = tq, tau = tau)
    
    c(cover = cov, length = len, time = (t1 - t0),
      lower = L, upper = U,
      center_bias = if (is.finite(ctr) && is.finite(tq)) ctr - tq else NA_real_,
      qhat_in = qin, relpos = relpos,
      lower_miss = lower_miss, upper_miss = upper_miss,
      is = is_val, pinball = pinb)
  }
  
  # ridge 处理器：允许 NULL / 数值 / 函数
  ridgify <- function(n){
    if (is.null(snqesa_ridge)) {
      # 更稳健的默认值（比 n^(-3/4) 稍大一点，改善覆盖）
      return(max(n^(-2/3), 3/n))
    }
    if (is.function(snqesa_ridge)) return(as.numeric(snqesa_ridge(n)))
    as.numeric(snqesa_ridge)
  }
  
  runner <- function(r){
    s <- if (is.null(seed)) NULL else seed + r
    if (!is.null(s)) set.seed(s)
    x <- rsamp(n, dist)
    
    qhat_t8 <- as.numeric(quantile(x, probs = tau, type = 8))
    qhat_HD <- tryCatch(.hd_qhat(x, tau), error = function(e) NA_real_)
    smry <- c(qhat_t8 = qhat_t8, qhat_HD = qhat_HD,
              mean = mean(x), sd = sd(x), min = min(x), max = max(x))
    
    out <- setNames(vector("list", length(METHODS)), METHODS)
    
    t0 <- proc.time()[3]; ci <- tryCatch(snqesa_ci_discrete(x, tau=tau, alpha=alpha, tail_split="equal"),
                                         error=function(e) NULL); t1 <- proc.time()[3]
    out$SNQESA_disc <- if (is.null(ci)) make_na_vec() else pack(c(ci$lower,ci$upper), t0, t1, qhat_t8)
    
    t0 <- proc.time()[3]; ci <- tryCatch(snqesa_ci_discrete(x, tau=tau, alpha=alpha, tail_split="minlength"),
                                         error=function(e) NULL); t1 <- proc.time()[3]
    out$SNQESA_minL <- if (is.null(ci)) make_na_vec() else pack(c(ci$lower,ci$upper), t0, t1, qhat_t8)
    
    
    t0 <- proc.time()[3]; ci <- tryCatch(snqesa_ci_discrete_rand(x, tau=tau, alpha=alpha, tail_split="minlength"),
                                         error=function(e) NULL); t1 <- proc.time()[3]
    out$SNQESA_rand <- if (is.null(ci)) make_na_vec() else pack(c(ci$lower, ci$upper), t0, t1, qhat_t8)
    
    t0 <- proc.time()[3]
    ci <- tryCatch(snqesa_ci_auto(x, tau = tau, alpha = alpha), error = function(e) NULL)
    t1 <- proc.time()[3]
    out$SNQESA_auto <- if (is.null(ci)) make_na_vec() else pack(c(ci$lower, ci$upper), t0, t1, qhat_t8)
    if (snqesa_enable) {
      # ===== 关键：把求解设置传给 snqesa_ci =====
      t0 <- proc.time()[3]
      ci <- tryCatch(
        snqesa_ci(x,
                  tau = tau, alpha = alpha,
                  engine = snqesa_engine,
                  calib  = snqesa_calib,
                  tail_split = snqesa_tail_split,
                  cc = snqesa_cc,
                  ridge = ridgify(n)),
        error = function(e) NULL
      )
      t1 <- proc.time()[3]
      out$SNQESA_cont <- if (is.null(ci)) make_na_vec() else pack(c(ci$lower,ci$upper), t0, t1, qhat_t8)
    } else {
      out$SNQESA_cont <- make_na_vec()
    }
    
    t0 <- proc.time()[3]; ci <- tryCatch(wald_kde_ci(x, tau=tau, alpha=alpha)[c("lower","upper")],
                                         error=function(e) c(NA,NA)); t1 <- proc.time()[3]
    out$WaldKDE <- pack(ci, t0, t1, qhat_t8)
    
    t0 <- proc.time()[3]; ci <- tryCatch(exact_order_ci(x, tau=tau, alpha=alpha)[c("lower","upper")],
                                         error=function(e) c(NA,NA)); t1 <- proc.time()[3]
    out$Exact <- pack(ci, t0, t1, qhat_t8)
    
    t0 <- proc.time()[3]; ci <- tryCatch(bca_ci(x, tau=tau, alpha=alpha, B=B, type="percentile")[c("lower","upper")],
                                         error=function(e) c(NA,NA)); t1 <- proc.time()[3]
    out$PctBoot <- pack(ci, t0, t1, qhat_t8)
    
    t0 <- proc.time()[3]; ci <- tryCatch(bca_ci(x, tau=tau, alpha=alpha, B=B, type="bca")[c("lower","upper")],
                                         error=function(e) c(NA,NA)); t1 <- proc.time()[3]
    out$BCa <- pack(ci, t0, t1, qhat_t8)
    
    t0 <- proc.time()[3]; ci <- tryCatch(ci_hd_boot(x, tau, alpha, B=B, seed = s)[c("lower","upper")],
                                         error=function(e) c(NA,NA)); t1 <- proc.time()[3]
    out$HD_Boot <- pack(ci, t0, t1, qhat_t8)
    
    t0 <- proc.time()[3]; ci <- tryCatch(ci_mj(x, tau, alpha)[c("lower","upper")],
                                         error=function(e) c(NA,NA)); t1 <- proc.time()[3]
    out$MaritzJar <- pack(ci, t0, t1, qhat_t8)
    
    t0 <- proc.time()[3]; ci <- tryCatch(ci_smoothed_boot(x, tau, alpha, B=B, seed = s)[c("lower","upper")],
                                         error=function(e) c(NA,NA)); t1 <- proc.time()[3]
    out$SmBoot <- pack(ci, t0, t1, qhat_t8)
    
    t0 <- proc.time()[3]; ci <- tryCatch(ci_subsampling(x, tau, alpha, B=B, seed = s)[c("lower","upper")],
                                         error=function(e) c(NA,NA)); t1 <- proc.time()[3]
    out$Subsample <- pack(ci, t0, t1, qhat_t8)
    
    t0 <- proc.time()[3]; ci <- tryCatch(ci_m_out_n(x, tau, alpha, B=B, seed = s)[c("lower","upper")],
                                         error=function(e) c(NA,NA)); t1 <- proc.time()[3]
    out$mOutOfn <- pack(ci, t0, t1, qhat_t8)
    
    t0 <- proc.time()[3]; ci <- tryCatch(ci_hsnyblom(x, tau, alpha)[c("lower","upper")],
                                         error=function(e) c(NA,NA)); t1 <- proc.time()[3]
    out$HS_Nyblom <- pack(ci, t0, t1, qhat_t8)
    
    list(measures = out, sample_info = smry)
  }
  
  if (parallel && .Platform$OS.type == "windows"){
    warning("parallel=TRUE 在 Windows 上回退为串行。"); parallel <- FALSE
  }
  
  out_list <- if (parallel) {
    if (!requireNamespace("pbmcapply", quietly = TRUE)) stop("需要包 pbmcapply 才能并行运行")
    pbmcapply::pbmclapply(seq_len(R), runner, mc.cores = 12)
  } else lapply(seq_len(R), runner)
  
  for (r in seq_len(R)){
    res <- out_list[[r]]
    samples[r, ] <- res$sample_info
    for (m in METHODS){
      vec <- res$measures[[m]]
      if (is.null(vec)) vec <- make_na_vec()
      acc[[m]][r,] <- vec
    }
  }
  
  summary <- do.call(rbind, lapply(METHODS, function(m) cbind(method = m, summarise_matrix(acc[[m]]))))
  rownames(summary) <- NULL
  summary <- summary[order(summary$mean_len), ]
  
  list(
    params  = list(R=R, n=n, tau=tau, alpha=alpha, dist=dist, B=B, seed=seed, parallel=parallel,
                   snqesa_engine=snqesa_engine, snqesa_calib=snqesa_calib,
                   snqesa_tail_split=snqesa_tail_split, snqesa_cc=snqesa_cc,
                   snqesa_ridge=if (is.null(snqesa_ridge)) "default:max(n^(-2/3),3/n)" else snqesa_ridge,
                   snqesa_enable=snqesa_enable),
    t_true  = tq,
    summary = summary,
    raw     = acc,
    samples = as.data.frame(samples)
  )
}




## ---------- 12. 示例 ----------


# x <- rsamp(200, 'mixnorm')
# x <- rsamp(100, 'cauchy')
# 
# set.seed(NULL); x <- rnorm(15)
# ans <- compare_ci_methods(x, tau = 0.95, alpha = 0.05, t_true = qnorm(0.95))
# #attr(ans, "sample_info")   # 查看样本概况与点估计
# ans[, c("lower","upper","length","center","qhat_t8","qhat_in_CI","hit_true","true_rel_pos")]

# qnorm(0.95)
# 
# 
# res_norm <- simulate_compare_methods(
#   R = 500, n = 150, tau = 0.05, alpha = 0.05,
#   seed = 1319901, B = 1000,
#   dist = "normal",
#   parallel = TRUE,
#   snqesa_engine = "2d",
#   snqesa_calib = "rstar",       
#   snqesa_tail_split = "minlength",
#   snqesa_cc = "jeffreys",
#   snqesa_ridge = NULL          # 使用默认 max(n^(-2/3), 3/n)
# )
# res_norm$summary
# 
# 
# 
# res_lognormal <- simulate_compare_methods(
#   R = 2000, n = 100, tau = 0.95, alpha = 0.05,
#   dist = 'lognormal',          # 同上：可换 "normal","lognormal","t3","cauchy"
#   B = 1000, seed = 1008611,
#   parallel = TRUE           # Linux/macOS 可设 TRUE（需要 pbmcapply 包）
# )
# 
# round(res_lognormal$summary[, sapply(res_lognormal$summary, is.numeric)], 4)
# res_lognormal$summary
# 
# 
# res_t <- simulate_compare_methods(
#   R = 2000, n = 100, tau = 0.95, alpha = 0.05,
#   dist = 't2',          # 同上：可换 "normal","lognormal","t3","cauchy"
#   B = 1000, seed = 1008611,
#   parallel = TRUE           # Linux/macOS 可设 TRUE（需要 pbmcapply 包）
# )
# res_t$summary
# 
# 
# 
# res_cauchy <- simulate_compare_methods(
#   R = 2000, n = 100, tau = 0.95, alpha = 0.05,
#   dist = 'cauchy',          # 同上：可换 "normal","lognormal","t3","cauchy"
#   B = 1000, seed = 1008611,
#   parallel = TRUE           # Linux/macOS 可设 TRUE（需要 pbmcapply 包）
# )
# res_cauchy$summary
# 
# 
# 
# res_mixnorm <- simulate_compare_methods(
#   R = 2000, n = 100, tau = 0.95, alpha = 0.05,
#   dist = 'mixnorm',          # 同上：可换 "normal","lognormal","t3","cauchy"
#   B = 1000, seed = 1008611,
#   parallel = TRUE           # Linux/macOS 可设 TRUE（需要 pbmcapply 包）
# )
# res_mixnorm$summary
# library(distributional)
# 
# 
# set.seed(NULL)
# n = 15
# x <- rcauchy(n)
# tau = 0.5
# dist = dist_pareto(3,5)
# x = generate(dist, 10)[[1]]
# x = rt(5, 3)
# x = rsamp_mixnorm(10)
# x = rtukey_gh(20)
# 

#qcauchy(tau)
# qt(tau, 3) 1.63
# snqesa_ci(x, tau, alpha,
#                 engine = "1d",
#                 calib  = "rstar",     # 或更保守 "Exact"
#                 tail_split = "equal",
#                 cc = "jeffreys",     # 这就是需要 k 的连续性更正
#                 ridge = max(n^(-2/3), 3/n))


#snqesa_ci_discrete(x, tau = tau, alpha = 0.05, tail_split = "minlength")

# ci_smoothed_boot(x, tau, 0.05, B = 1000)
# wald_kde_ci(x, tau, 0.05)
# c("normal","lognormal","t2","t3","cauchy","mixnorm")
# 
# res_normal_95 <- simulate_compare_methods(
#   R = 500, n = 100, tau = 0.95, alpha = 0.05,
#   seed = 10086123, B = 999,
#   dist = "normal",
#   parallel = TRUE,
#   snqesa_engine = "2d",
#   snqesa_calib = "rstar",       
#   snqesa_tail_split = "equal",
#   snqesa_cc = "none",
#   snqesa_ridge = 0.8          # 使用默认 max(n^(-2/3), 3/n)
# )
# res_normal_95$summary
# 
# summary_to_latex(res_normal_95,
#                  caption = "Coverage and length across methods (\\tau=0.95, \\n=100, \\alpha=0.05).",
#                  label = "tab:res_normal_95")


res_normal_99 <- simulate_compare_methods(
  R = 500, n = 30, tau = 0.7, alpha = 0.05,
  seed = 10086123, B = 999,
  dist = "normal",
  parallel = TRUE,
  snqesa_engine = "2d",
  snqesa_calib = "rstar",       
  snqesa_tail_split = "equal",
  snqesa_cc = "anscombe",
  snqesa_ridge = 0.0         # 使用默认 max(n^(-2/3), 3/n)
)
res_normal_99$summary
#qnorm(0.99) # 2.326348
qnorm(0.7) # 0.5244005
res_normal_99$raw$SNQESA_cont[3, ]

set.seed(10086123 + 3)

xx = rnorm(100)
hist(xx)
quantile(xx, 0.7)
xx[xx>quantile(xx, 0.7)]
#c("jeffreys","half","anscombe","none")

res = snqesa_ci(xx, 0.7, ridge = 0, engine = '2d', 
                calib = 'rstar', tol = 1e-10, 
                tail_split = 'minlength', cc = 'none')
res$lower; res$upper


summary_to_latex(res_normal_99,
                 caption = "Coverage and length across methods (\\tau=0.99, \\n=100, \\alpha=0.05).",
                 label = "tab:res_normal_99")



# res_normal_05 <- simulate_compare_methods(
#   R = 500, n = 100, tau = 0.05, alpha = 0.05,
#   seed = 10086123, B = 999,
#   dist = "normal",
#   parallel = TRUE,
#   snqesa_engine = "2d",
#   snqesa_calib = "rstar",       
#   snqesa_tail_split = "equal",
#   snqesa_cc = "none",
#   snqesa_ridge = 0.8        # 使用默认 max(n^(-2/3), 3/n)
# )
# #res_normal_05$raw$SNQESA_cont
# res_normal_05$summary
# 
# #c("jeffreys","half","anscombe","none")
# 
# summary_to_latex(res_normal_05,
#                  caption = "Coverage and length across methods (\\tau=0.05, \\n=100, \\alpha=0.05).",
#                  label = "tab:res_normal_05")


# res_lognormal_95 <- simulate_compare_methods(
#   R = 500, n = 100, tau = 0.95, alpha = 0.05,
#   seed = 10086123, B = 999,
#   dist = "lognormal",
#   parallel = TRUE,
#   snqesa_engine = "2d",
#   snqesa_calib = "rstar",       
#   snqesa_tail_split = "equal",
#   snqesa_cc = "none",
#   snqesa_ridge = 0.3          # 使用默认 max(n^(-2/3), 3/n)
# )
# res_lognormal_95$summary
# 
# summary_to_latex(res_lognormal_95,
#                  caption = "Coverage and length across methods (\\tau=0.9, \\alpha=0.05).",
#                  label = "tab:res_lognormal_95")




res_lognormal_50 <- simulate_compare_methods(
  R = 500, n = 80, tau = 0.5, alpha = 0.05,
  seed = 10086123, B = 999,
  dist = "lognormal",
  parallel = TRUE,
  snqesa_engine = "2d",
  snqesa_calib = "rstar",       
  snqesa_tail_split = "equal",
  snqesa_cc = "half",
  snqesa_ridge = 0          # 使用默认 max(n^(-2/3), 3/n)
)
res_lognormal_50$summary


qlnorm(0.05)
xx = rlnorm(200)
quantile(xx, 0.05)
#c("jeffreys","half","anscombe","none")
res = snqesa_ci(xx, 0.05, ridge = 0, engine = '2d', 
          calib = 'rstar', tol = 1e-10, 
          tail_split = 'equal', cc = 'none')
res$lower; res$upper


summary_to_latex(res_lognormal_05,
                 caption = "Coverage and length across methods (\\tau=0.05,\\n=100, \\alpha=0.05).",
                 label = "tab:res_lognormal_05")





res_t_95 <- simulate_compare_methods(
  R = 500, n = 100, tau = 0.95, alpha = 0.05,
  seed = 10086123, B = 999,
  dist = "t2",
  parallel = TRUE,
  snqesa_engine = "2d",
  snqesa_calib = "rstar",       
  snqesa_tail_split = "equal",
  snqesa_cc = "none",
  snqesa_ridge = 0.4          # 使用默认 max(n^(-2/3), 3/n)
)
res_t_95$summary

summary_to_latex(res_t_95,
                 caption = "Coverage and length across methods (\\tau=0.95, \\alpha=0.05).",
                 label = "tab:res_t_95")



res_t_50 <- simulate_compare_methods(
  R = 1000, n = 50, tau = 0.50, alpha = 0.05,
  seed = 10086123, B = 999,
  dist = "t2",
  parallel = TRUE,
  snqesa_engine = "2d",
  snqesa_calib = "rstar",       
  snqesa_tail_split = "equal",
  snqesa_cc = "half",
  snqesa_ridge = 0          # 使用默认 max(n^(-2/3), 3/n)
)
res_t_50$summary

xx = rt(50, 2)
hist(xx)

summary_to_latex(res_t_50,
                 caption = "Coverage and length across methods (\\tau=0.50,\\n = 100, \\alpha=0.05).",
                 label = "tab:res_t_50")




# res_cauchy_95 <- simulate_compare_methods(
#   R = 500, n = 100, tau = 0.95, alpha = 0.05,
#   seed = 10086123, B = 999,
#   dist = "cauchy",
#   parallel = TRUE,
#   snqesa_engine = "2d",
#   snqesa_calib = "rstar",       
#   snqesa_tail_split = "equal",
#   snqesa_cc = "none",
#   snqesa_ridge = 0.2          # 使用默认 max(n^(-2/3), 3/n)
# )
# res_cauchy_95$summary
# 
# summary_to_latex(res_cauchy_95,
#                  caption = "Coverage and length across methods (\\tau=0.9, \\alpha=0.05).",
#                  label = "tab:res_cauchy_95")


# res_cauchy_50 <- simulate_compare_methods(
#   R = 1000, n = 50, tau = 0.50, alpha = 0.05,
#   seed = 10086123, B = 999,
#   dist = "cauchy",
#   parallel = TRUE,
#   snqesa_engine = "2d",
#   snqesa_calib = "rstar",       
#   snqesa_tail_split = "equal",
#   snqesa_cc = "none",
#   snqesa_ridge = 0          # 使用默认 max(n^(-2/3), 3/n)
# )
# res_cauchy_50$summary
# 
# summary_to_latex(res_cauchy_50,
#                  caption = "Coverage and length across methods (\\tau=0.5, \\n = 50, \\alpha=0.05).",
#                  label = "tab:res_cauchy_50")

# res_beta_95 <- simulate_compare_methods(
#   R = 1000, n = 100, tau = 0.95, alpha = 0.05,
#   seed = 10086123, B = 999,
#   dist = "beta",
#   parallel = TRUE,
#   snqesa_engine = "2d",
#   snqesa_calib = "rstar",       
#   snqesa_tail_split = "equal",
#   snqesa_cc = "half",
#   snqesa_ridge = 0.5          # 使用默认 max(n^(-2/3), 3/n)
# )
# res_beta_95$summary
# 
# summary_to_latex(res_beta_95,
#                  caption = "Coverage and length across methods (\\tau=0.95,\\n = 100, \\alpha=0.05).",
#                  label = "tab:res_beta_95")


# res_mixnorm_95 <- simulate_compare_methods(
#   R = 500, n = 100, tau = 0.95, alpha = 0.05,
#   seed = 10086123, B = 999,
#   dist = "mixnorm",
#   parallel = TRUE,
#   snqesa_engine = "2d",
#   snqesa_calib = "rstar",       
#   snqesa_tail_split = "equal",
#   snqesa_cc = "half",
#   snqesa_ridge = 0.2          # 使用默认 max(n^(-2/3), 3/n)
# )
# res_mixnorm_95$summary
# 
# summary_to_latex(res_mixnorm_95,
#                  caption = "Coverage and length across methods (\\tau=0.9, \\alpha=0.05).",
#                  label = "tab:res_mixnorm_95")





res_exp_95 <- simulate_compare_methods(
  R = 1000, n = 100, tau = 0.95, alpha = 0.05,
  seed = 10086123, B = 999,
  dist = "exp",
  parallel = TRUE,
  snqesa_engine = "2d",
  snqesa_calib = "rstar",       
  snqesa_tail_split = "equal",
  snqesa_cc = "half",
  snqesa_ridge = 0          # 使用默认 max(n^(-2/3), 3/n)
)
res_exp_95$summary

summary_to_latex(res_exp_95,
                 caption = "Coverage and length across methods (\\tau=0.95,\\n = 100, \\alpha=0.05).",
                 label = "tab:res_exp_95")
















# 
# oracle_row <- snqesa_oracle_row(res_normal_95, target_coverage = 0.95,
#                                 baseline = "SNQESA_cont", anchor = "Exact",
#                                 count_anchor_time = FALSE)
# 
# cat(oracle_to_latex_row(oracle_row))
# 
# oracle_row = snqesa_oracle_row(res_t_95, target_coverage = 0.95,
#                   baseline = "SNQESA_cont", anchor = "Exact",
#                   count_anchor_time = FALSE)
# 
# 
# cat(oracle_to_latex_row(oracle_row))
