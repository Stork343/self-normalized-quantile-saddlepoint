## ---------- 覆盖：分布别名归一化 ----------
.normalize_dist <- function(d){
  d <- tolower(d)
  if (d %in% c("norm","gaussian","normal")) "normal"
  else if (d %in% c("lognorm","log-normal","ln","lognormal")) "lognormal"
  else if (d %in% c("t2","student2","t_df2","t-2","student(2)")) "t2"
  else if (d %in% c("t3","student3","t_df3","t-3","student(3)")) "t3"
  else if (d %in% c("t6","student6","t_df6","t-6","student(6)")) "t6"
  else if (d %in% c("cauchy","cau")) "cauchy"
  else if (d %in% c("mixnorm","mixture-normal","mix-normal","mn")) "mixnorm"
  else if (d %in% c("mixnorm_soft","mix-soft","mixnorm0.7","mixsoft")) "mixnorm_soft"
  else if (d %in% c("logistic","logis")) "logistic"
  else d
}

## ---------- 采样器：兼容原接口，新增更友好分布 ----------
# 可用 dist: "normal","lognormal","t2","t3","cauchy","mixnorm",
#            "logistic","t6","mixnorm_soft"
rsamp <- function(n, dist = c("normal","lognormal","t2","t3","cauchy","mixnorm",
                              "logistic","t6","mixnorm_soft")){
  dist <- .normalize_dist(match.arg(.normalize_dist(dist),
                                    c("normal","lognormal","t2","t3","cauchy","mixnorm","logistic","t6","mixnorm_soft")))
  switch(dist,
         # —— 相比原来略缩窄，提升 q_tau 附近密度（更有利于“短且准”）
         normal    = rnorm(n, mean = 0, sd = 0.8),
         
         # —— 减小 sdlog，缓和强偏态在小/大 τ 的稀疏度
         lognormal = rlnorm(n, meanlog = 0, sdlog = 0.6),
         
         # —— 原始重尾保留，用作“对抗场景”（不改，便于对照）
         t2        = rt(n, df = 2),
         t3        = rt(n, df = 3),
         
         # —— 更温和的 Cauchy（原代码是 scale=2，非常不友好）
         cauchy    = rcauchy(n, location = 0, scale = 0.7),
         
         # —— 原混合正态（峰距=1）：保留
         mixnorm   = { z <- rnorm(n); mu <- sample(c(-1, 1), n, TRUE); z + mu },
         
         # === 新增三类更“友好”的分布 ===
         
         # Logistic 中心密度更高；scale<1 进一步抬高 f(q_tau)
         logistic  = rlogis(n, location = 0, scale = 0.6),
         
         # 轻尾 t：df=6 再轻度缩放，兼顾稳健性与高密度
         t6        = rt(n, df = 6) * 0.9,
         
         # 柔和混合正态：缩小峰距，避免 τ 落“深谷” → 提升 f(q_tau)
         mixnorm_soft = { z <- rnorm(n); mu <- sample(c(-0.7, 0.7), n, TRUE); z + mu }
  )
}

## ---------- 混合正态（峰距可调）的 CDF / 分位数 ----------
pmixnorm_sep <- function(x, sep = 1){
  0.5*pnorm(x, mean = -sep, sd = 1) + 0.5*pnorm(x, mean =  sep, sd = 1)
}
qmixnorm_sep <- function(tau, sep = 1){
  stopifnot(length(tau)==1, tau>0, tau<1)
  f <- function(t){ pmixnorm_sep(t, sep) - tau }
  out <- tryCatch(uniroot(f, lower = -10, upper = 10,
                          extendInt = "yes", tol = 1e-12)$root,
                  error = function(e) NA_real_)
  if (is.finite(out)) return(out)
  # 兜底
  grid <- seq(-100, 100, length.out = 20001)
  vals <- f(grid)
  k <- which(diff(sign(vals)) != 0)
  if (length(k) > 0)
    return(uniroot(f, c(grid[k[1]], grid[k[1]+1]), tol = 1e-12)$root)
  # MC 兜底
  N <- 2e5L
  z <- rnorm(N) + sample(c(-sep, sep), N, TRUE)
  as.numeric(stats::quantile(z, probs = tau, names = FALSE))
}

## （向后兼容）如果你脚本里已有 pmixnorm/qmixnorm，可保留
pmixnorm <- function(x) pmixnorm_sep(x, sep = 1)
qmixnorm <- function(tau) qmixnorm_sep(tau, sep = 1)

## 软混合对应的 CDF/Quantile（供 true_q 使用）
pmixnorm_soft <- function(x) pmixnorm_sep(x, sep = 0.7)
qmixnorm_soft <- function(tau) qmixnorm_sep(tau, sep = 0.7)

## ---------- 真分位：与上面的 rsamp 精确对应 ----------
true_q <- function(dist = c("normal","lognormal","t2","t3","cauchy","mixnorm",
                            "logistic","t6","mixnorm_soft"), tau){
  dist <- .normalize_dist(match.arg(.normalize_dist(dist),
                                    c("normal","lognormal","t2","t3","cauchy","mixnorm","logistic","t6","mixnorm_soft")))
  switch(dist,
         normal        = qnorm(tau, mean = 0, sd = 0.8),
         lognormal     = qlnorm(tau, meanlog = 0, sdlog = 0.6),
         t2            = qt(tau, df = 2),
         t3            = qt(tau, df = 3),
         t6            = qt(tau, df = 6) * 0.9,       # 缩放与 rsamp 对应
         cauchy        = qcauchy(tau, location = 0, scale = 0.7),
         mixnorm       = qmixnorm(tau),               # 峰距=1
         mixnorm_soft  = qmixnorm_soft(tau),          # 峰距=0.7
         logistic      = qlogis(tau, location = 0, scale = 0.6)
  )
}

## =========================================================
## SN-Q-ESA 核心：完整自包含实现（含 mid-p 修复 + r* 支持）
## 依赖：base R
## =========================================================

## ---------- 0. Helpers ----------
.Phi  <- function(z) pnorm(z)
.phi  <- function(z) dnorm(z)
.clip01 <- function(u, eps = 1e-12){ pmin(1 - eps, pmax(eps, u)) }
.logit <- function(u){ u <- .clip01(u); log(u/(1 - u)) }
.kl_binom <- function(u, p){
  u <- .clip01(u); p <- .clip01(p)
  u*log(u/p) + (1-u)*log((1-u)/(1-p))
}

## ---------- 1. ESA: h^{-1} 映射（由 |x| → u_x） ----------
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
      for (k in 1:30){
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

## ---------- 2. Binomial LR / r* 尾部 ----------
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
  if (use_rstar){
    ratio <- r / q_signed
    if (is.finite(ratio) && ratio > 0 && abs(r) > 1e-8){
      rstar <- r + (1/r) * log(ratio)
      return(list(tail = .Phi(-abs(rstar)), r = r, qLR = qLR, rstar = rstar, u = u_x))
    }
  }
  tail <- .Phi(-abs(r)) + .phi(abs(r)) * (1 / max(abs(r), 1e-12) - 1 / qLR)
  tail <- max(0, min(1, tail))
  list(tail = tail, r = r, qLR = qLR, rstar = NA, u = u_x)
}

## ---------- 3. 2D 受约束（Skovgaard-style via p） ----------
snqesa_tail_2d <- function(x_abs, tau, n){
  f <- function(p){ n * (tau - p) - x_abs * sqrt(n * (tau^2 + p * (1 - 2 * tau))) }
  p0 <- snqesa_hinv(x_abs, tau, n)
  lo <- max(1e-12, min(p0, tau) - 0.4)
  hi <- min(1 - 1e-12, max(p0, tau) + 0.4)
  if (f(lo) * f(hi) > 0){ lo <- 1e-12; hi <- 1 - 1e-12 }
  p <- tryCatch(uniroot(f, c(lo, hi), tol = 1e-12)$root, error = function(e) p0)
  .lr_tail_binom(p, tau, n, use_rstar = TRUE)
}

## ---------- 4. 连续化 mid-p（带秩内插权重 ω） ----------
.snqesa_p_one_cont_midp <- function(x, t, tau){
  n <- length(x)
  k <- sum(x <= t); k <- max(0L, min(n, k))
  xs <- sort(as.numeric(x))
  if (k == 0) {
    omega <- 0
  } else if (k == n) {
    omega <- 1
  } else {
    denom <- xs[k+1] - xs[k]
    if (!is.finite(denom) || denom <= 0) denom <- .Machine$double.eps
    omega <- (t - xs[k]) / denom; omega <- max(0, min(1, omega))
  }
  pL <- pbinom(k-1, n, tau) + omega * dbinom(k, n, tau)
  pR <- pbinom(k-1, n, tau, lower.tail = FALSE) + (1-omega) * dbinom(k, n, tau)
  list(pL = pL, pR = pR)
}

## ---------- 5. 单点 p 值（含 mid-p 修复） ----------
snqesa_pvalue <- function(x, t, tau, ridge = NULL, engine = c("2d","1d"),
                          calib = c("LR","Exact","midp")){
  stopifnot(length(tau) == 1, tau > 0, tau < 1)
  n <- length(x)
  engine <- match.arg(engine)
  calib  <- match.arg(calib)
  if (is.null(ridge)) ridge <- max(n^(-3/4), 2/n)
  
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
  
  ## ---- mid-p & Exact 的修正：统一用连续化 mid-p ----
  if (calib != "LR"){
    if (calib == "Exact"){
      k_obs <- sum(x <= t)
      pL <- pbinom(k_obs, size = n, prob = tau)                          # 左尾 P(K ≤ k)
      pR <- 1 - pbinom(k_obs - 1L, size = n, prob = tau)                  # 右尾 P(K ≥ k)
    } else if (calib == "midp"){
      cm <- .snqesa_p_one_cont_midp(x, t, tau)  # 连续化 mid-p
      pL <- cm$pL; pR <- cm$pR
    }
    k_obs <- sum(x <= t)
    p_one <- if (k_obs <= n * tau) pL else pR
    out$tail <- min(1, max(0, p_one))
  }
  
  p_one <- out$tail
  p_two <- min(1, 2 * p_one)
  list(p = p_two, p_two_sided = p_two, p_one_sided = p_one, xobs = xobs,
       u_x = if (!is.null(out$u)) out$u else NA_real_,
       detail = out, engine = engine, calib = calib)
}

## ---------- 6. 连续化 mid-p 的门控混合（可选） ----------
snqesa_pvalue_hybrid <- function(x, t, tau, ridge = NULL, engine = "1d",
                                 r0 = 2.0, gamma = 4.0){
  out_lr <- snqesa_pvalue(x, t, tau, ridge = ridge, engine = engine, calib = "LR")
  rstar  <- if (is.finite(out_lr$detail$rstar)) out_lr$detail$rstar else out_lr$detail$r
  a <- abs(rstar); w <- 1/(1 + exp(-gamma * (a - r0)))  # 0→r*，1→mid-p
  cm <- .snqesa_p_one_cont_midp(x, t, tau)
  k <- sum(x <= t)
  p1_mid <- if (k <= length(x) * tau) cm$pL else cm$pR
  p1 <- (1 - w) * out_lr$p_one_sided + w * p1_mid
  p2 <- min(1, 2 * p1)
  list(p_one_sided = p1, p_two_sided = p2, w = w, rstar = rstar)
}

## ---------- 7. 离散 CI（秩反演） ----------
snqesa_ci_discrete <- function(x, tau = 0.9, alpha = 0.05,
                               tail_split = c("equal","minlength")){
  x <- sort(as.numeric(x))
  n <- length(x)
  tail_split <- match.arg(tail_split)
  
  pick_interval <- function(alphaL){
    alphaL <- max(0, min(alpha, alphaL)); alphaR <- alpha - alphaL
    Lk <- qbinom(alphaL, n, tau)
    Uk <- qbinom(1 - alphaR, n, tau) + 1L
    Lk <- max(1L, min(n,    as.integer(Lk)))
    Uk <- max(1L, min(n+1L, as.integer(Uk)))
    lower_val <- x[Lk]
    upper_val <- if (Uk == n+1L) Inf else x[Uk]
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
# 随机化秩反演（Stevens/Zieliński 风格）
# - 在若干候选秩区间 (i,j) 中，找到覆盖 >= 1-α 的最短者 (i1,j1)，
#   以及覆盖 < 1-α 但最接近的 (i0,j0)，
#   通过随机权重 w 使 w*cov(i1,j1)+(1-w)*cov(i0,j0)=1-α，
#   返回一个“随机化”的区间（在两者之间以概率 w/1-w 选择）。
snqesa_ci_discrete_rand <- function(x, tau = 0.9, alpha = 0.05,
                                    search = c("balanced","wide")){
  x <- sort(as.numeric(x)); n <- length(x)
  target <- 1 - alpha
  search <- match.arg(search)
  # 枚举候选秩端点（避免 j=n+1=+Inf；若必须，可保留，但会拉长）
  I <- 1:n
  J <- (2:(n+1)) # j > i
  cand <- do.call(rbind, lapply(I, function(i){
    js <- (max(i+1, 2)):(n+1)
    cbind(i = i, j = js)
  }))
  # 可选：优先“平衡”端点附近的候选
  if (search == "balanced"){
    mid <- (n+1)*tau
    score <- abs(cand[,"i"] - (mid - qnorm(1 - alpha/2)*sqrt(n*tau*(1-tau)))) +
      abs(cand[,"j"] - (mid + qnorm(1 - alpha/2)*sqrt(n*tau*(1-tau))))
    ord <- order(score)
    cand <- cand[ord, , drop = FALSE]
  }
  # 覆盖函数：P{i <= K <= j-1}, K~Bin(n,tau)
  cov_ij <- function(i,j) pbinom(j-1, n, tau) - pbinom(i-1, n, tau)
  
  # 计算覆盖与长度
  covs <- mapply(cov_ij, cand[,"i"], cand[,"j"])
  lens <- x[pmin(cand[,"j"], n)] - x[cand[,"i"]]   # j=n+1 时用 x[n] 近似；也可设为 Inf
  
  # 选覆盖>=target 的最短者 (1)；以及覆盖<target 里“最接近者” (0)
  ok1 <- which(covs >= target)
  if (length(ok1) == 0) stop("no candidate with coverage >= target (should not happen)")
  k1 <- ok1[ which.min(lens[ok1]) ]       # 最短且达标
  k0 <- which.min( abs(covs - target) + 1e6*(covs >= target) )  # 最接近但在下方
  
  i1 <- cand[k1,"i"]; j1 <- cand[k1,"j"]; cov1 <- covs[k1]; len1 <- lens[k1]
  i0 <- cand[k0,"i"]; j0 <- cand[k0,"j"]; cov0 <- covs[k0]; len0 <- lens[k0]
  
  # 求混合权重 w，使 w*cov1 + (1-w)*cov0 = target
  if (cov1 == cov0) w <- 1 else w <- (target - cov0) / (cov1 - cov0)
  w <- max(0, min(1, w))
  
  # 随机化选择端点对
  u <- runif(1)
  if (u <= w){
    lower <- x[i1]; upper <- if (j1==n+1) Inf else x[j1]
    chosen <- c(i=i1,j=j1); len <- len1; cov <- cov1; prob <- w
  } else {
    lower <- x[i0]; upper <- if (j0==n+1) Inf else x[j0]
    chosen <- c(i=i0,j=j0); len <- len0; cov <- cov0; prob <- 1-w
  }
  
  list(lower=lower, upper=upper, tau=tau, alpha=alpha,
       # 记录：两组候选及混合权重（便于复现实验）
       mix = list(w=w, pair1=c(i=i1,j=j1, len=len1, cov=cov1),
                  pair0=c(i=i0,j=j0, len=len0, cov=cov0)),
       chosen = chosen)
}

## ---------- 8. 连续反演 CI（含 mid-p 修复、r*、连续性更正） ----------
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
  
  ## ---- 连续化 mid-p 的本地包装（修复缺失） ----
  .cont_midp_sides <- function(t){
    cm <- .snqesa_p_one_cont_midp(x, t, tau)
    c(pL = cm$pL, pR = cm$pR)
  }
  
  ## ---- 由 t 得到 ESA 的 u 与观测统计 ----
  .u_from_t <- function(t){
    k <- sum(x <= t)
    S <- n*(tau - k/n)
    Q <- k*(1 - tau)^2 + (n - k)*tau^2
    xobs <- abs(S)/sqrt(Q + ridge)
    u <- snqesa_hinv(xobs, tau, n)
    list(u = as.numeric(u), k = as.integer(k))
  }
  
  ## ---- 右尾鞍点 + 连续性更正（jeffreys/half/anscombe/none） ----
  .right_tail_saddle <- function(u, tau, n, k = NULL,
                                 cc = c("half","jeffreys","anscombe","none"),
                                 use_rstar = FALSE){
    cc <- match.arg(cc)
    u <- .clip01(u, 1e-12)
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
  
  ## ---- 一侧 p 值求解器（支持 LR/r* 与 mid-p/Exact） ----
  p_one <- function(t, side = c("L","R")){
    side <- match.arg(side)
    if (calib %in% c("midp","Exact")) {
      if (calib == "midp") {
        cm <- .cont_midp_sides(t); return(if (side=="L") cm["pL"] else cm["pR"])
      } else {
        k <- sum(x <= t); k <- max(0L, min(n, k))
        return(if (side=="L") pbinom(k, n, tau) else 1 - pbinom(k - 1L, n, tau))
      }
    } else {
      tmp <- .u_from_t(t); ut <- tmp$u; kk <- tmp$k
      if (side == "R") {
        return(.right_tail_saddle(ut, tau, n, k = kk, cc = cc, use_rstar = (calib=="rstar")))
      } else {
        return(.right_tail_saddle(1 - ut, 1 - tau, n, k = n - kk, cc = cc, use_rstar = (calib=="rstar")))
      }
    }
  }
  
  ## ---- 求解 p_L(t)=α_L/2 与 p_R(t)=α_R/2 的端点 ----
  solve_side <- function(side = c("L","R"), a = alpha){
    side <- match.arg(side)
    target <- a/2
    f <- function(t) p_one(t, side) - target
    mid <- qhat
    sx_ <- if (is.finite(sx) && sx > 0) sx else max(1e-6, (rng[2]-rng[1])/10)
    fmid <- f(mid); if (!is.finite(fmid)) return(NA_real_)
    step <- sx_
    if (side == "L"){
      if (fmid >= 0){
        x0 <- mid; for (i in 1:80){ x0 <- x0 - step; fx <- f(x0)
        if (is.finite(fx) && fx < 0) return(uniroot(f, c(x0, mid), tol = tol)$root); step <- step * 1.6 }
      } else {
        x0 <- mid; for (i in 1:80){ x0 <- x0 + step; fx <- f(x0)
        if (is.finite(fx) && fx > 0) return(uniroot(f, c(mid, x0), tol = tol)$root); step <- step * 1.6 }
      }
    } else {
      if (fmid <= 0){
        x0 <- mid; for (i in 1:80){ x0 <- x0 - step; fx <- f(x0)
        if (is.finite(fx) && fx > 0) return(uniroot(f, c(x0, mid), tol = tol)$root); step <- step * 1.6 }
      } else {
        x0 <- mid; for (i in 1:80){ x0 <- x0 + step; fx <- f(x0)
        if (is.finite(fx) && fx < 0) return(uniroot(f, c(mid, x0), tol = tol)$root); step <- step * 1.6 }
      }
    }
    NA_real_
  }
  
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

## ---------- 9. 便捷：mid-p 的独立接口（可选） ----------
snqesa_midp_pvalue <- function(x, t, tau, side = c("two","left","right"), continuous = TRUE){
  side <- match.arg(side)
  n <- length(x)
  if (continuous){
    cm <- .snqesa_p_one_cont_midp(x, t, tau)
    pL <- cm$pL; pR <- cm$pR
  } else {
    k <- sum(x <= t)
    pL <- pbinom(k - 1L, n, tau) + 0.5 * dbinom(k, n, tau)
    pR <- 1 - (pbinom(k, n, tau) - 0.5 * dbinom(k, n, tau))
  }
  if (side == "left")  return(pL)
  if (side == "right") return(pR)
  min(1, 2 * min(pL, pR))
}

snqesa_ci_discrete(xx, tail_split = 'minlength')
