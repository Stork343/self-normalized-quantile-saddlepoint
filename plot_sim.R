# ===== 依赖 =====
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(forcats)
})

# 把 raw 矩阵摊平为 tidy long
as_tidy_runs <- function(res){
  ms <- names(res$raw)
  bind_rows(lapply(ms, function(m){
    as.data.frame(res$raw[[m]]) |>
      mutate(method = m, .before = 1)
  }))
}

# ------- 1) 覆盖率点图（带二项误差棒） -------
plot_coverage <- function(res, nominal = 1 - res$params$alpha){
  df <- res$summary |>
    mutate(
      lower95 = pmax(0, coverage - 1.96*cover_se),
      upper95 = pmin(1, coverage + 1.96*cover_se),
      method  = fct_reorder(method, coverage)
    )
  ggplot(df, aes(x = method, y = coverage)) +
    geom_hline(yintercept = nominal, linetype = "dashed") +
    geom_errorbar(aes(ymin = lower95, ymax = upper95), width = 0.2) +
    geom_point(size = 2) +
    coord_flip() +
    labs(
      title = sprintf("Coverage vs nominal (n=%d, tau=%.2f, dist=%s)",
                      res$params$n, res$params$tau, res$params$dist),
      x = NULL, y = "Empirical coverage"
    )
}

# ------- 2) 区间长度分布（箱线，log 轴） -------
plot_length_box <- function(res){
  df <- as_tidy_runs(res) |>
    filter(is.finite(length), length >= 0) |>
    mutate(method = fct_reorder(method, length, .fun = median, na.rm = TRUE))
  ggplot(df, aes(x = method, y = length)) +
    geom_boxplot(outlier.alpha = 0.2) +
    scale_y_continuous(trans = "log10") +
    coord_flip() +
    labs(title = "Distribution of CI lengths", x = NULL, y = "Length (log scale)")
}

# ------- 3) 失配方向分解 -------
plot_miss_breakdown <- function(res){
  df <- as_tidy_runs(res) |>
    summarise(
      lower = mean(lower_miss, na.rm = TRUE),
      upper = mean(upper_miss, na.rm = TRUE),
      .by = method
    ) |>
    pivot_longer(c(lower, upper), names_to = "side", values_to = "rate") |>
    mutate(
      side = factor(side, levels = c("lower","upper"), labels = c("Lower-miss","Upper-miss")),
      method = fct_reorder(method, rate, .fun = function(z) -sum(z)) # 上下合计降序
    )
  ggplot(df, aes(x = method, y = rate, fill = side)) +
    geom_col(position = "stack") +
    coord_flip() +
    labs(title = "Under-coverage decomposition", x = NULL, y = "Rate") +
    guides(fill = guide_legend(title = NULL))
}

# ------- 4) Interval Score（越低越好） -------
plot_is_box <- function(res){
  df <- as_tidy_runs(res) |>
    filter(is.finite(is)) |>
    mutate(method = fct_reorder(method, is, .fun = median, na.rm = TRUE))
  ggplot(df, aes(method, is)) +
    geom_boxplot(outlier.alpha = 0.2) +
    coord_flip() +
    labs(title = "Interval Score by method", x = NULL, y = "Interval Score")
}

# ------- 5) 效率前沿：长度 vs 欠覆盖 -------
plot_efficiency_frontier <- function(res, nominal = 1 - res$params$alpha){
  df <- res$summary |>
    transmute(
      method,
      mean_len,
      under = pmax(0, nominal - coverage),
      fail_rate
    )
  ggplot(df, aes(x = mean_len, y = under, label = method, size = pmax(fail_rate, 1e-3))) +
    geom_point(alpha = 0.85) +
    ggrepel::geom_text_repel(min.segment.length = 0, box.padding = 0.2, size = 3) +
    labs(
      title = "Efficiency frontier: shorter vs under-coverage",
      x = "Mean CI length", y = "Under-coverage (nominal − empirical)", size = "Fail rate"
    ) +
    theme(legend.position = "right")
}
# 若不想依赖 ggrepel，可删掉 geom_text_repel 并改用 geom_text(check_overlap=TRUE)

# ------- 6) 真分位相对位置的均匀性 -------
plot_true_relpos <- function(res){
  df <- as_tidy_runs(res) |>
    filter(is.finite(relpos), relpos >= 0, relpos <= 1)
  ggplot(df, aes(relpos)) +
    geom_histogram(bins = 20) +
    geom_hline(yintercept = nrow(df)/20, linetype = "dotted") +
    facet_wrap(~method, scales = "free_y") +
    labs(title = "Calibration via true-relative-position U ∼ U(0,1)", x = "U = (θ−L)/(U−L)", y = "Count")
}

# ------- 7) 运行时间与失败率 -------
plot_time_fail <- function(res){
  df <- res$summary |>
    mutate(method = fct_reorder(method, mean_time))
  p1 <- ggplot(df, aes(method, mean_time)) +
    geom_col() + coord_flip() +
    labs(title = "Mean runtime by method", x = NULL, y = "Seconds (mean)")
  p2 <- ggplot(df, aes(fct_reorder(method, fail_rate), fail_rate)) +
    geom_col() + coord_flip() +
    labs(title = "Failure rate by method", x = NULL, y = "Failure rate")
  list(time = p1, fail = p2)
}

# ------- 8) 多实验结果合并（跨分布/样本量） -------
bind_summaries <- function(named_list){
  bind_rows(lapply(names(named_list), function(nm){
    res <- named_list[[nm]]
    ss  <- res$summary
    ss$scenario <- nm
    ss$n   <- res$params$n
    ss$tau <- res$params$tau
    ss$dist<- res$params$dist
    ss
  }))
}



# 单一场景
p1 <- plot_coverage(res_cauchy_5)
p2 <- plot_length_box(res_cauchy_5)
p3 <- plot_miss_breakdown(res_cauchy_5)
p4 <- plot_is_box(res_cauchy_5)
p5 <- plot_efficiency_frontier(res_cauchy_5)
cal <- plot_true_relpos(res_cauchy_5)
ptf <- plot_time_fail(res_cauchy_5)  # 返回两个图：ptf$time 和 ptf$fail

# 多场景（比如你已经算了 res_norm, res_cauchy_5, res_mixnorm ...）
all_sum <- bind_summaries(list(
  lognormal  = res_lognormal_5,
  cauchy     = res_cauchy_5
  # mixnorm = res_mixnorm, ...
))
# 例如：跨分布的覆盖率对比（按方法分面）
ggplot(all_sum, aes(x = dist, y = coverage, group = method, color = method)) +
  geom_point(position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = pmax(0, coverage - 1.96*cover_se),
                    ymax = pmin(1, coverage + 1.96*cover_se)),
                width = 0.4, position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 1 - unique(all_sum$alpha), linetype = "dashed") +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(title = "Coverage across distributions", x = "Distribution", y = "Coverage")


res_cauchy_5$summary






## -------- [新] 实证数据绘图函数 (核心) --------
plot_empirical_data <- function(r, dates, main_title = "Empirical Data Analysis") {
  
  # 确保 r 是数值型
  # r <- as.numeric(r)
  
  # 记录并设置绘图参数：2x2 布局
  # par(mfrow = ...) 会被 run_empirical_plots 控制
  
  # --- 1. 收益率时序图 (Returns Time Series) ---
  plot(dates, r, type = 'l', col = "steelblue", 
       xlab = "Date", ylab = "Returns", main = '(a)')
  grid()
  
  # --- 2. 收益率直方图与密度 (Histogram & Density) ---
  r_mean <- mean(r, na.rm = TRUE)
  r_sd   <- sd(r, na.rm = TRUE)
  
  hist(r, breaks = 100, freq = FALSE, col = "gray", border = "white",
       xlab = "Returns", ylab = "Density", main = '(b)')
  
  # 叠加核密度估计
  tryCatch({
    lines(density(r, na.rm = TRUE), col = "red", lwd = 2)
  }, error = function(e) {})
  
  # 叠加正态分布曲线
  curve(dnorm(x, mean = r_mean, sd = r_sd), 
        add = TRUE, col = "blue", lty = 2, lwd = 2)
  legend("topright", legend = c("KDE", "Normal"), 
         col = c("red", "blue"), lty = c(1, 2), lwd = 2, bty = "n", cex = 0.8)
  
  # --- 3. 正态 Q-Q 图 (Normal Q-Q Plot) ---
  qqnorm(r, main = '(c)',
         xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
         pch = 19, cex = 0.5, col = rgb(0,0,0,0.3) )
  qqline(r, col = "red", lwd = 2)
  grid()
  
  # --- 4. 平方收益率的自相关 (ACF of Squared Returns) ---
  acf(r^2, na.action = na.pass, lag.max = 50, plot = TRUE, main = '(d)')
  
  # 添加总标题
  #title(main_title, outer = TRUE, cex.main = 1.5)
  
  invisible()
}



## -------- [新] 交互式绘图与保存 (控制器) --------
run_empirical_plots <- function(r, dates, 
                                main_title = "Empirical Data Analysis", 
                                filename = "empirical_analysis_plots.pdf") {
  
  # 1. 检查是否在交互式环境
  if (!interactive()) {
    cat("Not in interactive mode. Skipping plots display and save prompt.\n")
    return(invisible(NULL))
  }
  
  # 2. 准备绘图参数并显示
  cat("Displaying empirical data plots...\n")
  old_par <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 0, 0))
  on.exit(par(old_par), add = TRUE) # 确保函数退出时恢复
  
  tryCatch({
    plot_empirical_data(r, dates, main_title)
  }, error = function(e) {
    cat("Error displaying plots:", e$message, "\n")
    return(invisible(NULL))
  })
  
  # 3. 询问用户是否保存
  answer <- ""
  while (!answer %in% c("y", "n")) {
    prompt <- paste0("Plots displayed. Save plots to '", filename, "'? (y/n): ")
    answer <- readline(prompt)
    answer <- tolower(trimws(answer))
    if (answer == "") answer <- "n" # 默认不保存
  }
  
  # 4. 根据回答执行保存
  if (answer == "y") {
    cat("Saving plots to", filename, "...\n")
    ext <- tools::file_ext(filename)
    
    tryCatch({
      # 根据文件扩展名选择设备
      if (ext == "pdf") {
        pdf(filename, width = 10, height = 8)
      } else if (ext == "png") {
        png(filename, width = 1000, height = 800)
      } else {
        stop(paste("Unsupported file type:", ext, ". Please use .pdf or .png"))
      }
      
      # 再次调用绘图函数，这次是绘制到文件
      par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
      plot_empirical_data(r, dates, main_title)
      dev.off() # 关闭文件设备
      
      cat("Successfully saved to", filename, "\n")
      
    }, error = function(e) {
      cat("Error saving file:", e$message, "\n")
      if (names(dev.cur()) != "null device") {
        dev.off() # 确保关闭失败的设备
      }
    })
    
  } else {
    cat("Plots not saved.\n")
  }
  
  invisible(NULL)
}


run_empirical_plots(spx$r, spx$date, 
                    main_title = "S&P 500 Daily Returns Analysis",
                    filename = "spx_empirical_plots.pdf")

