library(dplyr); library(tidyr); library(ggplot2); library(zoo)

# 合并各模型：date, model, VaR, CI_low, CI_up, ES，并拼上 L_next 与超越指示 ex
tidy_var_data <- function(results_list, r_full, dates_full, nwin){
  stopifnot(exists("make_L_next"), exists("make_aligned_dates"))
  L_next <- make_L_next(r_full, nwin)
  aligned_dates <- as.Date(make_aligned_dates(dates_full, nwin))
  base_df <- data.frame(date = aligned_dates, L_next = as.numeric(L_next))
  
  combined <- base_df
  for (m in names(results_list)) {
    ser <- results_list[[m]]$series
    if (is.null(ser)) next
    ser <- ser %>% mutate(date = as.Date(date)) %>% select(date, VaR, CI_low, CI_up, ES)
    names(ser) <- c("date",
                    paste0("VaR_", m),
                    paste0("CI_low_", m),
                    paste0("CI_up_", m),
                    paste0("ES_", m))
    combined <- left_join(combined, ser, by = "date")
  }
  
  long <- pivot_longer(
    combined,
    cols = matches("^(VaR|CI_low|CI_up|ES)_"),
    names_to = c(".value","model"),
    names_pattern = "(VaR|CI_low|CI_up|ES)_(.*)"
  ) %>%
    arrange(model, date) %>%
    mutate(model = factor(model, levels = names(results_list)),
           width = CI_up - CI_low,
           ex = as.integer(L_next > VaR))
  long
}

df <- tidy_var_data(results_list, spx$r, spx$date, 250)



# title = "Exceedance Heatmap (L_next > VaR)"
ggplot(df, aes(x = date, y = model, fill = factor(ex))) +
  geom_tile(height = 0.9) +
  scale_fill_manual(values = c("0" = "#f0f0f0", "1" = "#D55E00"),
                    name = "Exceed", labels = c("No","Yes")) +
  labs(x = "Date", y = "Model") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())




tau <- 0.99; win <- 60
df_roll <- df %>% group_by(model) %>%
  arrange(date,.by_group=TRUE) %>%
  mutate(hit_rate = zoo::rollapply(ex, width = win, FUN = mean, align="right", fill = NA_real_))

ggplot(df_roll, aes(date, hit_rate, color = model)) +
  geom_hline(yintercept = 1 - tau, linetype = "dashed") +
  geom_line(linewidth = 0.7, na.rm = TRUE) +
  labs(title = sprintf("Rolling Hit Rate (window = %d)", win),
       y = "Exceedance Rate", x = "Date") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")




df_sev <- df %>% filter(ex == 1, is.finite(VaR), is.finite(L_next)) %>%
  mutate(severity = (L_next - VaR) / pmax(abs(VaR), 1e-12))

ggplot(df_sev, aes(date, severity)) +
  geom_segment(aes(xend = date, y = 0, yend = severity), linewidth = 0.5, color = "grey50") +
  geom_point(size = 1.2, color = "#D55E00") +
  facet_wrap(~ model, ncol = 2, scales = "free_y") +
  labs(title = "Exceedance Severity ( (L_next - VaR)/|VaR| )",
       y = "Severity", x = "Date") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())



ggplot(df %>% filter(is.finite(width)),
       aes(model, width, fill = model)) +
  geom_violin(trim = TRUE, alpha = 0.15, color = NA) +
  geom_boxplot(width = 0.75, outlier.size = 0.7, alpha = 0.7) +
  labs(x = "Model", y = "CI Width") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", panel.grid.minor = element_blank())



library(ggrepel)

# 取各模型总耗时（若没有就 NA）
time_tbl <- tibble(
  model = names(results_list),
  time_sec = sapply(results_list, function(x) x$time_sec %||% NA_real_)
)

sum_tbl <- df %>%
  group_by(model) %>%
  summarize(mean_width = mean(width, na.rm = TRUE),
            hit_rate   = mean(ex, na.rm = TRUE),
            .groups = "drop") %>%
  left_join(time_tbl, by = "model")

ggplot(sum_tbl, aes(mean_width, hit_rate)) +
  geom_hline(yintercept = 1 - tau, linetype = "dashed") +
  geom_point(aes(size = time_sec), alpha = 0.7) +
  ggrepel::geom_text_repel(aes(label = model), min.segment.length = 0) +
  scale_size_continuous(name = "Time (s)", range = c(2,8)) +
  labs(title = "Width–Hit Rate–Time Trade-off",
       x = "Mean CI Width", y = "Exceedance Rate") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())


df_es <- df %>% filter(ex == 1, is.finite(ES), ES > 0)
ggplot(df_es, aes((L_next/ES), fill = model)) +
  geom_histogram(position = "identity", alpha = 0.25, bins = 40) +
  facet_wrap(~ model, ncol = 2, scales = "free_y") +
  labs(title = "Exceedance: L_next / ES (Shouldn't Systematically > 1)",
       x = "L_next / ES", y = "Count") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")
