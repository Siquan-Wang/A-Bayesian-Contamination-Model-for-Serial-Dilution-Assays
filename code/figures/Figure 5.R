library(tidyverse)
library(rstan)
library(ggplot2)
library(scales)

# ===== 4-parameter function (same as Stan) =====
h_function_4para <- function(x, b1, b2, b3, b4) {
  b1 + b2/(1 + (x/b3)^(-b4))
}

# ===== Standard data (now with 3 extra points at 250, 125, 62.5) =====
standard_concentration <- 250
d_zm_standard <- seq(0, 1, 0.0001)

col1 <- c(2.303, 2.134, 2.471, 2.617, 1.668, 1.886, 1.261, 1.261)
col2 <- c(0.593, 0.652, 0.306, 0.312, 0.166, 0.168, 0.097, 0.096)
col3 <- c(0.071, 0.071, 0.062, 0.060, 0.055, 0.053, 0.050, 0.057)

# original 24 standard y + 3 “promoted” unknowns
y_standard <- c(col1, col2, col3, 2.105, 2.448, 1.952)

# original d for 125-curve
d_base <- c(
  1,1, 1/2,1/2, 1/4,1/4, 1/8,1/8,
  1/16,1/16, 1/32,1/32, 1/64,1/64,
  1/128,1/128, 1/256,1/256, 1/512,1/512,
  1/1024,1/1024, 1/2048,1/2048
)

# Stan used: c(d_base * 1/2, 1, 1/2, 1/4)
d_standard <- c(d_base * 0.5, 1, 1/2, 1/4)  # length 27

# concentration of each standard point (must match x_standard in Stan)
conc_points <- standard_concentration * d_standard

# ===== Read new Stan fit (250 top concentration) =====
fit <- readRDS("model_fit_newdata_no_initial_value_250.rds")
draws <- as.data.frame(fit, pars = "beta")
stopifnot(ncol(draws) == 4)
colnames(draws) <- c("beta_1","beta_2","beta_3","beta_4")

set.seed(1234)
n_draw <- min(4000, nrow(draws))
idx <- sample(seq_len(nrow(draws)), n_draw)

# grid on absolute concentration scale (0–250, exclude 0)
conc_grid <- standard_concentration * d_zm_standard
conc_grid_pos <- conc_grid[conc_grid > 0]

d_y_standard_more <- matrix(NA_real_, nrow = n_draw, ncol = length(conc_grid_pos))
for (i in seq_len(n_draw)) {
  b1 <- draws$beta_1[idx[i]]
  b2 <- draws$beta_2[idx[i]]
  b3 <- draws$beta_3[idx[i]]
  b4 <- draws$beta_4[idx[i]]
  d_y_standard_more[i, ] <- h_function_4para(conc_grid_pos, b1, b2, b3, b4)
}
med <- apply(d_y_standard_more, 2, median)
up  <- apply(d_y_standard_more, 2, quantile, prob = 0.975)
low <- apply(d_y_standard_more, 2, quantile, prob = 0.025)

# breaks from 250 down to ~0.06
breaks_conc <- standard_concentration * 2^-(0:12)

# ===== Unknown data =====
unknown_df <- tribble(
  ~UnknownSample, ~IsContaminated, ~Y,     ~X,
  1, "NO", 2.105, 250,
  2, "NO", 2.448, 125,
  3, "NO", 1.952, 62.5,
  4, "NO", 1.708, 31.25,
  5, "NO", 1.177, 15.625,
  6, "NO", 0.67,  7.8125,
  7, "NO", 0.356, 3.90625,
  8, "NO", 0.174, 1.953125,
  1, "NO", 1.473, 25,
  2, "NO", 1.064, 12.5,
  3, "NO", 0.441, 6.25,
  4, "NO", 0.361, 3.125,
  5, "NO", 0.155, 1.5625,
  6, "NO", 0.093, 0.78125,
  7, "NO", 0.066, 0.390625,
  8, "NO", 0.061, 0.1953125,
  1, "NO", 0.081, 0.25,
  2, "NO", 0.061, 0.125,
  3, "NO", 0.050, 0.0625,
  4, "NO", 0.056, 0.03125,
  5, "NO", 0.049, 0.015625,
  6, "NO", 0.050, 0.0078125,
  7, "NO", 0.049, 0.00390625,
  8, "NO", 0.048, 0.001953125,
  9,  "NO", 0.113, 0.9765625,
  10, "NO", 0.081, 0.48828125,
  11, "NO", 0.067, 0.244140625,
  12, "YES", 1.275, 125,
  13, "YES", 0.724, 125,
  14, "YES", 0.066, 125,
  15, "YES", 0.302, 15.625,
  16, "YES", 0.157, 15.625,
  9,  "NO", 0.055, 0.09765625,
  10, "NO", 0.053, 0.048828125,
  11, "NO", 0.052, 0.024414063,
  12, "YES", 0.692, 12.5,
  13, "YES", 0.524, 12.5,
  14, "YES", 0.049, 12.5,
  15, "YES", 0.120, 1.5625,
  16, "YES", 0.095, 1.5625,
  9,  "NO", 0.052, 0.000976563,
  10, "NO", 0.049, 0.000488281,
  11, "NO", 0.050, 0.000244141,
  12, "YES", 0.057, 0.125,
  13, "YES", 0.056, 0.125,
  14, "YES", 0.048, 0.125,
  15, "YES", 0.052, 0.015625,
  16, "YES", 0.051, 0.015625,
  17, "YES", 0.064, 15.625,
  18, "YES", 0.085, 1.95,
  19, "YES", 0.075, 1.95,
  20, "YES", 0.062, 1.95,
  21, "YES", 0.064, 0.24,
  22, "YES", 0.068, 0.24,
  23, "YES", 0.074, 0.24,
  17, "YES", 0.051, 1.5625,
  18, "YES", 0.060, 0.195,
  19, "YES", 0.057, 0.195,
  20, "YES", 0.051, 0.195,
  21, "YES", 0.055, 0.024,
  22, "YES", 0.052, 0.024,
  23, "YES", 0.052, 0.024,
  17, "YES", 0.062, 0.015625,
  18, "YES", 0.052, 0.00195,
  19, "YES", 0.054, 0.00195,
  20, "YES", 0.051, 0.00195,
  21, "YES", 0.055, 0.00024,
  22, "YES", 0.048, 0.00024,
  23, "YES", 0.049, 0.00024
) |>
  filter(!(UnknownSample %in% c(1, 2, 3))) |>
  mutate(
    group    = factor(if_else(IsContaminated == "YES",
                              "Contaminated", "Uncontaminated"),
                      levels = c("Uncontaminated","Contaminated")),
    sample_f = factor(UnknownSample)
  ) |>
  group_by(X, group) |>
  mutate(
    contam_rank = if_else(group == "Contaminated",
                          row_number(), NA_integer_)
  ) |>
  ungroup()

# legend 顺序
type_levels <- c(
  "Calibration sample",
  "Uncontaminated sample",
  "Contaminated sample, 5% contamination",
  "Contaminated sample, 10% contamination",
  "Contaminated sample, 20% contamination"
)

unknown_df <- unknown_df |>
  mutate(
    type = case_when(
      group == "Uncontaminated" ~ "Uncontaminated sample",
      group == "Contaminated" & contam_rank == 1 ~ "Contaminated sample, 5% contamination",
      group == "Contaminated" & contam_rank == 2 ~ "Contaminated sample, 10% contamination",
      group == "Contaminated" & contam_rank == 3 ~ "Contaminated sample, 20% contamination",
      TRUE ~ NA_character_
    ),
    type = factor(type, levels = type_levels)
  )

# 校准点
standard_df <- tibble(
  X    = conc_points,
  Y    = y_standard,
  type = factor("Calibration sample", levels = type_levels)
)

# ===== 把 250 / (125,2.448) / (62.5,1.952) 这三个校准点改成“Uncontaminated sample” =====
standard_df <- standard_df |>
  mutate(
    type = case_when(
      X == 250 ~ factor("Uncontaminated sample", levels = type_levels),
      X == 125  & abs(Y - 2.448) < 1e-6 ~ factor("Uncontaminated sample", levels = type_levels),
      X == 62.5 & abs(Y - 1.952) < 1e-6 ~ factor("Uncontaminated sample", levels = type_levels),
      TRUE ~ type
    )
  )

# ===== 形状：可填充 shape（21/22/23/24）实现“空心 vs 实心” =====
shape_vals <- c(
  "Calibration sample"                       = 21, # 圆
  "Uncontaminated sample"                    = 21, # 圆
  "Contaminated sample, 5% contamination"    = 24, # 三角
  "Contaminated sample, 10% contamination"   = 22, # 方块
  "Contaminated sample, 20% contamination"   = 23  # 菱形
)

# ===== 黑白灰：边框色 + 填充色 =====
# 你要求：5% contaminated 更“灰” => 这里用 grey60（可自行微调）
color_vals <- c(
  "Calibration sample"                       = "black",
  "Uncontaminated sample"                    = "grey45",
  "Contaminated sample, 5% contamination"    = "grey60",
  "Contaminated sample, 10% contamination"   = "grey65",
  "Contaminated sample, 20% contamination"   = "grey80"
)

fill_vals <- c(
  "Calibration sample"                       = "white",   # 空心
  "Uncontaminated sample"                    = "grey45",
  "Contaminated sample, 5% contamination"    = "grey60",
  "Contaminated sample, 10% contamination"   = "grey65",
  "Contaminated sample, 20% contamination"   = "grey80"
)

# y 上限
y_max <- max(up, unknown_df$Y, y_standard, na.rm = TRUE) * 1.1

# x 轴范围
x_upper  <- standard_concentration * 1.1
x_limits <- c(
  min(c(unknown_df$X[unknown_df$X > 0], conc_points), na.rm = TRUE),
  x_upper
)
x_breaks <- sort(unique(c(breaks_conc, standard_concentration)))

# ===== 画图 =====
p_all <- ggplot() +
  geom_ribbon(aes(x = conc_grid_pos, ymin = low, ymax = up),
              fill = "grey70", alpha = 0.25) +
  geom_line(aes(x = conc_grid_pos, y = med),
            color = "black", linewidth = 1) +
  
  # 1) 真实 Calibration standards：空心圆
  geom_point(
    data = standard_df |> filter(type == "Calibration sample"),
    aes(x = X, y = Y, shape = type, color = type),
    fill = "white",
    size = 2.2,
    stroke = 0.8
  ) +
  
  # 2) 被你“提升”为 Uncontaminated 的那 3 个 standards：按 Uncontaminated 的灰色实心画
  geom_point(
    data = standard_df |> filter(type != "Calibration sample"),
    aes(x = X, y = Y, shape = type, color = type, fill = type),
    size = 2.6,
    alpha = 0.95,
    stroke = 0.35
  ) +
  
  geom_path(
    data = unknown_df |>
      filter(group == "Contaminated") |>
      arrange(sample_f, X),
    aes(x = X, y = Y, group = sample_f, color = type),
    linewidth = 0.4, alpha = 0.6,
    show.legend = FALSE
  ) +
  
  # unknown 点：实心（fill 映射）
  geom_point(
    data = unknown_df,
    aes(x = X, y = Y, shape = type, color = type, fill = type),
    size = 2.6, alpha = 0.95, stroke = 0.35
  ) +
  
  # 三个 scale 都保留同一 breaks，且 drop=FALSE 防止 legend 丢 level
  scale_shape_manual(values = shape_vals, breaks = type_levels, drop = FALSE) +
  scale_color_manual(values = color_vals, breaks = type_levels, drop = FALSE, name = "Sample type") +
  scale_fill_manual(values = fill_vals, breaks = type_levels, drop = FALSE, name = "Sample type") +
  
  # 只保留一个 legend：用 color guide 承载，同时 override fill/shape，让 Calibration icon 显示为空心
  guides(
    fill  = "none",
    shape = "none",
    color = guide_legend(
      override.aes = list(
        shape = unname(shape_vals[type_levels]),
        fill  = unname(fill_vals[type_levels]),
        size  = 3,
        alpha = 1,
        stroke = 0.8
      )
    )
  ) +
  
  scale_x_log10(
    limits = x_limits, breaks = x_breaks,
    labels = label_number(accuracy = 0.01, trim = TRUE),
    expand = c(0, 0.02)
  ) +
  scale_y_continuous(
    limits = c(0, y_max),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  labs(x = "Concentration", y = "y") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)  # x 轴竖排
  )

print(p_all)

ggsave("standard_plus_unknown_samples_lines_250_bw.png",
       p_all, width = 8, height = 5, dpi = 300)
