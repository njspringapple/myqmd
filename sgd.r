# =============================================
#  改进版牛顿法可视化 — 切线清晰显示
# =============================================

# --- 1. 定义方程和导数 ---
f  <- function(x) x^6 - 5*x^4 + 3*x^3 - 7*x^2 + 2*x + 4
df <- function(x) 6*x^5 - 20*x^3 + 9*x^2 - 14*x + 2

# --- 2. 牛顿法 ---
newton <- function(x0, f, df, max_iter = 8, tol = 1e-10) {
  x_hist <- x0
  for (k in seq_len(max_iter)) {
    fx  <- f(x_hist[k])
    dfx <- df(x_hist[k])
    if (abs(dfx) < 1e-14) break
    x_new <- x_hist[k] - fx / dfx
    x_hist <- c(x_hist, x_new)
    if (abs(x_new - x_hist[k]) < tol) break
  }
  data.frame(k = 0:(length(x_hist)-1), x = x_hist, fx = f(x_hist))
}

# =============================================
# 方案A: 单初始点详细展示（最清晰）
# =============================================
x0_demo <- 2.8
d <- newton(x0_demo, f, df)

par(mfrow = c(1,1), mar = c(5, 5, 4, 2), bg = "white")
xr <- seq(-3.5, 3.8, length.out = 1000)

plot(xr, f(xr), type = "l", lwd = 3, col = "gray30",
     xlab = "x", ylab = "f(x)", cex.lab = 1.3,
     main = bquote("牛顿法迭代详解   " * x[0] == .(x0_demo)),
     ylim = c(-35, 50), cex.main = 1.4)
abline(h = 0, lwd = 1.5, col = "gray50")
grid(col = "gray90")

# 颜色渐变: 每一步不同颜色
step_cols <- c("red", "blue", "forestgreen", "darkorange",
               "purple", "brown", "deeppink", "cyan4")

n_steps <- min(nrow(d) - 1, 6)  # 最多画6步

for (j in 1:n_steps) {
  xk  <- d$x[j]
  fk  <- d$fx[j]
  dfk <- df(xk)
  xk1 <- d$x[j + 1]  # = xk - f(xk)/f'(xk)
  cc  <- step_cols[j]

  # ---- (1) 画切线: 在 xk 附近延伸一段 ----
  # 切线方程: y = f(xk) + f'(xk) * (t - xk)
  tangent_y <- function(t) fk + dfk * (t - xk)

  # 切线范围: 从 xk 左边一点 到 xk1 右边一点
  t_left  <- min(xk, xk1) - 0.4
  t_right <- max(xk, xk1) + 0.4
  t_seq   <- seq(t_left, t_right, length.out = 200)
  # 限制切线在画面范围内
  ty <- tangent_y(t_seq)
  keep <- (ty > -40) & (ty < 55)
  t_seq <- t_seq[keep]; ty <- ty[keep]

  lines(t_seq, ty, col = cc, lwd = 2.5, lty = 2)

  # ---- (2) 从 (xk, fk) 画一个大圆点 ----
  points(xk, fk, pch = 19, cex = 2.0, col = cc)

  # ---- (3) 从 (xk, fk) 到 (xk, 0) 的垂线 ----
  arrows(xk, fk, xk, 0, col = cc, lwd = 1.5, lty = 3,
         length = 0.1, angle = 20)

  # ---- (4) 从 (xk, 0) 沿 x 轴到 (xk1, 0) 的箭头 ----
  arrows(xk, 0, xk1, 0, col = cc, lwd = 2, length = 0.12,
         angle = 25, code = 2)

  # ---- (5) 在 xk1 处画 x 轴上的小点 ----
  points(xk1, 0, pch = 4, cex = 1.5, col = cc, lwd = 2)

  # ---- (6) 标注文字 ----
  # 标注 (xk, fk) 点
  label_pos <- if (fk > 0) 3 else 1  # 正值标上方, 负值标下方
  text(xk, fk,
       labels = bquote(bold("(")*x[.(j-1)] == .(round(xk, 2))*bold(")")),
       pos = label_pos, col = cc, cex = 0.9, font = 2, offset = 0.8)

  # 在 x 轴标注切线交点
  text(xk1, 0,
       labels = bquote(x[.(j)]),
       pos = 1, col = cc, cex = 1.0, font = 2, offset = 0.6)
}

# 最终收敛点
x_final <- tail(d$x, 1)
points(x_final, 0, pch = 8, cex = 2.5, col = "red", lwd = 3)
text(x_final, -3,
     labels = bquote(bold("根  ") * x^"*" == .(round(x_final, 6))),
     col = "red", cex = 1.1, font = 2)

# 图例
legend("topleft",
       legend = paste0("第", 1:n_steps, "步  x", 0:(n_steps-1),
                       " → x", 1:n_steps),
       col = step_cols[1:n_steps], lwd = 2.5, lty = 2, pch = 19,
       bg = "white", cex = 0.85, title = "迭代步骤",
       title.col = "black")

# 添加说明框
legend("bottomleft",
       legend = c("— 虚线 = 切线",
                  "↓ 垂线 = 从曲线到x轴",
                  "→ 箭头 = x轴上移动方向",
                  "✱ = 收敛根"),
       col = c("red", "red", "red", "red"),
       lty = c(2, 3, 1, NA), pch = c(NA, NA, NA, 8),
       bg = "white", cex = 0.8, title = "图例说明")


# =============================================
# 方案B: 多初始点，但每个单独一张子图
# =============================================
starts <- c(-2.8, -0.8, 1.0, 2.8)
cols   <- c("red", "dodgerblue", "forestgreen", "darkorange")

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), bg = "white")

for (i in seq_along(starts)) {
  d <- newton(starts[i], f, df)

  # 自适应 x 范围
  x_min <- min(d$x) - 1.0
  x_max <- max(d$x) + 1.0
  xr_local <- seq(x_min, x_max, length.out = 500)
  yr_local <- f(xr_local)
  y_min <- min(c(yr_local, d$fx)) * 1.1
  y_max <- max(c(yr_local, d$fx)) * 1.1

  plot(xr_local, yr_local, type = "l", lwd = 2.5, col = "gray30",
       xlab = "x", ylab = "f(x)",
       ylim = c(max(y_min, -40), min(y_max, 60)),
       main = bquote(x[0] == .(starts[i])))
  abline(h = 0, lwd = 1.2, col = "gray50")
  grid(col = "gray90")

  n_steps <- min(nrow(d) - 1, 5)

  for (j in 1:n_steps) {
    xk  <- d$x[j]; fk <- d$fx[j]; dfk <- df(xk)
    xk1 <- d$x[j + 1]
    tangent_y <- function(t) fk + dfk * (t - xk)

    # 切线
    t_left  <- min(xk, xk1) - 0.3
    t_right <- max(xk, xk1) + 0.3
    t_seq   <- seq(t_left, t_right, length.out = 200)
    ty <- tangent_y(t_seq)
    keep <- (ty > y_min * 1.2) & (ty < y_max * 1.2)
    if (sum(keep) > 2) {
      lines(t_seq[keep], ty[keep], col = cols[i], lwd = 2.5, lty = 2)
    }

    # 点
    points(xk, fk, pch = 19, cex = 1.8, col = cols[i])

    # 垂线 (xk, fk) → (xk, 0)
    segments(xk, fk, xk, 0, col = cols[i], lwd = 1.2, lty = 3)

    # 标注
    lab_pos <- if (fk >= 0) 3 else 1
    text(xk, fk, labels = bquote(x[.(j-1)]),
         pos = lab_pos, col = cols[i], cex = 0.9, font = 2)
  }

  # 收敛根
  x_final <- tail(d$x, 1)
  points(x_final, 0, pch = 8, cex = 2.0, col = "red3", lwd = 2.5)
  text(x_final, 0,
       labels = paste0("x*=", round(x_final, 4)),
       pos = 1, col = "red3", cex = 0.8, font = 2, offset = 0.8)
}


# =============================================
# 方案C: 动画式逐步展示 (单初始点)
# =============================================
x0_demo <- 2.8
d <- newton(x0_demo, f, df)
xr <- seq(-0.5, 3.8, length.out = 500)

n_steps <- min(nrow(d) - 1, 5)

for (step in 1:n_steps) {

  par(mfrow = c(1,1), mar = c(5, 5, 4, 2), bg = "white")

  plot(xr, f(xr), type = "l", lwd = 3, col = "gray30",
       xlab = "x", ylab = "f(x)", cex.lab = 1.3,
       main = bquote("牛顿法第 " * .(step) * " 步:  " *
                        x[.(step-1)] * " → " * x[.(step)]),
       ylim = c(-35, 50), cex.main = 1.5)
  abline(h = 0, lwd = 1.5, col = "gray50")
  grid(col = "gray90")

  # 画之前所有步（灰色）
  if (step > 1) {
    for (j in 1:(step - 1)) {
      xk <- d$x[j]; fk <- d$fx[j]; dfk <- df(xk); xk1 <- d$x[j+1]
      tangent_y <- function(t) fk + dfk * (t - xk)
      t_seq <- seq(min(xk, xk1) - 0.3, max(xk, xk1) + 0.3, len = 100)
      ty <- tangent_y(t_seq)
      lines(t_seq, ty, col = "gray75", lwd = 1.5, lty = 2)
      points(xk, fk, pch = 19, cex = 1.5, col = "gray60")
      segments(xk, fk, xk, 0, col = "gray75", lty = 3)
    }
  }

  # 当前步（高亮）
  xk  <- d$x[step]; fk <- d$fx[step]; dfk <- df(xk)
  xk1 <- d$x[step + 1]
  tangent_y <- function(t) fk + dfk * (t - xk)

  t_left  <- min(xk, xk1) - 0.5
  t_right <- max(xk, xk1) + 0.5
  t_seq   <- seq(t_left, t_right, length.out = 200)
  ty <- tangent_y(t_seq)
  keep <- (ty > -40) & (ty < 55)

  # ★ 切线（粗红色虚线）
  lines(t_seq[keep], ty[keep], col = "red", lwd = 3.5, lty = 2)

  # ★ 当前点（大圆）
  points(xk, fk, pch = 19, cex = 2.5, col = "red")

  # ★ 垂线
  arrows(xk, fk, xk, 0, col = "blue", lwd = 2, lty = 1,
         length = 0.12, angle = 20)

  # ★ x 轴箭头
  arrows(xk, 0, xk1, 0, col = "forestgreen", lwd = 2.5,
         length = 0.15, angle = 25)

  # ★ 新交点
  points(xk1, 0, pch = 17, cex = 2.0, col = "forestgreen")

  # 标注
  text(xk, fk,
       labels = bquote(bold("  ") * x[.(step-1)] == .(round(xk, 3))),
       pos = 3, col = "red", cex = 1.2, font = 2, offset = 1)
  text(xk, fk,
       labels = bquote("f(" * x[.(step-1)] * ") = " * .(round(fk, 2))),
       pos = 3, col = "red", cex = 0.9, offset = 2.2)
  text(xk1, 0,
       labels = bquote(bold("  ") * x[.(step)] == .(round(xk1, 4))),
       pos = 1, col = "forestgreen", cex = 1.2, font = 2, offset = 0.8)

  # 说明框
  legend("topleft",
         legend = c(
           bquote("切线斜率 f'(" * x[.(step-1)] * ") = " * .(round(dfk, 2))),
           bquote(x[.(step)] == x[.(step-1)] - frac(f(x[.(step-1)]),
                  f*"'"*(x[.(step-1)])) == .(round(xk1, 5)))),
         bty = "n", cex = 1.0, text.col = c("red", "forestgreen"))

  Sys.sleep(0.5)  # 如果想看动画效果
}

cat("\n✅ 三种方案均已绘制完成\n")
cat("   方案A: 单图完整展示所有步骤\n")
cat("   方案B: 四个初始点分面展示\n")
cat("   方案C: 逐步动画展示\n")