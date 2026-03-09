# 完整代码：对称 vs 非对称矩阵的特征分解过程对比

par(mfrow = c(2, 4), mar = c(3, 3, 3, 1))

# ========================================================
#                    对称矩阵（上排）
# ========================================================

A_sym <- matrix(c(3, 1, 
                  1, 3), 2, 2, byrow = TRUE)

eig_sym <- eigen(A_sym)
Q <- eig_sym$vectors  # 正交矩阵
v1 <- Q[, 1]
v2 <- Q[, 2]
lambda_sym <- eig_sym$values

# 输入向量
x <- c(2, 1)

# 计算各步骤
x_eig_sym <- t(Q) %*% x                          # Step1: Q^T x
x_stretched_sym <- diag(lambda_sym) %*% x_eig_sym  # Step2: Λ (Q^T x)
Ax_sym <- Q %*% x_stretched_sym                  # Step3: Q Λ Q^T x

# ---------- 对称：图1 原始 ----------
plot(NULL, xlim = c(-3, 4), ylim = c(-3, 4), asp = 1,
     main = "对称：原始向量 x", xlab = "", ylab = "")
abline(h = 0, v = 0, col = "gray70")

# 标准基
arrows(0, 0, 2.5, 0, col = "gray50", lwd = 1, length = 0.08)
arrows(0, 0, 0, 2.5, col = "gray50", lwd = 1, length = 0.08)
text(2.7, 0, expression(e[1]), col = "gray40", cex = 0.9)
text(0, 2.8, expression(e[2]), col = "gray40", cex = 0.9)

# 向量 x
arrows(0, 0, x[1], x[2], col = "blue", lwd = 3, length = 0.1)
points(x[1], x[2], pch = 19, col = "blue", cex = 1.5)
text(x[1] + 0.4, x[2] + 0.3, paste0("x=(", x[1], ",", x[2], ")"), col = "blue", cex = 0.9)

# ---------- 对称：图2 显示特征基（正交）----------
plot(NULL, xlim = c(-3, 4), ylim = c(-3, 4), asp = 1,
     main = "Step1: 转到特征坐标系", xlab = "", ylab = "")
abline(h = 0, v = 0, col = "gray80")

# 特征向量（正交）
arrows(0, 0, v1[1]*2.5, v1[2]*2.5, col = "red", lwd = 2, length = 0.08)
arrows(0, 0, v2[1]*2.5, v2[2]*2.5, col = "red", lwd = 2, length = 0.08)
text(v1[1]*2.8, v1[2]*2.8, expression(v[1]), col = "red", cex = 0.9)
text(v2[1]*2.8, v2[2]*2.8, expression(v[2]), col = "red", cex = 0.9)

# 直角符号
rect_size <- 0.2
segments(v2[1]*rect_size, v2[2]*rect_size, 
         v2[1]*rect_size + v1[1]*rect_size, v2[2]*rect_size + v1[2]*rect_size, 
         col = "red", lwd = 2)
segments(v1[1]*rect_size, v1[2]*rect_size, 
         v1[1]*rect_size + v2[1]*rect_size, v1[2]*rect_size + v2[2]*rect_size, 
         col = "red", lwd = 2)

# 向量 x 分解到特征方向
proj1_sym <- as.numeric(x_eig_sym[1]) * v1
proj2_sym <- as.numeric(x_eig_sym[2]) * v2

arrows(0, 0, proj1_sym[1], proj1_sym[2], col = "purple", lwd = 2, length = 0.08)
arrows(0, 0, proj2_sym[1], proj2_sym[2], col = "orange", lwd = 2, length = 0.08)
segments(proj1_sym[1], proj1_sym[2], x[1], x[2], col = "orange", lty = 2)
segments(proj2_sym[1], proj2_sym[2], x[1], x[2], col = "purple", lty = 2)

arrows(0, 0, x[1], x[2], col = "blue", lwd = 3, length = 0.1)
points(x[1], x[2], pch = 19, col = "blue", cex = 1.5)

text(0, -2.5, paste0("特征坐标: (", round(x_eig_sym[1], 2), ", ", round(x_eig_sym[2], 2), ")"), 
     col = "black", cex = 0.85)
text(0, -3, "v₁·v₂ = 0 (正交)", col = "darkgreen", cex = 0.85, font = 2)

# ---------- 对称：图3 拉伸 ----------
plot(NULL, xlim = c(-4, 10), ylim = c(-4, 10), asp = 1,
     main = paste0("Step2: 拉伸 (×", lambda_sym[1], ", ×", lambda_sym[2], ")"), 
     xlab = "", ylab = "")
abline(h = 0, v = 0, col = "gray80")

# 特征方向
arrows(0, 0, v1[1]*8, v1[2]*8, col = "red", lwd = 1, lty = 2, length = 0.08)
arrows(0, 0, v2[1]*5, v2[2]*5, col = "red", lwd = 1, lty = 2, length = 0.08)

# 拉伸前（虚线）
arrows(0, 0, proj1_sym[1], proj1_sym[2], col = "purple", lwd = 1, lty = 3, length = 0.08)
arrows(0, 0, proj2_sym[1], proj2_sym[2], col = "orange", lwd = 1, lty = 3, length = 0.08)

# 拉伸后
proj1_str_sym <- as.numeric(x_stretched_sym[1]) * v1
proj2_str_sym <- as.numeric(x_stretched_sym[2]) * v2

arrows(0, 0, proj1_str_sym[1], proj1_str_sym[2], col = "purple", lwd = 3, length = 0.1)
arrows(0, 0, proj2_str_sym[1], proj2_str_sym[2], col = "orange", lwd = 3, length = 0.1)

# 平行四边形
segments(proj1_str_sym[1], proj1_str_sym[2], Ax_sym[1], Ax_sym[2], col = "orange", lty = 2)
segments(proj2_str_sym[1], proj2_str_sym[2], Ax_sym[1], Ax_sym[2], col = "purple", lty = 2)

arrows(0, 0, Ax_sym[1], Ax_sym[2], col = "darkgreen", lwd = 3, length = 0.1)
points(Ax_sym[1], Ax_sym[2], pch = 19, col = "darkgreen", cex = 1.5)

text(proj1_str_sym[1]/2 - 0.5, proj1_str_sym[2]/2 + 0.5, 
     paste0(round(x_stretched_sym[1], 1)), col = "purple", cex = 0.9, font = 2)
text(proj2_str_sym[1]/2 + 0.5, proj2_str_sym[2]/2 - 0.3, 
     paste0(round(x_stretched_sym[2], 1)), col = "orange", cex = 0.9, font = 2)

# ---------- 对称：图4 结果 ----------
plot(NULL, xlim = c(-2, 10), ylim = c(-2, 8), asp = 1,
     main = "Step3: 结果 Ax", xlab = "", ylab = "")
abline(h = 0, v = 0, col = "gray70")

arrows(0, 0, 8, 0, col = "gray50", lwd = 1, length = 0.08)
arrows(0, 0, 0, 6, col = "gray50", lwd = 1, length = 0.08)

# 原向量
arrows(0, 0, x[1], x[2], col = "blue", lwd = 2, lty = 2, length = 0.08)
text(x[1] + 0.3, x[2] + 0.4, "x", col = "blue", cex = 1)

# 结果
arrows(0, 0, Ax_sym[1], Ax_sym[2], col = "darkgreen", lwd = 3, length = 0.1)
points(Ax_sym[1], Ax_sym[2], pch = 19, col = "darkgreen", cex = 1.5)
text(Ax_sym[1] + 0.3, Ax_sym[2] + 0.4, 
     paste0("Ax=(", round(Ax_sym[1], 1), ",", round(Ax_sym[2], 1), ")"), 
     col = "darkgreen", cex = 0.9)

text(4, -1.5, expression(Ax == Q * Lambda * Q^T * x), cex = 1)

# ========================================================
#                   非对称矩阵（下排）
# ========================================================

A_nonsym <- matrix(c(2, 1, 
                     0, 3), 2, 2, byrow = TRUE)

eig_nonsym <- eigen(A_nonsym)
P <- eig_nonsym$vectors  # 不正交！
u1 <- P[, 1]
u2 <- P[, 2]
lambda_nonsym <- eig_nonsym$values

# 计算各步骤
x_eig_nonsym <- solve(P) %*% x                              # Step1: P^{-1} x（不是转置！）
x_stretched_nonsym <- diag(lambda_nonsym) %*% x_eig_nonsym  # Step2: Λ (P^{-1} x)
Ax_nonsym <- P %*% x_stretched_nonsym                       # Step3: P Λ P^{-1} x

# ---------- 非对称：图1 原始 ----------
plot(NULL, xlim = c(-3, 4), ylim = c(-3, 4), asp = 1,
     main = "非对称：原始向量 x", xlab = "", ylab = "")
abline(h = 0, v = 0, col = "gray70")

arrows(0, 0, 2.5, 0, col = "gray50", lwd = 1, length = 0.08)
arrows(0, 0, 0, 2.5, col = "gray50", lwd = 1, length = 0.08)
text(2.7, 0, expression(e[1]), col = "gray40", cex = 0.9)
text(0, 2.8, expression(e[2]), col = "gray40", cex = 0.9)

arrows(0, 0, x[1], x[2], col = "blue", lwd = 3, length = 0.1)
points(x[1], x[2], pch = 19, col = "blue", cex = 1.5)
text(x[1] + 0.4, x[2] + 0.3, paste0("x=(", x[1], ",", x[2], ")"), col = "blue", cex = 0.9)

# ---------- 非对称：图2 显示特征基（不正交）----------
plot(NULL, xlim = c(-3, 4), ylim = c(-3, 4), asp = 1,
     main = "Step1: 转到特征基（斜基）", xlab = "", ylab = "")
abline(h = 0, v = 0, col = "gray80")

# 特征向量（不正交）
arrows(0, 0, u1[1]*2.5, u1[2]*2.5, col = "red", lwd = 2, length = 0.08)
arrows(0, 0, u2[1]*2.5, u2[2]*2.5, col = "red", lwd = 2, length = 0.08)
text(u1[1]*2.8, u1[2]*2.8, expression(u[1]), col = "red", cex = 0.9)
text(u2[1]*2.8 + 0.3, u2[2]*2.8, expression(u[2]), col = "red", cex = 0.9)

# 画夹角弧线（不是90°）
theta1 <- atan2(u1[2], u1[1])
theta2 <- atan2(u2[2], u2[1])
if (theta1 > theta2) { tmp <- theta1; theta1 <- theta2; theta2 <- tmp }
arc_theta <- seq(theta1, theta2, length.out = 30)
lines(0.4 * cos(arc_theta), 0.4 * sin(arc_theta), col = "orange", lwd = 2)

# 向量 x 分解到特征方向
proj1_nonsym <- as.numeric(x_eig_nonsym[1]) * u1
proj2_nonsym <- as.numeric(x_eig_nonsym[2]) * u2

arrows(0, 0, proj1_nonsym[1], proj1_nonsym[2], col = "purple", lwd = 2, length = 0.08)
arrows(0, 0, proj2_nonsym[1], proj2_nonsym[2], col = "orange", lwd = 2, length = 0.08)
segments(proj1_nonsym[1], proj1_nonsym[2], x[1], x[2], col = "orange", lty = 2)
segments(proj2_nonsym[1], proj2_nonsym[2], x[1], x[2], col = "purple", lty = 2)

arrows(0, 0, x[1], x[2], col = "blue", lwd = 3, length = 0.1)
points(x[1], x[2], pch = 19, col = "blue", cex = 1.5)

angle_deg <- acos(sum(u1 * u2) / (sqrt(sum(u1^2)) * sqrt(sum(u2^2)))) * 180 / pi
text(0, -2.5, paste0("特征坐标: (", round(x_eig_nonsym[1], 2), ", ", round(x_eig_nonsym[2], 2), ")"), 
     col = "black", cex = 0.85)
text(0, -3, paste0("夹角=", round(angle_deg, 1), "° (不正交)"), col = "red", cex = 0.85, font = 2)

# ---------- 非对称：图3 拉伸 ----------
plot(NULL, xlim = c(-3, 8), ylim = c(-3, 8), asp = 1,
     main = paste0("Step2: 拉伸 (×", lambda_nonsym[1], ", ×", lambda_nonsym[2], ")"), 
     xlab = "", ylab = "")
abline(h = 0, v = 0, col = "gray80")

# 特征方向
arrows(0, 0, u1[1]*6, u1[2]*6, col = "red", lwd = 1, lty = 2, length = 0.08)
arrows(0, 0, u2[1]*6, u2[2]*6, col = "red", lwd = 1, lty = 2, length = 0.08)

# 拉伸前（虚线）
arrows(0, 0, proj1_nonsym[1], proj1_nonsym[2], col = "purple", lwd = 1, lty = 3, length = 0.08)
arrows(0, 0, proj2_nonsym[1], proj2_nonsym[2], col = "orange", lwd = 1, lty = 3, length = 0.08)

# 拉伸后
proj1_str_nonsym <- as.numeric(x_stretched_nonsym[1]) * u1
proj2_str_nonsym <- as.numeric(x_stretched_nonsym[2]) * u2

arrows(0, 0, proj1_str_nonsym[1], proj1_str_nonsym[2], col = "purple", lwd = 3, length = 0.1)
arrows(0, 0, proj2_str_nonsym[1], proj2_str_nonsym[2], col = "orange", lwd = 3, length = 0.1)

# 平行四边形
segments(proj1_str_nonsym[1], proj1_str_nonsym[2], Ax_nonsym[1], Ax_nonsym[2], col = "orange", lty = 2)
segments(proj2_str_nonsym[1], proj2_str_nonsym[2], Ax_nonsym[1], Ax_nonsym[2], col = "purple", lty = 2)

arrows(0, 0, Ax_nonsym[1], Ax_nonsym[2], col = "darkgreen", lwd = 3, length = 0.1)
points(Ax_nonsym[1], Ax_nonsym[2], pch = 19, col = "darkgreen", cex = 1.5)

text(proj1_str_nonsym[1]/2 - 0.3, proj1_str_nonsym[2]/2 + 0.4, 
     paste0(round(x_stretched_nonsym[1], 1)), col = "purple", cex = 0.9, font = 2)
text(proj2_str_nonsym[1]/2 + 0.5, proj2_str_nonsym[2]/2 - 0.3, 
     paste0(round(x_stretched_nonsym[2], 1)), col = "orange", cex = 0.9, font = 2)

# ---------- 非对称：图4 结果 ----------
plot(NULL, xlim = c(-2, 8), ylim = c(-2, 6), asp = 1,
     main = "Step3: 结果 Ax", xlab = "", ylab = "")
abline(h = 0, v = 0, col = "gray70")

arrows(0, 0, 6, 0, col = "gray50", lwd = 1, length = 0.08)
arrows(0, 0, 0, 5, col = "gray50", lwd = 1, length = 0.08)

# 原向量
arrows(0, 0, x[1], x[2], col = "blue", lwd = 2, lty = 2, length = 0.08)
text(x[1] + 0.3, x[2] + 0.4, "x", col = "blue", cex = 1)

# 结果
arrows(0, 0, Ax_nonsym[1], Ax_nonsym[2], col = "darkgreen", lwd = 3, length = 0.1)
points(Ax_nonsym[1], Ax_nonsym[2], pch = 19, col = "darkgreen", cex = 1.5)
text(Ax_nonsym[1] + 0.3, Ax_nonsym[2] + 0.4, 
     paste0("Ax=(", round(Ax_nonsym[1], 1), ",", round(Ax_nonsym[2], 1), ")"), 
     col = "darkgreen", cex = 0.9)

text(3, -1.5, expression(Ax == P * Lambda * P^{-1} * x), cex = 1)

par(mfrow = c(1, 1))

# ==================== 打印对比 ====================
cat("\n================== 对比总结 ==================\n\n")

cat("【对称矩阵】 A = A^T\n")
print(A_sym)
cat("特征值:", lambda_sym, "\n")
cat("特征向量点积: v1·v2 =", round(sum(v1 * v2), 6), "→ 正交 ✓\n")
cat("分解: A = Q Λ Q^T，其中 Q^T = Q^{-1}\n")
cat("Step1 用 Q^T（转置 = 投影）\n\n")

cat("【非对称矩阵】 A ≠ A^T\n")
print(A_nonsym)
cat("特征值:", lambda_nonsym, "\n")
cat("特征向量点积: u1·u2 =", round(sum(u1 * u2), 4), "→ 不正交 ✗\n")
cat("分解: A = P Λ P^{-1}，其中 P^T ≠ P^{-1}\n")
cat("Step1 用 P^{-1}（求逆，不是投影）\n")






par(mfrow = c(1, 2), mar = c(4, 4, 4, 1))

x <- c(2, 1)

# ========== 左图：正交基，投影 ==========
A_sym <- matrix(c(3, 1, 1, 3), 2, 2, byrow = TRUE)
eig_sym <- eigen(A_sym)
Q <- eig_sym$vectors
v1 <- Q[, 1]
v2 <- Q[, 2]

plot(NULL, xlim = c(-1, 3), ylim = c(-1, 3), asp = 1,
     main = "正交基：投影得到坐标", xlab = "", ylab = "")
abline(h = 0, v = 0, col = "gray80")

# 特征向量
arrows(0, 0, v1[1]*2.5, v1[2]*2.5, col = "red", lwd = 2, length = 0.08)
arrows(0, 0, v2[1]*2.5, v2[2]*2.5, col = "red", lwd = 2, length = 0.08)
text(v1[1]*2.7, v1[2]*2.7, expression(v[1]), col = "red")
text(v2[1]*2.7, v2[2]*2.7, expression(v[2]), col = "red")

# 直角符号
rect_size <- 0.15
segments(v2[1]*rect_size, v2[2]*rect_size, 
         v2[1]*rect_size + v1[1]*rect_size, v2[2]*rect_size + v1[2]*rect_size, 
         col = "red", lwd = 2)
segments(v1[1]*rect_size, v1[2]*rect_size, 
         v1[1]*rect_size + v2[1]*rect_size, v1[2]*rect_size + v2[2]*rect_size, 
         col = "red", lwd = 2)

# 向量 x
arrows(0, 0, x[1], x[2], col = "blue", lwd = 3, length = 0.1)
points(x[1], x[2], pch = 19, col = "blue", cex = 1.5)
text(x[1] + 0.2, x[2] + 0.2, "x", col = "blue", cex = 1.2)

# 投影（垂直落下）
proj1 <- sum(x * v1) * v1  # x 在 v1 上的投影
proj2 <- sum(x * v2) * v2  # x 在 v2 上的投影

# 画垂直投影线（虚线）
segments(x[1], x[2], proj1[1], proj1[2], col = "purple", lty = 2, lwd = 2)
segments(x[1], x[2], proj2[1], proj2[2], col = "orange", lty = 2, lwd = 2)

# 投影点
points(proj1[1], proj1[2], pch = 19, col = "purple", cex = 1.2)
points(proj2[1], proj2[2], pch = 19, col = "orange", cex = 1.2)

# 投影分量
arrows(0, 0, proj1[1], proj1[2], col = "purple", lwd = 2, length = 0.08)
arrows(0, 0, proj2[1], proj2[2], col = "orange", lwd = 2, length = 0.08)

# 标注坐标值
c1 <- sum(x * v1)  # = (Q^T x)[1]
c2 <- sum(x * v2)  # = (Q^T x)[2]
text(proj1[1]/2 + 0.25, proj1[2]/2 - 0.2, round(c1, 2), col = "purple", font = 2)
text(proj2[1]/2 - 0.25, proj2[2]/2 + 0.2, round(c2, 2), col = "orange", font = 2)

text(1, -0.7, expression(c[1] == v[1]^T * x), col = "purple", cex = 0.9)
text(1, -0.95, expression(c[2] == v[2]^T * x), col = "orange", cex = 0.9)
text(1, -0.45, "投影（垂直落下）", col = "black", cex = 0.9)

# ========== 右图：非正交基，平行线 ==========
A_nonsym <- matrix(c(2, 1, 0, 3), 2, 2, byrow = TRUE)
eig_nonsym <- eigen(A_nonsym)
P <- eig_nonsym$vectors
u1 <- P[, 1]
u2 <- P[, 2]

plot(NULL, xlim = c(-1, 3), ylim = c(-1, 3), asp = 1,
     main = "非正交基：平行分解得坐标", xlab = "", ylab = "")
abline(h = 0, v = 0, col = "gray80")

# 特征向量
arrows(0, 0, u1[1]*2.5, u1[2]*2.5, col = "red", lwd = 2, length = 0.08)
arrows(0, 0, u2[1]*2.5, u2[2]*2.5, col = "red", lwd = 2, length = 0.08)
text(u1[1]*2.7, u1[2]*2.7, expression(u[1]), col = "red")
text(u2[1]*2.7 + 0.15, u2[2]*2.7, expression(u[2]), col = "red")

# 夹角弧线
theta1 <- atan2(u1[2], u1[1])
theta2 <- atan2(u2[2], u2[1])
if (theta1 > theta2) { tmp <- theta1; theta1 <- theta2; theta2 <- tmp }
arc_theta <- seq(theta1, theta2, length.out = 30)
lines(0.3 * cos(arc_theta), 0.3 * sin(arc_theta), col = "orange", lwd = 2)

# 向量 x
arrows(0, 0, x[1], x[2], col = "blue", lwd = 3, length = 0.1)
points(x[1], x[2], pch = 19, col = "blue", cex = 1.5)
text(x[1] + 0.2, x[2] + 0.2, "x", col = "blue", cex = 1.2)

# 非正交分解：x = c1*u1 + c2*u2，用 P^{-1} 求 c
c_nonsym <- solve(P) %*% x
c1_nonsym <- c_nonsym[1]
c2_nonsym <- c_nonsym[2]

comp1 <- as.numeric(c1_nonsym) * u1
comp2 <- as.numeric(c2_nonsym) * u2

# 画平行线（不是垂直的！）
# 从 x 画平行于 u2 的线到 u1 轴
segments(x[1], x[2], comp1[1], comp1[2], col = "purple", lty = 2, lwd = 2)
# 从 x 画平行于 u1 的线到 u2 轴  
segments(x[1], x[2], comp2[1], comp2[2], col = "orange", lty = 2, lwd = 2)

# 分量
arrows(0, 0, comp1[1], comp1[2], col = "purple", lwd = 2, length = 0.08)
arrows(0, 0, comp2[1], comp2[2], col = "orange", lwd = 2, length = 0.08)

points(comp1[1], comp1[2], pch = 19, col = "purple", cex = 1.2)
points(comp2[1], comp2[2], pch = 19, col = "orange", cex = 1.2)

# 标注坐标值
text(comp1[1]/2 + 0.2, comp1[2]/2 - 0.15, round(c1_nonsym, 2), col = "purple", font = 2)
text(comp2[1]/2 - 0.15, comp2[2]/2 + 0.2, round(c2_nonsym, 2), col = "orange", font = 2)

text(1, -0.7, expression(c == P^{-1} * x), col = "black", cex = 0.9)
text(1, -0.45, "平行分解（斜着过去）", col = "black", cex = 0.9)

par(mfrow = c(1, 1))

# ========== 数值验证 ==========
cat("\n========== 正交基（对称矩阵）==========\n")
cat("方法：Q^T x（投影）\n")
cat("c1 = v1·x =", round(sum(v1 * x), 4), "\n")
cat("c2 = v2·x =", round(sum(v2 * x), 4), "\n")
cat("验证：c1*v1 + c2*v2 =", round(sum(v1*x)*v1 + sum(v2*x)*v2, 4), "\n")
cat("原始 x =", x, "\n")

cat("\n========== 非正交基（非对称矩阵）==========\n")
cat("方法：P^{-1} x（解方程）\n")
cat("c1 =", round(c1_nonsym, 4), "\n")
cat("c2 =", round(c2_nonsym, 4), "\n")
cat("验证：c1*u1 + c2*u2 =", round(as.numeric(c1_nonsym)*u1 + as.numeric(c2_nonsym)*u2, 4), "\n")
cat("原始 x =", x, "\n")

cat("\n========== 关键区别 ==========\n")
cat("正交：虚线垂直于坐标轴（投影）\n")
cat("非正交：虚线平行于另一个轴（平行四边形法则）\n")



















# ========== 情况1：对称正定 → 奇异值 = 特征值 ==========
A_spd <- matrix(c(3, 1, 1, 3), 2, 2)

eig_vals <- eigen(A_spd)$values
sing_vals <- svd(A_spd)$d

cat("对称正定矩阵：\n")
print(A_spd)
cat("特征值：", sort(eig_vals, decreasing = TRUE), "\n")
cat("奇异值：", sing_vals, "\n")
cat("相等？", all.equal(sort(eig_vals, decreasing = TRUE), sing_vals), "\n\n")

# ========== 情况2：对称但有负特征值 → 奇异值 = |特征值| ==========
A_sym_neg <- matrix(c(2, 1, 1, -1), 2, 2)

eig_vals2 <- eigen(A_sym_neg)$values
sing_vals2 <- svd(A_sym_neg)$d

cat("对称但不正定：\n")
print(A_sym_neg)
cat("特征值：", eig_vals2, "\n")
cat("奇异值：", sing_vals2, "\n")
cat("奇异值 = |特征值|：", all.equal(sort(abs(eig_vals2), decreasing = TRUE), sing_vals2), "\n\n")

# ========== 情况3：非对称 → 完全不同 ==========
A_nonsym <- matrix(c(2, 1, 0, 3), 2, 2)

eig_vals3 <- eigen(A_nonsym)$values
sing_vals3 <- svd(A_nonsym)$d

cat("非对称矩阵：\n")
print(A_nonsym)
cat("特征值：", eig_vals3, "\n")
cat("奇异值：", sing_vals3, "\n")
cat("差别很大！\n\n")

# ========== 总结表格 ==========
cat("========== 总结 ==========\n")
cat("对称正定：   σ = λ\n")
cat("对称不定：   σ = |λ|\n")
cat("非对称：     σ 和 λ 没有简单关系\n")