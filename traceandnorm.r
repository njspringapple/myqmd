# ============================================================
# 矩阵范数与迹 - R语言交互式练习
# 目的：通过代码直观理解矩阵范数的概念和应用
# ============================================================

# 清空环境
rm(list = ls())

# ============================================================
# 第一部分：矩阵的迹（Trace）
# ============================================================

cat("=" |> rep(60) |> paste(collapse = ""), "\n")
cat("第一部分：矩阵的迹\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n\n")

# 创建一个简单的矩阵
A <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)
cat("矩阵 A:\n")
print(A)

# 手动计算迹：对角元素之和
trace_manual <- A[1,1] + A[2,2]
cat("\n手动计算迹 (a11 + a22):", trace_manual, "\n")

# 用R内置函数计算
trace_builtin <- sum(diag(A))
cat("用 sum(diag(A)) 计算:", trace_builtin, "\n")

# 练习1：验证迹等于特征值之和
cat("\n--- 练习1：迹 = 特征值之和 ---\n")
eigenvalues <- eigen(A)$values
cat("特征值:", eigenvalues, "\n")
cat("特征值之和:", sum(eigenvalues), "\n")
cat("迹:", trace_builtin, "\n")
cat("两者相等？", all.equal(sum(eigenvalues), trace_builtin), "\n")

# 练习2：验证循环性质 tr(AB) = tr(BA)
cat("\n--- 练习2：循环性质 tr(AB) = tr(BA) ---\n")
B <- matrix(c(5, 6, 7, 8), nrow = 2, byrow = TRUE)
cat("矩阵 B:\n")
print(B)
cat("\ntr(AB) =", sum(diag(A %*% B)), "\n")
cat("tr(BA) =", sum(diag(B %*% A)), "\n")
cat("相等？", sum(diag(A %*% B)) == sum(diag(B %*% A)), "\n")

# 注意：tr(AB) ≠ tr(A) * tr(B)
cat("\n注意：tr(AB) ≠ tr(A) × tr(B)\n")
cat("tr(AB) =", sum(diag(A %*% B)), "\n")
cat("tr(A) × tr(B) =", sum(diag(A)) * sum(diag(B)), "\n")

# ============================================================
# 第二部分：各种矩阵范数
# ============================================================

cat("\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n")
cat("第二部分：各种矩阵范数\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n\n")

# 创建测试矩阵
M <- matrix(c(2, 1, 0, 3), nrow = 2, byrow = TRUE)
cat("测试矩阵 M:\n")
print(M)

# --- Frobenius范数 ---
cat("\n--- Frobenius范数 ---\n")
cat("定义：所有元素平方和的平方根\n")
frob_manual <- sqrt(sum(M^2))
cat("手动计算: sqrt(", paste(M^2, collapse = " + "), ") =", frob_manual, "\n")
cat("用 norm() 函数:", norm(M, "F"), "\n")

# 验证：Frobenius范数² = tr(A'A)
cat("\n验证：||A||_F² = tr(A'A)\n")
cat("||M||_F² =", frob_manual^2, "\n")
cat("tr(M'M) =", sum(diag(t(M) %*% M)), "\n")

# --- 列和范数 (1-范数) ---
cat("\n--- 列和范数 (p=1) ---\n")
cat("定义：每列绝对值之和的最大值\n")
col_sums <- colSums(abs(M))
cat("各列绝对值之和:", col_sums, "\n")
cat("最大值 (列和范数):", max(col_sums), "\n")
cat("用 norm() 函数:", norm(M, "1"), "\n")

# --- 行和范数 (∞-范数) ---
cat("\n--- 行和范数 (p=∞) ---\n")
cat("定义：每行绝对值之和的最大值\n")
row_sums <- rowSums(abs(M))
cat("各行绝对值之和:", row_sums, "\n")
cat("最大值 (行和范数):", max(row_sums), "\n")
cat("用 norm() 函数:", norm(M, "I"), "\n")

# --- 谱范数 (2-范数) ---
cat("\n--- 谱范数 (p=2) ---\n")
cat("定义：最大奇异值 = sqrt(A'A的最大特征值)\n")
AtA <- t(M) %*% M
eigenvalues_AtA <- eigen(AtA)$values
cat("M'M 的特征值:", eigenvalues_AtA, "\n")
cat("最大特征值的平方根:", sqrt(max(eigenvalues_AtA)), "\n")
cat("用 norm() 函数:", norm(M, "2"), "\n")

# 也可以用SVD
svd_result <- svd(M)
cat("用 SVD 计算 (最大奇异值):", max(svd_result$d), "\n")

# ============================================================
# 第三部分：直观理解 - 矩阵对向量的拉伸
# ============================================================

cat("\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n")
cat("第三部分：矩阵拉伸向量的可视化\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n\n")

# 创建单位圆上的点
theta <- seq(0, 2*pi, length.out = 100)
unit_circle <- rbind(cos(theta), sin(theta))

# 用矩阵变换这些点
transformed <- M %*% unit_circle

# 计算每个变换后向量的长度
lengths_after <- sqrt(colSums(transformed^2))
cat("变换后向量长度的最大值:", max(lengths_after), "\n")
cat("谱范数:", norm(M, "2"), "\n")
cat("两者应该相等！\n")

# 绘图
par(mfrow = c(1, 2))

# 原始单位圆
plot(unit_circle[1,], unit_circle[2,], type = "l", 
     asp = 1, xlim = c(-4, 4), ylim = c(-4, 4),
     main = "原始：单位圆", xlab = "x", ylab = "y",
     col = "blue", lwd = 2)
abline(h = 0, v = 0, col = "gray", lty = 2)

# 变换后的椭圆
plot(transformed[1,], transformed[2,], type = "l",
     asp = 1, xlim = c(-4, 4), ylim = c(-4, 4),
     main = "变换后：椭圆", xlab = "x", ylab = "y",
     col = "red", lwd = 2)
abline(h = 0, v = 0, col = "gray", lty = 2)

# 标注最长的向量
max_idx <- which.max(lengths_after)
arrows(0, 0, transformed[1, max_idx], transformed[2, max_idx],
       col = "darkred", lwd = 3, length = 0.1)
text(transformed[1, max_idx]/2, transformed[2, max_idx]/2 + 0.3,
     paste("长度 =", round(max(lengths_after), 2)), col = "darkred")

par(mfrow = c(1, 1))

cat("\n图形说明：\n")
cat("- 左图：单位圆上所有长度为1的向量\n")
cat("- 右图：矩阵M把单位圆变成了椭圆\n")
cat("- 红色箭头：被拉伸最长的向量，其长度就是谱范数\n")

# ============================================================
# 第四部分：条件数
# ============================================================

cat("\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n")
cat("第四部分：条件数与病态性\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n\n")

# 良态矩阵
good_matrix <- matrix(c(2, 0, 0, 1.8), nrow = 2)
cat("良态矩阵:\n")
print(good_matrix)

svd_good <- svd(good_matrix)
cat("奇异值:", svd_good$d, "\n")
cat("条件数 (σ_max / σ_min):", max(svd_good$d) / min(svd_good$d), "\n")
cat("用 kappa() 函数:", kappa(good_matrix), "\n\n")

# 病态矩阵
bad_matrix <- matrix(c(2, 0, 0, 0.001), nrow = 2)
cat("病态矩阵:\n")
print(bad_matrix)

svd_bad <- svd(bad_matrix)
cat("奇异值:", svd_bad$d, "\n")
cat("条件数:", max(svd_bad$d) / min(svd_bad$d), "\n")
cat("用 kappa() 函数:", kappa(bad_matrix), "\n\n")

cat("解释：\n")
cat("- 条件数接近1：问题稳定，小误差不会被放大\n")
cat("- 条件数很大：问题病态，小误差会被放大很多倍\n")

# 可视化条件数的影响
par(mfrow = c(1, 2))

# 良态情况
unit_circle <- rbind(cos(theta), sin(theta))
good_transformed <- good_matrix %*% unit_circle
plot(good_transformed[1,], good_transformed[2,], type = "l",
     asp = 1, xlim = c(-3, 3), ylim = c(-3, 3),
     main = paste("良态 (κ ≈", round(kappa(good_matrix), 1), ")"),
     xlab = "x", ylab = "y", col = "green", lwd = 2)
abline(h = 0, v = 0, col = "gray", lty = 2)

# 病态情况
bad_transformed <- bad_matrix %*% unit_circle
plot(bad_transformed[1,], bad_transformed[2,], type = "l",
     asp = 1, xlim = c(-3, 3), ylim = c(-0.01, 0.01),
     main = paste("病态 (κ =", round(kappa(bad_matrix), 0), ")"),
     xlab = "x", ylab = "y", col = "red", lwd = 2)
abline(h = 0, v = 0, col = "gray", lty = 2)

par(mfrow = c(1, 1))

cat("\n图形说明：\n")
cat("- 良态：椭圆比较圆，各方向拉伸差不多\n")
cat("- 病态：椭圆极扁（几乎是一条线），某些方向被严重压缩\n")

# ============================================================
# 第五部分：Schatten范数家族
# ============================================================

cat("\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n")
cat("第五部分：Schatten范数 - 三种范数的统一\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n\n")

# 计算Schatten范数的函数
schatten_norm <- function(A, p) {
  singular_values <- svd(A)$d
  if (p == Inf) {
    return(max(singular_values))
  } else {
    return(sum(singular_values^p)^(1/p))
  }
}

cat("测试矩阵 M:\n")
print(M)

svd_M <- svd(M)
cat("\n奇异值:", svd_M$d, "\n\n")

cat("Schatten范数家族：\n")
cat("-" |> rep(40) |> paste(collapse = ""), "\n")

# p = ∞：谱范数
cat("p = ∞ (谱范数):\n")
cat("  定义：最大奇异值\n")
cat("  计算：max(", paste(round(svd_M$d, 3), collapse = ", "), ") =", 
    schatten_norm(M, Inf), "\n")
cat("  验证：norm(M, '2') =", norm(M, "2"), "\n\n")

# p = 2：Frobenius范数
cat("p = 2 (Frobenius范数):\n")
cat("  定义：奇异值平方和的平方根\n")
cat("  计算：sqrt(", paste(round(svd_M$d^2, 3), collapse = " + "), ") =",
    schatten_norm(M, 2), "\n")
cat("  验证：norm(M, 'F') =", norm(M, "F"), "\n\n")

# p = 1：核范数
cat("p = 1 (核范数):\n")
cat("  定义：奇异值之和\n")
cat("  计算：", paste(round(svd_M$d, 3), collapse = " + "), "=",
    schatten_norm(M, 1), "\n")
cat("  用途：逼迫矩阵低秩（矩阵补全问题）\n")

# ============================================================
# 第六部分：次可乘性验证
# ============================================================

cat("\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n")
cat("第六部分：次可乘性 ||AB|| ≤ ||A|| × ||B||\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n\n")

A <- matrix(c(1, 2, 3, 4), nrow = 2)
B <- matrix(c(5, 6, 7, 8), nrow = 2)

cat("矩阵 A:\n")
print(A)
cat("\n矩阵 B:\n")
print(B)

cat("\n矩阵 AB:\n")
print(A %*% B)

cat("\n验证次可乘性（谱范数）：\n")
norm_A <- norm(A, "2")
norm_B <- norm(B, "2")
norm_AB <- norm(A %*% B, "2")

cat("||A||_2 =", round(norm_A, 3), "\n")
cat("||B||_2 =", round(norm_B, 3), "\n")
cat("||AB||_2 =", round(norm_AB, 3), "\n")
cat("||A||_2 × ||B||_2 =", round(norm_A * norm_B, 3), "\n")
cat("||AB||_2 ≤ ||A||_2 × ||B||_2 ?", norm_AB <= norm_A * norm_B + 1e-10, "\n")

cat("\n意义：如果 ||A|| < 1，则 ||A^n|| ≤ ||A||^n → 0\n")
cat("这是分析迭代算法收敛性的基础！\n")

# ============================================================
# 第七部分：正交矩阵的特殊性质
# ============================================================

cat("\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n")
cat("第七部分：正交矩阵 - 不改变长度的变换\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n\n")

# 创建一个旋转矩阵（45度）
angle <- pi/4
Q <- matrix(c(cos(angle), -sin(angle), 
              sin(angle), cos(angle)), nrow = 2, byrow = TRUE)

cat("45度旋转矩阵 Q:\n")
print(round(Q, 3))

cat("\n验证正交性 Q'Q = I:\n")
print(round(t(Q) %*% Q, 10))

cat("\n正交矩阵的性质：\n")
cat("1. 条件数 =", kappa(Q), "(完美！)\n")
cat("2. 谱范数 =", norm(Q, "2"), "(不拉伸)\n")
cat("3. 所有奇异值:", svd(Q)$d, "(都是1)\n")

# 验证保长度
x <- c(3, 4)
Qx <- Q %*% x
cat("\n验证保长度：\n")
cat("原向量 x = (", paste(x, collapse = ", "), "), 长度 =", sqrt(sum(x^2)), "\n")
cat("Qx = (", paste(round(Qx, 3), collapse = ", "), "), 长度 =", sqrt(sum(Qx^2)), "\n")

# 可视化
par(mfrow = c(1, 1))
unit_circle <- rbind(cos(theta), sin(theta))
rotated <- Q %*% unit_circle

plot(unit_circle[1,], unit_circle[2,], type = "l",
     asp = 1, xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5),
     main = "正交矩阵：旋转但不变形",
     xlab = "x", ylab = "y", col = "blue", lwd = 2)
lines(rotated[1,], rotated[2,], col = "red", lwd = 2, lty = 2)
legend("topright", c("原始圆", "旋转后"), col = c("blue", "red"), lty = c(1, 2))
abline(h = 0, v = 0, col = "gray", lty = 2)

cat("\n图形说明：正交矩阵把圆旋转，但圆还是圆，不会变成椭圆！\n")

# ============================================================
# 第八部分：综合练习
# ============================================================

cat("\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n")
cat("第八部分：综合练习题\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n\n")

cat("请尝试回答以下问题，然后运行代码验证：\n\n")

# 练习矩阵
P <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE)
cat("给定矩阵 P (2×3):\n")
print(P)

cat("\n问题1：P的Frobenius范数是多少？\n")
cat("提示：所有元素平方和，再开根号\n")
cat("答案：", norm(P, "F"), "\n")

cat("\n问题2：P的列和范数是多少？\n")
cat("提示：每列绝对值之和的最大值\n")
cat("各列之和：", colSums(abs(P)), "\n")
cat("答案：", norm(P, "1"), "\n")

cat("\n问题3：P的行和范数是多少？\n")
cat("提示：每行绝对值之和的最大值\n")
cat("各行之和：", rowSums(abs(P)), "\n")
cat("答案：", norm(P, "I"), "\n")

cat("\n问题4：P'P的迹等于什么？\n")
cat("提示：应该等于Frobenius范数的平方\n")
cat("tr(P'P) =", sum(diag(t(P) %*% P)), "\n")
cat("||P||_F² =", norm(P, "F")^2, "\n")

# ============================================================
# 总结
# ============================================================

cat("\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n")
cat("总结：关键概念速查\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n\n")

cat("1. 迹 tr(A) = 对角元素之和 = 特征值之和\n\n")

cat("2. 常用范数：\n")
cat("   - Frobenius: ||A||_F = sqrt(元素平方和) = sqrt(tr(A'A))\n")
cat("   - 谱范数:    ||A||_2 = 最大奇异值 = 最大拉伸倍数\n")
cat("   - 列和范数:  ||A||_1 = max(每列绝对值之和)\n")
cat("   - 行和范数:  ||A||_∞ = max(每行绝对值之和)\n")
cat("   - 核范数:    ||A||_* = 奇异值之和\n\n")

cat("3. 条件数 κ(A) = σ_max / σ_min\n")
cat("   - κ ≈ 1: 良态，误差不放大\n")
cat("   - κ >> 1: 病态，误差被放大\n\n")

cat("4. 实际应用：\n")
cat("   - 岭回归: 惩罚 ||x||_2² (L2正则化)\n")
cat("   - LASSO:  惩罚 ||x||_1  (产生稀疏解)\n")
cat("   - 矩阵补全: 惩罚 ||X||_* (核范数，逼迫低秩)\n")