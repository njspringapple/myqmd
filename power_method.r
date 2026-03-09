power_method <- function(A, tol = 1e-8, max_iter = 1000) {
  v <- rnorm(nrow(A))  # 随机初始向量
  v <- v / sqrt(sum(v^2))  # 归一化
  
  for (i in 1:max_iter) {
    w <- A %*% v  # 乘以矩阵
    lambda <- sum(v * w)  # Rayleigh商估计特征值
    v_new <- w / sqrt(sum(w^2))  # 归一化
    
    if (sqrt(sum((v_new - v)^2)) < tol) break
    v <- v_new
  }
  
  list(value = lambda, vector = as.vector(v), iter = i)
}

# 测试
A <- matrix(c(4, 1, 1, 3), 2, 2)
power_method(A)
eigen(A)$values[1]  # 对比验证






# 原矩阵 A（旧坐标系下的变换）
A <- matrix(c(4, 1, 2, 3), nrow = 2, byrow = TRUE)

# 求特征值和特征向量
eig <- eigen(A)
eigenvalues <- eig$values
P <- eig$vectors

cat("特征值:", eigenvalues, "\n")  # [5, 2]

# 相似变换：B = P⁻¹AP
P_inv <- solve(P)
B <- P_inv %*% A %*% P

cat("\n原矩阵 A:\n")
print(A)

cat("\n过渡矩阵 P（特征向量）:\n")
print(P)

cat("\n相似变换后 B = P⁻¹AP:\n")
print(round(B, 10))  # 对角阵，对角线就是特征值

# 验证：A 和 B 的不变量相同
cat("\n验证不变量:\n")
cat(sprintf("det(A) = %.1f, det(B) = %.1f\n", det(A), det(B)))
cat(sprintf("trace(A) = %.1f, trace(B) = %.1f\n", sum(diag(A)), sum(diag(B))))