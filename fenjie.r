set.seed(42)
n <- 100 # 观测数量
p <- 3 # 变量数量

# y = 2 - x1 + 0.5x2
# 设计矩阵 100 * 2
X <- cbind(1,matrix(rnorm(n * (p - 1)), n, p - 1))
# 真实系数
beta_true <- c(2,-1,0.5)
# 响应变量（加噪声）
b <- X %*% beta_true + rnorm(n, sd = 0.5)

# ============================================================
# 方法1：QR分解
# ============================================================
qr_decomp <- qr(X)
Q <- qr.Q(qr_decomp)
R <- qr.R(qr_decomp)

cat("QR还原误差 ：" , norm(X - Q %*% R, "F") , "\n" )

# 求解 Rx = Q'b
# 可以用solve来求解，但是计算复杂度O(n^3)
# backsolve直接利用上三角结构求解，计算量O(n^2)，快了一个量级
# 利用QR分解已经把矩阵化成了上三角矩阵
beta_qr <- backsolve(R,t(Q) %*% b)
cat("QR解：", round(beta_qr,4),"\n")

# ============================================================
# 方法2：Cholesky分解
# ============================================================
# 先构造方程 X'X β = X'b
XtX <- t(X) %*% X
Xtb <- t(X) %*% b

L <- chol(XtX)
cat("Cholesky还原误差:", norm(XtX - t(L) %*% L, "F"), "\n")

# 两步求解：先解 L'y = X'b，再解 Ly = beta
y_tmp <- forwardsolve(t(L), Xtb)
beta_chol <- backsolve(L, y_tmp)
cat("Cholesky解:", round(beta_chol, 4), "\n")

# ============================================================
# 方法3：SVD分解
# ============================================================
sv <- svd(X)
U <- sv$u
D <- sv$d
V <- sv$v

# 验证 X = U diag(D) V'
# 还原误差取F范数，用F范数可以直接看出丢弃了多少"能量"
cat("SVD还原误差：", norm(X - U %*% diag(D) %*% t(V), "F"), "\n")

# 伪逆求解：beta = V * D^{-1} * U'b
beta_svd <- V %*% (t(U) %*% b / D)
cat("SVD解:", round(beta_svd, 4), "\n")

# ============================================================
# 对比验证
# ============================================================
beta_lm <- coef(lm(b ~ X - 1))   # lm()内置结果

cat("\n===== 对比 =====\n")
cat("真实值:    ", beta_true, "\n")
cat("QR解:      ", round(beta_qr, 4), "\n")
cat("Cholesky解:", round(beta_chol, 4), "\n")
cat("SVD解:     ", round(beta_svd, 4), "\n")
cat("lm()解:    ", round(beta_lm, 4), "\n")


A <- matrix(rnorm(100*100),nrow = 100)
qr(A)$rank
sv <- svd(A)
sv$d



set.seed(42)
# 4x4 正交矩阵
Q <- qr.Q(qr(matrix(rnorm(4*4),nrow=4)))
# rank 为 2 
sigma <- diag(c(5,3,0,0))
A <- Q %*% sigma %*% t(Q)
sum(svd(A)$d >  1e-10)
