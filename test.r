library(MASS)

# 练习1
# 非方阵，m * n，m < n, 欠定方程，n 维度 映射到 m 维度， 无穷多解，亏秩
# b 一定在 Ax 张成的空间, 因为行满秩
A <- matrix(rnorm(12),2)
A
b <- rnorm(2)
b
# ATA 的秩和 A 的秩一样，这里为最大3（行数），所以不可逆，得用伪逆求解
r <- qr(t(A) %*% A)$rank
r
# 求A的零空间，N 是 A 的零空间基向量
N <- MASS::Null(t(A))  # 零空间基向量
N
# 秩-零度定理
dim(MASS::Null(t(A)))
6 - r
# 求伪逆
Ap <- ginv(A)
# A 不可逆，所以方程没有唯一解
#solve(A,b)
x_hat <- Ap %*% b
x_hat
# 求 b_hat
b_hat <- A %*% x_hat
b_hat
# 求残差, 对于欠定方程，因为有无数解，所以
diff <- abs(b_hat - b)
diff
# 求L2范数，因为欠定方程，b 一定在 Ax 空间，所以残差 = 0
sum(diff^2)

# 练习2
# 非方阵，m * n，m > n, 超定方程，n 维度 映射到 m 维度， 无精确解，满秩
# 超定方程，列满秩的话，简单理解为就 N 个未知数，给了我一大堆（M）个条件，属于既要又要，所以，一般情况下无解
# 即 b 一般不在 Ax 张成的空间内
A <- matrix(rnorm(10),5)
A
b <- rnorm(5)
b
# 列满秩
r <- qr(t(A) %*% A)$rank
r
# 求A的零空间，N 是 A 的零空间基向量
N <- MASS::Null(t(A))  # 零空间基向量
N
2 - r
# 求伪逆
Ap <- ginv(A)
x_hat <- Ap %*% b
x_hat
# 求 b_hat
b_hat <- A %*% x_hat
b_hat
# 求残差, 对于超定方程，因为有无解，所以
diff <- abs(b_hat - b)
diff
# 求L2范数，因为超定方程，b 不一定在 Ax 空间，所以残差 > 0
sum(diff^2)


# 练习3
# 非方阵，超定方程，转为方阵
A <- matrix(rnorm(10),5)  # 5 x 2
A
b <- rnorm(5)
b
B <- t(A) %*% A
B
# 方法1：正规方程
# Ax = b ==> A^TA x = A^Tb => x = solve(ATA) * ATb
x2 <- solve(B) %*% t(A) %*% b 
x2
x3 <- solve(B,t(A) %*% b)  # 直接求解
x3
# 方法2：伪逆（通用）, 应该一致
Ap <- ginv(A)
x_hat <- Ap %*% b
x_hat



# 练习4
# 方阵，亏秩（不满秩）情况
A <- matrix(rnorm(12),3)  # 3 x 4
A
b <- rnorm(3)
b
# 验证亏秩
r <- qr(A)$rank
r  # 秩为 3, 不是 4

# 行列式为 0，不可逆
det(t(A) %*% A)

# 方法1：正规方程： Ax = b ==> A^TA x = A^Tb => x = solve(ATA) * ATb 
# 因为 ATA 不可逆，所以无法使用正规方程

# 方法2：只能用伪逆
Ap <- ginv(A)
x_hat <- Ap %*% b
x_hat
# 检查残差
b_hat <- A %*% x_hat
sum((b_hat - b)^2)  # 可能为 0（b 在列空间），也可能大于 0（b 不在列空间）

# 零空间
N <- MASS::Null(t(A))
N  # 1 维，因为 4 - 3 = 1


# 练习5
# 方阵但亏秩
A <- matrix(c(1, 2, 2, 4), 2, 2)  # 两列成比例
A
b <- rnorm(2)

qr(A)$rank     # 1，亏秩,不可逆
det(A)         # 0

# 不能用 solve
# solve(A, b)  # 报错

# 只能用伪逆
x_hat <- ginv(A) %*% b
b_hat <- A %*% x_hat
sum((b_hat - b)^2)  # 可能 > 0

# 零空间
N <- MASS::Null(t(A))
dim(N)  # 2 - 1 = 1






# 练习6 - SVD
# SVD图像压缩示例
library(jpeg)

# 1. 读取图像
img <- readJPEG("photo.jpg")
nrow(img)
ncol(img)
# 2. 分离RGB三个通道并做SVD
svd_R <- svd(img[,,1])
svd_G <- svd(img[,,2])
svd_B <- svd(img[,,3])
# 奇异值个数就等于宽度  数学上，m×n 矩阵的奇异值个数是 min(m, n)
length(svd_G$d)
# 3. 设置保留的奇异值数量
k <- 80

# 4. 重构压缩图像
# svd() 返回的奇异值是从大到小排好序的
compress <- function(s, k) {
  # 取矩阵 s$u 的前 k 列, U 的每一列对应一个奇异值, 同理 s$v[,1:k] 是取 V 的前 k 列
  s$u[,1:k] %*% diag(s$d[1:k]) %*% t(s$v[,1:k])
}

# 创建一个和原图同样大小的空数组，填满0
# 原图要存整个矩阵，m×n 个数
# SVD压缩后只需要存三样东西   
# U 的前 k 列：m×k 个数
# 前 k 个奇异值：k 个数
# V 的前 k 列：n×k 个数
# 总共 k×(m + n + 1) 个数，远小于 m×n
# 前k个奇异值 ，奇异向量，就是低秩近似
img_compressed <- array(0, dim = dim(img))
img_compressed[,,1] <- pmin(pmax(compress(svd_R, k), 0), 1)
img_compressed[,,2] <- pmin(pmax(compress(svd_G, k), 0), 1)
img_compressed[,,3] <- pmin(pmax(compress(svd_B, k), 0), 1)

# 5. 并排显示原图和压缩图
par(mfrow = c(1, 2), mar = c(2, 1, 2, 1))

# 创建一个空白画布 c(0, 1), c(0, 1) 设定 x 和 y 轴的范围都是 0 到 1 
# type = "n" 表示不画任何点，只建立坐标系
# axes = FALSE 不显示坐标轴
# xlab = "", ylab = "" 不显示轴标签
plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "", main = "原图")
rasterImage(img, 0, 0, 1, 1)

plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "", main = paste0("压缩 (k=", k, ")"))
rasterImage(img_compressed, 0, 0, 1, 1)

m <- dim(img)[1]
n <- dim(img)[2]
cat("压缩率:", round(m * n / (k * (m + n + 1)), 1), "倍\n")



# 用户-电影评分矩阵（NA表示没看过）
ratings <- matrix(c(
  5, 4, NA, 1, 1,
  4, 5, 4, 1, NA,
  NA, 4, 5, 1, 1,
  1, 1, 1, 5, 4,
  1, NA, 1, 4, 5
), nrow = 5, byrow = TRUE)

rownames(ratings) <- c("小明", "小红", "小华", "小李", "小张")
colnames(ratings) <- c("泰坦尼克", "罗马假日", "傲慢与偏见", "变形金刚", "速度与激情")

cat("原始评分矩阵：\n")
print(ratings)

# 用均值填充缺失值（简化处理）
ratings_filled <- ratings
ratings_filled[is.na(ratings_filled)] <- mean(ratings, na.rm = TRUE)

# SVD分解
s <- svd(ratings_filled)

# 只保留前2个奇异值（提取2个隐因子：爱情片 vs 动作片）
k <- 2
ratings_predicted <- s$u[,1:k] %*% diag(s$d[1:k]) %*% t(s$v[,1:k])
rownames(ratings_predicted) <- rownames(ratings)
colnames(ratings_predicted) <- colnames(ratings)

cat("\nSVD预测评分：\n")
print(round(ratings_predicted, 1))

# 找出原本缺失的位置，这些就是推荐预测
cat("\n预测结果：\n")
cat("小明对《傲慢与偏见》的预测评分:", round(ratings_predicted["小明", "傲慢与偏见"], 1), "\n")
cat("小红对《速度与激情》的预测评分:", round(ratings_predicted["小红", "速度与激情"], 1), "\n")
cat("小华对《泰坦尼克》的预测评分:", round(ratings_predicted["小华", "泰坦尼克"], 1), "\n")
cat("小李对《罗马假日》的预测评分:", round(ratings_predicted["小李", "罗马假日"], 1), "\n")

# 练习7 特征分解
# 随机矩阵几乎肯定能对角化，但特征值和特征向量可能是复数
# 能对角化不等于好用：
# 特征值是复数，解释起来不直观
# 特征向量也是复数，不一定垂直
# 但"能对角化"和"对角化好用"是两回事。实对称矩阵最好：特征值是实数，特征向量垂直，干干净净
# 这就是为什么 SVD 更实用——奇异值永远是非负实数
# U 和 V 永远是实正交矩阵（实矩阵的话）
A <- matrix(rnorm(16),4)
e <- eigen(A)
e

# 练习7 梯度下降对比

m <- 10000
n <- 1000
A <- matrix(rnorm(m*n), nrow=m)
b <- rnorm(m)

# 方法1：ginv 伪逆
cat("开始计算 ginv...\n")
t1 <- system.time({
  Ap <- ginv(A)
  x_ginv <- Ap %*% b
})
loss_ginv <- sum((A %*% x_ginv - b)^2)
cat("ginv 耗时:", t1[3], "秒\n")
cat("ginv 损失:", loss_ginv, "\n\n")

# 方法2：lm.fit QR分解
cat("开始计算 lm.fit...\n")
t2 <- system.time({
  x_lm <- lm.fit(A, b)$coefficients
})
loss_lm <- sum((A %*% x_lm - b)^2)
cat("lm.fit 耗时:", t2[3], "秒\n")
cat("lm.fit 损失:", loss_lm, "\n\n")

# 方法3：SGD
cat("开始计算 SGD...\n")
t3 <- system.time({
  x <- rep(0, n)
  lr <- 0.001
  epochs <- 200
  batch_size <- 512
  
  for (ep in 1:epochs) {
    idx <- sample(m)
    for (i in seq(1, m, batch_size)) {
      batch_idx <- idx[i:min(i+batch_size-1, m)]
      A_batch <- A[batch_idx, , drop=FALSE]
      b_batch <- b[batch_idx]
      residual <- A_batch %*% x - b_batch
      grad <- 2 * t(A_batch) %*% residual / length(batch_idx)
      x <- x - lr * grad
    }
    if (ep %% 10 == 0) {
      cat("Epoch", ep, "Loss:", sum((A %*% x - b)^2), "\n")
    }
  }
  x_sgd <- x
})
loss_sgd <- sum((A %*% x_sgd - b)^2)
cat("SGD 耗时:", t3[3], "秒\n")
cat("SGD 损失:", loss_sgd, "\n\n")

# 汇总对比
cat("===== 汇总对比 =====\n")
cat("方法\t\t耗时(秒)\t损失\n")
cat("ginv\t\t", t1[3], "\t\t", loss_ginv, "\n")
cat("lm.fit\t\t", t2[3], "\t\t", loss_lm, "\n")
cat("SGD\t\t", t3[3], "\t\t", loss_sgd, "\n")