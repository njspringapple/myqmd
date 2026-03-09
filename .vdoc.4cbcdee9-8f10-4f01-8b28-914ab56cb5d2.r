#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: data-gen
#| fig-cap: "模拟数据散点图"
#| fig-width: 6
#| fig-height: 4

set.seed(123)
n <- 100
x <- runif(n, -2, 2)
true_w <- 3.5
true_b <- 1.2
y <- true_w * x + true_b + rnorm(n, 0, 0.6)

library(ggplot2)
ggplot(data.frame(x = x, y = y), aes(x, y)) +
  geom_point(color = "steelblue", alpha = 0.7) +
  geom_abline(slope = true_w, intercept = true_b, color = "red", linetype = "dashed") +
  labs(title = "模拟数据（真实直线为红色虚线）", x = "x", y = "y") +
#
#
#
