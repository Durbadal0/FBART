# ============================================================================
# Test script for BART single-tree implementation
# ============================================================================
source("bart_single_tree.R")

cat("=== Test 1: Friedman function (non-linear regression) ===\n")
set.seed(42)
n <- 200
p <- 5
X <- matrix(runif(n * p), nrow = n, ncol = p)
# Friedman #1 (first 5 variables matter)
y_true <- 10 * sin(pi * X[, 1] * X[, 2]) + 20 * (X[, 3] - 0.5)^2 +
  10 * X[, 4] + 5 * X[, 5]
y <- y_true + rnorm(n, sd = 1)

fit <- bart_single_tree(X, y, n_iter = 500, n_burn = 100)

rmse <- sqrt(mean((fit$y_hat_mean - y_true)^2))
cor_val <- cor(fit$y_hat_mean, y_true)
cat(sprintf("  RMSE (vs true f):  %.3f\n", rmse))
cat(sprintf("  Correlation:       %.3f\n", cor_val))
cat(sprintf("  Mean sigma:        %.3f  (true = 1.0)\n", mean(fit$sigma_samples)))
cat(sprintf("  Prediction range:  [%.2f, %.2f]\n",
            min(fit$y_hat_mean), max(fit$y_hat_mean)))

stopifnot(cor_val > 0.5)
cat("  PASSED: correlation with truth > 0.5\n\n")

cat("=== Test 2: Simple linear signal y = 3*x + noise ===\n")
set.seed(123)
n2 <- 150
x2 <- matrix(runif(n2), ncol = 1)
y2_true <- 3 * x2[, 1]
y2 <- y2_true + rnorm(n2, sd = 0.3)

fit2 <- bart_single_tree(x2, y2, n_iter = 500, n_burn = 100)

rmse2 <- sqrt(mean((fit2$y_hat_mean - y2_true)^2))
cor2 <- cor(fit2$y_hat_mean, y2_true)
cat(sprintf("  RMSE (vs true f):  %.3f\n", rmse2))
cat(sprintf("  Correlation:       %.3f\n", cor2))
cat(sprintf("  Mean sigma:        %.3f  (true = 0.3)\n", mean(fit2$sigma_samples)))

stopifnot(cor2 > 0.8)
cat("  PASSED: correlation with truth > 0.8\n\n")

cat("=== Test 3: Constant signal (should recover the mean) ===\n")
set.seed(7)
n3 <- 100
x3 <- matrix(rnorm(n3 * 2), ncol = 2)
y3 <- rep(5, n3) + rnorm(n3, sd = 0.5)

fit3 <- bart_single_tree(x3, y3, n_iter = 400, n_burn = 100)

mean_pred <- mean(fit3$y_hat_mean)
cat(sprintf("  Mean of predictions: %.3f  (true mean = 5.0)\n", mean_pred))
cat(sprintf("  Pred std dev:        %.3f  (should be small)\n", sd(fit3$y_hat_mean)))

stopifnot(abs(mean_pred - 5) < 1.0)
cat("  PASSED: mean prediction within 1.0 of true mean\n\n")

cat("=== Test 4: Output dimensions are correct ===\n")
n_iter <- 300
n_burn <- 50
fit4 <- bart_single_tree(x3, y3, n_iter = n_iter, n_burn = n_burn)
stopifnot(nrow(fit4$y_hat_samples) == n_iter - n_burn)
stopifnot(ncol(fit4$y_hat_samples) == n3)
stopifnot(length(fit4$sigma_samples) == n_iter - n_burn)
stopifnot(length(fit4$tree_size) == n_iter)
cat("  PASSED: output dimensions are correct\n\n")

cat("All tests passed!\n")
