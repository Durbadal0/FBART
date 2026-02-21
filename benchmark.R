# ============================================================================
# Benchmark: Remote Server vs Local Machine
# Run this script on both machines and compare the timings.
# Usage: Rscript benchmark.R
# ============================================================================

cat("============================================\n")
cat("  Machine Benchmark\n")
cat("============================================\n\n")

# --- System info ---
cat(">> System Info:\n")
cat(sprintf("   OS:       %s\n", Sys.info()["sysname"]))
cat(sprintf("   Machine:  %s\n", Sys.info()["machine"]))
cat(sprintf("   Node:     %s\n", Sys.info()["nodename"]))
cat(sprintf("   R version: %s\n", R.version.string))
cat("\n")

# --- Test 1: Matrix multiplication (CPU + memory) ---
cat(">> Test 1: Matrix multiplication (1000x1000) ...\n")
set.seed(1)
A <- matrix(rnorm(1000 * 1000), 1000, 1000)
B <- matrix(rnorm(1000 * 1000), 1000, 1000)
t1 <- system.time({ C <- A %*% B })
cat(sprintf("   Time: %.3f seconds\n\n", t1["elapsed"]))

# --- Test 2: Sorting a large vector ---
cat(">> Test 2: Sorting 10 million random numbers ...\n")
set.seed(2)
x <- rnorm(1e7)
t2 <- system.time({ sort(x) })
cat(sprintf("   Time: %.3f seconds\n\n", t2["elapsed"]))

# --- Test 3: Loop-heavy computation (Fibonacci-style) ---
cat(">> Test 3: Loop-heavy computation (5 million iterations) ...\n")
t3 <- system.time({
  a <- 0; b <- 1
  for (i in seq_len(5e6)) {
    tmp <- a + b
    a <- b
    b <- tmp
  }
})
cat(sprintf("   Time: %.3f seconds\n\n", t3["elapsed"]))

# --- Test 4: Linear regression on large dataset ---
cat(">> Test 4: Linear regression (n=100000, p=50) ...\n")
set.seed(3)
n <- 1e5; p <- 50
X <- matrix(rnorm(n * p), n, p)
y <- X %*% rnorm(p) + rnorm(n)
t4 <- system.time({ fit <- lm.fit(X, y) })
cat(sprintf("   Time: %.3f seconds\n\n", t4["elapsed"]))

# --- Test 5: Monte Carlo simulation ---
cat(">> Test 5: Monte Carlo pi estimation (10 million samples) ...\n")
set.seed(4)
t5 <- system.time({
  u <- runif(1e7)
  v <- runif(1e7)
  pi_est <- 4 * mean(u^2 + v^2 <= 1)
})
cat(sprintf("   Time: %.3f seconds  (pi ~ %.5f)\n\n", t5["elapsed"], pi_est))

# --- Test 6: SVD of a matrix ---
cat(">> Test 6: SVD of 500x500 matrix ...\n")
set.seed(5)
M <- matrix(rnorm(500 * 500), 500, 500)
t6 <- system.time({ svd(M) })
cat(sprintf("   Time: %.3f seconds\n\n", t6["elapsed"]))

# --- Summary ---
total <- t1["elapsed"] + t2["elapsed"] + t3["elapsed"] +
         t4["elapsed"] + t5["elapsed"] + t6["elapsed"]

cat("============================================\n")
cat("  SUMMARY\n")
cat("============================================\n")
cat(sprintf("  1. Matrix multiply:    %.3f s\n", t1["elapsed"]))
cat(sprintf("  2. Sort 10M:           %.3f s\n", t2["elapsed"]))
cat(sprintf("  3. Loop (5M iters):    %.3f s\n", t3["elapsed"]))
cat(sprintf("  4. Linear regression:  %.3f s\n", t4["elapsed"]))
cat(sprintf("  5. Monte Carlo:        %.3f s\n", t5["elapsed"]))
cat(sprintf("  6. SVD 500x500:        %.3f s\n", t6["elapsed"]))
cat(sprintf("  ---------------------------------\n"))
cat(sprintf("  TOTAL:                 %.3f s\n", total))
cat("============================================\n")
