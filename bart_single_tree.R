# ============================================================================
# BART (Bayesian Additive Regression Trees) with a Single Tree
# ============================================================================
#
# A from-scratch implementation of BART restricted to one tree (m=1).
# Uses Metropolis-Hastings within Gibbs: at each MCMC iteration we
#   1. Propose a tree mutation (grow / prune) via the Bayesian back-fitting
#      scheme of Chipman, George & McCulloch (2010).
#   2. Accept or reject the proposal with the marginal likelihood ratio.
#   3. Draw the terminal-node (leaf) parameters from their conjugate
#      normal posterior.
#   4. Draw the residual variance sigma^2 from its conjugate inverse-gamma
#      posterior.
#
# The prior on tree structure follows CGM (1998):
#   P(node is internal) = alpha * (1 + depth)^(-beta)
# with defaults alpha = 0.95, beta = 2.
#
# Leaf values ~ N(mu_mu, sigma_mu^2)  (prior)
# sigma^2     ~ InvGamma(nu/2, nu*lambda/2) (prior)
# ============================================================================

# ---------------------------------------------------------------------------
# Tree data structure
# ---------------------------------------------------------------------------
# Each tree is stored as a list of nodes keyed by integer id.
# Node fields:
#   id        – integer node id (root = 1)
#   parent    – parent id (NA for root)
#   left      – left child id  (NA if leaf)
#   right     – right child id (NA if leaf)
#   is_leaf   – logical
#   var       – splitting variable index (NA if leaf)
#   val       – splitting value          (NA if leaf)
#   mu        – leaf mean                (NA if internal)
#   depth     – depth in the tree (root = 0)

new_leaf <- function(id, parent, depth) {
  list(id = id, parent = parent, left = NA, right = NA,
       is_leaf = TRUE, var = NA, val = NA, mu = 0, depth = depth)
}

new_tree <- function() {
  tree <- list()
  tree[["1"]] <- new_leaf(1L, NA, 0L)
  tree
}

get_leaves <- function(tree) {
  Filter(function(nd) nd$is_leaf, tree)
}

get_internals <- function(tree) {
  Filter(function(nd) !nd$is_leaf, tree)
}

next_id <- function(tree) {
  max(as.integer(names(tree))) + 1L
}

# ---------------------------------------------------------------------------
# Predict: route observations through the tree
# ---------------------------------------------------------------------------
predict_tree <- function(tree, X) {
  # X is a matrix (n x p)
  n <- nrow(X)
  preds <- numeric(n)
  for (i in seq_len(n)) {
    nd <- tree[["1"]]
    while (!nd$is_leaf) {
      if (X[i, nd$var] <= nd$val) {
        nd <- tree[[as.character(nd$left)]]
      } else {
        nd <- tree[[as.character(nd$right)]]
      }
    }
    preds[i] <- nd$mu
  }
  preds
}

# Assign each observation to a leaf id
assign_leaves <- function(tree, X) {
  n <- nrow(X)
  leaf_ids <- integer(n)
  for (i in seq_len(n)) {
    nd <- tree[["1"]]
    while (!nd$is_leaf) {
      if (X[i, nd$var] <= nd$val) {
        nd <- tree[[as.character(nd$left)]]
      } else {
        nd <- tree[[as.character(nd$right)]]
      }
    }
    leaf_ids[i] <- nd$id
  }
  leaf_ids
}

# ---------------------------------------------------------------------------
# Tree prior  P(node is non-terminal) = alpha / (1 + depth)^beta
# ---------------------------------------------------------------------------
psplit <- function(depth, alpha, beta) {
  alpha / (1 + depth)^beta
}

# Log of the tree structure prior ratio for a GROW move
# We grow leaf at depth d into an internal node with two children at d+1.
# Ratio = p_split(d) * (1-p_split(d+1))^2 / (1-p_split(d))
log_tree_ratio_grow <- function(depth, alpha, beta) {
  log(psplit(depth, alpha, beta)) +
    2 * log(1 - psplit(depth + 1, alpha, beta)) -
    log(1 - psplit(depth, alpha, beta))
}

# ---------------------------------------------------------------------------
# Marginal likelihood of data in a set of leaf observations
# Under the conjugate normal model:
#   y_i | mu ~ N(mu, sigma^2),  mu ~ N(mu_mu, sigma_mu^2)
# The marginal log-likelihood (integrated over mu) is standard.
# ---------------------------------------------------------------------------
log_marginal_leaf <- function(y_leaf, sigma2, sigma_mu2, mu_mu) {
  n_l <- length(y_leaf)
  if (n_l == 0) return(0)
  sum_y <- sum(y_leaf)
  tau <- 1 / sigma2
  tau_mu <- 1 / sigma_mu2
  post_var <- 1 / (n_l * tau + tau_mu)
  post_mean <- post_var * (tau * sum_y + tau_mu * mu_mu)
  # log marginal = -n/2 log(2pi) - (n-1)/2 log(sigma2)
  #                -0.5 log(n*sigma_mu2 + sigma2)
  #                -0.5*(SS - n_l*ybar^2*tau + mu_mu^2*tau_mu - post_mean^2/post_var)
  # We use the standard identity directly:
  ll <- -0.5 * n_l * log(2 * pi * sigma2) +
    0.5 * log(post_var) - 0.5 * log(sigma_mu2) -
    0.5 * tau * sum(y_leaf^2) -
    0.5 * tau_mu * mu_mu^2 +
    0.5 * post_mean^2 / post_var
  ll
}

# ---------------------------------------------------------------------------
# GROW proposal: pick a random leaf, split it
# ---------------------------------------------------------------------------
propose_grow <- function(tree, X, alpha, beta) {
  leaves <- get_leaves(tree)
  if (length(leaves) == 0) return(list(tree = tree, log_trans = -Inf))

  # Pick a random leaf
  leaf <- leaves[[sample.int(length(leaves), 1)]]
  depth <- leaf$depth

  p <- ncol(X)
  var_idx <- sample.int(p, 1)

  # Available split values: unique values of X[, var_idx] that fall in this

  # leaf's observations (we pass the full X and let the caller filter, but

  # here we just pick from all unique values for simplicity)
  vals <- sort(unique(X[, var_idx]))
  if (length(vals) < 2) return(list(tree = tree, log_trans = -Inf))
  # pick a split point uniformly among interior gaps
  cut_idx <- sample.int(length(vals) - 1, 1)
  cut_val <- (vals[cut_idx] + vals[cut_idx + 1]) / 2

  # Build new tree
  new_tree_obj <- tree
  lid <- next_id(new_tree_obj)
  rid <- lid + 1L

  # Update the chosen leaf -> internal
  new_tree_obj[[as.character(leaf$id)]]$is_leaf <- FALSE
  new_tree_obj[[as.character(leaf$id)]]$var <- var_idx
  new_tree_obj[[as.character(leaf$id)]]$val <- cut_val
  new_tree_obj[[as.character(leaf$id)]]$mu <- NA
  new_tree_obj[[as.character(leaf$id)]]$left <- lid
  new_tree_obj[[as.character(leaf$id)]]$right <- rid

  # Create two new leaves
  new_tree_obj[[as.character(lid)]] <- new_leaf(lid, leaf$id, depth + 1L)
  new_tree_obj[[as.character(rid)]] <- new_leaf(rid, leaf$id, depth + 1L)

  # Transition ratio: P(prune picks this node back) / P(grow picked this leaf)
  n_leaves_old <- length(leaves)
  # After grow, the number of "prunable" internal nodes (both children are leaves)
  internals_new <- get_internals(new_tree_obj)
  n_prunable <- sum(sapply(internals_new, function(nd) {
    lc <- new_tree_obj[[as.character(nd$left)]]
    rc <- new_tree_obj[[as.character(nd$right)]]
    lc$is_leaf && rc$is_leaf
  }))

  # log transition ratio  (q_prune / q_grow)
  # q_grow  = 0.5 * (1/n_leaves_old) * (1/p) * (1/(n_cuts))
  # q_prune = 0.5 * (1/n_prunable)
  n_cuts <- length(vals) - 1
  log_trans <- log(n_leaves_old) + log(p) + log(n_cuts) - log(n_prunable)

  # log tree structure prior ratio
  log_prior <- log_tree_ratio_grow(depth, alpha, beta)

  list(tree = new_tree_obj, log_trans = log_trans, log_prior = log_prior,
       old_leaf_id = leaf$id)
}

# ---------------------------------------------------------------------------
# PRUNE proposal: pick a random "prunable" internal node, collapse it
# ---------------------------------------------------------------------------
propose_prune <- function(tree, X, alpha, beta) {
  internals <- get_internals(tree)
  # Prunable = both children are leaves
  prunable <- Filter(function(nd) {
    lc <- tree[[as.character(nd$left)]]
    rc <- tree[[as.character(nd$right)]]
    lc$is_leaf && rc$is_leaf
  }, internals)

  if (length(prunable) == 0) return(list(tree = tree, log_trans = -Inf))

  node <- prunable[[sample.int(length(prunable), 1)]]
  depth <- node$depth

  new_tree_obj <- tree
  # Remove children
  new_tree_obj[[as.character(node$left)]] <- NULL
  new_tree_obj[[as.character(node$right)]] <- NULL
  # Convert back to leaf
  new_tree_obj[[as.character(node$id)]]$is_leaf <- TRUE
  new_tree_obj[[as.character(node$id)]]$left <- NA
  new_tree_obj[[as.character(node$id)]]$right <- NA
  new_tree_obj[[as.character(node$id)]]$var <- NA
  new_tree_obj[[as.character(node$id)]]$val <- NA
  new_tree_obj[[as.character(node$id)]]$mu <- 0

  n_prunable_old <- length(prunable)
  n_leaves_new <- length(get_leaves(new_tree_obj))

  p <- ncol(X)
  var_idx <- node$var
  vals <- sort(unique(X[, var_idx]))
  n_cuts <- max(length(vals) - 1, 1)

  # log transition ratio  (q_grow / q_prune)  – inverse of grow

  log_trans <- log(n_prunable_old) - log(n_leaves_new) - log(p) - log(n_cuts)

  log_prior <- -log_tree_ratio_grow(depth, alpha, beta)

  list(tree = new_tree_obj, log_trans = log_trans, log_prior = log_prior)
}

# ---------------------------------------------------------------------------
# Draw leaf parameters from their full conditional posterior
# ---------------------------------------------------------------------------
draw_leaf_params <- function(tree, X, R, sigma2, sigma_mu2, mu_mu) {
  leaf_ids <- assign_leaves(tree, X)
  leaves <- get_leaves(tree)
  for (leaf in leaves) {
    idx <- which(leaf_ids == leaf$id)
    n_l <- length(idx)
    if (n_l == 0) {
      tree[[as.character(leaf$id)]]$mu <- rnorm(1, mu_mu, sqrt(sigma_mu2))
      next
    }
    y_leaf <- R[idx]
    tau <- 1 / sigma2
    tau_mu <- 1 / sigma_mu2
    post_var <- 1 / (n_l * tau + tau_mu)
    post_mean <- post_var * (tau * sum(y_leaf) + tau_mu * mu_mu)
    tree[[as.character(leaf$id)]]$mu <- rnorm(1, post_mean, sqrt(post_var))
  }
  tree
}

# ---------------------------------------------------------------------------
# Draw sigma^2 from its full conditional (inverse-gamma)
# ---------------------------------------------------------------------------
draw_sigma2 <- function(R_vec, nu, lambda) {
  n <- length(R_vec)
  shape <- (nu + n) / 2
  rate  <- (nu * lambda + sum(R_vec^2)) / 2
  1 / rgamma(1, shape = shape, rate = rate)
}

# ---------------------------------------------------------------------------
# Main BART function (single tree, m = 1)
# ---------------------------------------------------------------------------
#' @param X       numeric matrix of predictors (n x p)
#' @param y       numeric response vector (length n)
#' @param n_iter  total MCMC iterations
#' @param n_burn  burn-in iterations to discard
#' @param alpha   tree prior parameter (default 0.95)
#' @param beta    tree prior parameter (default 2)
#' @param k       prior shrinkage: sigma_mu = (ymax-ymin) / (2*k)
#' @param nu      prior df for sigma^2
#' @param q       quantile of the prior on sigma^2 (used to set lambda)
#' @return list with posterior samples of predictions and sigma
bart_single_tree <- function(X, y,
                             n_iter = 1000, n_burn = 200,
                             alpha = 0.95, beta = 2,
                             k = 2, nu = 3, q = 0.90) {
  # Coerce to matrix
  if (!is.matrix(X)) X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  # Centre y
  y_mean <- mean(y)
  y_c <- y - y_mean

  # Prior calibration
  y_range <- max(y_c) - min(y_c)
  if (y_range == 0) y_range <- 1
  sigma_mu <- y_range / (2 * k)
  sigma_mu2 <- sigma_mu^2
  mu_mu <- 0  # centred

  # sigma^2 prior: choose lambda so that P(sigma < s_hat) = q
  # where s_hat = sd(y_c)
  s_hat <- sd(y_c)
  if (s_hat == 0) s_hat <- 1
  # Use qchisq to calibrate lambda
  lambda <- s_hat^2 * qchisq(1 - q, df = nu) / nu

  # Initialise
  tree <- new_tree()
  sigma2 <- s_hat^2

  # Storage
  n_keep <- n_iter - n_burn
  pred_samples <- matrix(NA, nrow = n_keep, ncol = n)
  sigma_samples <- numeric(n_keep)
  tree_size <- integer(n_iter)

  for (iter in seq_len(n_iter)) {
    # Current residual (with 1 tree R = y_c - tree_pred, but tree_pred IS
    # our only tree so partial residual = y_c)
    R <- y_c

    # --- Step 1: Propose tree structure change (grow or prune) ---
    move <- sample(c("grow", "prune"), 1)
    if (move == "grow") {
      proposal <- propose_grow(tree, X, alpha, beta)
    } else {
      proposal <- propose_prune(tree, X, alpha, beta)
    }

    if (is.finite(proposal$log_trans)) {
      # Compute log marginal likelihood ratio
      leaf_ids_old <- assign_leaves(tree, X)
      leaf_ids_new <- assign_leaves(proposal$tree, X)

      log_lik_old <- 0
      for (lf in get_leaves(tree)) {
        idx <- which(leaf_ids_old == lf$id)
        if (length(idx) > 0)
          log_lik_old <- log_lik_old +
            log_marginal_leaf(R[idx], sigma2, sigma_mu2, mu_mu)
      }

      log_lik_new <- 0
      for (lf in get_leaves(proposal$tree)) {
        idx <- which(leaf_ids_new == lf$id)
        if (length(idx) > 0)
          log_lik_new <- log_lik_new +
            log_marginal_leaf(R[idx], sigma2, sigma_mu2, mu_mu)
      }

      log_alpha <- log_lik_new - log_lik_old +
        proposal$log_prior + proposal$log_trans

      if (log(runif(1)) < log_alpha) {
        tree <- proposal$tree
      }
    }

    # --- Step 2: Draw leaf parameters ---
    tree <- draw_leaf_params(tree, X, R, sigma2, sigma_mu2, mu_mu)

    # --- Step 3: Draw sigma^2 ---
    fitted <- predict_tree(tree, X)
    resid <- y_c - fitted
    sigma2 <- draw_sigma2(resid, nu, lambda)

    tree_size[iter] <- length(get_leaves(tree))

    # Store post-burn-in
    if (iter > n_burn) {
      j <- iter - n_burn
      pred_samples[j, ] <- fitted + y_mean  # un-centre
      sigma_samples[j] <- sqrt(sigma2)
    }
  }

  list(
    y_hat_mean  = colMeans(pred_samples),
    y_hat_samples = pred_samples,
    sigma_samples = sigma_samples,
    tree_size     = tree_size,
    y_mean        = y_mean
  )
}
