expand_Y <- function(mY_int, ivItemcat) {
  # mY_int: N x H matrix of integer category values (0-based)
  # ivItemcat: length-H vector of number of categories per item
  N <- nrow(mY_int)
  H <- ncol(mY_int)
  out <- matrix(0, N, sum(ivItemcat))
  col_start <- 1L
  for (h in seq_len(H)) {
    K_h <- ivItemcat[h]
    for (k in seq_len(K_h)) {
      out[, col_start + k - 1L] <- as.integer(mY_int[, h] == (k - 1L))
    }
    col_start <- col_start + K_h
  }
  out
}

expand_Phi <- function(phi_mat, ivItemcat) {
  dichotomous <- ivItemcat == 2L
  result <- vector("list", length(ivItemcat))
  h_phi <- 1L
  for (h in seq_along(ivItemcat)) {
    if (dichotomous[h]) {
      result[[h]] <- rbind(1 - phi_mat[h_phi, ], phi_mat[h_phi, ])
      h_phi <- h_phi + 1L
    } else {
      K_h <- ivItemcat[h]
      result[[h]] <- phi_mat[h_phi:(h_phi + K_h - 1L), , drop = FALSE]
      h_phi <- h_phi + K_h
    }
  }
  do.call(rbind, result)
}

expand_Phi_free <- function(phi_free, ivItemcat) {
  result <- vector("list", length(ivItemcat))
  h_phi <- 1L
  for (h in seq_along(ivItemcat)) {
    K_h <- ivItemcat[h]
    if (K_h == 2L) {
      result[[h]] <- phi_free[h_phi, , drop = FALSE]
      h_phi <- h_phi + 1L
    } else {
      rows_h <- phi_free[h_phi:(h_phi + K_h - 2L), , drop = FALSE]
      result[[h]] <- rbind(1 - colSums(rows_h), rows_h)
      h_phi <- h_phi + K_h - 1L
    }
  }
  do.call(rbind, result)
}

log_lik_matrix <- function(Y, mPhi, mDesign = NULL) {
  if (is.null(mDesign)) {
    mDesign <- matrix(1L, nrow(Y), ncol(Y))
  }
  (mDesign * Y) %*% log(mPhi) #+ (mDesign * (1 - Y)) %*% log(1 - mPhi)
}

joint_log_lik <- function(Y, Z, mPhi, gamma.coefs, mDesign = NULL) {
  if (is.null(mDesign)) {
    mDesign <- matrix(1L, nrow(Y), ncol(Y))
  }

  log_P_Y_given_X <- log_lik_matrix(Y, mPhi, mDesign)

  eta <- Z %*% gamma.coefs
  eta_full <- cbind(0, eta)

  row_maxes_eta <- apply(eta_full, 1, max)
  log_denom_eta <- row_maxes_eta + log(rowSums(exp(eta_full - row_maxes_eta)))
  log_P_X_given_Z <- eta_full - log_denom_eta

  log_joint_prob <- log_P_Y_given_X + log_P_X_given_Z

  row_maxes_joint <- apply(log_joint_prob, 1, max)
  log_marg_prob <- row_maxes_joint +
    log(rowSums(exp(log_joint_prob - row_maxes_joint)))

  sum(log_marg_prob)
}

compute_posteriors <- function(Y, mDesign, theta1, ivItemcat, T) {
  vPi_free <- theta1[1:(T - 1L)]
  vPi <- c(1 - sum(vPi_free), vPi_free)
  n_free <- sum(ifelse(ivItemcat == 2L, 1L, ivItemcat - 1L))
  phi_free <- matrix(theta1[T:(T + n_free * T - 1L)], nrow = n_free, ncol = T)

  mPhi <- expand_Phi(expand_Phi_free(phi_free, ivItemcat), ivItemcat)

  log_joint <- sweep(log_lik_matrix(Y, mPhi, mDesign), 2, log(vPi), "+")
  row_max <- apply(log_joint, 1, max)
  log_denom <- row_max + log(rowSums(exp(log_joint - row_max)))
  exp(log_joint - log_denom)
}

compute_pwx_adj <- function(
  Y.obs,
  fit0,
  ivItemcat,
  mDesign = NULL,
  use.modal.assignment = TRUE,
  pi_adj = NULL # N x T covariate-adjusted class probs, or NULL for flat vPi
) {
  N <- nrow(Y.obs)
  T <- ncol(fit0$mPhi)

  mPhi_exp <- expand_Phi(fit0$mPhi, ivItemcat)
  phi_clamped <- pmax(pmin(mPhi_exp, 1 - 1e-10), 1e-10)

  log_p_it <- if (is.null(mDesign)) {
    Y.obs %*% log(phi_clamped)
  } else {
    (mDesign * Y.obs) %*% log(phi_clamped)
  }

  # use adjusted or flat priors
  if (!is.null(pi_adj)) {
    log_prior <- log(pi_adj) # N x T, person-specific
  } else {
    log_prior <- matrix(log(fit0$vPi), nrow = N, ncol = T, byrow = TRUE)
  }

  log_joint <- log_p_it + log_prior # N x T
  row_max <- apply(log_joint, 1, max)
  post <- exp(log_joint - row_max - log(rowSums(exp(log_joint - row_max)))) # N x T posteriors

  w.is <- if (use.modal.assignment) {
    w <- matrix(0L, N, T)
    w[cbind(seq_len(N), max.col(post))] <- 1L
    w
  } else {
    post
  }

  p.wx_joint <- (t(w.is) %*% post) / N
  p.wx_mat <- sweep(p.wx_joint, 2, colSums(p.wx_joint), "/")

  list(
    post = post,
    w.is = w.is,
    p.wx_mat = p.wx_mat
  )
}

# -- Step 2: Posteriors and classification-error matrix -----------------------
lca_step2 <- function(
  Y.obs,
  fit0,
  n_classes,
  use.modal.assignment,
  boundary.tol,
  use.simple.cov,
  ivItemcat,
  mDesign = NULL
) {
  if (is.null(mDesign)) {
    mDesign <- matrix(1L, nrow(Y.obs), ncol(Y.obs))
  }

  N <- nrow(Y.obs)
  T <- n_classes
  K <- ncol(Y.obs) #sum(K_h)

  # Number of free rows per item in mPhi
  n_free <- ivItemcat - 1L

  starts <- c(
    1L,
    cumsum(ifelse(ivItemcat == 2L, 1L, ivItemcat))[-length(ivItemcat)] + 1L
  )

  free_idx <- unlist(mapply(
    \(s, K_h) if (K_h == 2L) s else (s + 1L):(s + K_h - 1L),
    starts,
    ivItemcat,
    SIMPLIFY = FALSE
  ))

  phi_free <- fit0$mPhi[free_idx, ]

  theta1 <- c(
    fit0$vPi[2:T],
    phi_free
  )

  p.xy <- compute_posteriors(Y.obs, mDesign, theta1, ivItemcat, T)
  assignment <- max.col(p.xy)

  make_w <- function(posteriors) {
    if (use.modal.assignment) {
      w <- matrix(0L, nrow = N, ncol = T)
      w[cbind(seq_len(N), max.col(posteriors))] <- 1L
      w
    } else {
      posteriors
    }
  }

  w.is <- make_w(p.xy)

  compute_pwx <- function(t1) {
    post <- compute_posteriors(Y.obs, mDesign, t1, ivItemcat, T)
    w_local <- make_w(post)
    p.wx_joint <- (t(w_local) %*% post) / N
    sweep(p.wx_joint, 2, colSums(p.wx_joint), "/")
  }

  theta2_from_theta1 <- function(th1) {
    # rho <- c(1 - sum(th1[1:(T - 1)]), th1[1:(T - 1)])
    # phi <- matrix(th1[T:length(th1)], nrow = K, ncol = T)
    p.wx_mat <- compute_pwx(th1)
    log_ref <- log(diag(p.wx_mat))
    gamma_mat <- sweep(log(p.wx_mat), 2, log_ref, "-")
    gamma_mat[row(gamma_mat) != col(gamma_mat)]
  }

  gamma_vec_to_pwx <- function(gamma_vec) {
    gamma_mat <- matrix(0, nrow = T, ncol = T)
    gamma_mat[row(gamma_mat) != col(gamma_mat)] <- gamma_vec
    exp_mat <- exp(gamma_mat)
    sweep(exp_mat, 2, colSums(exp_mat), "/")
  }

  theta2 <- theta2_from_theta1(theta1)
  p.wx_mat <- gamma_vec_to_pwx(theta2)

  J.2 <- NULL
  if (!use.simple.cov) {
    rho_to_logit <- function(rho) log(rho[-1] / rho[1])
    logit_to_rho <- function(alpha) {
      d <- 1 + sum(exp(alpha))
      c(1 / d, exp(alpha) / d)
    }
    phi_to_logit <- function(phi) log(phi / (1 - phi))
    logit_to_phi <- function(beta) 1 / (1 + exp(-beta))

    # Jacobian of log-ratio transform for a K-simplex:
    # u_k = log(phi_k / phi_1), k = 2..K  =>  du_k/dphi_j
    # Returns (K-1) x (K-1) matrix
    simplex_jac_free <- function(phi_vec) {
      # phi_vec: full K_h-vector including reference phi_1
      K_h <- length(phi_vec)
      J <- matrix(1 / phi_vec[1L], K_h - 1L, K_h - 1L) # fill with 1/phi_1
      diag(J) <- diag(J) + 1 / phi_vec[-1L] # add 1/phi_{k+1} on diagonal
      J
    }

    theta1_to_u <- function(th1) {
      rho <- c(1 - sum(th1[1:(T - 1L)]), th1[1:(T - 1L)])
      phi_mat <- matrix(th1[T:length(th1)], nrow = nrow(phi_free), ncol = T)
      u_phi <- vector("list", length(ivItemcat))
      row_idx <- 1L
      for (h in seq_along(ivItemcat)) {
        K_h <- ivItemcat[h]
        if (K_h == 2L) {
          u_phi[[h]] <- phi_to_logit(phi_mat[row_idx, ])
          row_idx <- row_idx + 1L
        } else {
          rows_h <- row_idx:(row_idx + K_h - 2L)
          free_h <- phi_mat[rows_h, , drop = FALSE] # (K_h-1) x T
          # prepend reference category to get full simplex, then log-ratio transform
          full_h <- rbind(1 - colSums(free_h), free_h) # K_h x T
          u_phi[[h]] <- apply(full_h, 2, rho_to_logit) # (K_h-1) x T
          row_idx <- row_idx + K_h - 1L
        }
      }
      c(rho_to_logit(rho), unlist(u_phi))
    }

    u_to_theta1 <- function(u) {
      rho <- logit_to_rho(u[1:(T - 1L)])
      u_phi <- u[T:length(u)]
      phi_list <- vector("list", length(ivItemcat))
      idx <- 1L
      for (h in seq_along(ivItemcat)) {
        K_h <- ivItemcat[h]
        if (K_h == 2L) {
          phi_list[[h]] <- logit_to_phi(u_phi[idx:(idx + T - 1L)])
          idx <- idx + T
        } else {
          # (K_h-1)*T free params, reconstruct K_h x T block
          n_free_h <- (K_h - 1L) * T
          u_block <- matrix(u_phi[idx:(idx + n_free_h - 1L)], nrow = K_h - 1L)
          phi_block <- apply(u_block, 2, logit_to_rho) # K_h x T
          phi_list[[h]] <- phi_block[-1L, ] # drop ref row
          idx <- idx + n_free_h
        }
      }
      c(rho[-1L], do.call(rbind, phi_list))
    }

    J_theta1_to_u <- function(th1, eps = boundary.tol) {
      rho <- c(1 - sum(th1[1:(T - 1L)]), th1[1:(T - 1L)])
      phi_mat <- matrix(th1[T:length(th1)], nrow = nrow(phi_free), ncol = T)
      n <- length(th1)
      J <- matrix(0, n, n)

      # rho block: (T-1) x (T-1)
      rho1 <- rho[1L]
      for (k in seq_len(T - 1L)) {
        for (j in seq_len(T - 1L)) {
          J[k, j] <- (k == j) / rho[k + 1L] + 1 / rho1
        }
      }

      # phi blocks: one per item
      u_row <- T # current row in J (u-space), 1-indexed
      phi_row <- 1L # current row in phi_mat
      j_col <- T # current col in J (theta1-space)

      for (h in seq_along(ivItemcat)) {
        K_h <- ivItemcat[h]
        if (K_h == 2L) {
          # one free phi per class => T diagonal entries
          phi_h <- phi_mat[phi_row, ]
          is_bdry <- phi_h > (1 - eps) | phi_h < eps
          derivs <- ifelse(is_bdry, 0, 1 / (phi_h * (1 - phi_h)))
          for (t in seq_len(T)) {
            J[u_row + t - 1L, j_col + t - 1L] <- derivs[t]
          }
          u_row <- u_row + T
          phi_row <- phi_row + 1L
          j_col <- j_col + T
        } else {
          # K_h - 1 free phis per class => (K_h-1)*T block
          for (t in seq_len(T)) {
            phi_ht <- c(
              1 - sum(phi_mat[phi_row:(phi_row + K_h - 2L), t]),
              phi_mat[phi_row:(phi_row + K_h - 2L), t]
            )
            Jh <- simplex_jac_free(phi_ht) # (K_h-1) x (K_h-1)
            u_idx <- u_row + seq_len(K_h - 1L) - 1L
            j_idx <- j_col + seq_len(K_h - 1L) - 1L
            J[u_idx, j_idx] <- Jh
            u_row <- u_row + (K_h - 1L)
            j_col <- j_col + (K_h - 1L)
          }
          phi_row <- phi_row + K_h - 1L
        }
      }
      J
    }

    compute_J_unc_analytical <- function(
      p_ik,
      Y_obs,
      mDes,
      th1,
      ivItemcat,
      T_classes
    ) {
      N <- nrow(p_ik)

      # A_{st} = sum_i p_{is} p_{it}
      A <- t(p_ik) %*% p_ik
      A[A < 1e-12] <- 1e-12

      # Extract item probabilities to match the free parameter structure
      phi_mat <- matrix(
        th1[T_classes:length(th1)],
        nrow = sum(ivItemcat - 1L),
        ncol = T_classes
      )

      L_rho <- T_classes - 1L
      L_phi <- sum((ivItemcat - 1L) * T_classes)
      L <- L_rho + L_phi

      J <- matrix(0, nrow = T_classes * (T_classes - 1L), ncol = L)

      # Offsets for locating items and categories in the expanded matrices
      starts_Y <- c(1L, cumsum(ivItemcat)[-length(ivItemcat)] + 1L)
      item_offsets <- c(
        0L,
        cumsum((ivItemcat - 1L) * T_classes)[-length(ivItemcat)]
      )

      # Map (s, t) to the correct row in J (matching gamma_mat[row != col] column-major)
      st_idx <- 1L
      st_map <- matrix(0L, nrow = T_classes, ncol = T_classes)
      for (t in seq_len(T_classes)) {
        for (s in seq_len(T_classes)) {
          if (s != t) {
            st_map[s, t] <- st_idx
            st_idx <- st_idx + 1L
          }
        }
      }

      for (t in seq_len(T_classes)) {
        P_tt <- (p_ik[, t]^2) / A[t, t]

        for (s in seq_len(T_classes)) {
          if (s == t) {
            next
          }
          row_J <- st_map[s, t]
          P_st <- (p_ik[, s] * p_ik[, t]) / A[s, t]

          for (c_prime in seq_len(T_classes)) {
            I_s <- if (s == c_prime) 1.0 else 0.0
            I_t <- if (t == c_prime) 1.0 else 0.0

            # The core shared derivative component for class c_prime
            Q_stc <- P_st *
              (I_s + I_t - 2 * p_ik[, c_prime]) -
              2 * P_tt * (I_t - p_ik[, c_prime])

            sum_Q <- sum(Q_stc)

            # 1. Derivative for class prevalence u^rho (only c' >= 2)
            if (c_prime >= 2L) {
              col_rho <- c_prime - 1L
              J[row_J, col_rho] <- sum_Q
            }

            # 2. Derivative for item response u^phi
            for (h in seq_along(ivItemcat)) {
              K_h <- ivItemcat[h]
              n_free <- K_h - 1L
              Y_cols <- (starts_Y[h] + 1L):(starts_Y[h] + K_h - 1L)

              phi_row_start <- if (h == 1L) {
                1L
              } else {
                sum(ivItemcat[1:(h - 1L)] - 1L) + 1L
              }
              phi_vals <- phi_mat[
                phi_row_start:(phi_row_start + n_free - 1L),
                c_prime
              ]

              col_J_start <- L_rho + item_offsets[h] + (c_prime - 1L) * n_free

              for (k in seq_len(n_free)) {
                Y_col <- Y_cols[k]

                # Multiply by mDes to strictly handle missing observation partials
                val <- sum(Q_stc * Y_obs[, Y_col]) -
                  phi_vals[k] * sum(Q_stc * mDes[, Y_col])
                J[row_J, col_J_start + k] <- val
              }
            }
          }
        }
      }
      J
    }

    J_unc <- compute_J_unc_analytical(
      p.xy,
      Y.obs,
      if (is.null(mDesign)) matrix(1L, nrow(Y.obs), ncol(Y.obs)) else mDesign,
      theta1,
      ivItemcat,
      T
    )
    J.2 <- J_unc %*% J_theta1_to_u(theta1)
  }

  list(
    theta1 = theta1,
    theta2 = theta2,
    w.is = w.is,
    p.wx_mat = p.wx_mat,
    gamma_vec_to_pwx = gamma_vec_to_pwx,
    theta2_from_theta1 = theta2_from_theta1,
    J.2 = J.2
  )
}

ml_hessian_distal <- function(
  beta,
  p.zx,
  pi_mat,
  pwx,
  Z0_cc,
  w.is_cc,
  T,
  family,
  sigma2 = 1
) {
  N <- length(Z0_cc)
  pzx <- exp(pmax(p.zx(beta), -500))
  assignment_errors <- w.is_cc %*% pwx # N x T: P(W=s_i|X=t)
  q_i <- rowSums(pi_mat * pzx * assignment_errors) # N x 1
  r_it <- pi_mat * pzx * assignment_errors / q_i # N x T

  if (family == "gaussian") {
    mu <- beta
    g_it <- outer(Z0_cc, mu, "-") / sigma2 # N x T
    h_t <- rep(-1 / sigma2, T) # length T
  } else if (family == "poisson") {
    mu <- exp(beta)
    g_it <- outer(Z0_cc, rep(1, T)) -
      outer(rep(1, N), mu) # N x T: z_i - mu_t
    h_t <- -mu # length T
  } else if (family == "binomial") {
    mu <- 1 / (1 + exp(-beta))
    g_it <- outer(Z0_cc, rep(1, T)) -
      outer(rep(1, N), mu) # N x T: z_i - mu_t
    h_t <- -mu * (1 - mu) # length T
  }

  # diagonal: sum_i r_it * h_t + sum_i r_it * g_it^2 * (1 - r_it)
  diag_term <- colSums(r_it) * h_t + colSums(r_it * g_it^2 * (1 - r_it)) # length T

  # off-diagonal: -sum_i r_it * g_it * r_is * g_is
  rg <- r_it * g_it # N x T
  off_diag <- -crossprod(rg) # T x T, includes diagonal

  # combine: diagonal replaces the off_diag diagonal entries
  H_pos <- off_diag
  diag(H_pos) <- diag_term

  -H_pos # Hessian of neg.ll
}

lca_step3.distal <- function(
  neg.ll,
  beta_init,
  T,
  covariate.tol,
  use.bch = FALSE,
  Z0_cc = NULL,
  w.is_cc = NULL,
  pwx = NULL,
  em.maxIter = 200L,
  family = "gaussian",
  p.zx = NULL,
  vPi = NULL,
  pi_mat = NULL, # covariate-adjusted N x T, or NULL for flat vPi
  verbose = FALSE
) {
  N <- length(Z0_cc)
  # use covariate-adjusted pi if provided, otherwise flat vPi
  pi_s <- if (!is.null(pi_mat)) {
    pi_mat
  } else {
    matrix(vPi, ncol = T, nrow = N, byrow = TRUE)
  }

  if (use.bch) {
    D <- qr.solve(pwx)
    w.it <- w.is_cc %*% D # N x T

    # BCH score: w.it fixed, no pi_mat/pwx needed
    score_nt_bch <- function(mu) {
      if (family == "gaussian") {
        resid <- outer(Z0_cc, mu, "-")
        w.it * resid
      } else if (family == "poisson") {
        mu_val <- exp(mu)
        w.it * (outer(Z0_cc, rep(1, T)) - outer(rep(1, N), mu_val))
      } else if (family == "binomial") {
        mu_val <- 1 / (1 + exp(-mu))
        w.it * (outer(Z0_cc, rep(1, T)) - outer(rep(1, N), mu_val))
      }
    }

    w_colsums <- colSums(w.it)

    if (any(w_colsums < 0)) {
      warning(
        "BCH weights have negative column sums for at least one class. ",
        "The variance-covariance matrix will not be positive semi-definite. ",
        "SEs should be interpreted with caution. Consider use.bch = FALSE."
      )
    }

    if (family == "gaussian") {
      beta <- colSums(w.it * Z0_cc) / w_colsums # closed-form weighted mean
      resid <- outer(Z0_cc, beta, "-")
      sigma2 <- sum(w.it * resid^2) / sum(w.it)

      three_step.score <- function(params) {
        mu <- params[1:T]
        resid <- outer(Z0_cc, mu, "-")
        w.it * resid / sigma2
      }

      H.3.inv <- diag(sigma2 / w_colsums)

      res <- list(
        par = beta,
        value = neg.ll(beta),
        convergence = 0L,
        sigma2 = sigma2
      )
    } else {
      # unified NR loop for all families
      beta <- beta_init
      for (nr in seq_len(em.maxIter)) {
        grad_vec <- colSums(score_nt_bch(beta))

        if (family == "gaussian") {
          H_diag <- w_colsums
        } else if (family == "poisson") {
          H_diag <- exp(beta) * w_colsums
        } else if (family == "binomial") {
          mu_val <- 1 / (1 + exp(-beta))
          H_diag <- mu_val * (1 - mu_val) * w_colsums
        }

        direction <- grad_vec / H_diag

        step <- 1.0
        ll_cur <- -neg.ll(beta)
        for (ls in seq_len(20L)) {
          beta_new <- beta + step * direction
          ll_new <- tryCatch(-neg.ll(beta_new), error = function(e) -Inf)
          if (is.finite(ll_new) && ll_new > ll_cur) {
            break
          }
          step <- step * 0.5
        }

        delta <- step * direction
        beta <- beta + delta

        if (max(abs(delta)) < covariate.tol) {
          if (verbose) {
            message(sprintf("BCH NR converged in %d iterations.", nr))
          }
          break
        }
        if (nr == em.maxIter) warning("BCH NR reached maximum iterations.")
      }

      resid <- outer(Z0_cc, beta, "-")
      sigma2 <- sum(w.it * resid^2) / sum(w.it)

      three_step.score <- function(params) {
        mu <- params[1:T]
        if (family == "gaussian") {
          resid <- outer(Z0_cc, mu, "-")
          w.it * resid
        } else if (family == "poisson") {
          mu_val <- exp(mu)
          w.it * (outer(Z0_cc, rep(1, T)) - outer(rep(1, N), mu_val))
        } else if (family == "binomial") {
          mu_val <- 1 / (1 + exp(-mu))
          w.it * (outer(Z0_cc, rep(1, T)) - outer(rep(1, N), mu_val))
        }
      }

      H.3.inv <- tryCatch(
        diag(sigma2 / w_colsums),
        error = function(e) {
          warning("Hessian inversion failed. SEs will be NA.")
          matrix(NA_real_, T, T)
        }
      )
      res <- list(
        par = beta,
        value = neg.ll(beta),
        convergence = 0L,
        sigma2 = sigma2
      )
    }
  } else {
    if (is.null(p.zx)) {
      stop("p.zx must be provided for ML distal outcome estimation.")
    }

    # ML score helper
    score_nt_ml <- function(mu, pzx, sigma2 = 1) {
      assignment_errors <- w.is_cc %*% pwx # N x T: P(W=s_i|X=t)
      q_i <- rowSums(pi_s * pzx * assignment_errors) # N x 1

      if (family == "gaussian") {
        resid <- outer(Z0_cc, mu, "-")
        pi_s * pzx * assignment_errors * resid / (q_i * sigma2)
      } else if (family == "poisson") {
        mu_val <- exp(mu)
        score <- outer(Z0_cc, rep(1, T)) - outer(rep(1, N), mu_val)
        pi_s * pzx * assignment_errors * score / q_i
      } else if (family == "binomial") {
        mu_val <- 1 / (1 + exp(-mu))
        score <- outer(Z0_cc, rep(1, T)) - outer(rep(1, N), mu_val)
        pi_s * pzx * assignment_errors * score / q_i
      }
    }

    beta <- beta_init
    Z_long <- rep(Z0_cc, T)
    X_long <- factor(rep(seq_len(T), each = N))

    for (iter in seq_len(em.maxIter)) {
      pzx <- exp(pmax(p.zx(beta), -500))
      joint <- pi_s * pzx * (w.is_cc %*% pwx)
      w_tilde <- joint / rowSums(joint)
      w_long <- as.vector(w_tilde)

      fit <- glm(Z_long ~ X_long - 1, family = family, weights = w_long)
      beta_new <- coef(fit)

      if (abs(-neg.ll(beta_new) - (-neg.ll(beta))) < covariate.tol) {
        beta <- beta_new
        break
      }
      beta <- beta_new
      if (iter == em.maxIter) {
        warning("ML distal EM reached maximum iterations.")
      }
    }

    # estimate sigma2 after EM convergence
    sigma2 <- if (family == "gaussian") {
      pzx <- exp(pmax(p.zx(beta), -500))
      joint <- pi_s * pzx * (w.is_cc %*% pwx)
      w_tilde <- joint / rowSums(joint)
      resid <- outer(Z0_cc, beta, "-")
      sum(w_tilde * resid^2) / sum(w_tilde)
    } else {
      NULL
    }

    three_step.score <- function(params) {
      mu <- params[1:T]
      pzx <- exp(pmax(p.zx(params), -500))
      score_nt_ml(mu, pzx, sigma2 = 1)
    }

    H.3.inv <- tryCatch(
      qr.solve(ml_hessian_distal(
        beta,
        p.zx,
        pi_s,
        pwx,
        Z0_cc,
        w.is_cc,
        T,
        family
      )),
      error = function(e) {
        warning("Hessian inversion failed. SEs will be NA.")
        matrix(NA_real_, T, T)
      }
    )
    res <- list(
      par = beta,
      value = neg.ll(beta),
      convergence = 0L,
      sigma2 = sigma2
    )
  }

  return(list(
    res = res,
    H.3.inv = H.3.inv,
    three_step.score = three_step.score
  ))
}

lca_step3 <- function(
  neg.ll,
  gamma_init,
  Q,
  T,
  covariate.tol,
  use.bch = FALSE,
  gradient = NULL,
  Z_mat_cc = NULL,
  w.is_cc = NULL,
  p.xz = NULL,
  pwx = NULL,
  em.maxIter = 200L,
  verbose = FALSE,
  correct.spec = FALSE
) {
  N <- nrow(Z_mat_cc)
  beta <- matrix(gamma_init, nrow = Q, ncol = T - 1)
  ll_prev <- -neg.ll(c(beta))
  H <- NULL
  # print(ll_prev)
  if (use.bch) {
    D <- qr.solve(pwx)
    w.it <- w.is_cc %*% D # N x T, fixed BCH weights
    w.it_plus <- rowSums(w.it)

    for (nr in seq_len(em.maxIter)) {
      # print(beta)
      # print(-neg.ll(c(beta)))
      grad_vec <- -gradient(c(beta)) # gradient of pos. ll: Q*(T-1) vector

      H <- matrix(0, Q * (T - 1), Q * (T - 1))
      pi_ <- p.xz(beta)
      for (k in seq_len(T - 1)) {
        for (l in k:(T - 1)) {
          # upper triangle only
          w_kl <- w.it_plus * pi_[, k + 1L] * ((k == l) - pi_[, l + 1L])
          idx_k <- ((k - 1) * Q + 1):(k * Q)
          idx_l <- ((l - 1) * Q + 1):(l * Q)
          block <- -t(Z_mat_cc) %*% (w_kl * Z_mat_cc)
          H[idx_k, idx_l] <- block
          if (k != l) H[idx_l, idx_k] <- t(block) # Clairaut: H symmetric
        }
      }

      direction <- tryCatch(
        qr.solve(-H, grad_vec),
        error = function(e) rep(0, Q * (T - 1))
      )

      alpha <- 1
      current_ll <- -neg.ll(c(beta))

      beta_vec <- c(beta)

      while (alpha > (covariate.tol / 2)) {
        trial <- beta_vec + alpha * direction
        trial_ll <- tryCatch(-neg.ll(trial), error = function(e) -Inf)
        if (is.finite(trial_ll) && trial_ll > current_ll) {
          break
        }
        alpha <- alpha / 2
      }

      delta <- alpha * direction
      #print(max(abs(delta)))
      beta <- beta + matrix(delta, nrow = Q, ncol = T - 1)
      if (max(abs(delta)) < covariate.tol) {
        if (verbose) {
          message(sprintf("BCH NR converged in %d iterations.", nr))
        }
        break
      }
      if (nr == em.maxIter) warning("BCH NR reached maximum iterations.")
    }
    #print(H)
  } else {
    for (iter in seq_len(em.maxIter)) {
      #print(beta)
      #print(ll_prev)

      # E-step ######################################################################
      p <- p.xz(beta)

      q <- p %*% t(pwx)

      gamma <- matrix(0, nrow = N, ncol = T)
      for (t in seq_len(T)) {
        gamma[, t] <- p[, t] * rowSums(w.is_cc * outer(rep(1, N), pwx[, t]) / q)
      }
      ###############################################################################

      #M-step ########################################################################
      gamma_plus <- rowSums(gamma)
      Gamma_nr <- gamma[, -1, drop = FALSE]

      for (nr in seq_len(10L)) {
        p_nr <- p.xz(beta)
        p_nr1 <- p_nr[, -1, drop = FALSE]

        grad <- t(Z_mat_cc) %*% (Gamma_nr - p_nr1 * gamma_plus)

        H <- matrix(0, Q * (T - 1), Q * (T - 1))
        for (k in seq_len(T - 1)) {
          for (l in k:(T - 1)) {
            # upper triangle only
            w_kl <- gamma_plus * p_nr1[, k] * ((k == l) - p_nr1[, l])
            idx_k <- ((k - 1) * Q + 1):(k * Q)
            idx_l <- ((l - 1) * Q + 1):(l * Q)
            block <- -t(Z_mat_cc) %*% (w_kl * Z_mat_cc)
            H[idx_k, idx_l] <- block
            if (k != l) H[idx_l, idx_k] <- t(block) # Clairaut: H symmetric
          }
        }

        delta <- tryCatch(
          qr.solve(-H, as.vector(grad)),
          error = function(e) rep(0, Q * (T - 1))
        )
        beta <- beta + matrix(delta, nrow = Q, ncol = T - 1)
        if (max(abs(delta)) < covariate.tol) break
      }
      if (nr == 10L && max(abs(delta)) >= covariate.tol) {
        warning(sprintf(
          "M-step in EM algorithm did not converge at iteration %d",
          iter
        ))
      }

      #Alternatively, fit a multinomial logistic regression model ##########################
      # class_exp <- rep(seq_len(T), each = N)
      # Z_exp <- Z_mat_cc[rep(seq_len(N), T), , drop = FALSE]
      # w_exp <- as.vector(gamma)

      # fit <- nnet::multinom(
      #   class_exp ~ Z_exp - 1,
      #   weights = w_exp,
      #   trace = FALSE,
      #   maxit = 500L
      # )

      # beta <- t(coef(fit))
      # H <- qr.solve(-vcov(fit))
      ######################################################################################

      ll_curr <- -neg.ll(c(beta))
      if (abs(ll_curr - ll_prev) < covariate.tol && iter > 1L) {
        if (verbose) {
          message(sprintf("EM converged in %d iterations.", iter))
        }
        break
      }
      ll_prev <- ll_curr

      if (iter == em.maxIter) {
        warning("EM reached maximum iterations without converging.")
      }
    }
  }
  res <- list(par = c(beta), value = neg.ll(c(beta)), convergence = 0L)

  H.3.inv <- tryCatch(
    {
      if (!use.bch) {
        if (correct.spec) {
          matrix(NA_real_, Q * (T - 1), Q * (T - 1))
        } else {
          # -- Analytic observed-data Hessian of neg.ll --------------------------
          # neg.ll = -sum_i sum_s w_{is} * log(q_{is})
          # q_{is} = sum_t pi_{it}(beta) * pwx[s,t]
          # r_{is} = w_{is} / q_{is}
          #
          # H_{(q,k),(p,l)} = sum_i z_{iq}*z_{ip} * [
          #   pi_{i,k+1}*(I(k==l)-pi_{i,l+1}) * F_k           <- term A
          # - pi_{i,k+1} * pi_{i,l+1} * G_{kl}                 <- term B (from dF/dbeta)
          # ]
          # where:
          #   F_k  = sum_s w_{is}*pwx[s,k+1]/q_{is} - sum_s w_{is}
          #   G_{kl} = sum_s w_{is}*pwx[s,k+1]*(pwx[s,l+1]-q_{is})/q_{is}^2

          p_ <- p.xz(beta) # N x T
          q_mat <- p_ %*% t(pwx) # N x T: q[i,s]
          r_mat <- w.is_cc / q_mat # N x T: r[i,s]

          # F_k for each non-reference class k (N x (T-1))
          F_mat <- matrix(0, N, T - 1L)
          for (k in seq_len(T - 1L)) {
            F_mat[, k] <- rowSums(r_mat * pwx[, k + 1L][col(r_mat)]) -
              rowSums(w.is_cc)
          }

          # G_{kl} for each pair (k,l) -- precompute N x T(T-1)^2 is expensive;
          # compute inside the double loop instead
          H_obs <- matrix(0, Q * (T - 1L), Q * (T - 1L))
          for (k in seq_len(T - 1L)) {
            for (l in k:(T - 1L)) {
              # upper triangle only
              idx_k <- ((k - 1L) * Q + 1L):(k * Q)
              idx_l <- ((l - 1L) * Q + 1L):(l * Q)

              # term A: pi_{i,k+1}*(I(k==l)-pi_{i,l+1}) * F_k
              tA <- p_[, k + 1L] * ((k == l) - p_[, l + 1L]) * F_mat[, k]

              # term B: -pi_{i,k+1} * pi_{i,l+1} * G_{kl}
              # G_{kl}[i] = sum_s w_{is}*pwx[s,k+1]*(pwx[s,l+1]-q_{is})/q_{is}^2
              G_kl <- rowSums(
                w.is_cc *
                  pwx[, k + 1L][col(w.is_cc)] *
                  (pwx[, l + 1L][col(w.is_cc)] - q_mat) /
                  q_mat^2
              )
              tB <- -p_[, k + 1L] * p_[, l + 1L] * G_kl

              block <- -t(Z_mat_cc) %*% ((tA + tB) * Z_mat_cc) # Q x Q
              H_obs[idx_k, idx_l] <- block
              if (k != l) {
                H_obs[idx_l, idx_k] <- t(block)
              } # Clairaut: H is symmetric
            }
          }
          qr.solve(H_obs)
        }
      } else {
        qr.solve(-H)
      }
    },
    error = function(e) {
      if (!use.bch) {
        warning(
          "Exact Hessian inversion failed. Falling back to Hessian estimated from the cross-product of the case-wise log-likelihood gradients. This assumes a correct third-stage model specification."
        )
        matrix(NA_real_, Q * (T - 1), Q * (T - 1))
      } else {
        warning(
          "Hessian inversion failed after BCH. Try again with use.bch = FALSE."
        )
        matrix(NA_real_, Q * (T - 1), Q * (T - 1))
      }
    }
  )
  return(list(res = res, H.3.inv = H.3.inv))
}


# -- Variance estimation (Bakk et al., 2014) ----------------------------------
lca_step1_vcov <- function(
  fit0,
  theta1,
  boundary.tol
) {
  H.1.inv <- fit0$Varmat
  J.1 <- fit0$mScore

  is_bdry <- theta1 > (1 - boundary.tol) | theta1 < boundary.tol
  if (any(is_bdry)) {
    H.1.inv[is_bdry, ] <- 0
    H.1.inv[, is_bdry] <- 0
    J.1[, is_bdry] <- 0
  }

  Sigma.1 <- (H.1.inv %*% crossprod(J.1) %*% H.1.inv)

  Sigma.1
}

lca_vcov <- function(
  coefs,
  three_step.score,
  H.3.inv,
  fit0,
  theta1,
  theta2,
  J.2,
  p.wx_mat,
  w.is,
  Z_mat,
  n_classes,
  p.xz,
  s2,
  boundary.tol,
  use.simple.cov
) {
  J.3 <- three_step.score(c(coefs))

  Sigma.3.robust <- H.3.inv %*% crossprod(J.3) %*% H.3.inv
  if (!use.simple.cov) {
    Sigma.1 <- lca_step1_vcov(fit0, theta1, boundary.tol)

    # -- Analytic C_mat = d/d theta2 [colSums(score_3)] ----------------------------

    T_ <- n_classes
    pwx <- p.wx_mat
    pi_ <- p.xz(matrix(coefs, ncol = T_ - 1L)) # N x T
    N <- nrow(Z_mat)
    Q_ <- ncol(Z_mat)

    q <- pi_ %*% t(pwx) # N x T: q[i,s]
    V <- s2$w.is / q # N x T: w[i,s]/q[i,s]
    U <- s2$w.is / q^2 # N x T: w[i,s]/q[i,s]^2

    # VP[i, t0] = sum_s V[i,s] * pwx[s, t0]
    VP <- V %*% pwx

    n_theta2 <- T_ * (T_ - 1L)
    n_coef <- Q_ * (T_ - 1L)
    C_mat <- matrix(0, nrow = n_coef, ncol = n_theta2)

    # Map (s0, t0) mapping to the correct column index in C_mat
    theta2_idx <- matrix(0, T_, T_)
    theta2_idx[row(theta2_idx) != col(theta2_idx)] <- seq_len(n_theta2)

    for (t0 in seq_len(T_)) {
      for (s0 in seq_len(T_)) {
        if (s0 == t0) {
          next
        }
        idx_theta2 <- theta2_idx[s0, t0]

        # Store derivatives for all classes k for this specific (s0, t0)
        d_score_mat <- matrix(0, nrow = N, ncol = T_ - 1L)

        for (k in seq_len(T_ - 1L)) {
          # UP[i] = sum_s U[i,s] * pwx[s, k+1] * pwx[s, t0]
          UP_k_t0 <- rowSums(
            U *
              matrix(
                pwx[, k + 1] * pwx[, t0],
                nrow = N,
                ncol = T_,
                byrow = TRUE
              )
          )

          term1 <- 0
          if (k + 1L == t0) {
            term1 <- pi_[, k + 1] * (V[, s0] - VP[, t0])
          }

          term2 <- pi_[, k + 1] *
            pi_[, t0] *
            (U[, s0] * pwx[s0, k + 1] - UP_k_t0)

          d_score_mat[, k] <- pwx[s0, t0] * (term1 - term2)
        }

        # crossprod(Z_mat, d_score_mat) produces a Q x (T-1) matrix.
        # as.vector() unrolls it column by column, exactly matching the
        # flattened order of three_step.score(c(coefs), pwx)
        C_mat[, idx_theta2] <- as.vector(crossprod(Z_mat, d_score_mat))
      }
    }

    step1.uncertainty <- C_mat %*% J.2 %*% Sigma.1 %*% t(J.2) %*% t(C_mat)
    Sigma.3 <- H.3.inv %*% (crossprod(J.3) + step1.uncertainty) %*% H.3.inv
  } else {
    Sigma.3 <- Sigma.3.robust
    step1.uncertainty <- NULL
  }

  # Apply parameter names to vcov matrix
  param_names <- as.vector(outer(
    rownames(coefs),
    colnames(coefs),
    paste,
    sep = ":"
  ))
  rownames(Sigma.3) <- param_names
  colnames(Sigma.3) <- param_names

  Sigma.3
}

lca_vcov_distal <- function(
  mu_hat,
  three_step.score,
  pi_adj, # N x T class probabilities at converged gamma (or flat vPi)
  w.is, # N x T assignment indicator matrix
  p.wx_mat, # T x T misclassification matrix
  p.zx, # function: params -> N x T log-density matrix
  family, # "gaussian", "poisson", or "binomial"
  H.3.inv,
  fit0,
  s2,
  Sigma.3 = NULL, # covariate vcov, only needed when Zp present
  s3.par = NULL, # gamma estimates, only needed when Zp present
  p.xz.cov = NULL, # p.xz function for covariate model
  Z_mat_cov = NULL, # N x Q covariate design matrix
  T,
  boundary.tol,
  use.simple.cov,
  use.bch
) {
  J.3 <- three_step.score(mu_hat)
  meat <- crossprod(J.3)

  if (use.bch || use.simple.cov) {
    result <- H.3.inv %*% meat %*% H.3.inv
    class_names <- paste0("mu_C", seq_len(T))
    rownames(result) <- class_names
    colnames(result) <- class_names
    return(result)
  }

  Sigma.1 <- lca_step1_vcov(fit0, s2$theta1, boundary.tol)
  N <- nrow(w.is)
  pwx <- p.wx_mat # T x T

  # -- Quantities at converged mu_hat --------------------------------------------
  pzx_mat <- exp(pmax(p.zx(mu_hat), -500)) # N x T: f(z_i|mu_t)
  ae <- w.is %*% pwx # N x T: sum_s w_{is}*pwx[s,t]
  q_i <- rowSums(pi_adj * pzx_mat * ae) # N: marginal density
  r_it <- pi_adj * pzx_mat * ae / q_i # N x T: posterior r_{it}

  # g_it = unit score of log f(z_i|mu_t) wrt mu_t
  # three_step.score[i,t] = r_{it} * g_{it}  =>  g_it = score / r_it
  g_it <- J.3 / pmax(r_it, 1e-300) # N x T

  # -- C1_mat: d/d theta2 [colSums(score_distal)]  (T x T*(T-1)) ----------------
  # theta2 = off-diagonal log-ratios of pwx (column-softmax parameterzation).
  # d ae_{i,t0}/d g_{s0,t0} = pwx[s0,t0] * (w_{i,s0} - ae_{i,t0})
  # c_i = dae / ae_{i,t0}
  # t==t0: dr_{i,t0} = r_{i,t0} * c_i * (1 - r_{i,t0})
  # t!=t0: dr_{i,t}  = -r_{i,t} * r_{i,t0} * c_i
  # C1[t, col] = sum_i dr_{it} * g_{it}

  n_theta2 <- T * (T - 1L)
  C1_mat <- matrix(0, T, n_theta2)
  col_idx <- 0L

  for (t0 in seq_len(T)) {
    for (s0 in seq_len(T)) {
      if (s0 == t0) {
        next
      }
      col_idx <- col_idx + 1L
      v <- pwx[s0, t0]
      dae <- v * (w.is[, s0] - ae[, t0]) # N
      c_i <- dae / pmax(ae[, t0], 1e-300) # N
      for (t in seq_len(T)) {
        dr <- if (t == t0) {
          r_it[, t0] * c_i * (1 - r_it[, t0])
        } else {
          -r_it[, t] * r_it[, t0] * c_i
        }
        C1_mat[t, col_idx] <- sum(dr * g_it[, t])
      }
    }
  }

  step1.uncertainty <- C1_mat %*% s2$J.2 %*% Sigma.1 %*% t(s2$J.2) %*% t(C1_mat)

  # -- C_mat: d/d gamma [colSums(score_distal)]  (T x Q*(T-1)) -----------------
  # gamma enters through pi_adj = p.xz(gamma); pwx_adj treated as fixed.
  # m_{it} = pzx_{it}*ae_i[t] / q_i   (proportional to r_{it}/pi_{it})
  # A_{i,l} = pi_{i,l+1} * (m_{i,l+1} - sum_t' m_{it'}*pi_{it'})
  # d r_{it}/d gamma_{q,l} = z_{iq} * [m_{it}*pi_{it}*(I(t==l+1)-pi_{i,l+1}) - r_{it}*A_{il}]
  # C_mat[t, l*Q+q] = sum_i score_mu[it] * z_{iq} * (above)

  step2.uncertainty <- matrix(0, T, T)

  if (!is.null(s3.par) && !is.null(p.xz.cov) && !is.null(Z_mat_cov)) {
    Q_cov <- ncol(Z_mat_cov)
    C_mat <- matrix(0, T, (T - 1L) * Q_cov)

    m_it <- pzx_mat * ae / q_i # N x T
    pi_cov <- p.xz.cov(matrix(s3.par, ncol = T - 1L)) # N x T
    m_pi_sum <- rowSums(m_it * pi_cov) # N: sum_t m_{it}*pi_{it}

    for (l in seq_len(T - 1L)) {
      A_il <- pi_cov[, l + 1L] * (m_it[, l + 1L] - m_pi_sum) # N
      for (t in seq_len(T)) {
        delta_tl <- as.integer(t == l + 1L)
        inner_t <- m_it[, t] *
          pi_cov[, t] *
          (delta_tl - pi_cov[, l + 1L]) -
          r_it[, t] * A_il # N
        idx_l <- ((l - 1L) * Q_cov + 1L):(l * Q_cov)
        C_mat[t, idx_l] <- colSums(g_it[, t] * inner_t * Z_mat_cov)
      }
    }

    step2.uncertainty <- C_mat %*% Sigma.3 %*% t(C_mat)
  }

  result <- H.3.inv %*%
    (meat + step1.uncertainty + step2.uncertainty) %*%
    H.3.inv
  class_names <- paste0("mu_C", seq_len(T))
  rownames(result) <- class_names
  colnames(result) <- class_names
  result
}


#' Three-step LCA estimation with covariates and/or distal outcomes
#'
#' @param data A data.frame.
#' @param Y.names Character vector of indicator column names.
#' @param n_classes Integer. Number of latent classes.
#' @param Zp.names Character vector of covariate column names, or `NULL`.
#' @param Zo.name Single character name of the distal outcome column, or `NULL`.
#' @param step1 Pre-fitted Step-1 object (output of [tseLCA::lca_step1()]), or `NULL`.
#' @param use.two.step Logical. Use two-step starting values for Step 3.
#' @param use.modal.assignment Logical. Use modal (hard) class assignment.
#' @param include.intercept Logical. Include intercept in covariate model.
#' @param use.simple.cov Logical. Skip measurement-uncertainty correction.
#' @param incomplete Logical. FIML for missing indicators.
#' @param boundary.tol Boundary tolerance for phi parameters.
#' @param maxIter.measurement Maximum EM iterations for Step 1.
#' @param measurement.tol Convergence tolerance for Step 1.
#' @param covariate.tol Convergence tolerance for Step 3.
#' @param iter.measurement Random restarts when entropy R\eqn{^2} is low.
#' @param R2.threshold Entropy R\eqn{^2} restart threshold.
#' @param use.bch Logical. Use BCH weights instead of ML.
#' @param em.maxIter Maximum EM iterations for Step 3.
#' @param get.twostep.vcov Logical. If `TRUE`, obtain multilevLCA's
#'   bias-corrected variance-covariance matrix for the two-step gamma
#'   estimates and store it in `$two_step_vcov`. If the `fitZ` object passed
#'   via `step1` already contains a `Varmat_cor` (from a prior
#'   `fitZ_from_multiLCA` or plain `multiLCA` call), it is attached
#'   automatically even when `get.twostep.vcov = FALSE`. Default `FALSE`.
#' @param rebase Character (e.g. `"C1"`, `"C2"`) or integer specifying which
#'   latent class to use as the reference category in the multinomial logit
#'   for Steps 2 and 3. The measurement model is permuted so this class
#'   becomes column 1 before any structural estimation. Default `"C1"`.
#' @param family One of `"gaussian"` (default), `"poisson"`, `"binomial"`.
#' @param correct.spec Logical. Use model-robust (outer-product) Hessian.
#' @param verbose Logical. Print progress messages.
#'
#' @references
#' Bakk, Z., Tekle, F. B., & Vermunt, J. K. (2013). Estimating the association
#'   between latent class membership and external variables using bias-adjusted
#'   three-step approaches. *Sociological Methodology*, 43(1), 272--311.
#'   \doi{10.1177/0081175012470644}
#'
#' Bakk, Z., & Kuha, J. (2018). Two-step estimation of models between latent
#'   classes and external variables. *Psychometrika*, 83(4), 871--892.
#'   \doi{10.1007/s11336-017-9592-7}
#'
#' Bakk, Z., Pohle, M. J., & Kuha, J. (2025). Bias-adjusted three-step
#'   estimation of structural models for latent classes. *Multivariate
#'   Behavioral Research*. \doi{10.1080/00273171.2025.2473935}
#'
#' @return A list containing `$measurement_model`, `$covariate` (if
#'   `Zp.names` supplied), and/or `$distal` (if `Zo.name` supplied).
#' @examples
#' d <- generate_data(n = 200, separation = "high",
#'                    scenario = "covariate", seed = 1)
#'
#' # Measurement model only
#' fit_m <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3)
#' summary(fit_m)
#'
#'
#' # ML three-step with simple SEs (fast)
#' fit <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                   Zp.names = "Zp", use.simple.cov = TRUE)
#' summary(fit)
#' coef(fit)
#' vcov(fit)
#'
#' # Full measurement-uncertainty correction
#' fit_cor <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                       Zp.names = "Zp", use.simple.cov = FALSE,
#'                       use.modal.assignment = FALSE)
#' summary(fit_cor)
#'
#' # BCH estimator
#' fit_bch <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                       Zp.names = "Zp", use.bch = TRUE,
#'                       use.simple.cov = TRUE)
#' summary(fit_bch)
#'
#' # Change reference class
#' fit_c2 <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                      Zp.names = "Zp", use.simple.cov = TRUE,
#'                      rebase = "C2")
#' summary(fit_c2)
#'
#' # Gaussian distal outcome
#' d2 <- generate_data(200, "high", "distal", seed = 2)
#' fit_dis <- three_step(d2, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                       Zo.name = "Zo", family = "gaussian",
#'                       use.simple.cov = TRUE)
#' summary(fit_dis)
#'
#' # Pass a pre-fitted measurement model
#' fit_step1 <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3)
#' fit2 <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                    Zp.names = "Zp", step1 = fit_step1,
#'                    use.simple.cov = TRUE)
#' summary(fit2)
#'
#' # Plot item-response profiles
#' plot(fit)
#'
#' @export
three_step <- function(
  data,
  Y.names,
  n_classes,
  Zp.names = NULL,
  Zo.name = NULL,
  step1 = NULL,
  use.two.step = TRUE,
  use.modal.assignment = TRUE,
  include.intercept = TRUE,
  use.simple.cov = FALSE,
  incomplete = FALSE,
  boundary.tol = 1e-2,
  maxIter.measurement = 5000,
  measurement.tol = 1e-8,
  covariate.tol = 1e-6,
  iter.measurement = 10L,
  R2.threshold = 0.70,
  use.bch = FALSE,
  em.maxIter = 200L,
  get.twostep.vcov = FALSE,
  rebase = "C1",
  family = "gaussian",
  correct.spec = FALSE,
  verbose = FALSE
) {
  # -- Step 1: measurement model -----------------------------------------------
  ref_idx <- parse_rebase(rebase, n_classes)

  if (!is.null(step1)) {
    # Normalize: accept raw lca_step1() list or any tseLCA object
    s1 <- if (inherits(step1, "tseLCA")) step1$measurement_model else step1
    # Apply rebase permutation so the desired reference class is column 1.
    s1$fit0 <- permute_fit0_classes(s1$fit0, ref_idx)
    # Normalise fitZ names first (handles raw multiLCA output with
    # "gamma(X|C)" rownames), then rebase to the new reference class.
    if (!is.null(s1$fitZ)) {
      s1$fitZ <- normalise_fitZ_names(
        s1$fitZ,
        n_classes = n_classes
      )
      s1$fitZ <- permute_fitZ_classes(s1$fitZ, ref_idx)
    }
  } else {
    s1 <- lca_step1(
      data,
      Y.names,
      n_classes,
      Zp.names,
      maxIter.measurement,
      measurement.tol,
      covariate.tol,
      iter.measurement,
      R2.threshold,
      get.twostep.vcov,
      incomplete = incomplete,
      include.intercept = include.intercept,
      rebase = rebase,
      verbose = verbose
    )
  }
  fit0 <- s1$fit0

  # Compute two-step coefficients with EM algorithm holding step1 fixed.
  # This runs when use.two.step = TRUE and no fitZ was already computed
  # (e.g. when step1 was passed in from a measurement-only fit).
  if (use.two.step && is.null(s1$fitZ) && !is.null(Zp.names)) {
    s1$fitZ <- fitZ_from_fit0(
      fit0 = s1$fit0,
      data = data,
      Y.names = Y.names,
      Zp.names = Zp.names,
      tol = covariate.tol,
      maxIter = em.maxIter,
      incomplete = incomplete,
      include.intercept = include.intercept,
      rebase = rebase,
      verbose = verbose
    )
  }

  # Assign fitZ after the block above so it reflects any newly computed value.
  fitZ <- s1$fitZ

  # -- Data preparation --------------------------------------------------------
  cd <- clean_data(
    data = data,
    Y.names = Y.names,
    Zp.names = Zp.names,
    Zo.name = Zo.name,
    incomplete = incomplete,
    include.intercept = include.intercept,
    verbose = verbose
  )
  Y.obs <- cd$Y.obs
  mDesign <- cd$mDesign
  ivItemcat <- cd$ivItemcat
  keep_Y <- cd$keep_Y
  Z_mat <- cd$Z_mat
  keep_step3_Z_in_Y <- cd$keep_step3_Z_in_Y
  Z0_mat <- cd$Z0_mat
  keep_step3_Z0_in_Y <- cd$keep_step3_Z0_in_Y
  keep_step3_Z0 <- cd$keep_step3_Z0

  # -- Early return if no covariates -------------------------------------------
  if (is.null(Zp.names) && is.null(Zo.name)) {
    return(structure(
      list(
        measurement_model = s1,
        llik = s1$fit0$LLKSeries[nrow(s1$fit0$LLKSeries)],
        AIC = s1$fit0$AIC,
        BIC = s1$fit0$BIC,
        R2entr = s1$fit0$R2entr,
        n_classes = n_classes
      ),
      class = c("tseLCA_measurement", "tseLCA")
    ))
  }

  # -- Step 2 ------------------------------------------------------------------
  s2 <- lca_step2(
    Y.obs,
    fit0,
    n_classes,
    use.modal.assignment,
    boundary.tol,
    use.simple.cov || use.bch,
    ivItemcat = ivItemcat,
    mDesign = mDesign
  )

  Q <- if (!is.null(Z_mat)) ncol(Z_mat) else 0L
  T <- n_classes

  # Subset s2 outputs to the rows used in each Step 3 model.
  # s2 was estimated on all keep_Y rows; Step 3 models only use complete-Z rows.
  s2_for_cov <- if (!is.null(Z_mat)) {
    list(
      theta1 = s2$theta1,
      theta2 = s2$theta2,
      p.wx_mat = s2$p.wx_mat,
      gamma_vec_to_pwx = s2$gamma_vec_to_pwx,
      theta2_from_theta1 = s2$theta2_from_theta1,
      J.2 = s2$J.2,
      w.is = s2$w.is[keep_step3_Z_in_Y, , drop = FALSE],
      post = if (!is.null(s2$post)) {
        s2$post[keep_step3_Z_in_Y, , drop = FALSE]
      } else {
        NULL
      }
    )
  } else {
    NULL
  }

  s2_for_dis <- if (!is.null(Z0_mat)) {
    list(
      theta1 = s2$theta1,
      theta2 = s2$theta2,
      p.wx_mat = s2$p.wx_mat,
      gamma_vec_to_pwx = s2$gamma_vec_to_pwx,
      theta2_from_theta1 = s2$theta2_from_theta1,
      J.2 = s2$J.2,
      w.is = s2$w.is[keep_step3_Z0_in_Y, , drop = FALSE],
      post = if (!is.null(s2$post)) {
        s2$post[keep_step3_Z0_in_Y, , drop = FALSE]
      } else {
        NULL
      }
    )
  } else {
    NULL
  }

  #For covariate estimation
  if (!is.null(Z_mat)) {
    p.xz <- function(params) {
      eta_full <- cbind(0, Z_mat %*% params)
      row_max <- apply(eta_full, 1, max)
      exp_eta <- exp(eta_full - row_max)
      exp_eta / rowSums(exp_eta)
    }

    if (use.bch) {
      D <- qr.solve(s2_for_cov$p.wx_mat)
      w.it <- s2_for_cov$w.is %*% D

      .ll_bch <- function(params, pwx = NULL) {
        beta.cur <- matrix(params, ncol = T - 1)
        rowSums(w.it * log(p.xz(beta.cur)))
      }

      neg.ll <- function(params) {
        beta.cur <- matrix(params, ncol = T - 1)
        -sum(w.it * log(pmax(p.xz(beta.cur), 1e-6)))
      }

      .grad_bch <- function(params) {
        beta.cur <- matrix(params, ncol = T - 1)
        pi_ <- p.xz(beta.cur)
        resid <- w.it[, -1L, drop = FALSE] -
          pi_[, -1L, drop = FALSE] * rowSums(w.it)
        -as.vector(t(Z_mat) %*% resid)
      }

      .score_bch <- function(params) {
        beta.cur <- matrix(params, ncol = T - 1)
        pi_ <- p.xz(beta.cur)
        resid <- w.it[, -1L, drop = FALSE] -
          pi_[, -1L, drop = FALSE] * rowSums(w.it)
        resid[, rep(seq_len(T - 1L), each = Q)] *
          Z_mat[, rep(seq_len(Q), T - 1L)]
      }

      three_step.ll <- .ll_bch
      three_step.grad <- .grad_bch
      three_step.score <- .score_bch
    } else {
      .ll_ml <- function(params, pwx = s2_for_cov$p.wx_mat) {
        probs <- p.xz(matrix(params, ncol = T - 1))
        rowSums(s2_for_cov$w.is * log(probs %*% t(pwx)))
      }

      .grad_ml <- function(params, pwx = s2_for_cov$p.wx_mat) {
        beta <- matrix(params, ncol = T - 1)
        p <- p.xz(beta)
        q <- p %*% t(pwx)
        r <- s2_for_cov$w.is / q
        grad <- matrix(0, nrow = T - 1, ncol = Q)
        for (k in seq_len(T - 1L)) {
          score_i <- p[, k + 1L] *
            (r %*% pwx[, k + 1L] - rowSums(s2_for_cov$w.is))
          grad[k, ] <- t(Z_mat) %*% score_i
        }
        as.vector(t(grad))
      }

      .score_ml <- function(params, pwx = s2_for_cov$p.wx_mat) {
        beta <- matrix(params, ncol = T - 1)
        p <- p.xz(beta)
        q <- p %*% t(pwx)
        r <- s2_for_cov$w.is / q
        score_ik <- matrix(0, nrow = nrow(Z_mat), ncol = T - 1)
        for (k in seq_len(T - 1L)) {
          score_ik[, k] <- p[, k + 1L] *
            (r %*% pwx[, k + 1L] - rowSums(s2_for_cov$w.is))
        }
        score_ik[, rep(seq_len(T - 1L), each = Q)] *
          Z_mat[, rep(seq_len(Q), T - 1L)]
      }

      three_step.ll <- .ll_ml
      three_step.grad <- .grad_ml
      three_step.score <- .score_ml
      neg.ll <- function(params) -sum(three_step.ll(params))
    }

    gamma_init <- if (!is.null(fitZ$mGamma) && use.two.step) {
      c(fitZ$mGamma)
    } else {
      rep(0, Q * (T - 1))
    }

    # -- Step 3 ------------------------------------------------------------------
    s3 <- lca_step3(
      neg.ll,
      gamma_init,
      Q,
      T,
      covariate.tol,
      gradient = three_step.grad,
      use.bch = use.bch,
      Z_mat_cc = Z_mat,
      w.is_cc = s2_for_cov$w.is,
      p.xz = p.xz,
      pwx = s2_for_cov$p.wx_mat,
      em.maxIter = em.maxIter,
      verbose = verbose,
      correct.spec = correct.spec
    )
    if (
      (correct.spec && !use.bch) ||
        is.null(s3$H.3.inv) ||
        any(is.na(s3$H.3.inv))
    ) {
      s3$H.3.inv <- qr.solve(crossprod(three_step.score(s3$res$par)))
    }

    coefs <- matrix(s3$res$par, ncol = T - 1)
    ref_idx <- parse_rebase(rebase, T)
    non_ref_classes <- seq_len(T)[-ref_idx]
    colnames(coefs) <- paste0("C", non_ref_classes)
    rownames(coefs) <- c("Intercept", Zp.names)

    # -- Variance -----------------------------------------------------------------
    Sigma.3 <- lca_vcov(
      coefs = coefs,
      three_step.score = three_step.score,
      H.3.inv = s3$H.3.inv,
      fit0 = fit0,
      theta1 = s2_for_cov$theta1,
      theta2 = s2_for_cov$theta2,
      J.2 = s2_for_cov$J.2,
      p.wx_mat = s2_for_cov$p.wx_mat,
      w.is = s2_for_cov$w.is,
      Z_mat = Z_mat,
      n_classes = n_classes,
      p.xz = p.xz,
      s2 = s2_for_cov,
      boundary.tol = boundary.tol,
      use.simple.cov = use.simple.cov || use.bch
    )
    Sigma.3.covariate <- Sigma.3

    # -- Model fit ----------------------------------------------------------------
    # Compute log-likelihood on Step 3 complete-case rows (Y ? Z)
    Y_cc <- Y.obs[keep_step3_Z_in_Y, , drop = FALSE]
    mDes_cc <- if (!is.null(mDesign)) {
      mDesign[keep_step3_Z_in_Y, , drop = FALSE]
    } else {
      NULL
    }
    total.llik <- joint_log_lik(
      Y_cc,
      Z_mat,
      expand_Phi(fit0$mPhi, ivItemcat),
      coefs,
      mDes_cc
    )
    total.k <- (T * ncol(Y.obs)) + (Q * (T - 1))
    total.AIC <- -2 * total.llik + 2 * total.k
    total.BIC <- -2 * total.llik + total.k * log(nrow(Y_cc))

    # -- Optional two-step vcov from multiLCA ------------------------------------
    # get.twostep.vcov = TRUE requests multilevLCA's bias-corrected SEs.
    # We skip re-estimation if a Varmat_cor is already attached to fitZ --
    # this handles three cases:
    #   (a) fitZ_from_multiLCA output: Varmat_cor lives at fitZ$raw_fit$Varmat_cor
    #   (b) Plain multiLCA output passed directly: fitZ$Varmat_cor
    #   (c) fitZ_from_fit0 output: no Varmat_cor anywhere -> must re-estimate
    #
    # If a Varmat_cor is already present on fitZ (regardless of get.twostep.vcov),
    # we always attach it to the output -- the user shouldn't lose it just because
    # they didn't set get.twostep.vcov = TRUE.

    .extract_varmat <- function(fZ) {
      # Returns Varmat_cor from whichever location it lives in, or NULL.
      if (is.null(fZ)) {
        return(NULL)
      }
      if (!is.null(fZ$Varmat_cor)) {
        return(fZ$Varmat_cor)
      }
      if (!is.null(fZ$raw_fit$Varmat_cor)) {
        return(fZ$raw_fit$Varmat_cor)
      }
      if (!is.null(fZ$raw_fit$SEs_cor_gamma)) {
        return(diag(as.vector(fZ$raw_fit$SEs_cor_gamma)^2))
      }
      NULL
    }

    .name_varmat <- function(V, fZ) {
      if (is.null(V) || is.null(fZ$mGamma)) {
        return(V)
      }
      param_names <- as.vector(outer(
        rownames(fZ$mGamma),
        colnames(fZ$mGamma),
        paste,
        sep = ":"
      ))
      rownames(V) <- param_names
      colnames(V) <- param_names
      V
    }

    # Check for an already-present Varmat_cor before deciding whether to
    # call fitZ_from_multiLCA.
    existing_varmat <- .extract_varmat(fitZ)

    two_step_vcov <- if (!is.null(existing_varmat)) {
      # Already have it -- attach regardless of get.twostep.vcov
      .name_varmat(existing_varmat, fitZ)
    } else if (get.twostep.vcov) {
      # No existing vcov -- call fitZ_from_multiLCA
      fZ_ml <- fitZ_from_multiLCA(
        data = data,
        Y.names = Y.names,
        n_classes = n_classes,
        Zp.names = Zp.names,
        maxIter.measurement = maxIter.measurement,
        measurement.tol = measurement.tol,
        covariate.tol = covariate.tol,
        iter.measurement = iter.measurement,
        R2.threshold = R2.threshold,
        incomplete = incomplete,
        rebase = rebase,
        verbose = verbose
      )
      # Use the multiLCA result as fitZ if none existed
      if (is.null(fitZ)) {
        fitZ <- fZ_ml
        s1$fitZ <- fZ_ml
      }
      raw_varmat <- .extract_varmat(fZ_ml)
      if (!is.null(raw_varmat)) {
        .name_varmat(raw_varmat, fZ_ml)
      } else {
        warning(
          "get.twostep.vcov: neither Varmat_cor nor SEs_cor_gamma found."
        )
        NULL
      }
    } else {
      NULL
    }

    s3.covariate <- structure(
      list(
        measurement_model = s1,
        two_step = if (!is.null(fitZ)) fitZ$mGamma else NULL,
        two_step_vcov = two_step_vcov,
        three_step = coefs,
        three_step_vcov = Sigma.3,
        three_step.llik = -s3$res$value,
        neg.ll = neg.ll,
        llik = total.llik,
        AIC = total.AIC,
        BIC = total.BIC,
        n_classes = T,
        estimator = if (use.bch) "BCH" else "ML"
      ),
      class = c("tseLCA_covariate", "tseLCA")
    )
  }

  if (!is.null(Z0_mat)) {
    if (!(family %in% c("gaussian", "poisson", "binomial"))) {
      message(
        'Provided family is not one of "gaussian", "poisson", nor "binomial". Defaulting to family="gaussain".'
      )
    }

    if (!is.null(Zp.names)) {
      # Build Z_mat_dis: covariate design for the Z0-complete rows.
      # keep_step3_Z0 contains the original row indices of Z0-complete obs;
      # we re-index into the full raw Z columns from data.
      Z_mat_dis <- if (!is.null(Z_mat) && length(keep_step3_Z0) > 0L) {
        Z_full_raw <- if (include.intercept) {
          m <- cbind(1, as.matrix(data[, Zp.names, drop = FALSE]))
          colnames(m) <- c("Intercept", Zp.names)
          m
        } else {
          as.matrix(data[, Zp.names, drop = FALSE])
        }
        Z_full_raw[keep_step3_Z0, , drop = FALSE]
      } else {
        NULL
      }

      if (!is.null(Z_mat_dis)) {
        pi_adj_full <- cbind(1, p.xz(matrix(s3$res$par, ncol = T - 1)))
        # p.xz was built on the Z step3 rows; rebuild for Z0 rows
        p.xz_dis <- function(params) {
          eta_full <- cbind(0, Z_mat_dis %*% params)
          row_max <- apply(eta_full, 1L, max)
          exp_eta <- exp(eta_full - row_max)
          exp_eta / rowSums(exp_eta)
        }
        pi_adj <- p.xz_dis(matrix(s3$res$par, ncol = T - 1))
      } else {
        pi_adj <- matrix(
          fit0$vPi,
          nrow = length(keep_step3_Z0_in_Y),
          ncol = T,
          byrow = TRUE
        )
      }

      res_adj <- compute_pwx_adj(
        Y.obs[keep_step3_Z0_in_Y, , drop = FALSE],
        fit0,
        ivItemcat,
        if (!is.null(mDesign)) {
          mDesign[keep_step3_Z0_in_Y, , drop = FALSE]
        } else {
          NULL
        },
        use.modal.assignment,
        pi_adj = pi_adj
      )
    } else {
      pi_adj <- matrix(
        fit0$vPi,
        nrow = length(keep_step3_Z0_in_Y),
        ncol = T,
        byrow = TRUE
      )
      res_adj <- list(
        w.is = s2_for_dis$w.is,
        p.wx_mat = s2_for_dis$p.wx_mat
      )
    }

    # w.is for starting-value GLMs -- use distal-subset rows
    w.is_dis <- res_adj$w.is

    # Create p.zx : The function that maps latent indicators to oberved distal outcome Zo (need a choice of likelihood)

    if (family == "poisson") {
      p.zx <- function(params) {
        log_mu <- params[1:T]
        mu <- exp(log_mu) # length T
        # P(Z=z_i | X=t) = mu_t^z_i * exp(-mu_t) / z_i!
        # log P = z_i * log(mu_t) - mu_t - log(z_i!)
        z <- Z0_mat[, 1L] # N vector
        outer(z, log_mu, "*") - # N x T: z_i * log(mu_t)
          outer(rep(1, nrow(Z0_mat)), mu, "*") - # N x T: mu_t
          lgamma(z + 1L) # N x 1, recycled
      }
      starting.lm <- glm(
        Z0_mat[, 1L] ~ -1 + as.factor(max.col(w.is_dis)),
        family = poisson()
      )
      beta_init <- coef(starting.lm)
    } else if (family == "binomial") {
      # params: logit(mu_t) -- logit link, ensures mu_t in (0,1)
      p.zx <- function(params) {
        logit_mu <- params[1:T]
        mu <- 1 / (1 + exp(-logit_mu)) # length T
        # P(Z=z_i | X=t) = mu_t^z_i * (1-mu_t)^(1-z_i)
        # log P = z_i * log(mu_t) + (1-z_i) * log(1-mu_t)
        z <- Z0_mat[, 1L] # N vector
        outer(z, log(mu), "*") + # N x T
          outer(1 - z, log(1 - mu), "*") # N x T, then exp
      }
      starting.lm <- glm(
        Z0_mat[, 1L] ~ -1 + as.factor(max.col(w.is_dis)),
        family = binomial()
      )
      beta_init <- coef(starting.lm) # already on logit scale
    } else {
      #(family == "gaussian")
      p.zx <- function(params) {
        mu <- params[1:T]
        resid <- outer(Z0_mat[, 1L], mu, "-")
        -0.5 * resid^2 - 0.5 * log(2 * pi)
      }
      starting.lm <- lm(Z0_mat[, 1L] ~ -1 + as.factor(max.col(w.is_dis)))
      beta_init <- coef(starting.lm)
    }

    if (use.bch) {
      D <- qr.solve(res_adj$p.wx_mat)
      w.it <- res_adj$w.is %*% D

      neg.ll <- function(params) {
        -sum(w.it * p.zx(params))
      }
    } else {
      neg.ll <- function(params) {
        pzx <- exp(pmax(p.zx(params), -500))
        # classification error probabilities for each person's assignment
        assignment_errors <- res_adj$w.is %*% res_adj$p.wx_mat
        -sum(log(rowSums(pi_adj * pzx * assignment_errors)))
      }
    }

    s3.distal <- lca_step3.distal(
      neg.ll = neg.ll,
      em.maxIter = em.maxIter,
      pwx = res_adj$p.wx_mat,
      w.is_cc = res_adj$w.is,
      Z0_cc = Z0_mat[, 1L],
      use.bch = use.bch,
      covariate.tol = covariate.tol,
      T = T,
      beta_init = beta_init,
      family = family,
      p.zx = p.zx,
      vPi = fit0$vPi,
      pi_mat = pi_adj
    )

    #Variance-covariance
    Sigma.3.distal <- lca_vcov_distal(
      mu_hat = s3.distal$res$par,
      three_step.score = s3.distal$three_step.score,
      pi_adj = pi_adj,
      w.is = res_adj$w.is,
      p.wx_mat = res_adj$p.wx_mat,
      p.zx = p.zx,
      family = family,
      H.3.inv = s3.distal$H.3.inv,
      fit0 = fit0,
      s2 = s2_for_dis,
      Sigma.3 = if (!is.null(Zp.names)) Sigma.3 else NULL,
      s3.par = if (!is.null(Zp.names)) s3$res$par else NULL,
      p.xz.cov = if (!is.null(Zp.names) && exists("p.xz_dis")) {
        p.xz_dis
      } else if (!is.null(Zp.names)) {
        p.xz
      } else {
        NULL
      },
      Z_mat_cov = if (!is.null(Zp.names) && !is.null(Z_mat_dis)) {
        Z_mat_dis
      } else if (!is.null(Zp.names)) {
        Z_mat
      } else {
        NULL
      },
      T = T,
      boundary.tol = boundary.tol,
      use.simple.cov = use.simple.cov,
      use.bch = use.bch
    )

    distal_par <- s3.distal$res$par
    names(distal_par) <- paste0("mu_C", seq_len(T))

    s3.distal.list <- list(
      three_step = distal_par,
      three_step_vcov = Sigma.3.distal
    )
  }

  if (!is.null(Z0_mat) && is.null(Z_mat)) {
    out <- s3.distal.list
    out$family <- family
    out$n_classes <- T
    out$estimator <- if (use.bch) "BCH" else "ML"
    class(out) <- c("tseLCA_distal", "tseLCA")
    return(out)
  }
  if (!is.null(Z_mat) && is.null(Z0_mat)) {
    class(s3.covariate) <- c("tseLCA_covariate", "tseLCA")
    return(s3.covariate)
  }

  out <- list(
    covariate = s3.covariate,
    distal = s3.distal.list,
    family = family,
    n_classes = T,
    estimator = if (use.bch) "BCH" else "ML"
  )
  class(out) <- c("tseLCA_both", "tseLCA")
  return(out)
}

# -- S3 methods for tseLCA objects ----------------------------------------------
#
# Four subclasses:
#   tseLCA_measurement  - measurement model only (no Zp, no Zo)
#   tseLCA_covariate    - covariate model only   (Zp present, no Zo)
#   tseLCA_distal       - distal outcome only    (Zo present, no Zp)
#   tseLCA_both         - both covariate and distal

# -- helpers -------------------------------------------------------------------

.covariate_table <- function(x) {
  # x is a tseLCA_covariate or x$covariate for tseLCA_both
  est <- as.vector(x$three_step)
  se <- sqrt(diag(x$three_step_vcov))
  zval <- est / se
  pval <- 2 * pnorm(-abs(zval))
  nms <- rownames(x$three_step_vcov)
  data.frame(
    Estimate = est,
    Std.Error = se,
    z.value = zval,
    p.value = pval,
    row.names = nms,
    check.names = FALSE
  )
}

.distal_table <- function(x, family) {
  # x is a tseLCA_distal or x$distal for tseLCA_both
  est <- x$three_step
  se <- sqrt(diag(x$three_step_vcov))
  zval <- est / se
  pval <- 2 * pnorm(-abs(zval))
  scale_label <- switch(
    family,
    gaussian = "(mean)",
    poisson = "(log mean)",
    binomial = "(logit)",
    "(parameter)"
  )
  data.frame(
    Estimate = est,
    Std.Error = se,
    z.value = zval,
    p.value = pval,
    row.names = paste0(names(est), " ", scale_label),
    check.names = FALSE
  )
}

.format_pval <- function(p) {
  # All entries formatted to the same width so decimal points align in print():
  #   "< 0.001 ***"   (special case, 11 chars)
  #   "0.0099  ** "   (6-char number + space + 4-char marker)
  #   "0.0526  *  "
  #   "0.0800  .  "
  #   "0.2275     "
  ifelse(
    p < 0.001,
    "< 0.001 ***",
    ifelse(
      p < 0.01,
      sprintf("%.4f  ** ", p),
      ifelse(
        p < 0.05,
        sprintf("%.4f  *  ", p),
        ifelse(p < 0.10, sprintf("%.4f  .  ", p), sprintf("%.4f     ", p))
      )
    )
  )
}

.print_table <- function(df, digits = 4) {
  df_fmt <- df
  df_fmt$Estimate <- round(df$Estimate, digits)
  df_fmt$Std.Error <- round(df$Std.Error, digits)
  df_fmt$z.value <- round(df$z.value, digits)
  df_fmt$p.value <- .format_pval(df$p.value)
  print(df_fmt, quote = FALSE, right = TRUE)
  cat("---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
}


# -- print ---------------------------------------------------------------------

#' Print a tseLCA model object
#'
#' Compact one-line or table summary printed to the console.
#'
#' @param x   A `tseLCA` object returned by [tseLCA::three_step()].
#' @param digits Integer. Number of decimal places for coefficient tables.
#' @param ... Further arguments.
#' @return Invisibly returns `x`.
#' @examples
#' d    <- generate_data(100, "high", "covariate", seed = 1)
#' fit_m <- three_step(d, paste0("Y", 1:6), n_classes = 3)
#' print(fit_m)
#' @export
print.tseLCA_measurement <- function(x, ...) {
  cat("tseLCA -- measurement model\n")
  cat(sprintf(
    "  Classes: %d   Log-lik: %.4f   AIC: %.2f   BIC: %.2f\n",
    x$n_classes,
    x$llik,
    x$AIC,
    x$BIC
  ))
  if (!is.null(x$R2entr)) {
    cat(sprintf("  Entropy R\u00b2: %.4f\n", x$R2entr))
  }
  invisible(x)
}

#' @rdname print.tseLCA_measurement
#' @examples
#' \donttest{
#' d   <- generate_data(200, "high", "covariate", seed = 1)
#' fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
#'                   Zp.names = "Zp", use.simple.cov = TRUE)
#' print(fit)
#' }
#' @export
print.tseLCA_covariate <- function(x, digits = 4, ...) {
  est <- if (!is.null(x$estimator)) x$estimator else "ML"
  cat("tseLCA -- three-step covariate model\n")
  cat(sprintf(
    "  Classes: %d   Estimator: %s   Log-lik: %.4f   AIC: %.2f   BIC: %.2f\n",
    x$n_classes,
    est,
    x$llik,
    x$AIC,
    x$BIC
  ))
  cat("\nCovariate coefficients (three-step):\n")
  .print_table(.covariate_table(x), digits = digits)
  invisible(x)
}

#' @rdname print.tseLCA_measurement
#' @examples
#' \donttest{
#' d   <- generate_data(200, "high", "distal", seed = 2)
#' fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
#'                   Zo.name = "Zo", use.simple.cov = TRUE)
#' print(fit)
#' }
#' @export
print.tseLCA_distal <- function(x, digits = 4, ...) {
  fam <- if (!is.null(x$family)) x$family else "gaussian"
  est <- if (!is.null(x$estimator)) x$estimator else "ML"
  cat("tseLCA -- three-step distal outcome model\n")
  cat(sprintf(
    "  Classes: %d   Estimator: %s   Family: %s\n",
    x$n_classes,
    est,
    fam
  ))
  cat("\nDistal outcome means by class:\n")
  .print_table(.distal_table(x, fam), digits = digits)
  invisible(x)
}

#' @rdname print.tseLCA_measurement
#' @export
print.tseLCA_both <- function(x, digits = 4, ...) {
  fam <- if (!is.null(x$family)) x$family else "gaussian"
  est <- if (!is.null(x$estimator)) x$estimator else "ML"
  cat("tseLCA -- three-step model with covariate and distal outcome\n")
  cat(sprintf(
    "  Classes: %d   Estimator: %s   Family: %s\n",
    x$n_classes,
    est,
    fam
  ))
  cat(sprintf(
    "  Log-lik: %.4f   AIC: %.2f   BIC: %.2f\n",
    x$covariate$llik,
    x$covariate$AIC,
    x$covariate$BIC
  ))
  cat("\nCovariate coefficients (three-step):\n")
  .print_table(.covariate_table(x$covariate), digits = digits)
  cat("\nDistal outcome means by class:\n")
  .print_table(.distal_table(x$distal, fam), digits = digits)
  invisible(x)
}


# -- summary -------------------------------------------------------------------

#' Summarize a tseLCA model object
#'
#' Verbose summary including model fit, class prevalences, item-response
#' probabilities, and coefficient tables with standard errors and p-values.
#'
#' @param object A `tseLCA` object returned by [tseLCA::three_step()].
#' @param digits Integer. Number of decimal places for coefficient tables.
#' @param ... Further arguments (currently unused).
#' @return Invisibly returns `object`.
#' @examples
#' d    <- generate_data(100, "high", "covariate", seed = 1)
#' fit_m <- three_step(d, paste0("Y", 1:6), n_classes = 3)
#' summary(fit_m)
#' @export
summary.tseLCA_measurement <- function(object, ...) {
  cat("-- tseLCA Measurement Model --------------------------------\n")
  cat(sprintf("Latent classes : %d\n", object$n_classes))
  cat(sprintf("Log-likelihood : %.4f\n", object$llik))
  cat(sprintf("AIC            : %.4f\n", object$AIC))
  cat(sprintf("BIC            : %.4f\n", object$BIC))
  if (!is.null(object$R2entr)) {
    cat(sprintf("Entropy R\u00b2     : %.4f\n", object$R2entr))
  }
  cat("\nClass prevalences:\n")
  vPi <- object$measurement_model$fit0$vPi
  names(vPi) <- paste0("C", seq_along(vPi))
  print(round(vPi, 4))
  cat("\nItem-response probabilities (P(Y=1|class)):\n")
  print(round(object$measurement_model$fit0$mPhi, 4))
  invisible(object)
}

#' @rdname summary.tseLCA_measurement
#' @examples
#' \donttest{
#' d   <- generate_data(200, "high", "covariate", seed = 1)
#' fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
#'                   Zp.names = "Zp", use.simple.cov = TRUE)
#' summary(fit)
#' }
#' @export
summary.tseLCA_covariate <- function(object, digits = 4, ...) {
  est <- if (!is.null(object$estimator)) object$estimator else "ML"
  cat("-- tseLCA Three-Step Covariate Model -----------------------\n")
  cat(sprintf("Latent classes : %d\n", object$n_classes))
  cat(sprintf("Estimator      : %s\n", est))
  cat(sprintf("Log-likelihood : %.4f\n", object$llik))
  cat(sprintf("AIC            : %.4f\n", object$AIC))
  cat(sprintf("BIC            : %.4f\n", object$BIC))

  if (!is.null(object$two_step)) {
    ts <- object$two_step
    # Guard: ensure intercept row is labelled even for pre-fix objects
    if (nrow(ts) > 0L && (is.null(rownames(ts)) || rownames(ts)[1L] == "")) {
      rownames(ts)[1L] <- "Intercept"
    }
    cat("\nTwo-step (starting) estimates:\n")
    print(round(ts, digits))
  }

  cat("\nThree-step estimates:\n")
  .print_table(.covariate_table(object), digits = digits)

  invisible(object)
}

#' @rdname summary.tseLCA_measurement
#' @examples
#' \donttest{
#' d   <- generate_data(200, "high", "distal", seed = 2)
#' fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
#'                   Zo.name = "Zo", use.simple.cov = TRUE)
#' summary(fit)
#' }
#' @export
summary.tseLCA_distal <- function(object, digits = 4, ...) {
  fam <- if (!is.null(object$family)) object$family else "gaussian"
  est <- if (!is.null(object$estimator)) object$estimator else "ML"
  cat("-- tseLCA Three-Step Distal Outcome Model -------------------\n")
  cat(sprintf("Latent classes : %d\n", object$n_classes))
  cat(sprintf("Estimator      : %s\n", est))
  cat(sprintf("Family         : %s\n", fam))
  cat("\nDistal outcome estimates by class:\n")
  .print_table(.distal_table(object, fam), digits = digits)
  invisible(object)
}

#' @rdname summary.tseLCA_measurement
#' @export
summary.tseLCA_both <- function(object, digits = 4, ...) {
  fam <- if (!is.null(object$family)) object$family else "gaussian"
  est <- if (!is.null(object$estimator)) object$estimator else "ML"
  cat("-- tseLCA Three-Step Model: Covariate + Distal Outcome -----\n")
  cat(sprintf("Latent classes : %d\n", object$n_classes))
  cat(sprintf("Estimator      : %s\n", est))
  cat(sprintf("Family         : %s\n", fam))
  cat(sprintf("Log-likelihood : %.4f\n", object$covariate$llik))
  cat(sprintf("AIC            : %.4f\n", object$covariate$AIC))
  cat(sprintf("BIC            : %.4f\n", object$covariate$BIC))

  if (!is.null(object$covariate$two_step)) {
    ts <- object$covariate$two_step
    if (nrow(ts) > 0L && (is.null(rownames(ts)) || rownames(ts)[1L] == "")) {
      rownames(ts)[1L] <- "Intercept"
    }
    cat("\nCovariate -- two-step (starting) estimates:\n")
    print(round(ts, digits))
  }

  cat("\nCovariate -- three-step estimates:\n")
  .print_table(.covariate_table(object$covariate), digits = digits)

  cat("\nDistal outcome -- three-step estimates:\n")
  .print_table(.distal_table(object$distal, fam), digits = digits)

  invisible(object)
}


# -- coef ----------------------------------------------------------------------

#' Extract coefficients from a tseLCA model object
#'
#' @param object A `tseLCA` object returned by [tseLCA::three_step()].
#' @param which Character. For covariate and both models: `"three_step"`
#'   (default) or `"two_step"`. For both models also accepts `"covariate"`,
#'   `"distal"`, or `"both"`.
#' @param step  Character. For `tseLCA_both`: `"three_step"` (default) or
#'   `"two_step"`.
#' @param ... Further arguments (currently unused).
#' @return The coefficient matrix (covariate models), named numeric vector
#'   (distal models), or a named list of both (measurement or both models).
#' @examples
#' d    <- generate_data(100, "high", "covariate", seed = 1)
#' fit_m <- three_step(d, paste0("Y", 1:6), n_classes = 3)
#' coef(fit_m)   # returns list with $prevalences and $item_probs
#' @export
coef.tseLCA_measurement <- function(object, ...) {
  list(
    prevalences = object$measurement_model$fit0$vPi,
    item_probs = object$measurement_model$fit0$mPhi
  )
}

#' @rdname coef.tseLCA_measurement
#' @examples
#' \donttest{
#' d   <- generate_data(200, "high", "covariate", seed = 1)
#' fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
#'                   Zp.names = "Zp", use.simple.cov = TRUE)
#' coef(fit)                      # three-step estimates
#' coef(fit, which = "two_step")  # two-step starting values
#' }
#' @export
coef.tseLCA_covariate <- function(
  object,
  which = c("three_step", "two_step"),
  ...
) {
  which <- match.arg(which)
  if (which == "two_step") {
    if (is.null(object$two_step)) {
      stop("No two-step estimates available in this object.", call. = FALSE)
    }
    return(object$two_step)
  }
  object$three_step
}

#' @rdname coef.tseLCA_measurement
#' @examples
#' \donttest{
#' d   <- generate_data(200, "high", "distal", seed = 2)
#' fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
#'                   Zo.name = "Zo", use.simple.cov = TRUE)
#' coef(fit)   # named vector of class means
#' }
#' @export
coef.tseLCA_distal <- function(object, ...) {
  object$three_step
}

#' @rdname coef.tseLCA_measurement
#' @examples
#' \donttest{
#' d   <- generate_data(200, "high", "covariate", seed = 1)
#' d$Zo <- rnorm(200, mean = c(-1, 0, 1)[d$X], sd = 1)
#' fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
#'                   Zp.names = "Zp", Zo.name = "Zo",
#'                   use.simple.cov = TRUE)
#' coef(fit, which = "covariate")
#' coef(fit, which = "distal")
#' coef(fit, which = "both")
#' }
#' @export
coef.tseLCA_both <- function(
  object,
  which = c("covariate", "distal", "both"),
  step = c("three_step", "two_step"),
  ...
) {
  which <- match.arg(which)
  step <- match.arg(step)
  cov_coef <- if (step == "two_step") {
    if (is.null(object$covariate$two_step)) {
      stop("No two-step covariate estimates available.", call. = FALSE)
    }
    object$covariate$two_step
  } else {
    object$covariate$three_step
  }
  if (which == "covariate") {
    return(cov_coef)
  }
  if (which == "distal") {
    return(object$distal$three_step)
  }
  list(covariate = cov_coef, distal = object$distal$three_step)
}


# -- vcov ----------------------------------------------------------------------

#' Extract the variance-covariance matrix from a tseLCA model object
#'
#' @param object A `tseLCA` object returned by [tseLCA::three_step()].
#' @param which Character. `"three_step"` (default) or `"two_step"` for
#'   covariate models; `"covariate"`, `"distal"`, or `"both"` for both models.
#' @param step  Character. For `tseLCA_both`: `"three_step"` (default) or
#'   `"two_step"`.
#' @param ... Further arguments (currently unused).
#' @return A named square matrix. The two-step vcov is only available when
#'   `get.twostep.vcov = TRUE` was set in [tseLCA::three_step()].
#' @examples
#' d    <- generate_data(100, "high", "covariate", seed = 1)
#' fit_m <- three_step(d, paste0("Y", 1:6), n_classes = 3)
#' vcov(fit_m)   # returns NULL with a message
#' @export
vcov.tseLCA_measurement <- function(object, ...) {
  message(
    "No variance-covariance matrix available for measurement-only models."
  )
  invisible(NULL)
}

#' @rdname vcov.tseLCA_measurement
#' @examples
#' \donttest{
#' d   <- generate_data(200, "high", "covariate", seed = 1)
#' fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
#'                   Zp.names = "Zp", use.simple.cov = TRUE)
#' vcov(fit)   # Q*(T-1) x Q*(T-1) vcov matrix with named rows/cols
#' }
#' @export
vcov.tseLCA_covariate <- function(
  object,
  which = c("three_step", "two_step"),
  ...
) {
  which <- match.arg(which)
  if (which == "two_step") {
    if (is.null(object$two_step_vcov)) {
      stop(
        "No two-step vcov available. Set get.twostep.vcov = TRUE in three_step().",
        call. = FALSE
      )
    }
    return(object$two_step_vcov)
  }
  object$three_step_vcov
}

#' @rdname vcov.tseLCA_measurement
#' @examples
#' \donttest{
#' d   <- generate_data(200, "high", "distal", seed = 2)
#' fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
#'                   Zo.name = "Zo", use.simple.cov = TRUE)
#' vcov(fit)   # T x T vcov matrix with mu_C1..mu_CT row/col names
#' }
#' @export
vcov.tseLCA_distal <- function(object, ...) {
  object$three_step_vcov
}

#' @rdname vcov.tseLCA_measurement
#' @examples
#' \donttest{
#' d   <- generate_data(200, "high", "covariate", seed = 1)
#' d$Zo <- rnorm(200, mean = c(-1, 0, 1)[d$X], sd = 0.5)
#' fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
#'                   Zp.names = "Zp", Zo.name = "Zo",
#'                   use.simple.cov = TRUE)
#' vcov(fit, which = "covariate")
#' vcov(fit, which = "distal")
#' }
#' @export
vcov.tseLCA_both <- function(
  object,
  which = c("covariate", "distal", "both"),
  step = c("three_step", "two_step"),
  ...
) {
  which <- match.arg(which)
  step <- match.arg(step)
  cov_vcov <- if (step == "two_step") {
    if (is.null(object$covariate$two_step_vcov)) {
      stop(
        "No two-step vcov available. Set get.twostep.vcov = TRUE in three_step().",
        call. = FALSE
      )
    }
    object$covariate$two_step_vcov
  } else {
    object$covariate$three_step_vcov
  }
  if (which == "covariate") {
    return(cov_vcov)
  }
  if (which == "distal") {
    return(object$distal$three_step_vcov)
  }
  list(covariate = cov_vcov, distal = object$distal$three_step_vcov)
}

# -- plot methods --------------------------------------------------------------
#
# All four subclasses delegate to multilevLCA's plot.multiLCA, which draws
# item-response probability profiles across classes from the Step-1 fit.
# Extra arguments (horiz, clab, ...) are passed straight through.

.get_fit0 <- function(x) {
  # Extract the raw multiLCA fit0 object from any tseLCA subclass.
  if (inherits(x, "tseLCA_both")) {
    x$covariate$measurement_model$fit0
  } else if (inherits(x, c("tseLCA_covariate", "tseLCA_distal"))) {
    x$measurement_model$fit0
  } else {
    # tseLCA_measurement
    x$measurement_model$fit0
  }
}

#' Plot item-response probability profiles for a tseLCA model
#'
#' Delegates to `plot.multiLCA` from \pkg{multilevLCA}, which draws the
#' class-specific item-response probability profiles from the Step-1
#' measurement model.
#'
#' @param x    A `tseLCA` object returned by [tseLCA::three_step()].
#' @param horiz Logical. If `TRUE`, item labels are drawn horizontally.
#' @param clab  Optional character vector of length T giving class labels.
#' @param ...  Further arguments passed to `plot.multiLCA`.
#'
#' @return Called for its side effect (a base-graphics plot). Invisibly
#'   returns `NULL`.
#' @examples
#' d    <- generate_data(100, "high", "covariate", seed = 1)
#' fit_m <- three_step(d, paste0("Y", 1:6), n_classes = 3)
#' plot(fit_m)
#'
#' \donttest{
#' # Custom class labels
#' plot(fit_m, clab = c("Low risk", "Mixed", "High risk"))
#' }
#' @export
plot.tseLCA_measurement <- function(x, horiz = FALSE, clab = NULL, ...) {
  plot(.get_fit0(x), horiz = horiz, clab = clab, ...)
  invisible(NULL)
}

#' @rdname plot.tseLCA_measurement
#' @export
plot.tseLCA_covariate <- function(x, horiz = FALSE, clab = NULL, ...) {
  plot(.get_fit0(x), horiz = horiz, clab = clab, ...)
  invisible(NULL)
}

#' @rdname plot.tseLCA_measurement
#' @export
plot.tseLCA_distal <- function(x, horiz = FALSE, clab = NULL, ...) {
  plot(.get_fit0(x), horiz = horiz, clab = clab, ...)
  invisible(NULL)
}

#' @rdname plot.tseLCA_measurement
#' @export
plot.tseLCA_both <- function(x, horiz = FALSE, clab = NULL, ...) {
  plot(.get_fit0(x), horiz = horiz, clab = clab, ...)
  invisible(NULL)
}
