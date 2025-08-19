#' Internal helpers
#' 
#' Utility functions for spline basis construction, simple trapezoid integrals,
#' inverse-gamma sampling, and a Metropolis-Hastings update for the baseline hazard.
#' The documentation text here avoids LaTeX commands and backslashes on purpose.
#' @keywords internal
#' @noRd
#' @importFrom splines bs
#' @importFrom stats qr rgamma rnorm runif
NULL

rinvgamma <- function(n, shape, rate) 1 / stats::rgamma(n, shape = shape, rate = rate)

make_diff_matrix <- function(K, order = 2) {
  if (order == 1) return(base::diff(diag(K), differences = 1))
  if (order == 2) return(base::diff(diag(K), differences = 2))
  stop("order must be 1 or 2")
}

setup_bs_logtime <- function(t, n_basis = 30, degree = 3) {
  tlog <- log(pmax(t, 1e-10))
  B0   <- splines::bs(tlog, df = n_basis, degree = degree, intercept = TRUE)
  list(
    knots    = attr(B0, "knots"),
    boundary = attr(B0, "Boundary.knots"),
    degree   = degree,
    n_basis  = ncol(B0)
  )
}

eval_bs_logtime <- function(t, bs_def) {
  tlog <- log(pmax(t, 1e-10))
  splines::bs(tlog,
              knots = bs_def$knots,
              Boundary.knots = bs_def$boundary,
              degree = bs_def$degree,
              intercept = TRUE)
}

# Compute o_{i,k} as the integral over time of B_k(log u)*exp(g * Z_i(u))*exp(f(x_i))
# evaluated on subject-specific grids using mid-point rule.
compute_Oik <- function(surv_df, bs_def, eg, exp_fx_by_id) {
  n  <- length(unique(surv_df$id))
  K  <- bs_def$n_basis
  O  <- matrix(0, nrow = n, ncol = K)
  ids <- sort(unique(surv_df$id))

  for (ii in seq_along(ids)) {
    s <- surv_df[surv_df$id == ids[ii], ]
    s <- s[order(s$t0), ]

    if (nrow(s) < 2) next
    dt   <- base::diff(s$t0)
    idx  <- which(dt > 0)
    if (!length(idx)) next

    mid  <- s$t0[idx] + dt[idx]/2
    Bmid <- eval_bs_logtime(mid, bs_def)
    wt   <- (eg^s$Z[idx]) * dt[idx]

    O[ii, ] <- exp_fx_by_id[ii] * colSums(Bmid * wt)
  }
  O
}

logpost_phi <- function(phi, B_event, Oik, delta, Q, sigma2) {
  exp_phi <- exp(phi)
  haz_ev  <- as.numeric(B_event %*% exp_phi)
  term1   <- sum(delta * log(pmax(haz_ev, 1e-12)))
  term2   <- sum(as.numeric(Oik %*% exp_phi))
  pen     <- 0.5 / sigma2 * drop(t(phi) %*% (Q %*% phi))
  term1 - term2 - pen
}

mh_update_phi_sigma <- function(
  phi, sigma2, B_event, Oik, delta, Q,
  prop_sd = 0.15, n_steps = 1, kappa1 = 1e-3, kappa2 = 1e-3
) {
  acc <- 0L
  lp  <- logpost_phi(phi, B_event, Oik, delta, Q, sigma2)

  for (s in seq_len(n_steps)) {
    prop <- phi + stats::rnorm(length(phi), sd = prop_sd)
    lp_p <- logpost_phi(prop, B_event, Oik, delta, Q, sigma2)
    if (log(stats::runif(1)) < (lp_p - lp)) { phi <- prop; lp <- lp_p; acc <- acc + 1L }
  }

  rQ      <- stats::qr(Q)$rank
  shape   <- kappa1 + 0.5 * rQ
  rate    <- kappa2 + 0.5 * drop(t(phi) %*% (Q %*% phi))
  sigma2  <- rinvgamma(1, shape = shape, rate = rate)

  list(phi = phi, sigma2 = sigma2, acc_rate = acc / n_steps)
}

cumtrapz_by_id <- function(t, f, id) {
  out <- numeric(length(t))
  split_idx <- split(seq_along(t), id)
  for (ix in split_idx) {
    tg <- t[ix]; fg <- f[ix]
    if (length(ix) > 1) {
      dt  <- base::diff(tg)
      avg <- (fg[-length(fg)] + fg[-1]) / 2
      integ <- c(0, cumsum(dt * avg))
      out[ix] <- integ
    } else {
      out[ix] <- 0
    }
  }
  out
}
