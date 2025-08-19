#' Main simulation and Gibbs sampler
#'
#' This function simulates the dataset and runs a Gibbs sampler that alternates
#' between a log-linear BART update for subject effects and a spline-based
#' baseline hazard updated through a Metropolis-Hastings step.
#'
#' @param N number of subjects.
#' @param g log hazard ratio for treatment.
#' @param a0 baseline family: "Weibull" or "Exponential".
#' @param BARTsetting choice of f(x) generator, 1 or 2.
#' @param Mi number of time grid steps per subject (excluding time zero).
#' @param iterations total number of Gibbs iterations.
#' @return A list with elements `chain` (samples of exp(g)) and `survival_data`.
#' @export
#' @importFrom R2BayesX bayesx sx
Simulation_Gibbs <- function(N, g, a0 = c("Weibull","Exponential"),
                             BARTsetting = 1, Mi = 5, iterations = 10000) {
  a0 <- match.arg(a0)
  nv <- if (a0 == "Weibull") 2 else 1
  n  <- N

  survival_data <- Simulate_data(N, g, a0 = a0, BARTsetting = BARTsetting, Mi = Mi)

  events <- data.frame(id = 1:n,
                       event = sapply(1:n, function(i) event_i(i, survival_data, Mi)))
  survival_data <- merge(survival_data, events, by = "id")

  for (i in 1:n) {
    if (survival_data[survival_data$id == i, ]$event[1] == Mi + 1) {
      survival_data[survival_data$id == i, ]$y_i <-
        rep(survival_data[survival_data$id == i, ]$t0[(Mi + 1)], (Mi + 1))
    }
  }

  compute_diff <- function(x) c(base::diff(x), NA)
  survival_data$Lam0_diff <- ave(survival_data$Lam0, survival_data$id, FUN = compute_diff)
  survival_data$Lam0_diff[is.na(survival_data$Lam0_diff)] <- 0

  results <- data.frame(
    id   = 1:n,
    z_yi = sapply(1:n, function(i) calculate_z_i(i, survival_data[survival_data$id == i, ]$y_i[1],
                                                 survival_data, Mi))
  )
  survival_data <- merge(survival_data, results, by = "id")

  survival_data <- do.call(rbind, lapply(1:n, function(x) {
    subject_data <- survival_data[survival_data$id == x, ]
    if (subject_data$event[1] != (Mi + 1)) {
      subject_data <- subject_data[2:(subject_data$event[1] + 1), ]
    } else { subject_data <- subject_data[2:(Mi + 1), ] }
    subject_data
  }))
  survival_data$log_fx0 <- survival_data$log_fx

  initial_g <- g
  chain_eg <- array(dim = c(iterations + 1, 3))
  colnames(chain_eg) <- c("eg", "A", "C")
  chain_eg[1, ] <- c(exp(initial_g), 0, 0)
  survival_data$eg0 <- chain_eg[1, 1]
  chain_eg <- as.data.frame(chain_eg)

  A <- sum(sapply(1:n, function(i) (survival_data[survival_data$id == i, ]$ind_c[1] *
                                      survival_data[survival_data$id == i, ]$z_yi[1])))

  run_ind <- 0
  for (i in 1:iterations) {
    if (run_ind == 0) {
      survival_data <- survival_data[, !colnames(survival_data) %in% "B_ij"]

      wi <- sapply(1:n, function(i) {
        subject_data <- survival_data[survival_data$id == i, ]
        exp(subject_data$ind_c[1] * (subject_data$a0[nrow(subject_data)] +
                                       log(subject_data$eg0[1]) * subject_data$z_yi[1]))
      })
      vi <- sapply(1:n, function(i) {
        subject_data <- survival_data[survival_data$id == i, ]
        as.numeric(t(exp(log(subject_data$eg0) * subject_data$Z)) %*% subject_data$Lam0_diff)
      })

      m <- 5
      aa0 <- sqrt(m) * sd(survival_data$log_fx0)
      aaa <- m / (aa0^2) + 0.5
      bbb <- m / (aa0^2)

      survival_data_fit <- unique(survival_data[, colnames(survival_data) %in%
                                                  c("id","x1","x2","x3","x4","log_fx0","ind_c")])
      X <- as.matrix(survival_data_fit[, c("x1","x2","x3","x4")])

      bart_fit <- loglinearBART::loglinearbart(
        X, survival_data_fit$log_fx0, survival_data_fit$ind_c, wi, vi,
        usequants = FALSE, ntree = 200, ndpost = 2, nskip = 0,
        printevery = 100, numcut = 100, aa = aaa, bb = bbb
      )

      if (sum(is.na(bart_fit$yhat.train)) > 0) {
        run_ind <- 1
        C <- mean(tail(chain_eg$C[!is.na(chain_eg$C)], 10))
      } else {
        log_fx0_df <- data.frame(id = 1:n, log_fx0 = apply(bart_fit$yhat.train, 2, mean))
        survival_data <- survival_data[, !colnames(survival_data) %in% c("log_fx0")]
        survival_data <- merge(survival_data, log_fx0_df, by = "id")

        survival_data <- survival_data[, !colnames(survival_data) %in%
                                         c("a0","lam0","Lam0","offsets","Lam0_diff")]
        survival_data$offsets <- survival_data$log_fx0 + log(survival_data$eg0) * survival_data$Z

        b1 <- R2BayesX::bayesx(
          log_hazard ~ R2BayesX::sx(logt0, knots = 30, degree = 3) + R2BayesX::sx(offsets, bs = "offset"),
          data = survival_data
        )
        survival_data$a0   <- predict(b1, newdata = survival_data) - survival_data$offsets
        survival_data$lam0 <- exp(survival_data$a0)

        K_basis <- 30
        bs_def  <- setup_bs_logtime(survival_data$t0, n_basis = K_basis, degree = 3)
        Dmat    <- make_diff_matrix(K_basis, order = 2)
        Qpen    <- t(Dmat) %*% Dmat

        if (!exists("phi_curr", inherits = FALSE))   phi_curr  <- rep(0, K_basis)
        if (!exists("sigma2_curr", inherits = FALSE)) sigma2_curr <- 1
        if (!exists("mh_sd", inherits = FALSE))       mh_sd     <- 0.15

        B_all   <- eval_bs_logtime(survival_data$t0, bs_def)

        ids     <- sort(unique(survival_data$id))
        n_subj  <- length(ids)
        B_event <- matrix(0, nrow = n_subj, ncol = K_basis)
        delta   <- numeric(n_subj)
        yvec    <- numeric(n_subj)
        exp_fx_by_id <- numeric(n_subj)

        for (ii in seq_along(ids)) {
          si <- survival_data[survival_data$id == ids[ii], ]
          yvec[ii]         <- si$y_i[1]
          B_event[ii, ]    <- eval_bs_logtime(si$y_i[1], bs_def)
          delta[ii]        <- si$ind_c[1]
          exp_fx_by_id[ii] <- exp(si$log_fx0[1])
        }

        eg_now <- survival_data$eg0[1]
        Oik    <- compute_Oik(survival_data, bs_def, eg = eg_now, exp_fx_by_id = exp_fx_by_id)

        mh_res <- mh_update_phi_sigma(
          phi = phi_curr, sigma2 = sigma2_curr, B_event  = B_event,
          Oik = Oik, delta = delta, Q = Qpen, prop_sd = mh_sd, n_steps = 1,
          kappa1 = 1e-3, kappa2 = 1e-3
        )
        phi_curr    <- mh_res$phi
        sigma2_curr <- mh_res$sigma2
        if (mh_res$acc_rate > 0.4) mh_sd <- mh_sd * 1.1
        if (mh_res$acc_rate < 0.2) mh_sd <- mh_sd * 0.9

        lam0_all <- as.numeric(B_all %*% exp(phi_curr))
        a0_all   <- log(pmax(lam0_all, 1e-12))
        survival_data$a0   <- a0_all
        survival_data$lam0 <- lam0_all

        survival_data$Lam0 <- cumtrapz_by_id(survival_data$t0, survival_data$lam0, survival_data$id)
        survival_data$Lam0_diff <- ave(survival_data$Lam0, survival_data$id,
                                       FUN = function(x) { z <- c(base::diff(x), 0); z[is.na(z)] <- 0; z })
        survival_data$B_ij <- survival_data$Lam0_diff * exp(survival_data$log_fx0)

        C <- sum(sapply(1:n, function(i) {
          subject_data <- survival_data[survival_data$id == i, ]
          as.numeric(t(subject_data$B_ij) %*% subject_data$Z)
        }))
      }
      chain_eg[i + 1, ] <- stats::rgamma(1, shape = A + exp(g) * (n/5), rate = C + 1 * (n/5))
      chain_eg[i + 1, "A"] <- A
      chain_eg[i + 1, "C"] <- C
      survival_data$eg0 <- chain_eg[i + 1, 1]
    } else {
      chain_eg[i + 1, ] <- stats::rgamma(1, shape = A + exp(g) * (n/5), rate = C + 1 * (n/5))
      chain_eg[i + 1, "A"] <- A
      chain_eg[i + 1, "C"] <- C
    }
    if (i %% 50 == 0) message("Iteration: ", i)
  }
  list(chain = chain_eg, survival_data = survival_data)
}
