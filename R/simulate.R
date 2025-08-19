#' Generate a survival time with a piecewise treatment process
#'
#' Draw a single event time for a subject given a time grid and a 0/1 treatment
#' path, a subject-specific log multiplier H_i, a treatment log effect g, and a
#' baseline shape parameter nv (1 for Exponential, 2 for Weibull).
#'
#' @param t_input numeric vector of times in ascending order.
#' @param z_input numeric (0/1) vector with the same length as t_input.
#' @param H_i scalar, subject-specific log-effect.
#' @param g scalar, log hazard ratio for treatment.
#' @param nv positive scalar; 1 means Exponential baseline, 2 means Weibull with shape 2.
#' @return A single positive event time.
#' @export
survival_time <- function(t_input, z_input, H_i, g, nv) {
  surv_df <- data.frame(t = t_input, z = z_input)
  first_one_index <- which(surv_df$z == 1)[1]
  if (is.na(first_one_index)) {
    l <- 1; t <- Inf
  } else {
    result_df <- surv_df[first_one_index:nrow(surv_df), ]
    change_indices <- c(1, which(base::diff(result_df$z) != 0) + 1)
    result_df2 <- result_df[change_indices, ]
    t <- as.numeric(result_df2$t)
    l <- length(t)
  }
  u <- -log(stats::runif(1, 0, 1))

  if (l == 5) {
    R <- c((exp(H_i)) * t[1]^nv,
           (exp(H_i)) * (t[1]^nv + exp(g) * (t[2]^nv - t[1]^nv)),
           (exp(H_i)) * (t[1]^nv + exp(g) * (t[2]^nv - t[1]^nv) + (t[3]^nv - t[2]^nv)),
           (exp(H_i)) * (t[1]^nv + exp(g) * (t[2]^nv - t[1]^nv) + (t[3]^nv - t[2]^nv) + exp(g) * (t[4]^nv - t[3]^nv)),
           (exp(H_i)) * (t[1]^nv + exp(g) * (t[2]^nv - t[1]^nv) + (t[3]^nv - t[2]^nv) + exp(g) * (t[4]^nv - t[3]^nv) + (t[5]^nv - t[4]^nv)))
    if      (u < R[1])                      L_inver <- (u/(exp(H_i)))^(1/nv)
    else if (u >= R[1] & u < R[2])          L_inver <- ((u - exp(H_i) * t[1]^nv + exp(H_i + g) * t[1]^nv)/(exp(H_i + g)))^(1/nv)
    else if (u >= R[2] & u < R[3])          L_inver <- ((u - exp(H_i) * t[1]^nv - exp(H_i + g) * (t[2]^nv - t[1]^nv) + exp(H_i) * t[2]^nv)/(exp(H_i)))^(1/nv)
    else if (u >= R[3] & u < R[4])          L_inver <- ((u - exp(H_i) * t[1]^nv - exp(H_i + g) * (t[2]^nv - t[1]^nv) - exp(H_i) * (t[3]^nv - t[2]^nv) + exp(H_i + g) * t[3]^nv)/(exp(H_i + g)))^(1/nv)
    else if (u >= R[4] & u < R[5])          L_inver <- ((u - exp(H_i) * t[1]^nv - exp(H_i + g) * (t[2]^nv - t[1]^nv) - exp(H_i) * (t[3]^nv - t[2]^nv) - exp(H_i + g) * (t[4]^nv - t[3]^nv) + exp(H_i) * t[4]^nv)/(exp(H_i)))^(1/nv)
    else                                     L_inver <- ((u - exp(H_i) * t[1]^nv - exp(H_i + g) * (t[2]^nv - t[1]^nv) - exp(H_i) * (t[3]^nv - t[2]^nv) - exp(H_i + g) * (t[4]^nv - t[3]^nv) - exp(H_i) * (t[5]^nv - t[4]^nv) + exp(H_i + g) * t[5]^nv)/(exp(H_i + g)))^(1/nv)
  } else if (l == 4) {
    R <- c((exp(H_i)) * t[1]^nv,
           (exp(H_i)) * (t[1]^nv + exp(g) * (t[2]^nv - t[1]^nv)),
           (exp(H_i)) * (t[1]^nv + exp(g) * (t[2]^nv - t[1]^nv) + (t[3]^nv - t[2]^nv)),
           (exp(H_i)) * (t[1]^nv + exp(g) * (t[2]^nv - t[1]^nv) + (t[3]^nv - t[2]^nv) + exp(g) * (t[4]^nv - t[3]^nv)))
    if      (u < R[1])                      L_inver <- (u/(exp(H_i)))^(1/nv)
    else if (u >= R[1] & u < R[2])          L_inver <- ((u - exp(H_i) * t[1]^nv + exp(H_i + g) * t[1]^nv)/(exp(H_i + g)))^(1/nv)
    else if (u >= R[2] & u < R[3])          L_inver <- ((u - exp(H_i) * t[1]^nv - exp(H_i + g) * (t[2]^nv - t[1]^nv) + exp(H_i) * t[2]^nv)/(exp(H_i)))^(1/nv)
    else if (u >= R[3] & u < R[4])          L_inver <- ((u - exp(H_i) * t[1]^nv - exp(H_i + g) * (t[2]^nv - t[1]^nv) - exp(H_i) * (t[3]^nv - t[2]^nv) + exp(H_i + g) * t[3]^nv)/(exp(H_i + g)))^(1/nv)
    else                                     L_inver <- ((u - exp(H_i) * t[1]^nv - exp(H_i + g) * (t[2]^nv - t[1]^nv) - exp(H_i) * (t[3]^nv - t[2]^nv) - exp(H_i + g) * (t[4]^nv - t[3]^nv) + exp(H_i) * t[4]^nv)/(exp(H_i)))^(1/nv)
  } else if (l == 3) {
    R <- c((exp(H_i)) * t[1]^nv,
           (exp(H_i)) * (t[1]^nv + exp(g) * (t[2]^nv - t[1]^nv)),
           (exp(H_i)) * (t[1]^nv + exp(g) * (t[2]^nv - t[1]^nv) + (t[3]^nv - t[2]^nv)))
    if      (u < R[1])                      L_inver <- (u/(exp(H_i)))^(1/nv)
    else if (u >= R[1] & u < R[2])          L_inver <- ((u - exp(H_i) * t[1]^nv + exp(H_i + g) * t[1]^nv)/(exp(H_i + g)))^(1/nv)
    else if (u >= R[2] & u < R[3])          L_inver <- ((u - exp(H_i) * t[1]^nv - exp(H_i + g) * (t[2]^nv - t[1]^nv) + exp(H_i) * t[2]^nv)/(exp(H_i)))^(1/nv)
    else                                     L_inver <- ((u - exp(H_i) * t[1]^nv - exp(H_i + g) * (t[2]^nv - t[1]^nv) - exp(H_i) * (t[3]^nv - t[2]^nv) + exp(H_i + g) * t[3]^nv)/(exp(H_i + g)))^(1/nv)
  } else if (l == 2) {
    R <- c((exp(H_i)) * t[1]^nv,
           (exp(H_i)) * (t[1]^nv + exp(g) * (t[2]^nv - t[1]^nv)))
    if      (u < R[1])                      L_inver <- (u/(exp(H_i)))^(1/nv)
    else if (u >= R[1] & u < R[2])          L_inver <- ((u - exp(H_i) * t[1]^nv + exp(H_i + g) * t[1]^nv)/(exp(H_i + g)))^(1/nv)
    else                                     L_inver <- ((u - exp(H_i) * t[1]^nv - exp(H_i + g) * (t[2]^nv - t[1]^nv) + exp(H_i) * t[2]^nv)/(exp(H_i)))^(1/nv)
  } else {
    R <- (exp(H_i)) * t[1]^nv
    if      (u < R[1])                      L_inver <- (u/(exp(H_i)))^(1/nv)
    else                                     L_inver <- ((u - exp(H_i) * t[1]^nv + exp(H_i + g) * t[1]^nv)/(exp(H_i + g)))^(1/nv)
  }
  L_inver
}

#' Simulate a dataset with a time-varying treatment
#'
#' @param N number of subjects.
#' @param g log hazard ratio for treatment.
#' @param a0 baseline family: "Weibull" (nv = 2) or "Exponential" (nv = 1).
#' @param BARTsetting either 1 or 2, controls how covariates enter f(x).
#' @param Mi number of grid points per subject (excluding time zero).
#' @return A data.frame in long format containing the stacked subject grids.
#' @export
Simulate_data <- function(N, g, a0 = c("Weibull","Exponential"), BARTsetting = 1, Mi = 5) {
  a0 <- match.arg(a0)
  nv <- if (a0 == "Weibull") 2 else 1

  n <- N
  x1 <- sample(c(0, 1), size = n, replace = TRUE)
  x2 <- stats::rnorm(n, 0, 1)
  x3 <- stats::rnorm(n, 10, 5)
  x4 <- stats::rnorm(n, 20, 8)

  if (BARTsetting == 1) {
    b0 <- -2; b1 <- 1; b2 <- 0.25
    log_fx <- 0.1 * (b0 + b1*x1 + b2*x2^2 + 0.5*sin(x3) - 0.5*cos(x4) + 1.25)
    log_fx <- 0.05 * scale(log_fx)
  } else {
    b0 <- -3; b1 <- 2; b2 <- -0.5; b3 <- 0.75
    log_fx <- b0 + b1*x1 + b2*x2^2 + b3*x3 + sin(x2) - cos(x3) + 2*sin(x4) - 1.5*cos(x2*x3)
    log_fx <- 0.05 * scale(log_fx)
  }

  dat <- data.frame(id = 1:n, x1 = x1, x2 = x2, x3 = x3, x4 = x4, log_fx = as.numeric(log_fx))
  dat$H <- dat$log_fx
  dat$g <- g

  treat_times <- data.frame(V = rep(NA, (Mi+1)*n), j = rep(0:Mi, n), id = rep(1:n, each = (Mi+1)))
  lambda_rate <- 3
  treat_times$V <- as.vector(sapply(1:n, function(x) { c(0, cumsum(stats::rexp(Mi, rate = lambda_rate))) }))
  treat_times$Z <- as.vector(sapply(1:n, function(x) { c(0, stats::rbinom(Mi, 1, 0.5)) }))

  dat_treat <- merge(dat, treat_times, by = "id")
  dat_treat2 <- dat_treat[dat_treat$j != 0, ]

  survivaltime <- sapply(1:n, function(x) {
    survival_time(dat_treat2[dat_treat2$id == x, ]$V,
                  dat_treat2[dat_treat2$id == x, ]$Z,
                  dat_treat2[dat_treat2$id == x, ]$H[1],
                  dat_treat2[dat_treat2$id == x, ]$g[1],
                  nv = nv)
  })
  survivaltimedf <- data.frame(id = 1:n, T_i = survivaltime)
  final_data <- merge(dat_treat, survivaltimedf, by = "id")

  c_df <- data.frame(C_i = stats::runif(n, min = 0, max = max(survivaltime)), id = 1:n)
  final_data2 <- merge(final_data, c_df, by = "id")
  final_data2$ind_c <- ifelse(final_data2$T_i <= final_data2$C_i, 1, 0)
  final_data2$y_i   <- ifelse(final_data2$T_i <= final_data2$C_i, final_data2$T_i, final_data2$C_i)

  final_data2$t0    <- ifelse(final_data2$y_i >= final_data2$V, final_data2$V, final_data2$y_i)
  final_data2$logt0 <- round(log(final_data2$t0), 5)
  final_data2$a0    <- log(nv * (final_data2$t0)^(nv-1))
  final_data2$lam0  <- nv * (final_data2$t0)^(nv-1)
  final_data2$Lam0  <- (final_data2$t0)^nv
  final_data2$log_hazard <- final_data2$a0 + final_data2$log_fx + final_data2$g * final_data2$Z

  final_data2
}

event_i <- function(i, data, Mi) {
  subject_data <- data[data$id == i, ]
  event_sub <- 0
  for (j in 0:Mi) {
    if (subject_data$V[j + 1] <= subject_data$y_i[1]) event_sub <- event_sub + 1
  }
  event_sub
}

calculate_z_i <- function(i, u, data, Mi) {
  subject_data <- data[data$id == i, ]
  zi_u <- 0
  for (j in 1:Mi) {
    if (u >= subject_data$V[j] && u < subject_data$V[j + 1]) zi_u <- zi_u + subject_data$Z[j]
  }
  if (u >= subject_data$V[(Mi + 1)]) zi_u <- subject_data$Z[(Mi + 1)]
  zi_u
}
