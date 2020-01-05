fit_half_life_all_df <- function(data, conc_col, lloq_col, id) {
  cat("Calculating for ", id, "\n", sep="")
  if (is.numeric(lloq_col)) {
    lloq_col_name <- paste0(max(names(data)), "X")
    if (lloq_col == 0) {
      lloq_col <- min(data[[conc_col]][data[[conc_col]] > 0])
    }
    data[[lloq_col_name]] <- lloq_col
    lloq_col <- lloq_col_name
  }
  fit_half_life_all(
    conc=data[[conc_col]],
    time=data$time,
    lloq=data[[lloq_col]],
    extravascular=all(data$Route_Text == "extravascular")
  )
}

fit_half_life_all <- function(conc, time, lloq, extravascular=FALSE) {
  tlast <- pk.calc.tlast(conc, time, check=FALSE)
  data_fit <-
    tibble(
      log_conc=log(conc),
      time=time,
      log_lloq=log(lloq),
      mask_blq=conc < lloq
    )
  if (extravascular) {
    tmax <- pk.calc.tmax(conc, time, check=FALSE)
    data_fit <-
      data_fit %>%
      filter(time > tmax)
  }
  data_fit_no_blq <-
    data_fit %>%
    filter(!mask_blq)
  ret_std <- list()
  idx_tlast_no_blq <- which(data_fit_no_blq$time == tlast)
  for (idx_first in seq_len(max(idx_tlast_no_blq-2, 0))) {
    tmp_fit_std <-
      PKNCA:::fit_half_life(
        data_fit_no_blq[idx_first:idx_tlast_no_blq, ],
        tlast=tlast
      ) %>%
      rename(
        lambda.z_std=lambda.z,
        clast.pred_std=clast.pred,
        half.life_std=half.life,
        span.ratio_std=span.ratio
      )
    ret_std <- append(ret_std, list(tmp_fit_std))
  }
  ret_tobit <- list()
  # Require at least 2 above LOQ points for a non-trivial solution
  idx_tlast <- which(data_fit$time == data_fit_no_blq$time[nrow(data_fit_no_blq)-1])
  for (idx_first in seq_len(max(idx_tlast, 0))) {
    tmp_fit_tobit <-
      fit_half_life_tobit(
        data_fit[idx_first:nrow(data_fit),],
        tlast
      ) %>%
      rename(
        lambda.z_tobit=lambda.z,
        clast.pred_tobit=clast.pred,
        half.life_tobit=half.life,
        span.ratio_tobit=span.ratio
      )
    ret_tobit <- append(ret_tobit, list(tmp_fit_tobit))
  }
  df_std <- bind_rows(ret_std)
  df_tobit <- bind_rows(ret_tobit)
  if (nrow(df_std) & nrow(df_tobit)) {
    full_join(df_std, df_tobit)
  } else if (nrow(df_std)) {
    df_std
  } else {
    df_tobit
  }
}

#' Perform the half-life fit given the data.  The function simply fits 
#' the data without any validation.  No selection of points or any other
#' components are done.
#' 
#' @param data The data to fit.  Must have columns named "log_conc", "time",
#'   "log_lloq", and "mask_blq"
#' @param tlast The time of last observed concentration above the limit
#'   of quantification.
#' @return A data.frame with one row and columns named "tobit_residual", 
#'   "adj_tobit_residual", "lambda.z", "clast.pred", 
#'   "lambda.z.n.points", "half.life", "span.ratio"
#' @seealso \code{\link{pk.calc.half.life}}
#' @importFrom stats optim
fit_half_life_tobit <- function(data, tlast) {
  fit <-
    stats::optim(
      par=
        c(
          log_c0=max(data$log_conc),
          lambda_z=
            -log(2)*
            diff(range(data$log_conc[!data$mask_blq]))/
            diff(range(data$time)),
          log_resid_error=log(sd(data$log_conc[!data$mask_blq])/5)
        ),
      fn=fit_half_life_tobit_LL,
      log_conc=data$log_conc,
      time=data$time,
      mask_blq=data$mask_blq,
      log_lloq=data$log_lloq
    )
  tobit_residual <- exp(fit$par[["log_resid_error"]])
  n_points_before_tlast <- sum(data$time <= tlast)
  ret <-
    data.frame(
      tobit_residual=tobit_residual,
      adj_tobit_residual=tobit_residual*(n_points_before_tlast - 2)/(n_points_before_tlast - 1),
      lambda.z=fit$par[["lambda_z"]],
      clast.pred=exp(sum(fit$par[1:2]*c(1, -tlast))),
      lambda.z.time.first=min(data$time, na.rm=TRUE),
      lambda.z.n.points_all=nrow(data),
      lambda.z.n.points=sum(!data$mask_blq)
    )
  ret$half.life <- log(2)/ret$lambda.z
  ret$span.ratio <- diff(range(data$time[!data$mask_blq]))/ret$half.life
  ret
}

#' Sum of log-likelihood helper function for fit_half_life_tobit
#'
#' @param par A 3-length numeric vector of c(log(C0), lambda_z,
#'   log(resid_error))
#' @param log_conc The natural logarithm of the concentration (may include -Inf
#'   for BLQ)
#' @param time Time
#' @param mask_blq A logical vector indicating if \code{log_conc} is below the
#'   limit of quantification
#' @param log_lloq The natural logarithm of the lower limit of quantification.
#' @return -sum(likelihood)
#' @seealso \code{\link{fit_half_life_tobit}}
fit_half_life_tobit_LL <- function(par, log_conc, time, mask_blq, log_lloq) {
  log_c0 <- par[[1]]
  lambda_z <- par[[2]]
  resid_error <- exp(par[[3]])
  est <- -lambda_z * time + log_c0
  ret <- rep(NA_real_, length(time))
  if (any(mask_blq)) {
    ret[mask_blq] <-
      pnorm(
        q=log_lloq[mask_blq],
        mean=est[mask_blq],
        sd=resid_error,
        log.p=TRUE
      )
  }
  ret[!mask_blq] <-
    dnorm(
      x=log_conc[!mask_blq],
      mean=est[!mask_blq],
      sd=resid_error,
      log=TRUE
    )
  -sum(ret)
}

make_d_nca_select <- function(data, mantissa) {
  data %>%
    mutate(
      best_std_set=lambda.z.n.points >= 3 & !is.na(adj.r.squared),
      best_tobit_set=lambda.z.n.points >= 3 & !is.na(adj_tobit_residual),
      mantissa=mantissa,
      tobit_residual_mantissa=tobit_residual*mantissa^(-lambda.z.n.points_all)
    ) %>%
    group_by(ID) %>%
    mutate(
      best_std=
        case_when(
          best_std_set & adj.r.squared == max(adj.r.squared[best_std_set], -Inf, na.rm=TRUE)~TRUE,
          TRUE~FALSE
        ),
      best_tobit=
        case_when(
          best_tobit_set & tobit_residual_mantissa == min(tobit_residual_mantissa[best_tobit_set], Inf, na.rm=TRUE)~TRUE,
          TRUE~FALSE
        ),
      best_tobit_adjusted=
        case_when(
          best_tobit_set & adj_tobit_residual == min(adj_tobit_residual[best_tobit_set], Inf, na.rm=TRUE)~TRUE,
          TRUE~FALSE
        )
    )
}

make_d_nca_select_summary <- function(data) {
  data %>%
    summarize(
      r.squared=c(r.squared[best_std], NA)[1],
      adj.r.squared=c(adj.r.squared[best_std], NA)[1],
      lambda.z_std=c(lambda.z_std[best_std], NA)[1],
      clast.pred_std=c(clast.pred_std[best_std], NA)[1],
      lambda.z.time.first_std=c(lambda.z.time.first[best_std], NA)[1],
      lambda.z.n.points_std=c(lambda.z.n.points[best_std], NA)[1],
      half.life_std=c(half.life_std[best_std], NA)[1],
      span.ratio_std=c(span.ratio_std[best_std], NA)[1],
      lambda.z.time.first_tobit=c(lambda.z.time.first[best_tobit], NA)[1],
      lambda.z.n.points_tobit=c(lambda.z.n.points[best_tobit], NA)[1],
      tobit_residual=c(tobit_residual[best_tobit], NA)[1],
      adj_tobit_residual=c(adj_tobit_residual[best_tobit], NA)[1],
      lambda.z_tobit=c(lambda.z_tobit[best_tobit], NA)[1],
      clast.pred_tobit=c(clast.pred_tobit[best_tobit], NA)[1],
      lambda.z.n.points_all=c(lambda.z.n.points_all[best_tobit], NA)[1],
      half.life_tobit=c(half.life_tobit[best_tobit], NA)[1],
      span.ratio_tobit=c(span.ratio_tobit[best_tobit], NA)[1]
    )
}

optimize_mantissa <- function(test, reference) {
  opt_fun <- function(mantissa) {
    test_selected <-
      left_join(
        make_d_nca_select_summary(
          make_d_nca_select(test, mantissa=mantissa)
        ),
        reference,
        by="ID"
      )
    test_hl <- test_selected$half.life_tobit
    test_hl[test_hl < 0] <- NA_real_
    ret <- sum(abs(log(test_hl/test_selected$thalf)), na.rm=TRUE)
    cat(mantissa, ": ", ret, "\n")
    ret
  }
  optim(par=1.5, fn=opt_fun, method="Brent", lower=0, upper=5, control=list(reltol=0.01))
}

make_blq_summary <- function(data) {
  ret <-
    tibble(
      tlast=pk.calc.tlast(conc=data$CP_LOQ20pct, time=data$time, check=FALSE),
      tmax=pk.calc.tmax(conc=data$CP_LOQ20pct, time=data$time, check=FALSE)
    )
  d_calc <-
    if (unique(data$Route_Text) == "intravascular bolus") {
      data
    } else {
      data[data$time > ret$tmax,]
    }
  ret <-
    ret %>%
    mutate(
      n_above_loq=sum(d_calc$CP_LOQ20pct > 0),
      n_below_loq=sum(d_calc$CP_LOQ20pct <= 0),
      n_below_loq_before_or_at_tlast=sum(d_calc$CP_LOQ20pct[d_calc$time <= ret$tlast] <= 0),
      n_below_loq_after_tlast=sum(d_calc$CP_LOQ20pct[d_calc$time > ret$tlast] <= 0)
    )
  ret
}
