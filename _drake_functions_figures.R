make_d_plot <- function(half_life, dose, blq, loq00, loq20, noerror) {
  result <-
    half_life %>%
    full_join(dose) %>%
    full_join(blq) %>%
    full_join(
      loq00 %>%
        rename_at(.vars=-1, .funs=function(x) paste0(x, "_loq00"))
    ) %>%
    full_join(
      loq20 %>%
        rename_at(.vars=-1, .funs=function(x) paste0(x, "_loq20"))
    ) %>%
    full_join(
      noerror %>%
        rename_at(.vars=-1, .funs=function(x) paste0(x, "_noerror"))
    ) %>%
    mutate(
      tobit_20_00_same=abs((half.life_tobit_loq20 - thalf)/thalf) < 0.00001,
      std_20_00_same=abs((half.life_std_loq20 - thalf)/thalf) < 0.00001
    )
  
  almost_zero <- c(result$half.life_tobit_loq20, result$half.life_std_loq20)
  almost_zero <- min(almost_zero[!is.na(almost_zero) & almost_zero > 0])/2
  ret <-
    result %>%
    mutate(
      half.life_tobit_loq20_noneg=
        case_when(
          half.life_tobit_loq20 < 0~0,
          TRUE~half.life_tobit_loq20
        ),
      half.life_std_loq20_noneg=
        case_when(
          half.life_std_loq20 < 0~0,
          TRUE~half.life_std_loq20
        ),
      half.life_tobit_loq20_slight_positive=
        case_when(
          half.life_tobit_loq20 < 0~almost_zero,
          TRUE~half.life_tobit_loq20
        ),
      half.life_std_loq20_slight_positive=
        case_when(
          half.life_std_loq20 < 0~almost_zero,
          TRUE~half.life_std_loq20
        )
    )
  ret
}

verify_and_count_missing <- function(d_plot) {
  d_plot %>%
    # Both or neither method is missing results on the same rows
    verify(!xor(is.na(half.life_std_loq20), is.na(half.life_std_loq20))) %>%
    # NA values are only when there are <= 2 points available
    verify(xor(n_above_loq > 2, is.na(half.life_std_loq20)))
  sum(is.na(d_plot$half.life_tobit_loq20))
}

make_count_negative <- function(d_plot) {
  tibble(
    tobit_count=sum(d_plot$half.life_tobit_loq20 < 0),
    std_count=sum(d_plot$half.life_std_loq20 < 0),
    n=nrow(d_plot),
    tobit_percent=100*tobit_count/n,
    std_percent=100*std_count/n
  )
}

make_figure1 <- function(d_plot) {
  d_plot_mod <-
    d_plot %>%
    mutate(
      n_below_loq_before_or_at_tlast_Text=
        paste0("`", n_below_loq_before_or_at_tlast, " BLQ before t`[last]"),
      n_below_loq_after_tlast_Text=
        paste0("`", n_below_loq_after_tlast, " BLQ after t`[last]")
    )
  ggplot(d_plot_mod) +
    geom_vline(xintercept=1, colour="gray", size=1) +
    geom_vline(xintercept=c(0.5, 2), colour="gray", size=1, linetype="63") +
    stat_ecdf(aes(x=half.life_tobit_loq20_slight_positive/thalf, colour="tobit"), size=1) +
    stat_ecdf(aes(x=half.life_std_loq20_slight_positive/thalf, colour="least-squares"), size=1) +
    scale_x_log10(breaks=c(0.1, 0.5, 1, 2, 10)) +
    coord_cartesian(xlim=c(0.1, 10)) +
    labs(
      x="Estimated/Theoretical Half-Life Ratio",
      y="Cumulative Distribution of Ratios",
      colour=NULL
    ) +
    theme(
      legend.position=c(0.95, 0.05),
      legend.justification=c(1, 0)
    )
}

make_figure1_tmdd <- function(p) {
  p +
    facet_grid(~TMDD_Text) +
    theme(
      legend.position="bottom",
      legend.justification=NULL
    )
}

make_figure1_early_blq <- function(p) {
  p_mod <- p
  # Insert the hline under all the other elements so that it doesn't obscure the
  # ecdf lines.
  p_mod$layers <-
    c(
      geom_hline(yintercept=0.5, colour="gray"),
      p_mod$layers
    )
  p_mod +
    facet_grid(
      n_below_loq_before_or_at_tlast_Text~.,
      labeller=label_parsed
    ) +
    theme(
      legend.position="bottom",
      legend.justification=NULL
    )
}

make_figure1_late_blq <- function(p) {
  p_mod <- p
  # Insert the hline under all the other elements so that it doesn't obscure the
  # ecdf lines.
  p_mod$layers <-
    c(
      geom_hline(yintercept=0.5, colour="gray"),
      p_mod$layers
    )
  p_mod +
    facet_grid(
      n_below_loq_after_tlast_Text~.,
      labeller=label_parsed
    ) +
    theme(
      legend.position="bottom",
      legend.justification=NULL
    )
}

make_figure2 <- function(d_plot) {
  ggplot(d_plot) +
    geom_bar(
      aes(
        x=lambda.z.n.points_std_loq20,
        y = 100*(..count..)/sum(..count..),
        #colour="least-squares",
        fill="least-squares"
      ),
      colour=NA,
      width=0.2, position=position_nudge(x=0.2)
    ) +
    geom_bar(
      aes(
        x=lambda.z.n.points_tobit_loq20,
        y = 100*(..count..)/sum(..count..),
        #colour="tobit (above LOQ)",
        fill="tobit (above LOQ)"
      ),
      colour=NA,
      width=0.2, position=position_nudge(x=0)
    ) +
    geom_bar(
      aes(
        x=lambda.z.n.points_all_loq20,
        y = 100*(..count..)/sum(..count..),
        #colour="tobit (total)",
        fill="tobit (total)"
      ),
      colour=NA,
      width=0.2, position=position_nudge(x=-0.2)
    ) +
    scale_x_continuous(
      breaks=
        unique(c(
          d_plot$lambda.z.n.points_std_loq20,
          d_plot$lambda.z.n.points_all_loq20,
          d_plot$lambda.z.n.points_tobit_loq20
        ))
    ) +
    labs(
      #colour=NULL,
      fill=NULL,
      x="Number of Points",
      y="Percent of Profiles"
    ) +
    theme(legend.position="bottom")
}

make_figure2_tmdd <- function(p) {
  p +
    facet_wrap(~TMDD_Text, scales="free_y", nrow=1)
}