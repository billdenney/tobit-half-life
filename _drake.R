library(microbenchmark)
library(assertr)
library(pmxTools)
library(PKNCA)
library(mrgsolve)
library(drake)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

source("_drake_functions_scenarios.R")
source("_drake_functions_nca.R")
source("_drake_functions_figures.R")

model <-
  mread_cache(model="generic_tmdd.Rmd", soloc=".mrgsolve/") %>%
  update(atol=1e-50)

plan_scenarios <-
  drake_plan(
    n_subjects_per_scenario=2,
    ratios=2^c(-2, 0, 2),
    d_scenario_ncompartments=make_d_scenario_ncompartments(ratios),
    d_scenarios_tmdd=make_d_scenarios_tmdd(ratios),
    d_scenarios_route=make_d_scenarios_route(ratios),
    d_scenarios_residual=make_d_scenarios_residual(),
    d_scenarios=
      make_d_scenarios(crossing(
        d_scenario_ncompartments,
        d_scenarios_tmdd,
        d_scenarios_route,
        d_scenarios_residual
      )),
    d_obs_prep=make_d_obs_prep(),
    d_dose=make_d_dose(d_scenarios, n=n_subjects_per_scenario),
    true_half_life=
      with(
        d_dose,
        calc_derived_3cpt(
          CL=CL, V1=Vcentral,
          V2=Vperiph1, Q2=Qcp1,
          V3=Vperiph2, Q3=Qcp2,
          ka=Ka,
          sigdig=Inf
        )
      ),
    d_true_half_life=
      tibble(
        ID=d_dose$ID,
        thalf=
          case_when(
            d_dose$Compartment_Text %in% "One compartment"~true_half_life$thalf_alpha,
            d_dose$Compartment_Text %in% "Two compartment"~true_half_life$thalf_beta,
            d_dose$Compartment_Text %in% "Three compartment"~true_half_life$thalf_gamma
          )
      ) %>%
      verify(!is.na(thalf)),
    d_obs=crossing(d_obs_prep, tibble(ID=unique(d_dose$ID))),
    d_sim_data=
      bind_rows(d_obs, d_dose[, names(d_obs)]) %>%
      arrange(ID, time, -evid),
    o_simresult=mrgsim(model, data=d_sim_data, idata=d_dose),
    d_simresult=make_d_simresult(o_simresult, d_dose, model)
  )

plan_nca <-
  drake_plan(
    d_nca_prep=
      d_simresult %>%
      group_by(ID) %>%
      nest(),
    d_blq_summary=
      d_nca_prep %>%
      mutate(
        blq_summary=
          pmap(
            .l=list(data=data),
            .f=make_blq_summary
          )
      ) %>%
      select(-data) %>%
      unnest(cols="blq_summary"),
    d_nca_initial_calc_loq20=
      d_nca_prep %>%
      mutate(
        LOQ20pct=
          pmap(
            .l=list(data=data, id=ID),
            .f=fit_half_life_all_df,
            conc_col="CP_LOQ20pct",
            lloq_col="LOQ20pct"
          )
      ) %>%
      select(-data),
    d_nca_initial_calc_loq00=
      d_nca_prep %>%
      mutate(
        LOQ00pct=
          pmap(
            .l=list(data=data, id=ID),
            .f=fit_half_life_all_df,
            conc_col="CP",
            lloq_col=0
          )
      ) %>%
      select(-data),
    d_nca_initial_calc_noerror=
      d_nca_prep %>%
      mutate(
        noerror=
          pmap(
            .l=list(data=data, id=ID),
            .f=fit_half_life_all_df,
            conc_col="CP_noerror",
            lloq_col=0
          )
      ) %>%
      select(-data),
    d_nca_unnest_loq20=
      d_nca_initial_calc_loq20 %>%
      unnest(col="LOQ20pct"),
    d_nca_unnest_loq00=
      d_nca_initial_calc_loq00 %>%
      unnest(col="LOQ00pct"),
    d_nca_unnest_noerror=
      d_nca_initial_calc_noerror %>%
      unnest(col="noerror"),
    d_nca_select_loq20=make_d_nca_select(d_nca_unnest_loq20, mantissa=1),
    d_nca_select_loq00=make_d_nca_select(d_nca_unnest_loq00, mantissa=1),
    d_nca_select_noerror=make_d_nca_select(d_nca_unnest_noerror, mantissa=1),
    d_nca_select_summary_loq20=make_d_nca_select_summary(d_nca_select_loq20),
    d_nca_select_summary_loq00=make_d_nca_select_summary(d_nca_select_loq00),
    d_nca_select_summary_noerror=make_d_nca_select_summary(d_nca_select_noerror),
    best_mantissa=optimize_mantissa(test=d_nca_unnest_loq20, reference=d_true_half_life),
    best_mantissa_no_tmdd=
      optimize_mantissa(
        test=
          d_nca_unnest_loq20 %>%
          filter(ID %in% (d_dose %>% filter(TMDD_Text %in% "without TMDD"))$ID),
        reference=
          d_true_half_life %>%
          filter(ID %in% (d_dose %>% filter(TMDD_Text %in% "without TMDD"))$ID)
      )
  )

plan_benchmark <-
  drake_plan(
    benchamark_least_squares=
      microbenchmark::microbenchmark(
        PKNCA:::fit_half_life(
          data=
            tibble(
              log_conc=log(d_nca_prep$data[[1]]$CP_LOQ20pct),
              time=d_nca_prep$data[[1]]$time,
              log_lloq=log(d_nca_prep$data[[1]]$LOQ20pct),
              mask_blq=log_conc < log_lloq
            ) %>%
            filter(!mask_blq),
          tlast=
            pk.calc.tlast(
              conc=d_nca_prep$data[[1]]$CP_LOQ20pct,
              time=d_nca_prep$data[[1]]$time
            )
        ),
        times=1000
      ),
    benchmark_tobit=
      microbenchmark::microbenchmark(
        fit_half_life_tobit(
          data=
            tibble(
              log_conc=log(d_nca_prep$data[[1]]$CP_LOQ20pct),
              time=d_nca_prep$data[[1]]$time,
              log_lloq=log(d_nca_prep$data[[1]]$LOQ20pct),
              mask_blq=log_conc < log_lloq
            ),
          tlast=
            pk.calc.tlast(
              conc=d_nca_prep$data[[1]]$CP_LOQ20pct,
              time=d_nca_prep$data[[1]]$time
            )
        ),
        times=1000
      ),
    benchmark_tobit_ls_ratio=
      median(benchmark_tobit$time)/median(benchamark_least_squares$time)
  )

plan_figures <-
  drake_plan(
    d_plot_prep=make_d_plot(
      half_life=d_true_half_life,
      dose=d_dose,
      blq=d_blq_summary,
      loq00=d_nca_select_summary_loq00,
      loq20=d_nca_select_summary_loq20,
      noerror=d_nca_select_summary_noerror
    ),
    count_missing=verify_and_count_missing(d_plot_prep),
    count_negative=make_count_negative(d_plot),
    d_plot=d_plot_prep[!is.na(d_plot_prep$half.life_tobit_loq20),],
    p_figure_1=make_figure1(d_plot),
    p_figure_1_tmdd=make_figure1_tmdd(p_figure_1),
    p_figure_s1_early_blq=make_figure1_early_blq(p_figure_1),
    p_figure_s2_late_blq=make_figure1_late_blq(p_figure_1),
    p_figure_2=make_figure2(d_plot),
    p_figure_2_tmdd=make_figure2_tmdd(p_figure_2),
    save_figure_1=
      ggsave(
        filename=file_out("../../Output/Figures/figure_1.eps"),
        plot=p_figure_1 + theme(legend.position="none"),
        width=4, height=4, units="in"
      ),
    save_figure_2=
      ggsave(
        filename=file_out("../../Output/Figures/figure_2.eps"),
        plot=p_figure_2 + theme(legend.position="none"),
        width=4, height=4, units="in"
      ),
    save_figure_s1=
      ggsave(
        filename=file_out("../../Output/Figures/figure_s1.eps"),
        plot=p_figure_s1_early_blq + theme(legend.position="none"),
        width=4, height=7, units="in"
      ),
    save_figure_s2=
      ggsave(
        filename=file_out("../../Output/Figures/figure_s2.eps"),
        plot=
          p_figure_s2_late_blq +
          theme(
            legend.position="none",
            strip.text=element_text(size=rel(0.7))
          ),
        width=4, height=7, units="in"
      )
  )

make(
  bind_plans(
    plan_scenarios,
    plan_nca,
    plan_benchmark,
    plan_figures
  ),
  seed=5
)
