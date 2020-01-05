define_loq <- function(x, fraction) {
  if (fraction > 0 & fraction < 1) {
    # the fraction is a fraction of nonzero values
    quantile(x=x[x > 0], probs=fraction)
  } else if (fraction == 0) {
    -Inf
  } else if (fraction == 1) {
    Inf
  } else {
    stop("Invalid fraction value")
  }
}

make_d_scenario_ncompartments <- function(ratios) {
  mutate_at(
    .tbl=
      bind_rows(
        crossing(
          Vcentral=1,
          Vperiph1_Vcentral_ratio=1,
          Vperiph2_Vcentral_ratio=1,
          Qcp1_Vcentral_ratio=0,
          Qcp2_Vcentral_ratio=0,
          CL_Vcentral_ratio=ratios,
          Compartment_Text="One compartment"
        ),
        crossing(
          Vcentral=1,
          Vperiph1_Vcentral_ratio=ratios,
          Vperiph2_Vcentral_ratio=1,
          Qcp1_Vcentral_ratio=ratios,
          Qcp2_Vcentral_ratio=0,
          CL_Vcentral_ratio=ratios,
          Compartment_Text="Two compartment"
        ),
        crossing(
          Vcentral=1,
          Vperiph1_Vcentral_ratio=ratios,
          Vperiph2_Vcentral_ratio=ratios,
          Qcp1_Vcentral_ratio=ratios,
          Qcp2_Vcentral_ratio=ratios,
          CL_Vcentral_ratio=ratios,
          Compartment_Text="Three compartment"
        )
      ),
    .vars="Compartment_Text", .funs=fct_inorder
  )
}

make_d_scenarios_tmdd <- function(ratios) {
  mutate_at(
    bind_rows(
      crossing(
        CLTMDD_CL_ratio=0,
        C50TMDD=1,
        TMDD_Text="without TMDD"
      ),
      crossing(
        CLTMDD_CL_ratio=ratios,
        C50TMDD=1,
        TMDD_Text="with TMDD"
      )
    ),
    .vars="TMDD_Text", .funs=fct_inorder
  )
}

make_d_scenarios_route <- function(ratios) {
  mutate_at(
    bind_rows(
      # Oral
      tibble(
        Ka_Kel_ratio=ratios,
        rate_Kel_ratio=0,
        dose_cmt="GUT",
        Route_Text="extravascular"
      ),
      # IV bolus
      tibble(
        Ka_Kel_ratio=0,
        rate_Kel_ratio=0,
        dose_cmt="CENT",
        Route_Text="intravascular bolus"
      )
    ),
    .vars="Route_Text", .funs=fct_inorder
  )
}

make_d_scenarios_residual <- function() {
  mutate_at(
    .tbl=
      crossing(
        PROP=c(0.0025, 0.01, 0.04),
        ADD=0
      ) %>%
      mutate(
        Variability_Text=sprintf("%g prop variance, %g add variance", PROP, ADD)
      ),
    .vars="Variability_Text", .funs=fct_inorder
  )
}

make_d_scenarios <- function(data) {
  data %>%
    mutate(
      CL=Vcentral*CL_Vcentral_ratio,
      Ka=CL/Vcentral*Ka_Kel_ratio,
      Vperiph1=Vcentral*Vperiph1_Vcentral_ratio,
      Vperiph2=Vcentral*Vperiph2_Vcentral_ratio,
      Qcp1=Vcentral*Qcp1_Vcentral_ratio,
      Qcp2=Vcentral*Qcp2_Vcentral_ratio,
      CLTMDD=CL*CLTMDD_CL_ratio,
      rate=CL/Vcentral*rate_Kel_ratio
    )
}

make_d_obs_prep <- function() {
  crossing(
    cmt="CENT",
    evid=0,
    time=c(0, 0.5, 1, 2, 3, 4, 8, 12, 16, 24, 36, 48),
    rate=0,
    amt=0
  )
}

make_d_dose <- function(data, n) {
  data %>%
    crossing(
      ID_sub=seq(from=1, to=n, by=1)
    ) %>%
    mutate(
      cmt=dose_cmt,
      evid=1,
      time=0,
      amt=1,
      ID=row_number()
    )
}

make_d_simresult <- function(o_simresult, d_dose, model) {
  left_join(
    as.data.frame(o_simresult),
    unique(
      d_dose[,
             intersect(
               names(d_dose),
               c("ID",
                 names(param(model)),
                 grep(
                   x=names(d_dose),
                   pattern="Text",
                   fixed=TRUE,
                   value=TRUE
                 )
               )
             )
             ]
    )
  ) %>%
    unique() %>%
    group_by_at(.vars=c("CL", "Ka", "Vperiph1", "Vperiph2", "Qcp1", "Qcp2", "CLTMDD")) %>%
    # Define the limits of quantification 
    mutate(
      CP=pmax(CP, 0),
      CP_noerror=pmax(CP_noerror, 0),
      LOQ05pct=define_loq(CP, fraction=0.05),
      LOQ10pct=define_loq(CP, fraction=0.10),
      LOQ20pct=define_loq(CP, fraction=0.20)
    ) %>%
    ungroup() %>%
    mutate(
      CP_LOQ05pct=
        case_when(
          CP <= LOQ05pct~0,
          TRUE~CP
        ),
      CP_LOQ10pct=
        case_when(
          CP <= LOQ10pct~0,
          TRUE~CP
        ),
      CP_LOQ20pct=
        case_when(
          CP <= LOQ20pct~0,
          TRUE~CP
        )
    )
}
