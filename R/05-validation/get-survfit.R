get_survfit <- function(predictions_data,
                        casesubcohort_participant_data,
                        casesubcohort_sampling_weights,
                        cutoffs = c(1/5, 2/5, 3/5, 4/5),
                        adjust = "age_at_study_date*region"){
  
  casesubcohort_sampling_weights <- casesubcohort_sampling_weights[[2]]
  
  ## function to make adjusted risk predictions  
  make_adjusted_risk <- function(x, adjust){
    x$adjusted_risk <- resid(lm(
      reformulate(adjust, "linear_pred"),
      data = left_join(x,
                       select(casesubcohort_participant_data,
                              csid, region, age_at_study_date, sex),
                       by = "csid")))
    return(x)
  }
  
  survival_data <- predictions_data |> 
    
    ## add adjusted risk predictions
    rowwise() |> 
    mutate(preds = list(make_adjusted_risk(preds, adjust = .env$adjust))) |> 
    
    ## create predicted risk groups (within each validation data set, and within region)
    unnest(preds) |> 
    left_join(select(casesubcohort_participant_data, csid, region), by = "csid") |> 
    group_by(rep, test_fold) |> 
    mutate(pred_group = cut(adjusted_risk,
                            c(-Inf, quantile(adjusted_risk, cutoffs), Inf),
                            labels = c("Low", "Low-Middle", "Middle", "High-Middle", "High"))) |> 
    ungroup() |>
    
    ## join sampling weights and follow up time
    inner_join(casesubcohort_sampling_weights, by = "csid", relationship = "many-to-many") |> 
    
    ## collapse to have one observation per participant and risk group
    ## multiply sampling weights by number of times participant appears in
    ## that risk group (zero to number of repeats of cross-validation)
    select(csid, rep, test_fold, pred_group, time_in, time_out, status, weight) |> 
    group_by(csid, pred_group, time_in, time_out, status, weight) |> 
    summarise(n = n(),
              .groups = "drop") |> 
    mutate(weight = weight*n) |> 
    arrange(time_out, csid)
  
  ## survfit stratified by risk group, up to 10 years
  broom::tidy(survfit(Surv(time_in, time_out, status) ~ strata(pred_group),
                      data = survival_data,
                      weights = weight,
                      id = csid,
                      robust = TRUE)) |>
    mutate(pred_group = stringr::str_extract(strata,
                                             "(strata\\(pred_group\\)=)(.+)",
                                             group = 2)) |> 
    dplyr::filter(time < 365.25*10)
}
