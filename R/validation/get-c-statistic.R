get_c_statistic <- function(predictions, reference_predictions, participant_data){
  predictions |> 
    mutate(c = purrr::map_dbl(preds, \(x) calculate_c_stat(x, participant_data))) |> 
    select(-preds) |> 
    left_join(
      reference_predictions |> 
        mutate(ref_c = purrr::map_dbl(preds, \(x) calculate_c_stat(x, participant_data))) |> 
        select(-preds),
      by = c("rep","test_fold")
    ) |> 
    mutate(delta_c = c - ref_c) |> 
    group_by(rep) |> 
    summarise(se = sd(c)/sqrt(n()),
              se_delta = sd(delta_c)/sqrt(n()),
              c = mean(c),
              delta_c = mean(delta_c),
              .groups = "drop") |>
    summarise(c = mean(c),
              se = sqrt(mean(se^2)),
              delta_c = mean(delta_c),
              se_delta = sqrt(mean(se_delta^2))) |> 
    mutate(lci = c - 1.96*se,
           uci = c + 1.96*se,
           lci_delta = delta_c - 1.96*se_delta,
           uci_delta = delta_c + 1.96*se_delta)
}

calculate_c_stat <- function(preds, participant_data, min_t = 0, max_t = 20){
  x <- preds |> 
    inner_join(select(participant_data, csid, region, time_in, time_out, status), by = "csid")
  c <- concordance(Surv(time_in, time_out, status) ~ linear_pred + strata(region),
                   data       = x,
                   cluster    = csid,
                   ymin       = 365.25*min_t,
                   ymax       = 365.25*max_t,
                   keepstrata = TRUE,
                   reverse    = TRUE)
  
  return(c$concordance)
}
