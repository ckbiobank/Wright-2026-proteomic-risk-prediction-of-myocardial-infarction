make_casesubcohort_sampling_weights <- function(all_participant_data,
                                                casesubcohort_participant_data){
  
  all_participant_data <- all_participant_data |> 
    mutate(status = acute_mi,
           endpoint_date = acute_mi_date,
           time_out = as.numeric(endpoint_date - study_date))
  
  ## 'cross-sectional' sampling weights
  status0_sampling_frac <- sum(casesubcohort_participant_data$status == 0) / sum(all_participant_data$status == 0)
  status1_sampling_frac <- sum(casesubcohort_participant_data$status == 1) / sum(all_participant_data$status == 1)
  
  casesubcohort_participant_data_with_weights <- casesubcohort_participant_data |> 
    select(csid, time_in, time_out, status, ascertainment_group)
  
  casesubcohort_participant_data_with_weights$weight[casesubcohort_participant_data_with_weights$status == 0] <- 1 / status0_sampling_frac
  casesubcohort_participant_data_with_weights$weight[casesubcohort_participant_data_with_weights$status == 1] <- 1 / status1_sampling_frac
  
  ## 'time-to-event data' sampling weights
  case_sampling_frac <- tibble(period_start = seq(0, 6000, 100),
                               period_end = period_start + 100) |> 
    rowwise() |> 
    mutate(total = all_participant_data |>
             filter(status == 1, time_out > period_start, time_out <= period_end) |>
             nrow(),
           selected = casesubcohort_participant_data |>
             filter(status == 1, time_out > period_start, time_out <= period_end) |>
             nrow()) |> 
    mutate(sampled_frac = selected/total,
           tstart = period_start)
  
  cases_with_weights <- casesubcohort_participant_data |> 
    select(csid, time_in, time_out, status) |> 
    filter(status == 1) |> 
    mutate(time_in = time_out - 0.001) |>
    mutate(case_sampling_frac = approx(case_sampling_frac$period_start+0.001,
                                       case_sampling_frac$sampled_frac,
                                       time_out, method = "constant")[[2]]) |> 
    mutate(weight = 1 / case_sampling_frac) |> 
    select(csid, time_in, time_out, status, weight)
  
  
  subcohort_weights <- tibble(period_start = seq(0, 6000, 100),
                              period_end = period_start + 100) |> 
    rowwise() |> 
    mutate(total = all_participant_data |>
             filter(time_out > period_start) |>
             nrow(),
           selected = casesubcohort_participant_data |>
             filter(ascertainment_group == "subcohort", time_out > period_start) |>
             nrow()) |> 
    mutate(sampled_frac = selected/total,
           weight = 1/sampled_frac)
  
  subcohort_with_weights <- casesubcohort_participant_data |>
    filter(ascertainment_group == "subcohort") |> 
    crossing(subcohort_weights) |> 
    filter(time_out > period_start) |> 
    mutate(time_in = period_start,
           time_out = pmin(time_out, period_end),
           status = 0) |> 
    select(csid, time_in, time_out, status, weight)
  
  
  casesubcohort_sampling_weights <- rbind(subcohort_with_weights,
                                          cases_with_weights)
  return(lst(casesubcohort_participant_data_with_weights,
             casesubcohort_sampling_weights))
}
