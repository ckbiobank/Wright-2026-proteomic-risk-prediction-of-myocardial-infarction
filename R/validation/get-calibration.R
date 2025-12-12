
get_calibration <- function(predictions, casesubcohort_sampling_weights){
  
  casesubcohort_sampling_weights <- casesubcohort_sampling_weights[[1]]
  
  get_obs_prop <- function(x){
    est <- sum(x$status*x$weight)/sum(x$weight)
    
    design <- survey::svydesign(ids = ~csid,
                                weight = ~ weight,
                                data = x)
    
    ci <- as.numeric((attr(survey::svyciprop(~ status, design), "ci")))
    
    return(data.frame(obs_prop = est,
                      obs_prop_lci = ci[[1]],
                      obs_prop_uci = ci[[2]]))  
  }
  
  status1_weight <- casesubcohort_sampling_weights[casesubcohort_sampling_weights$status == 1]$weight[[1]]
  status0_weight <- casesubcohort_sampling_weights[casesubcohort_sampling_weights$status == 0]$weight[[1]]
  
  predictions |> 
    unnest(preds) |> 
    mutate(linear_pred_calibrated = linear_pred + log(status1_weight / status0_weight)) |> 
    mutate(prob = exp(linear_pred_calibrated) / (1 + exp(linear_pred_calibrated))) |>
    group_by(rep, test_fold) |> 
    mutate(prob_group = cut(prob, c(0, quantile(prob, seq(0.1, 0.9, 0.1)), 1), labels = FALSE)) |> 
    left_join(casesubcohort_sampling_weights,
              by = "csid",
              relationship = "many-to-many") |> 
    filter(!is.na(weight)) |> 
    group_by(csid, prob, prob_group, status, weight) |> 
    summarise(n = n(),
              .groups = "drop") |> 
    mutate(weight = weight*n) |> 
    nest_by(prob_group) |> 
    mutate(obs_prob = list(get_obs_prop(data)),
           ave_pred_prob = mean(data$prob)) |>
    unnest_wider(obs_prob) |> 
    select(prob_group, obs_prop, obs_prop_lci, obs_prop_uci, ave_pred_prob)
  
}
