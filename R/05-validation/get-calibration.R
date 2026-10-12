
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
  
  calibration_data <- predictions |> 
    unnest(preds) |> 
    mutate(linear_pred_calibrated = linear_pred + log(status1_weight / status0_weight)) |> 
    mutate(prob = exp(linear_pred_calibrated) / (1 + exp(linear_pred_calibrated))) |>
    group_by(rep, test_fold) |> 
    mutate(prob_group = cut(prob, c(0, quantile(prob, seq(0.1, 0.9, 0.1)), 1), labels = FALSE)) |> 
    ungroup() |> 
    left_join(casesubcohort_sampling_weights,
              by = "csid",
              relationship = "many-to-many") |> 
    filter(!is.na(weight))
  
  # observed risk vs. mean predicted risk 
  group_results <- calibration_data |> 
    nest_by(prob_group) |> 
    mutate(obs_prob = list(get_obs_prop(data)),
           ave_pred_prob = plogis(mean(data$linear_pred_calibrated))) |>
    unnest_wider(obs_prob) |> 
    select(prob_group, obs_prop, obs_prop_lci, obs_prop_uci, ave_pred_prob)
  
  # calibration intercept and slope
  design <- survey::svydesign(ids = ~ csid,
                              weight = ~ weight,
                              data = calibration_data)
  reg_intercept <- survey::svyglm(status ~ offset(linear_pred_calibrated),
                                  design = design,
                                  family = "binomial")
  reg_slope <- survey::svyglm(status ~ linear_pred_calibrated,
                              design = design,
                              family = "binomial")
  intercept_slope <- bind_rows(
    data.frame(measure = "intercept",
               coef = coef(reg_intercept)[[1]],
               lci = confint(reg_intercept)[[1,1]],
               uci = confint(reg_intercept)[[1,2]],
               row.names = NULL),
    data.frame(measure = "slope",
               coef = coef(reg_slope)[[2]],
               lci = confint(reg_slope)[[2,1]],
               uci = confint(reg_slope)[[2,2]],
               row.names = NULL)
  )
  
  return(lst(group_results, intercept_slope))
}
