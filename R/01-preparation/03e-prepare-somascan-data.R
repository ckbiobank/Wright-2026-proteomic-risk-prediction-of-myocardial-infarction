prepare_somascan_data <- function(somalogic, somalogic_samples, somalogic_prediction){
  # create sample information (plateid, hybcontrolnormscale) columns
  sample_info <- somalogic_samples |> 
    select(csid, qc, plateid, hybcontrolnormscale) |>
    
    ## use first measurements where participants samples have been split
    filter(qc == 0) |>
    
    mutate(plateid = factor(plateid))
  
  # SomaLogic risk scores
  risk_scores <- somalogic_prediction |> 
    filter(qc == 0) |>
    select(csid,
           somalogic_cv_risk_score = predict_primary_cardiovascular_risk_4yr_per,
           somalogic_cv_risk_score_2 = predict_secondary_cardiovascular_risk_4yr_per)
  
  # prepare aptamer data
  somalogic <- somalogic |>  
    
    ## remove any non-participant data
    filter(!is.na(csid)) |>  
    
    ## use first measurements where participants samples have been split
    filter(qc == 0) |>
    
    ## log_2 transform protein values (RFU)
    mutate(value = log2(value)) |> 
    select(csid, variable, value) |> 
    pivot_wider(names_from = "variable",
                values_from = "value") |> 
    
    mutate(apob = 2^(`seq.2797.56`),
           apoa = 2^(`seq.2750.3`)) |> 
    mutate(apob_apoa_ratio = apob/apoa) |> 
    
    left_join(sample_info, by = "csid") |>
    
    left_join(risk_scores, by = "csid") |> 
    
    ## exclude participant samples with hybcontrolnormscale
    ## out of range (i.e. less than 0.4 or greater than 2.5)
    filter(hybcontrolnormscale >= 0.4,
           hybcontrolnormscale <= 2.5)
  
  return(somalogic)
}

