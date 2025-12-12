create_test_data <- function(rep, test_fold, folds_data, participant_data, ...){
  
  test_data <- folds_data |> 
    filter(rep == .env$rep, test_fold == .env$test_fold) |> 
    inner_join(participant_data, by = "csid")
  
  if (!is.null(...)){
    for (l in list(...)){
      test_data <- test_data |> 
        inner_join(l, by = "csid")
    }
  }
  
  test_data <- test_data |>   
    mutate(across(starts_with(c("oid", "seq", "somalogic_cv_risk_score")),
                  \(x) coalesce(x, mean(x, na.rm = TRUE))))
  
  return(test_data)
}
