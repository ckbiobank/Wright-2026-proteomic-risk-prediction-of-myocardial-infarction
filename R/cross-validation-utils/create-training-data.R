create_training_data <- function(rep, test_fold, folds_data, participant_data, ..., .impute = TRUE){
  
  train_data <- folds_data |> 
    filter(rep == .env$rep, test_fold != .env$test_fold) |> 
    mutate(inner_fold = .data[[glue::glue("inner_fold{test_fold}")]]) |> 
    select(csid, inner_fold) |> 
    inner_join(participant_data, by = "csid")
  
  if (!is.null(...)){
    for (l in list(...)){
      train_data <- train_data |> 
        inner_join(l, by = "csid")
    }
  }
  
  if (.impute){
    train_data <- train_data |>
      mutate(across(starts_with(c("oid", "seq", "somalogic_cv_risk_score")),
                    \(x) coalesce(x, mean(x, na.rm = TRUE))))
  }
  
  return(train_data)
}
