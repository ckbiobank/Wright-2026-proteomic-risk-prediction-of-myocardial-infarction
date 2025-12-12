individual_proteins_analysis <- function(analysis_data, meta_data){
  
  fit_model_raw <- function(data,
                            protein_id,
                            init_values = c(0)){
    adjust_for = c("strata(sex)",
                   "strata(region)",
                   "poly(age_at_study_date, 2)",
                   "poly(hours_since_last_ate, 2)",
                   "poly(mean_temp, 2)",
                   "highest_education",
                   "plateid",
                   "smoking_category",
                   "alcohol_category",
                   "met",
                   "sbp_mean",
                   "has_diabetes",
                   "bmi_calc")
    
    data <- mutate(data, value = .data[[protein_id]])
    
    ## scale values by standard deviation in subcohort + cases
    sd <- data |>
      pull(value) |>
      sd(na.rm = TRUE)
    data$value_scaled <- data$value / sd
    
    ## create init values
    formula_no_int <- reformulate(c("-1", "value_scaled", adjust_for),
                                  "Surv(time_in, time_out, status)")
    model_matrix <- model.matrix(formula_no_int, data = data)
    col_names <- colnames(model_matrix)
    col_names <- col_names[!grepl("^strata\\(", col_names)]
    init_values_aligned <- init_values[match(col_names, names(init_values))]
    init_values_aligned[is.na(init_values_aligned)] <- 0
    
    ## model formula
    formula <- reformulate(c("value_scaled", adjust_for),
                           "Surv(time_in, time_out, status)")
    
    ## fit the model
    fit <- coxph(formula,
                 data = data,
                 cluster = csid,
                 init = init_values_aligned,
                 y = FALSE,
                 iter.max = 40)
    est <- broom::tidy(fit) |> 
      filter(term == "value_scaled") |> 
      select(-term)
    return(est)
  }
  
  fit_model_safely <- function(data,
                               protein_id,
                               init_values = c(0)){
    x <- purrr::quietly(purrr::safely(fit_model_raw))(data,
                                                      protein_id,
                                                      init_values)
    r <- x[["result"]][["result"]]
    r$e <- x[["result"]][["error"]] |> length_zero_to_na()
    r$w <- x[["warnings"]] |> length_zero_to_na()
    r$m <- x[["messages"]] |> length_zero_to_na()
    
    return(r)
  }
  
  make_init_values_for_models <- function(participant_data){
    
    formula <- reformulate(c("poly(age_at_study_date, 2)",
                             "poly(hours_since_last_ate, 2)",
                             "poly(mean_temp, 2)",
                             "highest_education",
                             "plateid",
                             "smoking_category",
                             "alcohol_category",
                             "met",
                             "sbp_mean",
                             "has_diabetes",
                             "bmi_calc"),
                           "Surv(time_in, time_out, status)")
    
    fit0 <- survival::coxph(formula,
                            data = participant_data,
                            cluster = csid,
                            y = FALSE,
                            iter.max = 40)
    
    return(fit0$coefficients)
  }
  
  analysis_data <- analysis_data |> 
    filter(!is.na(mean_temp))
  
  init_values <- make_init_values_for_models(analysis_data)
  result <- meta_data |> 
    select(protein_id) |> 
    mutate(result = purrr::map(protein_id,
                               \(x) fit_model_safely(data = analysis_data,
                                                     protein_id = x,
                                                     init_values = init_values))) |> 
    unnest_wider(result)
  return(result)
}

length_zero_to_na <- function(obj){
  if (is.null(obj)){
    return(NA_character_)
  }
  
  obj <- as.character(obj)
  if(length(obj)==0){
    return(NA_character_)
  }else if(length(obj) == 1){
    return(obj)
  }else{
    return(paste(obj, collapse = "; "))
  }
}

