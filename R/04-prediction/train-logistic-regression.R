train_logistic_regression <- function(train_data, predictors){
  formula <- reformulate(predictors, "status")
  
  fit <- glm(formula,
             family = binomial(link = "logit"),
             data = train_data)
  
  return(fit)
}
