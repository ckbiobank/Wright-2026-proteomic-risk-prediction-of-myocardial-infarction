train_lasso_logistic <- function(train_data,
                                 predictors,
                                 zero_penalty_factor = "^(region|prs_|somalogic_cv_risk_score|ukb_olink_score)"){
  x <- model.matrix(reformulate(predictors), train_data)  ## create matrix of predictors
  x <- x[,colnames(x)[colnames(x) != "(Intercept)"]]      ## remove any intercept column
  
  penalty_factor <- as.numeric(!grepl(zero_penalty_factor, colnames(x)))
  if (sum(penalty_factor) <= 1){
    return(train_logistic_regression(train_data, predictors))
  }
  
  y <- train_data$status
  
  fit <- glmnet::cv.glmnet(x,
                           y,
                           family = "binomial",
                           data = train_data,
                           foldid = train_data$inner_fold,
                           type.measure = "auc",
                           intercept = (predictors[[1]] != "-1"), ## include intercept first predictor is not "-1"
                           penalty.factor = penalty_factor)
  
  return(fit)
}
