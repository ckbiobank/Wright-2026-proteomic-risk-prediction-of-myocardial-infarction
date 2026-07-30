get_predictions <- function(test_data, model, predictors){
  if (inherits(model, "cv.glmnet")){
    x <- model.matrix(reformulate(predictors), test_data)  ## create matrix of predictors
    x <- x[,colnames(x)[colnames(x) != "(Intercept)"]]     ## remove any intercept column
    
    preds <- data.frame(csid = test_data$csid,
                        linear_pred = predict(model,
                                              newx = x,
                                              s = "lambda.1se",
                                              type = "link"))
    names(preds) <- c("csid", "linear_pred")
  }
  if (inherits(model, "glm")){
    preds <- data.frame(csid = test_data$csid,
                        linear_pred = predict(model,
                                              newdata = test_data))
  }
  
  return(preds)
}
