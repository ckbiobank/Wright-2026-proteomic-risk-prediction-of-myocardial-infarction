
get_boruta_selected <- function(proteins,
                                analysis_data){
  
  analysis_data <- analysis_data |>
    select(status, all_of(proteins))
  
  boruta_res <- Boruta(status ~ .,
                       data = analysis_data,
                       ntree = 100,
                       maxRuns = 400)
  
  boruta_selected <- names(boruta_res$finalDecision)[boruta_res$finalDecision == "Confirmed"]
  attr(boruta_selected, "Tentative") <- names(boruta_res$finalDecision)[boruta_res$finalDecision == "Tentative"] 
  
  return(boruta_selected)
}
