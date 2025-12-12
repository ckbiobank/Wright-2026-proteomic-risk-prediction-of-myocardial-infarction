filter_by_rep_test_fold <- function(rep, test_fold, x, col){
  if (is.vector(x)) {
    return(x)
  }
  
  x |> 
    filter(rep == .env$rep, test_fold == .env$test_fold) |> 
    pull(.data[[col]]) |> 
    first()
}
