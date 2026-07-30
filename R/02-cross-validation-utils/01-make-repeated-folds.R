make_repeated_folds <- function(participant_data,
                                n_folds = 10,
                                n_reps = 5){
  folds_data <- participant_data |>
    select(csid, group) |> 
    crossing(rep = 1:n_reps) |> 
    group_by(rep, group) |>
    mutate(test_fold = sample(rep(seq.int(n_folds), length = n()))) |>
    ungroup() |> 
    select(csid, rep, test_fold, group)
  
  for (i in seq.int(n_folds)){
    non_fold <- filter(folds_data, test_fold != i)
    inner_fold_col_name <- glue::glue("inner_fold{i}")
    innerfolds <- non_fold |>
      group_by(rep, group) |>
      mutate("{inner_fold_col_name}" := sample(rep(seq.int(n_folds), length = n()))) |>
      ungroup() |>
      select(csid, rep, all_of(inner_fold_col_name))
    folds_data <- left_join(folds_data, innerfolds, by = c("csid", "rep"))
  }
  
  folds_data <- folds_data |>
    select(-group) |> 
    arrange(rep, test_fold)
  
  return(folds_data)
}
