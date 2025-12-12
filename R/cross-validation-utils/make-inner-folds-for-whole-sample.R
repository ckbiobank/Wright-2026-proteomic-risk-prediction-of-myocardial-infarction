make_inner_folds_for_whole_sample <- function(participant_data,
                                              n_folds = 10,
                                              n_reps = 25){
  participant_data |>
    crossing(rep = 1:n_reps) |> 
    select(csid, group, rep) |> 
    group_by(rep, group) |>
    mutate(inner_fold1 = sample(rep(seq.int(n_folds), length = n()))) |>
    ungroup() |> 
    select(csid, rep, inner_fold1) |> 
    arrange(rep, inner_fold1) |> 
    mutate(test_fold = 0)
}
