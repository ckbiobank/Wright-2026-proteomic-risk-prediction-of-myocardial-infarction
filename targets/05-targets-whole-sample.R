
models <- tibble::tribble(
  ~name, ~protein_data, ~predictors,
  "conventional_whole_sample",           "null_protein_data", "predictors_conventional",
  "conventional_sig_olink_whole_sample", "olink_data_scaled", "predictors_conventional_sig_olink_whole_sample"
) |> 
  dplyr::mutate(protein_data = rlang::syms(protein_data),
                predictors = rlang::syms(predictors))





targets_whole_sample <- tar_plan(
  
  # individual protein associations with MI
  tar_target(olink_individual_proteins_analysis_whole_sample,
             individual_proteins_analysis(inner_join(casesubcohort_participant_data,
                                                     olink_data,
                                                     by = "csid"),
                                          olink_meta_data)),
  tar_target(somascan_individual_proteins_analysis_whole_sample,
             individual_proteins_analysis(inner_join(casesubcohort_participant_data,
                                                     somascan_data,
                                                     by = "csid"),
                                          somascan_meta_data)),
  
  
  # prediction using whole sample
  inner_folds_whole_sample = make_inner_folds_for_whole_sample(casesubcohort_participant_data,
                                                               n_folds = 10,
                                                               n_reps = 25),
  
  tar_group_count(reps_whole_sample, distinct(inner_folds_whole_sample, rep), 5),
  
  tar_map(
    values = models,
    names = "name",
    
    ## train models
    tar_target(trained_model,
               mutate(reps_whole_sample,
                      model = purrr::map(
                        rep,
                        \(rep) {
                          training_data <- create_training_data(rep,
                                                                test_fold = 1,
                                                                inner_folds_whole_sample,
                                                                casesubcohort_participant_data,
                                                                protein_data)
                          train_lasso_logistic(training_data, c("-1", "region", predictors))
                        })),
               pattern = map(reps_whole_sample)),
    
    ## select median lambda across repeats
    ## (all models with this lambda will be the same)
    tar_target(
      selected_model,
      trained_model |>
        mutate(lambda.1se = map_dbl(model,
                                    \(x) if (!is.null(x$lambda.1se)) x$lambda.1se else NA_real_)) |> 
        arrange(lambda.1se) |> 
        filter(row_number() == floor((n() + 1) / 2)) |> 
        pull(model) |> 
        first()
    ),
    
    tar_target(
      predictions,
      get_predictions(inner_join(casesubcohort_participant_data, protein_data),
                      selected_model,
                      c("-1", "region", predictors))
    )
  )
)
