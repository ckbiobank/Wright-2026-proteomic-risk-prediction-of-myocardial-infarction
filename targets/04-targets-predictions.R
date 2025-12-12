
models <- tibble::tibble(
  name = c(
    "age_sex",
    "conventional",
    
    "all_olink",
    "sig_olink",
    "boruta_olink",
    "all_somascan",
    "sig_somascan",
    "boruta_somascan",
    "all_both",
    "sig_both",
    "boruta_both",
    
    "conventional_all_olink",
    "conventional_sig_olink",
    "conventional_boruta_olink",
    "conventional_all_somascan",
    "conventional_sig_somascan",
    "conventional_boruta_somascan",
    "conventional_all_both",
    "conventional_sig_both",
    "conventional_boruta_both",
    
    "ukb_olink_score",
    "conventional_ukb_olink_score"
  )
) |> 
  dplyr::mutate(protein_data = dplyr::case_when(
    grepl("_olink$", name)         ~ rlang::syms("olink_data_scaled"),
    grepl("_somascan$", name)      ~ rlang::syms("somascan_data_scaled"),
    grepl("_both$", name)          ~ rlang::syms("olink_somascan_data_scaled"),
    grepl("ukb_olink_score", name) ~ rlang::syms("ukb_olink_score"),
    .default                       = rlang::syms("null_protein_data"))) |>
  dplyr::mutate(protein_data = rlang::syms(protein_data),
                predictors =  rlang::syms(paste0("predictors_", name)))

targets_predictions <- tar_plan(
  
  
  folds_data = make_repeated_folds(casesubcohort_participant_data,
                                   n_folds = 10,
                                   n_reps = 5),
  
  tar_group_count(folds_and_reps, distinct(folds_data, rep, test_fold), 5),
  
  ## cross-validation
  tar_map(
    values = models,
    names = "name",
    
    tar_target(trained_model,
               mutate(folds_and_reps,
                      model = purrr::pmap(list(rep, test_fold),
                                          \(rep, test_fold) {
                                            training_data <- create_training_data(rep, test_fold, folds_data,
                                                                                  casesubcohort_participant_data, protein_data)
                                            train_lasso_logistic(training_data,
                                                                 c("-1", "region",
                                                                   filter_by_rep_test_fold(rep, test_fold, predictors, "predictors")))
                                          })),
               pattern = map(folds_and_reps)),
    
    tar_target(predictions,
               mutate(trained_model,
                      preds = purrr::pmap(list(model, rep, test_fold),
                                          \(model, rep, test_fold) {
                                            test_data <- create_test_data(rep, test_fold, folds_data,
                                                                          casesubcohort_participant_data, protein_data)
                                            get_predictions(test_data,
                                                            model,
                                                            c("-1", "region",
                                                              filter_by_rep_test_fold(rep, test_fold, predictors, "predictors")))
                                          })) |>
                 select(rep, test_fold, preds),
               pattern = map(trained_model)),
    
    ## validation
    tar_target(c_statistic, get_c_statistic(predictions, predictions_conventional, casesubcohort_participant_data)),
    tar_target(nri, get_nri(predictions, predictions_conventional, casesubcohort_participant_data)),
    tar_target(survfit, get_survfit(predictions, casesubcohort_participant_data, casesubcohort_sampling_weights)),
    tar_target(calibration, get_calibration(predictions, casesubcohort_sampling_weights))
  )
)
