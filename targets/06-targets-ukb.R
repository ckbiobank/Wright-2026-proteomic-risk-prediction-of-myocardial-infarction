targets_ukb <- tar_plan(
  tar_file_read(ukb_raw_data,
                file.path(ukb_data_source, "data20062023.RData"),
                {attach(!!.x); return(data20062023)}),
  tar_file_read(ukb_raw_data_plateid,
                file.path(ukb_data_source, "olink_batch_plateid.RData"),
                {attach(!!.x); return(batch_plateid)}),
  
  tar_file_read(ukb_raw_data_ii,
                file.path(ukb_data_source, "2ndBatch", "olink_proteins.csv"),
                readr::read_csv(!!.x)),
  tar_file_read(ukb_raw_data_plateid_ii,
                file.path(ukb_data_source, "2ndBatch", "olink_batch_plateid.csv"),
                readr::read_csv(!!.x)),
  
  
  ukb_assay_data = prepare_ukb_assay_data(ukb_raw_data, ukb_raw_data_ii, olink_meta_data),
  ukb_participant_data = prepare_ukb_participant_data(ukb_raw_data,
                                                      ukb_raw_data_plateid,
                                                      ukb_raw_data_plateid_ii,
                                                      ukb_assay_data),
  
  ukb_folds_data = make_repeated_folds(ukb_participant_data,
                                       n_folds = 10,
                                       n_reps = 5),
  
  tar_group_count(ukb_folds_and_reps, distinct(ukb_folds_data, rep, test_fold), 5),
  
  
  protein_score_weights_all_olink = create_weights_for_protein_score(selected_model_conventional_all_olink_whole_sample,
                                                                     olink_meta_data,
                                                                     somascan_meta_data),
  protein_score_weights_all_somascan = create_weights_for_protein_score(selected_model_conventional_all_somascan_whole_sample,
                                                                        olink_meta_data,
                                                                        somascan_meta_data,
                                                                        convert_somascan_to_olink = TRUE),
  protein_score_weights_all_both = create_weights_for_protein_score(selected_model_conventional_all_both_whole_sample,
                                                                    olink_meta_data,
                                                                    somascan_meta_data,
                                                                    convert_somascan_to_olink = TRUE),
  protein_score_weights_sig_olink = create_weights_for_protein_score(selected_model_conventional_sig_olink_whole_sample,
                                                                     olink_meta_data,
                                                                     somascan_meta_data),
  protein_score_weights_sig_somascan = create_weights_for_protein_score(selected_model_conventional_sig_somascan_whole_sample,
                                                                        olink_meta_data,
                                                                        somascan_meta_data,
                                                                        convert_somascan_to_olink = TRUE),
  protein_score_weights_sig_both = create_weights_for_protein_score(selected_model_conventional_sig_both_whole_sample,
                                                                    olink_meta_data,
                                                                    somascan_meta_data,
                                                                    convert_somascan_to_olink = TRUE),
  protein_score_weights_boruta_olink = create_weights_for_protein_score(selected_model_conventional_boruta_olink_whole_sample,
                                                                        olink_meta_data,
                                                                        somascan_meta_data),
  protein_score_weights_boruta_somascan = create_weights_for_protein_score(selected_model_conventional_boruta_somascan_whole_sample,
                                                                           olink_meta_data,
                                                                           somascan_meta_data,
                                                                           convert_somascan_to_olink = TRUE),
  protein_score_weights_boruta_both = create_weights_for_protein_score(selected_model_conventional_boruta_both_whole_sample,
                                                                       olink_meta_data,
                                                                       somascan_meta_data,
                                                                       convert_somascan_to_olink = TRUE),
  
  
  ukb_participant_data_conventional = select(ukb_participant_data,
                                             any_of(c("csid",
                                                      "status",
                                                      "sex",
                                                      "age_at_study_date_c",
                                                      "smoking_category",
                                                      "has_diabetes",
                                                      "sbp_mean_c",
                                                      "waist_mm_c"))),
  
  tar_target(trained_model_ukb_conventional,
             mutate(ukb_folds_and_reps,
                    model = purrr::pmap(list(rep, test_fold),
                                        \(rep, test_fold) {
                                          training_data <- create_training_data(rep, test_fold, ukb_folds_data,
                                                                                ukb_participant_data_conventional,
                                                                                NULL)
                                          train_lasso_logistic(training_data,
                                                               c("sex",
                                                                 "age_at_study_date_c",
                                                                 "smoking_category",
                                                                 "has_diabetes",
                                                                 "sbp_mean_c",
                                                                 "waist_mm_c"))
                                        })),
             pattern = map(ukb_folds_and_reps)),
  tar_target(predictions_ukb_conventional,
             mutate(trained_model_ukb_conventional,
                    preds = purrr::pmap(list(model, rep, test_fold),
                                        \(model, rep, test_fold) {
                                          test_data <- create_test_data(rep, test_fold, ukb_folds_data,
                                                                        ukb_participant_data_conventional, NULL)
                                          get_predictions(test_data,
                                                          model,
                                                          c("sex",
                                                            "age_at_study_date_c",
                                                            "smoking_category",
                                                            "has_diabetes",
                                                            "sbp_mean_c",
                                                            "waist_mm_c"))
                                        })) |>
               select(rep, test_fold, preds),
             pattern = map(trained_model_ukb_conventional)),
  tar_target(c_statistic_ukb_conventional, get_c_statistic(predictions_ukb_conventional, predictions_ukb_conventional, ukb_participant_data)),
  
  
  tar_target(trained_model_ukb_age_sex,
             mutate(ukb_folds_and_reps,
                    model = purrr::pmap(list(rep, test_fold),
                                        \(rep, test_fold) {
                                          training_data <- create_training_data(rep, test_fold, ukb_folds_data,
                                                                                ukb_participant_data_conventional,
                                                                                NULL)
                                          train_lasso_logistic(training_data,
                                                               c("sex",
                                                                 "age_at_study_date_c"))
                                        })),
             pattern = map(ukb_folds_and_reps)),
  tar_target(predictions_ukb_age_sex,
             mutate(trained_model_ukb_age_sex,
                    preds = purrr::pmap(list(model, rep, test_fold),
                                        \(model, rep, test_fold) {
                                          test_data <- create_test_data(rep, test_fold, ukb_folds_data,
                                                                        ukb_participant_data_conventional, NULL)
                                          get_predictions(test_data,
                                                          model,
                                                          c("sex",
                                                            "age_at_study_date_c"))
                                        })) |>
               select(rep, test_fold, preds),
             pattern = map(trained_model_ukb_age_sex)),
  tar_target(c_statistic_ukb_age_sex, get_c_statistic(predictions_ukb_age_sex, predictions_ukb_conventional, ukb_participant_data)),
  
  
  tar_map(
    values = tibble::tibble(
      name = c("all_olink",
               "all_somascan",
               "all_both",
               "sig_olink",
               "sig_somascan",
               "sig_both",
               "boruta_olink",
               "boruta_somascan",
               "boruta_both")) |>
      dplyr::mutate(protein_weights = rlang::syms(paste0("protein_score_weights_", name))),
    names = "name",
    
    tar_target(
      ukb_participant_data_with_score,
      {
        data <- select(ukb_participant_data, !starts_with("oid") | all_of(protein_weights$term))
        data$protein_risk_score <- as.numeric(as.matrix(data[,protein_weights$term]) %*% protein_weights$coef)
        data
      }
    ),
    
    ## cross-validation
    tar_target(trained_model_ukb,
               mutate(ukb_folds_and_reps,
                      model = purrr::pmap(list(rep, test_fold),
                                          \(rep, test_fold) {
                                            
                                            
                                            training_data <- create_training_data(rep, test_fold, ukb_folds_data,
                                                                                  ukb_participant_data_with_score,
                                                                                  NULL)
                                            train_lasso_logistic(training_data,
                                                                 c("sex",
                                                                   "age_at_study_date_c",
                                                                   "smoking_category",
                                                                   "has_diabetes",
                                                                   "sbp_mean_c",
                                                                   "waist_mm_c",
                                                                   "protein_risk_score"))
                                          })),
               pattern = map(ukb_folds_and_reps)),
    tar_target(predictions_ukb,
               mutate(trained_model_ukb,
                      preds = purrr::pmap(list(model, rep, test_fold),
                                          \(model, rep, test_fold) {
                                            test_data <- create_test_data(rep, test_fold, ukb_folds_data,
                                                                          ukb_participant_data_with_score, NULL)
                                            get_predictions(test_data,
                                                            model,
                                                            c("sex",
                                                              "age_at_study_date_c",
                                                              "smoking_category",
                                                              "has_diabetes",
                                                              "sbp_mean_c",
                                                              "waist_mm_c",
                                                              "protein_risk_score"))
                                          })) |>
                 select(rep, test_fold, preds),
               pattern = map(trained_model_ukb)),
    
    tar_target(c_statistic_ukb, get_c_statistic(predictions_ukb, predictions_ukb_conventional, ukb_participant_data)),
    tar_target(nri_ukb, get_nri(predictions_ukb, predictions_ukb_conventional, ukb_participant_data))
  )
)

