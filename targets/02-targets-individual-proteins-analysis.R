

targets_individual_proteins_analysis <- list(
  
  tar_target(folds_and_reps_ungrouped,
             distinct(folds_data, rep, test_fold)),
  
  ## olink
  tar_target(olink_individual_proteins_analysis,
             mutate(folds_and_reps_ungrouped,
                    results = purrr::pmap(list(rep, test_fold),
                                          \(rep, test_fold) {
                                            training_data <- create_training_data(rep, test_fold, folds_data,
                                                                                  casesubcohort_participant_data, olink_data,
                                                                                  .impute = FALSE)
                                            individual_proteins_analysis(training_data, olink_meta_data)
                                          })),
             pattern = map(folds_and_reps_ungrouped)),
  tar_target(somascan_individual_proteins_analysis,
             mutate(folds_and_reps_ungrouped,
                    results = purrr::pmap(list(rep, test_fold),
                                          \(rep, test_fold) {
                                            training_data <- create_training_data(rep, test_fold, folds_data,
                                                                                  casesubcohort_participant_data, somascan_data,
                                                                                  .impute = FALSE)
                                            individual_proteins_analysis(training_data, somascan_meta_data)
                                          })),
             pattern = map(folds_and_reps_ungrouped))
  
)
