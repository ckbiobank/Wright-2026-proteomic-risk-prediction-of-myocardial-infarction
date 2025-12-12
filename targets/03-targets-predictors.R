
targets_predictors <- tar_plan(
  
  predictors_age_sex = c("sex",
                         "age_at_study_date_c"),
  
  predictors_conventional = c("sex",
                              "age_at_study_date_c",
                              "smoking_category",
                              "has_diabetes",
                              "sbp_mean_c",
                              "waist_mm_c"),
  
  ## olink proteins
  predictors_all_olink = olink_meta_data$protein_id,
  predictors_sig_olink = olink_individual_proteins_analysis |> 
    mutate(predictors = map(results, \(x) select_sig_proteins(x, olink_meta_data))) |> 
    select(rep, test_fold, predictors),
  tar_target(
    predictors_boruta_olink,
    predictors_sig_olink |> 
      mutate(predictors = pmap(
        list(rep, test_fold, predictors),
        \(x, y, z) get_boruta_selected(z,
                                       create_training_data(x, y, folds_data,
                                                            casesubcohort_participant_data,
                                                            centre_olink_data_on_plates(olink_data, casesubcohort_participant_data) ))
      )),
    pattern = map(predictors_sig_olink)
  ),
  
  ## somascan proteins
  predictors_all_somascan = somascan_meta_data |> 
    filter(organism == "Human", type == "Protein") |> 
    pull(protein_id),
  predictors_sig_somascan = somascan_individual_proteins_analysis |> 
    mutate(predictors = map(results, \(x) select_sig_proteins(x, somascan_meta_data))) |> 
    select(rep, test_fold, predictors),
  tar_target(
    predictors_boruta_somascan,
    predictors_sig_somascan |> 
      mutate(predictors = pmap(
        list(rep, test_fold, predictors),
        \(x, y, z) get_boruta_selected(z,
                                       create_training_data(x, y, folds_data,
                                                            casesubcohort_participant_data,
                                                            somascan_data))
      )),
    pattern = map(predictors_sig_somascan)
  ),
  
  ## olink and somascan
  predictors_all_both = c(predictors_all_olink, predictors_all_somascan),
  predictors_sig_both = predictors_sig_olink |> 
    rename(olink_predictors = predictors) |> 
    left_join(rename(predictors_sig_somascan, somascan_predictors = predictors), by = c("rep", "test_fold")) |> 
    mutate(predictors = pmap(list(olink_predictors, somascan_predictors), \(x, y) c(x, y))) |> 
    select(rep, test_fold, predictors),
  predictors_boruta_both = predictors_boruta_olink |> 
    rename(olink_predictors = predictors) |> 
    left_join(rename(predictors_boruta_somascan, somascan_predictors = predictors), by = c("rep", "test_fold")) |> 
    mutate(predictors = pmap(list(olink_predictors, somascan_predictors), \(x, y) c(x, y))) |> 
    select(rep, test_fold, predictors),
  
  ## with conventional
  
  predictors_conventional_all_olink = c(predictors_conventional, predictors_all_olink), 
  predictors_conventional_sig_olink = mutate(predictors_sig_olink,
                                             predictors = map(predictors, \(x) c(predictors_conventional, x))),
  predictors_conventional_boruta_olink = mutate(predictors_boruta_olink,
                                                predictors = map(predictors, \(x) c(predictors_conventional, x))),
  
  predictors_conventional_all_somascan = c(predictors_conventional, predictors_all_somascan), 
  predictors_conventional_sig_somascan = mutate(predictors_sig_somascan,
                                                predictors = map(predictors, \(x) c(predictors_conventional, x))),
  predictors_conventional_boruta_somascan = mutate(predictors_boruta_somascan,
                                                   predictors = map(predictors, \(x) c(predictors_conventional, x))),
  
  predictors_conventional_sig_both = mutate(predictors_sig_both,
                                            predictors = map(predictors, \(x) c(predictors_conventional, x))),
  
  predictors_conventional_boruta_both = mutate(predictors_boruta_both,
                                               predictors = map(predictors, \(x) c(predictors_conventional, x))),
  
  
  ## UKB Olink score
  predictors_ukb_olink_score = "ukb_olink_score",
  predictors_conventional_ukb_olink_score = c(predictors_conventional, predictors_ukb_olink_score),
  
  
  ## whole sample as training data
  predictors_sig_olink_whole_sample = select_sig_proteins(olink_individual_proteins_analysis_whole_sample, olink_meta_data),
  predictors_boruta_olink_whole_sample = get_boruta_selected(predictors_sig_olink_whole_sample,
                                                             casesubcohort_participant_data |> 
                                                               inner_join(centre_olink_data_on_plates(olink_data, casesubcohort_participant_data), by = "csid") |> 
                                                               mutate(across(starts_with("oid"), \(x) coalesce(x, mean(x, na.rm = TRUE))))),
  predictors_sig_somascan_whole_sample = select_sig_proteins(somascan_individual_proteins_analysis_whole_sample, somascan_meta_data),
  predictors_boruta_somascan_whole_sample = get_boruta_selected(predictors_sig_somascan_whole_sample,
                                                                casesubcohort_participant_data |> 
                                                                  inner_join(somascan_data, by = "csid") |> 
                                                                  mutate(across(starts_with("seq"), \(x) coalesce(x, mean(x, na.rm = TRUE))))),
  predictors_sig_both_whole_sample = c(predictors_sig_olink_whole_sample,
                                       predictors_sig_somascan_whole_sample),
  predictors_boruta_both_whole_sample = c(predictors_boruta_olink_whole_sample,
                                          predictors_boruta_somascan_whole_sample),
  predictors_conventional_sig_olink_whole_sample = c(predictors_conventional,
                                                     predictors_sig_olink_whole_sample),
  predictors_conventional_sig_somascan_whole_sample = c(predictors_conventional,
                                                        predictors_sig_somascan_whole_sample),
  predictors_conventional_sig_both_whole_sample = c(predictors_conventional,
                                                    predictors_sig_olink_whole_sample,
                                                    predictors_sig_somascan_whole_sample),
  predictors_conventional_boruta_olink_whole_sample = c(predictors_conventional,
                                                        predictors_boruta_olink_whole_sample),
  predictors_conventional_boruta_somascan_whole_sample = c(predictors_conventional,
                                                           predictors_boruta_somascan_whole_sample),
  predictors_conventional_boruta_both_whole_sample = c(predictors_conventional,
                                                       predictors_boruta_olink_whole_sample,
                                                       predictors_boruta_somascan_whole_sample)
)
