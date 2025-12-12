targets_prep <- tar_plan(
  
  ## Import data ----
  tar_file_read(baseline_questionnaires,
                file.path(rawdatadir, "data_baseline_questionnaires.rds"),
                readRDS(!!.x)),
  tar_file_read(baseline_ascertainments,
                file.path(rawdatadir, "ascertainments.csv"),
                read_csv(!!.x)),
  tar_file_read(endpoints,
                file.path(rawdatadir, "endpoints.rds"),
                readRDS(!!.x)),
  tar_file_read(endpoint_definitions,
                file.path(rawdatadir, "endpoint_definitions.rds"),
                readRDS(!!.x)),
  tar_file_read(somalogic,
                file.path(rawdatadir, "data_baseline_somalogic.csv"),
                read_csv(!!.x)),
  tar_file_read(somalogic_meta,
                file.path(rawdatadir, "somalogic_meta.csv"),
                read_csv(!!.x)),
  tar_file_read(somalogic_samples,
                file.path(rawdatadir, "data_baseline_somalogic_samples.csv"),
                read_csv(!!.x)),
  tar_file_read(somalogic_noanml,
                file.path(rawdatadir, "data_baseline_somalogic_qc.csv"),
                read_csv(!!.x)),
  tar_file_read(somalogic_samples_noanml,
                file.path(rawdatadir, "data_baseline_somalogic_samples_qc.csv"),
                read_csv(!!.x)),
  tar_file_read(somalogic_prediction,
                file.path(rawdatadir, "data_baseline_somalogic_prediction.csv"),
                read_csv(!!.x)),
  
  tar_file_read(olink_assays,
                file.path(rawdatadir, "olink_assays.csv"),
                read_csv(!!.x)),
  tar_file_read(olink_explore,
                file.path(rawdatadir, "data_baseline_olink_explore.csv"),
                read_csv(!!.x)),
  tar_file_read(olink_plates,
                file.path(rawdatadir, "data_baseline_olink_plates.csv"),
                read_csv(!!.x)),
  
  
  
  # Prepare data ----
  
  all_participant_data = prepare_participant_data(
    baseline             = baseline_questionnaires,
    ascertainments       = baseline_ascertainments,
    endpoints            = endpoints,
    endpoint_definitions = endpoint_definitions),
  
  
  casesubcohort_participant_data = select_study_participants(
    all_participant_data,
    olink_data,
    somascan_data),
  
  
  ### Olink data ----
  olink_meta_data = prepare_olink_meta_data(olink_assays),
  olink_data = prepare_olink_data(olink_explore,
                                  olink_meta_data,
                                  olink_plates),
  
  
  ### for risk prediction using Olink:
  ### center protein values on plate
  ### keep only data for participants also with somascan data
  ### scale protein values
  olink_data_scaled = scale_protein_data(
    olink_data |> 
      centre_olink_data_on_plates(casesubcohort_participant_data) |> 
      inner_join(select(somascan_data, csid),by= "csid"),
    casesubcohort_participant_data),
  
  
  ### SomaScan data ----  
  somascan_meta_data = prepare_somascan_meta_data(somalogic_meta),
  somascan_data = prepare_somascan_data(somalogic,
                                        somalogic_samples,
                                        somalogic_prediction),
  
  somascan_noanml_data = prepare_somascan_data(somalogic_noanml,
                                               somalogic_samples_noanml,
                                               somalogic_prediction),
  
  ### for risk prediction using SomaScan:
  ### scale protein values
  somascan_data_scaled = scale_protein_data(somascan_data,
                                            casesubcohort_participant_data),
  somascan_noanml_data_scaled = scale_protein_data(somascan_noanml_data,
                                                   casesubcohort_participant_data),

    
  ### Olink and SomaScan data ----
  olink_somascan_data_scaled = inner_join(olink_data_scaled, somascan_data_scaled, by = "csid"),
  null_protein_data = select(olink_somascan_data_scaled, csid),

  ### Sampling weights ----
  casesubcohort_sampling_weights = make_casesubcohort_sampling_weights(all_participant_data,
                                                                       casesubcohort_participant_data),
  
  
  # External score ----
  tar_file_read(gadd_protein_score_coefs,
                file.path(rawdatadir2, "gadd_protein_score_coefs.csv"),
                read_csv(!!.x)),
  ukb_olink_score = prepare_ukb_olink_score(gadd_protein_score_coefs,
                                            olink_data_for_risk_prediction,
                                            casesubcohort_participant_data)
)
