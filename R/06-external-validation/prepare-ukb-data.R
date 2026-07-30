

prepare_ukb_assay_data <- function(ukb_raw_data, ukb_raw_data_ii, assay_meta_data){
  
  look_up_olinkid <- function(assay){
    assay_meta_data$protein_id[[match(assay, assay_meta_data$target)]]
  }
  
  ukb_assay_data <- ukb_raw_data |> 
    select(csid = eid, any_of(assay_meta_data$assay)) |> 
    rename_with(\(x) purrr::map_chr(x, look_up_olinkid), .cols = -csid)
  
  ukb_assay_data_ii <- ukb_raw_data_ii |> 
    rename_with(\(x) stringr::str_extract(x, "^NPX_(.*)_", group = 1), .cols = starts_with("NPX_")) |> 
    select(csid = eid, any_of(assay_meta_data$assay)) |> 
    rename_with(\(x) purrr::map_chr(x, look_up_olinkid), .cols = -csid)
  
  ukb_assay_data <- ukb_assay_data |> 
    full_join(ukb_assay_data_ii, by = "csid")
  
  return(ukb_assay_data)  
}



prepare_ukb_participant_data <- function(ukb_raw_data,
                                         ukb_raw_data_plateid,
                                         ukb_raw_data_plateid_ii,
                                         ukb_assay_data){
  
  ukb_plateid <- ukb_raw_data_plateid |> 
    mutate(csid = eid,
           plateid = factor(PlateID, labels = "a"),
           batch = Batch) |> 
    select(csid, plateid, batch) |> 
    bind_rows(
      ukb_raw_data_plateid_ii |> 
        mutate(csid = as.integer(PT_ID),
               plateid = factor(PlateID, labels = "b")) |> 
        select(csid, plateid)
    ) |> 
    group_by(csid) |> 
    summarise(plateid = paste(unique(plateid), collapse = " "),
              batch = first(batch)) |> 
    mutate(plateid = factor(plateid))
  
  ## check 632 different plates
  stopifnot(length(levels(ukb_plateid$plateid)) == 632)
  
  ukb_participant_data <- ukb_raw_data |> 
    select(csid           = eid,
           sex            = BaC_Sex,
           ethnicity      = Eth_Ethnicity,
           year_of_birth  = BaC_BirthYear,
           month_of_birth = BaC_BirthMonth,
           study_date     = Rec_DateAssess,
           lost_fu        = BaC_DateLostFU,
           smoking        = Smo_Status,
           waist_cm       = BSM_Waist,
           sbp_manual0    = BlP_SBPMan0,
           sbp_manual1    = BlP_SBPMan1,
           sbp_auto0      = BlP_SBPAuto0,
           sbp_auto1      = BlP_SBPAuto1,
           diabetes_diag  = HMH_Diabetes,
           blood_glucose  = BBC_GLU_Result,
           med_chol_bp_dm_0,
           med_chol_bp_dm_1,
           med_chol_bp_dm_2,
           mi_prev        = mi_prev,
           ihd_prev       = ihd_i20i25_prev,
           mi_inc,
           mi_inc_date) |> 
    mutate(
      dob               = make_date(year = year_of_birth, month = month_of_birth, day = 15L),
      age_at_study_date = interval(dob, study_date) / years(1),
      sex               = factor(sex,
                                 levels = c(0, 1),
                                 labels = c("Female", "Male")),
      white_british = as.numeric(ethnicity == 1001),
      smoking_category  = factor(smoking,
                                 levels = c(0, 1, 2, -3),
                                 labels = c("Never", "Previous", "Current", "Prefer not to answer")),
      has_diabetes = factor(coalesce(as.numeric(diabetes_diag == 1 |
                                                  med_chol_bp_dm_0 == 1 |
                                                  med_chol_bp_dm_1 == 1 |
                                                  med_chol_bp_dm_2 == 1 |
                                                  blood_glucose > 11.1), 0),
                            levels = c(0, 1),
                            labels = c("No", "Yes")),
      sbp_mean          = (coalesce(sbp_auto0, sbp_manual0) + coalesce(sbp_auto1, sbp_manual1))/2,
      waist_mm          = waist_cm * 10,
      ascertainment_group = "cohort",
      group               = "cohort",
      region              = "UK",
      endpoint_date = pmin(mi_inc_date, lost_fu, na.rm = TRUE)
    ) |> 
    
    mutate(
      status = as.numeric(mi_inc == 1 & mi_inc_date == endpoint_date),
      time_in = 0,
      time_out = interval(study_date, endpoint_date) / days(1),
    ) |> 
    
    left_join(ukb_plateid, by = "csid") |> 
    
    filter(time_out > time_in) |>
    
    filter(white_british == 1,
           mi_prev       == 0,
           ihd_prev      == 0,
           !batch %in% c(0, 7)) |>
    
    ## centre continuous predictor variables
    mutate(
      age_at_study_date_c = scale(age_at_study_date, center = TRUE, scale = FALSE)[,1],
      sbp_mean_c = scale(sbp_mean, center = TRUE, scale = FALSE)[,1],
      waist_mm_c = scale(waist_mm, center = TRUE, scale = FALSE)[,1]
    ) |> 
    
    select(csid, age_at_study_date_c, sex, smoking_category, sbp_mean_c, has_diabetes, waist_mm_c,
           plateid,
           time_in, time_out, status,
           ascertainment_group, group, region) |>
    na.omit() |> 
    
        inner_join(ukb_assay_data, by= "csid") |> 
    group_by(plateid) |>
    mutate(across(starts_with("oid"),
                  \(x) x - median(x, na.rm = TRUE))) |>
    ungroup() |> 
    
    mutate(across(starts_with("oid"),
                  \(x) x/sd(x, na.rm = TRUE))) |> 
    mutate(across(starts_with("oid"),
                  \(x) coalesce(x, mean(x, na.rm = TRUE))))
}

