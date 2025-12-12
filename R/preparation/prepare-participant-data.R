
prepare_participant_data <- function(baseline,
                                     ascertainments,
                                     endpoints,
                                     endpoint_definitions){
  ascertainments <- ascertainments |> 
    mutate(ascertainment_group = case_when(
      olinkexploreexpansion_chd_b1_subcohort == 1  ~ "subcohort",
      olinkexploreexpansion_chd_b1_case_mi == 1    ~ "mi_case",
      olinkexploreexpansion_chd_b1_case_oami == 1  ~ "mi_case",
      olinkexploreexpansion_chd_b1_case_fihd  == 1 ~ "fihd_case")) |> 
    select(csid, ascertainment_group)
  
  baseline <- baseline |> 
    select(csid,
           dob_anon,
           study_date,
           censoring_date,
           region_code,
           is_female,
           age_at_study_date_x100,
           highest_education,
           smoking_category,
           alcohol_category,
           met,
           bmi_calc,
           waist_mm,
           has_diabetes,
           hours_since_last_ate_x10,
           sbp_mean,
           mean_temp,
           chd_diag,
           stroke_or_tia_diag,
           rheum_heart_dis_diag
    ) |> 
    mutate_at(vars(study_date, censoring_date), as.Date) |>
    mutate(age_at_study_date = age_at_study_date_x100 / 100,
           hours_since_last_ate = hours_since_last_ate_x10 / 10,
           
           sex = factor(is_female,
                        levels = attr(.data$is_female, "labels"),
                        labels = names(attr(.data$is_female, "labels"))),
           region = factor(region_code,
                           levels = attr(.data$region_code, "labels"),
                           labels = names(attr(.data$region_code, "labels"))),
           highest_education = factor(highest_education,
                                      levels = attr(.data$highest_education, "labels"),
                                      labels = names(attr(.data$highest_education, "labels"))),
           smoking_category = factor(smoking_category,
                                     levels = attr(.data$smoking_category, "labels"),
                                     labels = names(attr(.data$smoking_category, "labels"))),
           alcohol_category = factor(alcohol_category,
                                     levels = attr(.data$alcohol_category, "labels"),
                                     labels = names(attr(.data$alcohol_category, "labels"))),
           has_diabetes = factor(has_diabetes,
                                 levels = attr(.data$has_diabetes, "labels"),
                                 labels = names(attr(.data$has_diabetes, "labels")))
    ) |> 
    mutate(
      age_at_study_date_c = scale(age_at_study_date, center = TRUE, scale = FALSE)[,1],
      sbp_mean_c = scale(sbp_mean, center = TRUE, scale = FALSE)[,1],
      waist_mm_c = scale(waist_mm, center = TRUE, scale = FALSE)[,1]
    )
  
  
  ## endpoints
  endpoint_definitions <- endpoint_definitions |> 
    distinct(diy_epid, description) |> 
    mutate(name = case_match(description,
                             "CKB0003 - Ischaemic heart disease (ICD-10: I20-I25)" 
                             ~ "ihd",
                             "CKB0004 - Acute myocardial infarction (ICD-10: I21)" 
                             ~ "acute_mi",
                             "CKB0068 - Any myocardial infarction (ICD-10: I21-I23)" 
                             ~ "any_mi",
                             "CKB0005 - Other ischaemic heart disease (endpoint 0003 but not 0004)"
                             ~ "other_ihd")) |> 
    filter(!is.na(name))
  
  endpoints <- endpoints |> 
    rename_with(\(x) stringi::stri_replace_all_regex(x,
                                                     pattern = paste0("ep_", endpoint_definitions$diy_epid),
                                                     replacement = endpoint_definitions$name,
                                                     vectorize = FALSE)) |> 
    rename_with( ~ gsub("_ep", "", .x, fixed = TRUE)) |> 
    rename_with( ~ gsub("_datedeveloped", "_date", .x, fixed = TRUE)) |> 
    mutate_at(vars(ends_with("_date")), as.Date) |> 
    select(csid,
           ihd            = ihd_combined,
           ihd_date       = ihd_combined_date,
           acute_mi       = acute_mi_combined,
           acute_mi_date  = acute_mi_combined_date,
           fatal_ihd      = ihd_du,
           fatal_ihd_date = ihd_du_date,
           other_ihd      = other_ihd_combined,
           other_ihd_date = other_ihd_combined_date) |> 
    mutate(acute_mi_fatal_ihd      = pmax(acute_mi, fatal_ihd),
           acute_mi_fatal_ihd_date = pmin(acute_mi_date, fatal_ihd_date))
  
  ## join baseline and ascertainment
  participants <- full_join(baseline, ascertainments, by = "csid") |> 
    full_join(endpoints, by = "csid")
  
  ## participants ascertained as cases, but with no acute MI or fatal IHD endpoint
  ## not to be counted as "cases"
  participants <- participants |> 
    mutate(ascertainment_group = if_else(
      ascertainment_group %in% c("mi_case", "fihd_case") &
        acute_mi_fatal_ihd == 0,
      NA_character_,
      ascertainment_group
    ))
  
  
  ## in subcohort, censor MI + fatal IHD at other IHD
  participants <- participants |>  
    mutate(
      acute_mi_fatal_ihd = case_when(
        ascertainment_group == "subcohort" & acute_mi_fatal_ihd == 1 & acute_mi_fatal_ihd_date <= other_ihd_date
        ~ 1,
        ascertainment_group == "subcohort" ~ 0,
        TRUE ~ acute_mi_fatal_ihd
      ),
      acute_mi_fatal_ihd_date = case_when(
        ascertainment_group == "subcohort" ~ pmin(acute_mi_fatal_ihd_date, other_ihd_date),
        TRUE ~ acute_mi_fatal_ihd_date
      )
    )
  
  return(participants)
}
