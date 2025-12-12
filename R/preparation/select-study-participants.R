select_study_participants <- function(participants,
                                      olink_data,
                                      somascan_data){
  ## keep participants with ascertainment and protein data
  participants <- filter(participants) |>
    filter(!is.na(ascertainment_group)) |> 
    filter(csid %in% unique(c(somascan_data$csid, olink_data$csid)))
  
  ## define follow up times (in days)
  ### for Prentice case-subcohort method
  ### follow up time starts at zero for participants in subcohort
  ### but starts just before endpoint for cases outside subcohort
  participants <- participants |> 
    mutate(endpoint_date = acute_mi_fatal_ihd_date,
           status        = acute_mi_fatal_ihd) |> 
    mutate(time_out = as.numeric(endpoint_date - study_date),
           time_in  = if_else(ascertainment_group != "subcohort",
                              time_out - 0.001,
                              0)) |> 
    filter(ascertainment_group == "subcohort" | status == 1) |> 
    mutate(group = case_when(
      ascertainment_group == "subcohort" & status == 0 ~ "subcohort_not_case",
      ascertainment_group == "subcohort" & status == 1 ~ "subcohort_case",
      ascertainment_group != "subcohort" & status == 1 ~ "case"
    ))
  
  return(participants)
}
