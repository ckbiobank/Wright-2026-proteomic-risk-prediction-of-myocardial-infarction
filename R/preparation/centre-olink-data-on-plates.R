# applies Olink's "Intensity Normalisation" but using medians of subcohort values

centre_olink_data_on_plates <- function(
    olink_data,
    casesubcohort_participant_data
){
  casesubcohort_participant_data |> 
    
    ## join Olink data
    inner_join(olink_data, by = "csid") |> 
    
    ## median centring of assay values (using subcohort),
    group_by(plateid) |>
    mutate(across(starts_with("oid") & !any_of(c("oid20775", "oid20816", "oid21485")),
                  \(x) x - median(x[ascertainment_group == "subcohort"],
                                  na.rm = TRUE))) |>
    ungroup() |> 
    mutate(across(any_of(c("oid20775", "oid20816", "oid21485")),
                  \(x) x - median(x[ascertainment_group == "subcohort"],
                                  na.rm = TRUE))) |>
    select(csid, plateid, starts_with("oid"))
}
