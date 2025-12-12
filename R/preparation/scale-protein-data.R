# centres and scales protein values (using subcohort mean and SD) 
scale_protein_data <- function(protein_data, participant_data){
  participant_data |> 
    select(csid, ascertainment_group) |> 
    inner_join(protein_data, by = "csid") |> 
    mutate(across(starts_with(c("oid", "seq")),
                  \(x) (x - mean(x[ascertainment_group == "subcohort"], na.rm = TRUE))/sd(x[ascertainment_group == "subcohort"], na.rm = TRUE))) |> 
    select(-ascertainment_group)
}
