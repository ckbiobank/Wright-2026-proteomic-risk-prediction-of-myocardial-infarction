
rint <- function(x){
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

prepare_ukb_olink_score <- function(gadd_protein_score_coefs,
                                    olink_data,
                                    casesubcohort_participant_data){
  weights <- gadd_protein_score_coefs |> 
    mutate(protein =stringr::str_to_lower(stringr::str_extract(protein, "(.*\\.)(OID.*)(\\..*)", group = 2)))
  
  q <- olink_data|> 
    inner_join(select(casesubcohort_participant_data, csid), by = "csid") |> 
    mutate(across(starts_with("oid"),
                  \(x) coalesce(x, mean(x, na.rm = TRUE)))) |> 
    mutate(across(starts_with("oid"), rint)) |> 
    select(csid, all_of(weights$protein)) |> 
    mutate(across(starts_with("oid"),
                  \(x) x * weights$coefficient[match(cur_column(), weights$protein)])) |> 
    mutate(score = rowSums(across(starts_with("oid")))) |> 
    select(csid, ukb_olink_score = score)
  
  return(q)
}
