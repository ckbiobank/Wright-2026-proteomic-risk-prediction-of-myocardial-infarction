
create_weights_for_protein_score <- function(model,
                                             olink_meta_data,
                                             somascan_meta_data,
                                             convert_somascan_to_olink = FALSE){
  
  if (convert_somascan_to_olink){
    weights <- get_coefs(model) |>
      filter(grepl("^(oid)|(seq)", term),
             coef != 0) |> 
      left_join(select(somascan_meta_data, protein_id, uniprot), by = c("term" = "protein_id")) |> 
      left_join(select(olink_meta_data, uniprot, protein_id), by = "uniprot") |>
      mutate(protein_id = coalesce(protein_id, term)) |> 
      rename(protein_id = term, term = protein_id) |> 
      filter(grepl("^oid", term)) |> 
      group_by(term) |>
      summarise(coef = sum(coef),
                protein_id = paste(protein_id, collapse = ",")) |> 
      ungroup()
    
    stopifnot(length(unique(weights$term)) == length(weights$term))
    
    return(weights)
  }
  
  
  weights <- get_coefs(model) |>  
    filter(grepl("^(oid)|(seq)", term),
           coef != 0)
  
  return(weights)
}

