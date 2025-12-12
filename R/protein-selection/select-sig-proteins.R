select_sig_proteins <- function(result, meta_data){
  result |> 
    mutate(error_or_warning = !is.na(e) | !is.na(w)) |> 
    
    # add meta data
    inner_join(meta_data,
               by = "protein_id") |>
    
    # select one aptamer per protein target
    group_by(target) |> 
    arrange(p.value) |> 
    mutate(duplicate_target = row_number() > 1) |> 
    ungroup() |> 
    
    # identify human protein targets
    mutate(human_protein = (organism == "Human" & type == "Protein")) |> 
    
    # identify significant results
    mutate(p.value.orig = p.value,
           p.value = if_else(error_or_warning | duplicate_target | !human_protein,
                             NA,
                             p.value)) |> 
    mutate(p.value.bonferroni = p.adjust(p.value, method = "bonferroni", n = sum(!duplicate_target & human_protein)),
           p.value.fdr = p.adjust(p.value, method = "fdr", n = sum(!duplicate_target & human_protein))) |> 
    
    filter(p.value.fdr < 0.05) |>
    pull(protein_id)
}
