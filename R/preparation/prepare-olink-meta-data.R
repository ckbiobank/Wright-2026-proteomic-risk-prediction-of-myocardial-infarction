
prepare_olink_meta_data <- function(olink_assays){
  olink_assays |>
    filter(batch %in% c("explore", "expand")) |> 
    select(protein_id = olinkid, panel_full, assay, uniprot) |>
    distinct() |> 
    group_by(assay) |> 
    arrange(panel_full) |> 
    filter(row_number() == 1) |> 
    ungroup() |> 
    mutate(organism = "Human", type = "Protein", target = assay) |> 
    mutate(protein_id = tolower(protein_id))
}
