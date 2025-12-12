
prepare_olink_data <- function(olink_explore,
                               olink_meta_data,
                               olink_plates){
  # create plateid column
  plates <- olink_plates |> 
    group_by(csid) |> 
    summarise(plateid = paste(unique(plateid), collapse = " ")) |> 
    mutate(plateid = factor(plateid))
  
  ## check 46 different plates
  stopifnot(length(levels(plates$plateid)) == 46)
  
  # prepare assay data
  olink_explore |> 
    left_join(select(olink_meta_data, assay, protein_id, panel_full),
              by = c("assay", "panel_full")) |> 
    filter(!is.na(protein_id)) |> 
    mutate(npx = if_else(qc_warning == 1, NA_real_, npx)) |> 
    select(csid, protein_id, npx) |> 
    pivot_wider(id_cols = "csid",
                names_from = "protein_id",
                values_from = "npx") |> 
    left_join(plates, by = "csid") |> 
    relocate(plateid, .after = 1)
}
