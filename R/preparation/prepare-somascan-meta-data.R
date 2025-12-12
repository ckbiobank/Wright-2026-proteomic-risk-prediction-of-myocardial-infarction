prepare_somascan_meta_data <- function(meta){
  meta <- meta |> 
    select(protein_id = aptname, seqid, somaid,
           target_full_name, target, uniprot,
           organism, type, dilution2)
  
  return(meta)
}
