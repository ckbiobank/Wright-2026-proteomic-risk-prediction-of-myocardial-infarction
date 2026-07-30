get_nri <- function(predictions, reference_predictions, participant_data, ...){
  predictions |> 
    left_join(reference_predictions |> rename(preds_ref = preds), by = c("rep", "test_fold")) |> 
    mutate(nri = purrr::pmap(list(preds, preds_ref),
                             \(x, y) calculate_nri(x, y, participant_data))) |> 
    select(-preds, -preds_ref) |> 
    unnest(nri) |> 
    group_by(rep, test_fold) |> 
    summarise(nri = sum(nri),
              se = sqrt(sum(se^2)),
              .groups = "drop") |> 
    group_by(rep) |> 
    summarise(se = sd(nri)/sqrt(n()),
              nri = mean(nri),
              .groups = "drop") |>
    summarise(nri = mean(nri),
              se = sqrt(mean(se^2)))
}


# calculate continuous NRI
calculate_nri <- function(predictions, predictions_null, participant_data){
  if (!is.data.frame(predictions) | !is.data.frame(predictions_null)){return(tibble(statistic = c("nri.controls",
                                                                                                  "nri.cases"),
                                                                                    value = NA,
                                                                                    se = NA))}
  x <- select(predictions_null, csid, pred0 = linear_pred) |> 
    inner_join(select(predictions, csid, pred1 = linear_pred), by = "csid") |>
    inner_join(select(participant_data, csid, status), by = "csid") |> 
    mutate(pred0 = pred0,
           pred1 = pred1)
  
  x <- x |> 
    mutate(up = pred1 > pred0,
           down = pred1 < pred0,
           same = pred1 == pred0) |>
    group_by(status) |>
    summarise(n = n(),
              up = sum(up),
              down = sum(down)) |>
    mutate(pr_up = up / n,
           pr_down = down / n) |>
    mutate(nri = if_else(status == 1, pr_up - pr_down, pr_down - pr_up)) |>
    mutate(group = if_else(status == 0, "controls", "cases")) |> 
    select(group, nri, se)
  
  return(x)
}
