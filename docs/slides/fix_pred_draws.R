fix_draws = function(object, newdata, ..., func = tidybayes::predicted_draws) {
  draws = func(object, newdata, ...)

  n = names(draws)

  dplyr::full_join(
    draws |> dplyr::select(-.chain, -.iteration),
    tidybayes::tidy_draws(object) |>
      dplyr::select(.chain, .iteration, .draw),
    by = ".draw"
  ) |>
    dplyr::select(dplyr::all_of(n)) %>%
    dplyr::ungroup()
}

predicted_draws_fix = function(object, newdata, ...) {
  fix_draws(object, newdata, ..., func = tidybayes::predicted_draws)
}

epred_draws_fix = function(object, newdata, ...) {
  fix_draws(object, newdata, ..., func = tidybayes::epred_draws)
}

residual_draws_fix = function(object, newdata, ...) {
  fix_draws(object, newdata, ..., func = tidybayes::residual_draws)
}