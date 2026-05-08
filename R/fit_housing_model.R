#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param housing_2000
#' @param housing_2015
#' @return
#' @author geryan
#' @export
fit_housing_model <- function(
    housing_2000,
    housing_2015
  ) {

  h2000 <- values(housing_2000)
  h2015 <- values(housing_2015)

  notna <- which(!is.na(h2000))

  hdat <- tibble(
    h0 = h2000[notna],
    h1 = h2015[notna],
    t = 15
  ) |>
    mutate(
      hdelta = h1 - h0
    )


  hmod <- nls(
    #formula = h1 ~ h0 + SSlogis(asym, 1, mid, t),
    formula = h1 ~ h0 + 1 /(1  + exp(-r *t)),
    start = list(r = 0.5),
    data = hdat
  )

  summary(hmod)


  preds <- predict(
    hmod,
    newdata = tibble(
      h0 = 0,
      t = seq(0:100)
    )
  )

  plot(
    x = seq(0:100),
    y = preds
  )

}
