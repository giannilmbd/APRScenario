#' Exported version of full_scenarios_core
#' @inheritParams scenarios
#' @export
full_scenarios_core <- function(...) {
  .Call(`_APRScenario_full_scenarios_core`, ...)
}
