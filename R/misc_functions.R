

#' @returns UI and server code for shinydashboard and DT
#' @noRd
make_DT_table <- function(cox.fit.list) {
  uicodetop <- c("menuItem('Summary Tables',",
                 "         tabName = 'tables', icon = icon('table'),")
  uicodebottom <- c()

  servercode <- c()

  for(i in 1:length(cox.fit.list)) {
  if(i < length(cox.fit.list)) {
    uicodetop <- c(uicodetop,
                   paste0("menuSubItem('", names(cox.fit.list)[i], "', tabName = 'table", i, "'),")
    )
    uicodebottom <- c(uicodebottom,
                      paste0("tabItem('table", i, "',"),
                      "fluidRow(",
                      paste0("infoBoxOutput('subjects", i, "'),"),
                      paste0("infoBoxOutput('events", i, "')),"),
                      "h4('Hazard Ratio Summary Table'),",
                      paste0("DT::dataTableOutput(outputId = 'HR", i, "'),"),
                      "h4('Assessing the Proportional Hazards Assumption'),",
                      paste0("DT::dataTableOutput(outputId = 'PHA", i, "')"),
                      "),"
                      )
  }
    else {
    uicodetop <- c(uicodetop,
                   paste0("menuSubItem('", names(cox.fit.list)[i], "', tabName = 'table", i, "'))")
    )

    uicodebottom <- c(uicodebottom,
                      paste0("tabItem('table", i, "',"),
                      "fluidRow(",
                      paste0("infoBoxOutput('subjects", i, "'),"),
                      paste0("infoBoxOutput('events", i, "')),"),
                      "h4('Hazard Ratio Summary Table'),",
                      paste0("DT::dataTableOutput(outputId = 'HR", i, "'),"),
                      "h4('Assessing the Proportional Hazards Assumption'),",
                      paste0("DT::dataTableOutput(outputId = 'PHA", i, "')")
                   )

    }
    servercode <- c(servercode,
                    paste0("output$HR", i, "=DT::renderDataTable(DT::datatable(cox.fit.list[[", i, "]]$HR.table,"),
                           "                                     options = list(",
                           "                                     dom = 't'",
                           "                                     )) |>",
                           "                 # The formatRound() function is set to give four digits after the decimal, change 'digits' to alter this",
                           "                 DT::formatRound(columns = c('Hazard Ratio', 'Lower Bound', 'Upper Bound', 'p value'), digits = 4) |>",
                           "                 DT::formatStyle('p value',",
                           "                                 target = 'cell',",
                           "                                 # the fontweight argument will bold cells under a certain value, the default is 0.05",
                           "                                 fontweight = DT::styleInterval(0.05, c('bold', 'normal'))))",
                    paste0("output$PHA", i, "=DT::renderDataTable(DT::datatable(cox.fit.list[[", i, "]]$PHA.table$table,"),
                           "                                      options = list(",
                           "                                      dom = 't')))",
                           paste0("output$subjects", i, "=renderInfoBox({"),
                           paste0("infoBox('Subjects', cox.fit.list[[", i, "]]$nsample,"),
                           "       color = 'red', icon = icon('hospital-user'))})",
                           paste0("output$events", i, "=renderInfoBox({"),
                           paste0("infoBox('Events', cox.fit.list[[", i, "]]$nevents,"),
                    "       color = 'red', icon = icon('file-medical'))})"
                    )

  }
  code.res <- list(uitop = uicodetop,
                   uibottom = uicodebottom,
                   server.code = servercode)
  return(code.res)
}

#' Wrapper to create `survival::coxph()` object suitable for [shine_coxph()]
#'
#' Performs [survival::coxph()] with `model = TRUE` and `x = TRUE` as defaults.
#' Checks that Cox model is appropriate for use with [shine_coxph()].
#'
#' @param formula a formula object, with the response on the left of a `~`
#'   operator, and the terms on the right.  The response must be a survival
#'   object as returned by the `Surv` function.
#' @inheritParams survival::coxph
#' @param ... other arguments which will be passed to `coxph()`. Note that
#'  `x = TRUE` and `model = TRUE` are the default arguments (and required by
#'  [shine_coxph()]), you do not need to include them here.
#' @returns Object of class `"coxph"` representing the fit
#'
#' @examples
#' library(survival)
#' ovarianph <- make_coxph(Surv(futime, fustat) ~ age + strata(rx),
#' data = ovarian)
#' @export
make_coxph <- function(formula, data, ...) {
  Call <- match.call()
  if(missing(formula)) {
    stop("a formula argument is required")
  }
  if(missing(data)) {
    stop("a data argument is required")
  }
  Call$model <- TRUE
  Call$x <- TRUE
  Call[[1]] <- quote(coxph)
  cox_model <- eval(Call)
  if(!is.null(cox_model$call$tt)) {
    stop("tt terms cannot be used in the shine_cox() function")
  }
  if(inherits(cox_model, "coxphms")) {
    stop("This appears to be a multistate model, which is not appropiate for the shine_cox function.")
  }
  if(inherits(cox_model, "coxph.penal")) {
    stop("This appears to be a penalized Cox model, which is not currently
            supported for use with shine_coxph()")
  }
  if("fgwt" %in% names(data)) {
    warning("Your dataset appears to be made using the finegray() function. This is not an appropiate model to use with shine_coxph()")
  }
  return(cox_model)
}
