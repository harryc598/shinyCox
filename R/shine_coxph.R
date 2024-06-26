#' Generates a shiny app for predictions from Cox model(s)
#'
#' Writes a shiny app to visualize predicted survival curves from one or
#' multiple Cox models. One feature of this function is that
#' the shiny app, once created, will not contain any identifiable data,
#' containing only information necessary for predictions.
#'
#' @param ... Arbitrary number of Cox proportional hazard models, created by
#'   [survival::coxph()] or [make_coxph()], which automatically ensures the
#'   models are appropriate for `shine_coxph()`
#' @param app.dir Directory where shiny app is created. Specifically, a
#'   sub-folder will be made containing the `app.R` file as well as the `.Rdata`
#'   file within `app.dir`. If no directory is provided, execution will pause
#'   and the user will be asked to confirm whether this sub-folder may be
#'   created in the working directory or to stop the function and provide an
#'   input `app.dir`.
#' @param theme Theme of shiny app.
#'  * `"default`: default theme, requires only shiny
#'  * `"dashboard"`: requires `"shinydashboard"` and `"DT"` packages
#'
#' @section Notes:
#'
#'   There are some requirements in order for this function to run without
#'   error: in your original [survival::coxph()] function or functions,
#'   `model = TRUE` and `x = TRUE` are required arguments (used to create the
#'   simplified `"coxph"` object). Currently, this function does not support
#'   penalized models (e.g., as created by `ridge()` and `pspline()`). Multiple
#'   strata terms and strata by covariate interaction terms in the formula are
#'   also not currently supported, but workarounds are available by respectively
#'   using a new strata factor variable encompassing all combinations of desired
#'   stratum variable levels. Use of time-varying covariates (e.g. with `tt()`)
#'   and multi-state models is not supported in our function. The package is not
#'   intended to support Fine-gray models by [survival::finegray()] creating Cox
#'   models, but doing so will not result in an error.
#'
#' @section Guidelines:
#'
#'   This package is intended to visualize and present predicted survival
#'   functions for fitted Cox models. In regards to formula notation, the
#'   variable names used are ultimately what will be displayed in the
#'   application. Using functions in the formula will work, but with multiple
#'   nested functions it will fail. Using "." notation is not currently
#'   supported. The `na.action` is inherited from the Cox models, with `omit`
#'   being the only option with support at this time. For these reasons, we
#'   recommend creating all final variables (including suitable transformations)
#'   with meaningful names prior to using [survival::coxph()].
#'
#'
#' @returns A list containing Cox model information along with the shiny app
#'  code. The app is written to the directory while the function is operating.
#' @examplesIf interactive()
#'
#'   library(survival)
#'
#'   # Data used is from survival package, renamed for legibility
#'   names(leukemia)[names(leukemia) == "x"] <- "treatment"
#'   # Make Cox model, with x = TRUE and model = TRUE
#'   model1 <- coxph(Surv(time, status) ~ treatment,
#'   leukemia, x = TRUE, model = TRUE)
#'
#'   # Use shine_coxph() to create shiny app in temporary directory
#'   shine_coxph("Model 1" = model1)
#'
#'   # Get directory for shiny app (should be first, check file list if not)
#'   filedir <- list.files(tempdir())[1]
#'
#'   # Run shiny app from temporary directory
#'   shiny::runApp(paste0(tempdir(), "/", filedir))
#'   # Remove app from directory once finished
#'   unlink(paste0(tempdir(),"/",filedir), recursive = TRUE)
#'
#'
#' @export
#' @importFrom grDevices hcl.colors
#' @importFrom shiny showTab
#' @importFrom utils menu
shine_coxph <- function(..., app.dir = NULL, theme = c("default", "dashboard"))

{
  ########################
  # determine the class of each input argument
  input.list <- list(...)
  # strata by covariate error
  if(any(grepl(":strata", names(input.list[[1]]$coefficients), fixed = TRUE))) {
    stop("Your model containes strata by covariate interaction, which is not currently supported")
  }
  if(inherits(input.list[[1]], "coxph.penal")) {
    stop("shine_coxph does not currently support penalized Cox models")
  } else if(inherits(input.list[[1]], "coxms")) {
    stop("shine_coxph does not support multi-state models")
  } else
  n.list <- length(input.list)
  list.class <- rep("",n.list)
  for (i in 1:n.list) {
    list.class[i] <- class(input.list[[i]])
  }

  #######################
  # extract the coxph models from the input
  cox.list <- which(list.class == "coxph")
  n.model <- length(cox.list)
  model.list <- vector("list",n.model)
  for (i in 1:n.model) {
    model.list[[i]] <- input.list[[cox.list[i]]]
  }
  names(model.list) <- names(input.list)[cox.list]


  ##########################
  # get app directory
  if (is.null(app.dir)) {
    check <- utils::menu(c("Yes", "No, I'll provide a directory"),
                  title = "You have not specified an app directory. Would you like to use tempdir()?")
    if(check == 1) {
    app.dir <- tempdir()
    } else {
      return(invisible())
    }
  }
  date.time <- as.character(Sys.time())
  date.time <- gsub(":", "-", date.time, fixed=TRUE)
  date.time <- gsub(" ", "-", date.time, fixed=TRUE)
  app.dir <- paste0(app.dir, "/", date.time, "/")
  dir.create(app.dir)
  message(paste0("Shiny app will be written in directory ", app.dir, "."))

  #######################
  # convert the coxph models to coxfit objects
  coxfit.list <- vector("list",n.model)
  part1.out = vector("list", n.model)
  for (i in 1:n.model) {
    coxfit.list[[i]] <- prep_coxfit(model.list[[i]])
    part1.out[[i]] = surv_pred_info(model.list[[i]])
  }
  names(coxfit.list) <- names(model.list)

  #######################
  # Save app data into app directory
  cox.fit.list <- coxfit.list
  save(cox.fit.list = cox.fit.list, part1.out = part1.out, file = paste0(app.dir, "appData.Rdata"))

  #######################
  # get different code sections
  input.data.code <- write_coxfit_input_data_code(coxfit.list)
  KM.plot.code <- write_KM_plot_code(coxfit.list)
  theme <- match.arg(theme)
  if (theme == "default") {
  table.code <- prop_haz_tables(cox.fit.list)
  } else if (theme == "dashboard") {
     table.code <- make_DT_table(cox.fit.list)
   } else stop("Invalid theme, please choose 'default' or 'dashboard' as a theme")

  ########################
  # generate header code
  hdr.code <- c("#rm(list=ls())",
             "options(stringsAsFactors=FALSE)",
             "library(shinyCox)",
             paste0("load('", gsub("\\\\", "/", app.dir), "appData.Rdata')"))


  ########################
  # generate ui code
  if (theme == "default") {
  ui.code <- c("ui=navbarPage(",
               " # Change app title with the argument below",
            "             'Cox Model Survival Predictions',",
            "    tabPanel('Plot',",
            "             sidebarLayout(",
            "             sidebarPanel(",
            " # Inputs below; for any sliderInput with binary options, e.g. sex,",
            " # add the argument 'step = 1'",
            paste0("                          ",
                   input.data.code$ui.code,
                   ","),
            "textInput('predProbTimes','Times for predicted probabilities',placeholder='Enter values separated by a comma'),",
            "checkboxInput('confint', 'Confidence Intervals', value = FALSE),",
            "textInput('ylab', 'Y-axis Text'),",
            "actionButton(inputId = 'go', label = 'Generate Plot'),",
            "actionButton(inputId = 'reset', label = 'Reset'),",
            "actionButton(inputId = 'app.exit', label = 'Exit App'),",
            "selectInput('clrs', label = 'Choose Colors', choices = palette.pals(), selected = 'Okabe-Ito'),",
            "hr(),",
            "numericInput('height', 'Plot Height', value = 400),",
            "numericInput('width', 'Plot Width', value = 800),",
            "downloadButton('downloadPlot', 'Download Plot')",
            "),",
            "mainPanel(",
            "h3('Predicted Survival Curve'),",
            "plotOutput(outputId = 'KM'),",
            "h3('Predicted Probability at Fixed Times'),",
            "textOutput(outputId='noPredTimes'),",
            "tableOutput(outputId = 'cox.times')))),",
            "tabPanel('Summary Tables',",
            table.code$ui.code,
            "))")
  } else if (theme == "dashboard") {
    ui.code <- c("library(shinydashboard)",
                     "ui = dashboardPage(",
                     " # Change title of app with dashboardHeader argument",
                     "           dashboardHeader(title = 'Cox Model Survival Predictions'),",
                     "           dashboardSidebar(",
                     "             sidebarMenu(",
                     "             menuItem('Plot', tabName = 'plot'),",
                     table.code$uitop,
                     "             )",
                     "           ),",
                     "     dashboardBody(",
                     "       tabItems(",
                     "     tabItem(tabName = 'plot',",
                     "             sidebarLayout(",
                     "             sidebarPanel(",
                     " # Inputs below; for any sliderInput with binary options, e.g. sex,",
                     " # add the argument 'step = 1'",
                     paste0("                          ",
                            input.data.code$ui.code,
                            ","),
                     "textInput('predProbTimes','Times for predicted probabilities',placeholder='Enter values separated by a comma'),",
                     "checkboxInput('confint', 'Confidence Intervals', value = FALSE),",
                     "textInput('ylab', 'Y-axis Text'),",
                     "actionButton(inputId = 'go', label = 'Generate Plot'),",
                     "actionButton(inputId = 'reset', label = 'Reset'),",
                     "actionButton(inputId = 'app.exit', label = 'Exit App'),",
                     "selectInput('clrs', label = 'Choose Colors', choices = palette.pals(), selected = 'Okabe-Ito'),",
                     "hr(),",
                     "numericInput('height', 'Plot Height', value = 400),",
                     "numericInput('width', 'Plot Width', value = 800),",
                     "downloadButton('downloadPlot', 'Download Plot')",
                     "),",
                     "mainPanel(",
                     "h3('Predicted Survival Curve'),",
                     "plotOutput(outputId = 'KM'),",
                     "h3('Predicted Probability at Fixed Times'),",
                     "textOutput(outputId='noPredTimes'),",
                     "tableOutput(outputId = 'cox.times')))),",
                     table.code$uibottom,
                     "))))"
    )
  }

  #########################
  # generate server code
  server.code <- c("server = function(input,output)",
                "{",
                "          observeEvent(input$app.exit, {stopApp()}) # Exit when exit button is pressed",
                table.code$server.code,
                "          observeEvent(input$go, {",
                input.data.code$server.code,
                KM.plot.code$server.code,
                "predProbTable <- cox_times_table(KM.hat,input$predProbTimes)",
                "if (is.null(predProbTable)) output$noPredTimes <- renderText('No input times detected. If you provided times, check that you separated numbers with a single comma and you provided valid numbers.') else output$noPredTimes <- renderText(invisible())",
                "output$cox.times <- renderTable(predProbTable, rownames = TRUE)",
                "colors <- palette.colors(length(cox.fit.list), input$clrs)",
                "output$downloadPlot <- downloadHandler(",
                "  filename = function() { paste0('plot.png') },",
                "  content = function(file) {",
                "    png(file, input$width, input$height)",
                "    shinyCox::cox_KM_plots(KM.hat,clrs=colors, confint=input$confint, ylab=input$ylab)",
                "    dev.off()",
                "  })",
                "})",
                "          observeEvent(input$reset, {output$KM <- output$HR <- output$PHA <- output$cox.times <- NULL}) # Reset main area",
                "}")

  #--------------------------------------------
  # app code
  app.code <- c("cox.app = shinyApp(ui,server)")


  #--------------------------------------------
  # write the app file
  app.code.file <- paste0(app.dir, "app.R")
  write(unlist(hdr.code), file = app.code.file, sep = "\n")
  write(unlist(ui.code), file = app.code.file, sep = "\n", append = TRUE)
  write(unlist(server.code), file = app.code.file, sep = "\n", append = TRUE)
  write(unlist(app.code), file = app.code.file, sep = "\n", append = TRUE)



  #--------------------------------------------
  # return result

  res=list(app.dir = app.dir,
           app.data.file = paste0(app.dir, "/appData.Rdata"),
           app.code.file = app.code.file,
           coxfit.list = coxfit.list,
           hdr.code = hdr.code,
           ui.code = ui.code,
           server.code = server.code)

  return(res)

}
