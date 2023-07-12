#' Build a shiny app for coxph models
#'
#' Generates a shiny app to visualize Cox proportional hazard models created by
#' `coxph()`
#'
#' @param ... Any number of coxph models
#' @param app.dir directory where shiny app will be written
#' @section Details:
#'
#'   One of the important features of this function is that the shiny app, once
#'   created, will not contain any identifiable data, containing only
#'   information necessary for predictions. There are some requirements in order
#'   for this function to run without error: in your original `coxph()` function
#'   or functions, model = TRUE and x = TRUE are required arguments. Because of
#'   this requirement the use of `tt()` is not allowed in our function. As well
#'   as that, strata by covariate interaction terms are not allowed. Please do
#'   not use "." in the formula at this time either. There is a function
#'   included in the package that will handle these requirements and provide
#'   warnings if your model as constructed is unlikely to work with this
#'   function.
#'
#'
#' @returns A shiny app written into specified directory.
#' @examplesIf interactive()
#'
#'   library(survival)
#'   temp_app_dir <- tempdir()
#'   names(leukemia)[names(leukemia) == "x"] <- "treatment"
#'   model1 <- coxph(Surv(time, status) ~ treatment,
#'   leukemia, x=TRUE, model=TRUE)
#'
#'   shine_coxph("Model 1" = model1, app.dir =
#'   temp_app_dir)
#'
#'   filedir <- list.files(temp_app_dir)[1]
#'
#'   shiny::runApp(paste0(temp_app_dir, "/", filedir))
#'   files <- list.files(temp_app_dir)
#'   file.remove(files)
#'
#'
#' @export
#' @importFrom grDevices hcl.colors
#' @importFrom shiny showTab
shine_coxph=function(..., app.dir = NULL)

{
  ########################
  # determine the class of each input argument
  input.list=list(...)
  n.list=length(input.list)
  list.class=rep("",n.list)
  for (i in 1:n.list) {
    list.class[i]=class(input.list[[i]])
  }

  #######################
  # extract the coxph models from the input
  cox.list=which(list.class=="coxph")
  n.model=length(cox.list)
  model.list=vector("list",n.model)
  for (i in 1:n.model) {
    model.list[[i]]=input.list[[cox.list[i]]]
  }
  names(model.list)=names(input.list)[cox.list]

  #######################
  # get colors list
  # RColorBrewer brewer.pal
  clrs=input.list$clrs
  if (is.null(clrs)) {
    clrs=hcl.colors(n.model, "Dark 2", alpha = 1)
  }

  ##########################
  # get app directory
  #app.dir=input.list$app.dir
  if (is.null(app.dir)) {
    app.dir=getwd()
  }
  date.time=as.character(Sys.time())
  date.time=gsub(":","-",date.time,fixed=TRUE)
  date.time=gsub(" ","-",date.time,fixed=TRUE)
  app.dir=paste0(app.dir,"/",date.time,"/")
  dir.create(app.dir)
  message(paste0("Shiny app will be written in directory ",app.dir,"."))

  #######################
  # convert the coxph models to coxfit objects
  coxfit.list=vector("list",n.model)
  for (i in 1:n.model) {
    coxfit.list[[i]]=prep_coxfit(model.list[[i]])
  }
  names(coxfit.list)=names(model.list)

  #######################
  # Save app data into app directory
  cox.fit.list=coxfit.list
  save(cox.fit.list,clrs,file=paste0(app.dir,"appData.Rdata"))

  #######################
  # get different code sections
  input.data.code=write_coxfit_input_data_code(coxfit.list)
  KM.plot.code=write_KM_plot_code(coxfit.list,clrs)
  table.code=prop_haz_tables(cox.fit.list)

  ########################
  # generate header code
  hdr.code=c("#rm(list=ls())",
             "options(stringsAsFactors=FALSE)",
             paste0("load('",gsub("\\\\", "/", app.dir),"appData.Rdata')"))


  ########################
  # generate ui code
  ui.code=c("ui=navbarPage(",
            "             'Cox Model Survival Predictions',",
            "    tabPanel('Plot',",
            "             sidebarLayout(",
            "             sidebarPanel(",
            paste0("                          ", # SUBODH CHANGED FLOW
                   input.data.code$ui.code,
                   ","),
            "textInput('predProbTimes','Times for predicted probabilities',placeholder='Enter values separated by a comma'),", # SUBODH ADDITION
            "actionButton(inputId = 'go', label = 'Generate Plot'),",
            "actionButton(inputId = 'reset', label = 'Reset'),", # SUBODH ADDITION
            "actionButton(inputId = 'app.exit', label = 'Exit App'),",
            "selectInput('clrs', label = 'Choose Colors', choices = palette.pals(), selected = 'Okabe-Ito'),", # color choice
            "hr(),",
            "numericInput('height', 'Plot Height', value = 100),",
            "numericInput('width', 'Plot Width', value = 100),",
            "downloadButton('downloadPlot', 'Download Plot')",
            "),",
            "mainPanel(",
            "h3('Predicted Survival Curve'),",
            "plotOutput(outputId = 'KM'),",
            "h3('Predicted Probability at Fixed Times'),", # SUBODH ADDITION
            "textOutput(outputId='noPredTimes'),", # SUBODH ADDITION
            "tableOutput(outputId = 'cox.times')))),",# SUBODH ADDITION
            "tabPanel('Summary Tables',",
            "mainPanel(",
            # "h3('Hazard Ratio Summary Table'),", # SUBODH ADDITION
            # "tableOutput(outputId = 'HR'),", # SUBODH ADDITION
            # "h3('Assessing the Proportional Hazards Assumption'),", # SUBODH ADDITION
            # "tableOutput(outputId = 'PHA')", # SUBODH ADDITION
            table.code$ui.code, ### NEW
            ")))") # SUBODH CHANGED THE FLOW

  #########################
  # generate server code
  server.code=c("server=function(input,output)",
                "{",
                "          observeEvent(input$app.exit, {stopApp()}) # Exit when exit button is pressed",
                "          observeEvent(input$go, {",

                input.data.code$server.code,
                KM.plot.code$server.code,
                table.code$server.code, ### NEW
                "predProbTable <- cox_times_table(KM.hat,input$predProbTimes)", # SUBODH ADDITION
                "if (is.null(predProbTable)) output$noPredTimes <- renderText('No input times detected. If you provided times, check that you separated numbers with a single comma and you provided valid numbers.') else output$noPredTimes <- renderText(invisible())", # SUBODH ADDITION
                "output$cox.times=renderTable(predProbTable,rownames=TRUE)", # SUBODH ADDITION
                "output$HR=renderTable(cox.fit.list[[1]]$HR.table,rownames=TRUE)", # SUBODH ADDITION
                "output$PHA=renderTable(cox.fit.list[[1]]$PHA.table$table,rownames=TRUE)", # SUBODH ADDITION
                "colors = palette.colors(length(cox.fit.list), input$clrs)", # colors
                "output$downloadPlot <- downloadHandler(",
                "  filename = function() { paste0('plot.png') },",
                "  content = function(file) {",
                "    png(file, input$width, input$height)",
                "    cox.KM.plots(KM.hat,clrs=colors)",
                "    dev.off()",
                "  })",
                "})",
                "          observeEvent(input$reset, {output$KM <- output$HR <- output$PHA <- output$cox.times <- NULL}) # Reset main area", # SUBODH ADDITION
                "}") # SUBODH CHANGED THE FLOW

  ##########################
  # app code
  app.code=c("cox.app=shinyApp(ui,server)")
             #"runApp(cox.app)")

  ###########################
  # write the app file
  app.code.file=paste0(app.dir,"app.R")
  write(unlist(hdr.code),file=app.code.file,sep="\n")
  write(unlist(ui.code),file=app.code.file,sep="\n",append=TRUE)
  write(unlist(server.code),file=app.code.file,sep="\n",append=TRUE)
  write(unlist(app.code),file=app.code.file,sep="\n",append=TRUE)



  ##########################
  # return result

  res=list(app.dir=app.dir,
           app.data.file=paste0(app.dir,"/appData.Rdata"),
           app.code.file=app.code.file,
           coxfit.list=coxfit.list,
           hdr.code=hdr.code,
           ui.code=ui.code,
           server.code=server.code)

  return(res)

}
