# Generate ui and server code to pick one choice from a dropdown
#' Generate code for `selectInput()`
#'
#' Creates a `selectInput()` function in shiny ui and pairs it with a server
#' input.
#' @param vname variable name, or inputId
#' @param label visual label of select
#' @param choices vector of choices to be made available
#' @returns Vector containing ui and server code
#' @noRd
ez_pickone=function(vname,    # Variable name
                    label,    # Variable label in the GUI
                    choices)  # possible values of the variable


{
  ui.code=paste0("selectInput(inputId = '",vname,
                 "', label = '",label,"'",
                 ", choices = ",
                 paste0(choices,
                        collapse=", "),
                 ")")

  server.code=paste0("`",vname,"` = input$`",vname,"`")
  return(list(ui.code=ui.code,
              server.code=server.code))
}


#----------------------------------------------
# Generate ui and server code for a slider bar (numeric data input)
#' Generate code for `sliderInput()`
#'
#' Creates shiny ui code for `sliderInput()` and creates input for vname in
#' server
#' @param vname Name of variable
#' @param label Label to display in shiny ui
#' @param min minimum value of slider
#' @param max maximum value of slider
#' @param value default slider value
#' @returns Vector containing shiny ui and server code for `sliderInput()`
#' @noRd
ez_slider=function(vname,  # variable name
                   label,  # variable label in the GUI
                   min,    # minimum possible value of the variable
                   max,    # maximum possible value of the variable
                   value)  # default value of the variable

{
  ui.code=paste0("sliderInput(inputId = '",vname,"'",
                 ", label = '",label,"'",
                 ", min = ",min,
                 ", max = ",max,
                 ", value = ",value,")")
  server.code=paste0("`",vname,"` = input$`",vname,"`")
  return(list(ui.code=ui.code,
              server.code=server.code))
}


#' Generate code for `selectInput()`
#'
#' Creates a `selectInput()` function in shiny ui and
#' pairs it with a server input. Only for logical variables
#' @param vname variable name, or inputId
#' @param label visual label of select
#' @param choices vector of choices, must be TRUE and FALSE
#' @returns vector containing ui and server code for `selectInput()`
#' @noRd
ez_pickone_logic=function(vname,    # Variable name
                          label,    # Variable label in the GUI
                          choices)  # possible values of the variable
{
  ui.code=paste0("radioButtons(inputId = '",vname,
                 "', label = '",label,"'",
                 ", choices = ",
                 paste0(choices,
                        collapse=", "),
                 ")")
  server.code=paste0("`",vname,"` = as.logical(input$`",vname, "`)")
  return(list(ui.code=ui.code,
              server.code=server.code))
}





