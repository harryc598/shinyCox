########## SUBODH CHANGED ez.pickone

#####################################################################
# Generate ui and server code for to pick one choice from a dropdown
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
                 ")") # SUBODH CHANGED FLOW AND REMOVED THE BELOW:

  # ", choices = c(",  # this adds an extra c()
  # paste0(choices,
  # paste0("'",choices,"'"), # this causes quotes around everything
  #   collapse=", "),"))") # SUBODH CHANGED FLOW
  server.code=paste0("`",vname,"` = input$`",vname,"`")
  return(list(ui.code=ui.code,
              server.code=server.code))
}


#####################################################################
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


###################################################################################
# Generate ui and server code for an action button;
# embeds server.lines into the action of the button
# NOT USED IN CODE
ez_button=function(id,              # identifier of the button
                   label,           # label to show on the button
                   server.lines="") # lines to be executed in server when button is pushed

{
  ui.code=paste0("actionButton(inputId = '",id,
                 "', label = '",label,"')")
  server.code=c(paste0('observeEvent(input$',id,', {'),
                server.lines,'})')
  server.code=gsub('{"','{',server.code,fixed=T)
  res=list(ui.code=ui.code,
           server.code=server.code)
  return(res)
}


###################################
# Generate an exit app button
# NOT USED IN CODE
ez_exit=function(id="app.exit",
                 label="Exit App")
{
  exit.code="stopApp()"
  res=ez_button(id,label,exit.code)
  res$server.code=paste0(res$server.code,collapse="")
  res$server.code=paste0(res$server.code,"# Exit when exit button is pressed")
  return(res)
}

####################################
# Create a sidebar panel layout
# NOT USED IN CODE
ez_sidebarLayout=function(ui.sidebar,  # ui code for the sidebar
                          ui.main)     # ui code for the main panel
{
  # Add commas after all but the last line of ui.sidebar
  n.sidebar=length(ui.sidebar)
  ui.sidebar[-n.sidebar]=paste0(ui.sidebar[-n.sidebar],",")

  # Add some indention space to enhance readability
  ui.sidebar=paste0("    ",ui.sidebar)

  # similar operations for main panel ui code
  n.main=length(ui.main)
  ui.main[-n.main]=paste0(ui.main[-n.main],",")
  ui.main=paste0("    ",ui.main)

  # write the code before and after
  res=c("sidebarLayout(",
        "  sidebarPanel(",ui.sidebar,"    ), ",
        "  mainPanel(",ui.main,"    ))")
  return(res)
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





