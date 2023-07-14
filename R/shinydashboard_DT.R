



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
                      paste0("h5(textOutput(outputId = 'title", i, "')),"),
                      "h4('Hazard Ratio Summary Table'),",
                      paste0("DT::dataTableOutput(outputId = 'HR", i, "'),"),
                      "h4('Assesing the Proportional Hazards Assumption'),",
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
                      paste0("h5(textOutput(outputId = 'title", i, "')),"),
                      "h4('Hazard Ratio Summary Table'),",
                      paste0("DT::dataTableOutput(outputId = 'HR", i, "'),"),
                      "h4('Assesing the Proportional Hazards Assumption'),",
                      paste0("DT::dataTableOutput(outputId = 'PHA", i, "')")
                   )

    }
    servercode <- c(servercode,
                    paste0("output$HR", i, "=DT::renderDataTable(DT::datatable(cox.fit.list[[", i, "]]$HR.table,"),
                           "                                     options = list(",
                           "                                     dom = 't'",
                           "                                     )) |>",
                           "                 DT::formatRound(columns = c('Hazard Ratio', 'Lower Bound', 'Upper Bound', 'p value'), digits = 4) |>",
                           "                 DT::formatStyle('p value',",
                           "                                 target = 'cell',",
                           "                                 fontweight = DT::styleInterval(0.05, c('bold', 'normal'))))",
                    paste0("output$PHA", i, "=DT::renderDataTable(DT::datatable(cox.fit.list[[", i, "]]$PHA.table$table,"),
                           "                                      options = list(",
                           "                                      dom = 't')))",
                           paste0("output$title", i, "=renderText(paste(cox.fit.list[[",i,"]]$nsample, 'subjects,', cox.fit.list[[",i, "]]$nevents, 'events'))")
                    )

  }
  code.res <- list(uitop = uicodetop,
                   uibottom = uicodebottom,
                   server.code = servercode)
}


