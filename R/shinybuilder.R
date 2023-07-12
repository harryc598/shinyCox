################################
# write the shiny code to generate plots

write_KM_plot_code=function(cox.fit.list,clrs)
{
  n.models=length(cox.fit.list)
  ui.code=server.code=NULL

  ############
  # server code to initialize KM.hat object

  initialize.KM.hat=c("n.models=length(cox.fit.list)",
                      "KM.hat=vector('list',n.models)",
                      "lp=rep(NA,n.models)",
                      "names(KM.hat)=names(cox.fit.list)") # here
  server.code=c(server.code,
                initialize.KM.hat)

  ##########
  # server code to compute KM.hat object

  compute.KM.hat=c("for (i in 1:n.models)",
                   "{",
                   "   km.hat=predict_one_coxfit(cox.fit.list[[i]],new.data)",
                   "   lp[i]=attr(km.hat,'lp')",
                   "   sfit=list(time=km.hat$time,surv=km.hat$surv)",
                   "   class(sfit)='survfit'",
                   "   KM.hat[[i]]=sfit",
                   "}")
  server.code=c(server.code,
                compute.KM.hat)

  #########
  # server and ui code to display KM plots

  display.KM.server=c("output$KM=renderPlot({cox_KM_plots(KM.hat,clrs=colors)})")
  display.KM.ui=c("plotOutput(outputId = 'KM')")

  ui.code=c(ui.code,
            display.KM.ui)
  server.code=c(server.code,
                display.KM.server)

  res=list(ui.code=ui.code,
           server.code=server.code)

  return(res)
}


#############################
#' Generate Cox predicted KM plots
#'
#' @param KM.hat Kaplan-Meier estimate
#' @param clrs color of lines, consider for removal
#' @returns Plot of predicted survival curves
#' @export
#' @importFrom grDevices rainbow
#' @importFrom graphics lines
#' @importFrom graphics legend
cox_KM_plots=function(KM.hat,clrs=NULL)

{
  n.models=length(KM.hat)
  if (is.null(clrs)) {
    clrs=rainbow(n.models)
  }

  if (is.null(names(KM.hat))) {
    names(KM.hat)=paste0("model ",1:n.models)
  }

  max.time=0
  for (i in 1:n.models) {
    max.time=max(max.time,
                 max(KM.hat[[i]]$time,na.rm=TRUE))
  }

  plot(c(0,1.2*max.time),
       c(0,1),xlab="Time",las=1,
       ylab="Prob",type="n")

  for (i in 1:n.models) {
    lines(KM.hat[[i]],col=clrs[i], lwd = 2)
  }

  legend(1.05*max.time,1,
         col=clrs,lwd=1,
         legend=names(KM.hat),
         cex=1)
}

#############################
# Generate Cox predicted times table: SUBODH NEW ADDITION

predSurvTime <- function(kmIn,timeIn) { # expects a data frame with columns of time and surv
  kmIn$surv[max(which(kmIn$time <= timeIn))]
}

#' Create table of Cox predicted probabilities
#'
#' @param KM.hat Survival estimator
#' @param fixTimes vector of times
#' @returns Table of predicted probabilities for each time given
#' @export
cox_times_table=function(KM.hat,fixTimes=NULL)

{
  n.models=length(KM.hat)

  if (is.null(names(KM.hat))) {
    names(KM.hat)=paste0("model ",1:n.models)
  }

  if (is.null(fixTimes) | fixTimes=="") {
    return(NULL)
  } else {
    predTimes <- as.numeric(unlist(strsplit(fixTimes,split=","))) # expects an input character string of numbers each separated by a comma
    if (any(is.na(predTimes))) return(NULL)
  }

  tabOut <- matrix(nrow=n.models,ncol=length(predTimes))
  rownames(tabOut) <- names(KM.hat)
  colnames(tabOut) <- paste("Time:",predTimes)

  for (i in 1:n.models){
    for (j in 1:length(predTimes)){
      tabOut[i,j] <- predSurvTime(KM.hat[[i]],predTimes[j])
    }
  }
  return(tabOut)
}


#################################
# write the shiny code to obtain user inputs

write_coxfit_input_data_code=function(cox.fit.list)

{

  if(is.null(names(cox.fit.list))) {
    names(cox.fit.list)=paste0("model ",1:length(cox.fit.list))
  }

  ###############
  # Get the set of input variables across all models
  vnames=get_vnames_cox_fits(cox.fit.list)

  ############
  # Get the range of numeric predictor variables
  num.x.rng.mtx=get_xrng_cox_fits(cox.fit.list,
                                  vnames)

  ###########
  # Get the levels of categorical predictor variables
  cat.lvls=get_levels_cox_fits(cox.fit.list,
                               vnames)
  # NEW
  #######################################################################
  # Get levels of logic predictors
  logic.lvls <- get_logic_cox_fits(cox.fit.list, vnames)
  ########################################################################
  #############
  # Generate shiny code for each variable

  ui.code=server.code=NULL # initialize shiny code for ui and server

  if(!is.null(cat.lvls))
  {
    for (i in 1:length(cat.lvls))
    {
      cat.pick=ez_pickone(names(cat.lvls)[i],
                          tools::toTitleCase(names(cat.lvls)[i]),
                          cat.lvls[i])
      ui.code=c(ui.code,cat.pick$ui.code)
      server.code=c(server.code,
                    cat.pick$server.code)
    }
  }
  #NEW
  ######################################################################
  if(!is.null(logic.lvls))
  {
    for (i in 1:length(logic.lvls))
    {
      logic.pick=ez_pickone_logic(names(logic.lvls)[i],
                                  tools::toTitleCase(names(logic.lvls)[i]),
                                  logic.lvls[i])
      ui.code=c(ui.code,logic.pick$ui.code)
      server.code=c(server.code,
                    logic.pick$server.code)
    }
  }
  ##########################################################################

  if(!is.null(num.x.rng.mtx))
  {
    for (i in 1:ncol(num.x.rng.mtx))
    {
      x.slider=ez_slider(colnames(num.x.rng.mtx)[i],
                         colnames(num.x.rng.mtx)[i],
                         num.x.rng.mtx[1,i],
                         num.x.rng.mtx[2,i],
                         mean(num.x.rng.mtx[,i])) # Harrison fixed index
      ui.code=c(ui.code,x.slider$ui.code)
      server.code=c(server.code,
                    x.slider$server.code)
    }
  }

  new.data.code=paste0("new.data = cbind.data.frame(",
                       paste0(server.code,collapse=","),")")

  code.res=list(ui.code=ui.code,
                server.code=new.data.code)

  return(code.res)

}

#####################################
# get the set of unique predictor
# variable names
# from a list of cox.fit objects

#' @importFrom stats na.omit
get_vnames_cox_fits=function(cox.fit.list)
{
  n.models=length(cox.fit.list)
  var.name=NULL
  var.type=NULL
  for (i in 1:n.models) # loop over models
  {
    var.name=c(var.name,
               names(cox.fit.list[[i]]$types))
    var.type=c(var.type,
               cox.fit.list[[i]]$types)
  }
  dup.name=duplicated(var.name)
  var.name=var.name[!dup.name]
  var.type=var.type[!dup.name]
  # NEW
  #############################
  var.name=na.omit(var.name)
  var.type=na.omit(var.type)
  ###############################
  res=cbind.data.frame(var.name=var.name,
                       var.type=var.type)

  return(res)
}

###############################
# Get the levels for the categorical variables

get_levels_cox_fits=function(cox.fit.list,vnames)
{
  # how to deal with logical
  cat.vars=which(vnames[,"var.type"]!="numeric" & vnames[,"var.type"]!="logical")
  if (length(cat.vars)==0)
    return(NULL)

  n.vars=length(cat.vars)
  n.models=length(cox.fit.list)
  cat.lvls=vector("list",n.vars)
  cat.names=vnames[cat.vars,"var.name"]
  names(cat.lvls)=cat.names

  for (i in 1:n.models)
  {
    mod.lvls=cox.fit.list[[i]]$xlevels
    mod.vars=names(mod.lvls)
    for (j in mod.vars)
    {
      cat.lvls[[j]]=c(cat.lvls[[j]],
                      mod.lvls[[j]])
    }
  }

  for (j in 1:length(cat.lvls)) {
    cat.lvls[[j]]=unique(cat.lvls[[j]])
  }

  # cat.names=gsub("strata(","",cat.names,fixed=TRUE)
  # cat.names=gsub(")","",cat.names,fixed=TRUE)
  cat.names <- gsub("strata\\((\\w*)\\)", "\\1", cat.names)
  names(cat.lvls)=cat.names

  return(cat.lvls)

}

####################################
# Get the range of the numeric variables
# across a list of cox.fit objects
# Error, rng.mtx[1,x.name]=min(x.rng[1,j],rng.mtx[1,x.name],na.rm=T)
#        rng.mtx[2,x.name]=min(x.rng[2,j],rng.mtx[2,x.name],na.rm=T) fixed
get_xrng_cox_fits=function(cox.fit.list,vnames)

{
  num.vars=which(vnames[,"var.type"]=="numeric")
  if (length(num.vars)==0) {
    return(NULL)
  }

  n.models=length(cox.fit.list)
  rng.mtx=matrix(NA,2,length(num.vars))
  rng.mtx[1,]=NA
  rng.mtx[2,]=NA
  colnames(rng.mtx)=vnames[num.vars,"var.name"]


  for (i in 1:n.models)
  {
    x.rng=cox.fit.list[[i]]$num.x.rng
    for (j in 1:ncol(rng.mtx))
    {
      x.name=colnames(rng.mtx)[j]
      rng.mtx[1,x.name] <- min(x.rng[1,x.name],rng.mtx[1,x.name],na.rm = TRUE)
      rng.mtx[2,x.name] <- max(x.rng[2,x.name],rng.mtx[2,x.name],na.rm = TRUE) #changed from min()
    }
  }
  return(rng.mtx)
}


# add tab with plot hazard ratio table
# NEW function
########################################
get_logic_cox_fits <- function(cox.fit.list, vnames) {
  logic.vars=which(vnames[,"var.type"]=="logical")
  if (length(logic.vars)==0) {
    return(NULL)
  }

  n.vars <- length(logic.vars)

  logic.levels <- vector("list", n.vars)
  logic.names <- vnames[logic.vars, "var.name"]
  names(logic.levels) <- logic.names

  for (i in 1:n.vars) {
    logic.levels[[i]] <- c("TRUE", "FALSE")
  }

  return(logic.levels)
}

# NEW function, proportional hazard for each
##################################################################
prop_haz_tables <- function(cox.fit.list) {
  ui.code=c("tabsetPanel(")
  server.code=c()
  for (i in 1:length(cox.fit.list)) {

    server.code = c(server.code,paste0("output$HR", i, "=renderTable(cox.fit.list[[",i, "]]$HR.table,rownames=TRUE)"),
                    paste0("output$PHA", i, "=renderTable(cox.fit.list[[",i, "]]$PHA.table$table,rownames=TRUE)"),
                    paste0("output$title", i, "=renderText(paste(cox.fit.list[[",i,"]]$nsample, 'subjects,', cox.fit.list[[",i, "]]$nevents, 'events'))"))

    if(i < length(cox.fit.list)) {
      ui.code = c(ui.code,paste0("tabPanel('",names(cox.fit.list)[i],"',"),
                  paste0("h5(textOutput(outputId = 'title", i, "')),"),
                  "'Hazard Ratio Summary Table',",
                  paste0("column(12, align = 'center', tableOutput(outputId = 'HR", i, "')),"),
                  "h4('Assesing the Proportional Hazards Assumption'),",
                  paste0("column(12, align = 'center', tableOutput(outputId = 'PHA", i, "')),"),
                  paste0("textOutput(outputId = 'nevents", i, "')"),
                  "),") } else
                    ui.code = c(ui.code,paste0("tabPanel('",names(cox.fit.list)[i],"',"),
                                paste0("h5(textOutput(outputId = 'title", i, "')),"),
                                "'Hazard Ratio Summary Table',",
                                paste0("column(12, align = 'center', tableOutput(outputId = 'HR", i, "')),"),
                                "h4('Assesing the Proportional Hazards Assumption'),",
                                paste0("column(12, align = 'center', tableOutput(outputId = 'PHA", i, "')),"),
                                paste0("textOutput(outputId = 'nevents", i, "')"),
                                ")")
  }
  ui.code=c(ui.code, ")")
  code.res=list(ui.code=ui.code,
                server.code=server.code)
  return(code.res)
}
