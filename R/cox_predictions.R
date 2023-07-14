######################################################################
# Create a simplified representation of a coxph model fit
# with minimum information necessary to compute model predictions

#' Create simplified `coxph()` object
#'
#' Simplifies `coxph()` output and checks that predictions
#' match those of the full model
#' @param coxph.result Result returned by `coxph()`
#' @param tol numerical tolerance for prediction differences
#' @returns Simplified `coxph()` result, in list form, include what's in it
#' @export
prep_coxfit=function(coxph.result,  # result of coxph
                     tol=1e-7)      # numerical tolerance of prediction differences

{
  check_coxph(coxph.result)               # check that coxph.result is adequate for our purposes
  cox.fit=simplify_coxph(coxph.result)    # simplify coxph.result to cox.fit
  check_coxfit(cox.fit,coxph.result,tol)  # check that predictions of cox.fit match those of coxph.result
  return(cox.fit)                         # return cox.fit
}


###############################################################################################
# Simplifies a coxph.result object for visualization in the app
# Strips the predictor data and model data to best protect patient privacy
# Keeps the model coefficients and the assignment of input data to model data

#' @importFrom stats predict
#' @importFrom survival survfit
#' @importFrom stats coef
#' @importFrom stats confint
#' @importFrom survival cox.zph
simplify_coxph=function(coxph.result)

{
  ########
  # extract survival and linear predictor data from coxph.result
  sdat=coxph.result$y
  lp=stats::predict(coxph.result,type="lp")
  lp.surv=cbind.data.frame(sdat=sdat,lp=lp)

  ########
  # compute the baseline Cox survival function in a data.frame format
  bl.cox=survival::survfit(coxph.result)
  bl.surv=cbind.data.frame(time=bl.cox$time,  # baseline survival function corresponding to baseline hazard function
                           surv=bl.cox$surv)

  ########
  # add strata column if necessary
  if (!is.null(bl.cox$strata))
  {
    strt=rep(names(bl.cox$strata),bl.cox$strata)
    if(any(grepl("=", strt, fixed = TRUE))) {
      eq.pos=regexpr("=",strt,fixed=TRUE)
      strt.vname=substring(strt[1],1,eq.pos[1]-1)
    } else {
      where.strata <- which(grepl("strata(", names(coxph.result$model), fixed = TRUE))
      strt.vname <- names(coxph.result$model)[where.strata]
      strt.vname <- gsub("strata(", "", strt.vname, fixed = TRUE)
      strt.vname <- gsub(")", "", strt.vname, fixed=TRUE)
    }
    bl.surv=cbind.data.frame(strata=strt,bl.surv)
    names(bl.surv)=c(strt.vname,"time","surv")
  }

  ########
  # extract tables of estimates and evaluation of proportional hazards assumption
  cox.smry=stats::coef(summary(coxph.result))
  cox.CIs=stats::confint(coxph.result)
  cox.pha=survival::cox.zph(coxph.result)
  HR.tbl=cbind.data.frame("Hazard Ratio"=cox.smry[,"exp(coef)"],
               "Lower Bound"=exp(cox.CIs[,1]),
               "Upper Bound"=exp(cox.CIs[,2]),
               "p value"=cox.smry[,"Pr(>|z|)"])
  pha.tbl=cox.pha
  nevents=coxph.result$nevent
  nsample=coxph.result$n


  ########
  # extract minimal sufficient information to compute survival estimates

  cox.terms=attr(coxph.result$terms,"term.labels")
  cox.types=attr(coxph.result$terms,"dataClasses")
  cox.types=cox.types[cox.terms]
  #cox.coefs=coef(coxph.result)

  # EXPERIMENTING
  ##############################################################
  coefs <- stats::coef(coxph.result)
  if(any(is.na(coefs))) {
    stop("One or more of your coefficients is NA")
  }
  pattern <- "\\("
  badnames <- which(grepl(pattern, names(coefs)))
  if(any(grepl("\\)\\w+$", names(coefs)))) {
    new_coef_names <- gsub("\\)", "\\)\\`", names(coefs)[badnames])
    new_coef_names <- paste0("`", new_coef_names)
    attr(coefs, "names")[badnames] = new_coef_names
  } else {
    newcoefnames <- paste0("`", names(coefs)[badnames], "`")
    attr(coefs, "names")[badnames] = newcoefnames
  }
  cox.coefs=coefs
  ##############################################################

  #cox.assign=coxph.result$assign
  cox.means=coxph.result$means
  xlevels=coxph.result$xlevels
  strata=coxph.result$strata
  form=coxph.result$formula
  x.rng=apply(coxph.result$x,2,range)
  rownames(x.rng)=c("minimum","maximum")

  res=list(bl.surv=bl.surv,    # baseline survival function estimate
           lp.surv=lp.surv,    # subject level linear predictor and survival data
           #assign=cox.assign,  # assignment of input data values to regression model matrix values
           types=cox.types,    # vector of data types (character, numeric, strata, etc)
           coefs=cox.coefs,    # regression coefficient estimates
           means=cox.means,    # means of the regression model matrix columns
           num.x.rng=x.rng,    # range of the regression model matrix columns
           xlevels=xlevels,    # levels of categorical predictors
           strata=strata,      # strata
           form=form,          # model formula
           HR.table=HR.tbl,    # hazard ratio estimates
           PHA.table=pha.tbl,  # tests of proportional hazards assumption
           nevents=nevents,
           nsample=nsample)    # number of events

  return(res)                  # returns list object with all the above information
}



#####################################################
# Computes predicted survival outcomes for one patient based on a Cox model fit
#' Compute predicted survival outcome
#'
#' Computes predicted survival outcome for one patient based
#' on a Cox model fit.
#' @param coxfit coxfit object result of prep.coxfit
#' @param newdata vector of new data
#' @returns data.frame of predicted survival outcomes
#' @export
predict_one_coxfit=function(coxfit,          # result of prep.coxfit
                            newdata)         # new input data vector
{

  ok=check_coxfit_newdata(coxfit,newdata)               # check input newdata
  x=compute_coxfit_xvector(coxfit,newdata)              # convert the input data into an x vector
  coef.names=names(coxfit$coefs)                        # extract the names of the coeffients
  if(length(coef.names) == 1){
    x.beta=sum(coxfit$coefs*(x-coxfit$means))
  } else
    x.beta=sum(coxfit$coefs*(x[coef.names]-coxfit$means)) # see help(predict.coxph) for notes on centering
  res=coxfit$bl.surv                                    # initialize result as baseline survival function
  res[,"surv"]=coxfit$bl.surv$surv^exp(x.beta)          # formula for predicted survival of one new patient

  ###########
  # if stratified model, then limit to the stratum of newdata
  if (!is.null(coxfit$strata))
  {
    #strt.var=coxfit$strata[1]
    strt.col.name <- names(coxfit$bl.surv[,!names(coxfit$bl.surv) %in% c("time", "surv"), drop = FALSE])
    strt.var <- newdata[, strt.col.name]
    # Add note
    if(grepl("=", strt.var, fixed = TRUE)) {
    eq.pos=regexpr("=",strt.var,fixed=TRUE)
    strt.var=substring(strt.var,1,eq.pos-1)
    strt.mtch <- (newdata[, strt.var] == res[, strt.var])
    } else if(is.character(strt.var)) {
    strt.column <- which(grepl(strt.var, newdata, fixed = TRUE))
    strt.mtch=(newdata[,strt.column]==res[,1])
    } else {
      strt.mtch <- (newdata[, strt.col.name] == res[, 1])
    }
    if (!any(strt.mtch)) {
      stop(paste0("Unable to match strata variable ",
                  strt.var," in newdata."))
    }

    res=res[strt.mtch,]
  }

  res=as.data.frame(res)
  attr(res,"lp")=x.beta
  attr(res,"newdata")=newdata
  return(res)
}


#############################################
# compute the x vector for a newdata observation for a coxfit object
#' @importFrom stats as.formula
#' @importFrom stats model.matrix
compute_coxfit_xvector=function(coxfit,
                                newdata)

{
  # create right-hand side of formula
  #form=deparse(coxfit$form)
  # This works with interactions x*y already
  # form <- Reduce(paste, deparse(coxfit$form))
  # tilde.pos=regexpr("~",form,fixed=T)
  # form=substring(form,tilde.pos)
  #
  # # use model matrix to create x vector
  # form=as.formula(form)
  #EXPERIMENTING
  ####################################################################
  testform <- Reduce(paste, deparse(coxfit$form))
  splitform <- strsplit(testform, "~")[[1]][2]
  fixed.functions <- gsub("([a-zA-Z0-9\\.]+\\(.*?\\))", "\\`\\1\\`", splitform)
  fixedstrata <- gsub("\\`(strata\\(\\w*\\))\\`", "\\1", fixed.functions)
  fixedstrata <- gsub("\\`(\\w+\\:strata\\(.*\\))\\`", "\\1", fixedstrata)
  readyform <- paste0("~", fixedstrata)
  formasform <- stats::as.formula(readyform, env = parent.frame())
  mtx <- stats::model.matrix(formasform, data=newdata, xlev = coxfit$xlevels)
  ######################################################

  #mtx=model.matrix(form,data=newdata,xlev=coxfit$xlevels)
  if(length(mtx[1,]) == 1) {
    x=mtx[1,names(coxfit$coefs), drop = FALSE]
  } else x=mtx[1,names(coxfit$coefs)]
  return(x)
}


##############################################################
# Check that predictions from simplified coxph model object
# match those of the input coxph model for each subject

check_coxfit=function(cox.fit,coxph.result,tol=1e-7)

{
  message("Double checking calculations from minimal coxph fit information to those from full coxph information.")
  cox.input=coxph_input_data(coxph.result)
  n=nrow(cox.input)
  for (i in 1:n)
  {
    if(length(cox.input[i, -1]) == 1) {
      new.data=as.data.frame(cox.input[i,-1, drop = FALSE])          # get data for subject i
    } else
      new.data=as.data.frame(cox.input[i,-1])
    # --------------------------------------------------------------------
    # This section is made to handle certain variable transformations
    # such as log(x). Due to model.matrix's behavior, the unaltered
    # 'x' from the original data frame needs to be added to our newdata

    vars <- all.vars(coxph.result$formula)
    if(!exists(coxph.result$call$data)) {
      stop("Dataset(s) must be in environment while the function is running.
           Afterwords the original data are no longer needed.")
    }
    original.data <- eval(coxph.result$call$data)
    original.data <- as.data.frame(original.data)
    notvars <- setdiff(names(original.data), vars)
    original.data <- original.data[,!names(original.data) %in% notvars]
    stringform <- Reduce(paste, deparse(coxph.result$formula))
    output <- strsplit(stringform, " ~ ")[[1]][1]        # remove response from formula
    isSurv <- grepl("Surv\\(", output)                   # Checks form of response
    if(isSurv) {
      lhsadj1 <- strsplit(output, "Surv\\(")[[1]][2]       # strips Surv( from response
      lhsadj2 <- strsplit(lhsadj1, "\\,")[[1]]             # Splits string by comma separation
      lhsadj3 <- gsub("\\(", "", lhsadj2)                  # Removes parenthesis
      lhsadj4 <- gsub("\\)", "", lhsadj3)                  # Removes paranthesis
      lhsadj5 <- strsplit(lhsadj4, "\\$")                  # If data$variable, removes data$
      lhsadj6 <- c()
      for (j in 1:length(lhsadj5)) {                       # makes vector of names in response
        if(is.na(lhsadj5[[j]][2])) {
          lhsadj6 <- c(lhsadj6, lhsadj5[[j]][1])
        } else {
          lhsadj6 <- c(lhsadj6, lhsadj5[[j]][2])
        }
        if(grepl("\\s", lhsadj6[j])) {
          lhsadj6[j] <- gsub("\\s*(\\w)?\\s*", "\\1", lhsadj6[j])
        }
      }
      original.data <- original.data[, !names(original.data) %in% lhsadj6, drop = FALSE]          # removes response values if in original dataset
    } else
    {original.data <- original.data[, !names(original.data) %in% output, drop = FALSE]}
    unused <- setdiff(names(original.data), names(new.data))            # Looks for any variables in original data not in newdata
    if (inherits(coxph.result$na.action, "omit")) {
    original.data <- as.data.frame(na.omit(original.data))
    }
    new.data <- cbind.data.frame(new.data, original.data[i, unused, drop = FALSE])
    ########################################################################
    cox.pred1=predict_one_coxfit(cox.fit,new.data)                 # prediction with new object
    cox.pred2=one_survfit(coxph.result,newdata=new.data)           # prediction by survival package
    pred.diff=abs(cox.pred1[,"surv"]-cox.pred2$surv)               # absolute value of difference
    ok=(max(pred.diff)<tol)                                        # check if within numerical error
    if (!ok) {                                                       # stop if error detected
      stop("Error encountered in double-checking survival prediction calculations from reduced coxph result object.")
    }
  }
  return(invisible())
}

###################################
# Compute predicted survival outcomes using survival package

#' @importFrom survival survfit
one_survfit=function(coxph.result,newdata)

{

  sf.res=survival::survfit(coxph.result,newdata=newdata)
  res=cbind.data.frame(time=sf.res$time,
                       surv=sf.res$surv)
  # NOT NEEDED AFTER SURVIVAL 3.3
  # if (!is.null(coxph.result$strata))
  # {
  #   strt=rep(names(sf.res$strata),
  #            sf.res$strata)
  #   res=cbind.data.frame(strata=strt,res)
  #   vname=res$strata[1]
  #   eq.pos=regexpr("=",vname,fixed=T)
  #   vname=substring(vname,1,eq.pos-1)
  #   mtch=which(res$strata==newdata[,vname])
  #   res=res[mtch,]
  # }
  return(res)
}



#######################################
# Recreate input data from coxph.result object

coxph_input_data=function(coxph.result)

{
  res=coxph.result$model
  if (!is.null(coxph.result$strata))
  {
    orig.vname=coxph.result$strata[1]
    eq.pos=regexpr("=",orig.vname,fixed=TRUE)
    orig.vname=substring(orig.vname,1,eq.pos-1)

    orig.value=substring(coxph.result$strata,eq.pos+1)

    new.vname=paste0("strata(",orig.vname,")")
    #res[,new.vname]=orig.value
    new.clm=grep(new.vname,colnames(res),fixed=T)
    colnames(res)[new.clm]=orig.vname
  }
  res=as.data.frame(res)

  return(res)
}


###############################################
# Check that a coxph result object is suitable for generating a shiny app

check_coxph=function(coxph.result)

{
  #########
  # stop for invalid input
  if (!any(class(coxph.result)=="coxph")) {
    stop("The argument coxph.object must be the result object of the coxph function.")
  }

  cox.names=names(coxph.result)
  if (!all(is.element(c("x","model"),cox.names))) {
    stop("To perform this calculation, coxph.result must have components 'model' and 'x'.  Rerun coxph with model=TRUE and x=TRUE.")
  }


  #########
  # get the terms of the Cox model and their data types
  cox.terms=attr(coxph.result$terms,
                 "dataClasses")       # gets all the terms of the model
  cox.terms=cox.terms[-1]             # drop the survival response variable to get predictor terms

  ########
  # Stop if there are no predictor terms
  if (length(cox.terms)==0)         # no predictor terms
  {
    stop("The coxph.result object has no predictor terms.")
  }

  return(invisible())

}


############################################
# Check input newdata for predict.one.coxfit

check_coxfit_newdata=function(coxfit,   # result of prep.coxfit
                              newdata)  # data.frame with data for one patient

{

  ##############
  # Check that input data has only one row
  if (is.data.frame(newdata))
  {
    if (nrow(newdata)>1) {
      stop("Input data MUST have only ONE row.")
    }
    newdata=as.vector(newdata[1,])
  }

  ##############
  # check that newdata contains essential names
  new.names=names(newdata)
  cox.names=names(coxfit$assign)
  if (length(setdiff(cox.names,new.names))>0) {
    stop("newdata missing some essential elements for prediction.")
  }

  return(invisible())
}

