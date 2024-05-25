#' Obtains information for standard errors of predictions
#'
#' @description Computes necessary information to calculate standard errors and
#'   confidence intervals in shiny app. This is adapted from parts of
#'   'survfit.coxph'. This function is meant to be used in conjuction with
#'   'predict_se'.
#'
#' @param model a 'coxph' object
#' @param ctype whether the cumulative hazard computation should have a
#'   correction for ties, 1=no, 2=yes.
#' @param individual deprecated argument, replaced by 'id'
#' @param id optional variable name of subject identifiers. Not supported in app
#' @param se.fit a logical value indicating whether standard errors should be
#'   computed. Default is TRUE for standard models, FALSE for multi-state (code
#'   not yet present for that case.)
#' @param stype computation of the survival curve, 1=direct, 2= exponenial of
#'   the cumulative hazard. Default is 2.
#' @returns A list of information needed for computing predicted standard
#'   errors.
#' @examplesIf interactive()
#' library(survival)
#'
#' colondeaths <- colon[colon$etype == 2, ]
#' split_colon <- split(colondeaths, colondeaths$rx)
#' colon_arm1 <- split_colon$Obs
#'
#' colon1ph <- coxph(Surv(time, status) ~ factor(extent) + nodes + strata(surg)
#'                    + factor(differ),
#'                    colon_arm1,
#'                    x = TRUE, model = TRUE)
#' surv_pred_info(colon1ph)
#'
#' @import stats
#' @import survival
#' @export
surv_pred_info = function(model, ctype, individual = FALSE, id, se.fit = TRUE, stype = 2) {
  object <- model

  Terms  <- terms(object)
  robust <- !is.null(object$naive.var)

  if (!is.null(attr(object$terms, "specials")$tt)) {
    stop("The survfit function can not process coxph models with a tt term")
  }


  if (missing(ctype)) {
    # Use the appropriate one from the model
    temp1 <- match(object$method, c("exact", "breslow", "efron"))
    ctype <- c(1,1,2)[temp1]
  } else if (!(ctype %in% 1:2)) {stop ("ctype must be 1 or 2")}
  if (!(stype %in% 1:2)) {stop("stype must be 1 or 2")}

  tfac <- attr(Terms, 'factors')
  temp <- attr(Terms, 'specials')$strata
  has.strata <- !is.null(temp)
  if (has.strata) {
    stangle = untangle.specials(Terms, "strata")  #used multiple times, later
    # Toss out strata terms in tfac before doing the test 1 line below, as
    #  strata end up in the model with age:strat(grp) terms or *strata() terms
    #  (There might be more than one strata term)
    for (i in temp) tfac <- tfac[,tfac[i,] ==0]  # toss out strata terms
  } else stangle = NULL
  if (any(tfac >1)) {
    stop("not able to create a curve for models that contain an interaction without the lower order effect")
  }

  Terms <- object$terms
  n <- object$n[1]
  if (!has.strata) {strata <- NULL} else {strata <- object$strata}

  if (!missing(individual)) {warning("the `id' option supersedes `individual'")}
  # This part will need addressing, for now missid will just be set to true
  missid <- missing(id) # I need this later, and setting id below makes
  # "missing(id)" always false
  #missid = TRUE

  if (!missid) {individual <- TRUE
  } else if (missid && individual) {id <- rep(0L,n)
  }else {id <- NULL}
  #id = NULL
  # Maybe this is needed? Don't know yet
  # if (individual & missing(newdata)) {
  #   stop("the id option only makes sense with new data")
  # }

  if (has.strata) {
    temp <- attr(Terms, "specials")$strata
    factors <- attr(Terms, "factors")[temp,]
    strata.interaction <- any(t(factors)*attr(Terms, "order") >1)
  }

  # We don't allow multistate at this point, so maybe remove this
  coxms <- inherits(object, "coxphms")
  if (coxms || is.null(object$y) || is.null(object[['x']]) ||
      !is.null(object$call$weights) || !is.null(object$call$id) ||
      (has.strata && is.null(object$strata)) ||
      !is.null(attr(object$terms, 'offset'))) {

    mf <- stats::model.frame(object)
  } else mf <- NULL

  position <- NULL
  Y <- object[['y']]
  if (is.null(mf)) {
    weights <- object$weights  # let offsets/weights be NULL until needed
    offset <- NULL
    offset.mean <- 0
    X <- object[['x']]
  } else {
    weights <- model.weights(mf)
    offset <- model.offset(mf)
    if (is.null(offset)) offset.mean <- 0
    else {
      if (is.null(weights)) offset.mean <- mean(offset)
      else offset.mean <- sum(offset * (weights/sum(weights)))
    }
    X <- model.matrix.coxph(object, data=mf)
    if (is.null(Y) || coxms) {
      Y <- model.response(mf)
      if (is.null(object$timefix) || object$timefix) Y <- aeqSurv(Y)
    }
    oldid <- model.extract(mf, "id")
    if (length(oldid) && ncol(Y)==3) position <- survflag(Y, oldid)
    else position <- NULL
    if (!coxms && (nrow(Y) != object$n[1]))
      stop("Failed to reconstruct the original data set")
    if (has.strata) {
      if (length(strata)==0) {
        if (length(stangle$vars) ==1) strata <- mf[[stangle$vars]]
        else strata <- strata(mf[, stangle$vars], shortlabel=TRUE)
      }
    }

  }

  varmat <- object$var
  beta <- ifelse(is.na(object$coefficients), 0, object$coefficients)
  xcenter <- sum(object$means * beta) + offset.mean
  if (!is.null(object$frail)) {
    keep <- !grepl("frailty(", dimnames(X)[[2]], fixed=TRUE)
    X <- X[,keep, drop=F]
  }

  if (is.null(offset)) {risk <- c(exp(X%*% beta - xcenter))
  }else   {risk <- c(exp(X%*% beta + offset - xcenter))}

  wt <- weights
  if (missing(strata) || length(strata)==0) {strata <- rep(0L, nrow(Y))}

  if (is.factor(strata)) {ustrata <- levels(strata)
  } else {ustrata <- sort(unique(strata))}
  nstrata <- length(ustrata)
  survlist <- vector('list', nstrata)
  names(survlist) <- ustrata
  survtype <- if (stype==1) {1
  }else {ctype+1}
  vartype <- survtype
  if (is.null(wt)) wt <- rep(1.0, nrow(Y))
  if (is.null(strata)) strata <- rep(1L, nrow(Y))
  for (i in 1:nstrata) {
    indx <- which(strata== ustrata[i])
    survlist[[i]] <- agsurv(Y[indx,,drop=F], X[indx,,drop=F],
                                       wt[indx], risk[indx],
                                       survtype, vartype)
  }

  stuff_for_later <- list(survlist = survlist, Terms = Terms, has.strata = has.strata,
                         stangle = stangle, xlevels = object$xlevels, means = object$means,
                         beta = beta, xcenter = xcenter, se.fit = se.fit, varmat = varmat,
                         survtype = survtype)

  return(stuff_for_later)
}

#' Creates predicted survival and standard errors for confidence intervals
#'
#' @description Adapted from parts of 'survfit.coxph', computes predictions for
#'   standard errors based on 'surv_pred_info' output and 'newdata' from the
#'   shiny app.
#'
#' @param listsurv Output from 'surv_pred_info' function
#' @param coxfit coxfit object created for predictions. Used to find strata
#' @param newdata Data used to make predicted standard errors
#'
#' @returns a list of number of subjects for each curve, times at which the
#'   curve has a step, number at risk for each time, number of events at each
#'   time, number censored at each time (no event but exit risk set), estimated
#'   survival, cumulative hazard at each transition, and standard error of the
#'   cumulative hazard.
#'
#' @examplesIf interactive()
#'   library(survival)
#'   library(shinyCox)
#'   colondeaths <- colon[colon$etype == 2, ]
#'   split_colon <- split(colondeaths, colondeaths$rx)
#'
#'   colon_arm1 <- split_colon$Obs
#'   colon1ph <- coxph(Surv(time, status) ~
#'   factor(extent) + nodes + strata(surg) + factor(differ), colon_arm1, x =
#'   TRUE, model = TRUE)
#'
#'   new.data = cbind.data.frame(`factor(extent)` = 3, `surg` =
#'   "surg=0",`factor(differ)` = 2,`nodes` = 5)
#'
#'
#'   coxfit = prep_coxfit(colon1ph)
#'   coxlist = surv_pred_info(colon1ph)
#'
#'   predict_se(coxlist, coxfit, new.data)
#'
#' @import survival
#' @import stats
#' @export
predict_se = function(listsurv, coxfit, newdata) {
  object <- list()
  survlist <- listsurv[[1]]
  Terms <- listsurv[[2]]
  has.strata <- listsurv[[3]]
  stangle <- listsurv[[4]]
  object$xlevels <- listsurv[[5]]
  object$means <- listsurv[[6]]
  beta <- listsurv[[7]]
  xcenter <- listsurv[[8]]
  se.fit <- listsurv[[9]]
  varmat <- listsurv[[10]]
  #strata <- listsurv[[11]]
  survtype <- listsurv[[11]]
  found.strata <- FALSE
  ustrata <- NULL

  strt.col.name <- names(coxfit$bl.surv[, !names(coxfit$bl.surv) %in% c("time", "surv"), drop = FALSE])

  strt.var <- newdata[, strt.col.name]
  #strata.new = paste0(strt.col.name, '=', strt.var)

  survlist <- list(survlist[[strt.var]])

  step_1 <- gsub("factor", "", names(newdata))
  step_2 <- gsub("\\(", "", step_1)
  step_3 <- gsub("\\)", "", step_2)
  names(newdata) <- step_3
  strt.val <- strsplit(newdata[[strt.col.name]], "=")[[1]][2]
  newdata[[strt.col.name]] <- as.numeric(strt.val)

  Call <- match.call()
  Call[[1]] <- as.name("survfit")
  individual = FALSE
  if (!individual)  {
    Terms2 <- delete.response(Terms)
    y2 <- NULL  # a dummy to carry along, for the call to coxsurv.fit
  }

  if (is.vector(newdata, "numeric")) {
    if (individual) stop("newdata must be a data frame")
    if (is.null(names(newdata))) {
      stop("Newdata argument must be a data frame")
    }
    newdata <- data.frame(as.list(newdata), stringsAsFactors=FALSE)
  }  else if (is.list(newdata)) newdata <- as.data.frame(newdata)
  if (has.strata) {
    found.strata <- TRUE
    tempenv <- new.env(, parent=emptyenv())
    assign("strata", function(..., na.group, shortlabel, sep)
      list(...), envir=tempenv)
    assign("list", list, envir=tempenv)
    for (svar in stangle$vars) {
      temp <- try(eval(parse(text=svar), newdata, tempenv),
                  silent=TRUE)
      if (!is.list(temp) ||
          any(unlist(lapply(temp, class))== "function"))
        found.strata <- FALSE
    }

    if (!found.strata) {
      ss <- untangle.specials(Terms2, "strata")
      Terms2 <- Terms2[-ss$terms]
    }
  }

  tcall <- Call[c(1, match(c('id', "na.action"),
                           names(Call), nomatch=0))]
  tcall$data <- newdata
  tcall$formula <- Terms2
  tcall$xlev <- object$xlevels[match(attr(Terms2,'term.labels'),
                                     names(object$xlevels), nomatch=0)]
  tcall[[1L]] <- quote(stats::model.frame)
  mf2 <- eval(tcall)

  if (has.strata && found.strata) { #pull them off
    temp <- untangle.specials(Terms2, 'strata')
    # strata2 <- strata(mf2[temp$vars], shortlabel=TRUE)
    # strata2 <- factor(strata2, levels=levels(strata))
    # if (any(is.na(strata2)))
    #   stop("New data set has strata levels not found in the original")
    # An expression like age:strata(sex) will have temp$vars= "strata(sex)"
    #  and temp$terms = integer(0).  This does not work as a subscript
    if (length(temp$terms) >0) Terms2 <- Terms2[-temp$terms]
  } else strata2 <- factor(rep(0, nrow(mf2)))

  offset2 <- model.offset(mf2)
  if (length(offset2)==0 ) offset2 <- 0
  # a model with only an offset, but newdata containing a value for it
  if (length(object$means)==0) {x2 <- 0
  } else x2 <- model.matrix(Terms2, mf2)[,-1, drop=FALSE]

  if (length(object$means)) {
    risk2 <- exp(c(x2 %*% beta) + offset2 - xcenter)
  } else risk2 <- exp(offset2 - xcenter)

  expand <- function(fit, x2, varmat, se.fit) {
    if (survtype==1) {
      surv <- cumprod(fit$surv)
    } else surv <- exp(-fit$cumhaz)

    if (is.matrix(x2) && nrow(x2) >1) {  #more than 1 row in newdata
      fit$surv <- outer(surv, risk2, '^')
      dimnames(fit$surv) <- list(NULL, row.names(x2))
      if (se.fit) {
        varh <- matrix(0., nrow=length(fit$varhaz), ncol=nrow(x2))
        for (i in 1:nrow(x2)) {
          dt <- outer(fit$cumhaz, x2[i,], '*') - fit$xbar
          varh[,i] <- (cumsum(fit$varhaz) + rowSums((dt %*% varmat)* dt))*
            risk2[i]^2
        }
        fit$std.err <- sqrt(varh)
      }
      fit$cumhaz <- outer(fit$cumhaz, risk2, '*')
    }
    else {
      fit$surv <- surv^risk2
      if (se.fit) {
        dt <-  outer(fit$cumhaz, c(x2)) - fit$xbar
        varh <- (cumsum(fit$varhaz) + rowSums((dt %*% varmat)* dt)) *
          risk2^2
        fit$std.err <- sqrt(varh)
      }
      fit$cumhaz <- fit$cumhaz * risk2
    }
    fit
  }
  result <- lapply(survlist, expand, x2, varmat, se.fit)
  unlist = TRUE
  if (unlist) {
    if (length(result)==1) { # the no strata case
      if (se.fit) {
        result = result[[1]][c("n", "time", "n.risk", "n.event", "n.censor",
                               "surv", "cumhaz", "std.err")]
        return(result)
      } else {result[[1]][c("n", "time", "n.risk", "n.event", "n.censor",
                            "surv", "cumhaz")]}
    }
    else {
      temp <-list(n   =    unlist(lapply(result, function(x) x$n),
                                  use.names=FALSE),
                  time=    unlist(lapply(result, function(x) x$time),
                                  use.names=FALSE),
                  n.risk=  unlist(lapply(result, function(x) x$n.risk),
                                  use.names=FALSE),
                  n.event= unlist(lapply(result, function(x) x$n.event),
                                  use.names=FALSE),
                  n.censor=unlist(lapply(result, function(x) x$n.censor),
                                  use.names=FALSE),
                  strata = sapply(result, function(x) length(x$time)))
      names(temp$strata) <- names(result)


        temp$surv <- t(matrix(unlist(lapply(result,
                                            function(x) t(x$surv)), use.names=FALSE),
                              nrow= nrow(x2)))
        dimnames(temp$surv) <- list(NULL, row.names(x2))
        temp$cumhaz <- t(matrix(unlist(lapply(result,
                                              function(x) t(x$cumhaz)), use.names=FALSE),
                                nrow= nrow(x2)))
        if (se.fit)
          temp$std.err <- t(matrix(unlist(lapply(result,
                                                 function(x) t(x$std.err)), use.names=FALSE),
                                   nrow= nrow(x2)))

      return(temp)
    }
  }
  else {
    names(result) <- ustrata
    return(result)
  }
}

#' Get confidence intervals for predicted survival curves
#'
#' Creates confidence levels for plotting predicted survival curves.
#'
#' @param p Vector of survival probabilities
#' @param se Vector of standard errors
#' @param conf.type Type of confidence interval, includes 'plain', 'log',
#' 'log-log', 'logit', and 'arcsin'.
#' @param conf.int The level for two-sided confidence interval on the predicted
#' survival curve, default is 0.95.
#' @param ulimit Should upper bound be limited to 1, default is TRUE
#' @returns list of length two, containing the lower and upper confidence levels
#'
#' @examplesIf interactive()
#' library(survival)
#' library(shinyCox)
#' colondeaths <- colon[colon$etype == 2, ]
#' split_colon <- split(colondeaths, colondeaths$rx)
#'
#' colon_arm1 <- split_colon$Obs
#' colon1ph <- coxph(Surv(time, status) ~ factor(extent) + nodes + strata(surg)
#'                   + factor(differ),
#'                   colon_arm1,
#'                   x = TRUE, model = TRUE)
#'
#' new.data = cbind.data.frame(`factor(extent)` = 3,
#'                          `surg` = "surg=0",`factor(differ)` = 2,`nodes` = 5)
#'
#'
#' coxfit = prep_coxfit(colon1ph)
#' coxlist = surv_pred_info(colon1ph)
#'
#' for_ci = predict_se(coxlist, coxfit, new.data)
#'
#' get_confint(for_ci$surv, for_ci$std.err, conf.int = 0.95,
#'             conf.type = "log-log")
#'
#' @export
get_confint <- function(p, se, conf.type, conf.int, ulimit=TRUE) {
  zval <- qnorm(1- (1-conf.int)/2, 0,1)
  scale <- 1.0

  if (conf.type=='plain') {
    se2 <- se* p * zval  # matches equation 4.3.1 in Klein & Moeschberger
    if (ulimit) list(lower= pmax(p -se2*scale, 0), upper = pmin(p + se2, 1))
    else  list(lower= pmax(p -se2*scale, 0), upper = p + se2)
  }
  else if (conf.type=='log') {
    #avoid some "log(0)" messages
    xx <- ifelse(p==0, NA, p)
    se2 <- zval* se
    temp1 <- exp(log(xx) - se2*scale)
    temp2 <- exp(log(xx) + se2)
    if (ulimit) list(lower= temp1, upper= pmin(temp2, 1))
    else  list(lower= temp1, upper= temp2)
  }
  else if (conf.type=='log-log') {
    xx <- ifelse(p==0 | p==1, NA, p)
    se2 <- zval * se/log(xx)
    temp1 <- exp(-exp(log(-log(xx)) - se2*scale))
    temp2 <- exp(-exp(log(-log(xx)) + se2))
    list(lower = temp1 , upper = temp2)
  }
  else if (conf.type=='logit') {
    xx <- ifelse(p==0, NA, p)  # avoid log(0) messages
    se2 <- zval * se *(1 + xx/(1-xx))

    temp1 <- 1- 1/(1+exp(log(p/(1-p)) - se2*scale))
    temp2 <- 1- 1/(1+exp(log(p/(1-p)) + se2))
    list(lower = temp1, upper=temp2)
  }
  else if (conf.type=="arcsin") {
    xx <- ifelse(p==0, NA, p)
    se2 <- .5 *zval*se * sqrt(xx/(1-xx))
    list(lower= (sin(pmax(0, asin(sqrt(xx)) - se2*scale)))^2,
         upper= (sin(pmin(pi/2, asin(sqrt(xx)) + se2)))^2)
  }
  else stop("invalid conf.int type")
}

agsurv <- utils::getFromNamespace("agsurv", "survival")

#' For  use in internal surv_pred_info function, borrowed from survival package
#'
#' @param y Surv object
#' @param id vector of ids coinciding with y.
survflag <- function(y, id) {
  if (!inherits(y, "Surv")) stop("y must be a Surv object")
  if (nrow(y) != length(id)) stop("length mismatch")
  if (ncol(y) != 3) stop("y needs to be of (tstart, tstop) form")

  n <- nrow(y)
  indx <- order(id, y[,2])  # sort the data by time within id
  y2 <- y[indx,]
  id2 <- id[indx]

  newid <- (id2[-n] != id2[-1])
  gap <-  (y2[-n,2] < y2[-1,1])

  flag <- 1L*c(TRUE, newid | gap) + 2L*c(newid | gap, TRUE)
  flag[indx] <- flag   # return it to data order
  flag
}



#' Model.matrix method for coxph models
#'
#' Reconstruct the model matrix for a cox model
#' @noRd
model.matrix.coxph <- function(object, data=NULL,
                               contrast.arg=object$contrasts, ...) {
  #
  # If the object has an "x" component, return it, unless a new
  #   data set is given
  if (is.null(data) && !is.null(object[['x']]))
    return(object[['x']]) #don't match "xlevels"

  Terms <- delete.response(object$terms)
  if (is.null(data)) mf <- stats::model.frame(object)
  else {
    if (is.null(attr(data, "terms")))
      mf <- stats::model.frame(Terms, data, xlev=object$xlevels)
    else mf <- data  #assume "data" is already a model frame
  }

  cluster <- attr(Terms, "specials")$cluster
  if (length(cluster)) {
    temp <- untangle.specials(Terms, "cluster")
    dropterms <- temp$terms
  }
  else dropterms <- NULL

  strats <- attr(Terms, "specials")$strata
  hasinteractions <- FALSE
  if (length(strats)) {
    stemp <- untangle.specials(Terms, 'strata', 1)
    if (length(stemp$vars)==1) strata.keep <- mf[[stemp$vars]]
    else strata.keep <- strata(mf[,stemp$vars], shortlabel=TRUE)
    istrat <- as.integer(strata.keep)

    for (i in stemp$vars) {  #multiple strata terms are allowed
      # The factors attr has one row for each variable in the frame, one
      #   col for each term in the model.  Pick rows for each strata
      #   var, and find if it participates in any interactions.
      if (any(attr(Terms, 'order')[attr(Terms, "factors")[i,] >0] >1))
        hasinteractions <- TRUE
    }
    if (!hasinteractions) dropterms <- c(dropterms, stemp$terms)
  } else istrat <- NULL


  if (length(dropterms)) {
    Terms2 <- Terms[-dropterms]
    X <- model.matrix(Terms2, mf, constrasts.arg=contrast.arg)
    # we want to number the terms wrt the original model matrix
    temp <- attr(X, "assign")
    shift <- sort(dropterms)
    for (i in seq(along.with=shift))
      temp <- temp + 1*(shift[i] <= temp)
    attr(X, "assign") <- temp
  }
  else X <- model.matrix(Terms, mf, contrasts.arg=contrast.arg)

  # drop the intercept after the fact, and also drop strata if necessary
  Xatt <- attributes(X)
  if (hasinteractions) adrop <- c(0, untangle.specials(Terms, "strata")$terms)
  else adrop <- 0
  xdrop <- Xatt$assign %in% adrop  #columns to drop (always the intercept)
  X <- X[, !xdrop, drop=FALSE]
  attr(X, "assign") <- Xatt$assign[!xdrop]
  attr(X, "contrasts") <- Xatt$contrasts
  X
}
