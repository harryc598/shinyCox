#' Gets the bits for part2
#' @import stats
#' @import survival
#' @noRd
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
                         strata = strata, survtype = survtype)

  return(stuff_for_later) # need surv as well, can get from earlier stuff.
}

#' Creates predicted survival and standard errors for confidence intervals
#'
#' @param listsurv Output from 'surv_pred_info' function
#' @param coxfit coxfit object created for predictions. Used to find strata
#' @param newdata Data used to make predicted standard errors
#' @import survival
#' @import stats
part2 = function(listsurv, coxfit, newdata) {
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
  strata <- listsurv[[11]]
  survtype <- listsurv[[12]]
  found.strata <- FALSE

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
        return(list(result,
                    #strata2,
                    mf2, has.strata, found.strata, se.fit))
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

     # if ((missing(id2) || is.null(id2)) && nrow(x2)>1) {
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
      # }
      # else {
      #   temp$surv <- unlist(lapply(result, function(x) x$surv),
      #                       use.names=FALSE)
      #   temp$cumhaz <- unlist(lapply(result, function(x) x$cumhaz),
      #                         use.names=FALSE)
      #   if (se.fit)
      #     temp$std.err <- unlist(lapply(result,
      #                                   function(x) x$std.err), use.names=FALSE)
      # }
      return(list(temp,
                  #strata2,
                  mf2, has.strata, found.strata, se.fit))
    }
  }
  else {
    names(result) <- ustrata
    return(list(result,
                #strata2,
                mf2, has.strata, found.strata, se.fit))
  }
}


#' Creates confidence levels for plotting predicted survival curves.
#' @references survival authort
get_confint <- utils::getFromNamespace("survfit_confint", "survival")


agsurv <- utils::getFromNamespace("agsurv", "survival")
