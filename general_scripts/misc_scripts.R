#' Extension of abcrf:::abcrf.formula
#' allows to manually set importance measure as add. option 'impinp'
#' new default: "impurity_corrected" (Sandri & Zuccolotto (2008)) 

abcrf_imp <- function (formula, data, group = list(), lda = TRUE, ntree = 500, 
                       sampsize = min(1e+05, nrow(data)), paral = FALSE,
                       impinp = "impurity_corrected",
                       ncores = if (paral) max(detectCores() - 1, 1) else 1, ...) {
  if (!inherits(formula, "formula")) 
    stop("abcrf.formula is only for formula objects")
  if (!inherits(data, "data.frame")) 
    stop("data needs to be a data.frame object")
  if (is.na(ncores)) {
    warning("Unable to automatically detect the number of CPU cores, \n1 CPU core will be used or please specify ncores.")
    ncores <- 1
  }
  if ((!is.logical(paral)) || (length(paral) != 1L)) 
    stop("paral should be TRUE or FALSE")
  if (!is.list(group)) 
    stop("group needs to be a list")
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf))
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  if (!is.factor(model.response(mf))) 
    stop("response should be a factor containing the model indexes")
  if (nrow(data) == 0L || is.null(nrow(data))) 
    stop("no simulation in the reference table (response, sumstat)")
  if ((!is.logical(lda)) && (length(lda) != 1L)) 
    stop("lda should be TRUE or FALSE")
  if (length(group) != 0) {
    ngroup <- length(group)
    varn <- formula[[2]]
    data[[as.character(varn)]] <- as.vector(data[[as.character(varn)]])
    allmod <- unique(data[[as.character(varn)]])
    for (k in 1:ngroup) for (l in 1:length(group[[k]])) data[[as.character(varn)]][which(data[[as.character(varn)]] == 
                                                                                           group[[k]][l])] <- paste("g", k, sep = "")
    if (!setequal(allmod, unlist(group))) {
      diffe <- setdiff(allmod, unlist(group))
      for (l in 1:length(diffe)) data <- data[-which(data[[as.character(varn)]] == 
                                                       diffe[l]), ]
    }
    data[[as.character(varn)]] <- as.factor(data[[as.character(varn)]])
  }
  if (lda) {
    model.lda <- lda(formula, data)
    data <- cbind(data, as.matrix(predict(model.lda, data)$x))
  }
  else {
    model.lda <- NULL
  }
  m <- names(match.call(expand.dots = TRUE))
  if ((!"sampsize" %in% m) && (nrow(data) <= 15)) 
    sampsize <- as.integer(sampsize/10)
  if (sampsize > nrow(data)) 
    stop("sampsize too large")
  model.rf <- ranger(formula, data, num.trees = ntree, sample.fraction = sampsize/nrow(data), 
                     num.threads = ncores, keep.inbag = TRUE, importance = impinp, 
                     ...)
  class.error = vector()
  for (i in 1:nrow(model.rf$confusion.matrix)) {
    rowSum <- sum(model.rf$confusion.matrix[i, ])
    accurate <- diag(model.rf$confusion.matrix)[i]
    error <- rowSum - accurate
    class.error[i] <- error/rowSum
  }
  model.rf$confusion.matrix <- cbind(model.rf$confusion.matrix, 
                                     class.error)
  colnames(model.rf$confusion.matrix) <- c(paste(model.rf$forest$levels), 
                                           "class.error")
  model.rf$model.rf
  cl <- match.call()
  cl[[1]] <- as.name("abcrf")
  x <- list(call = cl, lda = lda, formula = formula, group = group, 
            model.rf = model.rf, model.lda = model.lda, prior.err = model.rf$prediction.error)
  class(x) <- "abcrf"
  x
}

environment(abcrf_imp) <- asNamespace('abcrf')

#' Function to perturbate theta
#' Input: A scaled mutation rate theta. Range1, s.t. the perturbed theta comes from [theta/range1,theta*range1],
#' steps: How many log-equidistant steps between the bounds
#' Output: A perturbed theta as the binomial draw (fair coin) from the grid of potential values (between which the steps are, both borders are 
#' possible).

log_smear <- function(theta1,range1=5,steps1=10){logth <- log(theta1)
range1 <- log(range1)
pos1 <- rbinom(1,steps1,0.5)
out1 <- logth - range1 + 2*pos1*range1/steps1
return(exp(out1))}
