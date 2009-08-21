##
## these are OBSOLETE:  see xyplot (not traceplot) method in coda
##
## new functions for operating on Bayesian stuff
## traceplot <- function (x, smooth=TRUE, ylab="", sys="lattice", ...) {
##     UseMethod("traceplot")
##   }

## traceplot.mcmc.list <- function(x,smooth=TRUE,ylab="", sys="lattice", ...) {
##   require(reshape)
##   require(lattice)
##   nchains = length(x)
##   niter = nrow(x[[1]])
##   x2 = melt(data.frame(do.call(rbind,x),
##     iter=rep(1:niter,nchains),
##     chain=rep(1:nchains,each=niter)),
##     id.var=c("chain","iter"))
##   type <- "l"
##   if (smooth) type <- c(type,"smooth")
##   xyplot(value~iter|variable,groups=chain,type=type,data=x2,
##          scales=list(y=list(relation="free")),
##          ylab=ylab,xlab="iteration")
## }

## traceplot.mcmc <- function(x,smooth=TRUE,ylab="", sys="lattice", ...) {
##   require(reshape)
##   require(lattice)
##   niter = nrow(x)
##   x2 = melt(data.frame(x,
##     iter=1:niter),
##     id.var="iter")
##   type <- "l"
##   if (smooth) type <- c(type,"smooth")
##   xyplot(value~iter|variable,type=type,data=x2,
##          scales=list(y=list(relation="free")),
##          ylab=ylab,xlab="iteration")
## }

## copied from traceplot in coda package, but uses lattice
## traceplot.mcmc <- function (x, data = NULL, outer, aspect = "xy",
##                             default.scales = list(y=list(relation = "free")), 
##                             start = 1, thin = 1,
##                             main = attr(x, "title"), xlab = "",
##                             ylab="",
##                             plot.points = "rug", ...,
##                             subset = coda:::thinned.indices(x, start = start, 
##                               thin = thin))
##   ## would like to supply subset argument that would also
##   ## select specific variables
##   ## consider aspect="fill" as an alternative
## {
##     require(lattice)
##     if (!is.R()) {
##         stop("This function is not yet available in S-PLUS")
##     }
##     if (!missing(outer)) 
##         warning("specification of outer ignored")
##     data <- as.data.frame(x)
##     v <- seq(nrow(x))
##     form <- as.formula(paste(paste(lapply(names(data), as.name),
##         collapse = "+"), "~v"))
##     xyplot(form, data = data[subset, ], outer = TRUE, aspect = aspect, 
##            default.scales = default.scales, main = main, xlab = xlab,
##            ylab=ylab,
##            type="l",
##            plot.points = plot.points, ...)
## }

HPDregionplot <- function(x,vars=1:2,h=c(1,1),n=50,lump=TRUE,prob=0.95,
                          xlab=NULL,ylab=NULL,...) {
  require("MASS") ## for kde2d
  parnames <- if (class(x)=="mcmc.list") colnames(x[[1]]) else colnames(x)
  if (is.character(vars)) {
    vars <- match(vars,parnames)
    if (any(is.na(vars))) stop("some variable names do not match")
  }
  varnames <- parnames[vars]
  mult <- (class(x)=="mcmc.list" && !lump)
  if (mult) stop("multiple chains without lumping not yet implemented")
  if (class(x)=="mcmc.list") {
    if (lump) var1 <- c(sapply(x,function(z)z[,vars[1]]))
    else var1 <- lapply(x,function(z)z[,vars[1]])
  } else var1 <- x[,vars[1]]
  if (class(x)=="mcmc.list") {
    if (lump) var2 <- c(sapply(x,function(z)z[,vars[2]]))
    else var2 <- lapply(x,function(z)z[,vars[2]])
  } else var2 <- x[,vars[2]]
  if (!mult) {
    post1 = kde2d(var1,var2,n=n,h=h)
    ## post0 = post1
  } else {
    post1 = mapply(kde2d,var1,var2,MoreArgs=list(n=n))
  }
  dx = diff(post1$x[1:2])
  dy = diff(post1$y[1:2])
  sz = sort(post1$z)
  c1 = cumsum(sz)*dx*dy
  levels = sapply(prob,function(x) {approx(c1,sz,xout=1-x)$y})
  ## meanvec <- c(mean(var1),mean(var2))
  if (is.null(xlab)) xlab <- varnames[1]
  if (is.null(ylab)) ylab <- varnames[2]
  contour(post1$x,post1$y,post1$z,level=levels,
          xlab=xlab,ylab=ylab,drawlabels=FALSE,...)
  invisible(contourLines(post1$x,post1$y,post1$z,level=levels))
}

## make an mcmc object out of an mcmc.list
lump.mcmc.list <- function(x) {
  x2 <- do.call("rbind",x)
  mcpars <- sapply(x,attr,"mcpar")
  class(x2) <- "mcmc"
  if (var(mcpars[1,])>0 || var(mcpars[3,])>0)
    stop("can't combine chains with unequal start/thin")
  attr(x2,"mcpar") <- c(mcpars[1,1],sum(mcpars[2,]),mcpars[3,1])
  x2
}

perturb.params = function(base,alt,which,
  mult=FALSE,use.base=TRUE) {
  if (!missing(alt)) {
    chlist = mapply(
      function(name,vals) {
        x = lapply(vals,
               function(z2) {
                 if (mult) {
                   base[[name]] = base[[name]]*z2}
                 else base[[name]]=z2
                 base
               })
        if (is.list(x)) x else list(x)
      },
      names(alt),alt)
    chlist = unlist(chlist,recursive=FALSE)
  }
  if (use.base) {
    chlist=append(list(base),chlist)
  }
  chlist
}

