lambertW = function(z,b=0,maxiter=10,eps=.Machine$double.eps,
    min.imag=1e-9) {
  badz = is.na(z)
  z.old = z
  z = z[!badz]
  if (any(round(Re(b)) != b))
    stop("branch number for W must be an integer")
  if (!is.complex(z) && any(z<0)) z=as.complex(z)
  ## series expansion about -1/e
  ##
  ## p = (1 - 2*abs(b)).*sqrt(2*e*z + 2);
  ## w = (11/72)*p;
  ## w = (w - 1/3).*p;
  ## w = (w + 1).*p - 1
  ##
  ## first-order version suffices:
  ##
  w = (1 - 2*abs(b))*sqrt(2*exp(1)*z + 2) - 1
  ## asymptotic expansion at 0 and Inf
  ##
  v = log(z + as.numeric(z==0 & b==0)) + 2*pi*b*1i;
  v = v - log(v + as.numeric(v==0))
  ## choose strategy for initial guess
  ##
  c = abs(z + exp(-1));
  c = (c > 1.45 - 1.1*abs(b));
  c = c | (b*Im(z) > 0) | (!Im(z) & (b == 1))
  w = (1 - c)*w + c*v
  ## Halley iteration
  ##
  for (n in 1:maxiter) {
    p = exp(w)
    t = w*p - z
    f = (w != -1)
    t = f*t/(p*(w + f) - 0.5*(w + 2.0)*t/(w + f))
    w = w - t
    if (any(is.na(t) | is.na(w))) {
      print(t)
      print(w)
    }
    if (abs(Re(t)) < (2.48*eps)*(1.0 + abs(Re(w)))
        && abs(Im(t)) < (2.48*eps)*(1.0 + abs(Im(w))))
      break
  }
  if (n==maxiter) warning(paste("iteration limit (",maxiter,
        ") reached, result of W may be inaccurate",sep=""))
  if (all(Im(w[!is.na(w)])<min.imag)) w = as.numeric(w)
  if (sum(badz)>1) {
    w.new = numeric(length(z.old))
    w.new[!badz] = w
    w.new[badz] = NA
    w = w.new
  }
  return(w)
}

apply2d = function(fun,x,y,...) {
  if (is.character(fun)) fun <- get(fun)
  return(matrix(apply(expand.grid(x,y),1,function(z) { fun(z[1],z[2],...)}),
                nrow=length(x)))
}

## TO DO: log scales
curve3d <- function (expr, from=c(0,0), to=c(1,1), n = c(41,41), add = FALSE, 
                     zlab = NULL, log = NULL, 
                     sys3d = c("persp","wireframe","rgl","contour","image",
                       "none"),
                     ...) 
{
  sys3d <- match.arg(sys3d)
  if (add && !(sys3d %in% c("contour","rgl")))
    stop("can only add contour or rgl to 3D plots")
  n = rep(n,length.out=2)
  sexpr <- substitute(expr)
  if (is.name(sexpr)) {
    fcall <- paste(sexpr, "(x,y)")
    expr <- parse(text = fcall)
    if (is.null(zlab)) 
      zlab <- fcall
  } else {
    if (!(is.call(sexpr) && 
          (match("x", all.vars(sexpr), nomatch = 0) || 
           match("y", all.vars(sexpr), nomatch = 0))))
      stop("'expr' must be a function or an expression containing 'x' and 'y'")
    expr <- sexpr
    if (is.null(zlab)) 
      zlab <- deparse(sexpr)
  } 
  lg <- if (length(log)) 
    log
  if (length(lg) == 0) 
    lg <- ""
  x <- if (lg != "" && "x" %in% strsplit(lg, NULL)[[1]]) {
        if (any(c(from[1], to[1]) <= 0)) 
            stop("'from[1]' and 'to[1]' must be > 0 with log=\"x\"")
        exp(seq(log(from[1]), log(to[1]), length = n))
      }  else seq(from[1], to[1], length = n[1])
    y <- if (lg != "" && "y" %in% strsplit(lg, NULL)[[1]]) {
        if (any(c(from[2], to[2]) <= 0)) 
            stop("'from[2]' and 'to[2]' must be > 0 with log=\"y\"")
        exp(seq(log(from[2]), log(to[2]), length = n))
    } else seq(from[2], to[2], length = n[2])
    tmpfun <- function(x,y) {
      eval(expr, envir = list(x = x, y=y), enclos = parent.frame())
    }
    z <- apply2d(tmpfun,x,y)
    switch(sys3d,
           persp=persp(x,y,z,zlab=zlab,...),
           contour=contour(x,y,z,...),
           image=image(x,y,z,...),
           none=NA,
           wireframe={require("lattice"); 
                      dimnames(z) <- list(x=x,y=y);
                      print(wireframe(z))},
           rgl={require("rgl"); rgl::persp3d(x,y,z,zlab=zlab,add=add,...)})
  invisible(list(x=x,y=y,z=z))
}

get.emdbook.packages <- function() {
   pkglist = c("adapt","bbmle","chron","coda","ellipse","gplots","gtools","gdata",
     "MCMCpack","odesolve","plotrix","R2WinBUGS","reshape","rgl",
     "scatterplot3d")
   inst.pkgs = rownames(installed.packages())
   newpkgs <- pkglist[!pkglist %in% inst.pkgs]
   if (length(newpkgs)>0)
     sapply(pkglist,install.packages)
 }


## expr: expression to evaluate (raw form)
## meanval: values of the mean (possibly named)
## vars: names
## Sigma: var-cov function
deltavar <- function(fun,meanval=NULL,vars,Sigma,verbose=FALSE) {
  expr <- as.expression(substitute(fun))
  if (missing(vars)) {
    if (missing(meanval) || is.null(names(meanval)))
      stop("must specify either variable names or named values for means")
    vars <- names(meanval)
  }
  derivs <- try(lapply(vars,D,expr=expr),silent=TRUE)
  if (inherits(derivs,"try-error")) {
    if (length(grep("is not in the derivatives table",derivs))) {
      ## take numeric derivative
      nderivs <- with(as.list(meanval),
                      numericDeriv(expr[[1]],theta=vars))
      nderivs <- attr(nderivs,"gradient")
    } else {
      stop(paste("Error within derivs:",derivs))
    }
  } else {
    if (verbose) print(derivs)
    nderivs <- sapply(derivs,eval,envir=as.list(meanval))
    if (verbose) cat(nderivs,"\n")
  }
  if (!is.matrix(Sigma) && length(Sigma)>1) Sigma <- diag(Sigma)
  ## if (!is.matrix(Sigma)) sum(Sigma*nderivs^2) else
  if (is.matrix(nderivs)) {
    apply(nderivs,1,function(z) c(z %*% Sigma %*% matrix(z)))
  } else c(nderivs %*% Sigma %*% matrix(nderivs))
  ## really only want diagonal
}

deltamethod <- function(fun,z,var="x",params=NULL,max.order=2) {
  d0 <- as.expression(substitute(fun))
  dvals <- list()
  dvals[[1]] <- d0
  for (i in 2:(max.order+1)) {
    dvals[[i]] <- D(dvals[[i-1]],var)
  }
  mvals <- numeric(max.order-1)
  m <- mean(z)
  for (i in 1:(max.order-1)) {
    mvals[i] <- mean((z-m)^(i+1))
  }
  mvals[1] <- var(z)  ## kluge
  ev1 <- c(as.list(params),list(m))
  ev2 <- c(as.list(params),list(z))
  names(ev1)[length(ev1)] <-   names(ev2)[length(ev2)] <- var
  r0 <- mean(eval(d0,ev2)) ## true value
  evals <- sapply(dvals,eval,ev1)  ## evaluated derivatives
  gvals <- gamma(c(1,3:(max.order+1)))
  deltavals <- cumsum(c(1,mvals)*c(evals[-(2)])/gvals)
  results <- c(r0,deltavals)
  names(results) = c("delta","E(f(x))",paste("delta",2:max.order,sep=""))
  results
}

## zero-inflated negative binomial

dzinbinom = function(x,mu,size,zprob,log=FALSE) {
  logv = log(1-zprob) + dnbinom(x,mu=mu,size=size,log=TRUE)
  logv = ifelse(x==0,log(zprob+exp(logv)),logv)
  if (log) logv else exp(logv)
}

rzinbinom = function(n,mu,size,zprob) {
  ifelse(runif(n)<zprob,0,rnbinom(n,mu=mu,size=size))
}

## add q, p functions for zinbinom?  where else is zinbinom
## implemented?

lseq <- function(from,to,length.out) {
  exp(seq(log(from),log(to),length.out=length.out))
}

## utility function for formatting
scinot <- function(x,format=c("latex","expression"),delim="$",
                   pref="",...) {
  format <- match.arg(format)
  y <- strsplit(as.character(formatC(x,format="e",...)),"e")[[1]]
  y[1] <- gsub("^0+","",y[1])
  y[2] <- ifelse(length(grep("^\\+",y[2]))>0,
                 gsub("^\\+0+","",y[2]),
                 gsub("^-0+","-",y[2]))
  if (format=="latex") { 
    v <- paste(delim,y[1],"\\\\times 10^{",y[2],"}",delim,sep="")
  } else if (format=="expression") {
    if (as.numeric(y[1])==1) {
      v <- substitute(expression(paste(pref,10^b)),list(pref=pref,b=as.numeric(y[2])))
    } else {
      v <- substitute(expression(paste(pref,a %*% 10^b)),
                      list(pref=pref,a=as.numeric(y[1]),b=as.numeric(y[2])))
    }
  }
  v
}
 
## convert R2WinBUGS output to coda/mcmc
as.mcmc.bugs <- function(x) {
  if (!require("coda")) stop("coda is required to use as.mcmc.bugs")
  if (x$n.chains>1) {
    z <- list()
    for (i in 1:x$n.chains) {
      z[[i]] <- mcmc(x$sims.array[,i,],start=1,thin=x$n.thin)
    }
    class(z) <- "mcmc.list"
  } else {
    z <- mcmc(x$sims.matrix,start=1,thin=x$n.thin)
  } 
  return(z)
}

## credible interval for a theoretical distribution
tcredint = function(dist,parlist,ranges,level=0.95,eps=1e-5,verbose=FALSE) {
  qfun = function(x) do.call(paste("q",dist,sep=""),c(list(x),parlist))
  dfun = function(x) do.call(paste("d",dist,sep=""),c(list(x),parlist))
  pfun = function(x) do.call(paste("p",dist,sep=""),c(list(x),parlist))
  if (missing(ranges))  ## set upper/lower limits for search by quantiles
    ranges <- qfun(c(eps,0.5,1-eps))
  lims <- function(pdens) { ## find lower and upper values for which prob dens = target value
    lower <- uniroot(function(p) {dfun(p)-pdens},
                     interval=ranges[1:2])$root
    upper <- uniroot(function(p) {dfun(p)-pdens},
                     interval=ranges[2:3])$root
    c(lower,upper)
  }
  limarea <- function(pdens) { ## find area between target values
     intlim <- lims(pdens)
     d <- diff(pfun(intlim))
     ## cat(pdens,intlim,d,"\n")
     d
   }
  ## these limits must be within limits set above
  v1 <- qfun(c(0.6,0.9999)) # quantiles
  v2 <- dfun(v1)  ## bracketing densities
  u <- uniroot(function(x) {limarea(x)-level},
                 interval=v2)
  intlim <- lims(u$root)
  r = c(intlim,dfun(intlim[1]),limarea(u$root))
  names(r) = c("lower","upper","p","area")
  if (verbose) r else r[1:2]
}

## credible interval for (1D) posterior distribution stored in a numeric
## vector.  assumed unimodal!  pvec (vector of parameter values),
## npost (vector of posterior densities), level, tolerance
ncredint <- function(pvec,npost,level=0.95,tol=0.01,verbose=FALSE) {
  dx = diff(pvec)[1]
  cumdist <- cumsum(npost)*dx
  midpt <- which.min(abs(cumdist-0.5))
  lims <- function(pdens) { ## find lower and upper values for which
    ## prob dens is closest to target value
    lower <- which.min(abs(npost[1:midpt]-pdens))
    upper <- which.min(abs(npost[(midpt+1):length(npost)]-pdens))+midpt
    c(lower,upper)
  }
  limarea <- function(pdens) {
    intlim <- lims(pdens)
    d <- sum(npost[intlim[1]:intlim[2]])*dx
    ##    cat(pdens,intlim,d,"\n")
    d
  }
  ## find credible interval
  v2 <- seq(0,max(npost),by=tol)
  vals <- sapply(v2,limarea)
  w <- which.min(abs(vals-level))
  r = c(pvec[lims(v2[w])],v2[w],limarea(v2[w]))
  names(r) = c("lower","upper","p","area")
  ## credible intervals; posterior density, area
  if (verbose) return(r) else return(r[1:2])
}

## credible interval for a sample (e.g. Markov chain)
##
## OBSOLETE: HPDinterval [coda] is better
##
## qcredint <- function(x,level=0.95,tol=0.01,verbose=FALSE,...) {
##   d = density(x,...)
##   pvec = d$x
##   npost = d$y
##   dx = diff(pvec)[1]
##   cumdist <- cumsum(npost)*dx
##   midpt <- which.min(abs(cumdist-0.5))
##   lims <- function(pdens) { ## find lower and upper values for which
##     ## prob dens is closest to target value
##     lower <- which.min(abs(npost[1:midpt]-pdens))
##     upper <- which.min(abs(npost[(midpt+1):length(npost)]-pdens))+midpt
##     c(lower,upper)
##   }
##   limarea <- function(pdens) {
##     intlim <- lims(pdens)
##     d <- sum(npost[intlim[1]:intlim[2]])*dx
##     ##    cat(pdens,intlim,d,"\n")
##     d
##   }
##   ## find credible interval
##   v2 <- seq(0,max(npost),by=tol)
##   vals <- sapply(v2,limarea)
##   w <- which.min(abs(vals-level))
##   r = c(pvec[lims(v2[w])],v2[w],limarea(v2[w]))
##   names(r) = c("lower","upper","p","area")
##   ## credible intervals; posterior density, area
##   if (verbose) return(r) else return(r[1:2])
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

calcslice <- function(fit1,fit2,fn=fit1@minuslogl,range=c(-0.1,1.1),
                  n=400) {
  slicep = seq(range[1],range[2],length=n)
  slicepars = t(sapply(slicep,function(x) (1-x)*coef(fit1)+x*coef(fit2)))
  v = apply(slicepars,1,function(x) do.call("fn",as.list(x)))
  list(x=slicep,y=v)
}

rchibarsq <- function(n,df=1,mix=0.5) {
  ifelse(runif(n)>mix,
         rchisq(n,df),
         if (df==1) 0 else rchisq(n,df-1))
}

dchibarsq <- function(x,df=1,mix=0.5,log=FALSE) {
  df <- rep(df,length.out=length(x))
  mix <- rep(mix,length.out=length(x))
  c1 <- ifelse(df==1,0,dchisq(x,df-1))
  c2 <- dchisq(x,df)
  r <- mix*c1+(1-mix)*c2
  zeros <- (x==0 & df==1)
  if (any(zeros)) {
    r[zeros] <- Inf
  }
  if (log) log(r) else r
}

pchibarsq <- function(p,df=1,mix=0.5,lower.tail=TRUE,log.p=FALSE) {
  df <- rep(df,length.out=length(p))
  mix <- rep(mix,length.out=length(p))
  c1 <- ifelse(df==1,if (lower.tail) 1 else 0,
               pchisq(p,df-1,lower.tail=lower.tail))
  c2 <- pchisq(p,df,lower.tail=lower.tail)
  r <- mix*c1+(1-mix)*c2
  if (log.p) log(r) else r
}

qchibarsq <- function(q,df=1,mix=0.5) {
  n <- max(length(q),length(df),length(mix))
  df <- rep(df,length.out=n)
  mix <- rep(mix,length.out=n)
  q <- rep(q,length.out=n)
  tmpf2 <- function(q,df,mix) {
    if (df>1) {
      tmpf <- function(x) {
        pchibarsq(x,df,mix)-q
      }
      uniroot(tmpf,lower=qchisq(q,df-1),upper=qchisq(q,df))$root
    } else {
      newq <- (q-mix)/(1-mix)
      ifelse(newq<0,0,qchisq(newq,df=1))
    }
  }
  mapply(tmpf2,q,df,mix)
  ##  if (any(df>1)) stop("df>1 not implemented yet")
  ## ?? need uniroot() solution to find quantiles of
  ##   mixtures?
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

dmvnorm <- function (x, mu, Sigma, log = FALSE, tol = 1e-06) {
    if (is.vector(x)) 
        x = t(as.matrix(x))
    n = length(mu)
    if (is.vector(mu)) {
        p <- length(mu)
        if (is.matrix(x)) {
            mu <- matrix(rep(mu, nrow(x)), ncol = p, byrow = TRUE)
        }
    }
    else {
        p <- ncol(mu)
    }
    if (!all(dim(Sigma) == c(p, p)) || nrow(x) != nrow(mu)) 
        stop("incompatible arguments")
    eS <- eigen(Sigma, sym = TRUE, EISPACK = TRUE)
    ev <- eS$values
    if (!all(ev >= -tol * abs(ev[1]))) 
        stop("Sigma is not positive definite")
    z = t(x - mu)
    logdetS = try(determinant(Sigma, logarithm = TRUE)$modulus)
    attributes(logdetS) <- NULL
    iS = try(solve(Sigma))
    if (class(iS) == "try-error" || class(logdetS) == "try-error") {
        warning("difficulty inverting/taking determinant of Var-Cov matrix")
        return(NA)
    }
    ssq = diag(t(z) %*% iS %*% z)
    loglik = -(n * (log(2*pi)) +  logdetS + ssq)/2
    if (log) loglik else exp(loglik)
  }


dbetabinom <- function(x,prob,size,theta,shape1,shape2,log=FALSE) {
  if (missing(prob) && !missing(shape1) && !missing(shape2)) {
    prob = shape1/(shape1+shape2)
    theta = shape1+shape2
  }
  v <- lchoose(size,x)-lbeta(theta*(1-prob),theta*prob)+lbeta(size-x+theta*(1-prob),x+theta*prob)
  if (log) v else exp(v)
}

rbetabinom <- function(n,prob,size,theta,shape1,shape2) {
  if (!missing(prob) && !missing(size) && missing(shape1) && missing(shape2)) {
    shape1 <- theta*prob
    shape2 <- theta*(1-prob)
  }
  rbinom(n,size=size,prob=rbeta(n,shape1,shape2))
}

## could implement distribution function as:
## D(x) = 1 - (n B(b+n-x-1,a+x+1) Gamma(n) 3_F_2(1,a+x+1,-n+x+1;
##   x+2, -b-n+x+2;1))/(B(a,b) B(n-x,x+2) Gamma(n+2))
## (from Mathworld) -- would need to get 3F2 somewhere?
## not in GSL!
## http://tolstoy.newcastle.edu.au/R/e2/help/07/02/10988.html
##  suggests Davies package could be adapted
##  or ?? SuppDists ??

gridsearch2d <- function(fun,v1min,v2min,v1max,v2max,n1=20,n2=20,
                         logz=FALSE,sys3d=c("both","contour","image"),...) {
  sys3d=match.arg(sys3d)
  redraw = function(v1vec,v2vec,m) {
    if (sys3d=="contour") {
      contour(v1vec,v2vec,m,xlab="",ylab="",...)
    } else image(v1vec,v2vec,m,xlab="",ylab="",...)
    if (sys3d=="both") contour(v1vec,v2vec,m,add=TRUE,...)
  }
  recalc = function(v1vec,v2vec,logz) {
    m = apply2d(fun,v1vec,v2vec)
    mindm = pmax(diff(sort(m)),1e-10)
    if (logz) m = log10(m-min(m)+mindm)
    m
  }
  v1vec = seq(v1min,v1max,length=n1)
  v2vec = seq(v2min,v2max,length=n2)
  m = recalc(v1vec,v2vec,logz)
  stop1 = FALSE
  first = TRUE
  while (!stop1) {
    redraw(v1vec,v2vec,m)
    if (!first) {
      resp = readline("Continue (y|n)? ")
      stop1 = (toupper(resp)=="N")
    }
    first = FALSE
    if (!stop1) {
      cat("click on box corners\n")
      z = lapply(locator(2),sort)
      rect(z$x[1],z$y[1],z$x[2],z$y[2])
      resp = readline("OK (y|n)? ")
      if (toupper(resp)=="Y") {
        v1min=z$x[1]
        v1max=z$x[2]
        v2min=z$y[1]
        v2max=z$y[2]
        v1vec = seq(v1min,v1max,length=n1)
        v2vec = seq(v2min,v2max,length=n2)
        m = recalc(v1vec,v2vec,logz)
      }
    }
  }
  resp = readline("get point (y|n)? ")
  if (toupper(resp)=="Y") {
    cat("click on point\n")
    z=locator(1)
    return(z)
  } else invisible(NULL)
}

metropSB <- function(fn,start,deltap=NULL,
                      scale=1,rptfreq=-1,
                      acceptscale=1.01,rejectscale=0.99,
                      nmax=10000,retvals=FALSE,retfreq=100,
                      verbose=FALSE,...) {
  ## initialization
  ndim <- length(start)  # number of parameters/dimensions
  p <- start             # current parameter vector = starting vector
  minp <- start          # parameter vector corresponding to minimum value so far
  val <- fn(p,...)     # current function value
  minval <- val          # minimum value so far
  if (retvals) {
    info <- matrix(nrow=round(nmax/retfreq),ncol=3*ndim+3)
    dimnames(info) <- list(NULL,c(paste("p",1:ndim,sep=""),
                                  paste("minp",1:ndim,sep=""),
                                  paste("deltap",1:ndim,sep=""),
                                        "val","minval","accept"))
    # save info on parameter values, function value, whether accepted or not, current best parameters,
    # current best value, current jump size
  }
  it <- 1                # iteration count
  if (is.null(deltap)) deltap <- p*0.05   # default starting jump size(s)
  while (it<=nmax) {
    oldp <- p                             # save current parameters and function value
    oldval <- val
    p <- p + runif(ndim,-1,1)*deltap      # perturb current values
    val <- fn(p,...)                    # new function value
    dval <- val-oldval                    # change
    saveinfo <- (it %% retfreq == 0)
    savecount <- it %/% retfreq
    if (saveinfo) {
      info[savecount,-ncol(info)] <- c(p,minp,deltap,val,minval)
    }
    if (verbose)
      if (it %% 100 == 0)
         cat(it,"\n")
    if (dval<0 || (exp(-dval*scale)>runif(1))) {
      ## accept new value
      if (saveinfo) info[savecount,ncol(info)] <- 1
      if (val<minval) {  # if better than best so far, record new minimum
	minval <- val
	minp <- p
      }
      deltap <- deltap*acceptscale   # increase jump size
    } else {
      ## replace old values
      if (saveinfo) info[savecount,ncol(info)] <- 0
      val <- oldval
      p <- oldp
      deltap <- deltap*rejectscale   # decrease jump size
    }
    it <- it+1
    if ((rptfreq>0) && (it %% rptfreq)==1)
      cat("it=",it," value=",val," min. val.=",minval,"\n")
  }
  if (retvals)
    return(list(minimum=minval,estimate=minp,funcalls=it,retvals=info))
  else
    return(list(minimum=minval,estimate=minp,funcalls=it))
}

trcoef = function(x,inverse=FALSE) {
  ## n.b. logit should come before log
  vals=list(pattern=c("logit","log","sqrt"),
    fun=c(plogis,exp,function(x) {x^2}),
    ifun=c(qlogis,log,sqrt))
  if (is.null(attr(x,"transf"))) attr(x,"transf") = character(length(x))
  for (p in 1:length(vals[[1]])) {
    pat = vals$pattern[[p]]
    fun = vals$fun[[p]]
    ifun = vals$ifun[[p]]
    if (!inverse) {
      w = grep(paste("^",pat,sep=""),names(x))
      if (length(w)>0) {
        attr(x,"transf")[w]=pat
        names(x)[w] = gsub(pat,"",names(x)[w])      
        x[w] = sapply(x[w],fun)
      }
    } else {
      w = which(attr(x,"transf")==pat)
      if (length(w)>0) {
        attr(x,"transf")[w] = ""
        names(x)[w] = paste(pat,names(x)[w],sep="")
        x[w] = sapply(x[w],ifun)
      }
    }
  }
  if (all(attr(x,"transf")=="")) attr(x,"transf") = NULL
  x
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
