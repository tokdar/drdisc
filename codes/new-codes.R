library(BayesLogit)
library(splines)
library(parallel)
library(dplyr)
library(MASS)

tcol <- Vectorize(function(col, alpha = 1) {x <- col2rgb(col)/255; return(rgb(x[1],x[2],x[3],alpha = alpha))}, "col")
extreme <- function(x) return(c(which.min(x), which.max(x)))
extract <- function(lo, varname) return(lo[[varname]])
threshold <- function(x, cutoff) return(x * (abs(x) > cutoff))
half.kern.erode <- function(x, bw=0.16) return(-(x < 0)*exp(-0.5*x^2/bw^2))
half.kern.bunch <- function(x, bw=0.16) return((x >= 0)*exp(-0.5*x^2/bw^2))
half.kern <- ifelse(bunch, half.kern.bunch, half.kern.erode)

pospart <- function(x) return(x*(x>0))
logsum <- function(lx) return(max(lx) + log(sum(exp(lx-max(lx)))))

rtt <- function(n, center=0,scale=1,df=1,lo=-Inf,hi=Inf){
  phi <- pt((hi-center)/scale,df)
  plo <- pt((lo-center)/scale,df)
  u <- runif(n)
  return(center+scale*qt((phi-plo)*u + plo,df))
}
dtt <- Vectorize(function(x, center=0, scale=1,df=1,lo=-Inf,hi=Inf,log=FALSE){
  lden <- -Inf  
  if(lo < x & x < hi){
    phi <- pt((hi-center)/scale,df)
    plo <- pt((lo-center)/scale,df)
    lden <- dt((x-center)/scale,df,log=TRUE) - log(scale) - log(phi-plo)
  }  
  if(!log) lden <- exp(lden)
  return(lden)
}, "x")

###### basis functions ######

type <- "kern"
lpoly <- list(
  #  Vectorize(function(x) return(1/sqrt(2)), "x"), ## k = 0
  function(x) return(sqrt(3/2)*x), ## k = 1
  function(x) return(sqrt(5/8)*(3*x^2 - 1)), ## k = 2
  function(x) return(sqrt(7/8)*(5*x^3 - 3*x)) ## k = 3
)
ypoly <- function(y) return(sapply(lpoly, function(f) f(y)))

kbw <- 0.5
knots <- seq(-1,1,kbw)
bsd.sq <- (0.67*kbw)^2
gausskern <- lapply(knots, function(mu) return(function(x) return(exp(-0.5*(x-mu)^2/bsd.sq))))
ykern <- function(y) return(sapply(gausskern, function(f) f(y)))

nbw <- 0.25
knots <- seq(-1, 1, nbw)
ns.basis <- ns(knots, knots=knots[-extreme(knots)])
yns <- function(y) return(predict(ns.basis, y))

fn.types <- list(poly=ypoly, kern=ykern, ns=yns)
yFn <- fn.types[[type]]

order <- length(as.numeric(yFn(0)))
wFn <- function(y, beta, alpha=0, bw=0.16) return(c(matrix(yFn(y), nrow=length(y)) %*% beta) + alpha*half.kern(y,bw))
Phi <- function(x) plogis(x)
fFn <- function(y, shapes=c(1,1), ...) return(0.5*dbeta((y+1)/2,shapes[1],shapes[2])*Phi(wFn(y, ...)))

y.grid <- seq(-1,1,.01)

get.poly.mat <- yFn

###### get.f and simulate.fx ######

get.f <- function(b, a=0, sh=c(1,1), jbw=0.16) {
  f.grid <- fFn(y.grid, beta=b, alpha=a, shapes=sh, bw=jbw)
  normalization.const <- (integrate(fFn, -1, 0, beta=b, alpha=a, shapes=sh, bw=jbw)$value
                          + integrate(fFn, 0, 1, beta=b, alpha=a, shapes=sh, bw=jbw)$value)
  return(f.grid/normalization.const)
}

simulate.fx <- function(n, b, x, a, shapes=c(1,1), jbw=0.16, positive=FALSE){
  x <- matrix(x, nrow=n)
  p <- ncol(x)
  b <- matrix(b, ncol=p)
  y <- rep(NA, n)
  xa <- c(x %*% a)
  if(positive) xa <- pmax(0, xa)
  
  ix.remain <- 1:n
  n.remain <- n
  while(n.remain > 0){
    z <- 2*rbeta(n.remain, shapes[1],shapes[2])-1
    pm <- matrix(yFn(z), nrow=n.remain)
    hk <- half.kern(z,jbw)
    x.remain <- x[ix.remain,,drop=FALSE]
    
    xa.remain <- xa[ix.remain]
    w <- rowSums(x.remain * (pm %*% b)) + xa.remain*hk
    u <- (runif(n.remain) < Phi(w))
    y[ix.remain[u]] <- z[u]
    ix.remain <- ix.remain[!u]
    n.remain <- length(ix.remain)
  }
  return(y)
}

### get.f0 ###

p0 <- c(.5, .4, .1)
shape1 <- c(10, 2, 1)
shape2 <- c(10, 1, 1)
get.f0 <- function(y) return((p0[1]*dbeta((y+1)/2,shape1[1],shape2[1])
                              + p0[2]*dbeta((y+1)/2,shape1[2],shape2[2])
                              + p0[3]*dbeta((y+1)/2,shape1[3],shape2[2]))/2)
f0 <- get.f0(y.grid)
h2.loss <- function(b){
  f.nc <- integrate(fFn, -1, 1, beta=b)$value
  return(integrate(function(y) (sqrt(fFn(y, beta=b)/f.nc) - sqrt(get.f0(y)))^2,-1,1)$value)
}
best.fit <- optim(rep(0,order), h2.loss, method="BFGS")
b0 <- best.fit$par
f0b <- get.f(b0)

###### get.miss and bde ######

get.miss.thin <- function(b, poiRate, ...){
  n.miss.cap <- rpois(1, poiRate)
  z.miss.cap <- runif(n.miss.cap,-1,1)
  Poly.miss.cap <- get.poly.mat(z.miss.cap)
  w.miss.cap <- rowSums(Poly.miss.cap %*% b)
  u.miss.cap <- (runif(n.miss.cap) < Phi(w.miss.cap))
  z.miss <- z.miss.cap[!u.miss.cap]
  n.miss <- length(z.miss)
  Poly.miss <- Poly.miss.cap[!u.miss.cap,]
  w.miss <- w.miss.cap[!u.miss.cap]
  return(list(Poly.miss=Poly.miss, w.miss=w.miss, n.miss=n.miss))
}
get.miss.reject <- function(b, n, ...){
  n.hits <- 0
  n.miss <- 0
  Poly.miss <- NULL
  w.miss <- NULL
  n.remain <- n
  while(n.remain > 0){
    z <- runif(n.remain, -1,1)
    pm <- matrix(yFn(z), nrow=n.remain)
    w <- rowSums(pm %*% b)
    u <- (runif(n.remain) < Phi(w))
    n.reject <- sum(!u)
    if(n.reject > 0){
      Poly.miss <- rbind(Poly.miss, pm[!u,])
      w.miss <- c(w.miss, w[!u])
    }
    n.remain <- n.reject
    n.miss <- n.miss + n.reject
  }
  return(list(Poly.miss=as.matrix(Poly.miss), w.miss=w.miss, n.miss=n.miss))
}

bde <- function(y, method=c("thin","reject")[1], b=rep(0,order), 
                nsamp=100, thin=1, burn=100, poiRate=4*length(y), 
                gam.pars = c(.01,.01), prec=1){
  n <- length(y)
  
  Poly.obs <- get.poly.mat(y)
  b.store <- matrix(NA, order, nsamp)
  Ntot.store <- rep(NA, nsamp)
  count.store <- 0
  
  get.miss <- ifelse(method=="thin", get.miss.thin, get.miss.reject)
  
  time.stamp.0 <- proc.time()
  for(iter in -(burn-1):(nsamp*thin)){
    w.obs <- c(Poly.obs %*% b)
    missing.stuff <- get.miss(b, n=n, poiRate=poiRate)
    n.miss <- missing.stuff$n.miss
    Poly.miss <- missing.stuff$Poly.miss
    w.miss <- missing.stuff$w.miss
    
    Ntot <- n + n.miss  
    poiRate <- rgamma(1, Ntot + gam.pars[1], 1 + gam.pars[2])
    
    Poly.design <- rbind(Poly.obs, Poly.miss)
    w.design <- c(w.obs, w.miss)
    u.design <- rep(c(1,0), c(n, n.miss))
    pg.draws <- pmax(1e-12, rpg.devroye(Ntot, h=1, z=w.design))
    pg.resp <- (u.design - 0.5)/pg.draws
    
    Xw <- Poly.design * sqrt(pg.draws)
    Yw <- pg.resp * sqrt(pg.draws)
    Xw.qr <- qr(Xw)
    Qw <- qr.Q(Xw.qr)
    Rw <- qr.R(Xw.qr)
    #b <- c(backsolve(Rw, crossprod(Qw, Yw) + rnorm(order)))
    Uw <- chol(crossprod(Rw) + diag(prec,order))
    b <- c(backsolve(Uw, backsolve(Uw, crossprod(Xw,Yw), transpose=TRUE) + rnorm(order)))
    if(any(is.na(b))) stop("NAs in b")
    
    if(iter > 0 & (iter %% thin == 0)){
      count.store <- count.store + 1
      b.store[,count.store] <- b
      Ntot.store[count.store] <- Ntot
    }
  }
  run.time <- proc.time() - time.stamp.0
  return(list(b = b.store, Ntot=Ntot.store, runtime=run.time))
}

###### get.miss.reject.x and bdreg ######

get.miss.reject.x <- function(b, n, x, a, shapes=c(1,1), jbw=0.16, positive=FALSE){
  n.hits <- 0
  n.miss <- 0
  Poly.miss <- NULL
  Half.miss <- NULL
  w.miss <- NULL
  x.miss <- NULL
  y.miss <- NULL
  
  xa <- c(x %*% a)
  if(positive) xa <- pmax(0, xa)
  
  n.remain <- n
  ix.remain <- 1:n
  while(n.remain > 0){
    z <- 2*rbeta(n.remain, shapes[1],shapes[2])-1
    pm <- matrix(yFn(z), nrow=n.remain)
    hk <- half.kern(z,jbw)
    x.remain <- x[ix.remain,,drop=FALSE]
    xa.remain <- xa[ix.remain]
    w <- rowSums(x.remain * (pm %*% b)) + xa.remain*hk
    u <- (runif(n.remain) < Phi(w))
    n.reject <- sum(!u)
    if(n.reject > 0){
      Poly.miss <- rbind(Poly.miss, pm[!u,])
      Half.miss <- c(Half.miss, hk[!u])
      w.miss <- c(w.miss, w[!u])
      x.miss <- rbind(x.miss, x.remain[!u,,drop=FALSE])
      y.miss <- c(y.miss, z[!u])
    }
    n.remain <- n.reject
    ix.remain <- ix.remain[!u]
    n.miss <- n.miss + n.reject
  }
  return(list(Poly.miss=as.matrix(Poly.miss), Half.miss=Half.miss, 
              w.miss=w.miss, n.miss=n.miss, x.miss=x.miss, y.miss=y.miss))
}

bdreg <- function(y, x=1, b=NULL, nsamp=100, thin=1, burn=100, prec=1){
  
  n <- length(y)
  if(is.null(dim(x))) x <- matrix(x, nrow=n)
  p <- ncol(x)
  if(is.null(b)) b <- replicate(p, rep(0,order))
  b <- matrix(b, order, p)
  if(length(prec) < p) prec <- rep(prec, p)[1:p]
  
  Poly.obs <- get.poly.mat(y)
  b.store <- array(NA, dim=c(order, p, nsamp))
  Ntot.store <- rep(NA, nsamp)
  count.store <- 0
  zeros <- rep(0,p)
  
  time.stamp.0 <- proc.time()
  for(iter in -(burn-1):(nsamp*thin)){
    w.obs <- rowSums(x * (Poly.obs %*% b))
    missing.stuff <- get.miss.reject.x(b, n, x, zeros)
    n.miss <- missing.stuff$n.miss
    Poly.miss <- missing.stuff$Poly.miss
    w.miss <- missing.stuff$w.miss
    x.miss <- missing.stuff$x.miss
    
    Ntot <- n + n.miss  
    
    Poly.comb <- rbind(Poly.obs, Poly.miss)
    w.comb <- c(w.obs, w.miss)
    x.comb <- rbind(x, x.miss)
    u.comb <- rep(c(1,0), c(n, n.miss))
    pg.draws <- pmax(1e-12, rpg.devroye(Ntot, h=1, z=w.comb))
    pg.resp <- (u.comb - 0.5)/pg.draws
    pg.residual <- pg.resp - w.comb
    
    for(j in 1:p){
      xPoly.j <- x.comb[,j] * Poly.comb
      pg.residual <- pg.residual + c(xPoly.j %*% b[,j])
      
      Yw <- pg.residual * sqrt(pg.draws)
      Xw <- xPoly.j * sqrt(pg.draws)
      Xw.qr <- qr(Xw)
      Qw <- qr.Q(Xw.qr)
      Rw <- qr.R(Xw.qr)
      Uw <- chol(crossprod(Rw) + diag(prec[j],order))
      b[,j] <- c(backsolve(Uw, backsolve(Uw, crossprod(Xw,Yw), transpose=TRUE) + rnorm(order)))
      if(any(is.na(b[,j]))) stop("NAs in b")
      pg.residual <- pg.residual - c(xPoly.j %*% b[,j])
    }
    
    if(iter > 0 & (iter %% thin == 0)){
      count.store <- count.store + 1
      b.store[,,count.store] <- b
      Ntot.store[count.store] <- Ntot
    }
  }
  run.time <- proc.time() - time.stamp.0
  return(list(b = b.store, Ntot=Ntot.store, runtime=run.time))
}

###### bdregjump ######

beta.loglik <- function(shapes, yy) return(sum(dbeta(yy, shapes[1], shapes[2], log=TRUE)))
beta.loglik2 <- function(bp, yy) beta.loglik(bp[2]*plogis(bp[1]*c(1,-1)), yy)
betapost.norm.approx <- function(y){
  mu.hat <- mean(y)
  mu.vec <- c(mu.hat, 1-mu.hat)
  g.prime.mu <- 1/prod(mu.vec)
  bt.hat <- qlogis(mu.hat)
  prec.hat <- max(0.5, prod(mu.vec)/var(y) - 1)
  tri.vec <- trigamma(prec.hat*mu.vec)
  w <- prec.hat * sum(tri.vec) / g.prime.mu^2
  cee <- -diff(tri.vec*mu.vec)*prec.hat
  dee <- sum(tri.vec*mu.vec^2) - trigamma(prec.hat)
  K.bb <- prec.hat*w
  K.pp <- dee
  K.bp <- cee/g.prime.mu
  K <- length(y) * matrix(c(K.bb,K.bp,K.bp,K.pp),2,2)
  R <- chol(solve(K))
  return(list(mean=c(bt.hat,prec.hat),half.var=R))
}
rprior.persistence <- function(n) return(rbeta(n,1,2))
dprior.persistence <- function(pers,...) return(dbeta(pers,1,2,...))
az.df <- 6
bp.df <- 6

### algorithms ###

bdregjump_neal8 <- function(y, x=1, b=NULL, nsamp=100, thin=1, burn=100,
                            shapes=c(1,1), prec=1,
                            jump=list(a=NULL, persistence=0.5, positive=TRUE,
                                      prec=1, ncand=10, update.jbw=TRUE, neal8=TRUE)){
  
  n <- length(y)
  if(is.null(dim(x))) x <- matrix(x, nrow=n)
  p <- ncol(x)
  if(is.null(b)) b <- replicate(p, rep(0,order))
  b <- matrix(b, order, p)
  if(length(prec) < p) prec <- rep(prec, p)[1:p]
  shapes.bp <- c(qlogis(shapes[1]/sum(shapes)), sum(shapes))
  
  
  a <- jump$a
  if(is.null(a)) a <- rep(0, p)
  persistence <- jump$persistence
  if(is.null(persistence)) persistence <- 0.5
  jbw <- 0.16*2*persistence
  pers.m <- 0.5; pers.s <- 3
  positive <- jump$positive
  if(is.null(positive)) positive <- FALSE
  hprec <- jump$prec
  if(is.null(hprec)) hprec <- 1
  ncand <- jump$ncand
  if(is.null(ncand)) ncand <- 10
  update.jbw <- jump$update.jbw
  if(is.null(update.jbw)) update.jbw <- TRUE
  neal8 <- jump$neal8
  if(is.null(neal8)) neal8 <- TRUE
  
  nu <- 3
  acand.mean <- rep(0,p)
  acand.var <- diag(1,p)
  acand.R <- chol(acand.var)
  
  Poly.obs <- get.poly.mat(y)
  b.store <- array(NA, dim=c(order, p, nsamp))
  a.store <- matrix(NA, p, nsamp)
  jbw.store <- rep(NA, nsamp)
  Ntot.store <- rep(NA, nsamp)
  shapes.store <- matrix(NA, 2, nsamp)
  count.store <- 0
  nacpt <- 0
  nacpt.bp <- 0
  
  time.stamp.0 <- proc.time()
  samp.count <- 1
  pb <- txtProgressBar(min = 0, max = burn+nsamp*thin, style = 3)
  for(iter in -(burn-1):(nsamp*thin)){
    samp.count <- samp.count + 1
    setTxtProgressBar(pb, samp.count)
    Half.obs <- half.kern(y, jbw)
    xa.obs <- c(x %*% a)
    if(positive) xa.obs <- pmax(0, xa.obs)
    
    w.obs <- rowSums(x * (Poly.obs %*% b)) + xa.obs * Half.obs
    
    missing.stuff <- get.miss.reject.x(b, n, x, a, shapes, jbw, positive)
    n.miss <- missing.stuff$n.miss
    Poly.miss <- missing.stuff$Poly.miss
    Half.miss <- missing.stuff$Half.miss
    w.miss <- missing.stuff$w.miss
    x.miss <- missing.stuff$x.miss
    y.miss <- missing.stuff$y.miss
    
    if(sum(y.miss==-1) > 0) y.miss[y.miss==-1] <- (-1 + 1e-10)
    if(sum(y.miss==1) > 0) y.miss[y.miss==1] <- (1 - 1e-10)
    
    Ntot <- n + n.miss  
    
    Poly.comb <- rbind(Poly.obs, Poly.miss)
    Half.comb <- c(Half.obs, Half.miss)
    w.comb <- c(w.obs, w.miss)
    x.comb <- rbind(x, x.miss)
    y.comb <- c(y, y.miss)
    
    u.comb <- rep(c(1,0), c(n, n.miss))
    pg.draws <- pmax(1e-12, rpg.devroye(Ntot, h=1, z=w.comb))
    pg.resp <- (u.comb - 0.5)/pg.draws
    pg.residual <- pg.resp - w.comb
    
    ## update b
    for(j in 1:p){
      xPoly.j <- x.comb[,j] * Poly.comb
      pg.residual <- pg.residual + c(xPoly.j %*% b[,j])
      
      Yw <- pg.residual * sqrt(pg.draws)
      Xw <- xPoly.j * sqrt(pg.draws)
      Uw <- chol(crossprod(Xw) + diag(prec[j],order))
      b[,j] <- c(backsolve(Uw, backsolve(Uw, crossprod(Xw,Yw), transpose=TRUE) + rnorm(order)))
      if(any(is.na(b[,j]))) stop("NAs in b")
      pg.residual <- pg.residual - c(xPoly.j %*% b[,j])
    }   
    ## update a
    xa <- c(x.comb %*% a)
    if(positive) xa <- pmax(0, xa)
    pg.residual <- pg.residual + xa * Half.comb
    
    if(neal8){
      az.cand <- replicate(ncand, rnorm(p, 0, sqrt(1/rgamma(1,nu/2,nu/2))))
      a.cand <- acand.mean + crossprod(acand.R, az.cand)
      a.cand[,1] <- a
      xa.cand <- x.comb %*% a.cand
      if(positive) xa.cand <- matrix(pmax(0, xa.cand), ncol=ncand)
      res.cand <- pg.residual - xa.cand * Half.comb

      lp.cand <- colSums(-0.5*pg.draws*res.cand^2) + 
        colSums(-0.5*hprec*a.cand^2) + 
        0.5*(nu+p)*log1p(colSums(az.cand^2)/nu)
      
      cand.pick <- sample(ncand, size=1, prob=exp(lp.cand-logsum(lp.cand)))
      nacpt <- nacpt+(cand.pick > 1)
      a <- a.cand[,cand.pick]
      xa <- xa.cand[,cand.pick]
      
      adapt.w <- 1/samp.count^(2/3)
      a.res <- a - acand.mean
      acand.mean <- acand.mean + adapt.w * a.res
      acand.var <- acand.var + adapt.w * (tcrossprod(a.res)-acand.var)
      acand.R <- chol(acand.var) 
      
    } else {
      xHalf <- x.comb*Half.comb ## nxp matrix
      Yw <- pg.residual * sqrt(pg.draws)
      Xw <- xHalf * sqrt(pg.draws)
      Uw <- chol(crossprod(Xw) + diag(hprec,p))
      mu.proto <- backsolve(Uw, crossprod(Xw,Yw), transpose=TRUE)
      if(positive){
        mu <- backsolve(Uw, mu.proto)
        xmu <- c(x%*%mu)
        ipos <- (xmu >= 0)
        Uw <- chol(crossprod(Xw[ipos,,drop=FALSE]) + diag(hprec,p))
        mu.proto <- backsolve(Uw, crossprod(Xw[ipos,,drop=FALSE],Yw[ipos]), transpose=TRUE)
      }      
      
      az.new <- rnorm(p) / sqrt(rgamma(1,az.df/2,az.df/2))
      a.new <- c(backsolve(Uw, mu.proto + az.new))
      xa.new <- c(x.comb %*% a.new)
      if(positive) xa.new <- pmax(0, xa.new)
      
      log.hastings.ratio <- 0
      if(positive){
        az <- c(Uw %*% a - mu.proto)
        
        res.new <- pg.residual - xa.new * Half.comb
        res <- pg.residual - xa * Half.comb
        log.hastings.ratio <- (sum(dnorm(res.new,0,sqrt(1/pg.draws),log=TRUE)
                                   - dnorm(res,0,sqrt(1/pg.draws),log=TRUE))
                               + sum(dnorm(a.new,0,sqrt(1/hprec),log=TRUE)
                                     - dnorm(a,0,sqrt(1/hprec),log=TRUE))
                               + 0.5*(az.df+p)*(log1p(sum(az.new^2)/az.df) - log1p(sum(az^2)/az.df)))
      }
      if(log(runif(1)) < log.hastings.ratio){
        a <- a.new
        xa <- xa.new
        nacpt <- nacpt+1
      }
    }
    pg.residual <- pg.residual - xa * Half.comb
    
    ## update jbw
    if(update.jbw){
      #persistence.cand <- rprior.persistence(ncand)
      persistence.cand <- rtt(ncand, center=pers.m, scale=pers.s, df=6, lo=0, hi=1)
      persistence.cand[1] <- persistence
      jbw.cand <- 2*0.16*persistence.cand
      Half.comb.cand <- sapply(jbw.cand, function(bw) half.kern(y.comb, bw))
      res.cand <- pg.residual + xa * Half.comb - xa * Half.comb.cand
      lp0 <- dnorm(pg.residual,0,sqrt(1/pg.draws),log=TRUE)
      lp.cand <- (
        apply(res.cand, 2, 
              function(res) sum(dnorm(res,0,sqrt(1/pg.draws),log=TRUE) - lp0))
        + dprior.persistence(persistence.cand,log=TRUE)
        - dtt(persistence.cand,center=pers.m,scale=pers.s,df=6,lo=0,hi=1,log=TRUE))
      cand.draw <- sample(ncand,1,prob=exp(lp.cand - logsum(lp.cand)))
      persistence <- persistence.cand[cand.draw]
      jbw <- jbw.cand[cand.draw]

      adapt.wt <- 1/samp.count^(2/3)
      pers.m <- (1-adapt.wt)*pers.m + adapt.wt*persistence
      pers.s <- sqrt((1-adapt.wt)*pers.s^2 + adapt.wt*(persistence-pers.m)^2)
      
    }    
    
    ## update beta shapes
    beta.y <- (1+y.comb)/2
    oo <- betapost.norm.approx(beta.y)
    bp.hat <- oo$mean
    bp.R <- oo$half.var
    bp.z <- backsolve(bp.R, shapes.bp - bp.hat, transpose=TRUE)
    cont <- TRUE
    while(cont){
      bp.znew <- rnorm(2)/sqrt(rgamma(1,bp.df/2,bp.df/2))
      bp.new <- bp.hat+c(crossprod(bp.R, bp.znew))
      cont <- bp.new[2] < 0
    }
    ll.diff <- beta.loglik2(bp.new, beta.y) - beta.loglik2(shapes.bp, beta.y)
    log.hastings <- ll.diff + 0.5*(2+bp.df)*(log1p(sum(bp.znew^2)/bp.df) - log1p(sum(bp.z^2)/bp.df))
    if(log(runif(1)) < log.hastings){
      shapes.bp <- bp.new
      nacpt.bp <- nacpt.bp+1
    }
    shapes <- shapes.bp[2]*plogis(shapes.bp[1]*c(-1,1))
    
    if(iter > 0 & (iter %% thin == 0)){
      count.store <- count.store + 1
      b.store[,,count.store] <- b
      a.store[,count.store] <- a
      jbw.store[count.store] <- jbw
      Ntot.store[count.store] <- Ntot
      shapes.store[,count.store] <- shapes
    }
  }
  run.time <- proc.time() - time.stamp.0
  return(list(b=b.store, a=a.store, jbw=jbw.store, Ntot=Ntot.store, 
              shapes=shapes.store, runtime=run.time, 
              acpt=c(a=nacpt, store=nacpt.bp)/(burn+nsamp*thin)))
}

bdregjump_adapt <- function(y, x=1, b=NULL, nsamp=100, thin=1, burn=100, 
                            shapes=c(1,1), prec=1, order = 5,
                            jump=list(a=NULL, persistence=0.5, positive=TRUE, 
                                      prec=1, ncand=10, update.jbw=TRUE),
                            adap=list(mu_adap=NULL, Sigma_adap=NULL,
                                      alpha_star = 0.234, r_adap = 2/3,
                                      chol.pivot=FALSE),
                            print.process = FALSE, print.i = FALSE){
  
  kbw <- 2/(order-1)
  knots <- seq(-1,1,kbw)
  bsd.sq <- (0.67*kbw)^2
  gausskern <- lapply(knots, function(mu) return(function(x) return(exp(-0.5*(x-mu)^2/bsd.sq))))
  get.poly.mat <- function(y) return(sapply(gausskern, function(f) f(y)))
  
  n <- length(y)
  if(is.null(dim(x))) x <- matrix(x, nrow=n)
  p <- ncol(x)
  if(is.null(b)) b <- replicate(p, rep(0,order))
  b <- matrix(b, order, p)
  if(length(prec) < p) prec <- rep(prec, p)[1:p]
  shapes.bp <- c(qlogis(shapes[1]/sum(shapes)), sum(shapes))
  
  a <- jump$a
  if(is.null(a)) a <- rep(0, p)
  persistence <- jump$persistence
  if(is.null(persistence)) persistence <- 0.5
  jbw <- 0.16*2*persistence
  positive <- jump$positive
  if(is.null(positive)) positive <- FALSE
  hprec <- jump$prec
  if(is.null(hprec)) hprec <- 1
  ncand <- jump$ncand
  if(is.null(ncand)) ncand <- 10
  update.jbw <- jump$update.jbw
  if(is.null(update.jbw)) update.jbw <- TRUE
  
  loglambda_adap <- log(2.38^2/p)
  
  mu_adap <- adap$mu_adap
  if(is.null(mu_adap)) mu_adap <- matrix(rep(0,p), ncol = 1)
  Sigma_adap <- adap$Sigma_adap
  if(is.null(Sigma_adap)) Sigma_adap <- diag(p)
  alpha_star <- adap$alpha_star
  if(is.null(alpha_star)) alpha_star <- 0.234
  r_adap <- adap$r_adap
  if(is.null(r_adap)) r_adap <- 2/3
  chol.pivot <- adap$chol.pivot
  if(is.null(chol.pivot)) chol.pivot <- FALSE
  
  Poly.obs <- get.poly.mat(y)
  b.store <- array(NA, dim=c(order, p, nsamp))
  a.store <- matrix(NA, p, nsamp)
  az.new.store <- matrix(NA, p, nsamp)
  jbw.store <- rep(NA, nsamp)
  Ntot.store <- rep(NA, nsamp)
  shapes.store <- matrix(NA, 2, nsamp)
  Sigma_adap.store <- array(NA, dim=c(p, p, nsamp))
  loglambda_adap.store <- rep(NA, nsamp)
  sample_i_count.store <- rep(NA, nsamp)
  p_accept.store <- rep(NA, nsamp)
  count.store <- 0
  nacpt <- 0
  nacpt.bp <- 0
  bp.df <- 6
  az.new <- rep(0,p)
  
  diff_prev <- burn_i_count <- sample_i_count <- 0
  
  ## get.miss.reject.x ##
  yFn <- get.poly.mat
  
  get.miss.reject.x <- function(b, n, x, a, shapes=c(1,1), jbw=0.16, positive=FALSE){
    n.hits <- 0
    n.miss <- 0
    Poly.miss <- NULL
    Half.miss <- NULL
    w.miss <- NULL
    x.miss <- NULL
    y.miss <- NULL
    
    xa <- c(x %*% a)
    if(positive) xa <- pmax(0, xa)
    
    n.remain <- n
    ix.remain <- 1:n
    while(n.remain > 0){
      z <- 2*rbeta(n.remain, shapes[1],shapes[2])-1
      pm <- matrix(yFn(z), nrow=n.remain)
      hk <- half.kern(z,jbw)
      x.remain <- x[ix.remain,,drop=FALSE]
      xa.remain <- xa[ix.remain]
      w <- rowSums(x.remain * (pm %*% b)) + xa.remain*hk
      u <- (runif(n.remain) < Phi(w))
      n.reject <- sum(!u)
      if(n.reject > 0){
        Poly.miss <- rbind(Poly.miss, pm[!u,])
        Half.miss <- c(Half.miss, hk[!u])
        w.miss <- c(w.miss, w[!u])
        x.miss <- rbind(x.miss, x.remain[!u,,drop=FALSE])
        y.miss <- c(y.miss, z[!u])
      }
      n.remain <- n.reject
      ix.remain <- ix.remain[!u]
      n.miss <- n.miss + n.reject
    }
    return(list(Poly.miss=as.matrix(Poly.miss), Half.miss=Half.miss, 
                w.miss=w.miss, n.miss=n.miss, x.miss=x.miss, y.miss=y.miss))
  }
  
  ##
  
  time.stamp.0 <- proc.time()
  for(iter in -(burn-1):(nsamp*thin)){
    
    Half.obs <- half.kern(y, jbw)
    xa.obs <- c(x %*% a)
    if(positive) xa.obs <- pmax(0, xa.obs)
    
    w.obs <- rowSums(x * (Poly.obs %*% b)) + xa.obs * Half.obs
    
    missing.stuff <- get.miss.reject.x(b, n, x, a, shapes, jbw, positive)
    n.miss <- missing.stuff$n.miss
    Poly.miss <- missing.stuff$Poly.miss
    Half.miss <- missing.stuff$Half.miss
    w.miss <- missing.stuff$w.miss
    x.miss <- missing.stuff$x.miss
    y.miss <- missing.stuff$y.miss
    
    if(sum(y.miss==-1) > 0) y.miss[y.miss==-1] <- (-1 + 1e-10)
    if(sum(y.miss==1) > 0) y.miss[y.miss==1] <- (1 - 1e-10)
    
    Ntot <- n + n.miss  
    
    Poly.comb <- rbind(Poly.obs, Poly.miss)
    Half.comb <- c(Half.obs, Half.miss)
    w.comb <- c(w.obs, w.miss)
    x.comb <- rbind(x, x.miss)
    y.comb <- c(y, y.miss)
    
    u.comb <- rep(c(1,0), c(n, n.miss))
    pg.draws <- pmax(1e-12, rpg.devroye(Ntot, h=1, z=w.comb))
    pg.resp <- (u.comb - 0.5)/pg.draws
    pg.residual <- pg.resp - w.comb
    
    ## update b
    for(j in 1:p){
      xPoly.j <- x.comb[,j] * Poly.comb
      pg.residual <- pg.residual + c(xPoly.j %*% b[,j])
      
      Yw <- pg.residual * sqrt(pg.draws)
      Xw <- xPoly.j * sqrt(pg.draws)
      Uw <- chol(crossprod(Xw) + diag(prec[j],order))
      b[,j] <- c(backsolve(Uw, backsolve(Uw, crossprod(Xw,Yw), transpose=TRUE) + rnorm(order)))
      if(any(is.na(b[,j]))) stop("NAs in b")
      pg.residual <- pg.residual - c(xPoly.j %*% b[,j])
    }
    
    ## update a
    xa <- c(x.comb %*% a)
    if(positive) xa <- pmax(0, xa)
    
    pg.residual <- pg.residual + xa * Half.comb
    
    Sigma_chol <- exp(loglambda_adap/2) * chol(Sigma_adap + diag(1e-10, p), pivot = chol.pivot)
    if(chol.pivot == TRUE){
      var_order <- attr(Sigma_chol, "pivot")
      az.new[var_order] <- crossprod(Sigma_chol, rnorm(p))
    } else{
      az.new <- crossprod(Sigma_chol, rnorm(p))
    }
    
    a.new <- a + az.new
    xa.new <- c(x.comb %*% a.new)
    if(positive) xa.new <- pmax(0, xa.new)

    res.new <- pg.residual - xa.new * Half.comb
    res <- pg.residual - xa * Half.comb
    log.hastings.ratio <- (sum(dnorm(res.new,0,sqrt(1/pg.draws),log=TRUE)
                               - dnorm(res,0,sqrt(1/pg.draws),log=TRUE))
                           + sum(dnorm(a.new,0,sqrt(1/hprec),log=TRUE)
                                 - dnorm(a,0,sqrt(1/hprec),log=TRUE)))
    
    if(log(runif(1)) < log.hastings.ratio){
      a <- a.new
      xa <- xa.new
      nacpt <- nacpt+1
    }
    p_accept <- min(1, exp(log.hastings.ratio))
    
    if(iter <= 0) {
      diff_curr <- p_accept - alpha_star
      burn_i_count <- burn_i_count + ((diff_prev * diff_curr) <= 0)
      gamma_adap <- min(0.5, 1/(burn_i_count^r_adap))
      diff_prev <- diff_curr
      
      if(print.i) print(burn_i_count)
    }
    else {
      diff_curr <- p_accept - alpha_star
      sample_i_count <- sample_i_count + ((diff_prev * diff_curr) <= 0)
      gamma_adap <- min(0.5, 1/(sample_i_count^r_adap))
      diff_prev <- diff_curr
      
      if(print.i) print(sample_i_count)
    }
    
    loglambda_adap <- loglambda_adap + gamma_adap * (p_accept - alpha_star)
    a_mu_diff <- a - mu_adap
    mu_adap <- mu_adap + gamma_adap * a_mu_diff
    Sigma_adap <- Sigma_adap + gamma_adap*(tcrossprod(a_mu_diff, a_mu_diff) - Sigma_adap)
    
    pg.residual <- pg.residual - xa * Half.comb
    
    ## update jbw
    if(update.jbw){
      persistence.cand <- rprior.persistence(ncand)
      jbw.cand <- 2*0.16*persistence.cand
      Half.comb.cand <- sapply(jbw.cand, function(bw) half.kern(y.comb, bw))
      res.cand <- pg.residual + xa * Half.comb - xa * Half.comb.cand
      lp0 <- dnorm(pg.residual,0,sqrt(1/pg.draws),log=TRUE)
      lp.cand <- c(apply(res.cand, 
                         2, 
                         function(res) sum(dnorm(res,0,sqrt(1/pg.draws),log=TRUE) - lp0)), 
                   0)
      cand.draw <- sample(ncand+1,1,prob=exp(lp.cand - logsum(lp.cand)))
      if(cand.draw < ncand+1) jbw <- jbw.cand[cand.draw]
    }    
    
    ## update beta shapes
    beta.y <- (1+y.comb)/2
    oo <- betapost.norm.approx(beta.y)
    bp.hat <- oo$mean
    bp.R <- oo$half.var
    bp.z <- backsolve(bp.R, shapes.bp - bp.hat, transpose=TRUE)
    cont <- TRUE
    while(cont){
      bp.znew <- rnorm(2)/sqrt(rgamma(1,bp.df/2,bp.df/2))
      bp.new <- bp.hat+c(crossprod(bp.R, bp.znew))
      cont <- bp.new[2] < 0
    }
    ll.diff <- beta.loglik2(bp.new, beta.y) - beta.loglik2(shapes.bp, beta.y)
    log.hastings <- ll.diff + 0.5*(2+bp.df)*(log1p(sum(bp.znew^2)/bp.df) - log1p(sum(bp.z^2)/bp.df))
    
    if(log(runif(1)) < log.hastings){
      shapes.bp <- bp.new
      nacpt.bp <- nacpt.bp+1
    }
    shapes <- shapes.bp[2]*plogis(shapes.bp[1]*c(-1,1))
    
    if(print.process) print(iter)
    
    if(iter > 0 & (iter %% thin == 0)){
      count.store <- count.store + 1
      b.store[,,count.store] <- b
      a.store[,count.store] <- a
      az.new.store[,count.store] <- az.new
      jbw.store[count.store] <- jbw
      Ntot.store[count.store] <- Ntot
      shapes.store[,count.store] <- shapes
      Sigma_adap.store[,,count.store] <- Sigma_adap
      loglambda_adap.store[count.store] <- loglambda_adap
      sample_i_count.store[count.store] <- sample_i_count
      p_accept.store[count.store] <- p_accept
    }
  }
  run.time <- proc.time() - time.stamp.0
  return(list(last_missing = missing.stuff,
              b=b.store, a=a.store, az.new = az.new.store, jbw=jbw.store, Ntot=Ntot.store,
              Sigma_adap = Sigma_adap.store, loglambda_adap = loglambda_adap.store,
              shapes=shapes.store, sample_i_count = sample_i_count.store,
              p_accept = p_accept.store,
              runtime=run.time, 
              acpt=c(a=nacpt, store=nacpt.bp)/(burn+nsamp*thin))
  )
}

bdregjump_adapt_all <- function(y, x=1, b=NULL, nsamp=100, thin=1, burn=100, 
                                shapes=c(1,1), prec=1, order=5,
                                jump=list(a=NULL, persistence=0.5, positive=TRUE, 
                                          prec=1, ncand=10),
                                adap=list(mu_adap=NULL, Sigma_adap=NULL,
                                          alpha_star = 0.234, r_adap = 2/3,
                                          chol.pivot=FALSE),
                                print.process = FALSE, print.i = FALSE){
  
  kbw <- 2/(order-1)
  knots <- seq(-1,1,kbw)
  bsd.sq <- (0.67*kbw)^2
  gausskern <- lapply(knots, function(mu) return(function(x) return(exp(-0.5*(x-mu)^2/bsd.sq))))
  get.poly.mat <- function(y) return(sapply(gausskern, function(f) f(y)))
  
  n <- length(y)
  if(is.null(dim(x))) x <- matrix(x, nrow=n)
  p <- ncol(x)
  if(is.null(b)) b <- replicate(p, rep(0,order))
  b <- matrix(b, order, p)
  if(length(prec) < p) prec <- rep(prec, p)[1:p]
  shapes.bp <- c(qlogis(shapes[1]/sum(shapes)), sum(shapes))
  
  a <- jump$a
  if(is.null(a)) a <- rep(0, p)
  persistence <- jump$persistence
  if(is.null(persistence)) persistence <- 0.5
  jbw <- 0.16*2*persistence
  log_jbw <- log(jbw)
  positive <- jump$positive
  if(is.null(positive)) positive <- FALSE
  hprec <- jump$prec
  if(is.null(hprec)) hprec <- 1
  ncand <- jump$ncand
  if(is.null(ncand)) ncand <- 10
  update.jbw <- jump$update.jbw
  if(is.null(update.jbw)) update.jbw <- TRUE
  a_logjbw <- c(a, log_jbw)
  
  loglambda_adap <- log(2.38^2/(p+1))
  
  mu_adap <- adap$mu_adap
  if(is.null(mu_adap)) mu_adap <- matrix(rep(0,p+1), ncol = 1)
  Sigma_adap <- adap$Sigma_adap
  if(is.null(Sigma_adap)) Sigma_adap <- diag(p+1)
  alpha_star <- adap$alpha_star
  if(is.null(alpha_star)) alpha_star <- 0.234
  r_adap <- adap$r_adap
  if(is.null(r_adap)) r_adap <- 2/3
  chol.pivot <- adap$chol.pivot
  if(is.null(chol.pivot)) chol.pivot <- FALSE
  
  Poly.obs <- get.poly.mat(y)
  b.store <- array(NA, dim=c(order, p, nsamp))
  a.store <- matrix(NA, p, nsamp)
  az.new.store <- matrix(NA, p, nsamp)
  jbw.store <- rep(NA, nsamp)
  Ntot.store <- rep(NA, nsamp)
  shapes.store <- matrix(NA, 2, nsamp)
  Sigma_adap.store <- array(NA, dim=c(p+1, p+1, nsamp))
  loglambda_adap.store <- rep(NA, nsamp)
  sample_i_count.store <- rep(NA, nsamp)
  p_accept.store <- rep(NA, nsamp)
  count.store <- 0
  nacpt <- 0
  nacpt.bp <- 0
  a_logjbw_z.new <- rep(0,p+1)
  
  diff_prev <- burn_i_count <- sample_i_count <- 0
  
  ## get.miss.reject.x ##
  yFn <- get.poly.mat
  
  get.miss.reject.x <- function(b, n, x, a, shapes=c(1,1), jbw=0.16, positive=FALSE){
    n.hits <- 0
    n.miss <- 0
    Poly.miss <- NULL
    Half.miss <- NULL
    w.miss <- NULL
    x.miss <- NULL
    y.miss <- NULL
    
    xa <- c(x %*% a)
    if(positive) xa <- pmax(0, xa)
    
    n.remain <- n
    ix.remain <- 1:n
    while(n.remain > 0){
      z <- 2*rbeta(n.remain, shapes[1],shapes[2])-1
      pm <- matrix(yFn(z), nrow=n.remain)
      hk <- half.kern(z,jbw)
      x.remain <- x[ix.remain,,drop=FALSE]
      xa.remain <- xa[ix.remain]
      w <- rowSums(x.remain * (pm %*% b)) + xa.remain*hk
      u <- (runif(n.remain) < Phi(w))
      n.reject <- sum(!u)
      if(n.reject > 0){
        Poly.miss <- rbind(Poly.miss, pm[!u,])
        Half.miss <- c(Half.miss, hk[!u])
        w.miss <- c(w.miss, w[!u])
        x.miss <- rbind(x.miss, x.remain[!u,,drop=FALSE])
        y.miss <- c(y.miss, z[!u])
      }
      n.remain <- n.reject
      ix.remain <- ix.remain[!u]
      n.miss <- n.miss + n.reject
    }
    return(list(Poly.miss=as.matrix(Poly.miss), Half.miss=Half.miss, 
                w.miss=w.miss, n.miss=n.miss, x.miss=x.miss, y.miss=y.miss))
  }
  
  ##
  
  time.stamp.0 <- proc.time()
  for(iter in -(burn-1):(nsamp*thin)){
    
    Half.obs <- half.kern(y, jbw)
    xa.obs <- c(x %*% a)
    if(positive) xa.obs <- pmax(0, xa.obs)
    
    w.obs <- rowSums(x * (Poly.obs %*% b)) + xa.obs * Half.obs
    
    missing.stuff <- get.miss.reject.x(b, n, x, a, shapes, jbw, positive)
    n.miss <- missing.stuff$n.miss
    Poly.miss <- missing.stuff$Poly.miss
    Half.miss <- missing.stuff$Half.miss
    w.miss <- missing.stuff$w.miss
    x.miss <- missing.stuff$x.miss
    y.miss <- missing.stuff$y.miss
    
    if(sum(y.miss==-1) > 0) y.miss[y.miss==-1] <- (-1 + 1e-10)
    if(sum(y.miss==1) > 0) y.miss[y.miss==1] <- (1 - 1e-10)
    
    Ntot <- n + n.miss  
    
    Poly.comb <- rbind(Poly.obs, Poly.miss)
    Half.comb <- c(Half.obs, Half.miss)
    w.comb <- c(w.obs, w.miss)
    x.comb <- rbind(x, x.miss)
    y.comb <- c(y, y.miss)
    
    u.comb <- rep(c(1,0), c(n, n.miss))
    pg.draws <- pmax(1e-12, rpg.devroye(Ntot, h=1, z=w.comb))
    pg.resp <- (u.comb - 0.5)/pg.draws
    pg.residual <- pg.resp - w.comb
    
    ## update b
    for(j in 1:p){
      xPoly.j <- x.comb[,j] * Poly.comb
      pg.residual <- pg.residual + c(xPoly.j %*% b[,j])
      
      Yw <- pg.residual * sqrt(pg.draws)
      Xw <- xPoly.j * sqrt(pg.draws)
      Uw <- chol(crossprod(Xw) + diag(prec[j],order))
      b[,j] <- c(backsolve(Uw, backsolve(Uw, crossprod(Xw,Yw), transpose=TRUE) + rnorm(order)))
      if(any(is.na(b[,j]))) stop("NAs in b")
      pg.residual <- pg.residual - c(xPoly.j %*% b[,j])
    }
    
    ## update a and log_jbw
    xa <- c(x.comb %*% a)
    if(positive) xa <- pmax(0, xa)
    
    pg.residual <- pg.residual + xa * Half.comb
    
    Sigma_chol <- exp(loglambda_adap/2) * chol(Sigma_adap + diag(1e-10, p+1), pivot = chol.pivot)
    if(chol.pivot == TRUE){
      var_order <- attr(Sigma_chol, "pivot")
      a_logjbw_z.new[var_order] <- crossprod(Sigma_chol, rnorm(p+1))
    } else{
      a_logjbw_z.new <- crossprod(Sigma_chol, rnorm(p+1))
    }
    
    az.new <- a_logjbw_z.new[-(p+1)]
    logjbw_z.new <- a_logjbw_z.new[p+1]
    
    a.new <- a + az.new
    log_jbw.new <- log_jbw + logjbw_z.new
    jbw.new <- exp(log_jbw.new)
    xa.new <- c(x.comb %*% a.new)
    if(positive) xa.new <- pmax(0, xa.new)
    
    log.hastings.ratio <- 0
    if(positive){
      Half.comb.new <- half.kern(y.comb, jbw.new)
      res.new <- pg.residual - xa.new * Half.comb.new
      res <- pg.residual - xa * Half.comb
      log.hastings.ratio <- (sum(dnorm(res.new,0,sqrt(1/pg.draws),log=TRUE)
                                 - dnorm(res,0,sqrt(1/pg.draws),log=TRUE))
                             + sum(dnorm(a.new,0,sqrt(1/hprec),log=TRUE)
                                   - dnorm(a,0,sqrt(1/hprec),log=TRUE))
                             + (log_jbw - log_jbw.new))
    }
    if(log(runif(1)) < log.hastings.ratio){
      a <- a.new
      xa <- xa.new
      jbw <- jbw.new
      nacpt <- nacpt+1
    }
    log_jbw <- log(jbw)
    a_logjbw <- c(a, log_jbw)
    p_accept <- min(1, exp(log.hastings.ratio))
    
    if(iter <= 0) {
      diff_curr <- p_accept - alpha_star
      burn_i_count <- burn_i_count + ((diff_prev * diff_curr) <= 0)
      gamma_adap <- min(0.5, 1/(burn_i_count^r_adap))
      diff_prev <- diff_curr
      
      if(print.i) print(burn_i_count)
    }
    else {
      diff_curr <- p_accept - alpha_star
      sample_i_count <- sample_i_count + ((diff_prev * diff_curr) <= 0)
      gamma_adap <- min(0.5, 1/(sample_i_count^r_adap))
      diff_prev <- diff_curr
      
      if(print.i) print(sample_i_count)
    }
    
    loglambda_adap <- loglambda_adap + gamma_adap * (p_accept - alpha_star)
    a_logjbw_mu_diff <- a_logjbw - mu_adap
    mu_adap <- mu_adap + gamma_adap * a_logjbw_mu_diff
    Sigma_adap <- Sigma_adap + gamma_adap*(tcrossprod(a_logjbw_mu_diff, a_logjbw_mu_diff) - Sigma_adap)
    
    ## update beta shapes
    beta.y <- (1+y.comb)/2
    oo <- betapost.norm.approx(beta.y)
    bp.hat <- oo$mean
    bp.R <- oo$half.var
    bp.z <- backsolve(bp.R, shapes.bp - bp.hat, transpose=TRUE)
    cont <- TRUE
    while(cont){
      bp.znew <- rnorm(2)/sqrt(rgamma(1,bp.df/2,bp.df/2))
      bp.new <- bp.hat+c(crossprod(bp.R, bp.znew))
      cont <- bp.new[2] < 0
    }
    ll.diff <- beta.loglik2(bp.new, beta.y) - beta.loglik2(shapes.bp, beta.y)
    log.hastings <- ll.diff + 0.5*(2+bp.df)*(log1p(sum(bp.znew^2)/bp.df) - log1p(sum(bp.z^2)/bp.df))
    if(log(runif(1)) < log.hastings){
      shapes.bp <- bp.new
      nacpt.bp <- nacpt.bp+1
    }
    shapes <- shapes.bp[2]*plogis(shapes.bp[1]*c(-1,1))
    
    if(print.process) print(iter)
    
    if(iter > 0 & (iter %% thin == 0)){
      count.store <- count.store + 1
      b.store[,,count.store] <- b
      a.store[,count.store] <- a
      az.new.store[,count.store] <- az.new
      jbw.store[count.store] <- jbw
      Ntot.store[count.store] <- Ntot
      shapes.store[,count.store] <- shapes
      Sigma_adap.store[,,count.store] <- Sigma_adap
      loglambda_adap.store[count.store] <- loglambda_adap
      sample_i_count.store[count.store] <- sample_i_count
      p_accept.store[count.store] <- p_accept
    }
  }
  run.time <- proc.time() - time.stamp.0
  return(list(last_missing = missing.stuff,
              b=b.store, a=a.store, az.new = az.new.store, jbw=jbw.store, Ntot=Ntot.store,
              Sigma_adap = Sigma_adap.store, loglambda_adap = loglambda_adap.store,
              shapes=shapes.store, sample_i_count = sample_i_count.store,
              p_accept = p_accept.store,
              runtime=run.time, 
              acpt=c(a=nacpt, store=nacpt.bp)/(burn+nsamp*thin))
  )
}



bdregjump_adapt_comp <- function(y, x=1, b=NULL, nsamp=100, thin=1, burn=100, 
                            shapes=c(1,1), prec=1, order = 5,
                            jump=list(a=NULL, persistence=0.5, positive=TRUE, 
                                      prec=1, ncand=10, update.jbw=TRUE),
                            adap=list(mu_adap=NULL, Sigma_adap=NULL,
                                      alpha_star = 0.234, r_adap = 2/3,
                                      chol.pivot=FALSE),
                            print.process = FALSE){
  
  kbw <- 2/(order-1)
  knots <- seq(-1,1,kbw)
  bsd.sq <- (0.67*kbw)^2
  gausskern <- lapply(knots, function(mu) return(function(x) return(exp(-0.5*(x-mu)^2/bsd.sq))))
  get.poly.mat <- function(y) return(sapply(gausskern, function(f) f(y)))
  
  n <- length(y)
  if(is.null(dim(x))) x <- matrix(x, nrow=n)
  p <- ncol(x)
  if(is.null(b)) b <- replicate(p, rep(0,order))
  b <- matrix(b, order, p)
  if(length(prec) < p) prec <- rep(prec, p)[1:p]
  shapes.bp <- c(qlogis(shapes[1]/sum(shapes)), sum(shapes))
  
  a <- jump$a
  if(is.null(a)) a <- rep(0, p)
  persistence <- jump$persistence
  if(is.null(persistence)) persistence <- 0.5
  jbw <- 0.16*2*persistence
  positive <- jump$positive
  if(is.null(positive)) positive <- FALSE
  hprec <- jump$prec
  if(is.null(hprec)) hprec <- 1
  ncand <- jump$ncand
  if(is.null(ncand)) ncand <- 10
  update.jbw <- jump$update.jbw
  if(is.null(update.jbw)) update.jbw <- TRUE
  
  loglambda_adap <- rep(log(2.38^2/p),p)
  
  mu_adap <- adap$mu_adap
  if(is.null(mu_adap)) mu_adap <- matrix(rep(0,p), ncol = 1)
  Sigma_adap <- adap$Sigma_adap
  if(is.null(Sigma_adap)) Sigma_adap <- diag(p)
  alpha_star <- adap$alpha_star
  if(is.null(alpha_star)) alpha_star <- 0.234
  r_adap <- adap$r_adap
  if(is.null(r_adap)) r_adap <- 2/3
  chol.pivot <- adap$chol.pivot
  if(is.null(chol.pivot)) chol.pivot <- FALSE
  
  Poly.obs <- get.poly.mat(y)
  b.store <- array(NA, dim=c(order, p, nsamp))
  a.store <- matrix(NA, p, nsamp)
  az.new.store <- matrix(NA, p, nsamp)
  jbw.store <- rep(NA, nsamp)
  Ntot.store <- rep(NA, nsamp)
  shapes.store <- matrix(NA, 2, nsamp)
  Sigma_adap.store <- array(NA, dim=c(p, p, nsamp))
  loglambda_adap.store <- matrix(NA, p, nsamp)
  p_accept.store <- rep(NA, nsamp)
  p_accept_corrd.store <- matrix(NA, p, nsamp)
  p_accept_corrd <- rep(NA, p)
  count.store <- 0
  nacpt <- 0
  nacpt.bp <- 0
  bp.df <- 6
  az.new <- rep(0,p)
  
  diff_prev <- burn_i_count <- sample_i_count <- 0
  
  ## get.miss.reject.x ##
  yFn <- get.poly.mat
  
  get.miss.reject.x <- function(b, n, x, a, shapes=c(1,1), jbw=0.16, positive=FALSE){
    n.hits <- 0
    n.miss <- 0
    Poly.miss <- NULL
    Half.miss <- NULL
    w.miss <- NULL
    x.miss <- NULL
    y.miss <- NULL
    
    xa <- c(x %*% a)
    if(positive) xa <- pmax(0, xa)
    
    n.remain <- n
    ix.remain <- 1:n
    while(n.remain > 0){
      z <- 2*rbeta(n.remain, shapes[1],shapes[2])-1
      pm <- matrix(yFn(z), nrow=n.remain)
      hk <- half.kern(z,jbw)
      x.remain <- x[ix.remain,,drop=FALSE]
      xa.remain <- xa[ix.remain]
      w <- rowSums(x.remain * (pm %*% b)) + xa.remain*hk
      u <- (runif(n.remain) < Phi(w))
      n.reject <- sum(!u)
      if(n.reject > 0){
        Poly.miss <- rbind(Poly.miss, pm[!u,])
        Half.miss <- c(Half.miss, hk[!u])
        w.miss <- c(w.miss, w[!u])
        x.miss <- rbind(x.miss, x.remain[!u,,drop=FALSE])
        y.miss <- c(y.miss, z[!u])
      }
      n.remain <- n.reject
      ix.remain <- ix.remain[!u]
      n.miss <- n.miss + n.reject
    }
    return(list(Poly.miss=as.matrix(Poly.miss), Half.miss=Half.miss, 
                w.miss=w.miss, n.miss=n.miss, x.miss=x.miss, y.miss=y.miss))
  }
  
  ##
  
  time.stamp.0 <- proc.time()
  for(iter in -(burn-1):(nsamp*thin)){
    iter_pos <- iter + burn
    
    Half.obs <- half.kern(y, jbw)
    xa.obs <- c(x %*% a)
    if(positive) xa.obs <- pmax(0, xa.obs)
    
    w.obs <- rowSums(x * (Poly.obs %*% b)) + xa.obs * Half.obs
    
    missing.stuff <- get.miss.reject.x(b, n, x, a, shapes, jbw, positive)
    n.miss <- missing.stuff$n.miss
    Poly.miss <- missing.stuff$Poly.miss
    Half.miss <- missing.stuff$Half.miss
    w.miss <- missing.stuff$w.miss
    x.miss <- missing.stuff$x.miss
    y.miss <- missing.stuff$y.miss
    
    if(sum(y.miss==-1) > 0) y.miss[y.miss==-1] <- (-1 + 1e-10)
    if(sum(y.miss==1) > 0) y.miss[y.miss==1] <- (1 - 1e-10)
    
    Ntot <- n + n.miss  
    
    Poly.comb <- rbind(Poly.obs, Poly.miss)
    Half.comb <- c(Half.obs, Half.miss)
    w.comb <- c(w.obs, w.miss)
    x.comb <- rbind(x, x.miss)
    y.comb <- c(y, y.miss)
    
    u.comb <- rep(c(1,0), c(n, n.miss))
    pg.draws <- pmax(1e-12, rpg.devroye(Ntot, h=1, z=w.comb))
    pg.resp <- (u.comb - 0.5)/pg.draws
    pg.residual <- pg.resp - w.comb
    
    ## update b
    for(j in 1:p){
      xPoly.j <- x.comb[,j] * Poly.comb
      pg.residual <- pg.residual + c(xPoly.j %*% b[,j])
      
      Yw <- pg.residual * sqrt(pg.draws)
      Xw <- xPoly.j * sqrt(pg.draws)
      Uw <- chol(crossprod(Xw) + diag(prec[j],order))
      b[,j] <- c(backsolve(Uw, backsolve(Uw, crossprod(Xw,Yw), transpose=TRUE) + rnorm(order)))
      if(any(is.na(b[,j]))) stop("NAs in b")
      pg.residual <- pg.residual - c(xPoly.j %*% b[,j])
    }
    
    ## update a
    gamma_adap <- min(0.5, 1/(iter_pos^r_adap))
    
    xa <- c(x.comb %*% a)
    if(positive) xa <- pmax(0, xa)
    
    pg.residual <- pg.residual + xa * Half.comb
    
    
    Sigma_chol <- chol(Sigma_adap + diag(1e-10, p), pivot = chol.pivot)
    if(chol.pivot == TRUE){
      var_order <- attr(Sigma_chol, "pivot")
      az.new <- crossprod(Sigma_chol[,order(var_order)] %*% diag(sqrt(exp(loglambda_adap))), rnorm(p))
    } else{
      az.new <- crossprod(Sigma_chol %*% diag(sqrt(exp(loglambda_adap))), rnorm(p))
    }
    
    a.new <- a + az.new
    xa.new <- c(x.comb %*% a.new)
    if(positive) xa.new <- pmax(0, xa.new)
    
    res <- pg.residual - xa * Half.comb
    
    for(adap_k in 1:p){
      az.new_k <- rep(0,p)
      az.new_k[adap_k] <- az.new[adap_k]
      a.new_k <- a + az.new_k
      xa.new_k <- c(x.comb %*% a.new_k)
      if(positive) xa.new_k <- pmax(0, xa.new_k)
      
      res.new_k <- pg.residual - xa.new_k * Half.comb
      
      log.hastings.ratio_k <- (sum(dnorm(res.new_k,0,sqrt(1/pg.draws),log=TRUE)
                                   - dnorm(res,0,sqrt(1/pg.draws),log=TRUE))
                               + sum(dnorm(a.new_k,0,sqrt(1/hprec),log=TRUE)
                                     - dnorm(a,0,sqrt(1/hprec),log=TRUE)))
      
      p_accept_k <- min(1, exp(log.hastings.ratio_k))
      loglambda_adap[adap_k] <- loglambda_adap[adap_k] + gamma_adap*(p_accept_k - alpha_star)
      p_accept_corrd[adap_k] <- p_accept_k
    }
    
    log.hastings.ratio <- 0
    if(positive){
      res.new <- pg.residual - xa.new * Half.comb
      log.hastings.ratio <- (sum(dnorm(res.new,0,sqrt(1/pg.draws),log=TRUE)
                                 - dnorm(res,0,sqrt(1/pg.draws),log=TRUE))
                             + sum(dnorm(a.new,0,sqrt(1/hprec),log=TRUE)
                                   - dnorm(a,0,sqrt(1/hprec),log=TRUE)))
    }
    if(log(runif(1)) < log.hastings.ratio){
      a <- a.new
      xa <- xa.new
      nacpt <- nacpt+1
    }
    p_accept <- min(1, exp(log.hastings.ratio))
    
    a_mu_diff <- a - mu_adap
    mu_adap <- mu_adap + gamma_adap * a_mu_diff
    Sigma_adap <- Sigma_adap + gamma_adap*(tcrossprod(a_mu_diff, a_mu_diff) - Sigma_adap)
    
    pg.residual <- pg.residual - xa * Half.comb
    
    ## update jbw
    if(update.jbw){
      persistence.cand <- rprior.persistence(ncand)
      jbw.cand <- 2*0.16*persistence.cand
      Half.comb.cand <- sapply(jbw.cand, function(bw) half.kern(y.comb, bw))
      res.cand <- pg.residual + xa * Half.comb - xa * Half.comb.cand
      lp0 <- dnorm(pg.residual,0,sqrt(1/pg.draws),log=TRUE)
      lp.cand <- c(apply(res.cand, 
                         2, 
                         function(res) sum(dnorm(res,0,sqrt(1/pg.draws),log=TRUE) - lp0)), 
                   0)
      cand.draw <- sample(ncand+1,1,prob=exp(lp.cand - logsum(lp.cand)))
      if(cand.draw < ncand+1) jbw <- jbw.cand[cand.draw]
    }    
    
    ## update beta shapes
    beta.y <- (1+y.comb)/2
    oo <- betapost.norm.approx(beta.y)
    bp.hat <- oo$mean
    bp.R <- oo$half.var
    bp.z <- backsolve(bp.R, shapes.bp - bp.hat, transpose=TRUE)
    cont <- TRUE
    while(cont){
      bp.znew <- rnorm(2)/sqrt(rgamma(1,bp.df/2,bp.df/2))
      bp.new <- bp.hat+c(crossprod(bp.R, bp.znew))
      cont <- bp.new[2] < 0
    }
    ll.diff <- beta.loglik2(bp.new, beta.y) - beta.loglik2(shapes.bp, beta.y)
    log.hastings <- ll.diff + 0.5*(2+bp.df)*(log1p(sum(bp.znew^2)/bp.df) - log1p(sum(bp.z^2)/bp.df))
    
    if(log(runif(1)) < log.hastings){
      shapes.bp <- bp.new
      nacpt.bp <- nacpt.bp+1
    }
    shapes <- shapes.bp[2]*plogis(shapes.bp[1]*c(-1,1))
    
    if(print.process) print(iter)
    
    if(iter > 0 & (iter %% thin == 0)){
      count.store <- count.store + 1
      b.store[,,count.store] <- b
      a.store[,count.store] <- a
      az.new.store[,count.store] <- az.new
      jbw.store[count.store] <- jbw
      Ntot.store[count.store] <- Ntot
      shapes.store[,count.store] <- shapes
      Sigma_adap.store[,,count.store] <- Sigma_adap
      loglambda_adap.store[,count.store] <- loglambda_adap
      p_accept.store[count.store] <- p_accept
      p_accept_corrd.store[,count.store] <- p_accept_corrd
    }
  }
  run.time <- proc.time() - time.stamp.0
  return(list(last_missing = missing.stuff,
              b=b.store, a=a.store, az.new = az.new.store, jbw=jbw.store, Ntot=Ntot.store,
              Sigma_adap = Sigma_adap.store, loglambda_adap = loglambda_adap.store,
              shapes=shapes.store,
              p_accept = p_accept.store, p_accept_corrd = p_accept_corrd.store,
              runtime=run.time, 
              acpt=c(a=nacpt, store=nacpt.bp)/(burn+nsamp*thin))
  )
}

bdregjump <- ifelse(adapt, bdregjump_adapt, bdregjump_neal8)


############################# Tempered Model ##################################

bdregjump_adapt_tempered <- function(y, x=1, b=NULL, nsamp=100, thin=1, burn=100, 
                                     shapes=c(1,1), prec=1, order = 15, 
                                     jump=list(a=NULL, persistence=0.5, positive=TRUE,
                                               prec=1, ncand=10, update.jbw=TRUE),
                                     adap=list(mu_adap=NULL, Sigma_adap=NULL,
                                               rho_adap = 1,
                                               alpha_star = 0.234, r_adap = 2/3,
                                               chol.pivot=FALSE, 
                                               parallel = 5, alpha_star_PT = 0.234),
                                     print.process = FALSE, print.i = FALSE){
  
  L <- adap$parallel
  if(is.null(L)) L <- 5
  
  kbw <- 2/(order-1)
  knots <- seq(-1,1,kbw)
  bsd.sq <- (0.67*kbw)^2
  gausskern <- lapply(knots, function(mu) return(function(x) return(exp(-0.5*(x-mu)^2/bsd.sq))))
  get.poly.mat <- function(y) return(sapply(gausskern, function(f) f(y)))
  
  n <- length(y)
  if(is.null(dim(x))) x <- matrix(x, nrow=n)
  p <- ncol(x)
  if(is.null(b)) b <- replicate(p, rep(0,order))
  b <- matrix(b, order, p)
  b_PT <- replicate(L, b)
  if(length(prec) < p) prec <- rep(prec, p)[1:p]
  shapes.bp <- c(qlogis(shapes[1]/sum(shapes)), sum(shapes))
  shapes_PT <- replicate(L, shapes)
  
  a <- jump$a
  a_PT <- replicate(L, a)
  if(is.null(a)){
    a <- rep(0, p)
    a_PT <- replicate(L, rep(0, p))
  }
  persistence <- jump$persistence
  if(is.null(persistence)) persistence <- 0.5
  jbw <- 0.16*2*persistence
  jbw_PT <- rep(jbw, L)
  positive <- jump$positive
  if(is.null(positive)) positive <- TRUE
  hprec <- jump$prec
  if(is.null(hprec)) hprec <- 1
  ncand <- jump$ncand
  if(is.null(ncand)) ncand <- 10
  update.jbw <- jump$update.jbw
  if(is.null(update.jbw)) update.jbw <- TRUE
  
  loglambda_adap_PT <- rep(log(2.38^2/p), L)
  
  mu_adap <- adap$mu_adap
  mu_adap_PT <- replicate(L, mu_adap)
  if(is.null(mu_adap)){
    mu_adap <- rep(0,p)
    mu_adap_PT <- replicate(L, rep(0,p))
  }
  Sigma_adap <- adap$Sigma_adap
  Sigma_adap_PT <- replicate(L, Sigma_adap)
  if(is.null(Sigma_adap)){
    Sigma_adap <- diag(p)
    Sigma_adap_PT <- replicate(L, diag(p))
  } 
  
  rho_adap <- adap$rho_adap
  rho_adap_PT <- rep(rho_adap, L-1)
  if(is.null(rho_adap)){
    rho_adap <- 1
    rho_adap_PT <- rep(1, L-1)
  }
  
  alpha_star <- adap$alpha_star
  if(is.null(alpha_star)) alpha_star <- 0.234
  alpha_star_PT <- adap$alpha_star_PT
  if(is.null(alpha_star_PT)) alpha_star_PT <- 0.234
  r_adap <- adap$r_adap
  if(is.null(r_adap)) r_adap <- 2/3
  chol.pivot <- adap$chol.pivot
  if(is.null(chol.pivot)) chol.pivot <- FALSE
  
  log.likelihood_PT <- rep(NA, L)
  
  Poly.obs <- get.poly.mat(y)
  b.store <- array(NA, dim=c(order, p, nsamp))
  a_PT.store <- array(NA, dim=c(p, L, nsamp))
  jbw.store <- rep(NA, nsamp)
  Ntot.store <- rep(NA, nsamp)
  shapes.store <- matrix(NA, 2, nsamp)
  Sigma_adap.store <- array(NA, dim=c(p, p, nsamp))
  loglambda_adap.store <- rep(NA, nsamp)
  p_accept.store <- rep(NA, nsamp)
  temp_PT.store <- matrix(NA, L, nsamp)
  log.likelihood_PT.store <- matrix(NA, L, nsamp)
  H_l_all <- rep(NA, L-1)
  H_l.store <- matrix(NA, L-1, nsamp)
  rho_adap_PT.store <- matrix(NA, L-1, nsamp)
  
  swap_info.store <- matrix(NA, 3, nsamp)
  
  count.store <- 0
  nacpt <- 0
  nacpt.bp <- 0
  nswap <- 0
  swap_pos <- NULL
  swap_accept <- NULL
  bp.df <- 6
  az.new <- rep(0,p)
  
  config_info_all <- list()
  swap_accept_prob <- rep(NA, L-1)
  swap_accept_prob.store <- matrix(NA, L-1, nsamp)
  
  ################### get.miss.reject.x ###################
  yFn <- get.poly.mat
  
  get.miss.reject.x <- function(b, n, x, a, shapes=c(1,1), jbw=0.16, positive=TRUE,
                                temp){
    n.hits <- 0
    n.miss <- 0
    Poly.miss <- NULL
    Half.miss <- NULL
    w.miss <- NULL
    x.miss <- NULL
    y.miss <- NULL
    
    xa <- c(x %*% a)
    if(positive) xa <- pmax(0, xa)
    
    n.remain <- n
    ix.remain <- 1:n
    while(n.remain > 0){
      z <- 2*rbeta(n.remain, shapes[1],shapes[2])-1
      pm <- matrix(yFn(z), nrow=n.remain)
      hk <- half.kern(z,jbw)
      x.remain <- x[ix.remain,,drop=FALSE]
      xa.remain <- xa[ix.remain]
      w <- temp * (rowSums(x.remain * (pm %*% b)) + xa.remain*hk)
      u <- (runif(n.remain) < Phi(w))
      n.reject <- sum(!u)
      if(n.reject > 0){
        Poly.miss <- rbind(Poly.miss, pm[!u,])
        Half.miss <- c(Half.miss, hk[!u])
        w.miss <- c(w.miss, w[!u])
        x.miss <- rbind(x.miss, x.remain[!u,,drop=FALSE])
        y.miss <- c(y.miss, z[!u])
      }
      n.remain <- n.reject
      ix.remain <- ix.remain[!u]
      n.miss <- n.miss + n.reject
    }
    return(list(Poly.miss=as.matrix(Poly.miss), Half.miss=Half.miss, 
                w.miss=w.miss, n.miss=n.miss, x.miss=x.miss, y.miss=y.miss))
  }
  
  ################### random_walk_step ###################
  
  random_walk_step <- function(y, x, prec, hprec, n, ncand, temp, Poly.obs,
                               positive = TRUE, alpha_star, r_adap, chol.pivot, 
                               b, a, shapes, jbw, update.jbw,
                               gamma_adap, loglambda_adap, mu_adap, Sigma_adap){
    nacpt <- 0
    nacpt.bp <- 0
    
    Half.obs <- half.kern(y, jbw)
    xa.obs <- c(x %*% a)
    if(positive) xa.obs <- pmax(0, xa.obs)
    
    w.obs <- temp * (rowSums(x * (Poly.obs %*% b)) + xa.obs * Half.obs)
    
    missing.stuff <- get.miss.reject.x(b, n, x, a, shapes, jbw, positive, temp)
    n.miss <- missing.stuff$n.miss
    Poly.miss <- missing.stuff$Poly.miss
    Half.miss <- missing.stuff$Half.miss
    w.miss <- missing.stuff$w.miss
    x.miss <- missing.stuff$x.miss
    y.miss <- missing.stuff$y.miss
    
    if(sum(y.miss==-1) > 0) y.miss[y.miss==-1] <- (-1 + 1e-10)
    if(sum(y.miss==1) > 0) y.miss[y.miss==1] <- (1 - 1e-10)
    
    Ntot <- n + n.miss  
    
    Poly.comb <- rbind(Poly.obs, Poly.miss)
    Half.comb <- c(Half.obs, Half.miss)
    w.comb <- c(w.obs, w.miss)
    x.comb <- rbind(x, x.miss)
    y.comb <- c(y, y.miss)
    
    u.comb <- rep(c(1,0), c(n, n.miss))
    pg.draws <- pmax(1e-12, rpg.devroye(Ntot, h=1, z=w.comb))
    pg.resp <- (u.comb - 0.5)/pg.draws
    pg.residual <- pg.resp - w.comb
    
    ## update b
    for(j in 1:p){
      xPoly.j <- temp * x.comb[,j] * Poly.comb
      pg.residual <- pg.residual + c(xPoly.j %*% b[,j])
      
      Yw <- pg.residual * sqrt(pg.draws)
      Xw <- xPoly.j * sqrt(pg.draws)
      Uw <- chol(crossprod(Xw) + diag(prec[j],order))
      b[,j] <- c(backsolve(Uw, backsolve(Uw, crossprod(Xw,Yw), transpose=TRUE) + rnorm(order)))
      if(any(is.na(b[,j]))) stop("NAs in b")
      pg.residual <- pg.residual - c(xPoly.j %*% b[,j])
    }
    
    ## update a
    xa <- temp * c(x.comb %*% a)
    if(positive) xa <- pmax(0, xa)
    
    pg.residual <- pg.residual + xa * Half.comb
    
    Sigma_chol <- exp(loglambda_adap/2) * chol(Sigma_adap + diag(1e-10, p), pivot = chol.pivot)
    if(chol.pivot == TRUE){
      var_order <- attr(Sigma_chol, "pivot")
      az.new[var_order] <- crossprod(Sigma_chol, rnorm(p))
    } else{
      az.new <- crossprod(Sigma_chol, rnorm(p))
    }
    
    a.new <- a + az.new
    xa.new <- temp * c(x.comb %*% a.new)
    if(positive) xa.new <- pmax(0, xa.new)
    
    res.new <- pg.residual - xa.new * Half.comb
    res <- pg.residual - xa * Half.comb
    log.hastings.ratio <- ( sum(dnorm(res.new,0,sqrt(1/pg.draws),log=TRUE)
                                - dnorm(res,0,sqrt(1/pg.draws),log=TRUE))
                            + sum(dnorm(a.new,0,sqrt(1/hprec),log=TRUE)
                                  - dnorm(a,0,sqrt(1/hprec),log=TRUE)) )
    
    if(log(runif(1)) < log.hastings.ratio){
      a <- a.new
      xa <- xa.new
      nacpt <- 1
      log.likelihood <- sum(dnorm(res.new,0,sqrt(1/pg.draws),log=TRUE))
    } else{
      log.likelihood <- sum(dnorm(res,0,sqrt(1/pg.draws),log=TRUE))
    }
    p_accept <- min(1, exp(log.hastings.ratio))
    
    loglambda_adap <- loglambda_adap + gamma_adap * (p_accept - alpha_star)
    a_mu_diff <- a - mu_adap
    mu_adap <- mu_adap + gamma_adap * a_mu_diff
    Sigma_adap <- Sigma_adap + gamma_adap*(tcrossprod(a_mu_diff, a_mu_diff) - Sigma_adap)
    
    pg.residual <- pg.residual - xa * Half.comb
    
    ## update jbw
    if(update.jbw){
      persistence.cand <- rprior.persistence(ncand)
      jbw.cand <- 2*0.16*persistence.cand
      Half.comb.cand <- sapply(jbw.cand, function(bw) half.kern(y.comb, bw))
      res.cand <- pg.residual + temp * (xa * Half.comb - xa * Half.comb.cand)
      lp0 <- dnorm(pg.residual,0,sqrt(1/pg.draws),log=TRUE)
      lp.cand <- c(apply(res.cand, 
                         2, 
                         function(res) sum(dnorm(res,0,sqrt(1/pg.draws),log=TRUE) - lp0)), 
                   0)
      cand.draw <- sample(ncand+1,1,prob=exp(lp.cand - logsum(lp.cand)))
      if(cand.draw < ncand+1) jbw <- jbw.cand[cand.draw]
    }    
    
    ## update beta shapes
    beta.y <- (1+y.comb)/2
    oo <- betapost.norm.approx(beta.y)
    bp.hat <- oo$mean
    bp.R <- oo$half.var
    bp.z <- backsolve(bp.R, shapes.bp - bp.hat, transpose=TRUE)
    cont <- TRUE
    while(cont){
      bp.znew <- rnorm(2)/sqrt(rgamma(1,bp.df/2,bp.df/2))
      bp.new <- bp.hat+c(crossprod(bp.R, bp.znew))
      cont <- bp.new[2] < 0
    }
    ll.diff <- beta.loglik2(bp.new, beta.y) - beta.loglik2(shapes.bp, beta.y)
    log.hastings <- ll.diff + 0.5*(2+bp.df)*(log1p(sum(bp.znew^2)/bp.df) - log1p(sum(bp.z^2)/bp.df))
    
    if(log(runif(1)) < log.hastings){
      shapes.bp <- bp.new
      nacpt.bp <- 1
    }
    shapes <- shapes.bp[2]*plogis(shapes.bp[1]*c(-1,1))
    
    config_info <- list(x.miss = x.miss, y.miss = y.miss)
    
    return(list(b = b, a = a, shapes = shapes, jbw = jbw, log.likelihood = log.likelihood,
                loglambda_adap = loglambda_adap, mu_adap = mu_adap, Sigma_adap = Sigma_adap,
                Ntot = Ntot, p_accept = p_accept, nacpt = nacpt, nacpt.bp = nacpt.bp,
                config_info = config_info))
  }
  
  ################### swap likelihood ###################
  
  get.swap.prob <-  function(temp1, temp2, x, y, x.miss1, y.miss1, x.miss2, y.miss2,
                             jbw1, jbw2, b1, b2, a1, a2, shapes1, shapes2, positive=TRUE){
    
    get.log.likelihood <- function(temp, x, y, x.miss, y.miss,
                                   jbw, b, a, shapes, positive = TRUE){
      y_b.miss <- (y.miss+1)/2
      
      Half.obs <- half.kern(y, jbw)
      Half.miss <- half.kern(y.miss, jbw)
      
      Poly.obs <- get.poly.mat(y)
      Poly.miss <- get.poly.mat(y.miss)
      
      xa.obs <- c(x %*% a)
      xa.miss <- c(x.miss %*% a)
      if(positive){
        xa.obs <- pmax(0, xa.obs)
        xa.miss <- pmax(0, xa.miss)
      }
      w.obs <- temp * (rowSums(x * (Poly.obs %*% b)) + xa.obs * Half.obs)
      w.miss <- temp * (rowSums(x.miss * (Poly.miss %*% b)) + xa.miss * Half.miss)
      
      log.likelihood <- sum(dbeta(y_b.miss, shapes[1], shapes[2], log = TRUE)) +
        sum(plogis(w.obs, log.p = TRUE)) +
        sum(plogis(w.miss, lower.tail = FALSE, log.p = TRUE))
      return(log.likelihood)
    }
    
    loglike12 <- get.log.likelihood(temp=temp1, x=x, y=y, x.miss=x.miss2, y.miss=y.miss2,
                                    jbw=jbw2, b=b2, a=a2, shapes=shapes2, positive = positive)
    loglike21 <- get.log.likelihood(temp=temp2, x=x, y=y, x.miss=x.miss1, y.miss=y.miss1,
                                    jbw=jbw1, b=b1, a=a1, shapes=shapes1, positive = positive)
    loglike11 <- get.log.likelihood(temp=temp1, x=x, y=y, x.miss=x.miss1, y.miss=y.miss1,
                                    jbw=jbw1, b=b1, a=a1, shapes=shapes1, positive = positive)
    loglike22 <- get.log.likelihood(temp=temp2, x=x, y=y, x.miss=x.miss2, y.miss=y.miss2,
                                    jbw=jbw2, b=b2, a=a2, shapes=shapes2, positive = positive)
    
    swap.prob <- min(1, exp(loglike12 + loglike21 - loglike11 - loglike22))
    return(swap.prob)
  }
  
  ########################################################################
  
  time.stamp.0 <- proc.time()
  for(iter in -(burn-1):(nsamp*thin)){
    
    iter_pos <- iter + burn
    gamma_adap <- min(0.5, 1/(iter_pos^r_adap))
    temp_PT <- 1 / cumsum(exp(c(0, rho_adap_PT)))
    
    for(chain_ind in 1:L){
      random_walk_result <- random_walk_step(y=y, x=x, prec=prec, hprec=hprec, Poly.obs = Poly.obs,
                                             n=n, ncand=ncand, temp = temp_PT[chain_ind],
                                             positive=positive, alpha_star=alpha_star,
                                             r_adap=r_adap, chol.pivot=chol.pivot, 
                                             b=b_PT[,,chain_ind], a=a_PT[,chain_ind],
                                             shapes=shapes_PT[,chain_ind], 
                                             jbw=jbw_PT[chain_ind], update.jbw=update.jbw,
                                             gamma_adap=gamma_adap,
                                             loglambda_adap=loglambda_adap_PT[chain_ind],
                                             mu_adap=mu_adap_PT[,chain_ind], 
                                             Sigma_adap=Sigma_adap_PT[,,chain_ind])
      b_PT[,,chain_ind] <- random_walk_result$b
      a_PT[,chain_ind] <- random_walk_result$a
      shapes_PT[,chain_ind] <- random_walk_result$shapes
      jbw_PT[chain_ind] <- random_walk_result$jbw
      log.likelihood_PT[chain_ind] <- random_walk_result$log.likelihood
      loglambda_adap_PT[chain_ind] <- random_walk_result$loglambda_adap
      mu_adap_PT[,chain_ind] <- random_walk_result$mu_adap
      Sigma_adap_PT[,,chain_ind] <- random_walk_result$Sigma_adap
      
      config_info_all[[chain_ind]] <- random_walk_result$config_info
      
      if(chain_ind == 1){
        Ntot <- random_walk_result$Ntot
        p_accept <- random_walk_result$p_accept
        nacpt <- nacpt + random_walk_result$nacpt
        nacpt.bp <- nacpt.bp + random_walk_result$nacpt.bp
      }
    }
    
    if(iter > 0 & (iter %% thin == 0)){
      count.store <- count.store + 1
      b.store[,,count.store] <- b_PT[,,1]
      a_PT.store[,,count.store] <- a_PT
      jbw.store[count.store] <- jbw_PT[1]
      Ntot.store[count.store] <- Ntot
      shapes.store[,count.store] <- shapes_PT[,1]
      Sigma_adap.store[,,count.store] <- Sigma_adap_PT[,,1]
      loglambda_adap.store[count.store] <- loglambda_adap_PT[1]
      p_accept.store[count.store] <- p_accept
      temp_PT.store[,count.store] <- temp_PT
      log.likelihood_PT.store[,count.store] <- log.likelihood_PT
    }
    
    ## likelihood ratio
    
    for(l in 1:(L-1)){
      config_info1 <- config_info_all[[l]]
      config_info2 <- config_info_all[[l+1]]
      swap_accept_prob[l] <- get.swap.prob(temp1=temp_PT[l], temp2=temp_PT[l+1],
                                           x=x, y=y, x.miss1=config_info1$x.miss,
                                           y.miss1=config_info1$y.miss,
                                           x.miss2=config_info2$x.miss,
                                           y.miss2=config_info2$y.miss,
                                           jbw1=jbw_PT[l], jbw2=jbw_PT[l+1],
                                           b1=b_PT[,,l], b2=b_PT[,,(l+1)],
                                           a1=a_PT[,l], a2=a_PT[,(l+1)],
                                           shapes1=shapes_PT[,l],
                                           shapes2=shapes_PT[,(l+1)], positive)
    }
    
    ## swap
    
    swap_ind <- sample(1:(L-1), 1)
    swap_tf <- FALSE
    swap_accept <- swap_accept_prob[swap_ind]
    
    if(runif(1) < swap_accept){
      swap_tf <- TRUE
      
      b_PT_s <- b_PT[,,swap_ind]
      a_PT_s <- a_PT[,swap_ind]
      jbw_PT_s <- jbw_PT[swap_ind]
      shapes_PT_s <- shapes_PT[,swap_ind]
      
      b_PT[,,swap_ind] <- b_PT[,,swap_ind+1]
      a_PT[,swap_ind] <- a_PT[,swap_ind+1]
      jbw_PT[swap_ind] <- jbw_PT[swap_ind+1]
      shapes_PT[,swap_ind] <- shapes_PT[,swap_ind+1]
      
      b_PT[,,swap_ind+1] <- b_PT_s
      a_PT[,swap_ind+1] <- a_PT_s
      jbw_PT[swap_ind+1] <- jbw_PT_s
      shapes_PT[,swap_ind+1] <- shapes_PT_s
      
      nswap <- nswap+1
      if(swap_ind == 1){
        swap_pos <- c(swap_pos, iter)
      }
    }
    
    ## update rho
    
    for(l in 1:(L-1)){
      H_l <- swap_accept_prob[l] - alpha_star_PT
      H_l_all[l] <- H_l
      rho_adap_PT[l] <- rho_adap_PT[l] + gamma_adap * H_l
    }
    
    if(iter > 0 & (iter %% thin == 0)){
      swap_info.store[1,count.store] <- swap_ind
      swap_info.store[2,count.store] <- swap_accept
      swap_info.store[3,count.store] <- swap_tf
      H_l.store[,count.store] <- H_l_all
      rho_adap_PT.store[,count.store] <- rho_adap_PT
      swap_accept_prob.store[,count.store] <- swap_accept_prob
    }
    
    if(print.process) print(iter)
    
  }
  
  run.time <- proc.time() - time.stamp.0
  return(list(b=b.store, a_PT=a_PT.store, jbw=jbw.store, Ntot=Ntot.store,
              Sigma_adap = Sigma_adap.store, loglambda_adap = loglambda_adap.store,
              shapes=shapes.store, p_accept = p_accept.store,
              temp = temp_PT.store, swap_info = swap_info.store,
              log.likelihood_PT = log.likelihood_PT.store, H_l = H_l.store,
              rho_adap_PT = rho_adap_PT.store, swap_accept_prob = swap_accept_prob.store,
              runtime=run.time, 
              acpt=c(a=nacpt, store=nacpt.bp, swap=nswap)/(burn+nsamp*thin))
  )
}
