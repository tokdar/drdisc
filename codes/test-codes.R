# fix cov, change lambda

bdregjump_adapt_handtune <- function(y, x=1, b=NULL, nsamp=100, thin=1, burn=100, 
                                    shapes=c(1,1), prec=1, cov_tune, lambda_tune,
                                    jump=list(a=NULL, persistence=0.5, positive=TRUE, 
                                              prec=1, ncand=10, update.jbw=TRUE),
                                    print.process = FALSE, print.i = FALSE){
  
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

  Poly.obs <- get.poly.mat(y)
  b.store <- array(NA, dim=c(order, p, nsamp))
  a.store <- matrix(NA, p, nsamp)
  az.new.store <- matrix(NA, p, nsamp)
  jbw.store <- rep(NA, nsamp)
  Ntot.store <- rep(NA, nsamp)
  shapes.store <- matrix(NA, 2, nsamp)
  sample_i_count.store <- rep(NA, nsamp)
  p_accept.store <- rep(NA, nsamp)
  count.store <- 0
  nacpt <- 0
  nacpt.bp <- 0
  bp.df <- 6
  az.new <- rep(0,p)
  
  diff_prev <- burn_i_count <- sample_i_count <- 0
  
  Sigma_tune <- chol(lambda_tune*cov_tune)
  
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
    
    az.new <- crossprod(Sigma_tune, rnorm(p))
    a.new <- a + az.new
    xa.new <- c(x.comb %*% a.new)
    if(positive) xa.new <- pmax(0, xa.new)
    
    log.hastings.ratio <- 0
    if(positive){
      res.new <- pg.residual - xa.new * Half.comb
      res <- pg.residual - xa * Half.comb
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
    }
  }
  run.time <- proc.time() - time.stamp.0
  return(list(last_missing = missing.stuff,
              b=b.store, a=a.store, az.new = az.new.store, jbw=jbw.store, Ntot=Ntot.store,
              shapes=shapes.store, 
              runtime=run.time, 
              acpt=c(a=nacpt, store=nacpt.bp)/(burn+nsamp*thin))
  )
}

###### no a ######

bdregjump_noa <- function(y, x=1, b=NULL, nsamp=100, thin=1, burn=100, 
                            shapes=c(1,1), prec=1, order = 5,
                            print.process = FALSE, print.i = FALSE){
  
  n <- length(y)
  if(is.null(dim(x))) x <- matrix(x, nrow=n)
  p <- ncol(x)
  if(is.null(b)) b <- replicate(p, rep(0,order))
  b <- matrix(b, order, p)
  if(length(prec) < p) prec <- rep(prec, p)[1:p]
  shapes.bp <- c(qlogis(shapes[1]/sum(shapes)), sum(shapes))
  
  kbw <- 2/(order-1)
  knots <- seq(-1,1,kbw)
  bsd.sq <- (0.67*kbw)^2
  gausskern <- lapply(knots, function(mu) return(function(x) return(exp(-0.5*(x-mu)^2/bsd.sq))))
  get.poly.mat <- function(y) return(sapply(gausskern, function(f) f(y)))
  
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
  
  Poly.obs <- get.poly.mat(y)
  b.store <- array(NA, dim=c(order, p, nsamp))
  Ntot.store <- rep(NA, nsamp)
  shapes.store <- matrix(NA, 2, nsamp)
  count.store <- 0
  nacpt <- 0
  nacpt.bp <- 0
  bp.df <- 6
  
  a <- rep(0, p)
  
  time.stamp.0 <- proc.time()
  for(iter in -(burn-1):(nsamp*thin)){
    
    w.obs <- rowSums(x * (Poly.obs %*% b))
    
    missing.stuff <- get.miss.reject.x(b, n, x, a, shapes)
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
      Ntot.store[count.store] <- Ntot
      shapes.store[,count.store] <- shapes
    }
  }
  run.time <- proc.time() - time.stamp.0
  return(list(last_missing = missing.stuff,
              b=b.store, Ntot=Ntot.store, shapes=shapes.store, runtime=run.time, 
              acpt=c(a=nacpt, store=nacpt.bp)/(burn+nsamp*thin))
  )
}
