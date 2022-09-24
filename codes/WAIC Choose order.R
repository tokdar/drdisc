######################## WAIC ##########################################
##          simulate data with order 2 (no trimming data/trimming data)
##          run model with order 2 and 3
##          use two type of WAICs to compare order
##          Now no trimming data makes sense, order 2 is a little better
##          than order 2. For trimming dta, there is a small problem maybe.
#########################################################################

library(abind)
library(knitr)
library(latex2exp)
require(magrittr)
require(plyr)
library(tidyverse)
library(coda)
library(doParallel)
library(doMC) 
library(foreach)
library(doRNG)
library(gridExtra)

adapt <- FALSE
bunch <- FALSE
N_para <- 2

registerDoMC(2)

source("./codes/new-codes.R")
source("./codes/test-codes.R")
source("./codes/other-codes.R")

simulate <- TRUE
## Density regression with jump
if(simulate){
  set.seed(123)
  p <- 2
  n.obs <- 6e4
  b0x <- cbind(b0, threshold(rt(order, df=6),1))
  while(ncol(b0x) < p) b0x <- cbind(b0x, 0)
  a0 <- c(1, 4, rep(0,p-2))
  shapes0 <- c(4,1)
  pers0 <- 1/2
  jbw0 <- 0.16*2*pers0
  pos.synth <- TRUE
  pos.est <- TRUE
  
  x.obs <- cbind(1, matrix(rnorm(n.obs*(p-1)), n.obs, p-1))
  y.obs <- simulate.fx(n.obs, b0x, x.obs, a0, shapes0, jbw0, positive=pos.synth)
}

N_para <- 2
set.seed(123)
b_init <- replicate(N_para, matrix(rnorm(n = order*p, sd = 3), order, p))
a_init <- replicate(N_para, rnorm(n = p, sd = 3))

set.seed(123)
order <- 2
b_init <- replicate(N_para, matrix(rnorm(n = order*p, sd=1), order, p))
a_init <- replicate(N_para, rnorm(n = p, sd = 1.5))

# sim_notrim_order2 <- foreach(i = 1:N_para) %dorng%
#   bdregjump_adapt_poly_trim_alphanoPG(y=y.obs, x=x.obs,
#                                       b=b_init[,,i], burn=15000, nsamp=10000,
#                                       thin=1, trim = 1, order = 2,
#                                       jump=list(a=a_init[,i], prec = 1, positive=pos.est,
#                                                 persistence=0.5, update.jbw=TRUE))
# save(sim_notrim_order2, file = "./report/report24/sim_notrim_order2.RData")
#load("./report/report24/sim_notrim_order2.RData")

set.seed(123)
order <- 3
b_init <- replicate(N_para, matrix(rnorm(n = order*p, sd=1), order, p))
a_init <- replicate(N_para, rnorm(n = p, sd = 1.5))

# sim_notrim_order3 <- foreach(i = 1:N_para) %dorng%
#   bdregjump_adapt_poly_trim_alphanoPG(y=y.obs, x=x.obs,
#                                       b=b_init[,,i], burn=15000, nsamp=10000,
#                                       thin=1, trim = 1, order = 3,
#                                       jump=list(a=a_init[,i], prec = 1, positive=pos.est,
#                                                 persistence=0.5, update.jbw=TRUE))
# save(sim_notrim_order3, file = "./report/report24/sim_notrim_order3.RData")
#load("./report/report24/sim_notrim_order3.RData")



WAIC1 <- function(y.obs, x.obs, chain, positive=T, usen, trim){
  
  order = dim(chain$b)[1]
  nsamp = length(chain$jbw)
  aa = chain$a
  bb = chain$b
  shsh = chain$shapes
  jbwjbw = chain$jbw
  
  lpoly <- list(
    #  Vectorize(function(x) return(1/sqrt(2)), "x"), ## k = 0
    function(x) return(sqrt(3/2)*x), ## k = 1
    function(x) return(sqrt(5/8)*(3*x^2 - 1)), ## k = 2
    function(x) return(sqrt(7/8)*(5*x^3 - 3*x)) ## k = 3
  )
  
  lpoly <- lpoly[1:order]
  ypoly <- function(y) return(sapply(lpoly, function(f) f(y)))
  fn.types <- list(poly=ypoly, kern=ykern, ns=yns)
  yFn <- fn.types[[type]]
  
  get.l.norm <- function(tt, sh, jbw, trim=1) {
    wFn.norm <- function(y, bw) tcrossprod(xb[tt,], yFn(y)) + xa[tt] * half.kern(y, jbw)
    Phi <- function(x) plogis(x)
    fFn.norm <- function(y, shapes, ...) return(dbeta((y+1)/2, shapes[1], shapes[2])*Phi(wFn.norm(y, ...)))
    norm.const <- integrate(fFn.norm, -trim, 0,  shapes=sh, bw=jbw)$value + 
      integrate(fFn.norm, 0, trim, shapes=sh, bw=jbw)$value
    return(log(norm.const))
  }
  n.obs = length(y.obs)
  pm = yFn(y.obs)
  y_b = (y.obs+1)/2
  log.lik.m = matrix(NA, nrow = n.obs, ncol = usen)
  for(i in 1:usen){
    a = aa[, nsamp-i+1]
    b = bb[,,nsamp-i+1]
    sh = shsh[, nsamp-i+1]
    jbw = jbwjbw[nsamp-i+1]
    xa <- c(x.obs %*% a)
    if(positive) xa <- pmax(0, xa)
    xb <- tcrossprod(x.obs, b)
    w <- rowSums(xb * pm) + xa * half.kern(y.obs, jbw)
    log.lik.m[,i] = dbeta(y_b, sh[1], sh[2], log = T) + plogis(w, log = T) -
      sapply(1:n.obs, function(t) get.l.norm(t, sh, jbw, trim=trim))
  }
  WAIC <- -2 * sum(rowMeans(log.lik.m))
  return(WAIC)
}

WAIC2 <- function(y.obs, x.obs, chain, positive=T, usen, trim){
  
  order = dim(chain$b)[1]
  nsamp = length(chain$jbw)
  aa = chain$a
  bb = chain$b
  shsh = chain$shapes
  jbwjbw = chain$jbw
  
  lpoly <- list(
    #  Vectorize(function(x) return(1/sqrt(2)), "x"), ## k = 0
    function(x) return(sqrt(3/2)*x), ## k = 1
    function(x) return(sqrt(5/8)*(3*x^2 - 1)), ## k = 2
    function(x) return(sqrt(7/8)*(5*x^3 - 3*x)) ## k = 3
  )
  
  lpoly <- lpoly[1:order]
  ypoly <- function(y) return(sapply(lpoly, function(f) f(y)))
  fn.types <- list(poly=ypoly, kern=ykern, ns=yns)
  yFn <- fn.types[[type]]
  
  get.l.norm <- function(tt, sh, jbw, trim=1) {
    wFn.norm <- function(y, bw) tcrossprod(xb[tt,], yFn(y)) + xa[tt] * half.kern(y, jbw)
    Phi <- function(x) plogis(x)
    fFn.norm <- function(y, shapes, ...) return(dbeta((y+1)/2, shapes[1], shapes[2])*Phi(wFn.norm(y, ...)))
    norm.const <- integrate(fFn.norm, -trim, 0,  shapes=sh, bw=jbw)$value +
      integrate(fFn.norm, 0, trim, shapes=sh, bw=jbw)$value
    return(log(norm.const))
  }
  n.obs = length(y.obs)
  pm = yFn(y.obs)
  y_b = (y.obs+1)/2
  log.lik.m = matrix(NA, nrow = n.obs, ncol = usen)
  for(i in 1:usen){
    a = aa[, nsamp-i+1]
    b = bb[,,nsamp-i+1]
    sh = shsh[, nsamp-i+1]
    jbw = jbwjbw[nsamp-i+1]
    xa <- c(x.obs %*% a)
    if(positive) xa <- pmax(0, xa)
    xb <- tcrossprod(x.obs, b)
    w <- rowSums(xb * pm) + xa * half.kern(y.obs, jbw)
    log.lik.m[,i] = dbeta(y_b, sh[1], sh[2], log = T) + plogis(w, log = T) -
      sapply(1:n.obs, function(t) get.l.norm(t, sh, jbw, trim=trim))
  }
  lik.m = exp(log.lik.m)
  WAIC <- -2 * sum(log(rowMeans(lik.m))) + 2 * sum(apply(log.lik.m, 1, var))
  return(WAIC)
}

#################### no trimming ####################

## order 2 = true order

time.stamp.0 <- proc.time()
WAIC1(y.obs, x.obs, sim_notrim_order2[[1]], positive=T, 1000, 1)
run.time <- proc.time() - time.stamp.0
run.time
# 5141.406  3.09h

time.stamp.0 <- proc.time()
WAIC2(y.obs, x.obs, sim_notrim_order2[[1]], positive=T, 1000, 1)
run.time <- proc.time() - time.stamp.0
run.time
# 5149.791  3.08h
                         
time.stamp.0 <- proc.time()
WAIC1(y.obs, x.obs, sim_notrim_order2[[2]], positive=T, 1000, 1)
run.time <- proc.time() - time.stamp.0
run.time
# 5141.045  3.27h

time.stamp.0 <- proc.time()
WAIC2(y.obs, x.obs, sim_notrim_order2[[2]], positive=T, 1000, 1)
run.time <- proc.time() - time.stamp.0
run.time
# 5149.08  3.30h
             


## order 3

time.stamp.0 <- proc.time()
WAIC1(y.obs, x.obs, sim_notrim_order3[[1]], positive=T, 1000, 1)
run.time <- proc.time() - time.stamp.0
run.time
# 5144.017  3.57h
             
time.stamp.0 <- proc.time()
WAIC2(y.obs, x.obs, sim_notrim_order3[[1]], positive=T, 1000, 1)
run.time <- proc.time() - time.stamp.0
run.time
# 5155.026 3.60h

time.stamp.0 <- proc.time()
WAIC1(y.obs, x.obs, sim_notrim_order3[[2]], positive=T, 1000, 1)
run.time <- proc.time() - time.stamp.0
run.time
# 5143.613 3.57h

time.stamp.0 <- proc.time()
WAIC2(y.obs, x.obs, sim_notrim_order3[[2]], positive=T, 1000, 1)
run.time <- proc.time() - time.stamp.0
run.time
# 5154.565 3.60h

#################### trimming 0.5 ####################

## order 2 = true order

trim05_ind <- (abs(y.obs) < 0.5)
sim_data_trim05_y <- y.obs[trim05_ind]
sim_data_trim05_x <- x.obs[trim05_ind,]

load("~/Desktop/Density Estimation/Github code 9-14/report/report24/sim_trim05_order2_new.RData")
time.stamp.0 <- proc.time()
WAIC1(sim_data_trim05_y, sim_data_trim05_x, sim_trim05_order2_new[[1]], positive=T, 1000, .5)
run.time <- proc.time() - time.stamp.0
run.time
# -19898.84 0.39h

time.stamp.0 <- proc.time()
WAIC2(sim_data_trim05_y, sim_data_trim05_x, sim_trim05_order2_new[[1]], positive=T, 1000, .5)
run.time <- proc.time() - time.stamp.0
run.time
# -19892.03 0.40h

time.stamp.0 <- proc.time()
WAIC1(sim_data_trim05_y, sim_data_trim05_x, sim_trim05_order2_new[[2]], positive=T, 1000, .5)
run.time <- proc.time() - time.stamp.0
run.time
# -19901.76 0.40h


time.stamp.0 <- proc.time()
WAIC2(sim_data_trim05_y, sim_data_trim05_x, sim_trim05_order2_new[[2]], positive=T, 1000, .5)
run.time <- proc.time() - time.stamp.0
run.time
# -19893.61  0.39h


## order 3

load("~/Desktop/Density Estimation/Github code 9-14/report/report24/sim_trim05_order3_new1.RData")

time.stamp.0 <- proc.time()
WAIC1(sim_data_trim05_y, sim_data_trim05_x, sim_trim05_order3_new1, positive=T, 1000, .5)
run.time <- proc.time() - time.stamp.0
run.time
# -19902.01  0.42h

time.stamp.0 <- proc.time()
WAIC2(sim_data_trim05_y, sim_data_trim05_x, sim_trim05_order3_new1, positive=T, 1000, .5)
run.time <- proc.time() - time.stamp.0
run.time
# -19893.44 0.42h

load("~/Desktop/Density Estimation/Github code 9-14/report/report24/sim_trim05_order3_new2.RData")

time.stamp.0 <- proc.time()
WAIC1(sim_data_trim05_y, sim_data_trim05_x, sim_trim05_order3_new2, positive=T, 1000, .5)
run.time <- proc.time() - time.stamp.0
run.time
# -19899.54  0.42h

time.stamp.0 <- proc.time()
WAIC2(sim_data_trim05_y, sim_data_trim05_x, sim_trim05_order3_new2, positive=T, 1000, .5)
run.time <- proc.time() - time.stamp.0
run.time
# -19891.09 0.40h

