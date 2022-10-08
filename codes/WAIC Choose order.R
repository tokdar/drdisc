######################### WAIC ##########################################
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
setwd("~/Desktop/Density Estimation/Github code 9-14")
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

set.seed(123)
WAIC <- function(y.obs, x.obs, chain, positive=T, usen, trim){
  
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
  WAIC1 <- -2 * sum(rowMeans(log.lik.m))
  lik.m = exp(log.lik.m)
  WAIC2 <- -2 * sum(log(rowMeans(lik.m))) + 2 * sum(apply(log.lik.m, 1, var))
  return(c(WAIC1, WAIC2))
}


#################### no trimming order = 2 = true order ####################

load("~/Desktop/Density Estimation/Github code 9-14/report/report24/sim_notrim_order2.RData")
n_chain = length(sim_notrim_order2)
no_cores <- detectCores() 
registerDoMC(no_cores)
res <- foreach(i = 1:n_chain) %dorng% WAIC(y.obs, x.obs, sim_notrim_order2[[i]], positive=T, usen = 1000, trim = 1)
colMeans(matrix(unlist(res), ncol = 2, byrow = T))
# 5142.700 5151.393

#################### no trimming order = 3 > true order=2 ####################

load("~/Desktop/Density Estimation/Github code 9-14/report/report24/sim_notrim_order3.RData")
n_chain = length(sim_notrim_order3)
no_cores <- detectCores()
registerDoMC(no_cores)
res1 <- foreach(i = 1:n_chain) %dorng% WAIC(y.obs, x.obs, sim_notrim_order3[[i]], positive=T, usen = 1000, trim = 1)
colMeans(matrix(unlist(res1), ncol = 2, byrow = T))
# 5144.021 5154.260

#################### trimming 0.5 order = 2 ####################

trim05_ind <- (abs(y.obs) < 0.5)
sim_data_trim05_y <- y.obs[trim05_ind]
sim_data_trim05_x <- x.obs[trim05_ind,]


load("~/Desktop/Density Estimation/Github code 9-14/report/report24/sim_trim05_order2_new.RData")
n_chain = length(sim_trim05_order2_new)
no_cores <- detectCores()
registerDoMC(no_cores)
res2 <- foreach(i = 1:n_chain) %dorng% WAIC(sim_data_trim05_y, sim_data_trim05_x, sim_trim05_order2_new[[i]], positive=T, usen = 1000, trim = 0.5)
colMeans(matrix(unlist(res2), ncol = 2, byrow = T))
# -19899.56 -19892.22

#################### trimming 0.5 order = 2 ####################

load("~/Desktop/Density Estimation/Github code 9-14/report/report24/sim_trim05_order3_new.RData")
n_chain = length(sim_trim05_order3_new)
no_cores <- detectCores()
registerDoMC(no_cores)
res3 <- foreach(i = 1:n_chain) %dorng% WAIC(sim_data_trim05_y, sim_data_trim05_x, sim_trim05_order3_new[[i]], positive=T, usen = 1000, trim = 0.5)
colMeans(matrix(unlist(res3), ncol = 2, byrow = T))
# -19897.73 -19888.63

######################### WAIC ##########################################
##          If we want to gfit no-trimming model and  evaluate the WAIC 
##          of the same trimmed data, try following code
#########################################################################

registerDoMC(no_cores)
res4 <- foreach(i = 1:n_chain) %dorng% WAIC(sim_data_trim05_y, sim_data_trim05_x, sim_notrim_order2[[i]], positive=T, usen = 1000, trim = 1)
colMeans(matrix(unlist(res4), ncol = 2, byrow = T))


registerDoMC(no_cores)
res5 <- foreach(i = 1:n_chain) %dorng% WAIC(sim_data_trim05_y, sim_data_trim05_x, sim_notrim_order3[[i]], positive=T, usen = 1000, trim = 1)
colMeans(matrix(unlist(res5), ncol = 2, byrow = T))
