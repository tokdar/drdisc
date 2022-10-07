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

order = 2
time.stamp.0 <- proc.time()
res = WAIC(y.obs, x.obs, sim_notrim_order2[[1]], positive=T, 1000, 1)
res
run.time <- proc.time() - time.stamp.0
run.time
# 5143.273 5152.375.  3.04h

time.stamp.0 <- proc.time()
res2 = WAIC(y.obs, x.obs, sim_notrim_order2[[2]], positive=T, 1000, 1)
res2
run.time <- proc.time() - time.stamp.0
run.time
# 5142.829 5151.546

time.stamp.0 <- proc.time()
res3 = WAIC(y.obs, x.obs, sim_notrim_order2[[3]], positive=T, 1000, 1)
res3
run.time <- proc.time() - time.stamp.0
run.time
# 5140.928 5148.544

time.stamp.0 <- proc.time()
res4 = WAIC(y.obs, x.obs, sim_notrim_order2[[4]], positive=T, 1000, 1)
res4
run.time <- proc.time() - time.stamp.0
run.time
# 5142.225 5150.903

time.stamp.0 <- proc.time()
res5 = WAIC(y.obs, x.obs, sim_notrim_order2[[5]], positive=T, 1000, 1)
res5
run.time <- proc.time() - time.stamp.0
run.time
# 5141.499 5149.395


time.stamp.0 <- proc.time()
res6 = WAIC(y.obs, x.obs, sim_notrim_order2[[6]], positive=T, 1000, 1)
res6
run.time <- proc.time() - time.stamp.0
run.time
# 5143.010 5151.806

#################### no trimming order = 3 > true order=2 ####################


load("~/Desktop/Density Estimation/Github code 9-14/report/report24/sim_notrim_order3.RData")

time.stamp.0 <- proc.time()
res = WAIC(y.obs, x.obs, sim_notrim_order3[[1]], positive=T, 1000, 1)
res
run.time <- proc.time() - time.stamp.0
run.time
# 5144.440 5155.418

time.stamp.0 <- proc.time()
res2 = WAIC(y.obs, x.obs, sim_notrim_order3[[2]], positive=T, 1000, 1)
res2
run.time <- proc.time() - time.stamp.0
run.time
# 5144.482 5154.484

time.stamp.0 <- proc.time()
res3 = WAIC(y.obs, x.obs, sim_notrim_order3[[3]], positive=T, 1000, 1)
res3
run.time <- proc.time() - time.stamp.0
run.time
# 5144.563 5155.208

time.stamp.0 <- proc.time()
res4 = WAIC(y.obs, x.obs, sim_notrim_order3[[4]], positive=T, 1000, 1)
res4
run.time <- proc.time() - time.stamp.0
run.time
# 5143.959 5154.276

time.stamp.0 <- proc.time()
res5 = WAIC(y.obs, x.obs, sim_notrim_order3[[5]], positive=T, 1000, 1)
res5
run.time <- proc.time() - time.stamp.0
run.time
# 5144.346 5154.914


time.stamp.0 <- proc.time()
res6 = WAIC(y.obs, x.obs, sim_notrim_order3[[6]], positive=T, 1000, 1)
res6
run.time <- proc.time() - time.stamp.0
run.time
# 5145.010 5155.733



#################### trimming 0.5 ####################

trim05_ind <- (abs(y.obs) < 0.5)
sim_data_trim05_y <- y.obs[trim05_ind]
sim_data_trim05_x <- x.obs[trim05_ind,]


# load("~/Desktop/Density Estimation/Github code 9-14/report/report24/sim_trim05_order2_new.RData")
# time.stamp.0 <- proc.time()
# WAIC1(sim_data_trim05_y, sim_data_trim05_x, sim_trim05_order2_new[[1]], positive=T, 1000, .5)
# run.time <- proc.time() - time.stamp.0
# run.time
# # -19898.84 0.39h
# 
# time.stamp.0 <- proc.time()
# WAIC2(sim_data_trim05_y, sim_data_trim05_x, sim_trim05_order2_new[[1]], positive=T, 1000, .5)
# run.time <- proc.time() - time.stamp.0
# run.time
# # -19892.03 0.40h
# 
# time.stamp.0 <- proc.time()
# WAIC1(sim_data_trim05_y, sim_data_trim05_x, sim_trim05_order2_new[[2]], positive=T, 1000, .5)
# run.time <- proc.time() - time.stamp.0
# run.time
# # -19901.76 0.40h
# 
# 
# time.stamp.0 <- proc.time()
# WAIC2(sim_data_trim05_y, sim_data_trim05_x, sim_trim05_order2_new[[2]], positive=T, 1000, .5)
# run.time <- proc.time() - time.stamp.0
# run.time
# # -19893.61  0.39h
# 
# 
# ## order 3
# 
# load("~/Desktop/Density Estimation/Github code 9-14/report/report24/sim_trim05_order3_new1.RData")
# 
# time.stamp.0 <- proc.time()
# WAIC1(sim_data_trim05_y, sim_data_trim05_x, sim_trim05_order3_new1, positive=T, 1000, .5)
# run.time <- proc.time() - time.stamp.0
# run.time
# # -19902.01  0.42h
# 
# time.stamp.0 <- proc.time()
# WAIC2(sim_data_trim05_y, sim_data_trim05_x, sim_trim05_order3_new1, positive=T, 1000, .5)
# run.time <- proc.time() - time.stamp.0
# run.time
# # -19893.44 0.42h
# 
# load("~/Desktop/Density Estimation/Github code 9-14/report/report24/sim_trim05_order3_new2.RData")
# 
# time.stamp.0 <- proc.time()
# WAIC1(sim_data_trim05_y, sim_data_trim05_x, sim_trim05_order3_new2, positive=T, 1000, .5)
# run.time <- proc.time() - time.stamp.0
# run.time
# # -19899.54  0.42h
# 
# time.stamp.0 <- proc.time()
# WAIC2(sim_data_trim05_y, sim_data_trim05_x, sim_trim05_order3_new2, positive=T, 1000, .5)
# run.time <- proc.time() - time.stamp.0
# run.time
# # -19891.09 0.40h
# 
