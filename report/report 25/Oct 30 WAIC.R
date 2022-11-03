rm(list = ls())

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
N_para <- 3

registerDoMC(2)
setwd("~/Desktop/Density Estimation/Github code 10-10")
source("./codes/new-codes.R")
source("./codes/test-codes.R")
source("./codes/other-codes.R")

load("~/Desktop/Density Estimation/Github code 10-10/Oct 30/Oct 30 code setup.RData")
trim05_ind <- (abs(y.obs) < 0.5)
sim_data_trim05_y <- y.obs[trim05_ind]  # 16513
sim_data_trim05_x <- x.obs[trim05_ind,]
length(sim_data_trim05_y) # 16407
hist(sim_data_trim05_y, main = 'Trimmed Data', xlim = c(-0.5, 0.5))

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
    function(x) return(sqrt(7/8)*(5*x^3 - 3*x)), ## k = 3
    function(x) return(sqrt(9/8)*(35*x^4 - 30*x^2 + 3)) ## k = 4
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


#################### no trimming order 2 ####################

load("~/Desktop/Density Estimation/Github code 10-10/Oct 30/Oct 30 sim_notrim_order2_20k.RData")
n_chain = length(sim_notrim_order2_20k)
registerDoMC(n_chain)
res <- foreach(i = 1:n_chain) %dorng% WAIC(y.obs, x.obs, sim_notrim_order2_20k[[i]], 
                                           positive=T, usen = 1000, trim = 1)
colMeans(matrix(unlist(res), ncol = 2, byrow = T))
#  -3219.456 -3214.052

#################### no trimming order 3 ####################

load("~/Desktop/Density Estimation/Github code 10-10/Oct 30/Oct 30 sim_notrim_order3_20k.RData")
n_chain = length(sim_notrim_order3_20k)
registerDoMC(n_chain)
not3 <- foreach(i = 1:n_chain) %dorng% WAIC(y.obs, x.obs, sim_notrim_order3_20k[[i]], 
                                            positive=T, usen = 1000, trim = 1)
colMeans(matrix(unlist(not3), ncol = 2, byrow = T))
#  -3784.375 -3777.419

#################### no trimming order 4 ####################

load("~/Desktop/Density Estimation/Github code 10-10/Oct 30/Oct 30 sim_notrim_order4_20k.RData")
n_chain = length(sim_notrim_order4_20k)
registerDoMC(n_chain)
not4 <- foreach(i = 1:n_chain) %dorng% WAIC(y.obs, x.obs, sim_notrim_order4_20k[[i]], 
                                            positive=T, usen = 1000, trim = 1)
colMeans(matrix(unlist(not4), ncol = 2, byrow = T))
#  -5002.221 -4994.198

#################### trimming order 2 ####################
load("~/Desktop/Density Estimation/Github code 10-10/Oct 30/Oct 30 sim_trim05_order2_20k.RData")
n_chain = length(sim_trim05_order2_20k)
registerDoMC(n_chain)
t2 <- foreach(i = 1:n_chain) %dorng% WAIC(sim_data_trim05_y, sim_data_trim05_x, sim_trim05_order2_20k[[i]], 
                                           positive=T, usen = 1000, trim = 0.5)
colMeans(matrix(unlist(t2), ncol = 2, byrow = T))
# -9947.736 -9943.592 

not2ont <- foreach(i = 1:n_chain) %dorng% WAIC(sim_data_trim05_y, sim_data_trim05_x, sim_notrim_order2_20k[[i]], 
                                          positive=T, usen = 1000, trim = 0.5)
colMeans(matrix(unlist(not2ont), ncol = 2, byrow = T))
# -8880.975 -8877.686

#################### trimming order 3 ####################
load("~/Desktop/Density Estimation/Github code 10-10/Oct 30/Oct 30 sim_trim05_order3_20k.RData")
n_chain = length(sim_trim05_order3_20k)
registerDoMC(n_chain)
t3 <- foreach(i = 1:n_chain) %dorng% WAIC(sim_data_trim05_y, sim_data_trim05_x, sim_trim05_order3_20k[[i]], 
                                          positive=T, usen = 1000, trim = 0.5)
colMeans(matrix(unlist(t3), ncol = 2, byrow = T))
# -9946.819 -9942.160

not3ont <- foreach(i = 1:n_chain) %dorng% WAIC(sim_data_trim05_y, sim_data_trim05_x, sim_notrim_order3_20k[[i]], 
                                               positive=T, usen = 1000, trim = 0.5)
colMeans(matrix(unlist(not3ont), ncol = 2, byrow = T))
# -9248.291 -9243.652

#################### trimming order 4 ####################
load("~/Desktop/Density Estimation/Github code 10-10/Oct 30/Oct 30 sim_trim05_order4_20k.RData")
n_chain = length(sim_trim05_order4_20k)
registerDoMC(n_chain)
t4 <- foreach(i = 1:n_chain) %dorng% WAIC(sim_data_trim05_y, sim_data_trim05_x, sim_trim05_order4_20k[[i]], 
                                          positive=T, usen = 1000, trim = 0.5)
colMeans(matrix(unlist(t4), ncol = 2, byrow = T))
#  -9975.652 -9970.131

not4ont <- foreach(i = 1:n_chain) %dorng% WAIC(sim_data_trim05_y, sim_data_trim05_x, sim_notrim_order4_20k[[i]], 
                                               positive=T, usen = 1000, trim = 0.5)
colMeans(matrix(unlist(not4ont), ncol = 2, byrow = T))
# -9642.399 -9636.573




