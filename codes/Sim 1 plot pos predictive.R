rm(list = ls())
setwd("~/Desktop/Density Estimation/Github code 3-16/drdisc-main/")

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

source("./codes/new-codes.R")
source("./codes/test-codes.R")
source("./codes/other-codes.R")

seed = 2
###################### load simulated data ###################### 
load(paste0("~/Desktop/Density Estimation/Github code 3-16/drdisc-main/sim_20k_case_b/sim_20k_data/sim_notrim_order2_20k_data_seed_", seed, ".RData"))
###################### load estimated parameters ###################### 
load(paste0("~/Desktop/Density Estimation/Github code 3-16/drdisc-main/sim_20k_case_b/sim_20k_res/sim_notrim_order2_20k_seed_", seed, ".RData"))

###################### simulation truth setup ###################### 
b0x = matrix(0, 2, 2)
b0x[1,1] <-  1.2598
b0x[1,2] <-  0
b0x[2,1] <-  -1.6498
b0x[2,2] <-  -1.3296
a0 = c(1, 4)
pos.synth = T
p = ncol(b0x)
true.order = 2
shapes0 <- c(4,1)
pers0 <- 1/2
jbw0 <- 0.16*2*pers0
pos.synth <- TRUE
pos.est <- TRUE

###################### hyperparameter setup ######################
order <- 2 # fitted order
trim = 1
xnew <- rep(1,4)
for(i in 1:(p-1)) xnew <- cbind(xnew, seq(-2,2,1.25))
nnew <- nrow(xnew)
xb0 <- tcrossprod(xnew, b0x)
xa0 <- c(xnew %*% a0)
if(pos.synth) xa0 <- pmax(0, xa0)
lpoly <- lpoly_all[1:true.order]
get.poly.mat <- function(y) return(sapply(lpoly, function(f) f(y)))
yFn <- get.poly.mat
wFn <- function(y, beta, alpha=0, bw=0.16) return(c(matrix(yFn(y), nrow=length(y)) %*% beta) + alpha*half.kern(y,bw))
fFn <- function(y, shapes=c(1,1), trim=1, ...) return((1/2)*dtrunc((y+1)/2, "beta", a = (1-trim)/2, b = (1+trim)/2, shape1 = shapes[1], shape2 = shapes[2])*Phi(wFn(y, ...)))
get.f <- function(b, a=0, sh=c(1,1), jbw=0.16, trim=1) {
  y.grid <- seq(-trim, trim, .01)
  f.grid <- fFn(y.grid, trim=trim, beta=b, alpha=a, shapes=sh, bw=jbw)
  normalization.const <- (integrate(fFn, -trim, 0, beta=b, alpha=a, shapes=sh, bw=jbw, trim = trim)$value
                          + integrate(fFn, 0, trim, beta=b, alpha=a, shapes=sh, bw=jbw, trim = trim)$value)
  
  return(f.grid/normalization.const)
}
f0x <- sapply(1:nnew, function(i) get.f(b=xb0[i,],a=xa0[i], sh=shapes0, jbw=jbw0,trim=trim))
a.est = sim_notrim_order2_60k$a
nsamp <- ncol(a.est)
b.est = sim_notrim_order2_60k$b
shapes.est = sim_notrim_order2_60k$shapes
jbw.est = sim_notrim_order2_60k$jbw

###################### plot posterior predictive ###################### 
par(mfrow = c(2,2))
for(cl in 1:3){
  ix.cl <- abs(x.obs[,2] - xnew[cl,2]) < 0.25
  hist(y.obs[ix.cl], freq=FALSE, main="", xlab="Y",
       col = 'white', cex.axis = 1.5, cex.lab = 1.5, xlim = c(-trim, trim))
  ylim <- par("usr")[3:4]
  xb <- apply(b.est, 3, function(b) tcrossprod(xnew[cl,,drop=FALSE], b))
  xa <- c(xnew[cl,,drop=FALSE] %*% a.est)
  if(pos.est) xa <- pmax(xa, 0)
  lpoly_all <- list(
    #  Vectorize(function(x) return(1/sqrt(2)), "x"), ## k = 0
    function(x) return(sqrt(3/2)*x), ## k = 1
    function(x) return(sqrt(5/8)*(3*x^2 - 1)), ## k = 2
    function(x) return(sqrt(7/8)*(5*x^3 - 3*x)), ## k = 3
    function(x) return(sqrt(9/128)*(35*x^4 - 30*x^2 + 3)) ## k = 4
  )
  
  lpoly <- lpoly_all[1:(dim(b.est)[1])]
  
  fx <- sapply(1:nsamp, function(s) get.f(b=xb[,s],a=xa[s],sh=shapes.est[,s],jbw=jbw.est[s],
                                          trim=trim))
  for(s in 1:nsamp) lines(y.grid, fx[,s], col=tcol(2,.5))
  lines(y.grid, rowMeans(fx), col=4, lwd=2)
  lines(y.grid, f0x[,cl], lwd=2)
  drop_pos <- trim/0.01
  dropx <- 1-fx[drop_pos,]/fx[drop_pos+1,]
  drop.est <- round(100*median(dropx))
  drop.ci <- round(100*quantile(dropx, pr=c(.025,.975)))
  text(-trim, ylim[2]*0.9, bquote('Drop%'==.(drop.est)[paste("[",.(drop.ci[1]),",",.(drop.ci[2]),"]")]), pos=4, cex = 1.5)
}


cl = 4
ix.cl <- abs(x.obs[,2] - xnew[cl,2]) < 0.25
hist(y.obs[ix.cl], freq=FALSE, main="", xlab="Y",
     col = 'white', cex.axis = 1.5, cex.lab = 1.5, xlim = c(-trim, trim), ylim  = c(0, 1.6))
ylim <- par("usr")[3:4]
xb <- apply(b.est, 3, function(b) tcrossprod(xnew[cl,,drop=FALSE], b))
xa <- c(xnew[cl,,drop=FALSE] %*% a.est)
if(pos.est) xa <- pmax(xa, 0)
lpoly_all <- list(
  #  Vectorize(function(x) return(1/sqrt(2)), "x"), ## k = 0
  function(x) return(sqrt(3/2)*x), ## k = 1
  function(x) return(sqrt(5/8)*(3*x^2 - 1)), ## k = 2
  function(x) return(sqrt(7/8)*(5*x^3 - 3*x)), ## k = 3
  function(x) return(sqrt(9/128)*(35*x^4 - 30*x^2 + 3)) ## k = 4
)
lpoly <- lpoly_all[1:(dim(b.est)[1])]
fx <- sapply(1:nsamp, function(s) get.f(b=xb[,s],a=xa[s],sh=shapes.est[,s],jbw=jbw.est[s],
                                        trim=trim))
for(s in 1:nsamp) lines(y.grid, fx[,s], col=tcol(2,.5))
lines(y.grid, rowMeans(fx), col=4, lwd=2)
lines(y.grid, f0x[,cl], lwd=2)
drop_pos <- trim/0.01
dropx <- 1-fx[drop_pos,]/fx[drop_pos+1,]
drop.est <- round(100*median(dropx))
drop.ci <- round(100*quantile(dropx, pr=c(.025,.975)))
text(-trim, ylim[2]*0.9, bquote('Drop%'==.(drop.est)[paste("[",.(drop.ci[1]),",",.(drop.ci[2]),"]")]), pos=4, cex = 1.5)
