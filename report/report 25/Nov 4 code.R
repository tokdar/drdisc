rm(list = ls())
setwd("~/Desktop/Density Estimation/Github code 11-4")

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

order <- 4
lpoly_all <- list(
  #  Vectorize(function(x) return(1/sqrt(2)), "x"), ## k = 0
  function(x) return(sqrt(3/2)*x), ## k = 1
  function(x) return(sqrt(5/8)*(3*x^2 - 1)), ## k = 2
  function(x) return(sqrt(7/8)*(5*x^3 - 3*x)), ## k = 3
  function(x) return(sqrt(9/128)*(35*x^4 - 30*x^2 + 3)) ## k = 4
)

lpoly <- lpoly_all[1:order]
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
wFn <- function(y, beta, alpha=0, bw=0.16) return(c(matrix(yFn(y), nrow=length(y)) %*% beta) + alpha*half.kern(y,bw))
Phi <- function(x) plogis(x)
fFn <- function(y, shapes=c(1,1), trim=1, ...) return((1/2)*dtrunc((y+1)/2, "beta", a = (1-trim)/2, b = (1+trim)/2,
                                                                   shape1 = shapes[1], shape2 = shapes[2])*Phi(wFn(y, ...)))

y.grid <- seq(-1,1,.01)
p0 <- c(.4, .5, .1)
shape1 <- c(10, .1, .1)
shape2 <- c(10, 100, 100)
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

seed = 1

simulate <- TRUE
## Density regression with jump
if(simulate){
  set.seed(seed)
  p <- 2
  n.obs <- 1e4
  b0x <- cbind(b0, threshold(rt(order, df=6),1))
  b0x[1,1] <-  -0.9706636
  b0x[1,2] <-  0
  b0x[2,1] <-  -15
  b0x[2,2] <-  -1.329144
  b0x[3,1] <-  -0.9773
  b0x[3,2] <-  0
  b0x[4,1] <-  0.5
  b0x[4,2] <-  -1.770323
  while(ncol(b0x) < p) b0x <- cbind(b0x, 0)
  a0 <- c(2, 4, rep(0,p-2))
  shapes0 <- c(4,1)
  pers0 <- 1/2
  jbw0 <- 0.16*2*pers0
  pos.synth <- TRUE
  pos.est <- TRUE

  x.obs <- cbind(1, matrix(rnorm(n.obs*(p-1)), n.obs, p-1))
  y.obs <- simulate.fx(n.obs, b0x, x.obs, a0, shapes0, jbw0, positive=pos.synth)
}

hist(y.obs, xlim = c(-1, 1))
b0x

trim05_ind <- (abs(y.obs) < 0.5)
sim_data_trim05_y <- y.obs[trim05_ind]  # 8318
sim_data_trim05_x <- x.obs[trim05_ind,]
hist(sim_data_trim05_y, main = 'Trimmed Data')
length(sim_data_trim05_y)

names = paste0("./", "Nov 4 setup seed", seed, ".RData")
save.image(names)
load(names)

###################### no trimming order 2,3,4 ####################
# N_para <- 10
# set.seed(seed)
# registerDoMC(N_para)
# order = 2
# p = dim(x.obs)[2]
# b_init <- replicate(N_para, matrix(rnorm(n = order*p, sd = 3), order, p))
# a_init <- replicate(N_para, rnorm(n = p, sd = 3))
# sim_notrim_order2_10k<- foreach(i = 1:N_para) %dorng% 
#   bdregjump_adapt_poly_trim_alphanoPG(y=y.obs, x=x.obs,
#                                       b=b_init[,,i], burn=20000, nsamp=10000,
#                                       thin=1, trim = 1, order = order, prec=1,
#                                       jump=list(a=a_init[,i], prec = 1, positive=T,
#                                                 persistence=0.5, update.jbw=TRUE))
# names = paste0("./", "sim_notrim_order2_10k seed ", seed, ".RData")
# save(sim_notrim_order2_10k, file = names)

set.seed(seed)
order <- 3
N_para <- 10
registerDoMC(N_para)
p = dim(x.obs)[2]
b_init <- replicate(N_para, matrix(rnorm(n = order*p, sd=1), order, p))
a_init <- replicate(N_para, rnorm(n = p, sd = 1.5))
sim_notrim_order3_10k <- foreach(i = 1:N_para) %dorng% bdregjump_adapt_poly_trim_alphanoPG(y=y.obs, x=x.obs,
                                                                                           b=b_init[,,i], burn=20000, nsamp=10000,
                                                                                           thin=1, trim = 1, order = order, prec = 1,
                                                                                           jump=list(a=a_init[,i], prec = 1, positive=T,
                                                                                                     persistence=0.5, update.jbw=TRUE))
names = paste0("./", "sim_notrim_order3_10k seed ", seed, ".RData")
save(sim_notrim_order3_10k, file = names)


set.seed(seed)
order = 4
N_para <- 10
registerDoMC(N_para)
p = dim(x.obs)[2]
b_init <- replicate(N_para, matrix(rnorm(n = order*p, sd = 3), order, p))
a_init <- replicate(N_para, rnorm(n = p, sd = 3))
sim_notrim_order4_10k<- foreach(i = 1:N_para) %dorng%
  bdregjump_adapt_poly_trim_alphanoPG(y=y.obs, x=x.obs,
                                      b=b_init[,,i], burn=20000, nsamp=10000,
                                      thin=1, trim = 1, order = order, prec = 1,
                                      jump=list(a=a_init[,i], prec = 1, positive=T,
                                                persistence=0.5, update.jbw=TRUE))
names = paste0("./", "sim_notrim_order4_10k seed ", seed, ".RData")
save(sim_notrim_order4_10k, file = names)

################################# trimming #####################################

# order = 2
# N_para <- 10
# registerDoMC(N_para)
# p = dim(x.obs)[2]
# b_init <- replicate(N_para, matrix(rnorm(n = order*p, sd=1), order, p))
# a_init <- replicate(N_para, rnorm(n = p, sd = 1.5))
# sim_trim05_order2_10k <- foreach(i = 1:N_para) %dorng% 
#   bdregjump_adapt_poly_trim_alphanoPG(y=sim_data_trim05_y, x=sim_data_trim05_x, 
#                                       b=b_init[,,i], prec = 1, burn=25000, nsamp=15000, 
#                                       thin=1, trim = 0.5, shapes = c(2,2), order = order,
#                                       jump=list(a=a_init[,i], prec = 1, positive=T,
#                                                 persistence=0.5,update.jbw=TRUE))
# names = paste0("./", "sim_trim05_order2_10k seed ", seed, ".RData")
# save(sim_trim05_order2_10k, file = names)


order = 3
N_para <- 10
registerDoMC(N_para)
p = dim(x.obs)[2]
b_init <- replicate(N_para, matrix(rnorm(n = order*p, sd=1), order, p))
a_init <- replicate(N_para, rnorm(n = p, sd = 1.5))
sim_trim05_order3_10k <- foreach(i = 1:N_para) %dorng% 
  bdregjump_adapt_poly_trim_alphanoPG(y=sim_data_trim05_y, x=sim_data_trim05_x, 
                                      b=b_init[,,i], prec = 0.1,burn=25000, nsamp=15000,
                                      thin=1, trim = 0.5, shapes = c(2,2), order = order,
                                      jump=list(a=a_init[,i], prec = 1, positive=T,
                                                persistence=0.5,update.jbw=TRUE))
names = paste0("./", "sim_trim05_order3_10k seed ", seed, ".RData")
save(sim_trim05_order3_10k, file = names)

order = 4
N_para <- 10
registerDoMC(N_para)
p = dim(x.obs)[2]
b_init <- replicate(N_para, matrix(rnorm(n = order*p, sd=1), order, p))
a_init <- replicate(N_para, rnorm(n = p, sd = 1.5))
sim_trim05_order4_10k <- foreach(i = 1:N_para) %dorng% 
  bdregjump_adapt_poly_trim_alphanoPG(y=sim_data_trim05_y, x=sim_data_trim05_x, 
                                      b=b_init[,,i], prec = 0.1,burn=25000, nsamp=15000,
                                      thin=1, trim = 0.5, shapes = c(2,2), order = order,
                                      jump=list(a=a_init[,i], prec = 1, positive=T,
                                                persistence=0.5,update.jbw=TRUE))
names = paste0("./", "sim_trim05_order4_10k seed ", seed, ".RData")
save(sim_trim05_order4_10k, file = names)
