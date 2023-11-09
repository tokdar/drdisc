###################### Load/Define functions setup ####################
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

order <- 2
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

my_summary <- function(x){
  c(quantile(x, 0.025), mean(x), quantile(x, 0.975))
}

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
    function(x) return(sqrt(9/128)*(35*x^4 - 30*x^2 + 3)) ## k = 4
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

###################### Parameter setup ####################

# y.grid <- seq(-1, 1, .01)
# p0 <- c(.5, .4, .1)
# shape1 <- c(10, 2, 1)
# shape2 <- c(10, 1, 1)
# get.f0 <- function(y) return((p0[1]*dbeta((y+1)/2,shape1[1],shape2[1])
#                               + p0[2]*dbeta((y+1)/2,shape1[2],shape2[2])
#                               + p0[3]*dbeta((y+1)/2,shape1[3],shape2[2]))/2)
# f0 <- get.f0(y.grid)
# h2.loss <- function(b){
#   f.nc <- integrate(fFn, -1, 1, beta=b)$value
#   return(integrate(function(y) (sqrt(fFn(y, beta=b)/f.nc) - sqrt(get.f0(y)))^2,-1,1)$value)
# }
# best.fit <- optim(rep(0,order), h2.loss, method="BFGS")
# b0 <- best.fit$par
# f0b <- get.f(b0)


###################### replicate function ####################

N_rep = 100
N_core = 8

rep_sim <- function(seed){
  simulate <- TRUE
  
  ## Simulate Data
  if(simulate){
    set.seed(seed)
    p <- 2
    n.obs <- 2e4
    b0x <- cbind(b0, threshold(rt(order, df=6),1))
    b0x[1,1] <-  1.2598
    b0x[1,2] <-  0
    b0x[2,1] <-  -1.6498
    b0x[2,2] <-  -1.3296
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
  trim05_ind <- (abs(y.obs) < 0.5)
  sim_data_trim05_y <- y.obs[trim05_ind]
  sim_data_trim05_x <- x.obs[trim05_ind,]
  names = paste0("./sim_20k_case_b/sim_20k_data/sim_notrim_order2_20k_data_seed_", seed, ".RData")
  save(sim_data_trim05_y, sim_data_trim05_x, x.obs, y.obs, file = names)
  
  ## hyperparameter/start point setup
  set.seed(seed)
  order <- 2
  p = dim(x.obs)[2]
  b_init <- matrix(rnorm(n = order*p, sd=1), order, p)
  a_init <- rnorm(n = p, sd = 1.5)
  
  ## fit BDRD model
  sim_notrim_order2_60k <- bdregjump_adapt_poly_trim_alphanoPG(y=y.obs, x=x.obs, 
                                                               b=b_init, burn=30000, nsamp=20000,
                                                               thin=1, trim = 1, order = order, prec = 1,
                                                               jump=list(a=a_init, prec = 1, positive=T,
                                                                         persistence=0.5, update.jbw=TRUE))
  names = paste0("./sim_20k_case_b/sim_20k_res/sim_notrim_order2_20k_seed_", seed, ".RData")
  save(sim_notrim_order2_60k, file = names)
  
  ## calculate WAIC and save result
  WAIC_res <- WAIC(sim_data_trim05_y, sim_data_trim05_x, sim_notrim_order2_60k, positive=T, usen = 1000, trim = 0.5)
  res <- data.frame(a1 = my_summary(sim_notrim_order2_60k$a[1,]),
                    a2 = my_summary(sim_notrim_order2_60k$a[2,]), 
                    b11 = my_summary(sim_notrim_order2_60k$b[1,1,]),
                    b12 = my_summary(sim_notrim_order2_60k$b[1,2,]), 
                    b21 = my_summary(sim_notrim_order2_60k$b[2,1,]), 
                    b22 = my_summary(sim_notrim_order2_60k$b[2,2,]), 
                    jbw = my_summary(sim_notrim_order2_60k$jbw), 
                    shapes1 = my_summary(sim_notrim_order2_60k$shapes[1, ]), 
                    shapes2 = my_summary(sim_notrim_order2_60k$shapes[2, ]))
  res <- cbind(res, WAIC1 = WAIC_res[1])
  res <- cbind(res, WAIC2 = WAIC_res[2])
  res <- as.matrix(res)
  names = paste0("./sim_20k_case_b/sim_20k_res/sim_notrim_order2_20k_model_seed_", seed, ".RData")
  save(res, file = names)
  return(res)
}


###################### replicate 100 times ####################
registerDoMC(N_core)
sim_order2_20k <- foreach(i = 1:N_rep) %dorng% rep_sim(i)
sim_order2_20k_notrim_order2 = array(NA, dim=c(3, 11, N_rep))
for(i in 1:N_rep){
  sim_order2_20k_notrim_order2[,,i] = sim_order2_20k[[i]]
}
sim_order2_20k_notrim_order2_summary = rowMeans(sim_order2_20k_notrim_order2, dims = 2)
colnames(sim_order2_20k_notrim_order2_summary) = colnames(sim_order2_20k[[i]])
sim_order2_20k_notrim_order2_summary
save.image("./sim_20k_case_b/sim_20k_res//sim_notrim_order2_20k_overall.RData")










