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


n_chain <- 10

################### load ###################

load("./real_data_result/Feb_19/realdata_notrim_order2_40k_seed_3.RData")
load("./real_data_result/Feb_19/realdata_trim05_order2_40k_seed_3.RData")

load("./real_data_result/Feb_20/realdata_trim09_order2_40k_seed_3.RData")
load("./real_data_result/Feb_20/realdata_trim08_order2_40k_seed_3.RData")
load("./real_data_result/Feb_20/realdata_trim07_order2_40k_seed_3.RData")
load("./real_data_result/Feb_20/realdata_trim06_order2_40k_seed_3.RData")

load("./real_data_result/Feb_21/realdata_notrim_order3_40k_seed_3.RData")
load("./real_data_result/Feb_21/realdata_trim09_order3_40k_seed_3.RData")
load("./real_data_result/Feb_21/realdata_trim08_order3_40k_seed_3.RData")
load("./real_data_result/Feb_21/realdata_trim07_order3_40k_seed_3.RData")
load("./real_data_result/Feb_21/realdata_trim06_order3_40k_seed_3.RData")
load("./real_data_result/Feb_21/realdata_trim05_order3_40k_seed_3.RData")


load("./real_data_result/Feb_25/realdata_notrim_order4_65k_seed_3.RData")
load("./real_data_result/Feb_23/realdata_trim09_order4_40k_seed_3.RData")
load("./real_data_result/Feb_23/realdata_trim08_order4_40k_seed_3.RData")
load("./real_data_result/Feb_23/realdata_trim07_order4_40k_seed_3.RData")
load("./real_data_result/Feb_25/realdata_trim06_order4_65k_seed_3.RData")
load("./real_data_result/Feb_25/realdata_trim05_order4_65k_seed_sel.RData")


################### data ###################

data_notrim <- get.trim.data(1)
data_trim09 <- get.trim.data(0.9)
data_trim08 <- get.trim.data(0.8)
data_trim07 <- get.trim.data(0.7)
data_trim06 <- get.trim.data(0.6)
data_trim05 <- get.trim.data(0.5)


#################### order = 2 ####################

# no trimming #

comp_notrim_order2 <- foreach(i = 1:n_chain) %dorng% WAIC(data_trim05$y, data_trim05$x, realdata_notrim_order2[[i]], positive=T, usen = 1000, trim = 1)
waic_notrim_order2 <- colMeans(matrix(unlist(comp_notrim_order2), ncol = 2, byrow = T))

save(comp_notrim_order2, file = "./real_data_result/Mar_4/comp_notrim_order2.RData")

# trim = 0.9 #

comp_trim09_order2 <- foreach(i = 1:n_chain) %dorng% WAIC(data_trim05$y, data_trim05$x, realdata_trim09_order2[[i]], positive=T, usen = 1000, trim = 0.9)
waic_trim09_order2 <- colMeans(matrix(unlist(comp_trim09_order2), ncol = 2, byrow = T))

save(comp_trim09_order2, file = "./real_data_result/Mar_4/comp_trim09_order2.RData")

# trim = 0.8 #

comp_trim08_order2 <- foreach(i = 1:n_chain) %dorng% WAIC(data_trim05$y, data_trim05$x, realdata_trim08_order2[[i]], positive=T, usen = 1000, trim = 0.8)
waic_trim08_order2 <- colMeans(matrix(unlist(comp_trim08_order2), ncol = 2, byrow = T))

save(comp_trim08_order2, file = "./real_data_result/Mar_4/comp_trim08_order2.RData")

# trim = 0.7 #

comp_trim07_order2 <- foreach(i = 1:n_chain) %dorng% WAIC(data_trim05$y, data_trim05$x, realdata_trim07_order2[[i]], positive=T, usen = 1000, trim = 0.7)
waic_trim07_order2 <- colMeans(matrix(unlist(comp_trim07_order2), ncol = 2, byrow = T))

save(comp_trim07_order2, file = "./real_data_result/Mar_4/comp_trim07_order2.RData")

# trim = 0.6 #

comp_trim06_order2 <- foreach(i = 1:n_chain) %dorng% WAIC(data_trim05$y, data_trim05$x, realdata_trim06_order2[[i]], positive=T, usen = 1000, trim = 0.6)
waic_trim06_order2 <- colMeans(matrix(unlist(comp_trim06_order2), ncol = 2, byrow = T))

save(comp_trim06_order2, file = "./real_data_result/Mar_4/comp_trim06_order2.RData")



#################### order = 3 ####################

# no trimming #

comp_notrim_order3 <- foreach(i = 1:n_chain) %dorng% WAIC(data_trim05$y, data_trim05$x, realdata_notrim_order3[[i]], positive=T, usen = 1000, trim = 1)
waic_notrim_order3 <- colMeans(matrix(unlist(comp_notrim_order3), ncol = 2, byrow = T))

save(comp_notrim_order3, file = "./real_data_result/Mar_4/comp_notrim_order3.RData")

# trim = 0.9 #

comp_trim09_order3 <- foreach(i = 1:n_chain) %dorng% WAIC(data_trim05$y, data_trim05$x, realdata_trim09_order3[[i]], positive=T, usen = 1000, trim = 0.9)
waic_trim09_order3 <- colMeans(matrix(unlist(comp_trim09_order3), ncol = 2, byrow = T))

save(comp_trim09_order3, file = "./real_data_result/Mar_4/comp_trim09_order3.RData")

# trim = 0.8 #

comp_trim08_order3 <- foreach(i = 1:n_chain) %dorng% WAIC(data_trim05$y, data_trim05$x, realdata_trim08_order3[[i]], positive=T, usen = 1000, trim = 0.8)
waic_trim08_order3 <- colMeans(matrix(unlist(comp_trim08_order3), ncol = 2, byrow = T))

save(comp_trim08_order3, file = "./real_data_result/Mar_4/comp_trim08_order3.RData")

# trim = 0.7 #

comp_trim07_order3 <- foreach(i = 1:n_chain) %dorng% WAIC(data_trim05$y, data_trim05$x, realdata_trim07_order3[[i]], positive=T, usen = 1000, trim = 0.7)
waic_trim07_order3 <- colMeans(matrix(unlist(comp_trim07_order3), ncol = 2, byrow = T))

save(comp_trim07_order3, file = "./real_data_result/Mar_4/comp_trim07_order3.RData")

# trim = 0.6 #

comp_trim06_order3 <- foreach(i = 1:n_chain) %dorng% WAIC(data_trim05$y, data_trim05$x, realdata_trim06_order3[[i]], positive=T, usen = 1000, trim = 0.6)
waic_trim06_order3 <- colMeans(matrix(unlist(comp_trim06_order3), ncol = 2, byrow = T))

save(comp_trim06_order3, file = "./real_data_result/Mar_4/comp_trim06_order3.RData")



#################### order = 4 ####################

# no trimming #

comp_notrim_order4 <- foreach(i = 1:n_chain) %dorng% WAIC(data_trim05$y, data_trim05$x, realdata_notrim_order4_more[[i]], positive=T, usen = 1000, trim = 1)
waic_notrim_order4 <- colMeans(matrix(unlist(comp_notrim_order4), ncol = 2, byrow = T))

save(comp_notrim_order4, file = "./real_data_result/Mar_4/comp_notrim_order4.RData")

# trim = 0.9 #

comp_trim09_order4 <- foreach(i = 1:n_chain) %dorng% WAIC(data_trim05$y, data_trim05$x, realdata_trim09_order4[[i]], positive=T, usen = 1000, trim = 0.9)
waic_trim09_order4 <- colMeans(matrix(unlist(comp_trim09_order4), ncol = 2, byrow = T))

save(comp_trim09_order4, file = "./real_data_result/Mar_4/comp_trim09_order4.RData")

# trim = 0.8 #

comp_trim08_order4 <- foreach(i = 1:n_chain) %dorng% WAIC(data_trim05$y, data_trim05$x, realdata_trim08_order4[[i]], positive=T, usen = 1000, trim = 0.8)
waic_trim08_order4 <- colMeans(matrix(unlist(comp_trim08_order4), ncol = 2, byrow = T))

save(comp_trim08_order4, file = "./real_data_result/Mar_4/comp_trim08_order4.RData")

# trim = 0.7 #

comp_trim07_order4 <- foreach(i = 1:n_chain) %dorng% WAIC(data_trim05$y, data_trim05$x, realdata_trim07_order4[[i]], positive=T, usen = 1000, trim = 0.7)
waic_trim07_order4 <- colMeans(matrix(unlist(comp_trim07_order4), ncol = 2, byrow = T))

save(comp_trim07_order4, file = "./real_data_result/Mar_4/comp_trim07_order4.RData")

# trim = 0.6 #

comp_trim06_order4 <- foreach(i = 1:n_chain) %dorng% WAIC(data_trim05$y, data_trim05$x, realdata_trim06_order4_more[[i]], positive=T, usen = 1000, trim = 0.6)
waic_trim06_order4 <- colMeans(matrix(unlist(comp_trim06_order4), ncol = 2, byrow = T))

save(comp_trim06_order4, file = "./real_data_result/Mar_4/comp_trim06_order4.RData")



stopCluster(cl)

