# package

library(abind)
library(knitr)
library(latex2exp)
require(magrittr)
require(plyr)
library(tidyverse)
library(coda)
library(doParallel)
library(foreach)
library(doRNG)
library(gridExtra)
library(BayesLogit)
par(mar=c(2,2,2,2))

# source

adapt <- FALSE
bunch <- FALSE
N_para <- 10
seed <- 3

# stopCluster(cl)
cl <- makeCluster(10)
registerDoParallel(cl)

source("./codes/new-codes.R")
source("./codes/test-codes.R")
source("./codes/other-codes.R")

# helper function

get.trim.data <- function(trim){
  vote <- subset(vote, abs(from_requirement_threshold) < trim/2)
  preds <- model.matrix(~ #major.requirement 
                          + ISS.against.recommendation 
                        #+ shares.oustanding.base 
                        #+ special.meeting 
                        + analyst.coverage 
                        #+ institutional.ownership 
                        + past.stock.return 
                        + q 
                        + firm.size, 
                        data = vote)[,-1]
  preds.all <- scale(preds)
  x <- cbind(1, preds.all)
  dimnames(x)[[2]][1] <- "Intercept"
  x.names <- dimnames(x)[[2]]
  y <- 2*vote[,"from_requirement_threshold"]
  return(list(x=x, y=y))
}

get.result.CI <- function(result_trim07_order2, N_para = 6, p = 6){
  result_CI_t7o2 <- data.frame()
  for(i in 1:N_para){
    one_set_r <- apply(result_trim07_order2[[i]]$a, 1, function(x) quantile(x, c(0.025, 0.975)))
    for(j in 1:p){
      result_CI_t7o2 <- rbind(result_CI_t7o2,
                              data.frame(lower = one_set_r[1,j],
                                         upper = one_set_r[2,j],
                                         variable = paste0("a", j),
                                         setting = as.character(i)))
    }
  }
  rownames(result_CI_t7o2) <- NULL
  return(result_CI_t7o2)
}


beta3.CI <- function(realdata_trim05_order3, N_para = 10, p = 6){
  b3_CI <- data.frame()
  for(i in 1:N_para){
    b3_cur <- realdata_trim05_order3[[i]]$b[3,,]
    one_set_r <- apply(b3_cur, 1, function(x) quantile(x, c(0.025, 0.975)))
    for(j in 1:p){
      b3_CI <- rbind(b3_CI,
                     data.frame(lower = one_set_r[1,j],
                                upper = one_set_r[2,j],
                                variable = paste0("b3", j),
                                setting = as.character(i)))
    }
  }
  rownames(b3_CI) <- NULL
  return(b3_CI)
}


# real data

vote <- read.csv("./codes/vote_data.csv")
str(vote)
vote <- subset(vote, abs(from_requirement_threshold) < 0.5)
bunch <- FALSE
pos.est <- TRUE
p <- 6


################################# Model #######################################

################## order = 3 ##################

set.seed(seed)
order <- 3
b_init <- replicate(N_para, matrix(rnorm(n = order*p, sd=1), order, p))
a_init <- replicate(N_para, rnorm(n = p, sd = 1.5))

## trimming = 1

data_notrim <- get.trim.data(1)

realdata_notrim_order3 <- foreach(i = 1:N_para,
                                  .packages='BayesLogit') %dorng% 
  bdregjump_adapt_poly_trim_alphanoPG(y=data_notrim$y, x=data_notrim$x,
                                      b=b_init[,,i], burn=15000, nsamp=25000,
                                      thin=1, trim = 1, order = order, prec = 1,
                                      jump=list(a=a_init[,i], prec = 1, positive=pos.est,
                                                persistence=0.5, update.jbw=TRUE))


names = paste0("./real_data_result/Feb_21/realdata_notrim_order3_40k_seed_", seed, ".RData")
save(realdata_notrim_order3, file = names)


realdata_notrim_order3_coda <- multi_chain_coda(realdata_notrim_order3,
                                                burn = 15000, N_para, p, order = order)
plot(realdata_notrim_order3_coda)
gelman.diag(realdata_notrim_order3_coda)
summary(realdata_notrim_order3_coda)

realdata_CI_nto3 <- get.result.CI(realdata_notrim_order3, N_para, p)

## trimming = 0.9

data_trim09 <- get.trim.data(0.9)

realdata_trim09_order3 <- foreach(i = 1:N_para,
                                  .packages='BayesLogit') %dorng% 
  bdregjump_adapt_poly_trim_alphanoPG(y=data_trim09$y, x=data_trim09$x,
                                      b=b_init[,,i], burn=15000, nsamp=25000,
                                      thin=1, trim = 0.9, order = order, prec = 1,
                                      jump=list(a=a_init[,i], prec = 1, positive=pos.est,
                                                persistence=0.5, update.jbw=TRUE))


names = paste0("./real_data_result/Feb_21/realdata_trim09_order3_40k_seed_", seed, ".RData")
save(realdata_trim09_order3, file = names)

realdata_trim09_order3_coda <- multi_chain_coda(realdata_trim09_order3,
                                                burn = 15000, N_para, p, order = order)
plot(realdata_trim09_order3_coda)
gelman.diag(realdata_trim09_order3_coda)
summary(realdata_trim09_order3_coda)

realdata_CI_t9o3 <- get.result.CI(realdata_trim09_order3, N_para, p)



## trimming = 0.8

data_trim08 <- get.trim.data(0.8)

realdata_trim08_order3 <- foreach(i = 1:N_para,
                                  .packages='BayesLogit') %dorng% 
  bdregjump_adapt_poly_trim_alphanoPG(y=data_trim08$y, x=data_trim08$x,
                                      b=b_init[,,i], burn=15000, nsamp=25000,
                                      thin=1, trim = 0.8, order = order, prec = 1,
                                      jump=list(a=a_init[,i], prec = 1, positive=pos.est,
                                                persistence=0.5, update.jbw=TRUE))


names = paste0("./real_data_result/Feb_21/realdata_trim08_order3_40k_seed_", seed, ".RData")
save(realdata_trim08_order3, file = names)

realdata_trim08_order3_coda <- multi_chain_coda(realdata_trim08_order3,
                                                burn = 15000, N_para, p, order = order)
plot(realdata_trim08_order3_coda)
gelman.diag(realdata_trim08_order3_coda)
summary(realdata_trim08_order3_coda)

realdata_CI_t8o3 <- get.result.CI(realdata_trim08_order3, N_para, p)



## trimming = 0.7

data_trim07 <- get.trim.data(0.7)

realdata_trim07_order3 <- foreach(i = 1:N_para,
                                  .packages='BayesLogit') %dorng% 
  bdregjump_adapt_poly_trim_alphanoPG(y=data_trim07$y, x=data_trim07$x,
                                      b=b_init[,,i], burn=15000, nsamp=25000,
                                      thin=1, trim = 0.7, order = order, prec = 1,
                                      jump=list(a=a_init[,i], prec = 1, positive=pos.est,
                                                persistence=0.5, update.jbw=TRUE))


names = paste0("./real_data_result/Feb_21/realdata_trim07_order3_40k_seed_", seed, ".RData")
save(realdata_trim07_order3, file = names)

realdata_trim07_order3_coda <- multi_chain_coda(realdata_trim07_order3,
                                                burn = 15000, N_para, p, order = order)
plot(realdata_trim07_order3_coda)
gelman.diag(realdata_trim07_order3_coda)
summary(realdata_trim07_order3_coda)

realdata_CI_t7o3 <- get.result.CI(realdata_trim07_order3, N_para, p)



## trimming = 0.6

data_trim06 <- get.trim.data(0.6)

realdata_trim06_order3 <- foreach(i = 1:N_para,
                                  .packages='BayesLogit') %dorng% 
  bdregjump_adapt_poly_trim_alphanoPG(y=data_trim06$y, x=data_trim06$x,
                                      b=b_init[,,i], burn=15000, nsamp=25000,
                                      thin=1, trim = 0.6, order = order, prec = 1,
                                      jump=list(a=a_init[,i], prec = 1, positive=pos.est,
                                                persistence=0.5, update.jbw=TRUE))


names = paste0("./real_data_result/Feb_21/realdata_trim06_order3_40k_seed_", seed, ".RData")
save(realdata_trim06_order3, file = names)

realdata_trim06_order3_coda <- multi_chain_coda(realdata_trim06_order3,
                                                burn = 15000, N_para, p, order = order)
plot(realdata_trim06_order3_coda)
gelman.diag(realdata_trim06_order3_coda)
summary(realdata_trim06_order3_coda)

realdata_CI_t6o3 <- get.result.CI(realdata_trim06_order3, N_para, p)


## trimming = 0.5

data_trim05 <- get.trim.data(0.5)

realdata_trim05_order3 <- foreach(i = 1:N_para,
                                  .packages='BayesLogit') %dorng% 
  bdregjump_adapt_poly_trim_alphanoPG(y=data_trim05$y, x=data_trim05$x,
                                      b=b_init[,,i], burn=15000, nsamp=25000,
                                      thin=1, trim = 0.5, order = order, prec = 1,
                                      jump=list(a=a_init[,i], prec = 1, positive=pos.est,
                                                persistence=0.5, update.jbw=TRUE))


names = paste0("./real_data_result/Feb_21/realdata_trim05_order3_40k_seed_", seed, ".RData")
save(realdata_trim05_order3, file = names)

realdata_trim05_order3_coda <- multi_chain_coda(realdata_trim05_order3,
                                                burn = 15000, N_para, p, order = order)
plot(realdata_trim05_order3_coda)
gelman.diag(realdata_trim05_order3_coda)
summary(realdata_trim05_order3_coda)

realdata_CI_t5o3 <- get.result.CI(realdata_trim05_order3, N_para, p)
beta3_CI_t5o3 <- beta3.CI(realdata_trim05_order3, N_para, p)





stopCluster(cl)

############### load ###############

load("./real_data_result/Feb_19/realdata_notrim_order2_40k_seed_3.RData")
load("./real_data_result/Feb_19/realdata_trim05_order2_40k_seed_3.RData")

load("./real_data_result/Feb_21/realdata_trim05_order3_40k_seed_3.RData")

############### plot ###############

realdata_CI_nto3$type = "nto3"
realdata_CI_t9o3$type = "t9o3"
realdata_CI_t8o3$type = "t8o3"
realdata_CI_t7o3$type = "t7o3"
realdata_CI_t6o3$type = "t6o3"
realdata_CI_t5o3$type = "t5o3"

result_CI_o3 <- rbind(realdata_CI_nto3,
                   realdata_CI_t9o3,
                   realdata_CI_t8o3,
                   realdata_CI_t7o3,
                   realdata_CI_t6o3,
                   realdata_CI_t5o3,
                   realdata_CI_t5o3)

result_CI_o3 %>%
  mutate(across(type, factor,
                levels=c("nto3", "t9o3", "t8o3", "t7o3", "t6o3", "t5o3"))) %>%
  ggplot(aes(x = variable, colour = setting)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  facet_wrap(vars(type))

beta3_CI_t5o3 %>%
  ggplot(aes(x = variable, colour = setting)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2)
