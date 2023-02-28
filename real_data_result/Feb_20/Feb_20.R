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

# real data

vote <- read.csv("./codes/vote_data.csv")
str(vote)
vote <- subset(vote, abs(from_requirement_threshold) < 0.5)
bunch <- FALSE
pos.est <- TRUE
p <- 6


################################# Model #######################################

################## order = 2 ##################


## trimming = 0.9

set.seed(seed)
order <- 2
b_init <- replicate(N_para, matrix(rnorm(n = order*p, sd=1), order, p))
a_init <- replicate(N_para, rnorm(n = p, sd = 1.5))

data_trim09 <- get.trim.data(0.9)

realdata_trim09_order2 <- foreach(i = 1:N_para,
                                  .packages='BayesLogit') %dorng% 
  bdregjump_adapt_poly_trim_alphanoPG(y=data_trim09$y, x=data_trim09$x,
                                      b=b_init[,,i], burn=15000, nsamp=25000,
                                      thin=1, trim = 0.9, order = order, prec = 1,
                                      jump=list(a=a_init[,i], prec = 1, positive=pos.est,
                                                persistence=0.5, update.jbw=TRUE))


names = paste0("./real_data_result/Feb_20/realdata_trim09_order2_40k_seed_", seed, ".RData")
save(realdata_trim09_order2, file = names)

realdata_trim09_order2_coda <- multi_chain_coda(realdata_trim09_order2,
                                                burn = 15000, N_para, p, order = 2)
plot(realdata_trim09_order2_coda)
gelman.diag(realdata_trim09_order2_coda)
summary(realdata_trim09_order2_coda)

realdata_CI_t9o2 <- get.result.CI(realdata_trim09_order2, N_para, p)



## trimming = 0.8

set.seed(seed)
order <- 2
b_init <- replicate(N_para, matrix(rnorm(n = order*p, sd=1), order, p))
a_init <- replicate(N_para, rnorm(n = p, sd = 1.5))

data_trim08 <- get.trim.data(0.8)

realdata_trim08_order2 <- foreach(i = 1:N_para,
                                  .packages='BayesLogit') %dorng% 
  bdregjump_adapt_poly_trim_alphanoPG(y=data_trim08$y, x=data_trim08$x,
                                      b=b_init[,,i], burn=15000, nsamp=25000,
                                      thin=1, trim = 0.8, order = order, prec = 1,
                                      jump=list(a=a_init[,i], prec = 1, positive=pos.est,
                                                persistence=0.5, update.jbw=TRUE))


names = paste0("./real_data_result/Feb_20/realdata_trim08_order2_40k_seed_", seed, ".RData")
save(realdata_trim08_order2, file = names)

realdata_trim08_order2_coda <- multi_chain_coda(realdata_trim08_order2,
                                                burn = 15000, N_para, p, order = 2)
plot(realdata_trim08_order2_coda)
gelman.diag(realdata_trim08_order2_coda)
summary(realdata_trim08_order2_coda)

realdata_CI_t8o2 <- get.result.CI(realdata_trim08_order2, N_para, p)



## trimming = 0.7

set.seed(seed)
order <- 2
b_init <- replicate(N_para, matrix(rnorm(n = order*p, sd=1), order, p))
a_init <- replicate(N_para, rnorm(n = p, sd = 1.5))

data_trim07 <- get.trim.data(0.7)

realdata_trim07_order2 <- foreach(i = 1:N_para,
                                  .packages='BayesLogit') %dorng% 
  bdregjump_adapt_poly_trim_alphanoPG(y=data_trim07$y, x=data_trim07$x,
                                      b=b_init[,,i], burn=15000, nsamp=25000,
                                      thin=1, trim = 0.7, order = order, prec = 1,
                                      jump=list(a=a_init[,i], prec = 1, positive=pos.est,
                                                persistence=0.5, update.jbw=TRUE))


names = paste0("./real_data_result/Feb_20/realdata_trim07_order2_40k_seed_", seed, ".RData")
save(realdata_trim07_order2, file = names)

realdata_trim07_order2_coda <- multi_chain_coda(realdata_trim07_order2,
                                                burn = 15000, N_para, p, order = 2)
plot(realdata_trim07_order2_coda)
gelman.diag(realdata_trim07_order2_coda)
summary(realdata_trim07_order2_coda)

realdata_CI_t7o2 <- get.result.CI(realdata_trim07_order2, N_para, p)



## trimming = 0.6

set.seed(seed)
order <- 2
b_init <- replicate(N_para, matrix(rnorm(n = order*p, sd=1), order, p))
a_init <- replicate(N_para, rnorm(n = p, sd = 1.5))

data_trim06 <- get.trim.data(0.6)

realdata_trim06_order2 <- foreach(i = 1:N_para,
                                  .packages='BayesLogit') %dorng% 
  bdregjump_adapt_poly_trim_alphanoPG(y=data_trim06$y, x=data_trim06$x,
                                      b=b_init[,,i], burn=15000, nsamp=25000,
                                      thin=1, trim = 0.6, order = order, prec = 1,
                                      jump=list(a=a_init[,i], prec = 1, positive=pos.est,
                                                persistence=0.5, update.jbw=TRUE))


names = paste0("./real_data_result/Feb_20/realdata_trim06_order2_40k_seed_", seed, ".RData")
save(realdata_trim06_order2, file = names)

realdata_trim06_order2_coda <- multi_chain_coda(realdata_trim06_order2,
                                                burn = 15000, N_para, p, order = 2)
plot(realdata_trim06_order2_coda)
gelman.diag(realdata_trim06_order2_coda)
summary(realdata_trim06_order2_coda)

realdata_CI_t6o2 <- get.result.CI(realdata_trim06_order2, N_para, p)





stopCluster(cl)

############### load ###############

load("./real_data_result/Feb_19/realdata_notrim_order2_40k_seed_3.RData")
load("./real_data_result/Feb_19/realdata_trim05_order2_40k_seed_3.RData")

load("./real_data_result/Feb_20/realdata_trim09_order2_40k_seed_3.RData")
load("./real_data_result/Feb_20/realdata_trim08_order2_40k_seed_3.RData")
load("./real_data_result/Feb_20/realdata_trim07_order2_40k_seed_3.RData")
load("./real_data_result/Feb_20/realdata_trim06_order2_40k_seed_3.RData")

############### plot ###############

realdata_CI_nto2$type = "nto2"
realdata_CI_t9o2$type = "t9o2"
realdata_CI_t8o2$type = "t8o2"
realdata_CI_t7o2$type = "t7o2"
realdata_CI_t6o2$type = "t6o2"
realdata_CI_t5o2$type = "t5o2"

avg_len_nto2 <- realdata_CI_nto2 %>%
  mutate(length = upper - lower) %>%
  group_by(variable) %>%
  summarize(avg_len = mean(length))

cal_avg_len <- function(realdata_CI_nto2){
  avg_len_nto2 <- realdata_CI_nto2 %>%
    mutate(length = upper - lower) %>%
    group_by(variable) %>%
    summarize(avg_len = mean(length))
  
  return(avg_len_nto2[,2])
}

avg_len_o2 <- cbind(
cal_avg_len(realdata_CI_nto2),
cal_avg_len(realdata_CI_t9o2),
cal_avg_len(realdata_CI_t8o2),
cal_avg_len(realdata_CI_t7o2),
cal_avg_len(realdata_CI_t6o2),
cal_avg_len(realdata_CI_t5o2))


rownames(avg_len_o2) <- paste0("a", 1:6)
colnames(avg_len_o2) <- paste0("t=", seq(1,0.5, by = -0.1))


result_CI_o2 <- rbind(realdata_CI_nto2,
                   realdata_CI_t9o2,
                   realdata_CI_t8o2,
                   realdata_CI_t7o2,
                   realdata_CI_t6o2,
                   realdata_CI_t5o2)


result_CI_o2 %>%
  mutate(across(type, factor,
                levels=c("nto2", "t9o2", "t8o2", "t7o2", "t6o2", "t5o2"))) %>%
  ggplot(aes(x = variable, colour = setting)) +
  geom_hline(yintercept = 0, colour = 'red') +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  facet_wrap(vars(type))
