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

################## order = 4, more sample ##################

set.seed(seed)
order <- 4
b_init <- replicate(N_para, matrix(rnorm(n = order*p, sd=1), order, p))
a_init <- replicate(N_para, rnorm(n = p, sd = 1.5))

## trimming = 1

data_notrim <- get.trim.data(1)

realdata_notrim_order4_more <- foreach(i = 1:N_para,
                                  .packages='BayesLogit') %dorng% 
  bdregjump_adapt_poly_trim_alphanoPG(y=data_notrim$y, x=data_notrim$x,
                                      b=b_init[,,i], burn=40000, nsamp=25000,
                                      thin=1, trim = 1, order = order, prec = 1,
                                      jump=list(a=a_init[,i], prec = 1, positive=pos.est,
                                                persistence=0.5, update.jbw=TRUE))


names = paste0("./real_data_result/Feb_25/realdata_notrim_order4_65k_seed_", seed, ".RData")
save(realdata_notrim_order4_more, file = names)


realdata_notrim_order4_more_coda <- multi_chain_coda(realdata_notrim_order4_more,
                                                burn = 40000, N_para, p, order = order)
plot(realdata_notrim_order4_more_coda)
gelman.diag(realdata_notrim_order4_more_coda)
summary(realdata_notrim_order4_more_coda)

realdata_CI_nto4 <- get.result.CI(realdata_notrim_order4_more, N_para, p)



## trimming = 0.6

data_trim06 <- get.trim.data(0.6)

realdata_trim06_order4_more <- foreach(i = 1:N_para,
                                  .packages='BayesLogit') %dorng% 
  bdregjump_adapt_poly_trim_alphanoPG(y=data_trim06$y, x=data_trim06$x,
                                      b=b_init[,,i], burn=40000, nsamp=25000,
                                      thin=1, trim = 0.6, order = order, prec = 1,
                                      jump=list(a=a_init[,i], prec = 1, positive=pos.est,
                                                persistence=0.5, update.jbw=TRUE))


names = paste0("./real_data_result/Feb_25/realdata_trim06_order4_65k_seed_", seed, ".RData")
save(realdata_trim06_order4_more, file = names)

realdata_trim06_order4_more_coda <- multi_chain_coda(realdata_trim06_order4_more,
                                                burn = 40000, N_para, p, order = order)
plot(realdata_trim06_order4_more_coda)
gelman.diag(realdata_trim06_order4_more_coda)
summary(realdata_trim06_order4_more_coda)

realdata_CI_t6o4 <- get.result.CI(realdata_trim06_order4_more, N_para, p)


## trimming = 0.5

data_trim05 <- get.trim.data(0.5)

realdata_trim05_order4_less <- foreach(i = 1:N_para,
                                  .packages='BayesLogit') %dorng% 
  bdregjump_adapt_poly_trim_alphanoPG(y=data_trim05$y, x=data_trim05$x,
                                      b=b_init[,,i], burn=0, nsamp=2000,
                                      thin=1, trim = 0.5, order = order, prec = 1,
                                      jump=list(a=a_init[,i], prec = 1, positive=pos.est,
                                                persistence=0.5, update.jbw=TRUE))


names = paste0("./real_data_result/Feb_25/realdata_trim05_order4_2k_seed_", seed, ".RData")
save(realdata_trim05_order4_less, file = names)

realdata_trim05_order4_less_coda <- multi_chain_coda(realdata_trim05_order4_less,
                                                burn = 0, N_para, p, order = order)
plot(realdata_trim05_order4_less_coda)
gelman.diag(realdata_trim05_order4_less_coda)
summary(realdata_trim05_order4_less_coda)


realdata_trim05_order4_more <- list()
realdata_trim05_order4_more[[1]] <- bdregjump_adapt_poly_trim_alphanoPG(y=data_trim05$y, x=data_trim05$x,
                                                                        b=b_init[,,1], burn=40000, nsamp=25000,
                                                                        thin=1, trim = 0.5, order = order, prec = 1,
                                                                        jump=list(a=a_init[,1], prec = 1, positive=pos.est,
                                                                                  persistence=0.5, update.jbw=TRUE),
                                                                        print.process = TRUE)

set.seed(5)

realdata_trim05_order4_more[[2]] <- bdregjump_adapt_poly_trim_alphanoPG(y=data_trim05$y, x=data_trim05$x,
                                                                        b=b_init[,,10], burn=40000, nsamp=25000,
                                                                        thin=1, trim = 0.5, order = order, prec = 1,
                                                                        jump=list(a=a_init[,10], prec = 1, positive=pos.est,
                                                                                  persistence=0.5, update.jbw=TRUE),
                                                                        print.process = TRUE)

set.seed(233)
realdata_trim05_order4_more[[3]] <- bdregjump_adapt_poly_trim_alphanoPG(y=data_trim05$y, x=data_trim05$x,
                                                                        b=b_init[,,1], burn=40000, nsamp=25000,
                                                                        thin=1, trim = 0.5, order = order, prec = 1,
                                                                        jump=list(a=a_init[,1], prec = 1, positive=pos.est,
                                                                                  persistence=0.5, update.jbw=TRUE),
                                                                        print.process = TRUE)


names = paste0("./real_data_result/Feb_25/realdata_trim05_order4_65k_seed_", "sel", ".RData")
save(realdata_trim05_order4_more, file = names)


realdata_trim05_order4_more_coda <- multi_chain_coda(realdata_trim05_order4_more,
                                                     burn = 40000, 2, p, order = order)
plot(realdata_trim05_order4_more_coda)
gelman.diag(realdata_trim05_order4_more_coda)
summary(realdata_trim05_order4_more_coda)


stopCluster(cl)

############### load ###############

load("./real_data_result/Feb_19/realdata_notrim_order2_40k_seed_3.RData")
load("./real_data_result/Feb_19/realdata_trim05_order2_40k_seed_3.RData")

load("./real_data_result/Feb_21/realdata_trim05_order3_40k_seed_3.RData")

load("./real_data_result/Feb_25/realdata_trim05_order4_65k_seed_sel.RData")
load("./real_data_result/Feb_25/realdata_trim06_order4_65k_seed_3.RData")
load("./real_data_result/Feb_25/realdata_notrim_order4_65k_seed_3.RData")

load("./real_data_result/Feb_23/realdata_trim07_order4_40k_seed_3.RData")
load("./real_data_result/Feb_23/realdata_trim08_order4_40k_seed_3.RData")
load("./real_data_result/Feb_23/realdata_trim09_order4_40k_seed_3.RData")

############### plot ###############

realdata_CI_nto4 <- get.result.CI(realdata_notrim_order4_more, N_para, p)
realdata_CI_t9o4 <- get.result.CI(realdata_trim09_order4, N_para, p)
realdata_CI_t8o4 <- get.result.CI(realdata_trim08_order4, N_para, p)
realdata_CI_t7o4 <- get.result.CI(realdata_trim07_order4, N_para, p)
realdata_CI_t6o4 <- get.result.CI(realdata_trim06_order4_more, N_para, p)
realdata_CI_t5o4 <- get.result.CI(realdata_trim05_order4_more, 2, p)


realdata_CI_nto4$type = "nto4"
realdata_CI_t9o4$type = "t9o4"
realdata_CI_t8o4$type = "t8o4"
realdata_CI_t7o4$type = "t7o4"
realdata_CI_t6o4$type = "t6o4"
realdata_CI_t5o4$type = "t5o4"

###

avg_len_o4 <- cbind(
  cal_avg_len(realdata_CI_nto4),
  cal_avg_len(realdata_CI_t9o4),
  cal_avg_len(realdata_CI_t8o4),
  cal_avg_len(realdata_CI_t7o4),
  cal_avg_len(realdata_CI_t6o4),
  cal_avg_len(realdata_CI_t5o4))


rownames(avg_len_o4) <- paste0("a", 1:6)
colnames(avg_len_o4) <- paste0("t=", seq(1,0.5, by = -0.1))

###


result_CI_o4 <- rbind(realdata_CI_nto4,
                   realdata_CI_t9o4,
                   realdata_CI_t8o4,
                   realdata_CI_t7o4,
                   realdata_CI_t6o4,
                   realdata_CI_t5o4)

result_CI_o4 %>%
  mutate(across(type, factor,
                levels=c("nto4", "t9o4", "t8o4", "t7o4", "t6o4", "t5o4"))) %>%
  ggplot(aes(x = variable, colour = setting)) +
  geom_hline(yintercept = 0, colour = 'red') +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  facet_wrap(vars(type))



beta3_CI_nto4 <-beta.CI(realdata_notrim_order4_more, order = 3, N_para, p)
beta3_CI_t9o4 <-beta.CI(realdata_trim09_order4, order = 3, N_para, p)
beta3_CI_t8o4 <-beta.CI(realdata_trim08_order4, order = 3, N_para, p)
beta3_CI_t7o4 <-beta.CI(realdata_trim07_order4, order = 3, N_para, p)
beta3_CI_t6o4 <-beta.CI(realdata_trim06_order4_more, order = 3, N_para, p)
beta3_CI_t5o4 <-beta.CI(realdata_trim05_order4_more, order = 3, 2, p)

beta3_CI_nto4$type = "nto4"
beta3_CI_t9o4$type = "t9o4"
beta3_CI_t8o4$type = "t8o4"
beta3_CI_t7o4$type = "t7o4"
beta3_CI_t6o4$type = "t6o4"
beta3_CI_t5o4$type = "t5o4"

beta3_CI_o4 <- rbind(beta3_CI_nto4,
                     beta3_CI_t9o4,
                     beta3_CI_t8o4,
                     beta3_CI_t7o4,
                     beta3_CI_t6o4,
                     beta3_CI_t5o4)

beta3_CI_o4 %>%
  mutate(across(type, factor,
                levels=c("nto4", "t9o4", "t8o4", "t7o4", "t6o4", "t5o4"))) %>%
  ggplot(aes(x = variable, colour = setting)) +
  geom_hline(yintercept = 0, colour = 'red') +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  facet_wrap(vars(type))



## beta4 o4

beta4_CI_nto4 <-beta.CI(realdata_notrim_order4_more, order = 4, N_para, p)
beta4_CI_t9o4 <-beta.CI(realdata_trim09_order4, order = 4, N_para, p)
beta4_CI_t8o4 <-beta.CI(realdata_trim08_order4, order = 4, N_para, p)
beta4_CI_t7o4 <-beta.CI(realdata_trim07_order4, order = 4, N_para, p)
beta4_CI_t6o4 <-beta.CI(realdata_trim06_order4_more, order = 4, N_para, p)
beta4_CI_t5o4 <-beta.CI(realdata_trim05_order4_more, order = 4, 2, p)

beta4_CI_nto4$type = "nto4"
beta4_CI_t9o4$type = "t9o4"
beta4_CI_t8o4$type = "t8o4"
beta4_CI_t7o4$type = "t7o4"
beta4_CI_t6o4$type = "t6o4"
beta4_CI_t5o4$type = "t5o4"

beta4_CI_o4 <- rbind(beta4_CI_nto4,
                     beta4_CI_t9o4,
                     beta4_CI_t8o4,
                     beta4_CI_t7o4,
                     beta4_CI_t6o4,
                     beta4_CI_t5o4)

beta4_CI_o4 %>%
  mutate(across(type, factor,
                levels=c("nto4", "t9o4", "t8o4", "t7o4", "t6o4", "t5o4"))) %>%
  ggplot(aes(x = variable, colour = setting)) +
  geom_hline(yintercept = 0, colour = 'red') +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  facet_wrap(vars(type))
