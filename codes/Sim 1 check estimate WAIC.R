rm(list = ls())

###################### load data ######################
load("~/Desktop/Density Estimation/Github code 3-16/drdisc-main/sim_20k_case_b/sim_20k_res/sim_notrim_order2_20k_overall.RData")
order2_colnames = colnames(sim_order2_20k_notrim_order2_summary)
rm(list=setdiff(ls(), c("sim_order2_20k_notrim_order2", "order2_colnames")))
load("~/Desktop/Density Estimation/Github code 3-16/drdisc-main/sim_20k_case_b/sim_20k_res/sim_notrim_order3_20k_overall.RData")
order3_colnames = colnames(sim_order2_20k_notrim_order3_summary)
rm(list=setdiff(ls(), c("sim_order2_20k_notrim_order2", "sim_order2_20k_notrim_order3","order2_colnames", "order3_colnames")))
load("~/Desktop/Density Estimation/Github code 3-16/drdisc-main/sim_20k_case_b/sim_20k_res/sim_notrim_order4_20k_overall.RData")

###################### obtain WAIC for 100 replicates ######################
WAIC_order2 = sim_order2_20k_notrim_order2[1,11,]
WAIC_order3 = sim_order2_20k_notrim_order3[1,13,]
WAIC_order4 = sim_order2_20k_notrim_order4[1,15,]
WAIC_order234 = cbind(WAIC_order2, WAIC_order3)
WAIC_order234 = cbind(WAIC_order234, WAIC_order4)

table(apply(WAIC_order234, 1, which.min))/100
# 63%, 24%, 13%
table(apply(WAIC_order234, 1, which.max))/100
# 19%, 21%, 60%

###################### obtain adaptive WAIC estimates ######################

a.WAIC.adpt <- array(0, c(3, 2, 100))
#apply(WAIC_order234, 1, which.min)
for(j in 1:100){
  ttt = (apply(WAIC_order234, 1, which.min))[j]
  if(ttt == 1){
    a.WAIC.adpt[,,j] = sim_order2_20k_notrim_order2[, 1:2, j]
  } else if(ttt == 2){
    a.WAIC.adpt[,,j] = sim_order2_20k_notrim_order3[, 1:2, j]
  } else {
    a.WAIC.adpt[,,j] = sim_order2_20k_notrim_order4[, 1:2, j]
  }
}

j = 1
length_95 = c()
MSE_mean = c()
coverage = c()
truth_2  = c(1, 4, 1.2598, 0, -1.6498, -1.3296, 0.16, 4, 1) 
for(i in 1:100){
  lower_par = a.WAIC.adpt[1,j,i]
  mean_par = a.WAIC.adpt[2,j,i]
  upper_par = a.WAIC.adpt[3,j,i]
  
  length_95[i] = upper_par - lower_par
  MSE_mean[i] = ( mean_par - truth_2[j] )^2
  coverage[i] = (truth_2[j] > lower_par) & (truth_2[j] < upper_par)
}
mean(length_95)
mean(MSE_mean)
mean(coverage)


j = 2
length_95 = c()
MSE_mean = c()
coverage = c()
truth_2  = c(1, 4, 1.2598, 0, -1.6498, -1.3296, 0.16, 4, 1) 
for(i in 1:100){
  lower_par = a.WAIC.adpt[1,j,i]
  mean_par = a.WAIC.adpt[2,j,i]
  upper_par = a.WAIC.adpt[3,j,i]
  
  length_95[i] = upper_par - lower_par
  MSE_mean[i] = ( mean_par - truth_2[j] )^2
  coverage[i] = (truth_2[j] > lower_par) & (truth_2[j] < upper_par)
}
mean(length_95)
mean(MSE_mean)
mean(coverage)

###################### plot WAIC curve ######################
plot(1:100, WAIC_order2, 'l', ylab = 'WAIC', xlab = 'Replicates')
lines(1:100, WAIC_order3, 'l', col = 'red', lty = 2)
lines(1:100, WAIC_order4, 'l', col = 'green', lty = 3)
legend("topleft", legend=c("Order 2", "Order 3", "Order 4"),
       col=c("black", "red", "green"), lty=1:3, cex=1)


###################### Obtain numbers for parameters in the table ######################

truth_2  = c(1, 4, 1.2598, 0, -1.6498, -1.3296, 0.16, 4, 1) 
truth_3  = c(1, 4, 1.2598, 0, -1.6498, -1.3296, 0, 0, 0.16, 4, 1) 
truth_4  = c(1, 4, 1.2598, 0, -1.6498, -1.3296, 0, 0, 0, 0, 0.16, 4, 1) 

order2_output = matrix(NA, 9, 3)
for(j in 1:9){
  length_95 = c()
  MSE_mean = c()
  coverage = c()
  for(i in 1:100){
    lower_par = sim_order2_20k_notrim_order2[1,j,i]
    mean_par = sim_order2_20k_notrim_order2[2,j,i]
    upper_par = sim_order2_20k_notrim_order2[3,j,i]
    
    length_95[i] = upper_par - lower_par
    MSE_mean[i] = ( mean_par - truth_2[j] )^2
    coverage[i] = (truth_2[j] > lower_par) & (truth_2[j] < upper_par)
  }
  
  order2_output[j, 1] = mean(length_95)
  order2_output[j, 2] = mean(MSE_mean)
  order2_output[j, 3] = mean(coverage)
}
rownames(order2_output) = order2_colnames[1:9]
colnames(order2_output) = c("Avg Length", "MSE", "Coverage Rate")
round(order2_output, 3)

order3_output = matrix(NA, 11, 3)
for(j in 1:11){
  length_95 = c()
  MSE_mean = c()
  coverage = c()
  for(i in 1:100){
    lower_par = sim_order2_20k_notrim_order3[1,j,i]
    mean_par = sim_order2_20k_notrim_order3[2,j,i]
    upper_par = sim_order2_20k_notrim_order3[3,j,i]
    
    length_95[i] = upper_par - lower_par
    MSE_mean[i] = ( mean_par - truth_3[j] )^2
    coverage[i] = (truth_3[j] > lower_par) & (truth_3[j] < upper_par)
  }
  
  order3_output[j, 1] = mean(length_95)
  order3_output[j, 2] = mean(MSE_mean)
  order3_output[j, 3] = mean(coverage)
}
rownames(order3_output) = order3_colnames[1:11]
colnames(order3_output) = c("Avg Length", "MSE", "Coverage Rate")
round(order3_output, 3)
knitr::kable(round(order3_output, 3), 'latex')

order4_output = matrix(NA, 13, 3)
for(j in 1:13){
  length_95 = c()
  MSE_mean = c()
  coverage = c()
  for(i in 1:100){
    lower_par = sim_order2_20k_notrim_order4[1,j,i]
    mean_par = sim_order2_20k_notrim_order4[2,j,i]
    upper_par = sim_order2_20k_notrim_order4[3,j,i]
    
    length_95[i] = upper_par - lower_par
    MSE_mean[i] = ( mean_par - truth_4[j] )^2
    coverage[i] = (truth_4[j] > lower_par) & (truth_4[j] < upper_par)
  }
  
  order4_output[j, 1] = mean(length_95)
  order4_output[j, 2] = mean(MSE_mean)
  order4_output[j, 3] = mean(coverage)
}
rownames(order4_output) = order4_colnames[1:13]
colnames(order4_output) = c("Avg Length", "MSE", "Coverage Rate")
round(order4_output, 3)

knitr::kable(round(order4_output, 3), 'latex')

