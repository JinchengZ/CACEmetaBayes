rm(list=ls())
library(MASS)
library(mvtnorm)
library(coda)
library(rjags)
library(VGAM)
library(parallel)


source('fun_coda_dic.R')
expit = function(a){
  m = exp(a)/(1 + exp(a))
  return(m)
}

I <- 20  # number of studies in each set
R <- 2000 # number of simulation

# assign true values
alpha_n <- -0.4
alpha_a <- -0.6
alpha_s <- 0.5
alpha_b <- -0.5
alpha_u <- -0.5
alpha_v <- 0.5
sigma_n <- sigma_a <- sigma_u <- sigma_v <- sigma_s <- sigma_b <- 0.5
beta0 <- 0  
beta1 <- 0  # 10  # 50
N0 <- N1 <- rep(175, I)

# assign empty
delta_n <- delta_a <- delta_u <- delta_v <- delta_s <- delta_b <- pi_c <- pi_n <- pi_a <- rep(NA,I)
n <- a <- u1 <- v1 <- s1 <- b1 <- p0 <- p1 <- miss_r0 <- miss_r1 <- rep(NA, I)
prob0 <- prob1 <- R0_4 <- R1_4 <- matrix(NA, nrow = I, ncol = 4)
R0_2 <- R1_2 <- matrix(0, nrow = I, ncol = 2)

# assign empty matrix to store final estimates
select0_1_pD <- select0_1_DIC <- rep(NA, 4)
select0_1_CACE <- rep(NA, 37)

### set parameters to run models ##########################
seed <- 2018
n.chains <- 3
n.adapt <- 1000
n.burnin <- 10000
n.iter <- 100000
n.thin <- 5
set.seed <- seed

wrap=function(j)
{
  # r <- rep(0.5, I)  # prob of Randomized trt=1
  N0 <- N1 <- rep(175, I) # number of subjects in each arm or each study
  
  # Generate parameters and data, 1) random delta_n, delta_a
  for (i in 1:I) {
    delta_n[i] <- rnorm(1, 0, sigma_n)
    delta_a[i] <- rnorm(1, 0, sigma_a)
    n[i] <- alpha_n + delta_n[i]
    a[i] <- alpha_a + delta_a[i]
    
    pi_n[i] = exp(n[i])/(1+exp(n[i])+exp(a[i]))
    pi_a[i] = exp(a[i])/(1+exp(n[i])+exp(a[i]))
    pi_c[i] = 1-pi_a[i]-pi_n[i]

    #delta_u[i] <- rnorm(1, 0, sigma_u)
    
    u1[i] <- probit(alpha_u, inverse=T) 
    v1[i] <- probit(alpha_v, inverse=T)
    s1[i] <- logit(alpha_s, inverse=T) 
    b1[i] <- logit(alpha_b, inverse=T) 
    
    # generate multinomial data
    prob0[i, ] <- c((pi_n[i]*(1-s1[i]) + pi_c[i]*(1-v1[i])), 
                   (pi_n[i]*s1[i] + pi_c[i]*v1[i]), 
                   (pi_a[i]*(1-b1[i])), 
                   (pi_a[i]*b1[i])    )
    prob1[i, ] <- c((pi_n[i]*(1-s1[i])), 
                    (pi_n[i]*s1[i]), 
                    (pi_c[i]*(1-u1[i])+pi_a[i]*(1-b1[i])), 
                    (pi_c[i]*u1[i]+pi_a[i]*b1[i])    )
    
    R0_4[i, ] <- t(rmultinom(1, N0[i], prob0[i, ]))
    R1_4[i, ] <- t(rmultinom(1, N1[i], prob1[i, ]))
    
    if (i<= 10) {
      miss_r0[i]=1
      R0_2[i, 1] <- R0_4[i, 1] + R0_4[i, 3]
      R0_2[i, 2] <- R0_4[i, 2] + R0_4[i, 4]
    }
    else {
      miss_r0[i]=0
      R0_2[i, 1] <- R0_2[i, 2] <- 0
    }
    if (i >= 6 & i<= 15) {
      miss_r1[i]=1
      R1_2[i, 1] <- R1_4[i, 1] + R1_4[i, 3]
      R1_2[i, 2] <- R1_4[i, 2] + R1_4[i, 4]
    }
    else {
      miss_r1[i]=0
      R1_2[i, 1] <- R1_2[i, 2] <- 0
    }
  }
  
  simdata <- cbind(N0, N1, R0_4, R1_4, R0_2, R1_2, miss_r0, miss_r1)
  simdata <- as.data.frame(simdata)
  
  select0_1_CACE[1] <- probit(alpha_u, inverse=T) - probit(alpha_v, inverse=T)
  
  newdata <- simdata[order(simdata[["miss_r0"]], simdata[["miss_r1"]]),]
  n1 <- sum(newdata[["miss_r0"]]==0 & newdata[["miss_r1"]]==0) #  4+4
  n2 <- sum(newdata[["miss_r0"]]==0 & newdata[["miss_r1"]]==1) #  4+2
  n3 <- sum(newdata[["miss_r0"]]==1 & newdata[["miss_r1"]]==0) #  2+4
  n4 <- sum(newdata[["miss_r0"]]==1 & newdata[["miss_r1"]]==1) #  4+4
  n <- length(newdata$miss_r0)
  
  newdata$R0_4 <- cbind(newdata$V3, newdata$V4, newdata$V5, newdata$V6)
  newdata$R1_4 <- cbind(newdata$V7, newdata$V8, newdata$V9, newdata$V10)
  newdata$R0_2 <- cbind(newdata$V11, newdata$V12)
  newdata$R1_2 <- cbind(newdata$V13, newdata$V14)
  
  data=list(N0=newdata$N0, N1=newdata$N1,
            R0_4=newdata$R0_4, R0_2=newdata$R0_2, R1_4=newdata$R1_4, R1_2=newdata$R1_2, 
            n1=n1, n2=n2, n3=n3, n4=n4)

  
  ### set initials ##########################
  init.seeds <- sample(1:1000000, n.chains)
  
  Inits <- vector("list", n.chains)
  for(i in 1:n.chains){
    Inits[[i]] <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = init.seeds[i]+j*R)
  }
  params <- c("CACE", "u1out", "v1out", "s1out", "b1out", "pic", "pin", "pia")  

  ############# Model 0) random none ###################################################
  params <- c("CACE", "u1out", "v1out", "s1out", "b1out", "pic", "pin", "pia")  
  jags.m0 <- jags.model(file="sim_m0.txt", data=data, inits=Inits, n.chains=n.chains, n.adapt=n.adapt)
  update(jags.m0, n.iter=n.burnin) # burn in
  samps.m0 <- coda.samples.dic(jags.m0, variable.names=params, n.iter=n.iter, thin=n.thin)
   
  select0_1_pD[1] <- samps.m0$dic[[2]]
  select0_1_DIC[1] <- samps.m0$dic[[1]] + samps.m0$dic[[2]]
  
  select0_1_CACE[2] <- summary(samps.m0$samples)[[1]][1,1]  # mean CACE
  select0_1_CACE[3:5] <- c(summary(samps.m0$samples)[[2]][1,1], 
                           summary(samps.m0$samples)[[2]][1,3], 
                           summary(samps.m0$samples)[[2]][1,5])  # 2.5% 50% 97.5% CACE
  select0_1_CACE[6] <- ifelse((select0_1_CACE[3]<=select0_1_CACE[1]) & (select0_1_CACE[5]>=select0_1_CACE[1]), 1, 0)
  select0_1_CACE[7] <- select0_1_CACE[5] - select0_1_CACE[3]  # CACE length
  select0_1_CACE[8] <- (select0_1_CACE[2] - select0_1_CACE[1])/select0_1_CACE[1]  # relative bias
  select0_1_CACE[9] <- select0_1_CACE[2] - select0_1_CACE[1]  # bias
  select0_1_CACE[10] <- (select0_1_CACE[2] - select0_1_CACE[1])^2  # bias^2
  
  ############# Model 1). delta_n, delta_a ###################################################
  params <- c("CACE", "u1out", "v1out", "s1out", "b1out", "pic", "pin", "pia")  
  jags.m1 <- jags.model(file="sim_m1.txt", data=data, inits=Inits, n.chains=n.chains, n.adapt=n.adapt)
  update(jags.m1, n.iter=n.burnin) # burn in
  samps.m1 <- coda.samples.dic(jags.m1, variable.names=params, n.iter=n.iter, thin=n.thin)
  
  select0_1_pD[2] <- samps.m1$dic[[2]]
  select0_1_DIC[2] <- samps.m1$dic[[1]] + samps.m1$dic[[2]]
  
  select0_1_CACE[2+9] <- summary(samps.m1$samples)[[1]][1,1]  # mean CACE
  select0_1_CACE[(3+9):(5+9)] <- c(summary(samps.m1$samples)[[2]][1,1], 
                                   summary(samps.m1$samples)[[2]][1,3], 
                                   summary(samps.m1$samples)[[2]][1,5])  # 2.5% 50% 97.5% CACE
  select0_1_CACE[6+9] <- ifelse((select0_1_CACE[3+9]<=select0_1_CACE[1]) & (select0_1_CACE[5+9]>=select0_1_CACE[1]), 1, 0)
  select0_1_CACE[7+9] <- select0_1_CACE[5+9] - select0_1_CACE[3+9]  # CACE length
  select0_1_CACE[8+9] <- (select0_1_CACE[2+9] - select0_1_CACE[1])/select0_1_CACE[1]  # relative bias
  select0_1_CACE[9+9] <- select0_1_CACE[2+9] - select0_1_CACE[1]  # bias
  select0_1_CACE[10+9] <- (select0_1_CACE[2+9] - select0_1_CACE[1])^2  # bias^2
  
  ############# Model 2). delta_n, delta_a, delta_s ###################################################
  params <- c("CACE", "u1out", "v1out", "s1out", "b1out", "pic", "pin", "pia")  
  jags.m2 <- jags.model(file="sim_m2.txt", data=data, inits=Inits, n.chains=n.chains, n.adapt=n.adapt)
  update(jags.m2, n.iter=n.burnin) # burn in
  samps.m2 <- coda.samples.dic(jags.m2, variable.names=params, n.iter=n.iter, thin=n.thin)
    
  select0_1_pD[3] <- samps.m2$dic[[2]]
  select0_1_DIC[3] <- samps.m2$dic[[1]] + samps.m2$dic[[2]]
  
  select0_1_CACE[2+9*2] <- summary(samps.m2$samples)[[1]][1,1]  # mean CACE
  select0_1_CACE[(3+9*2):(5+9*2)] <- c(summary(samps.m2$samples)[[2]][1,1], 
                                       summary(samps.m2$samples)[[2]][1,3], 
                                       summary(samps.m2$samples)[[2]][1,5])  # 2.5% 50% 97.5% CACE
  select0_1_CACE[6+9*2] <- ifelse((select0_1_CACE[3+9*2]<=select0_1_CACE[1]) & (select0_1_CACE[5+9*2]>=select0_1_CACE[1]), 1, 0)
  select0_1_CACE[7+9*2] <- select0_1_CACE[5+9*2] - select0_1_CACE[3+9*2]  # CACE length
  select0_1_CACE[8+9*2] <- (select0_1_CACE[2+9*2] - select0_1_CACE[1])/select0_1_CACE[1]  # relative bias
  select0_1_CACE[9+9*2] <- select0_1_CACE[2+9*2] - select0_1_CACE[1]  # bias
  select0_1_CACE[10+9*2] <- (select0_1_CACE[2+9*2] - select0_1_CACE[1])^2  # bias^2

  ############# Model 3). delta_n, delta_a, delta_s, delta_u ###################################################
  params <- c("CACE", "u1out", "v1out", "s1out", "b1out", "pic", "pin", "pia")  
  jags.m3 <- jags.model(file="sim_m3.txt", data=data, inits=Inits, n.chains=n.chains, n.adapt=n.adapt)
  update(jags.m3, n.iter=n.burnin) # burn in
  samps.m3 <- coda.samples.dic(jags.m3, variable.names=params, n.iter=n.iter, thin=n.thin)
  
  select0_1_pD[4] <- samps.m3$dic[[2]]
  select0_1_DIC[4] <- samps.m3$dic[[1]] + samps.m3$dic[[2]]
  
  select0_1_CACE[2+9*3] <- summary(samps.m3$samples)[[1]][1,1]  # mean CACE
  select0_1_CACE[(3+9*3):(5+9*3)] <- c(summary(samps.m3$samples)[[2]][1,1], 
                                       summary(samps.m3$samples)[[2]][1,3], 
                                       summary(samps.m3$samples)[[2]][1,5])  # 2.5% 50% 97.5% CACE
  select0_1_CACE[6+9*3] <- ifelse((select0_1_CACE[3+9*3]<=select0_1_CACE[1]) & (select0_1_CACE[5+9*3]>=select0_1_CACE[1]), 1, 0)
  select0_1_CACE[7+9*3] <- select0_1_CACE[5+9*3] - select0_1_CACE[3+9*3]  # CACE length
  select0_1_CACE[8+9*3] <- (select0_1_CACE[2+9*3] - select0_1_CACE[1])/select0_1_CACE[1]  # relative bias
  select0_1_CACE[9+9*3] <- select0_1_CACE[2+9*3] - select0_1_CACE[1]  # bias
  select0_1_CACE[10+9*3] <- (select0_1_CACE[2+9*3] - select0_1_CACE[1])^2  # bias^2
  
  write.table(select0_1_CACE, paste("sim_CACE/select0_1_CACE_", j, ".txt", sep=""), sep="\t")
  write.table(select0_1_pD, paste("sim_pD/select0_1_pD_", j, ".txt", sep=""), sep="\t")
  write.table(select0_1_DIC, paste("sim_DIC/select0_1_DIC_", j, ".txt", sep=""), sep="\t")

  output <- list(select0_1_pD, select0_1_DIC)
  return(output)
}

runsim <- mclapply(1:R, wrap, mc.cores=24)


######### summarize the results ################
select0_1_DIC <- matrix(NA, R, 4)
select0_1_pD <- matrix(NA, R, 4)
select0_1_loc <- rep(NA, R)
select0_1_CACE <- matrix(NA, R, 37)
colnames(select0_1_CACE) <- 
  c("True", "Mean", "Q2.5", "Q50", "Q97.5", "Coverage", "CILength", "Rel_Bias", "Bias", "Bias_sq", 
    "Mean_1", "Q2.5_1", "Q50_1", "Q97.5_1", "Coverage_1", "CILength_1", "Rel_Bias_1", "Bias_1", "Bias_sq_1", 
    "Mean_2", "Q2.5_2", "Q50_2", "Q97.5_2", "Coverage_2", "CILength_2", "Rel_Bias_2", "Bias_2", "Bias_sq_2", 
    "Mean_3", "Q2.5_3", "Q50_3", "Q97.5v", "Coverage_3", "CILength_3", "Rel_Bias_3", "Bias_3", "Bias_sq_3")

for (i in 1:R){
  select0_1_DIC[i, ] <- read.table(paste("sim_DIC/select0_1_DIC_",i,".txt",sep=""))$x
  select0_1_loc[i] <- which.min(select0_1_DIC[i,])
  select0_1_pD[i, ] <- read.table(paste("sim_pD/select0_1_pD_",i,".txt",sep=""))$x
  
  select0_1_CACE[i, ] <- read.table(paste("sim_CACE/select0_1_CACE_",i,".txt",sep=""))$x
}

table0_1 <- matrix(NA, 1, 4)
for (i in 1:4) {
  table0_1[i] <- length(which(select0_1_loc == i)) 
}

selectCACE0_1 <- colMeans(select0_1_CACE)

# Save the data
write.table(table0_1, "simout/table0_1.txt", sep="\t", row.names = FALSE)
write.table(selectCACE0_1, "simout/selectCACE0_1.txt", sep="\t", row.names = FALSE)