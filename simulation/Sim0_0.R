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
sigma_n <- sigma_a <-sigma_u <- sigma_v <- 0.5
beta0 <- beta00 <- beta01 <- 0  # 8.712  # 43.56
beta1 <- 0  # 10  # 50
N0 <- N1 <- rep(175, I)

# assign empty
delta_n <- delta_a <- delta_u <- delta_v <- pi_c <- pi_n <- pi_a <- rep(NA, I)
n <- a <- u1 <- v1 <- s1 <- b1 <- p0 <- p1 <- miss_r0 <- miss_r1 <- rep(NA, I)
prob0 <- prob1 <- R0_4 <- R1_4 <- matrix(NA, nrow = I, ncol = 4)
R0_2 <- R1_2 <- matrix(0, nrow = I, ncol = 2)

# assign empty matrix to store final estimates
sim0_0_pD <- sim0_0_DIC <- sim0_0c_pD <- sim0_0c_DIC <- rep(NA, 1)
sim0_0_CACE <- sim0_0c_CACE <- rep(NA, 10)

seed <- 2018
n.chains <- 3
n.adapt <- 1000
n.burnin <- 10000
n.iter <- 100000
n.thin <- 5
set.seed(seed)

wrap=function(j)
{
  # r <- rep(0.5, I)  # prob of Randomized trt=1
  N0 <- N1 <- rep(175, I) # number of subjects in each arm or each study
  
  # Generate parameters and data, 0) random none
  for (i in 1:I) {
    #delta_n[i] <- rnorm(1, 0, sigma_n)
    #delta_a[i] <- rnorm(1, 0, sigma_a)
    n[i] <- alpha_n
    a[i] <- alpha_a
    
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
  
  sim0_0_CACE[1] <- probit(alpha_u, inverse=T) - probit(alpha_v, inverse=T)
  #sim0_0_CACE[50] <- mean(u1-v1)
  
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
  
  write.table(newdata, paste("simdata0/data0_0_", j,".txt",sep=""), sep="\t")
  
  
  ### set initials ##########################
  set.seed(seed+j*R)
  init.seeds <- sample(1:1000000, n.chains)
  
  Inits <- vector("list", n.chains)
  for(i in 1:n.chains){
    Inits[[i]] <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = init.seeds[i]+j*R)
  }
  
  ############# Sec 1. Model 0) random none ###################################################
  params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
              "pic", "pin", "pia")  
  
  jags.m0_0 <- jags.model(file="sim_m0.txt", data=data, inits=Inits, n.chains=n.chains, n.adapt=n.adapt)
  update(jags.m0_0, n.iter=n.burnin) # burn in
  samps.m0_0 <- coda.samples.dic(jags.m0_0, variable.names=params, n.iter=n.iter, thin=n.thin)
  
  sim0_0_CACE[2] <- summary(samps.m0_0$samples)[[1]][1,1]  # mean CACE
  sim0_0_CACE[3:5] <- c(summary(samps.m0_0$samples)[[2]][1,1], 
                        summary(samps.m0_0$samples)[[2]][1,3], 
                        summary(samps.m0_0$samples)[[2]][1,5])  # 2.5% 50% 97.5% CACE
  sim0_0_CACE[6] <- ifelse((sim0_0_CACE[3]<=sim0_0_CACE[1]) & (sim0_0_CACE[5]>=sim0_0_CACE[1]), 1, 0)
  sim0_0_CACE[7] <- sim0_0_CACE[5] - sim0_0_CACE[3]  # CACE length
  sim0_0_CACE[8] <- (sim0_0_CACE[2] - sim0_0_CACE[1])/sim0_0_CACE[1]  # relative bias
  sim0_0_CACE[9] <- sim0_0_CACE[2] - sim0_0_CACE[1]  # bias
  sim0_0_CACE[10] <- (sim0_0_CACE[2] - sim0_0_CACE[1])^2  # bias^2
  
  sim0_0_pD[1] <- samps.m0_0$dic[[2]]
  sim0_0_DIC[1] <- samps.m0_0$dic[[1]] + samps.m0_0$dic[[2]]
  
  write.table(sim0_0_CACE, paste("sim_CACE/sim0_0_CACE_", j, ".txt", sep=""), sep="\t")
  write.table(sim0_0_pD, paste("sim_pD/sim0_0_pD_", j, ".txt", sep=""), sep="\t")
  write.table(sim0_0_DIC, paste("sim_DIC/sim0_0_DIC_", j, ".txt", sep=""), sep="\t")
  
  
  ############# Sec 2. Model 0) on complete data ###################################################
  cdata <- newdata[1:n1, ]  # only keep the complete data
  
  datac=list(N0=cdata$N0, N1=cdata$N1, R0_4=cdata$R0_4, R1_4=cdata$R1_4, n1=n1)
  
  params <- c("CACE", "u1out", "v1out", "s1out", "b1out", "pic", "pin", "pia")  
  
  jags.m0_0c <- jags.model(file="sim_c0.txt", data=datac, inits=Inits, n.chains=n.chains, n.adapt=n.adapt)
  update(jags.m0_0c, n.iter=n.burnin) # burn in
  samps.m0_0c <- coda.samples.dic(jags.m0_0c, variable.names=params, n.iter=n.iter, thin=n.thin)
  
  sim0_0c_CACE[1] <- sim0_0_CACE[1]
  sim0_0c_CACE[2] <- summary(samps.m0_0c$samples)[[1]][1,1]  # mean CACE
  sim0_0c_CACE[3:5] <- c(summary(samps.m0_0c$samples)[[2]][1,1], 
                        summary(samps.m0_0c$samples)[[2]][1,3], 
                        summary(samps.m0_0c$samples)[[2]][1,5])  # 2.5% 50% 97.5% CACE
  sim0_0c_CACE[6] <- ifelse((sim0_0c_CACE[3]<=sim0_0c_CACE[1]) & (sim0_0c_CACE[5]>=sim0_0c_CACE[1]), 1, 0)
  sim0_0c_CACE[7] <- sim0_0c_CACE[5] - sim0_0c_CACE[3]  # CACE length
  sim0_0c_CACE[8] <- (sim0_0c_CACE[2] - sim0_0c_CACE[1])/sim0_0c_CACE[1]  # relative bias
  sim0_0c_CACE[9] <- sim0_0c_CACE[2] - sim0_0c_CACE[1]  # bias
  sim0_0c_CACE[10] <- (sim0_0c_CACE[2] - sim0_0c_CACE[1])^2  # bias^2
  
  sim0_0c_pD[1] <- samps.m0_0c$dic[[2]]
  sim0_0c_DIC[1] <- samps.m0_0c$dic[[1]] + samps.m0_0c$dic[[2]]
  
  write.table(sim0_0c_CACE, paste("sim_CACE/sim0_0c_CACE_", j, ".txt", sep=""), sep="\t")
  write.table(sim0_0c_pD, paste("sim_pD/sim0_0c_pD_", j, ".txt", sep=""), sep="\t")
  write.table(sim0_0c_DIC, paste("sim_DIC/sim0_0c_DIC_", j, ".txt", sep=""), sep="\t")

  output <- list(sim0_0_CACE, sim0_0_pD, sim0_0_DIC, 
                 sim0_0c_CACE, sim0_0c_pD, sim0_0c_DIC)
  return(output)
}

runsim <- mclapply(1:R, wrap, mc.cores=24)



######### summarize the results ################
sim0_0_CACE <- sim0_0c_CACE <- matrix(NA, R, 10)
sim0_0_DIC <- sim0_0c_DIC <- matrix(NA, R, 1)
colnames(sim0_0_CACE) <- colnames(sim0_0c_CACE) <-
  c("True", "Mean", "Q2.5", "Q50", "Q97.5", "Coverage", "CILength", "Rel_Bias", "Bias", "Bias_sq")

for (i in 1:R){
  sim0_0_CACE[i, ] <- read.table(paste("sim_CACE/sim0_0_CACE_",i,".txt",sep=""))$x
  sim0_0c_CACE[i, ] <- read.table(paste("sim_CACE/sim0_0c_CACE_",i,".txt",sep=""))$x
  
  sim0_0_DIC[i, ] <- read.table(paste("sim_DIC/sim0_0_DIC_",i,".txt",sep=""))$x
  sim0_0c_DIC[i, ] <- read.table(paste("sim_DIC/sim0_0c_DIC_",i,".txt",sep=""))$x
}

write.table(sim0_0_CACE, "simout/CACE0_0.txt", sep="\t", row.names = FALSE)
write.table(sim0_0_DIC, "simout/DIC0_0.txt", sep="\t", row.names = FALSE)
write.table(sim0_0c_CACE, "simout/CACE0_0c.txt", sep="\t", row.names = FALSE)
write.table(sim0_0c_DIC, "simout/DIC0_0c.txt", sep="\t", row.names = FALSE)


