############### Figure 4. CACE MNAR Sens plot 

########## Sec 1. Save the CI_CACE data
rm(list=ls())
library(MASS)
library(mvtnorm)
library(coda)
library(rjags)
library(HDInterval)
library(parallel)


source('fun_coda_dic.R')
obs <- read.table("data/data.txt",header=T)
attach(obs)

n1 <- sum(obs$miss_r0==0 & obs$miss_r1==0) #  4+4
n2 <- sum(obs$miss_r0==0 & obs$miss_r1==1) #  4+2
n3 <- sum(obs$miss_r0==1 & obs$miss_r1==0) #  2+4
n4 <- sum(obs$miss_r0==1 & obs$miss_r1==1) #  4+4
n <- length(obs$Author)
n1+n2+n3+n4 == n  # check if all studies included

Ntol <- n000+n001+n010+n011+n100+n101+n110+n111+n1s1+n1s0+n0s1+n0s0
N0_4 <- n000+n001+n010+n011
N0_2 <- n0s1+n0s0
N1_4 <- n100+n101+n110+n111
N1_2 <- n1s1+n1s0

R <- cbind(n000,n001,n010,n011, n100,n101,n110,n111, n0s0,n0s1,n1s0,n1s1)
R0_4 <- cbind(n000,n001,n010,n011)
R0_2 <- cbind(n0s0,n0s1)
R1_4 <- cbind(n100,n101,n110,n111)
R1_2 <- cbind(n1s0,n1s1)
r <- (n100+n101+n110+n111+n1s1+n1s0)/Ntol
sum(Ntol)/n

R1 <- cbind(R0_4, R1_4)
R2 <- cbind(R0_4, R1_2)
R3 <- cbind(R0_2, R1_4)
R4 <- cbind(R0_2, R1_2)

data0 <- list(N0_4=N0_4, N0_2=N0_2, N1_4=N1_4, N1_2=N1_2, 
              R0_4=R0_4, R0_2=R0_2, R1_4=R1_4, R1_2=R1_2, 
              n1=n1, n2=n2, n3=n3, n4=n4,
              miss_r0=miss_r0, miss_r1=miss_r1)

seed <- 2017
n.chains <- 3
n.adapt <- 1000
n.burnin <- 10000
n.iter <- 100000
n.thin <- 1
set.seed(seed)
init.seeds <- sample(1:1000000, n.chains)

### set initials ##########################
init.jags <- vector("list", n.chains)
for(i in 1:n.chains){
  init.jags[[i]] <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = init.seeds[i])
}

### set sequence ###################
seq1 <- seq(-2.5, 2.5, by = 0.05)
seq0 <- 0
num <- length(seq1)


wrap = function(i)
{  
  gamma_10 <- seq0
  gamma_11 <- seq1[i]
  
  data=list(N0_4=N0_4, N0_2=N0_2, N1_4=N1_4, N1_2=N1_2, 
            R0_4=R0_4, R0_2=R0_2, R1_4=R1_4, R1_2=R1_2, 
            n1=n1, n2=n2, n3=n3, n4=n4,
            miss_r0=miss_r0, miss_r1=miss_r1,
            gamma_10=gamma_10, gamma_11=gamma_11)
  params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
              "pic", "pin", "pia", "sigma_n", "sigma_a", "sigma_s", "sigma_u", 
              "alpha_n", "alpha_a") 
  
  jags.m5c_mnar = jags.model(file='m5c_mnar.txt', data=data, 
                             inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
  update(jags.m5c_mnar, n.iter=n.burnin) # burn in
  samps.m5c_mnar <- coda.samples.dic(jags.m5c_mnar, variable.names=params, 
                                     n.iter=n.iter, thin=n.thin)
  
  write.table(summary(samps.m5c_mnar$samples)[[2]], 
              paste("output/MNARout2_5c_", i, ".txt", sep=""), sep="\t")
  write.table(t(hdi(samps.m5c_mnar$samples)), 
              paste("output/MNARHPD_5c_", i, ".txt", sep=""), sep="\t")
  
  return(0)
}

runsim <- mclapply(1:num, wrap, mc.cores=24)


# save the estimated CACE from sensitivity analysis
seq1 <- seq(-2.5, 2.5, by = 0.05)
seq0 <- 0
R <- length(seq1)
CI_CACE <- matrix(NA, R, 7) 
for (i in 1:R){
  value <- -2.5 + 0.05*(i-1)
  CI_CACE[i, 1] <- 0
  CI_CACE[i, 2] <- value
  CI_CACE[i, 3:5] <- t(read.table(paste("output/MNARout2_5c_",i,".txt",sep=""))
                       [c("CACE"), c("X2.5.", "X50.", "X97.5.")])
  CI_CACE[i, 6] <- read.table(paste("output/MNARHPD_5c_",i,".txt",sep=""))[c("CACE"), c("lower")]
  CI_CACE[i, 7] <- read.table(paste("output/MNARHPD_5c_",i,".txt",sep=""))[c("CACE"), c("upper")]
}
colnames(CI_CACE) <- c("gamma10", "gamma11", "lower", "median", "upper", "HPD_lower", "HPD_upper")
write.table(CI_CACE, "output/CI_CACE.txt", sep="\t")


########### Sec 2. Plot Figure 4. MNAR

CACE <- read.table("output/CI_CACE.txt", sep="\t")

value <- CACE$gamma11
est <- CACE$median
lower <- CACE$lower
upper <- CACE$upper
HPD_lower <- CACE$HPD_lower
HPD_upper <- CACE$HPD_upper

par(mar=c(4.3, 5.2, 1.5, 1.5))
n <- length(est)
plot(value, est, type="l",xlab=expression(gamma[11]), cex.lab=2, cex.axis=1.5,
     ylab=expression(paste(theta^"CACE")), ylim=c(-.03, 0.18), xaxs="i")  #,xaxt="n"
lines(value, est, lwd = 2)
lines(value, lower)
lines(value, upper)
lines(value, HPD_lower, lty = 2)
lines(value, HPD_upper, lty = 2)
abline(h = 0, lty = 3, col="blue")