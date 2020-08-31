rm(list=ls())
library(MASS)
library(mvtnorm)
library(coda)
library(rjags)


########### Sec 0. Process data
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

data <- list(N0_4=N0_4, N0_2=N0_2, N1_4=N1_4, N1_2=N1_2, 
             R0_4=R0_4, R0_2=R0_2, R1_4=R1_4, R1_2=R1_2, 
             n1=n1, n2=n2, n3=n3, n4=n4)

seed <- 1234
n.chains <- 3
n.adapt <- 1000
n.burnin <- 10000
n.iter <- 100000
n.thin <- 1
set.seed(seed)
init.seeds <- sample(1:1000000, n.chains)

### set initials 
init.jags <- vector("list", n.chains)
for(i in 1:n.chains){
  init.jags[[i]] <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = init.seeds[i])
}

########### Sec 1. Run the final model each forward selection step and save posterior CACE samples
### model 1 
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
            "pic", "pin", "pia", "alpha_n", "alpha_a", "alpha_s1", "alpha_b1") 
jags.m1 <- jags.model(file="model/m1.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.m1, n.iter=n.burnin) # burn in
samps.m1 <- coda.samples.dic(jags.m1, variable.names=params, n.iter=n.iter, thin=n.thin)

samps1 <- samps.m1$samples
samps1_CACE <- rbind(samps1[[1]], samps1[[2]], samps1[[3]])[, c("CACE")] # put three chains together
write.table(samps1_CACE, "output/samps1_CACE.txt", sep="\t")


### Model 2 f) 
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", "pic", "pin", "pia", "sigma_a", 
            "alpha_n", "alpha_a") 
jags.m2f <- jags.model(file="model/m2f.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.m2f, n.iter=n.burnin) # burn in
samps.m2f <- coda.samples.dic(jags.m2f, variable.names=params, n.iter=n.iter, thin=n.thin)

samps2 <- samps.m2f$samples
samps2_CACE <- rbind(samps2[[1]], samps2[[2]], samps2[[3]])[, c("CACE")] # put three chains together
write.table(samps2_CACE, "output/samps2_CACE.txt", sep="\t")


### Model 3 e)
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
            "pic", "pin", "pia", "sigma_n", "sigma_a", 
            "alpha_n", "alpha_a")  
jags.m3e <- jags.model(file="model/m3e.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.m3e, n.iter=n.burnin) # burn in
samps.m3e <- coda.samples.dic(jags.m3e, variable.names=params, n.iter=n.iter, thin=n.thin)

samps3 <- samps.m3e$samples
samps3_CACE <- rbind(samps3[[1]], samps3[[2]], samps3[[3]])[, c("CACE")] # put three chains together
write.table(samps3_CACE, "output/samps3_CACE.txt", sep="\t")


#### Model 4 a) 
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
            "pic", "pin", "pia", "sigma_n", "sigma_a", "sigma_s", 
            "alpha_n", "alpha_a")  
jags.m4a <- jags.model(file="model/m4a.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.m4a, n.iter=n.burnin) # burn in
samps.m4a <- coda.samples.dic(jags.m4a, variable.names=params, n.iter=n.iter, thin=n.thin)

samps4 <- samps.m4a$samples
samps4_CACE <- rbind(samps4[[1]], samps4[[2]], samps4[[3]])[, c("CACE")] # put three chains together
write.table(samps4_CACE, "output/samps4_CACE.txt", sep="\t")


#### Model 5 c) 
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
            "pic", "pin", "pia", "sigma_n", "sigma_a", "sigma_s", "sigma_u", 
            "alpha_n", "alpha_a") 
jags.m5c <- jags.model(file="model/m5c.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.m5c, n.iter=n.burnin) # burn in
samps.m5c <- coda.samples.dic(jags.m5c, variable.names=params, n.iter=n.iter, thin=n.thin)

samps5 <- samps.m5c$samples
samps5_CACE <- rbind(samps5[[1]], samps5[[2]], samps5[[3]])[, c("CACE")] # put three chains together
write.table(samps5_CACE, "output/samps5_CACE.txt", sep="\t")


############ Sec 2. Make Figure 2. Posterior Density

samps1 <- read.table("output/samps1_CACE.txt")
samps2 <- read.table("output/samps2_CACE.txt")
samps3 <- read.table("output/samps3_CACE.txt")
samps4 <- read.table("output/samps4_CACE.txt")
samps5 <- read.table("output/samps5_CACE.txt")
density(samps1$x)

plot(density(samps1$x), col="blue", main="", 
     ylim=c(0, 45), xlim=c(-0.02, 0.11), #ylab="Density", # yaxt='n', 
     xlab=expression(paste("Posterior ", theta^"CACE")))
lines(density(samps2$x), col="red", lty=2)
lines(density(samps3$x), col="brown", lty=3)
lines(density(samps4$x), col="green", lty=4)
lines(density(samps5$x), col="black", lty=5)
legend(0.05, 41, 
       legend=c("Model I (None)", 
                expression(paste("Model IIf (", delta[ia], ")")), 
                expression(paste("Model IIIe (", delta[i][n], ", ", 
                                 delta[ia], ")" )), 
                expression(paste("Model IVa (", delta[i][n], ", ", 
                                 delta[ia], ", ", delta[is], ")" )), 
                expression(paste("Model Vc (", delta[i][n], ", ", 
                                 delta[ia], ", ", delta[is],  ", ", delta[iu], ")" )) ),
       col=c("blue", "red", "brown", "green", "black"), 
       lty=1:5, box.lty=0)

