rm(list=ls())
library(MASS)
library(mvtnorm)
library(coda)
library(rjags)


########### Sec 0. Load data
source('fun_coda_dic.R')
obs13_final <- read.table("data/obs13_final.txt", sep="\t", header=T)
M <- length(unique(obs13_final$id))
N <- length(obs13_final$id)
data <- list(id=obs13_final$id, r=obs13_final$rr, y=obs13_final$y, n=obs13_final$n,
             M = M, N = N)

################# Sec 1. CR analysis 
seed <- 2017
n.chains <- 3
n.adapt <- 1000
n.burnin <- 10000
n.iter <- 40000
n.thin <- 1
set.seed(seed)
init.seeds <- sample(1:100000, n.chains)

### set initials 
init.jags <- vector("list", n.chains)
for(i in 1:n.chains){
  init.jags[[i]] <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = init.seeds[i])
}

params <- c("cp", "mu[1]", "mu[2]", "Sigma", "uout", "vout", "diff") 
jags.m_CR <- jags.model(file="model/m_CR_long.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.m_CR, n.iter=n.burnin) # burn in
samps.m_CR <- coda.samples.dic(jags.m_CR, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.m_CR$samples)[[1]], "output/out1_CR_long.txt", sep="\t")
write.table(summary(samps.m_CR$samples)[[2]], "output/out2_CR_long.txt", sep="\t")


############## Sec 2. ER analysis 
obs <- read.table("data/data_ITT.txt", header=T)
attach(obs)

n <- length(obs$Author_Year)
data <- list(n1s1=n1s1, n1ss=n1ss, n0s1=n0s1, n0ss=n0ss, n=n)

seed <- 2017
n.chains <- 3
n.adapt <- 1000
n.burnin <- 10000
n.iter <- 40000
n.thin <- 1
set.seed(seed)
init.seeds <- sample(1:100000, n.chains)

### set initials 
init.jags <- vector("list", n.chains)
for(i in 1:n.chains){
  init.jags[[i]] <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = init.seeds[i])
}

### model ER 
params <- c("p1s1", "p0s1", "u", "v", "Sigma", "uout", "vout", "diff") 
jags.m_ER <- jags.model(file="model/m_ER.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.m_ER, n.iter=n.burnin) # burn in
samps.m_ER <- coda.samples.dic(jags.m_ER, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.m_ER$samples)[[1]], "output/out1_ER.txt", sep="\t")
write.table(summary(samps.m_ER$samples)[[2]], "output/out2_ER.txt", sep="\t")



############# Sec 3. Data ERCR for plot (go to local computer ) 
out2_ER <- read.table("output/out2_ER.txt", sep="\t")
ER_out1 <- out2_ER[c("p1s1[1]", "p1s1[2]", "p1s1[3]", "p1s1[4]", "p1s1[5]", 
                     "p1s1[6]", "p1s1[7]", "p1s1[8]", "p1s1[9]", "p1s1[10]", 
                     "p1s1[11]", "p1s1[12]", "p1s1[13]", "p1s1[14]", "p1s1[15]", "p1s1[16]",  
                     "p1s1[17]", "p1s1[18]", "p1s1[19]", "p1s1[20]", "p1s1[21]", "p1s1[22]", 
                     "p1s1[23]", "p1s1[24]", "p1s1[25]", "p1s1[26]", "p1s1[27]"), 
                   c("X2.5.", "X50.", "X97.5.")]

ER_out0 <- out2_ER[c("p0s1[1]", "p0s1[2]", "p0s1[3]", "p0s1[4]", "p0s1[5]", 
                     "p0s1[6]", "p0s1[7]", "p0s1[8]", "p0s1[9]", "p0s1[10]", 
                     "p0s1[11]", "p0s1[12]", "p0s1[13]", "p0s1[14]", "p0s1[15]", "p0s1[16]",  
                     "p0s1[17]", "p0s1[18]", "p0s1[19]", "p0s1[20]", "p0s1[21]", "p0s1[22]", 
                     "p0s1[23]", "p0s1[24]", "p0s1[25]", "p0s1[26]", "p0s1[27]"), 
                   c("X2.5.", "X50.", "X97.5.")]
data_ER <- cbind(ER_out1, ER_out0)
names(data_ER) <- c("p1s1L", "p1s1M", "p1s1U", "p0s1L", "p0s1M", "p0s1U")


out2_CR <- read.table("output/out2_CR_long.txt", sep="\t")
rownames(out2_CR)
temp <- out2_CR[c(5:27), c("X2.5.", "X50.", "X97.5.")]
obs13_final <- read.table("data/obs13_final.txt", sep="\t")
CR1 <- cbind(obs13_final[c("id", "r")], temp)

obs_temp <- read.table("data/data_plot13.txt", header=T)
obs13 <- obs_temp[c("Author_Year", "n01s", "n0ss", "n10s", "n1ss", "miss_r0", "miss_r1", "study")]
obs13$id <- 1:nrow(obs13) 
obs_r0 <- obs13[c("id", "Author_Year", "n01s", "n0ss", "miss_r0", "study")]
names(obs_r0) <- c("id", "Author_Year", "y", "n", "mis", "study")
obs_r0$r <- 0
obs_r1 <- obs13[c("id", "Author_Year", "n10s", "n1ss", "miss_r1", "study")]
names(obs_r1) <- c("id", "Author_Year", "y", "n", "mis", "study")
obs_r1$r <- 1

CR2 <- rbind(obs_r0, obs_r1)
CR2 <- CR2[order(CR2[["id"]], CR2[["r"]]),]

CR3 <- merge(CR2, CR1, by=c("id","r"), all.x = T)
out_CR <- CR3[order(CR3[["r"]], CR3[["id"]]),]

CR_out1 <- out_CR[c(14:26), c("X2.5.", "X50.", "X97.5.")]
CR_out0 <- out_CR[c(1:13), c("X2.5.", "X50.", "X97.5.")]

data_CR <- cbind(CR_out1, CR_out0)
names(data_CR) <- c("c10sL", "c10sM", "c10sU", "c01sL", "c01sM", "c01sU")


data_ITT <- read.table("data/data_ITT.txt", sep='\t')
ER <- cbind(data_ITT, data_ER)
overall_p1s1 <- paste(formatC(round(out2_ER["uout", "X50."]*100, 3), format='f', digits=1 ), " (",
                      formatC(round(out2_ER["uout", "X2.5."]*100, 3), format='f', digits=1 ),", ", 
                      formatC(round(out2_ER["uout", "X97.5."]*100, 3), format='f', digits=1 ), ")", sep="")
overall_p0s1 <- paste(formatC(round(out2_ER["vout", "X50."]*100, 3), format='f', digits=1 ), " (",
                      formatC(round(out2_ER["vout", "X2.5."]*100, 3), format='f', digits=1 ),", ", 
                      formatC(round(out2_ER["vout", "X97.5."]*100, 3), format='f', digits=1 ), ")", sep="")
data_plot <- read.table("data/data_plot13.txt",header=T)
CR <- cbind(data_plot, data_CR)
overall_c10s <- paste(formatC(round(out2_CR["uout", "X50."]*100, 3), format='f', digits=1 ), " (",
                      formatC(round(out2_CR["uout", "X2.5."]*100, 3), format='f', digits=1 ),", ", 
                      formatC(round(out2_CR["uout", "X97.5."]*100, 3), format='f', digits=1 ), ")", sep="")
overall_c01s <- paste(formatC(round(out2_CR["vout", "X50."]*100, 3), format='f', digits=1 ), " (",
                      formatC(round(out2_CR["vout", "X2.5."]*100, 3), format='f', digits=1 ),", ", 
                      formatC(round(out2_CR["vout", "X97.5."]*100, 3), format='f', digits=1 ), ")", sep="")

ERCR <- merge(ER, CR, by="Author_Year", all.x = T)
write.table(ERCR, "output/ERCR13.txt", sep="\t")


########### Sec 4. Make Figure 1.
ERCR <- read.table("output/ERCR13.txt", sep="\t")
x1M <- ERCR[ , "p1s1M"]
x1L <- ERCR[ , "p1s1L"]
x1U <- ERCR[ , "p1s1U"]
y1M <- ERCR[ , "c10sM"]
y1L <- ERCR[ , "c10sL"]
y1U <- ERCR[ , "c10sU"]

x0M <- ERCR[ , "p0s1M"]
x0L <- ERCR[ , "p0s1L"]
x0U <- ERCR[ , "p0s1U"]
y0M <- ERCR[ , "c01sM"]
y0L <- ERCR[ , "c01sL"]
y0U <- ERCR[ , "c01sU"]


# plot in log scale x, dots in the middle
plot(x1U, y1M, pch = ".", col="red", log='x',
     xlim=c(0.01, 0.3),  ylim=c(0, 0.65), xaxt="n",yaxt="n", 
     xlab = "Event Rates (%)", ylab = "Noncompliance Rates (%)")
# labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
axis(1, at=c(0.01, 0.02, 0.05, 0.10, 0.20), labels=c(1, 2, 5, 10, 20))
axis(2, at=c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6), labels=c(0, 10, 20, 30, 40, 50, 60))
points(x1L, y1M, col="red", pch = ".")
points(x1M, y1L, col="red", pch = ".")
points(x1M, y1U, col="red", pch = ".")
points(x1M, y1M, col="red", pch=20)
segments(x1L, y1M, x1U, y1M, col="red", lty="dashed")
segments(x1M, y1L, x1M, y1U, col="red", lty="dashed")

points(x0U, y0M, pch = ".", col="blue")
points(x0L, y0M, pch = ".", col="blue")
points(x0M, y0L, pch = ".", col="blue")
points(x0M, y0U, pch = ".", col="blue")
points(x0M, y0M, col="blue", pch=18)
segments(x0L, y0M, x0U, y0M, col="blue", lty="solid")
segments(x0M, y0L, x0M, y0U, col="blue", lty="solid")
legend("topleft", legend = c("r=1", "r=0"), col = c("red", "blue"), pch = c(20,18), horiz = F, 
       lty = c("dashed","solid"), pt.cex = 1, cex = 1, text.col = "black", inset = c(0.1, 0.05))

