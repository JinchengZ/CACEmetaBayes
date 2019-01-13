data <- read.csv("data/epidural_orig.csv", header=T)

#### Sec 1.  Create data_plot13 and obs13_final for noncompliance analysis (Figure 1)
data$Author_Year <- paste(data$Author, data$Year, sep = "_")
data$study <- 1:nrow(data) 

data$n110 <- data$n11s-data$n111
data$n100 <- data$n10s-data$n101
data$n010 <- data$n01s-data$n011
data$n000 <- data$n00s-data$n001
data$n1s0 <- data$n1ss-data$n1s1
data$n0s0 <- data$n0ss-data$n0s1
data$n1s1 # double check, because if using attach data, this variable is masked 

data2 <- data
x <- cbind(data2$n111, data2$n110, data2$n101, data2$n100)
data2$n1ss <- rowSums(x)
y <- cbind(data2$n011, data2$n010, data2$n001, data2$n000)
data2$n0ss <- rowSums(y)

data3 <- data2
for(i in 1:length(data2$n0ss)) { 
  if (is.na(data2$n1ss[i])) {
    data3$n111[i] <-data3$n110[i] <- data3$n101[i] <- data3$n100[i] <- 0
    data3$n10s[i] <- data3$n1ss[i] <- NA
    data3$miss_r1[i] <- 1
  }
  else {
    data3$n1s0[i] <- data3$n1s1[i] <- 0
    data3$miss_r1[i] <- 0
  }
  if (is.na(data2$n0ss[i])) {
    data3$n011[i] <- data3$n010[i] <- data3$n001[i] <- data3$n000[i] <- 0
    data3$n01s[i] <- data3$n0ss[i] <- NA
    data3$miss_r0[i] <- 1
  }
  else {
    data3$n0s0[i] <- data3$n0s1[i] <- 0
    data3$miss_r0[i] <- 0
  }
  if (data3$miss_r1[i] == 0 & data3$miss_r0[i] == 0) 
    data3$miss[i] <- 0  
  else data3$miss[i] <- 1  
}


data_plot <- data3[ which(data3$miss_r0==0|data3$miss_r1==0), ] # 13 studies
write.table(data_plot, "data/data_plot13.txt", sep="\t")


obs13 <- data_plot[c("Author_Year", "n01s", "n0ss", "n10s", "n1ss", "miss_r0", "miss_r1", "study")]
obs13$id <- 1:nrow(obs13) 

obs_r0 <- obs13[c("id", "Author_Year", "n01s", "n0ss", "miss_r0", "study")]
names(obs_r0) <- c("id", "Author_Year", "y", "n", "mis", "study")
obs_r0$r <- 0
obs_r1 <- obs13[c("id", "Author_Year", "n10s", "n1ss", "miss_r1", "study")]
names(obs_r1) <- c("id", "Author_Year", "y", "n", "mis", "study")
obs_r1$r <- 1
obs13_long <- rbind(obs_r0, obs_r1)
obs13_long <- obs13_long[order(obs13_long[["id"]], obs13_long[["r"]]),]

obs13_final = obs13_long[obs13_long$mis==0,]
obs13_final$rr <- obs13_final$r + 1
write.table(obs13_final, "data/obs13_final.txt", sep="\t")



#### Sec 2. Create data_ITT and data.txt from Figure 2, 3, 4, analysis. 
data_final <- subset(data3, select=-c(n1ss, n0ss))
sum(data_final$miss) # there are 17 studies with missing data

data_ITT <- cbind(subset(data, select=c(n1s1, n1s0, n0s1, n0s0, n1ss, n0ss)), 
                  data_final$Author_Year, data_final$study, data_final$miss)
colnames(data_ITT) <- c("n1s1", "n1s0", "n0s1", "n0s0", "n1ss",	"n0ss",	
                        "Author_Year",	"study",	"miss")
write.table(data_ITT, "data/data_ITT.txt", sep="\t")

newdata <- data_final[order(data_final$miss_r0, data_final$miss_r1),]
write.table(newdata, "data/data.txt", sep="\t")
