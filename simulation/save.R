R <- 2000

setwd('simout')

CACE0_0 <- colMeans(cbind(read.table("CACE0_0.txt", sep="\t", header=T), 
                          read.table("CACE0_0c.txt", sep="\t", header=T) ) )
CACE0_1 <- colMeans(cbind(read.table("CACE0_1.txt", sep="\t", header=T), 
                          read.table("CACE0_1c.txt", sep="\t", header=T) ) )
CACE0_2 <- colMeans(cbind(read.table("CACE0_2.txt", sep="\t", header=T), 
                          read.table("CACE0_2c.txt", sep="\t", header=T) ) )
CACE0_3 <- colMeans(cbind(read.table("CACE0_3.txt", sep="\t", header=T),
                          read.table("CACE0_3c.txt", sep="\t", header=T) ) )
CACE1_1 <- colMeans(cbind(read.table("CACE1_1.txt", sep="\t", header=T), 
                          read.table("CACE1_1c.txt", sep="\t", header=T) ) )
CACE1_2 <- colMeans(cbind(read.table("CACE1_2.txt", sep="\t", header=T), 
                          read.table("CACE1_2c.txt", sep="\t", header=T) ) )
CACE1_3 <- colMeans(cbind(read.table("CACE1_3.txt", sep="\t", header=T), 
                          read.table("CACE1_3c.txt", sep="\t", header=T) ) )
CACE4_4_2 <- colMeans(cbind(read.table("CACE4_4_2.txt", sep="\t", header=T), 
                            read.table("CACE4_4_2c.txt", sep="\t", header=T) ) )
a <- rbind(CACE0_0, CACE0_1, CACE0_2, CACE0_3, 
           CACE1_1, CACE1_2, CACE1_3, 
           CACE4_4_2)
write.table(a, "results.txt", sep="\t")




b <- rbind(read.table("table0_0.txt", sep="\t", header=T),
           read.table("table0_1.txt", sep="\t", header=T),
           read.table("table0_2.txt", sep="\t", header=T),
           read.table("table0_3.txt", sep="\t", header=T),
           read.table("table1_1.txt", sep="\t", header=T),
           read.table("table1_2.txt", sep="\t", header=T),
           read.table("table1_3.txt", sep="\t", header=T),
           read.table("table4_4_2.txt", sep="\t", header=T))
rownames(b) <- c("select0_0", "select0_1", "select0_2", "select0_3", 
                "select1_1", "select1_2", "select1_3",  
                "select4_4_2")
write.table(b, "selectresults.txt", sep="\t")




cc <- cbind(read.table("selectCACE0_0.txt", sep="\t", header=T)$x,
           read.table("selectCACE0_1.txt", sep="\t", header=T)$x,
           read.table("selectCACE0_2.txt", sep="\t", header=T)$x,
           read.table("selectCACE0_3.txt", sep="\t", header=T)$x,
           read.table("selectCACE1_1.txt", sep="\t", header=T)$x,
           read.table("selectCACE1_2.txt", sep="\t", header=T)$x,
           read.table("selectCACE1_3.txt", sep="\t", header=T)$x,
           read.table("selectCACE4_4_2.txt", sep="\t", header=T)$x)
rownames(cc) <- 
  c("True", "Mean", "Q2.5", "Q50", "Q97.5", "Coverage", "CILength", "Rel_Bias", "Bias", "Bias_sq", 
    "Mean_1", "Q2.5_1", "Q50_1", "Q97.5_1", "Coverage_1", "CILength_1", "Rel_Bias_1", "Bias_1", "Bias_sq_1", 
    "Mean_2", "Q2.5_2", "Q50_2", "Q97.5_2", "Coverage_2", "CILength_2", "Rel_Bias_2", "Bias_2", "Bias_sq_2", 
    "Mean_3", "Q2.5_3", "Q50_3", "Q97.5v", "Coverage_3", "CILength_3", "Rel_Bias_3", "Bias_3", "Bias_sq_3")

table_all <- matrix(NA, 32, 4)

for (i in 1:8) {
  table_all[((i-1)*4)+1, ] <- c(cc["Rel_Bias", i], cc["Rel_Bias_1", i], cc["Rel_Bias_2", i], cc["Rel_Bias_3", i])
  table_all[((i-1)*4)+2, ] <- c(cc["Bias_sq", i], cc["Bias_sq_1", i], cc["Bias_sq_2", i], cc["Bias_sq_3", i])
  table_all[((i-1)*4)+3, ] <- c(cc["Coverage", i], cc["Coverage_1", i], cc["Coverage_2", i], cc["Coverage_3", i])
  table_all[((i-1)*4)+4, ] <- c(cc["CILength", i], cc["CILength_1", i], cc["CILength_2", i], cc["CILength_3", i])
}

rownames(table_all) <- c("Rel_bias_0_0", "MSE_0_0", "CP_0_0", "CIL_0_0", 
                         "Rel_bias_0_1", "MSE_0_1", "CP_0_1", "CIL_0_1", 
                         "Rel_bias_0_2", "MSE_0_2", "CP_0_2", "CIL_0_2", 
                         "Rel_bias_0_3", "MSE_0_3", "CP_0_3", "CIL_0_3", 
                         "Rel_bias_1_1", "MSE_1_1", "CP_1_1", "CIL_1_1", 
                         "Rel_bias_1_2", "MSE_1_2", "CP_1_2", "CIL_1_2", 
                         "Rel_bias_1_3", "MSE_1_3", "CP_1_3", "CIL_1_3",  
                         "Rel_bias_4_4_2", "MSE_4_4_2", "CP_4_4_2", "CIL_4_4_2")
colnames(table_all) <- c("Model_0", "Model_1", "Model_2", "Model_3")

write.table(table_all, "table_all.txt", sep="\t")
