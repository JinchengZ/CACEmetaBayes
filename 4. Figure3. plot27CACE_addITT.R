library(forestplot)
library(metafor)
setwd("submitted")
data <- read.table("data/data.txt", sep="", header = TRUE)
data$author_year <- paste(data$Author, data$Year, sep=", ")
data$author_year[1:10] <- paste(data$author_year[1:10], " ???", sep="")

data0 <- read.table("output/CACEi.txt", sep="", header = TRUE)
data1 <- cbind(data$author_year, data0)
colnames(data1)<- c("Author_year", "Lower", "Mean", "Median", "Upper", "SD", "HPD_l", "HPD_u")
data2 <- data1[order(data1$Author_year), ] 

result <- read.table("output/out2_5c_cacei.txt", sep="", header = TRUE)[c("CACE"), c("X2.5.", "X50.", "X97.5.")]

data_ITT <- read.table("data/data_ITT.txt", sep="", header = TRUE)
data_test <- data_ITT[data_ITT$miss==1, ]

ma_m <- rma(ai=n1s1, bi=n1s0, ci=n0s1, di=n0s0, data = data_ITT, measure = "RD", 
            method="FE", slab=study)
# ma_m <- rma(ai=n1s1, bi=n1s0, ci=n0s1, di=n0s0, data = data_temp, measure = "RR",
#             method="FE", slab=study)
# exp(ma_m$ci.lb)
# exp(ma_m$ci.ub)
# exp(ma_m$b)
# data_temp <- subset(data_ITT, data_ITT$Author_Year!="Nikkola_1997")


# add analysis for a single study, for example, Ramin 1995
data_R1995 <- subset(data_ITT, data_ITT$Author_Year=="Clark_1998")
ma_R1995 <- rma(ai=n1s1, bi=n1s0, ci=n0s1, di=n0s0, data = data_R1995, measure = "RD", 
            method="FE", slab=study)
ma_R1995


meandata <- data2$Mean
mediandata <- data2$Median
lowerdata <- data2$Lower
#lowerdata <- data2$HPD_l # HPD lower
upperdata <- data2$Upper
#upperdata <- data2$HPD_u # HPD upper

CACE_for_each <- 
  structure(list(
    median  = c(NA, NA, data2$Median, round(as.numeric(result[2]), 3), 
                round(as.numeric(ma_m$b), 3) ), 
    # mean  = c(NA, NA, data2$Mean, -0.026), 
    lower = c(NA, NA, lowerdata, round(as.numeric(result[1]), 3), 
              round(as.numeric(ma_m$ci.lb), 3) ), 
    upper = c(NA, NA, upperdata, round(as.numeric(result[3]), 3), 
              round(as.numeric(ma_m$ci.ub), 3) )),
    .Names = c("median", "lower", "upper"), 
    row.names = c(NA, -31L), 
    class = "data.frame")

tabletext<-cbind(
  c("", "Study (Author, Year)", paste(data2$Author_year), "Overall (Model Vc)",
    "Overall RD (ITT)"),
  
  c("", "CACE", paste(formatC(round(data2$Median, 3), format='f', digits=3 ), " (",
                      formatC(round(lowerdata, 3), format='f', digits=3 ),", ", 
                      formatC(round(upperdata, 3), format='f', digits=3 ), ")", sep="") , 
    paste(round(as.numeric(result[2]), 3), " (", 
          round(as.numeric(result[1]), 3), ", ", 
          round(as.numeric(result[3]), 3), ")", sep="" ) , 
    paste(round(ma_m$b, 3), " (", 
        round(ma_m$ci.lb, 3), ", ", 
        round(ma_m$ci.ub, 3), ")", sep="" ) )
)

own.f <- fpTxtGp(ticks = gpar(cex=0.85), xlab  = gpar(cex = 0.85))

xticks <- seq(from = -.3, to = 0.45, by = 0.1)
xtlab <- rep(c(TRUE, FALSE), length.out = length(xticks))
attr(xticks, "labels") <- xtlab

CACE_complete <- CACE_for_each
CACE_marginal <- CACE_for_each
CACE_marginal[c(3, 4, 10, 11, 14, 20, 21, 23, 24, 29),]<- NA
CACE_complete[-c(3, 4, 10, 11, 14, 20, 21, 23, 24, 29),]<- NA

forestplot(tabletext, 
           boxsize = .2, 
           hrzl_lines = gpar(lwd=1, col="#444444"),
           txt_gp = own.f,
           mean = cbind(CACE_complete[, "median"], CACE_marginal[, "median"]),
           lower = cbind(CACE_complete[, "lower"], CACE_marginal[, "lower"]),
           upper = cbind(CACE_complete[, "upper"], CACE_marginal[, "upper"]),
           new_page = TRUE, 
           is.summary=c(TRUE, TRUE, rep(FALSE,27), TRUE, TRUE, TRUE),
           clip=c(-1.0, 1.0), 
           lty.ci = c(1, 5),
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"), 
           xticks = xticks) 
