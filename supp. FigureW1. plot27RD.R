library(forestplot)
library(metafor)
setwd("submitted")
data <- read.table("data/data.txt", sep="", header = TRUE)
data$author_year <- paste(data$Author, data$Year, sep=", ")
data$author_year[1:10] <- paste(data$author_year[1:10], " *", sep="")

data0 <- read.table("output/CACEi.txt", sep="", header = TRUE)
data1 <- cbind(data$author_year, data0)
colnames(data1)<- c("Author_year", "Lower", "Mean", "Median", "Upper", "SD", "HPD_l", "HPD_u")
data2 <- data1[order(data1$Author_year), ] 

result <- read.table("output/out2_5c_cacei.txt", sep="", header = TRUE)[c("CACE"), c("X2.5.", "X50.", "X97.5.")]

data_ITT <- read.table("data/data_ITT.txt", sep="", header = TRUE)
data_test <- data_ITT[data_ITT$miss==1, ]

ma_m <- rma(ai=n1s1, bi=n1s0, ci=n0s1, di=n0s0, data = data_ITT, measure = "RD", 
            method="FE", slab=Author_Year, digits=3)
# # built in forest plot function:
# forest(ma_m, header=TRUE, digits=3)
ma_m$yi
ci.lb <- ma_m$yi - qnorm(0.025, lower.tail=FALSE) * sqrt(ma_m$vi)
ci.ub <- ma_m$yi + qnorm(0.025, lower.tail=FALSE) * sqrt(ma_m$vi)
meandata <- ma_m$yi[1:27]
mediandata <- ma_m$yi[1:27]
lowerdata <- ci.lb[1:27]
upperdata <- ci.ub[1:27]

RD_for_each <- 
  structure(list(
    mean  = c(NA, NA, meandata, round(as.numeric(ma_m$b), 3)), 
    lower = c(NA, NA, lowerdata, round(as.numeric(ma_m$ci.lb), 3) ), 
    upper = c(NA, NA, upperdata, round(as.numeric(ma_m$ci.ub), 3) )),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -30L), 
    class = "data.frame")

tabletext<-cbind(
  c("", "Study (Author, Year)", paste(data2$Author_year), "Overall"),
  
  c("", "Estimated", paste(formatC(round(meandata, 3), format='f', digits=3 ), " (",
                      formatC(round(lowerdata, 3), format='f', digits=3 ),", ", 
                      formatC(round(upperdata, 3), format='f', digits=3 ), ")", sep="") , 
    # paste(round(as.numeric(result[2]), 3), " (", 
    #       round(as.numeric(result[1]), 3), ", ", 
    #       round(as.numeric(result[3]), 3), ")", sep="" ) , 
    paste(round(ma_m$b, 3), " (", 
        round(ma_m$ci.lb, 3), ", ", 
        round(ma_m$ci.ub, 3), ")", sep="" ) )
)

own.f <- fpTxtGp(ticks = gpar(cex=0.85), xlab  = gpar(cex = 0.85))

xticks <- seq(from = -.3, to = 0.45, by = 0.1)
xtlab <- rep(c(TRUE, FALSE), length.out = length(xticks))
attr(xticks, "labels") <- xtlab

# CACE_complete <- CACE_for_each
# CACE_marginal <- CACE_for_each
# CACE_marginal[c(3, 4, 10, 11, 14, 20, 21, 23, 24, 29),]<- NA
# CACE_complete[-c(3, 4, 10, 11, 14, 20, 21, 23, 24, 29),]<- NA

CACE_for_each <- 
  structure(list(
    median  = c(NA, NA, data2$Median, round(as.numeric(result[2]), 3) ), 
    lower = c(NA, NA, data2$Lower, round(as.numeric(result[1]), 3)), 
    upper = c(NA, NA, data2$Upper, round(as.numeric(result[3]), 3))),
    .Names = c("median", "lower", "upper"), 
    row.names = c(NA, -30L), 
    class = "data.frame")


forestplot(tabletext[,1], 
           boxsize = .2, 
           legend_args = fpLegend(pos = list(x = .9, y = 0.45), 
                                  gp = gpar(col = "#CCCCCC", fill = "#F9F9F9")),
           legend = c("CACE", "RD (ITT)"),
           fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
           hrzl_lines = gpar(lwd=1, col="#444444"),
           txt_gp = own.f,
           mean = cbind(RD_for_each[, "mean"], CACE_for_each[, "median"]), 
           lower = cbind(RD_for_each[, "lower"], CACE_for_each[, "lower"]),
           upper = cbind(RD_for_each[, "upper"], CACE_for_each[, "upper"]),
           new_page = TRUE, 
           is.summary=c(TRUE, TRUE, rep(FALSE,27), TRUE, TRUE),
           clip=c(-1.0, 1.0), 
           # lty.ci = c(1, 5),
           col=fpColors(box=c("royalblue", "darkred"),
                        line=c("darkblue", "darkred"), 
                        summary=c("royalblue", "darkred") ), 
           xticks = xticks) 
