args <- commandArgs(trailingOnly = TRUE)

mpt.data<-read.table(args[1],header = F,sep = ",")
summary(mpt.data$V2)
print(c("Standard deviation:", round(sd(mpt.data$V2),digits = 2)))