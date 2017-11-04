#Gene dependency data was downloaded from https://portals.broadinstitute.org/achilles/datasets/12/download named "Gene solutions"
# then breast cancer cell lines data were extracted, data were splitted into two groups, one group contains 11392 genes' dependency in 32 breast cancer cell lines,
# another group includes 5707 genes' dependency in 22 breast cancer cell lines.
# Pearson correlation implemented in cor function was utilized to calculate correlation coefficiency.

#aim: to calculate TEAD2 dependency-correlated gene
#loading packages
library("Hmisc")
library("corrplot")
library("pheatmap")
library("scatterplot3d")
library("ggplot2")
#loading input data & make up new data frame
data=read.delim("~/Desktop/demo.txt",head=T,sep="\t")
data2=data[,2:ncol(data)]
rownames(data2)=data[,1:1]
#calculate pangene-correlation
res=cor(data2)
result=res[,1:1]
write.csv(result,"~/Desktop/not_intact_2_result.csv")

#calculate single P value & correlation
x=data$SND1
y=data$TEAD2_Dependency
plot(y,x,xlim=c(-5,5),ylim=c(-5,5))
results2=lm(data$SND~data$TEAD2_Dependency,data=data)
summary(results2)

#plot distribution of correlation value
#intact --- correlated gene (11392 genes have dependency score in 34 breast cancer cell lines)
data=read.delim("~/Desktop/intact_correlate.txt",head=T,sep="\t")
value=data$Correlate
hist(value,breaks=20,freq=FALSE,xlim=c(-1,1),ylim=c(0,2),main="correlation distribution of 11392 genes in 34 breast cancer cell lines") #plot histogram 
sd=sd(value)
mean=mean(value)
x=seq(min(value),max(value),by=0.001)
y=dnorm(x,mean,sd)
lines(x,y,col="green",lwd=2)
abline(v=0,col="blue",lwd=2)
#not intact ----- correlated gene (5707 genes have dependency score in 21 breast cancer cell lines)
data2=read.delim("~/Desktop/not_intact_correlate.txt",head=T,sep="\t")
value=data2$Correlate
hist(value,breaks=20,freq=FALSE,xlim=c(-1,1),ylim=c(0,2),main="correlation distribution of 5707 genes in 21 breast cancer cell lines") #plot histogram 
sd=sd(value)
mean=mean(value)
x=seq(min(value),max(value),by=0.001)
y=dnorm(x,mean,sd)
lines(x,y,col="green",lwd=2)
abline(v=0,col="blue",lwd=2)




