library(gmodels)
library(ggplot2)

group <- c("NF10_A","NF10_B","NF18_A","NF18_B","NF26_A","NF26_B","NF34_A","NF34_B","NR10_A","NR10_B","NR18_A","NR18_B","NR26_A","NR26_B","NR34_A","NR34_B","PF10_A","PF10_B","PF18_A","PF18_B","PF26_A","PF26_B","PF34_A","PF34_B","PF42_A","PF42_B","PF50_A","PF50_B")

expr <- read.table("miRNA_tpm.txt", header=T, row.names=1)
a <- expr[,-1]
data <- t(as.matrix(a))

data.pca <- fast.prcomp(data,retx=T,scale=F,center=T)

a <- summary(data.pca)
tmp <- a[4]$importance
pro1 <- as.numeric(sprintf("%.3f",tmp[2,1]))*100
pro2 <- as.numeric(sprintf("%.3f",tmp[2,2]))*100

pc = as.data.frame(a$x)
pc$group = group
pc$names = rownames(pc)

xlab=paste("PC1(",pro1,"%)",sep="") 
ylab=paste("PC2(",pro2,"%)",sep="")
pca=ggplot(pc,aes(PC1,PC2)) + geom_point(size=5,aes(color=group)) + geom_text(aes(label=names),size=2)+labs(x=xlab,y=ylab,title="PCA",size=5) + geom_hline(yintercept=0,linetype=4,color="grey") + geom_vline(xintercept=0,linetype=4,color="grey") + theme_bw()
ggsave("PCA_5tpm_miRNA.png",pca,width=10,height=8)