setwd("Your/working/directory")
load("GLbio_DR_pca.Rdata")
#install.packages("kernlab")
library(kernlab)
#install.packages("devtools")
library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)

pca <- prcomp(t(expression_mat))

plot(pca,type = "l")

#pdf("PCA.pdf")
g = ggbiplot(pca,choices = 1:2,var.axes = 0,groups = as.factor(rep(c(0,1,2,3,5,9,18),each= 17)))+ 
  theme_bw() + scale_color_brewer(palette="Reds") + labs(color='time(h)')  
#dev.off()
kpca = kpca(t(expression_mat),kernel="rbfdot",kpar = list(sigma = 0.00002))

#pdf("kpca.pdf")
plot(rotated(kpca),
     xlab="1st Principal Component",ylab="2nd Principal Component",
     col = rep(c("#FEE5D9","#FCBBA1","#FC9272","#FB6A4A","#EF3B2C","#CB181D","#99000D"),each= 17),
     pch=16)
#dev.off()
