#Diffusion map portion of this code adapted from Kye Taylor's eigenmap code
#https://www.mathworks.com/matlabcentral/fileexchange/36141-laplacian-eigenmap-diffusion-map-manifold-learning 
load('../Data/single_cell_processed_normal.RData')
library('Matrix')
library('pracma')
data=t(smat)
m=dim(data)[1]
knn=ceiling(0.03*m);
sigma2=100

#pairwise distance of cells based on scRNA-seq expression profiles
pdist=dist(data,method = "euclidean")
dt=as.matrix(pdist)

#KNN graph
srtdDt=sapply(split((dt), col((dt))), sort, index.return=TRUE);
sortedDtMat=matrix(0,nrow=m,ncol=m)
sortedDtIdx=matrix(0,nrow=m,ncol=m)
for (i in 1:m) {
sortedDtMat[,i]=t(srtdDt[1,i][[1]])
sortedDtIdx[,i]=t(srtdDt[2,i][[1]])
}
dt=sortedDtMat[1:(knn+1),]
nidx=sortedDtIdx[1:(knn+1),]

#compute weights - gaussian kernel 
tempW=exp(-dt^2/sigma2)
i=repmat(1:m,knn+1,1)
W = sparseMatrix(i=as.vector(i),j=as.vector(nidx),x=as.vector(tempW))
W2=pmax(W,t(W))
ld = diag(rowSums(W2)^(-1/2))
#symmetric normalized laplacian 
DO = ld %*% W2 %*% ld
DO = pmax(DO,t(DO))
ev <- eigen(DO)
k=5
v=ev$vectors[,1:k]
d=diag(ev$values[1:k])
eigVecIdx=rbind(c(2,3),c(2,4),c(3,4))

#make the lower dimensional embedding figures
for (i in 1:3) {
	jpeg(paste("lowDimSCrnaSeq",i,".jpg",sep=""), width = 600, height = 600)
	rbc <- rainbow(m)
#Rrb <- t(col2rgb(rbc))
	plot(v[,eigVecIdx[i,1]],v[,eigVecIdx[i,2]],col = rbc)
	dev.off()
}

#Continue with spectral clustering w/ k=5 (number of eigenvecs)
EV  = v;
for (j in 1:m){
        EV[j,] = v[j,]/Norm(v[j,]);
}
cidx = kmeans(EV,k)
cidx=t(cidx);
indx=sort(cidx[[1]],index.return=TRUE)
smat_sorted=data[indx[[2]],]

#Here I am organizing the cells by cluster and plotting the expression as a heatmap
library('devtools')
library('gplots')
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
jpeg(paste("heatmap2.jpg",sep=""), width = 600, height = 600)
colors=vector(mode="character", length=m)
for (i in 1:length(indx[[1]])){
if (indx[[1]][i]==1) {
colors[i]="darkred"
}
if (indx[[1]][i]==2) {
colors[i]="darkorchid"
}
if (indx[[1]][i]==3) {
colors[i]="green"
}
if (indx[[1]][i]==4) {
colors[i]="blue"
}
if (indx[[1]][i]==5) {
colors[i]="pink"
}
}
names(colors)="Cluster"
heatmap.3(as.matrix(smat_sorted),dendrogram='none', trace="none",RowSideColors=t(colors), labCol=FALSE, labRow=FALSE,breaks=seq(0,0.01,0.0001),col=redblue(100),symm=F,symkey=F,symbreaks=T,scale="none",Rowv=FALSE, Colv=TRUE)
dev.off()
