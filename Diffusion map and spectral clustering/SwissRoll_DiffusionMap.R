#Adapted from Kye Taylor's MATLAB eigenmap code
#https://www.mathworks.com/matlabcentral/fileexchange/36141-laplacian-eigenmap-diffusion-map-manifold-learning
#create the swiss roll data first 
N=2^11;
t=runif(N, min = 0, max = 1)
t=t(sort(4*pi*sqrt(t)))
h=runif(N, min = 0, max = 1)
z = 8*pi*h;
x = (t+0.1)*cos(t);
y = (t+0.1)*sin(t);
data = data.frame(x=t(x),y=t(y),z=z); 

#plot roll
library("scatterplot3d")
jpeg("swissroll.jpg", width = 600, height = 600)
rbc <- rainbow(N)
scatterplot3d(data$x, data$y, data$z,angle = 90,color = rbc)
dev.off()

#Changing these values will lead to different nonlinear embeddings
knn=ceiling(0.03*N);
sigma2=100
m=dim(data)[1]

#Distance matrix
library('pracma')
pdist=dist(data,method = "euclidean")
dt=as.matrix(pdist)

#KNN graph
srtdDt=sapply(split((dt), col((dt))), sort, index.return=TRUE)
sortedDtMat=matrix(0,nrow=N,ncol=N)
sortedDtIdx=matrix(0,nrow=N,ncol=N)
for (i in 1:N) {
sortedDtMat[,i]=t(srtdDt[1,i][[1]])
sortedDtIdx[,i]=t(srtdDt[2,i][[1]])
}
dt=sortedDtMat[1:(knn+1),]
nidx=sortedDtIdx[1:(knn+1),]

#compute weights - gaussian kernel 
tempW=exp(-dt^2/sigma2)
i=repmat(1:m,knn+1,1)
library(Matrix)
W = sparseMatrix(i=as.vector(i),j=as.vector(nidx),x=as.vector(tempW))
W2=pmax(W,t(W))
ld = diag(rowSums(W2)^(-1/2))
#symmetric normalized laplacian
DO = ld %*% W2 %*% ld
DO = pmax(DO,t(DO))
ev <- eigen(DO)
k=10
v=ev$vectors[,1:k]
d=diag(ev$values[1:k])
eigVecIdx=rbind(c(2,3),c(2,4),c(3,4))

#make the lower dimensional embedding figures from different eigenvectors
for (i in 1:3) {
	jpeg(paste("lowDim",i,".jpg",sep=""), width = 600, height = 600)
	rbc <- rainbow(N)
	plot(v[,eigVecIdx[i,1]],v[,eigVecIdx[i,2]],col = rbc)
	dev.off()
}
