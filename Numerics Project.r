library(ggplot2)
library(lattice)
library(scales)

nu=0.1
T=1
deltat=0.00001
NT=T/deltat

N=1000 #gives a total of N+1 x grid points
L=pi
xgrid=0:N*(L/N)
deltax=L/N

initialphi=sin(2*xgrid)^2
B<-2*nu*deltat/deltax^2
C<-U*deltat/deltax

#L=10
#xgrid=0:N*(L/N)
#deltax=L/N
#initialphi=1/(1+(xgrid-5)^2)

pl <- plot(xgrid,initialphi,type="l")
pl
currentphi=initialphi
for (k in 1:NT) {
  newphi0=currentphi[1]+nu*deltat/deltax^2*(currentphi[2]-2*currentphi[1]+currentphi[N+1])-deltat/deltax*currentphi[1]*(currentphi[1]-currentphi[N+1])
  newphiN=currentphi[N+1]+nu*deltat/deltax^2*(currentphi[1]-2*currentphi[N+1]+currentphi[N])-deltat/deltax*currentphi[N+1]*(currentphi[N+1]-currentphi[N])
  newphi=c(newphi0,1:(N-1)*0,newphiN)
  for (j in 2:N){
    newphi[j]=currentphi[j]+nu*deltat/deltax^2*(currentphi[j+1]-2*currentphi[j]+currentphi[j-1])-deltat/deltax*currentphi[j]*(currentphi[j]-currentphi[j-1])
  }
  currentphi=newphi
  #plot(xgrid,newphi)
  if (k%%floor(NT/10)==0){
    lines(xgrid,newphi,type="l",col="blue")
    }
}
#stability analysis
U<-1
modA2<-function(b,c,kdx){
  return ((1+(b-c)*(cos(kdx)-1))^2+(c*sin(kdx))^2)
}
maxmodA2<-function(b,c){
  return (max(modA2(b,c,1:313/100)))
}
BN<-800
CN<-400
division<-400
bvec<-0:BN/division
cvec<-0:CN/division
heatmat<-matrix(0,BN+1,CN+1)
for (i in 1:(BN+1)){
  for (j in 1:(CN+1)){
    heatmat[i,j]<-maxmodA2(bvec[i],cvec[j])
  }
}
#levelplot(heatmat,contour=TRUE,labels=T,cuts=8)
colnames(heatmat)<-cvec
rownames(heatmat)<-bvec
contourplot(heatmat,region=TRUE,xlab="b",ylab="c",pretty=T,cuts=20)
contour(bvec,cvec,heatmat,levels=c(0,1,2,3,4,5,6,7,8),xlab="b",ylab="c")
#exact solution


