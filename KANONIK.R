#INPUT DATA# 
data=read.csv(file.choose(),header=T,sep=";") 
data 
attach(data) 
#PENGUJIAN ASUMSI# 
Y=cbind(Y1,Y2) 
X=cbind(X1,X2,X3,X4,X5) 
 
#LINEARITAS 
library(lmtest) 
resettest(Y~X) 
 
#NORMALITAS MULTIVARIAT#
y=data 
mu=colMeans(y) 
n=row(y)
p=ncol(y) 
cov=cov(y) 
d=sort(mahalanobis(y,mu,cov)) 
j=qchisq(ppoints(n),df=p) 
ks.test(d,j,df=p) 
 
#QQ-PLOT# 
qqplot(j,d,main="QQ-Plot Normal Multivariat",ylab="Jarak Mahalanlobis")  
abline(0,1) 
 
#NON MULTIKOLINEARITAS 
r=cor(y) 
diag(solve(r)) 

#ANALISIS KORELASI KANONIK# 
#ANALISIS KORELASI KANONIK#
kanonik=function(data,p,q,type=c("kovarians","korelasi"),alpha){
if (type=="korelasi"){
##MATRIX KORELASI
cat("KORELASI KANONIK DENGAN MATRIX KORELASI\n")
cat("\n")
n=nrow(data)
matrix.correlation=cor(data)
cat("MATRIX KORELASI\n")
print(matrix.correlation)
cat("\n")
YY=matrix.correlation[1:q,1:q]
YX=matrix.correlation[1:q,(q+1):(q+p)]
XX=matrix.correlation[(q+1):(q+p),(q+1):(q+p)]
XY=matrix.correlation[(q+1):(q+p),1:q]
k=min(p,q)
vektor.eigen.YY.U=eigen(YY)$vectors[1:q,1:k]
cat("VARIABEL KANONIK U (kolom)\n")
print(vektor.eigen.YY.U)
cat("\n")
vektor.eigen.XX.V=eigen(XX)$vectors[1:p,1:k]
cat("VARIABEL KANONIK V (kolom)\n")
print(vektor.eigen.XX.V)
cat("\n")
cor2=solve(YY)%*%YX%*%solve(XX)%*%XY
eigen(cor2)
eigen.values.cor2=eigen(cor2)$values
corellation.UV=rep(k,1)
for (i in 1:k){
corellation.UV[i]=sqrt(eigen.values.cor2[i])
}
cat("KORELASI ANTAR VARIABEL KANONIK\n")
print(data.frame(Korelasi.UV=corellation.UV))
cat("\n")
cat("KORELASI VARIABEL ZY DENGAN VARIABEL U \n")
eigen.values.YY=eigen(YY)$values
vektor.eigen.YY=eigen(YY)$vectors
YU=matrix(nrow=q,ncol=q)
for (i in 1:q){
for (j in 1:q){
YU[i,j]=(vektor.eigen.YY[j,i]*(sqrt(eigen.values.YY[i])))/(sqrt(YY[j,j])
)
}
}
print(t(YU[1:k,1:q]))
cat("\n")
cat("KORELASI VARIABEL ZX DENGAN VARIABEL V \n")
eigen.values.XX=eigen(XX)$values
vektor.eigen.XX=eigen(XX)$vectors
XV=matrix(nrow=p,ncol=p)
for (i in 1:p){
for (j in 1:p){
XV[i,j]=(vektor.eigen.XX[j,i]*(sqrt(eigen.values.XX[i])))/(sqrt(XX[j,j])
)
}
}
print(t(XV[1:k,1:p]))
cat("\n")
cat("UJI SIGNIFIKANSI\n")
cat("Hipotesis\n")
cat("H0 : Korelasi Kanonik yang dihasilkan TIDAK SIGNIFIKAN\n")
cat("H1 : Korelasi Kanonik yang dihasilkan SIGNIFIKAN\n")
cat("Taraf Signifikan = ")
print(alpha)
cat("\n")
A=(det(matrix.correlation))/(det(YY)*det(XX))
chisquare.hitung=-(n-(0.5*(p+q+2)))*log(A)
chisquare.tabel=qchisq((1-alpha),p+q)
k<- if (chisquare.hitung>=chisquare.tabel)"H0 DITOLAK" else"H0 
DITERIMA"
hasil=data.frame(chisquare.hitung,chisquare.tabel,Kesimpulan=k)
print(hasil)
}
else {
##MATRIX KOVARIANS
cat("ANALISIS KORELASI KANONIK DENGAN MATRIX 
KOVARIANS\n")
cat("\n")
n=nrow(data)
matrix.covarians=cov(data)
cat("Matrix Covarians\n")
print(matrix.covarians)
cat("\n")
YY=matrix.covarians[1:q,1:q]
YX=matrix.covarians[1:q,(q+1):(q+p)]
XX=matrix.covarians[(q+1):(q+p),(q+1):(q+p)]
XY=matrix.covarians[(q+1):(q+p),1:q]
k=min(p,q)
vektor.YY.U=eigen(YY)$vectors[1:q,1:k]
cat("VARIABEL KANONIK U (kolom)\n")
print(vektor.YY.U)
cat("\n")
vektor.XX.V=eigen(XX)$vectors[1:p,1:k]
cat("VARIABEL KANONIK V (kolom)\n")
print(vektor.XX.V)
cat("\n")
cor2=solve(YY)%*%YX%*%solve(XX)%*%XY
eigen(cor2)
eigen.values.cor2=eigen(cor2)$values
cor.UV=rep(k,1)
for (i in 1:k){
cor.UV[i]=sqrt(eigen.values.cor2[i])
}
cat("KORELASI ANTAR VARIABEL KANONIK\n")
print(data.frame(KORELASI.UV=cor.UV))
cat("\n")
cat("KORELASI VARIABEL Y DENGAN VARIABEL U \n")
eigen.values.YY=eigen(YY)$values
eigen.vektor.YY=eigen(YY)$vectors
YU=matrix(nrow=q,ncol=q)
for (i in 1:q){
for (j in 1:q){
YU[i,j]=(eigen.vektor.YY[j,i]*(sqrt(eigen.values.YY[i])))/(sqrt(YY[j,j])
)
}
}
print(t(YU[1:k,1:q]))
cat("\n")
cat("KORELASI VARIABEL X DENGAN VARIABEL V \n")
eigen.values.XX=eigen(XX)$values
eigen.vektor.XX=eigen(XX)$vectors
XV=matrix(nrow=p,ncol=p)
for (i in 1:p){
for (j in 1:p){
XV[i,j]=(eigen.vektor.XX[j,i]*(sqrt(eigen.values.XX[i])))/(sqrt(XX[j,j])
)
}
}
print(t(XV[1:k,1:p]))
cat("\n")
cat("UJI SIGNIFIKANSI KORELASI KANONIK \n")
cat("Hipotesis\n")
cat("H0 : Korelasi Kanonik yang dihasilkan TIDAK SIGNIFIKAN\n")
cat("H1 : Korelasi Kanonik yang dihasilkan SIGNIFIKAN\n")
cat("Taraf Signifikan = ")
print(alpha)
cat("\n")
A=(det(matrix.covarians))/(det(YY)*det(XX))
chisquare.hitung=-(n-(0.5*(p+q+2)))*log(A)
chisquare.tabel=qchisq((1-alpha),p+q)
k <- if (chisquare.hitung>=chisquare.tabel)"H0 DITOLAK" else"H0 
DITERIMA"
hasil=data.frame(chisquare.hitung,chisquare.tabel,Kesimpulan=k)
print(hasil)
}
}
kanonik(data,5,2,type="kovarians",0.05)
