rm(list=ls(all=TRUE))

myfunction <- function(DATA){

d<-c(1:nrow(DATA))  

# Check for validity, 1 - ok, 0 - bad
which(DATA[,"validity"]<=0,arr.ind=TRUE) 

# Check for low intensities (<100 cnts), should not be interpreted
alpha<-which( colnames(DATA)=="Mg" )
omega<-which( colnames(DATA)=="Pb" )
x<-rep(NA,omega-alpha)
for (i in alpha:omega){
  x[i+1-alpha]<-mean(DATA[,i])
}
y<-which(x<=100,arr.ind=TRUE)
lowint<-names(DATA[y+alpha])

# Mg, P, Ar, Al, Mo, W cannot be measured correctly
unmeas<-c("Mg","P","Ar","Al","Mo","W")

# S, Cl, Br, Ba, Se, Ga and Y are mostly close to background level
backlvl<-c("S","Cl","Br","Ba","Se","Ga","Y")

# Eliminate DATA
no.DATA<-DATA[,!(names(DATA) %in% backlvl)]
no.DATA<-no.DATA[,!(names(no.DATA) %in% unmeas)]
coh<-no.DATA[,ncol(no.DATA)]
S<-DATA[,"S"]/coh
CaFe<-(DATA[,"Ca"]/coh)/(DATA["Fe"]/coh)
SrCa<-(DATA[,"Sr"]/coh)/(DATA["Ca"]/coh)
MnTi<-(DATA[,"Mn"]/coh)/(DATA["Ti"]/coh)
BaTi<-(DATA[,"Ba"]/coh)/(DATA["Ti"]/coh)

#no.DATA<-no.DATA[,!(names(no.DATA) %in% lowint)]
names(no.DATA)
mynames <- c("Si", "K","Ca","Sc","Ti", "V","Cr", "Mn", "Fe", "Ni", "Cu", "Zn", "Ge", "Kr","Rb", 
             "Sr", "Zr", "Cd", "Sn", "Gd","Pb","S","CaFe","SrCa","MnTi","BaTi")

# Normalization against Mo coh in order to account for matrix differences
MSE<-which( colnames(no.DATA)=="MSE" )
Fe.a.a<-which( colnames(no.DATA)=="Fe.a.a" )
D<-no.DATA[,(MSE+1):(Fe.a.a-1)]/no.DATA[,"Mo.coh"]
D<-cbind(D, S,CaFe,SrCa,MnTi,BaTi)

names(D) <- mynames
names(D)

which(D==0,arr.ind=TRUE)
which(D=="Inf",arr.ind=TRUE)

D<-log(D+100/mean(coh))

###############
# compute PCA #
# to z-standardize data (equal weight for all variables, note that the variables are dimensionally different)
zD<-scale(D,center=TRUE,scale=TRUE) # scale is command for standardization by standard deviation

pca<-prcomp(zD,retx=T,center=F,scale.=F) # retx - factor holdings TRUE (scores), scale.=F bc. already done this

pca$sdev # the standard deviations of the PCA axes (their squares are the eigenvalues)

pca$sdev^2
sum(pca$sdev^2)
pca$sdev^2/sum(pca$sdev^2)

pca$x # the site scores on all PCA-axes (the coordinates of observations in the ordination space defined by the PCs)
apply(pca$x,2,sd) # once again their standard deviations
scores<-pca$x

# the variable loadings (the coefficients of the linear combinations defining the PCs)
# here each column represents one PC as an additive linear combination of the coefficients given for each underlying variable 
pca$rotation # factor loadings, can be used for predictions in fellow-up study
head(scores)
zD %*% pca$rotation # to manually compute scores from variables and loadings, MATRIX COMPUTATION
loadings<-pca$rotation

summary(pca)

return(list(pca,d,D))
}