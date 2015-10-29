################# General information
# TEG source: microRFA scans by GEOPOLAR Bremen done for four sediment cores
# sediment cores were taken at Lake Tegel on 30.09.15
# TEG-1 deepest lake site, 2 near PEP, 3 Tegelort, 4 outer range of main basin

rm(list=ls(all=TRUE))

# loading of additional packages for plotting
library(shape) 
library(vegan)
library(ggplot2)

# loading of data sets for each sediment core
setwd("D:/work/data/Hupfer") # work directory, change to own one if it differs from default
source("D:/work/data/Hupfer/PCA/github/SubLakeTegel_PCA.R", echo=T)

TEG1<-read.table(file="PCA/TEG1_1mm.txt",header=T,sep="",skip=2) 
TEG2<-read.table(file="PCA/TEG2_1mm.txt",header=T,sep="",skip=2) 
TEG3<-read.table(file="PCA/TEG3_1mm.txt",header=T,sep="",skip=2) 
TEG4<-read.table(file="PCA/TEG4_1mm.txt",header=T,sep="",skip=2) 

# calls sub function and calculates pca for all cores
a<-myfunction(TEG1)
pca1<-a[[1]]
scores1<-pca1$x
loadings1<-pca1$rotation
d1<-a[[2]]
T1<-a[[3]]

a<-myfunction(TEG2)
pca2<-a[[1]]
scores2<-pca2$x
loadings2<-pca2$rotation
d2<-a[[2]]
T2<-a[[3]]

a<-myfunction(TEG3)
pca3<-a[[1]]
scores3<-pca3$x
loadings3<-pca3$rotation
d3<-a[[2]]
T3<-a[[3]]

a<-myfunction(TEG4)
pca4<-a[[1]]
scores4<-pca4$x
loadings4<-pca4$rotation
d4<-a[[2]]
T4<-a[[3]]

# for plotting as color, creating color ramp to have a transition from green (upper part) to red (lower part)
my.colorRamp.fct<-colorRamp(c("green","yellow","red")) # custom-built color ramp
my.colorRamp.fct(0.5)
rgb(my.colorRamp.fct(0.5),maxColorValue=255)

###############
# Check importance of PCs

layout(matrix(1:4,2,2))

pca1.pct<-100*round(summary(pca1)$importance[2,],3)          
barplot(pca1.pct, main="Deepest Site")
pca2.pct<-100*round(summary(pca2)$importance[2,],3)          
barplot(pca2.pct, main="P Elimination Plant")
pca3.pct<-100*round(summary(pca3)$importance[2,],3)          
barplot(pca3.pct, main="River")
pca4.pct<-100*round(summary(pca4)$importance[2,],3)          
barplot(pca4.pct, main = "Main Basin")
# impression that only for TEG1 both PC1 and PC2 are important

###############
# PCA biplots #
# for a DISTANCE BIPLOT (focus is on sites, "scaling 1")
# each principal component has variance given by eigenvalue, loadings remain unscaled

# in this plot:
# 1) distances among sites are approximating true Euclidean distances in multivariate space
# 2) angles between arrows do not reflect correlations among variables
# 3) projecting site on descriptor at right angle gives its appr. descriptor value

#layout(matrix(c(1,2,3,4,5,5), nrow = 3, ncol = 2, byrow=TRUE), heights = c(0.4, 0.4, 0.2))
layout(matrix(1:4,2,2))
xsi<-c(-10,10)
ysi<-c(-5,5)
ext<-7
arleng<-0.15

as.rgb.channels<-my.colorRamp.fct(decostand(d1,method="range"))
peak.colors<-rgb(as.rgb.channels,maxColorValue=255)
plot(scores1[,1]*(-1), scores1[,2]*(-1),asp=1,pch=21,bg=peak.colors, xlim = xsi,ylim=ysi, main="Deepest Site", xlab = "PC1", ylab ="PC2") # asp=1 x and y are equally scaled
arrows<- -loadings1*ext # with extension factor to get nice graphs
plotcircle(r=7*sqrt(2/ncol(T1)),lcol="blue") # circle of equilibrium contribution (equal contribution to all PCA-dimensions)
Arrows(x0=0,y0=0,x1=arrows[,1],y1=arrows[,2],col="darkgreen",arr.length=arleng)
# circle shows equilibrium, arrows that are longer than circle are more important for PCA
text(x=arrows[,1]*1.3,y=arrows[,2]*1.2,labels=names(T1),cex=1.5)
legend("topleft",c("Top","Middle","Bottom"),lty=c(1,1),
      lwd=c(2.5,2.5), col=c("green","yellow","red"))

as.rgb.channels<-my.colorRamp.fct(decostand(d2,method="range"))
peak.colors<-rgb(as.rgb.channels,maxColorValue=255)
plot(scores2[,1:2],asp=1,pch=21,bg=peak.colors, xlim = xsi,ylim=ysi, main="P Elimination Plant")
arrows<-loadings2*ext 
plotcircle(r=7*sqrt(2/ncol(T2)),lcol="blue") 
Arrows(x0=0,y0=0,x1=arrows[,1],y1=arrows[,2],col="darkgreen",arr.length=arleng)
text(x=arrows[,1]*1.3,y=arrows[,2]*1.2,labels=names(T2),cex=1.5)

as.rgb.channels<-my.colorRamp.fct(decostand(d3,method="range"))
peak.colors<-rgb(as.rgb.channels,maxColorValue=255)
plot(scores3[,1]*(-1), scores3[,2]*(-1),asp=1,pch=21,bg=peak.colors, xlim = xsi,ylim=ysi, main="River", xlab = "PC1", ylab ="PC2") 
arrows<- -loadings3*ext 
plotcircle(r=7*sqrt(2/ncol(T3)),lcol="blue") 
Arrows(x0=0,y0=0,x1=arrows[,1],y1=arrows[,2],col="darkgreen",arr.length=arleng)
text(x=arrows[,1]*1.3,y=arrows[,2]*1.2,labels=names(T1),cex=1.5)

as.rgb.channels<-my.colorRamp.fct(decostand(d4,method="range"))
peak.colors<-rgb(as.rgb.channels,maxColorValue=255)
plot(scores4[,1:2],asp=1,pch=21,bg=peak.colors, xlim = xsi,ylim=ysi, main = "Main Basin") 
arrows<-loadings4*ext 
plotcircle(r=7*sqrt(2/ncol(T4)),lcol="blue")
Arrows(x0=0,y0=0,x1=arrows[,1],y1=arrows[,2],col="darkgreen",arr.length=arleng)
text(x=arrows[,1]*1.3,y=arrows[,2]*1.2,labels=names(T1),cex=1.5)

# plot(1, type = "n", axes=FALSE, xlab="", ylab="")
# legend("center", pch=21, pt.bg = c("green", "yellow", "red"), legend=c("Top", "Mid", "Bottom"))

###############
# CORRELATION BIPLOT
# arrows which reach out of the circle contribute more than on average
# biplot(pca1,scale=0)
# for a CORRELATION BIPLOT (focus is on variables, "scaling 2") # also for paper: CORRELATION/DISTANCE BIPLOT
# each principal component is weighted by 1/sqrt(eigenvalue), so it has variance 1

# in this plot
# 1) distances among sites are not approximating true Euclidean distances in multivariate space
# 2) angles between arrows reflect correlations among variables (NOT proximity of arrow heads)
# 3) projecting site on descriptor at right angle gives its appr. descriptor value

layout(matrix(1:4,2,2))

var(scores1[,1]/pca1$sdev[1]) # just demo
plot(scores1[,1]/pca1$sdev[1],scores1[,2]/pca1$sdev[2],pch=21,bg=peak.colors,asp=1, main="Deepest Site")
# loadings are weighted by sqrt(eigenvalues) (multiplied by sqrt(eigenvalues))
arrows<-loadings1*matrix(pca1$sdev,nrow=nrow(loadings1),ncol=ncol(loadings1),byrow=TRUE) # scaling by STD
arrows<-arrows*2 # choose extension factor
Arrows(x0=0,y0=0,x1=arrows[,1],y1=arrows[,2],col="purple")
# as alternative just compute correlation of scores with original data ("structure coefficients")
(structure<-cor(T1,scores1))
structure<-2*structure
Arrows(x0=0,y0=0,x1=structure[,1],y1=structure[,2],col="red")
text(x=arrows[,1]*1.3,y=arrows[,2]*1.2,labels=names(T1),cex=1.5)

var(scores2[,1]/pca2$sdev[1]) 
plot(scores2[,1]/pca2$sdev[1],scores2[,2]/pca2$sdev[2],pch=21,bg=peak.colors,asp=1, main="P Elimination Plant")
arrows<-loadings2*matrix(pca2$sdev,nrow=nrow(loadings2),ncol=ncol(loadings2),byrow=TRUE) 
arrows<-arrows*2 
Arrows(x0=0,y0=0,x1=arrows[,1],y1=arrows[,2],col="purple")
(structure<-cor(T2,scores2))
structure<-2*structure
Arrows(x0=0,y0=0,x1=structure[,1],y1=structure[,2],col="red")
text(x=arrows[,1]*1.3,y=arrows[,2]*1.2,labels=names(T2),cex=1.5)

var(scores3[,1]/pca3$sdev[1]) 
plot(scores3[,1]/pca3$sdev[1],scores3[,2]/pca3$sdev[2],pch=21,bg=peak.colors,asp=1, main="River")
arrows<-loadings3*matrix(pca3$sdev,nrow=nrow(loadings3),ncol=ncol(loadings3),byrow=TRUE) 
arrows<-arrows*2 
Arrows(x0=0,y0=0,x1=arrows[,1],y1=arrows[,2],col="purple")
(structure<-cor(T3,scores3))
structure<-2*structure
Arrows(x0=0,y0=0,x1=structure[,1],y1=structure[,2],col="red")
text(x=arrows[,1]*1.3,y=arrows[,2]*1.2,labels=names(T3),cex=1.5)

var(scores4[,1]/pca4$sdev[1])
plot(scores4[,1]/pca4$sdev[1],scores4[,2]/pca4$sdev[2],pch=21,bg=peak.colors,asp=1, main="Main Basin")
arrows<-loadings4*matrix(pca4$sdev,nrow=nrow(loadings4),ncol=ncol(loadings4),byrow=TRUE) 
arrows<-arrows*2 
Arrows(x0=0,y0=0,x1=arrows[,1],y1=arrows[,2],col="purple")
(structure<-cor(T4,scores4))
structure<-2*structure
Arrows(x0=0,y0=0,x1=structure[,1],y1=structure[,2],col="red")
text(x=arrows[,1]*1.3,y=arrows[,2]*1.2,labels=names(T4),cex=1.5)

#biplot(pca1,scale=1) # alternative command

###############
# Comparison of PC1 over the depth for all cores
layout(matrix(1:4,2,2))
plot(scores1[,1],d1,ylim = rev(c(0,length(d1))), main="Deepest Site", xlab = "PC1", ylab = "Depth [mm]")
lines(scores1[,1],d1,ylim = rev(c(0,length(d1))))
plot(scores2[,1],d2,ylim = rev(c(0,length(d2))), main="P Elimination Plant", xlab = "PC1", ylab = "Depth [mm]")
lines(scores2[,1],d2,ylim = rev(c(0,length(d2))))
plot(-scores3[,1],d3,ylim = rev(c(0,length(d3))), main="River", xlab = "PC1", ylab = "Depth [mm]")
lines(-scores3[,1],d3,ylim = rev(c(0,length(d3))))
plot(scores4[,1],d4,ylim = rev(c(0,length(d4))), main="Main Basin", xlab = "PC1", ylab = "Depth [mm]")
lines(scores4[,1],d4,ylim = rev(c(0,length(d4))))

# Comparison of PC1 over all elements for all cores
layout(matrix(1:4,2,2))
barplot(loadings1[,1],ylim=c(-0.3,0.3), main="Deepest Site")
barplot(loadings2[,1],ylim=c(-0.3,0.3), main="P Elimination Plant")
barplot(loadings3[,1],ylim=c(-0.3,0.3), main="River")
barplot(loadings4[,1],ylim=c(-0.3,0.3), main="Main Basin")
# PEP + basin look similar

# Comparison of PC2 over all elements for all cores
layout(matrix(1:4,2,2))
barplot(loadings1[,2],ylim=c(-0.3,0.3), main="Deepest Site")
barplot(loadings2[,2],ylim=c(-0.3,0.3), main="P Elimination Plant")
barplot(loadings3[,2],ylim=c(-0.3,0.3), main="River")
barplot(loadings4[,2],ylim=c(-0.3,0.3), main="Main Basin")
# deepest + river look similar

# Comparison of PC1 over the depth for all cores in one figure
layout(1)
plot(scores1[,1],d1,ylim = rev(c(0,length(d1))),col="red")
lines(scores1[,1],d1,ylim = rev(c(0,length(d1))),col="red")
points(scores2[,1],d2,ylim = rev(c(0,length(d2))),col="blue")
lines(scores2[,1],d2,ylim = rev(c(0,length(d2))),col="blue")
points(scores3[,1],d3,ylim = rev(c(0,length(d3))),col="green")
lines(scores3[,1],d3,ylim = rev(c(0,length(d3))),col="green")
points(scores4[,1],d4,ylim = rev(c(0,length(d4))),col="yellow")
lines(scores4[,1],d4,ylim = rev(c(0,length(d4))),col="yellow")
legend(-9,600,c("Deepest Site","P Elimination Plant","River","Main Basin"),lty=c(1,1),
       lwd=c(2.5,2.5), col=c("red","blue","green","yellow"))
# correlation between PEP + main basin

# Comparison of PC2 over the depth for all cores in on figure
layout(1)
plot(scores1[,2],d1,ylim = rev(c(0,length(d1))),col="red")
lines(scores1[,2],d1,ylim = rev(c(0,length(d1))),col="red")
points(scores2[,2],d2,ylim = rev(c(0,length(d2))),col="blue")
lines(scores2[,2],d2,ylim = rev(c(0,length(d2))),col="blue")
points(scores3[,2],d3,ylim = rev(c(0,length(d3))),col="green")
lines(scores3[,2],d3,ylim = rev(c(0,length(d3))),col="green")
points(scores4[,2],d4,ylim = rev(c(0,length(d4))),col="yellow")
lines(scores4[,2],d4,ylim = rev(c(0,length(d4))),col="yellow")
legend(-9,600,c("Deepest Site","P Elimination Plant","River","Main Basin"),lty=c(1,1),
       lwd=c(2.5,2.5), col=c("red","blue","green","yellow"))
# correlation between deepest site + PEP, also river + main basin

# Comparison of PC3 over the depth for all cores in on figure
# NONSENSE!!! just a mess
layout(1)
plot(scores1[,3],d1,ylim = rev(c(0,length(d1))),col="red")
lines(scores1[,3],d1,ylim = rev(c(0,length(d1))),col="red")
points(scores2[,3],d2,ylim = rev(c(0,length(d2))),col="blue")
lines(scores2[,3],d2,ylim = rev(c(0,length(d2))),col="blue")
points(scores3[,3],d3,ylim = rev(c(0,length(d3))),col="green")
lines(scores3[,3],d3,ylim = rev(c(0,length(d3))),col="green")
points(scores4[,3],d4,ylim = rev(c(0,length(d4))),col="yellow")
lines(scores4[,3],d4,ylim = rev(c(0,length(d4))),col="yellow")
legend(-9,600,c("Deepest Site","P Elimination Plant","River","Main Basin"),lty=c(1,1),
       lwd=c(2.5,2.5), col=c("red","blue","green","yellow"))

# need some interpretation which PC explains how much for each core

###############
# Interpretation

# PC 1 shows depth of sediment core, also distinguishes different env. regime (shallow side, river, deepest site)
# PC 2 shows contamination history, deepest+PEP are heavily influenced by wastewaters, river+basin by river nutrients

###############
# ONE PCA FOR ALL CORES
# the fun starts here

T1$identifier<-rep(2,nrow(T1))
T2$identifier<-rep(3,nrow(T2))
T3$identifier<-rep(4,nrow(T3))
T4$identifier<-rep(1,nrow(T4))
# TotalT$identifier<-rep(1,nrow(T1)#c(1:nrow(T1)),c(nrow(T1)+1:nrow(T2)),c(nrow(T2)+1:nrow(T3)),c(nrow(T3)+1:nrow(T4)))
TotalT<-rbind(T1,T2,T3,T4)

TotalT$identifier <- factor(TotalT$identifier, colors=c("red","green","blue","black"))

zT<-scale(TotalT[,-ncol(TotalT)],center=TRUE,scale=TRUE) # scale is command for standardization by standard deviation
pca<-prcomp(zT,retx=T,center=F,scale.=F) # retx - factor holdings TRUE (scores), scale.=F bc. already done this
pca$sdev # the standard deviations of the PCA axes (their squares are the eigenvalues)

pca$sdev^2
sum(pca$sdev^2)
pca$sdev^2/sum(pca$sdev^2)
pca$x 
apply(pca$x,2,sd) 
scores<-pca$x
pca$rotation 
head(scores)
zD %*% pca$rotation 
loadings<-pca$rotation
summary(pca)

# for plotting as color, creating color ramp to have a transition from 
my.colorRamp.fct1<-colorRamp(c("darkred", "white")) # custom-built color ramp
my.colorRamp.fct1(0.5)
rgb(my.colorRamp.fct1(0.5),maxColorValue=255)

my.colorRamp.fct2<-colorRamp(c("darkgreen", "white")) # custom-built color ramp
my.colorRamp.fct2(0.5)
rgb(my.colorRamp.fct2(0.5),maxColorValue=255)

my.colorRamp.fct3<-colorRamp(c("darkblue", "white")) # custom-built color ramp
my.colorRamp.fct3(0.5)
rgb(my.colorRamp.fct3(0.5),maxColorValue=255)

my.colorRamp.fct4<-colorRamp(c("black", "white")) # custom-built color ramp
my.colorRamp.fct4(0.5)
rgb(my.colorRamp.fct4(0.5),maxColorValue=255)

###############
# Check importance of PCs

layout(1)

pca.pct<-100*round(summary(pca)$importance[2,],3)          
barplot(pca.pct, main="Lake Tegel")

###############
# PCA biplots #
# for a DISTANCE BIPLOT (focus is on sites, "scaling 1")
# each principal component has variance given by eigenvalue, loadings remain unscaled

# in this plot:
# 1) distances among sites are approximating true Euclidean distances in multivariate space
# 2) angles between arrows do not reflect correlations among variables
# 3) projecting site on descriptor at right angle gives its appr. descriptor value

layout(1)
plot(scores[,1:2],asp=1,pch=21, col = TotalT$identifier,
     xlim = xsi,ylim=ysi, main="Lake Tegel") # asp=1 x and y are equally scaled, bg=peak.colors, "1"=my.colorRamp.fct1, "2"=my.colorRamp.fct2, "3"=my.colorRamp.fct3 "4" =my.colorRamp.fct4
arrows<-loadings*ext
plotcircle(r=7*sqrt(2/ncol(TotalT[,-ncol(TotalT)])),lcol="blue") 
Arrows(x0=0,y0=0,x1=arrows[,1],y1=arrows[,2],col="darkgreen",arr.length=arleng)
text(x=arrows[,1]*1.3,y=arrows[,2]*1.2,labels=names(TotalT[,-ncol(TotalT)]),cex=1.5)
# legend("topleft",legend=levels(TotalT$identifier),col=1:4)

###############
# Comparison of PC1 over the depth for all cores
layout(1)
plot(scores[1:max(d1),1],d1,type="l",xlim = c(-7,7), ylim = rev(c(0,length(d1))), main="Lake Tegel", xlab = "PC1", ylab = "Depth [mm]", col="red")
lines(scores[(max(d1)+1):(max(d1)+max(d2)),1],d2,ylim = rev(c(0,length(d2))), col="green")
lines(scores[(max(d2)+max(d1)+1):(max(d1)+max(d2)+max(d3)),1],d3,ylim = rev(c(0,length(d3))), col="blue")
lines(scores[(max(d3)+max(d1)+max(d2)+1):(max(d1)+max(d2)+max(d3)+max(d4)),1],d4,ylim = rev(c(0,length(d4))), col="black")
legend("topleft",c("Deepest Site","P Elimination Plant","River","Main Basin"),lty=c(1,1),
       lwd=c(2.5,2.5), col=c("red","green","blue", "black"))

# Comparison of PC2 over the depth for all cores
layout(1)
plot(scores[1:max(d1),2],d1,type="l",xlim = c(-7,7), ylim = rev(c(0,length(d1))), main="Lake Tegel", xlab = "PC2", ylab = "Depth [mm]", col="red")
lines(scores[(max(d1)+1):(max(d1)+max(d2)),2],d2,ylim = rev(c(0,length(d2))), col="green")
lines(scores[(max(d2)+max(d1)+1):(max(d1)+max(d2)+max(d3)),2],d3,ylim = rev(c(0,length(d3))), col="blue")
lines(scores[(max(d3)+max(d1)+max(d2)+1):(max(d1)+max(d2)+max(d3)+max(d4)),2],d4,ylim = rev(c(0,length(d4))), col="black")
legend("topleft",c("Deepest Site","P Elimination Plant","River","Main Basin"),lty=c(1,1),
       lwd=c(2.5,2.5), col=c("red","green","blue", "black"))

# Comparison of PC3 over the depth for all cores
layout(1)
plot(scores[1:max(d1),3],d1,type="l",xlim = c(-7,7), ylim = rev(c(0,length(d1))), main="Lake Tegel", xlab = "PC3", ylab = "Depth [mm]", col="red")
lines(scores[(max(d1)+1):(max(d1)+max(d2)),3],d2,ylim = rev(c(0,length(d2))), col="green")
lines(scores[(max(d2)+max(d1)+1):(max(d1)+max(d2)+max(d3)),3],d3,ylim = rev(c(0,length(d3))), col="blue")
lines(scores[(max(d3)+max(d1)+max(d2)+1):(max(d1)+max(d2)+max(d3)+max(d4)),3],d4,ylim = rev(c(0,length(d4))), col="black")
legend("topleft",c("Deepest Site","P Elimination Plant","River","Main Basin"),lty=c(1,1),
       lwd=c(2.5,2.5), col=c("red","green","blue", "black"))

# Comparison of PC1 over all elements for all cores
layout(matrix(1:3,3,1))
barplot(loadings[,1],ylim=c(-0.4,0.7), main="PC1")
barplot(loadings[,2],ylim=c(-0.4,0.7), main="PC2")
barplot((-1)*loadings[,3],ylim=c(-0.4,0.7), main="PC3")

# Comparison of PC1/PC1 for all cores
layout(matrix(1:4,2,2))
plot(scores[1:max(d1),1],d1,type="l",xlim = c(-7,7), ylim = rev(c(0,length(d1))), main= "Deepest Site", xlab = "PC1", ylab = "Depth [mm]", col="red")
lines(scores1[,1],d1,ylim = rev(c(0,length(d1))))
plot(scores[(max(d1)+1):(max(d1)+max(d2)),1],d2,type="l",ylim = rev(c(0,length(d2))), main= "P Elimination Plant", xlab = "PC1", ylab = "Depth [mm]", col="green")
lines(scores2[,1],d2,ylim = rev(c(0,length(d2))))
plot(scores[(max(d2)+max(d1)+1):(max(d1)+max(d2)+max(d3)),1],d3,type="l",ylim = rev(c(0,length(d3))), main= "River", xlab = "PC1", ylab = "Depth [mm]", col="blue")
lines(-scores3[,1],d3,ylim = rev(c(0,length(d3))))
plot(scores[(max(d3)+max(d1)+max(d2)+1):(max(d1)+max(d2)+max(d3)+max(d4)),1],d4,type="l",ylim = rev(c(0,length(d4))),main= "Main Basin", xlab = "PC1", ylab = "Depth [mm]",  col="darkorange")
lines(scores4[,1],d4,ylim = rev(c(0,length(d4))))

# Surpirse
layout(matrix(1:4,2,2))
plot(scores[1:max(d1),1],d1,type="l",xlim = c(-7,7), ylim = rev(c(0,length(d1))), main= "Deepest Site", xlab = "PC1", ylab = "Depth [mm]", col="red")
lines(-scores1[,2],d1,ylim = rev(c(0,length(d1))))
plot(scores[(max(d1)+1):(max(d1)+max(d2)),1],d2,type="l",ylim = rev(c(0,length(d2))), main= "P Elimination Plant", xlab = "PC1", ylab = "Depth [mm]", col="green")
lines(scores2[,1],d2,ylim = rev(c(0,length(d2))))
plot(scores[(max(d2)+max(d1)+1):(max(d1)+max(d2)+max(d3)),1],d3,type="l",ylim = rev(c(0,length(d3))), main= "River", xlab = "PC1", ylab = "Depth [mm]", col="blue")
lines(-scores3[,1],d3,ylim = rev(c(0,length(d3))))
plot(scores[(max(d3)+max(d1)+max(d2)+1):(max(d1)+max(d2)+max(d3)+max(d4)),1],d4,type="l",ylim = rev(c(0,length(d4))),main= "Main Basin", xlab = "PC1", ylab = "Depth [mm]",  col="darkorange")
lines(scores4[,1],d4,ylim = rev(c(0,length(d4))))

# Comparison of PC1/PC1 for all cores
layout(matrix(1:4,2,2))
plot(-scores[1:max(d1),3],d1,type="l",xlim = c(-7,7), ylim = rev(c(0,length(d1))), main= "Deepest Site", xlab = "PC1", ylab = "Depth [mm]", col="red")
lines(scores1[,1],d1,ylim = rev(c(0,length(d1))))
plot(scores[(max(d1)+1):(max(d1)+max(d2)),1],d2,type="l",ylim = rev(c(0,length(d2))), main= "P Elimination Plant", xlab = "PC1", ylab = "Depth [mm]", col="green")
lines(scores2[,1],d2,ylim = rev(c(0,length(d2))))
plot(scores[(max(d2)+max(d1)+1):(max(d1)+max(d2)+max(d3)),1],d3,type="l",ylim = rev(c(0,length(d3))), main= "River", xlab = "PC1", ylab = "Depth [mm]", col="blue")
lines(-scores3[,1],d3,ylim = rev(c(0,length(d3))))
plot(scores[(max(d3)+max(d1)+max(d2)+1):(max(d1)+max(d2)+max(d3)+max(d4)),1],d4,type="l",ylim = rev(c(0,length(d4))),main= "Main Basin", xlab = "PC1", ylab = "Depth [mm]",  col="darkorange")
lines(scores4[,1],d4,ylim = rev(c(0,length(d4))))
