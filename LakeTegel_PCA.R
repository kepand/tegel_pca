################# TEG1 information
# TEG1 source: microRFA scans by GEOPOLAR Bremen done for four sediment cores
# sediment cores were taken at Lake Tegel on 30.09.15
# TEG-1 deepest lake site, 2 near PEP, 3 Tegelort, 4 outer range of main basin

<<<<<<< Updated upstream
# loading of additional packages for plotting
=======
# Testing testing, 1, 2, 3...

>>>>>>> Stashed changes
library(shape) # just plots nice arrows and circles
library(vegan)

# loading of data sets for each sediment core
setwd("D:/work/data/Hupfer") # work directory
TEG1<-read.table(file="PCA/TEG1_1mm.txt",header=T,sep="",skip=2) 
TEG2<-read.table(file="PCA/TEG2_1mm.txt",header=T,sep="",skip=2) 
TEG3<-read.table(file="PCA/TEG3_1mm.txt",header=T,sep="",skip=2) 
TEG4<-read.table(file="PCA/TEG4_1mm.txt",header=T,sep="",skip=2) 

# calls sub function and calculates pca
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

as.rgb.channels<-my.colorRamp.fct(decostand(d,method="range"))
d1.colors<-rgb(as.rgb.channels,maxColorValue=255)

###############
# PCA biplots #
# for a DISTANCE BIPLOT (focus is on sites, "scaling 1")
# each principal component has variance given by eigenvalue, loadings remain unscaled
layout(matrix(1:4,2,2))
xsi<-c(-10,10)
ysi<-c(-5,5)
ext<-7
arleng<-0.15

as.rgb.channels<-my.colorRamp.fct(decostand(d1,method="range"))
peak.colors<-rgb(as.rgb.channels,maxColorValue=255)
plot(scores1[,1:2],asp=1,pch=21,bg=peak.colors, xlim = xsi,ylim=ysi)#, xlim = c(-10,10),ylim=c(-8,8)) # asp=1 x and y are equally scaled
arrows<-loadings1*ext # with extension factor to get nice graphs
plotcircle(r=7*sqrt(2/ncol(T1)),lcol="blue") # circle of equilibrium contribution (equal contribution to all PCA-dimensions)
Arrows(x0=0,y0=0,x1=arrows[,1],y1=arrows[,2],col="darkgreen",arr.length=arleng)
# circle shows equilibrium, arrows that are longer than circle are more important for PCA
text(x=arrows[,1]*1.3,y=arrows[,2]*1.2,labels=names(T1),cex=1.5)

as.rgb.channels<-my.colorRamp.fct(decostand(d2,method="range"))
peak.colors<-rgb(as.rgb.channels,maxColorValue=255)
plot(scores2[,1:2],asp=1,pch=21,bg=peak.colors, xlim = xsi,ylim=ysi)#, xlim = c(-10,10),ylim=c(-8,8)) # asp=1 x and y are equally scaled
arrows<-loadings2*ext # with extension factor to get nice graphs
plotcircle(r=7*sqrt(2/ncol(T2)),lcol="blue") # circle of equilibrium contribution (equal contribution to all PCA-dimensions)
Arrows(x0=0,y0=0,x1=arrows[,1],y1=arrows[,2],col="darkgreen",arr.length=arleng)
# circle shows equilibrium, arrows that are longer than circle are more important for PCA
text(x=arrows[,1]*1.3,y=arrows[,2]*1.2,labels=names(T2),cex=1.5)

as.rgb.channels<-my.colorRamp.fct(decostand(d3,method="range"))
peak.colors<-rgb(as.rgb.channels,maxColorValue=255)
plot(scores3[,1:2],asp=1,pch=21,bg=peak.colors, xlim = xsi,ylim=ysi) #, xlim = c(-10,10),ylim=c(-8,8)) # asp=1 x and y are equally scaled
arrows<-loadings3*ext # with extension factor to get nice graphs
plotcircle(r=7*sqrt(2/ncol(T3)),lcol="blue") # circle of equilibrium contribution (equal contribution to all PCA-dimensions)
Arrows(x0=0,y0=0,x1=arrows[,1],y1=arrows[,2],col="darkgreen",arr.length=arleng)
# circle shows equilibrium, arrows that are longer than circle are more important for PCA
text(x=arrows[,1]*1.3,y=arrows[,2]*1.2,labels=names(T1),cex=1.5)

as.rgb.channels<-my.colorRamp.fct(decostand(d4,method="range"))
peak.colors<-rgb(as.rgb.channels,maxColorValue=255)
plot(scores4[,1:2],asp=1,pch=21,bg=peak.colors, xlim = xsi,ylim=ysi) # asp=1 x and y are equally scaled
arrows<-loadings4*ext # with extension factor to get nice graphs
plotcircle(r=7*sqrt(2/ncol(T4)),lcol="blue") # circle of equilibrium contribution (equal contribution to all PCA-dimensions)
Arrows(x0=0,y0=0,x1=arrows[,1],y1=arrows[,2],col="darkgreen",arr.length=arleng)
# circle shows equilibrium, arrows that are longer than circle are more important for PCA
text(x=arrows[,1]*1.3,y=arrows[,2]*1.2,labels=names(T1),cex=1.5)

###############
# CORRELATION BIPLOT
layout(matrix(1:4,2,2))

# arrows which reach out of the circle contribute more than on average
# biplot(pca1,scale=0)
# for a CORRELATION BIPLOT (focus is on variables, "scaling 2") # also for paper: CORRELATION/DISTANCE BIPLOT
# each principal component is weighted by 1/sqrt(eigenvalue), so it has variance 1
var(scores1[,1]/pca1$sdev[1]) # just demo
plot(scores1[,1]/pca1$sdev[1],scores1[,2]/pca1$sdev[2],pch=21,bg=peak.colors,asp=1)
# loadings are weighted by sqrt(eigenvalues) (multiplied by sqrt(eigenvalues))
arrows<-loadings1*matrix(pca1$sdev,nrow=nrow(loadings1),ncol=ncol(loadings1),byrow=TRUE) # scaling by STD
arrows<-arrows*2 # choose extension factor
Arrows(x0=0,y0=0,x1=arrows[,1],y1=arrows[,2],col="purple")
# as alternative just compute correlation of scores with original data ("structure coefficients")
(structure<-cor(T1,scores1))
structure<-2*structure
Arrows(x0=0,y0=0,x1=structure[,1],y1=structure[,2],col="red")
text(x=arrows[,1]*1.3,y=arrows[,2]*1.2,labels=names(T1),cex=1.5)

var(scores2[,1]/pca2$sdev[1]) # just demo
plot(scores2[,1]/pca2$sdev[1],scores2[,2]/pca2$sdev[2],pch=21,bg=peak.colors,asp=1)
# loadings are weighted by sqrt(eigenvalues) (multiplied by sqrt(eigenvalues))
arrows<-loadings2*matrix(pca2$sdev,nrow=nrow(loadings2),ncol=ncol(loadings2),byrow=TRUE) # scaling by STD
arrows<-arrows*2 # choose extension factor
Arrows(x0=0,y0=0,x1=arrows[,1],y1=arrows[,2],col="purple")
# as alternative just compute correlation of scores with original data ("structure coefficients")
(structure<-cor(T2,scores2))
structure<-2*structure
Arrows(x0=0,y0=0,x1=structure[,1],y1=structure[,2],col="red")
text(x=arrows[,1]*1.3,y=arrows[,2]*1.2,labels=names(T2),cex=1.5)

var(scores3[,1]/pca3$sdev[1]) # just demo
plot(scores3[,1]/pca3$sdev[1],scores3[,2]/pca3$sdev[2],pch=21,bg=peak.colors,asp=1)
# loadings are weighted by sqrt(eigenvalues) (multiplied by sqrt(eigenvalues))
arrows<-loadings3*matrix(pca3$sdev,nrow=nrow(loadings3),ncol=ncol(loadings3),byrow=TRUE) # scaling by STD
arrows<-arrows*2 # choose extension factor
Arrows(x0=0,y0=0,x1=arrows[,1],y1=arrows[,2],col="purple")
# as alternative just compute correlation of scores with original data ("structure coefficients")
(structure<-cor(T3,scores3))
structure<-2*structure
Arrows(x0=0,y0=0,x1=structure[,1],y1=structure[,2],col="red")
text(x=arrows[,1]*1.3,y=arrows[,2]*1.2,labels=names(T3),cex=1.5)

var(scores4[,1]/pca4$sdev[1]) # just demo
plot(scores4[,1]/pca4$sdev[1],scores4[,2]/pca4$sdev[2],pch=21,bg=peak.colors,asp=1)
# loadings are weighted by sqrt(eigenvalues) (multiplied by sqrt(eigenvalues))
arrows<-loadings4*matrix(pca4$sdev,nrow=nrow(loadings4),ncol=ncol(loadings4),byrow=TRUE) # scaling by STD
arrows<-arrows*2 # choose extension factor
Arrows(x0=0,y0=0,x1=arrows[,1],y1=arrows[,2],col="purple")
# as alternative just compute correlation of scores with original data ("structure coefficients")
(structure<-cor(T4,scores4))
structure<-2*structure
Arrows(x0=0,y0=0,x1=structure[,1],y1=structure[,2],col="red")
text(x=arrows[,1]*1.3,y=arrows[,2]*1.2,labels=names(T4),cex=1.5)

#biplot(pca1,scale=1)

###############
# Comparison of PC1 over the depth for all cores
layout(matrix(1:4,2,2))
plot(scores1[,1],d1,ylim = rev(c(0,length(d1))))
lines(scores1[,1],d1,ylim = rev(c(0,length(d1))))
plot(scores2[,1],d2,ylim = rev(c(0,length(d2))))
lines(scores2[,1],d2,ylim = rev(c(0,length(d2))))
plot(scores3[,1],d3,ylim = rev(c(0,length(d3))))
lines(scores3[,1],d3,ylim = rev(c(0,length(d3))))
plot(scores4[,1],d4,ylim = rev(c(0,length(d4))))
lines(scores4[,1],d4,ylim = rev(c(0,length(d4))))

# Comparison of PC1 over all elements for all cores
layout(matrix(1:4,2,2))
barplot(loadings1[,1],ylim=c(-0.3,0.3))
barplot(loadings2[,1],ylim=c(-0.3,0.3))
barplot(loadings3[,1],ylim=c(-0.3,0.3))
barplot(loadings4[,1],ylim=c(-0.3,0.3))
