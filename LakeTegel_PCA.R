################# General information
# TEG source: microRFA scans by GEOPOLAR Bremen done for four sediment cores
# sediment cores were taken at Lake Tegel on 30.09.15
# TEG-1 deepest lake site, 2 near PEP, 3 Tegelort, 4 outer range of main basin

# loading of additional packages for plotting
library(shape) 
library(vegan)

# loading of data sets for each sediment core
setwd("D:/work/data/Hupfer") # work directory, change to own one if it differs from default
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
plot(scores1[,1:2],asp=1,pch=21,bg=peak.colors, xlim = xsi,ylim=ysi, main="Deepest Site") # asp=1 x and y are equally scaled
arrows<-loadings1*ext # with extension factor to get nice graphs
plotcircle(r=7*sqrt(2/ncol(T1)),lcol="blue") # circle of equilibrium contribution (equal contribution to all PCA-dimensions)
Arrows(x0=0,y0=0,x1=arrows[,1],y1=arrows[,2],col="darkgreen",arr.length=arleng)
# circle shows equilibrium, arrows that are longer than circle are more important for PCA
text(x=arrows[,1]*1.3,y=arrows[,2]*1.2,labels=names(T1),cex=1.5)

as.rgb.channels<-my.colorRamp.fct(decostand(d2,method="range"))
peak.colors<-rgb(as.rgb.channels,maxColorValue=255)
plot(scores2[,1:2],asp=1,pch=21,bg=peak.colors, xlim = xsi,ylim=ysi, main="P Elimination Plant")
arrows<-loadings2*ext 
plotcircle(r=7*sqrt(2/ncol(T2)),lcol="blue") 
Arrows(x0=0,y0=0,x1=arrows[,1],y1=arrows[,2],col="darkgreen",arr.length=arleng)
text(x=arrows[,1]*1.3,y=arrows[,2]*1.2,labels=names(T2),cex=1.5)

as.rgb.channels<-my.colorRamp.fct(decostand(d3,method="range"))
peak.colors<-rgb(as.rgb.channels,maxColorValue=255)
plot(scores3[,1:2],asp=1,pch=21,bg=peak.colors, xlim = xsi,ylim=ysi, main="River") 
arrows<-loadings3*ext 
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

var(scores2[,1]/pca2$sdev[1]) 
plot(scores2[,1]/pca2$sdev[1],scores2[,2]/pca2$sdev[2],pch=21,bg=peak.colors,asp=1)
arrows<-loadings2*matrix(pca2$sdev,nrow=nrow(loadings2),ncol=ncol(loadings2),byrow=TRUE) 
arrows<-arrows*2 
Arrows(x0=0,y0=0,x1=arrows[,1],y1=arrows[,2],col="purple")
(structure<-cor(T2,scores2))
structure<-2*structure
Arrows(x0=0,y0=0,x1=structure[,1],y1=structure[,2],col="red")
text(x=arrows[,1]*1.3,y=arrows[,2]*1.2,labels=names(T2),cex=1.5)

var(scores3[,1]/pca3$sdev[1]) 
plot(scores3[,1]/pca3$sdev[1],scores3[,2]/pca3$sdev[2],pch=21,bg=peak.colors,asp=1)
arrows<-loadings3*matrix(pca3$sdev,nrow=nrow(loadings3),ncol=ncol(loadings3),byrow=TRUE) 
arrows<-arrows*2 
Arrows(x0=0,y0=0,x1=arrows[,1],y1=arrows[,2],col="purple")
(structure<-cor(T3,scores3))
structure<-2*structure
Arrows(x0=0,y0=0,x1=structure[,1],y1=structure[,2],col="red")
text(x=arrows[,1]*1.3,y=arrows[,2]*1.2,labels=names(T3),cex=1.5)

var(scores4[,1]/pca4$sdev[1])
plot(scores4[,1]/pca4$sdev[1],scores4[,2]/pca4$sdev[2],pch=21,bg=peak.colors,asp=1)
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