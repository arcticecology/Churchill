library(riojaPlot)
library(readxl)
library(dplyr)
#library(rioja)
library(analogue)
library(vegan)
library(ggpalaeo)
library(tidyverse)
library(mgcv)
library(draw)
#setwd("C:\\Users\\andre\\OneDrive\\Desktop\\2024\\Connor\\Churchill")
coreinput<-read.csv(file="left_2023.csv", row.names=1)


BCI<-coreinput

chron <- cbind(BCI[,1] , BCI[,2])
colnames(chron) <- c("Depth","Year")
chron=as.data.frame(chron)

BCI<-BCI[,-cbind(1:2)] #first two columns
core <- BCI / rowSums(BCI) * 100

BCI2<-round(BCI, digits = 0)
S <- specnumber(BCI2) # observed number of species
(raremax <- min(rowSums(BCI2)))
Srare <- rarefy(BCI2, raremax)

# Fortin et al. 2015 used filtered >=2%-in-2-lakes spp criteria in their analyses
# this code is to retain only taxa that are >=2% in 2 intervals in the training set
N <- 2
M <- 2
i <- colSums(core >= M) >= N
spp_red <- core[, i, drop = FALSE]
# 61 chironomid taxa remain in training set
ncol(spp_red)

#this is for spider, check for others
fos<-spp_red[,-cbind(20:23)] #remove undifferentiated left lake

nms <- colnames(fos)
nms3<-gsub("\\.", " ", nms)




depth<-as.numeric(chron$Depth)



#diversity

H<-diversity(core)
H<-as.data.frame(H)
J <- H/log(specnumber(core))
J<-as.data.frame(J)
#DCA and PCA
dca <- decorana(fos, iweigh=1)
sc1 <- scores(dca, display="sites", choices=1:2)
DC1<-as.data.frame(sc1)

H1.pca <- rda(fos)
sc.pca <-scores(H1.pca,display="sites", choices=1:2)
PC1<-as.data.frame(sc.pca)

div2<-data.frame(PC1,H, Srare)


isotopes_spider <-read.csv(file="left_2023_isotopes.csv", row.names=1)
recon_left <-read.csv(file="reconstruction_left.csv", row.names=1)

re04<-as.data.frame(recon_left$Comp01)
colnames(re04) <- c("Comp01")


err1<-as.data.frame(recon_left$err1)
colnames(err1) <- c("err1")


iso1<-isotopes_spider[,-cbind(1, 2, 3,4,5)] #d13c d15n d18o
iso2<-isotopes_spider[,-cbind(1, 2, 6,7,8)]# %C %N C:N
colnames(iso2) <- c("%C", "%N", "C:N")

colnames(iso1) <- c("δ13C", "δ15N", "δ18O")



diss <- dist(fos/100)
clust <- chclust(diss)
bstick(clust, 5)


tmp2 <- data.frame(x=re04$Comp01, y=chron$Year)
sdev<-err1$err1

fun.gam3 <- function(x, y, i, nm, style) {
  lo <- loess(re04$Comp01 ~ chron$Year)
  lines(lo$fitted, chron$Year, col = "darkgray", lwd=2)
  arrows(re04$Comp01-sdev, chron$Year, re04$Comp01+sdev, chron$Year, length=0.05, angle=180, code=3)
  
}

#gam1 <- mgcv::gam(chron$Year ~ s(re04$Comp02, k=5))
tiff("Left_2023a.tiff", width = 6, height = 4, units = 'in', res = 300)

rp1 <- riojaPlot(fos, chron, 
                 ymin = 1600, ymax=2020, yinterval = 20, 
                 yvar.name = "Year",
                 y.rev=FALSE,
                 yTop=0.6,
                 yBottom=0.20,
                 scale.percent=TRUE, plot.groups=FALSE, do.clust = TRUE,
                 plot.poly=FALSE, plot.line=FALSE, lwd.bar=3,
                 plot.clust=TRUE,  labels.italicise=TRUE,
                 plot.cumul=TRUE, cex.cumul=0.5, xlabels=nms3, labels.break.long=FALSE,
                 srt.xlabel=45, plot.bar=TRUE,
                 tcl=-0.1, cex.yaxis=0.4, cex.xlabel=0.4,  cex.ylabel=0.5, col.bar="black",
                 cex.xaxis=0.4, xRight = 0.65, plot.line=FALSE)

#par(fig=c(0.8,1,0,1), new=TRUE)


rp2 <- riojaPlot(iso2, isotopes_spider, yvar.name="Year",
                 riojaPlot=rp1,xGap = 0.01,
                 xRight=0.78, 
                 scale.minmax=FALSE, plot.bar=F, las.xaxis=2, cex.xlabel=0.4,
                 plot.line=T, cex.xaxis=0.4, plot.symb=TRUE, symb.cex=0.5)

rp3 <- riojaPlot(iso1, isotopes_spider, yvar.name="Year",
                 riojaPlot=rp2,xGap = 0.01,
                 xRight=0.92, las.xaxis=2,
                 scale.minmax=FALSE, plot.bar=F, cex.xlabel=0.4,
                 plot.line=T, cex.xaxis=0.4, plot.symb=TRUE, symb.cex=0.5)

rp4 <-riojaPlot(re04, chron, 
          riojaPlot=rp3, xGap = 0.01, yvar.name="Year",
          xRight=0.955, scale.minmax=FALSE, 
          plot.bar=FALSE, plot.line=FALSE, 
          plot.symb=TRUE, cex.xaxis=0.4, symb.cex=0.5, fun.xfront=fun.gam3)

drawBox(x=4.55, y = 1.88, width = 8.25, height = 0.12,  opacity=0.6, fillColor="grey")


#rp3a <- addRPClustZone(rp1, clust, nZone=2, xRight=0.70, col="#59595988", lwd=23, lty=1, alpha=0.5)
rp3b <- addRPClustZone(rp4, clust, nZone=2, xRight=0.99, col="#59595988", lty=2, lwd=2)
#
dev.off()
#
#
#esm2
tiff("Left_ESM2.tiff", width = 6, height = 4, units = 'in', res = 300)
rp1 <- riojaPlot(div2, chron, 
                 ymin = 1600, ymax=2020, yinterval = 20, 
                 yvar.name = "Year",
                 y.rev=FALSE,srt.xlabel=45,
                 yTop=0.6,   labels.italicise=TRUE,
                 yBottom=0.20,xRight = 0.60,
                 scale.minmax=FALSE, plot.bar=F, las.xaxis=2, cex.xlabel=0.4,
                 plot.line=T, cex.xaxis=0.4, plot.symb=TRUE, symb.cex=0.5)

rp2 <-riojaPlot(re04, chron, 
                riojaPlot=rp1, xGap = 0.01, yvar.name="Year",
                xRight=0.8, scale.minmax=FALSE, srt.xlabel=45,do.clust = TRUE,plot.clust=TRUE, 
                plot.bar=FALSE, plot.line=FALSE, 
                plot.symb=TRUE, cex.xaxis=0.4, symb.cex=0.5, fun.xfront=fun.gam3)

drawBox(x=2.55, y = 1.88, width = 3.8, height = 0.12,  opacity=0.6, fillColor="grey")


#rp3a <- addRPClustZone(rp1, clust, nZone=2, xRight=0.70, col="#59595988", lwd=23, lty=1, alpha=0.5)
rp3b <- addRPClustZone(rp2, clust, nZone=2, xRight=0.8, col="#59595988", lty=2, lwd=2)
#
dev.off()



#
#
#

coreinput<-read.csv(file="larch_2023.csv", row.names=1)


BCI<-coreinput

chron <- cbind(BCI[,1] , BCI[,2])
colnames(chron) <- c("Depth","Year")
chron=as.data.frame(chron)

BCI<-BCI[,-cbind(1:2)] #first two columns
core <- BCI / rowSums(BCI) * 100
BCI2<-round(BCI, digits = 0)
S <- specnumber(BCI2) # observed number of species
(raremax <- min(rowSums(BCI2)))
Srare <- rarefy(BCI2, raremax)
# Fortin et al. 2015 used filtered >=2%-in-2-lakes spp criteria in their analyses
# this code is to retain only taxa that are >=2% in 2 intervals in the training set
N <- 2
M <- 2
i <- colSums(core >= M) >= N
spp_red <- core[, i, drop = FALSE]
# 61 chironomid taxa remain in training set
ncol(spp_red)

#this is for spider, check for others
fos<-spp_red[,-cbind(20:23)] #remove undifferentiated left lake

nms <- colnames(fos)
nms3<-gsub("\\.", " ", nms)




depth<-as.numeric(chron$Depth)



#diversity

H<-diversity(core)
H<-as.data.frame(H)
J <- H/log(specnumber(core))
J<-as.data.frame(J)
#DCA and PCA
dca <- decorana(fos, iweigh=1)
sc1 <- scores(dca, display="sites", choices=1:2)
DC1<-as.data.frame(sc1)

H1.pca <- rda(fos)
sc.pca <-scores(H1.pca,display="sites", choices=1:2)
PC1<-as.data.frame(sc.pca)

div2<-data.frame(PC1,H, Srare)


isotopes_spider <-read.csv(file="larch_2023_isotopes.csv", row.names=1)
recon_left <-read.csv(file="reconstruction_larch.csv", row.names=1)

re04<-as.data.frame(recon_left$Comp01)
colnames(re04) <- c("Comp01")


err1<-as.data.frame(recon_left$err1)
colnames(err1) <- c("err1")


iso1<-isotopes_spider[,-cbind(1, 2, 3,4,5)] #d13c d15n d18o
iso2<-isotopes_spider[,-cbind(1, 2, 6,7,8)]# %C %N C:N
colnames(iso2) <- c("%C", "%N", "C:N")

colnames(iso1) <- c("δ13C", "δ15N", "δ18O")



diss <- dist(fos/100)
clust <- chclust(diss)
bstick(clust, 5)


tmp2 <- data.frame(x=re04$Comp01, y=chron$Year)
sdev<-err1$err1

fun.gam3 <- function(x, y, i, nm, style) {
  lo <- loess(re04$Comp01 ~ chron$Year)
  lines(lo$fitted, chron$Year, col = "darkgray", lwd=2)
  arrows(re04$Comp01-sdev, chron$Year, re04$Comp01+sdev, chron$Year, length=0.05, angle=180, code=3)
  
}


#gam1 <- mgcv::gam(chron$Year ~ s(re04$Comp02, k=5))
tiff("Larch_2023a.tiff", width = 6, height = 4, units = 'in', res = 300)

rp1 <- riojaPlot(fos, chron, 
                 ymin = 1750, ymax=2020, yinterval = 20, 
                 yvar.name = "Year",
                 y.rev=FALSE,
                 yTop=0.6,
                 yBottom=0.20,
                 scale.percent=TRUE, plot.groups=FALSE, do.clust = TRUE,
                 plot.poly=FALSE, plot.line=FALSE, lwd.bar=3,
                 plot.clust=TRUE,  labels.italicise=TRUE,
                 plot.cumul=TRUE, cex.cumul=0.5, xlabels=nms3, labels.break.long=FALSE,
                 srt.xlabel=45, plot.bar=TRUE,
                 tcl=-0.1, cex.yaxis=0.4, cex.xlabel=0.4,  cex.ylabel=0.5, col.bar="black",
                 cex.xaxis=0.4, xRight = 0.65, plot.line=FALSE)

#par(fig=c(0.8,1,0,1), new=TRUE)


rp2 <- riojaPlot(iso2, isotopes_spider, yvar.name="Year",
                 riojaPlot=rp1,xGap = 0.01,
                 xRight=0.78, 
                 scale.minmax=FALSE, plot.bar=F, las.xaxis=2, cex.xlabel=0.4,
                 plot.line=T, cex.xaxis=0.4, plot.symb=TRUE, symb.cex=0.5)

rp3 <- riojaPlot(iso1, isotopes_spider, yvar.name="Year",
                 riojaPlot=rp2,xGap = 0.01,
                 xRight=0.92, las.xaxis=2,
                 scale.minmax=FALSE, plot.bar=F, cex.xlabel=0.4,
                 plot.line=T, cex.xaxis=0.4, plot.symb=TRUE, symb.cex=0.5)

rp4 <-riojaPlot(re04, chron, 
                riojaPlot=rp3, xGap = 0.01, yvar.name="Year",
                xRight=0.955, scale.minmax=FALSE, 
                plot.bar=FALSE, plot.line=FALSE, 
                plot.symb=TRUE, cex.xaxis=0.4, symb.cex=0.5, fun.xfront=fun.gam3)

drawBox(x=4.55, y = 1.6, width = 8.25, height = 0.3,  opacity=0.6, fillColor="grey")


#rp3a <- addRPClustZone(rp1, clust, nZone=2, xRight=0.70, col="#59595988", lwd=23, lty=1)
rp3b <- addRPClustZone(rp4, clust, nZone=2, xRight=0.99, col="darkgrey", lty=2, lwd=2)
#
dev.off()

#
#
#ESM2
tiff("Larch_ESM2.tiff", width = 6, height = 4, units = 'in', res = 300)
rp1 <- riojaPlot(div2, chron, 
                 ymin = 1750, ymax=2020, yinterval = 20, 
                 yvar.name = "Year",
                 y.rev=FALSE,srt.xlabel=45,
                 yTop=0.6,   labels.italicise=TRUE,
                 yBottom=0.20,xRight = 0.60,
                 scale.minmax=FALSE, plot.bar=F, las.xaxis=2, cex.xlabel=0.4,
                 plot.line=T, cex.xaxis=0.4, plot.symb=TRUE, symb.cex=0.5)

rp2 <-riojaPlot(re04, chron, 
                riojaPlot=rp1, xGap = 0.01, yvar.name="Year",
                xRight=0.8, scale.minmax=FALSE, srt.xlabel=45,do.clust = TRUE,plot.clust=TRUE, 
                plot.bar=FALSE, plot.line=FALSE, 
                plot.symb=TRUE, cex.xaxis=0.4, symb.cex=0.5, fun.xfront=fun.gam3)


drawBox(x=2.55, y = 1.6, width = 3.8, height = 0.3,  opacity=0.6, fillColor="grey")


#rp3a <- addRPClustZone(rp1, clust, nZone=2, xRight=0.70, col="#59595988", lwd=23, lty=1, alpha=0.5)
rp3b <- addRPClustZone(rp2, clust, nZone=2, xRight=0.8, col="#59595988", lty=2, lwd=2)
#
dev.off()

#
#
#

Larch <-read.csv(file="reconstruction_larch.csv", row.names=1)
Left <-read.csv(file="reconstruction_left.csv", row.names=1)
metstation <-read.csv(file="climate_churchill.csv", row.names=1)
library(ggplot2)

#ESM5

p4 <- ggplot()+
  ggtitle("Churchill Temperature 1930-2010")+ # for the main title
  xlab("Year")+ # for the x axis label
  ylab(expression("Temperature Anomaly ("*~degree*C*")"))+ # for the y axis label
  geom_point(data=Larch, aes(x=Year.CE, y=Tanon, colour="Larch"), shape=15) + xlim(1900,2010) +
  geom_smooth(data=Larch, aes(x=Year.CE, y=Tanon), method=lm , color="darkgreen", fill="#69b3a2", se=TRUE) +
  geom_point(data=Left, aes(x=Year.CE, y=Tanon, colour="Left"), shape=17) +
  geom_smooth(data=Left, aes(x=Year.CE, y=Tanon), method=lm , color="red", fill="#FF9999", se=TRUE) +
  geom_point(data=metstation, aes(x=Year, y=Tanon, colour="MET")) +
  geom_smooth(data=metstation, aes(x=Year, y=Tanon), method=lm , color="blue", fill="#99CCFF", se=TRUE) +
  labs(colour="Legend")

p4 + theme(legend.position = c(0.1, 0.8))


p5 <- ggplot()+
  ggtitle("Churchill Temperature 1800-2010")+ # for the main title
  xlab("Year") + # for the x axis label
  ylab(expression("Temperature Anomaly ("*~degree*C*")")) + # for the y axis label
  geom_point(data=Larch, aes(x=Year.CE, y=Tanon, colour="Larch"), shape=17) + + xlim(1800,2010) +
  geom_smooth(data=Larch, aes(x=Year.CE, y=Tanon), method = "gam", formula = y ~ s(x, k = 3), size = 1 , color="red", fill="#FF9999", se=TRUE) +
  geom_point(data=Left, aes(x=Year.CE, y=Tanon, colour="Left"), shape=15) + xlim(1800,2010) +
  geom_smooth(data=Left, aes(x=Year.CE, y=Tanon), method = "gam", formula = y ~ s(x, k = 3), size = 1 , color="darkgreen", fill="#69b3a2", se=TRUE) +
  geom_point(data=metstation, aes(x=Year, y=Tanon, colour="MET")) +
  geom_smooth(data=metstation, aes(x=Year, y=Tanon), method = "gam", formula = y ~ s(x, k = 3), size = 1 , color="blue", fill="#99CCFF", se=TRUE) +
  labs(lake="Legend")

p5 + theme(legend.position = c(0.1, 0.8))

#
tiff("anomaly.tiff", width = 5, height = 4, units = 'in', res = 300)
p6 <- ggplot()+
  ggtitle("Churchill Temperature 1900-2010")+ # for the main title
  xlab("Year")+ # for the x axis label
  ylab(expression("Temperature Anomaly ("*~degree*C*")"))+ # for the y axis label
  geom_point(data=Larch, aes(x=Year.CE, y=Tanon2, colour="Larch"), shape=15) + xlim(1900,2010) +
  geom_smooth(data=Larch, aes(x=Year.CE, y=Tanon), method=gam ,  formula = y ~ s(x, k = 3), size = 1, color="red", fill="#FF9999", se=TRUE) +
  geom_point(data=Left, aes(x=Year.CE, y=Tanon, colour="Left"), shape=17) +
  geom_smooth(data=Left, aes(x=Year.CE, y=Tanon), method=gam ,  formula = y ~ s(x, k = 3), size = 1 ,color="darkgreen", fill="#69b3a2", se=TRUE) +
  geom_point(data=metstation, aes(x=Year, y=Tanon, colour="MET")) +
  geom_smooth(data=metstation, aes(x=Year, y=Tanon), method="gam", formula = y ~ s(x, k = 3), size = 1 , color="blue", fill="#99CCFF", se=TRUE) +
  labs(colour="Legend")

p6 + theme(legend.position = c(0.1, 0.8))

dev.off()
