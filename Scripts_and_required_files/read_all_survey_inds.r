library(strucchange);
source("C:/Users/JR13/Desktop/TyL/Scripts_and_required_files/Reflevels.r")
#more complicated below
#install.packages("rioja")
#install.packages("ade4")
source("C:/Users/JR13/Desktop/TyL/Scripts_and_required_files/REGIME2.r")


#[1] "BBICnSpaOT4" "BBICPorOT4"  "BBICsSpaOT1" "BBICsSpaOT4" "CSBBFraOT4"  "CSEngBT3"    "CSIreOT4"    "CSNIrOT1"   
#[9] "CSNIrOT4"    "CSScoOT1"    "CSScoOT4"   
#"GNSEngBT3"   "GNSFraOT4"   "GNSGerBT3"   "GNSIntOT1"   "GNSIntOT3"   "GNSNetBT3"   
#"WAScoOT3"    "WASpaOT3"  "CSFraOT4"

setwd(paste(PACKAGEDIRECTORY,"Outputs/", sep = ""))
##### now read GNS indicators all surveys
CS<-T#GNS
surveys <- c("GNSEngBT3",   "GNSFraOT4",   "GNSGerBT3",   "GNSIntOT1","GNSIntOT3",   "GNSNetBT3"   )
PELSURY<- c(2,4,5)
##### now read CS indicators all surveysCS<-T
if(CS) surveys <- c("BBICPorOT4", "CSBBFraOT4",
             "CSEngBT3","CSIreOT4","CSNIrOT1","CSNIrOT4","CSScoOT1", "CSScoOT4",
             "WAScoOT3","CSFraOT4")
if(CS) PELSURY<- -c(6)

DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep="")
sur<-surveys[1]
YRS <-  1983:2017
YRS_GNS<- 2002:2015 ## all surveys
if(CS) YRS_GNS<- 2002:2015 ## all surveys

ind<-NULL
###############################################################################
#Surv_biotonnes_Year

SP <- "DEM"
BIO<- matrix(NA,ncol=length(surveys),nrow=length(YRS))
rownames(BIO) <-YRS
colnames(BIO) <-surveys
for(sur in surveys){
  print(sur)
  ifelse(sur %in% c("CSFraOT4","CSBBFraOT4"), DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep=""),DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep=""))
  if(sur %in% c("CSBBFraOT4","CSScoOT1"))  DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep="")
  if(sur %in% c("CSFraOT4"))  DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep="")
  ind[[sur]] <- 
    read.csv( paste(sur,"/",sur,DTE,SP,"_Surv_biotonnes_Year.csv",sep='') )
  BIO[rownames(BIO) %in% ind[[sur]]$Year,which(colnames(BIO)==sur)] <- ind[[sur]]$CatCatchWgtSwept 
}
BIO<-data.frame(BIO)
cor(BIO,use="pairwise.complete.obs")

BIOl <- BIO
names(BIOl) <- paste("Biomass",names(BIOl))
#REFLVL(BIOl[,-3],H=0.3)

#plot
windows()
par(mfrow=c(2,2))
plot(YRS, BIO[,1]/1000,col="white",ylim=range(BIO,na.rm=T)/1000,main="Demersal fish", ylab="Biomass (kt)",xlab="Year")
for(sur in surveys) lines(YRS, BIO[,sur]/1000,col=which(sur == surveys))
plot(YRS, BIO[,1]/1000,col="white",ylim=range(BIO,na.rm=T)/1000,main="", ylab="",axes=F,xlab="")
legend("left",surveys,lty=1,col= 1:length(surveys))




if(CS) surveys <- c("BBICPorOT4", "CSBBFraOT4",
                    "CSIreOT4","CSNIrOT1","CSNIrOT4","CSScoOT1", "CSScoOT4",
                    "WAScoOT3","CSFraOT4")
if(CS) PELSURY<- -c(6)
BIO<- matrix(NA,ncol=length(surveys),nrow=length(YRS))
rownames(BIO) <-YRS
colnames(BIO) <-surveys

SP <- "PEL"
BIOp<- matrix(NA,ncol=length(surveys),nrow=length(YRS))
rownames(BIOp) <-YRS
colnames(BIOp) <-surveys

for(sur in surveys[PELSURY]){
  print(sur)
  ifelse(sur %in% c("CSFraOT4","CSBBFraOT4"), DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep=""),DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep=""))
  if(sur %in% c("CSBBFraOT4","CSScoOT1"))  DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep="")
  if(sur %in% c("CSFraOT4"))  DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep="")
  
  ind[[sur]] <- read.csv( paste(sur,"/",sur,DTE,SP,"_Surv_biotonnes_Year.csv",sep='') )
  BIOp[rownames(BIOp) %in% ind[[sur]]$Year,which(colnames(BIOp)==sur)] <- ind[[sur]]$CatCatchWgtSwept 
}
cor(BIOp,use="pairwise.complete.obs")
#plot pelagic par(mfrow=c(1,2))
plot(YRS, BIOp[,1]/1000,col="white",ylim=range(BIOp,na.rm=T)/1000,main="Pelagic fish", ylab="Biomass (kt)",xlab="Year")
for(sur in surveys) lines(YRS, BIOp[,sur]/1000,col=which(sur == surveys))
if(!CS) savePlot(filename = "GNSbiomass.bmp",type="bmp")
if(CS) savePlot(filename = "CSbiomass.bmp",type="bmp")

BIOpl<-data.frame(BIOp)
names(BIOpl) <- paste("BiomassPelagic",names(BIO))
#REFLVL(BIOpl[,PELSURY],H=0.3)

BAR<-apply(X=BIOp,FUN=mean,2,na.rm=T)
SD<-apply(X=BIOp,FUN=sd,2,na.rm=T)
BIOps <- scale(BIOp)
#plot(YRS, BIOps[,1],col="white",ylim=range(BIOps,na.rm=T),main="GNS Pelagic fish", ylab="Biomass (scaled)",xlab="Year")
#for(sur in surveys) lines(YRS, BIOps[,sur],col=which(sur == surveys))
###########################################################################
#MML

ind<-NULL
IND <- "MaxL"
SP <- "DEM"
MML<- matrix(NA,ncol=length(surveys),nrow=length(YRS))
rownames(MML) <-YRS
colnames(MML) <-surveys
for(sur in surveys){
  ifelse(sur %in% c("CSFraOT4","CSBBFraOT4"), DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep=""),DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep=""))
  if(sur %in% c("CSBBFraOT4","CSScoOT1"))  DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep="")
  if(sur %in% c("CSFraOT4"))  DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep="")
  
  ind[[sur]] <- read.csv( paste(sur,"/",sur,DTE,SP,"_",IND,"_swept_sea.csv",sep='') )
  MML[rownames(MML) %in% ind[[sur]]$Year,which(colnames(MML)==sur)] <- ind[[sur]]$sea
}
MML<-data.frame(MML)
cor(MML,use="pairwise.complete.obs")
MMLbio <- rowSums(BIO*MML,na.rm=T)/ rowSums(BIO,na.rm=T)
library(Hmisc)
rcorr(as.matrix(MML))

#by subdiv

for(sur in surveys){
if(!is.null(ind[[sur]] )){
  CORS<- (rcorr(as.matrix( ind[[sur]] )))
  #write.table(c(colnames(CORS$r)),paste(sur,"MMLd_CORS.csv",sep=''),sep=',',col.names=F,row.names=F)
  write.table("r",paste(sur,"MMLd_CORS.csv",sep=''),sep=',',col.names=F,row.names=F,append = F)
  write.table(CORS$r,paste(sur,"MMLd_CORS.csv",sep=''),sep=',',append = T,col.names=F)
  write.table("n",paste(sur,"MMLd_CORS.csv",sep=''),sep=',',append = T,col.names=F,row.names=F)
  write.table(CORS$n,paste(sur,"MMLd_CORS.csv",sep=''),sep=',',append = T,col.names=F)
  write.table("P",paste(sur,"MMLd_CORS.csv",sep=''),sep=',',append = T,col.names=F,row.names=F)
  write.table(CORS$P,paste(sur,"MMLd_CORS.csv",sep=''),append = T,sep=',',col.names=F)
}
}

#plot
windows()
par(mfrow=c(2,2))
plot(YRS, MML[,1],col="white",ylim=range(MML,na.rm=T),main="Demersal fish", ylab="MML (cm)",xlab="Year")
for(sur in surveys) lines(YRS, MML[,sur],col=which(sur == surveys))
#lines(YRS_GNS, MMLbio[which(names(MMLbio) %in% YRS_GNS)],col=1,lwd=4)
plot(YRS, MML[,1],col="white",ylim=range(MML,na.rm=T),main="", ylab="",axes=F,xlab="")
legend("left",surveys,lty=1,col= 1:length(surveys))

MMLl<-MML
names(MMLl) <- paste("MMLdemersal",names(MML))
#REFLVL(MMLl[,-3],H=0.3)
biodpr<-MMLl[paste(1999:2015),-3]
#REGIME(dat=prcomp(biodpr), tit='PCA MMLdemersal 99-15 step change 30pc',H=0.3,BREAKS=-11,LINEAR=0,PLOT_F=F,PLOT_ALT=F,COMP=c(1,2,3),ADD=F)

ind<-NULL
SP <- "PEL"

MMLp<- matrix(NA,ncol=length(surveys),nrow=length(YRS))
rownames(MMLp) <-YRS
colnames(MMLp) <-surveys
for(sur in surveys[PELSURY]){
  ifelse(sur %in% c("CSFraOT4","CSBBFraOT4"), DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep=""),DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep=""))
  if(sur %in% c("CSBBFraOT4","CSScoOT1"))  DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep="")
  if(sur %in% c("CSFraOT4"))  DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep="")
  
  ind[[sur]] <- read.csv( paste(sur,"/",sur,DTE,SP,"_",IND,"_swept_sea.csv",sep='') )
  MMLp[rownames(MMLp) %in% ind[[sur]]$Year,which(colnames(MMLp)==sur)] <- ind[[sur]]$sea
}
MMLbio_p<-rowSums(BIOp*MMLp,na.rm=T)/ rowSums(BIOp,na.rm=T)

cor(MMLp,use="pairwise.complete.obs")
rcorr(as.matrix(MMLp))##REFLVLs

#by subdiv
for(sur in surveys){
if(!is.null(ind[[sur]] )){
  CORS<- (rcorr(as.matrix( ind[[sur]] )))
  #write.table(c(colnames(CORS$r)),paste(sur,"MMLp_CORS.csv",sep=''),sep=',',col.names=F,row.names=F)
  write.table("r",paste(sur,"MMLp_CORS.csv",sep=''),sep=',',col.names=F,row.names=F,append = F)
  write.table(CORS$r,paste(sur,"MMLp_CORS.csv",sep=''),sep=',',append = T,col.names=F)
  write.table("n",paste(sur,"MMLp_CORS.csv",sep=''),sep=',',append = T,col.names=F,row.names=F)
  write.table(CORS$n,paste(sur,"MMLp_CORS.csv",sep=''),sep=',',append = T,col.names=F)
  write.table("P",paste(sur,"MMLp_CORS.csv",sep=''),sep=',',append = T,col.names=F,row.names=F)
  write.table(CORS$P,paste(sur,"MMLp_CORS.csv",sep=''),append = T,sep=',',col.names=F)
}
}
MMLpl<-data.frame(MMLp)
names(MMLpl) <- paste("MMLpelagic",names(MMLpl))
#REFLVL(MMLpl[,c(2,4,5)],H=0.3)
biodpr<-MMLpl[paste(1999:2015),c(2,4,5)]
#REGIME(dat=prcomp(biodpr), tit='PCA MMLpelagic 99-15 step change 30pc',H=0.3,BREAKS=-11,LINEAR=0,PLOT_F=F,PLOT_ALT=F,COMP=c(1,2,3),ADD=F)

#plot pelagic par(mfrow=c(1,2))
plot(YRS, MMLp[,1],col="white",ylim=range(MMLp,na.rm=T),main="Pelagic fish", ylab="MMLp (cm)",xlab="Year")
for(sur in surveys) lines(YRS, MMLp[,sur],col=which(sur == surveys))
#lines(YRS_GNS, MMLbio_p[which(names(MMLbio_p) %in% YRS_GNS)],col=1,lwd=4)
if(CS) savePlot(filename = "CS_MML.bmp",type="bmp")
if(!CS) savePlot(filename = "GNS_MML.bmp",type="bmp")



###########################################################################
#TyL

ind<-NULL
IND <- "TyL"
SP <- "DEM"
TyL<- matrix(NA,ncol=length(surveys),nrow=length(YRS))
rownames(TyL) <-YRS
colnames(TyL) <-surveys
for(sur in surveys){
  ifelse(sur %in% c("CSFraOT4","CSBBFraOT4"), DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep=""),DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep=""))
  if(sur %in% c("CSBBFraOT4","CSScoOT1"))  DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep="")
  if(sur %in% c("CSFraOT4"))  DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep="")
  
  ind[[sur]] <- read.csv( paste(sur,"/",sur,DTE,SP,"_",IND,".cm_swept_subdiv_sea.csv",sep='') )
  TyL[rownames(TyL) %in% ind[[sur]]$Year,which(colnames(TyL)==sur)] <- ind[[sur]]$sea
  TyL[rownames(TyL) %in% ind[[sur]]$X,which(colnames(TyL)==sur)] <- ind[[sur]]$sea
}
TyL<-data.frame(TyL)
cor(TyL,use="pairwise.complete.obs")
rcorr(as.matrix(TyL))

#by subdiv

for(sur in surveys){
if(!is.null(ind[[sur]] )){
  CORS<- (rcorr(as.matrix( ind[[sur]] )))
  #write.table(c(colnames(CORS$r)),paste(sur,"TyL_CORS.csv",sep=''),sep=',',col.names=F,row.names=F)
  write.table("r",paste(sur,"TyL_CORS.csv",sep=''),sep=',',col.names=F,row.names=F,append = F)
  write.table(CORS$r,paste(sur,"TyL_CORS.csv",sep=''),sep=',',append = T,col.names=F)
  write.table("n",paste(sur,"TyL_CORS.csv",sep=''),sep=',',append = T,col.names=F,row.names=F)
  write.table(CORS$n,paste(sur,"TyL_CORS.csv",sep=''),sep=',',append = T,col.names=F)
  write.table("P",paste(sur,"TyL_CORS.csv",sep=''),sep=',',append = T,col.names=F,row.names=F)
  write.table(CORS$P,paste(sur,"TyL_CORS.csv",sep=''),append = T,sep=',',col.names=F)
}
}
TyLbio <- rowSums(BIO*TyL,na.rm=T)/ rowSums(BIO,na.rm=T)

#plot
windows()
par(mfrow=c(2,2))
plot(YRS, TyL[,1],col="white",ylim=range(TyL,na.rm=T),main="Demersal fish", ylab="TyL (cm)",xlab="Year")
for(sur in surveys) lines(YRS, TyL[,sur],col=which(sur == surveys))
lines(YRS_GNS, TyLbio[which(names(TyLbio) %in% YRS_GNS)],col=1,lwd=4)
plot(YRS, TyL[,1],col="white",ylim=range(TyL,na.rm=T),main="", ylab="",axes=F,xlab="")
legend("left",surveys,lty=1,col= 1:length(surveys))

SP <- "PEL"
TyLp<- matrix(NA,ncol=length(surveys),nrow=length(YRS))
rownames(TyLp) <-YRS
colnames(TyLp) <-surveys
for(sur in surveys[PELSURY]){
  ifelse(sur %in% c("CSFraOT4","CSBBFraOT4"), DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep=""),DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep=""))
  if(sur %in% c("CSBBFraOT4","CSScoOT1"))  DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep="")
  if(sur %in% c("CSFraOT4"))  DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep="")
  ind[[sur]] <- read.csv( paste(sur,"/",sur,DTE,SP,"_",IND,".cm_swept_subdiv_sea.csv",sep='') )
  TyLp[rownames(TyLp) %in% ind[[sur]]$Year,which(colnames(TyLp)==sur)] <- ind[[sur]]$sea
  TyLp[rownames(TyLp) %in% ind[[sur]]$X,which(colnames(TyLp)==sur)] <- ind[[sur]]$sea
}
TyLbiop <-rowSums(BIOp*TyLp,na.rm=T)/ rowSums(BIOp,na.rm=T) 

cor(TyLp,use="pairwise.complete.obs")
rcorr(as.matrix(TyLp))

#by subdiv

for(sur in surveys){
if(!is.null(ind[[sur]] )){
  CORS<- (rcorr(as.matrix( ind[[sur]] )))
  #write.table(c(colnames(CORS$r)),paste(sur,"TyLp_CORS.csv",sep=''),sep=',',col.names=F,row.names=F)
  write.table("r",paste(sur,"TyLp_CORS.csv",sep=''),sep=',',col.names=F,row.names=F,append = F)
  write.table(CORS$r,paste(sur,"TyLp_CORS.csv",sep=''),sep=',',append = T,col.names=F)
  write.table("n",paste(sur,"TyLp_CORS.csv",sep=''),sep=',',append = T,col.names=F,row.names=F)
  write.table(CORS$n,paste(sur,"TyLp_CORS.csv",sep=''),sep=',',append = T,col.names=F)
  write.table("P",paste(sur,"TyLp_CORS.csv",sep=''),sep=',',append = T,col.names=F,row.names=F)
  write.table(CORS$P,paste(sur,"TyLp_CORS.csv",sep=''),append = T,sep=',',col.names=F)
}
}

#plot pelagic par(mfrow=c(1,2))
plot(YRS, TyLp[,1],col="white",ylim=range(TyLp,na.rm=T),main="GNS Pelagic fish", ylab="TyLp (cm)",xlab="Year")
for(sur in surveys) lines(YRS, TyLp[,sur],col=which(sur == surveys))
lines(YRS_GNS, TyLbiop[which(names(TyLbiop) %in% YRS_GNS)],col=1,lwd=4)
if(!CS) savePlot(filename = "GNS_TyL.bmp",type="bmp")
if(CS) savePlot(filename = "CS_TyL.bmp",type="bmp")


##REFLVLs

TyLl<-TyL
names(TyLl) <- paste("TyLdemersal",names(TyL))
#REFLVL(TyLl[,-3],H=0.3)
biodpr<-TyLl[paste(1999:2015),-3]
#REGIME(dat=prcomp(biodpr), tit='PCA TyLdemersal 99-15 step change 30pc',H=0.3,BREAKS=-11,LINEAR=0,PLOT_F=F,PLOT_ALT=F,COMP=c(1,2,3),ADD=F)

TyLpl<-data.frame(TyLp)
names(TyLpl) <- paste("TyLpelagic",names(TyLpl))
#REFLVL(TyLpl[,PELSURV],H=0.3)
biodpr<-TyLpl[paste(1999:2015),c(2,4,5)]
#REGIME(dat=prcomp(biodpr), tit='PCA TyLpelagic 99-15 step change 30pc',H=0.3,BREAKS=-11,LINEAR=0,PLOT_F=F,PLOT_ALT=F,COMP=c(1,2,3),ADD=F)

###############################################################################

##### weighted average
MMLbio <- rowSums(BIO*MML,na.rm=T)/ rowSums(BIO,na.rm=T)
MMLbio_p<-rowSums(BIOp*MMLp,na.rm=T)/ rowSums(BIOp,na.rm=T)
TyLbio <- rowSums(BIO*TyL,na.rm=T)/ rowSums(BIO,na.rm=T)
TyLbiop <-rowSums(BIOp*TyLp,na.rm=T)/ rowSums(BIOp,na.rm=T) 
  
windows()
par(mfrow=c(2,2))
YRS_GNS<- 2002:2015
plot(YRS_GNS, MMLbio[which(names(MMLbio) %in% YRS_GNS)],xlab="Year",ylab="MML (cm)",main="demersal fish",type='l',lwd=2)
plot(YRS_GNS, MMLbio_p[which(names(MMLbio_p) %in% YRS_GNS)],xlab="Year",ylab="MML (cm)",main="pelagic fish",type='l',lwd=2)
plot(YRS_GNS, TyLbio[which(names(TyLbio) %in% YRS_GNS)],xlab="Year",ylab="TyL (cm)",main="demersal fish",type='l',lwd=2)
plot(YRS_GNS, TyLbiop[which(names(TyLbiop) %in% YRS_GNS)],xlab="Year",ylab="TyL (cm)",main="pelagic fish",type='l',lwd=2)
if(!CS) savePlot(filename = "GNS_TyLMML_weighted.bmp",type="bmp")
if(CS) savePlot(filename = "CS_TyLMML_weighted.bmp",type="bmp")






######
SciName2GroupEff <- read.csv("Z:/Eco Indicators/DATRASoutput/MarScot/INDscripts/SciName2GroupEff.csv")
# add catchability by length group
  Q <- read.csv("Z:/Eco Indicators/DATRASoutput/Walker_GearEfficiency_ICESJMarSci17_Supp/EfficiencyTab.csv")
  Q$Code <- as.character(Q$Code)
  Qsp <- read.csv("Z:/Eco Indicators/DATRASoutput/Walker_GearEfficiency_ICESJMarSci17_Supp/specieslist.csv");names(Qsp)[1] <- "SciName"
  Qsp$Code <- as.character(Qsp$SpCode)
  Qs <- merge(x=Q,y=Qsp,by.x="Code", by.y="Code",all=T); 
  Qs$SciName <- as.character(Qs$SciName)#spaces at the end!
  Qs$SciName <- substr(Qs$SciName, start=1,stop=(nchar(Qs$SciName))-1)
  rm(Q,Qsp)
  ##apply to surveys
###### SSS
#logbio vs logL
# CSBBFraOT4_09Jan2017.LD_tonnes_Year_W.by.subdiv
sur <- "GNSIntOT3"
ifelse(sur %in% c("CSFraOT4","CSBBFraOT4"), DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep=""),DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep=""))
indsurdat0 <- read.csv( paste(sur,"/",sur,DTE,".LD_tonnes_Year_W.by.subdiv.csv",sep='') )
indsurdat <- merge(indsurdat0, Qs, by.x="SpeciesSciName",by.y="SciName",all.x=TRUE)
# indsurdat <- indsurdat0[indsurdat0$DEMPEL=="PEL",]
# indsurdat <- indsurdat0[indsurdat0$DEMPEL=="DEM",]

indsurdat$CatCatchWgtSwept_Log <- log(indsurdat$CatCatchWgtSwept)
indsurdat$FishLength_Log <- log(indsurdat$FishLength_cm)
indsurdat$FishCat_Log <- indsurdat$FishLength_Log*indsurdat$CatCatchWgtSwept

#by year
indsur<- aggregate(indsurdat[,which(names(indsurdat) %in% c("FishCat_Log","CatCatchWgtSwept"))],
                   by=list(indsurdat$FishLength_cm,indsurdat$Year), FUN=sum)
  indsur$CatCatchWgtSwept_Log <- log(indsur$CatCatchWgtSwept)
  indsur$FishLength_Log <- log(indsur$Group.1)
  
  with(indsur,xyplot(CatCatchWgtSwept_Log ~ FishLength_Log  | as.factor(Group.2)))
  LMIN<-30;  LMAX<-80
  with(indsur[indsur$Group.1>LMIN & indsur$Group.1<LMAX,],xyplot(CatCatchWgtSwept_Log ~ FishLength_Log | 
                                                              as.factor(Group.2)) )
  sss<-NULL
  for( SURV in as.numeric(indsur$Group.2)){
  sss[paste(SURV)] <- summary(lm(CatCatchWgtSwept_Log ~ FishLength_Log, 
                          data=indsur[as.numeric(indsur$Group.2)==SURV & indsur$Group.1>LMIN & indsur$Group.1<LMAX,]))$coef[2]
  }
  
  #cf sss tyl
  indsur2<- aggregate(indsur[,which(names(indsur) %in% c("FishCat_Log","CatCatchWgtSwept"))],
                     by=list(indsur$Group.2), FUN=sum)
  indsur2$TyL <- exp( indsur2$FishCat_Log / indsur2$CatCatchWgtSwept )
  par(mfrow=c(2,1)); plot(indsur2$TyL)
  plot(sss)
  summary(lm(indsur2$TyL ~ sss))
  
  
#by area
  indsur<- aggregate(indsurdat[,which(names(indsurdat) %in% c("FishCat_Log","CatCatchWgtSwept"))],
                     by=list(indsurdat$FishLength_cm,indsurdat$SurvStratum), FUN=sum)
  indsur$CatCatchWgtSwept_Log <- log(indsur$CatCatchWgtSwept)
  indsur$FishLength_Log <- log(indsur$Group.1)
  
  with(indsur,xyplot(CatCatchWgtSwept_Log ~ FishLength_Log  | Group.2))
  
  with(indsur[indsur$Group.1>20 & indsur$Group.1<80,],xyplot(CatCatchWgtSwept_Log ~ FishLength_Log | 
                                                               as.factor(Group.2)) )
  sss<-NULL
  for( SURV in indsur$Group.2){
    sss[SURV] <- summary(lm(CatCatchWgtSwept_Log ~ FishLength_Log, 
                            data=indsur[indsur$Group.2==SURV & indsur$Group.1>20 & indsur$Group.1<80,]))$coef[2]
  }
  
  #cf sss tyl
  indsur2<- aggregate(indsur[,which(names(indsur) %in% c("FishCat_Log","CatCatchWgtSwept"))],
                      by=list(indsur$Group.2), FUN=sum)
  indsur2$TyL <- exp( indsur2$FishCat_Log / indsur2$CatCatchWgtSwept )
  plot(indsur2$TyL , sss)
  summary(lm(indsur2$TyL ~ sss))
  
  
#by year and area
indsur<- aggregate(indsurdat[,which(names(indsurdat) %in% c("FishCat_Log","CatCatchWgtSwept"))],
                   by=list(indsurdat$FishLength_cm,indsurdat$SurvStratum,indsurdat$Year), FUN=sum)
indsur$CatCatchWgtSwept_Log <- log(indsur$CatCatchWgtSwept)
indsur$FishLength_Log <- log(indsur$Group.1)

with(indsur,xyplot(CatCatchWgtSwept_Log ~ FishLength_Log  | Group.2))
with(indsur,xyplot(CatCatchWgtSwept_Log ~ Group.1 | as.factor(Group.2),groups= Group.3) )
with(indsur[indsur$Group.1>12 & indsur$Group.1<100,],xyplot(CatCatchWgtSwept_Log ~ FishLength_Log | 
                                                              as.factor(Group.2),groups= Group.3) )


#by year and area
indsuragg<- aggregate(indsur[,which(names(indsur) %in% c("FishCat_Log","CatCatchWgtSwept"))],
                      by=list(indsur$Group.3,indsur$Group.2), FUN=sum)
indsuragg$TyL <- exp( indsuragg$FishCat_Log / indsuragg$CatCatchWgtSwept )
indsuragg$CatCatchWgtSwept_Log <-log(indsuragg$CatCatchWgtSwept)
indsuragg$LogLen <- indsuragg$FishCat_Log/indsuragg$CatCatchWgtSwept

xyplot(indsuragg$CatCatchWgtSwept_Log ~ indsuragg$LogLen | as.factor(indsuragg$Group.2), pch=19, groups=indsuragg$Group.2 )
x11();xyplot(indsuragg$CatCatchWgtSwept_Log ~ indsuragg$LogLen, groups=indsuragg$Group.2, pch=19)


indsuraggsea<- aggregate(indsur[,which(names(indsur) %in% c("FishCat_Log","CatCatchWgtSwept"))],
                         by=list(indsur$Group.3), FUN=sum)
indsuraggsea$TyL <- exp( indsuraggsea$FishCat_Log / indsuraggsea$CatCatchWgtSwept )

#plot(indsuragg[indsuragg$Group.2=="Gn1",]$Group.1, indsuragg[indsuragg$Group.2=="Gn1",]$TyL)
plot(indsuragg[indsuragg$Group.2=="KS",]$Group.1, indsuragg[indsuragg$Group.2=="KS",]$TyL)
plot(indsuragg[indsuragg$Group.2=="NE",]$Group.1, indsuragg[indsuragg$Group.2=="NE",]$TyL)
plot(indsuragg[indsuragg$Group.2=="CW",]$Group.1, indsuragg[indsuragg$Group.2=="CW",]$TyL)
plot(indsuragg[indsuragg$Group.2=="SW",]$Group.1, indsuragg[indsuragg$Group.2=="SW",]$TyL)
plot(indsuragg[indsuragg$Group.2=="SE",]$Group.1, indsuragg[indsuragg$Group.2=="SE",]$TyL)

