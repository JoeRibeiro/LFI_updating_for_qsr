#  Before using, please install used packages


# with option of bootstapping by haul and subdiv and trendline fitting
#
# Author: C Lynam, Cefas
# Contact: chris.lynam@cefas.co.uk
# Version: 2.3 
# Date: 10 Feb 2017
# edits: now choose surveys to loop through upfront, set paths used throughout upfront, uses new dataproduct structure
# reads shapefiles for subdivisions and includes spatial area of strata in calcs
# if both SAMP_STRAT and BYSUBDIV are true now averages catch by SAMP_STRAT then multiples by AREA; then sums within SUBDIV and raises to whole SUBDIV AREA if not fully sampled, so combination of SUBDIV not subject to change in relative contribution of SUBDIV due to sampling 
# link to loess ploting code for subdivs complete

# to do -
# link to shapefile survey maps colours by status
# alter bootstrap code so works with survey strata not just rectangles
#  Error in MINUS$centlon : $ operator is invalid for atomic vectors



#TODAYSDATE = ""

LFIthresholds <- data.frame(BBICnSpaOT4=25,WASpaOT3=40,BBICsSpaOT4=45,BBICPorOT4=40,CSBBFraOT4=35,CSEngBT3=35,CSIreOT4=30,CSNIrOT1=45,CSNIrOT4=40,CSScoOT1=35,CSScoOT4=40,GNSEngBT3=50,GNSFraOT4=50,GNSGerBT3=30,GNSIntOT1=50,GNSIntOT3=50,GNSNetBT3=30,WAScoOT3=50,CSFraOT4=40)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Automatically (or manually, if desired...) set working directory
#PACKAGEDIRECTORY = "Y:/MA016FI_Cefmat/Working_Area/D4/TyL/" # If manually changing directory here, then comment out the second reference to PACKAGEDIRECTORY
library(rstudioapi)    
this_script_dir = rstudioapi::getActiveDocumentContext()$path # DO NOT COPY AND PASTE THIS LINE INTO THE CONSOLE!! THIS LINE ONLY WORKS WHEN EXECUTED FROM THE SCRIPT
this_script_dir2 = strsplit(this_script_dir,  "Scripts_and_required_files"); this_script_dir2 = this_script_dir2[[1]]
PACKAGEDIRECTORY = this_script_dir2[1]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#location of the data and subscripts
MAINDIR<- paste(PACKAGEDIRECTORY,"Scripts_and_required_files/Data_Product_V3/", sep = "")
#location of the subscripts
SUBSCRIPTS_TRAITS<- paste(PACKAGEDIRECTORY,"Scripts_and_required_files/", sep = "")
#location of the shapefiles and where is the 'attributes' folder
SHAPEPATH<- paste(PACKAGEDIRECTORY,"Scripts_and_required_files/Strata/", sep = "")
#where save output?
OUTPATHstem<- paste(PACKAGEDIRECTORY,"Outputs/", sep = "")



## choices for analyses upfront
setwd(MAINDIR)
SURVEYS<-list.files("Sampling_Info\\")
SURVEYS<-strsplit(SURVEYS,split="_")
SURVEY<-NULL
for(i in 1:length(SURVEYS)) SURVEY[i] <- SURVEYS[[i]][2]
SURVEY<-SURVEY[SURVEY!="AllSurveys" & !is.na(SURVEY)]; rm(SURVEYS)
SURVEY[length(SURVEY)+1]<-"CSFraOT4"

#SURVEY<- c("GNSEngBT3", "GNSFraOT4",  "GNSGerBT3",   "GNSIntOT1",   "GNSIntOT3",  "GNSNetBT3")
EHDS_PP<-F; 
#SURVEY<- c("GNSIntOT1")
#SURVEY<- c("GNSIntOT3","GNSGerBT3","GNSIntOT1","GNSNetBT3")
#SURVEY<-SURVEY[!SURVEY%in%c("GNSIntOT3","GNSGerBT3","GNSIntOT1","GNSNetBT3")]
#full list
# "BBICnSpaOT4" "BBICPorOT4"  "BBICsSpaOT1" "BBICsSpaOT4" "CSBBFraOT4"  "CSEngBT3"    "CSIreOT4"    "CSNIrOT1"    "CSNIrOT4"    "CSScoOT1"   
# "CSScoOT4"    "GNSEngBT3"   "GNSFraOT4"   "GNSGerBT3"   "GNSIntOT1"   "GNSIntOT3"   "GNSNetBT3"   "WAScoOT3"    "WASpaOT3"    "CSFraOT4"  

#do you want to write outputs along the way and save workspace
WRITE <- T #save csvs as we go?
WRITE_LDs <- T #write Length distributions by species/year/subdiv? 
BOOTSTRAP <- F # invoke code? if F next 3 lines redundant
NBOOT<-200
WRITE_BOOT <- F # every bootstrap dataset and indicator output
SAVE <- T # save workspace (after bootstrap)
#which indicators?
TyL_GeoM <- F # OSPAR FW3
MaxL <- F # OSPAR FC3
MEANTLs <- F #similar to OSPAR FW4 # not will not be calc'd for WAsurveys as no data file for TL
MeanL <- F #not OSPAR but simple
LFI <- T # similar to OSPAR FC2 
OVERWITE_SUBDIV<-F ##not run##01Feb2017# if true ignore subdivisions [as done for LFI] and use sampstrat as subdivisions instead
# catchability for general species groups
SPECIES_IN_MOD_ONLY <-F
CATCHABILITY_COR_WALKER<-F
CATCHABILITY_COR_MOD<-F
SPECIES_IN_MOD_ONLY<-F
survey <-SURVEY
RAISE_Q<-F #can override below if Either CATCHABILITY_COR_ TRUE
#change lfi thr to 50
##NOTE additional choices below for plotting
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#useful packages
#install.packages("spatstat")
#install.packages("rgeos")
library(Hmisc)# for errbars
library(mgcv) # to investigate indicator trends with gam
library(maps) # for eyeball of data
library(lattice)
require(rgdal) # to read shapefiles for subdivisions
#to add shapefile info to data in Lynam_OSPARsubdiv_rectmod.r
library(sp)
library(spatstat)
library(rgeos)    
library(maptools)
library(stringr)

#increase memory to max
memory.limit(size =4095)
memory.limit(size = NA)

#some shorthand
an <- as.numeric
af <- as.factor
ac <- as.character
h <- head
u <- unique
count<- function(DATA){length(unique(DATA))}
neg<-function(x) -x 

#where are the subscripts and traits?
setwd(SUBSCRIPTS_TRAITS)
#additional functions
source("required.funcs.r")              # using tapply.ID
source("Lynam_INDfn_Dec2016 - OSPARrectmod.r") # Nov2016 update to output total of meanCPUE by strata and total survey area
if(BOOTSTRAP) source("Lynam_IND_BOOTfn_July_2016 - OSPAR.r")  # bootstrap the hauls by STSQ and subdiv
source("Lynam_INDPLOTFN.r")             # plot options 
#subscripts for indicators 
if(LFI) source("Lynam_INDfn_Dec2016_LFI_rectmod_May2018.R")   # Large Fish Indicator # Nov2016 no longer fall over if no fish above LFI_THRESHOLD
if(TyL_GeoM) source("Lynam_INDfn_Dec2016_TyL_GeoM.r")#Geometric mean length (Typical Length cm)
if(MeanL) source("Lynam_INDfn_Dec2016_MeanL.r") # Mean Length cm
if(MaxL) source("Lynam_INDfn_Dec2016_MaxL.r")  # MaxL
if(MEANTLs) source("Lynam_INDfn_Dec2016_MeanTL.r")# TL output by rectangle and year 
# info on trophic level, demersal/pelagic,  maxLength, etc
#trait_file <- paste("fishspecieslifehistorytraits_v3.csv",sep='')
# should have some similar file in final dataproduct!
# dataset with traits by species for indicator script
#trait <- read.csv(trait_file)
#trait$SpeciesSciName <- paste(trait$Genus,trait$Species,sep='')

#max length observed in ospar dataproduct by species
#Hyperoplus immaculatus 'Greater sand-eel' were missing -> recorded in GNS IBTS Q1 but without length  treat as Ammodytidae 'sandeel'
MAXL_file <- paste("traits_by_species_Mar2018.csv",sep='')
trait_MAXL <- read.csv(MAXL_file)
#trophic level info
#cnrs mtl fw4
# TLceltic <- read.csv("TL_complete14102016-1.csv")
# names(TLceltic)[1] <- "SpeciesSciName"
# TLceltic <- TLceltic[,c("SpeciesSciName","TL")]

# west of scotland
#TLwest <- NULL

# north sea
# TLnorth <- read.csv("refspp_LWmerged_edit06Nov2015.csv")
# names(TLnorth)[2] <- "SpeciesSciName"
# ECOPATH_TL <- F; if(ECOPATH_TL)  TLnorth$TL <- TLnorth$TL_Ecopath #for comparison with ewe
# TLnorth <- TLnorth[,c("SpeciesSciName","TL")]

# catchability for general species groups
if(CATCHABILITY_COR_WALKER){
  SciName2GroupEff <- read.csv("SciName2GroupEff.csv")
  # add catchability by length group
  Q <- read.csv(paste(PACKAGEDIRECTORY,"Scripts_and_required_files/Walker_GearEfficiency_ICESJMarSci17_Supp/EfficiencyTab.csv",sep=""))
  Q$Code <- ac(Q$Code)
  Qsp <- read.csv(paste(PACKAGEDIRECTORY,"Scripts_and_required_files/Walker_GearEfficiency_ICESJMarSci17_Supp//specieslist.csv",sep=""));names(Qsp)[1] <- "SciName"
  Qsp$Code <- ac(Qsp$Code)
  Qs <- merge(x=Q,y=Qsp,by.x="Code", by.y="Code",all=T); 
  Qs$SciName <- ac(Qs$SciName)#spaces at the end!
  Qs$SciName <- substr(Qs$SciName, start=1,stop=(nchar(Qs$SciName))-1)
  rm(Q,Qsp)
  ##apply to GNS surveys
  RAISE_Q<-T
}
if(CATCHABILITY_COR_MOD | SPECIES_IN_MOD_ONLY){
  # add catchability by length group
  if(survey == "CSFraOT4"){
    #Q <- read.table("O:/C5716_FizzyFish/Working_Area/WP3 Food web indicators/lenmodels/CSeaModSpecies.csv",header=T,sep=',',as.is=T);         SPECIES_SUBSET <- Q[,1]
    Q <- read.table("//lowfilecds/CDP/C5716_FizzyFish/Working_Area/WP3 Food web indicators/lenmodels/CSea RATIO SURVEY_MODEL.csv",header=T,sep=',',as.is=T)        
    
    SPECIES_SUBSET <- names(Q[1,-c(1,2)])
    SPECIES_SUBSET <-sub("[.]", " ", SPECIES_SUBSET) #remove the R created .
    Q[1,-c(1,2)] <- SPECIES_SUBSET
  }
  if(survey == "GNSIntOT1"){
    Q <- read.table("//lowfilecds/CDP/C5716_FizzyFish/Working_Area/WP3 Food web indicators/lenmodels/NSea RATIO SURVEY_MODEL.csv",header=T,sep=',',as.is=T)
    SPECIES_SUBSET <- names(Q[1,-c(1,2)])
    SPECIES_SUBSET <-sub("[.]", " ", SPECIES_SUBSET) #remove the R created .
    Q[1,-c(1,2)] <- SPECIES_SUBSET
  }
  ##apply to GNS surveys
  #RAISE_Q<-F
}
#### indicators survey loop ####
#"BBICnSpaOT4" "BBICPorOT4"  "BBICsSpaOT1" "BBICsSpaOT4" "CSBBFraOT4"
#"CSEngBT3"    "CSIreOT4"    "CSNIrOT1"  "CSNIrOT4"    "CSScoOT1"    "CSScoOT4"    
#"GNSEngBT3"   "GNSFraOT4"  "GNSGerBT3"   "GNSIntOT1"   "GNSIntOT3"  "GNSNetBT3"   
#"WAScoOT3"    "WASpaOT3"  "CSFraOT4"
setwd(MAINDIR)
#SURVEY=SURVEY[10:14]
#for(survey in  c("GNSFraOT4")){ # survey <- "GNSGerBT3" #SURVEY[1]#"WAScoOT3" # #"CSFraOT4"#WASpaOT3#[4]# BBICsSpaOT1 and BBICsSpaOT4 missing datayear causes warnings in LFI fn so have suppressed these
for(survey in  SURVEY){ # survey <- "GNSGerBT3" #SURVEY[1]#"WAScoOT3" # #"CSFraOT4"#WASpaOT3#[4]# BBICsSpaOT1 and BBICsSpaOT4 missing datayear causes warnings in LFI fn so have suppressed these
  print(survey)
  #add a directory for files out
  if(!file.exists(paste(OUTPATHstem,survey,sep=''))) dir.create(file.path(OUTPATHstem, survey))
  OUTPATH <- paste(OUTPATHstem,survey,"/",sep='')
  LFI_THRESHOLD <- LFIthresholds[,survey]
  print(LFIthresholds[,survey])
  SPECIES <- c("DEM","PEL"); #
  if(survey %in% c("GNSGerBT3","GNSNetBT3","GNSEngBT3","GNSEngBT4","CSEngBT3","CSEngBT4")){ SPECIES <- "DEM"  } else {  SPECIES <- c("DEM","PEL")} #
  
  
  
  if(survey %in% c("GNSGerBT3","GNSNetBT3","GNSIntOT1","GNSIntOT3")){
    SAMP_STRAT <- T # average hauls by rect in north sea #if set to FALSE need to update Attibutes table as area only given by rect
    BYSUBDIV <-   T # average indicator by LFI-subdivision 
  } else { #CS, BB and WA
    SAMP_STRAT <- F # do not average hauls by rect
    BYSUBDIV <-   T # average indicator by survey strata from shapefiles
  }
  #exceptions
  if(survey=="GNSFraOT4"){
    SAMP_STRAT <- T # average hauls over surveyed minigrid0.25deg  #if set to FALSE need to update Attibutes table as area only given by minigrid
    BYSUBDIV <-   T # no higher survey strata - could use depth strata from GNSEngBT3 - not implemented
  }
  
  #Trophic Level FW4 check
  if(MEANTLs==T & !(substr(survey,1,2) %in% c("CS","BB", "GN")) ) { MEANTLs<-F; print("no data for TL for WA surveys, deleting MTL selection") }
  
  ##load data
  #sampling data 
  if(survey=="CSFraOT4"){
    samp_file <- paste("Sampling_Info\\SamplingInfo_CSBBFraOT4_SSASMP.csv",sep='')
  } else {
    samp_file <- paste("Sampling_Info\\SamplingInfo_",survey,"_SSASMP.csv",sep='')
  }
  samp <- read.csv(samp_file)
  
  # biological data
  if(survey=="CSFraOT4"){
    bio_file <- paste("Biological_Info\\BiologicalInfo_CSBBFraOT4_SSASMP_kNN.csv",sep='')
  } else {
    bio_file <- paste("Biological_Info\\BiologicalInfo_",survey,"_SSASMP_kNN.csv",sep='')
  }
  bio <- read.csv(bio_file)
  #names(bio)[which(names(bio) =="IndivFishWght_g")] <-"IndivFishWgt_g" #typo correct
  #names(bio)[which(names(bio) =="DensAbund_N_sqkm")] <-"DensAbund_N_Sqkm" #typo correct
  
  # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#quick eyeball of sampling data 
  # # quick check
  # #par(mfrow=c(1,1))
  # x11(); par(mfrow=c(4,4))
  # plot(samp$ShootLong_degdec, samp$ShootLat_degdec); map(add=T)
  # abline(v=c(4,8));  abline(h=55.5)
  # with(samp, hist(MonthShot))
  # with(samp, hist(HaulDur_min))
  # with(samp, hist(Depth_m))
  # with(samp, hist(WingSpread_m))
  # with(samp, hist(DoorSpread_m))
  # with(samp, hist(NetOpen_m))
  # with(samp, hist(WingSwpArea_sqkm) )
  # with(samp, hist(WingSwpVol_CorF ))
  # with(samp, hist(DoorSwptArea_CorF) )
  # with(samp, hist(DoorSwptVol_CorF))
  # with(samp, hist(Distance_km))
  # samp[(samp$Distance_km)==max(samp$Distance_km),]
  # #samp[(samp$Distance_km)>5,] # big value
  # ##x11(); par(mfrow=c(2,3))
  # with(bio, hist(FishLength_cm))
  # with(bio, hist(IndivFishWgt_g))
  # #bio[!is.na(bio$IndivFishWgt_g) & (bio$IndivFishWgt_g)==max(bio$IndivFishWgt_g,na.rm=T),]
  # with(bio, hist(DensAbund_N_Sqkm))
  # #bio[!is.na(bio$DensAbund_N_Sqkm) & (bio$DensAbund_N_Sqkm)==max(bio$DensAbund_N_Sqkm,na.rm=T),]
  # # samp[samp$HaulID == bio[(bio$DensAbund_N_Sqkm)==max(bio$DensAbund_N_Sqkm),'HaulID'], ] # sometimes get multiple matches e.g sprat on north sea IBTS Q1
  # with(bio, hist(DensBiom_kg_Sqkm))
  # #bio[!is.na(bio$DensBiom_kg_Sqkm) & (bio$DensBiom_kg_Sqkm)==max(bio$DensBiom_kg_Sqkm,na.rm=T),]
  # #with(bio, hist(Number))
  # #eyeball end
  # savePlot(filename= paste(OUTPATH,survey,"_",format(Sys.time(), "%d%b%Y"),"eyeball",".bmp",sep=''),type="bmp")
  # dev.off()
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #now merge datasets on haul and species for indicator script
  #bio <- read.csv(bio_file) #reduce memory usage
  bio <- bio[,c("HaulID","SpeciesSciName","FishLength_cm","DensBiom_kg_Sqkm")]
  #add efficiency of E=GOV gear
  #Gear efficiency is the probability that fish in the path of a trawl will be caught and retained
  if(RAISE_Q & CATCHABILITY_COR_WALKER){  
    bio <- merge(bio, Qs[Qs$Gear=="GOV",],  by.x=c("SpeciesSciName","FishLength_cm"), 
                 by.y=c("SciName", "Length"), all.x=T, all.y=F)
    bio$mult <- 1/bio$Efficiency
    # zeros for 
    u(bio[is.infinite(bio$mult),'SpeciesSciName'])
    #GNSEngBT3: Agonus cataphractus    Arnoglossus laterna    Callionymus maculatus  .
    #           Microchirus variegatus Syngnathus rostellatus
    #if efficiency is zero use absolute estimate last resort
    bio[is.infinite(bio$mult),'mult'] <- 1/bio[is.infinite(bio$mult),'Absolute']
    
    #if no absolute estimate use group
    u(bio[is.na(bio$mult),'SpeciesSciName'])# check species names not missing due to spelling mistakes
    bio[is.na(bio$mult),'mult'] <- 1/bio[is.na(bio$mult),'Efficiency']
  }
  samp <-samp[,c("HaulID","YearShot","ShootLat_degdec","ShootLong_degdec","ICESStSq","SurvStratum")]
  dhspp <- merge(bio,samp,by="HaulID")
  rm(samp,bio)
  
  if(SPECIES_IN_MOD_ONLY) dhspp <- dhspp[dhspp$SpeciesSciName %in% SPECIES_SUBSET,] #use only those from thorpe model
  if(RAISE_Q & CATCHABILITY_COR_MOD){  
    
    dhspp$from <- dhspp$to <- NA
    #from # an(Q[-1,1])
    #to   # an(Q[-1,2])
    for(i in 2:nrow(Q)) dhspp[dhspp$FishLength_cm >= an(Q[i,1]) & dhspp$FishLength_cm < an(Q[i,2]),"from"] <- an(Q[i,1])
    for(i in 2:nrow(Q)) dhspp[dhspp$FishLength_cm >= an(Q[i,1]) & dhspp$FishLength_cm < an(Q[i,2]),"to"]   <- an(Q[i,2])
    if(nrow(dhspp[is.na(dhspp$from),])>0) print("some species len cat without scaling factors")
    
    #species names in colnames have dot in, sciname in first row does not so matches    SpeciesSciName
    Qss <-NULL
    for(i in 1:length(SPECIES_SUBSET) ){
      Qs<-data.frame(matrix(NA,length(Q[-1,1]),3))
      names(Qs)<-c("SpeciesSciName","from","Efficiency")
      Qs[,"from"] <- an(Q[-1,1]) #from
      Qs[,"SpeciesSciName"] <- SPECIES_SUBSET[i]
      Qs["Efficiency"] <- Q[-1,Q[1,]==SPECIES_SUBSET[i]]
      Qss <- rbind(Qss,Qs)
    }
    Qss$Efficiency <- an(Qss$Efficiency) 
    Qss$mult <- 1/(Qss$Efficiency) # this is a scaled multiplier from surveyed tonnes to model currency of grams so community metrics match those from model output
    Qss[is.infinite(Qss$mult),"mult"]<-0
    #Qss$Efficiency <- Qss$Efficiency * (10^8) # to avoid numerical issues #biomass units now in 100g
    dhspp <- merge(dhspp, Qss,  by.x=c("SpeciesSciName","from"), 
                   by.y=c("SpeciesSciName", "from"), all.x=T, all.y=F)
    dhspp$DensBiom_kg_Sqkm <- dhspp$DensBiom_kg_Sqkm  * dhspp$mult/1000 #thousand since multipler based on tonnes output from script
    
  }
  #(dhspp[ac(dhspp$SpeciesSciName) %in% ac(SPECIES_SUBSET),"SpeciesSciName"])
  
  #simplify names
  names(dhspp)[which(names(dhspp) =="YearShot")] <-"Year"
  #for tyl
  dhspp$LogLngtClass <- log(dhspp$FishLength_cm)
  dhspp$LogLngtBio <- dhspp$LogLngtClass * dhspp$DensBiom_kg_Sqkm
  dhspp<- dhspp[,-which(names(dhspp) %in% "SurvStratum" )]
    #### link subdivisional shapefiles to data ####
  # and read attributes of shapefiles create table ATTRIB with names SurvStratum & KM2_LAM
  source(paste(PACKAGEDIRECTORY,"Scripts_and_required_files/Lynam_OSPARsubdiv_rm2.r",sep=""))  

  #MaxLength and demersal or pelagic
  #MaxLength observed for MML
#  if(OVERWITE_SUBDIV) 
  dhspp <- merge(dhspp,trait_MAXL[,c("SpeciesSciName","maximum.length", "DEMPEL")],by="SpeciesSciName",all.x=T)
  names(dhspp)[which( names(dhspp)=="maximum.length")] <- "MaxL"
  #DEM are ecotypes Demersal + Bathydemersal + Bathypelagic  (except Micromesistius poutassou Blue whiting) + Benthopelagic (except Clupea harengus herring)
  #PEL are Pelagic  #sprat, mac, hmx, her, whb etc
  dhspp$DEMPEL <-as.character(dhspp$DEMPEL)
  dhspp$DEMPEL[dhspp$DEMPEL=="Demersal"] <- "DEM"
  dhspp$DEMPEL[dhspp$DEMPEL=="Pelagic"] <- "PEL"
  dhspp$DEMPEL <-as.character(dhspp$DEMPEL)
#  if(OVERWITE_SUBDIV) {dhspp$SurvStratum<-dhspp$sampstrat} # FOR SOME REASON SAMPSTRAT AND SURVSTRAT ARE THE SAME
#  if(survey %in% c("GNSIntOT1","GNSIntOT3","GNSEngBT4","GNSNetBT3","GNSGerBT3","GNSEngBT3")){dhspp$sampstrat<-dhspp$ICESStSq}
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### Calc ALL indicators ####
  if(EHDS_PP) SAMP_STRAT<-F
  if(SAMP_STRAT){ # dependent on survey see Lynam_OSPARsubdiv.r above# if(survey=="GNSFraOT4"){ dhspp$sampstrat <- dhspp$SurvStratum # minigrid of sqs
    dhspp <- dhspp[!is.na(dhspp$SurvStratum),] #rm outside area
  } else { dhspp$sampstrat <- NA} 
  if(BYSUBDIV){   dhspp$subdiv <- dhspp$SurvStratum; #dhspp$subdiv <-substr(dhspp$subdiv,1,2); # dependent on survey
  dhspp$STRAT_DIV <- paste(dhspp$sampstrat, dhspp$subdiv,sep="_") 
  } else { dhspp$subdiv <- dhspp$STRAT_DIV <- NA} #could be hauls(sampstrat) within subdiv or only 'NA_subdiv'
  STRATA <- T #for output
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  YRS <- sort(unique(dhspp$Year))
  
  FILENAM<-paste(OUTPATH,survey,"_",format(Sys.time(), "%d%b%Y"),sep="")
  FILENAM2<-paste(OUTPATHstem,survey,"_",format(Sys.time(), "%d%b%Y"),sep="")
  #all indicators calc here and biomass by rectangle
  print("Calc Biomass and Indicators")
  print(BYSUBDIV)
  IND_OUT <- INDfn( DATA=dhspp, WRITE=T, SPECIES=SPECIES, LFI_THRESHOLD=LFI_THRESHOLD, 
                    SAMP_STRAT=SAMP_STRAT, BYSUBDIV=BYSUBDIV, FILENAM=FILENAM,FILENAM2=FILENAM2,
                    LFI=LFI, MEANTL=MEANTLs, MaxL=MaxL, MeanL=MeanL, TyL_GeoM=TyL_GeoM)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
} # next survey
print("script 1 complete")
  

library(strucchange);

source(paste(PACKAGEDIRECTORY,"Scripts_and_required_files/Lynam_INDPLOTFN.r", sep = ""))
# plot options
setwd(paste(PACKAGEDIRECTORY,"Outputs/", sep = ""))
#survey<-"BBICsSpaOT1"; SP<-"PEL"

surveys <- c("GNSEngBT3",   "GNSFraOT4",   "GNSGerBT3",   "GNSIntOT1","GNSIntOT3",   "GNSNetBT3",
             "BBICnSpaOT4","BBICPorOT4","BBICsSpaOT1","BBICsSpaOT4", "CSBBFraOT4",
             "CSEngBT3","CSIreOT4","CSNIrOT1","CSNIrOT4","CSScoOT1", "CSScoOT4",
             "WASpaOT3","CSFraOT4","WAScoOT3")
#note Datestamp is important see below

SPECIES<-c("DEM","PEL")

#H<- 6#0.25 # for N in period
#surveys <- c("BBICPorOT4") #("CSScoOT1","CSBBFraOT4","CSFraOT4")
H<- 3 # for BBICPorOT4

#IND<-"MML"#IND<-"TyL";# SP<-"DEM"; survey<-"CSFraOT4"
COLNAMESshifts<- cbind("survey","ACCEPT","CHANGE DIR","% change rel segment1","area","SpeciesGroup", "Indicator", "P",
                       "segment1","segment2","segment3","segment4","segment5","segment6")

write.table(COLNAMESshifts,"SIGshifts2017.csv",sep=',',col.names=F,row.names=F,append = F)
for(IND in c("LFI")){
  for(SP in SPECIES){
    for(survey in surveys){
      #IND<-"MML"; SP<-"DEM"; survey<-"CSFraOT4"
      ifelse(survey %in% c("CSFraOT4","CSBBFraOT4"), DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep=""),DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep=""))
      if(survey %in% c("CSScoOT1","CSBBFraOT4"))  DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep="")
      if(survey %in% c("CSFraOT4"))  DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep="")
      if(IND=="TyL") DATA2PLOT <-  paste(PACKAGEDIRECTORY,"Outputs/",survey,DTE,"_",SP,"_TyL.cm_swept_subdiv_sea.csv",sep="")
      if(IND=="MML") DATA2PLOT <-  paste(PACKAGEDIRECTORY,"Outputs/",survey,DTE,"_",SP,"_MaxL_swept_sea.csv",sep="")
      if(IND=="LFI") DATA2PLOT <-  paste(PACKAGEDIRECTORY,"Outputs/",survey,DTE,SP,"_",SP,"_LFI_subregional.csv",sep="")
      if(!file.exists(DATA2PLOT)){print(paste("no survey",SP," ",survey,sep=' ')); next;}
      DATA2PLOT <-  data.frame(read.csv(DATA2PLOT))
      ifelse(is.null(DATA2PLOT$Year), row.names(DATA2PLOT) <- DATA2PLOT$X,
             row.names(DATA2PLOT) <- DATA2PLOT$Year )
      if(survey!="GNSFraOT4"){ DATA2PLOT<-DATA2PLOT[,-1]}
      if(survey=="GNSFraOT4"){ DATA2PLOT<-DATA2PLOT[,c(ncol(DATA2PLOT),ncol(DATA2PLOT))] }

      #plot TyL LOES BYSUBDIV
      if(nrow(DATA2PLOT)>3){
        if(IND=="LFI") YLAB <- "LFI"
        if(IND=="TyL") YLAB <- "TyL (cm)"
        if(IND=="MML") YLAB <- "MML (cm)"
        YLIM<-c(floor(min(DATA2PLOT,na.rm=T)*.95), ceiling(max(DATA2PLOT,na.rm=T)*1.05) )
        YLIM<-NULL
        N<-ncol(DATA2PLOT)
        winwidth <- ceiling(sqrt(N))
        winhgt <- ifelse(N==2,1,winwidth)
        windows(width=winwidth*8, height=winhgt*8)
        par(mfrow=c(winhgt,winwidth),mar=c(2,4,4,2),oma=c(1,1,3,1))
        for(n in 1:N){
          TITLE <- colnames(DATA2PLOT)[n]
          if(TITLE=="sea.1") next #replicate
          YRSPLOT <- as.numeric(row.names(DATA2PLOT[!is.na(DATA2PLOT[,n]),]))
          INDPLOTFN(DATA2PLOT=DATA2PLOT[!is.na(DATA2PLOT[,n]),n], BOOTDATA2PLOT=NULL,YRS=YRSPLOT,YLIM=YLIM,TITLE=TITLE, GAMMOD=NA, YLAB=YLAB
                    ,BOOTSTRAP=F, ADDBOOTTREND=F, ADDBOOTTREND_CI=F, ADDBOOT_ERRBAR=F, ADDGAM=F, BEST_AND_BOOT=F)

          #shift
          fs<-NULL
          try(fs<- Fstats(DATA2PLOT[!is.na(DATA2PLOT[,n]),n] ~ 1, from = H),silent=T)
          if(!is.null(fs)){
            test<-sctest(fs, type = "sup");
            #title( paste('sup.F =',round(test$statistic,1),'  p =', round(test$p.value,3),sep=' '), line= -1)
            title( paste('sup.F p =', round(test$p.value,3),sep=' '), line= +0.6)

            #constant model
            fm0 <- fm1 <- NULL
            fm0 <- lm(DATA2PLOT[!is.na(DATA2PLOT[,n]),n] ~ 1);
            bp<-breakpoints(DATA2PLOT[!is.na(DATA2PLOT[,n]),n] ~ 1, h = H);

            if(test$p.value<0.05 & (length(bp$breakpoints)>1 || is.na(bp$breakpoints)==F)){
              #alternative model
              ci <-  try(confint(bp), silent =T)
              fm1 <- lm(DATA2PLOT[!is.na(DATA2PLOT[,n]),n] ~ breakfactor(bp) )
              ifelse(test$p.value<0.05, COL<-'dark red', COL<-'dark grey')
              lines(YRSPLOT,fitted(fm1), col = COL,lwd=2,lty=2)
              
              
              FILENAMn<-paste(OUTPATHstem,survey,"_",format(Sys.time(), "%d%b%Y"),sep="")
              JR_BPTout <- data.frame(Year_X=YRSPLOT,Breakpoint_Y=fitted(fm1))
              JR_BPTout$Colour=COL            
              print(paste(FILENAMn,SP,'region',str_replace_all(TITLE, "[^[:alnum:]]", " "),'_breakpoint.csv',sep="_"))
              write.csv(JR_BPTout,paste(FILENAMn,SP,'region',str_replace_all(TITLE, "[^[:alnum:]]", " "),'_breakpoint.csv',sep="_"),row.names=F)
              
              
              
              #change last vs first
              CHANGE <- ifelse(
                (fitted(fm1)[breakfactor(bp)==levels(breakfactor(bp))[length(levels(breakfactor(bp)))]][1]) >
                  fitted(fm1)[breakfactor(bp)=="segment1"][1], "Increase", "Decrease")
              #last vs min level
              ACCEPT <- ifelse(
                (fitted(fm1)[breakfactor(bp)==levels(breakfactor(bp))[length(levels(breakfactor(bp)))]][1]) >
                  min(fitted(fm1)[breakfactor(bp) %in% (levels(breakfactor(bp))[-length(levels(breakfactor(bp)))]) ]) , "OK", "Minimum")

              OUT<- cbind(survey,  ACCEPT, CHANGE,
                          (1-(fitted(fm1)[breakfactor(bp)=="segment1"][1]) /
                             fitted(fm1)[breakfactor(bp)==levels(breakfactor(bp))[length(levels(breakfactor(bp)))]][1]
                          )
                          ,
                          colnames(DATA2PLOT)[n], SP, IND,  test$p.value,
                          fitted(fm1)[breakfactor(bp)=="segment1"][1],
                          fitted(fm1)[breakfactor(bp)=="segment2"][1],
                          fitted(fm1)[breakfactor(bp)=="segment3"][1],
                          fitted(fm1)[breakfactor(bp)=="segment4"][1],
                          fitted(fm1)[breakfactor(bp)=="segment5"][1],
                          fitted(fm1)[breakfactor(bp)=="segment6"][1]
              )
              write.table(OUT,"SIGshifts2017.csv",append = T,sep=',',col.names=F,row.names=F,)

            } else {
              lines(YRSPLOT,fitted(fm0), col = 'dark grey',lwd=2,lty=2)
              OUT<- cbind(survey,"OK","No change","", colnames(DATA2PLOT)[n], SP, IND, test$p.value,"","","","","","")
              #Encoding(OUT[5]) <- "bytes"
              write.table(OUT,"SIGshifts2017.csv",sep=',',append = T,col.names=F,row.names=F,)
            }
          }
          #plus mean lines
          lines(YRSPLOT[1]:(YRSPLOT[1]+5), rep( mean( DATA2PLOT[,n][1:6] ,na.rm=T),6) )
          lines( (YRSPLOT[length(YRSPLOT)]-5):YRSPLOT[length(YRSPLOT)] ,
                 rep(
                   mean( DATA2PLOT[((nrow(DATA2PLOT)-5):nrow(DATA2PLOT)),n]
                         ,na.rm=T ),6) )

        }

        if(SP=="DEM") TITLE <-paste(survey, "Demersal fish", sep=" ")
        if(SP=="PEL") TITLE <-paste(survey, "Pelagic fish", sep=" ")
        if(SP=="ALL") TITLE <-paste(survey, "All fish", sep=" ")
        mtext(TITLE,line=0.5,outer=T)
        savePlot(filename= paste(survey,"_",SP,"_",IND,"_LVL2017.bmp",sep=''),type="bmp")
        dev.off()
      }
    }}
}





# 
# ########### read_all_survey_inds ##################
# 
# library(strucchange);
# source(paste(PACKAGEDIRECTORY,"Scripts_and_required_files/Reflevels.r",sep=""))
# #more complicated below
# #install.packages("rioja")
# #install.packages("ade4")
# source(paste(PACKAGEDIRECTORY,"Scripts_and_required_files/REGIME2.r",sep=""))
# 
# 
# #[1] "BBICnSpaOT4" "BBICPorOT4"  "BBICsSpaOT1" "BBICsSpaOT4" "CSBBFraOT4"  "CSEngBT3"    "CSIreOT4"    "CSNIrOT1"
# #[9] "CSNIrOT4"    "CSScoOT1"    "CSScoOT4"
# #"GNSEngBT3"   "GNSFraOT4"   "GNSGerBT3"   "GNSIntOT1"   "GNSIntOT3"   "GNSNetBT3"
# #"WAScoOT3"    "WASpaOT3"  "CSFraOT4"
# 
# setwd(paste(PACKAGEDIRECTORY,"Outputs/", sep = ""))
# ##### now read GNS indicators all surveys
# CSGNS = c(T,F)
# for (TruFa in CSGNS){
#   CS<-TruFa#GNS
#   surveys <- c("GNSEngBT3",   "GNSFraOT4",   "GNSGerBT3",   "GNSIntOT1","GNSIntOT3",   "GNSNetBT3"   )
#   PELSURY<- c(2,4,5)
#   ##### now read CS indicators all surveysCS<-T
#   if(CS) surveys <- c("BBICPorOT4", "CSBBFraOT4",
#                       "CSEngBT3","CSIreOT4","CSNIrOT1","CSNIrOT4","CSScoOT1", "CSScoOT4",
#                       "CSFraOT4")
#   if(CS) PELSURY<- -c(6)
# 
#   DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep="")
#   sur<-surveys[1]
#   YRS <-  1983:2017
#   YRS_GNS<- 2002:2015 ## all surveys
#   if(CS) YRS_GNS<- 2002:2015 ## all surveys
# 
#   ind<-NULL
#   ###############################################################################
#   #Surv_biotonnes_Year
# 
#   SP <- "DEM"
#   BIO<- matrix(NA,ncol=length(surveys),nrow=length(YRS))
#   rownames(BIO) <-YRS
#   colnames(BIO) <-surveys
#   for(sur in surveys){
#     print(sur)
#     ifelse(sur %in% c("CSFraOT4","CSBBFraOT4"), DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep=""),DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep=""))
#     if(sur %in% c("CSBBFraOT4","CSScoOT1"))  DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep="")
#     if(sur %in% c("CSFraOT4"))  DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep="")
#     ind[[sur]] <-
#       read.csv( paste(sur,"/",sur,DTE,SP,"_Surv_biotonnes_Year.csv",sep='') )
#     BIO[rownames(BIO) %in% ind[[sur]]$Year,which(colnames(BIO)==sur)] <- ind[[sur]]$CatCatchWgtSwept
#   }
#   BIO<-data.frame(BIO)
#   cor(BIO,use="pairwise.complete.obs")
# 
#   BIOl <- BIO
#   names(BIOl) <- paste("Biomass",names(BIOl))
#   #REFLVL(BIOl[,-3],H=0.3)
# 
#   #plot
#   windows()
#   par(mfrow=c(2,2))
#   plot(YRS, BIO[,1]/1000,col="white",ylim=range(BIO,na.rm=T)/1000,main="Demersal fish", ylab="Biomass (kt)",xlab="Year")
#   for(sur in surveys) lines(YRS, BIO[,sur]/1000,col=which(sur == surveys))
#   plot(YRS, BIO[,1]/1000,col="white",ylim=range(BIO,na.rm=T)/1000,main="", ylab="",axes=F,xlab="")
#   legend("left",surveys,lty=1,col= 1:length(surveys))
# 
# 
# 
# 
#   if(CS) surveys <- c("BBICPorOT4", "CSBBFraOT4",
#                       "CSIreOT4","CSNIrOT1","CSNIrOT4","CSScoOT1", "CSScoOT4",
#                       "CSFraOT4")
#   if(CS) PELSURY<- -c(6)
#   BIO<- matrix(NA,ncol=length(surveys),nrow=length(YRS))
#   rownames(BIO) <-YRS
#   colnames(BIO) <-surveys
# 
#   SP <- "PEL"
#   BIOp<- matrix(NA,ncol=length(surveys),nrow=length(YRS))
#   rownames(BIOp) <-YRS
#   colnames(BIOp) <-surveys
# 
#   for(sur in surveys[PELSURY]){
#     print(sur)
#     ifelse(sur %in% c("CSFraOT4","CSBBFraOT4"), DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep=""),DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep=""))
#     if(sur %in% c("CSBBFraOT4","CSScoOT1"))  DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep="")
#     if(sur %in% c("CSFraOT4"))  DTE <- paste("_",format(Sys.time(), "%d%b%Y"),sep="")
# 
#     ind[[sur]] <- read.csv( paste(sur,"/",sur,DTE,SP,"_Surv_biotonnes_Year.csv",sep='') )
#     BIOp[rownames(BIOp) %in% ind[[sur]]$Year,which(colnames(BIOp)==sur)] <- ind[[sur]]$CatCatchWgtSwept
#   }
#   cor(BIOp,use="pairwise.complete.obs")
#   #plot pelagic par(mfrow=c(1,2))
#   plot(YRS, BIOp[,1]/1000,col="white",ylim=range(BIOp,na.rm=T)/1000,main="Pelagic fish", ylab="Biomass (kt)",xlab="Year")
#   for(sur in surveys) lines(YRS, BIOp[,sur]/1000,col=which(sur == surveys))
# #  if(!CS) savePlot(filename = "GNSbiomass.bmp",type="bmp")
# #  if(CS) savePlot(filename = "CSbiomass.bmp",type="bmp")
# 
#   BIOpl<-data.frame(BIOp)
#   names(BIOpl) <- paste("BiomassPelagic",names(BIO))
#   #REFLVL(BIOpl[,PELSURY],H=0.3)
# 
#   BAR<-apply(X=BIOp,FUN=mean,2,na.rm=T)
#   SD<-apply(X=BIOp,FUN=sd,2,na.rm=T)
#   BIOps <- scale(BIOp)
#   #plot(YRS, BIOps[,1],col="white",ylim=range(BIOps,na.rm=T),main="GNS Pelagic fish", ylab="Biomass (scaled)",xlab="Year")
#   #for(sur in surveys) lines(YRS, BIOps[,sur],col=which(sur == surveys))
# }
# ###########################################################################
# 
# Delete folders created that contain outputs we don't need for TyL indicator
for(survey in  SURVEY){ # survey <- "GNSGerBT3" #SURVEY[1]#"WAScoOT3" # #"CSFraOT4"#WASpaOT3#[4]# BBICsSpaOT1 and BBICsSpaOT4 missing datayear causes warnings in LFI fn so have suppressed these
  OUTPATH <- paste(OUTPATHstem,survey,sep='')
  unlink(OUTPATH, recursive = T)
}



#Load the results Csv for matching with shapefiles. rename regions/columns that don't match that of the shapefiles before merging
sigshiftscsv = read.csv(paste(PACKAGEDIRECTORY,"Outputs/SIGshifts2017.csv", sep = ""),stringsAsFactors = FALSE)
names(sigshiftscsv)[which(names(sigshiftscsv) %in% "area")]<- "SurvStratum"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='FRcoast..25m']="FRcoast <25m"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='UKcoast..25m']="UKcoast <25m"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='Eastern.Irish.Sea...50m']="Eastern Irish Sea, <50m"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='Irish.Coast....50m']="Irish Coast, < 50m"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='Isle.of.Man..50...100m']="Isle of Man, 50 - 100m"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='St.George.s.Channel..100m']="St George's Channel <100m"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X0']="0"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X1']="1"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X2']="2"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X3']="3"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X4']="4"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X5']="5"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X6']="6"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X7']="7"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X8']="8"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X9']="9"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X10']="10"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X11']="11"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X12']="12"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X13']="13"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X14']="14"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X15']="15"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X16']="16"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X17']="17"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X18']="18"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X19']="19"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X20']="20"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X21']="21"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X22']="22"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X23']="23"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X24']="24"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X25']="25"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X26']="26"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X27']="27"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X28']="28"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X29']="29"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X30']="30"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X32']="32"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X33']="33"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X38']="38"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X37']="37"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X39']="39"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X40']="40"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X42']="42"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X1931']="1931"
sigshiftscsv["SurvStratum"][sigshiftscsv["SurvStratum"]=='X2236']="2236"

sigshiftscsv["CHANGE.DIR"][sigshiftscsv["ACCEPT"]=='Minimum']="Decrease to minimum"

#names(sigshiftscsv)[which(names(sigshiftscsv) %in% NAMsubdiv)] <- "SurvStratum"; SAMP_FACT <- c(SAMP_FACT, "SurvStratum")



setwd(paste(PACKAGEDIRECTORY,"Outputs/", sep = ""))
for(spcs in c("DEM")){
  for(survey in SURVEY){
    #,"CSFraOT4","GNSFraOT4"
    
    if(survey %in% c("GNSIntOT1","GNSIntOT3")){
      SUBDIVshapefile <- readShapeSpatial(paste(SHAPEPATH,"GNS_rectstrat/GNSIntOT/GNSstrat_Atlantis.shp",sep='') ) 
      names(SUBDIVshapefile)[which(names(SUBDIVshapefile) %in% "LFIregion")]<- "SurvStratum"
      names(SUBDIVshapefile)[which(names(SUBDIVshapefile) %in% "ICESNAME")]<- "sampstrat"}

    if(survey %in% c("GNSGerBT3")){
      SUBDIVshapefile <- readShapeSpatial(paste(SHAPEPATH,"GNS_rectstrat/GNSGerBT3/GNSstrat_Atlantis.shp",sep='') ) 
      names(SUBDIVshapefile)[which(names(SUBDIVshapefile) %in% "NAME")]<- "SurvStratum"
      names(SUBDIVshapefile)[which(names(SUBDIVshapefile) %in% "ICESNAME")]<- "sampstrat"}

    if(survey %in% c("GNSNetBT3")){
      SUBDIVshapefile <- readShapeSpatial(paste(SHAPEPATH,"GNS_rectstrat/GNSNetBT3/GNSstrat_Atlantis.shp",sep='') ) 
      names(SUBDIVshapefile)[which(names(SUBDIVshapefile) %in% "NAME")]<- "SurvStratum"
      names(SUBDIVshapefile)[which(names(SUBDIVshapefile) %in% "ICESNAME")]<- "sampstrat"}

    if(survey %in% c("GNSEngBT4")){
      SUBDIVshapefile <- readShapeSpatial(paste(SHAPEPATH,"WChanEngBeam/WChanBT4.shp",sep='') ) 
      names(SUBDIVshapefile)[which(names(SUBDIVshapefile) %in% "Stratum")]<- "SurvStratum"}
    
    if(survey %in% c("GNSEngBT3")){
      SUBDIVshapefile <- readShapeSpatial(paste(SHAPEPATH,"EChanEngBeamSimple/EChanBT3Simple.shp",sep='') ) 
      names(SUBDIVshapefile)[which(names(SUBDIVshapefile) %in% "name")]<- "SurvStratum"}
    
    if(survey %in% c("CSEngBT3")){
      SUBDIVshapefile <- readShapeSpatial(paste(SHAPEPATH,"irish_seaBT//NI_IBTS.WGS84.shp",sep='') ) 
      names(SUBDIVshapefile)[which(names(SUBDIVshapefile) %in% "Features")]<- "SurvStratum"}
    
    if(survey %in% c("CSEngBT3bc")){
      SUBDIVshapefile <- readShapeSpatial(paste(SHAPEPATH,"BChanEngBeam/BChanBT3.shp",sep='') ) 
      names(SUBDIVshapefile)[which(names(SUBDIVshapefile) %in% "Stratum")]<- "SurvStratum"}
    
    if(survey %in% c("CSNIrOT4","CSNIrOT1")){
      SUBDIVshapefile <- readShapeSpatial(paste(SHAPEPATH,"irish_seaGFS//NI_IBTS.WGS84.shp",sep='') ) 
      names(SUBDIVshapefile)[which(names(SUBDIVshapefile) %in% "Features")]<- "SurvStratum"}
    
    if(survey %in% c("CSScoOT1")){
      SUBDIVshapefile <- readShapeSpatial(paste(SHAPEPATH,"SWCQ1.WGS84//SWC_Q1.shp",sep='') ) 
      names(SUBDIVshapefile)[which(names(SUBDIVshapefile) %in% "Name")]<- "SurvStratum"}
    
    if(survey %in% c("CSScoOT4")){
      SUBDIVshapefile <- readShapeSpatial(paste(SHAPEPATH,"SWCQ4.WGS84//SWC_Q4.shp",sep='') ) 
      names(SUBDIVshapefile)[which(names(SUBDIVshapefile) %in% "Name")]<- "SurvStratum"}
    
    if(survey %in% c("WAScoOT3")){
      SUBDIVshapefile <- readShapeSpatial(paste(SHAPEPATH,"SWC-RockQ3.WGS84//SWC_Q3.shp",sep='') ) 
      names(SUBDIVshapefile)[which(names(SUBDIVshapefile) %in% "Name")]<- "SurvStratum"}
    
    if(survey %in% c("BBICnSpaOT1","BBICnSpaOT4")){
      SUBDIVshapefile <- readShapeSpatial(paste(SHAPEPATH,"Sp-NGFS.WGS84//Sp_North.WGS84.shp",sep='') ) 
      names(SUBDIVshapefile)[which(names(SUBDIVshapefile) %in% "Primary")]<- "SurvStratum"}
    
    if(survey %in% c("BBICsSpaOT1")){
      SUBDIVshapefile <- readShapeSpatial(paste(SHAPEPATH,"Sp-Cadiz.WGS84/Sp_Cadiz.WGS84.shp",sep='') ) 
      names(SUBDIVshapefile)[which(names(SUBDIVshapefile) %in% "ESTRATO")]<- "SurvStratum"}
    
    if(survey %in% c("BBICsSpaOT4")){
      SUBDIVshapefile <- readShapeSpatial(paste(SHAPEPATH,"Sp-Cadiz.WGS84/Sp_Cadiz.WGS84.shp",sep='') ) 
      names(SUBDIVshapefile)[which(names(SUBDIVshapefile) %in% "ESTRATO")]<- "SurvStratum"}
    
    if(survey %in% c("WASpaOT3")){
      SUBDIVshapefile <- readShapeSpatial(paste(SHAPEPATH,"Sp-PorcGFS.WGS84/Porcupine.WGS84.shp",sep='') ) 
      names(SUBDIVshapefile)[which(names(SUBDIVshapefile) %in% "Primary")]<- "SurvStratum"}
    
    if(survey %in% c("CSIreOT4")){
      SUBDIVshapefile <- readShapeSpatial(paste(SHAPEPATH,"IGFS.WGS84//IGFS.WGS84.shp",sep='') ) 
      names(SUBDIVshapefile)[which(names(SUBDIVshapefile) %in% "Primary")]<- "SurvStratum"}
    
    if(survey %in% c("CSBBFraOT4")){
      SUBDIVshapefile <- readShapeSpatial(paste(SHAPEPATH,"Fr-EVHOE.WGS84//EVHOE.WGS84.shp",sep='') ) 
      names(SUBDIVshapefile)[which(names(SUBDIVshapefile) %in% "STRATE")]<- "SurvStratum"}
    
    if(survey %in% c("CSFraOT4")){
      SUBDIVshapefile <- readShapeSpatial(paste(SHAPEPATH,"Fr-EVHOE.WGS84_original//EVHOE.WGS84.shp",sep='') ) 
      names(SUBDIVshapefile)[which(names(SUBDIVshapefile) %in% "STRATE")]<- "SurvStratum"}
#      names(SUBDIVshapefile)[which(names(SUBDIVshapefile) %in% "STRATE_BAT")]<- "sampstrat"}
  
    
    if(survey %in% c("BBICPorOT4")){
      SUBDIVshapefile <- readShapeSpatial(paste(SHAPEPATH,"BBICPorOT4/Contour3strata_sector.shp",sep='') ) 
      names(SUBDIVshapefile)[which(names(SUBDIVshapefile) %in% "ID_1")]<- "SurvStratum"}
    
    if(survey %in% c("GNSFraOT4")){
      SUBDIVshapefile <- readShapeSpatial(paste(SHAPEPATH,"GNSFraOT4/GNSFraOT4_EngBT3Simple.shp",sep='') ) 
      names(SUBDIVshapefile)[which(names(SUBDIVshapefile) %in% "name")]<- "SurvStratum"
      names(SUBDIVshapefile)[which(names(SUBDIVshapefile) %in% "FID_GNSFra")]<- "sampstrat"}

    
    sigshiftssurvey = sigshiftscsv[sigshiftscsv$survey==survey & sigshiftscsv$SpeciesGroup==spcs,]
    mergeddf = merge(SUBDIVshapefile,sigshiftssurvey,by="SurvStratum")
    
    # Add data for shapefiles with sample strata
    if(survey %in% c("GNSFraOT4","GNSNetBT3","GNSIntOT1","GNSIntOT3","GNSGerBT3")){if(spcs %in% c("DEM")){
      
      FILENAM2<-paste(OUTPATHstem,survey,"_",format(Sys.time(), "%d%b%Y"),sep="")
      sampstratcsv = t(read.csv(paste(FILENAM2,spcs,'_',spcs,'_LFI_sampstrat.csv',sep=""),check.names=FALSE,stringsAsFactors = FALSE))
      colnames(sampstratcsv) <- sampstratcsv[1, ]; sampstratcsv= sampstratcsv[-1, ]
      sampstratcsv <- as.data.frame(sampstratcsv);     sampstratcsv$sampstrat=rownames(sampstratcsv)
            #      names(sampstratcsv)[which(names(sampstratcsv) =="")] <-"year"
      mergeddf2=merge(mergeddf,sampstratcsv,by="sampstrat")
      
    }} else {

      FILENAM2<-paste(OUTPATHstem,survey,"_",format(Sys.time(), "%d%b%Y"),sep="")
      subregcsv = t(read.csv(paste(FILENAM2,spcs,'_',spcs,'_LFI_subregional.csv',sep=""),check.names=FALSE,stringsAsFactors = FALSE))
      colnames(subregcsv) <- subregcsv[1, ]; subregcsv= subregcsv[-1, ]
      subregcsv <- as.data.frame(subregcsv);     subregcsv$SurvStratum=rownames(subregcsv)
      mergeddf2=merge(mergeddf,subregcsv,by="SurvStratum")
    }
    nameoutstr = paste(survey,"_",spcs,"_results", sep = "")
    writeOGR(mergeddf2, ".", nameoutstr, driver="ESRI Shapefile") # paste(PACKAGEDIRECTORY,"Outputs/",survey,"_",spcs,"_results.shp", sep = "")
  }
}
# What about these GNS surveys? Are they not sampled by rect?
#"GNSNetBT3","GNSIntOT3"


# If error remove the -1000000 from Sys.time()-100000
#for(spcs in c("DEM")){
#  for(survey in c("CSFraOT4","GNSFraOT4","GNSIntOT1","GNSGerBT3")){

# year is X, and will have to be a result in the shapefile after merging


