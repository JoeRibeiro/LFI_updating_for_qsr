#20 dec 2016 # edit to avoid creating NAs in  c("BBICsSpaOT1", "BBICsSpaOT4","WASpaOT3") when missing data by year/strata
INDfn_LFI <- function(species_bio_by_area, numhaulsyr, numsampstrat_by_sea, SP, WRITE=F, FILENAM="",FILENAM2="",BYSUBDIV=T,BYSAMPSTRAT=F){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ LARGE fish Indicator prep  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~ LARGE fish -> species_bioL_by_area and by region and overall sea area  ~~~~~~~
    LFI_FACT<- c("Year","FishLength_cm","SpeciesSciName")
    if(BYSUBDIV) LFI_FACT<- c(LFI_FACT,"subdiv")
    if(BYSAMPSTRAT) LFI_FACT<- c(LFI_FACT,"sampstrat")
    # by length class  #species sum bio cpue at len by rect
    if(nrow(species_bio_by_area[species_bio_by_area$FishLength_cm >LFI_THRESHOLD,])>0){
      suppressWarnings( #NAs introduced due to missing years on following surveys: BBICsSpaOT1
         species_bioL_by_area <- tapply.ID(df=species_bio_by_area[species_bio_by_area$FishLength_cm >LFI_THRESHOLD,], 
                                      datacols=c("CatCatchWgtSwept"), 
                                      factorcols=LFI_FACT, 
                                      sum,c("CatCatchWgtSwept_Large"))
                        )
       
      if(survey %in% c("BBICsSpaOT1", "BBICsSpaOT4","WASpaOT3")) species_bioL_by_area <- species_bioL_by_area[!is.na(species_bioL_by_area$Year),]
    } else { species_bioL_by_area <- species_bio_by_area[species_bio_by_area$FishLength_cm >LFI_THRESHOLD,] }
    #large average over hauls by sampstrat (e.g. rectangle) #species mean bio cpue at len by sampstrat
    if(nrow(species_bioL_by_area)==0){ print(paste("no data above LFI threshold for ",SP,sep='')); LFIout<-LFI_by_sub<-NULL; } else {
    
    #~~~~~~~~~~~~~~ all and large fish by subdivision (already ave by #hauls in sampstrat (rect or other)
      if(BYSUBDIV){
        suppressWarnings( #NAs introduced due to missing years on following surveys: BBICsSpaOT1
          species_bio_by_subdiv <- tapply.ID(df=species_bio_by_area, datacols=c("CatCatchWgtSwept"), 
                                             factorcols=c("Year","FishLength_cm","SpeciesSciName","subdiv"), sum,c("CatCatchWgtSwept"))
        )
        if(survey %in% c("BBICsSpaOT1", "BBICsSpaOT4","WASpaOT3")) species_bio_by_subdiv <- species_bio_by_subdiv[!is.na(species_bio_by_subdiv$Year),]
        
        species_bioL_by_subdiv <- tapply.ID(df=species_bioL_by_area, datacols=c("CatCatchWgtSwept_Large"), 
                                            factorcols=c("Year","FishLength_cm","SpeciesSciName","subdiv"), sum,c("CatCatchWgtSwept_Large"))
      }
      if(BYSAMPSTRAT){
        suppressWarnings( #NAs introduced due to missing years on following surveys: BBICsSpaOT1
          species_bio_by_sampstrat <- tapply.ID(df=species_bio_by_area, datacols=c("CatCatchWgtSwept"), 
                                             factorcols=c("Year","FishLength_cm","SpeciesSciName","sampstrat"), sum,c("CatCatchWgtSwept"))
        )
        if(survey %in% c("BBICsSpaOT1", "BBICsSpaOT4","WASpaOT3")) species_bio_by_subdiv <- species_bio_by_subdiv[!is.na(species_bio_by_subdiv$Year),]
        
        species_bioL_by_sampstrat <- tapply.ID(df=species_bioL_by_area, datacols=c("CatCatchWgtSwept_Large"), 
                                            factorcols=c("Year","FishLength_cm","SpeciesSciName","sampstrat"), sum,c("CatCatchWgtSwept_Large"))
      }
      
    #~~~~~~~~~~~~~~ all and large fish by regional sea scale from sampstrat not subdiv
    # from rects if SAMP_STRAT/subdiv above or from subdiv only
    # corrected for any change in sampling between subdiv
    suppressWarnings( #NAs introduced due to missing years on following surveys: BBICsSpaOT1
      species_bio_by_sea <- tapply.ID(df=species_bio_by_area, datacols=c("CatCatchWgtSwept"), 
                                      factorcols=c("Year","FishLength_cm","SpeciesSciName"), sum,c("CatCatchWgtSwept"))
    )
    if(survey %in% c("BBICsSpaOT1", "BBICsSpaOT4","WASpaOT3")) species_bio_by_sea <- species_bio_by_sea[!is.na(species_bio_by_sea$Year),]
    
    
    species_bioL_by_sea <- tapply.ID(df=species_bioL_by_area, datacols=c("CatCatchWgtSwept_Large"), 
                                    factorcols=c("Year","FishLength_cm","SpeciesSciName"), sum,c("CatCatchWgtSwept_Large"))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  #  sum for LFI
  # sum numerator of LFI
  #browser()
  LFInum <- tapply.ID(df=species_bioL_by_sea, datacols=c("CatCatchWgtSwept_Large"), factorcols=c("Year"), sum,c("CatCatchWgtSwept_Large"))
  # denominator of LFI
  LFIden <- tapply.ID(df=species_bio_by_sea, datacols=c("CatCatchWgtSwept"), factorcols=c("Year"), sum,c("CatCatchWgtSwept"))
  #sometimes no large fish in num so be careful with years...
  rownames(LFIden)<- LFIden[,2];   
  LFIy <- LFIden; LFIy[,1]<-LFIy[,1]*0; #set up holder
  
  LFIy[rownames(LFIy)%in%LFInum[,2],1] <-  LFInum[,1] / LFIden[rownames(LFIy)%in%LFInum[,2],1]
  #plot(LFI[,2:1],type='b',xlab="",ylim=c(0,1.3e9)) #perfect!
  
  #lfi and out; LFI is num/den, so, numerator is lfi * denom
  LFIout <- data.frame(Year=LFIden[,2],numhauls=numhaulsyr[,1],numsampstrata=numsampstrat_by_sea[,1], 
                       numerator=LFIy[,1]*LFIden[,1], denominator=LFIden[,1], LFIregional=LFIy[,1])
  if(WRITE) write.csv(LFIout,paste(FILENAM2,SP,'LFI_sea.csv',sep="_"),row.names=F)
  
  # LFI
  if(BYSUBDIV){
    
    FACT <- c("Year","subdiv")
    # sum numerator of LFI by sub divisions
    LFInumreg <- tapply.ID(df=species_bioL_by_subdiv, datacols=c("CatCatchWgtSwept_Large"), factorcols=FACT, sum,c("CatCatchWgtSwept_Large"))
    # denominator of LFI by sub divisions
    LFIdenreg <- tapply.ID(df=species_bio_by_subdiv, datacols=c("CatCatchWgtSwept"), factorcols=FACT, sum,c("CatCatchWgtSwept"))
    LFIreg <- merge(LFInumreg,LFIdenreg,by=FACT,all.y=T)
    LFIreg$LFI <- LFIreg[,'CatCatchWgtSwept_Large']/ LFIreg[,'CatCatchWgtSwept']
    
    LFI_by_sub <- xtabs(LFI ~ Year + subdiv, LFIreg) # end up with missing years if no large fish
    if(WRITE) write.csv(LFI_by_sub,paste(FILENAM2,SP,'LFI_subregional.csv',sep="_"),row.names=T)
    LFI_by_suborig=LFI_by_sub
      } else { LFI_by_sub <- NULL }
  if(BYSAMPSTRAT){
    
    FACT <- c("Year","sampstrat")
    # sum numerator of LFI by sub divisions
    LFInumreg <- tapply.ID(df=species_bioL_by_sampstrat, datacols=c("CatCatchWgtSwept_Large"), factorcols=FACT, sum,c("CatCatchWgtSwept_Large"))
    # denominator of LFI by sub divisions
    LFIdenreg <- tapply.ID(df=species_bio_by_sampstrat, datacols=c("CatCatchWgtSwept"), factorcols=FACT, sum,c("CatCatchWgtSwept"))
    LFIreg <- merge(LFInumreg,LFIdenreg,by=FACT,all.y=T)
    LFIreg$LFI <- LFIreg[,'CatCatchWgtSwept_Large']/ LFIreg[,'CatCatchWgtSwept']
    
    LFI_by_samp <- xtabs(LFI ~ Year + sampstrat, LFIreg) # end up with missing years if no large fish
    if(WRITE) write.csv(LFI_by_samp,paste(FILENAM2,SP,'LFI_sampstrat.csv',sep="_"),row.names=T)
    } else { LFI_by_samp <- NULL }
  
  
  
  #plot code
  if(!BOOTSTRAP & BYSUBDIV & !is.null(LFI_by_sub)){
    if(nrow(LFI_by_sub)>3){
    YLAB <- "LFI"
    
    DATA2PLOT <- LFI_by_sub 
    YLIM<-c(floor(min(DATA2PLOT,na.rm=T)*.95), ceiling(max(DATA2PLOT,na.rm=T)*1.05) )
    N<-ncol(DATA2PLOT)
    winwidth <- ceiling(sqrt(N))
    winhgt <- ifelse(N==2,1,winwidth)
    windows(width=winwidth*8, height=winhgt*8)
    par(mfrow=c(winhgt,winwidth),mar=c(2,4,2,2),oma=c(1,1,3,1))
    for(n in 1:N){
      TITLE <- colnames(DATA2PLOT)[n]
      YRSPLOT <- an(names(DATA2PLOT[!is.na(DATA2PLOT[,n]),n]))
      INDPLOTFN(DATA2PLOT=DATA2PLOT[!is.na(DATA2PLOT[,n]),n], BOOTDATA2PLOT=NULL,YRS=YRSPLOT,YLIM=YLIM,TITLE=TITLE, GAMMOD=NA, YLAB=YLAB
                ,BOOTSTRAP=F, ADDBOOTTREND=F, ADDBOOTTREND_CI=F, ADDBOOT_ERRBAR=F, ADDGAM=F, BEST_AND_BOOT=F,FILENAM3=FILENAM2,SP=SP)
    }
    
    if(SP=="DEM") TITLE <-paste(survey, "Demersal fish", sep=" ")
    if(SP=="PEL") TITLE <-paste(survey, "Pelagic fish", sep=" ")
    if(SP=="ALL") TITLE <-paste(survey, "All fish", sep=" ")
    mtext(TITLE,line=0.5,outer=T)
    savePlot(filename= paste(FILENAM2,survey,"_",SP,"_",format(Sys.time(), "%d%b%Y"),"_LFI.bmp",sep=''),type="bmp")
    
    dev.off()  
    }
  }
    }
  return(list(LFIout, LFI_by_sub))
  
}