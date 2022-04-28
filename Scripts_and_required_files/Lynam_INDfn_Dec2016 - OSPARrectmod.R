#ChrisLynam@Cefas.co.uk
# update 20 Dec 2016 to give LD by species and total biomass by strata
# update 09 Jan 2017 to pass SP to indicators and enable plotting LOESS fns
# species_bio_by_area includes numhauls by subdivision if no sampling strata applied (i.e. no rectangles or minigrid)
# average over hauls then raise up by area of sampstrat or subdivision
##second level raising if needed, here weighting is correct for coverage of subdiv (sum of sampstrat areas should be subdiv area - but prone to issues with missing sampstrat)

INDfn <- function(DATA, WRITE=F, BOOT=F, LFI=T, LFI_THRESHOLD=40, FILENAM="",FILENAM2="", SAMP_STRAT=T, BYSUBDIV=T,
                   MEANTL=T, MaxL=T, MeanL=T, TyL_GeoM=F, SPECIES = c("DEM") ){
   
    #SAMP_STRAT<-T; BYSUBDIV<-T; SPECIES<-c("ALL","DEM","PEL"); BOOT<-F
    #DATA<-dhspp; BOOT<-F;  WRITE=T; LFI=T; MEANTL=MEANTLs; MaxL=T; MeanL=T; TyL_GeoM=T; LFI_THRESHOLD<-40; #possible outputs #
    #DATA<-dhspp; MEANTL<-MEANTLs
    LFIout    <- LFI_by_sub    <- FishLength_cmsea    <- MaxLsea    <- TLsea    <- TyL.cm.sea    <- NULL
    LFIout$LFIregional <- NULL
    LFIoutpel <- LFI_by_subpel <- FishLength_cmseapel <- MaxLseapel <- TLseapel <- TyL.cm.seapel <- LFIout$LFIregionalpel <- NULL
    LFIoutdem <- LFI_by_subdem <- FishLength_cmseadem <- MaxLseadem <- TLseadem <- TyL.cm.seadem <- LFIout$LFIregionaldem <- NULL
    
    TyL.cm.sea_pel<- FishLength_cmsea_pel <- MaxLsea_pel <- TLsea_pel <-  LFIbind_pel <- NULL  
    TyL.cm.sea_dem<- FishLength_cmsea_dem <- MaxLsea_dem <- TLsea_dem <-  LFIbind_dem <- NULL  
    TyL.cm.sea_all<- FishLength_cmsea_all <- MaxLsea_all <- TLsea_all <-  LFIbind_all <- NULL  
    

      # how many hauls?    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        numhauls <- DATA; # numhauls <- dhspp
        numhauls$ones <- 0 # 1 val as a marker to pivot around
        #HaulID and StNo are lost so now using HaulID
        FACTHAUL <-  c("Year","HaulID","ShootLat_degdec","ShootLong_degdec")
        if(SAMP_STRAT) FACTHAUL <-  c(FACTHAUL,"sampstrat")
        if(BYSUBDIV)   FACTHAUL <-  c(FACTHAUL,"subdiv","STRAT_DIV")
        
        numhauls <- tapply.ID(df=numhauls, datacols=c("ones"),factorcols=FACTHAUL,sum,c("ones"))
        numhauls$ones <- 1  # now 1 val per haul      
        if(WRITE) write.table(numhauls,paste(FILENAM,"numhauls.csv",sep="_"),row.names =F,sep=',')
        
        #add centlat centlon i.e. centre points of ICES stsqs
        numhauls$centlon <- floor(numhauls$ShootLong_degdec)+.5
        numhauls$centlat <- round(numhauls$ShootLat_degdec)
        cor <- ifelse(numhauls$centlat < numhauls$ShootLat_degdec, 0.25,  
                      ifelse(numhauls$centlat > numhauls$ShootLat_degdec, -0.25,  
                             + 0.25) #if x.00
                      )
        numhauls$centlat <- (numhauls$centlat + cor) 
        
        if(SAMP_STRAT){ #might be STSQ or minigrid see sampstrat
          FACTHAUL <-  c("Year","centlon","centlat","sampstrat")
          if(BYSUBDIV) FACTHAUL <-  c(FACTHAUL,"subdiv","STRAT_DIV")
          
          numhaulsBYsampstrat <- tapply.ID(df=numhauls, datacols=c("ones"), factorcols=FACTHAUL, sum,c("numhauls"));
          #if(WRITE) write.table(numhaulsBYsampstrat,paste(FILENAM,"numhaulsBYsampstrat.csv",sep="_"),row.names =F,sep=',')
          #and reshape since have one value per year and subdiv combination
          numhaulsBYsampstratout <- (tapply(numhaulsBYsampstrat$numhauls,list(numhaulsBYsampstrat$Year, numhaulsBYsampstrat$sampstrat), FUN=sum, na.rm=T))
          if(WRITE) write.table(numhaulsBYsampstratout,paste(FILENAM,"numhaulsBYsampstrat.csv",sep="_"),row.names =T,sep=',')
          rm(numhaulsBYsampstratout)
        } else { numhaulsBYsampstrat <- NULL }
        
        if(BYSUBDIV){#user_defined or survey poly
          numhaulsBYsubdiv <- tapply.ID(df=numhauls, datacols=c("ones"),
                                        factorcols=c("Year","subdiv"),
                                        sum,c("numhauls"));
          #and reshape since have one value per year and subdiv combination
          numhaulsBYsubdivout <- (tapply(numhaulsBYsubdiv$numhauls,list(numhaulsBYsubdiv$Year, numhaulsBYsubdiv$subdiv), FUN=sum, na.rm=T))
#          if(WRITE) write.table(numhaulsBYsubdivout,paste(FILENAM2,"numhaulsBYsubdiv.csv",sep="_"),row.names =T,sep=',')
# J Ribeiro mdified from C.Lynam's verson as the tables had an offset when using write.table
          if(WRITE) write.csv(numhaulsBYsubdivout,paste(FILENAM2,"numhaulsBYsubdiv.csv"),row.names =T,sep=',')
          rm(numhaulsBYsubdivout)
        } else { numhaulsBYsubdiv <- NULL }
        numhaulsyr <- tapply.ID(df=numhauls, datacols=c("ones"),factorcols=c("Year"),sum,c("numhauls"));  # now 1 val per STSQ
        numhauls<-numhauls[,-1]#now just a list of hauls
        #browser()
        #plot(numhaulsBYsubdiv[numhaulsBYsubdiv$BOX_ID==11,2:1])
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   species_bio_by_area    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # species bioLD by rect (strata1) and subdivision (strata2) from haul data
      
      # run tapply.id over these factorcols FACT
       FACT <- c("Year","FishLength_cm","SpeciesSciName")
       if(SAMP_STRAT)FACT <- c(FACT,"sampstrat")
       if(BYSUBDIV)  FACT <- c(FACT,"subdiv")
       #if(MaxL)      FACT <- c(FACT,"MaxL")
       #if(MEANTL)    FACT <- c(FACT,"TL")
      
      # sum by species by length cat by sampstrat (e.g ICES STSQ)
       # and subdivisional strata (e.g. NE North Sea or 'survstrata') #DATA<-dhspp
       suppressWarnings( #NAs introduced by coercion since some MaxL and TL are NA
         species_bio_by_area <- tapply.ID(df=DATA, datacols=c("DensBiom_kg_Sqkm"), factorcols=FACT, sum,c("CatCatchWgtSwept"))  
         )  
       # to average LD must add num hauls to bio data
       if(SAMP_STRAT){
         species_bio_by_area <- merge(x = species_bio_by_area,
                                      y = numhaulsBYsampstrat[,which(names(numhaulsBYsampstrat) != "STRAT_DIV" & names(numhaulsBYsampstrat) != "subdiv")], #avoid replicating names and creating .x .y
                                      by = c("Year","sampstrat"),all.x=T)
         
         #average species cpue over hauls by rectangle-strata for MaxL, TL , Len, TyL                                                                      
         species_bio_by_area$CatCatchWgtSwept <- species_bio_by_area$CatCatchWgtSwept / species_bio_by_area$numhauls
       } else {
         if(BYSUBDIV){ # only do here if not using rects/minigrid/etc
           species_bio_by_area <- merge(x = species_bio_by_area,
                                        y = numhaulsBYsubdiv[,which(names(numhaulsBYsubdiv) != "STRAT_DIV")], #avoid replicating names and creating .x .y
                                        by = c("Year","subdiv"),all.x=T)
           
           #average species cpue over hauls by rectangle-strata for MaxL, TL , Len, TyL                                                                      
           species_bio_by_area$CatCatchWgtSwept <- species_bio_by_area$CatCatchWgtSwept / species_bio_by_area$numhauls
         }
       }
       #if both SAMP_STRAT and BYSUBDIV are true will have to sum up catch by SAMP_STRAT within SUBDIV later to avoid change between years due to change in relative sampling of strata..
       
       #introduce MAXL, TL, DEMPEL as lost now! warning here is an opportunity for NAs to appear!
       species_bio_by_area <- merge(species_bio_by_area,trait_MAXL[,c("SpeciesSciName","maximum.length", "DEMPEL")],by="SpeciesSciName",all.x=T)
       # species_bio_by_area[is.na(species_bio_by_area$DEMPEL),] 
       names(species_bio_by_area)[which( names(species_bio_by_area)=="maximum.length")] <- "MaxL"
       species_bio_by_area$DEMPEL <-as.character(species_bio_by_area$DEMPEL)
       species_bio_by_area$DEMPEL[species_bio_by_area$DEMPEL=="Demersal"] <- "DEM"
       species_bio_by_area$DEMPEL[species_bio_by_area$DEMPEL=="Pelagic"] <- "PEL"
       species_bio_by_area$DEMPEL <-as.character(species_bio_by_area$DEMPEL)
       #Trophic Level FW4
       if(MEANTLs){
         if(substr(survey,1,2) == "GN") species_bio_by_area <- merge(x=species_bio_by_area,y=TLnorth,by="SpeciesSciName",all.x=T,all.y=F)
         if(substr(survey,1,2) %in% c("CS","BB")) species_bio_by_area <- merge(x=species_bio_by_area,y=TLceltic,by="SpeciesSciName",all.x=T,all.y=F)
         if(substr(survey,1,2) == "WA") MEANTLs <- F
       }
       #save a copy 'species_bio_by_area_DEMPEL' so can loop through DEM or PEL etc
       species_bio_by_area_DEMPEL <- species_bio_by_area
     
       #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       # record sampling effort for indicators
       # e.g. num rects sampled by subdiv  sumsampstrat_by_sub   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       if(SAMP_STRAT){
         FACT<-c("Year","sampstrat")
         if(BYSUBDIV) FACT<-c(FACT,"subdiv")
         numsampstrat_by_sea <- tapply.ID(df=species_bio_by_area, datacols=c("CatCatchWgtSwept"), 
                                          factorcols=FACT, sum,c("CatCatchWgtSwept"))
         numsampstrat_by_sea$numsampstrat <- 1
         if(BYSUBDIV){ 
           numsampstrat_by_sub <- tapply.ID(df=numsampstrat_by_sea, datacols=c("numsampstrat"), factorcols=c("Year","subdiv"), sum,c("numsampstrat"))
           sumsampstrat_by_sub <- xtabs(numsampstrat ~ Year + subdiv, numsampstrat_by_sub)
         }
         numsampstrat_by_sea <- tapply.ID(df=numsampstrat_by_sea, datacols=c("numsampstrat"), factorcols=c("Year"), sum,c("numsampstrat"))
         if(BYSUBDIV){ numsampstrat_by_sea <- cbind(sumsampstrat_by_sub,sea=numsampstrat_by_sea[,1])
         rm(sumsampstrat_by_sub)}
         if(WRITE) write.table(numsampstrat_by_sea,paste(FILENAM,"num_rects_sampled_BY_reg_yr.csv",sep="_"),row.names =T,sep=',')
       }
       if(BYSUBDIV){ #if no SAMP_STRAT and sampstrat=NA, then above gives same as this
         num_by_sub <- tapply.ID(df=species_bio_by_area, datacols=c("CatCatchWgtSwept"), 
                                 factorcols=c("Year","subdiv"), sum,c("CatCatchWgtSwept")) 
         num_by_sub$numsamp <- 1
         sum_by_sub <- xtabs(numsamp ~ Year + subdiv, num_by_sub)
         num_by_sea <- tapply.ID(df=num_by_sub, datacols=c("numsamp"), factorcols=c("Year"), sum,c("numsamp"))
         if(BYSUBDIV) num_by_sea <- cbind(sum_by_sub,sea=num_by_sea[,1])
         rm(sum_by_sub)
         if(WRITE) write.table(num_by_sea,paste(FILENAM,"num_subdiv_sampled_BY_yr.csv",sep="_"),row.names =T,sep=',')
         
         #correction for regional sea sampling area required if missing part of SUBDIV
          #area sampled
          num_by_sub <- merge(x=num_by_sub,y=ATTRIB,by.x="subdiv",by.y="SurvStratum",all=T)
          areasurveyed_by_sub <- tapply.ID(df=num_by_sub, datacols=c("KM2_LAM"), 
                                 factorcols=c("Year"), sum,c("KM2_LAM")) 
          #proportion of regional sea area sampled #ATTRIB_SUBDIV is same as totalarea for GNS 'SAMP_STRAT+BYSUBDIV'
          if(survey %in% c("GNSIntOT1","GNSIntOT3","GNSNetBT3","GNSGerBT3")) {
                   totalarea <- sum(ATTRIB_SUBDIV$KM2_LAM)
          } else { totalarea <- sum(ATTRIB$KM2_LAM) }
          areasurveyed_by_sub$scale <- totalarea/areasurveyed_by_sub$KM2_LAM
          if(length(areasurveyed_by_sub[areasurveyed_by_sub$scale>1,'Year']) >0) print( paste("survey area not fully covered in", areasurveyed_by_sub[areasurveyed_by_sub$scale>1,'Year'] ) )
          write.table(areasurveyed_by_sub,paste(FILENAM,"areasurveyed_by_sub.csv",sep="_"),row.names =F,sep=',')
          #not used further
       }
       
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #LD+catch raised by kmsq for species 
      #summarise by sampstrat+subdiv raised by spatial area in kmsq and then kg->tonnes
      if(BYSUBDIV | SAMP_STRAT){
         #find sampling areas from ATTRIB and merge with species data
         if(SAMP_STRAT){ meanLD_bio_by_area <- merge(x=ATTRIB, y=species_bio_by_area, all.y=TRUE, by=c("sampstrat") ) 
         } else {        meanLD_bio_by_area <- merge(x=ATTRIB, y=species_bio_by_area, all.y=TRUE, by.x=c("SurvStratum") , by.y=c("subdiv") ) 
                         if( any(names(meanLD_bio_by_area) %in% "SurvStratum") & !any(names(meanLD_bio_by_area) %in% "subdiv") ) meanLD_bio_by_area$subdiv <-meanLD_bio_by_area$SurvStratum
                         #meanLD_bio_by_area <- meanLD_bio_by_area[,-which(names(meanLD_bio_by_area) == "SurvStratum")] #rm duplicate col
         }
         #raise by area of lowest resolution of sampling strategy/subdiv and return in tonnes
         meanLD_bio_by_area$CatCatchWgtSwept <- meanLD_bio_by_area$CatCatchWgtSwept*meanLD_bio_by_area$KM2_LAM/1000
         
         #output
         if(SAMP_STRAT & WRITE_LDs & !BOOTSTRAP) write.csv(meanLD_bio_by_area,file = paste(paste(FILENAM,sep='_'),"LD_tonnes_Year_W.by.sampstrat.csv",sep=".") )
         if(BYSUBDIV & !SAMP_STRAT & WRITE_LDs & !BOOTSTRAP) write.csv(meanLD_bio_by_area,file = paste(paste(FILENAM,sep='_'),"LD_tonnes_Year_W.by.subdiv.csv",sep=".") )
          
         #second level raising if both levels applied i.e. c("GNSGerBT3","GNSNetBT3","GNSIntOT1","GNSIntOT3")
         if(survey %in% c("GNSIntOT1","GNSIntOT3","GNSNetBT3","GNSGerBT3")){ #if(BYSUBDIV & SAMP_STRAT){
           #use catches scaled by size of grid (rects not constant over sea area)
           # and scale to SUBDIV (beware GNSGerBT3 only sampled a small part of NE so should not do this)
           if(!survey %in% c("GNSIntOT1","GNSIntOT3","GNSNetBT3","GNSGerBT3")) print(paste(survey,"survey does not have two level stratification"))
           #work out value to scale up subdiv by
           #lose species and length (otherwise inflate sum of areas)
           area_by_subdiv <- tapply.ID(df=species_bio_by_area, datacols=c("CatCatchWgtSwept"), factorcols=c("sampstrat","Year","subdiv"), sum,c("CatCatchWgtSweptsum"))  
           #merge in area for sample coverage
           area_by_subdiv <- merge(x=ATTRIB, y=area_by_subdiv, all.y=TRUE, by=c("sampstrat") ) 
           #sum area by subdiv sampled
           area_by_subdiv <- tapply.ID(df=area_by_subdiv, datacols=c("KM2_LAM"), factorcols=c("Year","subdiv"), sum,c("KM2_LAMsum"))  
           #compare to area of subdivision for survey (all years)
           area_by_subdiv <- merge(x=ATTRIB_SUBDIV, y=area_by_subdiv, all.y=TRUE, by.x=c("SurvStratum") , by.y=c("subdiv") ) 
           #ratio to scale up to subdiv estimate
           area_by_subdiv$scale <- area_by_subdiv$KM2_LAM / area_by_subdiv$KM2_LAMsum
           write.csv(area_by_subdiv, file = paste(FILENAM,"_area_by_subdiv.csv",sep=''))
           #big raising factors?
           #area_by_subdiv[area_by_subdiv$scale>2,]
           #now sum catch by sampstrat to subdiv area and scale for missing area (each year)
           FACT <- c("Year","FishLength_cm","SpeciesSciName","subdiv","DEMPEL")
           meanLD_bio_by_subdiv <- tapply.ID(df=meanLD_bio_by_area, datacols=c("CatCatchWgtSwept"), factorcols=FACT, sum,c("CatCatchWgtSwept"))  
           meanLD_bio_by_subdiv <- merge(x=area_by_subdiv[,c("SurvStratum","Year","scale")], y=meanLD_bio_by_subdiv, all.y=TRUE, by.x=c("SurvStratum","Year") , by.y=c("subdiv","Year") )  
           meanLD_bio_by_subdiv$CatCatchWgtSwept <- meanLD_bio_by_subdiv$scale*meanLD_bio_by_subdiv$CatCatchWgtSwept
    
           if(WRITE_LDs & !BOOTSTRAP) write.csv(meanLD_bio_by_subdiv,file = paste(paste(FILENAM,sep='_'), "LD_tonnes_Year_W.by.subdiv.csv",sep=".") )
           #lost DEMPEL, MAXL and TL again
          }
          
         #if missing completely a survey stratum from sampling will have underestimate here
         #mean length distribution by species over all sampling strat
         meanLD_bio_by_areaNOYEAR <- tapply.ID(df=meanLD_bio_by_area, datacols=c("CatCatchWgtSwept"), factorcols= c("FishLength_cm","SpeciesSciName"), mean,c("CatCatchWgtSwept"))  
         if(WRITE_LDs & !BOOTSTRAP) write.csv(meanLD_bio_by_areaNOYEAR,file = paste(paste(FILENAM,sep='_'),"LD_tonnes_YEARave.csv",sep=".") )
   
         #third level biomass to regional sea (or survey extent i.e. coverage of 'subdiv' in year)
         for(SP in SPECIES){
           if(SP=="ALL") SP<-c("DEM","PEL")
           if(survey %in% c("GNSIntOT1","GNSIntOT3","GNSNetBT3","GNSGerBT3")){
             bio_spp_subdiv <- tapply.ID(df=meanLD_bio_by_subdiv[meanLD_bio_by_subdiv$DEMPEL %in% SP,], datacols=c("CatCatchWgtSwept"), factorcols=c("Year","SurvStratum","SpeciesSciName"), sum,c("CatCatchWgtSwept"))  
             bio_spp_area <- tapply.ID(df=meanLD_bio_by_subdiv[meanLD_bio_by_subdiv$DEMPEL %in% SP,], datacols=c("CatCatchWgtSwept"), factorcols=c("Year","SpeciesSciName"), sum,c("CatCatchWgtSwept"))  
             bio_by_subdiv <- tapply.ID(df=meanLD_bio_by_subdiv[meanLD_bio_by_subdiv$DEMPEL %in% SP,], datacols=c("CatCatchWgtSwept"), factorcols=c("Year","SurvStratum"), sum,c("CatCatchWgtSwept"))  
             bio_by_area <- tapply.ID(df=meanLD_bio_by_subdiv[meanLD_bio_by_subdiv$DEMPEL %in% SP,], datacols=c("CatCatchWgtSwept"), factorcols=c("Year"), sum,c("CatCatchWgtSwept"))  
           } else {#one or other of BYSUBDIV | SAMP_STRAT
             if(BYSUBDIV) bio_spp_subdiv <- tapply.ID(df=meanLD_bio_by_area[meanLD_bio_by_area$DEMPEL %in% SP,], datacols=c("CatCatchWgtSwept"), factorcols=c("Year","SurvStratum","SpeciesSciName"), sum,c("CatCatchWgtSwept"))  
             if(BYSUBDIV) bio_by_subdiv <- tapply.ID(df=meanLD_bio_by_area[meanLD_bio_by_area$DEMPEL %in% SP,], datacols=c("CatCatchWgtSwept"), factorcols=c("Year","SurvStratum"), sum,c("CatCatchWgtSwept"))  
             bio_spp_area <- tapply.ID(df=meanLD_bio_by_area[meanLD_bio_by_area$DEMPEL %in% SP,], datacols=c("CatCatchWgtSwept"), factorcols=c("Year","SpeciesSciName"), sum,c("CatCatchWgtSwept"))  
             bio_by_area <- tapply.ID(df=meanLD_bio_by_area[meanLD_bio_by_area$DEMPEL %in% SP,], datacols=c("CatCatchWgtSwept"), factorcols=c("Year"), sum,c("CatCatchWgtSwept"))  
           }
           if(length(SP)==2) SP<-"ALL"
           if(WRITE & !BOOTSTRAP){#all
             if(BYSUBDIV) write.csv(bio_spp_subdiv,file = paste(paste(FILENAM,SP,sep=''), "_Surv_biotonnes_AllSpecies_subdivYear.csv",sep="") )
             if(BYSUBDIV) write.csv(bio_by_subdiv,file = paste(paste(FILENAM,SP,sep=''), "_Surv_biotonnes_subdivYear.csv",sep="") )
             write.csv(bio_spp_area,file = paste(paste(FILENAM,SP,sep=''), "_Surv_biotonnes_AllSpecies_Year.csv",sep="") )
             write.csv(bio_by_area,file = paste(paste(FILENAM,SP,sep=''), "_Surv_biotonnes_Year.csv",sep="") )
             #plot biomass by area once
             windows(width=8, height=4)
             plot(bio_by_area$Year, bio_by_area$CatCatchWgtSwept/1000,type="l",
                  xlab="Year",ylab="Surveyed biomass (kt)",main=paste(survey,SP) )
             points(bio_by_area$Year, bio_by_area$CatCatchWgtSwept/1000,pch=19)
             savePlot(filename= paste(FILENAM,SP,"BIO.bmp",sep='_'),type="bmp")
             dev.off()
           }#plot(bio_by_area[bio_by_area$SpeciesSciName=="Clupea harengus",2:1])
         
           #plot(bio_by_area[bio_by_area$SpeciesSciName=="Clupea harengus",2], bio_by_area[bio_by_area$SpeciesSciName=="Clupea harengus",1]/2000,main="000t, assume half spawning female",type='b')
           
           if(BYSUBDIV){
           bio_by_subdiv$SurvStratumName <- as.factor(bio_by_subdiv$SurvStratum)
           NL<- nlevels(bio_by_subdiv$SurvStratumName)
           PLOTCOLN<- ceiling(sqrt(NL))
           PLOTROWNbio <- ifelse(NL==2,1,PLOTCOLN)
           PLOTROWNbio <- ifelse(NL==5 | NL==6,2,PLOTCOLN)
           #windows(width=8*PLOTCOLN, height=4*PLOTROWNbio)
           
           bmp(filename= paste(FILENAM,SP,"BIOstrata.bmp",sep='_'))
           xy<- xyplot(data=bio_by_subdiv, CatCatchWgtSwept~Year | SurvStratumName,type="b",
                xlab="Year",ylab="Surveyed biomass (tonnes)",main=paste(survey,SP) )
           print(xy)
           dev.off()
           bio_by_subdiv <- bio_by_subdiv[,-which(names(bio_by_subdiv) == "SurvStratumName")]
           }
           
         } #end species loop
       } #END if(BYSUBDIV | SAMP_STRAT)
       
     
        
     #indicators by species dem pel groups    h(species_bioL_by_area_DEMPEL[species_bioL_by_area_DEMPEL$DEMPEL=='PEL',])
     FILENAM_DEMPEL <- FILENAM # copy as overwrite later
     FILENAM_DEMPEL2 <- FILENAM2 # copy as overwrite later
     
     for(SP in SPECIES){ # SP<-"PEL" 
       print(SP)
       FILENAM <-  paste(FILENAM_DEMPEL,SP,sep='') 
       FILENAM2 <-  paste(FILENAM_DEMPEL2,SP,sep='') 
       #meanLD_bio_by_area# is species_bio_by_area_DEMPEL but raised by area to lowest sampstrat
       if(BYSUBDIV | SAMP_STRAT){
         if(SP == "ALL"){ species_bio_by_area <- meanLD_bio_by_area
         } else { species_bio_by_area <- meanLD_bio_by_area[meanLD_bio_by_area$DEMPEL==SP,] }
       } else { #not raised by area above
         if(SP == "ALL"){ species_bio_by_area <- species_bio_by_area_DEMPEL
        } else { species_bio_by_area <- species_bio_by_area_DEMPEL[species_bio_by_area_DEMPEL$DEMPEL==SP,]}
       }
       if(nrow(species_bio_by_area)==0){ print(paste("no species_bio_by_area data for ",SP,sep=''));  break}
       # include correction for area of strata here so have CPUE_estimates * area of sampstrat (or subdiv if lowest level)
       if(survey %in% c("GNSIntOT1","GNSIntOT3","GNSNetBT3","GNSGerBT3")){ # merge in scaling factor
         species_bio_by_area <- merge(x=species_bio_by_area, y=area_by_subdiv, by.x = c("Year","subdiv"),by.y = c("Year","SurvStratum"),all.x=T)
         species_bio_by_area$CatCatchWgtSwept <- species_bio_by_area$CatCatchWgtSwept*species_bio_by_area$scale
        }
         # corrected for any change in sampling between subdiv, but not scaled up to include missing subdiv
        
        #Large Fish Indicator
        if(!SAMP_STRAT) numsampstrat_by_sea <- 0*numhaulsyr
        if(LFI){   IND_LFI <- INDfn_LFI( species_bio_by_area=species_bio_by_area, SP=SP,
                                      numhaulsyr=numhaulsyr,numsampstrat_by_sea=numsampstrat_by_sea, WRITE=WRITE, FILENAM=FILENAM,FILENAM2=FILENAM2,BYSUBDIV=BYSUBDIV, BYSAMPSTRAT=SAMP_STRAT)
        } else { IND_LFI <-NULL }
        #Mean Length cm by sampstrata and year 
        if(MeanL) IND_MeanL <- INDfn_MeanL(species_bio_by_area=species_bio_by_area, WRITE=WRITE, FILENAM=FILENAM,SAMP_STRAT=SAMP_STRAT,BYSUBDIV=BYSUBDIV,SP=SP)
        
        #MaxL by rectangle and year 
        if(MaxL) IND_MaxL <- INDfn_MaxL(species_bio_by_area=species_bio_by_area, WRITE=WRITE, FILENAM=FILENAM,SAMP_STRAT=SAMP_STRAT,BYSUBDIV=BYSUBDIV,SP=SP)
       
        #TL by rectangle and year 
        if(MEANTL) IND_MeanTL <- INDfn_MeanTL(species_bio_by_area=species_bio_by_area, WRITE=WRITE, FILENAM=FILENAM,SAMP_STRAT=SAMP_STRAT,BYSUBDIV=BYSUBDIV,SP=SP)
        
        #Geometric mean length (Typical Length cm)
        # Weight in kg raised to 60 min haul 
        # log[ length(cm) ]
        if(TyL_GeoM) IND_TyL_GeoM <-INDfn_TyL_GeoM(species_bio_by_area=species_bio_by_area, WRITE=WRITE, FILENAM=FILENAM,FILENAM2=FILENAM2,SAMP_STRAT=SAMP_STRAT,BYSUBDIV=BYSUBDIV,SP=SP)
        
       #noddy way to make a list with all indicators for pelagic and demersals
        if(SP=="ALL"){      
          if(TyL_GeoM) TyL.cm.sea_all <- IND_TyL_GeoM;
          if(MeanL) FishLength_cmsea_all<-IND_MeanL; 
          if(MaxL)  MaxLsea_all<-IND_MaxL; 
          if(MEANTL) TLsea_all<-IND_MeanTL; 
          if(LFI & !is.null(IND_LFI[[1]]) ){
            LFIbind_all <-   IND_LFI[[1]]       #all years
            rownames(LFIbind_all) <- LFIbind_all$Year
            # add ncol(LFI_by_sub_all) cols
            if(!is.null(IND_LFI[[2]])){ 
              LFI_by_sub_all <- IND_LFI[[2]]; 
              for(i in 1:ncol(LFI_by_sub_all)) LFIbind_all <- cbind(LFIbind_all,NA)
              names(LFIbind_all)[(ncol(LFIbind_all)-ncol(LFI_by_sub_all)+1):ncol(LFIbind_all)]  <- colnames(LFI_by_sub_all) 
              
              LFIbind_all[ rownames(LFIbind_all) %in% rownames(LFI_by_sub_all),
                           (ncol(LFIbind_all)-ncol(LFI_by_sub_all)+1):ncol(LFIbind_all)] <- (LFI_by_sub_all)
            }
          }
        }
        if(SP=="DEM"){      
          if(TyL_GeoM) TyL.cm.sea_dem <- IND_TyL_GeoM;
          if(MeanL) FishLength_cmsea_dem<-IND_MeanL; 
          if(MaxL)  MaxLsea_dem<-IND_MaxL; 
          if(MEANTL) TLsea_dem<-IND_MeanTL; 
          if(LFI & !is.null(IND_LFI[[1]]) ){
            LFIbind_dem <-   IND_LFI[[1]]       #all years
            rownames(LFIbind_dem) <- LFIbind_dem$Year
            # add ncol(LFI_by_sub_dem) cols
            if(!is.null(IND_LFI[[2]])){ 
              LFI_by_sub_dem <- IND_LFI[[2]]; 
              for(i in 1:ncol(LFI_by_sub_dem)) LFIbind_dem <- cbind(LFIbind_dem,NA)
              names(LFIbind_dem)[(ncol(LFIbind_dem)-ncol(LFI_by_sub_dem)+1):ncol(LFIbind_dem)]  <- colnames(LFI_by_sub_dem) 
              
              LFIbind_dem[ rownames(LFIbind_dem) %in% rownames(LFI_by_sub_dem),
                           (ncol(LFIbind_dem)-ncol(LFI_by_sub_dem)+1):ncol(LFIbind_dem)] <- (LFI_by_sub_dem)
            }
          }
        }
        if(SP=="PEL"){
                    
        if(TyL_GeoM)    TyL.cm.sea_pel<-IND_TyL_GeoM
        if(MeanL)  FishLength_cmsea_pel<-IND_MeanL
        if(MaxL)   MaxLsea_pel <- IND_MaxL
        if(MEANTL) TLsea_pel   <- IND_MeanTL
        if(LFI & !is.null(IND_LFI[[1]]) ){
          LFIbind_pel <-  IND_LFI[[1]]     #all years
          rownames(LFIbind_pel) <- LFIbind_pel$Year
          # add ncol(LFI_by_sub_dem) cols
          if(!is.null(IND_LFI[[2]])){ 
            LFI_by_sub_pel<-IND_LFI[[2]]; 
            for(i in 1:ncol(LFI_by_sub_pel)) LFIbind_pel <- cbind(LFIbind_pel,NA)
            names(LFIbind_pel)[(ncol(LFIbind_pel)-ncol(LFI_by_sub_pel)+1):ncol(LFIbind_pel)]  <- colnames(LFI_by_sub_pel)              
            LFIbind_pel[rownames(LFIbind_pel) %in% rownames(LFI_by_sub_pel),(ncol(LFIbind_pel)-ncol(LFI_by_sub_pel)+1):ncol(LFIbind_pel)] <- (LFI_by_sub_pel)              
          }
        }
      }
     }#species set loop 
     
      if(!BOOT) return(list(                                               
                    LFI_by_sub_all = LFIbind_all, TyL.cm.sea_all = TyL.cm.sea_all, FishLength_cmsea_all =FishLength_cmsea_all, MaxLsea_all =MaxLsea_all, TLsea_all =TLsea_all 
                    ,
                    LFI_by_sub_dem = LFIbind_dem, TyL.cm.sea_dem = TyL.cm.sea_dem, FishLength_cmsea_dem =FishLength_cmsea_dem, MaxLsea_dem =MaxLsea_dem, TLsea_dem =TLsea_dem 
                    ,
                    LFI_by_sub_pel = LFIbind_pel, TyL.cm.sea_pel = TyL.cm.sea_pel, FishLength_cmsea_pel =FishLength_cmsea_pel, MaxLsea_pel =MaxLsea_pel, TLsea_pel =TLsea_pel
                    ,
                    species_bio_by_area=species_bio_by_area_DEMPEL, 
                    numhauls=numhauls, numhaulsyr=numhaulsyr, numhaulsBYsampstrat=numhaulsBYsampstrat, numhaulsBYsubdiv=numhaulsBYsubdiv) )
     #dont save everything in bootstrap
      if(BOOT) return(list(                    
                    LFI_by_sub_all = LFIbind_all, TyL.cm.sea_all = TyL.cm.sea_all, FishLength_cmsea_all =FishLength_cmsea_all, MaxLsea_all =MaxLsea_all, TLsea_all =TLsea_all 
                    ,
                    LFI_by_sub_dem = LFIbind_dem, TyL.cm.sea_dem = TyL.cm.sea_dem, FishLength_cmsea_dem =FishLength_cmsea_dem, MaxLsea_dem =MaxLsea_dem, TLsea_dem =TLsea_dem
                    ,
                    LFI_by_sub_pel = LFIbind_pel, TyL.cm.sea_pel = TyL.cm.sea_pel, FishLength_cmsea_pel =FishLength_cmsea_pel, MaxLsea_pel =MaxLsea_pel, TLsea_pel =TLsea_pel 
                    ) )
 }
