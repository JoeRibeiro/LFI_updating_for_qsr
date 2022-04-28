library(rioja); #library(palaeo); # rioja includes updated fns from palaeo package with bug fixes
library(strucchange); library(ade4);  library(vegan)

REGIME<-function(dat,H=0.3,tit,BREAKS=-11,LINEAR=2,PLOT_F=F,PLOT_ALT=F,COMP=NA,ADD=F){
  # dat may be raw ts or pca/prcomp
  # tit is title to add to plot
  # H is minimum number of points in a regime (or a proportion of ts) - will be altered and warning given if too small

  # funtion will fit segmented regression lines (const or linear) to data if breaks are identified
    # if breaks == -11, will auto select breakpoints using breaks fn for ALL fits
    # if linear == 2, will select between const or linear fits using BIC
    # if breaks == 0 or 1 then this is the maximum number of breaks allowed
    # if breaks == -1, will auto select breakpoints using breaks fn for chosen plot
  # if plot_ALT == T, will plot const/linear fit not selected manually on input or auto by BIC
  # note will not plot alternative if best fit has no breaks
  #if(!ADD) windows(width=25,height=15);par(mfrow=c(2,3),oma=c(3,1,1,1));
  fmL<-fmC<-NULL
  if(class(dat)!="data.frame"){ pca<-TRUE;} else {pca<-FALSE;} # catch prcomp or other pca types this way
  
  if(pca==T){ 
    THRESHDAT<-predict(dat); 
    sum_m<-summary(dat);# summary for titles      if(any(is.na(COMP)==F)){THRESHDAT<-THRESHDAT[,COMP];} 
    PLOTF<-PLOT_F; # plot F statistics
    if(ADD==F){
      if(PLOTF==F){windows(width=55,height=80); par(mfrow=c(3,1),oma=c(3,1,1,1));
      } else { windows(width=110,height=80);par(mfrow=c(3,2),oma=c(3,1,1,1)); }
    }
  } else {
    PLOTF<-FALSE; #too many plots if using raw data
    THRESHDAT<-na.trim(dat);   windows(width=110,height=80);par(mfrow=c(3,3),oma=c(3,1,1,1)); 
  }
  BREAKSALL<-NULL; allci<-NULL;
  Year<-as.numeric(rownames(THRESHDAT));
  if(pca==T){ n<-min(3,ncol(THRESHDAT))  } else { n<-ncol(THRESHDAT)  }

  #option to choose only a subset of components to model/plot
  if(any(is.na(COMP)==F)){RANGE<-COMP;} else {RANGE<-1:n;}      
  for(x in RANGE){# 1:n
    if(pca==F) print(names(dat)[x]);#  x<-1
    if(sum(x==seq(10,100,9))==1){
      title(tit,outer=TRUE,line=-.5);
      mtext(paste('h =',H,sep=' '),side=1,outer=T);
      savePlot(filename = paste( paste(tit,(x-1),'h',H, sep="_"),'jpg', sep="."), type = "jpg");
      windows(width=110,height=80);par(mfrow=c(3,3),oma=c(3,1,1,1));}
     dattest<-data.frame(THRESHDAT[,x],Year); dattest<-na.trim(dattest)

    ## F statistics test NULL hyp: no breaks vs Alt hyp: 1 breakpoint present
    tsTHRESHDAT <- ts(dattest$THRESHDAT,start=dattest$Year[1])
    if(LINEAR==2){
       ## TRY linear model first
        test1<-sctest(na.trim(THRESHDAT[,x]) ~ Year, type = "sup",data=dattest);
        fs <- Fstats(tsTHRESHDAT ~ Year, from = H, data=dattest);
        bp1<-breakpoints(tsTHRESHDAT ~ Year,data=dattest, h = H);
        if(length(bp1$breakpoints)>1 || is.na(bp1$breakpoints)==F) ci <- try(confint(bp1), silent =T)
        if(PLOTF==T) { ry1<-range(fs$Fstats); fs1<-fs;      }
       ## which is equivalent to segmenting the regression via
        fm10 <- lm(tsTHRESHDAT ~ Year,data=dattest);                  ## linear model
        if(length(bp1$breakpoints)>1 || is.na(bp1$breakpoints)==F){
         fm11 <- lm(tsTHRESHDAT ~ Year*breakfactor(bp1),data=dattest);## linear model
         anova(fm10, fm11)
        } else { fm11 <- fm10; }
       ## and constant model
        test<-sctest(na.trim(THRESHDAT[,x]) ~ 1, type = "sup",data=dattest);
        fs <- Fstats(tsTHRESHDAT ~ 1,from = H,data=dattest);
        bp<-breakpoints(tsTHRESHDAT ~ 1,data=dattest, h = H);
        if(length(bp$breakpoints)>1 || is.na(bp$breakpoints)==F) ci <- try(confint(bp), silent =T)
        if(PLOTF==T){  ry2<-range(fs$Fstats,ry1); #,xlim=c(Year[1],Year[length(Year)])
                         plot(fs1,ylim=c(ry2),main='',xlab='',col='orange');
                         title('linear (orange) and const (blue)',line=0.5)
                         lines(fs,col=4); #mtext("const",side=4,line=2,col=4); box();          
                     if(length(bp$breakpoints)>1 || is.na(bp$breakpoints)==F){lines(bp,col=4);}
                     if(length(bp1$breakpoints)>1 || is.na(bp1$breakpoints)==F){lines(bp1,col='orange');}
                    } 
       ## which is equivalent to segmenting the regression via
        fm0 <- lm(tsTHRESHDAT ~ 1,data=dattest);                     #constant model
        if(length(bp$breakpoints)>1 || is.na(bp$breakpoints)==F){    #constant model
         fm1 <- lm(tsTHRESHDAT ~ breakfactor(bp),data=dattest); anova(fm0, fm1)
         BICM1<-round(BIC(fm1),2);
        } else { fm1 <- fm0; BICM1<-NA;}
        BEST<-rownames(BIC(fm11,fm10,fm1,fm0)[ BIC(fm11,fm10,fm1,fm0)[,2]==min(BIC(fm11,fm10,fm1,fm0)[,2]),])
        if(length(BEST)>1) BEST<-BEST[1]
                                      
        if(BEST=='fm11'){fm1<-fm11; fm0<-fm10; bp<-bp1; test<-test1; LINEAR_BEST<-TRUE;# need to keep linear form
        } else {LINEAR_BEST<-FALSE;} # simply keep CONST form - no need to overwrite

    } else {
      if(LINEAR==TRUE){ LINEAR_BEST<-TRUE;   # LINEAR==1==TRUE
       ## sup test on NSea climate data
        test<-sctest(na.trim(THRESHDAT[,x]) ~ Year, type = "sup",data=dattest);
        fs <- Fstats(tsTHRESHDAT ~ Year,from = H,data=dattest);
        bp<-breakpoints(tsTHRESHDAT ~ Year,data=dattest, h = H);
        if(length(bp$breakpoints)>1 || is.na(bp$breakpoints)==F) ci <- try(confint(bp), silent =T)
        if(PLOTF==T){  plot(fs);
                       if(length(bp$breakpoints)>1 || is.na(bp$breakpoints)==F){ lines(bp); 
                        if(any(is.na(ci$confint))==F){ lines(ci);}
                       } 
                    }
        ## which is equivalent to segmenting the regression via
        fm0 <- lm(tsTHRESHDAT ~ Year,data=dattest);
          if(length(bp$breakpoints)>1 || is.na(bp$breakpoints)==F){
           fm1 <- lm(tsTHRESHDAT ~ Year*breakfactor(bp),data=dattest); anova(fm0, fm1)
            BICM1<-round(BIC(fm1),2);
          } else { BICM1<-NA;}
      } else { LINEAR_BEST<-FALSE;       # LINEAR==0
        test<-sctest(na.trim(THRESHDAT[,x]) ~ 1, type = "sup",data=dattest);
        fs <- Fstats(tsTHRESHDAT ~ 1,from = H,data=dattest);
        bp<-breakpoints(tsTHRESHDAT ~ 1,data=dattest, h = H);
        if(length(bp$breakpoints)>1 || is.na(bp$breakpoints)==F) ci <-  try(confint(bp), silent =T)
        if(PLOTF==T){  plot(fs);
                       if(length(bp$breakpoints)>1 || is.na(bp$breakpoints)==F){ lines(bp); 
                        if(any(is.na(ci$confint))==F){ lines(ci);}
                       } 
                    }
        ## which is equivalent to segmenting the regression via
         fm0 <- lm(tsTHRESHDAT ~ 1,data=dattest);
          if(length(bp$breakpoints)>1 || is.na(bp$breakpoints)==F){
           fm1 <- lm(tsTHRESHDAT ~ breakfactor(bp),data=dattest); anova(fm0, fm1)
           BICM1<-round(BIC(fm1),2)
          } else {BICM1<-NA;}
      }
    }

    #summary(fm0);  summary(fm1)  ## estimates  & ADD!=T
    if(pca==T){ YLAB=paste("PC",x,"var =",100*round(sum_m$importance[2,x],2)[1],'%',sep=' ' )
    } else { YLAB=names(dat)[x] }
    plot(tsTHRESHDAT,type='p',xlab='',main='',col='dark grey',ylab=YLAB); abline(h=0);
      #polygon(x=c(1996,1996:2003,2003),y=c(0,dattest$THRESHDAT[14:21],0),col=3)
      COL.MAIN<-1; #if(LINEAR==2 & LINEAR_BEST){ COL.MAIN<-4;} else {COL.MAIN<-1}
      if(LINEAR==2){ title(sub=paste('LINEAR BEST:',LINEAR_BEST));}
      title(col.main=COL.MAIN, paste('BIC M0 =', round(BIC(fm0),2),' M1 =', BICM1,
            '  sup.F =',round(test$statistic,1),'  p =', round(test$p.value,3),sep=' '), cex.main=1 ,line=1.0)#,col.main='dark grey'

      
  #  if(names(THRESHDAT)[x]=="StatlantCDab" || names(THRESHDAT)[x]=="StatlantCLong.rough.dab" || names(THRESHDAT)[x]=="StatlantCThornback...Spotted.ray") next;
  #  if(names(THRESHDAT)[x]=="StatlantCMackerel" || names(THRESHDAT)[x]=="StatlantCMiscellaneous.filter.feeding.pelagic.fish") next;


  ### below code for selecting breakpoints manually
    if(length(bp$breakpoints)>1 || is.na(bp$breakpoints)==F){ #if breaks found
        if(BREAKS==22 || BREAKS==-11){
        } else {
          ANSWER <- readline("Max number of breakpoints? (or enter -1 for auto) ")
          BREAKS<-as.numeric(ANSWER); BREAKSALL<-cbind(BREAKSALL,BREAKS)
        }

        REGIMEEND<- length(dattest$Year);
        lines(ts(fitted(fm1)[1:bp$breakpoints[1]], start = dattest$Year[1]), col = 'grey',lwd=1,lty=1)

        for(i in 1:length(bp$breakpoints) ){
         if(i==length(bp$breakpoints)){
          lines(ts(fitted(fm1)[(bp$breakpoints[i]+1):REGIMEEND], start=dattest$Year[1+bp$breakpoints[i]]),col='grey',lwd=1,lty=1)
         } else {
          lines(ts(fitted(fm1)[(bp$breakpoints[i]+1):bp$breakpoints[i+1]], start=dattest$Year[1+bp$breakpoints[i]]),col='grey',lwd=1,lty=1)
         }
        }
        if(BREAKS==0){
          lines(ts(fitted(fm0), start = dattest$Year[1]), col = 1,lwd=2,lty=1);
        } else {
          if(BREAKS==-1 || BREAKS==-11){
            bpL<-breakpoints(tsTHRESHDAT ~ Year,data=dattest,h=H);
            bpC<-breakpoints(tsTHRESHDAT ~ 1,data=dattest,h=H);
          } else {
            if(BREAKS==22){ BREAKSt<-2;} else {BREAKSt<-BREAKS}
            bpL<-breakpoints(tsTHRESHDAT ~ Year,data=dattest,h=H,breaks=BREAKSt);
            bpC<-breakpoints(tsTHRESHDAT ~ 1,data=dattest,h=H,breaks=BREAKSt);
          }
  
          if(length(bpL$breakpoints)>1 || is.na(bpL$breakpoints)==F) fmL <- lm(tsTHRESHDAT ~ Year*breakfactor(bpL),data=dattest)
          if(length(bpC$breakpoints)>1 || is.na(bpC$breakpoints)==F) fmC <- lm(tsTHRESHDAT ~ breakfactor(bpC),data=dattest)
          bp2<-NULL;
          if(LINEAR_BEST==TRUE){
            if(length(bpL$breakpoints)>1 || is.na(bpL$breakpoints)==F){ fm1 <- fmL; bp <- bpL; }
            if(length(bpC$breakpoints)>1 || is.na(bpC$breakpoints)==F){ fm2 <- fmC; bp2 <- bpC;}
          } else {
            if(length(bpL$breakpoints)>1 || is.na(bpC$breakpoints)==F){ fm1 <- fmC; bp <- bpC;} 
            if(length(bpC$breakpoints)>1 || is.na(bpL$breakpoints)==F){ fm2 <- fmL; bp2 <- bpL;}
          }
          if(is.null(bp2)==F && (length(bp2$breakpoints)>1 || is.na(bp2$breakpoints)==F)){ 
            if(PLOT_ALT==T){ lines(bp2,col=4,lty=3); ci2 <-  try(confint(bp2), silent =T);             
              if(any(is.na(ci2$confint))==F){ lines(ci2,col=3,lwd=1,lty=2);}
            }
          }
          lines(bp,col='grey',lty=3); #plot best fit
          ci <-  try(confint(bp), silent =T); 
          if(any(is.na(ci$confint))==F){ lines(ci,col=1,lwd=2,lty=1);}
          #if(BREAKS>-1) title( paste('BIC M1 =', round(BIC(fm1),2),sep=' '), line=1.5,col.main=1,font.main=3)
  
          REGIMEEND<- length(dattest$Year); #plot chosen/best fit
          lines(ts(fitted(fm1)[1:bp$breakpoints[1]], start = dattest$Year[1]), col = 1,lwd=2,lty=1)
          for(i in 1:length(bp$breakpoints) ){
           if(i==length(bp$breakpoints)){
            lines(ts(fitted(fm1)[(bp$breakpoints[i]+1):REGIMEEND], start=dattest$Year[1+bp$breakpoints[i]]),col=1,lwd=2,lty=1)
           } else {
            lines(ts(fitted(fm1)[(bp$breakpoints[i]+1):bp$breakpoints[i+1]], start=dattest$Year[1+bp$breakpoints[i]]),col=1,lwd=2,lty=1)
           }
          }
          #### alternative ####  bp2 & fm2 are for const or linear model (whichever not best)
           if(PLOT_ALT==T){                       
            if(is.null(bp2)==F && (length(bp2$breakpoints)>1 || is.na(bp2$breakpoints)==F)){
              lines(ts(fitted(fm2)[1:bp2$breakpoints[1]], start = dattest$Year[1]), col=3,lwd=1,lty=1)
              for(i in 1:length(bp2$breakpoints) ){
               if(i==length(bp2$breakpoints)){
                lines(ts(fitted(fm2)[(bp2$breakpoints[i]+1):REGIMEEND], start=dattest$Year[1+bp2$breakpoints[i]]),col=3,lwd=1,lty=1)
               } else {
                lines(ts(fitted(fm2)[(bp2$breakpoints[i]+1):bp2$breakpoints[i+1]], start=dattest$Year[1+bp2$breakpoints[i]]),col=3,lwd=1,lty=1)
               }
              }
            } else { # no breakpoints in alternative
              if(LINEAR_BEST==0) lines(Year,fitted(fm10), col=3,lwd=1,lty=1);              
            }# fm10
           }
          #####################     
       allci  <- rbind(allci,cbind( names(THRESHDAT)[x], dattest$Year[1]+ci$confint-1 ))
      }
    } else {    lines(ts(fitted(fm0), start = dattest$Year[1]), col = 1,lwd=2,lty=1);}
    points(tsTHRESHDAT,col='grey',pch=19,cex=1.5)
  }
  if(pca==T & (n>2)){if(ADD!=T) title(paste(tit, ":  PC1-3","var =",100*round(sum_m$importance[3,3],2)[1],'%' ,sep=' '),outer=TRUE,line=-1.5);          
    } else { title(tit,outer=TRUE,line=- .5);}  
  
  if(any(is.na(COMP))){ write.table( allci,file=paste(tit,'h',H,'ci.csv',sep="_") );
                   mtext(paste('h =',H,sep=' '),side=1,outer=T)
                   savePlot(filename = paste( paste(tit,'h',H, sep="_"),'jpg', sep="."), type = "jpg")
     } else { write.table( allci,file=paste(tit,'ci.csv',sep="_"),append=TRUE );  
              mtext(paste('h =',H,sep=' '),side=4,line=1.5)
              savePlot(filename = paste( paste(tit, sep="_"),'jpg', sep="."), type = "jpg")
  }
}

