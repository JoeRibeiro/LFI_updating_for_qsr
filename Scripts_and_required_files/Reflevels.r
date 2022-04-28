#require(strucchange)

REFLVL<-function(dat,H=0.3,COMP=NA){
  # dat is raw ts
  # H is minimum number of points in a regime (or a proportion of ts) - will be altered and warning given if too small
  # funtion will fit segmented regression lines (const) to data if breaks are identified

  windows()
    THRESHDAT<- dat;  BREAKSALL<-NULL; allci<-NULL;
  Year<-as.numeric(rownames(THRESHDAT));
  n<-ncol(THRESHDAT)

  #option to choose only a subset of components to model/plot
  if(any(is.na(COMP)==F)){RANGE<-COMP;} else {RANGE<-1:n;}      
  for(x in RANGE){# 1:n
    print(names(dat)[x]);#  x<-1

    dattest<-data.frame(THRESHDAT[,x],Year); dattest<-na.trim(dattest)

    ## F statistics test NULL hyp: no breaks vs Alt hyp: 1 breakpoint present
    tsTHRESHDAT <- ts(dattest$THRESHDAT,start=dattest$Year[1])

        test<-sctest(na.trim(THRESHDAT[,x]) ~ 1, type = "sup",data=dattest);
        fs <- Fstats(tsTHRESHDAT ~ 1,from = H,data=dattest);
        bp<-breakpoints(tsTHRESHDAT ~ 1,data=dattest, h = H);
        if(length(bp$breakpoints)>1 || is.na(bp$breakpoints)==F) ci <-  try(confint(bp), silent =T)

        ## which is equivalent to segmenting the regression via
         fm0 <- lm(tsTHRESHDAT ~ 1,data=dattest);
          if(length(bp$breakpoints)>1 || is.na(bp$breakpoints)==F){
           fm1 <- lm(tsTHRESHDAT ~ breakfactor(bp),data=dattest); anova(fm0, fm1)
          } else {fm1 <-NA;}

   YLAB=names(dat)[x]
    plot(tsTHRESHDAT,type='p',xlab='',main='',col='dark grey',ylab=YLAB); abline(h=0);
      COL.MAIN<-1;
      title(col.main=COL.MAIN, paste('sup.F =',round(test$statistic,1),'  p =', round(test$p.value,3),sep=' '), cex.main=1 ,line=1.0)
      if(length(bp$breakpoints)==0 || is.na(bp$breakpoints)==T){
        lines(ts(fitted(fm0), start = dattest$Year[1]), col = 'grey',lwd=2,lty=1)
          next;}
        REGIMEEND<- length(dattest$Year);
        lines(ts(fitted(fm1)[1:bp$breakpoints[1]], start = dattest$Year[1]), col = 'grey',lwd=1,lty=1)

        for(i in 1:length(bp$breakpoints) ){
         if(i==length(bp$breakpoints)){
          lines(ts(fitted(fm1)[(bp$breakpoints[i]+1):REGIMEEND], start=dattest$Year[1+bp$breakpoints[i]]),col='grey',lwd=1,lty=1)
         } else {
          lines(ts(fitted(fm1)[(bp$breakpoints[i]+1):bp$breakpoints[i+1]], start=dattest$Year[1+bp$breakpoints[i]]),col='grey',lwd=1,lty=1)
         }
        }
        bpL<-breakpoints(tsTHRESHDAT ~ Year,data=dattest,h=H);
        bpC<-breakpoints(tsTHRESHDAT ~ 1,data=dattest,h=H);

  
          if(length(bpL$breakpoints)>1 || is.na(bpL$breakpoints)==F) fmL <- lm(tsTHRESHDAT ~ Year*breakfactor(bpL),data=dattest)
          if(length(bpC$breakpoints)>1 || is.na(bpC$breakpoints)==F) fmC <- lm(tsTHRESHDAT ~ breakfactor(bpC),data=dattest)
          bp2<-NULL;
          
            if(length(bpL$breakpoints)>1 || is.na(bpC$breakpoints)==F){ fm1 <- fmC; bp <- bpC;} 
            if(length(bpC$breakpoints)>1 || is.na(bpL$breakpoints)==F){ fm2 <- fmL; bp2 <- bpL;}

          lines(bp,col='grey',lty=3); #plot best fit
          ci <-  try(confint(bp), silent =T); 
          if(any(is.na(ci$confint))==F){ lines(ci,col=1,lwd=2,lty=1);}

          REGIMEEND<- length(dattest$Year); #plot chosen/best fit
          lines(ts(fitted(fm1)[1:bp$breakpoints[1]], start = dattest$Year[1]), col = 4,lwd=2,lty=1)
          for(i in 1:length(bp$breakpoints) ){
           if(i==length(bp$breakpoints)){
            lines(ts(fitted(fm1)[(bp$breakpoints[i]+1):REGIMEEND], start=dattest$Year[1+bp$breakpoints[i]]),col=4,lwd=2,lty=1)
           } else {
            lines(ts(fitted(fm1)[(bp$breakpoints[i]+1):bp$breakpoints[i+1]], start=dattest$Year[1+bp$breakpoints[i]]),col=4,lwd=2,lty=1)
           }
          }
       allci  <- rbind(allci,cbind( names(THRESHDAT)[x], dattest$Year[1]+ci$confint-1 ))

  if(any(is.na(COMP))){ write.table( allci,file=paste(YLAB,'h',H,'ci.csv',sep="_") );
                   mtext(paste('h =',H,sep=' '),side=1,outer=T)
                   savePlot(filename = paste( paste(YLAB,'h',H, sep="_"),'png', sep="."), type = "png")
     } else { write.table( allci,file=paste(YLAB,'ci.csv',sep="_"),append=TRUE );
              mtext(paste('h =',H,sep=' '),side=4,line=1.5)
              savePlot(filename = paste( paste(YLAB, sep="_"),'png', sep="."), type = "png")
      }
  }
}

