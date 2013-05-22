##for examining results

setwd("C:/Users/dmcglinn/Documents/Grad Ecology/hump paper/diversity version/resubmit/online supplement")

load("simoutput.RData")

par(mfrow=c(2,3))
hump.plot(output,type='l',avg=T,ci=T,poly=T,cls=c(purples[1],NA,NA),cls.poly=c(purples[2],NA,NA),lwd=lwd,pch=pchs)
hump.plot(output,type='l',avg=T,ci=T,poly=T,cls=c(NA,blues[1],NA),cls.poly=c(NA,blues[2],NA),lwd=lwd,pch=pchs)
hump.plot(output,type='l',avg=T,ci=T,poly=T,cls=c(NA,NA,greens[2]),cls.poly=c(NA,NA,greens[1]),lwd=lwd,pch=pchs)

hump.plot2(output,type='l',avg=T,ci=T,poly=T,cls=c(purples[1],NA,NA),cls.poly=c(purples[2],NA,NA),lwd=lwd,pch=pchs)
hump.plot2(output,type='l',avg=T,ci=T,poly=T,cls=c(NA,blues[1],NA),cls.poly=c(NA,blues[2],NA),lwd=lwd,pch=pchs)
hump.plot2(output,type='l',avg=T,ci=T,poly=T,cls=c(NA,NA,greens[2]),cls.poly=c(NA,NA,greens[1]),lwd=lwd,pch=pchs)


##graphical script for the biomass-richness relationship##
hump.plot2<-function(holder,avg=TRUE,log='',add=FALSE,ci=FALSE,poly=FALSE,type='l',lwd=2,pch=19,cls=palette(),cls.poly=rep('grey',dim(holder)[2]),main='',xlab='biomass',ylab='species richness',allo=8/3){
 ##plots number of species vs. total biomass
 ##FUNCTION ARGUMENTS:
 ##'holder' is any output object from the function "hump.sim"
 ##'avg' if TRUE then the average across all permuations is plotted
 ##'log' may be either '','x','y', or 'xy' which results in a log transformation of neither axis, x only, y only, or both x and y axes respectively
 ##'add' if TRUE then the function plots the results on top of the previous active graphical display
 ##'ci' if TRUE then the 95% confidence intervals are plotted, only works if 'avg' is TRUE
 ##'poly', if TRUE then the confidence intervals are represented as a colored polygon and not simply as lines, only works if 'ci' is TRUE
 ##'type' may be either 'l','p','o', see ?plot.default
 ##'lwd' specifies line width any positive numeric
 ##'pch' specfies point style, see ?points
 ##'cls' specifes the color scheme for the points or lines of the results
 ##'cls.poly' specifies the color scheme for the confidence interval
 ##'main' specifies the main title for the graphic, default is blank
 ##'xlab' label of the x-axis
 ##'ylab' label of the y-axis
 ##'allo' the allometric exponent by which diameter is converted to above-ground biomass
 #####################
 xlims<-c(min(apply(apply(holder[,1,1,1,,,4]^allo,c(1,2),function(x){sum(x,na.rm=TRUE)}),2,mean)),
          max(apply(apply(holder[,dim(holder)[2],dim(holder)[3],dim(holder)[4],,,4]^allo,c(1,2),function(x){sum(x,na.rm=TRUE)}),2,mean)))
 ylims<-c(1,numsp)
 if(avg){ ##if avg is true then averages over the permutations are calculated and ploted
  if(add==FALSE){ ##if add is false then new plot is created
   plot(apply(apply(holder[,1,1,1,,,4]^allo,c(1,2),function(x){sum(x,na.rm=TRUE)}),2,mean),apply(apply(holder[,1,1,1,,,3],c(1,2),function(x){length(unique(x[!is.na(x)]))}),2,mean),type='n',
        xlim=xlims,ylim=ylims,log=log,xlab=xlab,ylab=ylab,main=main)
  }
  if(ci){ ##if ci is true then confidence intervals are added as lines
   for(s in 1:dim(holder)[2]){ ##different spatial distributions
    if(is.na(cls[s])) next 
    for(d in 1:dim(holder)[3]){ ##different starting diameters
     for(n in 1:dim(holder)[4]){ ##differnt starting individual densities
      b.low<-apply(apply(holder[,s,d,n,,,4]^allo,c(1,2),function(x){sum(x,na.rm=TRUE)}),2,function(x){quantile(x,probs=.025)})
      b.high<-apply(apply(holder[,s,d,n,,,4]^allo,c(1,2),function(x){sum(x,na.rm=TRUE)}),2,function(x){quantile(x,probs=.975)})
      s.low<-apply(apply(holder[,s,d,n,,,3],c(1,2),function(x){length(unique(x[!is.na(x)]))}),2,function(x){quantile(x,probs=.025)})
      s.high<-apply(apply(holder[,s,d,n,,,3],c(1,2),function(x){length(unique(x[!is.na(x)]))}),2,function(x){quantile(x,probs=.975)})
      ##necessary to work with logs...
      b.low<-ifelse(b.low==0,b.low+.00001,b.low)
      s.low<-ifelse(s.low==0,s.low+.00001,s.low)
      if(poly){ 
       polygon(c(b.low,rev(b.high)),c(s.low,rev(s.high)),col=cls.poly[s],border=NA)
      }
      else{
       points(b.low,s.low,type='l',col=cls[s],lty=2)
       points(b.high,s.high,type='l',col=cls[s],lty=2)
  }}}}}  
  for(s in 1:dim(holder)[2]){ ##different spatial distributions
   if(is.na(cls[s])) next
   for(d in 1:dim(holder)[3]){ ##different starting diameters
    for(n in 1:dim(holder)[4]){ ##differnt starting individual densities
     points(apply(apply(holder[,s,d,n,,,4]^allo,c(1,2),function(x){sum(x,na.rm=TRUE)}),2,mean),apply(apply(holder[,s,d,n,,,3],c(1,2),function(x){length(unique(x[!is.na(x)]))}),2,mean),type=type,lwd=lwd,col=cls[s],pch=pch)
  }}}
 }
 else{ #avg is FALSE and no averaging is performed
  if(add==FALSE){ ##if add is false then new plot is created
   plot(apply(holder[,1,1,1,,,4]^allo,c(1,2),function(x){sum(x,na.rm=TRUE)}),apply(holder[,1,1,1,,,3],c(1,2),function(x){length(unique(x[!is.na(x)]))}),type='n',
      xlim=xlims,ylim=ylims,log=log,xlab=xlab,ylab=ylab,main=main)
  }
  for(s in 1:dim(holder)[2]){ ##different spatial distributions
   if(is.na(cls[s])) next
   for(d in 1:dim(holder)[3]){ ##different starting diameters
    for(n in 1:dim(holder)[4]){ ##differnt starting individual densities
     points(apply(holder[,s,d,n,,,4]^allo,c(1,2),function(x){sum(x,na.rm=TRUE)}),apply(holder[,s,d,n,,,3],c(1,2),function(x){length(unique(x[!is.na(x)]))}),type=type,lwd=lwd,col=cls[s],pch=pch)
  }}}
 }
 abline(h=numsp,col='grey',lty=2)
}
#########################