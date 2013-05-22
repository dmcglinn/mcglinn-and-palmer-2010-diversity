##for each start diam, each individual start numb,each time val,each indiv, c(xcord,ycord,sp.id,diameter values)
hump.sim<-function(nperm,diams,ns,RAD='even',space='random',intensity=.8,all.space=FALSE,mort.buffer=1,expon=TRUE,overlap=FALSE){
 ##nperm is number of times to run the simulation
 ##diams is the input starting diameters
 ##ns is the input starting individual numbers
 ##space specifies spatial distribution: 'clustered','random','uniform', or 'grid'
 ##intensity specifies how strongly spatial clustering will occur (0-1)
 ##all.space means that all possible spatial distributions should be used
 ##mort.buff is a way to decrease the intensity of the mortaility (values less than one)
 ##expon is whether or not exponential decay of mortality shoudl be used, linear methods are available
 ##overlap specifies whether the degree of competition should be based on the area of overlap  with neighbors (TRUE) 
 ##or the summed diameter or neighbors (FALSE)
 RADs<-c('even','lnorm(0,1)','geometric','lnorm(0,4)','uneven')
 rad.num<-match(RAD,RADs)
 if(is.na(rad.num)){
  print(paste(c('Choose one of the appropriate RADs:',RADs)))
  stop()
 }
 spaces<-c('clustered','random','uniform','grid')
 if(all.space){ ##if all.space=TRUE then results generated for each spatial distribution
  space<-spaces
 }
 for(i in 1:length(space)){
  spa.num<-match(space[i],spaces)
  if(is.na(spa.num)){
   print(paste(c('Choose one or more of the appropriate spatial distributions:',spaces)))
   stop()
 }}
 len.space<-length(space) 
 len.diam<-length(diams)
 len.n<-length(ns)
 holder<-array(NA,dim=c(nperm,len.space,len.diam,len.n,Tim,max(ns),4))
 for(iperm in 1:nperm){
  for(s in 1:len.space){ #spatial distributions
   for(d in 1:len.diam){ #dimeter sizes
    for(n in 1:len.n){ #ninit is initial density
     #assign individuals to grid:
     if(space[s]=='random')
      holder[iperm,s,d,n,1,1:ns[n],1:2]<-runif(2*ns[n])#x&y
     if(space[s]=='clustered'){ 
      ##fraction of ns[n] will be the parents
      nparents<-round((1-intensity)*ns[n])##how many parent points
      parent.cords<-cbind(runif(nparents),runif(nparents))#x&y
      ##now distribute individuals clustered around parents
      icount<-1
      iparent<-1
      while(icount!=(ns[n]+1)){
       if(iparent %% nparents == 0) iparent <- 1
       holder[iperm,s,d,n,1,icount,1:2]<-c(rnorm(1,mean=parent.cords[iparent,1],sd=diams[d]),
                                                    rnorm(1,mean=parent.cords[iparent,2],sd=diams[d]))
       ##the following if statement checks to make sure the point did not fall outside the 1x1 window
       ##if it did fall outside the algo goes to the next parent and tries to place an inidividual there
       if(holder[iperm,s,d,n,1,icount,1]>1|holder[iperm,s,d,n,1,icount,1]<0|
          holder[iperm,s,d,n,1,icount,2]>1|holder[iperm,s,d,n,1,icount,2]<0){
        holder[iperm,s,d,n,1,icount,1:2]<-c(NA,NA)
        iparent<-iparent+1
       }
       else{
        icount<-icount+1
        iparent<-iparent+1
     }}}
     if(space[s]=='uniform'){
      require(spatstat)
      pts<-rSSI(sqrt((0.4*4)/(pi*ns[n])),ns[n])  ##Simple Sequential Inhibition,first arg is radius, 2nd arg is # of pts
      holder[iperm,s,d,n,1,1:ns[n],1:2]<-c(pts$x,pts$y)#x&y
     }
     if(space[s]=='grid'){
      sns<-sqrt(ns[n])
      if(sns==round(sqrt(ns[n]))){ ##first check that n[ns] is a perferct square
       pts<- c(0,seq(1/(sns-1),(sns-2)*(1/(sns-1)),1/(sns-1)),1)
       holder[iperm,s,d,n,1,1:ns[n],1:2]<-as.vector(as.matrix(expand.grid(pts,pts)))
      }
      else{
        stop('# of individuals does not completely fill grid, n = 2^x where x is an integer')
     }}
     ##assign species ids and initial diameters
     holder[iperm,s,d,n,1,1:ns[n],3]<-sample(x=1:numsp,size=ns[n],prob=P.mat[,rad.num],replace=TRUE) #sp.id
     holder[iperm,s,d,n,1,1:ns[n],4]<-rep(diams[d],ns[n]) #diameter
     for(t in 1:Tim){
      #growth phase
      #check for closeness
      hasdiam<-(1:ns[n])[!is.na(holder[iperm,s,d,n,t,1:ns[n],4])]
      ndiams<-length(hasdiam)
      if(ndiams>0){
       diamsum<-rep(0,ndiams)
       dists<-as.matrix(dist(cbind(holder[iperm,s,d,n,t,hasdiam,1],holder[iperm,s,d,n,t,hasdiam,2]),upper=T,diag=T))
       for(i in 1:ndiams){
        for(j in (1:ndiams)[-i]){
         if( dists[i,j] < (0.5 * (holder[iperm,s,d,n,t,hasdiam[i],4]+ holder[iperm,s,d,n,t,hasdiam[j],4])) & dists[i,j]>0){
          if(overlap){ ##area of overlap is used as a proxy for competition severity
           ##area of overlap for two intersecting circles
           ##http://mathworld.wolfram.com/Circle-CircleIntersection.html equ 14
           r<-0.5 *holder[iperm,s,d,n,t,hasdiam[i],4] ##radius  1
           R<-0.5 *holder[iperm,s,d,n,t,hasdiam[j],4] ##radius 2
           diamsum[i] <- diamsum[i] + r^2*acos((dists[i,j]^2+r^2-R^2)/(2*dists[i,j]*r)) + R^2*acos((dists[i,j]^2+R^2-r^2)/(2*dists[i,j]*R)) - .5* ((-dists[i,j]+r+R)*(dists[i,j]+r-R)*(dists[i,j]-r+R)*(dists[i,j]+r+R))^.5
          }
          else{ ##just use the summed diameters of neighbors
           diamsum[i] <- diamsum[i] + holder[iperm,s,d,n,t,hasdiam[j],4]
       }}}}
       deltadiam <- rep(0,ndiams)
       for(i in 1:ndiams){
        if(diamsum [i] > 0 ){ 
         if(expon){ ##negative exponential prob of mortality (default)
#          deltadiam[i] <- diamch * exp(-diamsum[i]  * mort.buffer[holder[iperm,s,d,n,t,hasdiam[i],3]]) ##simplest and more complex produce similar results
          deltadiam[i] <- diamch * exp(-diamsum[i]  * mort.buffer) ##simplest and more complex produce similar results
#          deltadiam[i] <- diamch * exp(-diamsum[i] * (1-diams[d]) * mort.buffer)
#          deltadiam[i] <- diamch * exp(-diamsum[i] * (1-holder[iperm,s,d,n,t,hasdiam[i],4]) * mort.buffer)
         }
         else{ ##piecewise linear decreased of mortality, piecwise b/c not negative
          if((diamsum[i] * mort.buffer) < diamch){
           deltadiam[i] <- diamch - (diamsum[i] * (1-diams[d]) * mort.buffer)
          }
          else{
           deltadiam[i] <- 0  
        }}}
        else{
         deltadiam[i] <- diamch
        }
       }
       #probabilistic mortality related to growth
       if(t < Tim){
        alive<-runif(ndiams) > 1 - (deltadiam / diamch)
        holder[iperm,s,d,n,t+1,hasdiam[alive],1:3] <- holder[iperm,s,d,n,t,hasdiam[alive],1:3]
        holder[iperm,s,d,n,t+1,hasdiam[alive],4] <- holder[iperm,s,d,n,t,hasdiam[alive],4] + deltadiam[alive]
  }}}}}}
  print(iperm)
 }
 holder
}

hump.plot<-function(holder,avg=TRUE,log='',add=FALSE,ci=FALSE,poly=FALSE,type='l',lwd=2,pch=19,cls=palette(),cls.poly=rep('grey',dim(holder)[2]),xlab='biomass(diam^bexp)',ylab='species richness',bexp=8/3){
 xlims<-c(min(apply(apply(holder[,1,1,1,,,4]^bexp,c(1,2),function(x){sum(x,na.rm=TRUE)}),2,mean)),
          max(apply(apply(holder[,dim(holder)[2],dim(holder)[3],dim(holder)[4],,,4]^bexp,c(1,2),function(x){sum(x,na.rm=TRUE)}),2,mean)))
 ylims<-c(1,numsp)
 if(avg){ ##if avg is true then averages over the permutations are calculated an ploted
  if(add==FALSE){ ##if add is false then new plot is created
   plot(apply(apply(holder[,1,1,1,,,4]^bexp,c(1,2),function(x){sum(x,na.rm=TRUE)}),2,mean),apply(apply(holder[,1,1,1,,,3],c(1,2),function(x){length(unique(x[!is.na(x)]))}),2,mean),type='n',
        xlim=xlims,ylim=ylims,log=log,xlab=xlab,ylab=ylab)
  }
  if(ci){ ##if ci is true then confidence intervals are added as lines
   for(s in 1:dim(holder)[2]){ ##different spatial distributions
    if(is.na(cls[s])) next
    for(d in 1:dim(holder)[3]){ ##different starting diameters
     for(n in 1:dim(holder)[4]){ ##differnt starting individual densities
#      bios<-apply(apply(holder[,s,d,n,,,4]^bexp,c(1,2),function(x){sum(x,na.rm=TRUE)}),2,mean)
      b.low<-apply(apply(holder[,s,d,n,,,4]^bexp,c(1,2),function(x){sum(x,na.rm=TRUE)}),2,function(x){quantile(x,probs=.025)})
      b.high<-apply(apply(holder[,s,d,n,,,4]^bexp,c(1,2),function(x){sum(x,na.rm=TRUE)}),2,function(x){quantile(x,probs=.975)})
      s.low<-apply(apply(holder[,s,d,n,,,3],c(1,2),function(x){length(unique(x[!is.na(x)]))}),2,function(x){quantile(x,probs=.025)})
      s.high<-apply(apply(holder[,s,d,n,,,3],c(1,2),function(x){length(unique(x[!is.na(x)]))}),2,function(x){quantile(x,probs=.975)})
      ##necessary to work with logs...
      b.low<-ifelse(b.low==0,b.low+.00001,b.low)
      s.low<-ifelse(s.low==0,s.low+.00001,s.low)
      if(poly){ 
       polygon(c(b.low,rev(b.high)),c(s.low,rev(s.high)),col=cls.poly[s],border=NA)
      }
      else{
#       points(bios,s.low,type='l',col=cls[s],lty=2)
#       points(bios,s.high,type='l',col=cls[s],lty=2)
       points(b.low,s.low,type='l',col=cls[s],lty=2)
       points(b.high,s.high,type='l',col=cls[s],lty=2)
  }}}}}  
  for(s in 1:dim(holder)[2]){ ##different spatial distributions
   if(is.na(cls[s])) next
   for(d in 1:dim(holder)[3]){ ##different starting diameters
    for(n in 1:dim(holder)[4]){ ##differnt starting individual densities
     points(apply(apply(holder[,s,d,n,,,4]^bexp,c(1,2),function(x){sum(x,na.rm=TRUE)}),2,mean),apply(apply(holder[,s,d,n,,,3],c(1,2),function(x){length(unique(x[!is.na(x)]))}),2,mean),type=type,lwd=lwd,col=cls[s],pch=pch)
  }}}
 }
 else{ #avg is FALSE and no averaging is performed
  if(add==FALSE){ ##if add is false then new plot is created
   plot(apply(holder[,1,1,1,,,4]^bexp,c(1,2),function(x){sum(x,na.rm=TRUE)}),apply(holder[,1,1,1,,,3],c(1,2),function(x){length(unique(x[!is.na(x)]))}),type='n',
      xlim=xlims,ylim=ylims,log=log,xlab=xlab,ylab=ylab)
  }
  for(s in 1:dim(holder)[2]){ ##different spatial distributions
   if(is.na(cls[s])) next
   for(d in 1:dim(holder)[3]){ ##different starting diameters
    for(n in 1:dim(holder)[4]){ ##differnt starting individual densities
     points(apply(holder[,s,d,n,,,4]^bexp,c(1,2),function(x){sum(x,na.rm=TRUE)}),apply(holder[,s,d,n,,,3],c(1,2),function(x){length(unique(x[!is.na(x)]))}),type=type,lwd=lwd,col=cls[s],pch=pch)
  }}}
 }
 abline(h=numsp,col='grey',lty=2)
}

indiv.plot<-function(holder,avg=TRUE,log='',add=FALSE,ci=FALSE,poly=FALSE,type='l',lwd=2,pch=19,cls=palette(),cls.poly=rep('grey',dim(holder)[2]),xlab='biomass(diam^bexp)',ylab='# of individuals',bexp=8/3){
 ##plots biomass vs num of indivs
 xlims<-c(min(apply(apply(holder[,1,1,1,,,4]^bexp,c(1,2),function(x){sum(x,na.rm=TRUE)}),2,mean)),
          max(apply(apply(holder[,dim(holder)[2],dim(holder)[3],dim(holder)[4],,,4]^bexp,c(1,2),function(x){sum(x,na.rm=TRUE)}),2,mean)))
 ylims<-c(1,max(ns))
 if(avg){ ##if avg is true then averages over the permutations are calculated an ploted
  if(add==FALSE){ ##if add is false then new plot is created
   plot(apply(apply(holder[,1,1,1,,,4]^bexp,c(1,2),function(x){sum(x,na.rm=TRUE)}),2,mean),apply(apply(holder[,1,1,1,,,3],c(1,2),function(x){length((x[!is.na(x)]))}),2,mean),type='n',
        xlim=xlims,ylim=ylims,log=log,xlab=xlab,ylab=ylab)
  }
  if(ci){ ##if ci is true then confidence intervals are added as lines
   for(s in 1:dim(holder)[2]){ ##different spatial distributions
    if(is.na(cls[s])) next
    for(d in 1:dim(holder)[3]){ ##different starting diameters
     for(n in 1:dim(holder)[4]){ ##differnt starting individual densities
#      bios<-apply(apply(holder[,s,d,n,,,4]^bexp,c(1,2),function(x){sum(x,na.rm=TRUE)}),2,mean)
      b.low<-apply(apply(holder[,s,d,n,,,4]^bexp,c(1,2),function(x){sum(x,na.rm=TRUE)}),2,function(x){quantile(x,probs=.025)})
      b.high<-apply(apply(holder[,s,d,n,,,4]^bexp,c(1,2),function(x){sum(x,na.rm=TRUE)}),2,function(x){quantile(x,probs=.975)})
      s.low<-apply(apply(holder[,s,d,n,,,3],c(1,2),function(x){length((x[!is.na(x)]))}),2,function(x){quantile(x,probs=.025)})
      s.high<-apply(apply(holder[,s,d,n,,,3],c(1,2),function(x){length((x[!is.na(x)]))}),2,function(x){quantile(x,probs=.975)})
      ##necessary to work with logs...
      b.low<-ifelse(b.low==0,b.low+.00001,b.low)
      s.low<-ifelse(s.low==0,s.low+.00001,s.low)
      if(poly){ 
       polygon(c(b.low,rev(b.high)),c(s.low,rev(s.high)),col=cls.poly[s],border=NA)
      }
      else{
#       points(bios,s.low,type='l',col=cls[s],lty=2)
#       points(bios,s.high,type='l',col=cls[s],lty=2)
       points(b.low,s.low,type='l',col=cls[s],lty=2)
       points(b.high,s.high,type='l',col=cls[s],lty=2)
  }}}}}  
  for(s in 1:dim(holder)[2]){ ##different spatial distributions
   if(is.na(cls[s])) next
   for(d in 1:dim(holder)[3]){ ##different starting diameters
    for(n in 1:dim(holder)[4]){ ##differnt starting individual densities
     points(apply(apply(holder[,s,d,n,,,4]^bexp,c(1,2),function(x){sum(x,na.rm=TRUE)}),2,mean),apply(apply(holder[,s,d,n,,,3],c(1,2),function(x){length((x[!is.na(x)]))}),2,mean),type=type,lwd=lwd,col=cls[s],pch=pch)
  }}}
 }
 else{ #avg is FALSE and no averaging is performed
  if(add==FALSE){ ##if add is false then new plot is created
   plot(apply(holder[,1,1,1,,,4]^bexp,c(1,2),function(x){sum(x,na.rm=TRUE)}),apply(holder[,1,1,1,,,3],c(1,2),function(x){length((x[!is.na(x)]))}),type='n',
      xlim=xlims,ylim=ylims,log=log,xlab=xlab,ylab=ylab)
  }
  for(s in 1:dim(holder)[2]){ ##different spatial distributions
   if(is.na(cls[s])) next
   for(d in 1:dim(holder)[3]){ ##different starting diameters
    for(n in 1:dim(holder)[4]){ ##differnt starting individual densities
     points(apply(holder[,s,d,n,,,4]^bexp,c(1,2),function(x){sum(x,na.rm=TRUE)}),apply(holder[,s,d,n,,,3],c(1,2),function(x){length((x[!is.na(x)]))}),type=type,lwd=lwd,col=cls[s],pch=pch)
  }}}
 }
}

time.plot<-function(holder,avg=TRUE,log='',add=FALSE,ci=FALSE,poly=FALSE,type='l',lwd=2,pch=19,cls=palette(),cls.poly=rep('grey',dim(holder)[2]),xlab='Time',ylab='# of individuals'){
 ##plots time vs num of indivs
 xlims<-c(1,Tim)
 ylims<-c(1,max(ns))
 if(avg){ ##if avg is true then averages over the permutations are calculated an ploted
  if(add==FALSE){ ##if add is false then new plot is created
   plot(1:Tim,apply(apply(holder[,1,1,1,,,3],c(1,2),function(x){length((x[!is.na(x)]))}),2,mean),type='n',
        xlim=xlims,ylim=ylims,log=log,xlab=xlab,ylab=ylab)
  }
  if(ci){ ##if ci is true then confidence intervals are added as lines
   for(s in 1:dim(holder)[2]){ ##different spatial distributions
    if(is.na(cls[s])) next
    for(d in 1:dim(holder)[3]){ ##different starting diameters
     for(n in 1:dim(holder)[4]){ ##differnt starting individual densities
      s.low<-apply(apply(holder[,s,d,n,,,3],c(1,2),function(x){length((x[!is.na(x)]))}),2,function(x){quantile(x,probs=.025)})
      s.high<-apply(apply(holder[,s,d,n,,,3],c(1,2),function(x){length((x[!is.na(x)]))}),2,function(x){quantile(x,probs=.975)})
      ##necessary to work with logs...
      s.low<-ifelse(s.low==0,s.low+.00001,s.low)
      if(poly){ 
       polygon(c(1:Tim,rev(1:Tim)),c(s.low,rev(s.high)),col=cls.poly[s],border=NA)
      }
      else{
       points(1:Tim,s.low,type='l',col=cls[s],lty=2)
       points(1:Tim,s.high,type='l',col=cls[s],lty=2)
  }}}}}  
  for(s in 1:dim(holder)[2]){ ##different spatial distributions
   if(is.na(cls[s])) next
   for(d in 1:dim(holder)[3]){ ##different starting diameters
    for(n in 1:dim(holder)[4]){ ##differnt starting individual densities
     points(1:Tim,apply(apply(holder[,s,d,n,,,3],c(1,2),function(x){length((x[!is.na(x)]))}),2,mean),type=type,lwd=lwd,col=cls[s],pch=pch)
  }}}
 }
 else{ #avg is FALSE and no averaging is performed
  if(add==FALSE){ ##if add is false then new plot is created
   plot(1:Tim,apply(holder[,1,1,1,,,3],c(1,2),function(x){length((x[!is.na(x)]))}),type='n',
      xlim=xlims,ylim=ylims,log=log,xlab=xlab,ylab=ylab)
  }
  for(s in 1:dim(holder)[2]){ ##different spatial distributions
   if(is.na(cls[s])) next
   for(d in 1:dim(holder)[3]){ ##different starting diameters
    for(n in 1:dim(holder)[4]){ ##differnt starting individual densities
     points(1:Tim,apply(holder[,s,d,n,,,3],c(1,2),function(x){length((x[!is.na(x)]))}),type=type,lwd=lwd,col=cls[s],pch=pch)
  }}}
 }
}



thin.plot<-function(holder,avg=TRUE,log='',add=FALSE,ci=FALSE,poly=FALSE,type='l',lwd=2,pch=19,cls=palette(),cls.poly=rep('grey',dim(holder)[2]),xlab='number of individuals',ylab='mean biomass/individual',bexp=8/3){
 ##plots mean biomass per tree vs # of individuals
 ylims<-c(min(apply(apply(holder[,dim(holder)[2],dim(holder)[3],dim(holder)[4],,,4]^bexp,c(1,2),function(x){mean(x/length(x[!is.na(x)]),na.rm=TRUE)}),2,mean),na.rm=T),
  max(apply(apply(holder[,dim(holder)[2],dim(holder)[3],dim(holder)[4],,,4]^bexp,c(1,2),function(x){mean(x/length(x[!is.na(x)]),na.rm=TRUE)}),2,mean),na.rm=T))
 xlims<-c(1,max(ns))
 if(avg){ ##if avg is true then averages over the permutations are calculated an ploted
  if(add==FALSE){ ##if add is false then new plot is created
   plot(apply(apply(holder[,1,1,1,,,3],c(1,2),function(x){length((x[!is.na(x)]))}),2,mean),apply(apply(holder[,1,1,1,,,4]^bexp,c(1,2),function(x){mean(x/length(x[!is.na(x)]),na.rm=TRUE)}),2,mean),type='n',
        xlim=xlims,ylim=ylims,log=log,xlab=xlab,ylab=ylab)
  }
  if(ci){ ##if ci is true then confidence intervals are added as lines
   for(s in 1:dim(holder)[2]){ ##different spatial distributions
    if(is.na(cls[s])) next
    for(d in 1:dim(holder)[3]){ ##different starting diameters
     for(n in 1:dim(holder)[4]){ ##differnt starting individual densities
      s.low<-apply(apply(holder[,s,d,n,,,4]^bexp,c(1,2),function(x){mean(x/length(x[!is.na(x)]),na.rm=TRUE)}),2,function(x){quantile(x,probs=.025)})
      s.high<-apply(apply(holder[,s,d,n,,,4]^bexp,c(1,2),function(x){mean(x/length(x[!is.na(x)]),na.rm=TRUE)}),2,function(x){quantile(x,probs=.975)})
      b.low<-apply(apply(holder[,s,d,n,,,3],c(1,2),function(x){length((x[!is.na(x)]))}),2,function(x){quantile(x,probs=.025)})
      b.high<-apply(apply(holder[,s,d,n,,,3],c(1,2),function(x){length((x[!is.na(x)]))}),2,function(x){quantile(x,probs=.975)})
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
     points(apply(apply(holder[,s,d,n,,,3],c(1,2),function(x){length((x[!is.na(x)]))}),2,mean),apply(apply(holder[,s,d,n,,,4]^bexp,c(1,2),function(x){mean(x/length(x[!is.na(x)]),na.rm=TRUE)}),2,mean),type=type,lwd=lwd,col=cls[s],pch=pch)
  }}}
 }
 else{ #avg is FALSE and no averaging is performed
  if(add==FALSE){ ##if add is false then new plot is created
   plot(apply(holder[,1,1,1,,,3],c(1,2),function(x){length((x[!is.na(x)]))}),apply(holder[,1,1,1,,,4]^bexp,c(1,2),function(x){mean(x/length(x[!is.na(x)]),na.rm=TRUE)}),type='n',
      xlim=xlims,ylim=ylims,log=log,xlab=xlab,ylab=ylab)
  }
  for(s in 1:dim(holder)[2]){ ##different spatial distributions
   if(is.na(cls[s])) next
   for(d in 1:dim(holder)[3]){ ##different starting diameters
    for(n in 1:dim(holder)[4]){ ##differnt starting individual densities
     points(apply(holder[,s,d,n,,,3],c(1,2),function(x){length((x[!is.na(x)]))}),apply(holder[,s,d,n,,,4]^bexp,c(1,2),function(x){mean(x/length(x[!is.na(x)]),na.rm=TRUE)}),type=type,lwd=lwd,col=cls[s],pch=pch)
  }}}
 }
 ##-3/2 thinning lines for comparison
 x<-seq(1,max(ns),1)
 c<-seq(.01,.5,.05)
 for(i in 1:length(c)){
  lines(x,c[i]*x^-(3/2),col='grey')
 }
 ##-4/3 thinning lines for comparison
 x<-seq(1,max(ns),1)
 c<-seq(.01,.5,.05)
 for(i in 1:length(c)){
  lines(x,c[i]*x^-(4/3),col='grey50',lty=2,lwd=2)
 }
 legend('bottomleft',c('-3/2 thinning lines','-4/3 thinning lines'),col=c('grey','grey50'),lty=c(1,2),lwd=c(1,2),bty='n')
}

raref.plot<-function(holder,avg=TRUE,log='',add=FALSE,ci=FALSE,poly=FALSE,type='l',lwd=2,pch=19,cls=palette(),cls.poly=rep('grey',dim(holder)[2]),ylab='species richness',xlab='# of individuals'){
 ##plots num of indivs vs richness
 xlims<-c(0,max(ns))
 ylims<-c(1,numsp)
 if(avg){ ##if avg is true then averages over the permutations are calculated an ploted
  if(add==FALSE){ ##if add is false then new plot is created
   plot(apply(apply(holder[,1,1,1,,,3],c(1,2),function(x){length((x[!is.na(x)]))}),2,mean),apply(apply(holder[,1,1,1,,,3],c(1,2),function(x){length(unique(x[!is.na(x)]))}),2,mean),type='n',
        xlim=xlims,ylim=ylims,log=log,xlab=xlab,ylab=ylab)
  }
  if(ci){ ##if ci is true then confidence intervals are added as lines
   for(s in 1:dim(holder)[2]){ ##different spatial distributions
    if(is.na(cls[s])) next
    for(d in 1:dim(holder)[3]){ ##different starting diameters
     for(n in 1:dim(holder)[4]){ ##differnt starting individual densities
      b.low<-apply(apply(holder[,s,d,n,,,3],c(1,2),function(x){length((x[!is.na(x)]))}),2,function(x){quantile(x,probs=.025)})
      b.high<-apply(apply(holder[,s,d,n,,,3],c(1,2),function(x){length((x[!is.na(x)]))}),2,function(x){quantile(x,probs=.975)})
      s.low<-apply(apply(holder[,s,d,n,,,3],c(1,2),function(x){length(unique(x[!is.na(x)]))}),2,function(x){quantile(x,probs=.025)})
      s.high<-apply(apply(holder[,s,d,n,,,3],c(1,2),function(x){length(unique(x[!is.na(x)]))}),2,function(x){quantile(x,probs=.975)})
      ##necessary to work with logs...
      b.low<-ifelse(b.low==0,b.low+.00001,b.low)
      s.low<-ifelse(s.low==0,s.low+.00001,s.low)
      if(poly){ 
       polygon(c(b.low,rev(b.high)),c(s.low,rev(s.high)),col=cls.poly[s],border=NA)
      }
      else{
#       points(bios,s.low,type='l',col=cls[s],lty=2)
#       points(bios,s.high,type='l',col=cls[s],lty=2)
       points(b.low,s.low,type='l',col=cls[s],lty=2)
       points(b.high,s.high,type='l',col=cls[s],lty=2)
  }}}}}  
  for(s in 1:dim(holder)[2]){ ##different spatial distributions
   if(is.na(cls[s])) next
   for(d in 1:dim(holder)[3]){ ##different starting diameters
    for(n in 1:dim(holder)[4]){ ##differnt starting individual densities
     points(apply(apply(holder[,s,d,n,,,3],c(1,2),function(x){length((x[!is.na(x)]))}),2,mean),apply(apply(holder[,s,d,n,,,3],c(1,2),function(x){length(unique(x[!is.na(x)]))}),2,mean),type=type,lwd=lwd,col=cls[s],pch=pch)
  }}}
 }
 else{ #avg is FALSE and no averaging is performed
  if(add==FALSE){ ##if add is false then new plot is created
   plot(apply(holder[,1,1,1,,,3],c(1,2),function(x){length((x[!is.na(x)]))}),apply(holder[,1,1,1,,,3],c(1,2),function(x){length(unique(x[!is.na(x)]))}),type='n',
      xlim=xlims,ylim=ylims,log=log,xlab=xlab,ylab=ylab)
  }
  for(s in 1:dim(holder)[2]){ ##different spatial distributions
   if(is.na(cls[s])) next
   for(d in 1:dim(holder)[3]){ ##different starting diameters
    for(n in 1:dim(holder)[4]){ ##differnt starting individual densities
     points(apply(holder[,s,d,n,,,3],c(1,2),function(x){length((x[!is.na(x)]))}),apply(holder[,s,d,n,,,3],c(1,2),function(x){length(unique(x[!is.na(x)]))}),type=type,lwd=lwd,col=cls[s],pch=pch)
  }}}
 }
 abline(h=numsp,col='grey',lty=2)
}

siml<-function(holder,abu=TRUE,ci=FALSE){
 ##calculates average similarity across permuations for a given point in time
 require(vegan) 
 dists<-array(NA,dim=dim(holder)[2:5])
 ci.int<-array(NA,dim=c(dim(holder)[2:5],2))
 for(s in 1:dim(holder)[2]){
  for(d in 1:dim(holder)[3]){
   for(n in 1:dim(holder)[4]){
    spxsite<-array(0,dim=c(dim(holder)[1],numsp,Tim))
    for(iperm in 1:dim(holder)[1]){
     for(t in 1:Tim){
      for(indiv in 1:ns[n]){
       if(!is.na(holder[iperm,s,d,n,t,indiv,4])){
        spxsite[iperm,holder[iperm,s,d,n,t,indiv,3],t] <- spxsite[iperm,holder[iperm,s,d,n,t,indiv,3],t] + holder[iperm,s,d,n,t,indiv,4]^(8/3)
    }}}}
    if(ci){
     if(abu){
      ci.int[s,d,n,,1]<-apply(spxsite,3,function(x){quantile(vegdist(x,method='horn'),probs=.025,na.rm = TRUE)})
      ci.int[s,d,n,,2]<-apply(spxsite,3,function(x){quantile(vegdist(x,method='horn'),probs=.975,na.rm = TRUE)})
     }
     else{
      ci.int[s,d,n,,1]<-apply(spxsite,3,function(x){quantile(vegdist(ifelse(x>0,1,0),method='jaccard'),probs=.025,na.rm = TRUE)})
      ci.int[s,d,n,,2]<-apply(spxsite,3,function(x){quantile(vegdist(ifelse(x>0,1,0),method='jaccard'),probs=.975,na.rm = TRUE)})
    }}
    else{
     if(abu)
      dists[s,d,n,]<-apply(spxsite,3,function(x){mean(vegdist(x,method='horn'))})
     else
      dists[s,d,n,]<-apply(spxsite,3,function(x){mean(vegdist(ifelse(x>0,1,0),method='jaccard'))})
 }}}}
 if(ci)
  ci.int
 else
  dists
}



#quantile(apply(holder[,s,d,n,1,,3],2,function(x){length(unique(x[!is.na(x)]))}),probs=c(.025,.957))

##non-parametric ci
#hist(apply(holder[,s,d,n,1,,3],2,function(x){length(unique(x[!is.na(x)]))}))
#points(quantile(apply(holder[,s,d,n,1,,3],2,function(x){length(unique(x[!is.na(x)]))}),probs=c(.025,.957)),c(0,0),col='red')
#sum(apply(holder[,s,d,n,1,,3],2,function(x){length(unique(x[!is.na(x)]))})>58)
#sum(apply(holder[,s,d,n,1,,3],2,function(x){length(unique(x[!is.na(x)]))})<=69)
#      

#     Bmeans<-apply(apply(holder[,s,d,n,,,4]^bexp,c(1,2),function(x){sum(x,na.rm=TRUE)}),2,mean)
#     Bsds<-apply(apply(holder[,s,d,n,,,4]^bexp,c(1,2),function(x){sum(x,na.rm=TRUE)}),2,sd)
#     Bzval<-ifelse(Bsds>0,(Bmeans-Bsds)/Bsds,0)
#     SRsds<-apply(apply(holder[,s,d,n,,,3],c(1,2),function(x){length(unique(x[!is.na(x)]))}),2,sd)
#     SRmeans<-apply(apply(holder[,s,d,n,,,3],c(1,2),function(x){length(unique(x[!is.na(x)]))}),2,mean)
#     SRzval<-(SRmeans-SRsds)/SRsds
