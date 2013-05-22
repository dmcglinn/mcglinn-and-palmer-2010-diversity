########################################################
# FILENAME: "hump.sim.R"
# AUTHOR: Dan McGlinn
# EMAIL: danmcglinn@gmail.com
# DATE CREATED: 11.15.2009
# LAST MODIFIED: 3.10.2010
# PURPOSE: To provide the code for generating the relative abundance distributions (RADs), 
# simulation model and secondary graphical programs needed to reproduce the results 
# discussed in the manuscript
########################################################
########################################################
##create the nine relative abundance distributions (RADs)##
S<-numsp #specify size of the species pool
##even RAD##
p1<-rep(1/S,S)
#############
##now generate the lognormals as averages of random lognormal distributions##
reps<-1000 ##number of replictates to average over
vars<-c(1,2,4) ##variances to use for the 3 LOGN RADs
X<-array(NA,dim=c(S,reps,length(vars))) #create holding array
for(k in 1:length(vars)){ #loop through # of variances
 for(i in 1:1000){ #loop through replicates
  X[,i,k]<-sort(rlnorm(S,0,vars[k]),decreasing=TRUE) #populate array with sorted values from largest to smallest
 }
}
X.m<-array(NA,dim=c(S,length(vars))) #create holding array
for(k in 1:length(vars)) 
 X.m[,k]<-sort(apply(X[,,k],1,mean),decreasing=TRUE)##calcualte mean and then order each distribution from largest to smallest
P.m<-array(NA,dim=c(S,length(vars))) #holding arrary for RAD
for(k in 1:length(vars))
 P.m[,k]<-X.m[,k]/sum(X.m[,k]) #now convert the abundance distribution into relative abundance distributions
##########################
##create the uneven RAD##
prop<-.99  #proportion of species pool dominated by one species
p5<-rep(0,S) #create holding vector
p5[1]<-prop #domiant species set
p5[-1]<-rep((1-prop)/(S-1),S-1) #the remaining S-1 species set
########################
##geometric series##
geom <- function(a1, k, S){
 ai <- a1*k^((1:S)-1)
 pi <- ai/sum(ai)
 pi
}
p6<-geom(10,.9,S)
##################
##broken stick##
brok <- function(abar, S){
 ai <- rep(NA,S)
 for (i in 1:S)
  ai[i] <- abar*sum(1/i:S)
 pi <- ai/sum(ai)
 pi
}
p7<-brok(10,S)
##################
##Zipf##
zipf <- function(a1,gamma,S){
 ai <- a1*(1:S)^-gamma
 pi <- ai / sum(ai)
 pi
}
p8<-zipf(10,1.3,S)
###################
##Zipf-Mandelbrot##
zipf.man <- function(a1,gamma,beta,S){
 ai <- a1*((1:S)+beta)^-gamma
 pi <- ai / sum(ai)
 pi
}
p9<-zipf.man(10,1.3,100,S)
#########################
P.mat<-cbind(p1,P.m[,1],p6,P.m[,3],p5) #bind some of the RAD vectors together into one matrix
colnames(P.mat)<-c('even','lnorm(0,1)','geometric','lnorm(0,4)','uneven')
##note that 'P.mat' only contains the RADs mentioned in the manuscript

########################################################
########the function that simulates the results#########
########################################################
hump.sim<-function(nperm,diams,ns,RAD='even',space='random',intensity=.8,all.space=FALSE,expon=TRUE,overlap=FALSE){
 ###the simulation model###
 ##FUNCTION ARGUMENTS:
 ##'nperm' is number of permutations a single set of starting parameters will be used
 ##'diams' is the input starting diameters
 ##'ns' is the input starting individual numbers
 ##'space' specifies spatial distribution: 'clustered','random','uniform'
 ##'intensity' specifies how strongly spatial clustering will occur (0-1)
 ##'all.space' if TRUE then all three spatial distributions will be simulated
 ##'expon' if TRUE then the probability of mortality will follow a negative exponential function, linear methods are available
 ##'overlap' if TRUE then the degree of competition is based on the area of overlap with neighbors rather than the the summed diameter or neighbors
 ##Note: specification of the allometric exponent is in the function 'hump.plot' not this function
 #####################
 ##Check that the RAD is parameterized correctly
 RADs<-c('even','lnorm(0,1)','geometric','lnorm(0,4)','uneven')
 rad.num<-match(RAD,RADs)
 if(is.na(rad.num)){
  print(paste(c('Choose one of the appropriate RADs:',RADs)))
  stop()
 }
 ##Check that the spatial distributions are parameterized correctly
 spaces<-c('clustered','random','uniform')
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
 #####################
 ##create holder object that will keep the output from the simulation
 holder<-array(NA,dim=c(nperm,len.space,len.diam,len.n,Tim,max(ns),4))
 ##loop through all parameter combinations
 for(iperm in 1:nperm){ ##permuations
  for(s in 1:len.space){ ##spatial distributions
   for(d in 1:len.diam){ ##dimeter sizes
    for(n in 1:len.n){ ##initial densities
     #######################
     ##1-ESTABLISHMENT PHASE
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
      sns<-sqrt(ns[n])
      if(sns==round(sqrt(ns[n]))){ ##first check that n[ns] is a perferct square
       pts<- c(0,seq(1/(sns-1),(sns-2)*(1/(sns-1)),1/(sns-1)),1)
       holder[iperm,s,d,n,1,1:ns[n],1:2]<-as.vector(as.matrix(expand.grid(pts,pts)))
      }
      else{
        stop('# of individuals does not completely fill grid, n = 2^x where x is an integer')
     }}
     ##assign species ids and initial diameters
     ##note that the object 'P.mat' which contains each species probability of being drawn from the 
     ##the species pool is defined outside of this function ('hump.sim').
     holder[iperm,s,d,n,1,1:ns[n],3]<-sample(x=1:numsp,size=ns[n],prob=P.mat[,rad.num],replace=TRUE) #sp.id
     holder[iperm,s,d,n,1,1:ns[n],4]<-rep(diams[d],ns[n]) #starting diameter
     for(t in 1:Tim){ ##loop through time steps
      ################
      ##2-GROWTH PHASE
      hasdiam<-(1:ns[n])[!is.na(holder[iperm,s,d,n,t,1:ns[n],4])]
      ndiams<-length(hasdiam)
      if(ndiams>0){
       diamsum<-rep(0,ndiams)
       ##check for closeness
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
         if(expon){ ##negative exponential prob of mortality (default), Equ.1
          deltadiam[i] <- diamch * exp(-diamsum[i]) ##simplest and more complex produce similar results
         }
         else{ ##piecewise linear decreased of mortality, piecwise b/c prob of mort cannot be a negative value
          if(diamsum[i] < diamch){
           deltadiam[i] <- diamch - (diamsum[i] * (1-diams[d]))
          }
          else{
           deltadiam[i] <- 0  
        }}}
        else{
         deltadiam[i] <- diamch
        }
       }
       ##################
       ##3-THINNING PHASE
       if(t < Tim){
        alive<-runif(ndiams) > 1 - (deltadiam / diamch) ##Equ.2
        holder[iperm,s,d,n,t+1,hasdiam[alive],1:3] <- holder[iperm,s,d,n,t,hasdiam[alive],1:3]
        holder[iperm,s,d,n,t+1,hasdiam[alive],4] <- holder[iperm,s,d,n,t,hasdiam[alive],4] + deltadiam[alive]
  }}}}}}
  print(iperm) ##let user know what permutation the simulation is on
 }
 holder ##output simulation result
}
########################################################
########################################################


##graphical script for the biomass-richness relationship##
hump.plot<-function(holder,avg=TRUE,log='',add=FALSE,ci=FALSE,poly=FALSE,type='l',lwd=2,pch=19,cls=palette(),cls.poly=rep('grey',dim(holder)[2]),main='',xlab='biomass',ylab='species richness',allo=8/3){
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
      b.low<-apply(apply(holder[,s,d,n,-1,,4]^allo,c(1,2),function(x){sum(x,na.rm=TRUE)}),2,function(x){quantile(x,probs=.025)})
      b.high<-apply(apply(holder[,s,d,n,-1,,4]^allo,c(1,2),function(x){sum(x,na.rm=TRUE)}),2,function(x){quantile(x,probs=.975)})
      s.low<-apply(apply(holder[,s,d,n,-1,,3],c(1,2),function(x){length(unique(x[!is.na(x)]))}),2,function(x){quantile(x,probs=.025)})
      s.high<-apply(apply(holder[,s,d,n,-1,,3],c(1,2),function(x){length(unique(x[!is.na(x)]))}),2,function(x){quantile(x,probs=.975)})
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
     points(apply(apply(holder[,s,d,n,-1,,4]^allo,c(1,2),function(x){sum(x,na.rm=TRUE)}),2,mean),apply(apply(holder[,s,d,n,-1,,3],c(1,2),function(x){length(unique(x[!is.na(x)]))}),2,mean),type=type,lwd=lwd,col=cls[s],pch=pch)
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
     points(apply(holder[,s,d,n,-1,,4]^allo,c(1,2),function(x){sum(x,na.rm=TRUE)}),apply(holder[,s,d,n,-1,,3],c(1,2),function(x){length(unique(x[!is.na(x)]))}),type=type,lwd=lwd,col=cls[s],pch=pch)
  }}}
 }
 abline(h=numsp,col='grey',lty=2)
}
##########################################################

indiv.plot<-function(holder,avg=TRUE,log='',add=FALSE,ci=FALSE,poly=FALSE,type='l',lwd=2,pch=19,cls=palette(),cls.poly=rep('grey',dim(holder)[2]),main='',xlab='biomass',ylab='# of individuals',allo=8/3){
 ##plots num of indivs vs. total biomass
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
 ylims<-c(1,max(ns))
 if(avg){ ##if avg is true then averages over the permutations are calculated and ploted
  if(add==FALSE){ ##if add is false then new plot is created
   plot(apply(apply(holder[,1,1,1,,,4]^allo,c(1,2),function(x){sum(x,na.rm=TRUE)}),2,mean),apply(apply(holder[,1,1,1,,,3],c(1,2),function(x){length((x[!is.na(x)]))}),2,mean),type='n',
        xlim=xlims,ylim=ylims,log=log,xlab=xlab,ylab=ylab,main=main)
  }
  if(ci){ ##if ci is true then confidence intervals are added as lines
   for(s in 1:dim(holder)[2]){ ##different spatial distributions
    if(is.na(cls[s])) next
    for(d in 1:dim(holder)[3]){ ##different starting diameters
     for(n in 1:dim(holder)[4]){ ##differnt starting individual densities
      b.low<-apply(apply(holder[,s,d,n,-1,,4]^allo,c(1,2),function(x){sum(x,na.rm=TRUE)}),2,function(x){quantile(x,probs=.025)})
      b.high<-apply(apply(holder[,s,d,n,-1,,4]^allo,c(1,2),function(x){sum(x,na.rm=TRUE)}),2,function(x){quantile(x,probs=.975)})
      s.low<-apply(apply(holder[,s,d,n,-1,,3],c(1,2),function(x){length((x[!is.na(x)]))}),2,function(x){quantile(x,probs=.025)})
      s.high<-apply(apply(holder[,s,d,n,-1,,3],c(1,2),function(x){length((x[!is.na(x)]))}),2,function(x){quantile(x,probs=.975)})
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
     points(apply(apply(holder[,s,d,n,-1,,4]^allo,c(1,2),function(x){sum(x,na.rm=TRUE)}),2,mean),apply(apply(holder[,s,d,n,,,3],c(1,2),function(x){length((x[!is.na(x)]))}),2,mean),type=type,lwd=lwd,col=cls[s],pch=pch)
  }}}
 }
 else{ #avg is FALSE and no averaging is performed
  if(add==FALSE){ ##if add is false then new plot is created
   plot(apply(holder[,1,1,1,,,4]^allo,c(1,2),function(x){sum(x,na.rm=TRUE)}),apply(holder[,1,1,1,-1,,3],c(1,2),function(x){length((x[!is.na(x)]))}),type='n',
      xlim=xlims,ylim=ylims,log=log,xlab=xlab,ylab=ylab,main=main)
  }
  for(s in 1:dim(holder)[2]){ ##different spatial distributions
   if(is.na(cls[s])) next
   for(d in 1:dim(holder)[3]){ ##different starting diameters
    for(n in 1:dim(holder)[4]){ ##differnt starting individual densities
     points(apply(holder[,s,d,n,-1,,4]^allo,c(1,2),function(x){sum(x,na.rm=TRUE)}),apply(holder[,s,d,n,-1,,3],c(1,2),function(x){length((x[!is.na(x)]))}),type=type,lwd=lwd,col=cls[s],pch=pch)
  }}}
 }
}

time.plot<-function(holder,avg=TRUE,log='',add=FALSE,ci=FALSE,poly=FALSE,type='l',lwd=2,pch=19,cls=palette(),cls.poly=rep('grey',dim(holder)[2]),main='',xlab='Time',ylab='# of individuals'){
 ##plots time vs. num of indivs
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
 #####################
 xlims<-c(1,Tim)
 ylims<-c(1,max(ns))
 if(avg){ ##if avg is true then averages over the permutations are calculated an ploted
  if(add==FALSE){ ##if add is false then new plot is created
   plot(1:Tim,apply(apply(holder[,1,1,1,,,3],c(1,2),function(x){length((x[!is.na(x)]))}),2,mean),type='n',
        xlim=xlims,ylim=ylims,log=log,xlab=xlab,ylab=ylab,main=main)
  }
  if(ci){ ##if ci is true then confidence intervals are added as lines
   for(s in 1:dim(holder)[2]){ ##different spatial distributions
    if(is.na(cls[s])) next
    for(d in 1:dim(holder)[3]){ ##different starting diameters
     for(n in 1:dim(holder)[4]){ ##differnt starting individual densities
      s.low<-apply(apply(holder[,s,d,n,-1,,3],c(1,2),function(x){length((x[!is.na(x)]))}),2,function(x){quantile(x,probs=.025)})
      s.high<-apply(apply(holder[,s,d,n,-1,,3],c(1,2),function(x){length((x[!is.na(x)]))}),2,function(x){quantile(x,probs=.975)})
      ##necessary to work with logs...
      s.low<-ifelse(s.low==0,s.low+.00001,s.low)
      if(poly){ 
       polygon(c(2:Tim,rev(2:Tim)),c(s.low,rev(s.high)),col=cls.poly[s],border=NA)
      }
      else{
       points(2:Tim,s.low,type='l',col=cls[s],lty=2)
       points(2:Tim,s.high,type='l',col=cls[s],lty=2)
  }}}}}  
  for(s in 1:dim(holder)[2]){ ##different spatial distributions
   if(is.na(cls[s])) next
   for(d in 1:dim(holder)[3]){ ##different starting diameters
    for(n in 1:dim(holder)[4]){ ##differnt starting individual densities
     points(2:Tim,apply(apply(holder[,s,d,n,-1,,3],c(1,2),function(x){length((x[!is.na(x)]))}),2,mean),type=type,lwd=lwd,col=cls[s],pch=pch)
  }}}
 }
 else{ #avg is FALSE and no averaging is performed
  if(add==FALSE){ ##if add is false then new plot is created
   plot(1:Tim,apply(holder[,1,1,1,,,3],c(1,2),function(x){length((x[!is.na(x)]))}),type='n',
      xlim=xlims,ylim=ylims,log=log,xlab=xlab,ylab=ylab,main=main)
  }
  for(s in 1:dim(holder)[2]){ ##different spatial distributions
   if(is.na(cls[s])) next
   for(d in 1:dim(holder)[3]){ ##different starting diameters
    for(n in 1:dim(holder)[4]){ ##differnt starting individual densities
     points(2:Tim,apply(holder[,s,d,n,-1,,3],c(1,2),function(x){length((x[!is.na(x)]))}),type=type,lwd=lwd,col=cls[s],pch=pch)
  }}}
 }
}

raref.plot<-function(holder,avg=TRUE,log='',add=FALSE,ci=FALSE,poly=FALSE,type='l',lwd=2,pch=19,cls=palette(),cls.poly=rep('grey',dim(holder)[2]),main='',xlab='# of individuals',ylab='species richness'){
 ##plots num of indivs vs richness
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
 #####################
 xlims<-c(0,max(ns))
 ylims<-c(1,numsp)
 if(avg){ ##if avg is true then averages over the permutations are calculated an ploted
  if(add==FALSE){ ##if add is false then new plot is created
   plot(apply(apply(holder[,1,1,1,,,3],c(1,2),function(x){length((x[!is.na(x)]))}),2,mean),apply(apply(holder[,1,1,1,,,3],c(1,2),function(x){length(unique(x[!is.na(x)]))}),2,mean),type='n',
        xlim=xlims,ylim=ylims,log=log,xlab=xlab,ylab=ylab,main=main)
  }
  if(ci){ ##if ci is true then confidence intervals are added as lines
   for(s in 1:dim(holder)[2]){ ##different spatial distributions
    if(is.na(cls[s])) next
    for(d in 1:dim(holder)[3]){ ##different starting diameters
     for(n in 1:dim(holder)[4]){ ##differnt starting individual densities
      b.low<-apply(apply(holder[,s,d,n,-1,,3],c(1,2),function(x){length((x[!is.na(x)]))}),2,function(x){quantile(x,probs=.025)})
      b.high<-apply(apply(holder[,s,d,n,-1,,3],c(1,2),function(x){length((x[!is.na(x)]))}),2,function(x){quantile(x,probs=.975)})
      s.low<-apply(apply(holder[,s,d,n,-1,,3],c(1,2),function(x){length(unique(x[!is.na(x)]))}),2,function(x){quantile(x,probs=.025)})
      s.high<-apply(apply(holder[,s,d,n,-1,,3],c(1,2),function(x){length(unique(x[!is.na(x)]))}),2,function(x){quantile(x,probs=.975)})
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
     points(apply(apply(holder[,s,d,n,-1,,3],c(1,2),function(x){length((x[!is.na(x)]))}),2,mean),apply(apply(holder[,s,d,n,-1,,3],c(1,2),function(x){length(unique(x[!is.na(x)]))}),2,mean),type=type,lwd=lwd,col=cls[s],pch=pch)
  }}}
 }
 else{ #avg is FALSE and no averaging is performed
  if(add==FALSE){ ##if add is false then new plot is created
   plot(apply(holder[,1,1,1,,,3],c(1,2),function(x){length((x[!is.na(x)]))}),apply(holder[,1,1,1,,,3],c(1,2),function(x){length(unique(x[!is.na(x)]))}),type='n',
      xlim=xlims,ylim=ylims,log=log,xlab=xlab,ylab=ylab,main=main)
  }
  for(s in 1:dim(holder)[2]){ ##different spatial distributions
   if(is.na(cls[s])) next
   for(d in 1:dim(holder)[3]){ ##different starting diameters
    for(n in 1:dim(holder)[4]){ ##differnt starting individual densities
     points(apply(holder[,s,d,n,-1,,3],c(1,2),function(x){length((x[!is.na(x)]))}),apply(holder[,s,d,n,-1,,3],c(1,2),function(x){length(unique(x[!is.na(x)]))}),type=type,lwd=lwd,col=cls[s],pch=pch)
  }}}
 }
 abline(h=numsp,col='grey',lty=2)
}
