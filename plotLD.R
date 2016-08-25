#https://fabiomarroni.wordpress.com/2011/08/09/estimate-decay-of-linkage-disequilibrium-with-distance/

chromLD<-function(chr){
  
  #subset table to only include one chromosome
  c<-ldtable[ldtable$CHR_A==chr,]
  distance<-c$distance
  LD.data<-c$R2
  
  n<-10 #sample size
  LD.st<-c(b0=12.9)
  distance.mb<-distance/1000000
  LD.nonlinear<-nls(LD.data~(1-distance.mb)^b0,start=LD.st,control=nls.control(minFactor=1/1000000000,maxiter=100,warnOnly=T))
  summ<-summary(LD.nonlinear)
  param<-summ$parameters
  beta0<-param["b0","Estimate"]
  fpoints<-(1-distance.mb)^beta0
  
  #Estimate halfpoint of decay
  df<-data.frame(distance,fpoints)
  maxld<-max(LD.data)
  #You could elucubrate if it's better to use the maximum ESTIMATED value of LD
  #In that case you just set: maxld<-max(fpoints)
  h.decay<-maxld/2
  half.decay.distance<-df$distance[which.min(abs(df$fpoints-h.decay))]
  number=length(LD.data)
  dec=paste0("Half decay distance for chromosome #",chr,": ",half.decay.distance," #SNPs: ",number)
  print(dec)
  
  #PLOT
  ld.df<-data.frame(distance,fpoints)
  ld.df<-ld.df[order(ld.df$distance),]
  maintitle=paste0("Linkage disequilibrium decay for chromosome ",chr,"\nHalf decay point at ",half.decay.distance, 
                    " bp distance\n16.08.16")
  plot(distance,LD.data,pch=19,cex=0.9,xlim=c(0,10000),ylab="Linkage disequilibrium (R^2)",xlab="Distance between SNPs (bp)",main=maintitle,col="darkgreen",cex.lab=.8,cex.main=.7)
  lines(ld.df$distance,ld.df$fpoints,lty=3,lwd=1.2)
  abline(v=half.decay.distance,col="red")
  
}

allLD<-function(maxchr){
  
  colormatrix<-c("red","blue","aquamarine4","forestgreen","purple","darkgoldenrod3","steelblue1","brown2","darkmagenta","darkorange")
  distance<-ldtable$distance
  LD.data<-ldtable$R2
  maintitle=paste0("Linkage disequilibrium decay\n16.08.16")
  plot(distance,LD.data,pch=20,xlim=c(0,5000),cex=0.9,ylab="Linkage disequilibrium (R^2)",xlab="Distance between SNPs (bp)",main=maintitle,col="darkgray",cex.lab=.8,cex.main=.7)
  legend(4000,.5,1:10,fill=colormatrix,title="CHR",cex=.7,box.col="white") # gives the legend lines the correct color and width
  
  #sub=c(1,2,5,6,7,8,9,10)
  for (chr in 1:10){
    
    #subset table to only include one chromosome
    colnow=colormatrix[chr]
    c<-ldtable[ldtable$CHR_A==chr,]
    distance<-c$distance
    LD.data<-c$R2
    
    n<-10 #sample size
    LD.st<-c(b0=12.9)
    distance.mb<-distance/1000000
    LD.nonlinear<-nls(LD.data~(1-distance.mb)^b0,start=LD.st,control=nls.control(minFactor=1/1000000000,maxiter=100,warnOnly=T))
    summ<-summary(LD.nonlinear)
    param<-summ$parameters
    beta0<-param["b0","Estimate"]
    fpoints<-(1-distance.mb)^beta0
    
    #Estimate halfpoint of decay
    df<-data.frame(distance,fpoints)
    maxld<-max(LD.data)
    #You could elucubrate if it's better to use the maximum ESTIMATED value of LD
    #In that case you just set: maxld<-max(fpoints)
    h.decay<-maxld/2
    half.decay.distance<-df$distance[which.min(abs(df$fpoints-h.decay))]
    number=length(LD.data)
    dec=paste0("Half decay distance for chromosome #",chr,": ",half.decay.distance," #SNPs: ",number)
    #dec=paste0("Half decay distance for chromosome #",chr,": ",half.decay.distance,"    ",colnow)
    print(dec)
    
    #PLOT
    ld.df<-data.frame(distance,fpoints)
    ld.df<-ld.df[order(ld.df$distance),]
    lines(ld.df$distance,ld.df$fpoints,lty=3,lwd=1.2,col=colnow)
    abline(v=half.decay.distance,col=colnow,lty=1)

  }
}

broken<-function(){
  
    df1 <-  transform(dr, group=cut(as.numeric(sub('[%]', '', distance)), breaks=dr.breaks,labels=dr.labels))
    res <- do.call(data.frame,aggregate(R2~group, df1,FUN=function(x) c(Count=length(x), Sum=sum(x))))
    dNew <- data.frame(group=levels(df1$group))
    m<-merge(res, dNew, all=TRUE)
    m$R2.Avg<-m$R2.Sum/m$R2.Count
    plot(m$group,m$R2.Avg)
    
    dr.breaks
    dr.breaks=c(0,100,200,300,400,600,1000,1500,2000,5000,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000)
    dr.labels=c('0-100','100-200','200-300','300-400','400-600','600-1000','1000-1500','1500-2000',5000,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000)

}