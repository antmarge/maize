#outlier removal based on SAS mixed model REML

getOutliers<-function(obs){
  
  library(lme4)
  colnames(obs)<-c("year","plot","tasselNm","tasselWt","tasselBr","bn2","secBN","L1","L2","BLtop","BLmid","BLlow","range","row","item","accession","notes")
  editedCSV="edited_raw1.csv"

  traits<-c("tasselWt","tasselBr","bn2","secBN","L1","L2","BLtop","BLmid","BLlow")
  
  # calculate threshold
  N=nrow(obs)
  #threshold=2
  threshold=qt((1-.1/(2*N)),N-5)
  print("Threshold")
  print(threshold)
  
  #To store indices of rows that will be altered
  allInd<-vector(mode="numeric", length=0)
  
  print ("For each trait, looking for outliers based on studentized residuals that do no satisfy the threshold")
  for (trait in traits){
    print(trait)
    before<-obs[,trait]
    model<-lmer(obs[,trait] ~ (1|item)+(1|year)+(1|year:plot)+(1|year:row)+(1|year:range),data=obs,REML=TRUE,na.action=na.exclude)
    #plot(model)
    stdev=sd(residuals(model),na.rm=TRUE)
    stud.resid<-intstudres(residuals(model))
    
    #Append internally studentized residuals in column of main data
    rescol=paste0(trait,".studres")
    obs[,rescol]<-stud.resid
    
    #How many outliers were identified for this trait?
    remInd<-which(abs(stud.resid)>threshold)
    print(length(remInd)) 
    #Change the trait values of those outliers to NA
    obs[[trait]][abs(obs[[rescol]])>threshold]<-NA
    
    #Box plot to show outlier removal
    ptitle<-paste0("Effect of outlier removal for ",trait)
    labs=c("Before","After")
    after<-obs[,trait]
    boxplot(before,after,use.columns=TRUE,col="grey",outcol="darkred",main=ptitle,names=labs)
    
    s<-summary(model)
    capture.output(s, file = "modelSummary.txt",append=TRUE)
    allInd<-c(allInd,remInd)
  }
  length(allInd)
  print("Removing outliers from data\n")
  #Now remove all unwanted observations (that didn't satisfy the threshold)
  outliers<-obs[allInd,]
  numOut=nrow(outliers)
  print("Outliers removed: ")
  print(numOut)
  keep<-obs[-allInd,]

  #Do boxplot
  total<-obs[,4:12]
  labs<-c("weight","nPB","nPBwXB","nSB","len1","len2(MR)","BLtop","BLmid","BLlow")
  boxplot(total,use.columns=TRUE,col="grey",outcol="darkred",las=2,names=labs,main="After outlier removal based on residuals")


  write.csv(obs,"outlierRem_raw2.csv")
  print ("Check output files for outliers and filtered.")

  return(obs)
}

intstudres<-function(residuals){
  #Internally studentized residuals quantify how large the residuals are 
  #in standard deviation units, and therefore can be easily used to identify outliers.
  
  #Store all new internal standard dev 
  allsd<-vector()
  final<-data.frame(matrix())
  len=length(residuals)
  
  for (i in 1:len){
    res=residuals[i]
    #remove that residual to calcualte sd on temp vector
    temp<-residuals[-i]
    indivsd<-sd(temp,na.rm=TRUE)
    allsd[i]<-indivsd
  }
  
  #Now calculate internal studentized residual
  stud.resid<-(residuals/allsd)
  return(stud.resid)
  
}
