library(lme4)

traitBlups<-function(){
    
    traits<-colnames(orig[,4:12])
    for (t in traits){
      print("********************************")
      print(t)
      getBlups(orig,t)
    }
}

getBlups<-function(d,trait){

    #remove missing values
    data<-d[-which(is.na(d[,trait])),]
    
    #Boxplot
    #boxplot(data[[trait]] ~ data$accession, xlab="Genotype (Accession#)",ylab=trait,main=paste0(trait," by Genotype"),col="red",las=2,xaxt='n')
    
    # Rename variables for ease of use
    data[[trait]] = as.numeric(data[[trait]])
    data$accession = as.factor(data$accession)
    data$year = as.factor(data$year)
    
    ## Calculate variance 
    
    # Load required lme4 package
    library(lme4)
    
    # Linear Model with random effects for variance components
    # <- lmer(tasselWt ~ accession + (1|year), data)
    model <- lmer(data[[trait]] ~ (1|accession) + (1|year),data,na.action=na.omit)
    
    #Extract variance components
      s<-summary(model)
      print(s)
      capture.output(s,file=paste0(trait,"_model.txt"))
    
    #Estimate BLUPS
      allblup<- ranef(mod2,na.action=na.omit)
    
    #Look at output structure
      #str(traitblup)
    
    #Extract blup for genotype (accession #/id)
    traitblup = allblup$accession

    #See the structure of the blup for each line
      #str(traitblup)
    
    #Create a numeric vector with the BLUP for each line
    blupnum = traitblup[,1]
    
    #Save the traitblup output to a separate .csv file
      col2<-paste0("intercept",".",trait)
      outblup<-cbind(rownames(traitblup),blupnum)
      colnames(outblup)<-c("genotype",col2)
      filename<-paste0(trait,"_BLUP.csv")
      write.csv(outblup, file=filename,row.names = FALSE,quote=FALSE)
    
    ## Creating plots with the BLUPs
  
    # Create a histogram with the BLUP for each line
    #hist(blupnum, col="darkgreen",main=paste0("BLUP for each genotype for ",trait))

    ## Compare BLUP to genotype trait averages on a scatterplot
    
    #Get average of trait for each genotype
    acc<-data$accession[complete.cases(data[[trait]])]
    genoAgg<-aggregate(data[[trait]],list(data$accession),mean)
    traitxlab=paste0("Genotype BLUPs for ",trait)
    traitylab=paste0("Average ",trait," for genotype")
    traitmain=paste0("Average trait value vs BLUPs per genotype for ",trait)
    #plot(blupnum, genoAgg$x, col="darkblue",xlab=traitxlab,ylab=traitylab,main=traitmain)
    
}

badblup <- function(data,trait) {
  lmeout1 <- lme(tasselWt ~ (1|year/range) + (1|year/plot), data)
  ped.hat1 <- lmeout1$coef$fixed
  ped.hat1[-1] <- ped.hat1[-1] + ped.hat1[1]
  names(ped.hat1)[1] = "B73xZ001"
  names(ped.hat1) <- gsub("Genotype2", "", names(ped.hat1))
  tped <- data.frame(Genotype = names(ped.hat1), trait = ped.hat1)
  names(tped)[2] <- trait
  return(tped)
}