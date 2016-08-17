

internalSD<-function(residuals){
  
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
  return(allsd)
  
  
  
}