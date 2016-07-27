getTree<-function(table,plottype){
  t<-as.matrix(read.table(table))
  t<-data.frame(t)
  #remove alternate and reference
  tt<-t[c(-1,-2)]
  tt<-tt[c(-1,-2),]
  tt<-(100-tt)
  View(tt)
  colnames(tt)<-c("DH10","DH12","DH14","DH16","DH3","DH8","Kui2007","Ki11","Ki3","Ki43")
  rownames(tt)<-c("DH10","DH12","DH14","DH16","DH3","DH8","Ki2007","Ki11","Ki3","Ki43")
  arbol<-nj(as.dist(tt))
  plotcolors=c("darkorange3","darkorange3","darkorange3","darkorange3","darkorange3","darkorange3","black","black","black","black")
  plot(arbol,type=plottype,tip.color=plotcolors)
}