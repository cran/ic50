preimage<-function(y0,x,y){
  n<-length(y)
  i<-1
  while(abs(y[i]-y0)>=0.005){
    i<-i+1
    if(i==n+1) return(NA)
  }
  return(x[i])
}

linearCurve<-function(x,y){ #Werte für stückweise lineare Kurve berechnen
  x<-round(x,2)
  xgr<-round(seq(from=min(x),to=max(x),by=.005),3)
  ylin<-numeric(0)
  n<-length(x)
  for(i in 2:n){
    m<-(y[i]-y[i-1])/(x[i]-x[i-1])
    section<-(x[i-1]<=xgr & xgr<=x[i])
    ylin[section]<-y[i-1]+m*(xgr[section]-min(xgr[section]))
  }
  return(list(xgr,ylin))  
}

expCurve<-function(x,y){
  xgr<-linearCurve(x,y)[[1]]
  ygr<-linearCurve(x,y)[[2]]

  a=1-min(ygr)
  x0=preimage(1-a/2,xgr,ygr)
  y0=min(ygr)
  b<-preimage(y0+a/(1+exp(-.5)),xgr,ygr)-preimage(y0+a/(1+exp(.5)),xgr,ygr)

  return(list(xgr,y0+(a/(1+exp(-(xgr-x0)/b)))))
}


defaultfiles.write<-function(){
  if(!file.exists(".last384_measure.txt")){
    measure384<-data(default384_measure)
    write.table(measure384,file=".last384_measure.txt",col.names=FALSE,sep="\t",row.names=FALSE)
  }
  if(!file.exists(".last384_control.txt")){
    control384<-data(default384_control)
    write.table(control384,file=".last384_control.txt",col.names=FALSE,sep="\t",row.names=FALSE)
  }
  if(!file.exists(".last384_dilution.txt")){
    dilution384<-data(default384_dilution)
    write.table(dilution384,file=".last384_dilution.txt",col.names=FALSE,sep="\t",row.names=FALSE)
  }
  if(!file.exists(".last96_measure.txt")){
    measure96<-data(default96_measure)
    write.table(measure96,file=".last96_measure.txt",col.names=FALSE,sep="\t",row.names=FALSE)
  }
  if(!file.exists(".last96_control.txt")){
    control96<-data(default96_control)
    write.table(control96,file=".last96_control.txt",col.names=FALSE,sep="\t",row.names=FALSE)
  }
  if(!file.exists(".last96_dilution.txt")){
    dilution96<-data(default96_dilution)
    write.table(dilution96,file=".last96_dilution.txt",col.names=FALSE,sep="\t",row.names=FALSE)
  }
}
