
backcalculate.xf <- function(lambda0,D,p,n){
  
  m=length(D);
  time=(m+1-length(lambda0)):m; #time是lambda0的每个元素对应的时间点，D的元素对应的时间点是1~m
  
  #用D1,...,Dm和p确定t1和t2
  for(i in 0:-1000){
    Q1=sum(p[(2+abs(i)):(m+1+abs(i))]);
    if(Q1<=(5/100)){break};#探测得到的信息不得少于5%（这个值自己设定）
    t1=i;
  }
  
  for(i in 1:m){
    Q2=sum(p[1:(m+1-i)]);
    if(Q2<=(5/100)){break};
    t2=i;
  }
  
  q=rep(0,t2-t1+1);
  for(j in 1:(2-t1)){q[j]=sum(p[(2-(t1+j-1)):(m+1-(t1+j-1))])}
  for(j in (2-t1+1):length(q)){q[j]=sum(p[1:(m-(j+t1-2))])}
  
  #迭代初值
  Dn=D;
  
  nt1=which(time==t1);nt2=which(time==t2);
  lambda=lambda0[nt1:nt2];
  
  for(i in 1:t2){Dn[i]=sum(lambda[1:(i+1-t1)]*p[(i+1-t1):1])}
  if((t2+1)<=m){
    for(i in (t2+1):m){Dn[i]=sum(lambda*p[(i+1-t1):(i+1-t2)])}
  }
  
  
  for(t in 1:n){
    
    a=lambda; b=Dn;
    
    y=rep(0,length(lambda));
    for(l in 1:(2-t1)){y[l]=sum(p[(2-(t1+l-1)):(m+1-(t1+l-1))]*D/b)}
    for(l in (3-t1):length(y)){y[l]=sum(c(rep(0,l-(2-t1)),p[1:(m+(2-t1)-l)])*D/b)}
    lambda=a/q*y
    
    for(i in 1:t2){Dn[i]=sum(lambda[1:(i+1-t1)]*p[(i+1-t1):1])}
    
    if((t2+1)<=m){for(i in (t2+1):m){Dn[i]=sum(lambda*p[(i+1-t1):(i+1-t2)])}}
    
    e=sum((Dn-D)^2/Dn)/length(D);
    
    
    plot(t1:m,t1:m,pch='',ylim=c(0,max(lambda*1.1)));lines(t1:t2,lambda);lines(1:m,D,col=2);lines(1:m,Dn,col=3);
    legend("topright", legend = c("infection","reported case", "estimated reported case"),lwd=1,col = c(1,2,3))
  }
  
  return(list(error=e,time=c(t1,t2),estimation=lambda));
}





