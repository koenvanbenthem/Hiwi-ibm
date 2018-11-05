#Basic IBM
rm(list=ls())

#fecundity

a<-0.5
b<-0.5

c1<-1
c2<-1
c3<-1
c4<-1
c5<-1
c6<-1
c7<-1
c8<-1
c9<-1
c10<-1


w<-function(a,b,z,N,Np){
  y=a+b*log(c1+c2*N+c3*z+c4*(0.5*N-Np)+c5*N^2+c6*z^2+c7*(0.5*N-Np)^2+c8*z*N+
              c9*z*(0.5*N-Np)+c10*N*(0.5*N-Np))
return(y)
  }

  
  #Patches
  N1<-round(rnorm(1, mean=50, sd=10))
 
  N2<-round(rnorm(1, mean=50, sd=10))
 
  
  patch<-c(rep(1,N1),rep(2,N2))
  
  trait<-c(rep(5,N1),rep(5,N2))
  
  survival<-c(rep(1,N1),rep(1,N2))
  
  pop<-data.frame(patch,trait,survival)
  
  #Time
  Time<-10000
  
  #Census->Offspring->Survival
###Generation-loop##Start##########  
  
  for(i in 1:Time){
    
    N1<-nrow(subset(old,old[,1]==1))
    
    N2<-nrow(subset(old,old[,1]==2))
    
    N<-nrow(old)
    
    newp<-c()
    newt<-c()
    newa<-c()
    
    
    #offspring
    for(z in 2:Time){
    new<- data.frame(nrow=0,ncol=3)
    newp<-c()
    newt<-c()
    newa<-c()
    }
    
  if(N>0){
    for(i in 1:N){
     # Nchild<-round(w(a,b,old[i,2],N,))
      newp<-c(newp,(rep(old$patch[i],Nchild)))
      newt<-c(newt,(rep(old$trait[i],Nchild)))
      newa<-c(newa,(rep(old$age[i],Nchild)))
      
    }
  }
   pop<-
     
  }#generation end
  
  