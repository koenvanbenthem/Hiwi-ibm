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
  
  trait<-c(rep(0.5,N1),rep(0.5,N2))
  
  survival<-c(rep(1,N1),rep(1,N2))
  
  pop<-data.frame(patch,trait,survival)
  
  #Time
  Nt<-10000
  
  #Census->Offspring->Survival
###Generation-loop##Start##########  
  
  for(i in 1:Nt){
    
    N1<-nrow(subset(pop,pop[,1]==1))
    
    N2<-nrow(subset(pop,pop[,1]==2))
    
    N<-nrow(pop)
   
    
    
    #offspring
 
   offspring<-c()
    
  if(N>0){
    for(i in 1:N){
     if(pop[i,1]<2){
       Nchild<-round(w(a,b,pop[i,2],1,N))
      offspring<-c(offspring,Nchild)
     }else{
       Nchild<-round(w(a,b,pop[i,2],2,N))
       offspring<-c(offspring,Nchild)
     }
     }
      }
    
    
  }
     
  }#generation end
  
  pop<-pop[c(1:nrow(pop),1,1,1,2,2,4,5),]
  
  #Death
  pop[1:N,3]<-pop[1:N,3]-1
  pop<-subset(pop,pop[,3]>0)
  