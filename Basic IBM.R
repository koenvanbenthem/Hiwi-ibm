#Basic IBM
rm(list=ls())

<<<<<<< HEAD

##### PARAMETERS #####
Time <- 10000 #generations to run

a<-0.5 #alpha
b<-0.5 #beta
=======
#fecundity

a<-0.5
b<-0.5

>>>>>>> Matthias
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


<<<<<<< HEAD
##### FUNCTIONS #####
w<-function(a,b,z,N,Np){ #Fitness-function
  y=a+b*log(c1+c2*N+c3*z+c4*(0.5*N-Np)+c5*N^2+c6*z^2+c7*(0.5*N-Np)^2+c8*z*N+c9*z*(0.5*N-Np)+c10*N*(0.5*N-Np)) # !!!!! #need to do 2 functions for Np1 and Np2
  return(y)
}

  
##### PATCHES #####
N1<-round(rnorm(1, mean=50, sd=10)) #patch 1 is drawn 
N2<-round(rnorm(1, mean=50, sd=10)) #patch 2 is drawn
 
patch<-c(rep(1,N1),rep(2,N2)) #vector patch: is filled with patch 1 (=1) and patch 2 (=2)
trait<-c(rep(5,N1),rep(5,N2)) #vector trait: is for all indicviduals from both patches set as 5
survival<-c(rep(1,N1),rep(1,N2)) #vector survival: is for all new individuals of both patches 1
  
pop<-data.frame(patch,trait,survival) #data frame including all individuals out of both patches with the columns: patch, trait & survival
  
  
##### GENERATION LOOP START #####  
for(i in 1:Time){
    
  N1<-nrow(subset(pop,pop[,1]==1)) #N1 is every generation overwritten to keep updated 
  N2<-nrow(subset(pop,pop[,1]==2)) #N2 is every generation overwritten to keep updated
  N<-nrow(pop) #how many individuals there are in total
    
  ##### OFFSPRING #####
  for(z in 2:Time){
    
  }
    
  if(N>0){
    for(i in 1:N){
     # Nchild<-round(w(a,b,old[i,2],N,))
      
    }
  }
}
  ##### GENERATION LOOP END #####
  
  
  
  pop[c(1:nrow(pop),1,1,1,2,2,3,4),] #takes the old data frame and adds 
=======
w<-function(a,b,z,N,Np){
  y=a+b*plogis(c1+c2*N+c3*z+c4*(0.5*N-Np)+c5*N^2+c6*z^2+c7*(0.5*N-Np)^2+c8*z*N+
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
>>>>>>> Matthias
  