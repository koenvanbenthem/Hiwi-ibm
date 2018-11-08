#Basic IBM
rm(list=ls())

##### PARAMETERS #####
Nt<-10000 #generations

<<<<<<< HEAD
a<-0.49649467

b<-1.47718931

c1<-0.72415095

c2<--0.24464625

c3<-0.99490196

c4<--1.31337296

c5<--0.06855583

c6<-0.32833236

c7<--20.88383990

c8<--0.66263785

c9<-2.39334027

c10<-0.11670283

=======
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

##### FUNCTIONS #####
>>>>>>> 12102607a97417bdd58941ca79fe039747574d42
w<-function(a,b,z,N,Np){
  y=a+b*plogis(c1+c2*N+c3*z+c4*(0.5*N-Np)+c5*N^2+c6*z^2+c7*(0.5*N-Np)^2+c8*z*N+
              c9*z*(0.5*N-Np)+c10*N*(0.5*N-Np))
  return(y)
  }

  
##### PATCHES #####
N1<-abs(round(rnorm(1, mean=50, sd=10)))
N2<-abs(round(rnorm(1, mean=50, sd=10)))
 
patch<-c(rep(1,N1),rep(2,N2))
trait<-c(rep(0.5,N1),rep(0.5,N2))
survival<-c(rep(1,N1),rep(1,N2))

pop<-data.frame(patch,trait,survival)
  
<<<<<<< HEAD
  patch<-c(rep(1,N1),rep(2,N2))
  
  trait<-c(rep(0.5,N1),rep(0.5,N2))
  
  survival<-c(rep(1,N1),rep(1,N2))
  
  pop<-data.frame(patch,trait,survival)
  
  #Time
  Nt<-3
  
  #Census->Offspring->Survival
###Generation-loop##Start##########  
  
  for(t in 1:Nt){
    
    N1<-nrow(subset(pop,pop[,1]==1))
    
    N2<-nrow(subset(pop,pop[,1]==2))
    
    N<-nrow(pop)
=======
  
##### GENERATION LOOP START #####  
for(t in 2:Nt){
  N1<-nrow(subset(pop,pop[,1]==1))
  N2<-nrow(subset(pop,pop[,1]==2))
  N<-nrow(pop)
>>>>>>> 12102607a97417bdd58941ca79fe039747574d42
   
  ##### OFFSPRING #####
  offspring<-c()
    
  if(N>0){
    for(i in 1:N){
<<<<<<< HEAD
     if(pop[i,1]<2){
       Nchild<-round(w(a,b,pop[i,2],1,N))
      offspring<-c(offspring,Nchild)
     }else{
       Nchild<-round(w(a,b,pop[i,2],2,N))
       offspring<-c(offspring,Nchild)
     }
      
    }
    Hera<-c()
    for(h in 1:N){
      Child<-c(rep(h,offspring[h]))
      Hera<-c(Hera,Child)
      }
    pop<-pop[c(1:nrow(pop),Hera),]
    
  }
   pop[1:N,3]<-pop[1:N,3]-1
   pop<-subset(pop,pop[,3]>0)
   
  }#generation end
  
  
=======
      if(pop[i,1]<2){
        Nchild<-round(w(a,b,pop[i,2],1,N))
        offspring<-c(offspring,Nchild)
      }else{
        Nchild<-round(w(a,b,pop[i,2],2,N)) #number of offspring the individual i becomes calculated with the fitness function
        offspring<-c(offspring,Nchild) #new number for the individual i is added to the already existing numbers of offspring from the individuals before
      }
    }
    Hera<-c() #empty vector
    for(h in 1:N){ 
      Child<-c(rep(h,offspring[h])) #replicates the number of the individual * the number of offspring it becomes
      Hera<-c(Hera,Child) #individual is so often named by its number, how many offspring it becomes
    }
    pop<-pop[c(1:nrow(pop),Hera),] #adds the clons of the individuals the the population data frame
  }
  ##### DEATH #####
  pop[1:N,3]<-pop[1:N,3]-1
  pop<-subset(pop,pop[,3]>0)
} ##### GENERATION LOOP END #####
>>>>>>> 12102607a97417bdd58941ca79fe039747574d42
