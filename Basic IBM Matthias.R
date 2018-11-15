#Basic IBM
rm(list=ls())

##### PARAMETERS #####
replic<-2
Nt<-10 #generations
mig <- 0.05 #migrationfactor

#fecundity
a <-  0.49649467
b  <-  1.47718931
c1    <-  0.72415095
c2    <- -0.24464625
c3    <-  0.99490196
c4    <- -1.31337296
c5    <- -0.06855583
c6    <-  0.32833236
c7    <--20.88383990
c8    <- -0.66263785
c9    <-  2.39334027
c10   <-  0.11670283

##### FUNCTIONS #####
w<-function(a,b,z,N,Np){
  y=a+b*plogis(c1+c2*N+c3*z+c4*(0.5*N-Np)+c5*N^2+c6*z^2+c7*(0.5*N-Np)^2+c8*z*N+c9*z*(0.5*N-Np)+c10*N*(0.5*N-Np))
  return(y)
  }

  
##### PATCHES #####
N1<-abs(rnorm(1, mean=250, sd=10)) #patch 1 is drawn 
N2<-abs(rnorm(1, mean=250, sd=10)) #patch 2 is drawn
 
patch<-c(rep(1,N1),rep(2,N2)) #vector patch: is filled with patch 1 (=1) and patch 2 (=2)
trait<-c(rep(0.5,N1),rep(0.5,N2)) #vector trait: is for all indicviduals from both patches set as 5
survival<-c(rep(1,N1),rep(1,N2)) #vector survival: is for all new individuals of both patches 1
  
pop<-data.frame(patch,trait,survival) #data frame including all individuals out of both patches with the columns: patch, trait & survival


##### VECTORS #####
pop.N1.vector <- rep(0,Nt) #empty vector for the populationsize of each generation in patch 1
pop.N2.vector <- rep(0,Nt) #empty vector for the populationsize of each generation in patch 2
trait.N1.vector <- rep(0,Nt) #empty vector for the average trait-value of each generation in patch 1
trait.N2.vector <- rep(0,Nt) #empty vector for the average trait-value of each generation in patch 2

pop.N1.vector[1] <- N1 #populationsize for the first generation of patch 1
pop.N2.vector[1] <- N2 #populationsize for the first generation of patch 2
trait.N1.vector <- mean(pop$trait[pop$patch==1]) #average trait-value for the first generation of patch 1
trait.N2.vector <- mean(pop$trait[pop$patch==2]) #average trait-value for the first generation of patch 2


##### REPLICATION LOOP START#####

for(r in 1:replic){
  

##### GENERATION LOOP START #####  
for(t in 2:Nt){
  N1<-nrow(subset(pop,pop$patch==1)) #N1 is every generation overwritten to keep updated 
  N2<-nrow(subset(pop,pop$patch==2)) #N2 is every generation overwritten to keep updated
  N<-c(nrow(pop)) #how many individuals there are in both patches


  ##### OFFSPRING #####
  N.0<-N/500
  N.l <- c(N1/500,N2/500) # vector of local population sizes
  
  if(N>0){
    Nchild <- rpois(nrow(pop),w(a,b,pop$trait,N.0,N.l[pop$patch])) #vectoe with number of offspring for each individual out of a poisson distribution
  }
    Hera <- rep(1:nrow(pop),Nchild) #replicates the individual to clone times the number it becomes offspring
    pop<-pop[c(1:nrow(pop),Hera),] #adds the clons of the individuals the the population data frame


  ##### DEATH #####
  pop$survival[1:N]<-pop$survival[1:N]-1 #survival set on 0
  pop <-subset(pop,pop$survival>0)
  ##### END DEATH #####


  ##### MIGRATION START #####
  mig.N1 <- runif(nrow((subset(pop,pop$patch==1))),0,1) #draws so many uniformmly distributed numbers as there are individuals in patch 1
  mig.N2 <- runif(nrow((subset(pop,pop$patch==2))),0,1) #draws so many uniformmly distributed numbers as there are individuals in patch 2

  mig.N1 <- ifelse(mig.N1>mig,1,2) #the individuals with a random number lower then the migration rate get the value 2 (migrates to patch 2) & and the ones higher as the migration rate get the value 1 (dont migrate, stay in patch 1)
  mig.N2 <- ifelse(mig.N2>mig,2,1) #the individuals with a random number lower then the migration rate get the value 1 (migrate to patch 1) & and the ones higher as the migration rate get the value 2 (dont migrate,stay in patch 2)
  
  migration<-c(mig.N1,mig.N2)
  pop$patch<-migration
  
  chaos<-order(pop$patch)
  pop<-pop[chaos,]
  
  
  #### !!! ##### all individuals with a 1 need to migrate in the other patch #### !!! ####
  ##### MIGRATION END #####
  
  
  pop.N1.vector[t] <-sum(pop$patch==1) #overwrites the populationsizes for each generation in the empty vector (patch 1)
  pop.N2.vector[t] <-sum(pop$patch==2) #overwrites the average trait-value for each generation in the empty vector (patch 2)
 
  trait.N1.vector[t] <- mean(pop$trait[pop$patch==1]) #overwrites the average trait-value for each generation in the empty vector (patch 1)
  trait.N2.vector[t] <- mean(pop$trait[pop$patch==2]) #overwrites the average trait-value for each generation in the empty vector (patch 2)

  
  rownames(pop) <- 1:nrow(pop)        #re-indexing the population to prevent 1.1.3.2.4.....
} 
##### GENERATION LOOP END #####

#pdf(paste("graph",r,".pdf",sep=""))
#plot(pop.N1.vector, main="populationsize over the generations",xlab="generations",ylab="populationsize",type="l",col="darkorange3") #plot populationsize
#lines(pop.N2.vector,type="l",col="green") #includes the populationsize of patch 2
#legend("topright",legend=c("patch 1","patch 2"),lty=1,col=c("darkorange3","green"))
#dev.off()
#  print(r)
}
##### REPLICATION LOOP END#####


##### PLOTS #####

plot(trait.N1.vector,main="average trait-value over the generations", xlab="generations",ylab="average trait-value",type="l",col="red") #plot traitvalue
lines(trait.N2.vector,type="l",col="blue") #includes the average trait-value of patch 2
legend("topright",legend=c("patch 1","patch 2"),lty=1,col=c("red","blue"))
