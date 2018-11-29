#Basic IBM
rm(list=ls())

##### PARAMETERS #####
replic<-2 #replicates
Nt<-10 #generations
mig <- 0.05 #migrationfactor
max.Age<-2 # age limit

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
N1<-abs(round(rnorm(1, mean=250, sd=10))) #patch 1 is drawn 
N2<-abs(round(rnorm(1, mean=250, sd=10))) #patch 2 is drawn
N1.m<-round(runif(1,N1/4,3*N1/4)) #males in patch 1
N2.m<-round(runif(1,N1/4,3*N1/4)) #males in patch 2

ID <- c(1:(N1+N2)) #vector ID: gives each individual an ID
patch<-c(rep(1,N1),rep(2,N2)) #vector patch: is filled with patch 1 (=1) and patch 2 (=2)
gender<-c(rep("male",N1.m),rep("female",N1-N1.m),rep("male",N2.m),rep("female",N2-N2.m)) #vector gender: is filled with males and females
trait<-c(rep(0.5,N1),rep(0.5,N2)) #vector trait: is for all individuals from both patches set as 5
survival<-c(rep(max.Age,N1),rep(max.Age,N2)) #vector survival: is for all new individuals of both patches 1

pop<-data.frame(ID,patch,gender,trait,survival) #data frame including all individuals out of both patches with the columns: patch, trait & survival


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
  
  population <- round(N1) + round(N2) #number of individuals
  loci <- matrix(NA,nrow=population,ncol=20+1) #empty matrix for the locis
  for(x in 1:population){ #for each individual
    loci[x,] <- round(runif(21,1,10)) #each individual has 20 random numbers (first 10:row //last 10:column)
    loci[x,21] <- x
  }
  
  values <- matrix(NA,nrow=population,ncol=10) #empty matrix for the trait values for each loci
  for(y in 1:population){ #for each individual
    for(z in 1:10){ 
      values[y,z] <- gen_phen_map[z,loci[y,z],loci[y,10+z]]
    }
    pop[y,4] <- abs(sum(values[y,])) ##### USE OF NR
  }
  
  
  ##### GENERATION LOOP START #####  
  for(t in 2:Nt){
    N1<-nrow(subset(pop,pop$patch==1)) #N1 is every generation overwritten to keep updated 
    N2<-nrow(subset(pop,pop$patch==2)) #N2 is every generation overwritten to keep updated
    N<-c(nrow(pop)) #how many individuals there are in both patches
    
    
    ##### MATRICES #####
    N.w <- subset(pop,pop$gender=="female") #female individuals in total
    N1.w <- subset(pop,pop$gender=="female"&pop$patch==1) #female individuals from patch 1
    N2.w <- subset(pop,pop$gender=="female"&pop$patch==2) #female individuals from patch 2
    N1.m <- subset(pop,pop$gender=="male"&pop$patch==1) #male individuals from patch 1
    N2.m <- subset(pop,pop$gender=="male"&pop$patch==2) #male individuals from patch 
    
    
    ##### OFFSPRING #####
    N.0<-N/500
    N.l <- c(N1/500,N2/500) # vector of local population sizes
    
    if(nrow(N.w)>0){ #number of offspring per female
      Nchild <- rpois(nrow(N.w),w(a,b,N.w$trait,N.0,N.l[N.w$patch])) #each female gets a random number of offspring
    }
    
    ID.children <- c(rep(0,nrow=sum(Nchild))) #empty vector for the ID
    patch.children <- c(rep(0,nrow=sum(Nchild))) #empty vector for the patch
    gender.children <- c(rep(0,nrow=sum(Nchild))) #empty vector for the gender
    trait.children <- c(rep(0,nrow=sum(Nchild))) #empty vector for the trait
    survival.children <- c(rep(max.Age,nrow=sum(Nchild))) #each child gets the survival of the maximum age
    pop.new <- data.frame(ID.children,patch.children,gender.children,trait.children,survival.children)
    
    loci.new <- matrix(NA,nrow=sum(Nchild),ncol=21) #empty matrix: children locis
    
    #### LOOP PARTNERFINDING #####
    for(u in 1:nrow(N.w)){ #start loop
      
      ##### FEMALES PATCH 1 #####
      if(N.w[u,2]<2){ ###==1
        father <- sample(N1.m,size=1) #samples one individual out of the patch
        
        #GENETICS PATCH 1:
        loci.mother <- loci[nrow(N1.m)+u,] #vector of locis of the mother
        loci.father <- loci[u,] #vector of locis of the father
        loci.child <- rep(0,ncol(loci)) #empty vector with fixed legth
        
        for(o in 1:Nchild[u]){ #for loop for the number of children per female
          for(p in 1:(10)){ #loop over the 10 locis
            if(runif(1,0,1)>0.5){ #if the random number is higher then 0.5:
              loci.child[p] <- loci.mother[p] #child gets the top allel (spot p) from mother
            } else{
              loci.child[p] <- loci.mother[10+p] #child gets the bottom allel (spot 10+p) from mother
            }
            if(runif(1,0,1)>0.5){ #if the random number is higher then 0.5:
              loci.child[10+p] <- loci.father[p] #child gets the top allel (spot p) from father
            } else{
              loci.child[10+p] <- loci.father[10+p] #child gets the bottom allel (spot 10+p) from mother
            }
          } #end loop 10 locis
          
          #FILLS CHILDREN MATRIX
          loci.new[o,] <- loci.child #Loci of the child are written into the matrix for the children loci
          loci.child[o,21] <- (nrow(pop))+o #fills in the ID in the children loci matrix
          pop.new[o,2] <- 1 #each child gets the patch 1
          pop.new[o,4] <- abs(sum(loci.child[1:20])) #each child gets their traitvalue
          if(runif(1,0,1)>0.5){ #if random number is higher als 0.5, child is female
            pop.new[o,3] <- "female"
          } else{
            pop.new[o,3] <- "male"
          }
          pop.new[o,1] <- (nrow(pop))+o #fills in the ID in the children pop matrix of the individual
        } #end loop number children
      } #end females patch 1
      
      
      #### FEMALES PATCH 2 #####
      if(N.w[u,2]>1){  ###==2
        father <- sample(N2.m,size=1) #samples one individual out of the patch 
        
        #GENETICS PATCH 2:
        loci.mother <- loci[nrow(N1.m)+nrow(N1.w)+nrow(N2.m)+u,] #vector of locis of the mother
        loci.father <- loci[nrow(N1.m)+nrow(N1.w)+u,] #vector of locis of the father
        loci.child <- rep(0,ncol(loci)) #empüty vector with fixed legth
        
        for(q in 1:Nchild[u]){ #for loop for the number of children per female
          for(s in 1:10){ #loop over the 10 locis
            if(runif(1,0,1)>0.5){ #if the random number is higher then 0.5:
              loci.child[q] <- loci.mother[q] #child gets the top allel (spot p) from mother
            } else{
              loci.child[q] <- loci.mother[10+q] #child gets the bottom allel (spot 10+p) from mother
            }
            if(runif(1,0,1)>0.5){ #if the random number is higher then 0.5:
              loci.child[10+q] <- loci.father[q] #child gets the top allel (spot p) from father
            } else{
              loci.child[10+q] <- loci.father[10+q] #child gets the bottom allel (spot 10+p) from mother
            }
            loci.child[21] <- (nrow(pop))+o
          } #end loop 10 locis
          
          #FILLS CHILDREN MATRIX
          loci.new[o+q,] <- loci.child #Loci of the child are written into the matrix for the children loci
          loci.child[o+q,21] <- (nrow(pop))+o+q #fills in the ID in the children loci matrix
          pop.new[o+q,2] <- 2 #each child gets the patch 2
          pop.new[o+q,4] <- abs(sum(loci.child[1:20])) #each child gets their traitvalue
          if(runif(1,0,1)>0.5){ #if random number is higher als 0.5, child is female
            pop.new[o+q,3] <- "female"
          } else{
            pop.new[o+q,3] <- "male"
          }
          pop.new[o+q,1] <- (nrow(pop))+o+q #fills in the ID in the children pop matrix of the individual
        } #end loop number children
      } #end females patch 2
      
    } #END LOOP PARTNERFINDING
    pop<-pop[c(1:nrow(pop),pop.new),] #adds the children to the population data frame
    
    
    ##### DEATH #####
    pop$survival[1:N]<-pop$survival[1:N]-1 #survival set on 0
    
    for(v in 1:nrow(pop)){ #for each individual
      if(pop$survival<=0){ #if the survival is 0, it replaces the first loci with -2
        locis[v,1] <- -2
      }
    }
    
    locis <- subset(locis,locis[,1]>-2 ) #all rows with a -2 in the beginning are deleted
    pop <-subset(pop,pop$survival>0) #Individuals which have a survival higher then 0 stay alive in the dataframe
    ##### END DEATH #####
    
    
    ##### MIGRATION START #####
    mig.N1 <- runif(nrow((subset(pop,pop$patch==1))),0,1) #draws so many uniformmly distributed numbers as there are individuals in patch 1
    mig.N2 <- runif(nrow((subset(pop,pop$patch==2))),0,1) #draws so many uniformmly distributed numbers as there are individuals in patch 2
    
    mig.N1 <- ifelse(mig.N1>mig,1,2) #the individuals with a random number lower then the migration rate get the value 2 (migrates to patch 2) & and the ones higher as the migration rate get the value 1 (dont migrate, stay in patch 1)
    mig.N2 <- ifelse(mig.N2>mig,2,1) #the individuals with a random number lower then the migration rate get the value 1 (migrate to patch 1) & and the ones higher as the migration rate get the value 2 (dont migrate,stay in patch 2)
    
    migration<-c(mig.N1,mig.N2)
    pop$patch<-migration
    
    chaos<-order(pop$patch) #orderd after patches
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

