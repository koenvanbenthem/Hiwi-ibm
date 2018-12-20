#Basic IBM
rm(list=ls())
switch(Sys.info()['user'],
       bio23user = {setwd("/home/bio23user/Documents/Projects/Hiwi-ibm/Hiwi-ibm/")},
       Leron = {setwd("C:/Users/Leron/Desktop/IBM_code/")})

source('Gene_generator.R')
#USED INDIC######
#   replication -r
#   loci/traitv.-x,y,z    
#   genertion   -t
#   partner     -u
#   gentics     -o,p
#   survival    -v
#   gender loop -g
#   patches     -k
#   values      -d,f
#   offspring   -ß,j
#   migration   -g
##### PARAMETERS #####
replic<-1 #replicates
Nt<-100 #generations
mig <- 0.05 #migrationfactor
max.Age<-2 # age limit
patches<-3 # Number of Patches
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
pop<-c()
for(k in 1:patches){

N_patchx<-abs(round(rnorm(1, mean=250, sd=10))) #Number of individuals in the patch 

patchx_m<-round(runif(1,N_patchx/4,3*N_patchx/4)) #Number of males in the patch

ID <- c(1:(N_patchx)) #vector ID: gives each individual an ID
patch<-c(rep(k,N_patchx)) #vector patch: gives each individual their patch Nr.
gender<-c(rep("male",patchx_m),rep("female",N_patchx-patchx_m)) #vector gender: is filled with males and females
trait<-c(rep(0.5,N_patchx)) #vector trait: is for all individuals from both patches set as 5
survival<-c(rep(max.Age,N_patchx)) #vector survival: is for all new individuals of both patches 1

patchx<-data.frame(ID,patch,gender,trait,survival)
pop<-rbind(pop,patchx)  #data frame including all individuals of all patches
}

pop$ID<-c(1:nrow(pop))#new ID for the population
home<-c(1:patches)#vector of patchnr.
##### VECTORS #####
Npop <- rep(0,Nt) #empty vector for the populationsize of each generation in patch 1
#pop.N2.vector <- rep(0,Nt) #empty vector for the populationsize of each generation in patch 2
#trait.N1.vector <- rep(0,Nt) #empty vector for the average trait-value of each generation in patch 1
#trait.N2.vector <- rep(0,Nt) #empty vector for the average trait-value of each generation in patch 2

Npop[1] <- nrow(pop) #populationsize for the first generation of patch 1
#pop.N2.vector[1] <- N2 #populationsize for the first generation of patch 2
#trait.N1.vector <- mean(pop$trait[pop$patch==1]) #average trait-value for the first generation of patch 1
#trait.N2.vector <- mean(pop$trait[pop$patch==2]) #average trait-value for the first generation of patch 2


##### REPLICATION LOOP START#####

for(r in 1:replic){

  population <- nrow(pop) #number of individuals
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
    pop[y,4] <- abs(sum(values[y,])) ##### USE OF COLUMN.NR
  }
  
  
  ##### GENERATION LOOP START #####  
  for(t in 1:Nt){
    ######IS ANYBODY THERE? START##############
    if(nrow(pop)>0){
      
      
    N<-c(nrow(pop)) #how many individuals there are in both patches
    N.x<-c()
    ##### MATRICES #####
    N.w <- subset(pop,pop$gender=="female") #female individuals in total
    N.m<- subset(pop,pop$gender=="male") #male individuals in total
    
    ##### OFFSPRING #####
    chance<-c()
    for(pls in 1:patches){###IS OFFSPRING POSSIBLE START
      chance<-c(chance,nlevels(subset(pop,pop$patch==pls)$gender))#
    }#create a Vector which shows how many different arguments(levels) are in a vector 
    if(max(chance)==2){#if one patch contains both genders then it has a level of 2
    
    N.0<-N/500
    for(j in 1:patches){
      N.x<-c(N.x,nrow(subset(pop,pop$patch==j))/500)
    }
    # vector of local population sizes
    
    if(nrow(N.w)>0){ #number of offspring per female
      Nchild <- 2*rpois(nrow(N.w),w(a,b,N.w$trait,N.0,N.x[N.w$patch])) #each female gets a random number of offspring
    }
    
    ID.children <- c() #empty vector for the ID
    patch.children <- c() #empty vector for the patch
    gender.children <- c() #empty vector for the gender
    trait.children <- c() #empty vector for the trait
    survival.children <- c() #each child gets the survival of the maximum age
    
    
    loci.new <- c() #empty vector: children locis
    
    #### START LOOP PARTNERFINDING #####
    patchbook <- c()
    gendergram <- c()
    if(nrow(N.w)>0){#ANY FEMALES START#####
    for(u in 1:nrow(N.w)){ #loop mother 
      if(Nchild[u]>0){ #just if the mother gets offspring
        mother<-N.w$ID[u] #gives the ID of the mother
        
        ###FATHER####
        if(nrow(subset(N.m,N.m$patch==N.w$patch[u]))>0){#ANY MALES IN THE PATCH OF THE MOTHER? START
        # sample the ID of one male which patchnr. is the same as the patchnr. of the mother
        father<-sample(subset(N.m$ID,N.m$patch==subset(N.w$patch,N.w$ID==mother)),1)
        #GENETICS:
        loci.mother <- subset(loci,loci[,21]==mother) #vector of locis of the mother
        loci.father <- subset(loci,loci[,21]==father) #vector of locis of the father
        loci.child <- rep(0,ncol(loci)) #empty vector with fixed length
        
        if(Nchild[u]>0){ #just if the mother gets offspring
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
          loci.new <-  rbind(loci.new,loci.child) #connects loci of the child to the matrix of the other children
          
          if(runif(1,0,1)>0.5){ #if random number is higher als 0.5, child is female
            gendergram <- c(gendergram,"female")  
          } else{ #it is male
            gendergram <- c(gendergram,"male")     
          }
        }#END LOOP NUMBER CHILDREN
        } 
        patchbook <- c(patchbook, rep(subset(pop,pop$ID==mother)[2],Nchild[u]),recursive=TRUE) #each kid gets the patch of the mother
        }#END ANY MALES?
        }#does the mother get offspring
    } #END LOOP PARTNERFINDING/mother
    
    ID.children <- c(rep(0,length(patchbook)))
    trait.children <- c(rep(0,length(patchbook))) 
    survival.children <- c(rep(max.Age,length(patchbook))) #each child gets the survival of the maximum age
    gender.children <- gendergram #gender of the children are written into the matrix
    patch.children <- patchbook #patches of children are written into the matrix
    pop.new <- data.frame(ID.children,patch.children,gender.children,trait.children,survival.children)
    colnames(pop.new)<-c("ID","patch","gender","trait","survival") #colum names
    
    values.new <- matrix(NA,nrow=sum(Nchild),ncol=10) #empty matrix for the trait values for each loci
    for(d in 1:sum(Nchild)){ #for each individual offspring
      for(f in 1:10){ 
        values.new[d,f] <- gen_phen_map[f,loci.new[d,f],loci.new[d,10+f]]
      }
      pop.new[d,4] <- abs(sum(values.new[d,])) ##### USE OF COLUMN.NR
    }
    

    
    pop<-rbind(pop,pop.new)
    rownames(pop) <- 1:nrow(pop)
    loci<-rbind(loci,loci.new)
    }#END ANY FEMALES?
    }#END IS OFFSPRING POSSIBLE
    ##### DEATH START #####
    pop$survival[1:N]<-pop$survival[1:N]-1 #every adult loses one survival counter
    for(v in 1:nrow(pop)){ #for each individual
      if(pop[v,5]==0){ #if the survival is 0, it replaces the first loci with -2
        loci[v,1] <- -2
      }
    }
    
    loci <- subset(loci,loci[,1]>(-2 )) #all rows with a -2 in the beginning are deleted
    pop <-subset(pop,pop$survival>0) #Individuals which have a survival higher then 0 stay alive in the dataframe
    ##### END DEATH #####
    
    
    ##### MIGRATION START #####
    if(nrow(pop)>0){
    wanderlust<-runif(nrow(pop),0,1)# draws one uniformmly distribued number for every individual
    new.patch<-c()
    for(g in 1:length(wanderlust)){#for every individual
      if(wanderlust[g]<mig){#if wanderlust is lower than mig sample one patchnumber out of a vector of patchnr. other than your own
        new.patch<-c(new.patch,sample(subset(home,home!=(pop$patch[g])),1),recursive=TRUE)#samples one patchnr. out of a vector of patchnr. with the exception of the individuals own patch
      }else{# if mig is higher than wanderlust just put in your own patchnr.
        new.patch<-c(new.patch,pop$patch[g],recursive=TRUE)
      }
          }
    pop$patch<-new.patch#override the vector in the population
    
    rownames(pop) <- 1:nrow(pop) #re-indexing the population
    pop$ID<-c(1:nrow(pop))#new ID for the population
    loci[,21]<-c(1:nrow(pop))#new ID for the loci
    }
    ##### MIGRATION END #####
    
    
    Npop[t] <-nrow(pop) #overwrites the populationsizes for each generation in the empty vector (patch 1)
    #pop.N2.vector[t] <-sum(pop$patch==2) #overwrites the average trait-value for each generation in the empty vector (patch 2)
    #trait.N1.vector[t] <- mean(pop$trait[pop$patch==1]) #overwrites the average trait-value for each generation in the empty vector (patch 1)
    #trait.N2.vector[t] <- mean(pop$trait[pop$patch==2]) #overwrites the average trait-value for each generation in the empty vector (patch 2)
    }### IS ANYBODY THERE? END ######
    } ##### GENERATION LOOP END #####
  
  #pdf(paste("graph",r,".pdf",sep=""))
  #plot(pop.N1.vector, main="populationsize over the generations",xlab="generations",ylab="populationsize",type="l",col="darkorange3") #plot populationsize
  #lines(pop.N2.vector,type="l",col="green") #includes the populationsize of patch 2
  #legend("topright",legend=c("patch 1","patch 2"),lty=1,col=c("darkorange3","green"))
  #dev.off()
  #  print(r)
}
##### REPLICATION LOOP END#####


##### PLOTS #####

plot(Npop,main="Population over time", xlab="generations",ylab="population",type="l",col="red") #plot traitvalue
#lines(trait.N2.vector,type="l",col="blue") #includes the average trait-value of patch 2
#legend("topright",legend=c("patch 1","patch 2"),lty=1,col=c("red","blue"))



#######Order######
pop$ID<-c(1:nrow(pop))#new ID for the population
rownames(pop) <- 1:nrow(pop) #re-indexing the population to prevent 1.1.3.2.4.....
loci[,21]<-c(1:nrow(pop))#new ID for the loci
chaos<-order(pop$patch) #vector of indices to orderd after patches
pop<-pop[chaos,] #order the pop matrix


sorting <- c() #
for(h in 1:nrow(pop)){
  sorting <- rbind(sorting, subset(loci,loci[21]==pop[h,1]))
}
loci <- sorting

