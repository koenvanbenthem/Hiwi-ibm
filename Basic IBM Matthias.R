#Basic IBM
rm(list=ls())

Morpheus<-function(replic=1,
                   Nt=1000, #generations
                   mig=0.05, #migrationfactor
                   max.Age=2, # age limit
                   patches=3, # Number of Patches
                   mutate=0.05, #mutationfactor
                   DIE=0.05, #Chance to die
                   
                   #fecundity
                   a=0.49649467,
                   b=1.47718931,
                   c1=0.72415095,
                   c2=-0.24464625,
                   c3=0.99490196,
                   c4=-1.31337296,
                   c5=-0.06855583,
                   c6 = 0.32833236,
                   c7=-20.88383990,
                   c8=-0.66263785,
                   c9=2.39334027,
                   c10=0.11670283
){
switch(Sys.info()['user'],
       bio23user = {setwd("/home/bio23user/Documents/Projects/Hiwi-ibm/Hiwi-ibm/")},
       Leron = {setwd("C:/Users/Leron/Desktop/IBM_code/")},
       Anwender = {setwd("C:/Users/Anwender/Desktop/")})

source('Gene_generator.R')

# library(profvis)
# profvis({
  
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
#   offspring   -pls,j
#   migration   -g

##### FUNCTIONS #####
w<-function(a,b,z,N,Np){
  y=a+b*plogis(c1+c2*N+c3*z+c4*(0.5*N-Np)+c5*N^2+c6*z^2+c7*(0.5*N-Np)^2+c8*z*N+c9*z*(0.5*N-Np)+c10*N*(0.5*N-Np))
  return(y)
}


Here_is_your_ID<-function(Nchild){ #ID function
  ID.children <-   ID_scan:(ID_scan+sum(Nchild)-1)
  ID_scan <<- ID_scan + sum(Nchild)
  return(ID.children)
}


f <- function(row,pop.matrix,value.matrix,loci.matrix){ #trait value function
  value.matrix <- matrix(NA,nrow=row,ncol=10) #empty matrix for the trait values for each loci
  for(y in 1:row){ #for each individual
    for(z in 1:10){ 
      value.matrix[y,z] <- gen_phen_map[z,loci.matrix[y,z],loci.matrix[y,10+z]]
    }
    pop.matrix[y,4] <- abs(sum(value.matrix[y,]))
  }
  return(pop.matrix)
}


stat.fun <- function(pop.matrix, Npatch){ #patch/statistic function
  tmp <- aggregate(pop.matrix$trait,by=list(patch = pop.matrix$patch),mean)
  traits <- tmp$x[match(1:Npatch,tmp$patch)] 
  
  cbind(table(factor(pop.matrix$patch,levels=1:Npatch)),
        table(factor(pop.matrix[pop.matrix$gender=='male',]$patch,levels=1:Npatch)),
        table(factor(pop.matrix[pop.matrix$gender=='female',]$patch,levels=1:Npatch)),
        as.numeric(traits))
}


##### REPLICATION LOOP START#####
for(r in 1:replic){
  
  
  ##### INITIALISATION PATCHES #####
  pop<-c()
  stats <- array(NA,dim=c(patches,4,Nt)) #empty array for the statistics
  
  for(k in 1:patches){
  
    N_patchx<-abs(round(rnorm(1, mean=250, sd=10))) #Number of individuals in the patch 
    patchx_m<-round(runif(1,N_patchx/4,3*N_patchx/4)) #Number of males in the patch
    
    ID <- c(1:(N_patchx)) #vector ID: gives each individual an ID
    patch<-c(rep(k,N_patchx)) #vector patch: gives each individual their patch Nr.
    gender<-c(rep("male",patchx_m),rep("female",N_patchx-patchx_m)) #vector gender: is filled with males and females
    trait<-c(rep(0.5,N_patchx)) #vector trait: is for all individuals from both patches set as 5
    survival<-c(rep(max.Age,N_patchx)) #vector survival: is for all new individuals of both patches 1
    ID.mother <- c(rep(NA,N_patchx))
    ID.father <- c(rep(NA,N_patchx))
    
    patchx<-data.frame(ID,patch,gender,trait,survival,ID.mother,ID.father)
    pop<-rbind(pop,patchx)  #data frame including all individuals of all patches
  }
  pop$ID<-c(1:nrow(pop))#new ID for the population
  home<-c(1:patches)#vector of patchnr.
 
  ID_scan<-nrow(pop)+1
  
  
  ##### Statistic #####
  Npop <- rep(0,Nt) #empty vector for the populationsize of each generation 
  Trait.pop <- rep(0,Nt) #empty vector for the mean traitvalue of each generation
  
  Npop[1] <- nrow(pop) #populationsize for the first generation 
  Trait.pop[1] <- mean(pop$trait) #mean traitvalue for the first generation 
  ########Statistic End#####
  
  
    population <- nrow(pop) #number of individuals
    loci <- matrix(NA,nrow=population,ncol=20+1) #empty matrix for the locis
    
    for(x in 1:population){ #for each individual
      loci[x,] <- ceiling(runif(21,1e-16,10)) #each individual has 20 random numbers (first 10:row //last 10:column)
      loci[x,21] <- x
    }
    
    pop <- f(population,pop,values,loci) #traitvalues for the population
    
  
    ##### GENERATION LOOP START #####  
    for(t in 1:Nt){
      N<- nrow(pop) #how many individuals there are in all patches

      
      ######IS ANYBODY THERE? START##############
      if(N>0) {
        N.x<-c()
        
        
        ##### MATRICES #####
        N.w <- subset(pop,pop$gender=="female") #female individuals in total
        N.m<- subset(pop,pop$gender=="male") #male individuals in total
        
        
        ##### OFFSPRING #####
        chance<-c()
        
        for(pls in 1:patches){#START IS OFFSPRING POSSIBLE?
          chance<-c(chance,nlevels(subset(pop,pop$patch==pls)$gender))#
        }#create a Vector which shows how many different arguments(levels) are in a vector 
        
        if(max(chance)==2){#if one patch contains both genders then it has a level of 2
          N.0<-N/500
          for(j in 1:patches){
            N.x<-c(N.x,nrow(subset(pop,pop$patch==j))/500) # vector of local population sizes
          }
          if(nrow(N.w)>0){ #number of offspring per female
            Nchild <- 2*rpois(nrow(N.w),w(a,b,N.w$trait,N.0,N.x[N.w$patch])) #each female gets a random number of offspring
          }
        
          ID.children <- c() #empty vector for the ID
          patch.children <- c() #empty vector for the patch
          gender.children <- c() #empty vector for the gender
          trait.children <- c() #empty vector for the trait
          survival.children <- c() #each child gets the survival of the maximum age
          ID.mother.children <- c() #empty vector for the mothers ID
          ID.father.children <- c() #empty vector for the fathers ID
          
          loci.new <- matrix(NA,nrow=sum(Nchild),ncol=21) #empty vector: children locis
          
          
          #### START LOOP PARTNERFINDING #####
          patchbook <- c()
          gendergram <- c()
          
          N.w.patch <- table(factor(N.w$patch,levels = 1:patches))# number of females in each patch (as a vector) (N.w.patch[k], number of females in patch k)
          N.m.patch <- table(factor(N.m$patch,levels = 1:patches))
          curr_child <- 1 # counter that keeps track of how many kids have emerged so far during the loop below
          
          if(nrow(N.w)>0){ #START ANY FEMALES?
            for(u in 1:nrow(N.w)){ #START LOOP PARTNERFINDING/mother 
              if(Nchild[u]>0){ #START GETS THE MOTHER OFFSPRING?
                mother<-N.w$ID[u] #gives the ID of the mother
                ID.mother.children <- c(ID.mother.children, rep(mother,Nchild[u])) #ID of the mother is written into the vector for all her offspring
                
              
                ###FATHER####
                if(N.m.patch[N.w$patch[u]]>0){ #START ANY MALES IN THE PATCH OF THE MOTHER?
                  father<-sample(N.m$ID[N.m$patch==N.w$patch[u]],1) #sample the ID of one male which patchnr. is the same as the patchnr. of the mother
                  ID.father.children <- c(ID.father.children,rep(father,Nchild[u])) #ID of the father is written into the vector
                  
                  #GENETICS:
                  loci.mother <- loci[loci[,21]==mother,] #vector of locis of the mother
                  loci.father <- loci[loci[,21]==father,] #vector of locis of the father
                  loci.child <- rep(0,ncol(loci)) #empty vector with fixed length
                  
                  if(Nchild[u]>0){ #START GETS THE MOTHER OFFSPRING?
                    for(o in 1:Nchild[u]){ #START LOOP NUMBER CHILDREN per female
                      loci.child[1:10] <- loci.mother[(1:10) +sample(c(0,10),10,replace=TRUE)]
                      loci.child[11:20] <- loci.father[(1:10) +sample(c(0,10),10,replace=TRUE)]
                      
                      #MUTATION
                      if(runif(1,0,1)<mutate){
                        loci.child[round(runif(1,1,20))]<-round(runif(1,1,10))
                      }
                      
                      loci.new[curr_child,] <-  loci.child #connects loci of the child to the matrix of the other children
                      curr_child <- curr_child + 1
            
                      if(runif(1,0,1)>0.5){ #if random number is higher als 0.5, child is female
                        gendergram <- c(gendergram,"female")  
                      } else{ #it is male
                        gendergram <- c(gendergram,"male")     
                      }
                      
                    } #END LOOP NUMBER CHILDREN
                  } #END GETS THE MOTHER OFFSPRING?
              } #END ANY MALES IN THE PATCH OF THE MOTHER?
            } #END GETS THE MOTHER OFFSPRING?
          } #END LOOP PARTNERFINDING/mother
          
          patchbook <- rep(N.w$patch,Nchild) #each kid gets the patch of the mother
          ID.children <- Here_is_your_ID(Nchild)
          trait.children <- c(rep(0,length(patchbook))) 
          survival.children <- c(rep(max.Age,length(patchbook))) #each child gets the survival of the maximum age
          gender.children <- gendergram #gender of the children are written into the matrix
          patch.children <- patchbook #patches of children are written into the matrix
          pop.new <- data.frame(ID.children,patch.children,gender.children,trait.children,survival.children,ID.mother.children,ID.father.children)
          colnames(pop.new)<-c("ID","patch","gender","trait","survival","ID.mother","ID.father") #colum names
          loci.new[,21]<-ID.children
          
          pop.new <- f(sum(Nchild),pop.new,values.new,loci.new) #trait value function for the offspring
          
          pop<-rbind(pop,pop.new)
          rownames(pop) <- 1:nrow(pop)
          loci<-rbind(loci,loci.new)
          
        }#END ANY FEMALES?
      }#END IS OFFSPRING POSSIBLE?
        
        
      ##### DEATH START #####
      pop$survival[1:N]<-pop$survival[1:N]-1 #every adult loses one survival counter
      loci[pop$survival==0,1] <- -2  #if the survival is 0, it replaces the first loci with -2
       
      #random Death
      die_ind <- runif(nrow(pop),0,1) < DIE
      pop$survival[die_ind]<-0
      loci[die_ind,1] <- -2
        
      loci <- subset(loci,loci[,1]>(-2 )) #all rows with a -2 in the beginning are deleted
      pop <-subset(pop,pop$survival>0) #Individuals which have a survival higher then 0 stay alive in the dataframe
      ##### END DEATH #####
        
        
      ##### MIGRATION START #####
      if(nrow(pop)>0){
        wanderers<-runif(nrow(pop),0,1) < mig# draws one uniformmly distribued number for every individual
        pop$patch[wanderers] <- (pop$patch[wanderers] - 1 + floor(runif(sum(wanderers),1,patches)))%%patches + 1
          
        rownames(pop) <- 1:nrow(pop) #re-indexing the population
      }
      ##### MIGRATION END #####
        
      
      ###Statistic 2##
      Npop[t] <-nrow(pop) #overwrites the populationsizes for each generation in the empty vector (patch 1)
      Trait.pop[t] <- mean(pop$trait)#trait.N1.vector[t] <- mean(pop$trait[pop$patch==1]) #overwrites the average trait-value for each generation in the empty vector (patch 1)
        
      stats[,,t] <- stat.fun(pop,patches) #fills the arry with the statistic: N-pop, m-pop, w-pop, mean trait
      ##### End Statistic 2#############
        
      
    }#END IS ANYBODY THERE? 
    print(t)
  }##### END GENERATION LOOP #####
}##### END REPLICATION LOOP #####
#})

colours <- c("turquoise","violet","orange")
par(mfrow=c(1,2))

plot(Npop,main="Population over time", xlab="generations",ylab="population",type="l",col="black",ylim =c(0,max(Npop))) #plot traitvalue
for(ink in 1:patches){
  lines(stats[ink,1,],type="l",col=colours[ink])
}

plot(Trait.pop,main="mean trait value over time", xlab="generations",ylab="trait value",type="l",col="black",ylim = c(min(Trait.pop),max(Trait.pop)))

par(mfrow=c(2,2))
for(paint in 1:patches){
  plot(stats[paint,3,],main="Frequency of the Sexes Patch x", xlab="generations",ylab="frequency",type="l",col="red")
  lines(stats[paint,2,],type="l",col="green")
}

}#END MORPHEUS
Morpheus()
