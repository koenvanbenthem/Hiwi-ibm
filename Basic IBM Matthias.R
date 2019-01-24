#Basic IBM
rm(list=ls())

Morpheus<-function(replicates=1, #number of replicates
                   time=1000, #number of generations
                   migrate=0.05, #migrationfactor
                   age=2, #age limit for an individual
                   patches=4, #number of Patches
                   mutate=0.05, #mutationfactor
                   die=0.05, #chance to die
                   
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
  

##### FUNCTIONS #####
w<-function(a,b,z,N,Np){ #FITNESS-FUNCTION
  y=a+b*plogis(c1+c2*N+c3*z+c4*(0.5*N-Np)+c5*N^2+c6*z^2+c7*(0.5*N-Np)^2+c8*z*N+c9*z*(0.5*N-Np)+c10*N*(0.5*N-Np))
  return(y)
}


Here_is_your_ID<-function(Nchild){ #ID-FUNCTION
  ID.children <-   ID_scan:(ID_scan+sum(Nchild)-1)
  ID_scan <<- ID_scan + sum(Nchild)
  return(ID.children)
}


f <- function(row,pop.matrix,value.matrix,loci.matrix){ #TRAIT-VALUE-FUNCTION
  value.matrix <- matrix(NA,nrow=row,ncol=10) #empty matrix for the trait values for each loci
  for(y in 1:row){ #for each individual
    for(z in 1:10){ 
      value.matrix[y,z] <- gen_phen_map[z,loci.matrix[y,z],loci.matrix[y,10+z]]
    }
    pop.matrix[y,4] <- abs(sum(value.matrix[y,]))
  }
  return(pop.matrix)
}


stat.fun <- function(pop.matrix, Npatch){ #PATCH/STATISTIC-FUNCTION
  tmp <- aggregate(pop.matrix$trait,by=list(patch = pop.matrix$patch),mean)
  traits <- tmp$x[match(1:Npatch,tmp$patch)] 
  
  cbind(table(factor(pop.matrix$patch,levels=1:Npatch)),
        table(factor(pop.matrix[pop.matrix$gender=='male',]$patch,levels=1:Npatch)),
        table(factor(pop.matrix[pop.matrix$gender=='female',]$patch,levels=1:Npatch)),
        as.numeric(traits))
}


##### REPLICATION LOOP START#####
for(r in 1:replicates){
  
  
  ##### INITIALISATION PATCHES #####
  pop <- c() #emptc vector for the population matrix
  stats <- array(NA,dim=c(patches,4,time)) #empty array for the statistics
  
  
  for(k in 1:patches){ #LOOP OVER PATCHES
    N_patchx <- abs(round(rnorm(1, mean=250, sd=10))) #Number of individuals in the patch 
    patchx_m <- round(runif(1,N_patchx/4,3*N_patchx/4)) #Number of males in the patch
    
    ID <- c(1:(N_patchx)) #vector ID: gives each individual an ID
    patch <- c(rep(k,N_patchx)) #vector patch: gives each individual their patch Nr.
    gender <- c(rep("male",patchx_m),rep("female",N_patchx-patchx_m)) #vector gender: is filled with males and females
    trait <- c(rep(0.5,N_patchx)) #vector trait: is for all individuals from both patches set as 0.5
    survival <- c(rep(age,N_patchx)) #vector survival: is for all new individuals of both patches the pre defined age limit 
    ID.mother <- c(rep(NA,N_patchx)) #the first generation has no mother and therefore no ID in the column for the mothers ID
    ID.father <- c(rep(NA,N_patchx)) #the first generation has no father and therefore no ID in the column for the fathers ID
    
    patchx <- data.frame(ID,patch,gender,trait,survival,ID.mother,ID.father) #the dataframe is constructed for each patch including all vectors which where defined just before
    pop <- rbind(pop,patchx)  #data frame including all individuals of all patches (the dataframe of a patch is included in the population matrix)
  }
  
  pop$ID <- c(1:nrow(pop)) #the first generation of the population becomes a new ID
  home <- c(1:patches) #vector of patchnumbers
 
  ID_scan <- nrow(pop)+1
  
  
  ##### STATISTIC START #####
  Npop <- rep(0,time) #empty vector for the populationsize of each generation 
  Trait.pop <- rep(0,time) #empty vector for the mean traitvalue of each generation
  
  Npop[1] <- nrow(pop) #the populationsize for the first generation is written into the vector
  Trait.pop[1] <- mean(pop$trait) #the mean traitvalue for the first generation is written into the vector
  ########STATISTIC END  #####
  
  
    population <- nrow(pop) #number of individuals
    loci <- matrix(NA,nrow=population,ncol=20+1) #empty matrix for the locis (20 numbers) and the ID of the individual (+1 number)
    
    for(x in 1:population){ #LOOP OVER THE INDIVIDUALS
      loci[x,] <- ceiling(runif(21,1e-16,10)) #each individual has 20 random numbers (first 10:row //last 10:column)
      loci[x,21] <- x #the last vector-spot is defined as x (the ID of the individual) for the first generation
    }
    
    pop <- f(population,pop,values,loci) #traitvalue-function: traitvalues for the population are included and overwrite the population matrix
    
  
    ##### GENERATION LOOP START #####  
    for(t in 1:time){
      N <- nrow(pop) #number of individuals in total (all patches included)

      
      if(N>0) { #START IS ANYBODY THERE-LOOP: if there are any individuals and the population is not extinct 
        N.x <- c() #empty vector for local populationsize
        N.w <- subset(pop,pop$gender=="female") #number of female individuals in total
        N.m <- subset(pop,pop$gender=="male") #number of male individuals in total
        chance <- c() #empty vector
        
        for(pls in 1:patches){ #START IS OFFSPRING POSSIBLE?
          chance <- c(chance,nlevels(subset(pop,pop$patch==pls)$gender)) #create a Vector which shows how many different arguments(levels) are in a vector 
        }
        
        if(max(chance)==2){ #if one patch contains both genders then it has a level of 2
          N.0 <- N/500
          
          for(j in 1:patches){ #loop over patches
            N.x <- c(N.x,nrow(subset(pop,pop$patch==j))/500) #vector of local population sizes
          }
          
          if(nrow(N.w)>0){ #number of offspring per femal
            Nchild <- 2*rpois(nrow(N.w),w(a,b,N.w$trait,N.0,N.x[N.w$patch])) #each female gets a random number of offspring based on the fitness-function
          }
        
          ID.children <- c() #empty vector for the ID of the offspring
          patch.children <- c() #empty vector for the patch of the offspring
          gender.children <- c() #empty vector for the gender of the offspring
          trait.children <- c() #empty vector for the trait of the offspring
          survival.children <- c() #each offspring gets the survival of the maximum age
          ID.mother.children <- c() #empty vector for the mothers ID of the offspring
          ID.father.children <- c() #empty vector for the fathers ID of the offspring
          
          loci.new <- matrix(NA,nrow=sum(Nchild),ncol=21) #empty vector for the locis of the offspring
          
          
          #### START LOOP PARTNERFINDING #####
          patchbook <- c() #empty vector for the patchnumber of the offspring
          gendergram <- c() #empty vector for the gender of the offspring
          
          N.w.patch <- table(factor(N.w$patch,levels = 1:patches)) #number of females in each patch (as a vector)
          N.m.patch <- table(factor(N.m$patch,levels = 1:patches)) #number of males in each patch (as a vector)
          curr_child <- 1 #counter that keeps track of how many kids have emerged so far during the loop below
          
          if(nrow(N.w)>0){ #START ANY FEMALES?: loop starts if there is at least one female individual
            for(u in 1:nrow(N.w)){ #START LOOP PARTNERFINDING/mother 
              if(Nchild[u]>0){ #START GETS THE MOTHER OFFSPRING?
                mother <- N.w$ID[u] #gives the ID of the mother
                ID.mother.children <- c(ID.mother.children, rep(mother,Nchild[u])) #ID of the mother is written into the vector for all her offspring
                
              
                ###FATHER####
                if(N.m.patch[N.w$patch[u]]>0){ #START ANY MALES IN THE PATCH OF THE MOTHER?: loop starts if there is at least one male individual in the mothers patch
                  father <- sample(N.m$ID[N.m$patch==N.w$patch[u]],1) #sample the ID of one male which patchnumber is the same as the patchnumber of the mother
                  ID.father.children <- c(ID.father.children,rep(father,Nchild[u])) #ID of the father is written into the vector as often as he becomes offspring with the mother
                  
                  #GENETICS:
                  loci.mother <- loci[loci[,21]==mother,] #vector of locis of the mother
                  loci.father <- loci[loci[,21]==father,] #vector of locis of the father
                  loci.child <- rep(0,ncol(loci)) #empty vector with fixed length for the locis of the offspring
                  
                  #if(Nchild[u]>0){ #START GETS THE MOTHER OFFSPRING?
                    for(o in 1:Nchild[u]){ #START LOOP NUMBER CHILDREN per female
                      loci.child[1:10] <- loci.mother[(1:10) +sample(c(0,10),10,replace=TRUE)] #the offspring becomes 10 locis sampled from the mother
                      loci.child[11:20] <- loci.father[(1:10) +sample(c(0,10),10,replace=TRUE)] #the offspring becomes 10 locis sampled from the father
                      
                      #MUTATION
                      if(runif(1,0,1) < mutate){ #if a random number is lower than the mutationrate the offspring becomes a random distributed loci
                        loci.child[round(runif(1,1,20))] <- round(runif(1,1,10))
                      }
                      
                      loci.new[curr_child,] <-  loci.child #connects loci of the offspring to the matrix of the other offspring in this generation
                      curr_child <- curr_child + 1
            
                      if(runif(1,0,1)>0.5){ #if random number is higher as 0.5, the offspring is female
                        gendergram <- c(gendergram,"female") #the gender is written in the gender vector for the offspring
                      } else{ #otherwise the offspring is male
                        gendergram <- c(gendergram,"male") #the gender is written in the gender vector for the offspring  
                      }
                      
                    } #END LOOP NUMBER CHILDREN
                  #} #END GETS THE MOTHER OFFSPRING?
              } #END ANY MALES IN THE PATCH OF THE MOTHER?
            } #END GETS THE MOTHER OFFSPRING?
          } #END LOOP PARTNERFINDING/mother
          
          patchbook <- rep(N.w$patch,Nchild) #each offspring becomes the patchnumber of the mother
          ID.children <- Here_is_your_ID(Nchild) #the ID of the offspring is calculated by the ID-function and written into the vector for their ID
          trait.children <- c(rep(0,length(patchbook))) #the traitvalue of the offspring is set to 0 for the moment
          survival.children <- c(rep(age,length(patchbook))) #each offspring gets the survival of the age limit pre defined
          gender.children <- gendergram #genders of the offspring are written into the matrix
          patch.children <- patchbook #patches of offspring are written into the matrix
          pop.new <- data.frame(ID.children,patch.children,gender.children,trait.children,survival.children,ID.mother.children,ID.father.children) #a new dataframe is made for the offspring of this generation
          colnames(pop.new) <- c("ID","patch","gender","trait","survival","ID.mother","ID.father") #column names of the dataframe
          loci.new[,21] <- ID.children #the ID of the offspring is written into the matrix of the locis of the offspring
          
          pop.new <- f(sum(Nchild),pop.new,values.new,loci.new) #the offspring matrix is overwritten including the traitvalues calculated by the traitvalue-function
          
          pop <- rbind(pop,pop.new) #the offspring population matrix is added to the general population matrix
          rownames(pop) <- 1:nrow(pop) #rownames are overwritten
          loci <- rbind(loci,loci.new) #the offspring loci matrix is added to the general loci matrix
          
        }#END ANY FEMALES?
      }#END IS OFFSPRING POSSIBLE?
        
        
      ##### DEATH START #####
      #death by age:
      pop$survival[1:N] <- pop$survival[1:N]-1 #every adult loses one survival counter per generation
      loci[pop$survival==0,1] <- -2  #if the survival is 0, it replaces the first loci with -2
       
      #random Death:
      die_ind <- runif(nrow(pop),0,1) < die #for each individual is a random number distributed. if the number is belo the deathrate the individual is written into a vector
      pop$survival[die_ind] <- 0 #the individuals that where written into the vactor below, become a 0 in their survival
      loci[die_ind,1] <- -2 #the individuals that where written into the vactor below, become a -2 in the first space of their loci row
        
      #erasing dead individuals:
      loci <- subset(loci,loci[,1]>(-2 )) #loci matrix: all rows with a -2 in the beginning are deleted
      pop <-subset(pop,pop$survival>0) #population matrix: Individuals which have a survival higher then 0 stay alive in the dataframe. the others are deleted
      ##### END DEATH #####
        
        
      ##### MIGRATION START #####
      if(nrow(pop)>0){ #loop starts if there is at least one individual in the population
        wanderers <- runif(nrow(pop),0,1) < migrate #draws one uniformly distribued number for every individual deciding if and where the individual migrates to
        pop$patch[wanderers] <- (pop$patch[wanderers] - 1 + floor(runif(sum(wanderers),1,patches)))%%patches + 1 #migration
          
        rownames(pop) <- 1:nrow(pop) #re-indexing the population
      }
      ##### MIGRATION END #####
        
      
      ###Statistic 2##
      Npop[t] <-nrow(pop) #overwrites the populationsizes for each generation in the empty vector
      Trait.pop[t] <- mean(pop$trait)#trait.N1.vector[t] <- mean(pop$trait[pop$patch==1]) #overwrites the average trait-value for each generation in the empty vector
        
      stats[,,t] <- stat.fun(pop,patches) #fills the arry with the statistic: N-pop, m-pop, w-pop, mean trait
      ##### End Statistic 2#############
        
      
    }#END IS ANYBODY THERE? 
    print(t)
  }##### END GENERATION LOOP #####
}##### END REPLICATION LOOP #####
#})


#PLOT 1 & 2
colours <- c("turquoise","violet","orange","blue","red4","seagreen4")
par(mfrow=c(1,2))

plot(Npop,main="Population over time", xlab="generations",ylab="population",type="l",col="black",ylim =c(0,max(Npop))) #plot traitvalue
for(ink in 1:patches){
  lines(stats[ink,1,],type="l",col=colours[ink])
}

plot(Trait.pop,main="mean trait value over time", xlab="generations",ylab="trait value",type="l",col="black",ylim = c(min(Trait.pop),max(Trait.pop)))


#PLOTS FOR EACH PATCH
par(mfrow=c(2,2))

for(paint in 1:patches){
  plot(stats[paint,3,],main="Frequency of the Sexes Patch x", xlab="generations",ylab="frequency",type="l",col="red")
  lines(stats[paint,2,],type="l",col="green")
}

}#END MORPHEUS
Morpheus()
