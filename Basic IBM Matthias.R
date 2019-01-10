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
       Leron = {setwd("C:/Users/Leron/Desktop/IBM_code/")})

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

Here_is_your_ID<-function(Nchild,ID_scan){
  for(hebe in 1:sum(Nchild)){
  ID.children <- c(ID.children,ID_scan)
ID_scan<-ID_scan+1
  }
}

##### REPLICATION LOOP START#####
for(r in 1:replic){
  ##### INITIALISATION PATCHES #####
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
 
   ID_scan<-nrow(pop)+1
  ##### Statistic #####
  #population
  Npop <- rep(0,Nt) #empty vector for the populationsize of each generation in patch 1
  Trait.pop <- rep(0,Nt) #empty vector for the mean traitvalue of each generation
  
  Npop[1] <- nrow(pop) #populationsize for the first generation of patch 1
  Trait.pop[1] <- mean(pop$trait) #populationsize for the first generation of patch 2
  
  #patches
  #population,average traitvalue and genderdistribution per patch
  Npatch.1 <- rep(0,Nt) 
  Trait.patch.1 <- rep(0,Nt)
  Females.patch.1<- rep(0,Nt)
  Males.patch.1<- rep(0,Nt)
  
  Npatch.1[1] <-  nrow(subset(pop,pop$patch==1))
  Trait.patch.1[1] <- mean(subset(pop,pop$patch==1)$trait)
  Females.patch.1[1] <-nrow(subset(subset(pop,pop$patch==1),gender=="female"))
  Males.patch.1[1] <- nrow(subset(subset(pop,pop$patch==1),gender=="male"))
  

  Npatch.2 <- rep(0,Nt) 
  Trait.patch.2 <- rep(0,Nt)
  Females.patch.2<- rep(0,Nt)
  Males.patch.2<- rep(0,Nt)
  
  Npatch.2[1] <-  nrow(subset(pop,pop$patch==2))
  Trait.patch.2[1] <- mean(subset(pop,pop$patch==2)$trait)
  Females.patch.2[1] <-nrow(subset(subset(pop,pop$patch==2),gender=="female"))
  Males.patch.2[1] <- nrow(subset(subset(pop,pop$patch==2),gender=="male"))
  
  
  Npatch.3 <- rep(0,Nt) 
  Trait.patch.3 <- rep(0,Nt)
  Females.patch.3<- rep(0,Nt)
  Males.patch.3<- rep(0,Nt)
  
  Npatch.3[1] <-  nrow(subset(pop,pop$patch==3))
  Trait.patch.3[1] <- mean(subset(pop,pop$patch==3)$trait)
  Females.patch.3[1] <-nrow(subset(subset(pop,pop$patch==3),gender=="female"))
  Males.patch.3[1] <- nrow(subset(subset(pop,pop$patch==3),gender=="male"))
  
  ########Statistic End#####
  
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
      pop$trait[y] <- abs(sum(values[y,])) ##### USE OF COLUMN.NR
    }
    
  
  ##### GENERATION LOOP START #####  
  for(t in 1:Nt){
    
    N<- nrow(pop) #how many individuals there are in both patches
    
    ######IS ANYBODY THERE? START##############
    if(N>0){

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
    
    
    loci.new <- matrix(NA,nrow=sum(Nchild),ncol=21) #empty vector: children locis
    
    #### START LOOP PARTNERFINDING #####
    patchbook <- c()
    gendergram <- c()
    
    N.w.patch <- table(factor(N.w$patch,levels = 1:patches))# number of females in each patch (as a vector) (N.w.patch[k], number of females in patch k)
    N.m.patch <- table(factor(N.m$patch,levels = 1:patches))
    curr_child <- 1 # counter that keeps track of how many kids have emerged so far during the loop below
    if(nrow(N.w)>0){#ANY FEMALES START#####
      
    for(u in 1:nrow(N.w)){ #loop mother 
      if(Nchild[u]>0){ #just if the mother gets offspring
        mother<-N.w$ID[u] #gives the ID of the mother
        
        ###FATHER####
        if(N.m.patch[N.w$patch[u]]>0){#ANY MALES IN THE PATCH OF THE MOTHER? START
        # sample the ID of one male which patchnr. is the same as the patchnr. of the mother
        father<-sample(N.m$ID[N.m$patch==N.w$patch[u]],1)
        #GENETICS:
        loci.mother <- loci[loci[,21]==mother,] #vector of locis of the mother
        loci.father <- loci[loci[,21]==father,] #vector of locis of the father
        loci.child <- rep(0,ncol(loci)) #empty vector with fixed length
        
        if(Nchild[u]>0){ #just if the mother gets offspring
        for(o in 1:Nchild[u]){ #for loop for the number of children per female


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
        }#END LOOP NUMBER CHILDREN
        } 
        
        }#END ANY MALES?
        }#does the mother get offspring
    } #END LOOP PARTNERFINDING/mother
    
    patchbook <- rep(N.w$patch,Nchild) #each kid gets the patch of the mother
    ID.children <- Here_is_your_ID(Nchild,ID_scan)??????????????
    trait.children <- c(rep(0,length(patchbook))) 
    survival.children <- c(rep(max.Age,length(patchbook))) #each child gets the survival of the maximum age
    gender.children <- gendergram #gender of the children are written into the matrix
    patch.children <- patchbook #patches of children are written into the matrix
    pop.new <- data.frame(ID.children,patch.children,gender.children,trait.children,survival.children)
    colnames(pop.new)<-c("ID","patch","gender","trait","survival") #colum names
    
    # gen_phen_map[locus,bla,bla](loci.new)
    values.new <- matrix(NA,nrow=sum(Nchild),ncol=10) #empty matrix for the trait values for each loci

    for(d in 1:sum(Nchild)){ #for each individual offspring
        values.new[d,] <- gen_phen_map[cbind(1:10,loci.new[d,1:10],loci.new[d,11:20])]
    }
    pop.new$trait <- abs(rowSums(values.new)) ##### USE OF COLUMN.NR
    
    pop<-rbind(pop,pop.new)
    rownames(pop) <- 1:nrow(pop)
    loci<-rbind(loci,loci.new)
    }#END ANY FEMALES?
    }#END IS OFFSPRING POSSIBLE
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
    
    Npatch.1[t] <-  nrow(subset(pop,pop$patch==1))
    Trait.patch.1[t] <- mean(subset(pop,pop$patch==1)$trait)
    Females.patch.1[t] <-nrow(subset(subset(pop,pop$patch==1),gender=="female"))
    Males.patch.1[t] <- nrow(subset(subset(pop,pop$patch==1),gender=="male"))
    
    Npatch.2[t] <-  nrow(subset(pop,pop$patch==2))
    Trait.patch.2[t] <- mean(subset(pop,pop$patch==2)$trait)
    Females.patch.2[t] <-nrow(subset(subset(pop,pop$patch==2),gender=="female"))
    Males.patch.2[t] <- nrow(subset(subset(pop,pop$patch==2),gender=="male"))
    
    Npatch.3[t] <-  nrow(subset(pop,pop$patch==3))
    Trait.patch.3[t] <- mean(subset(pop,pop$patch==3)$trait)
    Females.patch.3[t] <-nrow(subset(subset(pop,pop$patch==3),gender=="female"))
    Males.patch.3[t] <- nrow(subset(subset(pop,pop$patch==3),gender=="male"))
    
    ##### End Statistic 2#############
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
# })

}#END MORPHEUS

##### PLOTS #####

plot(Npop,main="Population over time", xlab="generations",ylab="population",type="l",col="black",ylim =c(0,max(Npop))) #plot traitvalue
lines(Npatch.1,type="l",col="turquoise")
lines(Npatch.2,type="l",col="violet")
lines(Npatch.3,type="l",col="orange")



plot(Trait.pop,main="mean trait value over time", xlab="generations",ylab="trait value",type="l",col="black",ylim = c(0.001,0.02))
lines(Trait.patch.1,type="l",col="turquoise")
lines(Trait.patch.2,type="l",col="violet")
lines(Trait.patch.3,type="l",col="orange")



plot(Females.patch.1,main="Frequency of the Sexes Patch 1", xlab="generations",ylab="frequency",type="l",col="red")
lines(Males.patch.1,type="l",col="green")


plot(Females.patch.2,main="Frequency of the Sexes Patch 2", xlab="generations",ylab="frequency",type="l",col="red")
lines(Males.patch.2,type="l",col="green")

plot(Females.patch.3,main="Frequency of the Sexes Patch 3", xlab="generations",ylab="frequency",type="l",col="red")
lines(Males.patch.3,type="l",col="green")
