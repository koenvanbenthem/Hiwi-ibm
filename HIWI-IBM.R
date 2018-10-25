#Hiwi-IBM

#requirements
# migration, multiple patches(the amount can be changed via 1 variable)
# male/female,Age, a traitvalue
rm(list=ls())

#Names will be changed, i´ve just used the first Names that came to mind


##### DEFINED PARAMETERS #####

P <-5 #number of patches
mig <- 0 #migrationrate
mut <- 0 #mutationrate
T <- 0 #traitvalue
d <- 0 #deathrate
off <- 0 #possible offspringnumber
cc <- 0 #carrying capacity
Pmax <- 0 #maximal size of patch
Pmin <- 0 #minimal size of patch


 ##### LOOP FOR PATCHES #####

  for(i in 1:P){ 
    Apollo <- round(rnorm(1, mean=50, sd=10)) #number male subjectsper patch
    Artemis <- round(rnorm(1, mean=50, sd=10)) #number female subjects per patch
    Leto <- rbind(matrix(rep("female",Artemis),ncol=1),matrix(rep("male",Apollo),ncol=1)) #matrix for the male and female subjects of the patch
    Nike <- cbind(round(runif(1,1,5)),round(runif(1,1,6))) #
    
    for(z in 2:(Apollo+Artemis)){
      Nike<-rbind(Nike,cbind(round(runif(1,1,5)),round(runif(1,1,6))))
    }
    
    Leto<-cbind(Leto,Nike)
    colnames(Leto)<-c("Gender","Age","Trait")
  }
