#### Hiwi-IBM #####

#Names will be changed, i�ve just used the first Names that came to mind


rm(list=ls())

##### DEFINED PARAMETERS #####

P <-5 #number of patches
#mig <- 0 #migrationrate
#mut <- 0 #mutationrate
T <- 6 #traitvalue
#d <- 0 #deathrate
#off <- 0 #possible offspringnumber
#cc <- 0 #carrying capacity
#Pmax <- 0 #maximal size of patch
#Pmin <- 0 #minimal size of patch
age <- 5 #maximal age of subject


 ##### LOOP FOR PATCHES #####

  for(i in 1:P){ 
    N.female <- round(rnorm(1, mean=50, sd=10)) #number female subjects per patch
    N.male <- round(rnorm(1, mean=50, sd=10)) #number male subjects per patch
    N.total <- N.female + N.male #total number of subjects

    if(i==1){ #creates an array in the loop for the first patch
        patch.matrix <- array(-1,dim=c(N.total,3,P)) #empty Array including a matrix for each patch containing gender, age and trait value for each subject
    }

    population <- rbind(matrix(rep("female",N.female),ncol=1),matrix(rep("male",N.male),ncol=1)) #matrix for the male and female subjects (one column)
    
    ##### LOOP FOR AGE AND TRAIT #####

    for(z in 2:N.total){
      age.trait <-rbind(age.trait,cbind(round(runif(1,1,age)),round(runif(1,1,T)))) #creates uniformly distributed age and trait value for the other subjects connected to the privious subject
    }
    
   population <-cbind(population,age.trait)    
   patch.matrix[,,i] <- population #fills the first column (gender) of the array for the i patch    age.trait <- cbind(round(runif(1,1,age)),round(runif(1,1,T))) #creates uniformly distributed age and trait value for the first subject

    #colnames(population)<-c("Gender","Age","Trait")
  }
