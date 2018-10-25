#Hiwi-IBM

#requirements
# migration, multiple patches(the amount can be changed via 1 variable)
# male/female,Age, a traitvalue
rm(list=ls())

#Names will be changed, i´ve just used the first Names that came to mind
patches <-5
  
  for(i in 1:patches){
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
