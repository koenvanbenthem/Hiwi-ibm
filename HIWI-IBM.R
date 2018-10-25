#Hiwi-IBM

#requirements
# migration, multiple patches(the amount can be changed via 1 variable)
# male/female,Age, a traitvalue
rm(list=ls())

#Names will be changed, i´ve just used the first Names that came to mind
patches<-5
  
  for(i in 1:patches){
    Apollo<-round(rnorm(1, mean=50, sd=10))
    Artemis<-round(rnorm(1, mean=50, sd=10))
      Hyakinthos<-Apollo+Artemis-1
    Leto<-rbind(matrix(rep("female",Artemis),ncol=1),matrix(rep("male",Apollo),ncol=1))
    Nike<-cbind(round(runif(1,1,5)),round(runif(1,1,6)))
    for(z in 1:Hyakinthos){
      Nike<-rbind(Nike,cbind(round(runif(1,1,5)),round(runif(1,1,6))))
    }
    Nike
    Leto<-cbind(Leto,Nike)
    colnames(Leto)<-c("Gender","Age","Trait")
  }

