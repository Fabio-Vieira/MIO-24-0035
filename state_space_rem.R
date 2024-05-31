#In this script we will generate a network with a seasonal effect

#Installing packages
packages <- c("remify", "remstimate", "remstats", "devtools")
installed <- which(packages %in% rownames(installed.packages()))
if(sum(installed) == 0){
  sapply(packages, install.packages)
} else if (sum(installed) > 0){
  not_installed <- which(1:3 != installed)
  sapply(packages[not_installed], install.packages)
}

actors <- 1:10
D <- 10 #Number of weeks
intercept <- -3
beta_inertia <- c(.3, .8)
#Creating the riskset
riskset <- expand.grid(actors, actors)
riskset <- riskset[!(riskset[,1] == riskset[,2]),]
names(riskset) <- c("sender", "receiver")
#Creating the statistics matrix
outdegree <- matrix(NA, nrow = nrow(riskset), ncol = 0)
#Initializing outdegree of the sender vector
outdegree_vector <- rep(0, nrow(riskset))
#Creating event matrix
event <- matrix(NA, nrow = 0, ncol = 3)
#Counter to increament the day count
counter <- 1
time <- 0
for(i in 1:D){
  #We will create a day variable that when it reaches a certain threshold,
  #we will change the outdegree values, because we assume that it is weekend
  day <- 1
  while(day < 8){
    #Computing the event rate
    if(day < 6){ #weekdays
      if(sd(outdegree_vector) > 0){
        outdegree_vector_standadirze <- (outdegree_vector - mean(outdegree_vector))/sd(outdegree_vector)
        lambda <- exp(intercept + beta_inertia[1] * outdegree_vector_standadirze)
      } else {
        lambda <- exp(intercept + beta_inertia[1] * outdegree_vector)
      }
    } else { #weekends
      if(sd(outdegree_vector) > 0){
        outdegree_vector_standadirze <- (outdegree_vector - mean(outdegree_vector))/sd(outdegree_vector)
        lambda <- exp(intercept + beta_inertia[2] * outdegree_vector_standadirze)
      } else {
        lambda <- exp(intercept + beta_inertia[2] * outdegree_vector)
      }
    }
    time <- time + rexp(1, sum(lambda))
    #Sampling dyad
    dyad <- sample(1:nrow(riskset), size = 1, prob = lambda/sum(lambda))
    event <- rbind(event, cbind(time, riskset[dyad,]))
    if((time / (24*counter)) > 1){
      day <- day + 1
      counter <- counter + 1
    }
    #Updating inertia
    indicators <- which(riskset[,1] == riskset[dyad,1])
    outdegree_vector[indicators] <- outdegree_vector[indicators] + 1
    outdegree <- cbind(outdegree, outdegree_vector)
  }
  print(paste("This is week", i, "out of", D))
}
row.names(event) <- NULL
saveRDS(event, "edgelist.rds")

#######################################################################################################
#######################################################################################################
#######################################################################################################
#rm(list=ls())
if(!("remx" %in% installed.packages())){
  devtools::install_github("TilburgNetworkGroup/remx")
}
library(remx)
#Computing the statistics
library(remstats)
edgelist <- readRDS("edgelist.rds")
#Checking how many events per day we have
start <- as.numeric(10)
stop <- as.numeric(10)
start[1] <- 0
stop[1] <- 24
for(i in 2:70){
  start[i] <- stop[i-1]
  stop[i] <- stop[i-1]+24
}
stop[70] <- stop[70] + 1
#Getting the number of events per day
events <- sapply(1:70, function(x) nrow(edgelist[edgelist[,1] > start[x] & edgelist[,1] < stop[x],]))
barplot(events, ylim = c(0, max(events)+10))
box(which = "plot")

#######################################################################################################
#######################################################################################################
#######################################################################################################

#Creating weekly indicator
events_weekly <- cumsum(c(1, events))
events_weekly[71] <- sum(events)

#Declaring which effects we want remstats to compute
effects <- ~ remstats::outdegreeSender(scaling = "std")

#Getting the remify object for the entire sequence
rehObj <- remify::remify(edgelist, model = "tie")

data <- vector("list")

for(i in 2:length(events_weekly)){
  #Computing statistics
  stats <- remstats::tomstats(effects, rehObj, start = events_weekly[i-1], stop = events_weekly[i]-1)
  edl <- edgelist[events_weekly[i-1]:(events_weekly[i]-1),]
  #Every piece needs to be stored a in a list with edgelist and statistics
  data[[i-1]] <- list(edgelist = edl,
                      reh = remify::remify(edl, model = "tie", 
                                           actors = stream$actors), #pre processing
                      statistics = stats)
}

fit <- strem(data)
plot(fit)

#######################################################################################################
#######################################################################################################
#######################################################################################################




