#########################################################################
#      branching process
#########################################################################
#   model definitions   # parms3 = c(p1=0.5, p3=0.5, m0d=40, md=60, m0p=40, mp=8) 

branching.model<- function(bparms){

  require(distr)
  ngen <- 20                             # initialising default parameters
  p1 = bparms[1]; p3= bparms[2] ;       
  m0d=bparms[3]; md=bparms[4]; m0p=bparms[5]; mp=bparms[6]; 
  N <- rep(0,ngen+1)
  N[1] <-  50000          
  t_pts <- seq(0,120,12)              
  
  # distribution of time to division or proliferation
  Gdiv <- function(n, t){ 
    if (n==0) { 
      return( dexp(t, rate = 1/m0p))
    }
    else{
      return( dexp(t, rate = 1/mp))
    }
  }
  
  # distribution of time to death
  Ddie <- function(n, t){ 
    if (n==0) {
      return( dexp(t, rate = 1/m0d))
      
    }
    else{
      return( dexp(t, rate = 1/md))
    }
  }
  
  # define a function for the convolution of N folds of exponential distributions
  require(distr)
  nfoldconv <- function(n,x,r){
    dist.sum <- 0
    for (i in 1:n) {
      dist.sum <- dist.sum + Exp(rate = 1/r)
    }
    return(d(dist.sum)(x))
  }
  
  # convolution of the distributions of G0 and n folds of G 
  # n folds of G, x - times, r- division parameter for gen>0, 
  # r0- division parameter for zeroth generation
  
  
  convG0Gn <- function(n,x,r,r0){
    dist.sum <- 0
    for (i in 1:n) {
      dist.sum <- dist.sum + Exp(rate = 1/r)
    }
    dist.sum <- Exp(rate = 1/r0) + dist.sum
    return(d(dist.sum)(x))
  }
  
  # convolution of the distributions of G0, nfoldG and D(distribution of deaths)
  # n folds of G, x - times, r- division parameter for gen>0, r0- div par for zeroth
  # dr- parameter for exponential rate distribution
  
  convG0GnD <- function(n,x,r,r0,dr){
    dist.sum <- 0
    for (i in 1:n) {
      dist.sum <- dist.sum + Exp(rate = 1/r)
    }
    dist.sum <- Exp(rate = 1/r0) + dist.sum + Exp(rate = 1/dr)
    return(d(dist.sum)(x))
  }
  ####################################################################################
  ###        Cell populations is different generations
  ####################################################################################
  
  # number of cells of type 1 in all generations
  N1 <- function(n,t){
    if (n == 0) {
      return(p1*N[1]*(1-Ddie(0,t)) )
    }
    else if( n > 0){
      return(p3*N[1]*(2*(1-p3))^n*(convG0Gn(n-1,t,mp,m0p)-convG0GnD(n-1,t,mp,m0p,md)))
    }
    else{
      print("Invalid")
      break
    }
  }
  # number of cells of type 3 in each generation
  N3  <- function(n,t){
    if (n == 0) {
      return(p3*N[1]*(1-Gdiv(0,t)) )
    }
    else if( n > 0){
      return(p3*N[1]*(2*p3)^n*(convG0Gn(n-1,t,mp,m0p)-convG0Gn(n,t,mp,m0p)))
    }
    else{
      print("Invalid")
      break
    }
  }
  
  # Expectation of the number of death cells 
  N_neg1 <- function(n,t){
    if (n == 0) {
      return(p1*N[1]*Ddie(0,t) )
    }
    else if( n > 0){
      return(p3*N[1] * (2*(1-p3))^n * convG0GnD(n-1,t,mp,m0p,md))
    }
    else{
      print("Invalid")
      break
    }
  }
  matN <- matrix( rep(0, (length(t_pts)*(ngen+1))), nrow = length(t_pts), ncol = ngen+1)
  #number of cells in each division class
  matN[1,1] <- N[1]
  
  for (n in 1:21 ) {  
    for (t in 2:length(t_pts)) {
      matN[t,n]<- N1(n,t_pts[t]) + N3(n,t_pts[t]) + N_neg1(n,t_pts[t])
    }
  }
  # print(round(matN))
  # branchingdata <-ifelse(matN<0, 0, round(matN) )
  branchingdata <-ifelse(matN<0, 0, matN)
  return(branchingdata)
}


##################################################################
# Optimization of smith martin model
##################################################################
require(distr)
ngen <- 20                          # initialising default parameters
N <- rep(0,ngen+1)
N[1] <-  50000          
t_pts <- seq(0,120,12)

###############   parameter for cyton model #######################
parms3 = c(p1=0.5, p3=0.5, m0d=40, md=60, m0p=40, mp=8)  

# read the simulated data generated from agent based model
df <- read.table("../datasets/agent_based_sim_output_1.txt",header = TRUE)

sum.cols <- function(idx, data) return (rowSums(data[,idx]))

# reshape data to have total number of cells in each division class in a single column
dataset <- matrix( rep(0, (length(t_pts)*(ngen+1))), nrow = length(t_pts))
dataset[1,1] <- N[1]
dataset[-1,] <- sapply(split(2:43,rep(1:21,each=2)),sum.cols,df)



sinsum.cols <- function(idx,data) return(data[,idx])

# generates a candidate dataset using the parameter to simulate model
f2.sm <- function(parms3){
  out <- branching.model(parms3)
  return (sapply(split(1:21,rep(1:21,each=1)),sinsum.cols,out))
}

# fits the simulated data from agent based model to model simulation
f1.sm <- function(parms3, dataset){
  sim.dataset <- f2.sm(parms3)
  return(sum((sim.dataset - dataset)^2))
}
# optim function determines the best parameter canditate
out.optim <- optim(parms3, f1.sm, gr=NULL, dataset)
print(out.optim)







#######################################################################
# plots branching process
#######################################################################
model.sim <- branching.model(parms3)
agent.sim <- dataset
time <- t_pts
######################## plots ########################################

png(file="gen0.png",width=400,height=450,res = 120)
plot(c(0,120), c(0,50000), type = "n", xlab = "time (hrs)", ylab = "cell count",
     main = "Gen. 0",font.main=1)
# set up the axis
lines(time, model.sim[,1] , col="red", lwd=1.5)
lines(time, agent.sim[,1] , col="blue", lwd=1.5)
legend(40,49000, c('B-P','A-B'),
       lwd=c(1.5,1.5),col=c("red","blue"),lty=c(1,1))

dev.off()

png(file="gen1.png",width=400,height=450,res = 120)
plot(c(0,120), c(0,50000), type = "n", xlab = "time (hrs)", ylab = "cell count",
     main = "Gen. 1",font.main=1)
# set up the axis
lines(time, model.sim[,2] , col="red", lwd=1.5)
lines(time, agent.sim[,2] , col="blue", lwd=1.5)
legend(40,49000,  c('B-P','A-B'),
       lwd=c(1.5,1.5),col=c("red","blue"),lty=c(1,1))

dev.off()
png(file="gen2.png",width=400,height=450,res = 120)
plot(c(0,120), c(0,50000), type = "n", xlab = "time (hrs)", ylab = "cell count",
     main = "Gen. 2",font.main=1)
# set up the axis
lines(time, model.sim[,3] , col="red", lwd=1.5)
lines(time, agent.sim[,3] , col="blue", lwd=1.5)
legend(40,49000,  c('B-P','A-B'),
       lwd=c(1.5,1.5),col=c("red","blue"),lty=c(1,1))

dev.off()

png(file="gen3.png",width=400,height=450,res = 120)
plot(c(0,120), c(0,50000), type = "n", xlab = "time (hrs)", ylab = "cell count",
     main = "Gen. 3",font.main=1)
# set up the axis
lines(time, model.sim[,4] , col="red", lwd=1.5)
lines(time, agent.sim[,4] , col="blue", lwd=1.5)
legend(40,49000,  c('B-P','A-B'),
       lwd=c(1.5,1.5),col=c("red","blue"),lty=c(1,1))

dev.off()

png(file="gen4.png",width=400,height=450,res = 120)
plot(c(0,120), c(0,50000), type = "n", xlab = "time (hrs)", ylab = "cell count",
     main = "Gen. 4",font.main=1)
# set up the axis
lines(time, model.sim[,5] , col="red", lwd=1.5)
lines(time, agent.sim[,5] , col="blue", lwd=1.5)
legend(40,49000,  c('B-P','A-B'),
       lwd=c(1.5,1.5),col=c("red","blue"),lty=c(1,1))

dev.off()

png(file="gen5.png",width=400,height=450,res = 120)
plot(c(0,120), c(0,50000), type = "n", xlab = "time (hrs)", ylab = "cell count",
     main = "Gen. 5",font.main=1)
# set up the axis
lines(time, model.sim[,6] , col="red", lwd=1.5)
lines(time, agent.sim[,6] , col="blue", lwd=1.5)
legend(40,49000, c('B-P','A-B'),
       lwd=c(1.5,1.5),col=c("red","blue"),lty=c(1,1))

dev.off()
png(file="gen6.png",width=400,height=450,res = 120)
plot(c(0,120), c(0,50000), type = "n", xlab = "time (hrs)", ylab = "cell count",
     main = "Gen. 6",font.main=1)
# set up the axis
lines(time, model.sim[,7] , col="red", lwd=1.5)
lines(time, agent.sim[,7] , col="blue", lwd=1.5)
legend(40,49000, c('B-P','A-B'),
       lwd=c(1.5,1.5),col=c("red","blue"),lty=c(1,1))

dev.off()
png(file="gen7.png",width=400,height=450,res = 120)
plot(c(0,120), c(0,50000), type = "n", xlab = "time (hrs)", ylab = "cell count",
     main = "Gen. 7",font.main=1)
# set up the axis
lines(time, model.sim[,8] , col="red", lwd=1.5)
lines(time, agent.sim[,8] , col="blue", lwd=1.5)
legend(40,49000,  c('B-P','A-B'),
       lwd=c(1.5,1.5),col=c("red","blue"),lty=c(1,1))

dev.off()
png(file="gen8.png",width=400,height=450,res = 120)
plot(c(0,120), c(0,50000), type = "n", xlab = "time (hrs)", ylab = "cell count",
     main = "Gen. 8",font.main=1)
# set up the axis
lines(time, model.sim[,9] , col="red", lwd=1.5)
lines(time, agent.sim[,9] , col="blue", lwd=1.5)
legend(40,49000, c('B-P','A-B'),
       lwd=c(1.5,1.5),col=c("red","blue"),lty=c(1,1))

dev.off()
#########################################################################################












