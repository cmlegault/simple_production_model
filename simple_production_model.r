# simple production model
# compare response surface to estimates of K and r from TMB
# assumes B1=K to allow easy response surface
# created 28 December 2015  Chris Legault
# last modified 31 December 2015

# run this section only once when first start program
first.run <- FALSE
if (first.run){
  library(TMB)
  library(Hmisc)
  
  # have to add RTools to path to allow cpp to be compiled
  # only run once per R session
  source("C:\\Users\\chris.legault\\Documents\\Working\\Rstuff\\TMB\\add_RTools_to_path.r")
  # now can compile
  
  my.dir <- "C:\\Users\\chris.legault\\Documents\\Working\\Rstuff\\TMB\\simple_production_model"
  setwd(my.dir)
}

# create a new set of data
create.new.data <- FALSE
if (create.new.data){
  true.K <- 20000
  true.r <- 0.4
  
  ny <- 36
  yr <- 1:ny
  # 4 blocks of F values randomly selected within range
  F.min <- c(0.01, 0.05, 0.30, 0.10) * true.r
  F.max <- c(0.05, 0.20, 0.60, 0.30) * true.r
  Fval <-  c( runif(ceiling(ny/4), F.min[1], F.max[1]),
              runif(ceiling(ny/4), F.min[2], F.max[2]),
              runif(ceiling(ny/4), F.min[3], F.max[3]),
              runif(ceiling(ny/4), F.min[4], F.max[4]))
  Fval <- round(Fval[1:ny],3)                  
  
  true.B <- rep(NA, ny)
  true.Y <- rep(NA, ny)
  true.B[1] <- true.K # assumes start at K in year 1
  for (i in 1:(ny-1)){
    true.Y[i] <- Fval[i] * true.B[i]
    true.B[i+1] <- true.B[i] + true.r * true.B[i] * (1 - true.B[i] / true.K) - true.Y[i] 
  }
  true.Y[ny] <- Fval[ny] * true.B[ny]
  
  B.sd <- 0.1
  obs.B <- true.B * exp(rnorm(ny) * B.sd)
  
  Y.sd <- 0.05
  obs.Y <- true.Y * exp(rnorm(ny) * Y.sd)

}
# or else read in data from rdat list
else{
  mydata <- dget("mydata.rdat")
  true.K <- mydata$true.K
  true.r <- mydata$true.r
  ny     <- mydata$ny
  Fval   <- mydata$Fval
  true.B <- mydata$true.B
  true.Y <- mydata$true.Y
  B.sd   <- mydata$B.sd
  Y.sd   <- mydata$Y.sd
  obs.B  <- mydata$obs.B
  obs.Y  <- mydata$obs.Y
  yr <- 1:ny
}

# save data
save.data <- FALSE
if(save.data){
  mydata <- list()
  mydata$true.K <- true.K
  mydata$true.r <- true.r
  mydata$ny     <- ny
  mydata$Fval   <- Fval
  mydata$true.B <- true.B
  mydata$true.Y <- true.Y
  mydata$B.sd   <- B.sd
  mydata$Y.sd   <- Y.sd
  mydata$obs.B  <- obs.B
  mydata$obs.Y  <- obs.Y
  dput(mydata,"mydata.rdat")
}

# plot the basic data
plot.basic <- FALSE
if (plot.basic){
  plot(yr,Fval)
  plot(yr,true.Y,ylim=c(0,max(c(true.Y,obs.Y))))
    points(yr,obs.Y,pch=16)
  plot(yr,true.B,ylim=c(0,max(c(true.B,obs.B))))
    points(yr,obs.B,pch=16)
}

#-------------------------------------------------------------------------------------------
# grid search of RSS for range of K and r values

# function to compute residuals between observed and predicted B given params and observed Y
calc.resid <- function(obs.B,K,r,Y){
  ny <- length(obs.B)
  B <- rep(NA, ny)
  B[1] <- K
  for(ii in 1:(ny-1)){
    B[ii+1] <- B[ii] + r * B[ii] * (1 - B[ii] / K) - Y[ii]  
  }
  if(min(B) < 0) resid <- 1e6
  else{
    resids <- (log(B) - log(obs.B))^2
    resid <- sum(resids)
  }
  return(resid)           
}
  
rf.K <- 0.5   # range factor for K
rf.r <- 0.9   # range factor for r
nbin <- 101   # number of bins in each direction for grid search
est.K <- seq((1-rf.K) * true.K, (1+rf.K) * true.K, length.out=nbin)
est.r <- seq((1-rf.r) * true.r, (1+rf.r) * true.r, length.out=nbin)

res <- matrix(NA, nrow=nbin, ncol=nbin)  # container for residuals

minres <- 1e6
for(i in 1:nbin){
  K <- est.K[i]
  for(j in 1:nbin){
    r <- est.r[j]
    resid <- calc.resid(obs.B,K,r,obs.Y)
    res[i,j] <- resid
    if (!is.na(resid) & resid < minres){
      minres <- resid
      minres.i <- K
      minres.j <- r
    }
  }
}
# plot the result of the grid search
contour(x=est.K,y=est.r,z=res,levels=c(2,4,6,8,10,20,50,100),col=c("red",rep("black",7)),xlab="K",ylab="r")
  points(minres.i,minres.j,pch=16,col="red")
  points(true.K,true.r,pch=15,col="blue")

simp.K <- minres.i
simp.r <- minres.j
simp.B <- rep(NA, ny)
simp.B[1] <- simp.K
for(ii in 1:(ny-1)){
  simp.B[ii+1] <- simp.B[ii] + simp.r * simp.B[ii] * (1 - simp.B[ii] / simp.K) - obs.Y[ii]  
}
simp.F <- obs.Y / simp.B

#-------------------------------------------------------------------------------------------
# now estimate K and r using TMB
# slightly different approach to avoid NaN - estimate F each year as a random walk random effect

# compile cpp code and load dll
compile.and.load <- FALSE
if (compile.and.load){
  compile("simple_production_model.cpp")
  dyn.load(dynlib("simple_production_model"))
}

# set up data and parameters
dat <- list(
  obs_B=obs.B,
  obs_Y=obs.Y
)

parameters <- list(
  logTMB_B=rep(0,ny),
  logitTMB_F=rep(0,(ny-1)),
  logTMB_K=log(true.K),
  logTMB_r=log(true.r),
  logsigmaF=0,
  logsigmaB1=log(0.01),
  logsigmaB2=log(B.sd),
  logsigmaY=log(Y.sd)
)

Map = list()
Map[["logsigmaF"]] = factor(NA)
Map[["logsigmaB1"]] = factor(NA)
Map[["logsigmaB2"]] = factor(NA)
Map[["logsigmaY"]] = factor(NA)

# now estimate surplus production model using TMB
obj <- MakeADFun(dat,parameters,DLL="simple_production_model", random=c("logTMB_B","logitTMB_F"), map=Map, silent=TRUE)#FALSE)
#obj <- MakeADFun(dat,parameters,DLL="simple_production_model", map=Map, silent=TRUE)#FALSE)

opt <- nlminb(obj$par, obj$fn, obj$gr, control=list(trace=1,iter.max=1000,eval.max=1000))
obj$report()
rep <- sdreport(obj)
srep <- summary(sdreport(obj))
srep
TMB_B <- srep[rownames(srep)=="TMB_B",]
TMB_F <- srep[rownames(srep)=="TMB_F",]

make.compare.plots <- FALSE
if(make.compare.plots){
  plot(1:ny,true.B,ylim=c(0,max(c(true.B,simp.B,TMB_B[,1]+2*TMB_B[,2]))))
    points(1:ny,simp.B,pch=4,col="red")
    errbar(1:ny,TMB_B[,1],TMB_B[,1]+2*TMB_B[,2],TMB_B[,1]-2*TMB_B[,2],pch=16,col="blue",add=T)

  plot(1:ny,Fval,ylim=c(0,max(c(Fval,simp.F,TMB_F[,1]+2*TMB_F[,2]))))
    points(1:ny,simp.F,pch=4,col="red")
    errbar(1:(ny-1),TMB_F[,1],TMB_F[,1]+2*TMB_F[,2],TMB_F[,1]-2*TMB_F[,2],pch=16,col="blue",add=T)
}
