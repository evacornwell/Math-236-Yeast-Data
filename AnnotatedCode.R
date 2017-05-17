## Final Project Annoted R Code
## Eva Cornwell and Mattie Johnson

library(manipulate)
library(deSolve)
library(RCurl)

######################################################

# The following code pulls our data directly from GitHub...
# There's no need to download the csv file manually
rawWebC <- getURL("https://raw.githubusercontent.com/evacornwell/Math-236-Yeast-Data/master/Yeast%20and%20Ethanol%20Data%20-%20Control.csv")
rawWeb1 <- getURL("https://raw.githubusercontent.com/evacornwell/Math-236-Yeast-Data/master/Yeast%20and%20Ethanol%20Data%20-%201per.csv")
rawWeb5 <- getURL("https://raw.githubusercontent.com/evacornwell/Math-236-Yeast-Data/master/Yeast%20and%20Ethanol%20Data%20-%205per.csv")

yeastcont <- read.csv(text=rawWebC)
yeast1 <- read.csv(text=rawWeb1)
yeast5 <- read.csv(text=rawWeb5)

####################################
# Results: Raw Data
####################################

time <- yeastcont$Hour

# Plot of yeast population data over time
plot(time,yeastcont$Cells,
     ylim=c(41220000,53580000),
     xlab="time (hours)", ylab="population (cells/mL)")
points(time,yeast1$Cells,
       col="blue",pch=1)
points(time,yeast5$Cells,
       col="red",pch=1)
legend("bottomright",c("CTRL","1% EtOH","5% EtOH"),
       pch=1,col=c("black","blue","red"),cex=0.8)

# Plot of unadjusted ethanol concentrations over time
plot(time,yeastcont$ABV, ylim=c(0,2), 
     xlab="Time (hours)", ylab="Ethanol Concentration (mL/mL)")
points(time,yeast1$ABV,col="blue")
points(time,yeast5$ABV,col="red")


####################################
# Hydrometer Model
####################################

expect <- c(0,1,5) #The ethanol concentrations we know to be true
obs<-c(yeastcont$ABV[1],yeast1$ABV[1],yeast5$ABV[1]) #Observed ethanol concentrations using hydrometer

# Fitting a quadratic function to the expected vs. observed Alcohol By Volume
plot(obs,expect,col="blue")
x=seq(0,1,by=0.01)
lines(x,6*x^2)

#Linear Model
fit<-lm(expect~obs+0)
plot(obs,expect,col="blue",
     xlab="Observed % ABV (mL/mL)",ylab="True % ABV (mL/mL)")
abline(fit)
coef(fit)[1]*yeast1$ABV
text(x=0.6,y=2,"y = 4.992x",cex=0.9)

# Ultimately, we chose use the linear model to adjust our abserved hydrometer readings
# Athough the quadratic model seems to fit the data very well, this model does not makes
# sense at high levels of observed % ABV. 

# Adding a column to each dataframes of the observed ABV adjusted using the linear model
yeastcont$trueABV=coef(fit)[1]*yeastcont$ABV
yeast1$trueABV=coef(fit)[1]*yeast1$ABV
yeast5$trueABV=coef(fit)[1]*yeast5$ABV


# Plot
plot(time,yeastcont$trueABV,ylim=c(0,10), xlab="Time (hours)",
     ylab="Concentration of ethanol (mL/mL)")
points(time,yeast1$trueABV,col="blue")
points(time,yeast5$trueABV,col="red")
legend("bottomright",c("CTRL","1% EtOH","5% EtOH"),pch=1,col=c("black","blue","red"),cex=0.8)



####################################
# Model
####################################

mod <- function(t,y,parms) {
  Y = y[1]; G=y[2];E = y[3];
  r=parms[1]
  gamma=parms[2]
  phi=parms[3]
  omega=parms[4]
  dY= r*Y*G-gamma*E*Y
  dG= -omega*r*G*Y
  dE= phi*G*Y
  return(list(c(dY,dG,dE)));
}


manipulate({
  y0=c(42540000,0.5681,0); times=seq(0,25,by=0.2)
  parms=c(r,gamma,phi,omega)
  out=ode(y0,times,mod,parms)
  par(mfrow=c(1,2))
  plot(out[,1],out[,2],type="l",xlab="Time (days)",ylab="Population",ylim=c(42540000,56540000),xlim=c(0,25))
  points(time,yeastcont$Cells,col="black")
  plot(out[,1],out[,4],lty=3,col="red",type="l",ylim=c(0,6))
  points(time,yeastcont$trueABV,col="red")
  lines(out[,1],out[,3],lty=1,col="blue")
  
},
r=slider(0,.5),gamma=slider(0,0.001),phi=slider(0,3.00e-8),omega=slider(0,.000000042)
)

## Fitting the Model
par(mfrow=c(1,1))

parms=c(0.1,0.001,.00000003,.000000042)

y0=c(42540000,0.5681,0); times=seq(0,25,by=0.2)
out=ode(y = y0, times, func = mod, parms=parms)

plot(out[,1],out[,2],type="l",xlab="Time (hours)",ylab="Population (cells/mL)",ylim=c(42540000,56540000),xlim=c(0,25))
points(time,yeastcont$Cells,col="black")
legend("bottomright", c("Obs","Model"),pch=c(1,NA),lty=c(NA,1),cex=0.8)


plot(out[,1],out[,4],lty=3,col="red",type="l",ylim=c(0,5),xlab="Time (hours)",
     ylab="Concentration (mL/mL or g/mL)")
points(time,yeastcont$trueABV,col="red")
lines(out[,1],out[,3],lty=1,col="blue")
legend("bottomright",c("Obs. [EtOH]","Model [EtOH]","[Sugar]"),
       col=c("red","red","blue"),pch=c(1,NA,NA),lty=c(NA,3,1),cex=0.8)




# Far into the future
y0=c(42540000,0.5681,0); times=seq(0,1000,by=0.2)
out=ode(y = y0, times, func = mod, parms=parms)
plot(out[,1],out[,2],type="l",xlab="Time (hours)",ylab="Population (cells/mL)",ylim=c(0,56540000),xlim=c(0,1000))

plot(out[,1],out[,4],lty=3,col="red",type="l",ylim=c(0,8),xlab="Time (hours)",
     ylab="Concentration (mL/mL and g/mL)")
lines(out[,1],out[,3],lty=1,col="blue")

####################################
# 1% Ethanol at start
####################################

manipulate({
  y0=c(yeast1$Cells[1],0.5681,1); times=seq(0,25,by=0.2)
  parms=c(r,gamma,phi,omega)
  out=ode(y0,times,mod,parms)
  par(mfrow=c(1,2))
  plot(out[,1],out[,2],type="l",xlab="Time (hours)",ylab="Population (cells/mL)",ylim=c(42540000,56540000),xlim=c(0,25))
  points(time,yeast1$Cells,col="black")
  plot(out[,1],out[,4],lty=3,col="red",type="l",ylim=c(0,6))
  points(time,yeast1$trueABV,col="red")
  lines(out[,1],out[,3],lty=1,col="blue")
  
},
r=slider(0,.5),gamma=slider(0,0.001),phi=slider(0,.00000003),omega=slider(0,.000000042))

## Fitting the Model
par(mfrow=c(1,1))

parms=c(.1,.001,.00000003,.000000042)

y0=c(yeast1$Cells[1],0.5681,1); times=seq(0,25,by=0.2)
out=ode(y = y0, times, func = mod, parms=parms)

plot(out[,1],out[,2],type="l",xlab="Time (hours)",ylab="Population (cells/mL)",ylim=c(yeast1$Cells[1],56540000),xlim=c(0,25))
points(time,yeast1$Cells,col="black")
legend("bottomright", c("Obs","Model"),pch=c(1,NA),lty=c(NA,1),cex=0.8)


plot(out[,1],out[,4],lty=3,col="red",type="l",ylim=c(0,6),xlab="Time (hours)",
     ylab="Concentration (mL/mL or g/mL)")
points(time,yeast1$trueABV,col="red")
lines(out[,1],out[,3],lty=1,col="blue")
legend("bottomright",c("Obs. [EtOH]","Model [EtOH]","[Sugar]"),
       col=c("red","red","blue"),pch=c(1,NA,NA),lty=c(NA,3,1),cex=0.8)




# Far into the future
y0=c(yeast1$Cells[1],0.5681,1); times=seq(0,700,by=0.2)
out=ode(y = y0, times, func = mod, parms=parms)
plot(out[,1],out[,2],type="l",xlab="Time (hours)",ylab="Population",ylim=c(0,56540000),xlim=c(0,700))
points(time,yeast1$Cells,col="black")
legend("topright", c("Obs","Model"),pch=c(1,NA),lty=c(NA,1),cex=0.8)

plot(out[,1],out[,4],lty=3,col="red",type="l",ylim=c(0,8),xlab="Time (hours)",
     ylab="% ABV (L/L)")
points(time,yeast1$trueABV,col="red")
lines(out[,1],out[,3],lty=1,col="blue")

####################################
# 5% Ethanol at start
####################################


manipulate({
  y0=c(yeast5$Cells[1],0.5681,5); times=seq(0,25,by=0.2)
  parms=c(r,gamma,phi,omega)
  out=ode(y0,times,mod,parms)
  par(mfrow=c(1,2))
  plot(out[,1],out[,2],type="l",xlab="Time (days)",ylab="Population",ylim=c(42540000,56540000),xlim=c(0,25))
  points(time,yeast5$Cells,col="black")
  plot(out[,1],out[,4],lty=3,col="red",type="l",ylim=c(0,10))
  points(time,yeast5$trueABV,col="red")
  lines(out[,1],out[,3],lty=1,col="blue")
  
},
r=slider(0,.5),gamma=slider(0,0.001),phi=slider(0,.00000003),omega=slider(0,.000000042))


## Fitting the Model
par(mfrow=c(1,1))

parms=c(0.1,0.001,.00000003,.000000042)

y0=c(yeast5$Cells[1],0.5681,5); times=seq(0,25,by=0.2)
out=ode(y = y0, times, func = mod, parms=parms)

plot(out[,1],out[,2],type="l",xlab="Time (hours)",ylab="Population (cells/mL)",ylim=c(yeast1$Cells[1],56540000),xlim=c(0,25))
points(time,yeast5$Cells,col="black")
legend("bottomright", c("Obs","Model"),pch=c(1,NA),lty=c(NA,1),cex=0.8)


plot(out[,1],out[,4],lty=3,col="red",type="l",ylim=c(0,10),xlab="Time (hours)",
     ylab="Concentration (mL/mL or g/mL)")
points(time,yeast5$trueABV,col="red")
lines(out[,1],out[,3],lty=1,col="blue")
legend("bottomright",c("Obs. [EtOH]","Model [EtOH]","[Sugar]"),
       col=c("red","red","blue"),pch=c(1,NA,NA),lty=c(NA,3,1),cex=0.8)

# Far into the future
y0=c(yeast5$Cells[1],0.5681,5); times=seq(0,700,by=0.2)
out=ode(y = y0, times, func = mod, parms=parms)
plot(out[,1],out[,2],type="l",xlab="Time (hours)",ylab="Population",ylim=c(0,56540000),xlim=c(0,700))
points(time,yeast5$Cells,col="black")
legend("topright", c("Obs","Model"),pch=c(1,NA),lty=c(NA,1),cex=0.8)

plot(out[,1],out[,4],lty=3,col="red",type="l",ylim=c(0,10),xlab="Time (hours)",
     ylab="% ABV (L/L)")
points(time,yeast5$trueABV,col="red")
lines(out[,1],out[,3],lty=1,col="blue")
