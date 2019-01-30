##Simulation2.R

#Second Simulation, lack of pure independence

require(BnpSysRel)
require(mvtnorm)

##Calculate the truth
tires<-rmvnorm(100000, c(100, 110), 10*diag(.8, 2)+.2)
frontTire=tires[,1]
backTire = tires[,2]
brakes<-rmvnorm(100000, c(105, 105), 20*diag(.6, 2)+.4)
frontBrake = brakes[,1]
backBrake = brakes[,2]
gears<-rexp(100000, 1/5)+100

Tires<-apply(tires, 1, min)
Brakes<-apply(brakes, 1, max)
Bike<-apply(cbind(Tires, Brakes, gears), 1,min)


mean(frontBrake==Bike)
mean(backBrake==Bike)
mean(frontTire==Bike)
mean(backTire==Bike)
mean(gears==Bike)
