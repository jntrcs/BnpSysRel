###EstimateSystemReliability() Tests
require(BnpSysRel)
priorList<-list(Generator1 =bsp(c(5, 8, 10), c(.2, .6, .8), precision = 2) ,
                Ignition = bsp(10, .3, 1))
priorList$Generator2=priorList$Generator1

dataList<-list(Generator1 = matrix(c(5.9, 1,
                                    6.8, 1,
                                    11, 0), byrow = T, nrow=3),
               Generator2 = matrix(c(5.9, 1,
                                     6.8, 1,
                                     11, 0), byrow = T, nrow=3),
               Ignition=matrix(c(9, 1,
                                 15, 1), byrow = T, nrow=2),#,
               PowerSupply= matrix(c(30, 1,
                                     40, 0), byrow = T, nrow=2),
               Starter=matrix(c(12, 1), byrow = T, nrow=1))

a=estimateSystemReliability(file="Test/StarterTestFile.txt", priorList, dataList)
a$Generator1
a$Ignition
a$PowerSupply
plot(a$PowerSupply, withConfInt = T, withPaths = T)
