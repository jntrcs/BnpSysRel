###EstimateSystemReliability() Tests

priorList<-list(Generator =bsp(c(5, 8, 10), c(.2, .6, .8), precision = 2) ,
                Ignition = bsp(10, .3, 1))

dataList<-list(Generator = matrix(c(5.9, 1,
                                    6.8, 1,
                                    11, 0), byrow = T, nrow=3),
               Ignition=matrix(c(9, 1,
                                 15, 1), byrow = T, nrow=2),#,
               PowerSupply= matrix(c(30, 1,
                                     40, 0), byrow = T, nrow=2),
               Starter=matrix(c(12, 1), byrow = T, nrow=1))

a=estimateSystemReliability("StarterTestFile.txt", priorList, dataList)
a$Generator
a$Ignition
a$PowerSupply
