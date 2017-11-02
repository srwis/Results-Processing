install.packages(jsonlite)
library(jsonlite)
library(ggplot2)

  max.sim = fromJSON(txt = G:/BRC/SimResults/max) #read in the values that were set in the simulation file for the Max sim rates simulation
  max.dat = read.delim(file=G:/BRC/SimResults/maxRates.SRV.results.txt, header=TRUE,sep=\t,as.is = T);
  attach(max.dat)
  
  plot(max.dat$BUSTED.SRV.omega1.MLE)
  abline(h = 1)
  plot(max.dat$BUSTED.SRV.omega1.prop)
  abline(h=1)
 
  points(max.dat$BUSTED.SRV.omega2.MLE, col = red)
  abline(h=0.6666666666666666, col = red)
  
  plot(max.dat$BUSTED.SRV.omega3.MLE, col = blue)
  abline(h= 1)
  
  plot(max.dat$BUSTED.SRV.omega2.prop)
p = ggplot(max.dat, aes(FILE,BUSTED.SRV.omega1.MLE))


bothSel.dat = read.delim(file="G:/BRC/SimResults/bothSel.results.txt", header=TRUE,sep="\t",as.is = T)
plot(bothSel.dat$BUSTED.SRV.P, bothSel.dat$BUSTED.P)

summary(bothSel.dat$BUSTED.P)
summary(bothSel.dat$BUSTED.SRV.P)


test1.dat=compile("G:/BRC/SimResults/test1/")
max.dat = compile("G:/BRC/SimResults/maxRates/")
test2.dat = compile("G:/BRC/SimResults/test2/jsons/")
test3.dat = compile("G:/BRC/SimResults/test3/")
#graphs of P values for BUSTED vs SRV

plot(test1.dat$BUSTED.P,test1.dat$BUSTED.SRV.P)
abline(h=0.05,v=0.05, col = 'red')
ggplot(test1.process.dat, aes(x="true.alpha1",true.alpha1))+geom_boxplot()+facet_grid(p.cat~p.srv.cat)

plot(test2.dat$BUSTED.P,test2.dat$BUSTED.SRV.P, xlim = c(0,0.053))
abline(h=0.05,v=0.05, col = 'red')
test2.process.dat = process.dat(test2.dat)
ggplot(test2.process.dat, aes(x="true.alpha1",true.alpha1))+geom_boxplot()+facet_grid(p.cat~p.srv.cat)

plot(test3.dat$BUSTED.P,test3.dat$BUSTED.SRV.P,xlim = c(0,0.053))
abline(h=0.05,v=0.05, col = 'red')
test3.process.dat = process.dat(test3.dat)
ggplot(test3.process.dat, aes(x="true.alpha1",true.alpha1))+geom_boxplot()+facet_grid(p.cat~p.srv.cat)

plot(max.dat$BUSTED.P,max.dat$BUSTED.SRV.P)
max.process.dat =process.dat(max.dat)
abline(h=0.05,v=0.05, col = 'red')
ggplot(max.process.dat, aes(x="true.alpha1",true.alpha1))+geom_boxplot()+facet_grid(p.cat~p.srv.cat)
