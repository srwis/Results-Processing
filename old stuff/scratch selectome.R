#scrat for selectome

Selectome.dat=read.csv("G:/Selectome/Selectome_results.csv",  as.is = T)
#Selectome.brlen = branch_length("G:/Selectome/jsons/jsons")
Selectome.process.dat = process.dat(Selectome.dat)
Selectome.melt.dat=melt(Selectome.process.dat[,c(1,31:35)])
Selectome.sum.stat = box.sum.stats(Selectome.process.dat)
Selectome.bias = Selectome.dat[,c(1,17:18,21:22,4:5)]

ggplot(Selectome.dat, )