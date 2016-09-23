#Check if Rates are ordered correctly

NoNo.process.dat$BUSTED.omega1.MLE < NoNo.process.dat$BUSTED.omega2.MLE
NoNo.process.dat$BUSTED.omega2.MLE < NoNo.process.dat$BUSTED.omega3.MLE
NoNo.process.dat$BUSTED.SRV.omega1.MLE < NoNo.process.dat$BUSTED.omega2.MLE
NoNo.process.dat$BUSTED.SRV.omega2.MLE < NoNo.process.dat$BUSTED.SRV.omega3.MLE
NoNo.process.dat$true.alpha1 < NoNo.process.dat$true.alpha2
NoNo.process.dat$true.alpha2 < NoNo.process.dat$true.alpha3

NoYes.process.dat$BUSTED.omega1.MLE < NoYes.process.dat$BUSTED.omega2.MLE
NoYes.process.dat$BUSTED.omega2.MLE < NoYes.process.dat$BUSTED.omega3.MLE
NoYes.process.dat$BUSTED.SRV.omega1.MLE < NoYes.process.dat$BUSTED.omega2.MLE
NoYes.process.dat$BUSTED.SRV.omega2.MLE < NoYes.process.dat$BUSTED.SRV.omega3.MLE
NoYes.process.dat$true.alpha1 < NoYes.process.dat$true.alpha2
NoYes.process.dat$true.alpha2 < NoYes.process.dat$true.alpha3

YesYes.process.dat$BUSTED.omega1.MLE <= YesYes.process.dat$BUSTED.omega2.MLE
YesYes.process.dat$BUSTED.omega2.MLE <= YesYes.process.dat$BUSTED.omega3.MLE
YesYes.process.dat$BUSTED.SRV.omega1.MLE <= YesYes.process.dat$BUSTED.omega2.MLE
YesYes.process.dat$BUSTED.SRV.omega2.MLE <= YesYes.process.dat$BUSTED.SRV.omega3.MLE
YesYes.process.dat$true.alpha1 <= YesYes.process.dat$true.alpha2
YesYes.process.dat$true.alpha2 <= YesYes.process.dat$true.alpha3


YesNo.process.dat$BUSTED.omega1.MLE < YesNo.process.dat$BUSTED.omega2.MLE
YesNo.process.dat$BUSTED.omega2.MLE < YesNo.process.dat$BUSTED.omega3.MLE
YesNo.process.dat$BUSTED.SRV.omega1.MLE < YesNo.process.dat$BUSTED.omega2.MLE
YesNo.process.dat$BUSTED.SRV.omega2.MLE < YesNo.process.dat$BUSTED.SRV.omega3.MLE
YesNo.process.dat$true.alpha1 < YesNo.process.dat$true.alpha2
YesNo.process.dat$true.alpha2 < YesNo.process.dat$true.alpha3

# use number of rows of selectome data minus sum of the logic statement to check how many don't follow that rule
# cause it should just be 0 if everything is in order
nrow(Selectome.process.dat)-sum(Selectome.process.dat$BUSTED.omega1.MLE < Selectome.process.dat$BUSTED.omega2.MLE) #310
nrow(Selectome.process.dat)-sum(Selectome.process.dat$BUSTED.omega2.MLE < Selectome.process.dat$BUSTED.omega3.MLE) #0
nrow(Selectome.process.dat)-sum(Selectome.process.dat$BUSTED.SRV.omega1.MLE < Selectome.process.dat$BUSTED.omega2.MLE) #362
nrow(Selectome.process.dat)-sum(Selectome.process.dat$BUSTED.SRV.omega2.MLE < Selectome.process.dat$BUSTED.SRV.omega3.MLE)#0
nrow(Selectome.process.dat)-sum(Selectome.process.dat$true.alpha1 < Selectome.process.dat$true.alpha2) #337
nrow(Selectome.process.dat)-sum(Selectome.process.dat$true.alpha2 < Selectome.process.dat$true.alpha3) #665

ggplot(Selectome.melt.dat, aes(value))+geom_histogram()+facet_wrap(~variable, scales = "free_x")
Selectome.bias %>% melt() %>% ggplot(aes(value))+geom_histogram() +facet_wrap(~variable, scales = "free", ncol =2)
sum.table.omegas =gen.summary.o.tabs(Selectome.process.dat)[,3:4]
sum.table.alphas = gen.summary.a.tabs(Selectome.process.dat)[,3:4]

#Removing Numerical INFs
apply(Selectome.process.dat, 2, max)

no.outs.Selectome =Selectome.process.dat[-c(which(Selectome.process.dat$BUSTED.omega3.MLE >9999), 
                                            which(Selectome.process.dat$BUSTED.SRV.omega3.MLE>9999)),]

Selectome.bias = no.outs.Selectome[,c(1,17:18,21:22,4:5)]



nrow(no.outs.Selectome)-sum(no.outs.Selectome$BUSTED.omega1.MLE <= no.outs.Selectome$BUSTED.omega2.MLE) #227
nrow(no.outs.Selectome)-sum(no.outs.Selectome$BUSTED.omega2.MLE <= no.outs.Selectome$BUSTED.omega3.MLE) #0
nrow(no.outs.Selectome)-sum(no.outs.Selectome$BUSTED.SRV.omega1.MLE <= no.outs.Selectome$BUSTED.omega2.MLE) #302
nrow(no.outs.Selectome)-sum(no.outs.Selectome$BUSTED.SRV.omega2.MLE <= no.outs.Selectome$BUSTED.SRV.omega3.MLE)#0
nrow(no.outs.Selectome)-sum(no.outs.Selectome$true.alpha1 <= no.outs.Selectome$true.alpha2) #315
nrow(no.outs.Selectome)-sum(no.outs.Selectome$true.alpha2 <= no.outs.Selectome$true.alpha3) #648



