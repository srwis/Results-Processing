
source("./R_scripts/useful_functions.R")

library(readr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyr)

dat <- read_csv("C:/Users/srwisots/Google Drive/Muse/BRC/Results Processing/test_11_7_17_v2.csv")
View(dat)

#check how many omega rates
num.srv.omega.rates <- dat %>% select(contains("srv.omega.rate")) %>% ncol()
num.busted.omega.rate <- dat %>% select(contains("busted.omega.rate")) %>% ncol()

 dat %>% select(contains("omega.rate.3")) %>% summarise_each(funs(max(.,na.rm=TRUE)))


temp <- dat %>% select(BUSTED.SRV.LR,CV.SRV,BUSTED.SRV.LR,BUSTED.SRV.AICc,
                                Sites,Sequences, BUSTED.LR, BUSTED.P, BUSTED.SRV.P,BUSTED.AICc,
                       BUSTED.SRV.treelength, BUSTED.treelength)
pairs(temp)

gen.sig.table(dat) #this function still work as far as I can tell

idk <- add_p_cats(dat)

box.stats <- box.sum.stats(idk)

gen.boxplots(idk)

gen.p.plot(idk, "Images/p_plot_testset_11_7.png")

dat %>% ggplot(aes(x=BUSTED.P,y=BUSTED.SRV.P))+geom_bin2d()
dat %>% ggplot(aes(x=BUSTED.P,y=BUSTED.SRV.P))+geom_bin2d() + xlim (0,0.10) + ylim (0,0.10)




dat %>%  select(contains("omega.rate")) %>%  melt() %>% separate(variable, c("test","rate", "idk", "cat")) %>% 
  ggplot(aes(value)) + geom_histogram() + facet_grid(test ~ cat, scales = "free") + labs(title ="Omega rates")

dat %>%  select(contains("omega.rate")) %>% filter() %>%  melt() %>% separate(variable, c("test","rate", "idk", "cat")) %>% 
  ggplot(aes(value)) + geom_histogram() + facet_grid(test ~ cat, scales = "free") + labs(title ="Omega rates") 

dat %>% select(contains("omega.rate.1")) %>% melt() %>% separate(variable, c("test","rate", "idk", "cat")) %>% 
    ggplot(aes(value)) + geom_histogram() + facet_grid(test ~ cat) 

dat %>% select(contains("omega.rate.2")) %>% melt() %>% separate(variable, c("test","rate", "idk", "cat")) %>% 
  ggplot(aes(value)) + geom_histogram() + facet_grid(test ~ cat)

dat %>% select(contains("omega.rate.3")) %>% melt() %>% separate(variable, c("test","rate", "idk", "cat")) %>% 
  ggplot(aes(value)) + geom_histogram() + facet_grid(test ~ cat) + xlim(0,5)


dat %>% select(contains("omega.rate.3")) %>% summarise_each(funs(min, max))

plot_hist()

dat %>%  select(contains("omega")) %>%  melt() %>% separate(variable, c("test","rate", "idk", "cat")) %>% 
  ggplot(aes(value)) + geom_histogram(aes(fill = idk)) + facet_grid(test+idk ~ cat, scales = "free") +
  labs(title ="Omega rates")

agree.1 <- idk %>% filter((BUSTED.P >= 0.05 & BUSTED.SRV.P >=0.05) | (BUSTED.SRV.P <=0.05 & BUSTED.P <= 0.05))
agrre.2 <- idk %>% filter(BUSTED.SRV.P <=0.05, BUSTED.P <= 0.05)

agree <- union(agree.1, agrre.2)
not.agree <- setdiff(idk, agree)


dat %>% ggplot(aes(x= BUSTED.P, y = BUSTED.SRV.P)) + geom_point(aes(color = CV.SRV))

dat %>% select(BUSTED.P,BUSTED.SRV.P,CV.SRV) %>% pairs()

not.agree %>% select(BUSTED.P,BUSTED.SRV.P,CV.SRV) %>% pairs

dat %>% select(BUSTED.P,BUSTED.SRV.P,CV.SRV) %>% melt(id.vars = "CV.SRV") %>% ggplot(aes(y=CV.SRV, x = value))+
  geom_boxplot(aes(fill = variable))

dat %>% select(BUSTED.P,BUSTED.SRV.P,CV.SRV) %>% melt(id.vars = "CV.SRV") %>% ggplot(aes(y=CV.SRV, x = value))+
  geom_boxplot() + facet_grid(.~variable)

dat %>% select(FILE,BUSTED.P,BUSTED.SRV.P,CV.SRV) %>% melt(id.vars = c("CV.SRV","FILE")) %>% ggplot(aes(y=CV.SRV, x = value))+
  geom_point() + facet_grid(.~variable)
dat %>% select(FILE,BUSTED.P,BUSTED.SRV.P,CV.SRV) %>% melt(id.vars = c("CV.SRV","FILE")) %>% ggplot(aes(x=CV.SRV, y = value))+
  geom_point() + facet_grid(.~variable) + coord_cartesian(xlim = c(0,2))

dat<-dat %>% mutate(diff.P = BUSTED.SRV.P-BUSTED.P)

dat %>% ggplot(aes(x=CV.SRV, y = diff.P))+  geom_point() + coord_cartesian(xlim = c(0,2))

dat %>% filter(BUSTED.P < 0.001) %>% ggplot(aes(x=CV.SRV, y = diff.P))+  geom_point() +geom_smooth() + coord_cartesian(xlim = c(0,2)) 


dat %>% select(BUSTED.P,BUSTED.SRV.P,CV.SRV) %>% ggplot(aes(BUSTED.P,BUSTED.SRV.P)) + geom_bin2d(aes(color = CV.SRV))

dat %>% select(BUSTED.P,BUSTED.SRV.P,CV.SRV) %>% ggplot(aes(BUSTED.P,BUSTED.SRV.P, fill = CV.SRV)) + geom_raster()

dat %>% select(FILE,BUSTED.P,BUSTED.SRV.P,CV.SRV) %>% melt(id.vars = c("CV.SRV","FILE")) %>% ggplot(aes(y=CV.SRV, x = value))+
  geom_point(aes(color = variable))

gen.frac.sel.graphs(dat)


###Fraction under selection verses Cv
#want to see if the fraction under selection has any relation to CV values

dat = cbind(dat, srv.cv.cat = cut2(
  dat$CV.SRV , m = 100, levels.mean = TRUE
))

# data = cbind(data, omega2.cat = cut2(data$BUSTED.omega2.MLE  , m = 100, levels.mean = TRUE))
#then find the fraction of data sets under positive selection according to
#BUSTED+SRV p <= 0.05 in each bin

frac.cv.srv = with(dat, tapply(BUSTED.SRV.P, srv.cv.cat, function(x)
  length(which(x <= 0.05)) / length(x)))
frac.cv = with(dat, tapply(BUSTED.P, srv.cv.cat, function(x)
  length(which(x <= 0.05)) / length(x)))

plot(
  x = as.numeric(levels(dat$srv.cv.cat)),y = frac.cv.srv, ylim = c(0,1),
  pch = 19,
  xlab = "CV of SRV", ylab = "Fraction under selection"
)
points(
  x = as.numeric(levels(dat$srv.cv.cat)), y = frac.cv, col = 'red', pch = 19
)
lines(lowess(x = as.numeric(levels(
  dat$srv.cv.cat
)),y = frac.cv.srv, f = 0.25))
lines(lowess(x = as.numeric(levels(dat$srv.cv.cat)), y = frac.cv, f = 0.25))


#just curious can I make the same plot in ggplot2?

dat %>% ggplot(aes(x = as.numeric(levels(dat$srv.cv.cat)),y = frac.cv.srv))+ geom_point()
