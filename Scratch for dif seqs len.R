###Scratch to compare the files with selection but different seq lengths

YesYes_10000_7 = read.csv("E:/BRC/SimResults/Five_seq/YesYes_results.csv")
YesYes_1000_3 = read.csv("E:/BRC/SimREsults/Five_seq/YesYes_1000_results.csv", as.is = T)
YesYes_6000_3 = read.csv("E:/BRC/SimResults/Five_seq/YesYes_6000_results.csv", as.is = T)
YesNo_10000_7 = read.csv("E:/BRC/SimResults/Five_seq/YesNo_results.csv", as.is = T)
YesNo_6000_3 = read.csv("E:/BRC/SimResults/Five_seq/YesNo_6000_results.csv", as.is = T)

library(dplyr)

YesYes_all.dat = bind_rows(YesYes_10000_7,YesYes_1000_3, YesYes_6000_3) 
YesNo_all.dat = bind_rows(YesNo_6000_3,YesNo_10000_7)

YesYes_all.dat %>% ggplot(aes(BUSTED.P,BUSTED.SRV.P))+geom_point(aes(shape = factor(Sites)))

YesNo_all.dat %>% ggplot(aes(BUSTED.P,BUSTED.SRV.P))+geom_point(aes(shape = factor(Sites)))

YesYes_all.dat %>% group_by(Sites) %>% select(ends_with(".P")) %>%  summarise_each(funs(mean))

YesNo_all.dat %>%  group_by(Sites) %>% select(ends_with(".P"))%>% summarise_each(funs(mean))



YesYes_all.dat %>% group_by(Sites) %>% mutate()




```{r summary stats table YesNo_1000 varying Omegas, echo = FALSE}
YesNo_all.dat =YesNo_all.dat %>% mutate(True.Rate.Omega.1 = 0.1, True.prop.Omega.1=0.5, True.Rate.Omega.2 = 0.5, True.prop.Omega.2 = 0.25, True.prop.Omega.3 = 0.25 )
```

```{r}
q = cbind("Truth"=as.numeric(sum.table[1:6,3]),variable=c("BUSTED.omega1.MLE","BUSTED.SRV.omega1.MLE","BUSTED.omega2.MLE","BUSTED.SRV.omega2.MLE","BUSTED.omega3.MLE", "BUSTED.SRV.omega3.MLE")) %>% as.data.frame()
q$Truth = as.numeric(sum.table[1:6,3])
YesNo_1000.bias %>% melt() %>% ggplot(aes(value))+geom_histogram() +facet_wrap(~variable, scales = "free", ncol = 2)+ geom_vline(data = q,aes(xintercept = as.numeric(Truth)),color = "red") + ggtitle("Distribution of Omega R")

true.vals = YesNo_all.dat %>% select(contains("omega")) %>% select(contains("Rate")) %>% glimpse()

YesNo_all.dat %>% group_by(True.Rate.Omega.3) %>% select(contains("omega")) %>% select(contains("3")) %>%  melt(id.vars = "True.Rate.Omega.3") %>% ggplot(aes(value))+geom_histogram(aes(fill=factor(True.Rate.Omega.3))) +facet_wrap(~variable, scales = "free", ncol = 2) + geom_vline(aes(xintercept =True.Rate.Omega.3),color = "red")


YesNo_all.dat %>% group_by(True.Rate.Omega.3) %>% ggplot(aes(BUSTED.SRV.omega3.MLE))+geom_histogram(aes(fill=factor(True.Rate.Omega.3))) +  geom_vline(aes(xintercept =True.Rate.Omega.3,color=factor(True.Rate.Omega.3)))+xlim(0.5,3.5)

```


YesYes_all.dat = bind_rows(YesYes_1000.dat,YesYes_1000_O3_1_1.dat,YesYes_1000_O3_1_2.dat,YesYes_1000_O3_1_3.dat,YesYes_1000_O3_1_4.dat,YesYes_1000_O3_1_5.dat,YesYes_1000_O3_1_6.dat,YesYes_1000_O3_1_7.dat,YesYes_1000_O3_1_8.dat,YesYes_1000_O3_1_9.dat)
YesYes_all.dat %>% ggplot(aes(BUSTED.P,BUSTED.SRV.P))+geom_point(aes(color = factor(True.Rate.Omega.3)))
YesYes_all.dat %>% group_by(True.Rate.Omega.3) %>% select(matches("*\\.P$|*omega3.MLE")) %>% summarise_each(funs(mean)) %>% kable("markdown")


YesYes_all.dat %>% group_by(True.Rate.Omega.3) %>% ggplot(aes(BUSTED.SRV.omega3.MLE))+geom_histogram(aes(fill=factor(True.Rate.Omega.3)), bins=50) + xlim(0,10) + geom_vline(aes(xintercept = True.Rate.Omega.3, color = factor(True.Rate.Omega.3)))

YesYes_all.dat %>% group_by(True.Rate.Omega.3) %>% ggplot(aes(BUSTED.SRV.omega3.MLE))+geom_density(aes(color=factor(True.Rate.Omega.3))) + xlim(0,10)