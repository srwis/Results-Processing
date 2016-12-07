#scratch for tables

library(jsonlite)
library(ggplot2)
library(knitr)
library(stringr)
library(reshape2)
library(dplyr)
library(xtable)

compile <- function(cur.dir,csv){
  jsons <- list.files(path = cur.dir,
                      pattern = '*.json', recursive = TRUE)
  
  #create empty data.frame with variable names
  names = c("FILE", "BUSTED.LR", "BUSTED.SRV.LR", "BUSTED.omega3.MLE", "BUSTED.SRV.omega3.MLE", "BUSTED.omega3.prop",
            "BUSTED.SRV.omega3.prop", 'CV.SRV', 'BUSTED.P', 'BUSTED.SRV.P','BUSTED.AICc','BUSTED.SRV.AICc',
            'BUSTED.treelength' ,'BUSTED.SRV.treelength', 'Sites', 'Sequences', 
            'BUSTED.omega1.MLE','BUSTED.SRV.omega1.MLE', 'BUSTED.omega1.prop','BUSTED.SRV.omega1.prop',
            'BUSTED.omega2.MLE','BUSTED.SRV.omega2.MLE', 'BUSTED.omega2.prop','BUSTED.SRV.omega2.prop', 'SRV.alpha3.MLE',
            'SRV.alpha3.prop','SRV.alpha1.MLE','SRV.alpha1.prop','SRV.alpha2.MLE','SRV.alpha2.prop')
  
  classes = c("character", "numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric",
              "numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric",
              "numeric","numeric","numeric")
  df =read.table(text="", col.names = names, colClasses = classes)
  for (i in  seq(from=1, to=length(jsons), by=2)){
    filepath = paste(cur.dir,jsons[i], sep="")
    
    test = filepath %>% readLines() %>% gsub(x=.,pattern="nan",replacement ='"NA"') %>% fromJSON()
    
    
    
    FILE = jsons[i]
    Sites = length(test$profiles$unconstrained)
    tree_string = test$fits$`Unconstrained model`$`tree string`
    x= tree_string %>% str_replace_all("\\(",":") %>% str_replace_all("\\)",":") %>%     str_replace_all(",",":") %>% str_split(":")
    x= unlist(x)
    x =x[x !=""]
    br_len = matrix(x,ncol = 2,  byrow = T)
    colnames(br_len) = c("Branch", "Length")
    
    Sequences = sum(grepl("Node*", br_len[,1]) == FALSE)
    
    if (grepl("BUSTED-SRV",jsons[i])){
      filepath = paste(cur.dir,jsons[i], sep="")
      
      test = filepath %>% readLines() %>% gsub(x=.,pattern="nan",replacement ='"NA"') %>% fromJSON()      
      
      
      BUSTED.SRV.P = test$`test results`$p
      BUSTED.SRV.LR =test$`test results`$LR
      BUSTED.SRV.AICc = test$fits$`Unconstrained model`$`AIC-c`
      BUSTED.SRV.treelength = test$fits$`Unconstrained model`$`tree length`
      
      #OMEGA values for BUSTED.SRV
      BUSTED.SRV.omega3.MLE = test$fits$`Unconstrained model`$`rate distributions`$FG[3,1]
      BUSTED.SRV.omega3.prop = test$fits$`Unconstrained model`$`rate distributions`$FG[3,2]
      BUSTED.SRV.omega2.MLE = test$fits$`Unconstrained model`$`rate distributions`$FG[2,1]
      BUSTED.SRV.omega2.prop = test$fits$`Unconstrained model`$`rate distributions`$FG[2,2]
      BUSTED.SRV.omega1.MLE = test$fits$`Unconstrained model`$`rate distributions`$FG[1,1]
      BUSTED.SRV.omega1.prop = test$fits$`Unconstrained model`$`rate distributions`$FG[1,2]
      #ALPHA values for BUSTED.SRV
      SRV.alpha3.MLE = test$fits$`Unconstrained model`$`rate distributions`$SRV[3,1]
      SRV.alpha3.prop = test$fits$`Unconstrained model`$`rate distributions`$SRV[3,2]
      SRV.alpha2.MLE = test$fits$`Unconstrained model`$`rate distributions`$SRV[2,1]
      SRV.alpha2.prop = test$fits$`Unconstrained model`$`rate distributions`$SRV[2,2]
      SRV.alpha1.MLE = test$fits$`Unconstrained model`$`rate distributions`$SRV[1,1]
      SRV.alpha1.prop = test$fits$`Unconstrained model`$`rate distributions`$SRV[1,2]
      
      mom2 = SRV.alpha3.MLE^2*SRV.alpha3.prop+ SRV.alpha1.MLE^2*SRV.alpha1.prop+ SRV.alpha2.MLE^2*SRV.alpha2.prop
      mean= SRV.alpha3.MLE*SRV.alpha3.prop+ SRV.alpha1.MLE*SRV.alpha1.prop+ SRV.alpha2.MLE*SRV.alpha2.prop
      CV.SRV = sqrt(mom2-mean^2)/mean
      
    }
    if (grepl("BUSTED-SRV",jsons[i+1])==FALSE){
      filepath = paste(cur.dir,jsons[i+1], sep="")
      
      test = filepath %>% readLines() %>% gsub(x=.,pattern="nan",replacement ='"NA"') %>% fromJSON()    
      BUSTED.P = test$`test results`$p
      BUSTED.LR = test$`test results`$LR
      BUSTED.AICc = test$fits$`Unconstrained model`$`AIC-c`
      BUSTED.treelength = test$fits$`Unconstrained model`$`tree length`
      
      #OMEGA values for BUSTED
      BUSTED.omega3.MLE = test$fits$`Unconstrained model`$`rate distributions`$FG[3,1]
      BUSTED.omega3.prop = test$fits$`Unconstrained model`$`rate distributions`$FG[3,2]
      BUSTED.omega2.MLE = test$fits$`Unconstrained model`$`rate distributions`$FG[2,1]
      BUSTED.omega2.prop = test$fits$`Unconstrained model`$`rate distributions`$FG[2,2]
      BUSTED.omega1.MLE = test$fits$`Unconstrained model`$`rate distributions`$FG[1,1]
      BUSTED.omega1.prop = test$fits$`Unconstrained model`$`rate distributions`$FG[1,2]
      #ALPHA values for BUSTED
      #       BUSTED.alpha3.MLE = test$fits$`Unconstrained model`$`rate distributions`$SRV[3,1]
      #       BUSTED.alpha3.prop = test$fits$`Unconstrained model`$`rate distributions`$SRV[3,2]
      #       BUSTED.alpha2.MLE = test$fits$`Unconstrained model`$`rate distributions`$SRV[2,1]
      #       BUSTED.alpha2.prop = test$fits$`Unconstrained model`$`rate distributions`$SRV[2,2]
      #       BUSTED.alpha1.MLE = test$fits$`Unconstrained model`$`rate distributions`$SRV[1,1]
      #       BUSTED.alpha1.prop = test$fits$`Unconstrained model`$`rate distributions`$SRV[1,2]
      
    }
    
    df[nrow(df)+1,] <- c(FILE, BUSTED.LR, BUSTED.SRV.LR, BUSTED.omega3.MLE, BUSTED.SRV.omega3.MLE, BUSTED.omega3.prop,
                         BUSTED.SRV.omega3.prop, CV.SRV, BUSTED.P, BUSTED.SRV.P,BUSTED.AICc,BUSTED.SRV.AICc,
                         BUSTED.treelength ,BUSTED.SRV.treelength, Sites, Sequences, 
                         BUSTED.omega1.MLE,BUSTED.SRV.omega1.MLE, BUSTED.omega1.prop,BUSTED.SRV.omega1.prop,
                         BUSTED.omega2.MLE,BUSTED.SRV.omega2.MLE, BUSTED.omega2.prop,BUSTED.SRV.omega2.prop, SRV.alpha3.MLE,
                         SRV.alpha3.prop,SRV.alpha1.MLE,SRV.alpha1.prop,SRV.alpha2.MLE,SRV.alpha2.prop)
    
  }
  df[,2:30]=as.numeric(unlist(df[,2:30]))
  write.csv(file = csv, x = df, row.names= F)
  # return(df)
}

simulation_inputs <- function(dir,csv){
  require("stringr")
  require("jsonlite")
  require("dplyr")
  list = list.files(path = dir, recursive = T, pattern ="^([^.]+)$")
  #set up the empty data frame
  names = c("Sites" ,"Cat", "Omega 1 Rate", "Omega 2 Rate", "Omega 3 Rate", "Omega 1 Weight", "Omega 2 Weight", "Omega 3 Weight",
            "Alpha 1 Rate", "Alpha 2 Rate", "Alpha 3 Rate","Alpha 1 Weight", "Alpha 2 Weight", "Alpha 3 Weight")
  setup.tab = read.table(text = "",col.names = names)
  #loop thru each file to get info in correct format
  for(i in seq(from = 1, to= length(list))){
    x=readLines(paste(dir,list[i], sep = "/"))
    #making this a readable json
    x1 = str_replace(x,":\\{", ":\\[")
    
    x1[c(5:6,12:13)] = str_replace(x1[c(5:6,12:13)],"\\{([-+]?[0-9]*\\.?[0-9]+,[-+]?[0-9]*\\.?[0-9]+)\\}", "\\[\\1\\]," )
    
    x1[c(7,14)] = str_replace(x1[c(7,14)],"\\{([-+]?[0-9]*\\.?[0-9]+,[-+]?[0-9]*\\.?[0-9]+)\\}", "\\[\\1\\]" )
    x1[c(8,15)] = str_replace(x1[c(8,15)], "\\}", "\\]")
    r= fromJSON(x1)
    omega_rates = r$`omega distribution`[,1]
    names(omega_rates)= c("Omega 1 Rate", "Omega 2 Rate", "Omega 3 Rate")
    omega_weights = r$`omega distribution`[,2]
    names(omega_weights) = c("Omega 1 Weight", "Omega 2 Weight", "Omega 3 Weight")
    
    Alpha_rates = r$`alpha distribution`[,1]
    names(Alpha_rates)= c("Alpha 1 Rate", "Alpha 2 Rate", "Alpha 3 Rate")
    Alpha_weights = r$`alpha distribution`[,2]
    names(Alpha_weights) = c("Alpha 1 Weight", "Alpha 2 Weight", "Alpha 3 Weight")
    
    
    setup.tab[nrow(setup.tab)+1,] = c(r$sites,list[i], omega_rates,omega_weights,Alpha_rates,Alpha_weights)
  }
  setup.tab =select(setup.tab, Sites, Cat, Omega.1.Rate, Omega.1.Weight, Omega.2.Rate, Omega.2.Weight,Omega.3.Rate, Omega.3.Weight,
                    Alpha.1.Rate, Alpha.1.Weight, Alpha.2.Rate, Alpha.2.Weight,Alpha.3.Rate, Alpha.3.Weight)
  setup.tab[,3:14] = as.numeric(unlist(setup.tab[,3:14]))
  mom2 = setup.tab$Alpha.1.Rate^2*setup.tab$Alpha.1.Weight+setup.tab$Alpha.2.Rate^2*setup.tab$Alpha.2.Weight + setup.tab$Alpha.3.Rate^2*setup.tab$Alpha.3.Weight
  
  mean = setup.tab$Alpha.1.Rate*setup.tab$Alpha.1.Weight+setup.tab$Alpha.2.Rate*setup.tab$Alpha.2.Weight + setup.tab$Alpha.3.Rate*setup.tab$Alpha.3.Weight
  
  setup.tab = mutate(setup.tab,CV.SRV = sqrt(mom2-mean^2)/mean)
  
  #return(setup.tab)
  write.csv(file = csv, x = setup.tab, row.names= F)
}


add_truth <- function(dat, truth){
  dat = dat %>% mutate(True.CV = 1337)
  dat = dat %>% mutate(True.omega3.rate = 8008)
  dat = dat %>% mutate(True.omega1.rate = 8008)
  dat = dat %>% mutate(True.omega2.rate = 8008)
  dat = dat %>% mutate(True.alpha3.rate = 8008)
  dat = dat %>% mutate(True.alpha2.rate = 8008)
  dat = dat %>% mutate(True.alpha1.rate = 8008)
  dat = dat %>% mutate(True.omega3.prop = 8008)
  dat = dat %>% mutate(True.omega1.prop = 8008)
  dat = dat %>% mutate(True.omega2.prop = 8008)
  dat = dat %>% mutate(True.alpha3.prop = 8008)
  dat = dat %>% mutate(True.alpha2.prop = 8008)
  dat = dat %>% mutate(True.alpha1.prop = 8008)
  for( i in seq(from = 1, to = nrow(truth), by = 1)){
    temp = which(str_detect(dat$FILE,truth$Cat[i]))
    dat$True.CV[temp] = truth$CV.SRV[i] 
    dat$True.omega3.rate[temp] = truth$Omega.3.Rate[i]
    dat$True.omega2.rate[temp] = truth$Omega.2.Rate[i]
    dat$True.omega1.rate[temp] = truth$Omega.1.Rate[i]
    dat$True.alpha1.rate[temp] = truth$Alpha.1.Rate[i]
    dat$True.alpha2.rate[temp] = truth$Alpha.2.Rate[i]
    dat$True.alpha3.rate[temp] = truth$Alpha.3.Rate[i]
    dat$True.omega3.prop[temp] = truth$Omega.3.Weight[i]
    dat$True.omega2.prop[temp] = truth$Omega.2.Weight[i]
    dat$True.omega1.prop[temp] = truth$Omega.1.Weight[i]
    dat$True.alpha1.prop[temp] = truth$Alpha.1.Weight[i]
    dat$True.alpha2.prop[temp] = truth$Alpha.2.Weight[i]
    dat$True.alpha3.prop[temp] = truth$Alpha.3.Weight[i]
  }
  return(dat)
}


dir= "G:/BRC/SimResults/Five_seq/"
basename = "Five_seq_all"
process_dat <- function(dir, basename){
  temp = paste(dir,basename,"_Truth.csv", sep = "")
  simulation_inputs(dir,temp)
  truth = read.csv(temp, as.is = T)
  temp = paste(dir,basename,"_results.csv", sep = "")
  compile(dir,temp)
  dat = read.csv(temp, as.is = T)
  dat = add_truth(dat, truth)
  dat = mutate(dat, Cat = str_extract(dat$FILE, "YesYes|YesNo|NoNo|NoYes"))
  dat$True.CV = round(dat$True.CV, 3)
  return(dat)
}


pwr_tab <- function(dat) {
  
  
  A1.dat = dat
  A1.basic = A1.dat %>% group_by(True.omega3.rate, True.CV) %>% summarise(num_reps = n())
  
  A1.BUSTED = A1.dat %>% group_by(True.omega3.rate,True.CV, Cat) %>% filter(BUSTED.P <
                                                                               0.05) %>% tally()
  A1.BUSTED = rename(A1.BUSTED, "BUSTED_PWR" = n)
  
  A1.SRV = A1.dat %>%  group_by(True.omega3.rate, True.CV, Cat) %>% filter(BUSTED.SRV.P <
                                                                              0.05) %>%  tally()
  A1.SRV = rename(A1.SRV, "SRV_PWR" = n)
  
  A1.pwr.tab = full_join(A1.BUSTED, A1.SRV, by = c("True.omega3.rate", "True.CV", "Cat")) %>% 
    full_join(., A1.basic, by = c("True.omega3.rate", "True.CV"))
  
  
  
  A1.means = A1.dat %>% group_by(True.omega3.rate) %>%   summarise(
    "$BUSTED \\omega_3$ MLE" = mean(BUSTED.omega3.MLE),
    "SRV $\\omega_3$ MLE" = mean(BUSTED.SRV.omega3.MLE),
    "Mean CV" = mean(CV.SRV)
  )
  A1.pwr.tab = full_join(A1.pwr.tab, A1.means, by = "True.omega3.rate")
  A1.pwr.tab = replace(A1.pwr.tab, is.na(A1.pwr.tab), 0)
  A1.pwr.tab$BUSTED_PWR = A1.pwr.tab$BUSTED_PWR / A1.pwr.tab$num_reps
  A1.pwr.tab$SRV_PWR = A1.pwr.tab$SRV_PWR / A1.pwr.tab$num_reps
  return(A1.pwr.tab)
}

all.dat = process_dat(dir, basename)

YesYes_1000 = all.dat %>% filter(Sites == 1000, True.omega3.rate > 1, True.CV >= 0)

sub = YesYes_1000  %>% group_by(True.CV) %>% tally %>% filter(n > 500)
a =YesYes_1000 %>% filter(True.CV %in% sub$True.CV)
a.pwr.tab = pwr_tab(a)
a.pwr.tab %>% select(one_of(
c("True.omega3.rate", "True.CV", "BUSTED_PWR", "SRV_PWR")
)) %>% melt(id.vars = c("True.omega3.rate", "True.CV")) %>% ggplot(aes(
x = True.omega3.rate,
y = value,
color = interaction(round(True.CV,3), variable, sep = " ")
)) +
geom_point() + geom_smooth() + labs(x = "Omega 3 Rate", y = "Power", color = "CV and Analysis")




#make tables???
a %>% glimpse

omega1.bias = a %>% group_by(True.CV) %>% summarise_at(vars(matches("omega1")), funs(mean,median))

 b= omega1.bias %>% gather(key = test, value = estimate, -True.CV, -True.omega1.rate_mean, -True.omega1.prop_mean, -True.omega1.rate_median, -True.omega1.prop_median) 

  b =b %>% mutate( analysis = str_extract(b$test,"BUSTED.SRV\\.|BUSTED\\.|True"))
 b = b %>% arrange(True.CV)
 b = mutate(b, 
            stat = str_extract(b$test,"mean|median"), idk = str_extract(b$test, "MLE|prop"),
            thing = interaction(stat,idk))

 glimpse(b)
c = spread(select(b, -test, - stat, -idk), key = thing, value = estimate)

#rename some shit

c =c %>% rename(True.Rate.Mean =True.omega1.rate_mean, True.Prop.Mean =True.omega1.prop_mean, 
             True.Rate.Median =True.omega1.rate_median,True.Prop.Median =True.omega1.prop_median)
c = c %>% select(True.CV,analysis, mean.MLE, True.Rate.Mean, median.MLE,mean.prop,True.Prop.Mean,
                 median.prop)


#Alpha tables


alpha1.bias = a %>% group_by(True.CV) %>% summarise_at(vars(matches("alpha1")), funs(mean,median))

b= alpha1.bias %>% gather(key = test, value = estimate, -True.CV, -True.alpha1.rate_mean, -True.alpha1.prop_mean, -True.alpha1.rate_median, -True.alpha1.prop_median) 

b =b %>% mutate( analysis = str_extract(b$test,"SRV\\.|BUSTED\\.|True"))
b = b %>% arrange(True.CV)
b = mutate(b, 
           stat = str_extract(b$test,"mean|median"), idk = str_extract(b$test, "MLE|prop"),
           thing = interaction(stat,idk))

glimpse(b)
c = spread(select(b, -test, - stat, -idk), key = thing, value = estimate)

#rename some shit

c =c %>% rename(True.Rate.Mean =True.alpha1.rate_mean, True.Prop.Mean =True.alpha1.prop_mean, 
                True.Rate.Median =True.alpha1.rate_median,True.Prop.Median =True.alpha1.prop_median)
c = c %>% select(True.CV,analysis, mean.MLE, True.Rate.Mean, median.MLE,mean.prop,True.Prop.Mean,
                 median.prop)

#calculate CV???
SRV.alpha3.MLE =0.3
SRV.alpha3.prop =0.1
SRV.alpha1.MLE =0.1
SRV.alpha1.prop = 0.5
SRV.alpha2.MLE =0.2
SRV.alpha2.prop =0.4
  
mom2 = SRV.alpha3.MLE^2*SRV.alpha3.prop+ SRV.alpha1.MLE^2*SRV.alpha1.prop+ SRV.alpha2.MLE^2*SRV.alpha2.prop
mean= SRV.alpha3.MLE*SRV.alpha3.prop+ SRV.alpha1.MLE*SRV.alpha1.prop+ SRV.alpha2.MLE*SRV.alpha2.prop
CV.SRV = sqrt(mom2-mean^2)/mean

temp=   a.pwr.tab %>% select(one_of(
  c("True.omega3.rate", "True.CV", "BUSTED_PWR", "SRV_PWR")
)) %>% filter(True.omega3.rate >= 1.1  && True.CV != 1.031) %>% melt(id.vars = c("True.omega3.rate", "True.CV")) 
temp %>% ggplot(aes(x=True.omega3.rate, y = True.CV))+ geom_tile(aes(fill = value))+ facet_grid(~variable)


#Treelength???
 glimpse(a)

 Tree.bias = a %>% group_by(True.CV, True.omega3.rate) %>% summarise_at(vars(matches("treelength")), funs(mean,median))
 
 b= Tree.bias %>% gather(key = test, value = estimate, -True.CV, - True.omega3.rate)
 
 b =ungroup(b) %>% mutate(analysis = str_extract(b$test,"BUSTED.SRV\\.|BUSTED\\.|True"))
 b = b %>% arrange(True.CV)
 b = mutate(b, 
            stat = str_extract(b$test,"mean|median"))
 
 glimpse(b)
 c = spread(select(b, -test), key = stat, value = estimate)
 c %>% kable
 
 
 #grappppph 
 
 
 omega1.bias = a %>% group_by(True.omega3.rate, True.CV) %>% summarise_at(vars(matches("omega1")), funs(mean,median))
 
 b= omega1.bias %>% gather(key = test, value = estimate, - True.omega3.rate, -True.CV, -True.omega1.rate_mean, -True.omega1.prop_mean, -True.omega1.rate_median, -True.omega1.prop_median) 
 
 b = ungroup(b) %>% mutate( analysis = str_extract(b$test,"BUSTED.SRV\\.|BUSTED\\.|True"))
 b = b %>% arrange(True.omega3.rate)
 b = mutate(b, 
            stat = str_extract(b$test,"mean|median"), idk = str_extract(b$test, "MLE|prop"),
            thing = interaction(stat,idk))
 
glimpse(b)
 c = spread(select(b, -test, - stat, -idk), key = thing, value = estimate)
 
 #rename some shit
 
 c =c %>% rename(True.Rate.Mean =True.omega1.rate_mean, True.Prop.Mean =True.omega1.prop_mean, 
                 True.Rate.Median =True.omega1.rate_median,True.Prop.Median =True.omega1.prop_median)
 c = c %>% select(True.omega3.rate, True.CV,analysis, mean.MLE, True.Rate.Mean, median.MLE,mean.prop,True.Prop.Mean,
                  median.prop)
 c %>% arrange(True.omega3.rate) %>% kable()
 
d = c %>% select(True.omega3.rate,True.CV, analysis,mean.MLE, True.Rate.Mean, median.MLE) %>% 
  melt(id.vars = c("True.omega3.rate","True.CV", "analysis"))

d  %>%  ggplot(aes(x = True.omega3.rate, y = value,color = variable )) + 
  geom_point() + geom_smooth()+ facet_grid(True.CV~analysis) 


#branch lengths

#From the LF file generated by BUSTED+SRV
#these would be our true values
srv.T_scaler = 4
srv.tree.Node1.t=0
srv.tree.2.t =srv.T_scaler*0.02289285117051735
srv.tree.3.t =srv.T_scaler*0.04231146330982721
srv.tree.Node5.t =srv.T_scaler*0.01750439883884518
srv.tree.4.t =srv.T_scaler*0.02419017369544371
srv.tree.5.t =srv.T_scaler*0.04335211232556986
srv.tree.1.t=srv.T_scaler*0.0953224280755873
srv.tree.Node4.t=srv.T_scaler*0.1641916117433053

#true values from BUSTED.LF (for 1000 Codons)
busted.T_scaler = 4

busted.tree.Node1.t=0;
busted.tree.2.t=busted.T_scaler*0.02289269236042673
busted.tree.3.t=busted.T_scaler*0.04231155126191058
busted.tree.Node5.t=busted.T_scaler*0.01750393389073822
busted.tree.4.t=busted.T_scaler*0.02419036132816281
busted.tree.5.t=busted.T_scaler*0.04335169428624834
busted.tree.1.t=busted.T_scaler*0.09532222084774838
busted.tree.Node4.t=busted.T_scaler*0.1641932190637932

#trying things out
#adding this function
add_truth <- function(dat, truth){
  dat = dat %>% mutate(True.CV = 1337)
  dat = dat %>% mutate(True.omega3.value = 8008)
  dat = dat %>% mutate(True.omega1.value = 8008)
  dat = dat %>% mutate(True.omega2.value = 8008)
  dat = dat %>% mutate(True.alpha3.value = 8008)
  dat = dat %>% mutate(True.alpha2.value = 8008)
  dat = dat %>% mutate(True.alpha1.value = 8008)
  dat = dat %>% mutate(True.omega3.prop = 8008)
  dat = dat %>% mutate(True.omega1.prop = 8008)
  dat = dat %>% mutate(True.omega2.prop = 8008)
  dat = dat %>% mutate(True.alpha3.prop = 8008)
  dat = dat %>% mutate(True.alpha2.prop = 8008)
  dat = dat %>% mutate(True.alpha1.prop = 8008)
  for( i in seq(from = 1, to = nrow(truth), by = 1)){
    temp = which(str_detect(dat$FILE,truth$Cat[i]))
    dat$True.CV[temp] = truth$CV.SRV[i] 
    dat$True.omega3.value[temp] = truth$Omega.3.value[i]
    dat$True.omega2.value[temp] = truth$Omega.2.value[i]
    dat$True.omega1.value[temp] = truth$Omega.1.value[i]
    dat$True.alpha1.value[temp] = truth$Alpha.1.value[i]
    dat$True.alpha2.value[temp] = truth$Alpha.2.value[i]
    dat$True.alpha3.value[temp] = truth$Alpha.3.value[i]
    dat$True.omega3.prop[temp] = truth$Omega.3.prop[i]
    dat$True.omega2.prop[temp] = truth$Omega.2.prop[i]
    dat$True.omega1.prop[temp] = truth$Omega.1.prop[i]
    dat$True.alpha1.prop[temp] = truth$Alpha.1.prop[i]
    dat$True.alpha2.prop[temp] = truth$Alpha.2.prop[i]
    dat$True.alpha3.prop[temp] = truth$Alpha.3.prop[i]
  }
  return(dat)
}
simulation_inputs <- function(dir,csv){
  require("stringr")
  require("jsonlite")
  require("dplyr")
  list = list.files(path = dir, recursive = T, pattern ="^([^.]+)$")
  #set up the empty data frame
  names = c("Sites" ,"Cat", "Omega 1 value", "Omega 2 value", "Omega 3 value", "Omega 1 prop", "Omega 2 prop", "Omega 3 prop",
            "Alpha 1 value", "Alpha 2 value", "Alpha 3 value","Alpha 1 prop", "Alpha 2 prop", "Alpha 3 prop")
  setup.tab = read.table(text = "",col.names = names)
  #loop thru each file to get info in correct format
  for(i in seq(from = 1, to= length(list))){
    x=readLines(paste(dir,list[i], sep = "/"))
    #making this a readable json
    x1 = str_replace(x,":\\{", ":\\[")
    
    x1[c(5:6,12:13)] = str_replace(x1[c(5:6,12:13)],"\\{([-+]?[0-9]*\\.?[0-9]+,[-+]?[0-9]*\\.?[0-9]+)\\}", "\\[\\1\\]," )
    
    x1[c(7,14)] = str_replace(x1[c(7,14)],"\\{([-+]?[0-9]*\\.?[0-9]+,[-+]?[0-9]*\\.?[0-9]+)\\}", "\\[\\1\\]" )
    x1[c(8,15)] = str_replace(x1[c(8,15)], "\\}", "\\]")
    r= fromJSON(x1)
    omega_rates = r$`omega distribution`[,1]
    names(omega_rates)= c("Omega 1 value", "Omega 2 value", "Omega 3 value")
    omega_weights = r$`omega distribution`[,2]
    names(omega_weights) = c("Omega 1 prop", "Omega 2 prop", "Omega 3 prop")
    
    Alpha_rates = r$`alpha distribution`[,1]
    names(Alpha_rates)= c("Alpha 1 value", "Alpha 2 value", "Alpha 3 value")
    Alpha_weights = r$`alpha distribution`[,2]
    names(Alpha_weights) = c("Alpha 1 prop", "Alpha 2 prop", "Alpha 3 prop")
    
    
    setup.tab[nrow(setup.tab)+1,] = c(r$sites,list[i], omega_rates,omega_weights,Alpha_rates,Alpha_weights)
  }
  setup.tab =select(setup.tab, Sites, Cat, Omega.1.value, Omega.1.prop, Omega.2.value, Omega.2.prop,Omega.3.value, Omega.3.prop,
                    Alpha.1.value, Alpha.1.prop, Alpha.2.value, Alpha.2.prop,Alpha.3.value, Alpha.3.prop)
  setup.tab[,3:14] = as.numeric(unlist(setup.tab[,3:14]))
  mom2 = setup.tab$Alpha.1.value^2*setup.tab$Alpha.1.prop+setup.tab$Alpha.2.value^2*setup.tab$Alpha.2.prop + setup.tab$Alpha.3.value^2*setup.tab$Alpha.3.prop
  
  mean = setup.tab$Alpha.1.value*setup.tab$Alpha.1.prop+setup.tab$Alpha.2.value*setup.tab$Alpha.2.prop + setup.tab$Alpha.3.value*setup.tab$Alpha.3.prop
  
  setup.tab = mutate(setup.tab,CV.SRV = sqrt(mom2-mean^2)/mean)
  
  return(setup.tab)
  #write.csv(file = csv, x = setup.tab, row.names= F)
}

#this function from required Funs.R
br = branch_length("C:/Users/srwisots/Drive/Muse/BRC/SimResults/Five_seq/1000_Codons/VaryO_lowCV/")
truth = simulation_inputs("C:/Users/srwisots/Drive/Muse/BRC/SimResults/Five_seq/1000_Codons/VaryO_lowCV/", csv = "truth_test.csv")
#truth = read.csv("C:/Users/srwisots/Drive/Muse/BRC/SimResults/Five_seq/All_Truth.csv", as.is = TRUE)
br.truth = add_truth(br,truth)
br.truth$length.SRV = as.numeric(br.truth$length.SRV)
br.truth$length.BUSTED = as.numeric(br.truth$length.BUSTED)

br.1 = br.truth %>% filter(Branch == 1) %>% group_by(True.omega3.value) %>% summarise_at(vars(matches("length")), funs(mean, median))

t =br.1  %>%
  melt(id.vars = c("True.omega3.value", "length.BUSTED_median", "length.SRV_median"), variable.name = "mean_method", value.name = "mean") %>% 
  melt(id.vars = c("True.omega3.value", "mean_method", "mean"), variable.name = "median_method", value.name = "median")

t %>% ggplot(aes(x = True.omega3.value))+ geom_point(aes(y = mean, color = mean_method))+ geom_point(aes(y=median, color = mean_method))+
  geom_hline(yintercept = srv.tree.1.t, color = 'red')+
  geom_vline(xintercept = busted.tree.1.t, color = 'red')

t = br.1 %>% melt(id.vars = c("True.omega3.value")) %>% mutate(method = str_extract(variable,"BUSTED|SRV|truth"), 
                                                               measure = str_extract(variable, "mean|median|truth")
                                                                )

t %>% ggplot(aes(x = True.omega3.value, y = value, color = method, shape = measure))+ geom_point()+
  geom_hline(yintercept = srv.tree.1.t, color = 'red')


br.1 %>% ggplot(aes(x = length.BUSTED, y = length.SRV))+ geom_point(aes(color= True.omega3.value))+ geom_hline(yintercept = srv.tree.1.t, color = 'red')+
  geom_vline(xintercept = busted.tree.1.t, color = 'red')

br.truth %>% filter(Branch == 2) %>% ggplot(aes(x = length.BUSTED, y = length.SRV))+ geom_point(aes(color= True.omega3.value))+ 
  geom_hline(yintercept = srv.tree.2.t, color = 'red')+
  geom_vline(xintercept = busted.tree.2.t, color = 'red')

br.truth %>% filter(Branch == 3) %>% ggplot(aes(x = length.BUSTED, y = length.SRV))+ geom_point(aes(color= True.omega3.value))+ 
  geom_hline(yintercept = srv.tree.3.t, color = 'red')+
  geom_vline(xintercept = busted.tree.3.t, color = 'red')

#compiles all on the cluster. (see branch_length.R on the BRC cluster)
br.all = read.csv("C:/Users/srwisots/Drive/Muse/BRC/SimResults/Five_seq/Branch_lengths.csv", as.is = TRUE) #reading in the csv that was generated

#adding sites cause this was missing and it makes things easier to filter
br.all= br.all %>% mutate(Sites = str_extract(FILE, "1000|5000|10000")) #technically 5000 should be 6000 maybe look at way to change this
                                                                        #could use Sites from the Truth part as they're there

#get just the branch lengths for 1000 Codons with selection 
#trying to get same subset I showed spencer
idk  = br.all %>% filter(Sites == 1000, True.omega3.value > 1, True.CV >= 0) 

#filter out the ones with weird CV values 
#they exsisted cause I di vary CV while holding omega3 = 2 at one point

sub = idk  %>% group_by(True.CV) %>% tally %>% filter(n > 8400)
#a is going to be what we actually end up using
#it's correctly filtered and everything
a =idk %>% filter(True.CV %in% sub$True.CV)

#now we're just taking branch 1 and getting the mean and median length as grouped by true omega3 and true CV

br.1 = a %>% filter(Branch == 1)%>% group_by(True.omega3.value, True.CV) %>% summarise_at(vars(matches("length")), funs(mean, median))

#melting things to allow us to color and facet as we need
t = br.1 %>% melt(id.vars = c("True.omega3.value", "True.CV")) %>% mutate(method = str_extract(variable,"BUSTED|SRV|truth"), 
                                                               measure = str_extract(variable, "mean|median|truth")
)

#plot it all
#facet by method and by True CV 
#color by measure
#red hline is the true value for branch 1 according to what I pulled out of the LF 
t %>% ggplot(aes(x = True.omega3.value, y = value, color = measure))+ geom_point()+ facet_grid(True.CV~method)+
  geom_hline(yintercept = srv.tree.1.t, color = 'red')


##########explore bias with the longer sequence lengths
#subset the data

all.dat  = read.csv("C:/Users/srwisots/Drive/Muse/BRC/SimResults/Five_seq/All_12_7_16_processed.csv", as.is = TRUE)

#subset for 10,0000
sub.10 = all.dat %>% filter(Sites == 10000, True.omega3.value > 1, True.CV >= 0, True.CV != 1.906, True.omega3.value <= 2)


#subset for 6,000

sub.6 = all.dat %>% filter(Sites==6000, True.omega3.value > 1, True.CV >= 0, True.CV != 1.906, True.omega3.value <= 2)

#omega1 for 10,000

omega1.bias = sub.10 %>% group_by(True.omega3.value, True.CV) %>% summarise_at(vars(matches("omega1")), funs(mean,median))

b= omega1.bias %>% gather(key = test, value = estimate, - True.omega3.value, -True.CV, -True.omega1.value_mean, -True.omega1.prop_mean, -True.omega1.value_median, -True.omega1.prop_median) 

b = ungroup(b) %>% mutate( analysis = str_extract(b$test,"BUSTED.SRV\\.|BUSTED\\.|True"))
b = b %>% arrange(True.omega3.value)
b = mutate(b, 
           stat = str_extract(b$test,"mean|median"), idk = str_extract(b$test, "MLE|prop"),
           thing = interaction(stat,idk))

# glimpse(b)
c = spread(select(b, -test, - stat, -idk), key = thing, value = estimate)

#rename some shit

c =c %>% rename(True.value.Mean =True.omega1.value_mean, True.Prop.Mean =True.omega1.prop_mean, 
                True.value.Median =True.omega1.value_median,True.Prop.Median =True.omega1.prop_median)
c = c %>% select(True.omega3.value, True.CV,analysis, mean.MLE, True.value.Mean, median.MLE,mean.prop,True.Prop.Mean,
                 median.prop)

#that works let us compare between seq length

all.sub = rbind(sub.1,sub.6,sub.10)
 all.sub = all.sub %>%  filter(True.CV %in% c(0, 0.415), True.omega3.value <=2)

 #just using the median as our stat
omega1.bias = all.sub %>% group_by(True.omega3.value, True.CV, Sites) %>% summarise_at(vars(matches("omega1")), median)

b= omega1.bias %>% gather(key = test, value = estimate, - True.omega3.value, -True.CV, -Sites, -True.omega1.value, -True.omega1.prop) 

b = mutate(ungroup(b), stat = str_extract(b$test, "MLE|prop"), analysis = str_extract(b$test,"BUSTED.SRV\\.|BUSTED\\.|True"))

b %>%  filter( stat == "MLE") %>% ggplot(aes(x =True.omega3.value, y = estimate, color = factor(Sites))) + geom_point() + geom_smooth(se =F)+
  facet_grid(True.CV ~ analysis) + geom_hline(yintercept = b$True.omega1.value, color = "red", show.legend = TRUE))
