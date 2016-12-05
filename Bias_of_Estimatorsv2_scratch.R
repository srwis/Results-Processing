
library(jsonlite)

library(stringr)

library(dplyr)



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
  setup.tab =select(setup.tab, Sites, Cat, Omega.1.value, Omega.1.prop, Omega.2.value, Omega.2.prop,Omega.3.value, Omega.3.prop,
                    Alpha.1.value, Alpha.1.prop, Alpha.2.value, Alpha.2.prop,Alpha.3.value, Alpha.3.prop)
  setup.tab[,3:14] = as.numeric(unlist(setup.tab[,3:14]))
  mom2 = setup.tab$Alpha.1.value^2*setup.tab$Alpha.1.prop+setup.tab$Alpha.2.value^2*setup.tab$Alpha.2.prop + setup.tab$Alpha.3.value^2*setup.tab$Alpha.3.prop
  
  mean = setup.tab$Alpha.1.value*setup.tab$Alpha.1.prop+setup.tab$Alpha.2.value*setup.tab$Alpha.2.prop + setup.tab$Alpha.3.value*setup.tab$Alpha.3.prop
  
  setup.tab = mutate(setup.tab,CV.SRV = sqrt(mom2-mean^2)/mean)
  
  #return(setup.tab)
  write.csv(file = csv, x = setup.tab, row.names= F)
}


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



process_dat <- function(dir, basename){
  temp = paste(dir,basename,"_Truth.csv", sep = "")
  #simulation_inputs(dir,temp)
  truth = read.csv(temp, as.is = T)
  temp = paste(dir,basename,"_results.csv", sep = "")
  # compile(dir,temp)
  dat = read.csv(temp, as.is = T)
  dat = add_truth(dat, truth)
  dat = mutate(dat, Cat = str_extract(dat$FILE, "YesYes|YesNo|NoNo|NoYes"))
  dat$True.CV = round(dat$True.CV, 3)
  return(dat)
}




















```{r graph for omega1, echo = FALSE, fig.cap="Bias of $\\omega_1$ value across increasing $\\omega_3$ values and facted by the analysis type used and the CV of SRV simulated with. "}

d = c %>% select(True.omega3.value,True.CV, analysis,mean.MLE, True.value.Mean, median.MLE) %>%
  melt(id.vars = c("True.omega3.value","True.CV", "analysis"))

d  %>%  ggplot(aes(x = True.omega3.value, y = value,color = variable )) +
  geom_point() + geom_smooth()+ facet_grid(True.CV~analysis, scales = "free") + 
  labs(x = expression("True "*omega[3]), y = expression("Expected "*omega[1]))+
  scale_color_discrete(name = "Measure", breaks = c("mean.MLE", "True.value.Mean", "median.MLE"), labels = c("Mean MLE", expression("True "*omega[1]), "Median MLE"))
```

```{r, echo =FALSE}
d = c %>% select(True.omega3.value,True.CV, analysis,mean.prop, True.Prop.Mean, median.prop) %>%
  melt(id.vars = c("True.omega3.value","True.CV", "analysis"))

d  %>%  ggplot(aes(x = True.omega3.value, y = value,color = variable )) +
  geom_point() + geom_smooth()+ facet_grid(True.CV~analysis, scales = "free") +
  labs(x = expression("True "*omega[3]), y = expression("Expected "*omega[1]))+
  scale_color_discrete(name = "Measure", breaks = c("mean.prop", "True.Prop.Mean", "median.prop"), labels = c("Mean prop", expression("True "*omega[1]*" prop"), "Median prop"))


```