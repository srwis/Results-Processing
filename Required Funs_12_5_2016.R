#Required functions to process results
#and loading libraries

library(jsonlite)
library(ggplot2)
library(knitr)
library(stringr)
library(reshape2)
library(dplyr)

#Right now these are just he BUSTED-SRV branch lengths

cur.dir = "C:/Users/srwisots/Drive/Muse/BRC/SimResults/Five_seq/1000_Codons/VaryO_lowCV/"
branch_length <- function(cur.dir){
  jsons <- list.files(path = cur.dir,
                      pattern = '*.json', recursive = TRUE)
  
  #create empty data.frame with variable names
  
  df = read.table(text="", col.names = c("Branch", "FILE", "length.SRV","length.BUSTED"))
  for (i in  seq(from=1, to=length(jsons), by=2)){
    
    
    FILE = jsons[i]
    
    
    if(grepl("SRV", jsons[i])){
      method = "SRV"
      filepath = paste(cur.dir,jsons[i], sep="")
      
      test = fromJSON(filepath)
      tree_string = test$fits$`Unconstrained model`$`tree string`
      
      x= tree_string %>% str_replace_all("\\(",":") %>% str_replace_all("\\)",":") %>%
        str_replace_all(",",":") %>% str_split(":")
      x= unlist(x)
      x =x[x !=""]
      br_len = matrix(x,ncol = 2,  byrow = T)
      
    }
    length.srv = br_len[,2]
    
    if(grepl("SRV",jsons[i+1])==FALSE){
      method = "BUSTED"   
      filepath = paste(cur.dir,jsons[i+1], sep="")
      
      test = fromJSON(filepath)
      tree_string = test$fits$`Unconstrained model`$`tree string`
      
      x= tree_string %>% str_replace_all("\\(",":") %>% str_replace_all("\\)",":") %>% 
        str_replace_all(",",":") %>% str_split(":")
      x= unlist(x)
      x =x[x !=""]
      br_len = matrix(x,ncol = 2,  byrow = T)
      
    }
    length.BUSTED = br_len[,2]
    for(k in seq(from=1,to=length(br_len)/2)){
      df[nrow(df)+1,]=c(br_len[k,1], FILE, length.srv[k],length.BUSTED[k])
    }
  }
  #df_wide = spread(df,File, length)
  return(df)
}

#adds significant categories and normalizes alphas
process.dat <- function(dat, redo.alphas = TRUE){
  ### Add sig categories
  dat= cbind(dat,p.cat=cut(dat$BUSTED.P,breaks=c(-Inf,0.05,Inf),
                           labels = c("S V B", "NS V B")))
  dat= cbind(dat,p.srv.cat=cut(dat$BUSTED.SRV.P,breaks=c(-Inf,0.05,Inf),
                               labels = c("S V BSRV", "NS V BSRV")))
  ###rearrange alphas
  if(redo.alphas==TRUE){
    normalize = dat$SRV.alpha1.MLE*dat$SRV.alpha1.prop + 
      dat$SRV.alpha2.MLE*dat$SRV.alpha2.prop+ 
      dat$SRV.alpha3.MLE*dat$SRV.alpha3.prop
    
    
    dat = cbind(dat,true.alpha1 = dat$SRV.alpha1.MLE/normalize,
                true.alpha2 = dat$SRV.alpha2.MLE/normalize,
                true.alpha3 = dat$SRV.alpha3.MLE/normalize)
  }
  else{
    #sets the true.alpha values to just the standard MLE 
    #makes the rest of the stream line analysis easier. 
    #it's either this or make you specify alpha variables in later steps...
    dat = cbind(dat,true.alpha1 = dat$SRV.alpha1.MLE,
                true.alpha2 = dat$SRV.alpha2.MLE,
                true.alpha3 = dat$SRV.alpha3.MLE)  }
  
}

# here's where you'd have to start specifying other alphas 
# returns a list of statistics generated from the boxplot function 
# note the omega3 value here is from BUSTED SRV 
# not the omega3 from BUSTED for the same file
box.sum.stats <- function(dat, subset = NULL){
  omega3 = boxplot(
    dat$BUSTED.SRV.omega3.MLE ~ dat$p.cat + dat$p.srv.cat, data = dat, 
    subset = subset, plot = FALSE
  )
  
  alpha1 = boxplot(
    dat$true.alpha1 ~ dat$p.cat + dat$p.srv.cat, data = dat,
    subset = subset, plot = FALSE
  )
  
  alpha2 = boxplot(
    dat$true.alpha2 ~ dat$p.cat + dat$p.srv.cat, data = dat,
    subset = subset, plot = FALSE
  )
  
  alpha3 = boxplot(
    dat$true.alpha3 ~ dat$p.cat + dat$p.srv.cat, data = dat,
    subset = subset,plot = FALSE
  )
  results <- list(omega3 =omega3, alpha1 = alpha1,alpha2= alpha2, alpha3 = alpha3)
  return(results)
}

gen.tables <- function(box.stats){
  catnames = c("Sel for both", "No Sel for BUSTED", "No Sel for SRV", "No Sel for both")
  statlist = c("lower whisker","lower hinge","median","upper hinge","upper whisker", "number of files")
  omega3 = box.stats$omega3$stats
  omega3 = rbind(omega3,box.stats$omega3$n)
  colnames(omega3) = catnames
  rownames(omega3) = statlist
  
  alpha1 = box.stats$alpha1$stats
  alpha1 = rbind(alpha1,box.stats$alpha1$n)
  colnames(alpha1) = catnames
  rownames(alpha1) = statlist
  
  alpha2 = box.stats$alpha2$stats
  alpha2 = rbind(alpha2,box.stats$alpha2$n)
  colnames(alpha2) = catnames
  rownames(alpha2) = statlist
  
  alpha3 = box.stats$alpha3$stats
  alpha3 = rbind(alpha3,box.stats$alpha3$n)
  colnames(alpha3) = catnames
  rownames(alpha3) = statlist
  
  results <- list(omega3 =omega3, alpha1 = alpha1,alpha2= alpha2, alpha3 = alpha3)
  return(results)
}

gen.summary.o.tabs <- function(dat){
  omega1.busted =  summary(dat$BUSTED.omega1.MLE)
  omega1.srv = summary(dat$BUSTED.SRV.omega1.MLE)
  omega2.busted =summary(dat$BUSTED.omega2.MLE)
  omega2.srv = summary(dat$BUSTED.SRV.omega2.MLE)
  omega3.busted =summary(dat$BUSTED.omega3.MLE)
  omega3.srv = summary(dat$BUSTED.SRV.omega3.MLE)
  results = rbind(omega1.busted,omega1.srv,omega2.busted,omega2.srv,omega3.busted,omega3.srv)
  return(results)
}
gen.summary.a.tabs <- function(dat){
  
  alpha1.srv = summary(dat$true.alpha1)
  alpha2.srv = summary(dat$true.alpha2)
  alpha3.srv = summary(dat$true.alpha3)
  results = rbind(alpha1.srv, alpha2.srv,alpha3.srv)
  return(results)
}

gen.sig.table <- function(dat){
  require("xtable")
  library("xtable")
  
  under.sel.busted = which(dat$BUSTED.P<=0.05)
  under.sel.srv = which(dat$BUSTED.SRV.P<=0.05)
  if( length(under.sel.srv) == 0){
    
    mat = matrix( rep(0),nrow=length(dat$BUSTED.P),ncol = 2,
                  dimnames= list(1:length(dat$BUSTED.P),
                                 c("BUSTED", "BUSTED-SRV")))
    mat[under.sel.busted,1] = 1
    mat[under.sel.srv,2] = 1 
    sel.tab = table(mat[,1],mat[,2], dnn = colnames(mat))
    row.names(sel.tab) = c("Selection:BUSTED")
    colnames(sel.tab) = c("No Selection:BUSTED+SRV")
    Sel.Prop =prop.table(sel.tab)
    test.table = xtable(Sel.Prop)
    #print.xtable(test.table, type = "html")
    #print(kable(test.table))
    #return("All alignments are not significant according to BUSTED+SRV")
  }
  else if( length(under.sel.busted) == 0){
    
    mat = matrix( rep(0),nrow=length(dat$BUSTED.P),ncol = 2,
                  dimnames= list(1:length(dat$BUSTED.P),
                                 c("BUSTED", "BUSTED-SRV")))
    mat[under.sel.busted,1] = 1
    mat[under.sel.srv,2] = 1 
    sel.tab = table(mat[,1],mat[,2], dnn = colnames(mat))
    row.names(sel.tab) = c("No Selection:BUSTED")
    colnames(sel.tab) = c("Selection:BUSTED+SRV")
    Sel.Prop =prop.table(sel.tab)
    test.table = xtable(Sel.Prop)
    #print.xtable(test.table, type = "html", file="sel-table.html")
    #return("All alignments are not significant according to BUSTED")
  }
  else {
    mat = matrix( rep(0),nrow=length(dat$BUSTED.P),ncol = 2,
                  dimnames= list(1:length(dat$BUSTED.P),
                                 c("BUSTED", "BUSTED-SRV")))
    mat[under.sel.busted,1] = 1
    mat[under.sel.srv,2] = 1 
    sel.tab = table(mat[,1],mat[,2], dnn = colnames(mat))
    row.names(sel.tab) = c("No Selection", "Selection")
    colnames(sel.tab) = c("No Selection", "Selection")
    Sel.Prop =prop.table(sel.tab)
    test.table = xtable(Sel.Prop)
    #print.xtable(test.table, type = "html", file="sel-table.html")
  }
  return(Sel.Prop)
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
  
  #return(setup.tab)
  write.csv(file = csv, x = setup.tab, row.names= F)
}



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
    
    if (grepl("BUSTED-SRV",jsons[i+1])){
      filepath = paste(cur.dir,jsons[i+1], sep="")
      
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
    if (grepl("BUSTED.json",jsons[i])){
      filepath = paste(cur.dir,jsons[i], sep="")
      
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



pwr_tab_omegas<- function(csv_data, csv_truth){
  A1.dat = read.csv(csv_data, as.is = TRUE)
  A1.dat = str_extract(A1.dat$FILE, "2_[0-9]*|1_[0-9]*") %>% str_replace("_",".") %>% as.numeric() %>% mutate(A1.dat, True.Rate.Omega.3 = .)
  
  A1.True.vals = read.csv(csv_truth, as.is = TRUE)
  A1.basic = A1.dat %>% group_by(True.Rate.Omega.3) %>% summarise(num_reps = n())
  
  A1.BUSTED = A1.dat %>% group_by(True.Rate.Omega.3) %>% filter(BUSTED.P<0.05) %>% tally()
  colnames(A1.BUSTED)[2] ="BUSTED_PWR"
  
  A1.SRV = A1.dat %>%  group_by(True.Rate.Omega.3) %>% filter(BUSTED.SRV.P<0.05) %>%  tally()
  colnames(A1.SRV)[2]= "SRV_PWR"
  
  A1.pwr.tab = full_join(A1.BUSTED,A1.SRV, by = "True.Rate.Omega.3") %>% full_join(.,A1.basic, by = "True.Rate.Omega.3") %>% 
    mutate(True.SRV.CV = A1.True.vals$CV.SRV[1])
  
  A1.means = A1.dat %>% group_by(True.Rate.Omega.3) %>%   summarise("$BUSTED \\omega_3$ MLE" = mean(BUSTED.omega3.MLE),
                                                                    "SRV $\\omega_3$ MLE" = mean(BUSTED.SRV.omega3.MLE),
                                                                    "Mean CV" = mean(CV.SRV))
  A1.pwr.tab = full_join(A1.pwr.tab, A1.means, by = "True.Rate.Omega.3")
  A1.pwr.tab = replace(A1.pwr.tab, is.na(A1.pwr.tab), 0)
  A1.pwr.tab$BUSTED_PWR = A1.pwr.tab$BUSTED_PWR/A1.pwr.tab$num_reps
  A1.pwr.tab$SRV_PWR = A1.pwr.tab$SRV_PWR/A1.pwr.tab$num_reps
  return(A1.pwr.tab)
}

pwr_tab_alpha <- function(NoYes.O2.dat, NoYes.O2.Truth) {
  NoYes.O2.dat = read.csv(
    "E:/BRC/SimResults/Five_seq/1000_Codons/noSel_yesSRV/NoYes_VaryA_O2.csv", as.is = TRUE)
  NoYes.O2.Truth = read.csv(
    "E:/BRC/SimResults/Five_seq/1000_Codons/noSel_yesSRV/NoYes_VaryA_O2_Truth.csv", as.is = TRUE)
  #add CV column for easier manipulation
  NoYes.O2.dat = NoYes.O2.dat %>% mutate(True.CV = 1337)
  #add omega 3 for same reason
  NoYes.O2.dat = NoYes.O2.dat %>% mutate(True.Rate.Omega.3 = 8008)
  for( i in seq(from = 1, to = length(NoYes.O2.Truth), by = 1)){
    temp = which(str_detect(NoYes.O2.dat$FILE,NoYes.O2.Truth$Cat[i]))
    NoYes.O2.dat$True.CV[temp] = NoYes.O2.Truth$CV.SRV[i] 
    NoYes.O2.dat$True.Rate.Omega.3[temp] = NoYes.O2.Truth$Omega.3.Rate[i]
    
  }
  
  O2.basic = NoYes.O2.dat %>% group_by(True.CV) %>% summarise(num_reps = n())
  
  O2.BUSTED = NoYes.O2.dat %>% group_by(True.CV) %>% filter(BUSTED.P<0.05) %>% tally()
  colnames(O2.BUSTED)[2] ="BUSTED_PWR"
  
  O2.SRV = NoYes.O2.dat %>%  group_by(True.CV) %>% filter(BUSTED.SRV.P<0.05) %>%  tally()
  colnames(O2.SRV)[2]= "SRV_PWR"
  
  O2.pwr.tab = full_join(O2.BUSTED,O2.SRV, by = "True.CV") %>% full_join(.,O2.basic, by = "True.CV") %>% 
    
    
    O2.means = NoYes.O2.dat %>% group_by(True.CV) %>%   summarise("$BUSTED \\omega_3$ MLE" = mean(BUSTED.omega3.MLE),
                                                                  "SRV $\\omega_3$ MLE" = mean(BUSTED.SRV.omega3.MLE),
                                                                  "Mean CV" = mean(CV.SRV),
                                                                  "True.Rate.Omega3" = mean(True.Rate.Omega.3))
  O2.pwr.tab = full_join(O2.pwr.tab, O2.means, by = "True.CV")
  O2.pwr.tab = replace(O2.pwr.tab, is.na(O2.pwr.tab), 0)
  O2.pwr.tab$BUSTED_PWR = O2.pwr.tab$BUSTED_PWR/O2.pwr.tab$num_reps
  O2.pwr.tab$SRV_PWR = O2.pwr.tab$SRV_PWR/O2.pwr.tab$num_reps
  return(O2.pwr.tab)
}





type_I_error <- function(dat){
  
  
  basic.noSel = dat %>% filter(cat == "NoNo"| cat == "NoYes") %>% group_by(True.Rate.Omega.3,Sites, cat,True.SRV.CV) %>% summarise(num_reps = n())
  
  #Type I error
  #number of reps that report selection when there is none
  #aka they incorrectly rejest the null hypothesis
  type.I.BUSTED = dat %>% filter(cat == "NoNo"| cat == "NoYes") %>% 
    group_by(True.Rate.Omega.3, Sites, cat, True.SRV.CV) %>% filter(BUSTED.P < 0.05 ) %>% tally()
  type.I.BUSTED = rename(type.I.BUSTED, BUSTED.Type.I.Error = n)
  
  type.I.SRV  =  dat %>% filter(cat == "NoNo"| cat == "NoYes") %>% 
    group_by(True.Rate.Omega.3,Sites, cat,True.SRV.CV) %>% filter(BUSTED.SRV.P < 0.05 ) %>% tally()
  type.I.SRV= rename(type.I.SRV,"SRV.Type.I.Error" = n)
  
  type.I.error = full_join(type.I.SRV, type.I.BUSTED, by = c("cat","Sites","True.Rate.Omega.3","True.SRV.CV")) %>% 
    full_join(., basic.noSel, by = c("cat", "Sites","True.Rate.Omega.3", "True.SRV.CV"))
  type.I.error = replace(type.I.error, is.na(type.I.error),0)
  return(type.I.error)
  
}


type_II_error <- function(dat){
  #Type II error
  #the number of reps that report there is no selection when there is
  #akak they fail to reject the null hypothesis when they should
  
  basic.YesSel = dat %>% filter(cat == "YesNo"| cat == "YesYes") %>% group_by(True.Rate.Omega.3,Sites, cat,True.SRV.CV) %>% summarise(num_reps = n())
  
  type.II.BUSTED = dat %>% filter(cat == "YesNo"| cat == "YesYes") %>% 
    group_by(True.Rate.Omega.3,Sites, cat,True.SRV.CV) %>% filter(BUSTED.P > 0.05 ) %>% tally()
  type.II.BUSTED = rename(type.II.BUSTED, BUSTED.Type.II.Error = n)
  
  type.II.SRV  =  dat %>% filter(cat == "YesNo"| cat == "YesYes") %>% 
    group_by(True.Rate.Omega.3,Sites, cat,True.SRV.CV) %>% filter(BUSTED.SRV.P > 0.05 ) %>% tally()
  type.II.SRV = rename(type.II.SRV, SRV.Type.II.Error = n)
  
  type.II.error = full_join(type.II.SRV, type.II.BUSTED, by = c("cat","Sites","True.Rate.Omega.3","True.SRV.CV")) %>% 
    full_join(., basic.YesSel, by = c("cat", "Sites", "True.Rate.Omega.3", "True.SRV.CV"))
  type.II.error= replace(type.II.error, is.na(type.II.error), 0)
  return(type.II.error)
  
}