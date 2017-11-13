#Assorted function from old stuff that might be useful???
#variable names from the JSONs may have changed so expect errors




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