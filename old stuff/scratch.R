
library(jsonlite)
library(ggplot2)
library(knitr)
library(stringr)
library(reshape2)
library(dplyr)
library(xtable)



A2.dat =  read.csv("G:/BRC/SimResults/Five_seq/YesYes_A2_VaryO.csv", as.is = TRUE)
A2.dat =str_extract(dat$FILE, "2_[0-9]*") %>% str_replace("_",".") %>% as.numeric() %>% mutate(dat, True.Rate.Omega.3 = .)

A2.True.vals = read.csv("G:/BRC/SimResults/Five_seq/YesYes_A2_VaryO_Truth.csv",as.is = T)


A6.basic = all_omega_vary %>% group_by(True.Rate.Omega.3, cat) %>% summarise(num_reps =n()) %>% filter(cat == "YesYes")
A2.basic = dat %>% group_by(True.Rate.Omega.3) %>% summarise(num_reps = n())

simulation_inputs("G:/BRC/SimResults/Five_seq/1000_Codons/yesSel_yesSRV/VaryO_A6/","G:/BRC/SimResults/Five_seq/YesYes_A6_VaryO_Truth.csv" )
A6.True.vals = read.csv("G:/BRC/SimResults/Five_seq/YesYes_A6_VaryO_Truth.csv", as.is = T)

#calculate power
A6_correct_BUSTED = all_omega_vary %>% group_by(True.Rate.Omega.3, cat)  %>% filter(cat == "YesYes") %>% filter(BUSTED.P<0.05) %>% tally() 
colnames(A6_correct_BUSTED )[3] = "BUSTED_PWR"

A6_correct = all_omega_vary %>% group_by(True.Rate.Omega.3, cat)  %>% filter(cat == "YesYes") %>% filter(BUSTED.SRV.P<0.05) %>% tally() 
colnames(A6_correct)[3] = "SRV_PWR"

A6.pwr.tab = full_join(A6_correct_BUSTED,A6_correct, by = c("True.Rate.Omega.3", "cat")) %>%
    full_join(.,A6.basic, by = c("True.Rate.Omega.3", "cat"))

A6.pwr.tab = A6.pwr.tab %>% mutate(True.SRV.CV = A6.True.vals$CV.SRV[1]) %>% select(-cat)

A6.means = all_omega_vary %>% group_by(True.Rate.Omega.3, cat)  %>% filter(cat == "YesYes") %>%  
  summarise("$BUSTED \\omega_3$ MLE" = mean(BUSTED.omega3.MLE),
            "SRV $\\omega_3$ MLE" = mean(BUSTED.SRV.omega3.MLE),
            "Mean CV" = mean(CV.SRV))
A6.pwr.tab = full_join(A6.pwr.tab, A6.means, by = "True.Rate.Omega.3")

A2_correct_SRV = A2.dat %>% group_by(True.Rate.Omega.3) %>% filter(BUSTED.SRV.P<0.05) %>% tally()
colnames(A2_correct_SRV)[2]="SRV_PWR"

A2_correct_BUSTED = A2.dat %>% group_by(True.Rate.Omega.3) %>% filter(BUSTED.P<0.05) %>% tally()
colnames(A2_correct_BUSTED)[2]="BUSTED_PWR"
A2.pwr.tab = full_join(A2_correct_BUSTED, A2_correct_SRV, by = "True.Rate.Omega.3") %>%
  full_join(., A2.basic, by = "True.Rate.Omega.3")
A2.pwr.tab = A2.pwr.tab %>% mutate(True.SRV.CV = A2.True.vals$CV.SRV[1])
A2.means = A2.dat %>% group_by(True.Rate.Omega.3) %>%   summarise("$BUSTED \\omega_3$ MLE" = mean(BUSTED.omega3.MLE),
                                                                  "SRV $\\omega_3$ MLE" = mean(BUSTED.SRV.omega3.MLE),
                                                                  "Mean CV" = mean(CV.SRV))
A2.pwr.tab = full_join(A2.pwr.tab, A2.means, by = "True.Rate.Omega.3")


####A1
simulation_inputs("G:/BRC/SimResults/Five_seq/1000_Codons/VaryO_A1/","G:/BRC/SimResults/Five_seq/YesYes_A1_VaryO_Truth.csv" )

A1.dat = read.csv("G:/BRC/SimResults/Five_seq/YesYes_A1_VaryO.csv", as.is = TRUE)
A1.dat = str_extract(A1.dat$FILE, "O3_[0-9]_[0-9]*") %>% str_replace("_",".") %>% as.numeric() %>% mutate(A1.dat, True.Rate.Omega.3 = .)


O2.dat =read.csv("G:/BRC/SimResults/Five_seq/1000_Codons/noSel_yesSRV/NoYes_VaryA_O2.csv",as.is = TRUE)

O2.True.vals = read.csv("G:/BRC/SimResults/Five_seq/1000_Codons/noSel_yesSRV/NoYes_VaryA_O2_Truth.csv", as.is = TRUE)

str_extract(A1.dat$FILE, "_[0-9]_[0-9]*") %>% sub("_", "", x = .)%>% str_replace("_",".") %>% as.numeric()



A1.True.vals = read.csv("G:/BRC/SimResults/Five_seq/YesYes_A1_VaryO_Truth.csv", as.is = TRUE)
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

pwr.tab = bind_rows(A6.pwr.tab,A2.pwr.tab, A1.pwr.tab)
pwr.tab = replace(pwr.tab, is.na(pwr.tab), 0)
pwr.tab$BUSTED_PWR = pwr.tab$BUSTED_PWR/pwr.tab$num_reps
pwr.tab$SRV_PWR = pwr.tab$SRV_PWR/pwr.tab$num_reps
pwr.tab
pwr.tab =pwr.tab %>% arrange(True.Rate.Omega.3) %>% select(True.Rate.Omega.3,`$BUSTED \\omega_3$ MLE`, `SRV $\\omega_3$ MLE`, True.SRV.CV, `Mean CV`, BUSTED_PWR,SRV_PWR)

kable(pwr.tab)
pwr.tab %>% select(one_of(c("True.Rate.Omega.3", "True.SRV.CV","BUSTED_PWR", "SRV_PWR"))) %>% melt(id.vars= c("True.Rate.Omega.3", "True.SRV.CV")) %>% ggplot(aes(x= True.Rate.Omega.3, 
                                                                            y = value,color =interaction(True.SRV.CV,variable))) + 
  geom_point() + geom_smooth() +labs( x = "Omega 3 Rate", y = "Power", color = "CV and Analysis") 




#function to output power table
#requires csv of data and csv of simulation inputs
pwr_tab <- function(csv_data, csv_truth){
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





new.cats = c(
  "'Sel'^'+'* 'SRV'^'+'",
  "'Sel'^'+' * 'SRV'^'-'",
  "'Sel'^'-'*'SRV'^'-'",
  "'Sel'^'-'*'SRV'^'+'")





parse(text =levels(a))



SRV.alpha2.MLE = 0.9999999999390909
SRV.alpha3.MLE = 1.090909090842645
SRV.alpha1.MLE = 0.9090909090355371
SRV.alpha3.prop = 0.33333333367
SRV.alpha2.prop = 0.33333333333
SRV.alpha1.prop = 0.333333333
mom2 = SRV.alpha3.MLE^2*SRV.alpha3.prop+ SRV.alpha1.MLE^2*SRV.alpha1.prop+ SRV.alpha2.MLE^2*SRV.alpha2.prop
mean= SRV.alpha3.MLE*SRV.alpha3.prop+ SRV.alpha1.MLE*SRV.alpha1.prop+ SRV.alpha2.MLE*SRV.alpha2.prop
CV.SRV = sqrt(mom2-mean^2)/mean








#figure out how to redo fig 3

all_omega_vary = rbind( mutate(YesNo_Omega_all.dat, cat = "YesNo"),mutate(YesYes_all.dat, cat = "YesYes"))

basic = all_omega_vary %>% group_by(True.Rate.Omega.3, cat) %>% summarise(num_reps =n())

#both cases have selection (cuz oemga_3 >1)

YesYes_correct = all_omega_vary %>% group_by(True.Rate.Omega.3, cat)  %>% filter(cat == "YesYes") %>% filter(BUSTED.SRV.P<0.05) %>% tally() 
colnames(YesYes_correct)[3] = "SRV"

YesNo_correct = all_omega_vary %>% group_by(True.Rate.Omega.3, cat)  %>% filter(cat == "YesNo") %>% filter(BUSTED.SRV.P<0.05) %>% tally() 
colnames(YesNo_correct)[3] = "SRV"
SRV_w_sel = rbind(YesNo_correct,YesYes_correct)

YesYes_correct_BUSTED = all_omega_vary %>% group_by(True.Rate.Omega.3, cat)  %>% filter(cat == "YesYes") %>% filter(BUSTED.P<0.05) %>% tally() 
colnames(YesYes_correct_BUSTED )[3] = "BUSTED"

YesNo_correct_BUSTED = all_omega_vary %>% group_by(True.Rate.Omega.3, cat)  %>% filter(cat == "YesNo") %>% filter(BUSTED.P<0.05) %>% tally() 
colnames(YesNo_correct_BUSTED )[3] = "BUSTED"

BUSTED_w_sel = rbind(YesNo_correct_BUSTED,YesYes_correct_BUSTED)

w_sel = full_join(BUSTED_w_sel, SRV_w_sel, by = c("True.Rate.Omega.3","cat")) %>% full_join(.,basic, by= c("True.Rate.Omega.3","cat")) %>% arrange(True.Rate.Omega.3)
w_sel = replace(w_sel, is.na(w_sel), 0)
w_sel[,3:4] = w_sel[,3:4]/w_sel$num_reps
w_sel 
pdf(file = "G:/BRC/Latex report/figure-latex/pwrVaryO3.pdf")
w_sel %>% melt(id.vars= c("True.Rate.Omega.3","cat","num_reps")) %>% ggplot(aes(x= True.Rate.Omega.3, y = value,color =interaction(cat,variable))) + geom_point() + geom_smooth() +labs( x = "Omega 3 Rate", y = "Power", color = "Category and Analysis") 
dev.off()


compile("G:/BRC/SimResults/Five_seq/1000_Codons/", "G:/BRC/SimResults/Five_seq/1000_Codons/BigOldData.csv")
#compile the right files

compile("G:/BRC/SimResults/Five_seq/1000_Codons/yesSel_noSRV/",
        "G:/BRC/SimResults/Five_seq/1000_Codons/yesSel_noSRV/YesNo_A6_VaryO.csv")
simulation_inputs("G:/BRC/SimResults/Five_seq/1000_Codons/yesSel_noSRV/",
                  "G:/BRC/SimResults/Five_seq/1000_Codons/yesSel_noSRV/YesNo_A6_VaryO_Truth.csv")

YesYes.A6.dat = pwr_tab("G:/BRC/SimResults/Five_seq/YesYes_A6_VaryO.csv",
                                 "G:/BRC/SimResults/Five_seq/YesYes_A6_VaryO_Truth.csv")
YesNo.A6.dat = pwr_tab("G:/BRC/SimResults/Five_seq/1000_Codons/yesSel_noSRV/YesNo_A6_VaryO.csv",
                       "G:/BRC/SimResults/Five_seq/1000_Codons/yesSel_noSRV/YesNo_A6_VaryO_Truth.csv" )
all.pwr.tab = bind_rows(YesYes.A6.dat,YesNo)

all.pwr.tab = all.pwr.tab %>% arrange(True.Rate.Omega.3) %>% select(True.Rate.Omega.3,`$BUSTED \\omega_3$ MLE`, `SRV $\\omega_3$ MLE`, True.SRV.CV, `Mean CV`, BUSTED_PWR,SRV_PWR)

all.pwr.tab %>% filter(True.Rate.Omega.3 == 1.1|True.Rate.Omega.3==1.5 | True.Rate.Omega.3 == 2 |
                         True.Rate.Omega.3 ==2.5 | True.Rate.Omega.3 == 2.9 )

kable(all.pwr.tab %>% filter(True.Rate.Omega.3 == 1.1|True.Rate.Omega.3==1.5 | True.Rate.Omega.3 == 2 |
                               True.Rate.Omega.3 ==2.5 | True.Rate.Omega.3 == 2.9 ))

NoYes.O2.dat = read.csv(
  "G:/BRC/SimResults/Five_seq/1000_Codons/noSel_yesSRV/NoYes_VaryA_O2.csv", as.is = TRUE)
NoYes.O2.Truth = read.csv(
  "G:/BRC/SimResults/Five_seq/1000_Codons/noSel_yesSRV/NoYes_VaryA_O2_Truth.csv", as.is = TRUE)

pwr_tab_alpha <- function(NoYes.O2.dat, NoYes.O2.Truth) {
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
glimpse(O2.pwr.tab)



dir= "G:/BRC/SimResults/Five_seq/"
basename = "Five_seq_all"
process_dat <- function(dir, basename){
  temp = paste(dir,basename,"_Truth.csv", sep = "")
  simulation_inputs(dir,temp)
  truth = read.csv(temp, as.is = T)
  temp = paste(dir,basename,"_results.csv", sep = "")
  compile(dir,temp)
  dat = read.csv(temp, as.is = T)
  dat = dat %>% mutate(True.CV = 1337)
  #add omega 3 for same reason
  dat = dat %>% mutate(True.Rate.Omega.3 = 8008)
  for( i in seq(from = 1, to = nrow(truth), by = 1)){
    temp = which(str_detect(dat$FILE,truth$Cat[i]))
    dat$True.CV[temp] = truth$CV.SRV[i] 
    dat$True.Rate.Omega.3[temp] = truth$Omega.3.Rate[i]
  }
  dat = mutate(dat, cat = str_extract(dat$FILE, "YesYes|YesNo|NoNo|NoYes"))
  return(dat)
}


add_truth <- function(dat, truth){
  dat = dat %>% mutate(True.CV = 1337)
  dat = dat %>% mutate(True.Rate.Omega.3 = 8008)
  for( i in seq(from = 1, to = nrow(truth), by = 1)){
    temp = which(str_detect(dat$FILE,truth$Cat[i]))
    dat$True.CV[temp] = truth$CV.SRV[i] 
    dat$True.Rate.Omega.3[temp] = truth$Omega.3.Rate[i]
  }
  return(dat)
}

seq_len_all = x %>% filter(True.Rate.Omega.3 == 3|True.Rate.Omega.3 == 1) %>% glimpse


lowCV.dat = read.csv("E:/BRC/SimResults/Five_seq/YesYes_lowCV_VaryO.csv", as.is = TRUE)
lowCV.truth =read.csv("E:/BRC/SimResults/Five_seq/YesYes_lowCV_VaryO_Truth.csv", as.is = TRUE)


lowCV.dat = add_truth("E:/BRC/SimResults/Five_seq/YesYes_lowCV_VaryO.csv,lowCV.truth")

pwr_tab("E:/BRC/SimResults/Five_seq/YesYes_lowCV_VaryO.csv", "E:/BRC/SimResults/Five_seq/YesYes_lowCV_VaryO_Truth.csv")

lowCV.dat %>% group_by()


NoNo_1000 = all.dat %>% filter(Sites == 1000, True.Rate.Omega.3 <=1, True.CV == 0) 
YesYes_1000 = all.dat %>% filter(Sites == 1000, True.Rate.Omega.3 > 1, True.CV > 0)
YesNo_1000 = all.dat %>% filter(Sites == 1000, True.Rate.Omega.3 >1, True.CV == 0)
NoYes_1000 = all.dat %>% filter(Sites == 1000, True.Rate.Omega.3 <=1, True.CV>0)

#reworking the power table function
pwr_tab <- function(dat) {
  
  
  A1.dat = dat
  A1.basic = A1.dat %>% group_by(True.Rate.Omega.3) %>% summarise(num_reps = n())
  
  A1.BUSTED = A1.dat %>% group_by(True.Rate.Omega.3,True.CV, Cat) %>% filter(BUSTED.P <
                                                                               0.05) %>% tally()
  A1.BUSTED = rename(A1.BUSTED, "BUSTED_PWR" = n)
  
  A1.SRV = A1.dat %>%  group_by(True.Rate.Omega.3, True.CV, Cat) %>% filter(BUSTED.SRV.P <
                                                                              0.05) %>%  tally()
  A1.SRV = rename(A1.SRV, "SRV_PWR" = n)
  
  A1.pwr.tab = full_join(A1.BUSTED, A1.SRV, by = c("True.Rate.Omega.3", "True.CV", "Cat")) %>% full_join(., A1.basic, by = c("True.Rate.Omega.3"))
  
  
  
  A1.means = A1.dat %>% group_by(True.Rate.Omega.3) %>%   summarise(
    "$BUSTED \\omega_3$ MLE" = mean(BUSTED.omega3.MLE),
    "SRV $\\omega_3$ MLE" = mean(BUSTED.SRV.omega3.MLE),
    "Mean CV" = mean(CV.SRV)
  )
  A1.pwr.tab = full_join(A1.pwr.tab, A1.means, by = "True.Rate.Omega.3")
  A1.pwr.tab = replace(A1.pwr.tab, is.na(A1.pwr.tab), 0)
  A1.pwr.tab$BUSTED_PWR = A1.pwr.tab$BUSTED_PWR / A1.pwr.tab$num_reps
  A1.pwr.tab$SRV_PWR = A1.pwr.tab$SRV_PWR / A1.pwr.tab$num_reps
  return(A1.pwr.tab)
}


pwr_tab(NoYes_1000)


NoYes_1000 %>% ggplot(aes(x=True.CV,y =True.Rate.Omega.3, color = Cat))+ geom_point()


NoNo_1000 %>%  ggplot(aes(x = BUSTED.treelength, y = BUSTED.SRV.treelength))+ geom_point()


cur.dir = "C:/Users/srwisots/Google Drive/Muse/BRC/SimResults/Five_seq/1000_Codons/"
tree_length <- function(cur.dir){
  jsons <- list.files(path = cur.dir,
                      pattern = '*.json', recursive = T)
  
  #create empty data.frame with variable names
  
  df = read.table(text="", col.names = c("File", "length.Con.SRV", "length.Uncon.SRV",
                                         "length.Con.BUSTED", "length.Uncon.BUSTED"))
  for (i in  seq(from=1, to=length(jsons), by=2)){
    
    
    FILE = jsons[i]
    
    
    if(grepl("SRV", jsons[i])){
      method = "SRV"
      filepath = paste(cur.dir,jsons[i], sep="")
      
      test = fromJSON(filepath)
      tree_length_Uncon = test$fits$`Unconstrained model`$`tree length`
      if (test$timers$constrained > 0){

      tree_length_Con = test$fits$`Constrained model`$`tree length`
      }
      else{
        tree_length_Con  = "NA"
      }
 
      
    }
    length.Uncon.srv = tree_length_Uncon
    length.Con.srv = tree_length_Con 
    
    if(grepl("SRV",jsons[i+1])==FALSE){
      method = "BUSTED"   
      filepath = paste(cur.dir,jsons[i+1], sep="")
      
      test = fromJSON(filepath)
      tree_length_Uncon = test$fits$`Unconstrained model`$`tree length`
      if (test$timers$constrained > 0){
        
        tree_length_Con = test$fits$`Constrained model`$`tree length`
      }
      else{
        tree_length_Con  = "NA"
      }
      
    }
    length.Uncon.BUSTED = tree_length_Uncon 
    length.Con.BUSTED = tree_length_Con 
    
   
      df[nrow(df)+1,]=c( FILE, length.Con.srv, length.Uncon.srv,length.Con.BUSTED, length.Uncon.BUSTED)
    
  }
  #df_wide = spread(df,File, length)
  return(df)
}

tree = tree_length(dir)




