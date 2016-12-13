#script that loops through folder
#creates list of json files
#pulls info out from jsons into data.frame
#for easier manipulation in R

#try to loop thru jsons to get info needed

cur.dir = "E:/Selectome/jsons/jsons/"
csv ="E:/BRC/SimResults/test1/test1_results.csv"
install.packages("jsonlite")

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


###generate CSVs for first 3 tests

compile("E:/BRC/SimResults/test1/", "E:/BRC/SimResults/test1/test1_results.csv")
compile("E:/BRC/SimResults/test2/jsons/","E:/BRC/SimResults/test2/test2_results.csv")
compile("E:/BRC/SimResults/test2/jsons/", "E:/BRC/SimResults/test3/test3_results.csv")
compile("E:/BRC/SimResults/noSRV/","E:/BRC/SimResults/noSRV/noSRV_results.csv")
compile("E:/BRC/SimResults/highOlowA/","E:/BRC/SimResults/highOlowA/highOlowA_results.csv")
compile("E:/Selectome/jsons/","E:/Selectome/Selectome_results.csv")
compile("E:/BRC/SimResults/Five_seq/noSel_noSRV/", "E:/BRC/SimResults/Five_seq/NoNo_results.csv")
compile("E:/BRC/SimResults/Five_seq/noSel_yesSRV/", "E:/BRC/SimResults/Five_seq/NoYes_results.csv")

#1000 Codons
compile("E:/BRC/SimResults/Five_seq/1000_Codons/noSel_noSRV/", "E:/BRC/SimResults/Five_seq/NoNo_1000_results.csv")

compile("E:/BRC/SimResults/Five_seq/1000_Codons/noSel_yesSRV/", "E:/BRC/SimResults/Five_seq/NoYes_1000_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_noSRV/", "E:/BRC/SimResults/Five_seq/YesNo_1000_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_yesSRV/", "E:/BRC/SimResults/Five_seq/YesYes_1000_results.csv")

simulation_inputs("E:/BRC/SimResults/Five_seq/1000_Codons/noSel_noSRV/", "E:/BRC/SimResults/Five_seq/NoNo_1000_Truth.csv")
simulation_inputs("E:/BRC/SimResults/Five_seq/1000_Codons/noSel_yesSRV/", "E:/BRC/SimResults/Five_seq/NoYes_1000_Truth.csv")
simulation_inputs("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_noSRV/", "E:/BRC/SimResults/Five_seq/YesNo_1000_Truth.csv")
simulation_inputs("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_yesSRV/", "E:/BRC/SimResults/Five_seq/YesYes_1000_Truth.csv")



compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_yesSRV/O3_1_1/", "E:/BRC/SimResults/Five_seq/YesYes_1000_O3_1_1_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_yesSRV/O3_1_2/", "E:/BRC/SimResults/Five_seq/YesYes_1000_O3_1_2_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_yesSRV/O3_1_3/", "E:/BRC/SimResults/Five_seq/YesYes_1000_O3_1_3_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_yesSRV/O3_1_4/", "E:/BRC/SimResults/Five_seq/YesYes_1000_O3_1_4_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_yesSRV/O3_1_5/", "E:/BRC/SimResults/Five_seq/YesYes_1000_O3_1_5_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_yesSRV/O3_1_6/", "E:/BRC/SimResults/Five_seq/YesYes_1000_O3_1_6_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_yesSRV/O3_1_7/", "E:/BRC/SimResults/Five_seq/YesYes_1000_O3_1_7_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_yesSRV/O3_1_8/", "E:/BRC/SimResults/Five_seq/YesYes_1000_O3_1_8_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_yesSRV/O3_1_9/", "E:/BRC/SimResults/Five_seq/YesYes_1000_O3_1_9_results.csv")

compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_noSRV/O3_1_1/", "E:/BRC/SimResults/Five_seq/YesNo_1000_O3_1_1_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_noSRV/O3_1_2/", "E:/BRC/SimResults/Five_seq/YesNo_1000_O3_1_2_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_noSRV/O3_1_3/", "E:/BRC/SimResults/Five_seq/YesNo_1000_O3_1_3_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_noSRV/O3_1_4/", "E:/BRC/SimResults/Five_seq/YesNo_1000_O3_1_4_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_noSRV/O3_1_5/", "E:/BRC/SimResults/Five_seq/YesNo_1000_O3_1_5_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_noSRV/O3_1_6/", "E:/BRC/SimResults/Five_seq/YesNo_1000_O3_1_6_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_noSRV/O3_1_7/", "E:/BRC/SimResults/Five_seq/YesNo_1000_O3_1_7_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_noSRV/O3_1_8/", "E:/BRC/SimResults/Five_seq/YesNo_1000_O3_1_8_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_noSRV/O3_1_9/", "E:/BRC/SimResults/Five_seq/YesNo_1000_O3_1_9_results.csv")

compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_noSRV/O3_2/", "E:/BRC/SimResults/Five_seq/YesNo_1000_O3_2_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_noSRV/O3_2_1/", "E:/BRC/SimResults/Five_seq/YesNo_1000_O3_2_1_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_noSRV/O3_2_2/", "E:/BRC/SimResults/Five_seq/YesNo_1000_O3_2_2_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_noSRV/O3_2_3/", "E:/BRC/SimResults/Five_seq/YesNo_1000_O3_2_3_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_noSRV/O3_2_4/", "E:/BRC/SimResults/Five_seq/YesNo_1000_O3_2_4_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_noSRV/O3_2_5/", "E:/BRC/SimResults/Five_seq/YesNo_1000_O3_2_5_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_noSRV/O3_2_6/", "E:/BRC/SimResults/Five_seq/YesNo_1000_O3_2_6_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_noSRV/O3_2_7/", "E:/BRC/SimResults/Five_seq/YesNo_1000_O3_2_7_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_noSRV/O3_2_8/", "E:/BRC/SimResults/Five_seq/YesNo_1000_O3_2_8_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_noSRV/O3_2_9/", "E:/BRC/SimResults/Five_seq/YesNo_1000_O3_2_9_results.csv")




compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_yesSRV/A3_1/", "E:/BRC/SimResults/Five_seq/YesYes_1000_A3_1_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_yesSRV/A3_1_1/", "E:/BRC/SimResults/Five_seq/YesYes_1000_A3_1_1_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_yesSRV/A3_1_2/", "E:/BRC/SimResults/Five_seq/YesYes_1000_A3_1_2_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_yesSRV/A3_1_3/", "E:/BRC/SimResults/Five_seq/YesYes_1000_A3_1_3_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_yesSRV/A3_1_4/", "E:/BRC/SimResults/Five_seq/YesYes_1000_A3_1_4_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_yesSRV/A3_1_5/", "E:/BRC/SimResults/Five_seq/YesYes_1000_A3_1_5_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_yesSRV/A3_1_6/", "E:/BRC/SimResults/Five_seq/YesYes_1000_A3_1_6_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_yesSRV/A3_1_7/", "E:/BRC/SimResults/Five_seq/YesYes_1000_A3_1_7_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_yesSRV/A3_1_8/", "E:/BRC/SimResults/Five_seq/YesYes_1000_A3_1_8_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_yesSRV/A3_1_9/", "E:/BRC/SimResults/Five_seq/YesYes_1000_A3_1_9_results.csv")


compile("E:/BRC/SimResults/Five_seq/1000_Codons/NoSel_YesSRV/A3_1/", "E:/BRC/SimResults/Five_seq/NoYes_1000_A3_1_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/NoSel_YesSRV/A3_1_1/", "E:/BRC/SimResults/Five_seq/NoYes_1000_A3_1_1_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/NoSel_YesSRV/A3_1_2/", "E:/BRC/SimResults/Five_seq/NoYes_1000_A3_1_2_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/NoSel_YesSRV/A3_1_3/", "E:/BRC/SimResults/Five_seq/NoYes_1000_A3_1_3_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/NoSel_YesSRV/A3_1_4/", "E:/BRC/SimResults/Five_seq/NoYes_1000_A3_1_4_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/NoSel_YesSRV/A3_1_5/", "E:/BRC/SimResults/Five_seq/NoYes_1000_A3_1_5_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/NoSel_YesSRV/A3_1_6/", "E:/BRC/SimResults/Five_seq/NoYes_1000_A3_1_6_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/NoSel_YesSRV/A3_1_7/", "E:/BRC/SimResults/Five_seq/NoYes_1000_A3_1_7_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/NoSel_YesSRV/A3_1_8/", "E:/BRC/SimResults/Five_seq/NoYes_1000_A3_1_8_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/NoSel_YesSRV/A3_1_9/", "E:/BRC/SimResults/Five_seq/NoYes_1000_A3_1_9_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/NoSel_YesSRV/A3_2/", "E:/BRC/SimResults/Five_seq/NoYes_1000_A3_2_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/NoSel_YesSRV/A3_2_5/", "E:/BRC/SimResults/Five_seq/NoYes_1000_A3_2_5_results.csv")
compile("E:/BRC/SimResults/Five_seq/1000_Codons/NoSel_YesSRV/A3_3/", "E:/BRC/SimResults/Five_seq/NoYes_1000_A3_3_results.csv")



#6000 Codons
compile("E:/BRC/SimResults/Five_seq/5000_Codons/noSel_noSRV/", "E:/BRC/SimResults/Five_seq/NoNo_6000_results.csv")
simulation_inputs("E:/BRC/SimResults/Five_seq/5000_Codons/noSel_noSRV/", "E:/BRC/SimResults/Five_seq/NoNo_6000_Truth.csv")
compile("E:/BRC/SimResults/Five_seq/5000_Codons/noSel_yesSRV/", "E:/BRC/SimResults/Five_seq/NoYes_6000_results.csv")
simulation_inputs("E:/BRC/SimResults/Five_seq/5000_Codons/noSel_yesSRV/", "E:/BRC/SimResults/Five_seq/NoYes_6000_Truth.csv")
compile("E:/BRC/SimResults/Five_seq/5000_Codons/yesSel_noSRV/", "E:/BRC/SimResults/Five_seq/YesNo_6000_results.csv")
simulation_inputs("E:/BRC/SimResults/Five_seq/5000_Codons/yesSel_noSRV/", "E:/BRC/SimResults/Five_seq/YesNo_6000_Truth.csv")
compile("E:/BRC/SimResults/Five_seq/5000_Codons/yesSel_yesSRV/", "E:/BRC/SimResults/Five_seq/NoYes_6000_results.csv")
simulation_inputs("E:/BRC/SimResults/Five_seq/5000_Codons/yesSel_yesSRV/", "E:/BRC/SimResults/Five_seq/NoYes_6000_Truth.csv")

#Lower CV with varying Omega (Alpha 3 = 2)
compile("E:/BRC/SimResults/Five_seq/1000_Codons/VaryO_A2/", "E:/BRC/SimResults/Five_seq/YesYes_A2_VaryO.csv")
simulation_inputs("E:/BRC/SimResults/Five_seq/1000_Codons/VaryO_A2/", "E:/BRC/SimResults/Five_seq/YesYes_A2_VaryO_Truth.csv")

compile("E:/BRC/SimResults/Five_seq/1000_Codons/VaryO_A1/", "E:/BRC/SimResults/Five_seq/YesYes_A1_VaryO.csv")
simulation_inputs("E:/BRC/SimResults/Five_seq/1000_Codons/VaryO_A1/","E:/BRC/SimResults/Five_seq/YesYes_A1_VaryO_Truth.csv" )

compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_yesSRV/VaryO_A6/","E:/BRC/SimResults/Five_seq/YesYes_A6_VaryO.csv" )

compile("E:/BRC/SimResults/Five_seq/1000_Codons/VaryO_lowCV/", "E:/BRC/SimResults/Five_seq/YesYes_lowCV_VaryO.csv")
simulation_inputs("E:/BRC/SimResults/Five_seq/1000_Codons/VaryO_lowCV/","E:/BRC/SimResults/Five_seq/YesYes_lowCV_VaryO_Truth.csv" )

compile("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_yesSRV/VaryA_O2/", 
        "E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_yesSRV/YesYes_VaryA_O2.csv")
simulation_inputs("E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_yesSRV/VaryA_O2/", 
                  "E:/BRC/SimResults/Five_seq/1000_Codons/yesSel_yesSRV/YesYes_VaryA_O2_Truth.csv")

compile("E:/BRC/SimResults/Five_seq/1000_Codons/noSel_yesSRV/", 
        "E:/BRC/SimResults/Five_seq/1000_Codons/noSel_yesSRV/NoYes_VaryA_O2.csv")
simulation_inputs("E:/BRC/SimResults/Five_seq/1000_Codons/noSel_yesSRV/", 
                  "E:/BRC/SimResults/Five_seq/1000_Codons/noSel_yesSRV/NoYes_VaryA_O2_Truth.csv")
