#Required functions to process results
#and loading libraries

library(jsonlite)
library(ggplot2)
library(knitr)
library(stringr)
library(reshape2)
library(dplyr)

#Right now these are just he BUSTED-SRV branch lengths
branch_length <- function(cur.dir){
  jsons <- list.files(path = cur.dir,
                      pattern = '*.json')
  
  #create empty data.frame with variable names
  
  df = read.table(text="", col.names = c("Branch", "File", "length.SRV","length.BUSTED"))
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