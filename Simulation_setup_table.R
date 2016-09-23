#make a giant table showing all the simulation settings




dir = "G:/BRC/SimResults/Five_seq" #sets directory to look thru

list = list.files(path = "E:/BRC/SimResults/Five_seq", recursive = T, pattern ="^([^.]+)$") #lists only files with out an extension


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
  
  #return(setup.tab)
  write.csv(file = csv, x = setup.tab, row.names= F)
}

write.csv()


#reorder columns
setup.tab =select(setup.tab, Sites, Cat, Omega.1.Rate, Omega.1.Weight, Omega.2.Rate, Omega.2.Weight,Omega.3.Rate, Omega.3.Weight,
       Alpha.1.Rate, Alpha.1.Weight, Alpha.2.Rate, Alpha.2.Weight,Alpha.3.Rate, Alpha.3.Weight)


setup.tab$Cat = str_extract(setup.tab$Cat,"NoNo|YesYes|YesNo|NoYes") # make cats only YesYes, NoNo, etc



setup.tab[,3:14] = setup.tab[,3:14]  %>%  unlist() %>% as.numeric() #change strings into numeric for rates and weights


#make a table for only the varying sequence lengths
seq.len.tab=setup.tab %>% filter(Omega.3.Rate==3|Omega.3.Rate==7|Omega.3.Rate==1)
#input latex to format table and group by sites
seq.len.tab[,1] = paste0("\\multirow{4}{*}{", seq.len.tab[,1], "}")
#change column names so the grouping with multicol make sense
colnames(seq.len.tab)[3:14]=str_extract(colnames(setup.tab)[3:14],"Rate|Weight")
#make an xtable with a caption
seq.len.xtab = xtable(seq.len.tab, caption = "Set up values for sequence length comparison")

#align table and set fixed col width
align(seq.len.xtab) <- c( 'l', 'l',rep('p{1.25cm}',13))

#set multicol names
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <-paste0(paste0('& & \\multicolumn{2}{c}{',
                                 c("$\\omega_1$","$\\omega_2$","$\\omega_3$",
                                   "$\\alpha_1$", "$\\alpha_2$",
                                   "$\\alpha_3$"), '}', collapse=''),'\\\\')
print(xtab,  
      include.rownames = FALSE, add.to.row = addtorow, scalebox = 0.7, sanitize.text.function = force)
#print table
print(seq.len.xtab,
      include.rownames = FALSE, add.to.row = addtorow, scalebox = 0.7, sanitize.text.function = force) #force necessary for multirows to work



