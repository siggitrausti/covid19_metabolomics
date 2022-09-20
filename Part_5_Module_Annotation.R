# PART 5: Run the mummichog algorithm on SELECTED MODULES. Note that this originally 
# took all modules and performed automatic annotations on them all. But after 
# installing MetaboAnalyst on my computer, the package for MetaboAnalyst had not 
# been updated, so the annotation procedure could not take place in R. As a temporary 
# solution, I just picked a module of interest, created the spreadsheet that is 
# required by Mummichog and went to the Functional Analysis on www.metaboanalyst.ca
# and ran it there for my modules of interest (obtained from the modelling procedure
# performed in the next step in Jupyter). 

# NOTE: I have not made this into a function yet. This can be done in either Python
# or R, depending on the platform of interest. 

# Load workspace from part 1:
load('final_data_workspace_15SEP2022_PART2.Rdata')

# Set the module of interest:
module_of_interest = 'red'


signed_network = F # Set to T if the network uses signed correlation. 
name_to_check <- colnames(FinalData)
origin_of_dat <- substring(name_to_check,1,3)
masses_to_check <- as.numeric(substring(name_to_check,8))
o = which(colnames(datKME2) %in% paste0('MM.',module_of_interest))

# Use the empirical distribution function to generate p-values:
if (signed_network){
  t.scores <- datKME2[,o]
  cdf_fun <- ecdf(datKME2[,o])
} else {
  t.scores <- abs(datKME2[,o])
  cdf_fun <- ecdf(abs(datKME2[,o]))
}
p.values <- 1-cdf_fun(t.scores)

# ID the mode:
outp_dat1_mode <- c()
for (i in 1:length(name_to_check)){
  type1 = origin_of_dat[i]
  if (type1 == 'neg'){
    outp_dat1_mode <- c(outp_dat1_mode,'negative')
  } else if (type1 == 'pos'){
    outp_dat1_mode <- c(outp_dat1_mode,'positive')
  } else if (type1 == 'bas'){
    outp_dat1_mode <- c(outp_dat1_mode,'negative')
  } else {
    outp_dat1_mode <- c(outp_dat1_mode,NA)
  }
}

outp_dataf <- data.frame(cbind(m.z = masses_to_check,p.value = p.values,t.score = t.scores,mode = outp_dat1_mode))
outp_dataf$t.score <- as.numeric(as.character(outp_dataf$t.score))
outp_dataf$p.value <- as.numeric(as.character(outp_dataf$p.value))
outp_dataf <- outp_dataf[complete.cases(outp_dataf),]
outp_dataf <- outp_dataf[order(outp_dataf$t.score,decreasing = T),]

# Write the spreadsheet as a .csv to import into MetaboAnalyst online:
write.csv(outp_dataf,paste0(module_of_interest,'.csv'), row.names = FALSE)
