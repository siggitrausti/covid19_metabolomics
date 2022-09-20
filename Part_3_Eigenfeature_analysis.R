library(siggitRausti)
library(WGCNA)
library(readxl)
library(ggplot2)

# Load workspace from part 2:
load('final_data_workspace_15SEP2022.Rdata')

# Take out grey module (non-assigned metabolites):
MEs_orig <- MEs
MEs <- MEs[,-which(colnames(MEs) %in% 'MEgrey')]

# Pair clinical data so that it matches the module eigengene matrix:
id_vec = c()
for (patient in 1:length(patients)){
  tmp_patient = patients[patient]
  if (any(clin_dat$Subject.Identifier.for.the.Study %in% tmp_patient)){
    id_vec = append(id_vec,which(clin_dat$Subject.Identifier.for.the.Study %in% tmp_patient))
  }
}

# Test if matches:
siggitRausti::test_match_order(clin_dat$Subject.Identifier.for.the.Study[id_vec],patients)

if (length(id_vec) != 0){
  clin_dat = clin_dat[id_vec,]
}

# Remove clinical variables with more than 20% missing values
na_no <- rep(NA,ncol(clin_dat))
for (i in 1:ncol(clin_dat)){
  na_no[i] <- length(which(is.na(clin_dat[,i])))
}
id_remove2 <- which(na_no > 0.2*nrow(clin_dat)) 
clin_dat <- clin_dat[,-id_remove2]
clin_dat <- data.frame(clin_dat)


################################################################################
# PART 1: Trait vs module matrix. Shows correlation and associated significance of traits and module eigenfeature values:
traits <- clin_dat[,c(2,3,5,6,10,13,14,15,17,20,26,27,32)] # Traits can be manually altered according to what 
# one wants to show...

nSamples = nrow(MEs)
moduleTraitCor = WGCNA::cor(MEs, traits, use = "p",method = 'spearman')
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples) 
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "")

png("traits_vs_modules_15SEP2022.png",width = 1200,height = 2000,units = 'px',res = 150)
par(mar = c(7, 7, 1, 1));
labeledHeatmap(Matrix = moduleTraitCor,xLabels = names(traits),
               yLabels = names(MEs),ySymbols = names(MEs),colorLabels = FALSE,
               colors = blueWhiteRed(10),textMatrix = textMatrix,
               setStdMargins = FALSE,cex.text = 0.5,zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

# Save module membership values:
datKME=signedKME(FinalData, MEs_orig, outputColumnName="MM.")
datKME2 <- datKME[,which(substring(colnames(datKME),4) %in% substring(colnames(MEs),3))]

# Save workspace:
save.image("final_data_workspace_15SEP2022_PART2.RData")

################################################################################
# Part 2: Save all datasets to be used in Jupyter for Random Forest modeling. 
write.csv(MEs,"MEs_dataframe_15SEP2022.csv", row.names = FALSE)
write.csv(patients,"patients_vector_15SEP2022.csv", row.names = FALSE)
write.csv(clin_dat,"clinical_dataframe_15SEP2022.csv", row.names = FALSE)
write.csv(FinalData,"Final_dataframe_15SEP2022.csv", row.names = FALSE)
write.csv(mergedColors2,"module_assignments_15SEP2022.csv", row.names = FALSE)
write.csv(datKME2,"datkme_15SEP2022.csv", row.names = FALSE)

################################################################################

