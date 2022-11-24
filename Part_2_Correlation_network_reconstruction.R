library(siggitRausti)
library(WGCNA)
library(impute)
library(parallel)
library(readxl)
library(BiocManager)
library(GO.db)


# Perform Network analysis - extract module eigengenes and perform 
# correlation analysis. 

dat <- read.csv('All_data_normalised_24NOV2022.csv')

clin_dat = read.csv('Database.csv',sep = ';')

# Take out samples that do not have clinical information:
clinical_patients = clin_dat$Subject.Identifier.for.the.Study

id_vec = c()
for (patient in 1:length(clinical_patients)){
  tmp_patient = clinical_patients[patient]
  if (any(dat$X %in% tmp_patient)){
    id_vec = append(id_vec,which(dat$X %in% tmp_patient))
  } else {
   print(paste('MISSING VALUE FOR',tmp_patient)) 
  }
}

# Keep only samples with clinical information:
dat = dat[id_vec,]

# Take out patient ids from the dataframe into a separate file:
patients = dat$X
FinalData = dat[-1]

# Scan through the features and remove "bad features" (This is a pre-processing step that 
# is regularly performed in WGCNA. I have already done some pre-processing so this should 
# not remove any features... Keeping this for consistency if code will be used for other projects)

set.seed(1234)
gsg = goodSamplesGenes(FinalData, verbose = 3)
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(FinalData)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(FinalData)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  FinalData = FinalData[gsg$goodSamples, gsg$goodGenes]
}


# Define a threshold to use for the network correlation analysis (which power value of the 
# correlations will result in the most scale-free topology of the correlation network):
sampleTree = hclust(dist(FinalData,method = 'euclidean'),method = 'complete')
plot(sampleTree, main = "Sample clustering of trauma patients", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(FinalData, powerVector = powers, verbose = 5,networkType = 'unsigned',corFnc = 'bicor')

# Plot the results:
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
abline(h=0.9,col="red")


# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# Set the main correlation function (cor) as the WGCNA correlation function:
cor = WGCNA::cor

# Lets have the value at 5, since the scale-free topology kind of converges there for the next few values:
net = blockwiseModules(FinalData, power = 5,networkType = 'unsigned',TOMType = 'signed Nowick',corFnc = 'bicor',
                       deepSplit = 4,maxBlockSize = 8000)


# Convert labels to colors for plotting
mergedColors = net$colors

# Plot the dendrogram and the module colors underneath (for a single block)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],"Module colors",
                    dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)

# Calculate eigenfeatures
MEList = moduleEigengenes(FinalData, colors = mergedColors)
MEs = MEList$eigengenes


library(Hmisc)
res2 <- rcorr(as.matrix(MEs),type = 'spearman')
MEDiss = 1-res2$r

# Cluster module eigengenes
METree = hclust(as.dist(MEDiss));

# Plot the result
plot(METree, main = "Clustering of module eigenmetabolites",
     xlab = "", sub = "")

# Select a distance threshold (how similar modules have to be to be merged). 
# This is arbitrary, and can be changed. 

MEDissThres = 0.5

# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")

# Call an automatic merging function
merge = mergeCloseModules(FinalData, mergedColors, cutHeight = MEDissThres, verbose = 3)

# The merged module colors
mergedColors2 = merge$colors;

# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

plotDendroAndColors(net$dendrograms[[1]], cbind(mergedColors[net$blockGenes[[1]]], mergedColors2[net$blockGenes[[1]]]),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

MEList = moduleEigengenes(FinalData, colors = mergedColors2)
MEs = MEList$eigengenes
res2 <- rcorr(as.matrix(MEs),type = 'spearman')
MEDiss = 1-res2$r

# Cluster merged module eigenfeatures
METree = hclust(as.dist(MEDiss), method = "average");

# Plot them
plot(METree, main = "Clustering of module eigenmetabolites",
     xlab = "", sub = "")

# Save entire workspace:
save.image("final_data_workspace_24NOV2022.RData")
