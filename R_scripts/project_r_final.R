#Installing and importing libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("minfi")
library(minfi)
options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages("rstudioapi")  # library that recover the path where the R project is
library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))#seth the active path where the R project is, automatically
library(knitr)
library(IlluminaHumanMethylation450kmanifest)
##Load raw data with minfi and create an object called RGset storing the RGChannelSet object  
rm(list=ls())

#set directory of input data
baseDir <- (file.path(getwd(), "Input_Data"))
targets <- read.metharray.sheet(baseDir)

#1. Create an object of class RGChannelSet using the function read.metharray.ex
RGset <- read.metharray.exp(targets = targets)
save(RGset,file="RGset.RData")

# 2. Create the dataframes Red and Green to store the red and green fluorescences respectively.
Red <- data.frame(getRed(RGset))
#Show the dimension, the internal structure and the first six row of Red object
dim(Red) 

#Repeat the process for the Green dataset
Green <- data.frame(getGreen(RGset))
dim(Green)

# 3. Fill the following table: what are the Red and Green fluorescences for the address assigned to your group? Optional: check in the manifest file if the address corresponds to a Type I or a Type II probe and, in case of Type I probe, report its color.
#Retrieve Red and Green fluorescence for the group address
address<-"13673406"
probe_red<-Red[rownames(Red)==address,]
probe_green<-Green[rownames(Green)==address,]

# Load the cleaned manifest file to check the type of the probe
load("Illumina450Manifest_clean.RData")
paste("Type of address A: ",Illumina450Manifest_clean[Illumina450Manifest_clean$AddressA_ID==address,'Infinium_Design_Type'])
paste("Type of address B: ",Illumina450Manifest_clean[Illumina450Manifest_clean$AddressB_ID==address,'Infinium_Design_Type'])
# Load SampleSheet for later
SampleSheet <- read.csv("Input_Data\\SampleSheet_Report_II.csv",header=T)

#in our case is type II, so the probe can be either Red or Green. Infact AddressB is empty, and in AddressA is Type II.
#extract all data related to sample, red fluorescence, green fluorescence and type and color to fill the table
df_address <- data.frame(Sample = colnames(probe_green),'Number of probes' = unlist(probe_red, use.names = FALSE) + unlist(probe_green, use.names = FALSE),'Type' = Illumina450Manifest_clean[Illumina450Manifest_clean$AddressA_ID == address, 'Infinium_Design_Type'])
knitr::kable(df_address, col.names=c("Sample","Number or probes", "Infinium design type"))

#4. Create the object MSet.raw

# Use the preprocessRaw function to create the MSet.raw object from the RGChannelSet.
MSet.raw <- preprocessRaw(RGset)
save(MSet.raw,file="MSet_raw.RData")
dim(MSet.raw)

#5. Perform the following quality checks and provide a brief comment to each step: • QCplot • check the intensity of negative controls using minfi • calculate detection pValues; for each sample, how many probes have a detection p-value higher than the threshold assigned to each group?

#a. Generating QC plot
qc <- getQC(MSet.raw)
qc
plotQC(qc)
#b. Negative control intensity check
#only two samples are above the threshold line. However, the majority of the samples are distributed close to the line, and more importantly are all grouped together, except for two samples.
controlStripPlot(RGset, controls="NEGATIVE")

#c. Count Probes Exceeding Threshold 
threshold<-0.05
detP <- detectionP(RGset)
exceed <- detP>0.05
num_exceed = colSums(exceed)
df_exceed <- data.frame(Sample=colnames(probe_green), Exceed_Positions=num_exceed,perc_exceed_probes=(means_of_columns <- colMeans(exceed))*100, row.names = NULL)
knitr::kable(df_exceed, col.names=c("Sample","Number of exceed probes", "Percentage of exceed probes"))

# 6. Calculate raw beta and M values and plot the densities of mean methylation values, dividing the samples in CTRL and DIS (suggestion: subset the beta and M values matrices in order to retain CTRL or DIS subjects and apply the function mean to the 2 subsets). Do you see any difference between the two groups?

beta <- getBeta(MSet.raw)
M <- getM(MSet.raw)

ctrl_raw<- basename(targets[targets$Group == "CTRL", "Basename"])
dis_raw<- basename(targets[targets$Group == "DIS", "Basename"])
ctrlSet <- MSet.raw[, colnames(MSet.raw) %in% ctrl_raw]
disSet <- MSet.raw[, colnames(MSet.raw) %in% dis_raw]
ctrlBeta <- getBeta(ctrlSet)
ctrlM <- getM(ctrlSet)
disBeta <- getBeta(disSet)
disM <- getM(disSet)
mean_ctrlBeta <- apply(ctrlBeta, MARGIN=1, mean, na.rm=TRUE)
mean_disBeta <- apply(disBeta, MARGIN=1, mean, na.rm=TRUE)
mean_ctrlM <- apply(ctrlM, MARGIN=1, mean, na.rm=TRUE)
mean_disM <- apply(disM, MARGIN=1, mean, na.rm=TRUE)
d_mean_ctrlBeta <- density(mean_ctrlBeta, na.rm=TRUE)
d_mean_disBeta <- density(mean_disBeta, na.rm=TRUE)
d_mean_ctrlM <- density(mean_ctrlM, na.rm=TRUE)
d_mean_disM <- density(mean_disM, na.rm=TRUE)

#6b Plot the density obtained in point 6a
# Plot density of mean Beta and M values
par(mfrow=c(1,2))  # Set up the plotting area to have 1 row and 2 columns
# Plot density of mean Beta values
plot(d_mean_ctrlBeta, main="Density of Beta Values", col="#0067E6", lwd=2.5, xlab="Beta Values", ylab="Density")
lines(d_mean_disBeta, col="#E50068", lwd=2.5)
legend('topright', legend=c("CTRL", "DIS"), fill=c("#0067E6", "#E50068"), cex=0.5)
# Plot density of mean M values
plot(d_mean_ctrlM, main="Density of M Values", col="#0067E6", lwd=2.5, xlab="M Values", ylab="Density")
lines(d_mean_disM, col="#E50068", lwd=2.5)
legend('topright', legend=c("CTRL", "DIS"), fill=c("#0067E6", "#E50068"), cex=0.5)
# Reset plotting layout for subsequent plots
par(mfrow=c(1,1))

#7. Normalize the data using the function assigned to each group and compare raw data and normalized data. Produce a plot with 6 panels in which, for both raw and normalized data, you show the density plots of beta mean values according to the chemistry of the probes, the density plot of beta standard deviation values according to the chemistry of the probes and the boxplot of beta values. Provide a short comment about the changes you observe. Optional: do you think that the normalization approach that you used is appropriate considering this specific dataset? Try to color the boxplots according to the group (WT and MUT) and check whether the distribution of methylation values is different between the two groups, before and after normalization

# Subset the manifest into Type I (dfI) and Type II (dfII) probes based on Infinium_Design_Type
dfI <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type == "I", ]
dfI <- droplevels(dfI)
dfII <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type == "II", ]
dfII <- droplevels(dfII)

# Normalize using preprocessNoob
preprocessNoob_results <- preprocessNoob(RGset,dyeMethod="single")

# Separation of type I and type II probes.
beta_I <- beta[rownames(beta) %in% dfI$IlmnID, ]
beta_II <- beta[rownames(beta) %in% dfII$IlmnID, ]

# Calculate raw mean and sd
mean_of_beta_I <- apply(beta_I, 1, mean, na.rm = TRUE)
mean_of_beta_II <- apply(beta_II, 1, mean, na.rm = TRUE)
d_mean_of_beta_I <- density(mean_of_beta_I, na.rm = TRUE)
d_mean_of_beta_II <- density(mean_of_beta_II, na.rm = TRUE)
sd_of_beta_I <- apply(beta_I, 1, sd, na.rm = TRUE)
sd_of_beta_II <- apply(beta_II, 1, sd, na.rm = TRUE)
d_sd_of_beta_I <- density(sd_of_beta_I, na.rm = TRUE)
d_sd_of_beta_II <- density(sd_of_beta_II, na.rm = TRUE)

# Extract normalized beta values
beta_preprocessNoob <- getBeta(preprocessNoob_results)
beta_preprocessNoob_I <- beta_preprocessNoob[rownames(beta_preprocessNoob) %in% dfI$IlmnID, ]
beta_preprocessNoob_II <- beta_preprocessNoob[rownames(beta_preprocessNoob) %in% dfII$IlmnID, ]

# Calculate normalized mean and sd
mean_of_beta_preprocessNoob_I <- apply(beta_preprocessNoob_I, 1, mean, na.rm = TRUE)
mean_of_beta_preprocessNoob_II <- apply(beta_preprocessNoob_II, 1, mean, na.rm = TRUE)
d_mean_of_beta_preprocessNoob_I <- density(mean_of_beta_preprocessNoob_I, na.rm = TRUE)
d_mean_of_beta_preprocessNoob_II <- density(mean_of_beta_preprocessNoob_II, na.rm = TRUE)
sd_of_beta_preprocessNoob_I <- apply(beta_preprocessNoob_I, 1, sd, na.rm = TRUE)
sd_of_beta_preprocessNoob_II <- apply(beta_preprocessNoob_II, 1, sd, na.rm = TRUE)
d_sd_of_beta_preprocessNoob_I <- density(sd_of_beta_preprocessNoob_I, na.rm = TRUE)
d_sd_of_beta_preprocessNoob_II <- density(sd_of_beta_preprocessNoob_II, na.rm = TRUE)

# Set up the layout for the plots: 2 rows, 3 columns
par(mfrow = c(2, 3))

# Top Row: Raw Data
# Plot densities of mean beta values (raw)
plot(d_mean_of_beta_I, col = "blue", main = "Raw Beta Mean", xlim = c(0, 1), ylim = c(0, 5))
lines(d_mean_of_beta_II, col = "red")
legend("top", legend = c("Type I", "Type II"), col = c("blue", "red"), lty = 1, cex = 0.8)

# Plot densities of standard deviation beta values (raw)
plot(d_sd_of_beta_I, col = "blue", main = "Raw Beta SD", xlim = c(0, 0.6), ylim = c(0, 60))
lines(d_sd_of_beta_II, col = "red")
legend("topright", legend = c("Type I", "Type II"), col = c("blue", "red"), lty = 1, cex = 0.8)

# Boxplot of raw beta values
targets$Group <- as.factor(targets$Group)
palette(c("#E50068", "#0067E6"))
boxplot(beta, col = targets$Group)
title('Boxplot of raw Beta values')

# Bottom Row: Normalized Data
# Plot densities of mean beta values (preprocessNoob)
plot(d_mean_of_beta_preprocessNoob_I, col = "blue", main = "preprocessNoob Beta Mean", xlim = c(0, 1), ylim = c(0, 5))
lines(d_mean_of_beta_preprocessNoob_II, col = "red")
legend("top", legend = c("Type I", "Type II"), col = c("blue", "red"), lty = 1, cex = 0.8)

# Plot densities of standard deviation beta values (preprocessNoob)
plot(d_sd_of_beta_preprocessNoob_I, col = "blue", main = "preprocessNoob Beta SD", xlim = c(0, 0.6), ylim = c(0, 60))
lines(d_sd_of_beta_preprocessNoob_II, col = "red")
legend("topright", legend = c("Type I", "Type II"), col = c("blue", "red"), lty = 1, cex = 0.8)

# Boxplot of normalized beta values
boxplot(beta_preprocessNoob, col = targets$Group)
title('Boxplot of normalized Beta values')

# Reset par settings to default (1 row, 1 column) to avoid affecting subsequent plots
par(mfrow = c(1, 1))

# 8. Perform a PCA on the matrix of normalized beta values generated in step 7. Comment the plot (Do you see any outlier? Do the samples divide according to the group? Do they divide according to the sex of the samples? Do they divide according to the batch, that is the column Sentrix_ID?).

pca_results <- prcomp(t(beta_preprocessNoob), scale = TRUE)
install.packages("factoextra")
library(factoextra)
fviz_eig(pca_results, addlabels = TRUE, xlab = 'PC number', ylab ='% of variance', barfill = "#0063A6", barcolor = "black")

# PCA plot per Group
targets$Group <- as.factor(targets$Group)
palette(c("#E50068", "#0067E6"))
plot(pca_results$x[, 1], pca_results$x[, 2], cex = 1, pch = 19, col = targets$Group, xlab = "PC1", ylab = "PC2", xlim = c(-700, 700), ylim = c(-700, 700), main = 'PCA (Groups)')
text(pca_results$x[, 1], pca_results$x[, 2], labels = rownames(pca_results$x), cex = 0.4, pos = 3)
legend("bottomright", legend = levels(targets$Group), col = c(1, 2), pch = 19, cex = 1.0)

# PCA plot per Sex
targets$Sex <- as.factor(targets$Sex)
palette(c("#C364CA", "#8BC4F9"))
plot(pca_results$x[, 1], pca_results$x[, 2], cex = 1, pch = 19, col = targets$Sex, xlab = "PC1", ylab = "PC2", xlim = c(-700, 700), ylim = c(-700, 700), main = 'PCA (Sex)')
text(pca_results$x[, 1], pca_results$x[, 2], labels = rownames(pca_results$x), cex = 0.4, pos = 3)
legend("bottomright", legend = levels(targets$Sex), col = c(1:nlevels(targets$Sex)), pch = 19, cex = 1.0)

# PCA plot colored by batch
targets$Slide <- as.factor(targets$Slide)
palette(rainbow(length(levels(targets$Slide)))) 
plot(pca_results$x[, 1], pca_results$x[, 2], 
     cex = 1, pch = 19, col = as.numeric(targets$Slide),
     xlab = "PC1", ylab = "PC2", 
     xlim = c(-700, 700), ylim = c(-700, 700),
     main = "PCA (Batch)")
text(pca_results$x[, 1], pca_results$x[, 2], 
     labels = rownames(pca_results$x), cex = 0.4, pos = 3)
legend("bottomright", 
       legend = levels(targets$Slide), 
       col = 1:length(levels(targets$Slide)), 
       pch = 19, cex = 1.0)

#9. Using the matrix of normalized beta values generated in step 7, identify differentially methylated probes between CTRL and DIS groups using the function assigned to each group.

pheno = SampleSheet
MannWhitney_function <- function(x) {
  u_test <- wilcox.test(x ~ pheno$Group, exact=FALSE, correct = FALSE)
  return(u_test$p.value)
}

#execute the function
pValues_u_test <- apply(beta_preprocessNoob,1, MannWhitney_function)
final_u_test <- data.frame(beta_preprocessNoob, pValues_u_test)

# Count the number of probes with a p-value <= 0.05
final_utest_th <- final_u_test[final_u_test$pValues_u_test <= 0.05, ]
hist(final_u_test$pValues_u_test, main = "P-value distribution (u-test)", xlab ='p-val')
abline(v = 0.05, col = "red")
paste("Number of probes with a p-value ≤ 0.05: ",  dim(final_utest_th)[1])

#10. Apply multiple test correction and set a significant threshold of 0.05. How many probes do you identify as differentially methylated considering nominal pValues? How many after Bonferroni correction? How many after BH correction?

corrected_pValues_BH <- p.adjust(final_u_test$pValues_u_test,"BH")
corrected_pValues_Bonf <- p.adjust(final_u_test$pValues_u_test,"bonferroni")

final_u_test_corrected <- data.frame(final_u_test, corrected_pValues_BH, corrected_pValues_Bonf)

#How many probes survive the multiple test correction?
dim(final_u_test[final_u_test$pValues_u_test<=0.05,])
dim(final_u_test[final_u_test$corrected_pValues_BH<=0.05,])
dim(final_u_test_corrected[final_u_test_corrected$corrected_pValues_BH<=0.05,])
dim(final_u_test_corrected[final_u_test_corrected$corrected_pValues_Bonf<=0.05,])

# plotting
par(mfrow=c(1,1))
boxplot(final_u_test_corrected [,9:11], col = c("#7038FF", "#C7FF38", "seashell2"),names=NA, main='p-val before/after corrections')
legend("topright", legend=c("raw", "BH", "Bonferroni"),col=c("#7038FF", "#C7FF38", "seashell2"),pch=19, cex=0.5, xpd=TRUE)

#11. Produce a volcano plot and a Manhattan plot of the results of differential methylation analysis

#volcano plot
volcano_u_test <- data.frame(beta_preprocessNoob, pValues_u_test)
beta_volcano <- volcano_u_test[,1:8]
beta_volcano_groupCTRL <- beta_volcano[,pheno$Group=="CTRL"]
mean_beta_volcano_groupCTRL <- apply(beta_volcano_groupCTRL,1,mean)
beta_volcano_groupDIS <- beta_volcano[,pheno$Group=="DIS"]
mean_beta_volcano_groupDIS <- apply(beta_volcano_groupDIS,1,mean)

#Now we can calculate the difference between average values:
delta_volcano <- mean_beta_volcano_groupCTRL-mean_beta_volcano_groupDIS

#Now we create a dataframe with two columns, one containing the delta values and the other with the -log10 of p-values
toVolcPlot_volcano <- data.frame(delta_volcano, -log10(volcano_u_test$pValues_u_test))

#This line is useful to format in a good way the graph to avoid overlap of the legend with our graph, hiding data.
par(mar = c(5, 4, 4, 9), xpd = TRUE)

#create Volcano Plot
plot(toVolcPlot_volcano[,1], toVolcPlot_volcano[,2],
     pch = 16, cex = 0.5, col = "grey70",
     xlab = "delta beta (ctrl - dis)", ylab = "-log10(p-value)",
     xlim = c(-0.7, 0.7),
     ylim = c(0, max(toVolcPlot_volcano[,2])),
     asp = 0.5)

#add significance line at p-value threshold 0.05 
abline(h = -log10(0.05), col = "red", lty = 2, xpd = FALSE)

#add 10% differences in methylation line
abline(v = 0.1, col = "red", lty = 2, xpd = FALSE)
abline(v = -0.1, col = "red", lty = 2, xpd = FALSE)

# Highlight significant results
bh_sig <- toVolcPlot_volcano[
  volcano_u_test$pValues_u_test < 0.05 &
    abs(toVolcPlot_volcano[,1]) > 0.1, ]

points(bh_sig[,1], bh_sig[,2], pch = 16, cex = 0.8, col = "red")

# Legend
legend("bottomright",inset = c(-0.4, 0.12),
       legend = c("non significant",  "p < 0.05 & methylation difference ≥ 10%"),
       col = c("grey70", "red"),
       pch = 16, pt.cex = c(0.5, 0.8), cex = 0.7)

# manhattan plot
install.packages("qqman")
library(qqman)

#We prepare the annotation of the p-value per probe, recovering CHR, MAPINFO, IlmnID
final_u_test_manhattan <- data.frame(rownames(volcano_u_test),volcano_u_test)
colnames(final_u_test_manhattan)[1] <- "IlmnID"

final_u_test_manhattan_annotated <- merge(final_u_test_manhattan, Illumina450Manifest_clean,by="IlmnID")

# Now we can create the input for the Manhattan plot analysis. The input object should contain 4 info: probe, chromosome, position on the chromosome and p-value. We will select these columns in the final_u_test_manhattan_annotated object
input_Manhattan <- final_u_test_manhattan_annotated[colnames(final_u_test_manhattan_annotated) %in% c("IlmnID","CHR","MAPINFO","pValues_u_test")]

# It is better to reorder the levels of the CHR
order_chr <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
input_Manhattan$CHR <- factor(input_Manhattan$CHR,levels=order_chr )

# the column "CHR" should be numeric --> we will convert factors to numbers
input_Manhattan$CHR <- as.numeric(input_Manhattan$CHR)

# and finally we can produce our Manhattan plot
par(xaxs = "i")  
manhattan(input_Manhattan, snp="IlmnID",chr="CHR", bp="MAPINFO", p="pValues_u_test",col=rainbow(24) )
abline(h = -log10(0.05), col = "red", lty = 2)

# Extract the ids, chromosome and Mann-Whitney U test p-values of each probes with p-value above the threshold
significantly_probes <- input_Manhattan[input_Manhattan$pValues_u_test < 0.05, ]
#reorder these ids per chromosomes and insert them into the file.csv useful for possible further analysis
probes_table <- significantly_probes[, c("IlmnID", "pValues_u_test","CHR")]
probes_table <- probes_table[order(probes_table$CHR), ]
write.csv(probes_table, "significantly_probes.csv", row.names = FALSE)

# 12. Produce an heatmap of the top 100 significant, differentially methylated probes.

install.packages("gplots")
library(gplots)

# It wants as input a matrix. As an heatmap on several thousands of probes is computationally demanding, for our analysis we will use just the top 100 most significant CpG probes. We will extract only beta values for these top 100 probes and convert them in a matrix:
input_heatmap=as.matrix(volcano_u_test[1:100,1:8])

# Finally, we we create a bar of colors for Group A or Group B membership
colorbar_group <- c("green","orange","orange","green","orange","orange","green","green")

# We used average linkage mapping to find the distance relation between clusters
heatmap.2(input_heatmap,col=terrain.colors(100),Rowv=T,Colv=T,hclustfun = function(x) hclust(x,method = 'average'),dendrogram="both",key=T,ColSideColors=colorbar_group,density.info="none",trace="none",scale="none",symm=F,main="Average linkage - Groups")

#We performed also the analysis comparing sample using sex groups
colorbar_sex <- c("cyan","cyan","pink","cyan","pink","pink","pink","pink")

# Average
heatmap.2(input_heatmap,col=colorRampPalette(c("cyan","grey","pink"))(100),Rowv=T,Colv=T,hclustfun = function(x) hclust(x,method = 'average'),dendrogram="both",key=T,ColSideColors=colorbar_sex,density.info="none",trace="none",scale="none",symm=F,main="Average linkage - Sex")
