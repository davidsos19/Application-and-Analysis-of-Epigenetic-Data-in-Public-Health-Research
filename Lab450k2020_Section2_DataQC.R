#######################################################################
# EPIGENOMICS IN PUBLIC HEALTH, LAB CLASS
# Christine Ladd-Acosta, Andrew Jaffe, Kelly Bakulski, Brion Maher
# Johns Hopkins JUNE 1-4, 2020
#######################################################################


########################
# QUALITY CONTROL STEPS
########################

#################################
# VISUALIZE AND CLEAN RAW DATA
#################################

# MethylSet (Mset) contains metylated and unmethylated signals made using preprocessRaw()
rawMSet <- preprocessRaw( RGset )
rawMSet
#save( rawMSet, file = "rawMSet.rda" ) 
#M signal per probe, per sample
Meth <- getMeth( rawMSet )
Meth[ 1:5, 1:5 ]
#M signal per probe, per sample
Unmeth <- getUnmeth( rawMSet )
Unmeth[ 1:5, 1:5 ]

##################################################
#Overall intensity: M vs. U
pd$MQC <- log2( colMedians( Meth ) )
pd$UQC <- log2( colMedians( Unmeth ) )

pd$Array <- factor( pd$Array )
palette( graphColors )
pdf( "MvsUplot.pdf" )
plot( pd$UQC, pd$MQC, main = "M vs. U QC", pch = 16, xlab = "Log2 Median Unmethylated Intensity", ylab = "Log2 Median Methylated Intensity", cex.lab = 1.2, cex.main = 2 )
plot( pd$UQC, pd$MQC, col = as.factor( pd$Slide ), main = "M vs. U QC by Slide", pch = 16, xlab = "Log2 Median Unmethylated Intensity", ylab = "Log2 Median Methylated Intensity", cex.lab = 1.2, cex.main = 2 )
legend( "bottomright",levels(as.factor( pd$Slide ) ), fill = graphColors )
plot( pd$UQC, pd$MQC, col = pd$Array, main = "M vs. U QC by Position", pch = 16, xlab = "Log2 Median Unmethylated Intensity", ylab = "Log2 Median Methylated Intensity", cex.lab = 1.2, cex.main = 2 )
legend( "bottomright", levels( pd$Array ), fill = graphColors )
palette( myColors )
plot( pd$UQC, pd$MQC, col = pd$Batch, main = "M vs. U QC by Batch", pch = 16, xlab = "Log2 Median Unmethylated Intensity", ylab = "Log2 Median Methylated Intensity", cex.lab = 1.2, cex.main = 2 )
legend( "bottomright", c( "Batch 1", "Batch 2" ), fill = myColors )
dev.off()

#Drop (or if really small sample: watch out for): Samples with UQC<11 & MQC<11
length( which( pd$UQC < 11 ) )
length( which( pd$MQC < 11 ) )

##################################################
#Raw density plot
beta.raw <- getBeta( rawMSet )
mvalue.raw <- getM( rawMSet )
type <- getProbeType( rawMSet )
probe.type <- data.frame( Name = rownames( beta.raw ), Type = type )
pdf( "Density-plot-preproccessRaw.pdf" )
densityPlot( beta.raw, sampGroups = pd$casestatus, main = "Raw Beta by Case Status" )
densityPlot( beta.raw, sampGroups = pd$Batch, main = "Raw Beta by Batch" )
plotBetasByType( beta.raw[,1], probeTypes = probe.type, main = "Raw Beta by Probe Type, Sample 1" )
densityPlot( mvalue.raw, sampGroups = pd$Batch, main = "Raw M-value by Batch", xlab = "M-value" )
dev.off()

##################################################
#Raw PCA plots
beta.raw2 <- beta.raw
beta.raw2[ is.na( beta.raw2 ) ] <- 0

#makes the principal component object
prin <- prcomp( t( beta.raw2 ), center = T, scale. = F ) 

#pulls out proportion of variance explained by each PC
out.var <- prin$sdev^2 / sum( prin$sdev^2 ) #the percent variance at each PC
out.var[ 1:10 ]

sum( out.var )
pdf( "PCA-plots-Raw.pdf" )
screeplot( prin, col = "dodgerblue", xlab = "Principal Components of Raw Beta Values", main = "", cex.lab = 1.3 )
plot.new()
palette( myColors )
par( mar = c( 0, 0, 0, 0 ) )
legend( "bottom", levels( as.factor( pd$gender ) ), fill = myColors, title = "Principal Components by Sex" )
pairs( prin$x[,1:6], col = as.factor( pd$gender ), labels = c( "PC1", "PC2", "PC3", "PC4","PC5", "PC6" ), pch = 1, cex = 0.5 )
plot.new()
legend( "bottom", levels( as.factor( pd$Batch ) ), fill = myColors, title = "Principal Components by Batch" )
pairs( prin$x[,1:6], col = as.factor( pd$Batch ), labels = c( "PC1", "PC2", "PC3", "PC4","PC5", "PC6" ), pch = 1, cex = 0.5 )
plot.new()
palette( graphColors )
legend( "bottom", levels( as.factor( pd$Array ) ), fill = graphColors, title = "Principal Components by Array Position" )
pairs( prin$x[,1:6], col = as.factor( pd$Array ), labels = c( "PC1", "PC2", "PC3", "PC4","PC5", "PC6" ), pch = 1, cex = 0.5 )
plot.new()
legend( "bottom", levels( as.factor( pd$Slide ) ), fill = graphColors, title = "Principal Components by Slide" )
pairs( prin$x[,1:6], col = as.factor( pd$Slide ), labels = c( "PC1", "PC2", "PC3", "PC4","PC5", "PC6" ), pch = 1, cex = 0.5 )
dev.off()

##################################################
#Detection P 
detP <- detectionP( RGset )
dim( detP )
acorn( detP )
detP[ 1:5, 1:5 ]
#save(detP,file="detection-P.rda")
failed <- detP > 0.01
per.samp <- colMeans( failed ) # Fraction of failed positions per sample (length: 42)
summary( per.samp ) 
per.probe <- rowMeans( failed ) # Fraction of failed samples per position (length: 485512)
summary( per.probe )
sum( per.samp > 0.01 ) #How many samples had more than 1% of sites fail?
sum( per.probe > 0.1 ) # How many positions failed in >10% of samples? 

probe.fail <- failed[ per.probe > 0.1, ]
sample.fail <- per.samp[per.samp > 0.01 ]

pdf( "DetectionP-bySample.pdf" )
hist( per.samp, breaks = 40, col = "dodgerblue", cex.lab = 1.3, xlab = "Fraction of failed positions per sample" )
abline( v = 0.01, col = "red" )
hist( per.probe, breaks = 40, col = "dodgerblue", cex.lab = 1.3, xlab = "Fraction of failed samples per probe", ylim = c( 0,500 ) )
abline( v = 0.1, col = "red",lwd = 2 )
dev.off()

RGset.drop <- RGset[ ,!colnames( RGset ) %in% names( sample.fail ) ] #Drop samples that failed by detection p. 
dim( RGset.drop )
##################################################
#cross-reactive probes
load( "cross.probes.info.rda" )
dim( cross.probes.info )
head( cross.probes.info )
#rec is to use the 47 cross probe cutoff

##################################################
####VISUALIZE AND CLEAN NOOB DATA#################
##################################################

#Noob
noob <- preprocessNoob( RGset.drop, offset = 15, dyeCorr = TRUE, verbose = TRUE )

noob.dropP <- noob[ !rownames( noob ) %in% rownames( probe.fail ), ]
noob.dropCross <- noob.dropP[ !rownames( noob.dropP ) %in% cross.probes.info$TargetID, ]

noob <- noob.dropCross
#save(noob, file="noob.rda")
beta.n <- getBeta( noob )

pd <- pData( noob )
##################################################
##Sex Check
GmRawSet <- mapToGenome( rawMSet )
sex <- getSex( GmRawSet )
sex <- addSex( GmRawSet, sex = sex )
# #save(sex, file="Estimate-Sex.rda")
pd <- merge( pd, as.matrix( sex ), by = "row.names", sort = FALSE )
rownames( pd ) <- pd$Basename
pd <- pd[,-1]
table( pd$predictedSex, pd$gender )
     # F  M
  # F 22  0
  # M  0 18
pdf( "Sex-Plot.pdf" )
plotSex( sex )
dev.off()

##################################################
#Cell type
library( FlowSorted.Blood.450k )
pData( RGset )[ ,"Array" ] <- "Array.x"
pData( RGset )[ ,"Basename" ] <- "Basename.x"
cell <- estimateCellCounts( RGset.drop )
#save(cell, file="cell-type-estimates.rda")
dim( cell )
head( cell )
identical( rownames( cell ),colnames( noob ) ) # A quick sanity check to make sure our samples are in the right order. 
pd.n <- data.frame( pd,cell )

prin.cell <- prcomp( t( cell ), center = T, scale. = F ) 
out.var <- prin.cell$sdev^2 / sum( prin.cell$sdev^2 )
out.var


pdf( "Cell-Type-PCA.pdf" )
screeplot( prin.cell, col = "dodgerblue", xlab = "Principal Components of Estimated Cell Type", main = "", cex.lab = 1.3 )
myColors <- c( "seagreen3","dodgerblue", "darkorchid", "firebrick1", "darkorange", "khaki1" )
palette( myColors )
par( mar = c( 0, 0, 0, 0 ) )
plot.new()
legend( "bottom", c( "Bcell","CD4T", "CD8T","Gran","Mono", "NK" ), fill = myColors, title = "Principal Components by Estimated Cell Type")
pairs( prin.cell$x[,1:6], col = as.factor(rownames( prin.cell$x ) ), labels = c( "PC1", "PC2", "PC3", "PC4","PC5", "PC6" ), pch = 15 )
dev.off()

summary( pd.n$Gran )
pdf( "Cell-Type-Histogram.pdf" )
hist( pd.n$Gran, breaks = 20, col = "dodgerblue", xlab = "Estimated percent granulocytes", cex.lab = 1.3 )
hist( pd.n$Mono, breaks = 20, col = "dodgerblue", xlab = "Estimated percent monocytes", cex.lab = 1.3 )
hist( pd.n$NK, breaks = 20, col = "dodgerblue", xlab = "Estimated percent NK cells", cex.lab = 1.3 )
hist( pd.n$CD8T, breaks = 20, col = "dodgerblue", xlab = "Estimated percent CD8T", cex.lab = 1.3 )
hist( pd.n$CD4T, breaks = 20, col = "dodgerblue", xlab = "Estimated percent CD4T", cex.lab = 1.3 )
dev.off()

#Pick cutoffs for biological ranges for cell type estimates.
pd.n.drop <- pd.n[ pd.n$Mono < 0.2, ]
noob.drop <- noob[ ,colnames( noob ) %in% rownames( pd.n.drop ) ]
beta.noob <- getBeta( noob.drop )
dim( beta.noob )


##################################################
#Noob PCA plots
pd.n.drop$gran.quart <- cut2( pd.n.drop$Gran, g = 4 )

prin <- prcomp( t( beta.noob ), center = T, scale. = F )
out.var <- prin$sdev^2 / sum( prin$sdev^2 )
out.var[ 1:10 ]

ls()
rm( list = c( "RGset", "beta.raw" ) )

pdf( "PCA-plots-Noob.pdf" )
screeplot( prin, col="dodgerblue", xlab="Principal Components of Noob Beta Values", main="", cex.lab=1.3 )
plot.new()
palette( myColors )
par( mar = c( 0, 0, 0, 0 ) )
legend( "bottom", levels( as.factor( pd.n.drop$gender ) ), fill = myColors, title = "Principal Components by Sex" )
pairs(prin$x[,1:6], col = as.factor( pd.n.drop$gender ), labels = c( "PC1", "PC2", "PC3", "PC4","PC5", "PC6" ), pch = 1, cex = 0.5 )
plot.new()
legend( "bottom", levels( as.factor( pd.n.drop$Batch ) ), fill = myColors, title = "Principal Components by Batch" )
pairs(prin$x[,1:6], col = as.factor( pd.n.drop$Batch ), labels = c( "PC1", "PC2", "PC3", "PC4","PC5", "PC6" ), pch = 1, cex = 0.5 )
plot.new()
legend( "bottom", levels( as.factor( pd.n.drop$casestatus ) ), fill = myColors, title = "Principal Components by Case Status" )
pairs( prin$x[,1:6], col = as.factor( pd.n.drop$casestatus ), labels = c( "PC1", "PC2", "PC3", "PC4","PC5", "PC6" ), pch = 1, cex = 0.5 )
plot.new()
palette( graphColors )
legend( "bottom", levels( as.factor( pd.n.drop$Array ) ), fill = graphColors, title = "Principal Components by Array Position" )
pairs( prin$x[,1:6], col = as.factor( pd.n.drop$Array ), labels = c( "PC1", "PC2", "PC3", "PC4","PC5", "PC6" ), pch = 1, cex = 0.5 )
plot.new()
legend( "bottom", levels( as.factor( pd.n.drop$Slide ) ), fill = graphColors, title = "Principal Components by Slide" )
pairs( prin$x[,1:6], col = as.factor( pd.n.drop$Slide ), labels = c( "PC1", "PC2", "PC3", "PC4","PC5", "PC6" ), pch = 1, cex = 0.5 )
plot.new()
legend("bottom", levels( as.factor( pd.n.drop$gran.quart ) ), fill = graphColors, title = "Principal Components by Quartile of Gran" )
pairs( prin$x[,1:6], col = as.factor( pd.n.drop$gran.quart ), labels = c( "PC1", "PC2", "PC3", "PC4","PC5", "PC6"), pch = 1, cex = 0.5 )
dev.off()

##################################################
#Drop samples with missing data at key covariates
table( pd.n.drop$gender, useNA = "always" )
table( pd.n.drop$casestatus, useNA = "always" )
table( pd.n.drop$smoking, useNA = "always" )
table( pd.n.drop$Batch, useNA = "always" )
#No samples missing data, we're okay. But run this line anyway. 
pd.complete <- pd.n.drop[ !( is.na( pd.n.drop$smoking ) ), ] #Example of how to remove. 

##################################################
####VISUALIZE AND CLEAN COMBAT DATA###############
##################################################

beta.noob.complete <- beta.noob[ ,colnames( beta.noob ) %in% rownames( pd.complete ) ]
mod <- model.matrix( ~ pd.complete$gender + pd.complete$casestatus + pd.complete$smoking )

library( nlme )
library( sva )
combat.beta <- ComBat( dat = beta.noob.complete, batch = pd.complete$Batch, mod = mod )
#Optional, save the cleaned matrix.
#save(combat.beta, file="combat-beta.rda")
combat.beta[ 1:5,1:5 ]

prin <- prcomp( t( combat.beta ), center = T, scale. = F )
out.var <- prin$sdev^2 / sum( prin$sdev^2 )
out.var[ 1:10 ]

pd.n.drop$gran.quart <- cut2( pd.n.drop$Gran, g = 4 )

pdf( "PCA-plots-Noob-Combat.pdf" )
screeplot( prin, col = "dodgerblue", xlab = "Principal Components of Noob Beta Values", main = "", cex.lab = 1.3 )
plot.new()
palette( myColors )
par( mar = c( 0, 0, 0, 0 ) )
legend( "bottom", levels( as.factor( pd.n.drop$gender ) ), fill = myColors, title = "Principal Components by Sex" )
pairs( prin$x[,1:6], col = as.factor( pd.n.drop$gender ), labels = c( "PC1", "PC2", "PC3", "PC4","PC5", "PC6" ), pch=1, cex = 0.5 )
plot.new()
legend( "bottom", levels( as.factor( pd.n.drop$Batch ) ), fill = myColors, title = "Principal Components by Batch" )
pairs( prin$x[,1:6], col = as.factor( pd.n.drop$Batch ), labels = c( "PC1", "PC2", "PC3", "PC4","PC5", "PC6" ), pch = 1, cex = 0.5 )
plot.new()
legend ("bottom", levels( as.factor( pd.n.drop$casestatus ) ), fill = myColors, title = "Principal Components by Case Status" )
pairs( prin$x[,1:6], col = as.factor( pd.n.drop$casestatus ), labels = c( "PC1", "PC2", "PC3", "PC4","PC5", "PC6" ), pch = 1, cex = 0.5 )
plot.new()
palette( graphColors )
legend( "bottom", levels( as.factor( pd.n.drop$Array ) ), fill = graphColors, title = "Principal Components by Array Position" )
pairs( prin$x[,1:6], col = as.factor( pd.n.drop$Array ), labels = c( "PC1", "PC2", "PC3", "PC4","PC5", "PC6"), pch = 1, cex = 0.5 )
plot.new()
legend( "bottom", levels( as.factor( pd.n.drop$Slide ) ), fill = graphColors, title = "Principal Components by Slide" )
pairs( prin$x[,1:6], col = as.factor( pd.n.drop$Slide ), labels = c( "PC1", "PC2", "PC3", "PC4","PC5", "PC6" ), pch = 1, cex = 0.5 )
plot.new()
legend( "bottom", levels( as.factor( pd.n.drop$gran.quart ) ), fill = graphColors, title="Principal Components by Quartile of Gran" )
pairs( prin$x[,1:6], col = as.factor( pd.n.drop$gran.quart ), labels = c( "PC1", "PC2", "PC3", "PC4","PC5", "PC6"), pch = 1, cex = 0.5 )
dev.off()

table( pd.n.drop$gran.quart )

