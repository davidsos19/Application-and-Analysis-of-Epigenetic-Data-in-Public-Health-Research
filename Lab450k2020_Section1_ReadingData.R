#######################################################################
# EPIGENOMICS IN PUBLIC HEALTH, LAB CLASS
# Christine Ladd-Acosta, Andrew Jaffe, Kelly Bakulski, Brion Maher
# Johns Hopkins JUNE 1-4, 2020
#######################################################################


##################################################
### INSTALL RELEVANT PACKAGES (ONLY DO ONCE)
##################################################

if ( !requireNamespace( "BiocManager" ) ) install.packages( "BiocManager" )
BiocManager::install( "minfi" )
BiocManager::install( "limma" )
BiocManager::install( "IlluminaHumanMethylation450kmanifest" )
BiocManager::install( "IlluminaHumanMethylation450kanno.ilmn12.hg19" )
BiocManager::install( "sva" )
BiocManager::install( "FlowSorted.Blood.450k" )
install.packages( "RColorBrewer" )
install.packages( "ggplot2" )
install.packages( "matrixStats" )
install.packages( "MASS" )
install.packages( "abind" )
install.packages( "Hmisc" )



########################################################
###LOAD RELEVANT PACKAGES (EVERY TIME YOU DO ANALYSIS)
########################################################

library( minfi )
library( limma )
library( matrixStats )
library( MASS )
library( abind )
library( sva )
library( Hmisc )

#Set up color palettes for plotting
myColors <- c( "dodgerblue", "firebrick1", "seagreen3" )
graphColors = c( "#023FA5","#7D87B9","#BEC1D4","#D6BCC0","#BB7784", "#D33F6A", "#11C638","#8DD593","#C6DEC7","#EAD3C6",
                "#F0B98D","#EF9708", "#0FCFC0","#9CDED6","#D5EAE7","#F3E1EB","#F6C4E1","#F79CD4", "#4A6FE3","#8595E1",
                "#B5BBE3","#E6AFB9","#E07B91" )


####################################################
###SET YOUR WOKRING DIRECTORY AND READ IN RAW DATA
####################################################

#The rest of this script assumes that your
#working directory is the same as the folder
#into which you downloaded the raw .idat files
#from the Box link. Two options:
	#1) Place the file path to this folder between the quotation marks in the setwd command below
	#2) Manually change the working directory in RStudio: Session > Set Working Directory > Choose Directory
setwd( "C:/Users/David Sosnowski/Desktop/EpiLab2020" )  
getwd()
pheno <- read.csv( "samplesheet.csv", header = TRUE, stringsAsFactors = FALSE )
dim( pheno )
head( pheno )
RGset <- read.metharray.exp( file.path( "C:/Users/David Sosnowski/Desktop/EpiLab2020/idats" ), targets = pheno, verbose = TRUE )
dim( RGset )
manifest <- getManifest( RGset )
str( manifest )
annotation <- getAnnotation( RGset )
dim( annotation )
annotation[ 1:2, ]

typeof( annotation )
typeof( RGset )
getClass( RGset )
manifest
head( getProbeInfo( manifest ) )
dim( getProbeInfo( manifest ) )
table( getProbeInfo( manifest )$Color )
pd <- pData( RGset )
table( pd$casestatus )
table( pd$gender )
table( pd$casestatus, pd$gender )
table( pd$Batch)
table( pd$casestatus, pd$Batch )
dim( pd )
length( pd$GEOID )
summary( pd$age )
summary( pd$age[ pd$gender == "M" ] )
summary( pd$age[ pd$gender == "F" ] )
head( pd )

