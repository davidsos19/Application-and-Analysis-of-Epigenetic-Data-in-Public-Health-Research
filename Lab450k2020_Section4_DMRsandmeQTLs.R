#######################################################################
# EPIGENOMICS IN PUBLIC HEALTH, LAB CLASS
# Christine Ladd-Acosta, Andrew Jaffe, Kelly Bakulski, Brion Maher
# Johns Hopkins JUNE 1-2, 2020
#######################################################################

##################################################
# SECTION 4
##################################################

##################################################
# REGION-FINDING ASSOCIATION ANALYSIS
##################################################

mybumps <- bumphunter( combat.beta, mod, chr = chrnames, pos = pos,
	pickCutoff = TRUE, pickCutoffQ = 0.975, maxGap = 300, smooth = TRUE,
	smoothFunction = loessByCluster, B = 10, nullMethod = "bootstrap" )

#Annotate the bumphunter output with gene information
mybumps <- mybumps$table

#Plot the top 2 DMRs: here we are formatting the data in a certain way in order to fit
#the ggplot2 plotting function. We have also shifted the position of the controls by 10bp
#so that their methylation distribution is easier to see. 
pdf( "DMRexample.pdf" )  ##############################
for ( i in 1:2 ){
#for (i in 1:1){
	indexstart <- mybumps$indexStart[ i ]
	indexend <- mybumps$indexEnd[ i ]
	grabBeta <- pre.beta[ indexstart:indexend, ]
	formbeta <- c()
	for ( j in 1:nrow( grabBeta ) ){
		tempbeta <- grabBeta[ j, ]
		formbeta <- c( formbeta, tempbeta )
	}
	reppos <- rep( pos[ indexstart:indexend ], each = ncol( pre.beta ) )
	#reppos.real<-pos[reppos]
	xmin <- min( reppos ) - 100
	xmax <- max( reppos ) + 100
	status <- rep( pd$casestatus, nrow( grabBeta ) )
	reppos.real.shift <- ifelse( status == "Control", reppos + 10, reppos )
	toplot <- data.frame( ( formbeta ), reppos.real.shift, status )
	rownames( toplot )<-c( )
	colnames( toplot )<-c("Beta","Position","Status")
	#toplot$Status<-ifelse(toplot$Status==2,"normal","autism")
	p2 <- ggplot( data= toplot, aes( x = ( Position ),y = ( Beta ),color = factor( Status ) ) ) + geom_point( size = 0.75 ) + scale_colour_manual(values = c( "dodgerblue","black" ) ) +
	theme( legend.direction = "horizontal", legend.position = c( 0.5,0.95 ), legend.title = element_blank( ), panel.background = element_blank( ) ) +
	ylab( "Percent Methylation" ) + xlab( "Position" ) +
	scale_x_continuous( limits = c( xmin, xmax ) ) + stat_smooth( method = "loess", se = FALSE ) +
		scale_y_continuous( limits = c( 0,1 ), breaks = seq( 0, 1, by = 0.25 ), labels = c( "0","25","50","75","100" ) )
		#+stat_summary(fun.y=mean,geom="line",size=1)
	print( p2 )
}
dev.off()

##################################################
# meQTL analyses
##################################################

#Load the matrixEQTL package
library( MatrixEQTL )

#Load the genotype object
load( "genotypes.rda" )
load( "snp.pos.rda" )

#Subset methylation data and genomic positions to only our 
#chromosome of interest
B.mychr <- combat.beta[ which( chrnames == "chr22" ), ]
pos.mychr <- pos[ which( chrnames == "chr22" ) ]

#Subset the methylation data to probes in this particular region
B.LDblock <- B.mychr[ which( pos.mychr > 17583446 & pos.mychr < 17666446 ), ]   
pos.LDblock <- pos.mychr[ which( pos.mychr > 17583446 & pos.mychr < 17666446 ) ]

#Let's be sure to match to the samples in our methylation matrix
genotypes <- genotypes[, match (colnames( B.LDblock ), colnames( genotypes ) ) ]

#Format the genotypes and methylation objects for the package
genotypes.format <- SlicedData$new( genotypes )
meth.format <- SlicedData$new( B.LDblock )

#Call the function
results <- Matrix_eQTL_main( genotypes.format, meth.format, pvOutputThreshold = 0.05,
	snpspos = snp.pos,genepos = pos.LDblock, output_file_name = NULL, output_file_name.cis = NULL )

#Grab the meQTL results in a separate object for ease of plotting downstream
results.table <- results$all$eqtls

#Let's examine the spatial relationship between SNPs and CpG sites.
#First we make sure we grab the right positions
results.table$SNPpos <- snp.pos[ match( results.table$snps, rownames( genotypes ) ) ]
results.table$CGpos <- pos.LDblock[ match( results.table$gene, rownames( B.LDblock ) ) ]

#Transform the p-value to the -log10 scale
results.table$transP <- ( -1*log( results.table$pvalue, base = 10 ) )

results.table$Distance <- results.table$CGpos - results.table$SNPpos

pdf( "VolcanomeQTL.pdf" )
with( results.table, plot( Distance, transP, pch = 20, main = NULL, 
	xlab = "CpG Position - SNP Position (Kb)", ylab = "-log10 p-value" ) )
dev.off()
	

##################################################
# ADDENDUM: READ DATA FROM GEO
##################################################

#Use the 'GEOquery' package
library( GEOquery )
setwd( "WHERE YOU WANT TO DOWNLOAD THE RAW FILES" )

#Download the supplementary files attached to this 
#GEO ID. Raw .idat files are part of this group of
#supplementary files. 
getGEOSuppFiles( "GSE42861", makeDirectory = TRUE, baseDir = getwd( ) )

#General data available for this GEOID
mystudy <- getGEO( GEO = "GSE42861", destdir = getwd( ) )

#Phenotype/covariate information for these samples  
mypheno <- ( phenoData( mystudy$GSE42861_series_matrix.txt.gz ) )
variables <- varMetadata( mypheno )