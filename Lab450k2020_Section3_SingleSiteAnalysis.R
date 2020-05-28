#######################################################################
# EPIGENOMICS IN PUBLIC HEALTH, LAB CLASS
# Christine Ladd-Acosta, Andrew Jaffe, Kelly Bakulski, Brion Maher
# Johns Hopkins JUNE 1-2, 2020
#######################################################################

##################################################
# SECTION 3
##################################################

##################################################
# SINGLE-SITE ASSOCIATION ANALYSIS
##################################################

# noob <- noob[ ,colnames( noob ) %in% rownames( pd.complete ) ]
# pd <- pData( noob )

# cell.n <- cell[ rownames( cell ) %in% colnames( noob ),]

# pd <- data.frame( pd, cell.n )

#Make sure these output to TRUE
identical( colnames( combat.beta ), colnames( noob ) )
identical( rownames( pd.complete ), colnames( noob ) )
identical( rownames( noob ), rownames( combat.beta ) )

#Load the limma package
library( limma )

mapped <- mapToGenome( noob )
noob <- noob[ match( rownames( mapped ), rownames( noob ) ), ]
combat.beta <- combat.beta[ match( rownames( mapped ),rownames( combat.beta ) ), ]

pre.beta <- getBeta( noob )
chrnames <- as.character( seqnames( mapped ) )
pos <- as.numeric( start( mapped ) )

#Compute PCs for cell type

cellscores <- prcomp( ( pd.complete[ ,15:20 ] ), center = T, scale. = F )
round( ( cellscores$sdev^2 / sum( cellscores$sdev^2 ) )*100, digits = 2 ) #Variance explained by each PC

#Box Plots and t-statistics for each cell type
library( ggplot2 )
celltypes <- 15:20
pdf( "CompareCellComp.pdf" )
for ( i in celltypes ){
	myttest <- t.test( pd.complete[ which( pd.complete$casestatus == "RA" ), i ], pd.complete[ which( pd.complete$casestatus == "Control" ), i ] )
	thisplot <- data.frame( pd.complete[ ,i ], pd.complete$casestatus )
	colnames( thisplot ) <- c( "CellComp", "Status" )
	p1 <- ggplot( thisplot, aes( factor( Status ), CellComp ) ) + geom_boxplot( ) + xlab( "Case Status" ) +
	labs( title = paste( "Cell Type", colnames( pd.complete )[i], ", p-value of ", format( myttest$p.value, digits = 3, scientific = TRUE ), sep = "" ) ) +
	scale_y_continuous( limits = c( 0,1 ) )
	print( p1 )
}
dev.off()


#Construct the model matrix
mod <- model.matrix( ~ factor( pd.complete$casestatus ) + pd.complete$age + factor( pd.complete$gender ) + factor( pd.complete$smoking ) + cellscores$x[ ,1 ] )  

#Run the single site association model
out <- lmFit( combat.beta, mod )
out <- eBayes( out )
ss.hits <- topTable( out, coef = 2, number = nrow( combat.beta ) )

#Make a qq plot of our data
#install.packages( "qqman" )
library( qqman )
#png("qqplot.png",type="cairo")
qq( ss.hits$P.Value )
#dev.off()

#Calculate a lambda statistic
observed = -log10( sort( ss.hits$P.Value, decreasing = F ) )
expected = -log10( ppoints( length( ss.hits$P.Value ) ) )
lambda = median( observed ) / median( expected )

########Make a manhattan plot of our data########
#Grab the correct position and chromosome vectors
pos.ss <- pos[ match( rownames( ss.hits ), rownames( combat.beta ) ) ]
chr.ss <- chrnames[match(rownames(ss.hits), rownames( combat.beta ) ) ]
chr.ss.num <- as.numeric( unlist( lapply( chr.ss, function( x ){ strsplit( x, "chr", fixed = TRUE )[[1]][2] } ) ) )

#Construct data frame for function
formanhattan <- data.frame( pos.ss, chr.ss.num, ss.hits$P.Value )
colnames( formanhattan ) <- c( "BP","CHR","P" )

formanhattan <- formanhattan[ -which( is.na( formanhattan$CHR ) ), ]

#Call function
#png("manhattan.png",type="cairo")
manhattan( formanhattan )
#dev.off()

########Plot the top 6 hits and examine case vs control beta differences########
#Create a mini result list of the first 6 hits
mytrunc <- ss.hits[ 1:6, ]

#Create a results object to populate with the case-control differences
mymeans <- matrix( 0, nrow( mytrunc ), 2 )


#Plot the results: here we are formatting the data in a certain way in order to fit
#the ggplot2 plotting function. We have also shifted the position of the controls by 10bp
#so that their methylation distribution is easier to see. We are also populating our object
#"mymeans" with the case (first column) and control (second column) means at each site. 
pdf( "SingleSiteExamples.pdf" ) 
for ( k in 1:nrow( mytrunc ) ){
	probe <- which(rownames( pre.beta ) == rownames( mytrunc )[ k ] )
	grabBeta <- pre.beta[ probe, ]
	
		formbeta <- c()
		for ( leonard in 1:length( grabBeta ) ){
			tempbeta <- grabBeta[ leonard ]
			formbeta <- c( formbeta, tempbeta )
		}
		reppos <- rep( pos.ss[ k ], each = ncol( pre.beta ) )
	#reppos.real<-pos[reppos]
		xmin <- min( reppos ) - 10
		xmax <- max( reppos ) + 20
		status <- pd.complete$casestatus
		reppos.real.shift <- ifelse( status == "Control", reppos + 10, reppos )
		toplot <-data.frame( ( as.numeric( formbeta ) ), jitter( reppos.real.shift ), status )
		rownames( toplot ) <- c()
		#toplot<-data.frame(toplot)
		colnames( toplot ) <- c( "Beta","Position","Status" )
		
	casemean <- mean( as.numeric( grabBeta[ which( status == "RA" ) ] ) )
	contmean <- mean( as.numeric( grabBeta[ which( status == "Control" ) ] ) )
	casedat <- data.frame( x = ( pos.ss[k] - 2 ):( pos.ss[k] + 2 ), y = mean( as.numeric( grabBeta[ which( status == "RA" ) ] ) ), color = "red" )
	contdat <- data.frame( x = ( pos.ss[k] + 8 ):( pos.ss[k] + 12), y = mean( as.numeric( grabBeta[ which( status == "Control" ) ] ) ), color = "chartreuse" )

	mymeans[ k,1 ] <- mean( as.numeric( grabBeta[ which( status == "RA" ) ] ) )
	mymeans[ k,2 ] <- mean( as.numeric( grabBeta[ which( status == "Control" ) ] ) )
	
	meanDiff <- round( ( mean( as.numeric( grabBeta[ which( status == "RA" ) ] ) ) - mean( as.numeric( grabBeta[ which( status == "Control" ) ] ) ) )*100, digits = 1 )
	pVal <- format( mytrunc$P.Value[k], digits = 3 )
	
		p2 <- ggplot( data = toplot, aes( x = ( Position ), y = ( Beta ), color = factor( Status ) ) ) + geom_point( size = 2 ) + scale_colour_manual( values = c( "dodgerblue","black" ) ) +
		theme( legend.direction = "vertical", legend.position = c( 0.1,0.95 ), legend.title = element_blank( ), panel.background = element_blank( ) ) +
		ylab( "Percent Methylation" ) + xlab( paste0( rownames( mytrunc )[k], ":", " Delta M = ", meanDiff, " at p-val ", pVal ) ) +
		scale_x_continuous( limits = c( xmin, xmax ) ) + theme( axis.text.x = element_blank( ) ) +
		scale_y_continuous( limits = c( 0,1 ), breaks = seq( 0,1, by = 0.25 ), labels = c( "0","25","50","75","100" ) ) +
		geom_line( data = casedat, aes( x = x, y = y ), color = "red" ) + geom_line( data = contdat, aes( x = x, y = y ), color = "chartreuse" )

		#+stat_summary(fun.y=mean,geom="line",size=1)
		print( p2 )
}
dev.off()

##################################################
# GENE ONTOLOGY ANALYSES
##################################################

library( missMethyl )
library( BiasedUrn )
#Need to pick a significance cutoff for inclusion. Here is 1x10-5, but may need to be flexible with data.
gene.ontology <- gometh( as.character( rownames( ss.hits )[ ss.hits$P.Value < 1e-5 ] ), all.cpg = as.character( rownames( ss.hits ) ), plot.bias = FALSE, prior.prob = TRUE )
dim( gene.ontology )
summary( gene.ontology$P.DE )
summary( gene.ontology$FDR )
gene.ontology <- gene.ontology[ order( gene.ontology$P.DE ), ]
head( gene.ontology )

#Some functions that will be useful in making the plot
wrap.it <- function( x, len )
{ 
  sapply( x, function( y ) paste( strwrap( y, len ), 
                              collapse = "\n" ), 
         USE.NAMES = FALSE )
}
# Call this function with a list or vector
wrap.labels <- function( x, len )
{
  if ( is.list( x ) )
  {
    lapply( x, wrap.it, len )
  } else {
    wrap.it( x, len )
  }
}

pdf( file = "Enriched_Gene_Ontology.pdf", paper = "a4r" )
# # par(mar=c(4, 6, 2, 0)+0.1)
# par(pin=c(2,4))
# par(oma=c(0,0,0,0))
par( mai = c( 1,4,1,1 ) )
barplot( abs( log( as.numeric( gene.ontology$P.DE[ 1:10 ] ), base = 10 ) ), main = "Gene Ontology", horiz = TRUE, names.arg = wrap.labels( gene.ontology$TERM[ 1:10 ], 50 ), 
	xlab = "-log10(P value)", col = "dodgerblue", las = 1, cex.axis = 1.2, cex.main = 1.4, cex.lab = 1, cex.names = 1, space = 0.4, xlim = c( 0,5 ) )
dev.off()
