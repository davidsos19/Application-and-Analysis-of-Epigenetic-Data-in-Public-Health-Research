#Please run the following code to install various R packages before coming to class.

### Install BiocManager -- Allows user to install Bioconductor packages
if ( !requireNamespace( "BiocManager" ) ) install.packages( "BiocManager" )

### Install minfi -- 450k analysis package
BiocManager::install( "minfi" )

### Install RColorBrewer -- Useful package that provides color palettes to help visualize data
install.packages( "RColorBrewer" )

### Install sva -- Functions for batch effect correction
BiocManager::install( "sva" )

### Install limma -- Functions for single-site association analysis
BiocManager::install( "limma" )

### Install ggplot2 -- Easy plotting functions
install.packages( "ggplot2" )

### Install Go.db and topGo -- Functions for gene ontology queries
BiocManager::install( c( "GO.db","topGO","missMethyl", "GEOquery" ) )

### Install matrixEQTL --	For fast eQTL and meQTL analyses
install.packages( "MatrixEQTL" )

### Install qqman -- Functions for q-q and manhattan plots
install.packages( "qqman" )

### Install nlme -- Functions for non-linear mixed effects models
install.packages( "nlme" )

### Install other R packages
install.packages( c( "MaxtrixStats", "MASS", "abind", "Hmisc", "BiasedUrn" ) )

### Install packages for using 450k
BiocManager::install( c( "IlluminaHumanMethylation450kmanifest","IlluminaHumanMethylation450kanno.ilmn12.hg19", "FlowSorted.Blood.450k" ) )
