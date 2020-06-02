## via Hadley Wickham's R for Data Science
## https://r4ds.had.co.nz

## installing packages - you only need to do once
install.packages("tidyverse")

## using a package - you need to do in every new session
library(tidyverse)

## the pound symbol (#) is for writing comments
## the ?`function` searches for help

## example exploration
mpg
head(mpg)

## plots using ggplot2
ggplot(data = mpg) + 
  geom_point(mapping = aes(x = displ, y = hwy))

## add color
ggplot(data = mpg) + 
  geom_point(mapping = aes(x = displ, y = hwy, color = class))

## or shape
ggplot(data = mpg) + 
  geom_point(mapping = aes(x = displ, y = hwy, shape = drv))

## or both
ggplot(data = mpg) + 
  geom_point(mapping = aes(x = displ, y = hwy, shape = drv, color = class))

## and facets
ggplot(data = mpg) + 
  geom_point(mapping = aes(x = displ, y = hwy, color = drv)) + 
  facet_wrap(~ class, nrow = 2)

## different geometries
ggplot(data = mpg) + 
  geom_smooth(mapping = aes(x = displ, y = hwy, color = drv)) 

## plots using ggplot2 "quick" plot
qplot(x=displ, y=hwy, data=mpg, color=drv, facets = ~class)

## see more at ggplot2 homepage: https://ggplot2.tidyverse.org/
## more examples from our weeklong summer and winter class
## https://johnmuschelli.com/intro_to_r/Data_Visualization/Data_Visualization.html

## plots using base R
plot(x = mpg$displ, y= mpg$hwy)
plot(hwy ~ displ, data = mpg)

########################
## selecting columns ###
########################

## dplyr is a relatively new grammar of data manipulation
## and is generally easier than base R syntax

## selecting columns from mpg
select(mpg, manufacturer, model, displ)
mpg[,c("manufacturer", "model", "displ")]

## selecting all but certain columns
select(mpg, -drv)

## there are helpers
select(mpg, starts_with("c"))

## reorder columns
select(mpg, class, everything())

######################
### filtering rows ###
######################

## by categories
filter(mpg, class == "compact")

## by numerical variables
filter(mpg, cty > 20)

## by multiple columns
filter(mpg, model == "a4", drv == "f", hwy > 25)
filter(mpg, model == "a4" & drv == "f" & hwy > 25)
filter(mpg, model == "a4" | drv == "f" | hwy > 25)

## by multiple levels/values
filter(mpg, manufacturer %in% c("audi", "subaru"))
filter(mpg, manufacturer == "audi" | manufacturer =="subaru")

head(mpg[mpg$manufacturer == "audi",])

mpg$manufacturer[1] = "audi"

df= as.data.frame(mpg)
####################
# rows and columns #
####################

## nested
select(filter(mpg, model =="a4"), class, everything())

## "piping" is newer
mpg %>% filter(model == "a4") %>% select(class, everything())

## order matters here, not run cuz error
# mpg %>% select(class, drv, manufacturer) %>% filter(model == "a4")

#################
## assignment ###
#################

## all previous examples just displayed results
## you often want to retain filtering or selecting

## the equal sign is the assignment operator
audi = mpg %>% filter(manufacturer == "audi")

## <- also works, but is the same for 99.9% of uses
audi <- mpg %>% filter(manufacturer == "audi")

#################
### summarize ###
#################

## count elements, one way
mpg %>% count(manufacturer)

## count elements, two way
mpg %>% count(year, cyl)

## grouping of tibbles
mpg_group = mpg %>% group_by(class)

## tabulate here
mpg_group %>% tally()

## summaries
summarize(mpg_group, cty_avg = mean(cty))
mpg_group %>% summarize(cty_avg = mean(cty))

## one liner
mpg %>% group_by(class) %>% summarize(cty_avg=mean(cty))

####################################################
### read in results from differential expression  ##
####################################################

## importing data from your computer is important
## and is the starting point for any analysis
## we will use the `readr` package here, but there 
## are also a series of base R functions like read.csv

## lets read in the results from a recent preprint
## https://www.biorxiv.org/content/10.1101/426213v1
res = read_csv("data/SupplementaryTable2_sczd_gene_full.csv")

## filter on expression level and drop `gencodeTx` column
resExprsDlpfc = res %>% select(-gencodeTx) %>% filter(region == "DLPFC")
resExprsHippo = res %>% select(-gencodeTx) %>% filter(region == "HIPPO")

## volcano plot
ggplot(resExprsDlpfc) + 
  geom_point(mapping = aes(x = logFC, 
          y = -log10(P.Value), color = adj.P.Val < 0.05))

ggplot(resExprsHippo) + 
  geom_point(mapping = aes(x = logFC, 
                           y = -log10(P.Value), color = adj.P.Val < 0.05))

####################
## Bioconductor ####
####################

## installing bioconductor packages uses
##  different syntax than base R packages
##  read more: http://bioconductor.org/install/
## this has also changed relatively recently
## theres a ton of courses here: http://bioconductor.org/help/course-materials/

## you need the BiocManager CRAN package now
if (!requireNamespace("BiocManager"))  install.packages("BiocManager")

## more here in their vignette:
## https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html
BiocManager::install("GenomicRanges")
##this installs a bunch of packages and takes a bit of time

## we can check out a few of these course materials while that installs
# http://bioconductor.org/help/course-materials/2017/OSU/B1_Bioconductor_intro.html#bioconductor-sequencing-ecosystem 

### the rest is modeled after some of these courses
# http://bioconductor.org/help/course-materials/2017/OSU/B2_Common_Operations.html
# http://bioconductor.org/help/course-materials/2017/BioCAsia/A01-Bioc-Basics.html


#############################
## GenomicRanges / GRanges ##
#############################

# https://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html

library(GenomicRanges)

## ----example-GRanges-------------------------------------------------------
gr <- GRanges(
  seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
  ranges = IRanges(101:110, end = 111:120, names = head(letters, 10)),
  strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
  score = 1:10,
  GC = seq(1, 0, length=10))
gr


## ----example-overlaps---------------------------------------------------
gr1 <- GRanges(
  seqnames = "chr2",
  ranges = IRanges(103, 106),
  strand = "+",
  score = 5L, GC = 0.45)

## ----countOL---------------------------------------------------------------
countOverlaps(gr, gr1)

## ----subsetByOverlaps------------------------------------------------------
subsetByOverlaps(gr,gr1)

## ----findOverlaps------------------------------------------------------

findOverlaps(gr, gr1)
