[![Travis-CI Build Status](https://travis-ci.org/thibautjombart/treescape.png?branch=master)](https://travis-ci.org/thibautjombart/treescape)



```
## Error in parse(text = x, srcfile = src): <text>:11:32: unexpected '='
## 10: ## plot results
## 11: plotGroves(wm.res$pco, lab=show=
##                                    ^
```


*treescape*: exploration of landscapes of phylogenetic trees
=================================================
*treescape* implements new methods for the exploration and analysis of distributions of phylogenetic trees for a given set of taxa.


Installing *treescape*
-------------
To install the development version from github:

```r
library(devtools)
install_github("thibautjombart/treescape")
```

The stable version can be installed from CRAN using:

```r
install.packages("treescape")
```

Then, to load the package, use:

```r
library("treescape")
```


Content overview
-------------
The main functions implemented in *treescape* are:
* __`treescape`__: explore landscapes of phylogenetic trees
* __`treescapeServer`__: open up an application in a web browser for an interactive exploration of the diversity in a set of trees
* __`findGroves`__: identify clusters of similar trees
* __`plotGroves`__: scatterplot of groups of trees
* __`treeVec`__: characterise a tree by a vector
* __`treeDist`__: find the distance between two tree vectors
* __`multiDist`__: find the pairwise distances of a list of trees
* __`medTree`__: find geometric median tree(s) to summarise a group of trees


Distributed datasets include:
* __`woodmiceTrees`__: illustrative set of 201 trees built using the neighbour-joining and bootstrapping example from the *woodmice* dataset in the *ape* documentation.


Exploring trees with *treescape*
--------------

We first load *treescape*, and packages required for graphics:

```r
library("treescape")
library("ade4")
library("adegenet")
```

```
## 
##    /// adegenet 2.0.1 is loaded ////////////
## 
##    > overview: '?adegenet'
##    > tutorials/doc/questions: 'adegenetWeb()' 
##    > bug reports/feature resquests: adegenetIssues()
## 
## 
## 
## Attaching package: 'adegenet'
## 
## The following object is masked from 'package:treescape':
## 
##     .render.server.info
```

```r
library("adegraphics")
```

```
## 
## Attaching package: 'adegraphics'
## 
## The following objects are masked from 'package:ade4':
## 
##     kplotsepan.coa, s.arrow, s.class, s.corcircle, s.distri,
##     s.image, s.label, s.logo, s.match, s.traject, s.value,
##     table.value, triangle.class
## 
## The following object is masked from 'package:ape':
## 
##     zoom
```

```r
library("ggplot")
```

```
## Error in library("ggplot"): there is no package called 'ggplot'
```

The function __`treescape`__ defines typologies of phylogenetic trees using a two-steps approach:
1. perform pairwise comparisons of trees using various (Euclidean) metrics; by default, comparison uses the Kendall and Colijn metric (Kendall & Colijn, submitted) which is described in more detail below; other metrics rely on tips distances implemented in *adephylo* (Jombart *et al.* 2010).
2. use Metric Multidimensional Scaling (MDS, aka Principal Coordinates Analysis, PCoA) to summarise pairwise distances between the trees as well as possible into a few dimensions; output of MDS is typically visualised using scatterplots of the first few Principal Components (PCs); this step relies on the PCoA implemented in *ade4* (Dray & Dufour 2007).

The function `treescape` performs both tasks, returning both the matrix of pairwise tree comparisons (`$D`), and the PCoA (`$pco`).
This can be illustrated using randomly generated trees:

```r
## generate list of trees
x <- rmtree(10, 20)
names(x) <- paste("tree", 1:10, sep = "")

## use treescape
res <- treescape(x, nf=3)
names(res)
```

```
## [1] "D"   "pco"
```

```r
res
```

```
## $D
##        tree1 tree2 tree3 tree4 tree5 tree6 tree7 tree8 tree9
## tree2  33.00                                                
## tree3  41.40 34.01                                          
## tree4  39.08 30.40 36.21                                    
## tree5  33.24 32.80 38.22 34.70                              
## tree6  33.94 26.29 35.67 29.24 32.02                        
## tree7  35.69 26.78 37.82 32.11 31.35 29.36                  
## tree8  38.55 31.13 38.63 33.24 37.35 30.00 31.43            
## tree9  37.07 29.92 38.65 35.40 35.96 31.62 28.74 33.62      
## tree10 36.40 28.21 38.82 33.65 35.75 29.24 29.55 31.89 30.68
## 
## $pco
## Duality diagramm
## class: pco dudi
## $call: dudi.pco(d = D, scannf = is.null(nf), nf = nf)
## 
## $nf: 3 axis-components saved
## $rank: 9
## eigen values: 93.73 91.12 65.21 59.6 50.85 ...
##   vector length mode    content       
## 1 $cw    9      numeric column weights
## 2 $lw    10     numeric row weights   
## 3 $eig   9      numeric eigen values  
## 
##   data.frame nrow ncol content             
## 1 $tab       10   9    modified array      
## 2 $li        10   3    row coordinates     
## 3 $l1        10   3    row normed scores   
## 4 $co        9    3    column coordinates  
## 5 $c1        9    3    column normed scores
## other elements: NULL
```

Pairwise distances can be visualised using *adegraphics*:

```r
## table.image
table.image(res$D, nclass=30)
```

![plot of chunk distances](vignettes/figs/distances-1.png) 

```r
## table.value
table.value(res$D, nclass=10)
```

![plot of chunk distances](vignettes/figs/distances-2.png) 

```r
## with some customization
table.value(res$D, nclass=5, method="color" , symbol="circle", col=heat.colors(10))
```

![plot of chunk distances](vignettes/figs/distances-3.png) 

The best representation of these distances in a 2-dimensional space is given by the first 2 PCs of the MDS.
These can be visualised using *adegraphics*'s function `scatter`:

```r
scatter(res$pco)
```

![plot of chunk treescapescatter](vignettes/figs/treescapescatter-1.png) 
Alternatively, the function `plotGroves` can be used:

```r
plotGroves(res$pco, lab.show=TRUE, lab.cex=1.5)
```

![plot of chunk plotGroves](vignettes/figs/plotGroves-1.png) 


`treecsape` can be furthe illustrated using *ape*'s dataset *woodmouse*, from which we built the 201 trees supplied in __`woodmiceTrees`__ using the neighbour-joining and bootstrapping example from the *ape* documentation. 

