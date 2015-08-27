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
* __`medTree`__: find geometric median tree(s) to summarise a group of trees

Other functions are central to the computations of distances between trees:
* __`treeVec`__: characterise a tree by a vector
* __`treeDist`__: find the distance between two tree vectors
* __`multiDist`__: find the pairwise distances of a list of trees


Distributed datasets include:
* __`woodmiceTrees`__: illustrative set of 201 trees built using the neighbour-joining and bootstrapping example from the *woodmice* dataset in the *ape* documentation.


Exploring trees with *treescape*
--------------

We first load *treescape*, and packages required for graphics:

```r
library("treescape")
library("ade4")
library("adegenet")
library("adegraphics")
library("ggplot2")
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
## tree2  31.27                                                
## tree3  27.44 34.10                                          
## tree4  25.08 31.86 28.43                                    
## tree5  41.94 36.10 42.17 42.26                              
## tree6  28.64 30.85 29.55 28.76 41.84                        
## tree7  35.03 29.85 39.67 36.58 40.45 38.43                  
## tree8  25.51 30.38 25.96 25.34 36.66 29.48 33.91            
## tree9  22.54 27.20 27.55 23.22 40.71 27.13 34.28 25.51      
## tree10 30.30 29.77 32.70 31.00 39.66 30.46 31.51 30.15 27.20
## 
## $pco
## Duality diagramm
## class: pco dudi
## $call: dudi.pco(d = D, scannf = is.null(nf), nf = nf)
## 
## $nf: 3 axis-components saved
## $rank: 9
## eigen values: 132.3 86.92 53.57 44.51 41.1 ...
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

<img src="vignettes/figs/distances-1.png" title="plot of chunk distances" alt="plot of chunk distances" width="800px" />

```r
## table.value
table.value(res$D, nclass=10)
```

<img src="vignettes/figs/distances-2.png" title="plot of chunk distances" alt="plot of chunk distances" width="800px" />

```r
## with some customization
table.value(res$D, nclass=5, method="color" , symbol="circle", col=heat.colors(10))
```

<img src="vignettes/figs/distances-3.png" title="plot of chunk distances" alt="plot of chunk distances" width="800px" />

The best representation of these distances in a 2-dimensional space is given by the first 2 PCs of the MDS.
These can be visualised using *adegraphics*'s function `scatter`:

```r
scatter(res$pco)
```

<img src="vignettes/figs/treescapescatter-1.png" title="plot of chunk treescapescatter" alt="plot of chunk treescapescatter" width="800px" />

Alternatively, the function `plotGroves` can be used:

```r
plotGroves(res$pco, lab.show=TRUE, lab.cex=1.5)
```

<img src="vignettes/figs/plotGroves-1.png" title="plot of chunk plotGroves" alt="plot of chunk plotGroves" width="800px" />


`treecsape` can be furthe illustrated using *ape*'s dataset *woodmouse*, from which we built the 201 trees supplied in __`woodmiceTrees`__ using the neighbour-joining and bootstrapping example from the *ape* documentation. 

