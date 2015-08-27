[![Travis-CI Build Status](https://travis-ci.org/thibautjombart/treescape.png?branch=master)](https://travis-ci.org/thibautjombart/treescape)





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
## tree2  27.53                                                
## tree3  29.85 29.95                                          
## tree4  43.76 44.89 43.91                                    
## tree5  34.68 32.73 33.23 38.29                              
## tree6  31.27 34.67 34.51 42.13 36.59                        
## tree7  33.03 35.23 31.65 36.74 36.69 29.92                  
## tree8  26.19 27.60 28.09 42.77 32.88 32.19 34.68            
## tree9  28.09 31.22 32.16 40.52 36.50 30.81 31.59 27.29      
## tree10 30.59 33.88 32.23 39.46 37.24 32.86 31.42 27.64 21.93
## 
## $pco
## Duality diagramm
## class: pco dudi
## $call: dudi.pco(d = D, scannf = is.null(nf), nf = nf)
## 
## $nf: 3 axis-components saved
## $rank: 9
## eigen values: 136.7 91.33 68.85 54.05 45.94 ...
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

<img src="vignettes/figs/distances-1.png" title="plot of chunk distances" alt="plot of chunk distances" width="400px" />

```r
## table.value
table.value(res$D, nclass=10)
```

<img src="vignettes/figs/distances-2.png" title="plot of chunk distances" alt="plot of chunk distances" width="400px" />

```r
## with some customization
table.value(res$D, nclass=5, method="color" , symbol="circle", col=heat.colors(10))
```

<img src="vignettes/figs/distances-3.png" title="plot of chunk distances" alt="plot of chunk distances" width="400px" />

The best representation of these distances in a 2-dimensional space is given by the first 2 PCs of the MDS.
These can be visualised using *adegraphics*'s function `scatter`:

```r
scatter(res$pco)
```

<img src="vignettes/figs/treescapescatter-1.png" title="plot of chunk treescapescatter" alt="plot of chunk treescapescatter" width="400px" />

Alternatively, the function `plotGroves` can be used:

```r
plotGroves(res$pco, lab.show=TRUE, lab.cex=1.5)
```

<img src="vignettes/figs/plotGroves-1.png" title="plot of chunk plotGroves" alt="plot of chunk plotGroves" width="400px" />

`treecsape` can be further illustrated using *ape*'s dataset *woodmouse*, from which we built the 201 trees supplied in __`woodmiceTrees`__ using the neighbour-joining and bootstrapping example from the *ape* documentation. 

```r
data(woodmiceTrees)
wm.res <- treescape(woodmiceTrees,nf=3)

## this is the PCoA / MDS:
wm.res$pco
```

```
## Duality diagramm
## class: pco dudi
## $call: dudi.pco(d = D, scannf = is.null(nf), nf = nf)
## 
## $nf: 3 axis-components saved
## $rank: 54
## eigen values: 32.69 24.41 6.952 6.348 4.363 ...
##   vector length mode    content       
## 1 $cw    54     numeric column weights
## 2 $lw    201    numeric row weights   
## 3 $eig   54     numeric eigen values  
## 
##   data.frame nrow ncol content             
## 1 $tab       201  54   modified array      
## 2 $li        201  3    row coordinates     
## 3 $l1        201  3    row normed scores   
## 4 $co        54   3    column coordinates  
## 5 $c1        54   3    column normed scores
## other elements: NULL
```

```r
## PCs are stored in:
head(wm.res$pco$li)
```

```
##         A1     A2      A3
## 1  -0.9949 -1.363 -0.7918
## 2  -0.6137 -1.014 -0.6798
## 3   2.6667  4.219 -2.9293
## 4 -13.6081  1.854  1.0947
## 5   2.1980  4.176 -3.1960
## 6   3.6013  4.865  2.9853
```

```r
## plot results
plotGroves(wm.res$pco, lab.show=TRUE, lab.optim=FALSE)
```

<img src="vignettes/figs/woodmicePlots-1.png" title="plot of chunk woodmicePlots" alt="plot of chunk woodmicePlots" width="400px" />

```r
## visualising density of points
s.kde2d(wm.res$pco$li)
```

<img src="vignettes/figs/woodmicePlots-2.png" title="plot of chunk woodmicePlots" alt="plot of chunk woodmicePlots" width="400px" />

```r
## alternative visualisation
s.density(wm.res$pco$li, col=redpal(100), bandwidth=3)
```

<img src="vignettes/figs/woodmicePlots-3.png" title="plot of chunk woodmicePlots" alt="plot of chunk woodmicePlots" width="400px" />

```r
## same, other palette
s.density(wm.res$pco$li, col=rev(transp(spectral(100),.5)), bandwidth=3)
```

<img src="vignettes/figs/woodmicePlots-4.png" title="plot of chunk woodmicePlots" alt="plot of chunk woodmicePlots" width="400px" />

```r
## alternative using ggplot2
woodmiceplot <- ggplot(wm.res$pco$li, aes(x=A1, y=A2)) # create plot
woodmiceplot + geom_density2d(colour="gray80") + # contour lines
geom_point(size=6, shape=1, colour="gray50") + # grey edges
geom_point(size=6, alpha=0.2, colour="navy") + # transparent blue points
xlab("") + ylab("") + theme_bw(base_family="") # remove axis labels and grey background
```

<img src="vignettes/figs/woodmicePlots-5.png" title="plot of chunk woodmicePlots" alt="plot of chunk woodmicePlots" width="400px" />


Note that alternatively, the function __`multiDist`__ simply performs the pairwise comparison of trees and outputs a distance matrix. 
This function may be preferable for large datasets, and when principal co-ordinate analysis is not required. 
It includes an option to save memory at the expense of computation time.




Identifying clusters of trees
--------------
Once a typology of trees has been derived using the approach described above, one may want to formally identify clusters of similar trees.
One simple approach is:
1. select a few first PCs of the MDS (retaining signal but getting rid of random noise)
2. derive pairwise Euclidean distances between trees based on these PCs
3. use hierarchical clustering to obtain a dendrogram of these trees
4. cut the dendrogram to obtain clusters
 
In *treescape*, the function __`findGroves`__ implements this approach, offering various clustering options (see `?findGroves`):

```r
wm.groves <- findGroves(woodmiceTrees, nf=3, nclust=6)
names(wm.groves)
```

```
## [1] "groups"    "treescape"
```
Note that when the number of clusters (`nclust`) is not provided, the function will display a dendrogram and ask for a cut-off height. 

The results can be plotted directly using `plotGroves` (see `?plotGroves` for options):

```r
## basic plot
plotGroves(wm.groves)
```

<img src="vignettes/figs/plotgroves-1.png" title="plot of chunk plotgroves" alt="plot of chunk plotgroves" width="400px" />

```r
## alternative with inertia ellipses
plotGroves(wm.groves, type="ellipse")
```

<img src="vignettes/figs/plotgroves-2.png" title="plot of chunk plotgroves" alt="plot of chunk plotgroves" width="400px" />

```r
## plot axes 2-3
plotGroves(wm.groves, xax=2, yax=3)
```

<img src="vignettes/figs/plotgroves-3.png" title="plot of chunk plotgroves" alt="plot of chunk plotgroves" width="400px" />

```r
## customize graphics
plotGroves(wm.groves, bg="black", col.pal=lightseasun, lab.show=TRUE, lab.col="white", lab.cex=1.5)
```

<img src="vignettes/figs/plotgroves-4.png" title="plot of chunk plotgroves" alt="plot of chunk plotgroves" width="400px" /><img src="vignettes/figs/plotgroves-5.png" title="plot of chunk plotgroves" alt="plot of chunk plotgroves" width="400px" />



`treescapeServer`: a web application for *treescape*
--------------
The essential functionalities of `treescape` are also available via a user-friendly web interface, running locally on the default web browser.
It can be started using by simply typing `treescapeServer()`.
The interface allows one to import data, run `treescape` to explore the tree space, look for clusters using `findGroves`, customize MDS plots, visualise specific trees and saves results in various formats.
It is fully documented in the *help* tab.

![example of treescapeServer running](vignettes/figs/server.png) 




Finding median trees
--------------

When a set of trees have very similar structures, it makes sense to summarize them into as single tree, 'consensus' tree.
In `treescape`, this is achieved by finding the *median tree* is defined for a set of trees, as the tree which is closest to the centre of the set of trees in the tree landscape defined by `treescape`.
This procedure is implemented by the function **`medTree`**:

```r
## get first median tree
tre <- medTree(woodmiceTrees)$trees[[1]]

## plot tree
plot(tre,type="cladogram",edge.width=3, cex=0.8)
```

<img src="vignettes/figs/woodmiceMedian-1.png" title="plot of chunk woodmiceMedian" alt="plot of chunk woodmiceMedian" width="400px" />

However, a more complete and accurate summary of the data can be given by finding a summary tree from each cluster.
This is achieved using the `groups` argument of `medTree`:

```r
## identify 6 clusters
groves <- findGroves(woodmiceTrees, nf=3, nclust=6)

## find median trees
res <- medTree(woodmiceTrees, groves$groups)

## there isone output per cluster
names(res)
```

```
## [1] "1" "2" "3" "4" "5" "6"
```

```r
## get the first median of each
med.trees <- lapply(res, function(e) ladderize(e$trees[[1]]))

## plot trees
par(mfrow=c(2,3))
for(i in 1:length(med.trees)) plot(med.trees[[i]], main=paste("cluster",i))
```

<img src="vignettes/figs/woodmiceCluster1-1.png" title="plot of chunk woodmiceCluster1" alt="plot of chunk woodmiceCluster1" width="400px" />

These trees exhibit a number of topological differences, e.g. in the placement of the **(1007S,1208S,0909S)** clade. 
Performing this analysis enables the detection of distinct representative trees supported by data.



Method: characterising a tree by a vector
--------------
Kendall and Colijn proposed a [metric](http://arxiv.org/abs/1507.05211) for comparing rooted phylogenetic trees. Each tree is characterised by a vector which notes the placement of the most recent common ancestor (MRCA) of each pair of tips. Specifically, it records the distance between the MRCA of a pair of tips *(i,j)* and the root in two ways: the number of edges *m(i,j)*, and the path length *M(i,j)*. It also records the length *p(i)* of each 'pendant' edge between a tip *i* and its immediate ancestor. This procedure results in two vectors for a tree *T*:

*m(T) = (m(1,2), m(1,3),...,m(k-1,k),1,...,1)*

and

*M(T) = (M(1,2), M(1,3),...,M(k-1,k),p(1),...,p(k)).*

In *m(T)* we record the pendant lengths as 1, as each tip is 1 step from its immediate ancestor. We combine *m* and *M* with a parameter lambda between zero and one to weight the contribution of branch lengths, characterising each tree with a vector 

*v{lambda}(T) = (1-lambda)m(T) + lambda M(T)*.

This is implemented as the function __`treeVec`__. For example,

```r
## generate a random tree:
tree <- rtree(6)
## topological vector of mrca distances from root:
treeVec(tree)
```

```
##  [1] 1 1 1 0 0 2 3 0 0 2 0 0 0 0 1 1 1 1 1 1 1
```

```r
## vector of mrca distances from root when lambda=0.5:
treeVec(tree,0.5)
```

```
##  [1] 0.7473 0.7473 0.7473 0.0000 0.0000 1.3366 2.2140 0.0000 0.0000 1.3366
## [11] 0.0000 0.0000 0.0000 0.0000 0.8248 0.5023 0.6922 0.7611 0.9299 0.9121
## [21] 0.9191
```

```r
## vector of mrca distances as a function of lambda:
vecAsFunction <- treeVec(tree,return_lambda_function=TRUE)
```

```
## Error in treeVec(tree, return_lambda_function = TRUE): unused argument (return_lambda_function = TRUE)
```

```r
## evaluate the vector at lambda=0.5:
vecAsFunction(0.5)
```

```
## Error in eval(expr, envir, enclos): could not find function "vecAsFunction"
```

The metric -- the distance between two trees -- is the Euclidean distance between these vectors:

*d{lambda}(Ta, Tb) = || v{lambda}(Ta) - v{lambda}(Tb) ||.*


This can be found using __`treeDist`__:

```r
## generate random trees
tree_a <- rtree(6)
tree_b <- rtree(6)
## topological (lambda=0) distance:
treeDist(tree_a,tree_b) 
```

```
## [1] 5.196
```

```r
## branch-length focused (lambda=1) distance:
treeDist(tree_a,tree_b,1)
```

```
## [1] 2.197
```



References
--------------
* Kendall M & Colijn C (submitted) ...
* Dray S & Dufour AB (2007): The ade4 package: implementing the duality diagram for ecologists. Journal of Statistical Software 22(4): 1-20.
* Jombart R, Balloux F & Dray S (2010) adephylo: new tools for investigating the phylogenetic signal in biological traits. Bioinformatics 26: 1907-1909. Doi: 10.1093/bioinformatics/btq292




Authors / Contributors
--------------
Authors:
* [Thibaut Jombart](https://sites.google.com/site/thibautjombart/)
* [Michelle Kendall](http://www.imperial.ac.uk/people/m.kendall)
* [Jacob Almagro Garcia](http://www.well.ox.ac.uk/jacob-almagro-garcia)

Contributors:
* [Caroline Colijn](http://www.imperial.ac.uk/people/c.colijn)

Maintainer of the CRAN version:
* [Michelle Kendall](http://www.imperial.ac.uk/people/m.kendall)

