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

```
## Error in library("treescape"): there is no package called 'treescape'
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
```

```
## Error in library("treescape"): there is no package called 'treescape'
```

```r
library("ade4")
library("adegenet")
```

```
## 
##    /// adegenet 2.0.0 is loaded ////////////
## 
##    > overview: '?adegenet'
##    > tutorials/doc/questions: 'adegenetWeb()' 
##    > bug reports/feature resquests: adegenetIssues()
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
```

```r
library("ggplot2")
```

The function `treescape` defines typologies of phylogenetic trees using a two-step approach:

1. perform pairwise comparisons of trees using various (Euclidean) metrics; by default, the comparison uses the Kendall and Colijn metric (Kendall & Colijn, 2015) which is described in more detail below; other metrics rely on tips distances implemented in *adephylo* (Jombart *et al.* 2010).

2. use Metric Multidimensional Scaling (MDS, aka Principal Coordinates Analysis, PCoA) to summarise pairwise distances between the trees as well as possible into a few dimensions; the output of the MDS is typically visualised using scatterplots of the first few Principal Components (PCs); this step relies on the PCoA implemented in *ade4* (Dray & Dufour 2007).

The function `treescape` performs both tasks, returning both the matrix of pairwise tree comparisons (`$D`), and the PCoA (`$pco`).
This can be illustrated using randomly generated trees:

```r
## generate list of trees
set.seed(1)
x <- rmtree(10, 20)
```

```
## Error in eval(expr, envir, enclos): could not find function "rmtree"
```

```r
names(x) <- paste("tree", 1:10, sep = "")
```

```
## Error in names(x) <- paste("tree", 1:10, sep = ""): object 'x' not found
```

```r
## use treescape
res <- treescape(x, nf=3)
```

```
## Error in eval(expr, envir, enclos): could not find function "treescape"
```

```r
names(res)
```

```
## Error in eval(expr, envir, enclos): object 'res' not found
```

```r
res
```

```
## Error in eval(expr, envir, enclos): object 'res' not found
```

Pairwise distances can be visualised using *adegraphics*:

```r
## table.image
table.image(res$D, nclass=30)
```

```
## Error in eval(expr, envir, enclos): object 'res' not found
```

```r
## table.value with some customization
table.value(res$D, nclass=5, method="color", 
            symbol="circle", col=redpal(5))
```

```
## Error in eval(expr, envir, enclos): object 'res' not found
```

The best representation of these distances in a 2-dimensional space is given by the first 2 PCs of the MDS.
These can be visualised using *adegraphics*'s function `scatter`:

```r
scatter(res$pco)
```

```
## Error in scatter(res$pco): object 'res' not found
```

Alternatively, the function `plotGroves` can be used:

```r
plotGroves(res$pco, lab.show=TRUE, lab.cex=1.5)
```

```
## Error in eval(expr, envir, enclos): could not find function "plotGroves"
```

The functionality of `treecsape` can be further illustrated using *ape*'s dataset *woodmouse*, from which we built the 201 trees supplied in `woodmiceTrees` using the neighbour-joining and bootstrapping example from the *ape* documentation. 

```r
data(woodmiceTrees)
```

```
## Warning in data(woodmiceTrees): data set 'woodmiceTrees' not found
```

```r
wm.res <- treescape(woodmiceTrees,nf=3)
```

```
## Error in eval(expr, envir, enclos): could not find function "treescape"
```

```r
## this is the PCoA / MDS:
wm.res$pco
```

```
## Error in eval(expr, envir, enclos): object 'wm.res' not found
```

```r
## PCs are stored in:
head(wm.res$pco$li)
```

```
## Error in head(wm.res$pco$li): object 'wm.res' not found
```

```r
## plot results
plotGroves(wm.res$pco, lab.show=TRUE, lab.optim=FALSE)
```

```
## Error in eval(expr, envir, enclos): could not find function "plotGroves"
```

```r
## visualising density of points
s.kde2d(wm.res$pco$li)
```

```
## Error in data.frame(dfxy): object 'wm.res' not found
```

```r
## alternative visualisation
s.density(wm.res$pco$li, col=redpal(100), bandwidth=3)
```

```
## Error in s.density(wm.res$pco$li, col = redpal(100), bandwidth = 3): non convenient selection for dfxy (can not be converted to dataframe)
```

```r
## same, other palette
s.density(wm.res$pco$li, col=rev(transp(spectral(100),.5)), bandwidth=3)
```

```
## Error in s.density(wm.res$pco$li, col = rev(transp(spectral(100), 0.5)), : non convenient selection for dfxy (can not be converted to dataframe)
```

```r
## alternative using ggplot2
woodmiceplot <- ggplot(wm.res$pco$li, aes(x=A1, y=A2)) # create plot
```

```
## Error in ggplot(wm.res$pco$li, aes(x = A1, y = A2)): object 'wm.res' not found
```

```r
woodmiceplot + geom_density2d(colour="gray80") + # contour lines
geom_point(size=6, shape=1, colour="gray50") + # grey edges
geom_point(size=6, alpha=0.2, colour="navy") + # transparent blue points
xlab("") + ylab("") + theme_bw(base_family="") # remove axis labels and grey background
```

```
## Error in eval(expr, envir, enclos): object 'woodmiceplot' not found
```


Note that alternatively, the function `multiDist` simply performs the pairwise comparison of trees and outputs a distance matrix. 
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
 
In *treescape*, the function `findGroves` implements this approach, offering various clustering options (see `?findGroves`):

```r
wm.groves <- findGroves(woodmiceTrees, nf=3, nclust=6)
```

```
## Error in eval(expr, envir, enclos): could not find function "findGroves"
```

```r
names(wm.groves)
```

```
## Error in eval(expr, envir, enclos): object 'wm.groves' not found
```
Note that when the number of clusters (`nclust`) is not provided, the function will display a dendrogram and ask for a cut-off height. 

The results can be plotted directly using `plotGroves` (see `?plotGroves` for options):

```r
## basic plot
plotGroves(wm.groves)
```

```
## Error in eval(expr, envir, enclos): could not find function "plotGroves"
```

```r
## alternative with inertia ellipses
plotGroves(wm.groves, type="ellipse")
```

```
## Error in eval(expr, envir, enclos): could not find function "plotGroves"
```

```r
## plot axes 2-3
plotGroves(wm.groves, xax=2, yax=3)
```

```
## Error in eval(expr, envir, enclos): could not find function "plotGroves"
```

```r
## customize graphics
plotGroves(wm.groves, bg="black", col.pal=lightseasun, lab.show=TRUE, lab.col="white", lab.cex=1.5)
```

```
## Error in eval(expr, envir, enclos): could not find function "plotGroves"
```



`treescapeServer`: a web application for *treescape*
--------------
The essential functionalities of `treescape` are also available via a user-friendly web interface, running locally on the default web browser.
It can be started by simply typing `treescapeServer()`.
The interface allows one to import data, run `treescape` to explore the tree space, look for clusters using `findGroves`, customize MDS plots, visualise specific trees and save results in various formats.
It is fully documented in the *help* tab.

![example of treescapeServer running](vignettes/figs/server.png) 




Finding median trees
--------------

When a set of trees have very similar structures, it makes sense to summarize them into a single 'consensus' tree.
In `treescape`, this is achieved by finding the *median tree* for a set of trees according to the Kendall and Colijn metric.
That is, we find the tree which is closest to the centre of the set of trees in the tree landscape defined in `treescape`.
This procedure is implemented by the function **`medTree`**:

```r
## get first median tree
tre <- medTree(woodmiceTrees)$trees[[1]]
```

```
## Error in eval(expr, envir, enclos): could not find function "medTree"
```

```r
## plot tree
plot(tre,type="cladogram",edge.width=3, cex=0.8)
```

```
## Error in plot(tre, type = "cladogram", edge.width = 3, cex = 0.8): error in evaluating the argument 'x' in selecting a method for function 'plot': Error: object 'tre' not found
```

However, a more complete and accurate summary of the data can be given by finding a summary tree from each cluster.
This is achieved using the `groups` argument of `medTree`:

```r
## identify 6 clusters
groves <- findGroves(woodmiceTrees, nf=3, nclust=6)
```

```
## Error in eval(expr, envir, enclos): could not find function "findGroves"
```

```r
## find median trees
res <- medTree(woodmiceTrees, groves$groups)
```

```
## Error in eval(expr, envir, enclos): could not find function "medTree"
```

```r
## there isone output per cluster
names(res)
```

```
## Error in eval(expr, envir, enclos): object 'res' not found
```

```r
## get the first median of each
med.trees <- lapply(res, function(e) ladderize(e$trees[[1]]))
```

```
## Error in lapply(res, function(e) ladderize(e$trees[[1]])): object 'res' not found
```

```r
## plot trees
par(mfrow=c(2,3))
for(i in 1:length(med.trees)) plot(med.trees[[i]], main=paste("cluster",i),cex=1.5)
```

```
## Error in eval(expr, envir, enclos): object 'med.trees' not found
```

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
```

```
## Error in eval(expr, envir, enclos): could not find function "rtree"
```

```r
## topological vector of mrca distances from root:
treeVec(tree)
```

```
## Error in eval(expr, envir, enclos): could not find function "treeVec"
```

```r
## vector of mrca distances from root when lambda=0.5:
treeVec(tree,0.5)
```

```
## Error in eval(expr, envir, enclos): could not find function "treeVec"
```

```r
## vector of mrca distances as a function of lambda:
vecAsFunction <- treeVec(tree,return.lambda.function=TRUE)
```

```
## Error in eval(expr, envir, enclos): could not find function "treeVec"
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
```

```
## Error in eval(expr, envir, enclos): could not find function "rtree"
```

```r
tree_b <- rtree(6)
```

```
## Error in eval(expr, envir, enclos): could not find function "rtree"
```

```r
## topological (lambda=0) distance:
treeDist(tree_a,tree_b) 
```

```
## Error in eval(expr, envir, enclos): could not find function "treeDist"
```

```r
## branch-length focused (lambda=1) distance:
treeDist(tree_a,tree_b,1)
```

```
## Error in eval(expr, envir, enclos): could not find function "treeDist"
```



References
--------------
* Dray S & Dufour AB (2007): The ade4 package: implementing the duality diagram for ecologists. Journal of Statistical Software 22(4): 1-20.
* Jombart R, Balloux F & Dray S (2010) adephylo: new tools for investigating the phylogenetic signal in biological traits. Bioinformatics 26: 1907-1909. Doi: 10.1093/bioinformatics/btq292
* Kendall M & Colijn C (Preprint 2015) A tree metric using structure and length to capture distinct phylogenetic signals. arXiv 1507.05211




Authors / Contributors
--------------
Authors:
* [Thibaut Jombart](https://sites.google.com/site/thibautjombart/)
* [Michelle Kendall](http://www.imperial.ac.uk/people/m.kendall)
* [Jacob Almagro-Garcia](http://www.well.ox.ac.uk/jacob-almagro-garcia)

Contributors:
* [Caroline Colijn](http://www.imperial.ac.uk/people/c.colijn)

Maintainer of the CRAN version:
* [Michelle Kendall](http://www.imperial.ac.uk/people/m.kendall)

