---
title: "Statistical exploration of landscapes of phylogenetic trees"
author: "Thibaut Jombart, Michelle Kendall"
date: "2015-08-25"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{treescape: statistical exploration of landscapes of phylogenetic trees}
  \usepackage[utf8]{inputenc}
---




*treescape*: statistical exploration of landscapes of phylogenetic trees
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
## Loading required package: ape
## Loading required package: ade4
```


Content overview
-------------
The main functions implemented in *treescape* are:
* __`treescape`__: 
* __`treescapeServer`__: 
* __`findGroves`__: 
* __`plotGroves`__: 
* __`tree.dist`__: 
* __`med.tree`__: 


Distributed datasets include:
* __`woodmiceTrees`__: 



Authors / Contributors
--------------
Authors:
* [Thibaut Jombart](https://sites.google.com/site/thibautjombart/)
* [Michelle Kendall](http://www.imperial.ac.uk/people/m.kendall)

Contributors:
* [Jacob Almagro Garcia](http://www.well.ox.ac.uk/jacob-almagro-garcia)

Maintainer of the CRAN version:
* [Michelle Kendall](http://www.imperial.ac.uk/people/m.kendall)
