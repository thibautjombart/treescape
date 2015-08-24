---
title: "Statistical exploration of landscapes of phylogenetic trees"
author: "Thibaut Jombart, Michelle Kendall"
date: "2015-08-24"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{treescape: statistical exploration of landscapes of phylogenetic trees}
  \usepackage[utf8]{inputenc}
---




*treescape*: tatistical exploration of landscapes of phylogenetic trees
=================================================
*treescape* implements new methods for the exploration and analysis of distributions of phylogenetic trees for a given set of taxa.


Installing *apex*
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

