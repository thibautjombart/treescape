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
* __`plotGroves`__: scatterplot of groups of trees, and __`plotGrovesD3`__ which enables interactive plotting based on d3.js
* __`medTree`__: find geometric median tree(s) to summarise a group of trees

Other functions are central to the computations of distances between trees:
* __`treeVec`__: characterise a tree by a vector
* __`treeDist`__: find the distance between two tree vectors
* __`multiDist`__: find the pairwise distances of a list of trees
* __`refTreeDist`__: find the distances of a list of trees from a reference tree


Distributed datasets include:
* __`woodmiceTrees`__: illustrative set of 201 trees built using the neighbour-joining and bootstrapping example from the *woodmice* dataset in the *ape* documentation.
* __`DengueTrees`__: 500 trees sampled from a BEAST posterior set of trees from Drummond, A. J., and Rambaut, A. (2007) BEAST: Bayesian evolutionary analysis by sampling trees.



Exploring trees with *treescape*
--------------

We first load *treescape*, and the packages required for graphics:

```r
library("treescape")
library("ade4")
library("adegenet")
library("adegraphics")
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
## tree2  26.00                                                
## tree3  31.06 26.74                                          
## tree4  42.85 42.12 44.44                                    
## tree5  30.66 27.71 27.37 44.79                              
## tree6  36.50 31.18 30.18 41.81 31.59                        
## tree7  34.64 28.71 29.48 40.35 31.11 32.37                  
## tree8  28.97 26.29 24.45 43.74 23.47 30.41 29.00            
## tree9  29.63 27.42 27.48 45.61 26.31 30.89 29.77 24.60      
## tree10 34.87 30.00 29.44 44.97 34.06 31.05 34.41 31.54 32.59
## 
## $pco
## Duality diagramm
## class: pco dudi
## $call: dudi.pco(d = D, scannf = is.null(nf), nf = nf)
## 
## $nf: 3 axis-components saved
## $rank: 9
## eigen values: 142.1 76.52 62.69 49.88 41.07 ...
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
## table.value with some customization
table.value(res$D, nclass=5, method="color", 
            symbol="circle", col=redpal(5))
```

![plot of chunk distances](vignettes/figs/distances-2.png)

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

![plot of chunk plotgroves](vignettes/figs/plotgroves-1.png)

The functionality of `treecsape` can be further illustrated using *ape*'s dataset *woodmouse*, from which we built the 201 trees supplied in `woodmiceTrees` using the neighbour-joining and bootstrapping example from the *ape* documentation. 

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

![plot of chunk woodmicePlots](vignettes/figs/woodmicePlots-1.png)

```r
## visualising density of points
s.kde2d(wm.res$pco$li)
```

![plot of chunk woodmicePlots](vignettes/figs/woodmicePlots-2.png)

```r
## alternative visualisation
s.density(wm.res$pco$li, col=redpal(100), bandwidth=3)
```

![plot of chunk woodmicePlots](vignettes/figs/woodmicePlots-3.png)

```r
## same, other palette
s.density(wm.res$pco$li, col=rev(transp(spectral(100),.5)), bandwidth=3)
```

![plot of chunk woodmicePlots](vignettes/figs/woodmicePlots-4.png)

```r
## alternative using ggplot2
woodmiceplot <- ggplot(wm.res$pco$li, aes(x=A1, y=A2)) # create plot
woodmiceplot + geom_density2d(colour="gray80") + # contour lines
geom_point(size=6, shape=1, colour="gray50") + # grey edges
geom_point(size=6, alpha=0.2, colour="navy") + # transparent blue points
xlab("") + ylab("") + theme_bw(base_family="") # remove axis labels and grey background
```

![plot of chunk woodmicePlots](vignettes/figs/woodmicePlots-5.png)

```r
## alternative using plotGrovesD3
plotGrovesD3(wm.res$pco)
```

<!--html_preserve--><div id="htmlwidget-9880" style="width:504px;height:504px;" class="scatterD3 html-widget"></div>
<script type="application/json" data-for="htmlwidget-9880">{"x":{"data":{"x":[-0.994915614079129,-0.613678626917793,2.66670631307031,-13.6081153840016,2.19796507262781,3.60128203945996,-7.9346433024299,-0.601509902833136,-0.994915614079122,-0.601509902833136,-0.607150138418814,10.3943189332891,-8.22061660863031,-0.542706411478826,3.33541847941513,-0.0438280980188251,-0.16462933902331,2.79080622437025,10.5466992537476,-0.611741178555834,3.57629422667166,9.64506059773349,3.87015775914438,-0.926024164827251,-1.21908989674129,2.4723857807866,-8.86243680646928,-1.21908989674129,-14.5730155701834,-13.6024751484159,0.917863183050825,2.27894083741,2.28568544942178,-15.7748658745841,2.14920253942819,-0.790924629423664,12.7908135541976,0.944271832856848,-1.24839432652649,3.98230290310849,-1.01426319128175,2.01968513715778,2.67234654865598,3.60013557723756,-0.360155109563768,3.37184987522406,0.252595301625617,3.36472290920033,6.40461224612179,-0.994699490566317,-9.42122752144957,2.66670631307031,-9.4500778191334,-0.607150138418815,13.1790840072057,1.82079944018017,-1.24928257943979,-0.994915614079121,2.43713839500107,-13.044572686349,-15.9705489873276,-0.913373605366737,-13.6310955663374,3.64770068321681,0.6688659159456,8.1468271893242,13.4934309468058,4.04346690730909,2.2791569609228,-0.580632735597669,2.36612308170806,5.33313046919023,-0.994915614079121,11.0665629950681,-1.00840683218326,4.76538085706121,1.4438114263851,-0.171393662122702,2.02621362565675,-8.75385363504267,-0.858828154723405,-0.330268540830038,3.90067192022295,-0.0730717036886925,-16.5286674300406,-0.994915614079121,2.4723857807866,-9.65081856527717,-0.335908776415716,-0.796564865009344,-0.542706411478826,1.0492923610008,-1.21908989674129,-15.6745883154345,0.0979934237349436,-9.34792446440201,-0.437013152012128,2.27241234891103,-0.176892270892568,0.955121246501658,12.8975581653797,5.51089790555853,-10.4680320213257,-0.0399993481969892,9.52808714143766,3.58911331537531,4.42342442747863,0.04726565910757,-1.25417131433135,2.36112888973579,2.27894083741,-9.02693355729029,2.07898006954061,3.87668624764336,-6.4615827644661,2.03185386124243,1.82035414418588,-0.994915614079122,-1.01426319128175,3.37749011080974,0.199097799331354,-0.172281915035998,-1.00862295569607,4.15243988951517,-1.00840683218326,2.28458107299568,-0.607150138418815,-0.519052641259981,2.08128395756897,-16.5601837313279,2.46688265557635,-9.05280962299188,3.48076661737388,-0.995803866992419,-0.988171002067341,4.03546615887711,-0.994915614079122,-1.32470655220642,-0.989275378493443,3.56866136174658,-0.607150138418815,-7.66836564572728,-0.335908776415717,6.72470221778202,-0.27737989822945,-0.633026204120414,2.84248353506267,0.0979934237349434,2.92888779742133,3.20320883731572,-0.335908776415717,2.87234652407253,4.32578048491479,-1.25417131433135,2.83684329947699,-1.00751857926997,-0.915013358974151,3.57540597375836,-0.989275378493443,-0.707126558298265,-0.712766793883943,2.26633787221916,3.5888971918625,-1.01180587891182,-9.25652460111789,2.33076992578286,-0.165035101165792,3.30290284183796,-0.994699490566317,-16.2186567440298,3.5888971918625,-0.982530766481663,9.52808714143766,-0.606100942970156,2.73656237720091,3.46106473819439,2.27894083741,3.63680350037365,2.4723857807866,8.1468271893242,2.27894083741,-0.995803866992419,-9.25088436553222,0.859657106490925,0.774749830581256,0.04726565910757,-15.8890069786153,0.404970591334761,2.42557261667578,-0.613678626917791,-0.607013028043388,-0.919902093865712,2.28458107299568,-0.607150138418814,-0.765676710838192,4.1467996539295,-13.1809002016393,12.8975581653797,4.26288720150493,2.02621362565675,3.36367420397477],"y":[-1.36296405257431,-1.0143380748759,4.21882824322934,1.85395115102288,4.17637909087139,4.86461034763142,-16.009575686826,-1.02771387775047,-1.36296405257431,-1.02771387775047,-1.02766987185401,-0.50552388611458,-10.0515267421666,-1.07229905286387,4.15752547173998,2.76392139550187,2.59796872098458,4.93883384448699,-10.7625404554933,-2.21247979854658,4.87226502578944,-9.50887087450621,3.32153776588015,-0.842513536381463,-0.716141175163109,-7.16104731132083,5.75437287564331,-0.716141175163109,-1.7069313395939,1.85390714512642,4.8987315359459,3.88353406250903,5.61969323330277,0.344365695407378,-14.4171273048538,-0.501442487250261,-9.6000449103967,1.03726886777004,-0.33393376356846,3.46374535755798,-1.35535338031275,4.91950411344247,4.21878423733287,7.13159940553585,-0.00302341316305923,5.01212399682195,-13.0610240652768,3.77531806014533,7.01232473326826,0.386526915197534,3.75962479695718,4.21882824322933,-11.8666109777134,-1.02766987185401,-8.32290371040569,2.52086549795982,-0.320645972486813,-1.36296405257431,4.62733018857612,3.8619430148942,-1.58401124473219,-1.37092740080123,-1.87858215341395,3.08391000762073,-1.19130284875386,-0.939160305549966,-9.99053252948651,5.96476384086601,5.63302503028087,-1.03141128116272,3.87552670838565,7.5655093012157,-1.36296405257431,-11.1843323318248,0.394093581562642,6.98752111471013,1.37775884921857,1.22975877872164,4.90617231646436,6.55363566239471,-1.1898855504877,0.868592373260222,6.09901634763348,1.2290930557728,0.146376109093674,-1.36296405257431,-7.16104731132083,-11.8099177019773,0.868636379156678,-0.501398481353805,-1.07229905286387,1.43585335420713,-0.716141175163109,-1.76388587914128,-0.23437411826138,-3.78957669060228,0.645071817193465,3.89686585948713,2.27313496228741,1.1554629700066,-9.37652435432995,8.34207063093772,4.70275523659972,-5.51421132484612,-8.2299588005646,4.87798615050598,6.18643823118439,0.0191206331844138,-0.326994001640874,-7.3165426940803,3.88353406250902,4.08158718069937,-7.49629748614467,3.30820596890205,1.9385779426291,4.9061283105679,1.76288229528141,-1.36296405257431,-1.35535338031275,5.0120799909255,-0.010809556298168,1.24304656980329,-1.3553973862092,5.13648704660545,0.394093581562642,3.88349005661257,-1.02766987185401,-1.36622147874601,4.52661553061152,0.167362584229802,-7.16806857752342,4.10252964993904,-10.8783948870456,-1.34967626149267,0.37319511821943,3.41920418834103,-1.36296405257431,-0.871680563819035,-1.36300805847077,3.14939364607734,-1.02766987185401,-1.0051396096495,0.868636379156678,-2.10228430818339,0.0359583960557166,-1.00672740261434,5.89152592638035,-0.23437411826138,0.461240469186376,-9.87333975241743,0.868636379156678,-8.00759724325874,4.51980695910259,-0.326994001640875,5.8915699322768,0.380805790480995,-1.35124757466906,4.88555281687109,-1.36300805847077,-1.18325326640639,-1.18320926050994,5.62730390556433,3.12849518273413,-1.36839184423488,-11.4879553182351,4.29529384801837,-1.92883606503873,4.00946413995936,0.386526915197535,-0.0288756124843535,3.12849518273413,0.373151112322974,-8.2299588005646,-2.21252380444303,4.19946765707952,7.74251598858856,3.88353406250902,-10.4582816896664,-7.16104731132083,-0.939160305549964,3.88353406250903,-1.34967626149267,-11.4879993241316,0.931854402146932,10.0778445455006,0.0191206331844132,-1.59197459295911,1.73932530176124,-7.36117187509016,-1.01433807487591,-1.03473514395306,-1.35759560382313,3.88349005661257,-1.02766987185401,0.202085100647231,5.13653105250191,1.52547532730244,-9.37652435432995,8.19604761262773,4.90617231646436,7.45847251006686],"point_opacity":[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],"key_var":[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201]},"settings":{"labels_size":10,"point_size":64,"xlab":"Axis 1","ylab":"Axis 2","has_labels":false,"col_var":null,"col_lab":"groups","colors":null,"ellipses":false,"ellipses_data":[],"symbol_var":null,"symbol_lab":"symbol_var","size_var":null,"size_range":[10,300],"size_lab":"NULL","key_var":null,"type_var":null,"unit_circle":false,"has_color_var":false,"has_symbol_var":false,"has_size_var":false,"has_legend":false,"has_tooltips":true,"tooltip_text":null,"has_custom_tooltips":false,"fixed":false,"legend_width":150,"html_id":"scatterD3-KXCINBWK","xlim":null,"ylim":null,"lasso":false,"lasso_callback":null,"dom_id_reset_zoom":"scatterD3-reset-zoom","dom_id_svg_export":"scatterD3-svg-export","dom_id_lasso_toggle":"scatterD3-lasso-toggle","transitions":false,"hashes":[]}},"evals":[],"jsHooks":[]}</script><!--/html_preserve-->


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
 
In *treescape*, the function `findGroves` implements this approach, offering various clustering options (see `?findGroves`). Here we supply the function with our `treescape` output `wm.res` since we have already calculated it, but it is also possible to skip the steps above and directly supply `findGroves` with a multiPhylo list of trees.

```r
wm.groves <- findGroves(wm.res, nclust=6)
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

![plot of chunk plotgroves2](vignettes/figs/plotgroves2-1.png)

```r
## alternative with inertia ellipses
plotGroves(wm.groves, type="ellipse")
```

![plot of chunk plotgroves2](vignettes/figs/plotgroves2-2.png)

```r
## plot axes 2-3
plotGroves(wm.groves, xax=2, yax=3)
```

![plot of chunk plotgroves2](vignettes/figs/plotgroves2-3.png)

```r
## customize graphics
plotGroves(wm.groves, bg="black", col.pal=lightseasun, lab.show=TRUE, lab.col="white", lab.cex=1.5)
```

![plot of chunk plotgroves3](vignettes/figs/plotgroves3-1.png)


`treescapeServer`: a web application for *treescape*
--------------
The functionalities of `treescape` are also available via a user-friendly web interface, running locally on the default web browser.
It can be started by simply typing `treescapeServer()`.
The interface allows you to import trees and run `treescape` to view and explore the tree space in 2 or 3 dimensions.
It is then straightforward to analyse the tree space by varying lambda, looking for clusters using `findGroves` and saving results in various formats.
Individual trees can be easily viewed including median trees per cluster, and collections of trees can be seen together using `densiTree` from the package `phangorn`.
It is fully documented in the *help* tab.

![example of treescapeServer 3d](vignettes/figs/treescape3d.png) 

![example of treescapeServer tree](vignettes/figs/treescapeTree.png)

![example of treescapeServer densiTree](vignettes/figs/treescapeDensiTree.png)


Finding median trees
--------------

When a set of trees have very similar structures, it makes sense to summarize them into a single 'consensus' tree.
In `treescape`, this is achieved by finding the *median tree* for a set of trees according to the Kendall and Colijn metric.
That is, we find the tree which is closest to the centre of the set of trees in the tree landscape defined in `treescape`.
This procedure is implemented by the function **`medTree`**:

```r
## get first median tree
tre <- medTree(woodmiceTrees)$trees[[1]]

## plot tree
plot(tre,type="cladogram",edge.width=3, cex=0.8)
```

![plot of chunk woodmiceMedian](vignettes/figs/woodmiceMedian-1.png)

However, a more complete and accurate summary of the data can be given by finding a summary tree from each cluster.
This is achieved using the `groups` argument of `medTree`:

```r
## find median trees for the 6 clusters identified earlier:
res <- medTree(woodmiceTrees, wm.groves$groups)

## there is one output per cluster
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
for(i in 1:length(med.trees)) plot(med.trees[[i]], main=paste("cluster",i),cex=1.5)
```

<img src="vignettes/figs/woodmiceCluster1-1.png" title="plot of chunk woodmiceCluster1" alt="plot of chunk woodmiceCluster1" width="600px" />

These trees exhibit a number of topological differences, e.g. in the placement of the **(1007S,1208S,0909S)** clade. 
Performing this analysis enables the detection of distinct representative trees supported by data.

Note that we supplied the function `medTree` with the multiPhylo list of trees. A more computationally efficient process (at the expense of using more memory) is to use the option `return.tree.vectors` in the initial `treescape` call, and then supply these vectors directly to `medTree`.
In this case, the tree indices are returned by `medTree` but the trees are not (since they were not supplied).

Emphasising the placement of certain tips or clades
--------------

In some analyses it may be informative to emphasise the placement of particular tips or clades within a set of trees. This can be particularly useful in large trees where the study is focused on a smaller clade. Priority can be given to a list of tips using the argument `emphasise.tips`, whose corresponding values in the vector comparison will be given a weight of `emphasise.weight` times the others (the default is 2, i.e. twice the weight).

For example, if we wanted to emphasise where the woodmice trees agree and disagree on the placement of the **(1007S,1208S,0909S)** clade, we can simply emphasise that clade as follows: 

```r
wm3.res <- treescape(woodmiceTrees,nf=2,emphasise.tips=c("No1007S","No1208S","No0909S"),emphasise.weight=3)

## plot results
plotGroves(wm3.res$pco, lab.show=TRUE, lab.optim=FALSE)
```

![plot of chunk woodmice-tip-emphasis](vignettes/figs/woodmice-tip-emphasis-1.png)

It can be seen from the scale of the plot and the density of clustering that the trees are now separated into more distinct clusters.

```r
wm3.groves <- findGroves(woodmiceTrees,nf=3,nclust=6,emphasise.tips=c("No1007S","No1208S","No0909S"),emphasise.weight=3)
plotGroves(wm3.groves, type="ellipse")
```

![plot of chunk findgroves-with-emphasis](vignettes/figs/findgroves-with-emphasis-1.png)

Conversely, where the structure of a particular clade is not of interest (for example, lineages within an outgroup which was only included for rooting purposes), those tips can be given a weight less than 1 so as to give them less emphasis in the comparison. We note that although it is possible to give tips a weighting of 0, we advise caution with this as the underlying function will no longer be guaranteed to be a metric. That is, a distance of 0 between two trees will no longer necessarily imply that the trees are identical. In most cases it would be wiser to assign a very small weighting to tips which are not of interest.

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
##  [1] 0 0 1 1 1 1 0 0 0 0 0 0 3 2 2 1 1 1 1 1 1
```

```r
## vector of mrca distances from root when lambda=0.5:
treeVec(tree,0.5)
```

```
##  [1] 0.0000 0.0000 0.7979 0.7979 0.7979 0.5205 0.0000 0.0000 0.0000 0.0000
## [11] 0.0000 0.0000 2.2933 1.6387 1.6387 0.6866 0.7541 0.7458 0.6312 0.6345
## [21] 0.9824
```

```r
## vector of mrca distances as a function of lambda:
vecAsFunction <- treeVec(tree,return.lambda.function=TRUE)
## evaluate the vector at lambda=0.5:
vecAsFunction(0.5)
```

```
##  [1] 0.0000 0.0000 0.7979 0.7979 0.7979 0.5205 0.0000 0.0000 0.0000 0.0000
## [11] 0.0000 0.0000 2.2933 1.6387 1.6387 0.6866 0.7541 0.7458 0.6312 0.6345
## [21] 0.9824
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
## [1] 4.123
```

```r
## branch-length focused (lambda=1) distance:
treeDist(tree_a,tree_b,1)
```

```
## [1] 3.794
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

Contributors:
* [Jacob Almagro-Garcia](http://www.well.ox.ac.uk/jacob-almagro-garcia)
* [Caroline Colijn](http://www.imperial.ac.uk/people/c.colijn)

Maintainer of the CRAN version:
* [Michelle Kendall](http://www.imperial.ac.uk/people/m.kendall)

