## library(shiny)
## library(treescape)

source("helpers.R")

## DEFINE THE SERVER SIDE OF THE APPLICATION
shinyServer(function(input, output, session) {
  ## LOAD PACKAGES
  if(!require("ape")) stop("ape is required")
  if(!require("ade4")) stop("ade4 is required")
  if(!require("adegraphics")) stop("ade4 is required")
  if(!require("treescape")) stop("treescape is required")
  if(!require("adegenet")) stop("adegenet is required")
  if(!require("phangorn")) stop("phangorn is required")
  if(!require("shinyRGL")) stop("shinyRGL is required")
  
  # suppress warning messages from creating temporary directories when 3d plotting
  suppressWarnings(warning("dir.create(dir)"))
  
  # the following resets the DensiTree plot every time the number of clusters changes - it was really slow without this
  rvs <- reactiveValues(showDensiTree=NULL)
  observeEvent(input$nclust, {
    rvs$showDensiTree <- NULL
  })
  observeEvent(input$selectedDensiTree, {
    rvs$showDensiTree <- 1  
  })
  
  ## GET DYNAMIC ANNOTATION
  ## not used:
  #graphTitle <- reactive({
  #  xax <- getXax()
  #  yax <- getYax()
  #  zax <- getZax()
  #  paste0(input$dataset, ": MDS scatterplot, axes ", xax, yax, zax)
  #})
  
  ## DEFINE CAPTION
  ## not used:
  #output$caption <- renderText({
  #  graphTitle()
  #})
  
  ######################################
  ### Define main reactive functions
  ######################################
  
  getDataType <- reactive({
    input$datatype
  })
  
  getDataSet <- reactive({
    input$dataset
  })

  ## GET DATA ##
  getData <- reactive({
    out <- NULL
    dataType <- getDataType()
    dataSet <- getDataSet()
    
    ## data is a distributed dataset
    if(dataType=="expl"){
      if(dataSet=="woodmiceTrees") data("woodmiceTrees", package="treescape", envir=environment())
      out <- get(dataSet)
    }
    
    ## data is an input file
    if(dataType=="file" && !is.null(input$datafile)){
      ## need to rename input file
      oldName <- input$datafile$datapath
      extension <- adegenet::.readExt(input$datafile$name)
      newName <- paste(input$datafile$datapath, extension, sep=".")
      file.rename(oldName, newName)
      
      if(tolower(extension) %in% c("rdata","rda")){
        out <- get(load(newName))
      }
      if(tolower(extension) %in% c("nex", "nexus")){
        if(!require(ape)) stop("ape is required to read in NEXUS (.nex, .nexus) files")
        out <- read.nexus(file=newName)
      }
      
      ## fix potential bug with input of two trees
      validate(
        need(length(out)>2, "treescape expects at least three trees. The function treeDist is suitable for comparing two trees.")
      )
      
      
      ## warn of potential problems with input of too many / large trees
      if(length(out)*length(out[[1]]$tip.label)>=500) {
        warning("You have supplied many and/or large trees, so Shiny is likely to be slow. If no output appears soon, consider supplying fewer trees or comparing them using the standard R functions in treescape.")
      }
      
      ## fix potential bug with tip labels - they need to match
      tipLabelProblem <- FALSE
      for (i in 1:length(out)) {
        if (!setequal(out[[i]]$tip.label,out[[1]]$tip.label)) {
          tipLabelProblem <- TRUE
          validate(
            need(!tipLabelProblem, "Trees must have identical tip labels for the current version of treescape")
          )
        } 
      }
      
    }
    
    validate(
      need(!is.null(out), "Waiting for data")
    )
    
    ## fix potential bug with names - they need to be defined and unique
    if(is.null(names(out))) {names(out) <- 1:length(out)}
    if(length(unique(names(out)))!=length(out)){
      warning("duplicates detected in tree labels - using generic names")
      names(out) <- 1:length(out)
    }
    
    ## return data
    return(out)
  }) # end getData
  
  ## GET number of trees
  getLengthData <- reactive({
    x <- getData()
    validate(
      need(!is.null(x), "Loading data set")
    )
    return(length(x))
  })
  
  ## GET tree names
  getTreeNames <- reactive({
    x <- getData()
    validate(
      need(!is.null(x), "Loading data set")
    )
    return(names(x))
  })
  
  ## GET tip labels
  getTipLabels <- reactive({
    x <- getData()
    validate(
      need(!is.null(x), "Loading data set")
    )
    return(x[[1]]$tip.label)
  })
  
  ## GET tree method
  getTreemethod <- reactive({
    input$treemethod
  }) # end getTreemethod
  
  ## GET number of axes retained
  getNaxes <- reactive({
    if(is.null(input$naxes)){
      naxes <- 3
    }
    else {
      naxes <- as.numeric(input$naxes)
      # when naxes changes we update the options available for the axes 
      # unfortunately I think they have to reset to their original 1,2,3 values
      # but at least they now only do this when naxes changes; they used to also do it for lambda etc.
      
      updateNumericInput(session,"xax", "Indicate the x axis", value=1, min=1, max=naxes)
      updateNumericInput(session,"yax", "Indicate the y axis", value=2, min=1, max=naxes)
      
      # (if relevant, update z axis selector too)
      dim <- getPlotDim()
      if (dim==3){
       updateNumericInput(session,"zax", "Indicate the z axis", value=3, min=1, max=naxes)
      }
        
    } 
  return(naxes)  
  }) # end getNaxes
  
  ## GET lambda
  getLambda <- reactive({
    l <- input$lambda
    ## the following removes the lambda error messages:
    validate(
      need(!is.null(l), "Loading data set")
    )	
   return(l)
  }) # end getLambda

  # GET the tree vectors as functions of lambda
  getKCtreeVecs <- reactive({
    x <- getData()
    validate(
      need(!is.null(x), "Loading data set")
    )
    df <-sapply(x, function(i) treeVec(i, return.lambda.function=TRUE)) 
  })
  

  # GET the tree vectors evaluated at lambda
  getKCtreeVecsAtLambda <- reactive({
  vectors <- getKCtreeVecs()
  l <- getLambda()
  validate(
    need(!is.null(vectors), "Analysing data")
  )
  t(sapply(vectors, function(i) i(l)))
  })
    

  ## GET KC matrix, evaluated at lambda
  getKCmatrix <- reactive({
    vls <- getKCtreeVecsAtLambda()
    numtrees <- getLengthData()
    as.dist(sapply(1:numtrees, function(a) sapply(1:numtrees, function(b) if(a<b) {sqrt(sum((vls[a,]-vls[b,])^2))} else{0})))
  }) # end getKCmatrix
  
  ## GET medTrees for all clusters
  getMedTreesList <- reactive({
    mat <- getKCtreeVecsAtLambda()
    groves <- getClusters()
    if(!is.null(groves$groups)){ # if clusters have been picked
    numGroups <- length(unique(groves$groups))
    med <- medTree(mat,groves$groups)
    lapply(1:numGroups, function(x) med[[x]]$treenumbers[[1]])
    }
    else{
      medTree(mat)$treenumbers[[1]]
    }
  })
  
  getMedTree <- reactive({
    x <- getData()
    whichClust <- input$selectedMedTree
    medList <- getMedTreesList()
    if(whichClust=="all"){
    x[[medList[[1]]]]
    }
    else{
    x[[medList[[as.numeric(whichClust)]]]]
    }
  })
  
  ## GET PCO analysis ##
  getPCO <- reactive({
    D <- getKCmatrix()
    naxes <- getNaxes()
    validate(
      need(!is.null(D), "Analysing data")
    )
    validate(
      need(!is.null(naxes), "Analysing data")
    )
    dudi.pco(D,scannf=FALSE,nf=naxes)
  }) # end getPCO

  ## GET ANALYSIS ##
  getAnalysis <- reactive({
    x <- getData()
    validate(
      need(!is.null(x), "Loading data set")
    )
    
    naxes <- getNaxes()
    TM <- getTreemethod()

    ## select method used to summarise tree
    if(!is.null(TM)){
      if(TM %in% c("patristic","nNodes","Abouheif","sumDD")){
        ## run treescape
        res <- treescape(x, method=TM, nf=naxes)
      } else if(TM=="metric"){
        ## don't actually need to call treescape here, to save on recomputation for varying lambda
        D <- getKCmatrix()
        pco <- getPCO()
        res <- list(D=D, pco=pco) 
      } 
    }
    
    ## return results
    return(res)
  }) # end getAnalysis

#################################################
### Little "get" functions to support getClusters
#################################################

getNclust <- reactive({
  if(!is.null(input$nclust)) {
    input$nclust
  } else {
    2
  }
}) 
  
getClustmethod <- reactive({
  input$clustmethod
})
  

################
## GET CLUSTERS
################

getClusters <- reactive({
    ## stop if clusters not required
    if(!input$findGroves) return(NULL)
  
    ## reset the densiTree plot to accommodate number of clusters available
    choices <- getClustChoices()
    updateSelectInput(session, "selectedDensiTree", "Choose collection of trees to view in densiTree plot", 
                    choices=choices, selected="")
    
    ## reset the median tree choices to accommodate number of clusters available
    updateSelectInput(session, "selectedMedTree", "Median tree from:", 
                      choices=choices, selected="all")
  
    ## get dataset
    x <- getData()
    validate(
      need(!is.null(x), "Loading data set")
    )
    
    naxes <- getNaxes()
    TM <- getTreemethod()
    nclust <- getNclust()
    clustmethod <- getClustmethod()
 
    ## select method used to summarise tree
    if(!is.null(TM)){
      if(TM %in% c("patristic","nNodes","Abouheif","sumDD")){
        treeMethod <- function(x){return(adephylo::distTips(x, method=TM))}
        ## run findGroves
        res <- findGroves(x, method=treeMethod, nf=naxes,
                          nclust=nclust, clustering=clustmethod)
      } else if(TM=="metric"){
        res <- findGroves(getAnalysis(), nclust=nclust, clustering=clustmethod)
      } else {
        treeMethod <- adephylo::distTips
        ## run findGroves
        res <- findGroves(x, method=treeMethod, nf=naxes,
                          nclust=nclust, clustering=clustmethod)
      }
    }
    
    ## return results
    return(res)
  
  }) # end getClusters
  
  
  
  ## DYNAMIC UI COMPONENTS ##
  ## SELECTION OF MDS AXES
  output$naxes <- renderUI({
    if(!is.null(getLengthData())) {
      nmax <- getLengthData()
    } else {
      nmax <- 100
    }
    sliderInput("naxes", "Number of MDS axes retained:", min=2, max=nmax, value=3, step=1)
    })
  
  ## VALUE OF LAMBDA FOR METRIC
  output$lambda <- renderUI({
    ## if KC metric has been chosen
    TM <- getTreemethod()
    if(TM=="metric") {
      sliderInput("lambda", "Value of lambda", min=0, max=1, value=0, step=0.01)
    } else {
      NULL
    }
  })
  
  ## SELECTION OF NUMBER OF CLUSTERS
  output$nclust <- renderUI({
    if(!is.null(x <- getData())) {
      nmax <- length(x)
    } else {
      nmax <- 100
    }
    nmax <- min(20, nmax)
    sliderInput("nclust", "Number of clusters:", min=2, max=nmax, value=2, step=1)
  })
  

######################################################
### Little "get" functions to support getPlot
######################################################  

getPalette  <- reactive({
  get(input$palette)
})

getLabcol <- reactive({
  ifelse(!is.null(input$labcol), input$labcol, "black")
})

getBgcol <- reactive({
  ifelse(!is.null(input$bgcol), input$bgcol, "white")
})

getScattertype <- reactive({
input$scattertype
})

getXax <- reactive({
  input$xax
})  

getYax <- reactive({
  input$yax
})  

getZax <- reactive({
  input$zax
})  

getScreemds <- reactive({
  input$screemds
})  

getOptimlabels <- reactive({
  input$optimlabels
})  

getShowlabels <- reactive({
  input$showlabels
})

getLabelsize <- reactive({
  input$labelsize
})

getPointsize <- reactive({
  input$pointsize
})

##############  
## GET plot
##############

## GET whether plot is 2D (default) or 3D
getPlotDim <- reactive({
  if(input$plot3D==TRUE) {3}
  else {2}
})

## GET 2D plot
getPlot <- reactive({

    res <- getAnalysis()
    groves <- getClusters()
    
    ## get aesthetics
    pal <- getPalette()
    labcol <- getLabcol()
    bgcol <- getBgcol()
    scattertype <- getScattertype()
    xax <- getXax()
    yax <- getYax()
    screemds <- getScreemds()
    optimlabels <- getOptimlabels()
    showlabels <- getShowlabels()
    labelsize <- getLabelsize()
    pointsize <- getPointsize()
      
    ## plot without groups
    if(is.null(groves)){
      plotGroves(res$pco, type=scattertype, xax=xax, yax=yax,
                        scree.posi=screemds, lab.optim=optimlabels,
                        lab.show=showlabels, lab.cex=labelsize,
                        lab.col=labcol,
                        point.cex=pointsize, bg=bgcol)
    } else {
    ## plot with groups
      plotGroves(groves, type=scattertype, xax=xax, yax=yax,
                        scree.posi=screemds, lab.optim=optimlabels,
                        lab.show=showlabels, lab.cex=labelsize,
                        lab.col=labcol,
                        point.cex=pointsize, bg=bgcol, col.pal=pal)
    }
  })


## TREESCAPE IMAGE ##
output$treescapePlot <- renderUI({
  dim <- getPlotDim()
  if(dim==2){
  plotOutput("scatterplot", height = "800px"
             #, click="plot_click" # get this working later
  )}
  else{
    validate(
    need(packageVersion("rgl")<'0.95.1247',
    paste0("You are running version ",packageVersion("rgl")," of the package rgl. The Shiny wrapper for rgl is not supported for versions >=0.95.1247. We recommend deleting the current version of rgl from your library then installing this patch: devtools::install_github('trestletech/rgl@js-class')")
    ))
    webGLOutput("plot3D", height = "800px")
  }
})


output$scatterplot <- renderPlot({
    myplot <- getPlot()
    plot(myplot)
}, res=120)

getPlot3d <- reactive({
  res <- getAnalysis()
  xax <- getXax()
  yax <- getYax()
  zax <- getZax()
  
  # show clusters?
  clusts <- getClusters()
  if (!is.null(clusts)){
    pal <- getPalette()
    cols3d <- fac2col(clusts$groups,col.pal=pal)
  }
  else{cols3d <- "navy"}
  
  rgl::plot3d(res$pco$li[,xax],res$pco$li[,yax],res$pco$li[,zax], 
              type="s", size=getPointsize(),
              xlab="",ylab="",zlab="",
              col=cols3d)
})

output$plot3D <- renderWebGL({
  plot <- getPlot3d() # separated these out to enable easier save function
}) 
  
# get tree and aesthetics for plotting tree  
getTreeChoice <- reactive({
  input$treeChoice
})


getTree <- reactive({
  x <- getData()
  validate(
    need(!is.null(x), "Loading data set")
  )
  treechoice <- getTreeChoice()
  if(treechoice=="med"){
    tre <- getMedTree()
  }
  else{
    g <- input$selectedGenTree
    if(is.null(g)){tre <- NULL}  
    else{
    treeNum <- as.numeric(g)
    tre <- x[[treeNum]]
    }
  }
  
  #trelab <- input$selectedTree
  #if(trelab!=""){
  #  ## numeric label
  #  if(!is.na(as.numeric(trelab))){
  #    validate(
  #      need(as.numeric(trelab) %in% 1:length(x), paste0("Tree number must be in the range 1 to ",length(x)))
  #    )	
  #    tre <- x[[as.numeric(trelab)]]
  #  } else {
  #    ## text label
  #    validate(
  #      need(as.character(trelab) %in% names(x), "Tree name not recognised")
  #    )	
  #    tre <- x[[as.character(trelab)]]
  #  }
    
    # return tree
    if(!is.null(tre)){
      if(input$ladderize){
        tre <- ladderize(tre)
      }
    return(tre)   
    }
    else{
      NULL
    }
})  

## PHYLOGENY ##
output$tree <- renderPlot({
  tre <- getTree()
  if(!is.null(tre)){
  
  ## plot tree ##
  par(mar=rep(2,4), xpd=TRUE)
  plot(tre, type=input$treetype,
         use.edge.length=as.logical(input$edgelengths),
         show.tip.lab=input$showtiplabels, 
         font=as.numeric(input$tiplabelfont), 
         cex=input$tiplabelsize,
         direction=input$treedirection,
         edge.width=input$edgewidth,
         edge.color=input$edgecolor
         )
  }
})
  
## DENSITREE

# The slider bar is always at least 2 even when clusters haven't
# been requested, so we can't just use getNclust.



getNclustForDensiTree <- reactive({
  if(input$findGroves==FALSE){NULL}
  else{input$nclust}
}) 

getClustChoices <- reactive({
  nclust <- getNclustForDensiTree()
  if(is.null(nclust)){
    choices <- c("","all")
    names(choices) <- c("Choose one","All trees")
  }
  else{
    choices <- c("",1:nclust,"all")
    names(choices) <- c("Choose one",paste0("Cluster ",1:nclust),"All trees")
  }
  return(choices)
})

getDensiTree <- reactive({
  clusterNo <- input$selectedDensiTree
  if(clusterNo==""){
    NULL
  }
  else if(clusterNo=="all"){
    x <- getData()
    return(x)
  }
  else{
    x <- getData()
    clusts <- getClusters()
    clustTrees <- x[which(clusts$groups==as.numeric(clusterNo))]
    return(clustTrees)
  }
})  

output$densiTree <- renderPlot({
  if(is.null(rvs$showDensiTree)) {NULL}
  else{
  clustTrees <- getDensiTree()
  validate(
    need(!is.null(clustTrees), "Loading densiTree plot")
  )	
  densiTree(clustTrees, col=4, alpha=input$alpha, scaleX=input$scaleX)
  }
})


  ## EXPORT TREES ##
  output$exporttrees <- downloadHandler(
    filename = function() { paste(input$dataset, '.nex', sep='') },
    content = function(file) {
      if(!require(ape)) stop("ape is required to save trees into nexus file")
      x <- getData()
      if(!is.null(x) && inherits(x, "multiPhylo")) ape::write.nexus(x, file=file)
    })
  
  ## EXPORT ANALYSIS TO CSV ##
  output$exportrestocsv <- downloadHandler(
    filename = function() { paste(input$dataset, "-analysis", '.csv', sep='') },
    content = function(file) {
      x <- getData()
      res <- getClusters()
      if(!is.null(res)){
        tab <- cbind.data.frame(res$groups, res$treescape$pco$li)
        names(tab) <- c("cluster", paste("PC", 1:ncol(res$treescape$pco$li), sep="."))
        row.names(tab) <- names(x)
      } else{
        res <- getAnalysis()
        tab <- res$pco$li
        names(tab) <- paste("PC", 1:ncol(tab), sep=".")
        row.names(tab) <- names(x)
      }
      if(!is.null(res)) write.csv(tab, file=file)
    })
  
  
  ## EXPORT ANALYSIS TO RDATA ##
  output$exportrestordata <- downloadHandler(
    filename = function() { paste(input$dataset, "-analysis", '.RData', sep='') },
    content = function(file) {
      trees <- getData()
      analysis <- getClusters()
      if(is.null(analysis)) analysis <- getAnalysis()
      if(!is.null(analysis)) {
        save(trees, analysis, file=file)
      }
    })
  
  ## EXPORT MDS PLOT AS PNG ##
  output$downloadMDS <- downloadHandler(
    filename = function() { paste0(input$dataset,".png") },
    content = function(file) {
      dim <- getPlotDim()
      if (dim==2){
        myplot <- getPlot()
        png(file=file, width = 10, height = 10, units = 'in', res = 500)
        plot(myplot)
        dev.off()
      }
    else{
        validate(
          need(dim==2,"Sorry, saving is not yet enabled for 3D plots")
        )
      # Note rgl.postscript may be the way forward. I think it will require a separate download button.
        }
    contentType = 'image/png'  }
  )
  
  
  ## EXPORT TREE PLOT AS PNG ##
  output$downloadTree <- downloadHandler(
    filename = function() { paste(input$dataset, 'Tree',getTree(),'.png', sep='') },
    content = function(file) {
      tre <- getTree()
      png(file=file)
      plot(tre, type=input$treetype,
             show.tip.lab=input$showtiplabels, font=1, cex=input$tiplabelsize,
             direction=input$treedirection,
             edge.width=input$edgewidth)
      dev.off()
      contentType = 'image/png'
    }
    )
  
  ## EXPORT DENSITREE PLOT AS PNG ##
  output$downloadDensiTree <- downloadHandler(
    filename = function() { paste(input$dataset, 'DensiTree',input$selectedDensiTree,'.png', sep='') },
    content = function(file) {
      tre <- getDensiTree()
      png(file=file) 
      plot(tre, type=input$treetype, # CHANGE THIS!
           show.tip.lab=input$showtiplabels, font=1, cex=input$tiplabelsize,
           direction=input$treedirection,
           edge.width=input$edgewidth)
      dev.off()
      contentType = 'image/png'
    }
  )
  
output$selectedGenTree <- renderUI({
  numTrees <- getLengthData()
  treeNames <- getTreeNames()
  choices <- c("",1:numTrees)
  names(choices) <- c("Choose one",treeNames)
  selectInput("selectedGenTree", "Choose individual tree", 
              choices=choices, selected="")
  })
  
#  # Clicking
#  output$plot_click <- renderPrint({
#    res <- getAnalysis()
#    xax <- getXax()
#    yax <- getYax()
#    ## need to get the scaling right:
#    ## the input$plot_click points are in [0,1]x[0,1] (roughly but not quite!)
#    ## also need to know how scatter places the points...
#    minx <- min(res$pco$li[,xax])
#    miny <- min(res$pco$li[,yax])
#    res$pco$li[,xax] <- res$pco$li[,xax] - minx # since minx is negative
#    res$pco$li[,yax] <- res$pco$li[,yax] - miny
#    maxx <- max(res$pco$li[,xax])
#    maxy <- max(res$pco$li[,yax])
#    longestaxis <- max(maxx,maxy)
#    scaledCoords <- res$pco$li/longestaxis
#    if(!is.null(input$plot_click)){
#      test <- nearPoints(scaledCoords,input$plot_click,
#                         xvar=paste0("A",xax),yvar=paste0("A",yax),
#                         threshold=100,maxpoints=1)
#      nearTree <- intersect(which(res$pco$li[,xax]==test[[1]]),
#                            which(res$pco$li[,yax]==test[[2]]))
#    paste0(#"Click coordinates: x=",input$plot_click$x,", y=",input$plot_click$y,", giving tree ",
#      nearTree)
#           #intersect(which(res$pco$li[,as.numeric(xax)]==signif(input$plot_click$x,2)),
#          #           which(res$pco$li[,as.numeric(yax)]==signif(input$plot_click$y,2))))
#    }
#      })
  
  
  
  ## RENDER SYSTEM INFO ##
  output$systeminfo <- .render.server.info()
  
}) # end shinyServer