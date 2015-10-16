## library(shiny)
## library(treescape)

source("helpers.R")

## DEFINE THE SERVER SIDE OF THE APPLICATION
shinyServer(function(input, output) {
  ## LOAD PACKAGES
  if(!require("ape")) stop("ape is required")
  if(!require("ade4")) stop("ade4 is required")
  if(!require("adegraphics")) stop("ade4 is required")
  if(!require("treescape")) stop("treescape is required")
  if(!require("adegenet")) stop("adegenet is required")
  
  
  ## GET DYNAMIC ANNOTATION
  graphTitle <- reactive({
    xax <- getXax()
    yax <- getYax()
    paste(input$dataset, ": MDS scatterplot, axes ", xax,"-", yax, sep="")
  })
  
  ## DEFINE CAPTION
  output$caption <- renderText({
    graphTitle()
  })

  ######################################
  ### Define main reactive functions
  ######################################
  
  ## GET DATA ##
  getData <- reactive({
    out <- NULL
    
    ## data is a distributed dataset
    if(input$datatype=="expl"){
      if(input$dataset=="woodmiceTrees") data("woodmiceTrees", package="treescape", envir=environment())
      out <- get(input$dataset)
    }
    
    ## data is an input file
    if(input$datatype=="file" && !is.null(input$datafile)){
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
      
      ## fix potential bug with names - they need to be unique
      if(length(unique(names(out)))!=length(out)){
        warning("duplicates detected in tree labels - using generic names")
        names(out) <- 1:length(out)
      }
      
      ## fix potential bug with input of two trees
      if(length(out)<3) {
        stop("treescape expects at least three trees. The function treeDist is suitable for comparing two trees.")
      }
      
      
      ## warn of potential problems with input of too many / large trees
      if(length(out)*length(out[[1]]$tip.label)>=500) {
        warning("You have supplied many and/or large trees, so Shiny is likely to be slow. If no output appears soon, consider supplying fewer trees or comparing them using the standard R functions in treescape.")
      }
      
      ## fix potential bug with tip labels - they need to match
      for (i in 1:length(out)) {
        if (!setequal(out[[i]]$tip.label,out[[1]]$tip.label)) {
          stop("trees have different tip labels")
        } 
      }
      
    }
    
    ## return data
    return(out)
  }) # end getData
  
  ## GET tree method
  getTreemethod <- reactive({
    input$treemethod
  }) # end getTreemethod
  
  ## GET number of axes retained
  getNaxes <- reactive({
    if(!is.null(input$naxes)) {
      naxes <- input$naxes
    } else {
      naxes <- 2
    } 
  }) # end getNaxes
  
  ## GET lambda
  getLambda <- reactive({
    ## the following removes the lambda error messages:
    validate(
      need(input$lambda != "", "Loading data set")
    )	
   input$lambda
  }) # end getLambda
  
  # compute the tree vectors (as functions of lambda) only once per dataset to save on recomputation
  getKCmatrixfunction <- reactive({
   x <- getData()
    validate(
      need(!is.null(x), "Loading data set")
    )
   multiDist(x, return.lambda.function = TRUE)
  })
  
  
  ## GET KC matrix, evaluated at lambda
  getKCmatrix <- reactive({
    l <- getLambda()
    M <- getKCmatrixfunction()
    return(M(l))
  }) # end getKCmatrix
  
  ## GET PCO analysis ##
  getPCO <- reactive({
    D <- getKCmatrix()
    naxes <- getNaxes()
    dudi.pco(D,scannf=is.null(naxes),nf=naxes)
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
        treeMethod <- function(x){return(adephylo::distTips(x, method=TM))}
        ## run treescape
        res <- treescape(x, method=treeMethod, nf=naxes)
      } else if(TM=="metric"){
        ## don't actually need to call treescape here, to save on recomputation for varying lambda
        D <- getKCmatrix()
        pco <- getPCO()
        res <- list(D=D, pco=pco) 
      } else {
        treeMethod <- adephylo::distTips
        ## run treescape
        res <- treescape(x, method=treeMethod, nf=naxes)
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
    if(!is.null(x <- getData())) {
      nmax <- length(x)
    } else {
      nmax <- 100
    }
    sliderInput("naxes", "Number of MDS axes retained:", min=1, max=nmax, value=2, step=1)
  })
  
  ## SELECTION OF PLOTTED AXES
  output$xax <- renderUI({
    if(!is.null(x <- getAnalysis())){
      nmax <- x$pco$nf
    } else {
      nmax <- 100
    }
    numericInput("xax", "Indicate the x axis", value=1, min=1, max=nmax)
  })
  
  output$yax <- renderUI({
    if(!is.null(x <- getAnalysis())) {
      nmax <- x$pco$nf
    } else {
      nmax <- 100
    }
    numericInput("yax", "Indicate the y axis", value=2, min=1, max=nmax)
  })
  
  ## VALUE OF LAMBDA FOR METRIC
  output$lambda <- renderUI({
    ## if metric has been chosen
    if(input$treemethod=="metric") {
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

getPlot <- reactive({
    ## get dataset
    x <- getData()
    validate(
      need(!is.null(x), "Loading data set")
    )
    
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
output$scatterplot <- renderPlot({
    myplot <- getPlot()
    plot(myplot)
}, res=120)
  
# get tree and aesthetics for plotting tree  

getTree <- reactive({
  x <- getData()
  validate(
    need(!is.null(x), "Loading data set")
  )
  trelab <- input$selectedTree
  if(trelab!=""){
    ## numeric label
    if(!is.na(as.numeric(trelab))){
      validate(
        need(as.numeric(trelab) %in% 1:length(x), paste0("Tree number must be in the range 1 to ",length(x)))
      )	
      tre <- x[[as.numeric(trelab)]]
    } else {
      ## text label
      validate(
        need(as.character(trelab) %in% names(x), "Tree name not recognised")
      )	
      tre <- x[[as.character(trelab)]]
    }
    
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
  }
})  

## PHYLOGENY ##
output$tree <- renderPlot({
  tre <- getTree()
  if(!is.null(tre)){
  
  ## plot tree ##
  par(mar=rep(2,4), xpd=TRUE)
  plot(tre, type=input$treetype,
         show.tip.lab=input$showtiplabels, font=1, cex=input$tiplabelsize,
         direction=input$treedirection,
         edge.width=input$edgewidth)
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
      myplot <- getPlot()
      png(file=file, width = 10, height = 10, units = 'in', res = 500)
      plot(myplot)
      dev.off()
    contentType = 'image/png'
    })
  
  
  ## EXPORT TREE PLOT AS PNG ##
  output$downloadTree <- downloadHandler(
    filename = function() { paste(input$dataset, 'Tree',input$selectedTree,'.png', sep='') },
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