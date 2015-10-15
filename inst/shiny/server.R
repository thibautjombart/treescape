## library(shiny)
## library(treescape)


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
    paste(input$dataset, ": MDS scatterplot, axes ", input$xax,"-", input$yax, sep="")
  })
  
  ## DEFINE CAPTION
  output$caption <- renderText({
    graphTitle()
  })
  
  
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
  
  
  ## GET ANALYSIS ##
  getAnalysis <- reactive({
    
    ## get dataset
    x <- getData()
    
    ## the following removes the lambda error messages:
    validate(
      need(input$lambda != "", "Loading data set")
    )	
    
    ## Trying to remove the error message "only defined on a data frame with all numeric variables"
    ## This doesn't work, but presumably it's something like this?!
    #validate(
    #  need(input$treeMethod != "", "My test")
    #  )		
    
    ## stop if not data
    if(is.null(x)) return(NULL)
    
    ## get number of axes retained
    if(!is.null(input$naxes)) {
      naxes <- input$naxes
    } else {
      naxes <- 2
    }
    
    ## select method used to summarise tree
    if(!is.null(input$treemethod)){
      if(input$treemethod %in% c("patristic","nNodes","Abouheif","sumDD")){
        treeMethod <- function(x){return(adephylo::distTips(x, method=input$treemethod))}
      } else if(input$treemethod=="metric"){
        treeMethod <- function(x){return(treeVec(x, lambda=input$lambda))}
      } else {
        treeMethod <- adephylo::distTips
      }
    }
    
    ## run treescape
    res <- treescape(x, method=treeMethod, nf=naxes)
    
    ## return results
    return(res)
  }) # end getAnalysis
  
  
  ## GET CLUSTERS ##
  getClusters <- reactive({
    ## stop if clusters not required
    if(!input$findgroups) return(NULL)
    
    ## get dataset
    x <- getData()
    
    ## stop if not data
    if(is.null(x)) return(NULL)
    
    ## get number of axes retained
    if(!is.null(input$naxes)) {
      naxes <- input$naxes
    } else {
      naxes <- 2
    }
    
    ## get number of clusters
    if(!is.null(input$nclust)) {
      nclust <- input$nclust
    } else {
      nclust <- 2
    }
    
    
    ## select method used to summarise tree
    if(!is.null(input$treemethod)){
      if(input$treemethod %in% c("patristic","nNodes","Abouheif","sumDD")){
        treeMethod <- function(x){return(adephylo::distTips(x, method=input$treemethod))}
      } else if(input$treemethod=="metric"){
        treeMethod <- function(x){return(treeVec(x, lambda=input$lambda))}
      } else {
        treeMethod <- adephylo::distTips
      }
    }
    
    ## run findGroves
    res <- findGroves(x, method=treeMethod, nf=naxes,
                      nclust=nclust, clustering=input$clustmethod)
    
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
      sliderInput("lambda", "Value of lambda", min=0, max=1, value=0.5, step=0.01)
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
  
  ## ## SELECTION OF BACKGROUND COLOR
  ## output$bg <- renderUI({
  ##     selectInput("bg", "Background color:", colors(), selected="white")
  ## })
  
  ##  ## SELECTION OF LABELS COLOR
  ## output$labcol <- renderUI({
  ##     selectInput("labcol", "Label color:", colors(), selected="black")
  ## })
  
  
  ## ANALYSIS ##
  output$scatterplot <- renderPlot({
    ## get dataset
    x <- getData()
    
    if(!is.null(x)){
      ## get analysis
      res <- getAnalysis()
      
      ## get clusters
      groves <- getClusters()
      
      ## get palette
      pal <- get(input$palette)
      
      ## get colors
      labcol <- ifelse(!is.null(input$labcol), input$labcol, "black")
      bgcol <- ifelse(!is.null(input$bgcol), input$bgcol, "white")
      
      ## plot without groups
      if(is.null(groves)){
        y <- plotGroves(res$pco, type=input$scattertype, xax=input$xax, yax=input$yax,
                   scree.posi=input$screemds, lab.optim=input$optimlabels,
                   lab.show=input$showlabels, lab.cex=input$labelsize,
                   lab.col=labcol,
                   point.cex=input$pointsize, bg=bgcol)
      } else {
        ## plot with groups
        y <- plotGroves(groves, type=input$scattertype, xax=input$xax, yax=input$yax,
                   scree.posi=input$screemds, lab.optim=input$optimlabels,
                   lab.show=input$showlabels, lab.cex=input$labelsize,
                   lab.col=labcol,
                   point.cex=input$pointsize, bg=bgcol, col.pal=pal)
      }
    plot(y)
    }
  }, res=120)
  
  ## PHYLOGENY ##
  output$tree <- renderPlot({
    ## get dataset
    x <- getData()
    
    ## get right tree ##
    trelab <- input$selectedTree
    if(trelab!=""){
      ## numeric label
      if(!is.na(as.numeric(trelab))){
        tre <- x[[as.numeric(trelab)]]
      } else {
        ## text label
        tre <- x[[as.numeric(trelab)]]
      }
      
      if(input$ladderize){
        tre <- ladderize(tre)
      }
      
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
  # Note this seems inefficient - is there a more direct way to call the plot?
  output$downloadMDS <- downloadHandler(
    filename = function() { paste(input$dataset, '.png', sep='') },
    content = function(file) {
      x <- getData()
      
      if(!is.null(x)){
        ## get analysis
        res <- getAnalysis()
        
        ## get clusters
        groves <- getClusters()
        
        ## get palette
        pal <- get(input$palette)
        
        ## get colors
        labcol <- ifelse(!is.null(input$labcol), input$labcol, "black")
        bgcol <- ifelse(!is.null(input$bgcol), input$bgcol, "white")
      
        ## plot without groups
        if(is.null(groves)){
          y <- plotGroves(res$pco, type=input$scattertype, xax=input$xax, yax=input$yax,
                          scree.posi=input$screemds, lab.optim=input$optimlabels,
                          lab.show=input$showlabels, lab.cex=input$labelsize,
                          lab.col=labcol,
                          point.cex=input$pointsize, bg=bgcol)
        } else {
          ## plot with groups
          y <- plotGroves(groves, type=input$scattertype, xax=input$xax, yax=input$yax,
                          scree.posi=input$screemds, lab.optim=input$optimlabels,
                          lab.show=input$showlabels, lab.cex=input$labelsize,
                          lab.col=labcol,
                          point.cex=input$pointsize, bg=bgcol, col.pal=pal)
        }
      }
      png(file=file, width = 10, height = 10, units = 'in', res = 500)
      plot(y)
      dev.off()
    contentType = 'image/png'
    })
  
  
  ## EXPORT TREE PLOT AS PNG ##
  output$downloadTree <- downloadHandler(
    filename = function() { paste(input$dataset, 'Tree',input$selectedTree,'.png', sep='') },
    content = function(file) {
      png(file=file)
      x <- getData()
      
      ## get right tree ##
      trelab <- input$selectedTree
      if(trelab!=""){
        ## numeric label
        if(!is.na(as.numeric(trelab))){
          tre <- x[[as.numeric(trelab)]]
        } else {
          ## text label
          tre <- x[[as.numeric(trelab)]]
        }
        
        if(input$ladderize){
          tre <- ladderize(tre)
        }
        
        ## plot tree ##
        par(mar=rep(2,4), xpd=TRUE)
        plot(tre, type=input$treetype,
             show.tip.lab=input$showtiplabels, font=1, cex=input$tiplabelsize,
             direction=input$treedirection,
             edge.width=input$edgewidth)
      }
      dev.off()
      contentType = 'image/png'
    })
  
  ## RENDER SYSTEM INFO ##
  output$systeminfo <- .render.server.info()
  
}) # end shinyServer