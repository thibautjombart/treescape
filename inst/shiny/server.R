## library(shiny)
## library(treescape)


## DEFINE THE SERVER SIDE OF THE APPLICATION
shinyServer(function(input, output) {
    ## LOAD PACKAGES
    if(!require("ape")) stop("ape is required")
    if(!require("ade4")) stop("ade4 is required")
    if(!require("treescape")) stop("treescape is required")

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
        }

        ## return data
        return(out)
    })




    ## DYNAMIC UI COMPONENTS ##
    ## SELECTION OF MDS AXES
    output$naxes <- renderUI({
        if(!is.null(x <- getData())) {
            nmax <- length(x)
        } else {
            nmax <- 100
        }
        sliderInput("naxes", "Number of MDS axes retained:", min=1, max=nmax, value=10, step=1)
    })

    ## SELECTION OF PLOTTED AXES
    output$xax <- renderUI({
        if(!is.null(x <- getData())){
            nmax <- length(x)
        } else {
            nmax <- 100
        }
        numericInput("xax", "Indicate the x axis", value=1, min=1, max=nmax)
    })

    output$yax <- renderUI({
        if(!is.null(x <- getData())) {
            nmax <- length(x)
        } else {
            nmax <- 100
        }
        numericInput("yax", "Indicate the y axis", value=2, min=2, max=nmax)
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


    ## SCATTERPLOT ##
    output$scatterplot <- renderPlot({
        ## get dataset
        x <- getData()

        ## get number of axes retained
        if(!is.null(input$naxes)) {
            naxes <- input$naxes
        } else {
            naxes <- 2
        }

        if(!is.null(x)){
            ## select method used to summarise tree
            if(!is.null(input$treemethod)){
                if(input$treemethod %in% c("patristic","nNodes","Abouheif","sumDD")){
                    treeMethod <- function(x){return(adephylo::distTips(x, method=input$treemethod))}
                } else if(input$treemethod=="metric"){
                    treeMethod <- function(x){return(tree.vec(x, lambda=input$lambda))}
                } else {
                    treeMethod <- adephylo::distTips
                }
            }

            ## run treescape
            res <- treescape(x, method=treeMethod, nf=naxes)

            ## make scatterplot
            s.label(res$pco$li, xax=input$xax, yax=input$yax,
                    cpoint=input$pointsize, clabel=input$labelsize)

            ## add legend
            if(input$screemds!="none"){
                add.scatter.eig(res$pco$eig, res$pco$nf, xax=input$xax, yax=input$yax, posi=input$screemds)
            }
        } else {
            NULL
        }
    })

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
        }

        ## plot tree ##
        plot(tre)

    })

    ## RENDER SYSTEM INFO ##
    output$systeminfo <- .render.server.info()

}) # end shinyServer
