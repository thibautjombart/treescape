library(shiny)
library(exploratree)


## DEFINE THE SERVER SIDE OF THE APPLICATION
shinyServer(function(input, output) {

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

        if(!is.null(input$datafile)){
            ## need to rename input file
            oldName <- input$datafile$datapath
            extension <- .readExt(input$datafile$name)
            newName <- paste(input$datafile$datapath, extension, sep=".")
            file.rename(oldName, newName)

            if(extension %in% c("RData","Rdata","Rda","rda")){
                out <- get(load(newName))
            }
        }
        return(out)
    })




    ## DYNAMIC UI COMPONENTS ##
    ## SELECTION OF MDS AXES
    output$naxes <- renderUI({
        if(!is.null(x <- getData())) {
            nmax <- length(x@tab)
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

    output$xax <- renderUI({
        if(!is.null(x <- getData())) {
            nmax <- length(x)
        } else {
            nmax <- 100
        }
        numericInput("yax", "Indicate the y axis", value=2, min=2, max=nmax)
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
            res <- exploratree(x, nf=naxes)
            scatter(res$pco, xax=input$xax, yax=input$yax,
                    cex=input$pointsize, clabel=input$labelsize)
        } else {
            NULL
        }
    })

    ## RENDER SYSTEM INFO ##
    output$systeminfo <- .render.server.info()

}) # end shinyServer
