library(shiny)
##library(exploratree)


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
            extension <- adegenet::.readExt(input$datafile$name)
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
                    treeMethod <- function(x){return(distTips(x, method=input$treemethod))}
                } else{
                    treeMethod <- distTips
                }
            } else {
                treeMethod <- distTips
            }

            ## run exploratree
            res <- exploratree(x, method=treeMethod, nf=naxes)

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

    ## RENDER SYSTEM INFO ##
    output$systeminfo <- .render.server.info()

}) # end shinyServer
