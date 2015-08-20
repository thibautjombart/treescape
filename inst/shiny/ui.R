## CHECKS ##
if(!require(shiny)) stop("shiny is required")

## DEFINE UI ##
shinyUI(
    pageWithSidebar(
        ##  TITLE ##
        headerPanel(
            img(src="img/logo.png", height="160")
            ),

        ## SIDE PANEL CONTENT ##
        sidebarPanel(
            tags$head(tags$style(
                type = 'text/css',
                'form.well { max-height: 800px; overflow-y: auto; }'
                )),

            ## INPUT
            ## choice of type of data source
            conditionalPanel(condition = "$('li.active a').first().html()!= 'Help'",
                             img(src="img/line.png", width="100%"),
                             h2(HTML('<font color="#6C6CC4" size="6"> > Input / output </font>')),
                             radioButtons("datatype", HTML('<font size="4"> Choose data source:</font>'),
                                          list("Example from treescape"="expl","Input file"="file"))),

            ## choice of dataset if source is an example
            conditionalPanel(condition = "input.datatype=='expl'&& $('li.active a').first().html()!= 'Help'",
                             selectInput("dataset",
                                         HTML('<font size="4"> Select an example dataset:</font>'),
                                         choices=c("woodmiceTrees"))
                             ),

            ## choice of dataset if source is a file
            conditionalPanel(condition = "input.datatype=='file'&& $('li.active a').first().html()!= 'Help'",
                             fileInput("datafile", p(HTML(' <font size="4"> Choose input file:</font>'), br(),
                                                     strong("accepted formats:"), br(),
                                                     em("- multiphylo"), "saved from R (.RData/.rda)", br(),
                                                     em("- nexus"), "file (.nex/.nexus)")
                                       )
                             ),

            ## output data
            conditionalPanel(condition = " $('li.active a').first().html()!= 'Help'",
                             downloadButton('exporttrees', "Export trees (to nexus format)")
                             ),


            ## ANALYSIS
            conditionalPanel(
                ## condition
                "$('li.active a').first().html()!= 'Help'",
                img(src="img/line.png", width="100%"),
                h2(HTML('<font color="#6C6CC4" size="6"> > Analysis </font>')),
                selectInput("treemethod", "Choose a tree summary:",
                            choices=c(
                            "Metric" = "metric",
                            "Patristic distances" = "patristic",
                            "Number of nodes" = "nNodes",
                            "Abouheif's metric" = "Abouheif",
                            "Sum of direct descendents" = "sumDD"))
                ),

            conditionalPanel(
                ## condition
                "$('li.active a').first().html()!= 'Help'",
                uiOutput("lambda"),
                uiOutput("naxes")
                ),

            ## ## inputs specific of explorer tab
            conditionalPanel(
                ## condition
                "$('li.active a').first().html()=='Tree landscape explorer'",

                ## select first axis to plot
                ##numericInput("xax", "Indicate the x axis", value=1, min=1),
                uiOutput("xax"),

                ## select second axis to plot
                ##numericInput("yax", "Indicate the y axis", value=1, min=1),
                uiOutput("yax")
                ),

            conditionalPanel(
                ## condition
                "$('li.active a').first().html()==='Tree landscape explorer'",

                img(src="img/line.png", width="100%"),

                h2(HTML('<font color="#6C6CC4" size="6"> > Aesthetics </font>')),

                ## SCATTERPLOT AESTHETICS
                ## display labels
                checkboxInput("showlabels", label="Display labels?", value=TRUE),
                ## symbol size
                sliderInput("pointsize", "Size of the points", value=1, min=0, max=10, step=0.2),

                ## label size
                sliderInput("labelsize", "Size of the labels", value=1, min=0, max=10, step=0.2),

                ## add screeplot of MDS?
                selectInput("screemds", "Position of the MDS screeplot:",
                            choices=c("None" = "none",
                            "Bottom right" = "bottomright",
                            "Bottom left" = "bottomleft",
                            "Top right" = "topright",
                            "Top left" = "topleft"),
                            selected="bottomleft"),

                ## TREE AESTHETICS
                ## condition on tree being displayed
                conditionalPanel(condition = "input.selectedTree!=''",
                                 ## type of tree
                                 radioButtons("treetype", "Type of tree",
                                              choices=c("phylogram","cladogram", "fan", "unrooted", "radial"),
                                              selected="phylogram", width="100%"),

                                 ## tree direction
                                 radioButtons("treedirection", "Direction of the tree",
                                              choices=c("rightwards", "leftwards", "upwards", "downwards"),
                                              selected="rightwards", width="100%"),

                                 ## ladderize
                                 checkboxInput("ladderize", label="Ladderize the tree?", value=TRUE),

                                 ## tip labels
                                 checkboxInput("showtiplabels", label="Display tip labels?", value=TRUE),

                                 ## tip label size
                                 sliderInput("tiplabelsize", "Size of the tip labels", value=1, min=0, max=5, step=0.1),

                                 ## edge width
                                 sliderInput("edgewidth", "Width of the edges", value=1, min=0, max=20, step=0.2)

                                 )
                ),

            br(),br(),br(),br(),br(),br(),br(), # add some blank space at the end of side panel
            br(),br(),br(),br(),br(),br(),br(), # add some blank space at the end of side panel
            width=3), # end sidebarPanel; width is out of 12

        ## MAIN PANEL
        mainPanel(
            tabsetPanel(

                tabPanel("Tree landscape explorer",
                         plotOutput("scatterplot"),

                         ## add tree selector
                         textInput("selectedTree", "Choose tree (number or label)", value = ""),

                         ## conditional panel: plot tree if needed
                         conditionalPanel(
                             condition = "input.selectedTree!=''",
                             plotOutput("tree")
                             )
                         ),


                ## HELP SECTION
                tabPanel("Help",
                         HTML(paste(readLines("www/html/help.html"), collapse=" "))
                         ),

                ## SERVER INFO ##
                tabPanel("System info", verbatimTextOutput("systeminfo"))
                ) # end tabsetPanel
            ) # end mainPanel
        ) # end pageWithSidebar
    ) # end shinyUI
