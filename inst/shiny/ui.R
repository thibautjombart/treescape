library(shiny)
## library(treescape)

## DEFINE UI ##
shinyUI(
    pageWithSidebar(
        ##  TITLE
        headerPanel(
            img(src="img/logo.png", height="160")
            ),

        ## SIDE PANEL CONTENT
        sidebarPanel(
            tags$head(tags$style(
                type = 'text/css',
                'form.well { max-height: 800px; overflow-y: auto; }'
                )),

            ## choice of type of data source
            conditionalPanel(condition = "$('li.active a').first().html()!= 'Help'",
                             img(src="img/line.png", width="100%"),
                             h2(HTML('<font color="#6C6CC4" size="6"> > Input </font>')),
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
                                                     em("- multiphylo"), "saved from R (.RData/.Rda)", br(),
                                                     em("- nexus"), "file (.nex/.nexus)")
                                       )
                             ),


            conditionalPanel(
                ## condition
                "$('li.active a').first().html()!= 'Help'",
                img(src="img/line.png", width="100%"),
                h2(HTML('<font color="#6C6CC4" size="6"> > Options </font>')),
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

            ## ## inputs specific of scatterplot tab
            conditionalPanel(
                ## condition
                "$('li.active a').first().html()=='Scatterplot'",

                ## select first axis to plot
                ##numericInput("xax", "Indicate the x axis", value=1, min=1),
                uiOutput("xax"),

                ## select second axis to plot
                ##numericInput("yax", "Indicate the y axis", value=1, min=1),
                uiOutput("yax")
                ),

            conditionalPanel(
                ## condition
                "$('li.active a').first().html()==='Scatterplot'",

                img(src="img/line.png", width="100%"),

                h2(HTML('<font color="#6C6CC4" size="6"> > Aesthetics </font>')),

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
                            selected="bottomleft")
                ),

            br(),br(),br(),br(),br(),br(),br(), # add some blank space at the end of side panel
            width=3), # end sidebarPanel; width is out of 12

        ## MAIN PANEL
        mainPanel(
            tabsetPanel(

                tabPanel("Scatterplot",
                         plotOutput("scatterplot"),

                         ## add tree selector
                         textInput("selectedTree", "Choose tree (label)", value = "1")
                         ),


                ## HELP SECTION
                tabPanel("Help",
                         h3("Functions under development - do not use!")
                         ),

                ## SERVER INFO ##
                tabPanel("System info", verbatimTextOutput("systeminfo"))
                ) # end tabsetPanel
            ) # end mainPanel
        ) # end pageWithSidebar
    ) # end shinyUI
