library(shiny)
## library(exploratree)

## DEFINE UI ##
shinyUI(
    pageWithSidebar(
        ##  TITLE
        headerPanel("exploratree"),

        ## SIDE PANEL CONTENT
        sidebarPanel(

            ## choice of dataset if source is a file
            conditionalPanel(condition = "$('li.active a').first().html()!= 'Help'",
                             fileInput('datafile', 'Choose input file',
                                       accept=c('RData/Rdata/Rda/rda', 'R data (*.RData,*.Rda)')),
                             tags$hr()
                             ),

            conditionalPanel(
                ## condition
                "$('li.active a').first().html()!= 'Help'",
                selectInput("treemethod", "Choose a tree summary:",
                            choices=c("Patristic distances" = "patristic",
                            "Number of nodes" = "nNodes",
                            "Abouheif's metric" = "Abouheif",
                            "Sum of direct descendents" = "sumDD"))
                ),

            conditionalPanel(
                ## condition
                "$('li.active a').first().html()!= 'Help'",
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

                h3("Aesthetics"),

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
                            "Top left" = "topleft"))
                )

            ), # end sidebarPanel

        ## MAIN PANEL
        mainPanel(
            tabsetPanel(

                tabPanel("Scatterplot",plotOutput("scatterplot")),


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
