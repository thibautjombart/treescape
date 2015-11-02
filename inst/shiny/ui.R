## CHECKS ##
if(!require("shiny")) stop("shiny is required")
if(!require("RLumShiny")) stop("RLumShiny is required")
if(!require("shinyBS")) stop("shinyBS is required")
if(!require("shinyRGL")) stop("shinyRGL is required")

## DEFINE UI ##
shinyUI(
  navbarPage("",position="fixed-top", collapsible = TRUE,
             theme = "bootstrap.simplex.css", 
             
      tabPanel("Tree landscape explorer",
             tags$link(rel = 'stylesheet', type = 'text/css', href = 'styles.css'),
             tags$style(type="text/css", "body {padding-top: 40px;}"),
             pageWithSidebar(
               ##  TITLE ##
               headerPanel(
                 img(src="img/logo.png", height="160")
               ),
               
               ## SIDE PANEL CONTENT ##
               sidebarPanel(
                 tags$head(tags$style(
                   type = 'text/css',
                   'form.well { max-height: 1600px; overflow-y: auto; }'
                 )),
                 
      ## SPECIFIC TO TREE LANDSCAPE EXPLORER ##
                 conditionalPanel(condition = "$('li.active a').first().html()== 'Tree landscape explorer'",
                                  ## INPUT
                                  ## choice of type of data source
                                  img(src="img/line.png", width="100%"),
                                  h2(HTML('<font color="#6C6CC4" size="6"> > Input </font>')),
                                  radioButtons("datatype", HTML('<font size="4"> Choose data source:</font>'),
                                               list("Example: woodmice trees"="expl","Input file"="file")),
                                  
                                  ## choice of dataset if source is an example
                                  ## PUT THIS BACK IN IF WE INTRODUCE MORE EXAMPLES
                                  #conditionalPanel(condition = "input.datatype=='expl'",
                                  #                 selectInput("dataset",
                                  #                             HTML('<font size="4"> Select an example dataset:</font>'),
                                  #                             choices=c("woodmiceTrees"))
                                  #),
                                  
                                  ## choice of dataset if source is a file
                                  conditionalPanel(condition = "input.datatype=='file'",
                                                   fileInput("datafile", p(HTML(' <font size="4"> Choose input file:</font>'), br(),
                                                                           strong("accepted formats:"), br(),
                                                                           em("- multiphylo"), "saved from R (.RData/.rda)", br(),
                                                                           em("- nexus"), "file (.nex/.nexus)")
                                                   ),
                                                   checkboxInput("randSamp","Randomly sample from the trees?", value=TRUE),
                                                   bsTooltip("randSamp", "For large sets of trees and/or trees with many tips the app may be slow, so beginning the analysis with a small random sample is recommended.", 
                                                             placement = "right", trigger = "hover", options = NULL),
                                                   conditionalPanel(
                                                     condition="input.randSamp",
                                                     sliderInput("sampleSize", "Size of random sample:", value=10, min=10, max=300, step=10)
                                                   )
                                  ),

                                  
    ## ANALYSIS
                                  img(src="img/line.png", width="100%"),
                                  h4(HTML('<font color="#6C6CC4" size="6"> > Analysis </font>')),
                                  
                                  ## choose metric
                                  selectInput("treemethod", "Choose a tree summary:",
                                              choices=c(
                                                "Kendall Colijn Metric" = "metric",
                                                "Patristic distances" = "patristic",
                                                "Number of nodes" = "nNodes",
                                                "Abouheif's metric" = "Abouheif",
                                                "Sum of direct descendents" = "sumDD")),
                                  
                                  ## lambda, axes
                                  uiOutput("lambda"),
                                  bsTooltip("lambda","When lambda=0 trees are compared topologically; increasing lambda gives more focus to branch lengths"),
    
                                  conditionalPanel(
                                    condition="input.plotType==1",
                                      uiOutput("naxes")
                                      
                                      ## Future: highlight median trees (if plotType==1)
                                      #checkboxInput("showMedians", label=strong("Highlight median tree(s)?"), value=FALSE)
                                  ),
                                  
                                  
                                  ## find clusters?
                                  checkboxInput("findGroves", label=strong("Identify clusters?"), value=FALSE),
                                  bsTooltip("findGroves","Statistical tools for choosing an appropriate clustering method and number of clusters will be added to treescape soon.", placement="right"),
    
                                  
                                  conditionalPanel(
                                    ## condition
                                    condition="input.findGroves",
                                    
                                    ## clustering method
                                    selectInput("clustmethod", "Clustering method:",
                                                choices=c(
                                                  "Ward" = "ward.D2",
                                                  "Single" = "single",
                                                  "Complete" = "complete",
                                                  "UPGMA" = "average")),
                                    
                                    ## number of clusters
                                    uiOutput("nclust")
                                    ),
                                  
                                  ## relevant if method = KC metric, allow tip emphasis
                                  conditionalPanel(
                                    condition="input.treemethod=='metric'",
                                    ## Emphasise tips
                                    checkboxInput("emphTips", label=strong("Emphasise tips?"), value=FALSE),
                                    bsTooltip("emphTips","Choose tips to emphasise or de-emphasise: the vector elements corresponding to these tips are multiplied by the weight chosen below.", placement="right"),  
                                    ## if tip emphasis is chosen, provide options:
                                      conditionalPanel(
                                        condition="input.emphTips",
                                    
                                        uiOutput("whichTips"),
                                    
                                        sliderInput("emphWeight", "Weight of emphasis", value=2,min=0.1,max=100)
                                        )
                                  ),
                                  
              ## AESTHETICS
                                  img(src="img/line.png", width="100%"),
                                  
                                  h2(HTML('<font color="#6C6CC4" size="6"> > Aesthetics </font>')),
                                  
                                  ## tree landscape or compare to single reference tree
                                  radioButtons("plotType", "View",
                                               choices=c("Full tree landscape"=1,"Distances from a reference tree"=2),
                                               selected=1),
                                  bsTooltip("plotType", "Choose whether to view the relative distances between all trees, or a 1-dimensional plot of their distances from a fixed reference tree"),
                                  
                                  ## Dimensions (3D possible if 3 or more axes retained, and full tree landscape)
                                  conditionalPanel(condition="input.naxes>2",
                                         conditionalPanel(
                                           condition="input.plotType==1",
                                                   radioButtons("plot3D", "Plot dimensions",
                                                                choices=c("2D"=2,"3D"=3),
                                                                selected=2)
                                         )
                                  ),
                                  
                                  # if plotType=1, pick the axes to view:
                                  conditionalPanel(    
                                    condition="input.plotType==1",
                                  ## select first axis to plot
                                  numericInput("xax", "Indicate the x axis", value=1, min=1, max=3),
                                  
                                  ## select second axis to plot
                                  numericInput("yax", "Indicate the y axis", value=2, min=1, max=3),
                                  bsTooltip("yax", "If multiple MDS axes have been retained, any combination of axes can be viewed"),
                              
                                  ## if in 3D, need a z axis:
                                  conditionalPanel(condition="input.plot3D==3",
                                                   numericInput("zax", "Indicate the z axis", value=3, min=1, max=3)
                                                   )
                                  ),
    
                                  ## aesthetics for treel landscape view
                                  conditionalPanel(
                                    condition="input.plotType==1",
                                  ## type of graph (if clusters detected, and 2D)
                                  conditionalPanel(
                                    ## condition
                                    condition="input.findGroves",
                                    conditionalPanel(
                                      condition="input.plot3D==2",
                                    
                                      ## type of plot
                                      radioButtons("scattertype", "Type of scatterplot",
                                                 choices=c("chull","ellipse"),
                                                 selected="chull")
                                    )
                                  ),
                                  
                                  ## symbol size
                                  sliderInput("pointsize", "Size of the points", value=1, min=0, max=10, step=0.2),
                                  
                                  ## more 2D aesthetics
                                  conditionalPanel(
                                    condition="input.plot3D==2",
                                    ## display labels
                                    checkboxInput("showlabels", label="Display tree labels?", value=FALSE),
                                  
                                    ## optimize labels?
                                    conditionalPanel(
                                      ## condition
                                      condition="input.showlabels",
                                      checkboxInput("optimlabels", label="Optimize label position?", value=FALSE),
                                    
                                      ## label size
                                      sliderInput("labelsize", "Size of the labels", value=1, min=0, max=10, step=0.2)
                                      ),
                                  
                                  
                                      ## add screeplot of MDS
                                      ## Could add this to 3d plot? 
                                      selectInput("screemds", "Position of the MDS screeplot:",
                                              choices=c("None" = "none",
                                                        "Bottom right" = "bottomright",
                                                        "Bottom left" = "bottomleft",
                                                        "Top right" = "topright",
                                                        "Top left" = "topleft"),
                                              selected="bottomleft")
                                  
                                  )
                                  ),
    
                                  # if plotType=2, option to stretch
                                  conditionalPanel(
                                    condition="input.plotType==2",
                                    
                                    uiOutput("selectedRefTree"),
                                 
                                    sliderInput("stretch", "Height of plot (pixels)", value=1600, min=800, max=12800, step=200)
                                 
                                  ),
    
                                  
                                  ## choose color palette (if clusters detected)
                                  conditionalPanel(
                                    ## condition
                                    condition="input.findGroves",
                                    
                                    selectInput("palette", "Palette for the clusters",
                                                choices=c("funky", "spectral",
                                                          "seasun", "lightseasun", "deepseasun",
                                                          "rainbow", "azur", "wasp"),
                                                selected="funky")
                                  ),
                                  
                                conditionalPanel(
                                  condition="input.plotType==1",
                                  conditionalPanel(
                                    condition="input.plot3D==2",
                                    ## Could add this to 3d?
                                    ## choose background color
                                    jscolorInput("bgcol", "Background color", value="#FFFFFF", close=TRUE)
                                    ),
                                  ## choose label colors
                                  jscolorInput("labcol", "Label / point color", value="#1B2266", close=TRUE)
                                  
                                ),
    
                              
    
                                  
                                  br(),br(),br(),br(),br(),br(),br(), # add some blank space at the end of side panel
                                  br(),br(),br(),br(),br(),br(),br(), # add some blank space at the end of side panel
                                  width=4)), # end conditional panel and sidebarPanel; width is out of 12
               ## MAIN PANEL
               mainPanel("",
                         
                         # TITLE #
                         h2(HTML('<font color="#6C6CC4" size="6"> Tree landscape explorer </font>')),
                         br(),br(),
                         
                         ## function I was using for testing:
                         #verbatimTextOutput("plot_click"),
                         
                         # Removed:
                         #verbatimTextOutput("caption"),
                         
                         uiOutput("treescapePlot"),
                         
                         br(), br(),
                         
                         ## OUTPUT (save)
                         img(src="img/line.png", width="400px"),
                         h2(HTML('<font color="#6C6CC4" size="6"> > Output </font>')),
                         ## save MDS plot
                         downloadButton("downloadMDS", "Save treescape image"),
                         
                         ## save trees to nexus file
                         downloadButton('exporttrees', "Save trees to nexus file"),
                         
                         ## save results to csv
                         downloadButton('exportrestocsv', "Save results (MDS+clusters) to csv file"),
                         
                         ## save results to RData
                         downloadButton('exportrestordata', "Save results (MDS+clusters) to R object")
               ) # end mainPanel
             ) # end page with sidebar
             ), # end tabPanel
             tabPanel("Tree viewer",
                      tags$link(rel = 'stylesheet', type = 'text/css', href = 'styles.css'),
                      tags$style(type="text/css", "body {padding-top: 40px;}"),
                      pageWithSidebar(
                        ##  TITLE ##
                        headerPanel(
                          img(src="img/logo.png", height="160")
                        ),
                       
                        ## SIDE PANEL CONTENT ##
                        sidebarPanel(
                          tags$head(tags$style(
                            type = 'text/css',
                            'form.well { max-height: 1600px; overflow-y: auto; }'
                          )),
                          
                          ## INPUT
                          ## choice of tree type
                          img(src="img/line.png", width="100%"),
                          h2(HTML('<font color="#6C6CC4" size="6"> > Input </font>')),
                          
                          radioButtons("treeChoice", "Tree selection",
                                       choices=c("Median tree"="med","General tree selection"="gen"),
                                       selected="med", width="100%"),
                          bsTooltip("treeChoice", "A geometric median tree is plotted by default. If clusters have been identified, the median for each can be viewed. Alternatively, any individual tree can be plotted."),
                          
                          conditionalPanel(condition = "input.treeChoice=='med'",
                                         selectInput("selectedMedTree", "Median tree from:", 
                                                       choices=c("All trees"="all"))
                          ),
                          
                          conditionalPanel(condition = "input.treeChoice=='gen'",
                                           uiOutput("selectedGenTree")
                          ),
                          
                          ## old version of tree selector
                          #textInput("selectedTree", "Choose tree (number or label)", value = ""),
                          

                          ## TREE AESTHETICS
                          img(src="img/line.png", width="100%"),
                          h2(HTML('<font color="#6C6CC4" size="6"> > Aesthetics </font>')),
                          
                          ## condition on tree being displayed
                          conditionalPanel(condition = "input.selectedTree!=''",
                                           ## use edge lengths?
                                           checkboxInput("edgelengths", label="Use original branch lengths?", value=TRUE),
                                           
                                           ## ladderize
                                           checkboxInput("ladderize", label="Ladderize the tree?", value=TRUE),
                                                            
                                           ## type of tree
                                           radioButtons("treetype", "Type of tree",
                                                         choices=c("phylogram","cladogram", "fan", "unrooted", "radial"),
                                                         selected="phylogram", width="100%"),
                                                            
                                           ## tree direction
                                           radioButtons("treedirection", "Direction of the tree",
                                                         choices=c("rightwards", "leftwards", "upwards", "downwards"),
                                                         selected="rightwards", width="100%"),
                                           
                                           ## tip labels
                                           checkboxInput("showtiplabels", label="Display tip labels?", value=TRUE),
                                           
                                                         conditionalPanel(condition="input.showtiplabels",
                                                                          ## tip label font
                                                                          selectInput("tiplabelfont", "Tip label font", 
                                                                                      choices=c("Plain"="1","Bold"="2","Italic"="3","Bold italic"="4")),
                                                                               
                                                                          ## tip label size
                                                                          sliderInput("tiplabelsize", "Size of the tip labels", value=1, min=0, max=5, step=0.1)
                                                                                 
                                                                          ),
                                                            
                                                            
                                            ## edge width
                                            sliderInput("edgewidth", "Width of the edges", value=2, min=1, max=20, step=0.2),
                                                    
                                            ## edge colour
                                            jscolorInput("edgecolor", "Edge colour", value="#000000", close=TRUE),
                                                            
                                           
                                            br(),br(),br(),br(),br(),br(),br(), # add some blank space at the end of side panel
                                            br(),br(),br(),br(),br(),br(),br(), # add some blank space at the end of side panel
                                           width=4)), # end conditional panel and sidebarPanel; width is out of 12
                       
                        ## MAIN PANEL
                        mainPanel("",
                                  # TITLE #
                                  h2(HTML('<font color="#6C6CC4" size="6"> Tree viewer </font>')),
                                  br(),br(),
                                  
                                  ## conditional panel: plot tree if needed
                                  conditionalPanel(condition = "input.selectedTree!=''",
                                                   plotOutput("tree", height = "800px"),
                                                   
                                                   br(), br(),
                                                   
                                                   ## OUTPUT (save)
                                                   img(src="img/line.png", width="400px"),
                                                   h2(HTML('<font color="#6C6CC4" size="6"> > Output </font>')),
                                                   downloadButton("downloadTree", "Save tree image"),
                                                   
                                                   br(), br(), br(), br(), br(), br()
                                  ) # end single tree conditional panel  
                        ) # end mainPanel
                      ) # end page with sidebar
                      ), # end tabPanel "Tree Viewer"
                      tabPanel("densiTree viewer",
                               tags$link(rel = 'stylesheet', type = 'text/css', href = 'styles.css'),
                               tags$style(type="text/css", "body {padding-top: 40px;}"),
                               pageWithSidebar(
                                 ##  TITLE ##
                                 headerPanel(
                                   img(src="img/logo.png", height="160")
                                 ),
                                 
                                 ## SIDE PANEL CONTENT ##
                                 sidebarPanel(
                                   tags$head(tags$style(
                                     type = 'text/css',
                                     'form.well { max-height: 1600px; overflow-y: auto; }'
                                   )),
                                   
                                   ## INPUT
                                   ## choice of tree type
                                   img(src="img/line.png", width="100%"),
                                   h2(HTML('<font color="#6C6CC4" size="6"> > Input </font>')),
                                   
                                   ## add densiTree selector (gets updated to number of clusters by )
                                   selectInput("selectedDensiTree", "Choose collection of trees to view in densiTree plot", 
                                               choices=c("Choose one"="","All trees"="all"), width="100%"),
                                   #h2(HTML('<font color="#6C6CC4" size="2"> Note: this can be slow for large sets of trees </font>')),
                                   
                                   bsTooltip("selectedDensiTree", "View all trees together in a densiTree plot. If clusters have been identified, the set of trees from a single cluster can be plotted. Note this function can be slow if many trees are included.", placement="bottom"),
                                   
                                   
                                   ## DENSITREE AESTHETICS
                                   img(src="img/line.png", width="100%"),
                                   h2(HTML('<font color="#6C6CC4" size="6"> > Aesthetics </font>')),
                                   
                                   conditionalPanel(condition = "input.selectedDensiTree!=''",
                                                    
                                                    ## alpha (semitransparency of edges)
                                                    sliderInput("alpha", "Transparency of edges", value=0.5, min=0, max=1, step=0.05),
                                                    
                                                    
                                                    checkboxInput("scaleX", label="Scale trees to equal heights?", value=FALSE),                 
                                                    
                                   br(),br(),br(),br(),br(),br(),br(), # add some blank space at the end of side panel
                                   br(),br(),br(),br(),br(),br(),br(), # add some blank space at the end of side panel
                                   width=4)), # end conditional panel and sidebarPanel; width is out of 12
                               
                                 ## MAIN PANEL
                                 mainPanel("",
                                           
                                           # TITLE #
                                           h2(HTML('<font color="#6C6CC4" size="6"> densiTree viewer </font>')),
                                           br(),br(),
                                           
                                           ## conditional panel: plot tree if needed
                                           conditionalPanel(condition = "input.selectedDensiTree!=''",
                                                            
                                                            plotOutput("densiTree", height = "800px"),
                                                            
                                                            br(), br(),
                                                            
                                                            ## OUTPUT (save)
                                                            img(src="img/line.png", width="400px"),
                                                            h2(HTML('<font color="#6C6CC4" size="6"> > Output </font>')),
                                                            downloadButton("downloadDensiTree", "Save densiTree image"),
                                                            
                                                            br(), br(), br(), br(), br(), br()
                                                            ) # end densiTree conditional panel
                                 ) # end of main panel "Multi-tree viewer"
                               ) # end of page with sidebar
                               ),# end of tabPanel "densiTree viewer"
    ## HELP SECTION
    tabPanel("Help",
             tags$style(type="text/css", "body {padding-top: 40px;}"),
             HTML(paste(readLines("www/html/help.html"), collapse=" "))
    ),
    
    ## SERVER INFO ##
    tabPanel("System info", 
             tags$style(type="text/css", "body {padding-top: 40px;}"),
             verbatimTextOutput("systeminfo"))
                     ) # end of tabsetPanel
) # end of Shiny UI