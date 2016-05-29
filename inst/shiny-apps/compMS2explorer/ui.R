library(CompMS2miner)
library(shiny)
library(igraph)
library(rhandsontable)

shiny::shinyUI(shiny::fluidPage(
  shiny::titlePanel("Composite MS2 Explorer "), 
  shiny::fluidRow(shiny::column(width=2, 
                                shiny::tags$b(shiny::h4('Search CompMS2 results:')),
                                shiny::actionButton("goButton", "Search"),
                                shiny::tags$b(shiny::h4("Search Options:")),
                                shiny::checkboxInput("All_Features", "All features", value=TRUE), 
                                shiny::textInput("mass_to_charge", label=shiny::tags$b("mass to charge ratio : "), value=447.3468), 
                                shiny::numericInput("mass_accuracy", label=shiny::tags$b("mass accuracy (ppm) : "), value=10, min=0.1, max=10000, step=0.1), 
                                shiny::textInput("retentionTime", label=shiny::tags$b("retention time (seconds) : "), value=876.2210), 
                                shiny::numericInput("RTtolerance", label=shiny::tags$b("retention time tolerance (+/- seconds) : "), value=5, min=1, max=700, step=1), 
                                shiny::selectInput("SubStrTypes", label=shiny::tags$b("Substructure types to include (multiple combinations) :"), choices=c(" ", SubStrType.inputs), multiple=T), 
                                shiny::selectInput("NotSubStrTypes", label=shiny::tags$b("Substructure types to exclude (multiple combinations) :"), choices=c(" ", SubStrType.inputs), multiple=T), 
                                shiny::selectInput("subStrAnnoTypes", label=shiny::tags$b("substructure type annotation (multiple combinations) :"), choices=c(" ", subStrAnno.inputs), multiple=T), 
                                shiny::numericInput("subStrAnnoThresh", label=shiny::tags$b("summed relative intensity threshold (substructure type). "), value=10, min=0, max=1000, step=1)
                                
                                #selectInput("Possible_contaminants", label=tags$b("Possible contaminant/false positive annotations to exclude (name_frequency) :"), choices=c(" ", as.character(AnnoFreq$name)), multiple=T), 
                                
  ), 
  shiny::column(width=3,
                shiny::textInput("DB_match_name", "Or search for database matches by name", value="e.g. p-cresol"),
                shiny::radioButtons("DB_match_table", 'Search in:', choices=c('DB Annotations', 'Best Annotations')), 
                shiny::actionButton("DBbutton", "Search"), 
                # shiny::br(),
                # shiny::br(),
                # shiny::tags$b("Save comments and best candidates"),
                # shiny::actionButton("CommentButton", "SAVE"),
                # shiny::br(),
                
                #                br(), br(), sliderInput("Loess_smooth", "Select the Loess smoother span :", min=0.01, max=1, value=0.66), 
                #                numericInput("CorrCoefMin", "Minimum correlation coefficient :", min=0.3, max=0.99, value=0.3, step=0.01), 
                shiny::br(), 
                shiny::uiOutput("FeatureNames"), 
                shiny::br(), 
                shiny::h5(shiny::textOutput("matchSummaryText")), 
                shiny::tableOutput("matchSummary"), 
                shiny::uiOutput("tabbedPanelButton"), 
                shiny::br(), 
                shiny::uiOutput("Plot_select"),
                shiny::tags$b("metID comments"),
                shiny::br(),
                shiny::tags$b('(close app and save comments to .csv)'),
                shiny::br(),
                shiny::actionButton("CloseAppBtn", "Close"),
                shiny::br(),
                # metId comments
                
                rHandsontableOutput("hot")
  ), 
   shiny::column(width=7, 
                shiny::conditionalPanel(condition="input.tabbedPanelButton >'0' & input.FeatureNames!='No MS2 features found'", 
                                        shiny::tabsetPanel(
                                          shiny::tabPanel('Composite MS2 plot', 
                                                          shiny::verbatimTextOutput("compMS2Hover"),   
                                                          shiny::plotOutput("MS2_plot", width = "800px",  height = "600px", brush = 'compMS2_brush', hover='compMS2_hover'), 
                                                          verbatimTextOutput("compMS2tableInfo")),
                                          shiny::tabPanel('Overview plot', 
                                                          shiny::plotOutput("overview_plot", width = "1000px",  height = "600px", brush = 'overview_brush'), 
                                                          verbatimTextOutput("overviewtableInfo")),
                                          shiny::tabPanel('Correlation network plot',
                                                          shiny::verbatimTextOutput("nNodesEdges"),
                                                          shiny::plotOutput("network_plot", width = "1000px",  height = "800px", brush = 'network_brush')),                          
                                          # shiny::tabPanel("Composite MS2 plot", rCharts::chartOutput("Raw_data_plot", "PolyCharts")), #, #plotOutput("MS2_plot", width = "800px",  height = "600px")),  
                                          shiny::tabPanel("MS2 spectrum table", DT::dataTableOutput(outputId="MS2_data")), 
                                          shiny::tabPanel("MS1 MS2 match summary", DT::dataTableOutput(outputId="metadata")), 
                                          shiny::tabPanel("DB Annotations",  DT::dataTableOutput(outputId="DB.results")), 
                                          shiny::tabPanel("Best Annotations",  DT::dataTableOutput("BestCandidates")), 
                                          shiny::tabPanel("Substructure Annotations",  DT::dataTableOutput("BestSubStrAnno")), 
                                          #                                   tabPanel("Best candidates", dataTableOutput("BestCandidate")), 
                                          shiny::tabPanel("PubMed Word Cloud",  shiny::uiOutput("wordCloudSelect"),  shiny::uiOutput("nPMIDAbstracts"),  shiny::verbatimTextOutput("WordCloudText"), 
                                                          shiny::uiOutput("nRandomArticles"),  shiny::tableOutput(outputId="WordCloudTable"),  shiny::plotOutput("WordCloud", width = "800px",  height = "600px")), 
                                          shiny::tabPanel("MetFrag results",  shiny::tableOutput(outputId="MetFragTable"))#, 
                                          # 
                                          # shiny::tabPanel("Substructure type table"), #, tableOutput(outputId="SubStr_type")), 
                                          # shiny::tabPanel("InterFeature Correlation"), #, chartOutput("InterFeatureCorr", "morris")), #uiOutput("CorrCoefMin"), 
                                          # shiny::tabPanel("Chemical Similarity Scores"), #uiOutput("Compound"), uiOutput("IntraInter"), chartOutput("ChemicalSimilarity", "morris")), 
                                          # shiny::tabPanel("Comments"), #,  shiny::uiOutput("CommentSectionHeader"),  shiny::uiOutput("UserComments"))#, uiOutput("CommentButtonUI")
                                          # shiny::tabPanel("LogD RT plot")#, uiOutput("RTwindow"), chartOutput("LogD_RT_plot", "PolyCharts"))
                                        ))
  ))
 
  
  #         fluidRow(column(h4(paste0("Results summary ", "(nMatches=", TotalFeatures, ",  nCompositeSpectra=", TotalCompSpectra, ")")), width=6, 
  #                         dataTableOutput("MS2_features_detected")))
)) # end CompMS2 shiny UI