#' @import shiny
#' @import highcharter
#' @import shinycssloaders
#' @importFrom colourpicker colourInput updateColourInput
#' @importFrom DT renderDT DTOutput formatRound datatable
#' @importFrom bslib bs_theme
betASapp_ui <- function(){
  # :::: Variables ::::
  availabletools      <- c("vast-tools", "rMATS","whippet")
  yAxisStats          <- c("Pdiff (probability of differential splicing)", "F-statistic (median(|between|)/median(|within|))", "False discovery rate (FDR)")
  yAxisStats_multiple <- c("Pdiff (probability that |between| > |within|)", "F-statistic (median(|between|)/median(|within|))")
  eventTypesVT          <- c("Exon skipping (ES)"="EX", "Intron retention (IR)"="IR", "Alternative splice site (Altss)"="Altss")
  eventTypesWhippet     <- c("Core Exon (CE)"="CE",
                             "Alternative Acceptor splice site (AA)"="AA",
                             "Alternative Donor splice site (AD)"="AD",
                             "Intron retention (IR)"="RI",
                             "Tandem transcription start site (TS)"="TS",
                             "Tandem alternative polyadenylation site (TE)"="TE",
                             "Alternative First exon (AF)"="AF",
                             "Alternative Last exon (AL)"="AL",
                             "Circular back-splicing (BS)"="BS")
  exEventNames        <- c("HsaEX0007927", "HsaEX0032264", "HsaEX0039848", "HsaEX0029465", "HsaEX0026102", "HsaEX0056290", "HsaEX0035084", "HsaEX0065983", "HsaEX0036532", "HsaEX0049206")
  pastelColors        <- c("#FF9AA2", "#FFB7B2", "#FFDAC1", "#E2F0CB", "#B5EAD7", "#C7CEEA", "#FBE2FD", "#D9ECFE")

  # Themes
  # dark <- bs_theme(bootswatch = "darkly")
  # light <- bs_theme(bootswatch = "minty")

  # layout function: sets up the basic virtual structure of the page
  # defines user interface (ui)
  ui <- fluidPage(

    # tags$img(src="betasBanner.png", align = "right", height='80px', width='1000px'),

    # theme = bs_theme(bootswatch = "sketchy"),
    theme = bs_theme(bootswatch = "minty"),

    # theme = light,

    tags$head(
      tags$style(
        HTML(".shiny-notification {
             position:fixed;
             top: calc(10%);
             left: calc(50%);
             }"),
        'body {
      font-size: 16px}'
      )),

    tags$head(
      tags$style(HTML("
    .table.dataTable tbody td.active,
    .table.dataTable tbody tr.active td {
      background-color: #F3969A !important;
      color: white !important;  # Changing the text color to black for better contrast
    }
  "))
    ),



    # checkboxInput("dark_mode", "Dark mode", FALSE),

    titlePanel("betAS: intuitive visualisation of differential alternative splicing"),

    #h5(strong(div("(development version)", style = "color: #E0E0E0"))),
    h5(strong(div("v1.1.0", style = "color: #949494"))),

    tabsetPanel(

      tabPanel("Import inclusion levels",  icon = icon("play"),

               sidebarLayout(fluid = TRUE,

                             sidebarPanel(style = "max-height: 100%",
                                          # h3("Exploratory analysis of inclusion levels"),
                                          h4("Exploratory analysis of inclusion levels"),
                                          h5(icon("file-import", style = "color: #000000"),"Import data"),

                                          selectInput("sourcetool", "Select example dataset:", choices = c("Dataset 1 (obtained using vast-tools)" = "vast-tools1",
                                                                                                           "Dataset 2 (obtained using vast-tools)" = "vast-tools2",
                                                                                                           "Dataset 2 (obtained using rMATS)" = "rMATS",
                                                                                                           "Dataset 2 (obtained using whippet)" = "whippet")),
                                          "Alternatively, upload table with inclusion level quantification (e.g. PSI):",
                                          p("betAS currently supports inclusion level tables from vast-tools ",code("(INCLUSION_LEVELS_FULL*.tab.gz)", style = "font-size:12px; color: #AAAAAA"),
                                            ", rMATS ",code("(*.MATS.JC.txt)", style = "font-size:12px; color: #AAAAAA")," and whippet",code("(*.psi.gz).", style = "font-size:12px; color: #AAAAAA"),
                                            "Only a single file is supported when using rMATS or vast-tools results. Results from vast-tools module",code("tidy", style = "font-size:12px; color: #AAAAAA"),"not supported.
                                            For whippet, please upload one file per sample (at least two samples).", style = "font-size:12px; color: #AAAAAA"),


                                          fileInput("psitable", NULL, placeholder = "No file selected (max. 100MB)", multiple=T, accept = c(".tab",".txt",".psi",".gz")), #accept = c(".tab",".txt",".psi",".gz")

                                          hr(),

                                          h5(icon("file", style = "color: #000000"),"Your files:"),

                                          uiOutput("datasetInfo"),

                                          DTOutput("files"),

                                          hr(),

                                          h5(icon("filter", style = "color: #000000", class="Light"),"Filter events from loaded table(s):"),

                                          conditionalPanel(
                                            condition = "output.showchecks",
                                            checkboxGroupInput("types", label = "Event types to consider:", selected = c("EX"),  choices = eventTypesVT)
                                          ),
                                          sliderInput("psirange", "PSI values to consider:", value = c(1, 99), min = 0, max = 100, step=1),


                                          helpText((em("Consider only alternative splicing events with all PSI values within this range."))),

                                          p(" "),

                                          tags$div(numericInput("minNreads", "Minimum number of corrected reads per event:", 10, min = 0.5, width= "100%"),  style="display:inline-block"),
                                          tags$div(actionButton("buttonminreads", "Filter"),  style="display:inline-block"),
                                          verbatimTextOutput("test_box"),



                                          h6(textOutput("textTotalNumberEvents")),
                                          highchartOutput("eventsPiechart")

                             ),

                             mainPanel(
                               shinycssloaders::withSpinner(plotOutput("plot", height = "1000px"), type = 8, color = "#FF9AA2", size = 2),
                               #downloadButton("downloadPlot", "Download Plot"),
                               h6(textOutput("TextToolInfo")),
                               DTOutput("selected_file_table" )
                             )
               )
      ),

      tabPanel("Group definition", icon = icon("gears"),

               sidebarLayout(fluid = TRUE,

                             sidebarPanel(
                               selectInput("indbetaseventid", "Select alternative splicing event to plot:", choices = NULL),
                               h6(textOutput("EventInfoWhippet")),


                               uiOutput("url"),

                               h5(textOutput("samplesTabletext")),

                               DTOutput("sampleTable"),

                               h4("Organise samples into groups/conditions"),

                               helpText(tags$ul(
                                 tags$li("one beta distribution will be estimated per sample"),
                                 tags$li("groups defined can be used for differential splicing analyses")
                               )),

                               tags$ol(

                                 tags$li(h5("Add/delete groups by sample selection:")),

                                 textInput(inputId = "groupName", label = "Insert group name:", placeholder = "e.g. Control"),
                                 selectInput("groupSamples", "Choose samples in group:", choices = NULL, multiple = TRUE),
                                 colourInput("groupColor", "Select group color:", value = "#89C0AE"),
                                 actionButton("newGroup", "Create new group", icon = icon("layer-group"), class = "btn-secondary"),

                                 tags$li(h5("Automatic group creation:")),

                                 tags$ol(

                                   tags$li(h6("Group samples based on sample name similarities")),
                                   actionButton("findGroups", "Automatic group(s)", icon = icon("wand-sparkles"), class = "btn-info"),


                                   conditionalPanel(
                                     condition = 'output.showGroupingFeatureOption',
                                     tags$li(h6("Group samples based on a given feature")),
                                     selectInput("groupingFeature", label = "Select feature to group samples by:", choices = NULL),
                                     actionButton("findGroupsBasedSampleTable", label = "Feature-associated group(s)", icon = icon("robot"), class = "btn-info")
                                   )

                                 )),

                               h5("Current groups defined"),
                               DTOutput("groupsTable"),
                               helpText(em("Choose one table row for group deletion")),

                               actionButton("deleteGroups", "Delete group(s)", icon = icon("eraser"), class = "btn-danger")

                             ),

                             mainPanel(
                               shinycssloaders::withSpinner(plotOutput("densities", height = "1200px"), type = 8, color = "#FF9AA2", size = 2)
                             )
               )
      ),

      tabPanel("Differential alternative splicing", icon = icon("flag-checkered"),

               sidebarLayout(fluid = TRUE,

                             sidebarPanel(
                               h4("Differential splicing in groups/conditions"),
                               helpText(tags$ul(
                                 tags$li("All events from 'Import inclusion levels'"),
                                 tags$li("Groups created in 'Group definition'")
                               )),

                               selectInput("groupA", "Choose group A:", NULL),

                               selectInput("groupB", "Choose group B:", NULL),

                               radioButtons("volcanoYAxis", label = "Choose significance statistic (Y-axis) to consider:", choices = yAxisStats),

                               actionButton("rundiffbetas", "Run betAS", icon = icon("person-running"), class = "btn-info"),
                               # Other cool icons: play | play-circle

                               h5("Explore differential splicing at the event level:"),
                               helpText(em("Brush over the plot for event selection")),
                               DTOutput("brushed_data"),
                               h6(""),
                               selectInput("eventidtoplot", "Alternative splicing event to plot:" , choices = NULL),

                               uiOutput("urlplot"),
                               actionButton("plotEvent", "Plot considered event", icon = icon("mound"), class = "btn-secondary"),

                             ),

                             mainPanel(

                               shinycssloaders::withSpinner(plotOutput("volcano",
                                                                       height = "900px",
                                                                       brush = brushOpts(id = "plot_brush")),
                                                            type = 8,
                                                            color = "#FF9AA2",
                                                            size = 2),

                               fluidRow(
                                 shinycssloaders::withSpinner(plotOutput("densitiesSelectedEvent", height = "400px"), type = 8, color = "#FF9AA2", size = 2)
                               ),

                               conditionalPanel(
                                 condition = 'output.showIndividualEventPlots',

                                 fluidRow(

                                   column(4,
                                          h4("P(A - B) > x:"),
                                          h6("Probability that estimated PSI of one group is greater than the other by x.")),

                                   column(4,
                                          h4("F-statistic:"),
                                          h6("Ratio between absolute differences 'between' and 'within' groups.")),

                                   column(4,
                                          h4("False discovery rate (FDR):"),
                                          h6("Probability of getting an absolute difference in PSI between groups greater than the observed, under the null hypothesis that all PSIs come from the same distribution."))

                                 )),

                               fluidRow(
                                 column(4,

                                        shinycssloaders::withSpinner(plotOutput("PDiffEventPlot", height = "400px", width = "400px"), type = 8, color = "#FF9AA2", size = 2)),

                                 column(4,

                                        shinycssloaders::withSpinner(plotOutput("FstatEventPlot", height = "400px", width = "400px"), type = 8, color = "#FF9AA2", size = 2)),

                                 column(4,

                                        shinycssloaders::withSpinner(plotOutput("FDREventPlot", height = "400px", width = "400px"), type = 8, color = "#FF9AA2", size = 2))

                               )

                             ) # end mainPanel


               ) # end sidebarLayout
      ), # end tabPanel

      tabPanel("Differential alternative splicing (mutiple groups)", icon = icon("layer-group"),

               sidebarLayout(fluid = TRUE,

                             sidebarPanel(
                               h4("Generalised differential alternative splicing across groups/conditions"),
                               # "Perform differential alternative splicing analysis for all events",
                               helpText(tags$ul(
                                 tags$li("All events from 'Import inclusion levels'"),
                                 tags$li("Groups created in 'Group definition'")
                               )),

                               h6(""),

                               radioButtons("volcanoYAxis_mult", label = "Choose significance statistic (Y-axis) to consider:", choices = yAxisStats_multiple),
                               actionButton("rundiffbetas_mult", "Run betAS (multiple groups)", icon = icon("person-running"), class = "btn-info"),

                               h5("Current groups"),
                               helpText(em("Defined groups can be edited in 'Group definition' tab")),
                               DTOutput("groupsTableMultiple"),

                               h6(""),

                               h5("Explore differential splicing at the event level:"),
                               helpText(em("Brush over the plot for event selection")),
                               DTOutput("brushed_data_mult"),

                               selectInput("eventidtoplot_mult", "Select alternative splicing event to plot:", choices = NULL),
                               uiOutput("urlplot_mult"),

                               h6(""),
                               actionButton("plotEvent_mult", "Plot considered event", icon = icon("mound"), class = "btn-secondary"),

                             ),

                             mainPanel(

                               shinycssloaders::withSpinner(plotOutput("volcano_mult",
                                                                       height = "900px",
                                                                       brush = brushOpts(id = "plot_brush_mult")),
                                                            type = 8,
                                                            color = "#FF9AA2",
                                                            size = 2),
                               conditionalPanel(
                                 condition = 'output.showIndividualEventPlotsmult',
                                 fluidRow(

                                   column(3,
                                          h5("Absolute differences 'between' and 'within' groups.")),

                                   column(9,
                                          h5("PSI quantifications per sample across groups.")))),

                               fluidRow(

                                 column(3,

                                        plotOutput("FstatEventPlot_multgroup", height = "600px", width = "400px")

                                 ),

                                 column(9,
                                        shinycssloaders::withSpinner(plotOutput("violins_multgroup",
                                                                                height = "600px"),
                                                                     type = 8,
                                                                     color = "#FF9AA2",
                                                                     size = 2)))

                             )
               )
      )

    ))

  return(ui)

}

#' @importFrom thematic thematic_shiny
#' @importFrom utils read.table
betASapp_server <- function(){

  # file size limit = 10MB

  options(highcharter.theme = hc_theme_smpl(tooltip = list(valueDecimals = 2)), shiny.maxRequestSize = 100 * 1024^2)

  # :::: Variables ::::
  # file                  <- "test/INCLUSION_LEVELS_FULL-Hsa32-hg19_to_test.tab.gz"
  # testTable             <- read.table(gzfile(file), sep="\t", header=TRUE, quote="")
  # maxDevSimulationN100  <- readRDS(url("http://imm.medicina.ulisboa.pt/group/distrans/SharedFiles/Mariana/Splicing&SenescenceFLEX/xintercepts_100incr_100cov_100trials.R"))
  # maxDevSimulationN100  <- readRDS("test/xintercepts_100incr_100cov_100trials.rds")
  data("maxDevSimulationN100")
  pastelColors          <- c("#FF9AA2", "#FFB7B2", "#FFDAC1", "#E2F0CB", "#B5EAD7", "#C7CEEA", "#FBE2FD", "#D9ECFE")

  eventTypesVT          <- c("Exon skipping (ES)"="EX", "Intron retention (IR)"="IR", "Alternative splice site (Altss)"="Altss")
  eventTypesWhippet     <- c("Core Exon (CE)"="CE",
                             "Alternative Acceptor splice site (AA)"="AA",
                             "Alternative Donor splice site (AD)"="AD",
                             "Intron retention (IR)"="RI",
                             "Tandem transcription start site (TS)"="TS",
                             "Tandem alternative polyadenylation site (TE)"="TE",
                             "Alternative First exon (AF)"="AF",
                             "Alternative Last exon (AL)"="AL",
                             "Circular back-splicing (BS)"="BS")

  default_VT_events           <- c("EX")
  default_Whippet_events      <- c("CE")

  # Minimum number of reads for each event in each sample; events with less than one read in at least one sample will be filtered out
  defaultminNreads              <-  10


  # specifies the behaviour of the app by defining a server function
  server <- function(input, output, session){
    thematic_shiny()

    # not working, but based on slide 25 from https://talks.cpsievert.me/20201014/#25
    # observe(session$setCurrentTheme(
    #   if (input$dark_mode) dark else light
    # ))


    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # A. Import Inclusion Levels -----------------------------------------------
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


    ## Functions ---------------------------------------------------------------

    # auxiliary variable used to render the main plot of the shiny app
    # as the plot is rendered whenever a dataset is updated we need to ensure that
    # the density plots are only rendered after the default parameters are all defined
    # i.e., to not be rendered twice when the dataset & tool change

    dataset_updated <- reactiveVal(FALSE)

    # Loaded main original table
    dataset <- reactive({


      if(is.null(input$psitable)){

        if(input$sourcetool == "vast-tools1"){

          #testTable <- readRDS(file = "test/INCLUSION_LEVELS_FULL-hg19-98-v251.rds")
          data("VT1_data_human")
          testTable <- VT1_data_human

        } else {

          showNotification("Please note that importing data may take a few minutes.",
                           closeButton = TRUE,
                           duration = 10,
                           type = c("message"))


          if(input$sourcetool == "vast-tools2"){

            #testTable <- readRDS(file = "test/INCLUSION_LEVELS_FULL-mm10-8-v251.rds")
            testTable <- getDataset(pathTables=NULL, tool = "vast-tools")

          }else if(input$sourcetool == "rMATS"){

            #testTable <- read.delim(file = "test/SE.MATS.JC.txt")
            testTable <- getDataset(pathTables=NULL, tool = "rMATS")

          } else if (input$sourcetool == "whippet"){

            #testTable <- readRDS(file = "test/listdfs_WHippet.rds")
            testTable <- getDataset(pathTables=NULL, tool = "whippet")

          }

          return(testTable)

        }
      }else {



        if(length(input$psitable$datapath) > 1 & length(grep(pattern = "[.]psi", x = input$psitable$name)) == length(input$psitable$datapath) ){



          files <- as.list(input$psitable$datapath)
          files[grep("[.]gz",files)] <- lapply(files[grep("[.]gz",files)],gzfile)
          loadingFile <- lapply(files,read.delim)

          #loadingFile <- lapply(as.list(input$psitable$datapath),fread)
          names(loadingFile) <- sapply(input$psitable$name, function(file) gsub("\\..*","",gsub(".*/","",file)))

        } else if (length(input$psitable$datapath) == 1 & length(grep(pattern = "[.]psi", x = input$psitable$name)) == length(input$psitable$datapath) ){

          showNotification(HTML(paste("Number of files not supported. Please select at least two files when using whippet inclusion tables. <br/> Please refresh the page if you intend to use the default datasets.",collapse = "<br/>")),
                           closeButton = TRUE,
                           duration = 120,
                           type = c("error"))
          return(NULL)

        } else if (length(input$psitable$datapath) > 1 & length(grep(pattern = "[.]psi", x = input$psitable$name)) != length(input$psitable$datapath) ){

          showNotification(HTML(paste("Input files not supported. Please select files from only one tool. <br/> Please refresh the page if you intend to use the default datasets.",collapse = "<br/>")),
                           closeButton = TRUE,
                           duration = 120,
                           type = c("error"))
          return(NULL)

        } else {

          if(length(grep(pattern = "[.]gz", x = input$psitable$datapath)) == 0){

            loadingFile <- input$psitable$datapath
          } else {
            loadingFile <- gzfile(input$psitable$datapath)
          }

          loadingFile <- read.delim(loadingFile)

        }

        return(loadingFile)
      }

    })

    sampleTable <- reactive({

      req(dataset())
      req(input$sourcetool)

      if(is.null(input$psitable) & input$sourcetool == "vast-tools1"){

        #samplesTable <- readRDS(file = "test/samplesTable.rds")
        data("VT1_metadata_human")
        samplesTable <- VT1_metadata_human

        return(samplesTable)

      } else if(is.null(input$psitable) & input$sourcetool == "vast-tools2"){

        #samplesTable <- readRDS(file = "test/samplesTable_Whippet.rds")
        data("VT2_metadata_mouse")
        samplesTable <- VT2_metadata_mouse

        return(samplesTable)

      } else if(is.null(input$psitable) & input$sourcetool == "rMATS"){

        #samplesTable <- readRDS(file = "test/samplesTable_rMATS.rds")
        data("rMATS_metadata_mouse")
        samplesTable <- rMATS_metadata_mouse

        return(samplesTable)

      } else if(is.null(input$psitable) & input$sourcetool == "whippet"){

        #samplesTable <- readRDS(file = "test/samplesTable_Whippet.rds")
        data("whippet_metadata_mouse")
        samplesTable <- whippet_metadata_mouse

        return(samplesTable)

      }


    })




    #Check data format. Currently supports rMATS, vast-tools and whippet tables

    sourcetool <- reactive({

      req(dataset())

      # vast-tools inclusion tables have (at least) the following columns: GENE	EVENT	COORD	LENGTH	FullCO	COMPLEX	Sample1	Sample1.Q (...)
      requiredcols_vasttools <- c("GENE","EVENT","COORD","LENGTH","FullCO","COMPLEX")

      # rMATS inclusion tables have the following columns: ID	GeneID	geneSymbol	chr	strand IJC_SAMPLE_1 SJC_SAMPLE_1	IJC_SAMPLE_2	SJC_SAMPLE_2
      requiredcols_rMATS <- c("ID","GeneID","geneSymbol","chr","strand","IJC_SAMPLE_1","SJC_SAMPLE_1","IJC_SAMPLE_2","SJC_SAMPLE_2")

      # Whippet inclusion tables have the following columns: "Gene" "Node" "Coord" "Strand" "Type" "Psi" "CI_Width" "CI_Lo,Hi" "Total_Reads" "Complexity" "Entropy" "Inc_Paths" "Exc_Paths" "Edges"
      # Note: when using read.delim for importing the files, it changes all "," to "."; Meaning, "CI_Lo,Hi" changes to "CI_Lo.Hi"
      # For sake of simplicity, we are not considering "CI_Lo,Hi" as mandatory
      requiredcols_whippet <- c("Gene", "Node", "Coord", "Strand", "Type" ,"Psi", "CI_Width", "Total_Reads", "Complexity", "Entropy" ,"Inc_Paths", "Exc_Paths")

      if (class(dataset())=="list"){

        # for Whippet, the output for dataset() is a list of the original inclusion tables (one per sample), which should all have the same column names
        inputTablecols <- unique(lapply(dataset(), colnames))[[1]]

      } else {

        inputTablecols <- colnames(dataset())

      }

      if(length(inputTablecols[inputTablecols %in% requiredcols_vasttools]) == length(requiredcols_vasttools)){

        sourcetool <- "vast-tools"

      } else if( length(inputTablecols[inputTablecols %in% requiredcols_rMATS]) == length(requiredcols_rMATS)){

        sourcetool <- "rMATS"

      } else if(length(inputTablecols[inputTablecols %in% requiredcols_whippet]) == length(requiredcols_whippet)){

        sourcetool <- "whippet"

      }else{

        showNotification(HTML(paste( "The provided data is not supported by betAS. <br/> Please confirm that your data matches the input requirements. <br/> Please refresh the page and try again.",collapse = "<br/>")),
                         closeButton = FALSE,
                         duration = 120,
                         type = c("error"))

        return(NULL)

      }
      return(sourcetool)

    })



    GetTable <- reactive({

      req(dataset())
      req(sourcetool())

      FormattedTable <- getEvents(incTable=dataset(), tool=sourcetool())
      return(FormattedTable)

    })


    # auxiliary variable to save existing events in a dataset
    getExistingEvents <- reactive({

      req(sourcetool())
      req(GetTable())

      existingEvents <- names(GetTable()$EventsPerType)

      if(sourcetool()=="vast-tools"){

        temp <- c()
        if(any(existingEvents %in% c("C1", "C2", "C3", "S", "MIC"))) temp <- c(temp,"EX")
        if(any(existingEvents %in% c("IR-C", "IR-S", "IR"))) temp <- c(temp,"IR")
        if(any(existingEvents %in% c("Alt3", "Alt5"))) temp <- c(temp,"Altss")
        existingEvents <- temp

      }

      return(existingEvents)

    })


    # auxiliary variable to choose only one event type to filter as default
    getDefaultEvents <- reactive({

      req(sourcetool())
      req(GetTable())

      getDefaultEvents <- names(GetTable()$EventsPerType)

      if(sourcetool()=="rMATS"){
         DefaultEvents <- getDefaultEvents
      }

      if(sourcetool()=="whippet"){
        if ("CE" %in% getDefaultEvents) DefaultEvents <- "CE"
        else DefaultEvents <- getDefaultEvents[1]
      }

      if(sourcetool()=="vast-tools"){

        if(any(getDefaultEvents %in% c("C1", "C2", "C3", "S", "MIC"))) {

          DefaultEvents <- "EX"

        } else {

          if(any(getDefaultEvents %in% c("IR-C", "IR-S", "IR"))) DefaultEvents <- "IR"

          if(any(getDefaultEvents %in% c("Alt3", "Alt5"))) DefaultEvents <- "Altss"


        }
      }

      return(DefaultEvents)

    })

    # Update choice values for inputs when the dataset changes
    observeEvent(GetTable(),{

      req(sourcetool())
      req(dataset())
      req(GetTable())

      if(sourcetool()=="vast-tools"){

        #updateCheckboxGroupInput(session, "types", label = "Event types to consider:", selected = default_VT_events,  choices = eventTypesVT)
        updateCheckboxGroupInput(session, "types", label = "Event types to consider:", selected = getDefaultEvents(),  choices = eventTypesVT)

      }else if(sourcetool()=="whippet"){
        updateCheckboxGroupInput(session, "types", label = "Event types to consider:", selected = getDefaultEvents (),  choices = eventTypesWhippet)
        #updateCheckboxGroupInput(session, "types", label = "Event types to consider:", selected = default_Whippet_events,  choices = eventTypesWhippet)
      }

    })

    observeEvent(dataset(), {

      updateNumericInput(session, "minNreads", "Minimum number of corrected reads per event:", defaultminNreads, min = 1)

    })


    # Because the change using updates is not immediate, we have to initialize it in order to not use the values from the previous dataset
    # auxiliary variable that changes whenever the dataset changes OR the user clicks on the button "Filter"

    # to only render the plots when the user clicks on the Filter button OR the first time the dataset changes OR if the selected events change
    rendervar <- reactiveVal(1)
    # reactive value for minimum number of reads to be filtered
    minNreads <- reactiveVal(defaultminNreads)
    # reactive value for the default values for event types to consider
    input_types <- reactiveVal(default_VT_events)


    # defaults, when the dataset changes
    observeEvent(dataset(), {

      # check whether input$types actually changes in the dataset update observer.
      # If it doesn't change, we can immediately set dataset_updated(FALSE)

      new_input_types <- getDefaultEvents()

      # Check if input$types would actually change
      if(!identical(new_input_types, input$types)) {
        dataset_updated(TRUE)
      } else {
        dataset_updated(FALSE)
      }

      # Update input_types and rendervar
      input_types(new_input_types)
      rendervar(rendervar() + 1)


    })

    # update render variable when the user clicks on the Filter button
    # Also update the value of the minimum reads to be filtered to the one in the input$minNreads field
    observeEvent(input$buttonminreads, {

      render <- rendervar() + 1
      rendervar(render)

      minNreads(input$minNreads)

    })

    # update render variable when the selected input types change
    # Also update the value of the event types to be filtered to the ones in the input$types field
    observeEvent(input$types, {

      # if no event is selected, change to one of the default ones
      if(is.null(input$types)){
        if (sourcetool() != "rMATS"){
          showNotification("Please select at least one event type. Updating to default events.", duration = 10, type = c("error"))

          if(sourcetool()=="vast-tools"){
            updateCheckboxGroupInput(session, "types", label = "Event types to consider:", selected = getDefaultEvents (),  choices = eventTypesVT)
            #updateCheckboxGroupInput(session, "types", label = "Event types to consider:", selected = default_VT_events,  choices = eventTypesVT)
          }else if(sourcetool()=="whippet"){
            updateCheckboxGroupInput(session, "types", label = "Event types to consider:", selected = getDefaultEvents (),  choices = eventTypesWhippet)
            #updateCheckboxGroupInput(session, "types", label = "Event types to consider:", selected = default_Whippet_events,  choices = eventTypesWhippet)
          }
        }

        return(NULL)
      }

      # to prevent filterTable() to update twice, as the input types will change when changing tool (and, thus, dataset)
      if (dataset_updated()) {
        dataset_updated(FALSE)
        return(NULL) # Exit this observer early if dataset was just updated
      }

      input_types(input$types)
      render <- rendervar() + 1
      rendervar(render)

    }, ignoreNULL = FALSE)


    # 1. Filter table to remove events

    filterTable <- eventReactive(rendervar(), {

      req(sourcetool())
      req(GetTable())
      req(input_types())
      req(minNreads())


      if(sourcetool() == "vast-tools"){

        selectedEventTypes <- c()

        if("EX" %in% input_types()) selectedEventTypes <- c("C1", "C2", "C3", "S", "MIC")
        if("IR" %in% input_types()) selectedEventTypes <- c(selectedEventTypes, "IR-C", "IR-S", "IR") # different versions of VT handle intron retention events differently
        if("Altss" %in% input_types()) selectedEventTypes <- c(selectedEventTypes, "Alt3", "Alt5")

      } else if(sourcetool() == "rMATS"){

        selectedEventTypes <- NULL

      } else if(sourcetool() == "whippet"){

        selectedEventTypes <- input_types()

      }

      filteredList <- filterEvents(InputList=GetTable(), types = selectedEventTypes, N= minNreads())

      # if none of the checked events are in the filtered data, update the event to one that exists, in order to prevent the APP from crashing
      if (all(! selectedEventTypes %in% names(filteredList$EventsPerType)) & sourcetool() != "rMATS"){
        showNotification(HTML("There are no events. <br> Please confirm in the piechart below <br> if these events exist in your dataset. <br> Updating data to consider at least one event."),
                         closeButton = TRUE,
                         duration = 10,
                         type = c("error")  )



        if(sourcetool()=="vast-tools"){

          dataset_updated(FALSE)
          updateCheckboxGroupInput(session, inputId = "types", selected = getDefaultEvents(),  choices = eventTypesVT)

        }else if(sourcetool()=="whippet"){
          dataset_updated(FALSE)
          updateCheckboxGroupInput(session, inputId = "types", selected = getDefaultEvents(),  choices = eventTypesWhippet)
        }

        return(NULL)

      }


      return(filteredList)

    })



    selectAlternatives <- reactive({

      req(sourcetool())

      # This will make sure that selectAlternatives is updated whenever input$psirange changes
      input$psirange

      # This will make sure that selectAlternatives is also updated when filterTable changes
      isolate(filterTable())

      alternativeList <- alternativeEvents(filteredList=req(filterTable()), minPsi = input$psirange[1], maxPsi = input$psirange[2])

      if(nrow(alternativeList$PSI) == 0 & !is.null(filterTable)){

        showNotification("There are no PSI values for the events selected within such range. Updating PSI range to all possible values.",
                         closeButton = TRUE,
                         duration = 5,
                         type = c("error"))

        updateSliderInput(inputId = "psirange", value = c(0, 100))

        return(NULL)

      }

      return(alternativeList)

    })

    # create a reactive expression
    psidataset <- reactive({

      filterTable()$PSI

    })

    qualdataset <- reactive({

      filterTable()$Qual

    })

    # create a reactive expression
    psifiltdataset <- reactive({

      selectAlternatives()$PSI

    })

    qualfiltdataset <- reactive({

      req(selectAlternatives())
      return(selectAlternatives()$Qual)

    })

    eventNumber <- reactive({

      req(selectAlternatives())

      return(nrow(psifiltdataset()))
    })

    ## Outputs -----------------------------------------------------------------


    output$datasetInfo <- renderText({

      req(sourcetool())
      req(dataset())

      if(is.null(input$psitable)){

        if (input$sourcetool == "vast-tools1"){
          url <-" https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6814/"
          style <-  "style = 'font-size:14px; color: #000000' "
          HTML(paste0("<p",style,",>You are using a subset of the <i>'Human RNA-seq time-series of the development of seven major organs'</i> public human dataset (<a target='_blank' href='", url, "'>E-MTAB-6814</a>) as an example,
                where the inclusion tables have been obtained using vast-tools.</p>"))

        } else if (input$sourcetool == "vast-tools2"){
          url <- "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA185305/"
          style <-  "style = 'font-size:14px; color: #000000' "
          HTML(paste0("<p",style,",>You are using a subset of the <i>'Deep transcriptional profiling of longitudinal changes during neurogenesis and network
          maturation in vivo'</i> public mouse dataset (<a target='_blank' href='", url, "'>PRJNA185305</a>) as an example, where the inclusion tables have been obtained using vast-tools</p>"))


        } else if (input$sourcetool == "rMATS"){
          url <- "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA185305/"
          style <-  "style = 'font-size:14px; color: #000000' "
          HTML(paste0("<p",style,",>You are using a subset of the <i>'Deep transcriptional profiling of longitudinal changes during neurogenesis and network
          maturation in vivo'</i> public mouse dataset (<a target='_blank' href='", url, "'>PRJNA185305</a>) as an example, where the inclusion tables have been obtained using rMATS.</p>"))

        } else if (input$sourcetool == "whippet"){
          url <- "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA185305/"
          style <-  "style = 'font-size:14px; color: #000000' "
          HTML(paste0("<p",style,",>You are using a subset of the <i>'Deep transcriptional profiling of longitudinal changes during neurogenesis and network
          maturation in vivo'</i> public mouse dataset (<a target='_blank' href='", url, "'>PRJNA185305</a>) as an example, where the inclusion tables have been obtained using whippet. Browse through
          the original filtered tables: </p>"))

        }
      }


    })


    output$files <- renderDT({

      req(sourcetool())

      # Show input file names
      if(!is.null(input$psitable)){

        dt <- data.frame(FileName=input$psitable$name,
                         Size=paste0(input$psitable$size*10^(-6)," MB"))

      } else {
        if (input$sourcetool == "whippet"){

          dt <- data.frame(FileName=paste0("Example_Dataset2_Whippet_",GetTable()$Samples))

        } else if (input$sourcetool %in% c("vast-tools2","rMATS")){

          dt <- data.frame(FileName=paste0("Example_Dataset2_",sourcetool()))

        }else if (input$sourcetool == "vast-tools1"){

          dt <- data.frame(FileName=paste0("Example_Dataset1_",sourcetool()))

        }
      }

      datatable(dt, rownames=FALSE, selection = list(mode = 'single', selected = c(1)), options=list(pageLength = 10, scrollX = TRUE,dom = 't'))

    })



    output$selected_file_table <- renderDT(server = FALSE,{

      req(dataset())
      req(selectAlternatives())
      req(sourcetool())

      selectAlternatives()

      if (sourcetool()=="whippet"){
        # if(!is.null(input$psitable)){
        data <- dataset()[[input$files_rows_selected]]
        events_to_show <- paste0(data$Gene,":",data$Node,":",data$Type) %in% psifiltdataset()$EVENT
        data[events_to_show,]
        # }
      }else if(sourcetool() == "rMATS"){

        eventIDs <- isolate(psifiltdataset())$EVENT
        dataset()[dataset()$ID %in% eventIDs,]

      } else if (sourcetool() == "vast-tools"){

        eventIDs <- isolate(psifiltdataset())$EVENT
        dataset()[dataset()$EVENT %in% eventIDs,]

      }

    }, rownames = FALSE,

    extensions = 'Buttons',

    options = list(paging = TRUE,
                   pageLength = 10,
                   scrollX=TRUE,
                   searching = TRUE,
                   ordering = TRUE,
                   dom = 'Bfrtip',
                   buttons = list(list(extend='copy',
                                       filename = "betAS_SelectedEvents"),
                                  list(extend='csv',
                                       filename = "betAS_SelectedEvents"),
                                  list(extend='excel',
                                       filename= "betAS_SelectedEvents")),
                   lengthMenu=c(3,5,10)))



    # Auxiliary variable to show or not the event types to filter (rMATS tables cannot be filtered by event type, as they only have one)
    output$showchecks <- reactive({

      req(sourcetool())
      req(selectAlternatives())

      if (sourcetool() %in% c("vast-tools","whippet")){

        return(TRUE)

      }else{

        return(FALSE)

      }
    })
    outputOptions(output, 'showchecks', suspendWhenHidden=FALSE)


    # output_textTotalNumberEvents <- reactiveVal("")
    #
    # textTotalNumberEvents <- observeEvent(eventNumber(),{
    #
    #   if (isolate(sourcetool()) %in% c("vast-tools","whippet")){
    #
    #     output_textTotalNumberEvents(paste0("You have selected ", eventNumber(), " events"))
    #
    #   } else if (sourcetool() == "rMATS") {
    #
    #     rMATSEventType <- names(selectAlternatives()$EventsPerType)
    #
    #     if (rMATSEventType == "EX") {
    #       rMATSEventTypeText <- "Exon Skipping"
    #     } else if (rMATSEventType == "IR") {
    #       rMATSEventTypeText <- "Intron Retention"
    #     } else if (rMATSEventType == "Altss") {
    #       rMATSEventTypeText <- "Alternative splice site (Altss)"
    #     } else if (rMATSEventType == "MXE"){
    #       rMATSEventTypeText <- "Mutually Exclusive Exons"
    #     }
    #
    #     output_textTotalNumberEvents(paste0("You have selected ", eventNumber(), " ", rMATSEventTypeText," events"))
    #
    #   }
    #
    # })
    #
    # output$textTotalNumberEvents <- renderText({
    #   output_textTotalNumberEvents()
    # })

    # Step 1: Create a reactive value
    resultText <- reactiveVal("")

    # Step 2: Update the reactive value when eventNumber changes
    observeEvent(eventNumber(), {
      if (sourcetool() %in% c("vast-tools","whippet")){
        resultText(paste0("You have selected ", eventNumber(), " events"))
      }
      else if (sourcetool() == "rMATS") {
        rMATSEventType <- names(isolate(selectAlternatives())$EventsPerType)
        rMATSEventTypeText <- ""
        if (rMATSEventType == "EX") {
          rMATSEventTypeText <- "Exon Skipping"
        } else if (rMATSEventType == "IR") {
          rMATSEventTypeText <- "Intron Retention"
        } else if (rMATSEventType == "Altss") {
          rMATSEventTypeText <- "Alternative splice site (Altss)"
        } else if (rMATSEventType == "MXE"){
          rMATSEventTypeText <- "Mutually Exclusive Exons"
        }

        resultText(paste0("You have selected ", eventNumber(), " ", rMATSEventTypeText," events"))
      }
    })

    # Step 3: Render the reactive value with renderText
    output$textTotalNumberEvents <- renderText({
      resultText()
    })


    output$eventsPiechart <- renderHighchart({

      req(selectAlternatives())

      dataset()

      preparePieForVastToolsCOMPLEX(isolate(psifiltdataset()))

    })



    output$TextToolInfo <- renderText({

      req(sourcetool())

      if(sourcetool() %in% c("rMATS","whippet")){

        paste0("NOTE: ", sourcetool()," tables consider inclusion levels in the [0;1] interval")

      } else if (sourcetool() == "vast-tools"){

        "NOTE: vast-tools tables consider PSI values as percent spliced-in values, i.e., in the [0;100] interval"

      }

    })


    output$plot <- renderPlot({
      req(selectAlternatives())
      dataset()
      bigPicturePlot(isolate(psifiltdataset()))
    })

    # output$downloadPlot <- downloadHandler(
    #   filename = function() { "PSIdistributions.png" },
    #   content = function(file) {
    #     png(file, width = 1500, height = 1000)
    #     print(bigPicturePlot(isolate(psifiltdataset())))
    #     dev.off()
    #   }
    # )


    output$test_box <- renderText({
      validate(need(input$minNreads > 1, "Please use a value greater or equal than 1"))
    })



    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # B. Group Definition ------------------------------------------------------
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


    ## Functions ---------------------------------------------------------------


    selectedeventID <- reactive({
      return(input$indbetaseventid)
    })

    selectedeventIDDiff <- reactive({
      return(input$eventidtoplot)
    })

    selectedeventIDDiffMultiple <- reactive({
      return(input$eventidtoplot_mult)
    })


    selectedEventTable <- reactive({

      input$rundiffbetas

      isolate({

        groupA <- input$groupA
        groupB <- input$groupB

      })

      if(is.null(groupA) || is.null(groupB))return(NULL)

      samplesA    <- values$groups[[groupA]]$samples
      colsGroupA  <- match(samplesA, colnames(isolate(psifiltdataset())))

      samplesB    <- values$groups[[groupB]]$samples
      colsGroupB  <- match(samplesB, colnames(isolate(psifiltdataset())))

      eventList <- prepareTableEvent(eventID = selectedeventIDDiff(),
                                     psitable = psifiltdataset(),
                                     qualtable = qualfiltdataset(),
                                     npoints = 500,
                                     colsA = colsGroupA,
                                     colsB = colsGroupB,
                                     labA = groupA,
                                     labB = groupB,
                                     basalColor = "#89C0AE",
                                     interestColor = "#E69A9C",
                                     maxDevTable = maxDevSimulationN100,
                                     nsim = 1000)
      return(eventList)

    })








    observe({updateSelectInput(inputId = "groupSamples", choices = colnames(psifiltdataset())[-c(1:6)])})

    observe({updateSelectInput(inputId = "indbetaseventid", choices = psifiltdataset()$EVENT)})

    observe({updateSelectInput(inputId = "eventidtoplot_mult", choices = psifiltdataset()$EVENT)})

    values <- reactiveValues(groups = list())

    # plottingEvent <- reactiveValues(tableRow = list())

    # colors <- reactiveValues(used = list())

    observeEvent(input$newGroup, {
      currentNames <- names(values$groups)
      values$groups[[length(values$groups)+1]] <- list(name = input$groupName,
                                                       samples = input$groupSamples,
                                                       color = input$groupColor)
      names(values$groups) <- make.unique(c(currentNames, input$groupName))

      showNotification("A new group has been successfully created",
                       closeButton = TRUE,
                       duration = 5,
                       type = c("default"))

    })





    # nextSuggestedColor <- reactive({
    #   if(length(values$groups)<1){nextColor <- "#89C0AE"}else{
    #
    #     lastColor <- values$groups[[length(values$groups)]]$groupColor
    #     nextColor <- sample(pastelColors[-c(which(pastelColors == lastColor))], size = 1)
    #
    #   }
    #
    #   colors[[length(values$groups)+1]] <- list(color = nextColor)
    #   return(colors)
    #
    # })

    observe({updateSelectInput(inputId = "groupA", choices = names(values$groups))})

    observe({updateSelectInput(inputId = "groupB", choices = names(values$groups))})


    observeEvent(sourcetool(), {

      req(input$sourcetool)

      if(is.null(input$psitable) & input$sourcetool %in% c("rMATS","whippet", "vast-tools2")){

        updateSelectInput(inputId = "groupingFeature", choices = c("CellType","Div"))

      } else if (is.null(input$psitable) & input$sourcetool=="vast-tools1"){

        updateSelectInput(inputId = "groupingFeature", choices = c("organism_part", "developmental_stage", "sex"))


      }

    })

    observeEvent(input$deleteGroups, {

      if(length(values$groups)<1)return(NULL)

      pos <- input$groupsTable_rows_selected

      values$groups[pos] <- NULL

    })

    observeEvent(input$findGroups, {
      req(dataset())
      req(input$findGroups)

      names <- colnames(psifiltdataset())[-c(1:6)]


      not_grouped   <- names
      checked       <- c()
      random_colors <- pastelColors
      groups        <- LETTERS[seq(1, length = min(length(names), 26))]
      used_groups   <- groups

      while((length(checked) < length(names)) & length(used_groups) <= 26){

        groupNames  <- agrep(pattern = not_grouped[1], x = not_grouped, value = TRUE)

        # if(length(groupNames) == 1 | length(groupNames) == length(names)){
        #
        # }else{

        checked     <- c(checked, groupNames)
        not_grouped <- not_grouped[-c(match(groupNames, not_grouped))]

        # Assign new group
        currentNames <- names(values$groups)
        values$groups[[length(values$groups)+1]] <- list(name = used_groups[1],
                                                         samples = groupNames,
                                                         color = random_colors[1])
        names(values$groups) <- make.unique(c(currentNames, used_groups[1]))

        random_colors <- random_colors[-1]
        used_groups   <- used_groups[-1]
        #
        #         }

      }

      # return(NULL)

      showNotification("All found similarities in sample names have been used to define groups.",
                       closeButton = TRUE,
                       duration = 5,
                       type = c("default"))

    })

    observeEvent(input$findGroupsBasedSampleTable, {

      req(dataset())
      req(sampleTable())
      req(sourcetool())
      # req(input$findGroupsBasedSampleTable)
      # req(input$groupingFeatures)

      if(is.null(input$groupingFeature)){

        showNotification("Please select one sample feature to define groups",
                         closeButton = TRUE,
                         duration = 5,
                         type = c("error"))

        return(NULL)

      }

      sampleTable <- as.data.frame(sampleTable())
      columnToGroupBy <- match(input$groupingFeature, colnames(sampleTable))

      groups <- unique(sampleTable[,columnToGroupBy])
      random_colors <- pastelColors

      for(g in groups){

        #groupNames <- sampleTable$Run[which(sampleTable[,columnToGroupBy] == g)]
        if (sourcetool() %in% c("vast-tools","whippet")){
          groupNames <- sampleTable$Run[which(sampleTable[,columnToGroupBy] == g)]
        } else if (sourcetool()=="rMATS"){
          groupNames <- sampleTable$sampleID[which(sampleTable[,columnToGroupBy] == g)]
        }


        # Assign new group
        currentNames <- names(values$groups)
        values$groups[[length(values$groups)+1]] <- list(name = g,
                                                         samples = groupNames,
                                                         color = random_colors[1])
        names(values$groups) <- make.unique(c(currentNames, g))

        random_colors <- random_colors[-1]

      }


      showNotification("Groups based on the feature selected have been successfully created.",
                       closeButton = TRUE,
                       duration = 5,
                       type = c("default"))

    })





    observeEvent(dataset(), {

      req(dataset())
      req(sampleTable())
      req(sourcetool())

      values$groups <- NULL


    })





    ## Outputs -----------------------------------------------------------------




    output$samplesTabletext <- renderText({

      req(sourcetool())

      if(is.null(input$psitable)){


        print("Sample information")

      }

    })



    output$EventInfoWhippet <- renderText({
      req(sourcetool())

      if (sourcetool() == "whippet"){

        "Events should be in the format Gene:Node:Type"

      }

    })


    output$url <- renderUI({
      req(sourcetool())
      if(sourcetool() == "vast-tools"){

        event <- selectedeventID()
        url19 <- paste0("https://vastdb.crg.eu/event/", event, "@hg19")
        url38 <- paste0("https://vastdb.crg.eu/event/", event, "@hg38")

        HTML(paste0("<p>Check event details in VastDB (<a target='_blank' href='", url19, "'>hg19</a> or <a href='", url38, "'>hg38</a>)</p>"))

      }


    })

    output$sampleTable <- renderDT({

      if(is.null(input$psitable)){

        sampleTable()

      }
    }, options = list(pageLength = 10, scrollX = TRUE)
    #},
    #rownames = FALSE,
    #colnames = c("Run", "Age", "DevStage", "Organ", "Sex")
    )


    output$groupsTable <- renderDT({

      if(length(values$groups)<1)return(NULL)

      table <- c()
      for(i in 1:length(values$groups)){

        group <- c("Name" = values$groups[[i]]$name,
                   "Samples" = paste0(values$groups[[i]]$samples, collapse = " "),
                   "Color" = values$groups[[i]]$color)

        table <- rbind(table, group)

      }

      rownames(table) <- NULL

      return(table)

    },
    selection = c("multiple"),
    rownames = FALSE)


    output$densities <- renderPlot({

      if(length(values$groups)<1)return(NULL)

      # plotIndividualDensities(eventID = "HsaEX0055568",
      #                         npoints = 500,
      #                         psitable = psidataset(),
      #                         qualtable = qualdataset(),
      #                         colsA = values$groups[[1]]$samples,
      #                         colsB = values$groups[[2]]$samples,
      #                         labA = values$groups[[1]]$name,
      #                         labB = values$groups[[2]]$name,
      #                         colorA = values$groups[[1]]$color,
      #                         colorB =  values$groups[[2]]$color,
      #                         maxDevTable = maxDevSimulationN100)

      plotIndividualDensitiesList(eventID = selectedeventID(),
                                  npoints = 500,
                                  psitable = psifiltdataset(),
                                  qualtable = qualfiltdataset(),
                                  groupList = values$groups,
                                  maxDevTable = maxDevSimulationN100)

    })


    # auxiliary variable to show (or not) option to group based on variable (currently not supported by betAS for user-input data)

    output$showGroupingFeatureOption <- reactive({

      if (is.null(input$psitable)){

        return(TRUE)

      }else{

        return(FALSE)

      }
    })
    outputOptions(output, 'showGroupingFeatureOption', suspendWhenHidden=FALSE)

    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # C. Differential alternative splicing -------------------------------------
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


    ## Functions ---------------------------------------------------------------


    betasTableVolcano <- reactive({

      req(input$rundiffbetas)

      input$rundiffbetas

      isolate({

        groupA <- input$groupA
        groupB <- input$groupB
        yStat <- input$volcanoYAxis

      })


      if(is.null(groupA) || is.null(groupB) || is.null(yStat))return(NULL)else if(groupA == groupB){

        showNotification("Please select two distict groups for differential splicing analysis",
                         closeButton = TRUE,
                         duration = 5,
                         type = c("error"))
        return(NULL)

      }else if(

        yStat == "Pdiff (probability of differential splicing)"){

        samplesA <- values$groups[[groupA]]$samples
        colsGroupA <- match(samplesA, colnames(isolate(psifiltdataset())))

        samplesB <- values$groups[[groupB]]$samples
        colsGroupB <- match(samplesB, colnames(isolate(psifiltdataset())))

        table <- prepareTableVolcano(psitable = isolate(psifiltdataset()),
                                     qualtable = isolate(qualfiltdataset()),
                                     npoints = 500,
                                     colsA = colsGroupA,
                                     colsB = colsGroupB,
                                     labA = groupA,
                                     labB = groupB,
                                     basalColor = "#89C0AE",
                                     # interestColor = "#EE805B", #orange
                                     interestColor = "#E69A9C", #pink
                                     maxDevTable = maxDevSimulationN100)

      }else if(

        yStat == "F-statistic (median(|between|)/median(|within|))"){

        showNotification("Please note that F-statistic quantifications may take a few minutes.",
                         closeButton = TRUE,
                         duration = 10,
                         type = c("message"))

        samplesA <- values$groups[[groupA]]$samples
        colsGroupA <- match(samplesA, colnames(isolate(psifiltdataset())))

        samplesB <- values$groups[[groupB]]$samples
        colsGroupB <- match(samplesB, colnames(isolate(psifiltdataset())))

        table <- prepareTableVolcanoFstat(psitable = isolate(psifiltdataset()),
                                          qualtable = isolate(qualfiltdataset()),
                                          npoints = 500,
                                          colsA = colsGroupA,
                                          colsB = colsGroupB,
                                          labA = groupA,
                                          labB = groupB,
                                          basalColor = "#89C0AE",
                                          interestColor = "#E69A9C",
                                          maxDevTable = maxDevSimulationN100)

      }else if(

        yStat == "False discovery rate (FDR)"){

        samplesA <- values$groups[[groupA]]$samples
        colsGroupA <- match(samplesA, colnames(isolate(psifiltdataset())))

        samplesB <- values$groups[[groupB]]$samples
        colsGroupB <- match(samplesB, colnames(isolate(psifiltdataset())))

        table <- prepareTableVolcanoFDR(psitable = isolate(psifiltdataset()),
                                        qualtable = isolate(qualfiltdataset()),
                                        npoints = 500,
                                        colsA = colsGroupA,
                                        colsB = colsGroupB,
                                        labA = groupA,
                                        labB = groupB,
                                        basalColor = "#89C0AE",
                                        interestColor = "#E69A9C",
                                        maxDevTable = maxDevSimulationN100,
                                        nsim = 100)

      }

      return(table)

    })

    simplifiedTableVolcano <- reactive({

      req(betasTableVolcano())

      isolate({

        groupA <- input$groupA
        groupB <- input$groupB
        yStat <- input$volcanoYAxis

      })

      if(is.null(groupA) || is.null(groupB) || is.null(yStat))return(NULL)else if(

        yStat == "Pdiff (probability of differential splicing)"){
        column  <- 'Pdiff'

      }else if(

        yStat == "F-statistic (median(|between|)/median(|within|))"){
        column <- 'Fstat'


      }else if(

        yStat == "False discovery rate (FDR)"){
        column <- "invertedFDR"

      }

      simplifiedTable <- betasTableVolcano()
      simplifiedTable <- simplifiedTable[,match(c("EVENT", "GENE", "deltapsi", column), colnames(simplifiedTable))]

      return(simplifiedTable)

    })


    observe({updateSelectInput(inputId = "eventidtoplot", choices = psifiltdataset()$EVENT )})
    #
    # observeEvent(dataset(), {
    #
    #   req(dataset())
    #
    #   updateTextInput(session, "eventidtoplot", "Alternative splicing event to plot:", value = "")
    #
    # })




    ## Outputs -----------------------------------------------------------------


    output$brushed_data <- renderDT(server = FALSE,{

      req(input$plot_brush)
      req(colnameVector())

      datatable(
        brushedPoints(df = simplifiedTableVolcano(), input$plot_brush),
        rownames = FALSE,
        colnames = colnameVector(),
        extensions = 'Buttons',
        options = list(
          dom = 'Bfrtip',
          buttons = list(
            list(extend = 'copy', filename = "betAS_DifferentialAlternativeSplicing"),
            list(extend = 'csv', filename = "betAS_DifferentialAlternativeSplicing"),
            list(extend = 'excel', filename = "betAS_DifferentialAlternativeSplicing")
          ),
          pageLength = 10,
          lengthMenu = c(3, 5, 10),
          paging = TRUE,
          searching = TRUE,
          ordering = TRUE,
          scrollX = TRUE
        )
      ) %>% formatRound(columns = c(3:4), digits = 3)

    })


    output$urlplot <- renderUI({

      req(input$eventidtoplot)
      req(sourcetool())

      if(sourcetool() == "vast-tools"){

        event <- selectedeventIDDiff()

        url19 <- paste0("https://vastdb.crg.eu/event/", event, "@hg19")
        url38 <- paste0("https://vastdb.crg.eu/event/", event, "@hg38")

        HTML(paste0("<p>Check event details in VastDB (<a target='_blank' href='", url19, "'>hg19</a> or <a href='", url38, "'>hg38</a>)</p>"))

      }

    })

    # output$showMeEvent <- renderText({
    #
    # paste0("You have selected ", input$brushed_data_rows_selected$EVENT)
    #
    # })

    output$densitiesSelectedEvent <- renderPlot({

      input$rundiffbetas

      req(betasTableVolcano())
      req(selectedeventIDDiff())
      req(input$plotEvent)

      isolate({

        groupA <- input$groupA
        groupB <- input$groupB
        yStat <- input$volcanoYAxis

      })

      if(is.null(groupA) || is.null(groupB) || is.null(yStat))return(NULL)

      samplesA <- values$groups[[groupA]]$samples
      colsGroupA <- match(samplesA, colnames(isolate(psifiltdataset())))

      samplesB <- values$groups[[groupB]]$samples
      colsGroupB <- match(samplesB, colnames(isolate(psifiltdataset())))

      # plotIndividualDensities(eventID = "HsaEX0055568",
      #                         npoints = 500,
      #                         psitable = psifiltdataset(),
      #                         qualtable = qualfiltdataset(),
      #                         colsA = values$groups[[1]]$samples,
      #                         colsB = values$groups[[2]]$samples,
      #                         labA = values$groups[[1]]$name,
      #                         labB = values$groups[[2]]$name,
      #                         colorA = values$groups[[1]]$color,
      #                         colorB =  values$groups[[2]]$color,
      #                         maxDevTable = maxDevSimulationN100)

      plotIndividualDensitiesList(eventID = selectedeventIDDiff(),
                                  npoints = 500,
                                  psitable = psifiltdataset(),
                                  qualtable = qualfiltdataset(),
                                  groupList = values$groups[c(groupA, groupB)],
                                  maxDevTable = maxDevSimulationN100)

    })



    output$volcano <- renderPlot({

      input$rundiffbetas

      req(betasTableVolcano())
      req(simplifiedTableVolcano())

      isolate({

        groupA <- input$groupA
        groupB <- input$groupB
        yStat <- input$volcanoYAxis

      })


      if(is.null(groupA) || is.null(groupB) || is.null(yStat))return(NULL)else if(

        yStat == "Pdiff (probability of differential splicing)"){

        samplesA <- values$groups[[groupA]]$samples
        colsGroupA <- match(samplesA, colnames(isolate(psifiltdataset())))

        samplesB <- values$groups[[groupB]]$samples
        colsGroupB <- match(samplesB, colnames(isolate(psifiltdataset())))

        plotVolcano(betasTable = isolate(betasTableVolcano()),
                    labA = groupA,
                    labB = groupB,
                    basalColor = "#89C0AE",
                    # interestColor = "#EE805B", #orange
                    interestColor = "#E69A9C")

      }else if(

        yStat == "F-statistic (median(|between|)/median(|within|))"){

        samplesA <- values$groups[[groupA]]$samples
        colsGroupA <- match(samplesA, colnames(isolate(psifiltdataset())))

        samplesB <- values$groups[[groupB]]$samples
        colsGroupB <- match(samplesB, colnames(isolate(psifiltdataset())))

        plotVolcanoFstat(betasTable = isolate(betasTableVolcano()),
                         labA = groupA,
                         labB = groupB,
                         basalColor = "#89C0AE",
                         # interestColor = "#EE805B", #orange
                         interestColor = "#E69A9C")

      }else if(

        yStat == "False discovery rate (FDR)"){

        samplesA <- values$groups[[groupA]]$samples
        colsGroupA <- match(samplesA, colnames(isolate(psifiltdataset())))

        samplesB <- values$groups[[groupB]]$samples
        colsGroupB <- match(samplesB, colnames(isolate(psifiltdataset())))

        plotVolcanoFDR(betasTable = isolate(betasTableVolcano()),
                       labA = groupA,
                       labB = groupB,
                       basalColor = "#89C0AE",
                       # interestColor = "#EE805B", #orange
                       interestColor = "#E69A9C")

      }

    })


    output$PDiffEventPlot <- renderPlot({

      input$rundiffbetas

      req(betasTableVolcano())
      req(selectedeventIDDiff())
      req(selectedEventTable())
      req(input$plotEvent)


      isolate({

        groupA <- input$groupA
        groupB <- input$groupB

      })

      if(is.null(groupA) || is.null(groupB))return(NULL)

      samplesA <- values$groups[[groupA]]$samples
      colsGroupA <- match(samplesA, colnames(isolate(psifiltdataset())))

      samplesB <- values$groups[[groupB]]$samples
      colsGroupB <- match(samplesB, colnames(isolate(psifiltdataset())))

      plotPDiffFromEventObjList(eventObjList = isolate(selectedEventTable()))

    })

    output$FstatEventPlot <- renderPlot({

      input$rundiffbetas

      req(betasTableVolcano())
      req(selectedeventIDDiff())
      req(selectedEventTable())
      req(input$plotEvent)


      isolate({

        groupA <- input$groupA
        groupB <- input$groupB

      })

      if(is.null(groupA) || is.null(groupB))return(NULL)

      samplesA <- values$groups[[groupA]]$samples
      colsGroupA <- match(samplesA, colnames(isolate(psifiltdataset())))

      samplesB <- values$groups[[groupB]]$samples
      colsGroupB <- match(samplesB, colnames(isolate(psifiltdataset())))

      plotFstatFromEventObjList(eventObjList = isolate(selectedEventTable()))

    })

    output$FDREventPlot <- renderPlot({

      input$rundiffbetas

      req(betasTableVolcano())
      req(selectedeventIDDiff())
      req(selectedEventTable())
      req(input$plotEvent)


      isolate({

        groupA <- input$groupA
        groupB <- input$groupB

      })

      if(is.null(groupA) || is.null(groupB))return(NULL)

      samplesA <- values$groups[[groupA]]$samples
      colsGroupA <- match(samplesA, colnames(isolate(psifiltdataset())))

      samplesB <- values$groups[[groupB]]$samples
      colsGroupB <- match(samplesB, colnames(isolate(psifiltdataset())))

      plotFDRFromEventObjList(eventObjList = isolate(selectedEventTable()))

    })


    colnameVector <- reactive({

      req(input$plot_brush)
      # req(input$volcanoYAxis())
      req(betasTableVolcano())
      req(simplifiedTableVolcano())

      isolate({

        groupA <- input$groupA
        groupB <- input$groupB
        yStat   <- input$volcanoYAxis

      })

      if(yStat == "Pdiff (probability of differential splicing)"){

        column <- "Pdiff"
        colnameVector <- c("Event ID" = "EVENT", "Gene" = "GENE", "dPSI" = "deltapsi", "Pdiff" = "Pdiff")

      }else if(

        yStat == "F-statistic (median(|between|)/median(|within|))"){

        column <- "Fstat"
        colnameVector <- c("Event ID" = "EVENT", "Gene" = "GENE", "dPSI" = "deltapsi", "F-stat" = "Fstat")

      }else if(

        yStat == "False discovery rate (FDR)"){

        column <- "invertedFDR"
        colnameVector <- c("Event ID" = "EVENT", "Gene" = "GENE", "dPSI" = "deltapsi", "1-FDR" = "invertedFDR")

      }

      colnameVector

    })





    output$densitiesTest <- renderPlot({

      req(selectedEventTable())

      isolate({

        groupA <- input$groupA
        groupB <- input$groupB

      })

      if(is.null(groupA) || is.null(groupB))return(NULL)

      samplesA <- values$groups[[groupA]]$samples
      colsGroupA <- match(samplesA, colnames(isolate(psifiltdataset())))

      samplesB <- values$groups[[groupB]]$samples
      colsGroupB <- match(samplesB, colnames(isolate(psifiltdataset())))


      plotDensitiesFromEventObjList(eventObjList = isolate(selectedEventTable()))

    })



    # Auxiliary variable to show individual event plots only if a given event is selected
    output$showIndividualEventPlots <- reactive({

      req(betasTableVolcano())
      req(selectedeventIDDiff())
      req(selectedEventTable())
      req(input$plotEvent)

      if (!is.null(input$plotEvent)){

        return(TRUE)

      }else{

        return(FALSE)

      }
    })
    outputOptions(output, 'showIndividualEventPlots', suspendWhenHidden=FALSE)




    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # D. Differential alternative splicing (multiple groups) -------------------
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


    ## Functions ---------------------------------------------------------------


    betasTableVolcanoMultiple <- reactive({

      req(input$rundiffbetas_mult)

      input$rundiffbetas_mult

      isolate({

        yStat <- input$volcanoYAxis_mult

      })

      if(length(values$groups)<1){

        showNotification("Please define groups for differential splicing analysis (see 'Group definition' tab)",
                         closeButton = TRUE,
                         duration = 5,
                         type = c("error"))
        return(NULL)

      }

      showNotification("Please note that multiple-group quantifications may take some minutes.",
                       closeButton = TRUE,
                       duration = 10,
                       type = c("message"))

      table <- prepareTableVolcanoMultipleGroups(psitable = isolate(psifiltdataset()),
                                                 qualtable = isolate(qualfiltdataset()),
                                                 groupList = values$groups,
                                                 npoints = 500,
                                                 maxDevTable = maxDevSimulationN100)

      return(table)

    })

    simplifiedTableVolcanoMultiple <- reactive({

      req(betasTableVolcanoMultiple())

      isolate({

        yStat <- input$volcanoYAxis_mult

      })

      if(length(values$groups)<1){

        showNotification("Please define groups for differential splicing analysis (see 'Group definition' tab)",
                         closeButton = TRUE,
                         duration = 5,
                         type = c("error"))
        return(NULL)

      }else if(

        yStat == "Pdiff (probability that |between| > |within|)"){
        column  <- 'Pdiff'

      }else if(

        yStat == "F-statistic (median(|between|)/median(|within|))"){
        column <- 'Fstat'


      }

      simplifiedTable <- betasTableVolcanoMultiple()
      simplifiedTable <- simplifiedTable[,match(c("EVENT", "GENE", "deltaAbsolute", column), colnames(simplifiedTable))]

      return(simplifiedTable)

    })

    selectedEventTableMultiple <- reactive({

      input$rundiffbetas_mult

      eventList <- prepareTableEventMultiple(eventID = selectedeventIDDiffMultiple(),
                                             psitable = psifiltdataset(),
                                             qualtable = qualfiltdataset(),
                                             groupList = values$groups,
                                             npoints = 500,
                                             maxDevTable = maxDevSimulationN100)

      return(eventList)

    })



    colnameVector_mult <- reactive({

      req(input$plot_brush_mult)
      # req(input$volcanoYAxis())
      req(betasTableVolcanoMultiple())
      req(simplifiedTableVolcanoMultiple())

      isolate({

        yStat   <- input$volcanoYAxis_mult

      })

      if(yStat == "Pdiff (probability that |between| > |within|)"){

        column <- "Pdiff"
        colnameVector <- c("Event ID" = "EVENT", "Gene" = "GENE", "dPSI(between to within)" = "deltaAbsolute", "Pdiff" = "Pdiff")

      }else if(

        yStat == "F-statistic (median(|between|)/median(|within|))"){

        column <- "Fstat"
        colnameVector <- c("Event ID" = "EVENT", "Gene" = "GENE", "dPSI(between to within)" = "deltaAbsolute", "F-stat" = "Fstat")

      }

      colnameVector

    })


    ## Outputs -----------------------------------------------------------------



    output$groupsTableMultiple <- renderDT({

      if(length(values$groups)<1)return(NULL)

      table <- c()
      for(i in 1:length(values$groups)){

        group <- c("Name" = values$groups[[i]]$name,
                   "Samples" = paste0(values$groups[[i]]$samples, collapse = " "),
                   "Color" = values$groups[[i]]$color)

        table <- rbind(table, group)

      }

      rownames(table) <- NULL

      return(table)

    },
    selection = c("multiple"),
    rownames = FALSE)




    # Outputs from "multiple group" tab

    output$volcano_mult <- renderPlot({

      input$rundiffbetas_mult

      req(betasTableVolcanoMultiple())
      req(simplifiedTableVolcanoMultiple())

      isolate({

        yStat <- input$volcanoYAxis_mult

      })

      if(length(values$groups)<1){

        showNotification("Please define groups for differential splicing analysis (see 'Group definition' tab)",
                         closeButton = TRUE,
                         duration = 5,
                         type = c("error"))
        return(NULL)

      }else if(

        yStat == "Pdiff (probability that |between| > |within|)"){

        plotVolcano_MultipleGroups_Pdiff(betasTable = isolate(betasTableVolcanoMultiple()))

      }else if(

        yStat == "F-statistic (median(|between|)/median(|within|))"){

        plotVolcano_MultipleGroups_Fstat(betasTable = isolate(betasTableVolcanoMultiple()))

      }

    })

    output$brushed_data_mult <- renderDT(server = FALSE,{

      req(input$plot_brush_mult)
      req(colnameVector_mult())

      datatable(brushedPoints(df = simplifiedTableVolcanoMultiple(), input$plot_brush_mult),
                rownames = FALSE,
                colnames = colnameVector_mult(),
                extensions = 'Buttons',
                options = list(
                  dom = 'Bfrtip',
                  buttons = list(
                    list(extend = 'copy', filename = "betAS_DifferentialAlternativeSplicing_MultGroups"),
                    list(extend = 'csv', filename = "betAS_DifferentialAlternativeSplicing_MultGroups"),
                    list(extend = 'excel', filename = "betAS_DifferentialAlternativeSplicing_MultGroups")
                  ),
                  pageLength = 10,
                  lengthMenu = c(3, 5, 10),
                  paging = TRUE,
                  searching = TRUE,
                  ordering = TRUE,
                  scrollX = TRUE
                )) %>% formatRound(columns = c(3:4), digits = 3)

    })
    # Required if using renderDT from shiny but with DT (to allow formatRound()); rownames/colnames are replaced by the parameters inside DT::datatable
    # rownames = FALSE,
    # colnames = colnameVector())

    output$violins_multgroup <- renderPlot({

      req(input$plotEvent_mult)

      # req(input$eventidtoplot_mult)
      # req(input$plotEvent)

      if(length(values$groups)<1)return(NULL)

      plotIndividualViolinsList(eventID = selectedeventIDDiffMultiple(),
                                npoints = 500,
                                psitable = psifiltdataset(),
                                qualtable = qualfiltdataset(),
                                groupList = values$groups,
                                maxDevTable = maxDevSimulationN100)

    })

    output$FstatEventPlot_multgroup <- renderPlot({

      req(input$plotEvent_mult)
      # req(input$eventidtoplot_mult)
      # req(input$plotEvent_mult)

      if(length(values$groups)<1)return(NULL)

      plotFstatFromEventObjListMultiple(eventObjList = isolate(selectedEventTableMultiple()))

    })

    output$urlplot_mult <- renderUI({

      req(input$eventidtoplot_mult)
      req(sourcetool())

      if(sourcetool() == "vast-tools"){

        event <- selectedeventIDDiffMultiple()
        url19 <- paste0("https://vastdb.crg.eu/event/", event, "@hg19")
        url38 <- paste0("https://vastdb.crg.eu/event/", event, "@hg38")

        HTML(paste0("<p>Check event details in VastDB (<a target='_blank' href='", url19, "'>hg19</a> or <a href='", url38, "'>hg38</a>)</p>"))

      }
    })

    # Auxiliary variable to show individual event plots only if a given event is selected
    output$showIndividualEventPlotsmult <- reactive({

      req(betasTableVolcanoMultiple())
      req(selectedeventIDDiffMultiple())
      req(selectedEventTableMultiple())
      req(input$plotEvent_mult)

      if (!is.null(input$plotEvent_mult)){

        return(TRUE)

      }else{

        return(FALSE)

      }
    })
    outputOptions(output, 'showIndividualEventPlotsmult', suspendWhenHidden=FALSE)





  }

  return(server)

}

#' Run betAS graphical interface
#'
#' @param ... Arguments passed to shiny::runApp()
#'
#' @return betAS graphical interface (shiny app)
#' @importFrom shiny shinyApp runApp
#' @export
betASapp <- function(...){

  options(shiny.maxRequestSize = 100*1024^2, shiny.error=NULL)
  #options(shiny.maxRequestSize = 100*1024^2,  shiny.error = browser)
  # construct and start a Shiny application from UI and server
  app <- shinyApp(betASapp_ui(), betASapp_server())

  runApp(app, ...)

}


