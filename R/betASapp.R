#' @import shiny
#' @import highcharter
#' @import shinycssloaders
#' @importFrom colourpicker colourInput updateColourInput
#' @importFrom DT renderDT DTOutput formatRound datatable
#' @importFrom bslib bs_theme
betASapp_ui <- function(){
  # :::: Variables ::::
  # tools           <- c("vast-tools", "MISO", "SUPPA", "Other")
  availabletools      <- c("vast-tools", "rMATS")
  yAxisStats          <- c("Pdiff (probability of differential splicing)", "F-statistic (median(|between|)/median(|within|))", "False discovery rate (FDR)")
  yAxisStats_multiple <- c("Pdiff (probability that |between| > |within|)", "F-statistic (median(|between|)/median(|within|))")
  eventTypes          <- c("Exon skipping (ES)"="EX", "Intron retention (IR)"="IR", "Alternative splice site (Altss)"="Altss")
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

    # checkboxInput("dark_mode", "Dark mode", FALSE),

    titlePanel("betAS: intuitive visualisation of differential alternative splicing"),

    h5(strong(div("(development version)", style = "color: #E0E0E0"))),

    tabsetPanel(

      tabPanel("Import inclusion levels",  icon = icon("play"),

               sidebarLayout(fluid = TRUE,

                 sidebarPanel(style = "max-height: 100%",
                              # h3("Exploratory analysis of inclusion levels"),
                              h4("Exploratory analysis of inclusion levels"),


                              "Upload table with inclusion level quantification (e.g. PSI)",

                              tags$ul(
                                tags$li("each row is an alternative splicing event"),
                                tags$li("each column is a sample")
                              ),

                              fileInput("psitable", NULL),

                              HTML(paste0("<p>Alternatively, explore a subset of the publicly available dataset:",
                                          "<br>",
                                          "<a href='", "https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6814/",
                                          "'>Human RNA-seq time-series of the development of seven major organs</a></p>")),

                              helpText("(betAS currently supports inclusion level tables from: vast-tools and rMATS (*.MATS.JC.txt tables))"),
                              radioButtons("sourcetool", label = "Table source:", choices = availabletools),

                              h5("Filter events from loaded table:"),
                              checkboxGroupInput("types", label = "Event types to consider:", selected = c("EX", "IR"),  choices = eventTypes),
                              sliderInput("psirange", "PSI values to consider:", value = c(1, 99), min = 0, max = 100),
                              helpText((em("Consider only alternative splicing events with all PSI values within this range.")),
                              h6(textOutput("textTotalNumberEvents")),
                              # textOutput("textNumberEventsPerType"),
                              # actionButton("filter", "Filter original table", icon = icon("filter"), class = "btn-success"),
                              highchartOutput("eventsPiechart")

                 )),

                 mainPanel(
                   shinycssloaders::withSpinner(plotOutput("plot", height = "1200px"), type = 8, color = "#FF9AA2", size = 2),
                   # shinycssloaders::withSpinner(plotOutput("plot"), type = 8, color = "#FF9AA2", size = 2),
                   h6("NOTE: vast-tools tables consider PSI values as percent spliced-in values, i.e., in the [0;100] interval"),
                   DTOutput("table")
                 )
               )
      ),

      tabPanel("Group definition", icon = icon("gears"),

             sidebarLayout(fluid = TRUE,

               sidebarPanel(
                            selectInput("indbetaseventid", "Select alternative splicing event to plot:", choices = NULL),
                            uiOutput("url"),

                            h5("Sample information"),
                            DTOutput("sampleTable"),

                            h4("Organise samples into groups/conditions"),

                            helpText(tags$ul(
                              tags$li("one beta distribution will be estimated per sample"),
                              tags$li("groups defined can be used for differential splicing analyses")
                            )),


                            # h5("Group creation"),
                            #
                            # helpText(em("Choose one option for group creation and check the resulting grouping below")),

                            tags$ol(

                              tags$li(h5("Add/delete groups by sample selection:")),

                              textInput(inputId = "groupName", label = "Insert group name:", placeholder = "e.g. Control"),
                              selectInput("groupSamples", "Choose samples in group:", choices = NULL, multiple = TRUE),
                              colourInput("groupColor", "Select group color:", value = "#89C0AE"),
                              actionButton("newGroup", "Create new group", icon = icon("layer-group"), class = "btn-secondary"),

                              tags$li(h5("Automatic group creation:")),

                              tags$ol(

                                tags$li(h6("Group samples based on a given feature")),
                                selectInput("groupingFeature", label = "Select feature to group samples by:", choices = c("organism_part", "developmental_stage", "sex")),
                                actionButton("findGroupsBasedSampleTable", label = "Feature-associated group(s)", icon = icon("robot"), class = "btn-info"),

                                tags$li(h6("Group samples based on sample name similarities")),
                                actionButton("findGroups", "Automatic group(s)", icon = icon("wand-sparkles"), class = "btn-info")

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
                            # "Perform differential alternative splicing analysis for all events",
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
                 textInput("eventidtoplot", "Alternative splicing event to plot:", placeholder = "e.g. HsaEX0007927"),
                 uiOutput("urlplot"),
                 actionButton("plotEvent", "Plot considered event", icon = icon("mound"), class = "btn-secondary"),
                 # textOutput("showMeEvent"),
                 # plotOutput("densitiesSelectedEvent", height = "400px", width = "500px")

               ),

               mainPanel(

                   shinycssloaders::withSpinner(plotOutput("volcano",
                                                           height = "900px",
                                                           brush = brushOpts(id = "plot_brush")),
                                                type = 8,
                                                color = "#FF9AA2",
                                                size = 2),

                   fluidRow(
                     plotOutput("densitiesSelectedEvent", height = "400px")
                   ),

                   fluidRow(

                     column(4,
                            h4("P(A - B) > x:"),
                            h6("Probabilty that estimated PSI of one group is greater than the other by x.")),

                     column(4,
                            h4("F-statistic:"),
                            h6("Ratio between absolute differences 'between' and 'within' groups.")),

                     column(4,
                            h4("False discovery rate (FDR):"),
                            h6("Probability of getting an absolute difference in PSI between groups greater than the observed, under the null hypothesis that all PSIs come from the same distribution."))

                   ),

                   fluidRow(
                     column(4,

                            plotOutput("PDiffEventPlot", height = "400px", width = "400px")),

                     column(4,

                            plotOutput("FstatEventPlot", height = "400px", width = "400px")),

                     column(4,

                            plotOutput("FDREventPlot", height = "400px", width = "400px"))

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

                             fluidRow(

                               column(3,
                                      h5("Absolute differences 'between' and 'within' groups.")),

                               column(9,
                                      h5("PSI quantifications per sample across groups."))),

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

  options(highcharter.theme = hc_theme_smpl(tooltip = list(valueDecimals = 2)))

  # :::: Variables ::::
  # file                  <- "test/INCLUSION_LEVELS_FULL-Hsa32-hg19_to_test.tab.gz"
  # testTable             <- read.table(gzfile(file), sep="\t", header=TRUE, quote="")
  # maxDevSimulationN100  <- readRDS(url("http://imm.medicina.ulisboa.pt/group/distrans/SharedFiles/Mariana/Splicing&SenescenceFLEX/xintercepts_100incr_100cov_100trials.R"))
  maxDevSimulationN100  <- readRDS("test/xintercepts_100incr_100cov_100trials.rds")
  pastelColors          <- c("#FF9AA2", "#FFB7B2", "#FFDAC1", "#E2F0CB", "#B5EAD7", "#C7CEEA", "#FBE2FD", "#D9ECFE")

  # simplify test table
  # colnames(testTable)   <- gsub(x = colnames(testTable), pattern = "Sample_IMR90_", replacement = "")
  # samples               <- colnames(testTable)[grep(x = colnames(testTable), pattern = "-Q")-1]

  # specifies the behaviour of the app by defining a server function
  server <- function(input, output, session){
    thematic_shiny()

    # not working, but based on slide 25 from https://talks.cpsievert.me/20201014/#25
    # observe(session$setCurrentTheme(
    #   if (input$dark_mode) dark else light
    # ))

    # Loaded main original table
    dataset <- reactive({

      if(is.null(input$psitable)){

        testTable <- readRDS(file = "test/INCLUSION_LEVELS_FULL-hg19-98-v251.rds")

        if(input$sourcetool == "rMATS"){

          testTable <- read.delim(file = "test/SE.MATS.JC.txt")

        }

        return(testTable)

      }

      req(input$psitable)

      if(length(grep(pattern = "[.]gz", x = input$psitable$datapath)) == 0){

        loadingFile <- input$psitable$datapath

      }

      loadingFile <- gzfile(input$psitable$datapath)

      read.table(loadingFile, sep="\t", header=TRUE, quote="")

      })

    sampleTable <- reactive({

      req(dataset())

      if(is.null(input$psitable)){

        samplesTable <- readRDS(file = "test/samplesTable.rds")

        return(samplesTable)

      }

    })

    # 1. Filter table to remove events: with at least one NA, ANN and with at least one sample with less than VLOW quality
    filterVastToolsTable <- reactive({

      selectedEventTypes <- c()

      if("EX" %in% input$types) selectedEventTypes <- c("C1", "C2", "C3", "S", "MIC")
      if("IR" %in% input$types) selectedEventTypes <- c(selectedEventTypes, "IR-C", "IR-S")
      if("Altss" %in% input$types) selectedEventTypes <- c(selectedEventTypes, "Alt3", "Alt5")
      if(length(selectedEventTypes) == 0){
        showNotification("Please select at least one event type", duration = 5, type = c("error"))
        return(NULL)
      }

      filteredList <- filterVastTools(dataset(), types = selectedEventTypes)
      return(filteredList)

    })

    filterRMATSTable <- reactive({

      if(input$sourcetool == "rMATS"){

        filteredList <- filterrMATS(dataset())
        return(filteredList)

      }

    })

    selectAlternatives <- reactive({

      alternativeList <- alternativeVastTools(req(filterVastToolsTable()), minPsi = input$psirange[1], maxPsi = input$psirange[2])

      if(nrow(alternativeList$PSI) == 0){

        showNotification("There are no events with PSI values within such range.",
                         closeButton = TRUE,
                         duration = 5,
                         type = c("error"))
        return(NULL)

      }

      return(alternativeList)

    })

    selectAlternativesRM <- reactive({

      alternativeList <- alternativerMATS(req(filterRMATSTable()), minPsi = input$psirange[1], maxPsi = input$psirange[2])

      if(nrow(alternativeList$PSI) == 0){

        showNotification("There are no events with PSI values within such range.",
                         closeButton = TRUE,
                         duration = 5,
                         type = c("error"))
        return(NULL)

      }

      return(alternativeList)

    })

    # create a reactive expression
    psidataset <- reactive({

      if(input$sourcetool == "vast-tools"){

        return(filterVastToolsTable()$PSI)

      }

      if(input$sourcetool == "rMATS"){

        return(filterRMATSTable()$PSI)

      }

    })

    qualdataset <- reactive({

      if(input$sourcetool == "vast-tools"){

        return(filterVastToolsTable()$Qual)

      }

      if(input$sourcetool == "rMATS"){

        return(filterRMATSTable()$Qual)

      }

    })

    # create a reactive expression
    psifiltdataset <- reactive({

      if(input$sourcetool == "vast-tools"){

        req(selectAlternatives())
        return(selectAlternatives()$PSI)

      }

      if(input$sourcetool == "rMATS"){

        req(selectAlternativesRM())
        return(selectAlternativesRM()$PSI)

      }


    })

    qualfiltdataset <- reactive({

      if(input$sourcetool == "vast-tools"){

        req(selectAlternatives())
        return(selectAlternatives()$Qual)

      }

      if(input$sourcetool == "rMATS"){

        req(selectAlternativesRM())
        return(selectAlternativesRM()$Qual)

      }


    })

    eventNumber <- reactive({
      return(nrow(psifiltdataset()))
    })

    output$textTotalNumberEvents <- renderText({

      # req(selectAlternatives())

      paste0("You have selected ", eventNumber(), " events")

    })

    eventNumberPerType <- reactive({

      selTypes  <- input$types
      toPrint   <- selTypes
      printed   <- c()
      message   <- character()
      while(length(printed) < length(selTypes)){
        type  <- toPrint[1]

        if(type == "EX") selectedEventTypes <- c("C1", "C2", "C3", "S", "MIC")
        if(type == "IR") selectedEventTypes <- c("IR-C", "IR-S")
        if(type == "Altss") selectedEventTypes <- c("Alt3", "Alt5")

        count <- length(which(psifiltdataset()$COMPLEX %in% selectedEventTypes))

        message <- paste0(message, type, ": ", count, " ")
        printed <- c(printed, type)
        toPrint <- toPrint[-c(1)]
      }

      return(message)

    })

    output$textNumberEventsPerType <- renderText({

      req(selectAlternatives())

      eventNumberPerType()

    })

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

    # observe({updateColourInput(session,
    #                            inputID = "groupColor",
    #                            label = "Select group color:",
    #                            value = nextSuggestedColor())
    # })

    observeEvent(input$psitable, {

      req(dataset())
      req(sampleTable())
      # req(input$findGroupsBasedSampleTable)
      # req(input$groupingFeatures)

      if(is.null(input$psitable)){

        updateSelectInput(inputId = "groupingFeature", choices = NULL)

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

        if(length(groupNames) == 1 | length(groupNames) == length(names)){

        }else{

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

        }

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
      # req(input$findGroupsBasedSampleTable)
      # req(input$groupingFeatures)

      if(is.null(input$groupingFeature)){

        showNotification("Please select one sample feature to define groups",
                         closeButton = TRUE,
                         duration = 5,
                         type = c("error"))

        return(NULL)

      }

      sampleTable <- sampleTable()
      columnToGroupBy <- match(input$groupingFeature, colnames(sampleTable))

      groups <- unique(sampleTable[,columnToGroupBy])
      random_colors <- pastelColors

      for(g in groups){

        groupNames <- sampleTable$Run[which(sampleTable[,columnToGroupBy] == g)]

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

    output$sampleTable <- renderDT({

      if(is.null(input$psitable)){

       sampleTable()

      }

    },
    rownames = FALSE,
    colnames = c("Run", "Age", "DevStage", "Organ", "Sex"))

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

    output$url <- renderUI({
      event <- selectedeventID()
      url <- paste0("https://vastdb.crg.eu/event/", event, "@hg19")
      HTML(paste0("<p>Check event details in <a href='", url, "'>VastDB (hg19)</a></p>"))
    })

    output$eventsPiechart <- renderHighchart({

      req(selectAlternatives())

      preparePieForVastToolsCOMPLEX(psifiltdataset())

    })

    output$plot <- renderPlot({

      req(selectAlternatives())

      bigPicturePlot(psifiltdataset())

    })

    output$table <- renderDT({

      req(selectAlternatives())

      psifiltdataset()

    }, rownames = FALSE)

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

    output$brushed_data <- renderDT({

      req(input$plot_brush)
      req(colnameVector())

      # brushedPoints(df = simplifiedTableVolcano(), input$plot_brush, xvar = "deltapsi", yvar = column)
      datatable(brushedPoints(df = simplifiedTableVolcano(), input$plot_brush),
                    rownames = FALSE,
                    colnames = colnameVector()) %>% formatRound(columns = c(3:4), digits = 3)

    })
    # Required if using renderDT from shiny but with DT (to allow formatRound()); rownames/colnames are replaced by the parameters inside DT::datatable
    # rownames = FALSE,
    # colnames = colnameVector())

    output$urlplot <- renderUI({
      req(input$eventidtoplot)
      event <- selectedeventIDDiff()
      url <- paste0("https://vastdb.crg.eu/event/", event, "@hg19")
      HTML(paste0("<p>Check event details in <a href='", url, "'>VastDB (hg19)</a></p>"))
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

    output$brushed_data_mult <- renderDT({

      req(input$plot_brush_mult)
      req(colnameVector_mult())

      datatable(brushedPoints(df = simplifiedTableVolcanoMultiple(), input$plot_brush_mult),
                rownames = FALSE,
                colnames = colnameVector_mult()) %>% formatRound(columns = c(3:4), digits = 3)

    })
    # Required if using renderDT from shiny but with DT (to allow formatRound()); rownames/colnames are replaced by the parameters inside DT::datatable
    # rownames = FALSE,
    # colnames = colnameVector())

    output$violins_multgroup <- renderPlot({

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

      # req(input$eventidtoplot_mult)
      # req(input$plotEvent_mult)

      if(length(values$groups)<1)return(NULL)

      plotFstatFromEventObjListMultiple(eventObjList = isolate(selectedEventTableMultiple()))

    })

    output$urlplot_mult <- renderUI({
      req(input$eventidtoplot_mult)
      event <- selectedeventIDDiffMultiple()
      url <- paste0("https://vastdb.crg.eu/event/", event, "@hg19")
      HTML(paste0("<p>Check event details in <a href='", url, "'>VastDB (hg19)</a></p>"))
    })






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

  options(shiny.maxRequestSize = 200*1024^2)

  # construct and start a Shiny application from UI and server
  app <- shinyApp(betASapp_ui(), betASapp_server())

  runApp(app, ...)

}


