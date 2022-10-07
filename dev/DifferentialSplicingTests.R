
library(betAS)

# :::::::::::::::::::::::::::::::::::::::::::::::::::::
# 1. Simulate loading a (vast-tools) table into the app
# :::::::::::::::::::::::::::::::::::::::::::::::::::::

testTable   <- readRDS("test/INCLUSION_LEVELS_FULL-hg19-98-v251.rds")
sampleTable <- readRDS("test/samplesTable.rds")
# maxDevSimulationN100 <- readRDS(url("http://imm.medicina.ulisboa.pt/group/distrans/SharedFiles/Mariana/Splicing&SenescenceFLEX/xintercepts_100incr_100cov_100trials.R"))
maxDevSimulationN100  <- readRDS("test/xintercepts_100incr_100cov_100trials.rds")

cat(paste0("Initial events: ", nrow(testTable)))
table(testTable$COMPLEX)

pastelColors  <- c("#FF9AA2", "#FFB7B2", "#FFDAC1", "#E2F0CB", "#B5EAD7", "#C7CEEA", "#FBE2FD", "#D9ECFE")

# :::::::::::::::::::::::::::::::::::::::::::::::::::::
# 2. Select events to consider - filtered events (no NAs, minimal coverage, event type)
# :::::::::::::::::::::::::::::::::::::::::::::::::::::
testTable   <- filterVastTools(testTable, types = c("C1", "C2", "C3", "S", "MIC"))
cat(paste0("Filtered events: ", nrow(testTable$PSI)))
table(testTable$PSI$COMPLEX)

# :::::::::::::::::::::::::::::::::::::::::::::::::::::
# 3. Alternative events (PSI within selected values for all samples)
# :::::::::::::::::::::::::::::::::::::::::::::::::::::
testTable <- alternativeVastTools(testTable, minPsi = 1, maxPsi = 99)
cat(paste0("Alternative events: ", nrow(testTable$PSI)))

# :::::::::::::::::::::::::::::::::::::::::::::::::::::
# PSI and Qual table to be used in next sections
# :::::::::::::::::::::::::::::::::::::::::::::::::::::

psiObject <- testTable
psiTable  <- psiObject$PSI
qualTable <- psiObject$Qual

# ::: Plot pie chart
pieChart <- preparePieForVastToolsCOMPLEX(table = psiTable)

# ::: "Big picture" plot
bigPicturePlot <- bigPicturePlot(table = psiTable)

# :::::::::::::::::::::::::::::::::::::::::::::::::::::
# 4. Define groups based on sampleTable
# :::::::::::::::::::::::::::::::::::::::::::::::::::::

groups <- unique(sampleTable[,"organism_part"])
random_colors <- pastelColors
groupList <- list()

for(i in 1:length(groups$organism_part)){

  groupNames <- sampleTable$Run[which(sampleTable[,"organism_part"] == groups$organism_part[i])]

  # Assign new group
  currentNames <- names(groupList)
  groupList[[length(groupList)+1]] <- list(name = groups$organism_part[i],
                                           samples = groupNames,
                                           color = random_colors[1])
  names(groupList) <- make.unique(c(currentNames, groups$organism_part[i]))

  random_colors <- random_colors[-1]

}

# ::: Plot densities for defined groups/samples
tdensities <- plotIndividualDensitiesList(eventID = "HsaEX0019479",
                                          npoints = 500,
                                          psitable = psiTable,
                                          qualtable = qualTable,
                                          groupList = groupList,
                                          maxDevTable = maxDevSimulationN100)

# :::::::::::::::::::::::::::::::::::::::::::::::::::::
# Testing automatic group formation
# :::::::::::::::::::::::::::::::::::::::::::::::::::::

names <- colnames(psiTable)[-c(1:6)]

not_grouped   <- names
checked       <- c()
random_colors <- pastelColors
groups        <- LETTERS[seq(1, length = min(length(names), 26))]
used_groups   <- groups

groupList

while((length(checked) < length(names)) & length(used_groups) <= 26){

  groupNames  <- agrep(pattern = not_grouped[1], x = not_grouped, value = TRUE)

  if(length(groupNames) == 1 | length(groupNames) == length(names)){

  }else{

    checked     <- c(checked, groupNames)
    not_grouped <- not_grouped[-c(match(groupNames, not_grouped))]

    # Assign new groups (fake)
    currentNames <- names(groupList)
    groupList[[length(groupList)+1]] <- list(name = used_groups[1],
                                             samples = groupNames,
                                             color = random_colors[1])
    names(groupList) <- make.unique(c(currentNames, used_groups[1]))


    random_colors <- random_colors[-1]
    used_groups   <- used_groups[-1]

  }
}

groupList

# :::::::::::::::::::::::::::::::::::::::::::::::::::::
# 5. Perform betAS
# :::::::::::::::::::::::::::::::::::::::::::::::::::::

groupA  <- "heart"
groupB  <- "forebrain"
groupsAB <- c("heart", "forebrain")
yStat   <- "Pdiff (probability of differential splicing)"

samplesA    <- groupList[[groupA]]$samples
colsGroupA  <- match(samplesA, colnames(psiTable))

samplesB    <- groupList[[groupB]]$samples
colsGroupB  <- match(samplesB, colnames(psiTable))

#  Variable for volcano plot preparation
psitable = psiTable
qualtable = qualTable
npoints = 500
colsA = colsGroupA
colsB = colsGroupB
labA = groupA
labB = groupB
basalColor = "#89C0AE"
interestColor = "#E69A9C"
maxDevTable = maxDevSimulationN100

# Explore formation of individual betAS for one event

colsA    <- convertCols(psiTable, colsA)
samplesA <- names(colsA)

colsB    <- convertCols(psitable, colsB)
samplesB <- names(colsB)


x <- 10
eventPos <- 10
indBetasA <- individualBetas_nofitting_incr(table = qualTable[x,],
                               cols = colsA,
                               indpoints = npoints,
                               maxdevRefTable = maxDevTable)

indBetasB <- individualBetas_nofitting_incr(table = qualTable[x,],
                                            cols = colsB,
                                            indpoints = npoints,
                                            maxdevRefTable = maxDevTable)


x=c(1:10^6)
your.number=90000.43
which(abs(x-your.number)==min(abs(x-your.number)))


# ::: 5.1. Volcano plot with Pdiff as y-axis variable
volcanoTable <- prepareTableVolcano(psitable = psiTable,
                                    qualtable = qualTable,
                                    npoints = 500,
                                    colsA = colsGroupA,
                                    colsB = colsGroupB,
                                    labA = groupList$heart,
                                    labB = groupList$forebrain,
                                    basalColor = "#89C0AE",
                                    interestColor = "#E69A9C",
                                    maxDevTable = maxDevSimulationN100)

volcanoPdiff <- plotVolcano(betasTable = volcanoTable,
                            labA = groupList$heart,
                            labB = groupList$forebrain,
                            basalColor = "#89C0AE",
                            interestColor = "#E69A9C")

# ::: 5.2 Volcano plot with "F-stat." as y-axis variable
volcanoTableFstat <- prepareTableVolcanoFstat(psitable = test$PSI,
                                              qualtable = test$Qual,
                                              npoints = 500,
                                              colsA = colsGroupA,
                                              colsB = colsGroupB,
                                              labA = tgroupsList[[groupA]]$name,
                                              labB = tgroupsList[[groupB]]$name,
                                              basalColor = "#89C0AE",
                                              interestColor = "#E69A9C", #pink
                                              maxDevTable = maxDevSimulationN100)
volcanoFstat <- plotVolcanoFstat(betasTable = volcanoTableFstat,
                                 labA = tgroupsList[[groupA]]$name,
                                 labB = tgroupsList[[groupB]]$name,
                                 basalColor = "#89C0AE",
                                 interestColor = "#E69A9C")


# ::: 5.3. Volcano plot with FDR (Benilton's suggestion) as y-axis variable
volcanoTableFDR <- prepareTableVolcanoFDR(psitable = psiTable,
                                          qualtable = qualTable,
                                          npoints = 500,
                                          colsA = colsGroupA,
                                          colsB = colsGroupB,
                                          labA = groupList$heart,
                                          labB = groupList$forebrain,
                                          basalColor = "#89C0AE",
                                          interestColor = "#E69A9C", #pink
                                          maxDevTable = maxDevSimulationN100,
                                          nsim = 1000)

volcanoFDR <- plotVolcanoFstat(betasTable = volcanoTableFstat,
                               labA = tgroupsList[[groupA]]$name,
                               labB = tgroupsList[[groupB]]$name,
                               basalColor = "#89C0AE",
                               interestColor = "#E69A9C")


# :::::::::::::::::::::::::::::::::::::::::::::::::::::
# 6. Event-wise betAS
# :::::::::::::::::::::::::::::::::::::::::::::::::::::
eventID       <- "HsaEX0056290"
psitable      <- psi
qualtable     <- test$Qual
npoints       <- 500
colsA         <- colsGroupA
colsB         <- colsGroupB
labA          <- tgroupsList[[groupA]]$name
labB          <- tgroupsList[[groupB]]$name
basalColor    <- "#89C0AE"
interestColor <- "#E69A9C"
maxDevTable   <- maxDevSimulationN100
nsim          <- 1000


#Group betAS (A)
groupAbetAS <- individualBetas_nofitting_incr(table = qualtable[grep(eventID, qualtable$EVENT),],
                                              cols = colsA,
                                              indpoints = npoints,
                                              maxdevRefTable = maxDevTable)
#Group betAS (B)
groupBbetAS <- individualBetas_nofitting_incr(table = qualtable[grep(eventID, qualtable$EVENT),],
                                              cols = colsB,
                                              indpoints = npoints,
                                              maxdevRefTable = maxDevTable)

indBetasA <- groupAbetAS
indBetasB <- groupBbetAS
groupsAB  <- c(labA, labB)
nsim      <- 1000

tableEvent <- prepareTableEvent(eventID = "HsaEX0056290",
                                psitable = test$PSI,
                                qualtable = test$Qual,
                                npoints = 500,
                                colsA = colsGroupA,
                                colsB = colsGroupB,
                                labA = tgroupsList[[groupA]]$name,
                                labB = tgroupsList[[groupB]]$name,
                                basalColor = "#89C0AE",
                                interestColor = "#E69A9C", #pink
                                maxDevTable = maxDevSimulationN100,
                                nsim = 1000)

plotDensitiesFromEventObjList(eventObjList = tableEvent)
plotPDiffFromEventObjList(eventObjList = tableEvent)
plotFstatFromEventObjList(eventObjList = tableEvent)
plotFDRFromEventObjList(eventObjList = tableEvent)

eventObjList <- tableEvent
colorA <- "#89C0AE"
colorB <- "#E69A9C"

# :::::::::::::::::::::::::::::::::::::::::::::::::::::
# Explore FDR estimation
# :::::::::::::::::::::::::::::::::::::::::::::::::::::

plotSimulation <- function(i, color, sim, sampled){

  ggplot() +
    scale_x_continuous(breaks = seq(0, 1, 0.25),
                       limits = c(0,1)) +
    geom_density(aes(x = sim[[i]]),
                 fill = color,
                 show.legend = FALSE) +
    geom_point(aes(x = sampled[[i]],
                   y = runif(0, max(density(sim[[i]])$y), n = length(sampled[[i]]))),
               alpha = 0.5) +
    xlab(paste0(i))

}

simAPlots <- lapply(1:length(simulatedA),
                    function(x) plotSimulation(i = x, sim = simulatedA, sampled = sampledPointsA, color = groupList$heart$color))

simBPlots <- lapply(1:length(simulatedB),
                    function(x) plotSimulation(i = x, sim = simulatedB, sampled = sampledPointsB, color = groupList$forebrain$color))

plot_grid(plot_grid(plotlist = simAPlots, nrow = length(simulatedA)),
          plot_grid(plotlist = simBPlots, nrow = length(simulatedB)))

possiblefdrs <- seq(0, 1, 0.01)

plot(x = possiblefdrs, y = -log10(possiblefdrs))
