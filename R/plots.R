#' Convert between column names and indexes
#'
#' @param psitable PSI table
#' @param index (vector) column names or indexes
#'
#' @return named vector with column indexes and corresponding sample names (colnames)
#' @export
#'
#' @examples
#' testTable <- betAS:::testTable
#' psiTable   <- testTable$PSI
#' # get named vector with column indexes and corresponding sample names
#' convertCols(psiTable, index = c(10,11,13))
#' convertCols(psiTable, index = c("ERR2598269", "ERR2598270", "ERR2598271"))
convertCols <- function(psitable, index) {
  if (is.numeric(index)) {
    samples <- colnames(psitable[index])
  } else if (is.character(index)) {
    samples <- index
    index   <- match(samples, colnames(psitable))
  } else {
    stop("index needs to be a numeric or character vector")
  }
  names(index) <- samples
  return(index)
}

#' Density plots grid of all PSI values per sample
#'
#' @param table PSI table
#'
#' @return ggplot density plot grid
#' @export
#'
#' @importFrom reshape2 melt
#' @importFrom grDevices colorRampPalette
#' @import ggplot2
#'
#' @examples
#' testTable <- betAS:::testTable
#' psiTable   <- testTable$PSI
#' bigPicturePlot(table = psiTable)
bigPicturePlot <- function(table){

  transf <- melt(table[,-c(1,3:6)], id.vars = c("EVENT"))
  colors <- colorRampPalette(c("#F2969A","#55CC9D"))(length(unique(transf$variable)))

  plot <- ggplot(data = transf,
                 aes(x = value,
                     group = variable,
                     color = variable,
                     fill = variable)) +
    geom_density(show.legend = FALSE,
                 alpha = 0.8,
                 size = 2) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    xlab("Percent spliced-in (PSI)") +
    ylab("") +
    scale_x_continuous(breaks = seq(0,100, 25),
                       limits = c(0,100)) +
    facet_wrap(. ~ variable, ncol = 8) +
    theme_betAS() +
    theme(axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          strip.text = element_text(size = 14),
          plot.title = element_text(hjust = 0.5)) +
    ggtitle("PSI distribution per sample")

  return(plot)

}

#' Plot individual beta distributions per sample
#'
#' @param eventID
#' @param npoints
#' @param psitable
#' @param qualtable
#' @param colsA
#' @param colsB
#' @param labA
#' @param labB
#' @param colorA
#' @param colorB
#' @param seed (boolean) Boolean indicating if seed should be fixed. Default is TRUE
#' @param CoverageWeight (boolean)  Boolean indicating if number of points emitted in beta distributions should be weighted by the coverage of each sample. Default is FALSE
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom ggridges geom_density_ridges
#' @importFrom graphics title points
plotIndividualDensities <- function(eventID, npoints, psitable, qualtable, colsA, colsB, labA, labB, colorA, colorB, maxDevTable, seed=TRUE, CoverageWeight = FALSE){

  row <- which(psitable$EVENT == eventID)

  betasListA <- list()
  betasListB <- list()

  colsA    <- convertCols(psitable, colsA)
  samplesA <- names(colsA)

  colsB    <- convertCols(psitable, colsB)
  samplesB <- names(colsB)

  themeColors <- c("#89C0AE", "#E69A9C", "#76CAA0", "#EE805B", "#F7CF77", "#81C1D3")
  randomColorsA <- colorRampPalette(c(themeColors[1], themeColors[2]))(4)
  randomColorsB <- colorRampPalette(c(themeColors[3], themeColors[4]))(4)

  #Sample-wise betAS (A)
  for(a in 1:length(colsA)){

    indBetAS_A <- individualBetas_nofitting_incr(table = qualtable[row,],
                                                 cols = colsA[a],
                                                 indpoints = npoints,
                                                 maxdevRefTable = maxDevTable,
                                                 seed=seed,
                                                 CoverageWeight=CoverageWeight)

    betasListA[[a]] <- indBetAS_A

  }
  names(betasListA) <- samplesA

  #Sample-wise betAS (B)
  for(b in 1:length(colsB)){

    indBetAS_B <- individualBetas_nofitting_incr(table = qualtable[row,],
                                                 cols = colsB[b],
                                                 indpoints = npoints,
                                                 maxdevRefTable = maxDevTable,
                                                 seed=seed,
                                                 CoverageWeight=CoverageWeight)

    betasListB[[b]] <- indBetAS_B

  }
  names(betasListB) <- samplesB

  #Group betAS (A)
  groupAbetAS <- individualBetas_nofitting_incr(table = qualtable[row,],
                                                cols = colsA,
                                                indpoints = npoints,
                                                maxdevRefTable = maxDevTable,
                                                seed=seed,
                                                CoverageWeight=CoverageWeight)

  #Group betAS (B)
  groupBbetAS <- individualBetas_nofitting_incr(table = qualtable[row,],
                                                cols = colsB,
                                                indpoints = npoints,
                                                maxdevRefTable = maxDevTable,
                                                seed=seed,
                                                CoverageWeight=CoverageWeight)

  # Data-frame containing emitted points, sample names and group names
  densities_df  <- data.frame("points" = numeric(), "samples" = character())
  psis_df       <- data.frame("psis" = numeric(), "samples" = character())

  for(i in 1:length(colsA)){

    #prepare table with emitted points
    points_df_A   <- data.frame("points" = matrix(unlist(betasListA[[i]]$BetaPoints), nrow = npoints, byrow=TRUE), stringsAsFactors=FALSE)
    sample_df_A   <- cbind("points" = points_df_A, "samples" = rep(samplesA[i], times = length(points_df_A)))
    densities_df  <- rbind(densities_df, sample_df_A)

    #prepare table with psis
    psis_df_A       <- data.frame("psis" = matrix(unlist(betasListA[[i]]$PSI), nrow = 1, byrow=TRUE), stringsAsFactors=FALSE)
    samplepsis_df_A <- cbind("psis" = psis_df_A, "samples" = rep(samplesA[i], times = 1))
    psis_df <- rbind(psis_df, samplepsis_df_A)

  }

  for(j in 1:length(colsB)){

    #prepare table with emitted points
    points_df_B   <- data.frame("points" = matrix(unlist(betasListB[[j]]$BetaPoints), nrow = npoints, byrow=TRUE), stringsAsFactors=FALSE)
    sample_df_B   <- cbind("points" = points_df_B, "samples" = rep(samplesB[j], times = length(points_df_B)))
    densities_df  <- rbind(densities_df, sample_df_B)

    #prepare table with psis
    psis_df_B       <- data.frame("psis" = matrix(unlist(betasListB[[j]]$PSI), nrow = 1, byrow=TRUE), stringsAsFactors=FALSE)
    samplepsis_df_B <- cbind("psis" = psis_df_B, "samples" = rep(samplesB[j], times = 1))
    psis_df <- rbind(psis_df, samplepsis_df_B)

  }

  densities_df$group <- as.character(c(rep(labA, times = npoints*length(colsA)),
                                       rep(labB, times = npoints*length(colsB))))

  psis_df$group <- as.character(c(rep(labA, times = length(colsA)),
                                  rep(labB, times = length(colsB))))

  plot <- ggplot(densities_df,
         aes(x = points,
             y = samples,
             group = samples,
             color = group,
             fill = group)) +
    geom_density_ridges(size = 2,
                        alpha = 0.8,
                        show.legend = FALSE,
                        rel_min_height = 0.001) +
    geom_rug(data = psis_df,
             aes(x = psis,
                 color = group),
             size = 3,
             alpha = 0.8,
             show.legend = FALSE) +
    scale_color_manual(values = c(colorA, colorB)) +
    scale_fill_manual(values = c(colorA, colorB)) +
    xlab("Proportion spliced-in (PSI)") +
    ylab("") +
    labs(title ) +
    scale_x_continuous(breaks = seq(0,1, 0.25), limits = c(0,1)) +
    theme_betAS()
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line.x = element_line(color = "black"), axis.text.y = element_text(size = 16))

  return(plot)

}


#' Plot individual beta distributions per sample from groups list
#'
#' @param eventID (character) event ID
#' @param npoints number of points emitted by each beta distribution
#' @param psitable PSI table
#' @param qualtable Qual table
#' @param groupList group list
#' @param maxDevTable table with increments to be summed to zero-valued reads, per coverage
#' @param seed (boolean) Boolean indicating if seed should be fixed. Default is TRUE
#' @param CoverageWeight (boolean)  Boolean indicating if number of points emitted in beta distributions should be weighted by the coverage of each sample. Default is FALSE
#'
#' @return ggplot density plot grid
#' @export
#'
#' @importFrom ggridges geom_density_ridges
#' @examples
#' testTable <- betAS:::testTable
#' maxDevSimulationN100  <- betAS:::maxDevSimulationN100
#' psiTable   <- testTable$PSI
#' qualTable  <- testTable$Qual
#' testGroups <- list()
#' testGroups[[length(testGroups)+1]] <- list(name = "GroupA", samples = c("ERR2598266", "ERR2598267", "ERR2598268"), color = "#FF9AA2")
#' testGroups[[length(testGroups)+1]] <- list(name = "GroupB", samples = c("ERR2598270", "ERR2598307", "ERR2598351"), color = "#FFB7B2")
#' plotIndividualDensitiesList(eventID = "HsaEX0019479", npoints = 500, psitable = psiTable, qualtable = qualTable, groupList = testGroups, maxDevTable = maxDevSimulationN100)
plotIndividualDensitiesList <- function(eventID, npoints, psitable, qualtable, groupList, maxDevTable, seed=TRUE, CoverageWeight = FALSE){
  # Number of emmited points is fixed for visualization purposes, regardless of the value in CoverageWeight;
  # to be updated in the future if needed, npoints should also come as an output from individualBetas_nofitting_incr


   row <- which(psitable$EVENT == eventID)

  # Data-frame containing emitted points, sample names and group names
  densities_df  <- data.frame("points" = numeric(), "samples" = character())
  psis_df       <- data.frame("psis" = numeric(), "samples" = character())

  betasPerGroup           <- list()
  betasPerSamplePerGroup  <- list()

  definedColorsPerGroup <- c()
  definedGroups <- names(groupList)

  # for each defined group
  for(g in 1:length(groupList)){

    groupName     <- groupList[[g]]$name
    groupSamples  <- groupList[[g]]$samples
    groupColor    <- groupList[[g]]$color

    columns       <- convertCols(psitable, groupSamples)
    samples       <- names(columns)

    groupBetas    <- individualBetas_nofitting_incr(table = qualtable[row,],
                                                    cols = columns,
                                                    indpoints = npoints,
                                                    maxdevRefTable = maxDevTable,
                                                    seed=seed,
                                                    CoverageWeight = CoverageWeight)
    betasPerGroup[[g]] <- groupBetas

    groupIndSamplesList <- list()
    for(samp in 1:length(samples)){

      indBetasGroup <- individualBetas_nofitting_incr(table = qualtable[row,],
                                                      cols = columns[samp],
                                                      indpoints = npoints,
                                                      maxdevRefTable = maxDevTable,
                                                      seed=seed,
                                                      CoverageWeight = CoverageWeight)

      groupIndSamplesList[[samp]] <- indBetasGroup

      #prepare table with emitted points
      # Number of emmited points is fixed for visualization purposes, regardless of the value in CoverageWeight; to be updated in the future if needed, npoints should also come as an output from individualBetas_nofitting_incr
      points_df     <- data.frame("points" = matrix(unlist(indBetasGroup$BetaPoints), nrow = npoints, byrow=TRUE), stringsAsFactors=FALSE)
      sample_df     <- cbind("points" = points_df, "samples" = rep(samples[samp], times = length(points_df)), "group" = rep(groupName, times = length(points_df)))
      densities_df  <- rbind(densities_df, sample_df)

      #prepare table with psis
      psis          <- data.frame("psis" = matrix(unlist(indBetasGroup$PSI), nrow = 1, byrow=TRUE), stringsAsFactors=FALSE)
      samplepsis_df <- cbind("psis" = psis, "samples" = rep(samples[samp], times = 1), "group" = rep(groupName, times =  1))
      psis_df       <- rbind(psis_df, samplepsis_df)

    }

    names(groupIndSamplesList) <- samples
    betasPerSamplePerGroup[[g]] <- groupIndSamplesList
    definedColorsPerGroup <- c(definedColorsPerGroup, groupColor)

  }

  plot <- ggplot(densities_df,
                 aes(x = points,
                     y = samples,
                     # group = samples,
                     color = group,
                     fill = group)) +
    geom_density_ridges(size = 2,
                        alpha = 0.8,
                        show.legend = FALSE,
                        rel_min_height = 0.001) +
    geom_rug(data = psis_df,
             aes(x = psis,
                 color = group),
             size = 3,
             alpha = 0.8,
             show.legend = FALSE) +
    scale_color_manual(name = "", values = definedColorsPerGroup, breaks = definedGroups, limits =  definedGroups) +
    scale_fill_manual(name = "", values = definedColorsPerGroup, breaks = definedGroups, limits =  definedGroups) +
    xlab("Proportion spliced-in (PSI)") +
    ylab("") +
    labs(title ) +
    scale_x_continuous(breaks = seq(0,1, 0.1), limits = c(0,1)) +
    theme_betAS() +
    theme(legend.position = "bottom",
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 16)) +
    facet_grid(group ~ .,
               scales = "free_y")
  # theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line.x = element_line(color = "black"))

  return(plot)

}

#' Plot individual beta distributions as violin plots per sample from groups list
#'
#' @param eventID (character) event ID
#' @param npoints number of points emitted by each beta distribution
#' @param psitable PSI table
#' @param qualtable Qual table
#' @param groupList group list
#' @param maxDevTable table with increments to be summed to zero-valued reads, per coverage
#' @param seed (boolean) Boolean indicating if seed should be fixed. Default is TRUE
#' @param CoverageWeight (boolean)  Boolean indicating if number of points emitted in beta distributions should be weighted by the coverage of each sample. Default is FALSE
#'
#' @return ggplot density plot grid
#' @export
#'
#' @examples
plotIndividualViolinsList <- function(eventID, npoints, psitable, qualtable, groupList, maxDevTable, seed=TRUE, CoverageWeight = FALSE){

  row <- which(psitable$EVENT == eventID)

  # Data-frame containing emitted points, sample names and group names
  densities_df  <- data.frame("points" = numeric(), "samples" = character())
  psis_df       <- data.frame("psis" = numeric(), "samples" = character())

  betasPerGroup           <- list()
  betasPerSamplePerGroup  <- list()

  definedColorsPerGroup <- c()
  definedGroups <- names(groupList)

  # for each defined group
  for(g in 1:length(groupList)){

    groupName     <- groupList[[g]]$name
    groupSamples  <- groupList[[g]]$samples
    groupColor    <- groupList[[g]]$color

    columns       <- convertCols(psitable, groupSamples)
    samples       <- names(columns)

    groupBetas    <- individualBetas_nofitting_incr(table = qualtable[row,],
                                                    cols = columns,
                                                    indpoints = npoints,
                                                    maxdevRefTable = maxDevTable, seed=seed, CoverageWeight=CoverageWeight)
    betasPerGroup[[g]] <- groupBetas

    groupIndSamplesList <- list()
    for(samp in 1:length(samples)){

      indBetasGroup <- individualBetas_nofitting_incr(table = qualtable[row,],
                                                      cols = columns[samp],
                                                      indpoints = npoints,
                                                      maxdevRefTable = maxDevTable, seed=seed, CoverageWeight=CoverageWeight)

      groupIndSamplesList[[samp]] <- indBetasGroup

      #prepare table with emitted points
      points_df     <- data.frame("points" = matrix(unlist(indBetasGroup$BetaPoints), nrow = npoints, byrow=TRUE), stringsAsFactors=FALSE)
      sample_df     <- cbind("points" = points_df, "samples" = rep(samples[samp], times = length(points_df)), "group" = rep(groupName, times = length(points_df)))
      densities_df  <- rbind(densities_df, sample_df)

      #prepare table with psis
      psis          <- data.frame("psis" = matrix(unlist(indBetasGroup$PSI), nrow = 1, byrow=TRUE), stringsAsFactors=FALSE)
      samplepsis_df <- cbind("psis" = psis, "samples" = rep(samples[samp], times = 1), "group" = rep(groupName, times =  1))
      psis_df       <- rbind(psis_df, samplepsis_df)

    }

    names(groupIndSamplesList) <- samples
    betasPerSamplePerGroup[[g]] <- groupIndSamplesList
    definedColorsPerGroup <- c(definedColorsPerGroup, groupColor)

  }

  plot <- ggplot(densities_df,
                 aes(x = samples,
                     y = points,
                     # group = samples,
                     color = group,
                     fill = group)) +
    geom_violin(alpha = 0.6,
                show.legend = FALSE) +
    scale_color_manual(name = "",
                       values = definedColorsPerGroup,
                       breaks = definedGroups,
                       limits =  definedGroups) +
    scale_fill_manual(name = "",
                      values = definedColorsPerGroup,
                      breaks = definedGroups,
                      limits =  definedGroups) +
    xlab("Samples") +
    ylab("Proportion spliced-in (PSI)") +
    labs(title) +
    scale_y_continuous(breaks = seq(0,1, 0.25), limits = c(0,1)) +
    theme_betAS() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(size = 10, angle = 45)) +
    facet_grid(. ~ group,
               scales = "free_x")

  return(plot)

}


#' Prepare and return table used for betAS volcano plot (probability of splicing differences as y-axis)
#'
#' @param psitable PSI table
#' @param qualtable Qual table
#' @param npoints number of emitted points generated per beta distribution
#' @param colsA column indexes for samples in group A
#' @param colsB column indexes for samples in group B
#' @param labA group A label
#' @param labB group B label
#' @param basalColor general color for points (events)
#' @param interestColor color for highlighted points (events)
#' @param maxDevTable (data.frame) reference data frame with maximum increment values per coverage (inc+exc) to avoid emitted values to be artifitially beyond a certain threshold
#' @param seed (boolean) Boolean indicating if seed should be fixed. Default is TRUE
#' @param CoverageWeight (boolean)  Boolean indicating if number of points emitted in beta distributions should be weighted by the coverage of each sample. Default is FALSE
#'
#'
#' @return data table to be used as ggplot input for plotting
#' @export
#'
#' @examples
#' testTable <- betAS:::testTable
#' maxDevSimulationN100  <- betAS:::maxDevSimulationN100
#' psiTable   <- testTable$PSI
#' qualTable  <- testTable$Qual
#' testGroups <- list()
#' testGroups[[length(testGroups)+1]] <- list(name = "GroupA", samples = c("ERR2598266", "ERR2598267", "ERR2598268"), color = "#FF9AA2")
#' testGroups[[length(testGroups)+1]] <- list(name = "GroupB", samples = c("ERR2598270", "ERR2598307", "ERR2598351"), color = "#FFB7B2")
#' groupA    <- "GroupA"
#' groupB    <- "GroupB"
#' samplesA   <- testGroups[[1]]$samples
#' samplesB   <- testGroups[[2]]$samples
#' colsGroupA    <- convertCols(psiTable, samplesA)
#' colsGroupB    <- convertCols(psiTable, samplesB)
#' prepareTableVolcano(psitable = psiTable, qualtable = qualTable, npoints = 500, colsA = colsGroupA, colsB = colsGroupB, labA = groupA, labB = groupB, basalColor = "#89C0AE", interestColor = "#E69A9C", maxDevTable = maxDevSimulationN100)
prepareTableVolcano <- function(psitable, qualtable, npoints, colsA, colsB, labA, labB, basalColor, interestColor, maxDevTable, seed=TRUE, CoverageWeight=CoverageWeight){

  colsA    <- convertCols(psitable, colsA)
  samplesA <- names(colsA)

  colsB    <- convertCols(psitable, colsB)
  samplesB <- names(colsB)

  #Group betAS (A)
  groupAbetAS <- lapply(1:nrow(qualtable),
                        function(x)
                          individualBetas_nofitting_incr(table = qualtable[x,],
                                                         cols = colsA,
                                                         indpoints = npoints,
                                                         maxdevRefTable = maxDevTable,
                                                         seed=seed, CoverageWeight=CoverageWeight))

  #Group betAS (B)
  groupBbetAS <- lapply(1:nrow(qualtable),
                        function(x)
                          individualBetas_nofitting_incr(table = qualtable[x,],
                                                         cols = colsB,
                                                         indpoints = npoints,
                                                         maxdevRefTable = maxDevTable,
                                                         seed=seed, CoverageWeight=CoverageWeight))

  #Differential betAS (A vs. B)
  diffABbetAS <- lapply(1:nrow(qualtable),
                        function(x)
                          jointBetas_nofitting_Fast0(indBetasA = groupAbetAS[[x]],
                                                     indBetasB = groupBbetAS[[x]],
                                                     groupsAB = c(labA, labB)))

  names(groupAbetAS) <- qualtable$EVENT
  names(groupBbetAS) <- qualtable$EVENT
  names(diffABbetAS) <- qualtable$EVENT

  betasTable            <- psitable
  betasTable$Pdiff      <- as.numeric(as.vector(lapply(1:nrow(qualtable), function(x) diffABbetAS[[x]][[2]])))
  betasTable$betasPsiA  <- as.numeric(as.vector(lapply(1:nrow(qualtable), function(x) groupAbetAS[[x]]$MedianBeta)))
  betasTable$betasPsiB  <- as.numeric(as.vector(lapply(1:nrow(qualtable), function(x) groupBbetAS[[x]]$MedianBeta)))
  betasTable$deltapsi   <- betasTable$betasPsiB - betasTable$betasPsiA

  return(betasTable)

}

#' Plot betAS volcano plot
#'
#'
#' @param betasTable data table generated with prepareTableVolcano
#' @param labA group A label
#' @param labB group B label
#' @param basalColor general color for points (events)
#' @param interestColor color for highlighted points (events)
#'
#' @return ggplot scatterplot
#' @export
#'
#' @importFrom ggrepel geom_text_repel
#' @examples
#' testTable <- betAS:::testTable
#' maxDevSimulationN100  <- betAS:::maxDevSimulationN100
#' psiTable   <- testTable$PSI
#' qualTable  <- testTable$Qual
#' testGroups <- list()
#' testGroups[[length(testGroups)+1]] <- list(name = "GroupA", samples = c("ERR2598266", "ERR2598267", "ERR2598268"), color = "#FF9AA2")
#' testGroups[[length(testGroups)+1]] <- list(name = "GroupB", samples = c("ERR2598270", "ERR2598307", "ERR2598351"), color = "#FFB7B2")
#' groupA    <- "GroupA"
#' groupB    <- "GroupB"
#' samplesA   <- testGroups[[1]]$samples
#' samplesB   <- testGroups[[2]]$samples
#' colsGroupA    <- convertCols(psiTable, samplesA)
#' colsGroupB    <- convertCols(psiTable, samplesB)
#' volcanoTable <- prepareTableVolcano(psitable = psiTable, qualtable = qualTable, npoints = 500, colsA = colsGroupA, colsB = colsGroupB, labA = groupA, labB = groupB, basalColor = "#89C0AE", interestColor = "#E69A9C", maxDevTable = maxDevSimulationN100)
#' plotVolcano(betasTable = volcanoTable, labA = groupA, labB = groupB, basalColor = "#89C0AE", interestColor = "#E69A9C")
plotVolcano <- function(betasTable, labA, labB, basalColor, interestColor){

  refSeq      <- seq(0, 1.1, 0.1) # from 0 to 1.1, the latter needed to include deltaPSI=1
  maxDeltaPsi <- max(abs(betasTable$deltapsi))
  pos         <- which(abs(refSeq-maxDeltaPsi)==min(abs(refSeq-maxDeltaPsi)))
  refScale    <- refSeq[pos+1]

  ggplot(betasTable,
         aes(x = deltapsi,
             y = Pdiff)) +
    # y = -log10(1.001-Pdiff))) +
    geom_point(alpha = 0.3,
               color = basalColor,
               size = 4) +
    geom_point(data = betasTable[which(abs(betasTable$deltapsi) > 0.1),],
               size = 3,
               color = interestColor) +
    geom_text_repel(data = betasTable[which(abs(betasTable$deltapsi) > 0.1),],
                    color = interestColor,
                    aes(label = paste0(EVENT, "\n(", GENE, ")")),
                    size = 5) +
    scale_x_continuous(breaks = seq(-refScale, refScale, 0.1), limits = c(-refScale, refScale)) +
    # xlab(expression(PSI[B]-PSI[A])) +
    xlab(paste0("PSI(", labB, ") - PSI(", labA, ")")) +
    # ylab(expression(-log[10](paste("1"^"+")-P[diff.~splicing]))) +
    ylab(paste0("Probability of differential splicing between ", labA, " and ", labB)) +
    theme_betAS() +
    theme(axis.title.y = element_text(size = 24))

}


#' Prepare and return table used for betAS volcano plot with F-statistic for the y-axis
#'
#' @param psitable
#' @param qualtable
#' @param npoints
#' @param colsA
#' @param colsB
#' @param basalColor
#' @param interestColor
#' @param seed (boolean) Boolean indicating if seed should be fixed. Default is TRUE
#' @param CoverageWeight (boolean)  Boolean indicating if number of points emitted in beta distributions should be weighted by the coverage of each sample. Default is FALSE
#'
#' @return
#' @export
#'
#' @examples
prepareTableVolcanoFstat <- function(psitable, qualtable, npoints, colsA, colsB, labA, labB, basalColor, interestColor, maxDevTable, seed=TRUE, CoverageWeight=FALSE){

  colsA    <- convertCols(psitable, colsA)
  samplesA <- names(colsA)

  colsB    <- convertCols(psitable, colsB)
  samplesB <- names(colsB)

  #F-stat (2 groups) per event
  fstat2groups <- lapply(1:nrow(qualtable),
                         function(x)
                           fStatistic_2Groups_PerEvent(eventPos = x,
                                                       psitable = psitable,
                                                       qualtable = qualtable,
                                                       npoints = 500,
                                                       colsA = colsA,
                                                       colsB = colsB,
                                                       labA = labA,
                                                       labB = labB,
                                                       maxDevTable = maxDevTable,
                                                       seed=seed, CoverageWeight=CoverageWeight))


  names(fstat2groups) <- qualtable$EVENT

  betasTable            <- psitable
  betasTable$deltapsi   <- as.numeric(as.vector(lapply(1:nrow(qualtable), function(x) fstat2groups[[x]]$deltaPsi)))
  betasTable$Fstat      <- as.numeric(as.vector(lapply(1:nrow(qualtable), function(x) fstat2groups[[x]]$Fstat)))

  return(betasTable)

}

#' Plot betAS volcano plot with F-statistic for the y-axis
#'
#' @param psitable
#' @param qualtable
#' @param npoints
#' @param colsA
#' @param colsB
#' @param basalColor
#' @param interestColor
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom ggrepel geom_text_repel
plotVolcanoFstat <- function(betasTable, labA, labB, basalColor, interestColor){

  refSeq      <- seq(0,1.1, 0.1) # from 0 to 1.1, the latter needed to include deltaPSI=1
  maxDeltaPsi <- max(abs(betasTable$deltapsi))
  pos         <- which(abs(refSeq-maxDeltaPsi)==min(abs(refSeq-maxDeltaPsi)))
  refScale    <- refSeq[pos+1]

  ggplot(betasTable,
         aes(x = deltapsi,
             y = Fstat)) +
    geom_point(alpha = 0.3,
               color = basalColor,
               size = 4) +
    geom_point(data = betasTable[which(abs(betasTable$deltapsi) > 0.1),],
               size = 3,
               color = interestColor) +
    geom_text_repel(data = betasTable[which(abs(betasTable$deltapsi) > 0.1),],
                    color = interestColor,
                    aes(label = paste0(EVENT, "\n(", GENE, ")")),
                    size = 5) +
    scale_x_continuous(breaks = seq(-refScale,refScale,0.1), limits = c(-refScale, refScale)) +
    # xlab(expression(PSI[B]-PSI[A])) +
    xlab(paste0("PSI(", labB, ") - PSI(", labA, ")")) +
    ylab(expression("F-statistic")) +
    theme_betAS()


}


#' Prepare and return table used for betAS volcano plot with the estimated FPR
#'
#' @param psitable
#' @param qualtable
#' @param npoints
#' @param colsA
#' @param colsB
#' @param basalColor
#' @param interestColor
#' @param seed (boolean) Boolean indicating if seed should be fixed. Default is TRUE
#' @param CoverageWeight (boolean)  Boolean indicating if number of points emitted in beta distributions should be weighted by the coverage of each sample. Default is FALSE
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom ggrepel geom_text_repel
prepareTableVolcanoFDR <- function(psitable, qualtable, npoints, colsA, colsB, labA, labB, basalColor, interestColor, maxDevTable, nsim, seed=TRUE, CoverageWeight=FALSE){

  colsA    <- convertCols(psitable, colsA)
  samplesA <- names(colsA)

  colsB    <- convertCols(psitable, colsB)
  samplesB <- names(colsB)

  #Group betAS (A)
  groupAbetAS <- lapply(1:nrow(qualtable),
                        function(x)
                          individualBetas_nofitting_incr(table = qualtable[x,],
                                                         cols = colsA,
                                                         indpoints = npoints,
                                                         maxdevRefTable = maxDevTable, seed=seed, CoverageWeight=CoverageWeight))

  #Group betAS (B)
  groupBbetAS <- lapply(1:nrow(qualtable),
                        function(x)
                          individualBetas_nofitting_incr(table = qualtable[x,],
                                                         cols = colsB,
                                                         indpoints = npoints,
                                                         maxdevRefTable = maxDevTable, seed=seed, CoverageWeight=CoverageWeight))

  #Differential betAS (A vs. B)
  diffABbetAS <- lapply(1:nrow(qualtable),
                        function(x)
                          estimateFDR(indBetasA = groupAbetAS[[x]],
                                      indBetasB = groupBbetAS[[x]],
                                      groupsAB = c(labA, labB),
                                      nsim = nsim, seed=seed, CoverageWeight=CoverageWeight))

  names(groupAbetAS) <- qualtable$EVENT
  names(groupBbetAS) <- qualtable$EVENT
  names(diffABbetAS) <- qualtable$EVENT

  betasTable              <- psitable
  betasTable$FDR          <- as.numeric(as.vector(lapply(1:nrow(qualtable), function(x) diffABbetAS[[x]][[3]])))
  betasTable$invertedFDR  <- 1 - betasTable$FDR
  betasTable$betasPsiA    <- as.numeric(as.vector(lapply(1:nrow(qualtable), function(x) groupAbetAS[[x]]$MedianBeta)))
  betasTable$betasPsiB    <- as.numeric(as.vector(lapply(1:nrow(qualtable), function(x) groupBbetAS[[x]]$MedianBeta)))
  betasTable$deltapsi     <- betasTable$betasPsiB - betasTable$betasPsiA

  return(betasTable)

}

#' Plot betAS volcano plot with the estimated FPR
#'
#' @param psitable
#' @param qualtable
#' @param npoints
#' @param colsA
#' @param colsB
#' @param basalColor
#' @param interestColor
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom ggrepel geom_text_repel
plotVolcanoFDR <- function(betasTable, labA, labB, basalColor, interestColor){

  refSeq      <- seq(0,1.1, 0.1) # from 0 to 1.1, the latter needed to include deltaPSI=1
  maxDeltaPsi <- max(abs(betasTable$deltapsi))
  pos         <- which(abs(refSeq-maxDeltaPsi)==min(abs(refSeq-maxDeltaPsi)))
  refScale    <- refSeq[pos+1]

  # To be used if intented to use -log10 of FPR and still manage to see Inf points
  # betasTableTransf <- betasTable
  #
  # betasTableTransf$TransfFDR <- -log10(betasTableTransf$FDR)
  # betasTableTransf$TransfFDR[!is.finite(betasTableTransf$TransfFDR)] <- max(betasTableTransf$TransfFDR[is.finite(betasTableTransf$TransfFDR)])+0.0001
  # transfScale <- c(seq(from = 0, to = 2, by = 0.5), 2+0.0001)
  # transfScaleLabels <- c(seq(from = 0, to = 2, by = 0.5), Inf)

  ggplot(betasTable,
         aes(x = deltapsi,
             y = invertedFDR)) +
             # y = -log10(FDR))) +
    # y = -log10(1.001-Pdiff))) +
    geom_point(alpha = 0.3,
               color = basalColor,
               size = 4) +
    geom_point(data = betasTable[which(abs(betasTable$deltapsi) > 0.1),],
               size = 3,
               color = interestColor) +
    geom_text_repel(data = betasTable[which(abs(betasTable$deltapsi) > 0.1),],
                    color = interestColor,
                    aes(label = paste0(EVENT, "\n(", GENE, ")")),
                    size = 5) +
    scale_x_continuous(breaks = seq(-refScale, refScale, 0.1), limits = c(-refScale, refScale)) +
    # scale_y_continuous(breaks = transfScale, labels = transfScaleLabels) +
    # xlab(expression(PSI[B]-PSI[A])) +
    xlab(paste0("PSI(", labB, ") - PSI(", labA, ")")) +
    ylab("1 - False Discovery Rate") +
    # ylab(expression(-log[10](estimated~FDR+paste("0"^"+")))) +
    # ylab(expression(-log[10](estimated~FDR))) +
    theme_betAS()

}

#' Prepare and return table used for betAS volcano plot (multiple group analysis)
#'
#' @param psitable
#' @param qualtable
#' @param groupList
#' @param npoints
#' @param maxDevTable
#' @param seed (boolean) Boolean indicating if seed should be fixed. Default is TRUE
#' @param CoverageWeight (boolean)  Boolean indicating if number of points emitted in beta distributions should be weighted by the coverage of each sample. Default is FALSE
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom ggrepel geom_text_repel
prepareTableVolcanoMultipleGroups <- function(psitable, qualtable, groupList, npoints, maxDevTable, seed=TRUE,CoverageWeight=FALSE){

  # Prepare individual betAS object per group
  groupNames  <- names(groupList)
  listNames   <- groupNames
  #listNames   <- paste0(groupNames, "_indBetas")
  samplesPerGroupList <- list()
  indBetasList        <- list()

  for(i in 1:length(groupNames)){

    samples <- groupList[[i]]$samples
    pos <- match(samples, colnames(psitable))

    samplesPerGroupList[[i]] <- samples

    # Prepare individual betAS object for this group
    indbetas <- lapply(1:nrow(qualtable),
                         function(x)
                           individualBetas_nofitting_incr(table = qualtable[x,],
                                                          cols = pos,
                                                          indpoints = npoints,
                                                          maxdevRefTable = maxDevTable,
                                                          seed=seed, CoverageWeight=CoverageWeight))
    # Name each position in list after the event
    names(indbetas) <- qualtable$EVENT

    # Assign ind betas to global list
    indBetasList[[i]] <- indbetas

  }

  names(indBetasList) <- listNames
  names(samplesPerGroupList) <- listNames

  # Prepare directional ANOVA-like approach
  genGroupsBetas <- lapply(1:length(indBetasList[[1]]),
                             function(x)
                               generalisedGroupsBetas(eventPos = x,
                                                      indList = indBetasList,
                                                      samplesList = samplesPerGroupList,
                                                      groupNames = groupNames,
                                                      seed=seed,
                                                      CoverageWeight=CoverageWeight))

  # Name each position in list after the event
  names(genGroupsBetas) <- names(indBetasList[[1]])

  multipleGroupTable <- psitable

  multipleGroupTable$Fstat             <- as.numeric(as.vector(lapply(1:length(indBetasList[[1]]), function(x) genGroupsBetas[[x]][[3]])))
  multipleGroupTable$Pzero             <- as.numeric(as.vector(lapply(1:length(indBetasList[[1]]), function(x) genGroupsBetas[[x]][[4]])))
  multipleGroupTable$Pdiff             <- as.numeric(as.vector(lapply(1:length(indBetasList[[1]]), function(x) genGroupsBetas[[x]][[5]])))
  multipleGroupTable$medianBetweens    <- as.numeric(as.vector(lapply(1:length(indBetasList[[1]]), function(x) genGroupsBetas[[x]][[6]])))
  multipleGroupTable$deltaAbsolute     <- as.numeric(as.vector(lapply(1:length(indBetasList[[1]]), function(x) genGroupsBetas[[x]][[7]])))

  return(multipleGroupTable)

}


#' Plot betAS volcano plot (with Pdiff as y-axis metric)
#'
#'
#' @param betasTable data table generated with prepareTableVolcanoMultipleGroups
#' @param basalColor general color for points (events)
#' @param interestColor color for highlighted points (events)
#'
#' @return ggplot scatterplot
#' @export
#'
#' @importFrom ggrepel geom_text_repel
#' @examples
plotVolcano_MultipleGroups_Pdiff <- function(betasTable){

  ggplot(betasTable,
         aes(x = deltaAbsolute,
             y = Pdiff,
             color = Fstat)) +
    geom_point(size = 4,
               alpha = 0.8) +
    scale_color_gradient2(low = "white",
                          high = "#DC143C",
                          name = "F statistic") +
    # y = -log10(1.001-Pdiff))) +
    geom_point(data = betasTable[which(abs(betasTable$deltaAbsolute) > 0.1),],
               size = 4,
               shape = 21) +
    geom_text_repel(data = betasTable[which(abs(betasTable$deltaAbsolute) > 0.1),],
                    aes(label = paste0(EVENT, "\n(", GENE, ")")),
                    size = 5) +
    scale_x_continuous(breaks = seq(0,0.5,0.1), limits = c(-0.01,0.5)) +
    scale_y_continuous(breaks = seq(0.5,1,0.1), limits = c(0.5,1)) +
    xlab(expression(median["|between|"]-median["|within|"]~(PSI))) +
    ylab("Prob. |between| > |within|") +
    theme_betAS() +
    theme(legend.position = "top",
          legend.background = element_blank(),
          legend.key.width = unit(6, 'cm'),
          legend.direction = "horizontal")

}

#' Plot betAS volcano plot (with Fstat as y-axis metric)
#'
#'
#' @param betasTable data table generated with prepareTableVolcanoMultipleGroups
#' @param basalColor general color for points (events)
#' @param interestColor color for highlighted points (events)
#'
#' @return ggplot scatterplot
#' @export
#'
#' @importFrom ggrepel geom_text_repel
#' @examples
plotVolcano_MultipleGroups_Fstat <- function(betasTable){

  # maxFstat <- max(abs(betasTable$Fstat), na.rm = TRUE)

  ggplot(betasTable,
         aes(x = deltaAbsolute,
             y = Fstat,
             color = Pdiff)) +
    geom_point(size = 4,
               alpha = 0.8) +
    scale_color_gradient(low = "white",
                         high = "blue",
                         limits = c(0.5, 1),
                         name = "Prob. |between| > |within|") +
    # y = -log10(1.001-Pdiff))) +
    geom_point(data = betasTable[which(abs(betasTable$deltaAbsolute) > 0.1),],
               shape = 21,
               size = 4) +
    geom_text_repel(data = betasTable[which(abs(betasTable$deltaAbsolute) > 0.1),],
                    aes(label = paste0(EVENT, "\n(", GENE, ")")),
                    size = 5) +
    scale_x_continuous(breaks = seq(0,0.5,0.1), limits = c(-0.01,0.5)) +
    # scale_y_continuous(breaks = seq(0, round(maxFstat), 2), limits = c(0, maxFstat)) +
    xlab(expression(median["|between|"]-median["|within|"]~(PSI))) +
    ylab("F statistic") +
    theme_betAS() +
    theme(legend.position = "top",
          legend.background = element_blank(),
          legend.key.width = unit(6, 'cm'),
          legend.direction = "horizontal")

}



#' Prepare pie chart
#'
#' @param table PSI table
#'
#' @return highcharts pie chart
#' @export
#'
#' @import highcharter
#' @examples
#' testTable <- betAS:::testTable
#' preparePieForVastToolsCOMPLEX(table = testTable$PSI)
preparePieForVastToolsCOMPLEX <- function(table){

  # Set highcharter options
  options(highcharter.theme = hc_theme_smpl_tailored(tooltip = list(valueDecimals = 0),
                                                     chart = list(backgroundColor = "transparent")))

  df <- as.data.frame(table(table$COMPLEX))
  df <- df[,c(2,1)]
  df <- cbind(seq(0, nrow(df)-1), df)
  colnames(df) <- c("x", "y", "name")

  hc <- df %>%
    hchart(
      "pie", hcaes(x = name, y = y),
      name = "Number of events:"
    )

  hc

}

# Adapted from theme_minimal()
theme_betAS <- function(base_size = 30,
                        base_family = "",
                        base_line_size = base_size/22,
                        base_rect_size = base_size/22){
  theme_bw(base_size = base_size,
           base_family = base_family,
           base_line_size = base_line_size,
           base_rect_size = base_rect_size) %+replace%
    theme(axis.ticks = element_blank(),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          strip.background = element_blank(),
          plot.background = element_blank(),
          axis.text = element_text(size = 20), complete = TRUE)
}

# Adapted from theme_clean()
theme_clean20 <- function(base_size = 20, base_family = ""){
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.border     = element_blank(),
      axis.line        = element_line(colour = "black"),
      panel.grid.major = element_line(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.background = element_blank(),
      strip.background = element_blank()
    )
}

#' Adapted from highcharts themes: "simple theme for highcharts"
#'
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' @import highcharter
hc_theme_smpl_tailored <- function (...){
  theme <- hc_theme(colors = c("#FF9AA2", "#FFB7B2", "#FFDAC1", "#E2F0CB", "#B5EAD7", "#C7CEEA", "#FBE2FD", "#D9ECFE"),
                    chart = list(style = list(fontFamily = "Roboto", color = "#666666")),
                    title = list(align = "left", style = list(fontFamily = "Roboto Condensed", fontWeight = "bold")),
                    subtitle = list(align = "left", style = list(fontFamily = "Roboto Condensed")),
                    legend = list(align = "right", verticalAlign = "bottom"),
                    xAxis = list(gridLineWidth = 1, gridLineColor = "#F3F3F3", lineColor = "#F3F3F3", minorGridLineColor = "#F3F3F3", tickColor = "#F3F3F3", tickWidth = 1),
                    yAxis = list(gridLineColor = "#F3F3F3", lineColor = "#F3F3F3", minorGridLineColor = "#F3F3F3", tickColor = "#F3F3F3", tickWidth = 1),
                    plotOptions = list(line = list(marker = list(enabled = FALSE)),
                                       spline = list(marker = list(enabled = FALSE)), area = list(marker = list(enabled = FALSE)),
                                       areaspline = list(marker = list(enabled = FALSE)), arearange = list(marker = list(enabled = FALSE)),
                                       bubble = list(maxSize = "10%")))
  theme <- structure(theme, class = "hc_theme")
  if (length(list(...)) > 0) {
    theme <- hc_theme_merge(theme, hc_theme(...))
  }
  theme
}


#' Prepare and return table used for betAS event-wise plots
#'
#' @param eventID
#' @param psitable
#' @param qualtable
#' @param npoints
#' @param colsA
#' @param colsB
#' @param basalColor
#' @param interestColor
#' @param nsim
#' @param seed (boolean) Boolean indicating if seed should be fixed. Default is TRUE
#' @param CoverageWeight (boolean)  Boolean indicating if number of points emitted in beta distributions should be weighted by the coverage of each sample. Default is FALSE
#'
#' @return
#' @export
#'
#' @examples

prepareTableEvent <- function(eventID, psitable, qualtable, npoints, colsA, colsB, labA, labB, basalColor, interestColor, maxDevTable, nsim, seed=TRUE, CoverageWeight=FALSE){



  colsA    <- convertCols(psitable, colsA)
  samplesA <- names(colsA)

  colsB    <- convertCols(psitable, colsB)
  samplesB <- names(colsB)

  x <- match(eventID, qualtable$EVENT)

  #Group betAS (A)
  indBetasA <- individualBetas_nofitting_incr(table = qualtable[x,],
                                              cols = colsA,
                                              indpoints = npoints,
                                              maxdevRefTable = maxDevTable,
                                              seed=seed,
                                              CoverageWeight=CoverageWeight)

  #Group betAS (B)
  indBetasB <- individualBetas_nofitting_incr(table = qualtable[x,],
                                              cols = colsB,
                                              indpoints = npoints,
                                              maxdevRefTable = maxDevTable, seed=seed,
                                              CoverageWeight=CoverageWeight)

  eventPlotsObjs <- list()

  eventPlotsObjs[[1]] <- indBetasA
  eventPlotsObjs[[2]] <- indBetasB

  groupsAB <- c(labA, labB)

  # Get predefined objects A:
  artifA <- indBetasA$BetaPoints
  artiflabelsA <- names(indBetasA$BetaPoints)
  psisA <- indBetasA$PSI
  nA <- length(indBetasA$PSI)
  labelA <- labA
  covA <- indBetasA$inc + indBetasA$exc
  medianA <- indBetasA$MedianBeta

  # Get predefined objects B:
  artifB <- indBetasB$BetaPoints
  artiflabelsB <- names(indBetasB$BetaPoints)
  psisB <- indBetasB$PSI
  nB <- length(indBetasB$PSI)
  labelB <- labB
  covB <- indBetasB$inc + indBetasB$exc
  medianB <- indBetasB$MedianBeta


  if (seed){
    seed <- "21122023"
  } else {
    seed <- paste( sample( 0:9, 8, replace=TRUE ), collapse="" )
  }


  # :::::::::::::::::::::::::
  # 1. Pdiff approach
  # :::::::::::::::::::::::::

  # Number of points in generated Beta: (2 of them - first and last - are to remove to avoid Inf)
  # npoints <- jointpoints + 2
  xlength <- 1001

  indpoints <- length(artifA)/nA

  df_psis <- as.data.frame(rbind(cbind(psisA, names(psisA), rep(labelA, times = length(psisA))),
                                 cbind(psisB, names(psisB), rep(labelB, times = length(psisB)))))
  colnames(df_psis) <- c("PSI", "Sample", "Type")

  jntBetas <- list()

  prob3 <- c()
  seq3 <- seq(-1, 1, length.out = xlength)

  if(length(artifA) != length(artifB)){

    minlen <- min(length(artifA), length(artifB))

  }else{

    minlen <- length(artifB)

  }

  set.seed(seed)
  artifA <- sample(x = artifA, size = minlen)
  artifB <- sample(x = artifB, size = minlen)

  for(l in seq3){
    prob3 <- c(prob3, length(which(artifA - artifB > l))/length(artifA))
  }

  df3 <- as.data.frame(cbind(seq3, prob3))

  p_zero    <- df3$prob3[ceiling(xlength/2)]
  t_p_zero  <- abs(p_zero-0.5)+0.5

  eventPlotsObjs[[3]] <- df3
  eventPlotsObjs[[4]] <- p_zero
  eventPlotsObjs[[5]] <- t_p_zero

  # :::::::::::::::::::::::::
  # 2. F-statistic approach
  # :::::::::::::::::::::::::

  withins <- c()
  betweens <- c()
  pointsA <- indBetasA$BetaPoints
  pointsB <- indBetasB$BetaPoints

  samplesA <- unique(names(pointsA))
  samplesB <- unique(names(pointsB))


  Nbpoints_to_sample <- 1000 # number of differences to sample, will only affect the resolution of the result;
  # If CoverageWeight = T, samples with higher coverage will be more likely to be picked for a difference

  if (length(samplesA)!=1){


    # Matrix determining all possible combinations of sample pairs within group
    within_combs_A <- combn(samplesA, 2, simplify = TRUE)

    for(combA in 1:ncol(within_combs_A)){

      sample1 <- within_combs_A[1,combA]
      sample2 <- within_combs_A[2,combA]

      cov_sample1 <- indBetasA$tinc[sample1] + indBetasA$texc[sample1]
      cov_sample2 <- indBetasA$tinc[sample2] + indBetasA$texc[sample2]
      cov_total <- sum(indBetasA$tinc + indBetasA$texc)

      Prop_groupA <- (indBetasA$tinc+indBetasA$texc)/cov_total

      P12 <- as.numeric(Prop_groupA[sample1]*(Prop_groupA[sample2]/(sum(Prop_groupA[names(Prop_groupA)!=sample1]))) + Prop_groupA[sample2]*(Prop_groupA[sample1]/(sum(Prop_groupA[names(Prop_groupA)!=sample2]))))

      Nb_sampledpoints_12 <- round(Nbpoints_to_sample*P12)

      betas1  <- pointsA[grep(within_combs_A[1,combA], names(pointsA))] # define points of sample 1
      betas2  <- pointsA[grep(within_combs_A[2,combA], names(pointsA))] # define points of sample 2

      set.seed(seed)

      betas1  <- sample(betas1, Nb_sampledpoints_12, replace=T)
      betas2  <- sample(betas2, Nb_sampledpoints_12, replace=T)

      withins <- c(withins, betas2-betas1)
    }

  }



  if (length(samplesB)!=1){
    # Matrix determining all possible combinations of sample pairs within group
    within_combs_B <- combn(samplesB, 2, simplify = TRUE)

    for(combB in 1:ncol(within_combs_B)){

      sample1 <- within_combs_B[1,combB]
      sample2 <- within_combs_B[2,combB]

      cov_sample1 <- indBetasB$tinc[sample1] + indBetasB$texc[sample1]
      cov_sample2 <- indBetasB$tinc[sample2] + indBetasB$texc[sample2]
      cov_total <- sum(indBetasB$tinc + indBetasB$texc)

      Prop_groupB <- (indBetasB$tinc+indBetasB$texc)/cov_total

      P12 <- as.numeric(Prop_groupB[sample1]*(Prop_groupB[sample2]/(sum(Prop_groupB[names(Prop_groupB)!=sample1]))) + Prop_groupB[sample2]*(Prop_groupB[sample1]/(sum(Prop_groupB[names(Prop_groupB)!=sample2]))))

      Nb_sampledpoints_12 <- round(Nbpoints_to_sample*P12)

      betas1  <- pointsB[grep(within_combs_B[1,combB], names(pointsB))] # define points of sample 1
      betas2  <- pointsB[grep(within_combs_B[2,combB], names(pointsB))] # define points of sample 2

      set.seed(seed)
      betas1  <- sample(betas1, Nb_sampledpoints_12, replace=T)
      betas2  <- sample(betas2, Nb_sampledpoints_12, replace=T)

      withins <- c(withins, betas2-betas1)
    }

  }


  # calculate betweens
  set.seed(seed)
  betas1 <- sample(pointsA, Nbpoints_to_sample*2, replace=T)  # Nbpoints_to_sample*2 so that the number of points sampled in the numerator is the same as the denominator
  betas2 <- sample(pointsB, Nbpoints_to_sample*2, replace=T)
  betweens <-   betas2-betas1

  #
  # for (i in 1:Nbpoints_to_sample){
  #
  #   # Add to the withins one difference between samples from group A
  #   if (length(samplesA)!=1) { # if there is only one sample in group A, it shouldn't contribute to the calculation of the withins
  #     sampled_point_1 <- sample(pointsA,1)
  #     sampled_point_2 <- sample(pointsA[names(pointsA)!=names(sampled_point_1)],1)
  #     withins <- c(withins, sampled_point_2-sampled_point_1)
  #   }
  #
  #   if (length(samplesB)!=1) { # if there is only one sample in group B, it shouldn't contribute to the calculation of the withins
  #     # Add to the withins one difference between samples from group B
  #     sampled_point_1 <- sample(pointsB,1)
  #     sampled_point_2 <- sample(pointsB[names(pointsB)!=names(sampled_point_1)],1)
  #     withins <- c(withins, sampled_point_2-sampled_point_1)
  #   }
  #
  #   # Calculate betweens
  #   sampled_point_1 <- sample(pointsA,1)
  #   sampled_point_2 <- sample(pointsB,1)
  #   betweens <- c(betweens,sampled_point_2-sampled_point_1)
  #
  # }

  # Generate F-like statistic
  fstat  <- median(abs(betweens))/median(abs(withins))


  eventPlotsObjs[[6]] <- withins
  eventPlotsObjs[[7]] <- betweens
  eventPlotsObjs[[8]] <- fstat

  # :::::::::::::::::::::::::
  # 3. FDR approach
  # :::::::::::::::::::::::::

  # Under the null hypothesis of no difference in PSI across groups (all samples come from the same distribution)
  originalMedian  <- median(c(indBetasA$PSI, indBetasB$PSI))
  foundPsi        <- medianB - medianA

  simDeltaPsi <- c()

  npoints <- 10000
  #
  #   simulatedA  <- list()
  #   for(sample in unique(artiflabelsA)){
  #
  #     simReads    <- simulate_reads(cov = covA[which(names(covA) == sample)], psi = originalMedian)
  #     set.seed(seed)
  #     sampleDistA <- rbeta(npoints, shape1 = simReads$inc, shape2 = simReads$exc)
  #     simulatedA[[length(simulatedA)+1]] <- sampleDistA
  #
  #   }
  #
  #   simulatedB  <- list()
  #   for(sample in unique(artiflabelsB)){
  #
  #     simReads <- simulate_reads(cov = covB[which(names(covB) == sample)], psi = originalMedian)
  #     set.seed(seed)
  #     sampleDistB <- rbeta(npoints, shape1 = simReads$inc, shape2 = simReads$exc)
  #     simulatedB[[length(simulatedB)+1]] <- sampleDistB
  #
  #   }


  if (CoverageWeight) {

    NbPointsA <- round(nA * npoints * (indBetasA$inc + indBetasA$exc)/sum(covA))

  } else {

    NbPointsA <- rep(npoints, nA)
    names(NbPointsA) <- unique(artiflabelsA)

  }

  simulatedA  <- list()

  for(sample in unique(artiflabelsA)){

    simReads    <- simulate_reads(cov = covA[which(names(covA) == sample)], psi = originalMedian)
    set.seed(seed)
    sampleDistA <- rbeta(NbPointsA[sample], shape1 = simReads$inc, shape2 = simReads$exc)
    simulatedA[[length(simulatedA)+1]] <- sampleDistA

  }

  simulatedB  <- list()


  if (CoverageWeight) {

    NbPointsB <- round(nB * npoints * (indBetasB$inc + indBetasB$exc)/sum(covB))

  } else {

    NbPointsB <- rep(npoints, nB)
    names(NbPointsB) <- unique(artiflabelsB)

  }

  for(sample in unique(artiflabelsB)){

    simReads <- simulate_reads(cov = covB[which(names(covB) == sample)], psi = originalMedian)
    set.seed(seed)
    sampleDistB <- rbeta(NbPointsB[sample], shape1 = simReads$inc, shape2 = simReads$exc)
    simulatedB[[length(simulatedB)+1]] <- sampleDistB

  }


  sampledPointsA  <- lapply(1:length(simulatedA), function(x) sample(x = simulatedA[[x]], size = nsim))
  concatenatedA   <- as.numeric(as.vector(unlist(sampledPointsA)))
  psiA <- as.numeric(as.vector(lapply(1:nsim, function(x) median(concatenatedA[seq(from = x, to = nsim*nA, by = nsim)]))))

  sampledPointsB  <- lapply(1:length(simulatedB), function(x) sample(x = simulatedB[[x]], size = nsim))
  concatenatedB   <- as.numeric(as.vector(unlist(sampledPointsB)))
  psiB <- as.numeric(as.vector(lapply(1:nsim, function(x) median(concatenatedB[seq(from = x, to = nsim*nB, by = nsim)]))))

  simDeltaPsi <- psiB - psiA

  fdr <- length(which(abs(simDeltaPsi) >= abs(foundPsi)))/nsim

  eventPlotsObjs[[9]] <- simulatedA
  eventPlotsObjs[[10]] <- simulatedB
  eventPlotsObjs[[11]] <- simDeltaPsi
  eventPlotsObjs[[12]] <- fdr
  eventPlotsObjs[[13]] <- nsim

  eventPlotsObjs[[14]] <- foundPsi
  eventPlotsObjs[[15]] <- eventID
  eventPlotsObjs[[16]] <- groupsAB

  names(eventPlotsObjs) <- c("indBetasA", "indBetasB",
                             "DataFramePdiff", "Pzero", "Pdiff",
                             "Withins", "Betweens", "Fstat",
                             "simA", "simB", "simDelta", "FDR", "nsim",
                             "deltaPsi", "eventID", "groupsAB")

  return(eventPlotsObjs)

}
# prepareTableEvent <- function(eventID, psitable, qualtable, npoints, colsA, colsB, labA, labB, basalColor, interestColor, maxDevTable, nsim, seed=TRUE, CoverageWeight=FALSE){
#
#
#
#   colsA    <- convertCols(psitable, colsA)
#   samplesA <- names(colsA)
#
#   colsB    <- convertCols(psitable, colsB)
#   samplesB <- names(colsB)
#
#   x <- match(eventID, qualtable$EVENT)
#
#   #Group betAS (A)
#   indBetasA <- individualBetas_nofitting_incr(table = qualtable[x,],
#                                               cols = colsA,
#                                               indpoints = npoints,
#                                               maxdevRefTable = maxDevTable,
#                                               seed=seed,
#                                               CoverageWeight=CoverageWeight)
#
#   #Group betAS (B)
#   indBetasB <- individualBetas_nofitting_incr(table = qualtable[x,],
#                                               cols = colsB,
#                                               indpoints = npoints,
#                                               maxdevRefTable = maxDevTable, seed=seed,
#                                               CoverageWeight=CoverageWeight)
#
#   eventPlotsObjs <- list()
#
#   eventPlotsObjs[[1]] <- indBetasA
#   eventPlotsObjs[[2]] <- indBetasB
#
#   groupsAB <- c(labA, labB)
#
#   # Get predefined objects A:
#   artifA <- indBetasA$BetaPoints
#   artiflabelsA <- names(indBetasA$BetaPoints)
#   psisA <- indBetasA$PSI
#   nA <- length(indBetasA$PSI)
#   labelA <- labA
#   covA <- indBetasA$inc + indBetasA$exc
#   medianA <- indBetasA$MedianBeta
#
#   # Get predefined objects B:
#   artifB <- indBetasB$BetaPoints
#   artiflabelsB <- names(indBetasB$BetaPoints)
#   psisB <- indBetasB$PSI
#   nB <- length(indBetasB$PSI)
#   labelB <- labB
#   covB <- indBetasB$inc + indBetasB$exc
#   medianB <- indBetasB$MedianBeta
#
#
#   if (seed){
#     seed <- "21122023"
#   } else {
#     seed <- paste( sample( 0:9, 8, replace=TRUE ), collapse="" )
#   }
#
#
#   # :::::::::::::::::::::::::
#   # 1. Pdiff approach
#   # :::::::::::::::::::::::::
#
#   # Number of points in generated Beta: (2 of them - first and last - are to remove to avoid Inf)
#   # npoints <- jointpoints + 2
#   xlength <- 1001
#
#   indpoints <- length(artifA)/nA
#
#   df_psis <- as.data.frame(rbind(cbind(psisA, names(psisA), rep(labelA, times = length(psisA))),
#                                  cbind(psisB, names(psisB), rep(labelB, times = length(psisB)))))
#   colnames(df_psis) <- c("PSI", "Sample", "Type")
#
#   jntBetas <- list()
#
#   prob3 <- c()
#   seq3 <- seq(-1, 1, length.out = xlength)
#
#   if(length(artifA) != length(artifB)){
#
#     minlen <- min(length(artifA), length(artifB))
#
#   }else{
#
#     minlen <- length(artifB)
#
#   }
#   set.seed(seed)
#   artifA <- sample(x = artifA, size = minlen)
#   artifB <- sample(x = artifB, size = minlen)
#
#   for(l in seq3){
#     prob3 <- c(prob3, length(which(artifA - artifB > l))/length(artifA))
#   }
#
#   df3 <- as.data.frame(cbind(seq3, prob3))
#
#   p_zero    <- df3$prob3[ceiling(xlength/2)]
#   t_p_zero  <- abs(p_zero-0.5)+0.5
#
#   eventPlotsObjs[[3]] <- df3
#   eventPlotsObjs[[4]] <- p_zero
#   eventPlotsObjs[[5]] <- t_p_zero
#
#   # :::::::::::::::::::::::::
#   # 2. F-statistic approach
#   # :::::::::::::::::::::::::
#
#   #Calculate withins for each group (withins is the vector that will concatenate all within differences found)
#   withins <- c()
#
#   pointsA <- indBetasA$BetaPoints
#   pointsB <- indBetasB$BetaPoints
#
#   samplesA <- unique(names(pointsA))
#   samplesB <- unique(names(pointsB))
#
#   # Matrix determining all possible combinations of sample pairs within group
#   within_combs_A <- combn(samplesA, 2, simplify = TRUE)
#
#   # Define withins for each of A combinations
#   for(combA in 1:ncol(within_combs_A)){
#
#     betas1  <- pointsA[grep(within_combs_A[1,combA], names(pointsA))]
#     betas2  <- pointsA[grep(within_combs_A[2,combA], names(pointsA))]
#     withins <- c(withins, betas2-betas1)
#
#   }
#
#   # Matrix determining all possible combinations of sample pairs within group
#   within_combs_B <- combn(samplesB, 2, simplify = TRUE)
#
#   # Define withins for each of A combinations
#   for(combB in 1:ncol(within_combs_B)){
#
#     betas1  <- pointsB[grep(within_combs_B[1,combB], names(pointsB))]
#     betas2  <- pointsB[grep(within_combs_B[2,combB], names(pointsB))]
#     withins <- c(withins, betas2-betas1)
#
#   }
#
#   #Calculate betweens
#   groups <- c(labA, labB)
#
#   betweens <- c()
#
#   if(length(pointsA) != length(pointsB)){
#
#     minlen <- min(length(pointsA), length(pointsB))
#
#   }else{
#
#     minlen <- length(pointsB)
#
#   }
#   set.seed(seed)
#   pointsB <- sample(pointsB, size = minlen)
#   pointsA <- sample(pointsA, size = minlen)
#
#   betweens <- c(betweens, pointsB-pointsA)
#
#   # Generate F-like statistic
#   fstat  <- median(abs(betweens))/median(abs(withins))
#
#   # Calculate delta PSI
#   deltaPsi 	<- indBetasB$MedianBeta - indBetasA$MedianBeta
#
#   eventPlotsObjs[[6]] <- withins
#   eventPlotsObjs[[7]] <- betweens
#   eventPlotsObjs[[8]] <- fstat
#
#   # :::::::::::::::::::::::::
#   # 3. FDR approach
#   # :::::::::::::::::::::::::
#
#   # Under the null hypothesis of no difference in PSI across groups (all samples come from the same distribution)
#   originalMedian  <- median(c(indBetasA$PSI, indBetasB$PSI))
#   foundPsi        <- medianB - medianA
#
#   simDeltaPsi <- c()
#
#   npoints <- 10000
# #
# #   simulatedA  <- list()
# #   for(sample in unique(artiflabelsA)){
# #
# #     simReads    <- simulate_reads(cov = covA[which(names(covA) == sample)], psi = originalMedian)
# #     set.seed(seed)
# #     sampleDistA <- rbeta(npoints, shape1 = simReads$inc, shape2 = simReads$exc)
# #     simulatedA[[length(simulatedA)+1]] <- sampleDistA
# #
# #   }
# #
# #   simulatedB  <- list()
# #   for(sample in unique(artiflabelsB)){
# #
# #     simReads <- simulate_reads(cov = covB[which(names(covB) == sample)], psi = originalMedian)
# #     set.seed(seed)
# #     sampleDistB <- rbeta(npoints, shape1 = simReads$inc, shape2 = simReads$exc)
# #     simulatedB[[length(simulatedB)+1]] <- sampleDistB
# #
# #   }
#
#
#   if (CoverageWeight) {
#
#     NbPointsA <- round(nA * npoints * (indBetasA$inc + indBetasA$exc)/sum(covA))
#
#   } else {
#
#     NbPointsA <- rep(npoints, nA)
#     names(NbPointsA) <- unique(artiflabelsA)
#
#   }
#
#   simulatedA  <- list()
#
#   for(sample in unique(artiflabelsA)){
#
#     simReads    <- simulate_reads(cov = covA[which(names(covA) == sample)], psi = originalMedian)
#     set.seed(seed)
#     sampleDistA <- rbeta(NbPointsA[sample], shape1 = simReads$inc, shape2 = simReads$exc)
#     simulatedA[[length(simulatedA)+1]] <- sampleDistA
#
#   }
#
#   simulatedB  <- list()
#
#
#   if (CoverageWeight) {
#
#     NbPointsB <- round(nB * npoints * (indBetasB$inc + indBetasB$exc)/sum(covB))
#
#   } else {
#
#     NbPointsB <- rep(npoints, nB)
#     names(NbPointsB) <- unique(artiflabelsB)
#
#   }
#
#   for(sample in unique(artiflabelsB)){
#
#     simReads <- simulate_reads(cov = covB[which(names(covB) == sample)], psi = originalMedian)
#     set.seed(seed)
#     sampleDistB <- rbeta(NbPointsB[sample], shape1 = simReads$inc, shape2 = simReads$exc)
#     simulatedB[[length(simulatedB)+1]] <- sampleDistB
#
#   }
#
#   set.seed(seed)
#   sampledPointsA  <- lapply(1:length(simulatedA), function(x) sample(x = simulatedA[[x]], size = nsim))
#   concatenatedA   <- as.numeric(as.vector(unlist(sampledPointsA)))
#   psiA <- as.numeric(as.vector(lapply(1:nsim, function(x) median(concatenatedA[seq(from = x, to = nsim*nA, by = nsim)]))))
#   set.seed(seed)
#   sampledPointsB  <- lapply(1:length(simulatedB), function(x) sample(x = simulatedB[[x]], size = nsim))
#   concatenatedB   <- as.numeric(as.vector(unlist(sampledPointsB)))
#   psiB <- as.numeric(as.vector(lapply(1:nsim, function(x) median(concatenatedB[seq(from = x, to = nsim*nB, by = nsim)]))))
#
#   simDeltaPsi <- psiB - psiA
#
#   fdr <- length(which(abs(simDeltaPsi) >= abs(foundPsi)))/nsim
#
#   eventPlotsObjs[[9]] <- simulatedA
#   eventPlotsObjs[[10]] <- simulatedB
#   eventPlotsObjs[[11]] <- simDeltaPsi
#   eventPlotsObjs[[12]] <- fdr
#   eventPlotsObjs[[13]] <- nsim
#
#   eventPlotsObjs[[14]] <- foundPsi
#   eventPlotsObjs[[15]] <- eventID
#   eventPlotsObjs[[16]] <- groupsAB
#
#   names(eventPlotsObjs) <- c("indBetasA", "indBetasB",
#                              "DataFramePdiff", "Pzero", "Pdiff",
#                              "Withins", "Betweens", "Fstat",
#                              "simA", "simB", "simDelta", "FDR", "nsim",
#                              "deltaPsi", "eventID", "groupsAB")
#
#   return(eventPlotsObjs)
#
# }

#' Plot individual beta distributions per sample
#'
#' @param eventObjList
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom ggridges geom_density_ridges
plotDensitiesFromEventObjList <- function(eventObjList, colorA, colorB){

  eventID <- eventObjList$eventID
  indBetasA <- eventObjList$indBetasA
  indBetasB <- eventObjList$indBetasB
  groupsAB <- eventObjList$groupsAB

  # Get predefined objects A:
  artifA <- indBetasA$BetaPoints
  artiflabelsA <- names(indBetasA$BetaPoints)
  psisA <- indBetasA$PSI
  nA <- length(indBetasA$PSI)
  labelA <- groupsAB[1]
  covA <- indBetasA$inc + indBetasA$exc
  medianA <- indBetasA$MedianBeta

  # Get predefined objects B:
  artifB <- indBetasB$BetaPoints
  artiflabelsB <- names(indBetasB$BetaPoints)
  psisB <- indBetasB$PSI
  nB <- length(indBetasB$PSI)
  labelB <- groupsAB[2]
  covB <- indBetasB$inc + indBetasB$exc
  medianB <- indBetasB$MedianBeta

  # Extract sample names
  samplesA <- unique(artiflabelsA)
  samplesB <- unique(artiflabelsB)

  themeColors <- c("#89C0AE", "#E69A9C", "#76CAA0", "#EE805B", "#F7CF77", "#81C1D3")
  randomColorsA <- colorRampPalette(c(themeColors[1], themeColors[2]))(4)
  randomColorsB <- colorRampPalette(c(themeColors[3], themeColors[4]))(4)

  # Data-frame containing emitted points, sample names and group names
  densities_df  <- data.frame("points" = numeric(), "samples" = character())
  psis_df       <- data.frame("psis" = numeric(), "samples" = character())

  #prepare table with emitted points
  points_df_A   <- data.frame("points" = as.numeric(as.vector(unlist(artifA))))
  sample_df_A   <- cbind("points" = points_df_A, "samples" = names(artifA))
  densities_df  <- rbind(densities_df, sample_df_A)

  #prepare table with psis
  psis_df_A       <- data.frame("psis" = as.numeric(as.vector(unlist(psisA))))
  samplepsis_df_A <- cbind("psis" = psis_df_A, "samples" = samplesA)
  psis_df <- rbind(psis_df, samplepsis_df_A)

  #prepare table with emitted points
  points_df_B   <- data.frame("points" = as.numeric(as.vector(unlist(artifB))))
  sample_df_B   <- cbind("points" = points_df_B, "samples" = names(artifB))
  densities_df  <- rbind(densities_df, sample_df_B)

  #prepare table with psis
  psis_df_B       <- data.frame("psis" = as.numeric(as.vector(unlist(psisB))))
  samplepsis_df_B <- cbind("psis" = psis_df_B, "samples" = samplesB)
  psis_df <- rbind(psis_df, samplepsis_df_B)


  densities_df$group <- as.character(c(rep(labA, times = nrow(points_df_A)),
                                       rep(labB, times = nrow(points_df_B))))

  psis_df$group <- as.character(c(rep(labA, times = nA),
                                  rep(labB, times = nB)))


  themeColors <- c("#89C0AE", "#E69A9C", "#76CAA0", "#EE805B", "#F7CF77", "#81C1D3")
  colorA <- themeColors[1]
  colorB <- themeColors[6]

  plot <- ggplot(densities_df,
                 aes(x = points,
                     y = samples,
                     group = samples,
                     color = group,
                     fill = group)) +
    geom_density_ridges(size = 2,
                        alpha = 0.8,
                        show.legend = FALSE,
                        rel_min_height = 0.001) +
    geom_rug(data = psis_df,
             aes(x = psis,
                 color = group),
             size = 3,
             alpha = 0.8,
             show.legend = FALSE) +
    scale_color_manual(values = c(colorA, colorB)) +
    scale_fill_manual(values = c(colorA, colorB)) +
    xlab("PSI") +
    ylab("") +
    labs(title ) +
    scale_x_continuous(breaks = seq(0,1, 0.25), limits = c(0,1)) +
    theme_betAS()

  return(plot)

}

#' Plot Pdiff explanation plot per sample
#'
#' @param eventObjList
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom ggridges geom_density_ridges
plotPDiffFromEventObjList <- function(eventObjList){

  eventID <- eventObjList$eventID
  indBetasA <- eventObjList$indBetasA
  indBetasB <- eventObjList$indBetasB
  Pdiff_df <- eventObjList$DataFramePdiff
  pzero <- eventObjList$Pzero
  pdiff <- eventObjList$Pdiff
  groupsAB <- eventObjList$groupsAB

  # Get predefined objects A:
  artifA <- indBetasA$BetaPoints
  artiflabelsA <- names(indBetasA$BetaPoints)
  psisA <- indBetasA$PSI
  nA <- length(indBetasA$PSI)
  labelA <- groupsAB[1]
  covA <- indBetasA$inc + indBetasA$exc
  medianA <- indBetasA$MedianBeta

  # Get predefined objects B:
  artifB <- indBetasB$BetaPoints
  artiflabelsB <- names(indBetasB$BetaPoints)
  psisB <- indBetasB$PSI
  nB <- length(indBetasB$PSI)
  labelB <- groupsAB[2]
  covB <- indBetasB$inc + indBetasB$exc
  medianB <- indBetasB$MedianBeta

  # Extract sample names
  samplesA <- unique(artiflabelsA)
  samplesB <- unique(artiflabelsB)

  title <- paste0("P. diff = ", round(pdiff, digits = 3))

  plot <- ggplot(Pdiff_df,
                 aes(x = seq3,
                     y = prob3)) +
    geom_vline(xintercept = 0,
               linetype = 3,
               size = 1) +
    geom_line(size = 1,
              color = "#DC143C",
              show.legend = FALSE) +
    geom_hline(yintercept = pzero) +
    xlab("x") +
    ylab(paste0("P(", labelA, " - ", labelB, " > x)")) +
    labs(title = title) +
    scale_x_continuous(breaks = seq(-1,1, 0.5), limits = c(-1,1)) +
    theme_betAS() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 24),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA))

  return(plot)

}

#' Plot Fstat explanation plot per sample
#'
#' @param eventObjList
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom ggridges geom_density_ridges
plotFstatFromEventObjList <- function(eventObjList){

  eventID <- eventObjList$eventID
  indBetasA <- eventObjList$indBetasA
  indBetasB <- eventObjList$indBetasB
  withins <- eventObjList$Withins
  betweens <- eventObjList$Betweens
  fstat <- eventObjList$Fstat
  groupsAB <- eventObjList$groupsAB

  # Get predefined objects A:
  artifA <- indBetasA$BetaPoints
  artiflabelsA <- names(indBetasA$BetaPoints)
  psisA <- indBetasA$PSI
  nA <- length(indBetasA$PSI)
  labelA <- groupsAB[1]
  covA <- indBetasA$inc + indBetasA$exc
  medianA <- indBetasA$MedianBeta

  # Get predefined objects B:
  artifB <- indBetasB$BetaPoints
  artiflabelsB <- names(indBetasB$BetaPoints)
  psisB <- indBetasB$PSI
  nB <- length(indBetasB$PSI)
  labelB <- groupsAB[2]
  covB <- indBetasB$inc + indBetasB$exc
  medianB <- indBetasB$MedianBeta

  # Extract sample names
  samplesA <- unique(artiflabelsA)
  samplesB <- unique(artiflabelsB)

  title <- paste0("F-statistic = ", round(fstat, digits = 3))

  relHeight <- max(max(density(abs(withins))$y), max(density(abs(betweens))$y))

  plot <- ggplot() +
    geom_density(aes(x = abs(betweens)),
                 color = NA,
                 fill = "#DC143C",
                 alpha = 0.2) +
    geom_density(aes(x = abs(withins)),
                 color = NA,
                 fill = "darkgray",
                 alpha = 0.3) +
    geom_density(aes(x = abs(betweens)),
                 color = "#DC143C",
                 size = 1,
                 fill = NA) +
    geom_density(aes(x = abs(withins)),
                 color = "darkgray",
                 size = 1,
                 fill = NA) +
    annotate(x = 0.75,
             y = c(0.8, 0.7)*relHeight,
             size = 7,
             geom = "text",
             label = c("within groups", "between groups"),
             color = c("darkgray", "#DC143C")) +
    xlab("Absolute PSI difference") +
    ylab("Density") +
    labs(title = title) +
    scale_x_continuous(breaks = seq(0,1,0.25), limits = c(0,1)) +
    theme_betAS() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 24))

  return(plot)

}


#' Plot FPR explanation plot per sample
#'
#' @param eventObjList
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom ggridges geom_density_ridges
plotFDRFromEventObjList <- function(eventObjList){

  eventID <- eventObjList$eventID
  indBetasA <- eventObjList$indBetasA
  indBetasB <- eventObjList$indBetasB
  simulatedA <- eventObjList$simA
  simulatedB <- eventObjList$simB
  simDeltaPsi <- eventObjList$simDelta
  fdr <- eventObjList$FDR
  nsim <- eventObjList$nsim
  deltaPsi <- eventObjList$deltaPsi
  groupsAB <- eventObjList$groupsAB

  refSeq      <- seq(0,1, 0.1)
  maxDeltaPsi <- abs(deltaPsi)
  pos         <- which(abs(refSeq-maxDeltaPsi)==min(abs(refSeq-maxDeltaPsi)))
  refScale    <- refSeq[pos+1]

  # Get predefined objects A:
  artifA <- indBetasA$BetaPoints
  artiflabelsA <- names(indBetasA$BetaPoints)
  psisA <- indBetasA$PSI
  nA <- length(indBetasA$PSI)
  labelA <- groupsAB[1]
  covA <- indBetasA$inc + indBetasA$exc
  medianA <- indBetasA$MedianBeta

  # Get predefined objects B:
  artifB <- indBetasB$BetaPoints
  artiflabelsB <- names(indBetasB$BetaPoints)
  psisB <- indBetasB$PSI
  nB <- length(indBetasB$PSI)
  labelB <- groupsAB[2]
  covB <- indBetasB$inc + indBetasB$exc
  medianB <- indBetasB$MedianBeta

  # Extract sample names
  samplesA <- unique(artiflabelsA)
  samplesB <- unique(artiflabelsB)

  if(fdr == 0){

    res <- 1/nsim
    title <- paste0("FPR < ", res)

    }else{

    title <- paste0("FPR = ", round(fdr, digits = 3))

  }

  relHeight <- max(density(simDeltaPsi)$y)

  plot <- ggplot() +
    geom_density(aes(x = simDeltaPsi),
                 fill = "gray",
                 alpha = 0.4) +
    geom_vline(xintercept = medianB - medianA,
               color = "#DC143C",
               size = 2) +
    xlab(expression(Delta*PSI[simulated])) +
    annotate(x = deltaPsi,
             y = 0.75*relHeight,
             size = 6,
             geom = "label",
             label = paste0(round(deltaPsi, digits = 3)),
             color = "#DC143C") +
    scale_x_continuous(breaks = seq(-refScale,refScale,length.out = 5), limits = c(-refScale,refScale)) +
    labs(title = title) +
    ylab("Density") +
    theme_betAS() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 24))

  return(plot)

}

#' Prepare and return table used for betAS event-wise plots (multiple-group section)
#'
#' @param eventID
#' @param psitable
#' @param qualtable
#' @param groupList
#' @param npoints
#' @param maxDevTable
#' @param seed (boolean) Boolean indicating if seed should be fixed. Default is TRUE
#' @param CoverageWeight (boolean)  Boolean indicating if number of points emitted in beta distributions should be weighted by the coverage of each sample. Default is FALSE
#'
#' @return
#' @export
#'
#' @examples
prepareTableEventMultiple <- function(eventID, psitable, qualtable, groupList, npoints, maxDevTable, seed=TRUE, CoverageWeight=FALSE){

  x <- match(eventID, qualtable$EVENT)

  # Prepare individual betAS object per group
  groupNames  <- names(groupList)
  # listNames   <- paste0(groupNames, "_indBetas")
  listNames   <-  groupNames

  samplesPerGroupList <- list()
  indBetasList        <- list()

  for(i in 1:length(groupNames)){

    samples <- groupList[[i]]$samples
    pos <- match(samples, colnames(psitable))

    samplesPerGroupList[[i]] <- samples

    # Prepare individual betAS object for this group
    indbetas <- individualBetas_nofitting_incr(table = qualtable[x,],
                                               cols = pos,
                                               indpoints = npoints,
                                               maxdevRefTable = maxDevTable,
                                               seed=seed,
                                               CoverageWeight = CoverageWeight)
    # Name each position in list after the event
    # names(indbetas) <- groupNames[i]

    # Assign ind betas to global list
    indBetasList[[i]] <- indbetas

  }

  names(indBetasList) <- listNames
  names(samplesPerGroupList) <- listNames

  # Define groups

  samples <- c()
  groups  <- c()

  for(i in 1:length(samplesPerGroupList)){

    samples <- c(samples, samplesPerGroupList[[i]])
    groups  <- c(groups, rep(groupNames[i], times = length(samplesPerGroupList[[i]])))

  }


  # Store event-wise results
  allBetaPoints   <- c()
  allMedianPsi  <- c()
  withins     <- c()
  betweens    <- c()
  Nbpoints_to_sample <- 1000

  for(g in 1:length(unique(groups))){

    points <- indBetasList[[g]]$BetaPoints
    medians <- indBetasList[[g]]$MedianBeta
    groupSamples <- samplesPerGroupList[[g]]
    nsamp <- length(groupSamples)

    allBetaPoints <- c(allBetaPoints, points)
    allMedianPsi <- c(allMedianPsi, medians)

    # Matrix determining all possible combinations of sample pairs within group
    within_combs <- combn(groupSamples, 2, simplify = TRUE)

    # Define withins for each of those combinations
    for(c in 1:ncol(within_combs)){

      sample1 <- paste0(within_combs[1,c],".Q")
      sample2 <- paste0(within_combs[2,c],".Q")

      cov_sample1 <- indBetasList[[g]]$tinc[sample1] + indBetasList[[g]]$texc[sample1]
      cov_sample2 <- indBetasList[[g]]$tinc[sample2] + indBetasList[[g]]$texc[sample2]
      cov_total <- sum(indBetasList[[g]]$tinc + indBetasList[[g]]$texc)

      Prop_group <- (indBetasList[[g]]$tinc+indBetasList[[g]]$texc)/cov_total

      P12 <- as.numeric(Prop_group[sample1]*(Prop_group[sample2]/(sum(Prop_group[names(Prop_group)!=sample1]))) + Prop_group[sample2]*(Prop_group[sample1]/(sum(Prop_group[names(Prop_group)!=sample2]))))

      Nb_sampledpoints_12 <- round(Nbpoints_to_sample*P12)

      betas1  <- points[grep(sample1, names(points))] # define points of sample 1
      betas2  <- points[grep(sample2, names(points))] # define points of sample 2

      set.seed(seed)
      betas1  <- sample(betas1, Nb_sampledpoints_12, replace=T)
      betas2  <- sample(betas2, Nb_sampledpoints_12, replace=T)

      withins <- c(withins, betas2-betas1)

    }

  }

  between_combs <- combn(unique(groups), 2, simplify = TRUE)

  for(c in 1:ncol(between_combs)){

    posA <-   grep(between_combs[1,c], names(indBetasList))
    posB <-   grep(between_combs[2,c], names(indBetasList))

    betasA  <- indBetasList[[posA]]$BetaPoints
    betasB  <- indBetasList[[posB]]$BetaPoints

    set.seed(seed)
    betas1 <- sample(betasA, Nbpoints_to_sample*length(unique(groups)), replace=T) # Nbpoints_to_sample*length(unique(groups) so that the number of points sampled in the numerator is the same as the denominator
    betas2 <- sample(betasB, Nbpoints_to_sample*length(unique(groups)), replace=T)
    betweens <-   c(betweens,betas2-betas1)

  }

  # Generate F-like statistic
  fstat  <- median(abs(betweens))/median(abs(withins))

  # Generate Pdiff-like
  len_Fstat     <- min(length(betweens), length(withins))
  set.seed(seed)
  between_sample  <- sample(x = betweens, size = len_Fstat)
  within_sample   <- sample(x = withins, size = len_Fstat)

  p_zero_anova    <- length(which(abs(between_sample) - abs(within_sample) > 0))/len_Fstat
  pdiff_anova     <- abs(p_zero_anova-0.5)+0.5

  # Calculate a sort of delta PSI
  medianBetweens  <- median(betweens)
  deltaAbsolute   <- median(abs(betweens)) - median(abs(withins))


  # Store output as a list
  eventPlotsObjs <- list()

  eventPlotsObjs[[1]] <- p_zero_anova
  eventPlotsObjs[[2]] <- pdiff_anova

  eventPlotsObjs[[3]] <- withins
  eventPlotsObjs[[4]] <- betweens
  eventPlotsObjs[[5]] <- fstat

  eventPlotsObjs[[6]] <- medianBetweens
  eventPlotsObjs[[7]] <- deltaAbsolute

  eventPlotsObjs[[8]] <- eventID

  names(eventPlotsObjs) <- c("Pzero", "Pdiff",
                             "Withins", "Betweens", "Fstat",
                             "medianBetweens", "deltaAbsolute",
                             "eventID")

  return(eventPlotsObjs)

}

#' Plot Fstat explanation plot per sample
#'
#' @param eventObjList
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom ggridges geom_density_ridges
plotFstatFromEventObjListMultiple <- function(eventObjList){

  eventID <- eventObjList$eventID
  withins <- eventObjList$Withins
  betweens <- eventObjList$Betweens
  fstat <- eventObjList$Fstat

  title <- paste0("F-statistic = ", round(fstat, digits = 3))

  relHeight <- max(max(density(abs(withins))$y), max(density(abs(betweens))$y))

  plot <- ggplot() +
    geom_density(aes(x = abs(betweens)),
                 color = NA,
                 fill = "#DC143C",
                 alpha = 0.2) +
    geom_density(aes(x = abs(withins)),
                 color = NA,
                 fill = "darkgray",
                 alpha = 0.3) +
    geom_density(aes(x = abs(betweens)),
                 color = "#DC143C",
                 size = 1,
                 fill = NA) +
    geom_density(aes(x = abs(withins)),
                 color = "darkgray",
                 size = 1,
                 fill = NA) +
    annotate(x = 0.75,
             y = c(0.8, 0.7)*relHeight,
             size = 7,
             geom = "text",
             label = c("within groups", "between groups"),
             color = c("darkgray", "#DC143C")) +
    xlab("Absolute PSI difference") +
    ylab("Density") +
    labs(title = title) +
    scale_x_continuous(breaks = seq(0,1,0.25), limits = c(0,1)) +
    theme_betAS() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 24))

  return(plot)

}
