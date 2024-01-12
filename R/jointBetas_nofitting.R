#' Joint betAS (vast-tools)
#' Takes two indBetas objects and performs differential alternative splicing analysis considering the beta distributions between groups A and B. To be applied per row of Qual table,
#' i.e. per element of lists with indBetas event-wise results (list names = events IDs) using pblapply() or similar.
#' @param indBetasA (indBetas) individual betAS object of group A samples
#' @param indBetasB (indBetas) individual betAS object of group B samples
#' @param groupsAB (vector) vector providing labels for groups A and B (in this order)
#'
#' @return jntBetas object list with: 1) P(x = 0) and 2) Flipped P(x = 0), i.e. Pdiff
#' @export
#'
#' @examples
jointBetas_nofitting <- function(indBetasA, indBetasB, groupsAB){

  # Number of points in generated Beta: (2 of them - first and last - are to remove to avoid Inf)
  # npoints <- jointpoints + 2
  xlength <- 1001

  # Get predefined objects A:
  artifA <- indBetasA$BetaPoints
  artiflabelsA <- names(indBetasA$BetaPoints)
  psisA <- indBetasA$PSI
  nA <- length(indBetasA$PSI)
  labelA <- groupsAB[1]

  # Get predefined objects B:
  artifB <- indBetasB$BetaPoints
  artiflabelsB <- names(indBetasB$BetaPoints)
  psisB <- indBetasB$PSI
  nB <- length(indBetasB$PSI)
  labelB <- groupsAB[2]

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

  set.seed("21122023")

  artifA <- sample(x = artifA, size = minlen)
  artifB <- sample(x = artifB, size = minlen)

  for(l in seq3){
    prob3 <- c(prob3, length(which(artifA - artifB > l))/length(artifA))
  }

  df3 <- as.data.frame(cbind(seq3, prob3))

  p_zero <- df3$prob3[ceiling(xlength/2)]
  t_p_zero <- abs(p_zero-0.5)+0.5

  jntBetas[[1]] <- p_zero
  jntBetas[[2]] <- t_p_zero

  names(jntBetas) <- c("P(x=0)", "FlippedP(x=0)")
  return(jntBetas)

}

#' Joint betAS (vast-tools) - fast version
#' Takes two indBetas objects and performs differential alternative splicing analysis considering the beta distributions between groups A and B. To be applied per row of Qual table,
#' i.e. per element of lists with indBetas event-wise results (list names = events IDs) using pblapply() or similar.
#' @param indBetasA (indBetas) individual betAS object of group A samples
#' @param indBetasB (indBetas) individual betAS object of group B samples
#' @param groupsAB (vector) vector providing labels for groups A and B (in this order)
#'
#' @return jntBetas object list with: 1) P(x = 0) and 2) Flipped P(x = 0), i.e. Pdiff
#' @export
#'
#' @examples
jointBetas_nofitting_Fast0 <- function(indBetasA, indBetasB, groupsAB){

  # Get predefined objects A:
  artifA <- indBetasA$BetaPoints
  artiflabelsA <- names(indBetasA$BetaPoints)
  psisA <- indBetasA$PSI
  nA <- length(indBetasA$PSI)
  labelA <- groupsAB[1]

  # Get predefined objects B:
  artifB <- indBetasB$BetaPoints
  artiflabelsB <- names(indBetasB$BetaPoints)
  psisB <- indBetasB$PSI
  nB <- length(indBetasB$PSI)
  labelB <- groupsAB[2]

  indpoints <- length(artifA)/nA

  df_psis <- as.data.frame(rbind(cbind(psisA, names(psisA), rep(labelA, times = length(psisA))),
                                 cbind(psisB, names(psisB), rep(labelB, times = length(psisB)))))
  colnames(df_psis) <- c("PSI", "Sample", "Type")

  jntBetas <- list()

  if(length(artifA) != length(artifB)){

    minlen <- min(length(artifA), length(artifB))

    }else{

    minlen <- length(artifB)

  }

  set.seed("21122023")

  artifA <- sample(x = artifA, size = minlen)
  artifB <- sample(x = artifB, size = minlen)

  p_zero <- length(which(artifA - artifB > 0))/length(artifA)
  t_p_zero <- abs(p_zero-0.5)+0.5

  jntBetas[[1]] <- p_zero
  jntBetas[[2]] <- t_p_zero

  names(jntBetas) <- c("P(x=0)", "FlippedP(x=0)")
  return(jntBetas)

}


#' Joint betAS (vast-tools) - F-statistic
#' Takes samples from two groups, generates two indBetas objects and estimates the F-statistic using the ration between differences between beta distributions of groups A and B and within each group. To be applied per row of Qual table,
#' using pblapply() or similar.
#' @param eventPos
#' @param psitable
#' @param qualtable
#' @param npoints
#' @param colsA
#' @param colsB
#' @param labA
#' @param labB
#' @param seed (boolean) Boolean indicating if seed should be fixed. Default is TRUE
#' @param CoverageWeight (boolean)  Boolean indicating if number of points emitted in beta distributions should be weighted by the coverage of each sample. Default is FALSE
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom ggrepel geom_text_repel
#' @importFrom utils combn
#'
#'
fStatistic_2Groups_PerEvent <- function(eventPos, psitable, qualtable, npoints, colsA, colsB, labA, labB, maxDevTable, seed=TRUE,  CoverageWeight=FALSE){




  fStat_2Groups <- list()

  #Group betAS (A)
  indBetasA <- individualBetas_nofitting_incr(table = qualtable[eventPos,],
                                              cols = colsA,
                                              indpoints = npoints,
                                              maxdevRefTable = maxDevTable, seed=seed, CoverageWeight=CoverageWeight)

  #Group betAS (B)
  indBetasB <- individualBetas_nofitting_incr(table = qualtable[eventPos,],
                                              cols = colsB,
                                              indpoints = npoints,
                                              maxdevRefTable = maxDevTable,seed=seed,CoverageWeight=CoverageWeight)


  if (seed){
    seed <- "21122023"
  } else {
    seed <- paste( sample( 0:9, 8, replace=TRUE ), collapse="" )
  }


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
  betas1 <- sample(pointsA, Nbpoints_to_sample*2, replace=T) # Nbpoints_to_sample*2 so that the number of points sampled in the numerator is the same as the denominator
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

  # Calculate delta PSI
  deltaPsi 	<- indBetasB$MedianBeta - indBetasA$MedianBeta

  fStat_2Groups[[1]] <- fstat
  fStat_2Groups[[2]] <- deltaPsi

  names(fStat_2Groups) <- c("Fstat", "deltaPsi")
  class(fStat_2Groups) <- c(class(fStat_2Groups), "fStat_2Groups")
  return(fStat_2Groups)

}


# fStatistic_2Groups_PerEvent <- function(eventPos, psitable, qualtable, npoints, colsA, colsB, labA, labB, maxDevTable, seed=TRUE,  CoverageWeight=FALSE){
#
#   fStat_2Groups <- list()
#
#   #Group betAS (A)
#   indBetasA <- individualBetas_nofitting_incr(table = qualtable[eventPos,],
#                                               cols = colsA,
#                                               indpoints = npoints,
#                                               maxdevRefTable = maxDevTable, seed=seed, CoverageWeight=CoverageWeight)
#
#   #Group betAS (B)
#   indBetasB <- individualBetas_nofitting_incr(table = qualtable[eventPos,],
#                                               cols = colsB,
#                                               indpoints = npoints,
#                                               maxdevRefTable = maxDevTable,seed=seed,CoverageWeight=CoverageWeight)
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
#
#   set.seed("21122023")
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
#   fStat_2Groups[[1]] <- fstat
#   fStat_2Groups[[2]] <- deltaPsi
#
#   names(fStat_2Groups) <- c("Fstat", "deltaPsi")
#   class(fStat_2Groups) <- c(class(fStat_2Groups), "fStat_2Groups")
#   return(fStat_2Groups)
#
# }

#' Estimate FDR associated with the obtained delta PSI (vast-tools)
#' Takes samples from two groups, generates two indBetas objects and estimates the FDR of the obtained delta PSI by simulating nsim delta PSI under the null hypothesis that all samples have the same original PSI.
#' @param indBetasA (indBetas) individual betAS object of group A samples
#' @param indBetasB (indBetas) individual betAS object of group B samples
#' @param groupsAB (vector) vector providing labels for groups A and B (in this order)
#' @param nsim (numeric) number of simulations
#' @param seed (boolean) Boolean indicating if seed should be fixed. Default is TRUE
#' @param CoverageWeight (boolean)  Boolean indicating if number of points emitted in beta distributions should be weighted by the coverage of each sample. Default is FALSE
#'
#' @return jntBetas object list with: 1) vector of simulated delta PSIs; 2) original median PSI found, 3) FDR and 4) plot with the obtained FDR compared to the simulated
#' @export
#'
#' @examples
#' @importFrom stats density median rbeta sd
estimateFDR <- function(indBetasA, indBetasB, groupsAB, nsim, seed=TRUE, CoverageWeight=FALSE){

if (seed){
  seed <- "21122023"
} else {
  seed <- paste( sample( 0:9, 8, replace=TRUE ), collapse="" )
}

jntBetas <- list()

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

# Under the null hypothesis of no difference in PSI across groups (all samples come from the same distribution)
originalMedian  <- median(c(indBetasA$PSI, indBetasB$PSI))
foundPsi        <- medianB - medianA

simDeltaPsi <- c()

npoints <- 10000

simulatedA  <- list()


if (CoverageWeight) {

  NbPointsA <- round(nA * npoints * (indBetasA$inc + indBetasA$exc)/sum(covA))

} else {

  NbPointsA <- rep(npoints, nA)
  names(NbPointsA) <- unique(artiflabelsA)

}


for(sample in unique(artiflabelsA)){

  simReads    <- simulate_reads(cov = covA[which(names(covA) == sample)], psi = originalMedian)
  set.seed(seed)
  sampleDistA <- rbeta(NbPointsA[sample], shape1 = simReads$inc, shape2 = simReads$exc)
  simulatedA[[length(simulatedA)+1]] <- sampleDistA

}

simulatedB  <- list()


if (CoverageWeight) {

  NbPointsB <- round(nB * npoints * (indBetasB$exc + indBetasB$exc)/sum(covB))

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

set.seed(seed)

sampledPointsA  <- lapply(1:length(simulatedA), function(x) sample(x = simulatedA[[x]], size = nsim))
concatenatedA   <- as.numeric(as.vector(unlist(sampledPointsA)))
psiA <- as.numeric(as.vector(lapply(1:nsim, function(x) median(concatenatedA[seq(from = x, to = nsim*nA, by = nsim)]))))

set.seed(seed)

sampledPointsB  <- lapply(1:length(simulatedB), function(x) sample(x = simulatedB[[x]], size = nsim))
concatenatedB   <- as.numeric(as.vector(unlist(sampledPointsB)))
psiB <- as.numeric(as.vector(lapply(1:nsim, function(x) median(concatenatedB[seq(from = x, to = nsim*nB, by = nsim)]))))

simDeltaPsi <- psiB - psiA

fdr <- length(which(abs(simDeltaPsi) >= abs(foundPsi)))/nsim

jntBetas[[1]] <- simDeltaPsi
jntBetas[[2]] <- foundPsi
jntBetas[[3]] <- fdr

names(jntBetas) <- c("Simulated dPSI", "Found dPSI", "FDR")
return(jntBetas)

}

#' Estimate multiple-groups differential alternative splicing statistics
#' Takes group-wise beta distributions and estimates the emitted values differences between and within groups, from which statistics of differential alternative splicing across groups (Pdiff and F-statistic).
#' @param eventPos (numeric) position of event under analysis in provided list
#' @param indList (list) list of individual betAS objects per group
#' @param samplesList (list) list of vectors of samples corresponding to each group
#' @param groupNames (character vector) vector of group names in indList
#'
#' @return dirGroupBetas object list with: 1) samples considered; 2) group per sample, 3) F-stat, 4) P(|betweens| > |withins|), 5) Transformed Pdiff, 6) Median of between distances and 7) Difference between median(|betweens|) to median(|withins|)
#' @export
#'
#' @examples
#' @importFrom stats density median rbeta sd
generalisedGroupsBetas <- function(eventPos, indList, samplesList, groupNames, seed=TRUE, CoverageWeight=FALSE){

  dirGroupBetas <- list()

  if (seed){
    seed <- "21122023"
  } else {
    seed <- paste( sample( 0:9, 8, replace=TRUE ), collapse="" )
  }

  # indList is a list of objects coming from "individualBetas" function
  # groupNames is a vector of group names, each corresponding to one object in indList
  # samplesList is a list where each element is a vector of sample names,
  # 		that correspond to the samples in the same position in indBetas
  # groupColors are the colors to use in the plot of each group of samples

  samples <- c()
  groups  <- c()

  for(i in 1:length(samplesList)){

    samples <- c(samples, samplesList[[i]])
    groups 	<- c(groups, rep(groupNames[i], times = length(samplesList[[i]])))

  }

  # Assign samples and groups as part of function's output
  dirGroupBetas[[1]] <- samples
  dirGroupBetas[[2]] <- groups


  # Store event-wise results
  allBetaPoints 	<- c()
  allMedianPsi 	<- c()
  withins 		<- c()
  betweens 		<- c()

  for(g in 1:length(unique(groups))){

    points <- indList[[g]][[eventPos]]$BetaPoints
    medians <- indList[[g]][[eventPos]]$MedianBeta
    groupSamples <- samplesList[[g]]
    nsamp <- length(groupSamples)

    allBetaPoints <- c(allBetaPoints, points)
    allMedianPsi <- c(allMedianPsi, medians)

    # Matrix determining all possible combinations of sample pairs within group
    within_combs <- combn(groupSamples, 2, simplify = TRUE)

    Nbpoints_to_sample <- 1000

    # Define withins for each of those combinations
    for(c in 1:ncol(within_combs)){


      if (nsamp!=1){

        sample1 <- paste0(within_combs[1,c],".Q")
        sample2 <- paste0(within_combs[2,c],".Q")

        cov_sample1 <- indList[[g]][[eventPos]]$tinc[sample1] + indList[[g]][[eventPos]]$texc[sample1]
        cov_sample2 <- indList[[g]][[eventPos]]$tinc[sample2] + indList[[g]][[eventPos]]$texc[sample2]
        cov_total <- sum(indList[[g]][[eventPos]]$tinc + indList[[g]][[eventPos]]$texc)

        Prop_group <- (indList[[g]][[eventPos]]$tinc+indList[[g]][[eventPos]]$texc)/cov_total

        P12 <- as.numeric(Prop_group[sample1]*(Prop_group[sample2]/(sum(Prop_group[names(Prop_group)!=sample1]))) + Prop_group[sample2]*(Prop_group[sample1]/(sum(Prop_group[names(Prop_group)!=sample2]))))

        Nb_sampledpoints_12 <- round(Nbpoints_to_sample*P12)

        betas1  <- points[grep(within_combs[1,c], names(points))] # define points of sample 1
        betas2  <- points[grep(within_combs[2,c], names(points))] # define points of sample 2

        set.seed(seed)
        betas1  <- sample(betas1, Nb_sampledpoints_12, replace=T)
        betas2  <- sample(betas2, Nb_sampledpoints_12, replace=T)

        withins <- c(withins, betas2-betas1)

      }

    }

  }

  between_combs <- combn(unique(groups), 2, simplify = TRUE)

  for(c in 1:ncol(between_combs)){

    posA <-   grep(between_combs[1,c], names(indList))
    posB <-   grep(between_combs[2,c], names(indList))

    betasA  <- indList[[posA]][[eventPos]]$BetaPoints
    betasB  <- indList[[posB]][[eventPos]]$BetaPoints

    set.seed(seed)
    betas1 <- sample(betasA, Nbpoints_to_sample*length(unique(groups)), replace=T) # Nbpoints_to_sample*length(unique(groups) so that the number of points sampled in the numerator is the same as the denominator
    betas2 <- sample(betasB, Nbpoints_to_sample*length(unique(groups)), replace=T)
    betweens <-   c(betweens,betas2-betas1)

  }

  # Generate F-like statistic
  fstat  <- median(abs(betweens))/median(abs(withins))

  set.seed(seed)

  # Generate Pdiff-like
  len_Fstat 		<- min(length(betweens), length(withins))
  between_sample  <- sample(x = betweens, size = len_Fstat)
  within_sample   <- sample(x = withins, size = len_Fstat)

  p_zero_anova  	<- length(which(abs(between_sample) - abs(within_sample) > 0))/len_Fstat
  pdiff_anova   	<- abs(p_zero_anova-0.5)+0.5

  # Calculate a sort of delta PSI
  medianBetweens 	<- median(betweens)
  deltaAbsolute 	<- median(abs(betweens)) - median(abs(withins))

  dirGroupBetas[[3]] <- fstat
  dirGroupBetas[[4]] <- p_zero_anova
  dirGroupBetas[[5]] <- pdiff_anova
  dirGroupBetas[[6]] <- medianBetweens
  dirGroupBetas[[7]] <- deltaAbsolute

  names(dirGroupBetas) <- c("samples", "groups",
                            "Fstat", "Pzero", "Pdiff", "medianBetweens", "deltaAbsolute")
  class(dirGroupBetas) <- c(class(dirGroupBetas), "dirGroups")
  return(dirGroupBetas)

}
