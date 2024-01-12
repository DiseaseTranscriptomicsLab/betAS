IncExcReads <- data.frame(inc=c(10,100,56,93,910,400),
                          exc=c(95,980,450,11,102,43))
row.names(IncExcReads) <- paste0("Sample_",c(1:6))

nA <- 3
nB <- 3
indpoints_aux <- 500

groupA_samples <- paste0("Sample_",c(1:3))
groupB_samples <- paste0("Sample_",c(4:6))


tot_coveraveA <- sum(IncExcReads[groupA_samples,]$inc + IncExcReads[groupA_samples,]$exc) # number
NbPointsA <- round(nA * indpoints_aux * (IncExcReads[groupA_samples,]$inc + IncExcReads[groupA_samples,]$exc)/tot_coveraveA) # vector, each entry is the number of points per sample
names(NbPointsA) <- row.names(IncExcReads[groupA_samples,])
# NbPointsA <- rep(indpoints_aux, length(IncExcReads$inc))
dists <- list()
for (sample in groupA_samples){
  dists[[sample]] <-   rbeta(NbPointsA[sample],
                            shape1 = IncExcReads[sample,1],
                            shape2 = IncExcReads[sample,2])
}

tot_coveraveB <- sum(IncExcReads[groupB_samples,]$inc + IncExcReads[groupB_samples,]$exc) # number
NbPointsB <- round(nB * indpoints_aux * (IncExcReads[groupB_samples,]$inc + IncExcReads[groupB_samples,]$exc)/tot_coveraveB) # vector, each entry is the number of points per sample
names(NbPointsB) <- row.names(IncExcReads[groupB_samples,])
# NbPointsA <- rep(indpoints_aux, length(IncExcReads$inc))

for (sample in groupB_samples){
  dists[[sample]] <-   rbeta(NbPointsB[sample],
                             shape1 = IncExcReads[sample,1],
                             shape2 = IncExcReads[sample,2])
}


ggplot_densities <- melt(dists)
colnames(ggplot_densities) <- c("value","sample")
ggplot_densities <- merge(ggplot_densities, data.frame(sample=paste0("Sample_",c(1:6)),
                                                       group=c(rep("groupA",3), rep("groupB",3))))
ggplot(ggplot_densities, aes(x=value, y=sample, color=group, fill=group))+
  geom_density_ridges(size = 2,
                      alpha = 0.6,
                      show.legend = FALSE,
                      rel_min_height = 0.001) +
  xlab("Proportion spliced-in (PSI)") +
  ylab("") +
  labs(title ) +
  scale_x_continuous(breaks = seq(0,1, 0.1), limits = c(0,1)) +
  theme_betAS() +
  facet_grid(group ~ .,
             scales = "free_y")


# withins
dists_groupA <- dists[groupA_samples]
Nbpoints_to_sample <- sum( unlist(lapply(dists_groupA,length)))
pts_groupA <- melt(dists_groupA)
pts_groupA_names <- pts_groupA[,2]
pts_groupA <- pts_groupA[,1]
names(pts_groupA) <- pts_groupA_names
pts_groupA <- sample(pts_groupA)

frequency_pairs <- data.frame(x=NULL,y=NULL)

for (i in 1:Nbpoints_to_sample){

  sampled_point_1 <- sample(pts_groupA,1)
  sampled_point_2 <- sample(pts_groupA[names(pts_groupA)!=names(sampled_point_1)],1)

  frequency_pairs <- rbind(frequency_pairs, c(names(sampled_point_1),names(sampled_point_2)))

}


frequency_pairs <- apply(frequency_pairs, 1, function(row) paste(sort(row), collapse = "/"))
frequency_pairs <- table(frequency_pairs)

cov_groupA <- (IncExcReads[samples_groupA,1]+IncExcReads[samples_groupA,2])
names(cov_groupA) <- samples_groupA


props_cov_groupA_combs <- combn(samples_groupA, 2, simplify = TRUE)

# Calculate the mean for each pair
pair_means <- apply(props_cov_groupA_combs, 2, function(pair) {
  mean(c(cov_groupA[pair[1]], cov_groupA[pair[2]]))
})
names(pair_means) <- apply(props_cov_groupA_combs, 2, function(pair) paste(pair, collapse = "/"))





#dataset <- getDataset(pathTables=NULL, tool="vast-tools")
#dataset <- getEvents(dataset, tool = "vast-tools")
#dataset_filtered <- filterEvents(dataset, types=c("C1", "C2", "C3", "S", "MIC"), N=10)
#dataset_filtered <- alternativeEvents(dataset_filtered, minPsi=1, maxPsi=99)
# groupList <- list()
# for(i in 1:length(groups)){
#   groupNames <- samples[which(metadata[,groupingVariable] == groups[i])]
#   currentNames <- names(groupList)
#   groupList[[length(groupList)+1]] <- list(name = groups[i],
#                                            samples = groupNames,
#                                            color = random_colors[i])
#   names(groupList) <- make.unique(c(currentNames, groups[i]))
# }
# groupA    <- "Neuron"
# groupB    <- "ESC"
# samplesA    <- groupList[[groupA]]$samples
# samplesB    <- groupList[[groupB]]$samples
# colsGroupA    <- convertCols(dataset_filtered$PSI, samplesA)
# colsGroupB    <- convertCols(dataset_filtered$PSI, samplesB)
# data("maxDevSimulationN100")

eventPos <- 1
psitable <- dataset_filtered$PSI
qualtable <- dataset_filtered$Qual
npoints <-500
colsA <- colsGroupA
colsB <- colsGroupB
labA <- groupA
labB <- groupB
maxDevTable <- maxDevSimulationN100
seed <- T
CoverageWeight <- T

fStatistic_2Groups_PerEvent <- function(eventPos, psitable, qualtable, npoints, colsA, colsB, labA, labB, maxDevTable, seed=TRUE,  CoverageWeight=FALSE){


  if (seed){
    seed <- "21122023"
  } else {
    seed <- paste( sample( 0:9, 8, replace=TRUE ), collapse="" )
  }


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
  betas1 <- sample(pointsA, Nbpoints_to_sample, replace=T)
  betas2 <- sample(pointsB, Nbpoints_to_sample, replace=T)
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
  betas1 <- sample(pointsA, Nbpoints_to_sample, replace=T)
  betas2 <- sample(pointsB, Nbpoints_to_sample, replace=T)
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
