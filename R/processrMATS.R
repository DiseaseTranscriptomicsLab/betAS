# Format *.MATS.JC.txt table (rMATS) for quantified PSIs, for further analyses
# @param incTable rMATS' *.MATS.JC.txt table
#
# @return List with: 1) filtered table PSI columns, 2) filtered table "Qual" columns including inc and exc, 3) table with number of events per type and 4) Samples
# @export
#
# @examples

getrMATS <- function(incTable){


  filterRM <- list()

  colnames <- colnames(incTable)

  # Inspect rMATS colnames to identify type of event:
  # Each event type has its own separate file, and its own set of columns, as described in https://github.com/Xinglab/rmats-turbo/blob/v4.1.2/README.md

  if ("exonStart_0base" %in% colnames){
    eventType <- "EX"
    eventCols <- c("exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE")
    # The inclusion form includes the target exon (exonStart_0base, exonEnd)
    eventCoordinates <- paste0(incTable$chr,":",incTable$exonStart_0base,"-",incTable$exonEnd)

    # NOTE for "MXE" events:
    # reading tables using read.table(file, sep="\t", header=TRUE, quote=""), which is the case in the app, corrects colnames starting with numbers by adding an X
  } else if ("X1stExonStart_0base" %in% colnames){
    eventType <- "MXE"
    eventCols <- c("X1stExonStart_0base", "X1stExonEnd", "X2ndExonStart_0base", "X2ndExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE")
    # If the strand is +, then the inclusion form includes the 1st exon (1stExonStart_0base, 1stExonEnd) and skips the 2nd exon
    # If the strand is -, then the inclusion form includes the 2nd exon (2ndExonStart_0base, 2ndExonEnd) and skips the 1st exon
    eventCoordinates <- paste0(incTable$chr,":(1st)",incTable$X1stExonStart_0base,"-",incTable$X1stExonEnd,";(2nd)",incTable$X2ndExonStart_0base,"-",incTable$X2ndExonEnd)

  } else if ("longExonStart_0base"  %in% colnames){
    eventType <- "Altss"
    eventCols <- c("longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")
    # The inclusion form includes the long exon (longExonStart_0base, longExonEnd) instead of the short exon (shortES shortEE)
    eventCoordinates <- paste0(incTable$chr,":",incTable$longExonStart_0base,"-",incTable$longExonEnd)

  } else if ("riExonStart_0base"  %in% colnames){
    eventType <- "IR"
    eventCols <- c("riExonStart_0base", "riExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE")
    # The inclusion form includes (retains) the intron (upstreamEE, downstreamES)
    eventCoordinates <- paste0(incTable$chr,":",incTable$upstreamEE,"-",incTable$downstreamES)
  } else {
    print("The provided file is not supported")
  }

  # Mimicking vast-tools INCLUSION table structure (6 columns for event ID) to facilitate the compatibility with other betAS functions
  # "GENE" | "EVENT" | "COORD" | "LENGTH" | "FullCO" | "COMPLEX"
  commonCols <- cbind(incTable[,c("geneSymbol","ID")], eventCoordinates, rep(0,nrow(incTable)), eventCoordinates, rep(eventType,nrow(incTable)))

  # Remove "\"" from some gene symbols (probably not needed if read.table/read.delim is done with quote = "\"")
  commonCols$geneSymbol <- gsub(pattern = "\"", replacement = "", x = commonCols$geneSymbol)

  # Check if rMATS table was done with one or two groups
  if(all(is.na(c(incTable$IncLevel2, incTable$IJC_SAMPLE_2, incTable$SJC_SAMPLE_2)))){

    # If Group2 does not exist in rMATS table:
    # ::::::::::::::::::::::::::::::::::::::::

    # Infer number of samples from the first row in "SJC columns" by summing 1 to the number of ","
    Nsamples_Group1 <- length(gregexpr(",", incTable$SJC_SAMPLE_1[1], fixed = TRUE)[[1]])+1

    # Name samples from groups 1 and 2
    Samples_Group1 <- paste0("G1_S",1:Nsamples_Group1)

    # psiTable
    psiRM <- cbind(commonCols,
                   apply(matrix(unlist(strsplit(incTable$IncLevel1, ",")),ncol=Nsamples_Group1,byrow=T), 2, as.numeric))
    colnames(psiRM) <- c("GENE","EVENT","COORD","LENGTH","FullCO","COMPLEX",Samples_Group1)

    # qualTable
    inc <-  cbind(apply(matrix(unlist(strsplit(incTable$IJC_SAMPLE_1, ",")),ncol=Nsamples_Group1,byrow=T), 2, as.numeric)/incTable$IncFormLen)
    exc <-  cbind(apply(matrix(unlist(strsplit(incTable$SJC_SAMPLE_1, ",")),ncol=Nsamples_Group1,byrow=T), 2, as.numeric)/incTable$SkipFormLen)

    qualRM <- cbind(commonCols,
                    # Mimicking vast-tools INCLUSION table ".Q" columns to facilitate the compatibility with other betAS functions
                    matrix(paste0("A,A,0=0=0,A,",rep(eventType,nrow(incTable)),"@", inc , ",", exc ), nrow = nrow(inc)))
    colnames(qualRM) <- c("GENE","EVENT","COORD","LENGTH","FullCO","COMPLEX",paste0(Samples_Group1,".Q"))

  }else{

    # If Group2 exists in rMATS table:
    # ::::::::::::::::::::::::::::::::

    # Infer number of samples from the first row in "SJC columns" by summing 1 to the number of ","
    Nsamples_Group1 <- length(gregexpr(",", incTable$SJC_SAMPLE_1[1], fixed = TRUE)[[1]])+1
    Nsamples_Group2 <- length(gregexpr(",", incTable$SJC_SAMPLE_2[1], fixed = TRUE)[[1]])+1

    # Name samples from groups 1 and 2
    Samples_Group1 <- paste0("G1_S",1:Nsamples_Group1)
    Samples_Group2 <- paste0("G2_S",1:Nsamples_Group2)

    # psiTable
    psiRM <- cbind(commonCols,
                   apply(matrix(unlist(strsplit(incTable$IncLevel1, ",")),ncol=Nsamples_Group1,byrow=T), 2, as.numeric)*100,
                   apply(matrix(unlist(strsplit(incTable$IncLevel2, ",")),ncol=Nsamples_Group2,byrow=T), 2, as.numeric)*100)
    colnames(psiRM) <- c("GENE","EVENT","COORD","LENGTH","FullCO","COMPLEX",Samples_Group1,Samples_Group2)

    # qualTable
    inc <-  cbind(apply(matrix(unlist(strsplit(incTable$IJC_SAMPLE_1, ",")),ncol=Nsamples_Group1,byrow=T), 2, as.numeric)/incTable$IncFormLen,
                  apply(matrix(unlist(strsplit(incTable$IJC_SAMPLE_2, ",")),ncol=Nsamples_Group1,byrow=T), 2, as.numeric)/incTable$IncFormLen)
    exc <-  cbind(apply(matrix(unlist(strsplit(incTable$SJC_SAMPLE_1, ",")),ncol=Nsamples_Group1,byrow=T), 2, as.numeric)/incTable$SkipFormLen,
                  apply(matrix(unlist(strsplit(incTable$SJC_SAMPLE_2, ",")),ncol=Nsamples_Group1,byrow=T), 2, as.numeric)/incTable$SkipFormLen)

    qualRM <- cbind(commonCols,
                    # Mimicking vast-tools INCLUSION table ".Q" columns to facilitate the compatibility with other betAS functions
                    matrix(paste0("A,A,0=0=0,A,",rep(eventType,nrow(incTable)),"@", inc , ",", exc ), nrow = nrow(inc)))
    colnames(qualRM) <- c("GENE","EVENT","COORD","LENGTH","FullCO","COMPLEX",paste0(Samples_Group1,".Q"),paste0(Samples_Group2,".Q"))

  }


  filterRM[[1]] <- psiRM
  filterRM[[2]] <- qualRM
  filterRM[[3]] <- table(psiRM$COMPLEX)
  filterRM[[4]] <- colnames(psiRM)[-c(1:6)]

  names(filterRM) <- c("PSI", "Qual", "EventsPerType", "Samples")
  return(filterRM)



}


# Filter *.MATS.JC.txt table (rMATS) for quantified PSIs
# Filter table from rMATS to remove events containing NAs in at least one sample and split PSI and Qual tables
# Normalised inc and exc reads are obtained from "IJC" and "SJC" columns dividing junction counts by IncFormLen and SkipFormLen
# @param RMlist list containing PSI and Qual tables, as well as event and samples, obtained with getVastTools()
#
# @return List with: 1) filtered table PSI columns, 2) filtered table "Qual" columns including inc and exc, 3) table with number of events per type and 4) Samples
# @export
#
# @examples


filterrMATS <- function(RMlist){

  filterRM <- RMlist


  psiRM <- filterRM$PSI
  qualRM <- filterRM$Qual

  # Remove events containing at least one NA
  psiRM$AnyNA    <- apply(psiRM, 1, anyNA)
  psiRM          <- psiRM[which(psiRM$AnyNA == FALSE),]
  psiRM          <- psiRM[,-c(ncol(psiRM))]
  qualRM         <- qualRM[match(psiRM$EVENT, qualRM$EVENT),]

  filterRM[[1]] <- psiRM
  filterRM[[2]] <- qualRM
  filterRM[[3]] <- table(psiRM$COMPLEX)
  filterRM[[4]] <- colnames(psiRM)[-c(1:6)]

  names(filterRM) <- c("PSI", "Qual", "EventsPerType", "Samples")
  return(filterRM)

}

# Filter PSI table (rMATS) by alternativity
# Filter previously filtered PSI table from rMATS to consider only events with PSIs between (and including) minPsi and maxPSI.
# @param filteredRMList List containing PSI and Qual (rMATS) tables, obtained with filterrMATS()
# @param minPsi (numeric) Minimum PSI to consider
# @param maxPsi (numeric) Maximum PSI to consider
#
# @return List with: 1) filtered table PSI columns, 2) filtered table "Qual" columns, including inc and exc, 3) table with number of events per type and 4) Samples
# @export
#
# @examples


alternativerMATS <- function(filteredRMList, minPsi, maxPsi){
  alternativeRM <- list()

  psiTable <- filteredRMList$PSI
  qualTable <- filteredRMList$Qual

  originalColN <- ncol(psiTable)

  # Consider alternative events only
  psiTable$AllGreaterMin  <- apply(psiTable[,-c(1:6)], 1, all_grteq_row, minPsi)
  psiTable$AllLowerMax    <- apply(psiTable[,-c(1:6)], 1, all_lweq_row, maxPsi)
  psiTable                <- psiTable[which(psiTable$AllGreaterMin == TRUE & psiTable$AllLowerMax == TRUE),]
  qualTable               <- qualTable[match(psiTable$EVENT, qualTable$EVENT),]

  # Remove columns added
  psiTable                <- psiTable[,c(1:originalColN)]
  qualTable               <- qualTable[,c(1:originalColN)]

  alternativeRM[[1]] <- psiTable
  alternativeRM[[2]] <- qualTable
  alternativeRM[[3]] <- table(psiTable$COMPLEX)
  alternativeRM[[4]] <- colnames(psiTable)[-c(1:6)]

  names(alternativeRM) <- c("PSI", "Qual", "EventsPerType", "Samples")

  return(alternativeRM)

}
