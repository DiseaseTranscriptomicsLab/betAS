# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# A. Vast-tools ------------------------------------------------------------
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Tests for minimal coverage based (vast-tools')
# Checks if "Qual" columns from vast-tools' INCLUSION table are classified at least as "VLOW" based on score 2 (https://github.com/vastgroup/vast-tools)
# @param quals Set of "Qual" columns from vast-tools' INCLUSION table (one per sample) for a given event
#
# @return TRUE if all samples have minimal coverage.
# @export
#
# @examples
VT_all_minVLOW_tags <- function(quals){

  #Number of samples
  n <- length(quals)

  #Split string by "@" and keep the second element: "inc,exc"
  string <- as.character(unlist(quals))
  split <- unlist(strsplit(toString(unlist(string)), split = "[,@]"))
  #corrected reads' classifications are the 2nd score (of 7 fields after splitting by @) in each sample
  tags <- split[seq(from = 2, to = 7*n, by = 7)]
  test <- length(which(tags == "SOK" | tags == "OK" | tags == "VLOW" | tags == "LOW")) == n

  return(test)

}

# Tests for minimal coverage based (vast-tools')
# Checks if "Qual" columns from vast-tools' INCLUSION table have at least N reads in total (inc + exc)
# @param quals Set of "Qual" columns from vast-tools' INCLUSION table (one per sample) for a given event
# @param N minimum number of reads (inc + exc) for a given event in each sample
#
# @return TRUE if all samples have minimal coverage.
# @export
#
# @examples
VT_all_minReads <- function(quals, N){


  # Number of samples
  n <- length(quals)

  # Split string by "@" and keep the second element: "inc,exc"
  string <- as.character(unlist(quals))
  split <- unlist(strsplit(toString(unlist(string)), split = "[,@]"))

  # Inc reads are the 6th score in each sample; Exc reads the last and 7th
  inc <- as.numeric(as.vector(split[seq(from = 6, to = 7*n, by = 7)]))
  exc <- as.numeric(as.vector(split[seq(from = 7, to = 7*n, by = 7)]))

  test <- length(which(inc + exc >= N)) == n

  return(test)

}




# Format INCLUSION table (vast-tools) for further analyses
# @param incTable vast-tools' INCLUSION table
#
# @return List with: 1) table PSI columns, 2) table Qual columns, 3) table with number of events per type and 4) Samples (based on colnames)
# @export
#
# @examples



getVastTools <- function(incTable){

  filterVT <- list()

  # Structure of vast-tools "INCLUSION(...)" table is the following:
  # . 6 columns describing the event
  # . one column per sample with the PSI (between 0 and 100)
  # . one column per sample with the "Quality" information, including inc and exc
  # . each Quality column follows respective PSI column

  qual_cols <- grep("[.]Q", colnames(incTable))
  psi_cols  <- qual_cols-1

  psiVAST  <- incTable[,c(1:6,psi_cols)]
  qualVAST <- incTable[,c(1:6,qual_cols)]

  filterVT[[1]] <- psiVAST
  filterVT[[2]] <- qualVAST
  filterVT[[3]] <- table(psiVAST$COMPLEX)
  filterVT[[4]] <- colnames(psiVAST)[-c(1:6)]

  names(filterVT) <- c("PSI", "Qual", "EventsPerType", "Samples")
  return(filterVT)

}


# Filter INCLUSION table (vast-tools) for quality and event type
# Filter table from vast-tools to remove events containing NAs in at least one sample and those that do not have minimal coverage based on VT_all_minVLOW_tags() and split PSI and Qual tables
# @param VTlist list containing PSI and Qual tables, as well as event and samples, obtained with getVastTools()
# @param types (character) vast-tools' COMPLEX column letter code for alternative event types (Alt3, Alt5, ANN, C1, C2, C3, IR-C, IR-S, MIC or S) as described in https://github.com/vastgroup/vast-tools
# @param N minimum number of reads (inc + exc) for a given event in each sample
#
# @return List with: 1) filtered table PSI columns, 2) filtered table Qual columns, 3) table with number of events per type and 4) Samples (based on colnames)
# @export
#
# @examples

filterVastTools <- function(VTlist, types, N){

  filterVT <- VTlist

  psiVAST <- filterVT$PSI
  qualVAST <- filterVT$Qual

  originalColN <- ncol(psiVAST)

  # *** to be made REACTIVE!!! ***
  # Consider exon skipping events only
  # psiVAST <- psiVAST[which(psiVAST$COMPLEX == "S"),]
  psiVAST <- psiVAST[which(psiVAST$COMPLEX %in% types),]

  # Remove events containing at least one NA
  psiVAST$AnyNA    <- apply(psiVAST, 1, anyNA)
  psiVAST          <- psiVAST[which(psiVAST$AnyNA == FALSE),]
  psiVAST          <- psiVAST[,-c(ncol(psiVAST))]
  qualVAST         <- qualVAST[match(psiVAST$EVENT, qualVAST$EVENT),]

  # *** to be made REACTIVE!!! ***
  # should not be done while choosing only ES events
  # Filter out events from "ANN" module out
  # annot             <- which(psiVAST$COMPLEX == "ANN")
  # psiVAST           <- psiVAST[-c(annot),]
  # qualVAST          <- qualVAST[match(psiVAST$EVENT, qualVAST$EVENT),]

  # Calculate coverage/balance (use .Q columns)
  # qualVAST$AllminVLOW  <- apply(qualVAST[,grep("[.]Q", colnames(qualVAST))], 1, VT_all_minVLOW_tags)
  #   qualVAST             <- qualVAST[which(qualVAST$AllminVLOW == TRUE),]
  # psiVAST              <- psiVAST[match(qualVAST$EVENT, psiVAST$EVENT),]

  qualVAST$AllminReads <-  apply(qualVAST[,grep("[.]Q", colnames(qualVAST))], 1, FUN = function(X) VT_all_minReads(X,N))
  qualVAST             <- qualVAST[which(qualVAST$AllminReads == TRUE),]
  psiVAST              <- psiVAST[match(qualVAST$EVENT, psiVAST$EVENT),]

  # Remove columns added
  psiVAST               <- psiVAST[,c(1:originalColN)]
  qualVAST              <- qualVAST[,c(1:originalColN)]

  filterVT[[1]] <- psiVAST
  filterVT[[2]] <- qualVAST
  filterVT[[3]] <- table(psiVAST$COMPLEX)
  filterVT[[4]] <- colnames(psiVAST)[-c(1:6)]

  names(filterVT) <- c("PSI", "Qual", "EventsPerType", "Samples")
  return(filterVT)

}

#' Filter PSI table (vast-tools) by alternativity
#' Filter previously quality-filtered PSI table from vast-tools to consider only events with PSIs between (and including) minPsi and maxPSI.
#' @param filteredVTList List containing PSI and Qual (vast-tools) tables, obtained with filterVastTools()
#' @param minPsi (numeric) Minimum PSI to consider
#' @param maxPsi (numeric) Maximum PSI to consider
#'
#' @return List with: 1) filtered table PSI columns, 2) filtered table Qual columns, 3) table with number of events per type and 4) Samples (based on colnames)
#' @export
#'
#' @examples
#' testTable <- betAS:::testTable
#' # filter table to consider only events all PSIs between 1 and 99
#' alternativeVastTools(testTable, minPsi = 1, maxPsi = 99)
alternativeVastTools <- function(filteredVTList, minPsi, maxPsi){

  alternativeVT <- list()

  psiTable <- filteredVTList$PSI
  qualTable <- filteredVTList$Qual

  originalColN <- ncol(psiTable)

  # Consider alternative events only
  psiTable$AllGreaterMin  <- apply(psiTable[,-c(1:6)], 1, all_grteq_row, minPsi)
  psiTable$AllLowerMax    <- apply(psiTable[,-c(1:6)], 1, all_lweq_row, maxPsi)
  psiTable                <- psiTable[which(psiTable$AllGreaterMin == TRUE & psiTable$AllLowerMax == TRUE),]
  qualTable               <- qualTable[match(psiTable$EVENT, qualTable$EVENT),]

  # Remove columns added
  psiTable                <- psiTable[,c(1:originalColN)]
  qualTable               <- qualTable[,c(1:originalColN)]

  alternativeVT[[1]] <- psiTable
  alternativeVT[[2]] <- qualTable
  alternativeVT[[3]] <- table(psiTable$COMPLEX)
  alternativeVT[[4]] <- colnames(psiTable)[-c(1:6)]

  names(alternativeVT) <- c("PSI", "Qual", "EventsPerType", "Samples")

  return(alternativeVT)

}


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# B. Whippet ---------------------------------------------------------------
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# Tests for minimal coverage based (rMATS')
# Checks if "Qual" columns from rMATS' INCLUSION table have at least N reads in total (inc + exc)
# @param quals Set of "Qual" columns from rMATS' INCLUSION table (one per sample) for a given event
# @param N minimum number of reads (inc + exc) for a given event in each sample
#
# @return TRUE if all samples have minimal coverage.
# @export
#
# @examples
Whippet_all_minReads <- function(quals, N){


  # Number of samples
  n <- length(quals)

  # Split string by "@" and keep the second element: "inc,exc"
  string <- as.character(unlist(quals))
  split <- unlist(strsplit(toString(unlist(string)), split = "[,@]"))

  # Inc reads are the 6th score in each sample; Exc reads the last and 7th
  inc <- as.numeric(as.vector(split[seq(from = 6, to = 7*n, by = 7)]))
  exc <- as.numeric(as.vector(split[seq(from = 7, to = 7*n, by = 7)]))

  test <- length(which(inc + exc >= N)) == n

  return(test)

}



# Format *.psi.gz tables (whippet) for quantified PSIs, for further analyses
# @param incTable list of whippet's *.psi.gz tables, one per sample
#
# @return List with: 1) filtered table PSI columns, 2) filtered table "Qual" columns including inc and exc, 3) table with number of events per type and 4) Samples
# @export
#
# @examples

getWhippet <- function(listIncTables){

  filterWhip <- list()

  listfiles <- lapply(listIncTables, function(sample) {

    # Create column for eventID = GENE:NODE:TYPE
    sample$eventID <- paste0(sample$Gene, ":", sample$Node, ":", sample$Type)
    # Extract inclusion reads based on Ninc ~ Binomial(n=Ntotal, p=PSI) (See supplementary material from Whippet paper https://doi.org/10.1016/j.molcel.2018.08.018)
    sample$Ninc <- sample$Psi*sample$Total_Reads
    # Get exclusion reads based on Ntotal = Ninc + Nexc
    sample$Nexc <-sample$Total_Reads - sample$Ninc  # or sample$Ninc*(1/sample$Psi - 1)

    # Create auxiliary column that will be needed for the final table
    sample$length <- rep(0,nrow(sample))

    # Keep only necessary columns
    sample <- sample[,c(1,15,3,18,3,5,6,16,17)]
    colnames(sample) <- c("GENE","EVENT","COORD","LENGTH","FullCO","COMPLEX","PSI","INC","EXC")

    return(sample)

  })

  # Change column names so that PSI, INC and EXC are unique for each sample. Will be useful to merge
  newcolnames <- lapply(names(listfiles), function(x) paste0(c("PSI","INC","EXC"),".",x))
  for(i in 1:length(listfiles)){
    colnames(listfiles[[i]])[7:9] <- newcolnames[[i]]
  }

  # Merge all tables; might take a while, but ensures that the universe of events is the same
  mergeTable <- as.data.frame(Reduce(function(dtf1, dtf2) merge(dtf1, dtf2,
                                                                by.x=c("GENE","EVENT","COORD","LENGTH","FullCO","COMPLEX"),
                                                                by.y=c("GENE","EVENT","COORD","LENGTH","FullCO","COMPLEX"),
                                                                all = TRUE),
                                     listfiles))

  # psi matrix
  colpsis <- grep("PSI.",colnames(mergeTable))
  psiWhip <- mergeTable[,c(1:6,colpsis)]
  colpsis <- grep("PSI.",colnames(psiWhip))
  psiWhip[,colpsis] <- psiWhip[,colpsis]*100
  colnames(psiWhip)[colpsis] <- gsub("PSI.","",colnames(psiWhip)[colpsis])

  # qual matrix, which includes inc and exc counts

  inc <- as.matrix(mergeTable[,c(grep("INC.",colnames(mergeTable)))])
  exc <- as.matrix(mergeTable[,c(grep("EXC.",colnames(mergeTable)))])

  qualWhip <- cbind(mergeTable[,c(1:6)],
                    # Mimicking vast-tools INCLUSION table ".Q" columns to facilitate the compatibility with other betAS functions
                    matrix(paste0("A,A,0=0=0,A,",mergeTable$COMPLEX,"@", inc , ",", exc ), nrow = nrow(mergeTable[,c(1:6)])))
  colnames(qualWhip) <- c("GENE","EVENT","COORD","LENGTH","FullCO","COMPLEX",paste0(names(listfiles),".Q"))


  filterWhip[[1]] <- psiWhip
  filterWhip[[2]] <- qualWhip
  filterWhip[[3]] <- table(psiWhip$COMPLEX)
  filterWhip[[4]] <- colnames(psiWhip)[-c(1:6)]

  names(filterWhip) <- c("PSI", "Qual", "EventsPerType", "Samples")
  return(filterWhip)

}


# Filter whippet tables for quality (NAs) and event type
# Filter table from whippet to remove events containing NAs in at least one sample
# @param WhippetList list containing PSI and Qual tables, as well as event and samples, obtained with getWhippet()
# @param types (character) whippet's Type column letter code for alternative event types (CE,AA,AD,IR,TS,TE,AF,AL,BS) as described in https://github.com/timbitz/Whippet.jl
# @param N minimum number of reads (inc + exc) for a given event in each sample
#
# @return List with: 1) filtered table PSI columns, 2) filtered table Qual columns, 3) table with number of events per type and 4) Samples (based on colnames)
# @export
#
# @examples

filterWhippet <- function(WhippetList, types, N){

  filterWhippet <- WhippetList

  psiWhip <- filterWhippet$PSI
  qualWhip <- filterWhippet$Qual


  # Remove events containing at least one NA
  psiWhip$AnyNA    <- apply(psiWhip[,-c(1:6)], 1, anyNA)
  psiWhip          <- psiWhip[which(psiWhip$AnyNA == FALSE),]
  psiWhip          <- psiWhip[,-c(ncol(psiWhip))]
  qualWhip         <- qualWhip[match(psiWhip$EVENT, qualWhip$EVENT),]




  # *** to be made REACTIVE!!! ***
  # Consider exon skipping events only
  # psiVAST <- psiVAST[which(psiVAST$COMPLEX == "S"),]
  psiWhip <- psiWhip[which(psiWhip$COMPLEX %in% types),]

  # Calculate coverage/balance (use .Q columns)
  qualWhip$AllminReads <- apply(qualWhip[,grep("[.]Q", colnames(qualWhip))], 1, FUN = function(X) Whippet_all_minReads(X,N))
  qualWhip             <- qualWhip[which(qualWhip$AllminReads == TRUE),]
  psiWhip              <- psiWhip[match(qualWhip$EVENT, psiWhip$EVENT),]
  qualWhip$AllminReads <-  NULL

  filterWhippet[[1]] <- psiWhip
  filterWhippet[[2]] <- qualWhip
  filterWhippet[[3]] <- table(psiWhip$COMPLEX)
  filterWhippet[[4]] <- colnames(psiWhip)[-c(1:6)]

  names(filterWhippet) <- c("PSI", "Qual", "EventsPerType", "Samples")
  return(filterWhippet)

}

# Filter PSI table (whippet) by alternativity
# Filter previously filtered PSI table from whippet to consider only events with PSIs between (and including) minPsi and maxPSI.
# @param filteredWhippetList List containing PSI and Qual (whippet) tables, obtained with filterWhippet()
# @param minPsi (numeric) Minimum PSI to consider
# @param maxPsi (numeric) Maximum PSI to consider
#
# @return List with: 1) filtered table PSI columns, 2) filtered table "Qual" columns, including inc and exc, 3) table with number of events per type and 4) Samples
# @export
#
# @examples


alternativeWhippet <- function(filteredWhippetList, minPsi, maxPsi){

  alternativeWhippet <- list()

  psiTable <- filteredWhippetList$PSI
  qualTable <- filteredWhippetList$Qual

  originalColN <- ncol(psiTable)

  # Consider alternative events only
  psiTable$AllGreaterMin  <- apply(psiTable[,-c(1:6)], 1, all_grteq_row, minPsi)
  psiTable$AllLowerMax    <- apply(psiTable[,-c(1:6)], 1, all_lweq_row, maxPsi)
  psiTable                <- psiTable[which(psiTable$AllGreaterMin == TRUE & psiTable$AllLowerMax == TRUE),]
  qualTable               <- qualTable[match(psiTable$EVENT, qualTable$EVENT),]

  # Remove columns added
  psiTable                <- psiTable[,c(1:originalColN)]
  qualTable               <- qualTable[,c(1:originalColN)]

  alternativeWhippet[[1]] <- psiTable
  alternativeWhippet[[2]] <- qualTable
  alternativeWhippet[[3]] <- table(psiTable$COMPLEX)
  alternativeWhippet[[4]] <- colnames(psiTable)[-c(1:6)]

  names(alternativeWhippet) <- c("PSI", "Qual", "EventsPerType", "Samples")

  return(alternativeWhippet)

}


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# C. rMATS -----------------------------------------------------------------
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# Tests for minimal coverage based (rMATS')
# Checks if "Qual" columns from rMATS' INCLUSION table have at least N reads in total (inc + exc)
# @param quals Set of "Qual" columns from rMATS' INCLUSION table (one per sample) for a given event
# @param N minimum number of reads (inc + exc) for a given event in each sample
#
# @return TRUE if all samples have minimal coverage.
# @export
#
# @examples
RM_all_minReads <- function(quals, N){


  # Number of samples
  n <- length(quals)

  # Split string by "@" and keep the second element: "inc,exc"
  string <- as.character(unlist(quals))
  split <- unlist(strsplit(toString(unlist(string)), split = "[,@]"))

  # Inc reads are the 6th score in each sample; Exc reads the last and 7th
  inc <- as.numeric(as.vector(split[seq(from = 6, to = 7*n, by = 7)]))
  exc <- as.numeric(as.vector(split[seq(from = 7, to = 7*n, by = 7)]))

  test <- length(which(inc + exc >= N)) == n

  return(test)

}



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
    eventCoordinates <- paste0(incTable$chr,":",incTable$exonStart_0base + 1,"-",incTable$exonEnd)
    li_factor <- 2
    ls_factor <- 1

    # NOTE for "MXE" events:
    # reading tables using read.table(file, sep="\t", header=TRUE, quote=""), which is the case in the app, corrects colnames starting with numbers by adding an X
  } else if ("X1stExonStart_0base" %in% colnames){
    eventType <- "MXE"
    eventCols <- c("X1stExonStart_0base", "X1stExonEnd", "X2ndExonStart_0base", "X2ndExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE")
    # If the strand is +, then the inclusion form includes the 1st exon (1stExonStart_0base, 1stExonEnd) and skips the 2nd exon
    # If the strand is -, then the inclusion form includes the 2nd exon (2ndExonStart_0base, 2ndExonEnd) and skips the 1st exon
    eventCoordinates <- paste0(incTable$chr,":(1st)",incTable$X1stExonStart_0base + 1,"-",incTable$X1stExonEnd,";(2nd)",incTable$X2ndExonStart_0base + 1,"-",incTable$X2ndExonEnd)
    li_factor <- 2
    ls_factor <- 2

  } else if ("longExonStart_0base"  %in% colnames){
    eventType <- "Altss"
    eventCols <- c("longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")
    # The inclusion form includes the long exon (longExonStart_0base, longExonEnd) instead of the short exon (shortES shortEE)
    eventCoordinates <- paste0(incTable$chr,":",incTable$longExonStart_0base + 1,"-",incTable$longExonEnd)
    li_factor <- 2
    ls_factor <- 1

  } else if ("riExonStart_0base"  %in% colnames){
    eventType <- "IR"
    eventCols <- c("riExonStart_0base", "riExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE")
    # The inclusion form includes (retains) the intron (upstreamEE, downstreamES)
    eventCoordinates <- paste0(incTable$chr,":",incTable$upstreamEE,"-",incTable$downstreamES)
    li_factor <- 2
    ls_factor <- 1

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
    # inc <-  cbind(apply(matrix(unlist(strsplit(incTable$IJC_SAMPLE_1, ",")),ncol=Nsamples_Group1,byrow=T), 2, as.numeric)/incTable$IncFormLen)
    # exc <-  cbind(apply(matrix(unlist(strsplit(incTable$SJC_SAMPLE_1, ",")),ncol=Nsamples_Group1,byrow=T), 2, as.numeric)/incTable$SkipFormLen)
    # inc <-  cbind(apply(matrix(unlist(strsplit(incTable$IJC_SAMPLE_1, ",")),ncol=Nsamples_Group1,byrow=T), 2, as.numeric)*(incTable$SkipFormLen/incTable$IncFormLen))
    # exc <-  cbind(apply(matrix(unlist(strsplit(incTable$SJC_SAMPLE_1, ",")),ncol=Nsamples_Group1,byrow=T), 2, as.numeric))
    inc <-  cbind(apply(matrix(unlist(strsplit(incTable$IJC_SAMPLE_1, ",")),ncol=Nsamples_Group1,byrow=T), 2, as.numeric)/li_factor)
    exc <-  cbind(apply(matrix(unlist(strsplit(incTable$SJC_SAMPLE_1, ",")),ncol=Nsamples_Group1,byrow=T), 2, as.numeric)/ls_factor)


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
    inc <-  cbind(apply(matrix(unlist(strsplit(incTable$IJC_SAMPLE_1, ",")),ncol=Nsamples_Group1,byrow=T), 2, as.numeric)/li_factor,
                  apply(matrix(unlist(strsplit(incTable$IJC_SAMPLE_2, ",")),ncol=Nsamples_Group2,byrow=T), 2, as.numeric)/li_factor)
    exc <-  cbind(apply(matrix(unlist(strsplit(incTable$SJC_SAMPLE_1, ",")),ncol=Nsamples_Group1,byrow=T), 2, as.numeric)/ls_factor,
                  apply(matrix(unlist(strsplit(incTable$SJC_SAMPLE_2, ",")),ncol=Nsamples_Group2,byrow=T), 2, as.numeric)/ls_factor)


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
# @param N minimum number of reads (inc + exc) for a given event in each sample
#
# @return List with: 1) filtered table PSI columns, 2) filtered table "Qual" columns including inc and exc, 3) table with number of events per type and 4) Samples
# @export
#
# @examples


filterrMATS <- function(RMlist, N){

  filterRM <- RMlist


  psiRM <- filterRM$PSI
  qualRM <- filterRM$Qual

  # Remove events containing at least one NA
  psiRM$AnyNA    <- apply(psiRM, 1, anyNA)
  psiRM          <- psiRM[which(psiRM$AnyNA == FALSE),]
  psiRM          <- psiRM[,-c(ncol(psiRM))]
  qualRM         <- qualRM[match(psiRM$EVENT, qualRM$EVENT),]



  # Calculate coverage/balance (use .Q columns)
  qualRM$AllminReads <- apply(qualRM[,grep("[.]Q", colnames(qualRM))], 1, FUN = function(X) RM_all_minReads(X,N))
  qualRM             <- qualRM[which(qualRM$AllminReads == TRUE),]
  psiRM              <- psiRM[match(qualRM$EVENT, psiRM$EVENT),]
  qualRM$AllminReads <-  NULL

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



# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# D. Import data function --------------------------------------------------
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# get dataset tables. only useful for package, betASapp has a different way of dealing with it
# pathTables=NULL when using example datasets

# when running this script, one should be in the main folder

getDataset <- function(pathTables=NULL, tool){

  requiredcols_vasttools <- c("GENE","EVENT","COORD","LENGTH","FullCO","COMPLEX")
  requiredcols_rMATS <- c("ID","GeneID","geneSymbol","chr","strand","IJC_SAMPLE_1","SJC_SAMPLE_1","IJC_SAMPLE_2","SJC_SAMPLE_2")
  requiredcols_whippet <- c("Gene", "Node", "Coord", "Strand", "Type" ,"Psi", "CI_Width", "Total_Reads", "Complexity", "Entropy" ,"Inc_Paths", "Exc_Paths")

  # Use example datasets
  if (is.null(pathTables)){

    if (tool == "vast-tools"){

      incData <- readRDS(file = "test/INCLUSION_LEVELS_FULL-mm10-8-v251.rds")

    } else if (tool == "rMATS"){

      incData <- read.delim(file = "test/SE.MATS.JC.txt")

    } else if (tool == "whippet"){

      incData <- readRDS(file = "test/listdfs_WHippet.rds")

    } else {

      stop(paste0("Tool '", tool, "' is not supported by the current version of betAS."))

    }


    # Use user-input data
  } else {

    if (tool %in% c("vast-tools", "rMATS")){

      if(length(grep(pattern = "[.]gz", x = pathTables)) == 0){

        incData <- pathTables

      } else {

        incData <- gzfile(pathTables)

      }

      incData <- read.delim(incData)

    } else if (tool == "whippet"){

      if (length(pathTables) == 1){

        stop("Number of files not supported. Please select at least two files when using whippet inclusion tables.")

      } else {

        pathTables[grep("[.]gz",pathTables)] <- lapply(pathTables[grep("[.]gz",pathTables)],gzfile)
        incData <- lapply(pathTables,read.delim)

      }

    }

  }


  if (tool == "vast-tools" & length(incData[incData %in% requiredcols_vasttools]) == length(requiredcols_vasttools)){

    stop("The provided data is not supported by betAS. Please confirm that your data matches the input requirements.")

  } else if (tool == "rMATS" & length(incData[incData %in% requiredcols_rMATS]) == length(requiredcols_rMATS)){

    stop("The provided data is not supported by betAS. Please confirm that your data matches the input requirements.")

  } else if (tool == "whippet" & length( unique(lapply(incData, colnames))[[1]][unique(lapply(incData, colnames))[[1]] %in% requiredcols_rMATS]) == length(requiredcols_rMATS)) {

    stop("The provided data is not supported by betAS. Please confirm that your data matches the input requirements.")

  }

  return(incData)

}


# incTable is an inclusion table (or a list of inclusion tables). List of paths only works for Whippet files

getEvents <- function(incTable, tool){

  if (tool == "vast-tools"){

    ASList <- getVastTools(incTable)

  } else if (tool == "rMATS") {

    suppressWarnings({
      ASList <- getrMATS(incTable)
    })

  } else if (tool == "whippet") {

    ASList <- getWhippet(incTable)

  } else {
    stop(paste0("Tool '", tool, "' is not supported by the current version of betAS."))
  }

  return(ASList)

}

# output from getPSIs()
# types should follow the nomenclature of the tool being used
# types will not be used for rMATS results

filterEvents <- function(ASList, types=NULL, N=10, tool){

  if (is.null(types)){
    types <- names(ASList$EventsPerType)
  }

  if (tool == "vast-tools"){

    ASListFiltered <- filterVastTools(ASList, types, N)

  } else if (tool == "rMATS") {

    ASListFiltered <- filterrMATS(ASList, N)

  } else if (tool == "whippet") {

    ASListFiltered <- filterWhippet(ASList, types, N)

  } else {
    stop(paste0("Tool '", tool, "' is not supported by the current version of betAS."))
  }

  return(ASListFiltered)

}

# limits in percentage
alternativeEvents <- function(ASListFiltered, minPsi=1, maxPsi=99, tool ){


  if (tool == "vast-tools"){

    ASListAlternative <- alternativeVastTools(ASListFiltered, minPsi, maxPsi)

  } else if (tool == "whippet"){

    ASListAlternative <-  alternativeWhippet(ASListFiltered, minPsi, maxPsi)

  } else if (tool == "rMATS"){

    ASListAlternative <- alternativerMATS(ASListFiltered, minPsi, maxPsi)

  }

  return(ASListAlternative)

}






