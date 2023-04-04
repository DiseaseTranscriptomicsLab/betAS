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
#
# @return List with: 1) filtered table PSI columns, 2) filtered table Qual columns, 3) table with number of events per type and 4) Samples (based on colnames)
# @export
#
# @examples

filterVastTools <- function(VTlist, types){

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
  qualVAST$AllminVLOW  <- apply(qualVAST[,grep("[.]Q", colnames(qualVAST))], 1, VT_all_minVLOW_tags)
  qualVAST             <- qualVAST[which(qualVAST$AllminVLOW == TRUE),]
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

