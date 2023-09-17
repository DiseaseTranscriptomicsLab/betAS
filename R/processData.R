
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# A. General Functions ---------------------------------------------------------
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#' Retrieve Inclusion Tables for Specific Tools
#'
#' This function retrieves the inclusion tables specific to the chosen tool for further analyses.
#'
#' @param pathTables Character or list. Specifies the source of the dataset:
#'        \itemize{
#'          \item{\strong{vast-tools:}}{Should be a path pointing to a \code{*INCLUSION_LEVELS_FULL*.tab} file.}
#'          \item{\strong{rMATS:}}{Should be a path pointing to a \code{*MATS.JC.txt} file.}
#'          \item{\strong{whippet:}}{Should be a list of paths (minimum of 2) pointing to \code{*.psi.gz} files.}
#'        }
#'        If \code{NULL}, the function will load the default betAS test datasets for the chosen tool.
#' @param tool A character string specifying the desired tool for which the dataset
#'        should be fetched. Acceptable values are "vast-tools", "whippet", or "rMATS".
#'
#' @return Returns a data frame if the tool is "rMATS" or "vast-tools". If the tool is
#'        "whippet", the function returns a list of data frames.
#' @export
#'
#' @examples
#' # Example using the vast-tools tool with a single path
#' # getDataset(pathTables = "/path/to/vast_tools_file.txt", tool = "vast-tools")
#'
#' # Example using the whippet tool with multiple paths
#' # getDataset(pathTables = list("/path/to/whippet_file1.txt", "/path/to/whippet_file2.txt"), tool = "whippet")
#'
#' # Loading the betAS test datasets for rMATS
#' getDataset(pathTables = NULL, tool = "rMATS")
#'
getDataset <- function(pathTables=NULL, tool){

  requiredcols_vasttools <- c("GENE","EVENT","COORD","LENGTH","FullCO","COMPLEX")
  requiredcols_rMATS <- c("ID","GeneID","geneSymbol","chr","strand","IJC_SAMPLE_1","SJC_SAMPLE_1","IJC_SAMPLE_2","SJC_SAMPLE_2")
  requiredcols_whippet <- c("Gene", "Node", "Coord", "Strand", "Type" ,"Psi", "CI_Width", "Total_Reads", "Complexity", "Entropy" ,"Inc_Paths", "Exc_Paths")

  # Use example datasets
  if (is.null(pathTables)){

    if (tool == "vast-tools"){

      incData <- readRDS(file = "test/INCLUSION_LEVELS_FULL-mm10-8-v251.rds")
      colNames <- colnames(incData)

    } else if (tool == "rMATS"){

      incData <- read.delim(file = "test/SE.MATS.JC.txt")
      colNames <- colnames(incData)

    } else if (tool == "whippet"){

      incData <- readRDS(file = "test/listdfs_WHippet.rds")
      colNames <- Reduce(intersect, lapply(incData, colnames))

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

      colNames <- colnames(incData)

    } else if (tool == "whippet"){

      if (length(pathTables) == 1){

        stop("The chosen tool (", tool,") is not supported by betAS. Current tool supported by betAS are 'vast-tools', 'rMATS' or 'whippet'")

      } else {

        pathTables[grep("[.]gz",pathTables)] <- lapply(pathTables[grep("[.]gz",pathTables)],gzfile)
        incData <- lapply(pathTables,read.delim)

        colNames <- Reduce(intersect, lapply(incData, colnames))

      }

    } else {

    }

  }


  if (tool == "vast-tools" & length(intersect(colNames, requiredcols_vasttools)) != length(requiredcols_vasttools)){

    stop("The provided data is not supported by betAS. Please confirm that your data matches the input requirements.")

  }

  if (tool == "rMATS" & length(intersect(colNames, requiredcols_rMATS)) != length(requiredcols_rMATS)){

    stop("The provided data is not supported by betAS. Please confirm that your data matches the input requirements.")

  }

  if (tool == "whippet" & length(intersect(colNames, requiredcols_whippet)) != length(requiredcols_whippet)) {

    stop("The provided data is not supported by betAS. Please confirm that your data matches the input requirements.")

  }

  return(incData)

}




#' Format Inclusion Tables into a Standard betAS Format
#'
#' This function formats inclusion tables from one of the three specified tools into
#' a standard format that is compatible with betAS.
#'
#' @param incTable A data frame equivalent object (or list of data frames, when using tool=\code{whippet},
#'        representing the inclusion table for the chosen tool. This is typically the output of the `getDataset`
#'        function when not using the `betASapp()` shiny app.
#'
#' @param tool A character string specifying the tool from which the incTable is derived.
#'        It can take one of three values: "vast-tools", "whippet", or "rMATS".
#'
#' @return A list containing four data frames:
#'   \itemize{
#'     \item \strong{PSI}: A data frame with information regarding the splicing
#'           events and their inclusion levels.
#'     \item \strong{Qual}: A data frame detailing the events and information
#'           related to the coverage associated with each event. This representation
#'           follows the vast-tools representation of inclusion (inc) and
#'           exclusion (exc) reads.
#'     \item \strong{EventsPerType}: A named table summarising the number of
#'           events per type.
#'     \item \strong{Samples}: Vector containing the names of the samples considered.
#'   }
#'
#' @examples
#' # Assuming you have already loaded or created an incTable from rMATS using the getDataset() function:
#' formattedTable_rMATS <- getEvents(incTable = rMATS_incTable, tool = "rMATS")
#'
#' # Similarly, for vast-tools:
#' formattedTable_vasttools <- getEvents(incTable = vasttools_incTable, tool = "vast-tools")
#'
#' # And for whippet:
#' formattedTable_whippet <- getEvents(incTable = whippet_incTable, tool = "whippet")
#'
#' @export
#'
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


#' Verify Minimum Reads for a Given Event Across Samples
#'
#' This function checks if the coverage associated to a given event (inc+exc)
#' across samples for a specific event meet a certain threshold.
#'
#' @param quals A set of "Qual" columns (one per sample) for a given event.
#' @param N The minimum number of reads (inc + exc) for a given event in each sample.
#'
#' @return A boolean value indicating if all values of (inc + exc) across samples
#'         for a given event are greater than or equal to the specified minimum, \code{N}.
#'
#' @examples
#' # Assuming you have a set of "Qual" columns for a given event:
#' quals_data <- c(...) # Replace with your actual data
#' min_reads_check <- CheckMinReads(quals = quals_data, N = 5)
#' if (min_reads_check) {
#'   print("All samples meet the minimum reads criteria.")
#' } else {
#'   print("Not all samples meet the minimum reads criteria.")
#' }
#'
#' @export
CheckMinReads <- function(quals, N){

  # Number of samples
  n <- length(quals)

  # Split string by "@" and keep the second element with the reads "inc,exc"
  string <- as.character(quals)
  reads  <- unlist(strsplit(string, split = "@", fixed = TRUE))
  reads  <- reads[seq(2, length(reads), by = 2)]

  # Split reads by "," to separate inclusion from exclusion
  reads  <- strsplit(reads, split = ",", fixed = TRUE)
  reads  <- as.numeric(unlist(reads))
  inc <- reads[seq(1, length(reads), by = 2)]
  exc <- reads[seq(2, length(reads), by = 2)]

  # Test condition
  test <- length(which(inc + exc >= N)) == n
  return(test)

}






#' Filter Splicing Events by Quality, Event Type, and Minimum Reads
#'
#' This function filters splicing events based on quality (presence of NAs),
#' specified event types, and a minimum read count threshold.
#'
#' @param InputList A list, typically the output from \code{\link{getEvents}}.
#' @param types A character vector of column letter codes for alternative event types.
#'        For Whippet, valid types include c("CE", "AA", "AD", "IR", "TS", "TE", "AF", "AL", "BS")
#'        as described in \url{https://github.com/timbitz/Whippet.jl}.
#'        For vast-tools, valid types include c("Alt3", "Alt5", "ANN", "C1", "C2", "C3", "IR-C", "IR-S", "MIC", "S")
#'        as described in \url{https://github.com/vastgroup/vast-tools}.
#'        For rMATS, as each table only contains one event, this parameter should be set to NULL.
#' @param N An integer, representing the minimum number of reads (inc + exc) for a given event in each sample.
#'
#' @return A list containing:
#'        1. Psi filtered table.
#'        2. Qual filtered table.
#'        3. A table enumerating the number of events per type.
#'        4. A character vector, "Samples", derived from column names.
#'        Essentially, it mirrors the structure of \code{\link{getDataset}},
#'        but with events filtered by coverage and type.
#'
#' @examples
#' # Assuming you have an InputList from getEvents():
#' input_data <- getEvents(...) # Replace with your actual data
#' # Filtering whippet events, to isolate only CE and AA events with a minimum of 10 reads in each sample
#' filtered_data <- filterEvents(InputList = input_data, types = c("CE", "AA"), N = 10)
#'
#' @export
#'
filterEvents <- function(InputList, types=NULL, N=10){

  filterList <- InputList

  psi <- filterList$PSI
  qual <- filterList$Qual

  # Remove events containing at least one NA
  psi$AnyNA    <- apply(psi[,-c(1:6)], 1, anyNA)
  psi          <- psi[which(psi$AnyNA == FALSE),]
  psi$AnyNA    <- NULL
  qual         <- qual[match(psi$EVENT, qual$EVENT),]

  if (!is.null(types)){
    psi <- psi[which(psi$COMPLEX %in% types),]
  }


  # Calculate coverage/balance (use .Q columns)
  qual$AllminReads <- apply(qual[,grep("[.]Q", colnames(qual))], 1, FUN = function(X) CheckMinReads(X,N))
  qual             <- qual[which(qual$AllminReads == TRUE),]
  psi              <- psi[match(qual$EVENT, psi$EVENT),]
  qual$AllminReads <-  NULL

  filterList[[1]] <- psi
  filterList[[2]] <- qual
  filterList[[3]] <- table(psi$COMPLEX)
  filterList[[4]] <- colnames(psi)[-c(1:6)]

  names(filterList) <- c("PSI", "Qual", "EventsPerType", "Samples")
  return(filterList)


}



#' Filter Splicing Events by PSI Range
#'
#' This function filters splicing events based on a specified PSI range.
#' Events with a PSI value outside of the provided range are removed.
#'
#' @param filteredList A list, typically the output from \code{\link{filterEvents}} or
#'                     \code{\link{getEvents}}.
#' @param minPsi A numeric value between 0 and 100, representing the minimum PSI threshold.
#' @param maxPsi A numeric value between 0 and 100, representing the maximum PSI threshold.
#'
#' @details
#' The function filters out events whose PSI value (percentage spliced in) does not fall
#' within the [minPsi, maxPsi] range. Both bounds are inclusive.
#'
#' @return A list containing:
#'        1. PSI filtered table based on the provided PSI range.
#'        2. Corresponding Qual filtered table.
#'        3. A table enumerating the number of events per type.
#'        4. A character vector, "Samples", derived from column names.
#'        This mirrors the structure of \code{\link{getDataset}},
#'        but with additional filtering based on the PSI range.
#'
#' @examples
#' # Assuming you have a filteredList from filterEvents():
#' filtered_data <- filterEvents(...) # Replace with your actual data
#' alternative_data <- alternativeEvents(filteredList = filtered_data, minPsi = 20, maxPsi = 80)
#'
#' @export
alternativeEvents  <- function(filteredList, minPsi, maxPsi){

  alternative <- list()

  psiTable <- filteredList$PSI
  qualTable <- filteredList$Qual

  originalColN <- ncol(psiTable)

  # Consider alternative events only
  allGreaterMin  <- rowSums(psiTable[,-c(1:6)] < minPsi) > 0
  allLowerMax    <- rowSums(psiTable[,-c(1:6)] > maxPsi) > 0

  psiTable                <- psiTable[which(allGreaterMin == TRUE & allLowerMax == TRUE),]
  qualTable               <- qualTable[match(psiTable$EVENT, qualTable$EVENT),]

  # Remove columns added
  psiTable                <- psiTable[,c(1:originalColN)]
  qualTable               <- qualTable[,c(1:originalColN)]

  alternative[[1]] <- psiTable
  alternative[[2]] <- qualTable
  alternative[[3]] <- table(psiTable$COMPLEX)
  alternative[[4]] <- colnames(psiTable)[-c(1:6)]

  names(alternative) <- c("PSI", "Qual", "EventsPerType", "Samples")

  return(alternative)

}


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# B. Tool-specific functions ---------------------------------------------------
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#' Format vast-tools Inclusion Tables for betAS Compatibility
#'
#' This function reformats inclusion tables specifically from vast-tools
#' into a standardized format that is compatible with betAS.
#'
#' @param incTable A data frame or similar object representing the inclusion table
#'                 specifically from vast-tools. This is typically the output of the
#'                 `getDataset` function when not using the `betASapp()` shiny app.
#'                 Can also be a data frame obtained by reading a
#'                 \code{*INCLUSION_LEVELS_FULL*.tab} file.
#'
#' @details
#' vast-tools is one of several splicing tools that output inclusion tables.
#' For the sake of consistent analysis within betAS, these tables need to be
#' standardized. This function processes vast-tools output and structures it
#' in a way that betAS can work with, ensuring compatibility and ease of further analysis.
#'
#' @return A list containing:
#'        \itemize{
#'        \item{"PSI"}{A data frame containing information regarding the splicing events
#'                      and their inclusion level.}
#'        \item{"Qual"}{A data frame with information on the events as well as information
#'                       regarding the coverage associated with each, emulating vast-tools
#'                       way of representing inc and exc reads.}
#'        \item{"EventsPerType"}{A named table summarizing the number of events per type.}
#'        \item{"Samples"}{A character vector representing the names of the samples considered.}
#'        }
#'
#' @export
#'
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




#' Format Whippet Inclusion Tables for betAS Compatibility
#'
#' This function reformats inclusion tables specifically from Whippet
#' into a standardized format that is compatible with betAS.
#'
#' @param incTable A list of paths (minimum length of 2) pointing to
#'                 the inclusion tables specifically from Whippet. This is
#'                 typically the output of the `getDataset` function when not
#'                 using the `betASapp()` shiny app. Can also be a list of data frames
#'                 obtained by reading a list of \code{*.psi.gz} files.
#'
#' @details
#' Whippet is one of several splicing tools that output inclusion tables.
#' For the sake of consistent analysis within betAS, these tables need to be
#' standardized. This function processes Whippet output and structures it
#' in a way that betAS can work with, ensuring compatibility and ease of further analysis.
#'
#' @return A list containing:
#'        \itemize{
#'        \item{"PSI"}{A data frame containing information regarding the splicing events
#'                      and their inclusion level.}
#'        \item{"Qual"}{A data frame with information on the events as well as information
#'                       regarding the coverage associated with each, emulating Whippet
#'                       way of representing inc and exc reads.}
#'        \item{"EventsPerType"}{A named table summarizing the number of events per type.}
#'        \item{"Samples"}{A character vector representing the names of the samples considered.}
#'        }
#'
#' @export
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


#' Format rMATS Inclusion Tables for betAS Compatibility
#'
#' This function reformats inclusion tables specifically from rMATS
#' into a standardized format that is compatible with betAS.
#'
#' @param incTable A data frame or similar object representing
#'                 the inclusion table specifically from rMATS. This is
#'                 typically the output of the `getDataset` function when not
#'                 using the `betASapp()` shiny app. Can also be a data frame obtained
#'                 by reading a \code{*MATS.JC.txt} file
#'
#' @details
#' rMATS is one of several splicing tools that output inclusion tables.
#' For the sake of consistent analysis within betAS, these tables need to be
#' standardized. This function processes rMATS output and structures it
#' in a way that betAS can work with, ensuring compatibility and ease of further analysis.
#'
#' @return A list containing:
#'        \itemize{
#'        \item{"PSI"}{A data frame containing information regarding the splicing events
#'                      and their inclusion level.}
#'        \item{"Qual"}{A data frame with information on the events as well as information
#'                       regarding the coverage associated with each, emulating rMATS
#'                       way of representing inc and exc reads.}
#'        \item{"EventsPerType"}{A named table summarizing the number of events per type.}
#'        \item{"Samples"}{A character vector representing the names of the samples considered.}
#'        }
#'
#' @export
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


