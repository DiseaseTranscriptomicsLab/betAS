
getWhippet <- function(listIncTables){

  filterWhip <- list()

  listfiles <- lapply(listIncTables, function(sample) {
    sample$eventID <- paste0(sample$Gene, ":", sample$Node, ":", sample$Type)
    sample$Ninc <- sample$Psi*sample$Total_Reads
    sample$Nexc <-sample$Total_Reads - sample$Ninc  # or sample$Ninc*(1/sample$Psi - 1)
    return(sample)
  })


  # Use one of the samples to get the common columns
  commonCols <- as.data.frame(cbind(listfiles[[1]]$Gene, listfiles[[1]]$eventID, listfiles[[1]]$Coord, rep(0, nrow(listfiles[[1]])), listfiles[[1]]$Coord, listfiles[[1]]$Type))
  colnames(commonCols) <- c("GENE","EVENT","COORD","LENGTH","FullCO","COMPLEX")

  # Get only relevant columns for each data frame in the input list of data frames, and creat the psi, inc and exc data frames (one column per sample)
  # psis <- as.data.frame(matrix(NA,1))
  # colnames(psis) <- c("eventID")
  # inc <- as.data.frame(matrix(NA,1))
  # colnames(inc) <- c("eventID")
  # exc <- as.data.frame(matrix(NA,1))
  # colnames(exc) <- c("eventID")
  # for (sample in listfiles){
  #   psis <- merge(psis, sample[,c("eventID","Psi")], by="eventID", all = TRUE)
  #   inc <- merge(inc, sample[,c("eventID","Ninc")], by="eventID", all = TRUE)
  #   exc <- merge(exc, sample[,c("eventID","Nexc")], by="eventID", all = TRUE)
  # }


  # Explanation: You essentially want to do l[[1]][, my_names], l[[2]][,
  # my_names], ... lapply applies a function to every list element. In this
  # case, the function is [, which takes rows as its first argument (we leave it
  # blank to indicate all rows), columns as its second argument (we give it
  # my_names). It returns the results in a list.

  # lapply(listfiles, "[", 1: nrow(listfiles[[1]]) , c("Psi")) gets every Psi column for each data frame in the list; nrow(listfiles[[1]]) is used just to ensure
  # that we select all rows, and only the psi column
  # Assumptions: each sample has the same events, in the same order (which should be the case if Whippet is applied for all samples at the same time)

  #psis <- data.frame(eventID=commonCols$EVENT)
  #psis[,names(listfiles)] <- as.data.frame(lapply(lapply(listfiles, "[", 1:nrow(commonCols) , c("Psi")),cbind))

  psis <- as.data.frame(lapply(lapply(listfiles, "[", 1:nrow(listfiles[[1]])  , c("Psi")),cbind))*100
  # psis[,-1] <- psis[,-1]*100

  #inc <- data.frame(eventID=commonCols$EVENT)
  #inc[,names(listfiles)] <- as.data.frame(lapply(lapply(listfiles, "[", 1:nrow(commonCols) , c("Ninc")),cbind))
  inc <- as.data.frame(lapply(lapply(listfiles, "[", 1:nrow(listfiles[[1]]) , c("Ninc")),cbind))

  #exc <- data.frame(eventID=commonCols$EVENT)
  #exc[,names(listfiles)] <- as.data.frame(lapply(lapply(listfiles, "[", 1:nrow(commonCols) , c("Nexc")),cbind))
  exc <- as.data.frame(lapply(lapply(listfiles, "[",  1:nrow(listfiles[[1]])  , c("Nexc")),cbind))

  # qual matrix, which includes inc and exc counts
  qualWhip <- cbind(commonCols,
                    # Mimicking vast-tools INCLUSION table ".Q" columns to facilitate the compatibility with other betAS functions
                    matrix(paste0("A,A,0=0=0,A,",commonCols$COMPLEX,"@", as.matrix(inc) , ",", as.matrix(exc) ), nrow = nrow(commonCols)))
  colnames(qualWhip) <- c("GENE","EVENT","COORD","LENGTH","FullCO","COMPLEX",paste0(names(listfiles),".Q"))

  # psi matrix
  psiWhip <- cbind(commonCols,psis)
  colnames(psiWhip) <- c("GENE","EVENT","COORD","LENGTH","FullCO","COMPLEX",names(listfiles))


  filterWhip[[1]] <- psiWhip
  filterWhip[[2]] <- qualWhip
  filterWhip[[3]] <- table(psiWhip$COMPLEX)
  filterWhip[[4]] <- colnames(psiWhip)[-c(1:6)]

  names(filterWhip) <- c("PSI", "Qual", "EventsPerType", "Samples")
  return(filterWhip)

}


filterWhippet(whip1, names(whip1$EventsPerType))

# to filter based on event type
filterWhippet <- function(WhippetList, types){

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

  filterWhippet[[1]] <- psiWhip
  filterWhippet[[2]] <- qualWhip
  filterWhippet[[3]] <- table(psiWhip$COMPLEX)
  filterWhippet[[4]] <- colnames(psiWhip)[-c(1:6)]

  names(filterWhippet) <- c("PSI", "Qual", "EventsPerType", "Samples")
  return(filterWhippet)

}



