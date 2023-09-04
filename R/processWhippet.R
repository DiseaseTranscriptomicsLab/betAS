# # Tests for minimal coverage based (rMATS')
# # Checks if "Qual" columns from rMATS' INCLUSION table have at least N reads in total (inc + exc)
# # @param quals Set of "Qual" columns from rMATS' INCLUSION table (one per sample) for a given event
# # @param N minimum number of reads (inc + exc) for a given event in each sample
# #
# # @return TRUE if all samples have minimal coverage.
# # @export
# #
# # @examples
# Whippet_all_minReads <- function(quals, N){
#
#
#   # Number of samples
#   n <- length(quals)
#
#   # Split string by "@" and keep the second element: "inc,exc"
#   string <- as.character(unlist(quals))
#   split <- unlist(strsplit(toString(unlist(string)), split = "[,@]"))
#
#   # Inc reads are the 6th score in each sample; Exc reads the last and 7th
#   inc <- as.numeric(as.vector(split[seq(from = 6, to = 7*n, by = 7)]))
#   exc <- as.numeric(as.vector(split[seq(from = 7, to = 7*n, by = 7)]))
#
#   test <- length(which(inc + exc >= N)) == n
#
#   return(test)
#
# }
#
#
#
# # Format *.psi.gz tables (whippet) for quantified PSIs, for further analyses
# # @param incTable list of whippet's *.psi.gz tables, one per sample
# #
# # @return List with: 1) filtered table PSI columns, 2) filtered table "Qual" columns including inc and exc, 3) table with number of events per type and 4) Samples
# # @export
# #
# # @examples
#
# getWhippet <- function(listIncTables){
#
#   filterWhip <- list()
#
#   listfiles <- lapply(listIncTables, function(sample) {
#
#     # Create column for eventID = GENE:NODE:TYPE
#     sample$eventID <- paste0(sample$Gene, ":", sample$Node, ":", sample$Type)
#     # Extract inclusion reads based on Ninc ~ Binomial(n=Ntotal, p=PSI) (See supplementary material from Whippet paper https://doi.org/10.1016/j.molcel.2018.08.018)
#     sample$Ninc <- sample$Psi*sample$Total_Reads
#     # Get exclusion reads based on Ntotal = Ninc + Nexc
#     sample$Nexc <-sample$Total_Reads - sample$Ninc  # or sample$Ninc*(1/sample$Psi - 1)
#
#     # Create auxiliary column that will be needed for the final table
#     sample$length <- rep(0,nrow(sample))
#
#     # Keep only necessary columns
#     sample <- sample[,c(1,15,3,18,3,5,6,16,17)]
#     colnames(sample) <- c("GENE","EVENT","COORD","LENGTH","FullCO","COMPLEX","PSI","INC","EXC")
#
#     return(sample)
#
#   })
#
#   # Change column names so that PSI, INC and EXC are unique for each sample. Will be useful to merge
#   newcolnames <- lapply(names(listfiles), function(x) paste0(c("PSI","INC","EXC"),".",x))
#   for(i in 1:length(listfiles)){
#     colnames(listfiles[[i]])[7:9] <- newcolnames[[i]]
#   }
#
#   # Merge all tables; might take a while, but ensures that the universe of events is the same
#   mergeTable <- as.data.frame(Reduce(function(dtf1, dtf2) merge(dtf1, dtf2,
#                                                                 by.x=c("GENE","EVENT","COORD","LENGTH","FullCO","COMPLEX"),
#                                                                 by.y=c("GENE","EVENT","COORD","LENGTH","FullCO","COMPLEX"),
#                                                                 all = TRUE),
#                                      listfiles))
#
#   # psi matrix
#   colpsis <- grep("PSI.",colnames(mergeTable))
#   psiWhip <- mergeTable[,c(1:6,colpsis)]
#   colpsis <- grep("PSI.",colnames(psiWhip))
#   psiWhip[,colpsis] <- psiWhip[,colpsis]*100
#   colnames(psiWhip)[colpsis] <- gsub("PSI.","",colnames(psiWhip)[colpsis])
#
#   # qual matrix, which includes inc and exc counts
#
#   inc <- as.matrix(mergeTable[,c(grep("INC.",colnames(mergeTable)))])
#   exc <- as.matrix(mergeTable[,c(grep("EXC.",colnames(mergeTable)))])
#
#   qualWhip <- cbind(mergeTable[,c(1:6)],
#                     # Mimicking vast-tools INCLUSION table ".Q" columns to facilitate the compatibility with other betAS functions
#                     matrix(paste0("A,A,0=0=0,A,",mergeTable$COMPLEX,"@", inc , ",", exc ), nrow = nrow(mergeTable[,c(1:6)])))
#   colnames(qualWhip) <- c("GENE","EVENT","COORD","LENGTH","FullCO","COMPLEX",paste0(names(listfiles),".Q"))
#
#
#   filterWhip[[1]] <- psiWhip
#   filterWhip[[2]] <- qualWhip
#   filterWhip[[3]] <- table(psiWhip$COMPLEX)
#   filterWhip[[4]] <- colnames(psiWhip)[-c(1:6)]
#
#   names(filterWhip) <- c("PSI", "Qual", "EventsPerType", "Samples")
#   return(filterWhip)
#
# }
#
#
# # Filter whippet tables for quality (NAs) and event type
# # Filter table from whippet to remove events containing NAs in at least one sample
# # @param WhippetList list containing PSI and Qual tables, as well as event and samples, obtained with getWhippet()
# # @param types (character) whippet's Type column letter code for alternative event types (CE,AA,AD,IR,TS,TE,AF,AL,BS) as described in https://github.com/timbitz/Whippet.jl
# # @param N minimum number of reads (inc + exc) for a given event in each sample
# #
# # @return List with: 1) filtered table PSI columns, 2) filtered table Qual columns, 3) table with number of events per type and 4) Samples (based on colnames)
# # @export
# #
# # @examples
#
# filterWhippet <- function(WhippetList, types, N){
#
#   filterWhippet <- WhippetList
#
#   psiWhip <- filterWhippet$PSI
#   qualWhip <- filterWhippet$Qual
#
#
#   # Remove events containing at least one NA
#   psiWhip$AnyNA    <- apply(psiWhip[,-c(1:6)], 1, anyNA)
#   psiWhip          <- psiWhip[which(psiWhip$AnyNA == FALSE),]
#   psiWhip          <- psiWhip[,-c(ncol(psiWhip))]
#   qualWhip         <- qualWhip[match(psiWhip$EVENT, qualWhip$EVENT),]
#
#
#
#
#   # *** to be made REACTIVE!!! ***
#   # Consider exon skipping events only
#   # psiVAST <- psiVAST[which(psiVAST$COMPLEX == "S"),]
#   psiWhip <- psiWhip[which(psiWhip$COMPLEX %in% types),]
#
#   # Calculate coverage/balance (use .Q columns)
#   qualWhip$AllminReads <- apply(qualWhip[,grep("[.]Q", colnames(qualWhip))], 1, FUN = function(X) Whippet_all_minReads(X,N))
#   qualWhip             <- qualWhip[which(qualWhip$AllminReads == TRUE),]
#   psiWhip              <- psiWhip[match(qualWhip$EVENT, psiWhip$EVENT),]
#
#
#   filterWhippet[[1]] <- psiWhip
#   filterWhippet[[2]] <- qualWhip
#   filterWhippet[[3]] <- table(psiWhip$COMPLEX)
#   filterWhippet[[4]] <- colnames(psiWhip)[-c(1:6)]
#
#   names(filterWhippet) <- c("PSI", "Qual", "EventsPerType", "Samples")
#   return(filterWhippet)
#
# }
#
# # Filter PSI table (whippet) by alternativity
# # Filter previously filtered PSI table from whippet to consider only events with PSIs between (and including) minPsi and maxPSI.
# # @param filteredWhippetList List containing PSI and Qual (whippet) tables, obtained with filterWhippet()
# # @param minPsi (numeric) Minimum PSI to consider
# # @param maxPsi (numeric) Maximum PSI to consider
# #
# # @return List with: 1) filtered table PSI columns, 2) filtered table "Qual" columns, including inc and exc, 3) table with number of events per type and 4) Samples
# # @export
# #
# # @examples
#
#
# alternativeWhippet <- function(filteredWhippetList, minPsi, maxPsi){
#
#   alternativeWhippet <- list()
#
#   psiTable <- filteredWhippetList$PSI
#   qualTable <- filteredWhippetList$Qual
#
#   originalColN <- ncol(psiTable)
#
#   # Consider alternative events only
#   psiTable$AllGreaterMin  <- apply(psiTable[,-c(1:6)], 1, all_grteq_row, minPsi)
#   psiTable$AllLowerMax    <- apply(psiTable[,-c(1:6)], 1, all_lweq_row, maxPsi)
#   psiTable                <- psiTable[which(psiTable$AllGreaterMin == TRUE & psiTable$AllLowerMax == TRUE),]
#   qualTable               <- qualTable[match(psiTable$EVENT, qualTable$EVENT),]
#
#   # Remove columns added
#   psiTable                <- psiTable[,c(1:originalColN)]
#   qualTable               <- qualTable[,c(1:originalColN)]
#
#   alternativeWhippet[[1]] <- psiTable
#   alternativeWhippet[[2]] <- qualTable
#   alternativeWhippet[[3]] <- table(psiTable$COMPLEX)
#   alternativeWhippet[[4]] <- colnames(psiTable)[-c(1:6)]
#
#   names(alternativeWhippet) <- c("PSI", "Qual", "EventsPerType", "Samples")
#
#   return(alternativeWhippet)
#
# }
