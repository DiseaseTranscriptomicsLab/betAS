# Coverage (vast-tools)
# Gets coverage from Quality columns (vast-tools) table; requires a table with the first 6 columns and the “.Q” columns from vast-tools’ output table “INCLUSION”, an event name identifying the row and a sample name identifying the column from which to extract the inc, exc pair.
# @param eventname (character) name of the event (table row)
# @param samplename (factor) name of the sample (table column)
# @param vtqualTable (data.frame) vast-tools output INCLUSION table with the 6 identifier columns and “.Q” columns for each sample
#
# @return coverage
# @export
#
# @examples
getQualCov <- function(eventname, samplename, vtqualTable){

  qualCol <- list()

  row <- which(vtqualTable$EVENT == eventname)
  col <- match(samplename, gsub(colnames(vtqualTable), pattern = "[.]Q", replacement = ""))

  qual <- vtqualTable[row,col]

  #Split string by "@" and keep the second element: "inc,exc"
  string  <- as.character(unlist(qual))
  split   <- unlist(strsplit(toString(unlist(string)), split = "[,@]"))

  # Inc reads are the 6th score in each sample; Exc reads the last and 7th
  inc <- as.numeric(as.vector(split[6]))
  exc <- as.numeric(as.vector(split[7]))

  cov <- inc + exc

  qualCol[[1]] <- qual
  qualCol[[2]] <- inc
  qualCol[[3]] <- exc
  qualCol[[4]] <- cov

  names(qualCol) <- c("originalQual", "inc", "exc", "cov")

  return(qualCol)

}
