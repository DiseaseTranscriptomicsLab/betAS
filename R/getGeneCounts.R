# Gene counts (vast-tools)
# Get gene counts from Counts columns (vast-tools) table; requires a table with the first 2 columns and the “.Counts” columns from vast-tools’ output table “cRPKM_AND_COUNTS”, an event name identifying the row and a sample name identifying the column from which to extract the counts.
# @param genename (character) name of the gene (table row)
# @param samplename (factor) name of the sample (table column)
# @param vtcountsTable (data.frame) vast-tools output counts table with the 2 identifier columns and “.Counts” columns for each sample
#
# @return Gene counts
# @export
#
# @examples
getGeneCounts <- function(genename, samplename, vtcountsTable){

  row <- which(vtcountsTable$NAME == genename)
  col <- match(samplename, gsub(colnames(vtcountsTable), pattern = ".Counts", replacement = ""))

  genecounts <- vtcountsTable[row,col]

  if(length(genecounts) != 1){genecounts <- NA}

  return(genecounts)

}
