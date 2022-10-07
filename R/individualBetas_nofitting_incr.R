# Individual betAS
# Generates indpoints from a beta distribution for a group of samples, with shape parameters based on inc and exc from Qual columns (vast-tools) table; requires a table with the first 6 columns and the “.Q” columns from vast-tools’ output table “INCLUSION”, the columns (cols) to consider; requires a maximum deviation reference table maxdevRefTable to be used to increment inc or exc if any  is equal to zero.
# @param table (data.frame) row from vast-tools output INCLUSION table with the 6 identifier columns and “.Q” columns for each sample. Function is applied per row of INCLUSION table using pblapply() or similar.
# @param cols (integer) vector with columns of table to consider
# @param indpoints (integer) number of points/observations emitted from the beta distribution
# @param maxdevRefTable (data.frame) reference data frame with maximum increment values per coverage (inc+exc) to avoid emitted values to be artifitially beyond a certain threshold
#
# @return indBetas object list with, for event: 1) inc = number of inclusion reads, 2) exc = number of exclusion reads, 3) PSIs, 4 and 5) corrected inc and exc (using a small increment to avoid zeros), 6) corrected PSIs  (based on corrected inc and exc), 7) points emitted from beta distribution, 8) mean and 9) median of the distribution.
# @export
#
# @examples
individualBetas_nofitting_incr <- function(table, cols, indpoints, maxdevRefTable){

  indBetas <- list()

  # .Q columns to consider
  quals <- table[1,cols]
  # Number of samples
  n <- length(cols)
  # Number of points in individual Betas
  no_points <- indpoints
  # Sample names
  names <- colnames(table)[cols]

  # Split string by "@" and keep the second element: "inc,exc"
  string <- as.character(unlist(quals))
  split <- unlist(strsplit(toString(unlist(string)), split = "[,@]"))

  # Inc reads are the 6th score in each sample; Exc reads the last and 7th
  inc <- as.numeric(as.vector(split[seq(from = 6, to = 7*n, by = 7)]))
  exc <- as.numeric(as.vector(split[seq(from = 7, to = 7*n, by = 7)]))

  names(inc) <- names
  names(exc) <- names
  indBetas[[1]] <- inc
  indBetas[[2]] <- exc

  # Keep original PSI values
  psis <- inc/(inc+exc)
  names(psis) <- names
  indBetas[[3]] <- psis

  for(i in 1:n){
    # Sum incr only if read number is equal to zero
    if(inc[i] == 0){inc[i] <- givenCovfindIncr(cov = round(exc[i]), maxDevRef = maxdevRefTable)}
    if(exc[i] == 0){exc[i] <- givenCovfindIncr(cov = round(inc[i]), maxDevRef = maxdevRefTable)}
  }

  names(inc) <- names
  names(exc) <- names
  indBetas[[4]] <- inc
  indBetas[[5]] <- exc

  # Keep increment-corrected PSI values
  incrpsis <- inc/(inc+exc)
  names(incrpsis) <- names
  indBetas[[6]] <- incrpsis

  # Create individual artificial Betas and concatenate the n distributions
  artif <- c()
  artif_sample <- c()
  for(i in 1:length(inc)) {
    artif_sample <- c(artif_sample, rep(names[i], times = indpoints))
    artif <- c(artif, rbeta(no_points,
                            shape1 = inc[i],
                            shape2 = exc[i]))
  }

  names(artif) <- artif_sample

  indBetas[[7]] <- artif

  #Method of Moments
  mean_x <- mean(artif)
  variance_x <- sd(artif)^2
  median_x <- median(artif)

  indBetas[[8]] <- mean_x
  indBetas[[9]] <- median_x

  names(indBetas) <- c("inc", "exc", "PSI", "tinc", "texc", "incrPSI", "BetaPoints", "MeanBeta", "MedianBeta")
  class(indBetas) <- c(class(indBetas), "indBetas")
  return(indBetas)

}
