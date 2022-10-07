# Define increment for inclusion and exclusion reads
# Based on previous coverage/inclusion dependence study, use coverage to estimate a small increment to inc and exc so that PSI is not altered in more than 1%
# @param cov read coverage (inc + exc)
# @param maxDevRef table with previous dependence between coverage and PSI levels estimating the largest increment that ensures the PSI obtained with corrected reads (incremented) does notvary by more than 1%
#
# @return Increment
# @export
#
givenCovfindIncr <- function(cov, maxDevRef){

  if(cov > max(maxDevRef[,1])){

    incr <- 1

  }else{

    p <- which(cov == maxDevRef$cov)
    incr <- maxDevRef$maxIncr[p]

  }

  return(incr)

}
