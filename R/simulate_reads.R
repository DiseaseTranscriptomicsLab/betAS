# Simulate inclusion and exclusion reads
# Simulate random inclusion and exclusion reads using a binomial distribution where inclusion ~ success, with parameters coverage (number of trials) and  PSI (probability of success).
# @param cov (numeric) coverage (total number of reads)
# @param psi (numeric) percent spliced-in, reflecting the probability of success (i.e., inclusion)
#
# @return List with inclusion and exclusion reads
# @export
#
# @examples
#' @importFrom stats rbinom
simulate_reads <- function(cov, psi){

  simulation <- list()
  bin <- rbinom(n = cov, size = 1, prob = psi)
  inc <- sum(bin)
  exc <- cov - inc

  simulation[[1]] <- inc
  simulation[[2]] <- exc

  names(simulation) <- c("inc", "exc")
  return(simulation)

}
