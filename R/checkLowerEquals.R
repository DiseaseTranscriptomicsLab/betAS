# Lower or equal than k
# Checks if test is lower or equal than k.
# @param test Number to check if is lower or equal than k.
# @param value Value k to compare test with.
#
# @return TRUE if test is lower or equal than k.
#
# @examples
lweq <- function(test, value){
  r <- test <= value
  return(r)
}

# Lower or equal than (vector)
# Applies lweq() over a vector to check whether all its values are lower or equal than comparison value k.
# @param test_vec Vector to test if all entries are lower or equal than value.
# @param value Value to compare to.
#
# @return TRUE if all vector entries are lower or equal than value.
#
# @examples
all_lweq_row <- function(test_vec, value){
  r <- sapply(test_vec, lweq, value)
  test <- any(r == FALSE)
  bool <- !test
  return(bool)
}

# Greater or equal than k
# Checks if test is greater or equal than k.
# @param test Number to check if is greater or equal than k.
# @param value Value k to compare test with.
#
# @return TRUE if test is greater or equal than k.
#
# @examples
grteq <- function(test, value){
  r <- test >= value
  return(r)
}


# Greater or equal than (vector)
# Applies grteq() over a vector to check whether all its values are greater or equal than comparison value k.
# @param test_vec Vector to test if all entries are greater or equal than value.
# @param value  Value to compare to.
#
# @return TRUE if all vector entries are greater or equal than value.
#
# @examples
all_grteq_row <- function(test_vec, value){
  r <- sapply(test_vec, grteq, value)
  test <- any(r == FALSE)
  bool <- !test
  return(bool)
}


