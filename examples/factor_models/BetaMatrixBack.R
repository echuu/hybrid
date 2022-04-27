# a function that takes a vector of non zero entries
# and put it back into a m by k matrix

#vec is assumed to be a c(beta) of an original m by k beta matrix 
#, but without the zeroes
beta.back <- function(m, k, vec){
  
  mat  = matrix(0, nr = m,nc = k)
  vec.pos  = 0
  for(i in 1:k){
    mat[i:m ,i] = vec[(vec.pos+1):(vec.pos + m-i+1)]  
    vec.pos = vec.pos + m-i+1
  }
  return(mat)
}
