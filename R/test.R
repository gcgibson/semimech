

run_semi <- function(a,b,phi){
  i_t <- c()
  i_t[1] <- exp(a+b)
  for (i in 2:100){
    i_t[i] <- exp(a + b*i -  phi*sum(exp(a+b*(1:(i-1)))))
  }
  return (i_t[2:100])
}



plot((run_semi(0,.1,1e-2)),type='l')
