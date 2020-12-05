.Gradient_LogPitman_single_infinite <-
function(param,freq,k,n) {
  if (sigma < 0 | sigma >= 1 | theta <= -sigma | theta < 0) {
    return(NA)
  } else {
    d <- c(0,0)
    d[1] <- sum(1/(theta + c(1:(k-1))*sigma)) - sum(1/(theta+c(1:(n-1))))
    d[2] <- sum(c(1:(k-1))/(param[1] + c(1:(k-1))*param[2]))
    for(j in which(freq>1)) {
      d[2] <- d[2] + sum(1/(param[2]-c(1:(freq[j]-1))))
    }
    return(d)
  }
}


.LogPitman_single_infinite_multi <-
function(param, freq, k, n) {
  freq_2_occup <- function(freq){
    a <- as.data.frame(table(freq))
    colnames(a) <- c("Occup", "Freq")
    return (a)
  }

  theta = param[1]
  sigma = param[2]
  if (sigma < 0 | sigma >= 1 | theta <= -sigma | theta < 0) {
    return(NA)
  } else {
    nsample = ncol(freq)
    a = apply(freq, 2, freq_2_occup) # list X elements one per "column" in freq
    return(sum(sapply(a, function(occup) {
            lfactorial(n) + 
            sum(log(theta + (c(1:(k-1)*sigma)) )) -
            sum(log((theta + 1) + c(0:(n-2)))) +
            sum( apply(occup, 1, function(x) {
                        j = as.numeric(x[1])
                        aj = as.numeric(x[2])
                        if (j < 3) {
                          return( 0 - (lfactorial(j)*aj + lfactorial(aj)) )
                          } else {
                            return( sum(log(1-sigma + c(1:(j-2)))*aj -
                                        (lfactorial(j)*aj + lfactorial(aj))
                                       )
                                   )
                          } 
                        }
            ))
            })
          ))
  }
}


.pitman_ML_infinite_matrix <-
function(theta, sigma, freq){
        k = nrow(freq)
        n = colSums(freq)[1]
        param = optim(par=c(theta, sigma), 
                    fn = .LogPitman_single_infinite_multi, 
                    # gr = .Gradient_LogPitman_single_infinite, 
                    # method = "L-BFGS-B", 
                    # lower = c(0, Inf),
                    # upper = c(0,1), 
                    freq = freq, 
                    k = k, 
                    n = n, 
                    control = list(fnscale=-1))$par
        theta = param[1]
        sigma = param[2]
        p_new = pitman.probability_new_next(theta, sigma, k, n)
        return(c(theta=theta, sigma=sigma, p_new=p_new))
}
