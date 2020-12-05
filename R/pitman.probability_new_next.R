pitman.probability_new_next <-
function(theta, sigma, k, n){
    return( (theta + (k*sigma)) / (theta + n) )
}
