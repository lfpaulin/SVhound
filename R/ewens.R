.ewen.expected_k <-
function(theta, n){
    return(sum(theta/(theta+(1:n)-1)))
}


.ewen.theta.moments <-
function(frequencies=NULL, occupancies=NULL, k=NULL, n=NULL){
    # both are null, send error message
    if (is.null(k) == TRUE | is.null(n) == TRUE){
        if (is.null(frequencies) == TRUE & is.null(occupancies) == TRUE)
            stop('Either frequencies or occupancies should be provided')
        # both are available use frequencies
        if (!is.null(frequencies) == TRUE | !is.null(occupancies) == TRUE){
            k <- sum(frequencies>0)
            n <- sum(frequencies)
        }
        # if only frequencies is null use occupancies
        if (is.null(frequencies) == TRUE & !is.null(occupancies) == TRUE){
            k <- sum(occupancies)
            n <- sum(apply(X=cbind(1:length(occupancies)), 
                           MARGIN=1, 
                           FUN=function(x){occupancies[x]*x}))
        } else {
            k <- sum(frequencies>0)
            n <- sum(frequencies)
        }
    }
    if (k == 1){
        return (list(theta=0, k=k, n=n))
    }
    # ... implies that K is a sufficient statistic for theta
    # leads to the expression:
    # theta/theta + theta/(theta+1)...theta/(theta+n-1) ~ theta*log(n)
    # we cannot use theta*log(n) because here we need that theta << n and
    # we cannot assure that
    # so we estimate ir from the full formula
    n_operations=25
    epsilon=1e-9
    theta=k
    theta_Lower=0;
    upper_limit_n=1e3
    if (n < upper_limit_n){
        theta_Upper=n*n;
    } else {
        theta_Upper=upper_limit_n*upper_limit_n+n;
    }
    expected_k=0;
    while(abs(expected_k - k) > epsilon){ 
        expected_k <- .ewen.expected_k(theta, n)
        if (expected_k > k) {
            theta_Upper <- theta
            theta <- (theta + theta_Lower)/2
        } else { 
            theta_Lower <- theta
            theta <- (theta + theta_Upper)/2
        }
        n_operations<-n_operations-1; 
        if (n_operations==0){ break }
    }
    if (theta <= 0) theta <- 1e-3
    return (list(theta=theta, k=k, n=n))
}
