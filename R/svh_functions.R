# ##################################################################################################
# ### svhound #######################################################################################
# ##################################################################################################

apply_VERTICAL_MARGIN = 2
apply_HORIZONTAL_MARGIN = 1


sv_analysis_ESF <- function(sv_dataset=NULL, outprefix=NULL){
# ############################################## #
# Analysis of SV with the ESF                    #
#                                                #
# Genome-wide analysis to detect genomic regions #
# in which new alleles (SV) could be detected    #
# with high probability based on calculations    #
# with the ESF                                   #
# ############################################## #

    if(is.null(sv_dataset)) stop("the data-set to analyze is missing")
    if (is.null(outprefix)) outprefix="results"

    svhound_results <- list()
    nW <- nrow(sv_dataset)
    N <- ncol(sv_dataset)
    # object to save the results
    
    #  k and unique k
    k <- apply(sv_dataset, 
                apply_HORIZONTAL_MARGIN, 
                function(r) length(unique(r)) 
                )
    uniq_k <- sort(unique(as.vector(k)))

    #  theta and p new sv next individual
    theta <- sapply(uniq_k, function(x) .ewen.theta.moments(k=x, n=N)$theta)
    pnew  <- sapply(theta,  function(x) ewens.probability_new_next(x, n=N))

    # compile all the results
    sv_esf_pnew <- data.frame(k=k,
                    theta=sapply(k, function(x) theta[uniq_k == x]),
                    pnewEWENS=sapply(k, function(x) pnew[uniq_k == x])
    )
    svhound_results <- sv_esf_pnew$pnewEWENS
    names(svhound_results) <- rownames(sv_dataset)

    # save results to object
    outfile <- paste(outprefix, "-svhound_analysis.RData", sep="")
    svh_results <- svhound_results # new name everywhere in downstream analysis
    save(svh_results, file=outfile)
    return(c(results=outfile))
}

sv_analysis_ESF_window <- function(window=NULL, windowname=NULL, saved_results=NULL){
# ############################################## #
# Analysis of SV with the ESF                    #
#                                                #
# Genome-wide analysis to detect genomic regions #
# in which new alleles (SV) could be detected    #
# with high probability based on calculations    #
# with the ESF, for a single window              #
# ############################################## #

    if(is.null(window)) stop("the data to analyze is missing")
    if(is.null(windowname)) windowname <- "Window"

    svhound_results <- list()
    N <- length(window)
    # object to save the results
    
    # k
    k <- length(unique(window))

    prev_calc_pnew <- c()
    if (!is.null(saved_results)) {
        prev_calc_pnew <- saved_results$pnew[saved_results$k == k]
    } 
    
    if (length(prev_calc_pnew) != 0){
        pnew <- prev_calc_pnew
    } else {
        # theta and p new sv next individual
        theta <- .ewen.theta.moments(k=k, n=N)$theta
        pnew  <- ewens.probability_new_next(theta, n=N)
        saved_results$k <- c(saved_results$k, k)
        saved_results$pnew <- c(saved_results$pnew, pnew)
    }
    names(pnew) <- windowname

    # return results to object
    return(list(result=pnew, saved_results=saved_results))
}

sv_analysis_PSF <- function(sv_dataset=NULL){
# ############################################## #
# This is with the Pitman Sampling Formula       #
# to say that the results do not substantially   #
# improve to be worth the extra compiting time   #
# ############################################## #

    if(is.null(sv_dataset)) stop("the data-set to analyze is missing")

    nW <- nrow(sv_dataset)
    N <- ncol(sv_dataset)

    freq_vectors <- apply(sv_dataset, apply_HORIZONTAL_MARGIN, function(r) {
                        paste(sort(as.numeric(table(as.numeric(r)))), collapse=",")
                        })
    k <- apply(sv_dataset, 
                apply_HORIZONTAL_MARGIN, 
                function(r) length(unique(r)) 
                )
    uniq_k <- sort(unique(as.vector(k)))
    names(k) <- names(freq_vectors)

    # ordenar y solo usar los unicos frequency vectors
    # para los occupancy no importa cuales alelos son sino
    # cuantos aparecen
    uniq_freq_vectors <- unique(freq_vectors)
    uniq_freq_vectors_k <- sapply(strsplit(uniq_freq_vectors, ","), length)

    # calcular pnew con PSF, vemos un aumento en performance
    pitman_sv <- t(sapply(uniq_k, function(k) {
        if (k == 1) {
            c(k=k, theta=0, sigma=0, p_new = 0)
        } else {
            # fv is freq vector
            # fm is freq matrix
            fv = uniq_freq_vectors[uniq_freq_vectors_k == k]
            fm = matrix(as.integer(unlist(strsplit(fv, ","))), nrow=k)
            c(k=k, .pitman_ML_infinite_matrix(theta=1, sigma=0, freq=fm))
        }
    }))

    sv_psf_pnew <- data.frame(
        k=k,
        theta=sapply(k, function(x) pitman_sv[pitman_sv[,"k"] == x, "theta"]),
        sigma=sapply(k, function(x) pitman_sv[pitman_sv[,"k"] == x, "sigma"]),
        pnewPITMAN=sapply(k, function(x) pitman_sv[pitman_sv[,"k"] == x, "p_new"])
    )

    # save results to object
    svhound_results <- sv_psf_pnew$pnewPITMAN
    names(svhound_results) <- rownames(sv_dataset)
    svhound_metadata <- pitman_sv
    save(svhound_results, svhound_metadata, file="svhound_analysis_psf.RData")
    return(c(results="svhound_analysis_psf.RData"))
}


# ######### SUBSAMPLE ############################

subsample_analysis_ESF <- function(sv_dataset=NULL, subsample=NULL){
# ############################################## #
# Subsample analysis of SV with the ESF          #
#                                                #
# Genome-wide analysis to detect genomic regions #
# in which new alleles (SV) could be detected    #
# with high probability based on calculations    #
# with the ESF                                   #
# ############################################## #

    # number of individuals subsampled
    if(is.null(subsample)) stop("[ERRROR] subsample is not given")
    if(is.null(sv_dataset)) stop("[ERROR] data was not given")

    # object to save the results
    svhound_results <- list()
    nW <- nrow(sv_dataset)
    N <- ncol(sv_dataset)

    samp_ind <- lapply(subsample, function(s) {
                    sample(ncol(sv_dataset), s)
    })

    names(samp_ind) <- paste("s", as.character(subsample), sep="")

    k <- lapply(samp_ind, function(s) {
            apply(sv_dataset[, s], apply_HORIZONTAL_MARGIN, function(r) {
                    length(unique(r))
            })
        })
    uniq_k <- lapply(k, function(s){
        sort(unique(as.vector(s)))
    })

    pnewDATA <- lapply(samp_ind, function(s) {
                    apply(sv_dataset, apply_HORIZONTAL_MARGIN, function(r) {
                    S=r[s]; # alleles from sampled individuals
                    P=r[-s]; # alleles from NON-sampled individuals
                    sum(match(P, S, 0L) == 0L); # check those that NOT-match
                    }) / (N-length(s))
                })

    #  theta and p new sv next individual
    theta <- uniq_k
    pnew <- uniq_k
    for (i in seq(1, length(k))){
        theta[[i]] <- sapply(uniq_k[[i]], function(x) .ewen.theta.moments(k=x, n=subsample[i])$theta)
        pnew[[i]] <- sapply(theta[[i]], function(x) ewens.probability_new_next(x, n=subsample[i]))
    }

    # compile all the results
    results = list()
    for (i in seq(1, length(k))){
        results[[i]] <- data.frame(k=k[[i]],
                    theta=k[[i]],
                    pnewEWENS=k[[i]],
                    pnewDATA=pnewDATA[[i]]
        )
        for(uk in uniq_k[[i]]){
            results[[i]]$theta[ results[[i]]$k == uk ] = theta[[i]][ match(uk, uniq_k[[i]]) ]
            results[[i]]$pnewEWENS[ results[[i]]$k == uk ] = pnew[[i]][ match(uk, uniq_k[[i]]) ]
        }
    }
    names(results) <- paste("s", as.character(subsample), sep="")

    svhound_results = results 
    svhound_metadata = samp_ind
    save(svhound_results, svhound_metadata, file="svhound_analysis_subsample.RData")
    return(c(results="svhound_analysis_subsample.RData"))
}


subsample_analysis_PSF <- function(sv_dataset=NULL, subsample=NULL){
# ############################################## #
# This is with the Pitman Sampling Formula       #
# to say that the results do not substantially   #
# improve to be worth the extra compiting time   #
# ############################################## #

    # data-set    
    if(is.null(sv_dataset)) stop("[ERROR] data was not given")

    # object to save the results
    svhound_results <- list()
    nW <- nrow(sv_dataset)
    N <- ncol(sv_dataset)

    # number of individuals subsampled or previously used ones
    if(is.null(subsample)) stop("[ERRROR] subsample/esf_samp_ind is not given")
    samp_ind <- lapply(
        subsample, function(s) {
            sample(ncol(sv_dataset), s)
    })

    # sample individuals
    freq_vectors <- lapply(samp_ind, function(s) {
                        apply(sv_dataset, apply_HORIZONTAL_MARGIN, function(r) {
                        A=r[s]; # alleles from sampled individuals
                        # frequency vector
                        paste(sort(as.numeric(table(as.numeric(A)))), collapse=",")
                        })
                    })
    k <- lapply(freq_vectors, function(f) { sapply(strsplit(f, ","), length)  })
    # ordenar y solo usar los unicos frequency vectors
    # para los occupancy no importa cuales alelos son sino
    # cuantos aparecen
    uniq_freq_vectors <- lapply(freq_vectors, function(f) {
        ufv=unique(f)
        uk=sapply(strsplit(ufv, ","), length)
        data.frame(vector=ufv, k=uk)
    })

    # calcular pnew con PSF, vemos un aumento en performance
    pitman_sv <- lapply(uniq_freq_vectors, function(f){
        uniq_k = unique(f$k)
        t(sapply(uniq_k, function(k){
            if (k == 1) {
                c(k=k, theta=0, sigma=0, p_new = 0)
            # nore than one allele in the sample
            } else {
                c(k=k, .pitman_ML_infinite_matrix(theta=1, sigma=0, freq=
                    matrix(as.integer(unlist(strsplit(as.character(f[f$k==k,"vector"]), ","))), nrow=k)
                ))
            }
        } ))
    })
        
    # colectar los resultados
    results <- k
    for (nam in names(pitman_sv)){
        results[[nam]] <- t(
            sapply(k[[nam]], function(x) {
                pitman_sv[[nam]][pitman_sv[[nam]][,"k"] == x, ]
            })
        )
    }

    # save results to object
    svhound_results = results
    svhound_metadata = list(pitman = pitman_sv, subsample = samp_ind)
    save(svhound_results, svhound_metadata, file="svhound_analysis_subsample_psf.RData")
    return(c(results="svhound_analysis_subsample_psf.RData"))
}
