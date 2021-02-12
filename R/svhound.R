# ##################################################################################################
# ### SVHOUND ######################################################################################
# ##################################################################################################

svhound <- function(structuralVariantsDataFile=NULL, SVallelesTable=NULL, window_size=NULL, output_prefix=NULL, subsample=NULL, usePSF=FALSE, giveExampleData=FALSE, runExample=FALSE, onlypnew=TRUE){
# ############################################## #
# Analysis of SV with the ESF                    #
# Wreapper for the analysis that includes a      #
# dataset example.
# ############################################## #
    # first wants run-example
    if (runExample)  return(sv_analysis_ESF(sv_dataset = sv_data_example, outprefix = "Example_dataset"))

    # second wants example-data
    if (giveExampleData) return(sv_data_example)

    # no examples, check for mandatory parameters then
    is_data_provided <- is.null(structuralVariantsDataFile) & is.null(SVallelesTable)
    data_is_file <- !is.null(structuralVariantsDataFile)
    if (is_data_provided | is.null(window_size) ) stop("The input data and/or window size are missing. Set the 'giveExampleData' TRUE to get an example. Data can be: 1. the path of a VCF file (structuralVariantsDataFile), 2. the path of the a SV-allele table in plain text separated by 'spaces' or 'tabs' (structuralVariantsDataFile), or 3. the object of the SV-allele table (SVallelesTable). Window size is in kilobases.")
    window_size <- as.integer(window_size)

    # pre-cases
    # VCF is given give instructions
    if (data_is_file){
        if (isVCF(structuralVariantsDataFile)) {
            .actionVCF(structuralVariantsDataFile, window_size)
            return(0)
        }
        # SV-allele table is given
        SVallelesTable <- read.table(structuralVariantsDataFile, header = TRUE)
    }

    # check the data is like a matrix/data.frame
    cols <- ncol(SVallelesTable)
    rows <- nrow(SVallelesTable)
    if (is.null(ncol) | is.null(rows)) stop("The input data has the wrong format. Set the 'giveExampleData' TRUE to get an example.")

    # Optional parameters 4 cases:
    #   esf
    #   psf
    #   esf + subsaple
    #   psf + subsample
     
    # CASES 1 and 2: no subsample
    if (is.null(subsample)){
        # PTIMAN SAMPLING FORMULA
        if(usePSF){
            return(sv_analysis_PSF(sv_dataset = SVallelesTable, outprefix=output_prefix))
        # EWENS SAMPLING FORMULA
        } else {
            return(sv_analysis_ESF(sv_dataset = SVallelesTable, outprefix=output_prefix, only_pnew=onlypnew))
        }
    } else {
    # CASES 3 and 4: do subsample
        # subsample must be integer
        subsample = as.integer(subsample)
        # subsample must be smaller than total number of individuals
        for(s in subsample){
            if (s == 0) stop("The size of the subsample need to be larger than zero")
            if (ncol(SVallelesTable) < s) stop("The size of the subsample need to be smaller than the total of individuals in the sample")
        }
        # PITMAN SAMPLING FORMULA
        if(usePSF){
            return(subsample_analysis_PSF(sv_dataset = SVallelesTable, subsample = subsample))
        } else {
        # EWENS SAMPLING FORMULA
            return(subsample_analysis_ESF(sv_dataset = SVallelesTable, subsample = subsample))
        }

    }
}
