
isVCF <- function(this_file=NULL){
    return (grepl("vcf", this_file, fixed=TRUE) | grepl("VCF", this_file,fixed=TRUE))
}

.actionVCF <- function(vcf_file=NULL, window_size=NULL){
    outname_VCF <- strsplit(vcf_file, "[.]")[[1]][1]
    python_script <- paste("cat ", vcf_file, " | python vcf_parser_for_svhound.py ", window_size, " > ", outname_VCF, ".svalleles", sep="") 
    if(grepl("gz", vcf_file, fixed = TRUE)) python_script <- paste("gzip -dc ", vcf_file, " | python vcf_parser_for_svhound.py ", window_size, " > ", outname_VCF, ".svalleles", sep="") 
    cat("\nPlease use the 'vcf_parser_for_svhound.py' python script included in the package to parse the vcf file into a SV-allele table.\nRun as follows:\n")
    cat("  ", python_script, "\n\n")
    cat("Then use the ", outname_VCF, ".svalleles file as input\n\n", sep="")
}

GRCh37_SV_source <- function(){
    return(list(url=("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz"), readme="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/README_phase3_sv_callset_20150224", web="https://www.internationalgenome.org/phase-3-structural-variant-dataset", pub="https://www.nature.com/articles/nature15394"))
}