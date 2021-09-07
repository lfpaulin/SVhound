#!/usr/bin/env python

# Script parse vcf files to annotate the number of disitinct
# SV alleles for a defined window size (bp). It takes the data
# stream from the standard input.
# Getting user specified arguments:
# USAGE:  cat sv_file.vcf | python2 sv_vcf_parser.py <window size>
# EXAMPLES: zcat 1000genome_hg19.vcf.gz | python2 sv_vcf_parser.py  50000
#           zcat 1000genome_hg38.vcf.gz | python2 sv_vcf_parser.py  50000


# important imports
import sys
import gzip
from datetime import datetime
import string
import random


# error output
def error_message(text):
    sys.stderr.write("[ERROR] %s \n" % (text))

# error output
def log_message(text):
    print "[LOG] %s" % (text)

# incrementor class
class incrementor(object):
    def __init__(self, starting_value=None):
        if starting_value is None:
            self.value = 0
        else:
            self.value = starting_value
    def increment(self):
        self.value += 1
    def current_value(self):
        return self.value


# vcf line class
class vcf_line(object):
    # python 101 info:
    #   python array is [start:not_including]
    #                   [:until] from start
    #                   [from:]  until end

    # init
    def __init__(self, input_line):
        VCF_MANDATORY_FIELDS = 9
        tab_sep_fields = input_line.split("\t")
        # vcf files have 9 fields + the indiviuals data
        if len(tab_sep_fields) > VCF_MANDATORY_FIELDS:
            #extract mandatory fileds
            [self.CHROM, POS, self.ID, REF, ALT, QUAL, self.FILTER, INFO, FORMAT] = tab_sep_fields[:VCF_MANDATORY_FIELDS]

            # position is integer
            self.POS = int(POS)

            # pID is a chromosom-position unique id
            self.pID = "%s:%s" % (self.CHROM, self.POS)

            # quality is float or '.' when not reported
            # we change '.' for -1
            if (QUAL != "."):
                self.QUAL = float(QUAL)
            else:
                self.QUAL = -1.0

            # allele (genotype), unique ones and number of appearance
            # for unique allels
            self.GENOTYPE_UNIQ = {}
            self.GENOTYPE_COUNT = {}
            self.GENOTYPE_1000 = []
            # gives an unique numeric.id to the allele
            g = incrementor()

            # FORMAT index of GT
            split_format = FORMAT.split(":")
            use_format = False
            # Multiple fields
            if len(split_format) > 1:
                use_format = True
                gt_format_index = split_format.index("GT") # GT == genotype
                #ty_format_index = split_format.index("TY") # TY == type
                #co_format_index = split_format.index("CO") # CO == coordinates

            # indiviuals in the vcf file
            for gt in tab_sep_fields[VCF_MANDATORY_FIELDS:]:
                # use FORMAT duhh!!!
                if use_format:
                    sv_gt = gt.split(":")[gt_format_index]
                    if sv_gt == "./.":
                        # use ./. as 0/0
                        gt = "0/0"
                        # use ./. with type and coordinate for uniqueness
                        # gt = "|".join([
                        #     gt.split(":")[ty_format_index],
                        #     gt.split(":")[co_format_index]
                        # ])
                    else:
                        gt = sv_gt
                if gt not in self.GENOTYPE_UNIQ:
                    self.GENOTYPE_UNIQ[gt] = g.current_value()
                    self.GENOTYPE_COUNT[gt] = 1
                    self.GENOTYPE_1000.append(g.current_value())
                    g.increment()
                else:
                    self.GENOTYPE_1000.append(self.GENOTYPE_UNIQ[gt])
                    self.GENOTYPE_COUNT[gt] += 1
            self.ERROR = False
        else:
            # not in the vcf format, used to create empty objects
            # and get errors
            self.ERROR = True
            self.CHROM = ""
            self.POS = 0
            self.pID = ""
            self.ID = ""
            self.FILTER = ""
            self.QUAL = 0.0
            self.GENOTYPE_UNIQ = {}
            self.GENOTYPE_COUNT = {}
            self.GENOTYPE_1000 = []
        self.MERGE_SV = None

    # functions used to edit the values of each field of the object
    def not_error(self):
        self.ERROR = False

    def add_chromosome(self, chromosome):
        self.CHROM = chromosome

    def add_position(self, position):
        self.POS = position

    def add_pID(self, pID):
        self.pID = pID

    def add_id(self, id_):
        self.ID = id_

    def add_quality(self, qual):
        self.QUAL = qual

    def add_genotype_uniq(self, geno_uniq):
        self.GENOTYPE_UNIQ = geno_uniq

    def add_genotype_count(self, geno_count):
        self.GENOTYPE_COUNT = geno_count

    def add_genotype_list(self, geno_list):
        self.GENOTYPE_1000 = geno_list

    def add_genotype_merge(self, merge_list):
        self.MERGE_SV = merge_list

    # def add_filter(self, filter):
        # self.FILTER = filter

    # functions to print sv details
    def sv_print(self):
        print "%s\t%s" % (self.pID, "\t".join([str(g) for g in self.GENOTYPE_1000]))

    # functions to return sv details, for the case of write to file using open
    def sv_return(self):
        return "%s\t%s" % (self.pID, "\t".join([str(g) for g in self.GENOTYPE_1000]))

# merge the allele of two SV per individual
def merge_sv_alleles(sv_list):

    def mean(l):
        return(sum(l)/float(len(l)))

    # SV class properties
    #   self.CHROM is STRING
    #   self.POS is INT
    #   self.ID is STRING
    #   self.FILTER is STRING
    #   self.QUAL is FLOAT
    #   self.GENOTYPE_UNIQ is DICT
    #   self.GENOTYPE_COUNT is DICT
    #   self.GENOTYPE_1000 is LIST
    #   self.MERGE_SV is LIST
    #   self.ERROR is BOOL

    FIRST_ELEM = 0
    # create new sv object
    sv_merge = vcf_line("")

    # we expect/assume the same chromosome
    # we take all redundant from the first
    # object in the list
    sv_merge.add_chromosome(sv_list[FIRST_ELEM].CHROM)
    sv_merge.add_position(sv_list[FIRST_ELEM].CHROM)
    sv_merge.add_pID(sv_list[FIRST_ELEM].pID)
    sv_merge.add_id(",".join([sv.ID for sv in sv_list]))
    sv_merge.add_quality(mean([sv.QUAL for sv in sv_list]))
    l = len(sv_list[FIRST_ELEM].GENOTYPE_1000)
    sv_merge.add_genotype_merge([ ",".join(str(g.GENOTYPE_1000[i]) for g in sv_list) for i in range(l) ])
    # TODO: make it cleaner
    d = {}
    dc = {}
    svl = []
    g = -1
    for gt in sv_merge.MERGE_SV:
        if gt not in d:
            g += 1
            d[gt] = g
            dc[gt] = 1
            svl.append(g)
        else:
            dc[gt] += 1
            svl.append(g)
    # add important: alleles
    sv_merge.add_genotype_list(svl)
    sv_merge.add_genotype_uniq(d)
    sv_merge.add_genotype_count(dc)
    return sv_merge


def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))


def run_parser(filehander, analysis_window_size, samplefew=0, run_ID="", saveout=False):

    # go to start of file again (in case)
    filehander.seek(0)

    # Open stats file
    dateTimeObj = datetime.now()
    parser_stats = {} # id : number of SV-alleles
    if saveout:
        outfile_name = "svparser-results-%s-%s-%s.tsv" % (run_ID, analysis_window_size, dateTimeObj.strftime("%Y%m%d_%H%M%S"))
        outfile = open(outfile_name, "w")

    # temporary working variables
    rm_duplicates = {}
    sv_within_window = []
    current_chromosome = ""
    current_position = -1
    last_print = {}
    lineN = 0
    skip_sampling = True
    if samplefew != 0:
        skip_sampling = False
    # var used once
    print_header = True

    # get the input from stdin
    for line in filehander:
        # skip comments and in some cases lines of the file
        if not line.startswith("#") and (lineN < samplefew or skip_sampling):
            # use vcf class
            sv = vcf_line(line.rstrip("\n"))

            # check that the vcf line has the proper format
            if not sv.ERROR:
                # single usage to pront the table header
                if print_header:
                    if saveout:
                        outfile.write("\t" + "\t".join(["i%04d"%(i+1) for i in range(len(sv.GENOTYPE_1000))]) + "\n")
                    print_header = False

                # check chromosome and position
                # uninitialized
                if current_chromosome == "":
                    current_chromosome = sv.CHROM

                # different chromosome
                if current_chromosome != sv.CHROM:
                    # print current result
                    if len(sv_within_window) > 1:
                        sv_merge = merge_sv_alleles(sv_within_window)
                        if sv_merge.pID not in last_print:
                            if saveout:
                                outfile.write(sv_merge.sv_return() + "\n")
                            last_print[sv_merge.pID] = 1
                    else:
                        [sv_merge] = sv_within_window
                        if sv_merge not in last_print:
                            if saveout:
                                outfile.write(sv_merge.sv_return() + "\n")
                            last_print[sv_merge.pID] = 1
                    # init new chromosome
                    sv_within_window = []
                    current_chromosome = sv.CHROM
                    current_position = sv.POS
                    sv_within_window.append(sv)
                # same chromosome
                else:
                    # make a dictionary to remove duplicated lines
                    if sv.pID not in rm_duplicates:
                        rm_duplicates[sv.pID] = 1

                        # check the position of the SV
                        # if -1 (un-initialized) then the
                        # curren position is the the sv position
                        if current_position == -1:
                            current_position = sv.POS
                            sv_within_window.append(sv)

                        # in next sv if closer than the given
                        # window size, we merge both
                        else:
                            if sv.POS - current_position <= analysis_window_size:
                                sv_within_window.append(sv)

                            # next sv is father appart from window size
                            # print results and init
                            else:
                                if len(sv_within_window) > 1:
                                    sv_merge = merge_sv_alleles(sv_within_window)
                                    if sv_merge.pID not in last_print:
                                        if saveout:
                                            outfile.write(sv_merge.sv_return() + "\n")
                                        # number of SV-alleles per window in the individuales by window
                                        parser_stats[sv_merge.pID] = len(sv_merge.GENOTYPE_UNIQ)
                                        last_print[sv_merge.pID] = 1
                                else:
                                    [sv_merge] = sv_within_window
                                    if sv_merge.pID not in last_print:
                                        if saveout:
                                            outfile.write(sv_merge.sv_return() + "\n")
                                        # number of SV-alleles per window in the individuales by window
                                        parser_stats[sv_merge.pID] = len(sv_merge.GENOTYPE_UNIQ)
                                        last_print[sv_merge.pID] = 1
                                # init
                                sv_within_window = []
                                current_position = sv.POS
                                sv_within_window.append(sv)

        #add line counter
        lineN += 1

    # on final print for the last line
    if len(sv_within_window) > 1:
        sv_merge = merge_sv_alleles(sv_within_window)
        if sv_merge.pID not in last_print:
            if saveout:
                outfile.write(sv_merge.sv_return() + "\n")
            # number of SV-alleles per window in the individuales by window
            parser_stats[sv_merge.pID] = len(sv_merge.GENOTYPE_UNIQ)
            last_print[sv_merge.pID] = 1
    else:
        if len(sv_within_window) == 1:
            [sv_merge] = sv_within_window
            if sv_merge.pID not in last_print:
                if saveout:
                    outfile.write(sv_merge.sv_return() + "\n")
                # number of SV-alleles per window in the individuales by window
                parser_stats[sv_merge.pID] = len(sv_merge.GENOTYPE_UNIQ)
                last_print[sv_merge.pID] = 1

    if saveout:
        outfile.close()
        return outfile_name

    return parser_stats


def count_SVperWindow(parserStats):
    # parserStats is dict
    nsv = 0
    nline = 0
    nsv_all = 0
    for key in parserStats.keys():
        nsv = parserStats[key]
        nsv_all += nsv
        nline += 1
    return(float(nsv_all)/nline)


def autoparser(filehander, averageSVperWindow=10, minWindowSizeKb=10, maxWindowSizeKb=1000, espilonSVperWindow=0.2, sampleSize=0):

    # here we have the results
    SVperWindow = {}

    # run id
    run_ID = id_generator(8)

    # get n lines of fileV
    nlines = countlines(filehander)
    log_message("[%s] %s lines in file" % (run_ID, nlines))

    # sample some lines and make windows
    # samplefew = 1
    # if sampleSize != 0:
    #     samplefew = int(nlines/sampleSize)

    # window size start in small
    window_size = minWindowSizeKb
    current_averageSVwindow = 0

    testme = False
    if testme:
        # single
        window_size_bp = window_size * 1000
        parserStats = run_parser(filehander, window_size_bp, sampleSize, run_ID)
        current_averageSVwindow = count_SVperWindow(parserStats)
        log_message("[%s] Current window size = %s" % (run_ID, window_size))
        log_message("[%s] Current SVs per window = %s" % (run_ID, current_averageSVwindow))
    else:
        # auto parser SV analysis
        window_size_found = False
        while not window_size_found:
            log_message("########## %s ##########" % run_ID)
            log_message("[%s] Current window size = %s" % (run_ID, window_size))

            # compute number of sv-alleles per window
            window_size_bp = window_size * 1000
            parserStats = run_parser(filehander, window_size_bp, sampleSize, run_ID)
            current_averageSVwindow = count_SVperWindow(parserStats)

            # redo with diff window size if average sv-alleles per window != X +- er, with er = 0.2
            if (averageSVperWindow - espilonSVperWindow) < current_averageSVwindow < (averageSVperWindow + espilonSVperWindow):
                window_size_found = True
                break

            # TODO how to increase/decrease window size values
            if (averageSVperWindow - espilonSVperWindow) < current_averageSVwindow and not window_size_found:
                maxWindowSizeKb = window_size
                window_size = int((window_size + minWindowSizeKb)/2)

            if current_averageSVwindow < (averageSVperWindow + espilonSVperWindow) and not window_size_found:
                minWindowSizeKb = window_size
                window_size = int((window_size + maxWindowSizeKb)/2)

            log_message("[%s] Current SVs per window = %s" % (run_ID, current_averageSVwindow))

        # Done
        if window_size_found:
            log_message("[%s] run full" % run_ID)
            fileout = run_parser(filehander, window_size_bp, 0, run_ID, True)
            log_message("[%s] Results are in file: %s" % (run_ID, fileout))

    log_message("[%s] Done" % run_ID)


# number of lines of file
# source: https://stackoverflow.com/questions/845058/how-to-get-line-count-of-a-large-file-cheaply-in-python
# good for large files, bad for small ones
def bufcount(filehander):
    lines = 0
    buf_size = 1024 * 1024
    read_f = filehander.read # loop optimization

    buf = read_f(buf_size)
    while buf:
        lines += buf.count('\n')
        buf = read_f(buf_size)

    return lines

def simplecount(filehander):
    lines = 0
    for line in filehander:
        lines += 1
    return lines

def countlines(filehander, method="simple"):
    if method == "simple":
        return simplecount(filehander)
    if method == "buff":
        return bufcount(filehander)



usage = """
USAGE:   python2 sv_vcf_parser.py sv_file.vcf <average SV-alleles> <sample size>
EXAMPLE: python2 sv_vcf_parser.py 1000genome_sample_hg19.vcf.gz 10 1000
"""

# main
if __name__ == '__main__':

    # #########################
    # did the user specified a window size
    if len(sys.argv) != 4:
        error_message("Wrong number of parameters. %s" % usage)
        sys.exit("Program halt.")

    # Params
    #   filename     name of the input VCF file
    #   avSValleles  genome-wide average number of SV-alleles per window
    #   sampleSize   number of SV from the VCF file to use for the genome-wide average number of SV-alleles, 0 is all
    [script_name, filename, avSValleles, sampleSize] = sys.argv

    # is file compressed
    if "gz" in filename:
        filehander = gzip.open(filename, "r")
    else:
        filehander = open(filename, "r")

    # run parser
    autoparser(filehander=filehander, averageSVperWindow=int(avSValleles), sampleSize=int(sampleSize))
    filehander.close()
