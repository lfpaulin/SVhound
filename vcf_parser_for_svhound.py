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


# error output
def error_message(text):
    print "[ERROR] %s" % text


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
                ty_format_index = split_format.index("TY") # TY == type
                co_format_index = split_format.index("CO") # CO == coordinates

            # indiviuals in the vcf file
            for gt in tab_sep_fields[VCF_MANDATORY_FIELDS:]:
                # use FORMAT duhh!!! 
                if use_format:
                    tmp = gt.split(":")[gt_format_index]
                    if tmp == "./.":
                        gt = "|".join([
                            gt.split(":")[ty_format_index],
                            gt.split(":")[co_format_index]
                        ])
                    else:
                        gt = tmp
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
                        

usage = """
USAGE:  cat sv_file.vcf | python2 sv_vcf_parser.py <window size>
EXAMPLE: zcat 1000genome_hg19.vcf.gz | python2 sv_vcf_parser.py  50000
"""

# main
if __name__ == '__main__':
    # DEFAULT VALUES ##########
    
    #   Default is 10kb
    DEFAULT_WINDOW_SIZE = 10000
    
    # #########################
    
    # used once
    print_header = True
    
    # did the user specified a window size
    if len(sys.argv) != 2:
        print usage
        sys.exit("Program halt. Wrong number of parameters")
    [script_name, analysis_window_size] = sys.argv
    analysis_window_size = int(analysis_window_size)
    
    # temporary working variables
    rm_duplicates = {}
    sv_within_window = []
    current_chromosome = ""
    current_position = -1
    last_print = {}
    
    # get the input from stdin
    for line in sys.stdin:
        # skip comments
        if not line.startswith("#"):
            # use vcf class
            sv = vcf_line(line.rstrip("\n"))
            
            # check that the vcf line has the proper format
            if not sv.ERROR:
                # single usage to pront the table header
                if print_header:
                    print "\t" + "\t".join(["i%04d"%(i+1) for i in range(len(sv.GENOTYPE_1000))])
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
                            sv_merge.sv_print()
                            last_print[sv_merge.pID] = 1
                    else:
                        [sv_merge] = sv_within_window
                        if sv_merge not in last_print:
                            sv_merge.sv_print()
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
                                        sv_merge.sv_print()
                                        last_print[sv_merge.pID] = 1
                                else:
                                    [sv_merge] = sv_within_window
                                    if sv_merge.pID not in last_print:
                                        sv_merge.sv_print()
                                        last_print[sv_merge.pID] = 1
                                # init
                                sv_within_window = []
                                current_position = sv.POS
                                sv_within_window.append(sv)
    
    # on final print for the last line
    if len(sv_within_window) > 1:
        sv_merge = merge_sv_alleles(sv_within_window)
        if sv_merge.pID not in last_print:
            sv_merge.sv_print()
            last_print[sv_merge.pID] = 1
    else:
        if len(sv_within_window) == 1:
            [sv_merge] = sv_within_window
            if sv_merge.pID not in last_print:
                sv_merge.sv_print()
                last_print[sv_merge.pID] = 1
