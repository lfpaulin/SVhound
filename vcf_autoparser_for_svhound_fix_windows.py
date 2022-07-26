#!/usr/bin/env python3

# Script parse vcf files to annotate the number of disitinct
# SV alleles for a defined window size (bp). It takes the data
# stream from the standard input.
# Getting user specified arguments:
# USAGE:  cat sv_file.vcf | python2 sv_vcf_parser.py <window size>
# EXAMPLES: zcat 1000genome_hg19.vcf.gz | python3 sv_vcf_parser.py  50000
#           zcat 1000genome_hg38.vcf.gz | python3 sv_vcf_parser.py  50000


# important imports
import sys
import gzip
from datetime import datetime
import string
import random


# error output
class SimpleLogger(object):
    def __init__(self, use_as_dev=False):
        self.as_dev = use_as_dev

    @staticmethod
    def warn(msg):
        sys.stderr.write("[WARNING] %s\n" % msg)

    def dev(self, msg):
        if self.as_dev:
            sys.stderr.write("[DEV] %s\n" % msg)

    @staticmethod
    def error(msg):
        sys.stderr.write("[ERROR] %s\n" % msg)
        sys.exit("Program stopped.")

    @staticmethod
    def info(msg):
        sys.stdout.write("[INFO] %s\n" % msg)


# incrementor class
class Incrementor(object):
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
class VCFLine(object):
    # python 101 info:
    #   python array is [start:not_including]
    #                   [:until] from start
    #                   [from:]  until end

    # init
    def __init__(self, input_line):
        vcf_mandatory_fields = 9
        tab_sep_fields = input_line.split("\t")
        # vcf files have 9 fields + each genome data
        if len(tab_sep_fields) > vcf_mandatory_fields:
            # extract mandatory fields
            [self.CHROM, POS, self.ID, _, _, _, self.FILTER, self.INFO, FORMAT] = tab_sep_fields[:vcf_mandatory_fields]
            self.SAMPLES = tab_sep_fields[vcf_mandatory_fields:]
            # position is integer
            self.POS = int(POS)
            self.SVLEN = 0
            self.END = 0
            # pID is a chromosome-position unique id
            # allele (genotype), unique ones and number of appearance
            # for unique alleles
            self.GENOTYPE_UNIQ = {}
            self.GENOTYPE_COUNT = {}
            self.GENOTYPE_1000 = []
            self.get_gt(FORMAT, self.SAMPLES)
            self.parse_info()
            self.pID = f'{self.CHROM}:{self.POS}-{self.END}'
            self.genome_overlap = []
        else:
            # not in the vcf format, used to create empty objects
            self.CHROM = ""
            self.POS = 0
            self.pID = ""
            self.ID = ""
            self.FILTER = ""
            self.GENOTYPE_UNIQ = {}
            self.GENOTYPE_COUNT = {}
            self.GENOTYPE_1000 = []
        self.MERGE_SV = None

    def get_gt(self, _format, genotypes):
        # gives an unique numeric.id to the allele
        g = Incrementor()
        # each genome in the vcf file
        for each_gt in genotypes:
            for fmt, sv_gt in zip(_format.split(":"), each_gt.split(":")):
                if fmt == "GT":
                    if sv_gt == "./." or sv_gt == ".":
                        # use ./. as 0/0
                        sv_gt = "0/0"
                    if sv_gt not in self.GENOTYPE_UNIQ:
                        self.GENOTYPE_UNIQ[sv_gt] = g.current_value()
                        self.GENOTYPE_COUNT[sv_gt] = 1
                        self.GENOTYPE_1000.append(g.current_value())
                        g.increment()
                    else:
                        self.GENOTYPE_1000.append(self.GENOTYPE_UNIQ[sv_gt])
                        self.GENOTYPE_COUNT[sv_gt] += 1

    def parse_info(self):
        # INFO field extraction
        # ##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variation">
        # OR
        # ##INFO=<ID=END,Number=1,Type=Integer,Description="End position of structural variation">
        extract_info = ["SVLEN", "END"]
        if ";" in self.INFO:
            for each_info in self.INFO.split(";"):
                # DP2 unpack too many vals, SKIP IT
                if "DP2" not in each_info:
                    [info_key, info_val] = each_info.split("=") if "=" in each_info else [each_info, ""]
                    if info_key in extract_info:
                        self.SVLEN = int(info_val) if info_key == "SVLEN" else self.SVLEN
                        self.END = int(info_val) if info_key == "END" else self.END
        if self.SVLEN != 0:
            self.END = self.POS + self.SVLEN
        else:
            self.SVLEN = self.END - self.POS

    def add_chromosome(self, chromosome):
        self.CHROM = chromosome

    def add_position(self, position):
        self.POS = position

    def add_pid(self, pid):
        self.pID = pid

    def add_id(self, id_):
        self.ID = id_

    def add_genotype_uniq(self, geno_uniq):
        self.GENOTYPE_UNIQ = geno_uniq

    def add_genotype_count(self, geno_count):
        self.GENOTYPE_COUNT = geno_count

    def add_genotype_list(self, geno_list):
        self.GENOTYPE_1000 = geno_list

    def add_genotype_merge(self, merge_list):
        self.MERGE_SV = merge_list

    def get_genome_overlap(self, window_size):
        overlap_start = int(self.POS/window_size)*window_size + 1
        overlap_ends = max(overlap_start, int(self.END/window_size)*window_size)
        if overlap_start == overlap_ends:
            self.genome_overlap = [overlap_start]
        else:
            self.genome_overlap = range(overlap_start, overlap_ends, window_size)
            # self.genome_overlap = [x for x in range(overlap_start, overlap_ends, window_size)]

    # functions to print sv details
    def sv_print(self):
        sv_alleles = "\t".join([str(g) for g in self.GENOTYPE_1000])
        print(f'{self.pID}\t{sv_alleles}')

    # functions to return sv details, for the case of write to file using open
    def sv_return(self):
        sv_alleles = "\t".join([str(g) for g in self.GENOTYPE_1000])
        return f'{self.pID}\t{sv_alleles}'

    # functions to return sv details, for the case of write to file using open
    def sv_return_k(self):
        n_sv_allele = len(self.ID.split(","))
        return f'{self.pID}\t{len(self.GENOTYPE_COUNT)}\t{len(self.GENOTYPE_1000)}\t{n_sv_allele}\t{self.ID}'


def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))


def make_genome_windows(genome, window_size_bp):
    genome_windows = {}
    for each_chr in genome:
        genome_windows[each_chr] = {}
        for chr_pos in range(1, genome[each_chr], window_size_bp):
            genome_windows[each_chr][chr_pos] = []
    return genome_windows


# merge the allele of two SV per individual
def merge_sv_alleles(sv_list, this_window):
    # SV class properties
    #   self.CHROM is STRING
    #   self.POS is INT
    #   self.ID is STRING
    #   self.GENOTYPE_UNIQ is DICT
    #   self.GENOTYPE_COUNT is DICT
    #   self.GENOTYPE_1000 is LIST
    #   self.MERGE_SV is LIST
    first_elem = 0
    # create new sv object
    sv_merge = VCFLine("")
    # we expect/assume the same chromosome
    # we take all redundant from the first
    # object in the list
    this_chromosome = sv_list[first_elem].CHROM
    sv_merge.add_chromosome(this_chromosome)
    sv_merge.add_position(this_window)
    sv_merge.add_pid(f'{this_chromosome}:{this_window}')
    sv_merge.add_id(",".join([sv.ID for sv in sv_list]))
    n_sv_alleles = len(sv_list[first_elem].GENOTYPE_1000)
    sv_merge.add_genotype_merge([",".join(str(g.GENOTYPE_1000[i]) for g in sv_list) for i in range(n_sv_alleles)])
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


def run_parser(file_handler, analysis_window_size, run_id="", save_out=False, out_prefix=""):
    # go to start of file again (in case)
    file_handler.seek(0)
    # Open stats file
    date_time_obj = datetime.now()
    parser_stats = {}  # id : number of SV-alleles
    outfile_name = ""
    outfile_name_k = ""
    outfile = None
    outfile_k = None
    if save_out:
        window_kb = analysis_window_size/1000
        date_time_str = date_time_obj.strftime("%Y%m%d_%H%M%S")
        out_prefix1 = f'{out_prefix}-svparser-sv_alleles' if out_prefix != "" else "svparser-sv_alleles"
        outfile_name = f'{out_prefix1}-{run_id}-{window_kb}kb-{date_time_str}.tsv'
        out_prefix2 = f'{out_prefix}-svparser-k' if out_prefix != "" else "svparser-k"
        outfile_name_k = f'{out_prefix2}-{run_id}-{window_kb}kb-{date_time_str}.tsv'
        outfile = open(outfile_name, "w")
        outfile_k = open(outfile_name_k, "w")
    # make genome
    genome = {}
    current_chromosome = ""
    [chromosome, pos] = ["", 0]
    for line in file_handler:
        if not line.startswith("#"):
            [chromosome, pos] = line.split("\t")[:2]
            genome[chromosome] = int(pos)
    # for last chromosome
    genome[chromosome] = int(pos)
    genome_windows = make_genome_windows(genome, analysis_window_size)
    # temporary working variables
    sv_list = {}
    current_chromosome = ""
    # go to start
    file_handler.seek(0)
    for line in file_handler:
        # skip comments and in some cases lines of the file
        if line.startswith("#"):
            if save_out and not line.startswith("##"):
                genomes_start_at = 10-1
                outfile.write("\t" + "\t".join(line.rstrip("\n").split("\t")[genomes_start_at:]) + "\n")
        else:
            # use vcf class
            sv = VCFLine(line.rstrip("\n"))
            sv.get_genome_overlap(analysis_window_size)
            sv_list[sv.ID] = sv
            for overlap in sv.genome_overlap:
                genome_windows[sv.CHROM][overlap].append(sv.ID)

    for each_chr in genome_windows:
        for each_window in genome_windows[each_chr]:
            sv_list_merge = []
            for each_sv_id in genome_windows[each_chr][each_window]:
                sv_list_merge.append(sv_list[each_sv_id])
            if len(sv_list_merge) > 0:
                sv_merge = merge_sv_alleles(sv_list_merge, each_window)
                outfile.write(sv_merge.sv_return() + "\n")
                outfile_k.write(sv_merge.sv_return_k() + "\n")
    if save_out:
        outfile.close()
        outfile_k.close()
        return outfile_name

    return parser_stats


def count_sv_per_window(parser_stats):
    # parser_stats is dict
    n_line = 0
    nsv_all = 0
    for key in parser_stats.keys():
        nsv = parser_stats[key]
        nsv_all += nsv
        n_line += 1
    return float(nsv_all)/n_line


def autoparser(file_handler, average_sv_alleles_window=10, out_prefix="",
               fix_window=False, fix_window_size_kb=100,
               min_window_size_kb=10, max_window_size_kb=1000, epsilon_sv_window=0.2):
    # run id
    run_id = id_generator(8)
    # get n lines of fileV
    n_lines = count_lines(file_handler)
    SimpleLogger.info("[%s - %s] %s lines in file" % (run_id, out_prefix, n_lines))
    # window size start in small
    window_size = min_window_size_kb
    min_window_size_kb_update = min_window_size_kb
    max_window_size_kb_update = max_window_size_kb
    window_is_limit = 0
    limit_loops = 10  # failsafe
    window_size_bp = 0
    if fix_window:
        # single
        window_size = fix_window_size_kb  # in kb
        window_size_bp = window_size * 1000
        _ = run_parser(file_handler, window_size_bp, run_id, True, out_prefix)
        SimpleLogger.info("[%s - %s] Current window size = %skb" % (run_id, out_prefix, window_size))
    else:
        # auto parser SV analysis
        window_size_found = False
        window_size_is_limit = False
        while not window_size_found:
            SimpleLogger.info("########## %s ##########" % run_id)
            SimpleLogger.info("[%s - %s] Current window size = %s" % (run_id, out_prefix, window_size))

            # compute number of sv-alleles per window
            window_size_bp = window_size * 1000
            parser_stats = run_parser(file_handler, window_size_bp, run_id)
            current_average_sv_window = count_sv_per_window(parser_stats)

            # redo with diff window size if average sv-alleles per window != X +- er, with er = 0.2
            upp_pass = average_sv_alleles_window + epsilon_sv_window
            low_pass = average_sv_alleles_window - epsilon_sv_window
            if (low_pass < current_average_sv_window < upp_pass) or limit_loops == 0:
                if limit_loops != 0:
                    window_size_found = True
                else:
                    window_size_is_limit = True
                break
            # increase/decrease window size values and update limits.
            # Hard limits of min and max cannot go higher or lower
            if (average_sv_alleles_window - epsilon_sv_window) < current_average_sv_window and not window_size_found:
                max_window_size_kb_update = window_size
                window_size = int((window_size + min_window_size_kb_update)/2)

            if current_average_sv_window < (average_sv_alleles_window + epsilon_sv_window) and not window_size_found:
                min_window_size_kb_update = window_size
                window_size = int((window_size + max_window_size_kb_update)/2)

            # for cases where the window size will be one of the hard limits, if it tries 3 times, break
            # we still have a max of 10 loops to find the window size for the given SV density
            if window_size == max_window_size_kb or window_size == min_window_size_kb:
                window_is_limit += 1
            else:
                window_is_limit = 0

            if window_is_limit >= 3:
                SimpleLogger.info(f'[{run_id} - {out_prefix}] Window size hit the limit three times, '
                                  f'cannot go higher/lower')
                window_size_is_limit = True
                break

            SimpleLogger.info("[%s - %s] Current SVs per window = %s" % (run_id, out_prefix, current_average_sv_window))
            limit_loops -= 1
        # Done
        if window_size_found:
            SimpleLogger.info("[%s - %s] run full" % (run_id, out_prefix))
            fileout = run_parser(file_handler, window_size_bp, run_id, True, out_prefix)
            SimpleLogger.info("[%s - %s] Results are in file: %s" % (run_id, out_prefix, fileout))

        if window_size_is_limit:
            SimpleLogger.info("[%s - %s] run full (window size is in limit)" % (run_id, out_prefix))
            fileout = run_parser(file_handler, window_size_bp, run_id, True, out_prefix)
            SimpleLogger.info("[%s - %s] Results are in file: %s" % (run_id, out_prefix, fileout))

    SimpleLogger.info("[%s - %s] Done" % (run_id, out_prefix))


def buff_count(file_handler):
    lines = 0
    buf_size = 1024 * 1024
    read_f = file_handler.read  # loop optimization
    buf = read_f(buf_size)
    while buf:
        lines += buf.count('\n')
        buf = read_f(buf_size)
    return lines


def simple_count(file_handler):
    lines = 0
    for line in file_handler:
        lines += 1
    return lines


def count_lines(file_handler, method="simple"):
    # number of lines of file
    # source: https://stackoverflow.com/questions/845058/how-to-get-line-count-of-a-large-file-cheaply-in-python
    # good for large files, bad for small ones
    if method == "simple":
        return simple_count(file_handler)
    if method == "buff":
        return buff_count(file_handler)


def main():
    usage = "USAGE:   python3 sv_vcf_parser.py sv_file.vcf <average SV-alleles> <sample size> <out-pref>\n"
    usage += "EXAMPLE: python3 sv_vcf_parser.py 1000genome_sample_hg19.vcf.gz 10 example1"
    # #########################
    # did the user specified a window size
    if len(sys.argv) < 4:
        SimpleLogger.error(f'Wrong number of parameters {usage}')

    # Params
    #   input_file          name of the input VCF file
    #   average_sv_alleles  genome-wide average number of SV-alleles per window
    #   sample_size         number of SV from the VCF file to use for the genome-wide av. number of SV-alleles, 0 is all
    input_file, average_sv_alleles, outfile_prefix = "", "", ""
    dev_fix_the_window_size, dev_fixed_window_size_in_kb = False, 0
    if len(sys.argv) == 4:
        [_, input_file, average_sv_alleles, outfile_prefix] = sys.argv
    elif len(sys.argv) == 5:
        [_, input_file, average_sv_alleles, outfile_prefix, dev_fix_the_window_size] = sys.argv
        dev_fix_the_window_size = dev_fix_the_window_size != 0 and dev_fix_the_window_size != "False"
    elif len(sys.argv) == 6:
        [_, input_file, average_sv_alleles, outfile_prefix, dev_fix_the_window_size, dev_fixed_window_size_in_kb] = sys.argv
        dev_fix_the_window_size = dev_fix_the_window_size != 0 and dev_fix_the_window_size != "False"
        try:
            dev_fixed_window_size_in_kb = int(dev_fixed_window_size_in_kb)
        except ValueError:
            my_logger.warn(f'Could not convert value to integer, 100kb window is used instead')
            dev_fixed_window_size_in_kb = 100
    else:
        SimpleLogger.error(f'Wrong number of parameters {usage}')

    # is file compressed
    if "gz" in input_file:
        file_handler = gzip.open(input_file, "rt")
    else:
        file_handler = open(input_file, "r")
    # run parser
    autoparser(file_handler=file_handler, average_sv_alleles_window=int(average_sv_alleles), out_prefix=outfile_prefix,
               fix_window=dev_fix_the_window_size, fix_window_size_kb=dev_fixed_window_size_in_kb)
    file_handler.close()


# main
as_dev = True
my_logger = SimpleLogger(as_dev)
if __name__ == '__main__':
    main()
