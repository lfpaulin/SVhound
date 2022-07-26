#!/usr/bin/env python3

# Script parse vcf  and detect hotspots
# USAGE:  cat sv_file.vcf | python3 hotspot_detect_1kgp.py
# EXAMPLES: cat  sv-results.vcf    | python3 hotspot_detect_1kgp.py args
#           zcat sv-results.vcf.gz | python3 hotspot_detect_1kgp.py args


# important imports
import sys
import numpy as np


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


# vcf line class for multi individual/population file
class VCFLineSV1KGP(object):
    # init
    def __init__(self, input_line):
        # save the original
        self.svline = input_line
        # vcf files have 9 fields + the sample(s)/genome(s) data
        self.VCF_MANDATORY_FIELDS = 9
        tab_sep_fields = input_line.split("\t")
        self.ERROR = len(tab_sep_fields) < self.VCF_MANDATORY_FIELDS   # 9 VCF fields + genotypes
        if not self.ERROR:
            # extract mandatory fields
            # not used (_): REF, QUAL, FILTER
            [self.CHROM, POS, self.ID, _, _, _, _, self.INFO, self.FORMAT] = tab_sep_fields[:self.VCF_MANDATORY_FIELDS]
            self.SAMPLES = tab_sep_fields[self.VCF_MANDATORY_FIELDS:]
            self.POS = int(POS)
            self.SVLEN = 0
            self.END = 0
            self.hasSV = False
            self.maf = 0
            self.gt = []
            self.ac_af = []
            self.parse_gt()
            self.parse_info()
            rare_sv_maf_1percent = 0.01
            rare_sv_maf_5percent = 0.05
            rare_sv_maf_threshold = rare_sv_maf_1percent
            self.has_rare = self.maf < rare_sv_maf_threshold

    def parse_gt(self):
        for each_samp in self.SAMPLES:
            for fmt, gt in zip(self.FORMAT.split(":"), each_samp.split(":")):
                if fmt == "GT":
                    self.gt.append(gt)

    def parse_info(self):
        # INFO field extraction
        # ##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variation">
        # OR
        # ##INFO=<ID=END,Number=1,Type=Integer,Description="End position of structural variation">
        extract_info = ["SVLEN", "END", "AF", "AC"]
        if ";" in self.INFO:
            for each_info in self.INFO.split(";"):
                # DP2 unpack too many vals, SKIP IT
                [info_key, info_val] = each_info.split("=") if "=" in each_info else [each_info, ""]
                if info_key in extract_info:
                    self.SVLEN = int(info_val) if info_key == "SVLEN" else self.SVLEN
                    self.END = int(info_val) if info_key == "END" else self.END
                    if info_key == "AF":
                        self.ac_af.append(f'AF:{info_val}')
                        if "," in info_val:
                            af_list = [float(x) for x in info_val.split(",")]
                            self.maf = min(af_list)
                        else:
                            self.maf = float(info_val)
                    if info_key == "AC":
                        self.ac_af.append(f'AC:{info_val}')
        if self.SVLEN != 0:
            self.END = self.POS + self.SVLEN
        else:
            self.SVLEN = self.END - self.POS


def vcf_1kgp_rare_sv_detect():

    def get_genome_overlap(sv_info, window_size, the_chr):
        [start, end, has_raresv, ac_af, sv_id] = sv_info
        overlap_start = int(start/window_size)*window_size + 1
        overlap_ends = max(overlap_start, int(end/window_size)*window_size)
        if overlap_start == overlap_ends or overlap_start+window_size > overlap_ends:
            my_logger.dev(f'{the_chr}:{overlap_start}\t{the_chr}:{start}\t{has_raresv}\t{ac_af}\t{sv_id}')
            return [overlap_start], has_raresv
        else:
            for overlap in range(overlap_start, overlap_ends, window_size):
                my_logger.dev(f'{the_chr}:{overlap}\t{the_chr}:{start}\t{has_raresv}\t{ac_af}\t{sv_id}')
            return range(overlap_start, overlap_ends, window_size), has_raresv

    # loop file by line
    chr_lengths = {}
    current_chr = ""
    current_pos = 0
    sv_chromosomes = {}
    chr_rare_svs = {}
    vcf_entry = None
    # parse data
    for line in sys.stdin:
        if not line.startswith("#"):
            vcf_entry = VCFLineSV1KGP(line.rstrip("\n"))
            if vcf_entry.CHROM not in sv_chromosomes:
                sv_chromosomes[vcf_entry.CHROM] = [[vcf_entry.POS, vcf_entry.END, vcf_entry.has_rare,
                                                    vcf_entry.ac_af, vcf_entry.ID]]
            else:
                sv_chromosomes[vcf_entry.CHROM].append([vcf_entry.POS, vcf_entry.END, vcf_entry.has_rare,
                                                        vcf_entry.ac_af, vcf_entry.ID])
            if current_chr == "":
                current_chr = vcf_entry.CHROM
            elif current_chr != vcf_entry.CHROM:
                chr_lengths[current_chr] = current_pos
                current_chr = vcf_entry.CHROM
            else:
                pass
            current_pos = vcf_entry.POS
            if current_chr not in chr_rare_svs:
                chr_rare_svs[current_chr] = {}
    chr_lengths[current_chr] = current_pos
    window_size_analysis = 100000  # 100kb
    # init hotspot dict
    for each_chr in chr_lengths:
        for win_pos in range(1, chr_lengths[each_chr], window_size_analysis):
            chr_rare_svs[each_chr][win_pos] = False
    # check if has rare(s)
    for each_chr in chr_lengths:
        # loop 1: count SVs per window
        for each_sv in sv_chromosomes[each_chr]:
            sv_genome_window_overlap, has_rare_sv = get_genome_overlap(each_sv, window_size_analysis, each_chr)
            for win_pos in sv_genome_window_overlap:
                # if its rare, cannot be overwritten with not
                if not chr_rare_svs[each_chr][win_pos]:
                    chr_rare_svs[each_chr][win_pos] = has_rare_sv
    # # loop 3: call rare svs
    for each_chr in chr_lengths:
        for win_pos in chr_rare_svs[each_chr]:
            if chr_rare_svs[each_chr][win_pos]:
                print(f'{each_chr}:{win_pos}\thas_rare')
            else:
                print(f'{each_chr}:{win_pos}\tpass')


# main
if __name__ == '__main__':
    as_dev = True
    my_logger = SimpleLogger(as_dev)
    vcf_1kgp_rare_sv_detect()
