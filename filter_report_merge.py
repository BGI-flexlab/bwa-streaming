#!/usr/bin/env python
# encoding: utf-8
import os
import sys
import glob
import copy
from argparse import RawDescriptionHelpFormatter, ArgumentParser

__version__ = '1.0.0'
__date__ = '2017-06-28'
__updated__ = '2018-03-28'


def add(x, y): return x + y


def write_report(filepath, li):
    with open(filepath, 'w') as f:
        for l in li:
            f.write(l)
            f.write("\n")


class FastqInfo(object):
    def __init__(self):
        self.max_raw_read_len = 0
        self.max_clean_read_len = 0

        self.raw_base_num = 0
        self.clean_base_num = 0
        self.raw_base_A = 0
        self.clean_base_A = 0
        self.raw_base_C = 0
        self.clean_base_C = 0
        self.raw_base_G = 0
        self.clean_base_G = 0
        self.raw_base_T = 0
        self.clean_base_T = 0
        self.raw_base_N = 0
        self.clean_base_N = 0
        self.raw_q20 = 0
        self.clean_q20 = 0
        self.raw_q30 = 0
        self.clean_q30 = 0

        self.adapter_num = 0
        self.n_exceed_num = 0
        self.low_qual_num = 0
        self.low_mean_num = 0
        self.small_insert_num = 0
        self.polyA_num = 0

        self.max_quality_value = 0

        self.clean_read_len_distribution = []
        self.base = []
        self.clean_base = []
        self.q20q30 = []
        self.clean_q20q30 = []
        self.qual = []
        self.clean_qual = []

    def add_fqinfo(self, fqinfo):
        self.max_raw_read_len = fqinfo[0] if fqinfo[0] > self.max_raw_read_len else self.max_raw_read_len
        self.max_clean_read_len = fqinfo[1] if fqinfo[1] > self.max_clean_read_len else self.max_clean_read_len
        self.max_quality_value = fqinfo[24] if fqinfo[24] > self.max_quality_value else self.max_quality_value

        self.raw_base_num += fqinfo[2]
        self.clean_base_num += fqinfo[3]
        self.raw_base_A += fqinfo[4]
        self.clean_base_A += fqinfo[5]
        self.raw_base_C += fqinfo[6]
        self.clean_base_C += fqinfo[7]
        self.raw_base_G += fqinfo[8]
        self.clean_base_G += fqinfo[9]
        self.raw_base_T += fqinfo[10]
        self.clean_base_T += fqinfo[11]
        self.raw_base_N += fqinfo[12]
        self.clean_base_N += fqinfo[13]
        self.raw_q20 += fqinfo[14]
        self.clean_q20 += fqinfo[15]
        self.raw_q30 += fqinfo[16]
        self.clean_q30 += fqinfo[17]

        self.adapter_num += fqinfo[18]
        self.n_exceed_num += fqinfo[19]
        self.low_qual_num += fqinfo[20]
        self.low_mean_num += fqinfo[21]
        self.small_insert_num += fqinfo[22]
        self.polyA_num += fqinfo[23]

    def add_base_dis_info(self, base):
        if self.max_raw_read_len > len(self.base):
            extend_len = self.max_raw_read_len - len(self.base)
            i = 0
            while i < extend_len:
                self.base.append([0, 0, 0, 0, 0])
                self.clean_base.append([0, 0, 0, 0, 0])
                i += 1

        for pos, position_base in enumerate(base):
            i = 0
            while i < 5:
                self.base[pos][i] += position_base[i]
                i += 1

            i = 0
            while i < 5:
                self.clean_base[pos][i] += position_base[i + 5]
                i += 1

    def add_raw_qual(self, raw_qual):
        if self.max_raw_read_len > len(self.q20q30):
            extend_len = self.max_raw_read_len - len(self.q20q30)
            i = 0
            while i < extend_len:
                self.q20q30.append([0, 0])
                self.qual.append([0] * (len(raw_qual) - 2))
                # self.qual.append([0]*self.max_quality_value)
                i += 1

        for pos, qual in enumerate(raw_qual):
            self.q20q30[pos][0] += qual[0]
            self.q20q30[pos][1] += qual[1]

            i = 2
            while i < len(qual):
                self.qual[pos][i - 2] += qual[i]
                i += 1

    def add_clean_qual(self, clean_qual):
        if self.max_raw_read_len > len(self.clean_q20q30):
            extend_len = self.max_raw_read_len - len(self.clean_q20q30)
            i = 0
            while i < extend_len:
                self.clean_q20q30.append([0, 0])
                self.clean_qual.append([0] * (len(clean_qual) - 2))
                # self.qual.append([0]*self.max_quality_value)
                i += 1

        for pos, qual in enumerate(clean_qual):
            self.clean_q20q30[pos][0] += qual[0]
            self.clean_q20q30[pos][1] += qual[1]
            i = 2
            while i < len(qual):
                self.clean_qual[pos][i - 2] += qual[i]
                i += 1


class LaneReport(object):
    def __init__(self):
        self.total_raw_read_num = 0
        self.total_clean_read_num = 0
        self.total_short_length_num = 0
        self.total_n_exceed_num = 0
        self.total_low_qual_num = 0
        self.total_low_mean_num = 0
        self.total_adapter_num = 0
        self.total_cut_adapter_num = 0
        self.total_small_insert_num = 0
        self.total_polyA_num = 0
        self.read1_info = FastqInfo()
        self.read2_info = None

    def add(self, read_info):
        total_info = map(int, read_info[1].split())
        self.total_raw_read_num += total_info[0]
        self.total_clean_read_num += total_info[1]
        self.total_short_length_num += total_info[2]
        self.total_n_exceed_num += total_info[3]
        self.total_low_qual_num += total_info[4]
        self.total_low_mean_num += total_info[5]
        self.total_adapter_num += total_info[6]
        self.total_cut_adapter_num += total_info[7]
        self.total_small_insert_num += total_info[8]
        self.total_polyA_num += total_info[9]

        fq1_base = []
        fq1_raw_qual = []
        fq1_clean_qual = []

        fq2_base = []
        fq2_raw_qual = []
        fq2_clean_qual = []

        # Fq1_statistical_information
        fq1_info = map(int, read_info[3].split())
        self.read1_info.add_fqinfo(fq1_info)

        i = 5
        # Fq1 Base_distributions_by_read_position
        while not read_info[i].startswith('#'):
            fq1_base.append(map(int, read_info[i].split()))
            i += 1
        self.read1_info.add_base_dis_info(fq1_base)

        # Fq1 Raw_Base_quality_value_distribution_by_read_position
        i += 1
        while not read_info[i].startswith('#'):
            fq1_raw_qual.append(map(int, read_info[i].split()))
            i += 1
        self.read1_info.add_raw_qual(fq1_raw_qual)

        # Fq1 Clean_Base_quality_value_distribution_by_read_position
        i += 1
        while i < len(read_info) and not read_info[i].startswith('#'):
            fq1_clean_qual.append(map(int, read_info[i].split()))
            i += 1
        self.read1_info.add_clean_qual(fq1_clean_qual)

        if i >= len(read_info):
            return

        if self.read2_info is None:
            self.read2_info = FastqInfo()

        # Fq2_statistical_information
        i += 1
        self.read2_info.add_fqinfo(map(int, read_info[i].split()))

        # Fq2 Base_distributions_by_read_position
        i += 2
        while not read_info[i].startswith('#'):
            fq2_base.append(map(int, read_info[i].split()))
            i += 1
        self.read2_info.add_base_dis_info(fq2_base)

        # Fq2 Raw_Base_quality_value_distribution_by_read_position
        i += 1
        while not read_info[i].startswith('#'):
            fq2_raw_qual.append(map(int, read_info[i].split()))
            i += 1
        self.read2_info.add_raw_qual(fq2_raw_qual)

        # Fq2 Clean_Base_quality_value_distribution_by_read_position
        i += 1
        while i < len(read_info) and not read_info[i].startswith('#'):
            fq2_clean_qual.append(map(int, read_info[i].split()))
            i += 1
        self.read2_info.add_clean_qual(fq2_clean_qual)

    def print_basic_stat(self, report_prefix):
        # setw = len(str(self.total_raw_read_num)) + 10
        line = [''] * 14
        total_filter_read_num = self.total_raw_read_num - self.total_clean_read_num
        fq1_clean_base = self.read1_info.raw_base_num - self.read1_info.clean_base_num

        line[0] = "{:<65}\t{:<20}\t{:<20}".format('Item', 'raw reads(fq1)', 'clean reads(fq1)')
        line[1] = "{:<65}\t{:<20}\t{:<20}".format('Read length (max)',
                                                  self.read1_info.max_raw_read_len, self.read1_info.max_clean_read_len)
        line[2] = "{:<65}\t".format('Total number of reads')
        line[3] = "{:<65}\t".format('Number of filtered reads (%)')
        line[4] = "{:<65}\t".format('\nTotal number of bases')
        s = '{} ({:.2f}%)'.format(fq1_clean_base, 100.0 * fq1_clean_base / self.read1_info.raw_base_num)
        line[5] = "{:<65}\t{:<20}\t{:<20}\t".format('Number of filtered bases (%)', s, '-')
        line[6] = "{:<65}\t{:<20}\t{:<20}\t".format('\nReads related to Adapter and Trimmed (%)', '-', '-')
        raw_s = '{} ({:.2f}%)'.format(
            self.read1_info.raw_base_A, 100.0 * self.read1_info.raw_base_A / self.read1_info.raw_base_num)
        clean_s = '{} ({:.2f}%)'.format(
            self.read1_info.clean_base_A, 100.0 * self.read1_info.clean_base_A / self.read1_info.clean_base_num)
        line[7] = "{:<65}\t{:<20}\t{:<20}\t".format('\nNumber of base A (%)', raw_s, clean_s)
        raw_s = '{} ({:.2f}%)'.format(
            self.read1_info.raw_base_C, 100.0 * self.read1_info.raw_base_C / self.read1_info.raw_base_num)
        clean_s = '{} ({:.2f}%)'.format(
            self.read1_info.clean_base_C, 100.0 * self.read1_info.clean_base_C / self.read1_info.clean_base_num)
        line[8] = "{:<65}\t{:<20}\t{:<20}\t".format('Number of base C (%)', raw_s, clean_s)
        raw_s = '{} ({:.2f}%)'.format(
            self.read1_info.raw_base_G, 100.0 * self.read1_info.raw_base_G / self.read1_info.raw_base_num)
        clean_s = '{} ({:.2f}%)'.format(
            self.read1_info.clean_base_G, 100.0 * self.read1_info.clean_base_G / self.read1_info.clean_base_num)
        line[9] = "{:<65}\t{:<20}\t{:<20}\t".format('Number of base G (%)', raw_s, clean_s)
        raw_s = '{} ({:.2f}%)'.format(
            self.read1_info.raw_base_T, 100.0 * self.read1_info.raw_base_T / self.read1_info.raw_base_num)
        clean_s = '{} ({:.2f}%)'.format(
            self.read1_info.clean_base_T, 100.0 * self.read1_info.clean_base_T / self.read1_info.clean_base_num)
        line[10] = "{:<65}\t{:<20}\t{:<20}\t".format('Number of base T (%)', raw_s, clean_s)
        raw_s = '{} ({:.2f}%)'.format(
            self.read1_info.raw_base_N, 100.0 * self.read1_info.raw_base_N / self.read1_info.raw_base_num)
        clean_s = '{} ({:.2f}%)'.format(
            self.read1_info.clean_base_N, 100.0 * self.read1_info.clean_base_N / self.read1_info.clean_base_num)
        line[11] = "{:<65}\t{:<20}\t{:<20}\t".format('Number of base N (%)', raw_s, clean_s)
        raw_s = '{} ({:.2f}%)'.format(
            self.read1_info.raw_q20, 100.0 * self.read1_info.raw_q20 / self.read1_info.raw_base_num)
        clean_s = '{} ({:.2f}%)'.format(
            self.read1_info.clean_q20, 100.0 * self.read1_info.clean_q20 / self.read1_info.clean_base_num)
        line[12] = "{:<65}\t{:<20}\t{:<20}\t".format(
            '\nNumber of base calls with quality value of 20 or higher (Q20+) (%)', raw_s, clean_s)
        raw_s = '{} ({:.2f}%)'.format(
            self.read1_info.raw_q30, 100.0 * self.read1_info.raw_q30 / self.read1_info.raw_base_num)
        clean_s = '{} ({:.2f}%)'.format(
            self.read1_info.clean_q30, 100.0 * self.read1_info.clean_q30 / self.read1_info.clean_base_num)
        line[13] = "{:<65}\t{:<20}\t{:<20}\t".format(
            'Number of base calls with quality value of 30 or higher (Q30+) (%)', raw_s, clean_s)

        if self.read2_info is not None:
            line[0] += "\t{:<20}\t{:<20}".format('raw reads(fq2)', 'clean reads(fq2)')
            line[1] += "\t{:<20}\t{:<20}".format(self.read2_info.max_raw_read_len, self.read2_info.max_raw_read_len)
            line[2] += "{0:<20}\t{1:<20}\t{0:<20}\t{1:<20}".format(self.total_raw_read_num,
                                                                   self.total_clean_read_num)
            s = '{} ({:.2f}%)'.format(total_filter_read_num, 100.0 * total_filter_read_num / self.total_raw_read_num)
            line[3] += "{0:<20}\t{1:<20}\t{0:<20}\t{1:<20}".format(s, '-')
            line[4] += "{0:<20}\t{1:<20}\t{2:<20}\t{3:<20}".format(
                self.read1_info.raw_base_num, self.read1_info.clean_base_num,
                self.read2_info.raw_base_num, self.read2_info.clean_base_num)
            fq2_clean_base = self.read2_info.raw_base_num - self.read2_info.clean_base_num
            s = '{} ({:.2f}%)'.format(fq2_clean_base, 100.0 * fq2_clean_base / self.read2_info.raw_base_num)
            line[5] += "{0:<20}\t{1:<20}".format(s, '-')
            line[6] += "{0:<20}\t{0:<20}".format('-')
            raw_s = '{} ({:.2f}%)'.format(
                self.read2_info.raw_base_A, 100.0 * self.read2_info.raw_base_A / self.read2_info.raw_base_num)
            clean_s = '{} ({:.2f}%)'.format(
                self.read2_info.clean_base_A, 100.0 * self.read2_info.clean_base_A / self.read2_info.clean_base_num)
            line[7] += "{:<20}\t{:<20}".format(raw_s, clean_s)
            raw_s = '{} ({:.2f}%)'.format(
                self.read2_info.raw_base_C, 100.0 * self.read2_info.raw_base_C / self.read2_info.raw_base_num)
            clean_s = '{} ({:.2f}%)'.format(
                self.read2_info.clean_base_C, 100.0 * self.read2_info.clean_base_C / self.read2_info.clean_base_num)
            line[8] += "{:<20}\t{:<20}".format(raw_s, clean_s)
            raw_s = '{} ({:.2f}%)'.format(
                self.read2_info.raw_base_G, 100.0 * self.read2_info.raw_base_G / self.read2_info.raw_base_num)
            clean_s = '{} ({:.2f}%)'.format(
                self.read2_info.clean_base_G, 100.0 * self.read2_info.clean_base_G / self.read2_info.clean_base_num)
            line[9] += "{:<20}\t{:<20}".format(raw_s, clean_s)
            raw_s = '{} ({:.2f}%)'.format(
                self.read2_info.raw_base_T, 100.0 * self.read2_info.raw_base_T / self.read2_info.raw_base_num)
            clean_s = '{} ({:.2f}%)'.format(
                self.read2_info.clean_base_T, 100.0 * self.read2_info.clean_base_T / self.read2_info.clean_base_num)
            line[10] += "{:<20}\t{:<20}".format(raw_s, clean_s)
            raw_s = '{} ({:.2f}%)'.format(
                self.read2_info.raw_base_N, 100.0 * self.read2_info.raw_base_N / self.read2_info.raw_base_num)
            clean_s = '{} ({:.2f}%)'.format(
                self.read2_info.clean_base_N, 100.0 * self.read2_info.clean_base_N / self.read2_info.clean_base_num)
            line[11] += "{:<20}\t{:<20}".format(raw_s, clean_s)
            raw_s = '{} ({:.2f}%)'.format(
                self.read2_info.raw_q20, 100.0 * self.read2_info.raw_q20 / self.read1_info.raw_base_num)
            clean_s = '{} ({:.2f}%)'.format(
                self.read2_info.clean_q20, 100.0 * self.read2_info.clean_q20 / self.read1_info.clean_base_num)
            line[12] += "{:<20}\t{:<20}".format(raw_s, clean_s)
            raw_s = '{} ({:.2f}%)'.format(
                self.read2_info.raw_q30, 100.0 * self.read2_info.raw_q30 / self.read2_info.raw_base_num)
            clean_s = '{} ({:.2f}%)'.format(
                self.read2_info.clean_q30, 100.0 * self.read2_info.clean_q30 / self.read2_info.clean_base_num)
            line[13] += "{:<20}\t{:<20}".format(raw_s, clean_s)


        else:
            line[2] = "{0:<20}\t{1:<20}\t{2:<20}".format(
                'Total number of reads', self.total_raw_read_num, self.total_clean_read_num)
            line[3] += "{0:<20}\t{1:<20}".format(
                100.0 * total_filter_read_num / self.total_raw_read_num, '-')
            line[4] += "{0:<20}\t{1:<20}".format(self.read1_info.raw_base_num, self.read1_info.clean_base_num)

        write_report(report_prefix+'.txt', line)

    def print_filter_stat(self, report_prefix):
        if self.read2_info is not None:
            lines = ["{:<65}\t{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format(
                'Item', 'Total', 'Percentage', 'Counts(fq1)', 'Percentage', 'Counts(fq1)', 'Percentage')]
            filter_num = self.total_raw_read_num - self.total_clean_read_num
            fq1_filter_num = (self.total_raw_read_num - self.total_clean_read_num) / 2
            if filter_num == 0:
                filter_num = 1
                fq1_filter_num = 1
                lines.append("{:<65}\t{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format(
                    'Total filtered reads', 0, '{:.2%}'.format(0), 0, '{:.2%}'.format(0), 0, '{:.2%}'.format(0)))
            else:
                lines.append("{:<65}\t{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format(
                    'Total filtered reads', filter_num, '{:.2%}'.format(1), fq1_filter_num,
                    '{:.2%}'.format(1), fq1_filter_num, '{:.2%}'.format(1), ))
            p0 = '{:.2%}'.format(self.total_adapter_num * 1.0 / filter_num)
            p1 = '{:.2%}'.format(self.read1_info.adapter_num * 1.0 / fq1_filter_num)
            p2 = '{:.2%}'.format(self.read2_info.adapter_num * 1.0 / fq1_filter_num)
            lines.append("{:<65}\t{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format(
                'Reads with adapter', self.total_adapter_num, p0, self.read1_info.adapter_num, p1,
                self.read2_info.adapter_num, p2))
            p0 = '{:.2%}'.format(self.total_low_qual_num * 1.0 / filter_num)
            p1 = '{:.2%}'.format(self.read1_info.low_qual_num * 1.0 / fq1_filter_num)
            p2 = '{:.2%}'.format(self.read2_info.low_qual_num * 1.0 / fq1_filter_num)
            lines.append("{:<65}\t{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format(
                'Reads with low quality', self.total_low_qual_num, p0, self.read1_info.low_qual_num, p1,
                self.read2_info.low_qual_num, p2))
            p0 = '{:.2%}'.format(self.total_low_mean_num * 1.0 / filter_num)
            p1 = '{:.2%}'.format(self.read1_info.low_mean_num * 1.0 / fq1_filter_num)
            p2 = '{:.2%}'.format(self.read2_info.low_mean_num * 1.0 / fq1_filter_num)
            lines.append("{:<65}\t{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format(
                'Reads with low mean quality', self.total_low_mean_num, p0, self.read1_info.low_mean_num, p1,
                self.read2_info.low_mean_num, p2))
            # lines.append("{:<65}\t{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format(
            #     'Reads with duplications', 0, 0, 0, 0, 0, 0))
            p0 = '{:.2%}'.format(self.total_n_exceed_num * 1.0 / filter_num)
            p1 = '{:.2%}'.format(self.read1_info.n_exceed_num * 1.0 / fq1_filter_num)
            p2 = '{:.2%}'.format(self.read2_info.n_exceed_num * 1.0 / fq1_filter_num)
            lines.append("{:<65}\t{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format(
                'Read with n rate exceed', self.total_n_exceed_num, p0, self.read1_info.n_exceed_num, p1,
                self.read2_info.n_exceed_num, p2))
            p0 = '{:.2%}'.format(self.total_small_insert_num * 1.0 / filter_num)
            p1 = '{:.2%}'.format(self.read1_info.small_insert_num * 1.0 / fq1_filter_num)
            p2 = '{:.2%}'.format(self.read2_info.small_insert_num * 1.0 / fq1_filter_num)
            lines.append("{:<65}\t{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format(
                'Read with small insert size', self.total_small_insert_num, p0, self.read1_info.small_insert_num, p1,
                self.read2_info.small_insert_num, p2))
            p0 = '{:.2%}'.format(self.total_polyA_num * 1.0 / filter_num)
            p1 = '{:.2%}'.format(self.read1_info.polyA_num * 1.0 / fq1_filter_num)
            p2 = '{:.2%}'.format(self.read2_info.polyA_num * 1.0 / fq1_filter_num)
            lines.append("{:<65}\t{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format(
                'Reads with PolyA', self.total_polyA_num, p0, self.read1_info.polyA_num, p1,
                self.read2_info.polyA_num, p2))
        else:
            lines = ["{:<65}\t{:<20}\t{:<20}".format(
                'Item', 'Counts(fq1)', 'Percentage')]
            filter_num = self.total_raw_read_num - self.total_clean_read_num
            if filter_num == 0:
                filter_num = 1
                lines.append("{:<65}\t{:<20}\t{:<20}".format(
                    'Total filtered reads', 0, '{:.2%}'.format(0)))
            else:
                lines.append("{:<65}\t{:<20}\t{:<20}".format(
                    'Total filtered reads', filter_num, '{:.2%}'.format(1)))
            p0 = '{:.2%}'.format(self.total_adapter_num * 1.0 / filter_num)
            lines.append("{:<65}\t{:<20}\t{:<20}".format(
                'Reads with adapter', self.total_adapter_num, p0))
            p0 = '{:.2%}'.format(self.total_low_qual_num * 1.0 / filter_num)
            lines.append("{:<65}\t{:<20}\t{:<20}".format(
                'Reads with low quality', self.total_low_qual_num, p0))
            p0 = '{:.2%}'.format(self.total_low_mean_num * 1.0 / filter_num)
            lines.append("{:<65}\t{:<20}\t{:<20}".format(
                'Reads with low mean quality', self.total_low_mean_num, p0))
            # lines.append("{:<65}\t{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format(
            #     'Reads with duplications', 0, 0, 0, 0, 0, 0))
            p0 = '{:.2%}'.format(self.total_n_exceed_num * 1.0 / filter_num)
            lines.append("{:<65}\t{:<20}\t{:<20}".format(
                'Read with n rate exceed', self.total_n_exceed_num, p0))
            p0 = '{:.2%}'.format(self.total_small_insert_num * 1.0 / filter_num)
            lines.append("{:<65}\t{:<20}\t{:<20}".format(
                'Read with small insert size', self.total_small_insert_num, p0))
            p0 = '{:.2%}'.format(self.total_polyA_num * 1.0 / filter_num)
            lines.append("{:<65}\t{:<20}\t{:<20}".format('Reads with PolyA', self.total_polyA_num, p0))

        write_report(report_prefix + '.txt', lines)

    def print_base_dist(self, report_prefix):
        header = "{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format(
            'Pos', 'A', 'C', 'G', 'T', 'N', 'Clean A', 'Clean C', 'Clean G', 'Clean T', 'Clean N')
        i = 0
        r = self.read1_info
        lines = [header]
        while i < r.max_raw_read_len:
            n = (r.base[i][0] + r.base[i][1] + r.base[i][2] + r.base[i][3] + r.base[i][4])*1.0
            l = '{:<20}'.format(i+1)
            for num in r.base[i]:
                p = '{:.2%}'.format(num/n)
                l += '\t{:<20}'.format(p)

            for num in r.clean_base[i]:
                p = '{:.2%}'.format(num/n)
                l += '\t{:<20}'.format(p)
            i += 1
            lines.append(l)

        write_report(report_prefix + '_1.txt', lines)

        if self.read2_info is not None:
            i = 0
            r = self.read2_info
            lines = [header]
            while i < r.max_raw_read_len:
                n = (r.base[i][0] + r.base[i][1] + r.base[i][2] + r.base[i][3] + r.base[i][4]) * 1.0
                l = '{:<20}'.format(i + 1)
                for num in r.base[i]:
                    p = '{:.2%}'.format(num / n)
                    l += '\t{:<20}'.format(p)

                for num in r.clean_base[i]:
                    p = '{:.2%}'.format(num / n)
                    l += '\t{:<20}'.format(p)
                i += 1
                lines.append(l)

            write_report(report_prefix + '_2.txt', lines)

    def print_base_qual(self, report_prefix):
        i = 0
        h = ["{:<20}".format('Pos')]
        while i < self.read1_info.max_quality_value:
            h.append("Q{:<20}".format(i))
            i += 1

        # TODO add this stat
        # h.append('Mean')
        # h.append('Median')
        # h.append('Lower quartile')
        # h.append('Upper quartile')
        # h.append('10thpercentile')
        # h.append('90thpercentile')

        header = "\t".join(h)
        i = 0
        r = self.read1_info
        lines = [header]
        while i < r.max_raw_read_len:
            l = '{:<20}'.format(i + 1)
            j = 0
            while j < r.max_quality_value:
                l += '\t{:<20}'.format(r.qual[i][j])
                j += 1
            i += 1
            lines.append(l)

        lines.append("Clean Quality Value Distribute")
        lines.append(header)
        i = 0
        while i < r.max_raw_read_len:
            l = '{:<20}'.format(i + 1)
            j = 0
            while j < r.max_quality_value:
                l += '\t{:<20}'.format(r.clean_qual[i][j])
                j += 1
            i += 1
            lines.append(l)

        write_report(report_prefix + '_1.txt', lines)

        if self.read2_info is not None:
            header = "\t".join(h)
            i = 0
            r = self.read2_info
            lines = [header]
            while i < r.max_raw_read_len:
                l = '{:<20}'.format(i + 1)
                j = 0
                while j < r.max_quality_value:
                    l += '\t{:<20}'.format(r.qual[i][j])
                    j += 1

                i += 1
                lines.append(l)

            lines.append("Clean Quality Value Distribute")
            lines.append(header)
            i = 0
            while i < r.max_raw_read_len:
                l = '{:<20}'.format(i + 1)
                j = 0
                while j < r.max_quality_value:
                    l += '\t{:<20}'.format(r.clean_qual[i][j])
                    j += 1
                i += 1
                lines.append(l)

            write_report(report_prefix + '_2.txt', lines)

    def print_q20q30(self, report_prefix):
        h = ["Position in reads", "Percentage of Q20+ bases",
             "Percentage of Q30+ bases", "Percentage of Clean Q20+", "Percentage of Clean Q30+"]
        header = "\t".join(h)
        i = 0
        r = self.read1_info
        lines = [header]
        while i < r.max_raw_read_len:
            l = '{:<20}'.format(i + 1)
            pos_base_num = sum(r.qual[i])
            p20 = '{:.2%}'.format(float(r.q20q30[i][0])/ pos_base_num)
            p30 = '{:.2%}'.format(float(r.q20q30[i][1]) / pos_base_num)
            l += '{:<20}'.format(p20)
            l += '{:<20}'.format(p30)
            clean_pos_base_num = sum(r.clean_qual[i])
            p20 = '{:.2%}'.format(float(r.clean_q20q30[i][0]) / clean_pos_base_num)
            p30 = '{:.2%}'.format(float(r.clean_q20q30[i][1]) / clean_pos_base_num)
            l += '{:<20}'.format(p20)
            l += '{:<20}'.format(p30)
            i += 1
            lines.append(l)

        write_report(report_prefix + '_1.txt', lines)

        if self.read2_info is not None:
            i = 0
            r = self.read2_info
            lines = [header]
            while i < r.max_raw_read_len:
                l = '{:<20}'.format(i + 1)
                pos_base_num = sum(r.qual[i])
                p20 = '{:.2%}'.format(float(r.q20q30[i][0]) / pos_base_num)
                p30 = '{:.2%}'.format(float(r.q20q30[i][1]) / pos_base_num)
                l += '{:<20}'.format(p20)
                l += '{:<20}'.format(p30)
                clean_pos_base_num = sum(r.clean_qual[i])
                p20 = '{:.2%}'.format(float(r.clean_q20q30[i][0]) / clean_pos_base_num)
                p30 = '{:.2%}'.format(float(r.clean_q20q30[i][1]) / clean_pos_base_num)
                l += '{:<20}'.format(p20)
                l += '{:<20}'.format(p30)
                i += 1
                lines.append(l)

            write_report(report_prefix + '_2.txt', lines)


class Report(object):
    def __init__(self, outdir=None):
        self.reports = {}
        self.outdir = './'
        if outdir is not None:
            self.outdir = outdir

    def parse_report(self, report_file):
        lane_id, lane_name = (0, 'null')
        with open(report_file, 'r') as f:
            lines = (line.rstrip('#S\n') for line in f)
            for line in lines:
                if line.startswith('>'):
                    sample_line = line
                    lane_id, lane_name = sample_line.lstrip('>').rstrip().split()
                    if lane_id not in self.reports:
                        self.reports[lane_id] = {}
                        self.reports[lane_id]['lane_name'] = lane_name
                        self.reports[lane_id]['report'] = LaneReport()
                    self.reports[lane_id]['read_info'] = []
                else:
                    self.reports[lane_id]['read_info'].append(line.strip())

            for lane_id in self.reports:
                self.reports[lane_id]['report'].add(self.reports[lane_id]['read_info'])

    def add(self, report):
        self.parse_report(report)

    def print_quality_report(self):
        for lane_id in self.reports:
            prefix = self.reports[lane_id]['lane_name'] + '_'
            if len(self.reports) == 1:
                prefix = ''
            # quality_file = os.path.join(self.outdir, prefix + '.SEQUENCING_QUALITY.txt')
            f_basic_stat = os.path.join(self.outdir, prefix + 'Basic_Statistics_of_Sequencing_Quality')
            f_filter_stat = os.path.join(self.outdir, prefix + 'Statistics_of_Filtered_Reads')
            f_base_quality = os.path.join(self.outdir, prefix + 'Base_quality_value_distribution_by_read_position')
            f_base_dist = os.path.join(self.outdir, prefix + 'Base_distributions_by_read_position')
            f_q20q30 = os.path.join(self.outdir, prefix + 'Distribution_of_Q20_Q30_bases_by_read_position')
            report = self.reports[lane_id]['report']
            report.print_basic_stat(f_basic_stat)
            # report.print_quality_report()
            report.print_filter_stat(f_filter_stat)
            report.print_base_dist(f_base_dist)
            report.print_base_qual(f_base_quality)
            report.print_q20q30(f_q20q30)


def main():
    program_name = os.path.basename(sys.argv[0])
    program_license = '''{0}
      Created by huangzhibo on {1}.
      Last updated on {2}.
      Copyright 2017 BGI bigData. All rights reserved.
    USAGE'''.format(" v".join([program_name, __version__]), str(__date__), str(__updated__))

    parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", dest="input", help='input dir[requested] .')
    parser.add_argument("-o", "--outdir", dest="outdir", default='./', help='outdir [./] .')

    if len(sys.argv) == 1:
        parser.print_help()
        exit(0)
    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    report = Report(args.outdir)
    part_list = glob.glob(args.input + '/part*')
    for part in part_list:
        report.add(part)

    report.print_quality_report()


if __name__ == "__main__":
    sys.exit(main())
