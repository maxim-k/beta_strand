#!/usr/bin/env python

__author__ = 'Maxim K'
import string
import hist

def correspond_aa(thread_list, num):
    for line in thread_list:
        if line[0] == num:
            return line[1]
    return None

def count_pairwise(aa_list, dir_output):
    parallel = {}
    antiparallel = {}
    double_parallel = {}
    double_antiparallel ={}
    summary = {}

    for aa_line in aa_list:
        for sheet in aa_line:
            thread_list = []
            for thread in aa_line[sheet]:
                for line in aa_line[sheet][thread]:
                    if(thread == 's1')or(thread == 's2'):
                        if(len(line[0]) == 1):
                            line[0] = line[0][0]
                        if line[0][2] in string.ascii_lowercase:
                            line[0][2] = 'C'
                        thread_list.append([line[0][0], line[0][2], line[1]])

            for thread in aa_line[sheet]:
                for line in aa_line[sheet][thread]:
                    if(thread == 's1')or(thread == 's2'):
                        corr_aa = correspond_aa(thread_list, line[0][4] if thread == 's1' else line[0][5])
                        if corr_aa:
                            key = ' '.join(sorted([line[0][2], corr_aa]))
                            hist.if_add(summary, key)
                            if(line[1]):
                                hist.if_add(parallel, key)
                            else:
                                hist.if_add(antiparallel, key)

    for aa_line in aa_list:
        for sheet in aa_line:
            thread_list = []
            for thread in aa_line[sheet]:
                for line in aa_line[sheet][thread]:
                    if(thread == 's4')and(aa_line[sheet][thread]):
                        if(len(line[0]) == 1):
                            line[0] = line[0][0]
                        if line[0][2] in string.ascii_lowercase:
                            line[0][2] = 'C'
                        thread_list.append([line[0][0], line[0][2], line[1]])

            for thread in aa_line[sheet]:
                for line in aa_line[sheet][thread]:
                    if(thread == 's4')and(aa_line[sheet][thread]):
                        corr_aa = correspond_aa(thread_list, line[0][4] if line[0][5] == '0' else line[0][5])
                        if corr_aa:
                            key = ' '.join(sorted([line[0][2], corr_aa]))
                            hist.if_add(summary, key)
                            if(line[1]):
                                hist.if_add(double_parallel, key)
                            else:
                                hist.if_add(double_antiparallel, key)

    parallel_file = open(dir_output + '\\parallel.txt', 'w')
    hist.write_dic(parallel, parallel_file, 'aa1 aa2 quantity')
    antiparallel_file = open(dir_output + '\\antiparallel.txt', 'w')
    hist.write_dic(antiparallel, antiparallel_file, 'aa1 aa2 quantity')
    double_parallel_file = open(dir_output + '\\double_parallel.txt', 'w')
    hist.write_dic(parallel, double_parallel_file, 'aa1 aa2 quantity')
    double_antiparallel_file = open(dir_output + '\\double_antiparallel.txt', 'w')
    hist.write_dic(antiparallel, double_antiparallel_file, 'aa1 aa2 quantity')
    summary_file = open(dir_output + '\\summary.txt', 'w')
    hist.write_dic(summary, summary_file, 'aa1 aa2 quantity')