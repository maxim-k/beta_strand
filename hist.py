#!/usr/bin/env python

__author__ = 'Maxim K'
import string

def if_add(dic, val):
    if val in dic:
        dic[val] +=1
    else:
        dic[val] = 1

def write_dic(dic, file, header):
    file.write('%s\n' %(header))
    for key in dic:
        file.write('%s %d\n' % (key, dic[key]))

def count_aa(aa_list, dir_output):
    alone = {}
    single = {}
    double = {}
    double_parallel = {}
    double_antiparallel = {}
    parallel = {}
    antiparallel = {}
    summary = {}

    for aa_line in aa_list:
        for sheet in aa_line:
            for thread in aa_line[sheet]:
                for line in aa_line[sheet][thread]:
                    if(thread == 's1')or(thread == 's2'):
                        if(len(line[0]) == 1):
                            line[0] = line[0][0]
                        if line[0][2] in string.ascii_lowercase:
                            line[0][2] = 'C'
                        if_add(summary, line[0][2])
                        if_add(single, line[0][2])
                        if(line[1]):
                            if_add(parallel, line[0][2])
                        else:
                            if_add(antiparallel, line[0][2])
                    elif(thread == 's3'):
                        if line[2] in string.ascii_lowercase:
                            line[2] = 'C'
                        if_add(double, line[2])
                        if_add(summary, line[2])
                    elif(thread == 's4'):
                        if(len(line[0]) == 1):
                            line[0] = line[0][0]
                        if line[0][2] in string.ascii_lowercase:
                            line[0][2] = 'C'
                        if(line[1]):
                            if_add(double_parallel, line[0][2])
                        else:
                            if_add(double_antiparallel, line[0][2])
                    else:
                        if line[2] in string.ascii_lowercase:
                            line[2] = 'C'
                        if_add(alone, line[2])
                        if_add(summary, line[2])

    alone_file = open(dir_output + '\\alone.txt', 'w')
    write_dic(alone, alone_file, 'aa quantity')
    single_file = open(dir_output + '\\single.txt', 'w')
    write_dic(single, single_file, 'aa quantity')
    double_file = open(dir_output + '\\double.txt', 'w')
    write_dic(double, double_file, 'aa quantity')
    double_p_file = open(dir_output + '\\double_parallel.txt', 'w')
    write_dic(double_parallel, double_p_file, 'aa quantity')
    double_a_file = open(dir_output + '\\double_antiparallel.txt', 'w')
    write_dic(double_antiparallel, double_a_file, 'aa quantity')
    parallel_file = open(dir_output + '\\parallel.txt', 'w')
    write_dic(parallel, parallel_file, 'aa quantity')
    antiparallel_file = open(dir_output + '\\antiparallel.txt', 'w')
    write_dic(antiparallel, antiparallel_file, 'aa quantity')
    summary_file = open(dir_output + '\\summary.txt', 'w')
    write_dic(summary, summary_file, 'aa quantity')

