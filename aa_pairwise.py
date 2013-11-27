#!/usr/bin/env python

__author__ = 'Maxim K'
import string
import hist
from neighbours import generate_neib

def correspond_aa(thread_list, num):
    for line in thread_list:
        if line[0] == num:
            return line[1]
    return ''

def count_pairwise(aa_list, dir_output):
    parallel = {}
    antiparallel = {}
    double_parallel = {}
    double_antiparallel ={}
    summary = {}

    neib_single_parallel = []
    neib_single_antiparallel = []
    neib_double_parallel = []
    neib_double_antiparallel = []


    # заменяем цистеины с дисульфидыми мостиками, которые обозначены разными  буквами, на 'C'
    # создаём список с элементами вида '[позиция, аминокислота, параллельность]'
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

            # Из каждой первоначальной записи берется аминокислота и номер соседа.
            # Ищется номер соседа и определяется его аминокислота.
            # В словаре создаётся\инкрементируется ключ вида 'аминокислота аминокислота_соседа'.
            for thread in aa_line[sheet]:
                for line in aa_line[sheet][thread]:
                    if(thread == 's1')or(thread == 's2'):
                        corr_aa = correspond_aa(thread_list, line[0][4] if thread == 's1' else line[0][5])

                        if line[1]:
                            neib_single_parallel.append(generate_neib(thread_list,[line[0][4], line[0][0], line[0][5]]))
                        else:
                            neib_single_antiparallel.append(generate_neib(thread_list,[line[0][4], line[0][0], line[0][5]]))

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

                        if line[1]:
                            neib_double_parallel.append(generate_neib(thread_list,[line[0][4], line[0][0], line[0][5]]))
                        else:
                            neib_double_antiparallel.append(generate_neib(thread_list,[line[0][4], line[0][0], line[0][5]]))

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
    hist.write_dic(double_parallel, double_parallel_file, 'aa1 aa2 quantity')
    double_antiparallel_file = open(dir_output + '\\double_antiparallel.txt', 'w')
    hist.write_dic(double_antiparallel, double_antiparallel_file, 'aa1 aa2 quantity')
    summary_file = open(dir_output + '\\summary.txt', 'w')
    hist.write_dic(summary, summary_file, 'aa1 aa2 quantity')

    nsp_file = open(dir_output + '\\neib\\neib_single_parallel.txt', 'w')
    for thread in neib_single_parallel:
        for line in thread:
            nsp_file.write(';'.join(line) + ';')
        nsp_file.write('\n')

    nsa_file = open(dir_output + '\\neib\\neib_single_antiparallel.txt', 'w')
    for thread in neib_single_antiparallel:
        for line in thread:
            nsa_file.write(';'.join(line) + ';')
        nsa_file.write('\n')

    ndp_file = open(dir_output + '\\neib\\neib_double_parallel.txt', 'w')
    for thread in neib_double_parallel:
        for line in thread:
            ndp_file.write(';'.join(line) + ';')
        ndp_file.write('\n')

    nda_file = open(dir_output + '\\neib\\neib_double_antiparallel.txt', 'w')
    for thread in neib_double_antiparallel:
        for line in thread:
            nda_file.write(';'.join(line) + ';')
        nda_file.write('\n')
