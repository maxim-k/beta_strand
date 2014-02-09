#!/usr/bin/env python

__author__ = 'Maxim K'

import os
from hist import if_add, write_dic

def correspond_aa(thread_list, num):
    for line in thread_list:
        if line[0] == num:
            return line[1]
    return ''

def generate_neib(table, pos_neib):
    '''
    Для каждой записи генерируется список из трёх пятёрок:
    по два соседа на стренде вверх\вниз,
    такие же пять аминокислот на левом соседнем стренде,
    такие же пять аминокислот на правом соседнем стренде.
    '''
    neibs = []

    for pos in pos_neib:
        if(pos != '0'):
            pos = int(pos)
            neibs_thread = []
            for curr in range(pos - 2, pos + 3):
                neibs_thread.append(correspond_aa(table, str(curr)))
            neibs.append(neibs_thread)
        else:
            neibs.append(['', '', '', '', ''])

    return neibs



path = 'E:\\Science\\MG\\Marat\\server\\pairwise\\neib'
dir_output = 'E:\\Science\\MG\Marat\\pairwise\\neib_pairwise\\'

for dirname, dirnames, filenames in os.walk(path):
    for filename in filenames:
        file = open(os.path.join(dirname, filename), 'r')

        L_2, L_1, L_0, L1, L2 = {}, {}, {}, {}, {}
        R_2, R_1, R_0, R1, R2 = {}, {}, {}, {}, {}
        C_2, C_1, C1, C2 = {}, {}, {}, {}

        for line in file:
            aa = line.split(sep = ';')
            if aa[0]:
                tmp = [aa[0], aa[7]]
                if_add(L_2, ' '.join(tmp))
            elif aa[1]:
                tmp = [aa[1], aa[7]]
                if_add(L_1, ' '.join(tmp))
            elif aa[2]:
                tmp = [aa[2], aa[7]]
                if_add(L_0, ' '.join(tmp))
            elif aa[3]:
                tmp = [aa[3], aa[7]]
                if_add(L1, ' '.join(tmp))
            elif aa[4]:
                tmp = [aa[4], aa[7]]
                if_add(L2, ' '.join(tmp))
            elif aa[5]:
                tmp = [aa[5], aa[7]]
                if_add(C_2, ' '.join(tmp))
            elif aa[6]:
                tmp = [aa[6], aa[7]]
                if_add(C_1, ' '.join(tmp))
            elif aa[8]:
                tmp = [aa[8], aa[7]]
                if_add(C1, ' '.join(tmp))
            elif aa[9]:
                tmp = [aa[9], aa[7]]
                if_add(C2, ' '.join(tmp))
            elif aa[10]:
                tmp = [aa[10], aa[7]]
                if_add(R_2, ' '.join(tmp))
            elif aa[11]:
                tmp = [aa[11], aa[7]]
                if_add(R_1, ' '.join(tmp))
            elif aa[12]:
                tmp = [aa[12], aa[7]]
                if_add(R_0, ' '.join(tmp))
            elif aa[13]:
                tmp = [aa[13], aa[7]]
                if_add(R1, ' '.join(tmp))
            elif aa[14]:
                tmp = [aa[14], aa[7]]
                if_add(R2, ' '.join(tmp))

        L_2_file = open(dir_output + filename[:-4] + '_-2L.txt', 'w')
        write_dic(L_2, L_2_file, 'neighbour base quantity')
        L_1_file = open(dir_output + filename[:-4] + '_-1L.txt', 'w')
        write_dic(L_1, L_1_file, 'neighbour base quantity')
        L_0_file = open(dir_output + filename[:-4] + '_0L.txt', 'w')
        write_dic(L_0, L_0_file, 'neighbour base quantity')
        L1_file = open(dir_output + filename[:-4] + '_1L.txt', 'w')
        write_dic(L1, L1_file, 'neighbour base quantity')
        L2_file = open(dir_output + filename[:-4] + '_2L.txt', 'w')
        write_dic(L2, L2_file, 'neighbour base quantity')

        R_2_file = open(dir_output + filename[:-4] + '_-2R.txt', 'w')
        write_dic(R_2, R_2_file, 'neighbour base quantity')
        R_1_file = open(dir_output + filename[:-4] + '_-1R.txt', 'w')
        write_dic(R_1, R_1_file, 'neighbour base quantity')
        R_0_file = open(dir_output + filename[:-4] + '_0R.txt', 'w')
        write_dic(R_0, R_0_file, 'neighbour base quantity')
        R1_file = open(dir_output + filename[:-4] + '_1R.txt', 'w')
        write_dic(R1, R1_file, 'neighbour base quantity')
        R2_file = open(dir_output + filename[:-4] + '_2R.txt', 'w')
        write_dic(R2, R2_file, 'neighbour base quantity')

        C_2_file = open(dir_output + filename[:-4] + '_-2C.txt', 'w')
        write_dic(C_2, C_2_file, 'neighbour base quantity')
        C_1_file = open(dir_output + filename[:-4] + '_-1C.txt', 'w')
        write_dic(C_1, C_1_file, 'neighbour base quantity')
        C1_file = open(dir_output + filename[:-4] + '_1C.txt', 'w')
        write_dic(C1, C1_file, 'neighbour base quantity')
        C2_file = open(dir_output + filename[:-4] + '_2C.txt', 'w')
        write_dic(C2, C2_file, 'neighbour base quantity')

