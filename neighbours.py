#!/usr/bin/env python

__author__ = 'Maxim K'

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
        pos = int(pos)
        neibs_thread = []
        for curr in range(pos - 2, pos + 3):
            neibs_thread.append(correspond_aa(table, str(curr)))
        neibs.append(neibs_thread)

    return neibs