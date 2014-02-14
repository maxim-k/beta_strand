__author__ = 'Maxim K'

import os, math

path = 'E:\\Science\\MG\Marat\\pairwise\\neib_pairwise'


def get_strand_name(filename):
    '''
    Разделяет имя файла на части
    '''

    name = filename.replace('neib_','').split(sep='.')[0].split(sep='_')
    mult, parallel, strand = name
    return [mult, parallel, strand]


def set_normal(aa_count):
    '''
    Ставит в соответствие сочетаниям аминокислот вероятность
    '''
    count = 0

    norm_strands = {}

    aa = aa_count.split(sep='\n')[1:len(aa_count.split(sep='\n')) - 1]

    for line in aa:
        count += int(line.split()[2])

    for line in aa:
        if line.split()[0] in norm_strands.keys():
            norm_strands[line.split()[0]][line.split()[1]] = math.log2(int(line.split()[2]) / count)
        else:
            norm_strands[line.split()[0]] = {line.split()[1]: math.log2(int(line.split()[2]) / count)}

    return norm_strands

def set_max(aa_dic):
    '''
    Находит три самых вероятных аминокислоты для данной
    '''
    res = [0,0,0]
    val = list(aa_dic.values())
    val.sort()
    for i in range(-3,0):
        for key in aa_dic.keys():
            if aa_dic[key] == val[i]:
                res[abs(i)-1] = {key: aa_dic[key]}
    return res

def sum_strands(strand1, strand2):
    '''
    Складывает вероятности в двух стрендах
    '''
    sum = {}
    strands_keys = set(list(strand1.keys()) + list(strand2.keys()))
    for key in strands_keys:
        sum[key] = strand1.get(key, 0) + strand2.get(key, 0)
    return sum

def predict(line):
    sheet = [[dict(), dict(), dict(), dict(), dict()],
             [dict(), dict(), dict(), dict(), dict()],
             [dict(), dict(), dict(), dict(), dict()]]
    pattern = [['-2L', '-1L', '0L', '1L', '2L'],
               ['-2C', '-1C', '0C', '1C', '2С'],
               ['-2R', '-1R', '0R', '1R', '2R']]
    return None

for dirname, dirnames, filenames in os.walk(path):
    for filename in filenames:
        file = open(os.path.join(dirname, filename), 'r').read()
        strand_name = get_strand_name(filename)

        aa_norm = set_normal(file)
        n = set_max(sum_strands(aa_norm['A'], aa_norm['Y']))