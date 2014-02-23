__author__ = 'Maxim K'

import os, math

path = 'E:\\Science\\MG\Marat\\pairwise\\neib_pairwise'
path_0C = 'E:\\Science\\MG\\Marat\\server\\hist\\single.txt'

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

def set_normal_0(aa_count):
    '''
    Ставит в соответствие аминокислотам вероятность
    '''
    count = 0
    norm_strands = {}
    aa = aa_count.split(sep='\n')[1:len(aa_count.split(sep='\n')) - 1]

    for line in aa:
        count += int(line.split()[1])

    for line in aa:
        norm_strands[line.split()[0]] = math.log2(int(line.split()[1]) / count)

    return norm_strands

prob = dict()

for dirname, dirnames, filenames in os.walk(path):
    for filename in filenames:
        if 'C' in filename:
            file = open(os.path.join(dirname, filename), 'r').read()
            strand_name = get_strand_name(filename)
            aa_norm = set_normal(file)
            prob[strand_name[2]] = set_normal(file)

prob['0C'] = set_normal_0(open(path_0C, 'r').read())

print()