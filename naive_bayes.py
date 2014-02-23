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

def prob_2d(prob):
    '''
    Считает двумерную дискретную плотность распределения для соседей 0C
    '''
    for C in prob.keys():
        if(C != '0C'):
            for a in prob[C].keys():
                for aa in prob[C][a].keys():
                    prob[C][a][aa] += prob['0C'].get(a, 0)

def sum_prob(five):
    return None

prob = {'single': dict(), 'double': dict()}
for dirname, dirnames, filenames in os.walk(path):
    for filename in filenames:
        if 'C' in filename:
            if 'single' in filename:
                file = open(os.path.join(dirname, filename), 'r').read()
                strand_name = get_strand_name(filename)
                aa_norm = set_normal(file)
                prob['single'][strand_name[2]] = set_normal(file)
            elif 'double' in filename:
                file = open(os.path.join(dirname, filename), 'r').read()
                strand_name = get_strand_name(filename)
                aa_norm = set_normal(file)
                prob['double'][strand_name[2]] = set_normal(file)
prob['double']['0C'] = set_normal_0(open('E:\\Science\\MG\\Marat\\server\\hist\\single.txt', 'r').read())
prob['single']['0C'] = set_normal_0(open('E:\\Science\\MG\\Marat\\server\\hist\\double.txt', 'r').read())

for key in prob.keys():
    prob_2d(prob[key])

print()