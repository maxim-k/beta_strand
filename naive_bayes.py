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


def classify_five(five, prob):
    '''
    Классифицирует последовательность
    '''
    single = 0
    double = 0
    c = five[2]
    s = prob['single']
    d = prob['double']

    for i in range(-2, 3):
        a = five[i + 2]
        C = str(i) + 'C'
        if i and s[C].get(c, 0):
            single += s[C][c].get(a, 0)
        elif not(i):
            single += s[C].get(c, 0)

    for i in range(-2, 3):
        a = five[i + 2]
        C = str(i) + 'C'
        if i and d[C].get(c, 0):
            double += d[C][c].get(a, 0)
        elif not(i):
            double += d[C].get(c, 0)

    #print(abs(single-double))
    if single > double:
        return 'single'
    else:
        return 'double'

def read_strand(file):
    '''
    Читает 0C
    '''
    strand = ''
    for line in file:
        if len(line)>15:
            strand += line.strip().split(sep=';')[7]
    return strand

def split_five(strand):
    strand = strand[:len(strand)//5*5]
    return (strand[0+i:5+i] for i in range(0, len(strand), 5))

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

strand_path = 'E:\\Science\\MG\Marat\\server\\pairwise\\neib\\'
single_parallel = read_strand(open(strand_path + 'neib_single_parallel.txt','r').read().split(sep='\n')[1:])
single_parallel_5 = split_five(single_parallel)
count_all = 0
count_right = 0
for five in single_parallel_5:
    count_all += 1
    if classify_five(five, prob) == 'single':
        count_right += 1
print('strands: %s\nsingles: %s\npercent: %s' % (count_all, count_right, count_right/count_all))

#single_antiparallel = read_strand(open(strand_path + 'neib_single_antiparallel.txt','r').read().split(sep='\n')[1:])
#single_antiparallel_5 = split_five(single_antiparallel)
#double_parallel = read_strand(open(strand_path + 'neib_double_parallel.txt','r').read().split(sep='\n')[1:])
#double_parallel_5 = split_five(double_parallel)
#double_antiparallel = read_strand(open(strand_path + 'neib_double_antiparallel.txt','r').read().split(sep='\n')[1:])
#double_antiparallel_5 = split_five(double_antiparallel)