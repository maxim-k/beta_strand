#!/usr/bin/env python

__author__ = 'Maxim K'

import os
import hist
import aa_pairwise

def separate_struct_columns(line):
    '''
    Separate STRUCTURE, BP1, BP2 columns from table.
    '''
    num = line[0:5].strip()
    res_pos = line[6:11].strip()
    res_sheet_label = line[11:12].strip()
    aa = line[13:15].strip()
    secstruc = line[16:18].strip()
    helix3 = line[18:19].strip()
    helix4 = line[19:20].strip()
    helix5 = line[20:21].strip()
    geom_bend = line[21:22].strip()
    chiral = line[22:23].strip()
    beta_brige_lab1 = line[23:24].strip()
    beta_brige_lab2 = line[24:25].strip()
    bp1 = line[26:29].strip()
    bp2 = line[30:33].strip()
    beta_sheet_label = line[33:34].strip()
    acc = line[35:38].strip()
    nh_o1 = line[39:50].strip()
    o_hn1 = line[51:61].strip()
    nh_o2 = line[62:72].strip()
    o_hn2 = line[73:83].strip()
    tco = line[84:91].strip()
    kappa = line[92:97].strip()
    alpha = line[98:103].strip()
    phi = line[104:109].strip()
    psi = line[110:115].strip()
    x_ca = line[116:122].strip()
    y_ca = line[123:129].strip()
    z_ca = line[130:138].strip()
    return [num, res_sheet_label, aa, secstruc, bp1, bp2, beta_sheet_label]

def select_table(dssp):
    '''
    Select all lines down the table head line (#  RESIDUE AA STRUCTURE BP1 BP2 ...).
    '''
    return dssp[dssp.find('  #  RESIDUE') + 125:].split('\n')

def aggregate_dict(path):
    '''
    '''
    agg_dict = []
    for dirname, dirnames, filenames in os.walk(path):
        for filename in filenames:
            if '.dssp' in filename:
                agg_dict.append(os.path.join(dirname, filename))
    return agg_dict

def sel_neibor_strands(line):
    if (line[4] == '0')and(line[5] == '0'):
        return ['s0', [[line[0], '0', line[6]]]]
    elif (line[4] != '0')and(line[5] == '0'):
        return ['s1', [[line[0], line[4], line[6]]]]
    elif (line[4] == '0')and(line[5] != '0'):
        return ['s2', [[line[0], line[5], line[6]]]]
    else:
        return ['s3', [[line[0],line[4], line[6]],[line[0], line[5], line[6]]]]

def separate_beta(table):
    sheets = dict()
    sheets_list = []
    for line in table:
        if line[3] == 'E':
            if line[1] in sheets:
                if sel_neibor_strands(line)[0] in sheets[line[1]]:
                    sheets[line[1]][sel_neibor_strands(line)[0]].append(line)
                    sheets_list.append(line for line in sel_neibor_strands(line)[1])
                else:
                    sheets[line[1]][sel_neibor_strands(line)[0]] = [line]
                    sheets_list.append(line for line in sel_neibor_strands(line)[1])
            else:
                sheets[line[1]] = {sel_neibor_strands(line)[0] : [line]}
                sheets_list.append(sel_neibor_strands(line)[1])

    for sheet in sheets:
        if(not('s0' in sheets[sheet])): sheets[sheet]['s0'] = []
        if(not('s1' in sheets[sheet])): sheets[sheet]['s1'] = []
        if(not('s2' in sheets[sheet])): sheets[sheet]['s2'] = []
        if(not('s3' in sheets[sheet])): sheets[sheet]['s3'] = []

    # Pairwise-threads list for accordance check
    thread_list = []
    for tmp_line in sheets_list:
        for tmp_line2 in tmp_line:
            thread_list.append(tmp_line2)

    # Wrong accordance correction
    non_acc = check_accordance(thread_list)
    for na in non_acc:
        for sheet in sheets:
            for thread in sheets[sheet]:
                for line in sheets[sheet][thread]:
                    if (line[0] == na[0]) and (line[6] == na[2]) and ((line[4] == na[1]) or (line[5] == na[1])):
                        line_index = sheets[sheet][thread].index(line)
                        thread_index = sheets[sheet][thread][line_index].index(na[1])
                        sheets[sheet][thread][line_index][thread_index] = '0'
                        if (thread_index == 4) and (sheets[sheet][thread][line_index][5] != '0'):
                            sheets[sheet]['s2'].append(sheets[sheet][thread][line_index])
                            sheets[sheet]['s2'].sort()
                        elif (thread_index == 5) and (sheets[sheet][thread][line_index][4] != '0'):
                            sheets[sheet]['s1'].append(sheets[sheet][thread][line_index])
                            sheets[sheet]['s1'] = [[sheets[sheet][thread][line_index]]]
                        else:
                            sheets[sheet]['s0'].append(sheets[sheet][thread][line_index])
                            sheets[sheet]['s0'].sort()
                        sheets[sheet][thread].pop(line_index)
    print(non_acc)

    return sheets

def split_double(table):
    left = []
    right = []
    for line in table:
        left.append([line[0], line[1], line[2], line[3], line[4], '0', line[6]])
        right.append([line[0], line[1], line[2], line[3], '0', line[5], line[6]])
    left_parallel = [list(z) for z in (zip(left, set_parallel(left, 4, 6)))]
    right_parallel = [list(z) for z in (zip(right, set_parallel(right, 5, 6)))]
    left_parallel += right_parallel
    return left_parallel

def check_accordance(table):
    '''
    If residue #n is connected with residue #k, is residue #k is connected with residue #n?
    Return list of wrong accorded threads
    '''
    non_acc = [line for line in table if not([line[1], line[0], line[2]] in table)and(line[1] != '0')]
    return non_acc

def set_parallel(thread_list, col_num, sheet_label):
    '''
    Return 'True' if beta-strands are parallel. Else - 'False'.
    '''
    p = []
    if(len(thread_list) == 1):
        thread_list = thread_list[0]
        p.append(False)
    else:
        if thread_list:
            prev = thread_list[0]
            for line in thread_list:
                if(int(line[col_num]) - int(prev[col_num]) <= 0)and(int(line[col_num]) != 0)and(line[sheet_label] == prev[sheet_label]):
                    p.append(False)
                elif(line[sheet_label] == prev[sheet_label]):
                    p.append(True)
                prev = line

    if(len(thread_list) >= 3):
        for index in range(1, len(p) - 1):
            prev = p[index - 1]
            curr = p[index]
            next = p[index + 1]

            if(prev or curr)and(curr or next)and(not(prev ^ next)):
                p[index - 1] = p[index]
                p[index] = p[index + 1]

    return p

a = aggregate_dict('E:\\Science\\MG\\Marat\\data\\jp')
processed = []

for file in a:
    print(file)
    table = ([separate_struct_columns(line) for line in select_table(open(file).read())])
    s = separate_beta(table)
    for sheet in s:
        print(sheet)
        for thread in s[sheet]:
            print(thread)
            if(thread == 's1'):
                s[sheet][thread] = [list(z) for z in (zip(s[sheet][thread], set_parallel(s[sheet][thread], 4, 6)))]
            elif(thread == 's2'):
                s[sheet][thread] = [list(z) for z in zip(s[sheet][thread], set_parallel(s[sheet][thread], 5, 6))]
            elif(thread == 's3')and(s[sheet][thread]):
                s4 = split_double(s[sheet][thread])
            for line in s[sheet][thread]:
                print(line)
        if s4:
            s[sheet]['s4'] = s4
    if s:
        processed.append(s)

hist.count_aa(processed, 'E:\\Science\\MG\Marat\\hist')
aa_pairwise.count_pairwise(processed, 'E:\\Science\\MG\Marat\\pairwise')