# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
C:\Users\jeff\.spyder2\.temp.py
"""
import subprocess

from Bio.SeqUtils import ProtParam as PP

def calc_window_params(prot):
    seq = PP.ProteinAnalysis(prot)
    flex_list = PP.ProteinAnalysis.flexibility(seq)
    gravy = []
    iso = []
    arom = []
    ai = []
    
    for i in range(0, len(prot)):
        win = prot[i:i + 9]
        win_seq = PP.ProteinAnalysis(win)
        if len(win) == 9:

            tg = str(PP.ProteinAnalysis.gravy(win_seq))
            ti = str(PP.ProteinAnalysis.isoelectric_point(win_seq))
            tar = str(PP.ProteinAnalysis.aromaticity(win_seq))
            taa = PP.ProteinAnalysis.get_amino_acids_percent(win_seq)
            tai = str((100 * taa['A']) + (2.9 * 100 * taa['V']) + (3.9 * (100 * taa['I'] + 100 * taa['L'])))    

            gravy.append(tg)
            iso.append(ti)
            arom.append(tar)
            ai.append(tai)    
    
    return flex_list, gravy, iso, arom, ai


groups = {}    
with open('all_pfams.groups', 'r') as group_file:
    for line in group_file:
        line = line.rstrip()
        line = line.split()
        groups[line[0]] = line[1], line[2], line[3]
        ## prot = strain, group, pair

from Bio import SeqIO

clusters = ['p450_cluster_0', 'p450_cluster_1', 'Bac_luciferase_cluster_0', 'FA_desaturase_cluster_0', 'FA_desaturase_cluster_1', 'Pyr_redox_3_cluster_0']

for cluster in clusters:
    clustalo = subprocess.Popen('clustalo --force -i ' + cluster + '_pro.fasta -o ' + cluster + '_pro_aligned.fasta', shell = True)
    clustalo.communicate()
    
    with open(cluster + '_param_alignment.txt', 'w') as param_out, open(cluster + '_param_alignment.map', 'w') as map_out:    
        for record in SeqIO.parse(cluster + '_pro_aligned.fasta', 'fasta'):
                
            strain = groups[record.id][0]
            group = groups[record.id][1]
            pair = groups[record.id][2]
            
            seq = str(record.seq)
            degap_seq = seq.replace('-', '')
            flex, gravy, iso, arom, ai = calc_window_params(degap_seq)
            
            output_flex = []
            n = 0 ## number of nongap positions
            
            for i,p in enumerate(seq):
                if p != '-':
                    if (n - 4) >= 0:
                        try:
                            output_flex.append(str(flex[n - 4]))
                        except IndexError:
                            output_flex.append('NA')
                    n = n + 1
                    print >> map_out, record.id, n, i
                        
                else:
                    output_flex.append('NA')
            print len(output_flex)        
            output_flex = ','.join(output_flex)
            print >> param_out, record.id + ',' + strain + ',' + group + ',' + pair + ',' + 'flex' + ',' + output_flex
 