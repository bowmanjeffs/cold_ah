# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 10:50:36 2013

@author: Jeff

Takes as input a single aligned fasta file.  be sure to set appropriate gap
characters and the minimum number of positions.

Script takes some code from filter_seqs.py.  Main difference is that it sets
min and max to a select group of sequences and eliminates everything that
doesn't overlap in the alignment.

python filter_seqs_selective.py <file in>

"""
gap_character = ['-', '.']
min_position = 100
        
import re
import sys

try:
    file_in = sys.argv[1]
except IndexError:
    file_in = 'alignments/DUF900_pro_aligned.fasta'  ## use this for testing

name = re.sub('.fasta', '', file_in)
log = open(name+'.filter.log', 'w')

for character in gap_character:
    print >> log, 'gap = '+character
    
print >> log, '\n',

uninames = set()
with open('uniprot_alkane_catabolism_pfams_table.txt', 'r') as uniprot_in:
    for line in uniprot_in:
        line = line.split('\t')
        uniname = line[1]
        print 'uniname =', uniname
        uninames.add(uniname.rstrip())

## find start and end for each sequence
                    
start_end = {}
uni_start_end = {}
unistarts = set()
uniends = set()
starts = set()
ends = set()

print 'filtering to last start and first end'
print >> log, 'filtering to last start and first end'

with open(file_in, 'r') as fasta_in:
    l = 0
    seqname = ''
    lines = ''
    
    for line in fasta_in:
        if line.startswith('>') == True:
            if l != 0:
                rlines = lines[::-1]
                
                for i,p in enumerate(lines):
                    if p not in gap_character:
                        start = i - 1
                        break
                        
                for i,p in enumerate(rlines):
                    if p not in gap_character:
                        end = len(rlines) - i
                        break
                                
                print l, seqname, start, end
                print >> log, l, seqname, start, end
                if seqname.split()[0] in uninames:
                    uni_start_end[seqname] = start, end
                    unistarts.add(start)
                    uniends.add(end)
                start_end[seqname] = start, end
                starts.add(start)
                ends.add(end)
                lines = ''
            l = l + 1    
            seqname = line.rstrip()
            seqname = seqname.strip('>')
            
        else:
            line = line.rstrip()
            lines = lines + line


## make sure you get the last line!

for i,p in enumerate(lines):
    if p not in gap_character:
        start = i - 1
        break
        
for i,p in enumerate(rlines):
    if p not in gap_character:
        end = len(rlines) - i
        break 
                
print l, start, end
print >> log, l, seqname, start, end
if seqname.split()[0] in uninames:
    uni_start_end[seqname] = start, end
    unistarts.add(start)
    uniends.add(end)
start_end[seqname] = start, end
starts.add(start)
ends.add(end)

## now you've finished with the last line

for key in uni_start_end.keys():
    print >> log, 'uniname =', key

## set min allowable end point to be uniprot max start + 50
## set max allowable start point to be uniprot min end - 50
## this should guarantee a 100 residue overlap           
max_start = min(uniends) - min_position
min_end = max(unistarts) + min_position

print >> log, l, 'sequences evaluated'   
print >> log, 'max start / min end =', max_start, '/', min_end 

## eliminate all seqs that don't overlap with uniprot

bad = set()

for key in start_end.keys():
    if start_end[key][0] > max_start:
        bad.add(key)
        print >> log, key, 'is bad', start_end[key][0], start_end[key][1] 
    elif start_end[key][1] < min_end:
        bad.add(key)
        print >> log, key, 'is bad', start_end[key][0], start_end[key][1] 
 
## look for gap only columns among good sequences

gap = {} ## key = position, value = number of seqs position is gap
       
with open(file_in, 'r') as fasta_in:
    seqname = ''
    nseq = 0
    lines = ''
    for line in fasta_in:
        if line.startswith('>') == False:
            line = line.rstrip()
            lines = lines + line
        
        else:
            nseq = nseq + 1
            if nseq != 1:
                if seqname not in bad:
                    print nseq, seqname, 'is good'
                    
                    for i,p in enumerate(lines):
                        if p in gap_character:
                            try:
                                temp = gap[i]
                                temp = temp + 1
                                gap[i] = temp
                            except KeyError:
                                gap[i] = 1
            lines = ''
            seqname = line.rstrip()
            seqname = seqname.strip('>')


if seqname not in bad:
    print nseq, seqname
    
    for i,p in enumerate(lines):
        if p in gap_character:
            try:
                temp = gap[i]
                temp = temp + 1
                gap[i] = temp
            except KeyError:
                gap[i] = 1
                            
## determine which positions are used
                    
p_use = set()

for i in range(0, max(ends)):
        if i <= max_start:
            if i >= min_end:
                try:
                    if gap[i] < (len(start_end.keys()) - len(bad)):
                        p_use.add(i)
                except KeyError:
                    p_use.add(i)
            
with open(file_in, 'r') as fasta_in, open(name+'.filter.fasta', 'w') as fasta_out:
    seq = ''
    seq_out = ''
    keep = False
    
    for line in fasta_in:
        if line.startswith('>'):

            if keep == True:
                print seqname, 'generating filtered seq'
                for i,p in enumerate(seq):
                    if i in p_use:
                        seq_out = seq_out + p
                print >> fasta_out, seq_out
                seq = ''
                seq_out = ''
            
            seqname = line.rstrip()
            seqname = seqname.strip('>')
            
            if seqname in bad:
                print >> log, seqname, 'not added to filtered fasta'
                keep = False
            
            else:
                keep = True
                seq = ''
                print >> log, seqname, 'is added to filtered fasta'
                print >> fasta_out, line,
        
        elif keep == True:
            line = line.rstrip()
            seq = seq + line

    if keep == True:
        print seqname, 'generating filtered seq'
        for i,p in enumerate(seq):
            if i in p_use:
                seq_out = seq_out + p
        print >> fasta_out, seq_out
                
print >> log, '\n'
print >> log, 'filter:'
 
flter = ''           
for i in range(0, max(ends)):
    if i in p_use:
        flter = flter + '1'
    else:
        flter = flter + '0'

print >> log, flter

log.close()
                