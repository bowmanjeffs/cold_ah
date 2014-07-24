# -*- coding: utf-8 -*-
"""
Created on Wed Apr 09 10:41:27 2014

@author: Jeff

If you have done a partial run or have an old tre.dist file that you want to
write over, you will need to delete before you run the script!  Otherwise
script skips the long process of re-running.

"""
#### must set these parameters ####
dist_mat = '../FA_desaturase_mds_points.txt' ## existing distance matrix, must be gzipped
name = 'FA_desaturase_pro_aligned.filter.tre' ## name of tree to build distance matrix from
###################################

#dist_mat = 'DUF900_mds_points.txt'
#name = 'DUF900_pro_aligned.filter.tre'

from Bio import Phylo
import math
import gzip
import os

tree = Phylo.read(name, 'newick')   
t = ((len(tree.get_terminals()) ** 2) / 2) - len(tree.get_terminals())
tree_dist = {}
conf = {}

if name+'.dist.gz' not in os.listdir('.'):
      
    done = set()
    
    with gzip.open(name+'.dist.gz', 'wb') as output, open(name+'.conf', 'w') as conf_out:
        
        n = 0
        
        for clade1 in tree.get_terminals():
            name1 = clade1.name
            temp_conf = []
            path = tree.get_path(clade1)
            for clade in path:
                tc = clade.confidence
                if tc != None:
                    temp_conf.append(float(tc))
            try:
                conf[name1] = sum(temp_conf) / len(temp_conf)
            except ZeroDivisionError:
                conf[name1] = 0
            print >> conf_out, name1+'\t'+str(conf[name1])
            
            for clade2 in tree.get_terminals():
                name2 = clade2.name
                if name1 != name2:
                    if (name1, name2) not in done:
                        n = n + 1
                        dist = tree.distance(clade1, clade2)
                        print 'tree dist', n, 'out of', t, dist
                        print >> output, name1 + '\t' + name2 + '\t' + str(dist)               
                        done.add((name2, name1))
                        done.add((name1, name2))
                        tree_dist[name1+'\t'+name2] = dist
                        
else:
    with gzip.open(name+'.dist.gz', 'rb') as dist_in:
        for line in dist_in:
            line = line.rstrip()
            line = line.split()
            name1 = line[0]
            name2 = line[1]
            dist = line[2]
            tree_dist[name1+'\t'+name2] = dist
            
    with open(name+'.conf', 'r') as conf_in:
        for line in conf_in:
            line = line.rstrip()
            line = line.split()
            name1 = line[0]
            score = line[1]
            conf[name1] = score
        
mds_points = {}
mds_dist = {}
mds_done = set()
                
with open(dist_mat, 'r') as mds:
    for line in mds:
        if line.startswith('NMDS') == False:
            line = line.split()
            mds_points[line[0]] = (line[1], line[2])

n = 0            
for key1 in mds_points.keys():
    for key2 in mds_points.keys():
        if (key1, key2) not in mds_done:
            n = n + 1
            x1 = float(mds_points[key1][0])
            y1 = float(mds_points[key1][1])
            x2 = float(mds_points[key2][0])
            y2 = float(mds_points[key2][1])
            dist = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
            mds_dist[key1 + '\t' + key2] = dist
            mds_done.add(key1 + '\t' + key2)
            mds_done.add(key2 + '\t' + key1)
            print 'mds dist', n, 'out of', t * 2, dist

with gzip.open(name+'.compare.txt.gz', 'wb') as output:
    for key in tree_dist.keys():
        
        temp_key = key.split()
        temp_conf_1 = float(conf[temp_key[0]])
        temp_conf_2 = float(conf[temp_key[1]])
        local_conf = min((temp_conf_1, temp_conf_2))
        
        try:
            tdist = tree_dist[key]
            mdist = mds_dist[key]
            print >> output, key + '\t' + str(local_conf) + '\t' + str(tdist) + '\t' + str(mdist)
            print key + '\t' + str(tdist) + '\t' + str(mdist)
        except KeyError:
            tdist = tree_dist[key]
            temp = key.split()
            alt_key = key[1] + '\t' + key[0]
            try:
                mdist = mds_dist[alt_key]
                print >> output, key + '\t' + str(local_conf) + '\t' + str(tdist) + '\t' + str(mdist)
                print key + '\t' + str(tdist) + '\t' + str(mdist)              
            except KeyError:
                print 'error', key
        
    
                    