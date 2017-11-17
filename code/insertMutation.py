#! /usr/bin/env python
import os
import sys
from Bio import SeqIO
from optparse import OptionParser
import numpy as np


# function to read reference fasta file
def readFasta(fastafn):
    fastafh = open(fastafn, 'r')
    
    # parse the FASTA file
    seqs = {}
    for s in list(SeqIO.parse(fastafh, format='fasta')):
        seqs[s.name] = str(s.seq)

    return seqs


# parse mutation file
def parseMutation(mutation_file, ref_type = 'mt'):
    subs = dict()
    indels = dict()
    
    if ref_type != 'mt':
        chrs = map(str, range(1,23))
        chrs = list(chrs)
    
        chrs.append('X')
        chrs.append('Y')
        for ch in chrs:
            subs["chr" + ch] = dict()
            indels["chr" + ch] = dict()
    else:
        subs['chrM'] = dict()
        indels["chrM"] = dict()
        
        
    mutation_fh = open(mutation_file,'r')
    for line in mutation_fh:
        [ch, pos, var, freq] = [None] * 4
        if ref_type == 'mt':
            ch = 'chrM'
            [pos, var, freq] = line.split()
        else:
            [ch, pos, var, freq] = line.split()
            ch = 'chr' + ch
        
        if ch not in subs:
            continue
        pos, freq = int(pos), float(freq)
        [ref, new] = var.split('-')
        if len(new) == len(ref):
            subs[ch][pos] = [ch, pos, ref, new, freq]
        else:
            indels[ch][pos] = [ch, pos, ref, new, freq]

    return subs, indels


def insertSNPs(ref, mut, freq, force_insert):
    if force_insert == 'ins':
        freq = 1
    elif force_insert == 'no':
        freq =0
    
    if np.random.binomial(1,p=freq):
        return(mut)
    else:
        return(ref)


def insertIndels(ref, mut, freq, force_insert):
    # consider this later
    return ref


def main():
    usage = "usage: %prog [options] \n"
    parser = OptionParser(usage=usage)
    parser.add_option("--ref", default="/home/fs01/rz253/project/BTRY4840/NUMT/data/prepare.data/NUMTs.fa",help="input reference file")
    parser.add_option("--mutation", default="/home/fs01/rz253/project/BTRY4840/NUMT/data/prepare.data/numt.variants.frequency.txt",help="input mutation file")
    parser.add_option("--out", default="/home/fs01/rz253/project/BTRY4840/NUMT/data/prepare.data/test.fa",help="output fasta file")
    parser.add_option("--ref_type", type=str, default="numt", help="numt or mt? [numt]")
    parser.add_option("--insert_indel", action="store_false", help="insert indel? [F]")
    parser.add_option("--force_insert", type=str, default = 'file', help="force insert? [file|ins|no]")
    
    (o, args) = parser.parse_args()
    
    # parse mutation file
    [subs,indels]=parseMutation(o.mutation,ref_type=o.ref_type)
    
    #print(subs)
    # parse reference file
    ref_dict = readFasta(o.ref)
    out_fh = open(o.out,'w')
    
    # begin insert mutations
    for key, item in ref_dict.iteritems():
        [ch, ranges] = key.split(':')
        start = int(ranges.split('-')[0])
        string = '>' + key + '\n'
        for i, base in enumerate(item):
            pos = start + i
            new = base
            #print pos
            if pos in subs[ch]:
                #print(pos)
                [ref, mut, freq] = subs[ch][pos][2:5]
                new = insertSNPs(ref, mut, freq, o.force_insert)
            
            if o.insert_indel:
                if pos in indels[ch]:
                    [ref, mut, freq] = indels[ch][pos][2:5]
                    new = insertIndels(ref, mut, freq, o.force_insert)
            
            string = string + new
        
        string = string + '\n'
        out_fh.write(string) 
       
    
if __name__ == '__main__':
    main()