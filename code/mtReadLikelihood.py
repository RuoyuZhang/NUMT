#!/usr/bin/env python
import re
import math
import pysam
from optparse import OptionParser

# function to calculate Pr(a|g) given reference is g, observed a; only consider SNP
def cal_p_sub(ref_base, true_base, h=0.04, pts=2.0/3):
    p_sub = 0
    
    ts = ["AG","GA","CT","TC"]
    if ref_base == true_base:
        p_sub = 1 - h
    elif (ref_base + true_base) in ts:
        p_sub = h * pts
    else:
        p_sub = 1/2*(h - h*pts)
    
    return p_sub


# calculate Pr(x|a) given nucleotide is a, sequenced base is x
def cal_p_error(true_base, sequenced_base, Q):
    p_error = 10**(-Q/10.0)
    if true_base == sequenced_base:
        return 1-p_error
    else:
        return p_error/3
  

# function to calculate likehood for one site
def prob_site(ref_base, seq_base, Q, h, pts, ch, position, var_freq):
    #Q=100
    NT = ["A","T","G","C"]
    freq=dict()
    maf=0.0
    for a in NT:
        key = str(position) + "-" + ref_base + "-" + a
        if (ch in var_freq and key in var_freq[ch]):
            maf = maf + float(var_freq[ch][key])
            freq[ref_base + a]= float(var_freq[ch][key])
        else:
            freq[ref_base + a]=cal_p_sub(ref_base=ref_base, true_base=a, h=h,pts=pts)
    
    freq[ref_base + ref_base]= 1-maf
    
    prob_all = 0.0
    
    if (maf == 0):
        for a in NT:
            prob_all = prob_all + cal_p_sub(ref_base=ref_base, true_base=a, h=h,pts=pts) * cal_p_error(true_base=a,sequenced_base=seq_base,Q=Q)
    else:
        for a in NT:
            prob_all = prob_all + freq[ref_base+a] * cal_p_error(true_base=a,sequenced_base=seq_base,Q=Q)   
    
    return prob_all


# function to calculate log likelihood across entrire read
def prob_read(read, qual, mutation_matrix, h, pts, ch, st, var_freq):
    read_tides = list(read)
    #qual_tides = list(qual)
    
    loglike = 0.0
    i = 0
    indel = []
    for mu in mutation_matrix:
        if mu == 1:
            Q = qual[i]
            ref_seq = read_tides[i]
        elif mu.find('-') != -1 or mu.find('del') != -1:
            if mu.find('-') != -1:
                position = st + i
                mu = str(position) + '-' + mu
                indel.append(mu)
            continue # skip deletion
        elif mu.find('+') != -1 or mu.find('ins') != -1:
            if mu.find('+') != -1:
                position = st + i
                mu = str(position) + '-' + mu
                indel.append(mu)
            i = i + 1
            continue # skip insertion
        else:
            Q = qual[i]
            ref_seq = mu        
                
        position = st + i
        prob_this_site=prob_site(ref_base=ref_seq,seq_base=read_tides[i],Q=Q,h=h,pts=pts,ch=ch,position=position,var_freq=var_freq)
        if prob_this_site > 0:
            loglike = loglike + math.log(prob_this_site)
        
        i = i + 1
            
    return loglike


# Creat variants frequency dictionary
def cal_var_freq(mt_file,numt_file):
    var_dict = dict()
    var_dict["chrM"] = dict() # mitochondrial varaints
    with open (mt_file) as mt_fh:
        for i,line in enumerate(mt_fh):
            data=line.split()
            key = data[0] + "-" + data[1]
            var_dict["chrM"][key] = data[2]
                
    chrs = map(str, range(1,23))
    chrs = list(chrs)
    
    chrs.append('X')
    chrs.append('Y')
    for ch in chrs:
        var_dict["chr" + ch] = dict() 
    
    with open (numt_file) as numt_fh:
        for i,line in enumerate(numt_fh):
            data=line.split()
            if ("chr" + data[0]) in var_dict:
                key = data[1] + "-" + data[2]
                var_dict["chr"+str(data[0])][key] = data[3]
        
    return var_dict


# function to parse MD string
def parseMD(MD):
    mismatches, pos = [], 0
    for n in re.findall('(\d+|[a-zA-Z|^]+)', MD):
        if n.isalpha():
            mismatches.append(n)
            pos += 1
        #ignore deletions i.e. ^[A-Z] (isalpha() is False)
        elif not n.startswith('^'):
            pos += int(n)
            [mismatches.append(1) for i in range(int(n))]
                   
    #print(mismatches)
    return mismatches


# function to parse CIGAR and MD, output final mutation matrix
def parseCIGAR(seq, CIGAR, MD):    
    match = parseMD(MD)
    seqPos = 0
    
    for op, length in CIGAR:
        if op == 0: #match
            seqPos += length 
   
        elif op == 1: #insertion
            if length > 1:
                match = match[:seqPos] + ["+" + seq[seqPos:(seqPos+length)]] + ['ins_ext',] * (length - 1) + match[seqPos:]
            else:
                match = match[:seqPos] + [+1*length] + match[seqPos:]
            #match = match[:seqPos] + ["+" + seq[seqPos:(seqPos+length)]] + match[seqPos:]
            seqPos += length 
   
        elif op == 2: #deletion
            if length > 1:
                match = match[:seqPos] + [str(-length)] + ['del_ext',] * (length - 1) + match[seqPos:]
            else:
                match = match[:seqPos] + [str(-1*length)] + match[seqPos:]
            seqPos += length
        elif op == 4: #soft clip
            continue
    
    mutation = dict()
    for i,item in enumerate(match):
        if item !=1 and item !='del_ext':
            mutation[i+1] = item
             
    return match


# function to convert MD:Z tag to mutation matrix, assume no insertions and deletions
def convert_mutation_matrix(MDstring):
    mutations = re.findall('[ATGCN]',MDstring)
    numbers = re.sub('[ATGCN]',',',MDstring)
    numbers = numbers.split(',')
    
    position = 0
    mutation_matrix=dict()
    for i, num in enumerate(numbers):
        if num == '':
            position = position + 1
            mutation = mutations[0] 
            mutation_matrix[position] = mutation
        elif i == len(mutations):
            break
        else:
            position = position + int(num) + 1
            mutation = mutations[i] 
            mutation_matrix[position] = mutation

    return mutation_matrix


# function to parse read
def parse_read(r,h,pts_m,pts_n,var_freq):
    seq = r.query_alignment_sequence
    qual = r.query_alignment_qualities
    MD = r.get_tag('MD')
    ch = r.reference_name
    st = r.reference_start + 1
    pts = (pts_m if ch == "chrM" else pts_n)
    
    #mutation_matrix=convert_mutation_matrix(MD)
    mutation_matrix = parseCIGAR(seq=seq, CIGAR = r.cigartuples, MD=MD)
    loglike = prob_read(read=seq,qual=qual,mutation_matrix=mutation_matrix,h=h,pts=pts,ch=ch,st=st,var_freq=var_freq)
    return loglike


# function to summarize a read group
def summary_out(read1group,read2group,oh1,oh2,h,pts_m,pts_n,var_freq,mode):
    r1dict = dict()
    r2dict = dict()
    r1mt = -10000.0
    r2mt = -10000.0
    r1mtr = None
    r2mtr = None
    
    for r1 in read1group:
        loglike = parse_read(r=r1,h=h,pts_m=pts_m,pts_n=pts_n,var_freq=var_freq)
        r1dict[r1.reference_name + '-' + str(r1.reference_start)] = [loglike, r1]
        if r1.reference_name == 'chrM': r1mt = loglike; r1mtr = r1
            
    for r2 in read2group:
        loglike = parse_read(r=r2,h=h,pts_m=pts_m,pts_n=pts_n,var_freq=var_freq)
        r2dict[r2.reference_name + '-' + str(r2.reference_start)] = [loglike, r2]
        if r2.reference_name == 'chrM': r2mt = loglike; r2mtr = r2
    
    r1list = sorted(r1dict.iteritems(), key=lambda x : x[1],reverse=True)
    r2list = sorted(r2dict.iteritems(), key=lambda x : x[1],reverse=True)
    
    if mode == 'se':
        if r1list != []:
            
            readname = r1list[0][1][1].qname
            [map_type, chrM_like, NUMT_like, chrM_mis, NUMT_mis] = ['NA'] * 5
            for x in r1list[0:2]:
                if x[1][1].reference_name == 'chrM':
                    chrM_like, chrM_mis = x[1][0], x[1][1].get_tag('NM')
                else:
                    NUMT_like, NUMT_mis = x[1][0], x[1][1].get_tag('NM')
            
            map_type = 'multi' if chrM_like != 'NA' and NUMT_like != 'NA' else 'unique'
            
            oh1.write("%s %s %s %s %s %s\n" % (readname, map_type, str(chrM_like), str(chrM_mis), str(NUMT_like), str(NUMT_mis)))
            
        return ""       
    
    if r1list != [] and r2list != [] :
        
        readname = r1list[0][1][1].qname
        [map_type, chrM_like1, chrM_mis1, chrM_like2, chrM_mis2, NUMT_like1, NUMT_mis1, NUMT_like2, NUMT_mis2] = ['NA'] * 9
        for x in r1list[0:2]:
            if x[1][1].reference_name == 'chrM':
                chrM_like1, chrM_mis1 = x[1][0], x[1][1].get_tag('NM')
            else:
                NUMT_like1, NUMT_mis1 = x[1][0], x[1][1].get_tag('NM')
        
        for y in r2list[0:2]:
            if y[1][1].reference_name == 'chrM':
                chrM_like2, chrM_mis2 = y[1][0], y[1][1].get_tag('NM')
            else:
                NUMT_like2, NUMT_mis2 = y[1][0], y[1][1].get_tag('NM')
            
        map_type = 'multi' if chrM_like1 != 'NA' and NUMT_like1 != 'NA' else 'unique'
        oh1.write("%s %s %s %s %s %s %s %s %s %s\n" \
                  % (readname, map_type, str(chrM_like1), str(chrM_mis1), str(chrM_like2), str(chrM_mis2), str(NUMT_like1), str(NUMT_mis1), str(NUMT_like2), str(NUMT_mis2)))
        
        return ''


def main():
    usage = "usage: %prog [options] \n"
    parser = OptionParser(usage=usage)
    parser.add_option("--sam", default="/home/fs01/rz253/project/BTRY4840/NUMT/data/prepare.data/PCR.1000.sam",help="input sam file")
    parser.add_option("--stat", default="/home/fs01/rz253/project/BTRY4840/NUMT/data/prepare.data/PCR.stat",help="output likelihood file")
    #parser.add_option("--outsam", default="/home/fs01/rz253/project/BTRY4840/NUMT/data/prepare.data/PCR.rm.NUMT.sam",help="output sam file")
    parser.add_option("--mt_file", default="/home/fs01/rz253/project/BTRY4840/NUMT/data/prepare.data/mt.variants.frequency.txt",help="mtDNA variants file")
    parser.add_option("--numt_file", default="/home/fs01/rz253/project/BTRY4840/NUMT/data/prepare.data/numt.variants.frequency.txt",help="numt variants file")
    
    parser.add_option("-s", type=float, default=0.0005, help="DNA substitution rate")
    parser.add_option("--pts_m", type=float, default=23/24.0, help="mtDNA transition probability")
    parser.add_option("--pts_n", type=float, default=2.43/3.43, help="NUMT transition probability")
    parser.add_option("--mode", type=str, default="se", help="Singe end or pair end? [se]")

    (o, args) = parser.parse_args()
    
    oh1 = open(o.stat,"w+")
    if o.mode == 'se':
        oh1.write("read_name map chrM_like chrM_mis NUMT_like NUMT_mis\n")
    else:
        oh1.write("read_name map chrM_like1 chrM_mis1 chrM_like2 chrM_mis2 NUMT_like1 NUMT_mis1 NUMT_like2 NUMT_mis2\n")
        
    h = o.s
    
    mysam = o.sam
    
    pts_m = o.pts_m
    pts_n = o.pts_n
    var_freq = cal_var_freq(o.mt_file, o.numt_file) 
    
    samfile = pysam.AlignmentFile(mysam,"r")
    oh2 = pysam.AlignmentFile(o.outsam, "w", template=samfile)
    #oh3 = pysam.Samfile(o.outsam.replace(".sam","NUMTs.sam"), "w", template = samfile )
    
    # iterate reads
    read1group = []
    read2group = []
    readname = ""
    
    mt_total = 0
    numt_total = 0
    
    for r in samfile.fetch():
        #if re.findall('[IDS]',r.cigarstring) != []:
        #    continue
        if readname == "" or r.qname != readname:
            if readname == "":
                readname = r.qname
                if o.mode == 'se':
                    read1group.append(r)
                elif r.is_read1:
                    read1group.append(r)
                elif r.is_read2:
                    read2group.append(r)
                continue
            
            res = summary_out(read1group=read1group,read2group=read2group,oh1=oh1,oh2=oh2,h=h,pts_m=pts_m,pts_n=pts_n,var_freq=var_freq,mode = o.mode)
            if res == 'mt':
                mt_total = mt_total + 1
            elif res == 'numt':
                numt_total = numt_total + 1
        
            read1group = []
            read2group = []
            readname = r.qname
        
        if o.mode == 'se':
            read1group.append(r)
        elif r.is_read1:
            read1group.append(r)
        elif r.is_read2:
            read2group.append(r)
            
    res = summary_out(read1group=read1group,read2group=read2group,oh1=oh1,oh2=oh2,h=h,pts_m=pts_m,pts_n=pts_n,var_freq=var_freq, mode = o.mode)
    if res == 'mt':
        mt_total = mt_total + 1
    elif res == 'numt':
        numt_total = numt_total + 1
    
    #print "%i %i %f" %(mt_total, numt_total, float(numt_total)/(mt_total + numt_total))


if __name__ == '__main__':
    #test1()
    main()
    