#!/Users/emilydavis/anaconda/bin/python

# IMPORTS
import re
from pprint import pprint
import copy
import sys
import json
from operator import itemgetter

# GLOBALS
basepair = 60

codons = { 'ATT' : 'I',
           'ATC' : 'I',
           'ATA' : 'I',
           'CTT' : 'L',
           'CTC' : 'L',
           'CTA' : 'L',
           'CTG' : 'L',
           'TTA' : 'L',
           'TTG' : 'L',
           'GTT' : 'V',
           'GTC' : 'V',
           'GTA' : 'V',
           'GTG' : 'V',
           'TTT' : 'F',
           'TTC' : 'F',
           'ATG' : 'M',
           'TGT' : 'C',
           'TGC' : 'C',
           'GCT' : 'A',
           'GCC' : 'A',
           'GCA' : 'A',
           'GCG' : 'A',
           'GGT' : 'G',
           'GGC' : 'G',
           'GGA' : 'G',
           'GGG' : 'G',
           'CCT' : 'P',
           'CCC' : 'P',
           'CCA' : 'P',
           'CCG' : 'P',
           'ACT' : 'T',
           'ACC' : 'T',
           'ACA' : 'T',
           'ACG' : 'T',
           'TCT' : 'S',
           'TCC' : 'S',
           'TCA' : 'S',
           'TCG' : 'S',
           'AGT' : 'S',
           'AGC' : 'S',
           'TAT' : 'Y',
           'TAC' : 'Y',
           'TGG' : 'W',
           'CAA' : 'Q',
           'CAG' : 'Q',
           'AAT' : 'N',
           'AAC' : 'N',
           'CAT' : 'H',
           'CAC' : 'H',
           'GAA' : 'E',
           'GAG' : 'E',
           'GAT' : 'D',
           'GAC' : 'D',
           'AAA' : 'K',
           'AAG' : 'K',
           'CGT' : 'R',
           'CGC' : 'R',
           'CGA' : 'R',
           'CGG' : 'R',
           'AGA' : 'R',
           'AGG' : 'R',
           'TAA' : 'Stop',
           'TAG' : 'Stop',
           'TGA' : 'Stop'
        }


# FUNCTION DEFINITIONS
def validate_mutation(variant, cdna_pre, cdna_post):
    aa_position = variant[0][1:]
    aa_position = int(re.sub('[^0-9]','', aa_position)) - 1
    aa_pre      = variant[0][0]
    aa_post     = variant[0][-1]

    # build protein sequence pre-mutation
    cdna_pre_c  = copy.deepcopy(cdna_pre)
    protein_pre = ''
    while len(cdna_pre_c) != 0:
        subseq = codons["".join(cdna_pre_c[0:3])]
        protein_pre += subseq
        del cdna_pre_c[0:3]

    # build protein sequence post-mutation
    cdna_post_c = copy.deepcopy(cdna_post)
    protein_post= ''
    while len(cdna_post_c) != 0:
        subseq = codons["".join(cdna_post_c[0:3])]
        protein_post += subseq
        del cdna_post_c[0:3]
    
    result = False
    # validate based on variant action
    act = variant[1]
    if act == 'ins':
        if len(protein_pre) < len(protein_post):
            result = True
    elif act == 'del':
        if len(protein_pre) > len(protein_post):
            result = True
    elif act == 'mut':
        if protein_post[aa_position] == aa_post and protein_pre[aa_position] == aa_pre:
            result = True
        else:
            correction = aa_pre + str(aa_position) + protein_post[aa_position]
            print("warning: found mismatch %s is actually %s" % (variant[5].strip(), correction))
    return result

def reverse_primers(primers_scored):
    complements = { 'A' : 'T',
                    'T' : 'A',
                    'C' : 'G',
                    'G' : 'C',
                  }
    primers_reversed = {}
    for var, primers in primers_scored.items():
        primers_reversed[var] = []
        for p in primers:
            p_rev = copy.deepcopy(p)
            for i in range(len(p_rev)):
                p_rev[i] = complements[p_rev[i]]
            p_rev = p_rev[::-1]
            primer_couple = (p, p_rev)
            primers_reversed[var].append(primer_couple)
    return primers_reversed

def calc_gc(primer):
    g_count = 0
    c_count = 0 
    for p in primer:
        if p == 'G':
            g_count += 1
        if p == 'C':
            c_count += 1
    total_count = g_count + c_count
    my_score = (float(total_count) / float(len(primer))) * 100
    return my_score

def calc_tm(primer, gc_score):
    g_count = 0
    c_count = 0 
    a_count = 0 
    t_count = 0 
    for p in primer:
        if p == 'G':
            g_count += 1
        if p == 'C':
            c_count += 1
        if p == 'A':
            a_count += 1
        if p == 'T':
            t_count += 1
    # score = (81.5 + (.41*gc_score) - 675/len(primer)) - (float(1)/float(len(primer))) * 100
    # below formula is for basic Tm
    score = 64.9 + 41*(((g_count + c_count) - 16.4) / len(primer)) 
    return score

def count_gc(primer):
    gc_count = 0
    # only count G and C in the last 5 bases
    primer = primer[-5:]
    for p in primer:
        if p == 'G' or p == 'C':
            gc_count += 1
    return gc_count

def score_primers(primers):
    vars_final = {}
    for var, primer in primers.items():
         primers_final = []
         for p in primer:
            # calculate gc and tm scores of forward sequence
            gc_score   = calc_gc(p[0])
            tm_score   = calc_tm(p[0], gc_score)
            # calculate gc and tm scores of reverse sequence
            gc_score_r = calc_gc(p[1])
            tm_score_r = calc_tm(p[1], gc_score_r)
            # get the gc count for the forward sequence
            gc_count   = count_gc(p[0])
            # uncomment and comment lines depending if gc_count is needed
            # x = (p[0], tm_score, gc_score, gc_count)
            x = ("".join(p[0]), tm_score, gc_score)
            y = ("".join(p[1]), tm_score_r, gc_score_r)
            if tm_score < 65.0 or gc_count == 0 or p[0][-1] == 'A' or p[0][-1] == 'T':
                continue
            primers_final.append((x,y))
            vars_final[var] = primers_final
    return vars_final

# returns a sequence with the mutated base
def do_mutation(vars, cdna):
    vars_mutated = []
    for var in vars:
        cdna_copy = copy.deepcopy(cdna)
        act = var[1]
        idx = int(var[2])
        val = var[4]
        if act == "ins":
            for i in range(len(val)):
                cdna_copy.insert(idx+i, val[i])
        elif act == "del":
            for i in range(len(val)):
                del cdna_copy[idx-1]
        elif act == "mut":
            cdna_copy[idx - 1] = val
        result = validate_mutation(var, cdna, cdna_copy)
        subseq_pre = cdna[idx-basepair:idx+basepair]
        subseq_post = cdna_copy[idx-basepair:idx+basepair]
        # print var
        # print "".join(subseq_pre)
        # print "".join(subseq_post)
        vars_mutated.append((var, subseq_post))
    return vars_mutated
        
'''
Input: Takes a list of tuples(
Return: Returns candidate primers
Function: what it does
'''
def calc_primers(seqs):
    primerlen_start = 12
    primerlen_stop  = 30
    primer_candidates = {}
    for seq in seqs:
        primer_candidates[seq[0]] = []
        for i in range(primerlen_start, primerlen_stop):
            subseq = seq[1][basepair - 1 - i:basepair - 1 + i]
            primer_candidates[seq[0]].append(subseq)
        for i in range(primerlen_start, primerlen_stop):
            subseq = seq[1][basepair - 1 - i  - 1:basepair - 1 + i]
            primer_candidates[seq[0]].append(subseq)
    return primer_candidates

# generates 120 bp sequence with site to be mutated at position 60
def gen_subseqs(vars, cdna):
    list = []
    for var in vars:
        idx = int(var[2])
        subseq = cdna[idx-basepair:idx+basepair]
        # currently not supporting indexes at beginning of cdna
        if len(subseq) != 0:
            # set base to be mutated lowercase
            # subseq_original[basepair-1] = subseq_original[basepair-1].lower()
            list.append((var, subseq))
    return list

def readcdnas(fname):
    with open(fname, 'r') as f:
        contents = f.readlines()
    return list(contents[0].rstrip())

def reformatvariant(var):
    return formatted_variant

def readvariants(fname):
    with open(fname, 'r') as f:
        contents = f.readlines()
        list = []
        for line in contents:
            line_original = line
            var = line[:line.find("c.")].strip()
            line = line[line.find("c.") + 2:].rstrip()
            if "ins" in line:
                act = "ins"
                beg = line[:line.find("_")]
                end = line[line.find("_") + 1:line.find("ins")]
                val = line[line.find("ins") + 3:]
            elif "del" in line:
                act = "del"
                beg = line[:line.find("_")]
                end = line[line.find("_") + 1:line.find("del")]
                val = line[line.find("del") + 3:]
            elif ">" in line:
                act = "mut"
                beg = line[:-3].strip()
                end = line[-3]
                val = line[-1]
            else:
                print "error in file"
                exit(-1)
            list.append((var, act, beg, end, val, line_original))
    return list

def main(argv):
    # list of hgvs changes
    variants_file = argv[1]

    # reference file
    cdna_file     = argv[2]

    # vars is a list of tuples (var, act, beg, end, val) 
    vars = readvariants(variants_file)

    # cdna is a list of characters of the cdna where each character is {A, T, C, G}
    cdna = readcdnas(cdna_file)
    
    # seqs is a list of characters with the mutated cdna/base 
    seqs = do_mutation(vars, cdna)

    # sdm_seqs is a list of tuples(vars, seq)
    # sdm_seqs = gen_subseqs(vars, cdna)
    
    # primer_candidates is a dictionary of {variant : list of primer candidates)
    # where a primer candidate is of a specified length and terminates in either C or G
    primer_candidates = calc_primers(seqs)

    primers_reversed = reverse_primers(primer_candidates)

    primers_scored = score_primers(primers_reversed)

    with open("primers.txt", 'w') as f:
        for var, primer_info in primers_scored.items():
            f.write('-' * 50 + "\n")
            f.write((var[5]) + "\n")
            # for primer in primer_info:
                # primer[1] = float(primer[1])
            primer_info = sorted(primer_info, key=itemgetter(1))
            # iterate over all primer pairs
            for primer in primer_info:
                # iterate over the forward and reverse of the primer
                for p in primer:
                    f.write(str(len(p[0])) + ' ' +str(p) + '| ' + var[0] + ' SDM ' + ' R ' + p[0] + ' 25nm ' + 'STD\n')
                f.write('\n')

def usage():
    print "./makeprimers.py <variants_file> <cdna_file>"
    return 0

if __name__ == "__main__":
    if len(sys.argv) != 3:
        usage()
        exit(0)
    else:
        main(sys.argv)
