#!/usr/bin/env python
# coding: utf-8

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def dict_seqs(seqs):
    '''
    Indexing fasta sequences into seqrecords and writing id and sequence into a dictionary
    '''
    seq_dict = {}
    for key, SeqRecord in seqs.items():
        seq_dict[key] = str(SeqRecord.seq)
    return seq_dict

def sort_seqlens(seq_dict):
    '''
    Creates sorted list of tuples of seq id and sequence lengths from dictionary
    '''    
    seqlen_tuple = sorted((len(v),k) for (k,v) in seq_dict.items()) #list comprehension of tuple (k = id,v = length)
    return seqlen_tuple

def len_minmax(seqs_tuple):
    '''
    Returns max and min lengths of list of tuples
    '''
    longest = seqs_tuple[-1][0]
    shortest = seqs_tuple[0][0]
    return ('Shortest: ', shortest, 'Longest: ', longest)

def seq_dups(seq_tuples):
    '''
    Prints sequences with the same lengths from a list of tuples with id and sequence lengths
    '''
    #print(seq_tuples)
    len_dup = []
    for i in seq_tuples:  #checks for duplication of lengths i[0]
        if i[0] not in len_dup:
            len = i[0]
            len_dup.append(len)
        else:
            print('Same length, Id:', i[1], 'length:', i[0])
            continue
        #print(len_dup)

def rev(seq):
    '''
    Reverses sequence
    '''
    seq_rev = seq[::-1]
    return seq_rev

def complement(seq):  
    '''
    Uses dictionary for complementation partners simplified (only lowercase)
    '''
    seq = seq.upper()
    base_comp = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
    complement_seq = ''
    for n in seq:
        for key,value in base_comp.items():
            if n == key:
                complement_seq += value
    return complement_seq

def rev_complement(seq):
    seq = rev(seq)
    seq = complement(seq)
    return seq


def frame_dict(seq_dict):
    '''
    Returns a dictionary of each frame (all 6) of sequences in order: 1_F, 1_R, etc
    '''
    frames_dict = {}
    for key,seq in seq_dict.items():
        #print(key)
        revcomp_seq = rev_complement(seq_dict[key])
        #print(revcomp_seq)
        frames = []
        for j in range(0, 3):
            frame_f = []
            frame_r = []
            for n in range(j,len(seq_dict[key]),3):
                codon = seq_dict[key][n:n+3]
                frame_f.append(codon)
            frames.append(frame_f)
            for n in range(j,len(revcomp_seq),3):
                codon = revcomp_seq[n:n+3]
                #print(codon)
                frame_r.append(codon)
            frames.append(frame_r)
        frames_dict[key] = frames
        continue              
    return frames_dict

def find_start(l):
    '''
    Function returns index if a codon is id'd as start position in a list and returns a list of start indeces
    '''
    start = 'ATG'
    start_list = []
    for index,i in enumerate(l):  #order matters, index comes first then value of each item in list
        if i == start:  #if codon is start codon
            start_pos = index  #get starting index 
            start_list.append(start_pos)
        else:
            continue
    return start_list
    

def find_stop(l):
    '''
    Function finds each stop position in a list of string seq and returns a list of stop indeces
    '''
    stop = ['TAG', 'TGA', 'TAA']
    stop_list = []
    for index,i in enumerate(l):
        for codon in stop:
            if i == codon:  #for each stop codon in list
                stop_pos = index + 1  #get stop index (+1 to include the stop codon)
                #print('stop found', i)
                stop_list.append(stop_pos)
                #print(stop_list)
            else:
                continue
    return stop_list
            

def orfs(frames_dict):
    '''
    This function finds ORFs in each frame and adds them to a dictionary: key = frame, value = orf list
    '''
    for key, value in frames_dict.items():
        seqs = {}  #create dictionary for nesting orf dictionary
        #print(key,value)
        frames = value  #assigns all frames for a sequence to variable frames
        orf_dict = {}
        
        for index,frame in enumerate(frames):  #for each of the 6 frames, create frame label
            
            if index == 0:
                f = '1_F'
            if index == 1:
                f = '1_R'
            if index == 2:
                f = '2_F'
            if index == 3:
                f = '2_R'
            if index == 4:
                f = '3_F'
            if index == 5:
                f = '3_R'
              
            start = find_start(frame)
            stop = find_stop(frame)
            #print(start,stop)
            orf_dict[f] = []
            
            for i in start:
                orfs = []
                i = int(i)
                for j in stop:
                    j = int(j)
                    seq = ''.join(frame[i:j])
                    orfs.append(seq)
                #print(orfs)
                orf_dict[f] = orfs
        seqs[key] = orf_dict
        return seqs
       
    
def longest_orf(orf_dict):
    '''
    prints the longest orf of each sequence
    '''
    lens = []
    for key,value in orf_dict.items():
        for orf in value:
            lens.append(len(orf))
            sorted_lens = sorted(lens)
            #print(sorted_lens)
        longest_orf = sorted_lens[-1]
        return 'Longest orf for', key, 'is', longest_orf

def repeats(seqs_dict,n):
    for key,seq in seqs_dict.items():
        #print(seqs_dict[key])
        substrings = []
        repeat_dict = {}
        seq_list_by_n =[]
        
        for i in range(0,len(seq),n):  #default frame 1
            subseq = seq[i:i+n].upper()
            #print(subseq)
            if subseq not in substrings:
                substrings.append(subseq)
            #print(substrings)
            for s in substrings:
                repeat_dict[s] = 0
            #print(repeat_dict)
        
        for i in range(0,len(seq),n):  #default frame 1
            subseq = seq[i:i+n].upper()
            seq_list_by_n.append(subseq)
        
        #print(seq_list_by_n)
        for s in seq_list_by_n:
            repeat_dict[s] += 1
        return repeat_dict            
            
