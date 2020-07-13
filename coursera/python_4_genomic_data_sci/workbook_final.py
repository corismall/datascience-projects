#!/usr/bin/env python
# coding: utf-8

# In[1]:


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


# In[2]:


def dict_seqs(seqs):
    '''
    Indexing fasta sequences into seqrecords and writing id and sequence into a dictionary
    '''
    seq_dict = {}
    for key, SeqRecord in seqs.items():
        seq_dict[key] = str(SeqRecord.seq)
    return seq_dict


# In[124]:


def sort_seqlens(seq_dict):
    '''
    Creates sorted list of tuples of seq id and sequence lengths from dictionary
    '''    
    seqlen_tuple = sorted((len(v),k) for (k,v) in seq_dict.items()) #list comprehension of tuple (k = id,v = length)
    return seqlen_tuple


# In[149]:


def len_minmax(seqs_tuple):
    '''
    Returns max and min lengths of list of tuples
    '''
    longest = seqs_tuple[-1][0]
    shortest = seqs_tuple[0][0]
    return ('Shortest: ', shortest, 'Longest: ', longest)


# In[ ]:


def seq_dups(seq_tuples):
    '''
    Prints sequences with the same lengths from a list of tuples with id and sequence lengths
    '''
    #print(seq_tuples)
    len_dup = []
    for i in seq_tuples:  #checks for duplication of lengths i[0]
        if i[0] not in len_dup:
            length = i[0]
            len_dup.append(length)
            print('None')
        else:
            print('Same length, Id:', i[1], 'length:', i[0])
            continue
        #print(len_dup)


# In[ ]:


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


# In[10]:


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


# In[5]:


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


# In[6]:


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
            


# In[7]:


def orfs(frames_dict):
    '''
    Input: Dictionary of a list of codons by frame
    Output: ORFs in each frame, added to a dictionary by frame; nests frame dictionary into a dictionary by sequence
    '''
    seqs = {}  #create dictionary to hold orfs of all frames per sequence
    for key, value in frames_dict.items():
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
              
            start = find_start(frame)  #find start indeces
            stop = find_stop(frame)  #find stop indeces
            #print(f,'start',start,'stop',stop)
            orfs = []
            orf_lens = []
            for i in start:
                i = int(i)
                #print(i)
                for j in stop:
                    if j > i:
                        j = int(j)
                        #print(j)
                        seq = frame[i:j]
                        #print(f,seq,len(seq),'\n')
                        orf_tup = (key,f,len(seq))
                        orfs.append(''.join(seq))
                        orf_lens.append(orf_tup)
            
                #print(orfs)
                orf_dict[f] = orfs
                #print(orf_dict)
                s_orf_lengths = sorted(orf_lens, key = lambda tup : tup[2])
                print(s_orf_lengths)
                
        
            seqs[key] = orf_dict
            #print(seqs)
    return seqs


# In[8]:


def sorted_orfs(nested_dict):
    '''
    Input: nested dictionary of orfs by frame by identifier
    Output: A dictionary of sorted lists of orf lengths by identifier 
    '''
    orf_len = {}
    for key,value in nested_dict.items():
        #print(key,value)
        d = value
        #print(d)
        lengths = []
        for i,j in d.items():
            #print(i,j)
            for k in j:
                l = len(k)
                t = (i,l)
                lengths.append(t)
            #print(lengths)
            s = sorted(lengths, key = lambda tup : tup[1])
            #print(s)
        orf_len[key] = s
    return orf_len


# In[ ]:


def orf_len_by_frame(dictionary):
    l = sorted_orfs(dictionary)
    frame_list = []
    f_r = input('What frame do you want?')
    for k,v in l.items():
        #print(k,v)
        for t in v:
            if t[0] == f_r:
                frame_list.append(t[1])
            else:
                pass

    return frame_list




