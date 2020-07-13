#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#dnautil contains useful functions for dna sequence manipulation/stats


# In[ ]:


# Lecture 5.1 - functions

#computing GC content

def gc(seq):   # function counts gc content of dna string
    n_count = float(seq.count("n") + seq.count("N"))
    gc_total = float(seq.count("c") + seq.count("g"))
    gc_percent = (gc_total / (len(seq) - n_count)) * 100
    return gc_percent


# In[ ]:


#Lecture 5.2 - functions

#checking for stop codons
#addition: would like to identify which frame the stop codon is in

def stop_codon(seq):  # this function checks for codons in every frame
    frame = 1 
    stop_codon_found = False
    stop_codons = ["tga","tag","taa"]
    print(len(seq))
    for i in range(0,len(seq),3):  #default frame 1
        codon = seq[i:i+3].lower()
        if codon in stop_codons:
            stop_codon_found = True
            break
    for i in range(1,len(seq),3):  #default frame 0
        codon = seq[i:i+3].lower()
        frame = 2
        if codon in stop_codons:
            stop_codon_found = True
            break
    for i in range(2,len(seq),3):  #default frame 0
        codon = seq[i:i+3].lower()
        frame = 3
        if codon in stop_codons:
            stop_codon_found = True
            break
    return stop_codon_found, frame


# In[ ]:


#Lecture 5.3 - functions
#reverse complement; reversing a string, combining functions 
#my code first

def rev(dna):
    seq_rev = dna[::-1]
    return seq_rev

def complement(dna):  #use dictionary for complementation partners simplified (only lowercase)
    base_comp = {'a':'t','t':'a','c':'g','g':'c','n':'n'}
    complement_seq = ''
    for n in dna:
        for key,value in base_comp.items():
            if n == key:
                complement_seq += value
    return complement_seq

def rev_complement(dna):
    seq = rev(dna)
    seq = complement(seq)
    return seq
    


# In[4]:


#profs code EXAMPLE OF LIST COMPREHENSION

def rev_complement_lc(dna):
    letters = list(dna)
    base_comp = {'a':'t','t':'a','c':'g','g':'c','n':'n'}
    c_letters = [base_comp[base] for base in letters]  #list comp
    rev_c = ''
    return rev_c.join(c_letters[::-1])

seq = 'actagatag'
rev_complement(seq)

