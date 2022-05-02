# Author: Robin Pearce
# To report bugs and questions email: robpearc@umich.edu

import os
import torch
from torch.nn import functional as F
import numpy as np
import random

nucleotides = {
    'A': 1,
    'G': 2,
    'C': 3,
    'T': 4,
}

def read_hmm(hmm_file,L):
    count = 0
    seq_idx = 0
    hmm = np.zeros([L,5,3])
    f = open(hmm_file,'r')
    lines = f.readlines()
    for line in lines[int(-L*3-1):-1]:
        line = line.strip()
        line = line.split()
        if(count==0):
            try:
                hmm[seq_idx,0,0] = float(line[1])
            except:
                pass
            try:
                hmm[seq_idx,1,0] = float(line[2])
            except:
                pass
            try:
                hmm[seq_idx,2,0] = float(line[3])
            except:
                pass
            try:
                hmm[seq_idx,3,0] = float(line[4])
            except:
                pass
        elif(count==1):
            try:
                hmm[seq_idx,0,1] = float(line[0])
            except:
                pass
            try:
                hmm[seq_idx,1,1] = float(line[1])
            except:
                pass
            try:
                hmm[seq_idx,2,1] = float(line[2])
            except:
                pass
            try:
                hmm[seq_idx,3,1] = float(line[3])
            except:
                pass
        elif(count==2):
            try:
                hmm[seq_idx,0,2] = float(line[0])
            except:
                pass
            try:
                hmm[seq_idx,1,2] = float(line[1])
            except:
                pass
            try:
                hmm[seq_idx,2,2] = float(line[2])
            except:
                pass
            try:
                hmm[seq_idx,3,2] = float(line[3])
            except:
                pass
            try:
                hmm[seq_idx,4,0] = float(line[4])
            except:
                pass
            try:
                hmm[seq_idx,4,1] = float(line[5])
            except:
                pass
            try:
                hmm[seq_idx,4,2] = float(line[6])
            except:
                pass
        count+=1
        if(count>2):
            count = 0
            seq_idx+=1
    hmm = torch.from_numpy(hmm.reshape(L,15)).float()
    
    return hmm

def read_msa(msa_file):
    f = open(msa_file,'r')
    lines = f.read()
    f.close()

    lines = lines.split('\n')
    lines = lines[1::2]

    msa_ = np.array([list(s.strip()) for s in lines])
    msa = np.zeros_like(msa_,dtype=int)
    for akey in list(nucleotides.keys()):
        msa[msa_==akey]=nucleotides[akey]
    
    return msa

def read_seq(seq_file):
    sequence_lines=open(seq_file).readlines()
    sequence = sequence_lines[1].strip()
    L = len(sequence)
    seq_array = np.zeros(L)
    sequence_list = np.array(list(sequence))
    seq_array[sequence_list=='A']=1
    seq_array[sequence_list=='G']=2
    seq_array[sequence_list=='C']=3
    seq_array[sequence_list=='U']=4
    
    return seq_array

def randomly_sample_msa(unsampled_msa):
    msa_rand = unsampled_msa[1:,:]
    num_seqs,length = msa_rand.shape
    if num_seqs>=1:
        np.random.shuffle(msa_rand)
        num_sel = min(int(num_seqs), 128)
        idx = np.random.randint(num_seqs,size=num_sel)
        msa = msa_rand[idx,:]
    else:
        msa = np.int64(unsampled_msa)
    
    return msa

def parse_ss(ss_file):
    ss = np.load(ss_file)
    
    return ss

def collect_features(seq_file,msa_file,hmm_file,ss_file):
    seq = read_seq(seq_file)
    msa = read_msa(msa_file)
    msa = randomly_sample_msa(msa)
    msa[0] = seq[None,:]
    msa = torch.from_numpy(msa).long()
    num_seqs, length = msa.shape
    if num_seqs<1.9:
        msa = torch.cat([msa,msa],0)
    msa = F.one_hot(msa,6)

    hmm = None
    if(os.path.exists(hmm_file)):
        hmm = read_hmm(hmm_file,length)
    else:
        hmm = torch.zeros([length,15])
    
    ss = None
    if(os.path.exists(ss_file)):
        ss = np.loadtxt(ss_file,skiprows=1)[:-1,:]
    else:
        ss = np.zeros([length,length])
    #ss = ss + 0.0
    ss = torch.from_numpy(ss).float()[...,None]
    
    features = {}
    features["ss"] = ss
    features["msa"] = msa.to(torch.int64)
    features["seq"] = msa[0].float()
    features["hmm"] = hmm.float()
    
    return features

