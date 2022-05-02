# Author: Robin Pearce
# To report bugs and questions email: robpearc@umich.edu

import torch
from torch import nn
from torch.nn import functional as F


class Input_Embedder(nn.Module):
    def __init__(self,config,input_config,device):
        super(Input_Embedder,self).__init__()
        self.input_msa_dim = config["input_msa_dim"]
        self.input_seq_dim = config["input_seq_dim"]
        self.msa_dim = config["msa_dim"]
        self.pair_dim = config["pair_dim"]
        self.no_pos_bins_1 = input_config["no_pos_bins_1"]
        self.no_pos_bins_2 = input_config["no_pos_bins_2"]
        self.rel_pos_1d = int(input_config["rel_pos_1d"])
        self.min_rel_range = int(input_config["min_rel_range"])
        self.max_rel_range = int(input_config["max_rel_range"])
        self.max_seq_len = int(input_config["max_seq_length"])

        self.linear_msa_embed = nn.Linear(self.input_msa_dim,self.msa_dim)
        self.linear_seq_pair_i = nn.Linear(self.input_seq_dim,self.pair_dim)
        self.linear_seq_pair_j = nn.Linear(self.input_seq_dim,self.pair_dim)
        self.linear_seq_m = nn.Linear(self.input_seq_dim,self.msa_dim)
        self.linear_embed_pos_1 = nn.Linear(self.no_pos_bins_1,self.msa_dim)
        self.linear_embed_pos_2 = nn.Linear(self.no_pos_bins_2,self.pair_dim)
        self.pos_1 = self.compute_position_1(device)
        self.pos_2 = self.compute_position_2(device)

    def compute_position_1(self,device):
        pos = torch.arange(self.max_seq_len)
        rel_pos = ((pos[:,None] & (1 << torch.arange(self.rel_pos_1d)))) > 0
        
        return rel_pos.float().to(device)

    def compute_position_2(self,device):
        pos = torch.arange(self.max_seq_len)
        rel_pos = pos[None,:]-pos[:,None]
        rel_pos = rel_pos.clamp(self.min_rel_range,self.max_rel_range)
        rel_pos_encode = F.one_hot(rel_pos+self.max_rel_range,self.no_pos_bins_2)
        
        return rel_pos_encode.float().to(device)

    def forward(self,input_seq,input_msa):
        num_seqs, length, dim = input_msa.shape

        embedded_seq_m = self.linear_seq_m(input_seq)
        embedded_msa = self.linear_msa_embed(input_msa)
        rel_pos_1d = self.pos_1[:length]
        rel_pos_2d = self.pos_2[:length,:length]
        embedded_pos_1 = self.linear_embed_pos_1(rel_pos_1d)
        embedded_pos_2 = self.linear_embed_pos_2(rel_pos_2d)
        embedded_seq_pair_i = self.linear_seq_pair_i(input_seq)
        embedded_seq_pair_j = self.linear_seq_pair_j(input_seq)

        msa_embed = embedded_msa + embedded_seq_m[None,:,:] + embedded_pos_1[None,:,:]
        pair_embed = embedded_seq_pair_i[None,:,:] + embedded_seq_pair_j[:,None,:] + embedded_pos_2
        
        return msa_embed, pair_embed

class HMM_Embedder(nn.Module):
    def __init__(self,config):
        super(HMM_Embedder,self).__init__()
        self.msa_dim = config["msa_dim"]
        self.hmm_linear = nn.Linear(15,self.msa_dim)

    def forward(self,hmm):
        hmm_embed = self.hmm_linear(hmm)
        
        return hmm_embed

class Secondary_Structure_Embedder(nn.Module):
    def __init__(self,config):
        super(Secondary_Structure_Embedder,self).__init__()
        self.pair_dim = config["pair_dim"]
        self.ss_linear = nn.Linear(1,self.pair_dim)

    def forward(self,ss):
        pair_embed_ss = self.ss_linear(ss)
        
        return pair_embed_ss

class Recycling_Embedder_S(nn.Module):
    def __init__(self,config):
        super(Recycling_Embedder_S,self).__init__()
        self.seq_dim = config["seq_dim"]
        self.s_norm = nn.LayerNorm(self.seq_dim)
    
    def forward(self,seq_embed):
        seq_embed = self.s_norm(seq_embed)
            
        return seq_embed

class Recycling_Embedder(nn.Module):
    def __init__(self,config):
        super(Recycling_Embedder,self).__init__()
        self.msa_dim = config["msa_dim"]
        self.pair_dim = config["pair_dim"]
        self.m_norm = nn.LayerNorm(self.msa_dim)
        self.z_norm = nn.LayerNorm(self.pair_dim)
    
    def forward(self,msa_embed,pair_embed):
        msa_embed = self.m_norm(msa_embed)
        pair_embed = self.z_norm(pair_embed)
            
        return msa_embed, pair_embed
