# Author: Robin Pearce
# To report bugs and questions email: robpearc@umich.edu 

import torch
from torch import nn
from torch.nn import functional as F
import math


class MSA_Row_Att(nn.Module):
    def __init__(self,global_config,msa_row_att_config):
        super(MSA_Row_Att,self).__init__()
        self.msa_dim = global_config["msa_dim"]
        self.pair_dim = global_config["pair_dim"]
        self.no_heads = msa_row_att_config["no_heads"]
        self.hidden_dim = msa_row_att_config["hidden_dim"]
        
        self.m_norm = nn.LayerNorm(self.msa_dim)
        self.z_norm = nn.LayerNorm(self.pair_dim)
        self.bias_linear = nn.Linear(self.pair_dim,self.no_heads,bias=False)
        self.q_linear = nn.Linear(self.msa_dim,self.no_heads*self.hidden_dim,bias=False)
        self.k_linear = nn.Linear(self.msa_dim,self.no_heads*self.hidden_dim,bias=False)
        self.v_linear = nn.Linear(self.msa_dim,self.no_heads*self.hidden_dim,bias=False)
        self.gate_linear = nn.Linear(self.msa_dim,self.no_heads*self.hidden_dim)
        self.output_linear = nn.Linear(self.no_heads*self.hidden_dim,self.msa_dim)
        
        self.norm_factor = 1/math.sqrt(self.hidden_dim)
        self.sigmoid = nn.Sigmoid()

    def forward(self,msa_embed,pair_embed):
        num_seqs,length,embed_dim = msa_embed.shape
        
        msa_embed = self.m_norm(msa_embed)
        pair_embed = self.z_norm(pair_embed)

        query = self.q_linear(msa_embed) 
        key = self.k_linear(msa_embed)
        value = self.v_linear(msa_embed)
        
        query = query.contiguous().view(num_seqs,length,self.no_heads,self.hidden_dim)
        key = key.contiguous().view(num_seqs,length,self.no_heads,self.hidden_dim)
        value = value.contiguous().view(num_seqs,length,self.no_heads,self.hidden_dim)

        query = query * self.norm_factor

        pair_bias = self.bias_linear(pair_embed)
        pair_bias = pair_bias[None,:,:,:].permute(0,3,1,2)
        att_map = torch.einsum('rihd,rjhd->rhij',query,key) + pair_bias
        att_map = F.softmax(att_map,dim=-1)
        
        gate = self.gate_linear(msa_embed)
        gate = self.sigmoid(gate)
        gate = gate.contiguous().view(num_seqs,length,self.no_heads,self.hidden_dim)
        
        out = torch.einsum('rhij,rjhd->rihd',att_map,value)
        out = out * gate
        out = out.contiguous().view(num_seqs,length,-1)
        msa_embed_out = self.output_linear(out)
        
        return msa_embed_out

class MSA_Col_Att(nn.Module):
    def __init__(self,global_config,msa_col_att_config):
        super(MSA_Col_Att,self).__init__()
        self.msa_dim = global_config["msa_dim"]
        self.pair_dim = global_config["pair_dim"]
        self.no_heads = msa_col_att_config["no_heads"]
        self.hidden_dim = msa_col_att_config["hidden_dim"]

        self.m_norm = nn.LayerNorm(self.msa_dim)
        self.q_linear = nn.Linear(self.msa_dim,self.no_heads*self.hidden_dim,bias=False)
        self.k_linear = nn.Linear(self.msa_dim,self.no_heads*self.hidden_dim,bias=False)
        self.v_linear = nn.Linear(self.msa_dim,self.no_heads*self.hidden_dim,bias=False)
        self.gate_linear = nn.Linear(self.msa_dim,self.no_heads*self.hidden_dim)
        self.output_linear = nn.Linear(self.no_heads*self.hidden_dim,self.msa_dim)
        
        self.norm_factor = 1/math.sqrt(self.hidden_dim)
        self.sigmoid = nn.Sigmoid()

    def forward(self,msa_embed):
        num_seqs, length, embed_dim = msa_embed.shape
        
        msa_embed = self.m_norm(msa_embed)
        
        query = self.q_linear(msa_embed)
        key = self.k_linear(msa_embed)
        value = self.v_linear(msa_embed)

        query = query.contiguous().view(num_seqs,length,self.no_heads,self.hidden_dim)
        key = key.contiguous().view(num_seqs,length,self.no_heads,self.hidden_dim)
        value = value.contiguous().view(num_seqs,length,self.no_heads,self.hidden_dim)
        
        query = query * self.norm_factor
        
        attn_map = torch.einsum("ichd,jchd->hcij", query, key)
        attn_map = F.softmax(attn_map,dim=-1)

        gate = self.gate_linear(msa_embed) 
        gate = self.sigmoid(gate)
        gate = gate.contiguous().view(num_seqs,length,self.no_heads,self.hidden_dim)
        
        out = torch.einsum('hcij,jchd->ichd',attn_map, value) 
        out = out * gate
        out = out.contiguous().view(num_seqs,length,-1)
        m_embed_out = self.output_linear(out)
        
        return m_embed_out

class MSA_Transition(nn.Module):
    def __init__(self,global_config,msa_transition_config):
        super(MSA_Transition,self).__init__()
        self.msa_dim = global_config["msa_dim"]
        self.hidden_dim = msa_transition_config["hidden_dim"]
        
        self.m_norm = nn.LayerNorm(self.msa_dim)
        self.linear_trans_1 = nn.Linear(self.msa_dim,self.msa_dim*self.hidden_dim)
        self.linear_trans_2 = nn.Linear(self.msa_dim*self.hidden_dim,self.msa_dim)

        self.relu = nn.ReLU()

    def forward(self,msa_embed):
        msa_embed = self.m_norm(msa_embed)
        msa_embed = self.linear_trans_1(msa_embed)
        msa_embed = self.relu(msa_embed)
        msa_embed = self.linear_trans_2(msa_embed)
        
        return msa_embed

class MSA_Outer_Product_Mean(nn.Module):
    def __init__(self,global_config,msa_outer_product_mean_config):
        super(MSA_Outer_Product_Mean,self).__init__()
        self.msa_embed_dim = global_config["msa_dim"]
        self.pair_dim = global_config["pair_dim"]
        self.hidden_dim = msa_outer_product_mean_config["hidden_dim"]
        self.m_norm = nn.LayerNorm(self.msa_embed_dim)
        self.linear_opm_1 = nn.Linear(self.msa_embed_dim,self.hidden_dim)
        self.linear_opm_2 = nn.Linear(self.msa_embed_dim,self.hidden_dim)
        self.linear_opm_3 = nn.Linear(self.hidden_dim*self.hidden_dim,self.pair_dim)

    def forward(self,msa_embed):
        num_seqs, length, embed_dim = msa_embed.shape
        
        msa_embed = self.m_norm(msa_embed)

        msa_1 = self.linear_opm_1(msa_embed)
        msa_2 = self.linear_opm_2(msa_embed)
        pair_opm = torch.einsum('ria,rjb->rijab',msa_2,msa_1)
        pair_opm = torch.mean(pair_opm, dim=0)
        pair_opm = pair_opm.contiguous().view(length,length,-1)
        pair_opm = self.linear_opm_3(pair_opm)
        
        return pair_opm

