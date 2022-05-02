# Author: Robin Pearce
# To report bugs and questions email: robpearc@umich.edu 

import torch
from torch import nn
from torch.nn import functional as F
import math


class S_Att_Pair_Bias(nn.Module):
    def __init__(self,global_config,s_att_pair_bias_config):
        super(S_Att_Pair_Bias,self).__init__()
        self.seq_dim = global_config["seq_dim"]
        self.pair_dim = global_config["pair_dim"]
        self.no_heads = s_att_pair_bias_config["no_heads"]
        self.hidden_dim = s_att_pair_bias_config["hidden_dim"]
        
        self.s_norm = nn.LayerNorm(self.seq_dim)
        self.z_norm = nn.LayerNorm(self.pair_dim)
        self.bias_linear = nn.Linear(self.pair_dim,self.no_heads,bias=False)
        self.q_linear = nn.Linear(self.seq_dim,self.no_heads*self.hidden_dim,bias=False)
        self.k_linear = nn.Linear(self.seq_dim,self.no_heads*self.hidden_dim,bias=False)
        self.v_linear = nn.Linear(self.seq_dim,self.no_heads*self.hidden_dim,bias=False)
        self.gate_linear = nn.Linear(self.seq_dim,self.no_heads*self.hidden_dim)
        self.output_linear = nn.Linear(self.no_heads*self.hidden_dim,self.seq_dim)
        
        self.norm_factor = 1/math.sqrt(self.hidden_dim)
        self.sigmoid = nn.Sigmoid()

    def forward(self,seq_embed,pair_embed):
        length, embed_dim = seq_embed.shape
        
        seq_embed = self.s_norm(seq_embed)
        pair_embed = self.z_norm(pair_embed)
        
        query = self.q_linear(seq_embed)
        key = self.k_linear(seq_embed)
        value = self.v_linear(seq_embed)
       
        query = query.contiguous().view(1,length,self.no_heads,self.hidden_dim)
        key = key.contiguous().view(1,length,self.no_heads,self.hidden_dim)
        value = value.contiguous().view(1,length,self.no_heads,self.hidden_dim)
      
        query = query * self.norm_factor

        pair_bias = self.bias_linear(pair_embed)
        pair_bias = pair_bias[None,:,:,:].permute(0,3,1,2)
        att_map = torch.einsum('rihd,rjhd->rhij',query,key) + pair_bias
        att_map = F.softmax(att_map,dim=-1)

        gate = self.gate_linear(seq_embed)
        gate = self.sigmoid(gate)
        gate = gate.contiguous().view(1,length,self.no_heads,self.hidden_dim)

        out = torch.einsum('rhij,rjhd->rihd',att_map,value)
        out = out * gate
        out = out.contiguous().view(length,-1)
        seq_embed_out = self.output_linear(out)
        
        return seq_embed_out

class S_Att(nn.Module):
    def __init__(self,global_config,s_att_config):
        super(S_Att,self).__init__()
        self.seq_dim=global_config["seq_dim"]
        self.pair_dim=global_config["pair_dim"]
        self.no_heads=s_att_config["no_heads"]
        self.hidden_dim=s_att_config["hidden_dim"]

        self.s_norm = nn.LayerNorm(self.seq_dim)
        self.z_norm = nn.LayerNorm(self.pair_dim)
        self.q_linear = nn.Linear(self.seq_dim,self.no_heads*self.hidden_dim,bias=False)
        self.k_linear = nn.Linear(self.seq_dim,self.no_heads*self.hidden_dim,bias=False)
        self.v_linear = nn.Linear(self.seq_dim,self.no_heads*self.hidden_dim,bias=False)
        self.gate_linear = nn.Linear(self.seq_dim,self.no_heads*self.hidden_dim)
        self.output_linear = nn.Linear(self.no_heads*self.hidden_dim,self.seq_dim)
        
        self.norm_factor = 1/math.sqrt(self.hidden_dim)
        self.sigmoid = nn.Sigmoid()

    def forward(self,seq_embed):
        length, embed_dim = seq_embed.shape
        
        seq_embed = self.s_norm(seq_embed)
        
        query = self.q_linear(seq_embed)
        key = self.k_linear(seq_embed)
        value = self.v_linear(seq_embed)
       
        query = query.contiguous().view(1,length,self.no_heads,self.hidden_dim)
        key = key.contiguous().view(1,length,self.no_heads,self.hidden_dim)
        value = value.contiguous().view(1,length,self.no_heads,self.hidden_dim)
      
        query = query * self.norm_factor

        att_map = torch.einsum('rihd,rjhd->rhij',query,key)
        att_map = F.softmax(att_map,dim=-1)

        gate = self.gate_linear(seq_embed)
        gate = self.sigmoid(gate)
        gate = gate.contiguous().view(1,length,self.no_heads,self.hidden_dim)

        out = torch.einsum('rhij,rjhd->rihd',att_map,value)
        out = out * gate
        out = out.contiguous().view(length,-1)
        seq_embed_out = self.output_linear(out)
        
        return seq_embed_out

class S_Transition(nn.Module):
    def __init__(self, config):
        super(S_Transition, self).__init__()
        self.seq_dim = config["seq_dim"]

        self.s_norm = nn.LayerNorm(self.seq_dim)
        self.linear_1 = nn.Linear(self.seq_dim, self.seq_dim)
        self.linear_2 = nn.Linear(self.seq_dim, self.seq_dim)
        self.linear_3 = nn.Linear(self.seq_dim, self.seq_dim)

        self.relu = nn.ReLU()

    def forward(self, seq_embed):
        seq_embed = self.s_norm(seq_embed)
        seq_embed = self.linear_1(seq_embed)
        seq_embed = self.relu(seq_embed)
        seq_embed = self.linear_2(seq_embed)
        seq_embed = self.relu(seq_embed)
        seq_embed = self.linear_3(seq_embed)

        return seq_embed

