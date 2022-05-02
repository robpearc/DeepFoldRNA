# Author: Robin Pearce
# To report bugs and questions email: robpearc@umich.edu 

import torch
from torch import nn
from torch.nn import functional as F
import math


class Triangle_Att_Start(nn.Module):
    def __init__(self,global_config,tri_att_start_config):
        super(Triangle_Att_Start,self).__init__()
        self.pair_dim = global_config["pair_dim"]
        self.no_heads = tri_att_start_config["no_heads"]
        self.hidden_dim = tri_att_start_config["hidden_dim"]
        
        self.z_norm = nn.LayerNorm(self.pair_dim)
        self.q_linear = nn.Linear(self.pair_dim,self.hidden_dim*self.no_heads,bias=False)
        self.k_linear = nn.Linear(self.pair_dim,self.hidden_dim*self.no_heads,bias=False)
        self.v_linear = nn.Linear(self.pair_dim,self.hidden_dim*self.no_heads,bias=False)
        self.bias_linear = nn.Linear(self.pair_dim,self.no_heads,bias=False)
        self.gate_linear = nn.Linear(self.pair_dim,self.hidden_dim*self.no_heads)
        self.output_linear = nn.Linear(self.hidden_dim*self.no_heads,self.pair_dim)
        
        self.norm_factor = 1/math.sqrt(self.hidden_dim)
        self.sigmoid = nn.Sigmoid()

    def forward(self,pair_embed):
        length, _, embed_dim = pair_embed.shape
        
        pair_embed = self.z_norm(pair_embed)
        
        query = self.q_linear(pair_embed)
        key = self.k_linear(pair_embed)
        value = self.v_linear(pair_embed)

        query = query.contiguous().view(length,length,self.no_heads,self.hidden_dim)
        key = key.contiguous().view(length,length,self.no_heads,self.hidden_dim)
        value = value.contiguous().view(length,length,self.no_heads,self.hidden_dim)

        query = query * self.norm_factor

        pair_bias = self.bias_linear(pair_embed)
        pair_bias = pair_bias[None,:,:,:].permute(3,0,1,2)
        attn_map = torch.einsum("ijhd,ikhd->hijk", query, key) + pair_bias
        attn_map = F.softmax(attn_map,dim=-1)

        gate = self.gate_linear(pair_embed)
        gate = self.sigmoid(gate)
        gate = gate.contiguous().view(length,length,self.no_heads,self.hidden_dim)

        out = torch.einsum('hijk,ikhd->ijhd',attn_map, value)
        out = out * gate
        out = out.contiguous().view(length,length,-1)
        pair_embed_out = self.output_linear(out)

        return pair_embed_out

class Triangle_Att_End(nn.Module):
    def __init__(self,global_config,tri_att_end_config):
        super(Triangle_Att_End,self).__init__()
        self.pair_dim = global_config["pair_dim"]
        self.no_heads = tri_att_end_config["no_heads"]
        self.hidden_dim = tri_att_end_config["hidden_dim"]
        
        self.z_norm = nn.LayerNorm(self.pair_dim)
        self.q_linear = nn.Linear(self.pair_dim,self.hidden_dim*self.no_heads,bias=False)
        self.k_linear = nn.Linear(self.pair_dim,self.hidden_dim*self.no_heads,bias=False)
        self.v_linear = nn.Linear(self.pair_dim,self.hidden_dim*self.no_heads,bias=False)
        self.bias_linear = nn.Linear(self.pair_dim,self.no_heads,bias=False)
        self.gate_linear = nn.Linear(self.pair_dim,self.hidden_dim*self.no_heads)
        self.output_linear = nn.Linear(self.hidden_dim*self.no_heads,self.pair_dim)
       
        self.norm_factor = 1/math.sqrt(self.hidden_dim)
        self.sigmoid = nn.Sigmoid()
    
    def forward(self,pair_embed):
        length, _, embed_dim = pair_embed.shape
        
        pair_embed = self.z_norm(pair_embed)
        
        query = self.q_linear(pair_embed)
        key = self.k_linear(pair_embed)
        value = self.v_linear(pair_embed)

        query = query.contiguous().view(length,length,self.no_heads,self.hidden_dim)
        key = key.contiguous().view(length,length,self.no_heads,self.hidden_dim)
        value = value.contiguous().view(length,length,self.no_heads,self.hidden_dim)

        query = query * self.norm_factor

        pair_bias = self.bias_linear(pair_embed)
        pair_bias = pair_bias[None,:,:,:].permute(3,0,2,1)
        attn_map = torch.einsum("ijhd,kihd->hijk", query, key) + pair_bias
        attn_map = F.softmax(attn_map,dim=-1)

        gate = self.gate_linear(pair_embed)
        gate = self.sigmoid(gate)
        gate = gate.contiguous().view(length,length,self.no_heads,self.hidden_dim)

        out = torch.einsum('hijk,kjhd->ijhd',attn_map, value)
        out = out * gate
        out = out.contiguous().view(length,length,-1)
        pair_embed_out = self.output_linear(out)

        return pair_embed_out

class Triangle_Update_Outgoing(nn.Module):
    def __init__(self,global_config,tri_out_config):
        super(Triangle_Update_Outgoing,self).__init__()
        self.pair_dim = global_config["pair_dim"]
        self.hidden_dim = tri_out_config["hidden_dim"]

        self.z_norm = nn.LayerNorm(self.pair_dim)
        self.layer_norm = nn.LayerNorm(self.hidden_dim)
        self.a_linear = nn.Linear(self.pair_dim,self.hidden_dim)
        self.gate_a_linear = nn.Linear(self.pair_dim,self.hidden_dim)
        self.b_linear = nn.Linear(self.pair_dim,self.hidden_dim)
        self.gate_b_linear = nn.Linear(self.pair_dim,self.hidden_dim)
        self.output_linear = nn.Linear(self.hidden_dim,self.pair_dim)
        self.gate_linear = nn.Linear(self.pair_dim,self.pair_dim)

        self.sigmoid = nn.Sigmoid()

    def forward(self,pair_embed):
        pair_embed = self.z_norm(pair_embed)
        
        node_1 = self.a_linear(pair_embed) 
        gate_node_1 = self.gate_a_linear(pair_embed)
        gate_node_1 = self.sigmoid(gate_node_1)
        node_1 = node_1 * gate_node_1

        node_2 = self.b_linear(pair_embed)
        gate_node_2 = self.gate_b_linear(pair_embed)
        gate_node_2 = self.sigmoid(gate_node_2)
        node_2 = node_2 * gate_node_2

        pair_update = torch.einsum('ikh,jkh->ijh',node_1,node_2)
        pair_update = self.layer_norm(pair_update)
        pair_update = self.output_linear(pair_update)
        gate_out = self.gate_linear(pair_embed)
        gate_out = self.sigmoid(gate_out)
        pair_update = pair_update * gate_out
        
        return pair_update

class Triangle_Update_Incoming(nn.Module):
    def __init__(self,global_config,tri_in_config):
        super(Triangle_Update_Incoming,self).__init__()
        self.pair_dim = global_config["pair_dim"]
        self.hidden_dim = tri_in_config["hidden_dim"]
        
        self.layer_norm = nn.LayerNorm(self.hidden_dim)
        self.z_norm = nn.LayerNorm(self.pair_dim)
        self.a_linear = nn.Linear(self.pair_dim,self.hidden_dim)
        self.gate_a_linear = nn.Linear(self.pair_dim,self.hidden_dim)
        self.b_linear = nn.Linear(self.pair_dim,self.hidden_dim)
        self.gate_b_linear = nn.Linear(self.pair_dim,self.hidden_dim)
        self.output_linear = nn.Linear(self.hidden_dim,self.pair_dim)
        self.gate_linear = nn.Linear(self.pair_dim,self.pair_dim)

        self.sigmoid = nn.Sigmoid()

    def forward(self,pair_embed):
        pair_embed = self.z_norm(pair_embed)

        node_1 = self.a_linear(pair_embed) 
        gate_node_1 = self.gate_a_linear(pair_embed)
        gate_node_1 = self.sigmoid(gate_node_1)
        node_1 = node_1 * gate_node_1

        node_2 = self.b_linear(pair_embed)
        gate_node_2 = self.gate_b_linear(pair_embed)
        gate_node_2 = self.sigmoid(gate_node_2)
        node_2 = node_2 * gate_node_2

        pair_update = torch.einsum('ijh,ikh->jkh',node_1,node_2)
        pair_update = self.layer_norm(pair_update)
        pair_update = self.output_linear(pair_update)
        gate_out = self.gate_linear(pair_embed)
        gate_out = self.sigmoid(gate_out)
        pair_update = pair_update * gate_out
        
        return pair_update

class Pair_Transition(nn.Module):
    def __init__(self,global_config,pair_trans_config):
        super(Pair_Transition,self).__init__()
        self.pair_dim = global_config["pair_dim"]
        self.hidden_dim = pair_trans_config["hidden_dim"]
        
        self.z_norm = nn.LayerNorm(self.pair_dim)
        self.linear_trans_1 = nn.Linear(self.pair_dim,self.pair_dim*self.hidden_dim)
        self.linear_trans_2 = nn.Linear(self.pair_dim*self.hidden_dim,self.pair_dim)
        
        self.relu = nn.ReLU()

    def forward(self,pair_embed):
        pair_embed = self.z_norm(pair_embed)
        pair_embed = self.linear_trans_1(pair_embed)
        pair_embed = self.relu(pair_embed)
        pair_embed = self.linear_trans_2(pair_embed)
        
        return pair_embed
