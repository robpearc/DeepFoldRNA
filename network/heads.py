# Author: Robin Pearce
# To report bugs and questions email: robpearc@umich.edu

import torch
from torch import nn
from torch.nn import functional as F


class Angle_Block(nn.Module):
    def __init__(self, hidden_dim):
        super(Angle_Block, self).__init__()
        self.hidden_dim = hidden_dim
        self.linear_1 = nn.Linear(self.hidden_dim,self.hidden_dim)
        self.linear_2 = nn.Linear(self.hidden_dim,self.hidden_dim)

    def forward(self, seq_embed):
        seq_embed = self.linear_1(F.relu(seq_embed))
        seq_embed = self.linear_2(F.relu(seq_embed))
        
        return seq_embed

class Torsion_Head(nn.Module):
    def __init__(self,config,torsion_head_config):
        super(Torsion_Head,self).__init__()
        self.seq_dim = config["seq_dim"]
        self.hidden_dim = torsion_head_config["hidden_dim"]
        self.no_blocks = torsion_head_config["no_blocks"]
        self.no_angles = torsion_head_config["no_angles"]
        self.angle_bins = torsion_head_config["angle_bins"]
        
        self.linear_in = nn.Linear(self.seq_dim,self.hidden_dim)
        self.angle_layer = Angle_Block(self.hidden_dim)
        self.linear_out = nn.Linear(self.hidden_dim,self.no_angles*2)
        self.linear_out_dis = nn.Linear(self.hidden_dim,self.no_angles*self.angle_bins)

    def forward(self,seq_embed):
        length, embed_dim = seq_embed.shape
        seq_embed = self.linear_in(F.relu(seq_embed))

        for i in range(self.no_blocks):
            seq_embed = seq_embed + self.angle_layer(seq_embed)

        seq_embed = F.relu(seq_embed)

        angles_dis = self.linear_out_dis(seq_embed)
        angles_dis = F.log_softmax(angles_dis.contiguous().view(length,self.no_angles,self.angle_bins),dim=-1)

        output = {}
        output["angles_dis"] = angles_dis 

        return output

class Geometry_Head(nn.Module):
    def __init__(self,config,geometry_head_config):
        super(Geometry_Head,self).__init__()
        self.pair_dim = config["pair_dim"]
        self.dis_bins = geometry_head_config["dist_bins"]
        self.omg_bins = geometry_head_config["omega_bins"]
        self.theta_bins = geometry_head_config["theta_bins"]
        self.phi_bins = geometry_head_config["phi_bins"]
        
        self.linear_dis_n = nn.Linear(self.pair_dim,self.dis_bins)
        self.linear_dis_c4 = nn.Linear(self.pair_dim,self.dis_bins)
        self.linear_dis_p = nn.Linear(self.pair_dim,self.dis_bins)
        self.linear_omg = nn.Linear(self.pair_dim,self.omg_bins)
        self.linear_theta = nn.Linear(self.pair_dim,self.theta_bins)
        self.linear_phi = nn.Linear(self.pair_dim,self.phi_bins)

    def forward(self,pair_embed):
        pred_dis_n = self.linear_dis_n(pair_embed)
        pred_dis_n = F.log_softmax(pred_dis_n + pred_dis_n.permute(1,0,2),dim=-1)
        pred_dis_c4 = self.linear_dis_c4(pair_embed)
        pred_dis_c4 = F.log_softmax(pred_dis_c4 + pred_dis_c4.permute(1,0,2),dim=-1)
        pred_dis_p = self.linear_dis_p(pair_embed)
        pred_dis_p = F.log_softmax(pred_dis_p + pred_dis_p.permute(1,0,2),dim=-1)
        pred_omg = self.linear_omg(pair_embed)
        pred_omg = F.log_softmax(pred_omg + pred_omg.permute(1,0,2),dim=-1)
        pred_theta = F.log_softmax(self.linear_theta(pair_embed),dim=-1)
        pred_phi = F.log_softmax(self.linear_phi(pair_embed),dim=-1)

        output = {}
        output["pred_dis_n"] = pred_dis_n
        output["pred_dis_c4"] = pred_dis_c4
        output["pred_dis_p"] = pred_dis_p
        output["pred_omg"] = pred_omg
        output["pred_theta"] = pred_theta
        output["pred_phi"] = pred_phi

        return output

class MSA_Head(nn.Module):
    def __init__(self,config):
        super(MSA_Head,self).__init__()
        self.msa_dim = config["msa_dim"]
        self.in_msa_dim = config["input_msa_dim"]
        self.linear_msa = nn.Linear(self.msa_dim,self.in_msa_dim-1)

    def forward(self,msa_embed):
        pred_msa = F.log_softmax(self.linear_msa(msa_embed[:-1,:,:]),dim=-1)

        output = {}
        output["pred_msa"] = pred_msa

        return output

