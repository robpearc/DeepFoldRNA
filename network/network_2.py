# Author: Robin Pearce
# To report bugs and questions email: robpearc@umich.edu 

import torch
from torch import nn
import numpy as np
from torch.nn import functional as F

from network.config import network_config 
from network import embedders, msa_transformer, single_transformer, pair_transformer, heads


class Pred_Network(nn.Module):
    def __init__(self,device):
        super(Pred_Network,self).__init__()
        self.msa_dim = network_config["msa_dim"]
        self.seq_dim = network_config["seq_dim"]
        self.no_cycles = network_config["no_cycles"]
        
        self.input_embedder = embedders.Input_Embedder(network_config,network_config["input_embedder"],device) 
        self.recycling_embedder = embedders.Recycling_Embedder(network_config) 
        self.recycling_embedder_s = embedders.Recycling_Embedder_S(network_config)
        self.ss_embedder = embedders.Secondary_Structure_Embedder(network_config)
        self.hmm_embedder = embedders.HMM_Embedder(network_config)
        self.msa_transformer_layers = msa_transformer.MSA_Transformer(network_config)
        self.pair_transformer_layers = pair_transformer.Pair_Transformer(network_config)
        self.single_transformer = single_transformer.Single_Transformer(network_config)
        self.linear_s = nn.Linear(self.msa_dim,self.seq_dim)
        self.geometry_head = heads.Geometry_Head(network_config,network_config["geometry_head"])
        self.tor_head = heads.Torsion_Head(network_config,network_config["torsion_head"])
        self.msa_head = heads.MSA_Head(network_config)

    def iteration(self,features,msa_recycle,pair_recycle,seq_recycle,recycle_flag=False):
        num_seqs, length, embed_dim = features["msa"].shape
        
        seq = features["seq"]
        msa = features["msa"]
        ss = features["ss"]
        hmm = features["hmm"] 
        
        msa_mask = torch.zeros(num_seqs,length).to(msa.device)
        masked_msa = torch.cat([msa*(1-msa_mask[:,:,None]),msa_mask[:,:,None]],dim=-1)
        msa_embed, pair_embed = self.input_embedder(seq,masked_msa)
        hmm_embed = self.hmm_embedder(hmm)
        msa_embed = torch.cat([msa_embed,hmm_embed[None,:,:]],dim=0)
        
        pair_embed_ss = self.ss_embedder(ss) 
        pair_embed_ss = self.pair_transformer_layers(pair_embed_ss)
        pair_embed = pair_embed + pair_embed_ss

        if(recycle_flag):
            msa_recycle, pair_recycle = self.recycling_embedder(msa_recycle,pair_recycle) 
            msa_embed = msa_embed + msa_recycle
            pair_embed = pair_embed + pair_recycle

        msa_embed, pair_embed = self.msa_transformer_layers(msa_embed,pair_embed)
       
        sequence_embedding = msa_embed[0]
        seq_embed = self.linear_s(sequence_embedding)
        if(recycle_flag):
            seq_recycle = self.recycling_embedder_s(seq_recycle)
            seq_embed = seq_embed + seq_recycle 
        
        seq_embed = self.single_transformer(seq_embed,pair_embed)
        
        return msa_embed, pair_embed, seq_embed

    def forward(self,features):
        cycles=self.no_cycles

        msa_recycle = None
        pair_recycle = None
        seq_recycle = None
        recycle = False
        with torch.no_grad():
            for i in range(cycles-1):
                msa_recycle, pair_recycle, seq_recycle = self.iteration(features,msa_recycle,pair_recycle,seq_recycle,recycle)
                recycle = True
        
        msa_final, pair_final, seq_final = self.iteration(features,msa_recycle,pair_recycle,seq_recycle,recycle)
        
        output = {}
        output["geometry_head"] = self.geometry_head(pair_final)
        output["torsion_head"] = self.tor_head(seq_final)
        
        return output

