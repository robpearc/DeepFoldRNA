# Author: Robin Pearce
# To report bugs and questions email: robpearc@umich.edu 

import torch
from torch import nn

from network import msa, pair
from network.config import network_config


class MSA_Transformer_Block(nn.Module):
    def __init__(self,config):
        super(MSA_Transformer_Block,self).__init__()
        self.drop_msa = nn.Dropout(0.15)
        self.drop_pair_row = nn.Dropout(0.25)
        self.drop_pair_col = nn.Dropout(0.25)
        self.msa_row_att = msa.MSA_Row_Att(config,config["msa_row_att"])
        self.msa_col_att = msa.MSA_Col_Att(config,config["msa_col_att"])
        self.msa_transition = msa.MSA_Transition(config,config["msa_transition"])
        self.msa_outer_product_mean = msa.MSA_Outer_Product_Mean(config,config["msa_outer_product_mean"])
        self.triangle_update_outgoing = pair.Triangle_Update_Outgoing(config,config["tri_out"])
        self.triangle_update_incoming = pair.Triangle_Update_Incoming(config,config["tri_in"])
        self.triangle_att_start = pair.Triangle_Att_Start(config,config["tri_att_start"])
        self.triangle_att_end  = pair.Triangle_Att_End(config,config["tri_att_end"])
        self.pair_transition = pair.Pair_Transition(config,config["pair_transition"])
        
    def forward(self,msa_embed,pair_embed):
        msa_embed = msa_embed + self.drop_msa(self.msa_row_att(msa_embed,pair_embed))
        msa_embed = msa_embed + self.msa_col_att(msa_embed)
        msa_embed = msa_embed + self.msa_transition(msa_embed)
        pair_embed = pair_embed + self.msa_outer_product_mean(msa_embed)
        pair_embed = pair_embed + self.drop_pair_row(self.triangle_update_outgoing(pair_embed))
        pair_embed = pair_embed + self.drop_pair_row(self.triangle_update_incoming(pair_embed))
        pair_embed = pair_embed + self.drop_pair_row(self.triangle_att_start(pair_embed))
        pair_embed = pair_embed + self.drop_pair_col(self.triangle_att_end(pair_embed))
        pair_embed = pair_embed + self.pair_transition(pair_embed)
        return msa_embed, pair_embed

class MSA_Transformer(nn.Module):
    def __init__(self,config):
        super(MSA_Transformer,self).__init__()
        self.no_msa_trans_blocks = config["no_msa_transformer_blocks"]
        self.msa_transformer_blocks = nn.ModuleList()
        for _ in range(self.no_msa_trans_blocks):
            block = MSA_Transformer_Block(config)
            self.msa_transformer_blocks.append(block)
    
    def forward(self,msa_embed,pair_embed):
        for i in range(self.no_msa_trans_blocks):
            msa_embed, pair_embed = self.msa_transformer_blocks[i](msa_embed,pair_embed)
            
        return msa_embed, pair_embed

