# Author: Robin Pearce
# To report bugs and questions email: robpearc@umich.edu 

import torch
from torch import nn

from network import pair


class Pair_Transformer_Block(nn.Module):
    def __init__(self,config):
        super(Pair_Transformer_Block,self).__init__()
        self.drop_pair_row = nn.Dropout(0.25)
        self.drop_pair_col = nn.Dropout(0.25)
        self.triangle_update_outgoing = pair.Triangle_Update_Outgoing(config,config["tri_out"])
        self.triangle_update_incoming = pair.Triangle_Update_Incoming(config,config["tri_in"])
        self.triangle_att_start = pair.Triangle_Att_Start(config,config["tri_att_start"])
        self.triangle_att_end  = pair.Triangle_Att_End(config,config["tri_att_end"])
        self.pair_transition = pair.Pair_Transition(config,config["pair_transition"])
        
    def forward(self,pair_embed):
        pair_embed = pair_embed + self.drop_pair_row(self.triangle_update_outgoing(pair_embed))
        pair_embed = pair_embed + self.drop_pair_row(self.triangle_update_incoming(pair_embed))
        pair_embed = pair_embed + self.drop_pair_row(self.triangle_att_start(pair_embed))
        pair_embed = pair_embed + self.drop_pair_col(self.triangle_att_end(pair_embed))
        pair_embed = pair_embed + self.pair_transition(pair_embed)
        return pair_embed

class Pair_Transformer(nn.Module):
    def __init__(self,config):
        super(Pair_Transformer,self).__init__()
        self.no_pair_trans_blocks = config["no_pair_transformer_blocks"]
        self.pair_transformer_blocks = nn.ModuleList()
        for _ in range(self.no_pair_trans_blocks):
            block = Pair_Transformer_Block(config)
            self.pair_transformer_blocks.append(block)        
    
    def forward(self,pair_embed):
        for i in range(self.no_pair_trans_blocks):
            pair_embed = self.pair_transformer_blocks[i](pair_embed)
            
        return pair_embed

