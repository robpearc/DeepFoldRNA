# Author: Robin Pearce
# To report bugs and questions email: robpearc@umich.edu 

import torch
from torch import nn

from network import single, pair


class Single_Transformer_Block(nn.Module):
    def __init__(self,config):
        super(Single_Transformer_Block, self).__init__()
        self.drop_s = nn.Dropout(0.15)
        self.drop_pair_row = nn.Dropout(0.25)
        self.drop_pair_col = nn.Dropout(0.25)
        self.s_att_pair_bias = single.S_Att_Pair_Bias(config,config["s_att_pair_bias"])
        self.att_s = single.S_Att(config,config["s_att"])
        self.s_trans = single.S_Transition(config)

    def forward(self,seq_embed,pair_embed):
        seq_embed = seq_embed + self.drop_s(self.s_att_pair_bias(seq_embed,pair_embed))
        seq_embed = seq_embed + self.att_s(seq_embed)
        seq_embed = seq_embed + self.s_trans(seq_embed)

        return seq_embed

class Single_Transformer(nn.Module):
    def __init__(self,config):
        super(Single_Transformer,self).__init__()
        self.no_s_trans_blocks = config["no_single_transformer_blocks"]
        self.s_transformer_blocks = nn.ModuleList()
        for _ in range(self.no_s_trans_blocks):
            block = Single_Transformer_Block(config)
            self.s_transformer_blocks.append(block)

    def forward(self,seq_embed,pair_embed):
        for i in range(self.no_s_trans_blocks):
            seq_embed = self.s_transformer_blocks[i](seq_embed,pair_embed)

        return seq_embed
