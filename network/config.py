# Author: Robin Pearce
# To report bugs and questions email: robpearc@umich.edu 


network_config={
    "input_seq_dim": 6,
    "input_msa_dim": 7,
    "msa_dim": 32,
    "pair_dim": 32,
    "seq_dim": 64,
    "no_msa_transformer_blocks": 48,
    "no_single_transformer_blocks": 4,
    "no_pair_transformer_blocks": 2,
    "no_cycles": 4,

    "input_embedder": {
        "no_pos_bins_1": 14,
        "no_pos_bins_2": 65,
        "rel_pos_1d": 14,
        "min_rel_range": -32,
        "max_rel_range": 32,
        "max_seq_length": 2000,
    },

    "s_att_pair_bias": {
        "no_heads": 8,
        "hidden_dim": 16,
    },
    
    "s_att": {
        "no_heads": 8,
        "hidden_dim": 16,
    },

    "s_outer_product_mean": {
        "no_heads": 8,
        "hidden_dim": 16,
    },

    "msa_row_att": {
        "no_heads": 8,
        "hidden_dim": 16,
    },
    
    "msa_col_att": {
        "no_heads": 8,
        "hidden_dim": 16,
    },
    
    "msa_transition": {
        "no_heads": 8,
        "hidden_dim": 16,
    },

    "msa_outer_product_mean": {
        "no_heads": 8,
        "hidden_dim": 16,
    },
    
    "tri_out": {
        "hidden_dim": 32,
    },

    "tri_in": {
        "hidden_dim": 32,
    },

    "tri_att_start": {
        "no_heads": 8,
        "hidden_dim": 16,
    },
    
    "tri_att_end": {
        "no_heads": 8,
        "hidden_dim": 16,
    },
    
    "pair_transition": {
        "hidden_dim": 2,
    },

    "torsion_head": {
        "hidden_dim": 64,
        "no_blocks": 2,
        "no_angles": 9,
        "angle_bins": 24,
    },

    "geometry_head": {
        "min_dist": 2.0,
        "max_dist": 40.0,
        "dist_bins": 38+2,
        "omega_bins": 24+1,
        "theta_bins": 24+1,
        "phi_bins": 12+1,
    }
    
}
