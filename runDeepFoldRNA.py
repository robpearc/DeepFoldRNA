# Author: Robin Pearce
# To report bugs and questions email: robpearc@umich.edu

import os
import argparse
import subprocess
import torch
import numpy as np
import random
import math

random.seed(0)
np.random.seed(0)
torch.manual_seed(0)

from network import network_1, network_2, features 

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
fold_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),'fold')
network_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),'network')
bin_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),'bin')
home_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
os.environ['PETFOLDBIN'] = bin_dir

input_file_list = ['seq.fasta','eta.txt','theta.txt','chi.txt']

def generate_input(input_dir):
    if os.path.exists(input_dir+'/seq.afa') and os.path.exists(input_dir+'/seq.cm'):
        pass
    else:
        cmd = "perl " + bin_dir + "/rMSA/rMSA.pl " + input_dir + "/seq.fasta"
        os.system(cmd)

    if os.path.exists(input_dir+'/ss.txt'):
        pass
    else:
        cmd = "head -n 2000 " + input_dir + "/seq.afa > " + input_dir + "/msa_for_petfold.afa"
        os.system(cmd)

        cmd = bin_dir + "/PETfold" + " -f " + input_dir + "/msa_for_petfold.afa" + " -r " + input_dir + "/ss.txt"
        os.system(cmd)

    if os.path.exists(input_dir+'/eta.txt') and os.path.exists(input_dir+'theta.txt') and os.path.exists(input_dir+'/chi.txt'):
        pass
    else:
        cmd = "mkdir " + input_dir + "/angles"
        os.system(cmd)

        sequence_lines=open(input_dir+"/seq.fasta").readlines()
        sequence = sequence_lines[1].strip()
        fp = open(input_dir+"/angles/seq.fasta",'w')
        fp.write(">angles\n")
        fp.write(sequence)
        fp.close()
        
        cmd = '. ' + home_dir + '/conda_local/conda/etc/profile.d/conda.sh && conda activate venv && python3 ' + bin_dir + "/SPOT-RNA-1D/run.py --seq_file " + input_dir + "/angles/seq.fasta --save_outputs " + input_dir + "/angles"
        subprocess.call(cmd, shell=True, executable='/bin/bash')
        
        cmd = "perl " + bin_dir + "/parse_angles.pl " + input_dir + "/angles/angles.txt " + input_dir
        os.system(cmd)


def save_restraints(input_dir,pred_dis_n_ens,pred_dis_c4_ens,pred_dis_p_ens,pred_omg_ens,pred_theta_ens,
                    pred_phi_ens,pred_eta_bb_ens,pred_theta_bb_ens,pred_chi_ens,norm):
    if not os.path.exists(input_dir):
        cmd = "mkdir " + input_dir
        os.system(cmd)

    L,_,_ = pred_dis_n_ens.shape
    pred_dis_n_ens/=norm
    pred_dis_c4_ens/=norm
    pred_dis_p_ens/=norm
    pred_omg_ens/=norm
    pred_theta_ens/=norm
    pred_phi_ens/=norm
    pred_eta_bb_ens/=norm
    pred_theta_bb_ens/=norm
    pred_chi_ens/=norm

    pred_dis_n_ens = pred_dis_n_ens.cpu().data.numpy()
    pred_dis_c4_ens = pred_dis_c4_ens.cpu().data.numpy()
    pred_dis_p_ens = pred_dis_p_ens.cpu().data.numpy()
    pred_omg_ens = pred_omg_ens.cpu().data.numpy()
    pred_theta_ens = pred_theta_ens.cpu().data.numpy()
    pred_phi_ens = pred_phi_ens.cpu().data.numpy()
    pred_eta_bb_ens = pred_eta_bb_ens.cpu().data.numpy()
    pred_theta_bb_ens = pred_theta_bb_ens.cpu().data.numpy()

    ### save output restraints ###
    restraints = {}
    restraints['p_dist_n'] = np.zeros((L,L,39))
    restraints['p_dist_c4'] = np.zeros((L,L,39))
    restraints['p_dist_p'] = np.zeros((L,L,39))
    restraints['p_omega'] = np.zeros((L,L,25))
    restraints['p_theta'] = np.zeros((L,L,25))
    restraints['p_phi'] = np.zeros((L,L,13))
    restraints['p_eta_bb'] = np.zeros((L,24))
    restraints['p_theta_bb'] = np.zeros((L,24))
                
    restraints['p_eta_bb'] = pred_eta_bb_ens
    restraints['p_theta_bb'] = pred_theta_bb_ens
                
    for j in range(1,39):
        restraints['p_dist_n'][:,:,0] = pred_dis_n_ens[:,:,-1]
        restraints['p_dist_n'][:,:,j] = pred_dis_n_ens[:,:,j]
        restraints['p_dist_c4'][:,:,0] = pred_dis_c4_ens[:,:,-1]
        restraints['p_dist_c4'][:,:,j] = pred_dis_c4_ens[:,:,j]
        restraints['p_dist_p'][:,:,0] = pred_dis_p_ens[:,:,-1]
        restraints['p_dist_p'][:,:,j] = pred_dis_p_ens[:,:,j]
                
    for j in range(1,25):
        restraints['p_omega'][:,:,0] = pred_omg_ens[:,:,-1]
        restraints['p_omega'][:,:,j] = pred_omg_ens[:,:,j-1]
        restraints['p_theta'][:,:,0] = pred_theta_ens[:,:,-1]
        restraints['p_theta'][:,:,j] = pred_theta_ens[:,:,j-1]

    for j in range(1,13):
        restraints['p_phi'][:,:,0] = pred_phi_ens[:,:,-1]
        restraints['p_phi'][:,:,j] = pred_phi_ens[:,:,j-1]
                
    outfilename = input_dir+'/'+'DeepRNA_40'
    npz_file = outfilename+'.npz'
    np.savez_compressed(npz_file, dist_n=restraints['p_dist_n'], dist_c4=restraints['p_dist_c4'], dist_p=restraints['p_dist_p'], omega=restraints['p_omega'], theta=restraints['p_theta'], phi=restraints['p_phi'], eta_bb=restraints['p_eta_bb'], theta_bb=restraints['p_theta_bb'])

def run_prediction(input_dir,generate_all_models=True):
    num = 3
    if generate_all_models:
        num=6
    for i in range(num):
        name = "model_"+str(i+1)
        model_dir_name = os.path.join(input_dir,name)
        if os.path.exists(model_dir_name+'/DeepRNA_40.npz'):
            print("####### Restraints exist, skipping restraint generation. #######\n")
            return

    weight_file_location_1 = os.path.join(network_dir,'parameters','model_1')
    weight_file_location_2 = os.path.join(network_dir,'parameters','model_2')
    weight_file_location_1_2 = os.path.join(network_dir,'parameters','model_1_2')
    weight_file_location_2_2 = os.path.join(network_dir,'parameters','model_2_2')
    weights_model_1 = os.listdir(weight_file_location_1)
    weights_model_2 = os.listdir(weight_file_location_2)
    weights_model_1_2 = os.listdir(weight_file_location_1_2)
    weights_model_2_2 = os.listdir(weight_file_location_2_2)
    denominator_1 = len(weights_model_1) + len(weights_model_2)
    denominator_2 = len(weights_model_1_2) + len(weights_model_2_2)

    sequence_file = input_dir + "/seq.fasta"
    msa_file = input_dir + "/seq.afa"
    hmm_file = input_dir + "/seq.cm"
    ss_file = input_dir + "/ss.txt"
    with torch.no_grad():
        feature_dict = features.collect_features(sequence_file,msa_file,hmm_file,ss_file)
        for k, v in feature_dict.items():
            feature_dict[k] = v.to(device)
        
        N,L,_ = feature_dict["msa"].shape
        pred_dis_n_ens = torch.zeros([L,L,40]).to(device)
        pred_dis_c4_ens = torch.zeros([L,L,40]).to(device)
        pred_dis_p_ens = torch.zeros([L,L,40]).to(device)
        pred_omg_ens = torch.zeros([L,L,25]).to(device)
        pred_theta_ens = torch.zeros([L,L,25]).to(device)
        pred_phi_ens = torch.zeros([L,L,13]).to(device)
        pred_chi_ens = torch.zeros([L,24]).to(device)
        pred_eta_bb_ens = torch.zeros([L,24]).to(device)
        pred_theta_bb_ens = torch.zeros([L,24]).to(device)

        pred_dis_n_1 = torch.zeros([L,L,40]).to(device)
        pred_dis_c4_1 = torch.zeros([L,L,40]).to(device)
        pred_dis_p_1 = torch.zeros([L,L,40]).to(device)
        pred_omg_1 = torch.zeros([L,L,25]).to(device)
        pred_theta_1 = torch.zeros([L,L,25]).to(device)
        pred_phi_1 = torch.zeros([L,L,13]).to(device)
        pred_chi_1 = torch.zeros([L,24]).to(device)
        pred_eta_bb_1 = torch.zeros([L,24]).to(device)
        pred_theta_bb_1 = torch.zeros([L,24]).to(device)

        pred_dis_n_2 = torch.zeros([L,L,40]).to(device)
        pred_dis_c4_2 = torch.zeros([L,L,40]).to(device)
        pred_dis_p_2 = torch.zeros([L,L,40]).to(device)
        pred_omg_2 = torch.zeros([L,L,25]).to(device)
        pred_theta_2 = torch.zeros([L,L,25]).to(device)
        pred_phi_2 = torch.zeros([L,L,13]).to(device)
        pred_chi_2 = torch.zeros([L,24]).to(device)
        pred_eta_bb_2 = torch.zeros([L,24]).to(device)
        pred_theta_bb_2 = torch.zeros([L,24]).to(device)

        if True:
            for parameter_file in weights_model_1:
                print('####### Predicting Restraints for Network 1 Parameter File %s #######\n' %parameter_file)
                model = network_1.Pred_Network(device)
                model.eval()
                model.to(device)
                model_dict = torch.load(os.path.join(weight_file_location_1,parameter_file), map_location=device)
                model.load_state_dict(model_dict,strict=True)

                output = model(feature_dict)

                pred_dis_n = torch.exp(output["geometry_head"]["pred_dis_n"])
                pred_dis_c4 = torch.exp(output["geometry_head"]["pred_dis_c4"])
                pred_dis_p = torch.exp(output["geometry_head"]["pred_dis_p"])
                pred_omg = torch.exp(output["geometry_head"]["pred_omg"])
                pred_theta = torch.exp(output["geometry_head"]["pred_theta"])
                pred_phi = torch.exp(output["geometry_head"]["pred_phi"])
                pred_eta_bb = torch.exp(output["torsion_head"]["angles_dis"][:,7,:])
                pred_theta_bb = torch.exp(output["torsion_head"]["angles_dis"][:,8,:])
                pred_chi = torch.exp(output["torsion_head"]["angles_dis"][:,6,:])

                pred_dis_n_1+=pred_dis_n
                pred_dis_c4_1+=pred_dis_c4
                pred_dis_p_1+=pred_dis_p
                pred_omg_1+=pred_omg
                pred_theta_1+=pred_theta
                pred_phi_1+=pred_phi
                pred_eta_bb_1+=pred_eta_bb
                pred_theta_bb_1+=pred_theta_bb
                pred_chi_1+=pred_chi

            model_out_dir = os.path.join(input_dir,'model_2')
            norm = len(weights_model_1)
            save_restraints(model_out_dir,pred_dis_n_1,pred_dis_c4_1,pred_dis_p_1,pred_omg_1,pred_theta_1,
                            pred_phi_1,pred_eta_bb_1,pred_theta_bb_1,pred_chi_1,norm)

            for parameter_file in weights_model_2:
                print('####### Predicting Restraints for Network 2 Parameter File %s #######\n' %parameter_file)
                model = network_2.Pred_Network(device)
                model.eval()
                model.to(device)
                model_dict = torch.load(os.path.join(weight_file_location_2,parameter_file), map_location=device)
                model.load_state_dict(model_dict,strict=True)

                output = model(feature_dict)

                pred_dis_n = torch.exp(output["geometry_head"]["pred_dis_n"])
                pred_dis_c4 = torch.exp(output["geometry_head"]["pred_dis_c4"])
                pred_dis_p = torch.exp(output["geometry_head"]["pred_dis_p"])
                pred_omg = torch.exp(output["geometry_head"]["pred_omg"])
                pred_theta = torch.exp(output["geometry_head"]["pred_theta"])
                pred_phi = torch.exp(output["geometry_head"]["pred_phi"])
                pred_eta_bb = torch.exp(output["torsion_head"]["angles_dis"][:,7,:])
                pred_theta_bb = torch.exp(output["torsion_head"]["angles_dis"][:,8,:])
                pred_chi = torch.exp(output["torsion_head"]["angles_dis"][:,6,:])

                pred_dis_n_2+=pred_dis_n
                pred_dis_c4_2+=pred_dis_c4
                pred_dis_p_2+=pred_dis_p
                pred_omg_2+=pred_omg
                pred_theta_2+=pred_theta
                pred_phi_2+=pred_phi
                pred_eta_bb_2+=pred_eta_bb
                pred_theta_bb_2+=pred_theta_bb
                pred_chi_2+=pred_chi

            model_out_dir = os.path.join(input_dir,'model_3')
            norm = len(weights_model_2)
            save_restraints(model_out_dir,pred_dis_n_2,pred_dis_c4_2,pred_dis_p_2,pred_omg_2,pred_theta_2,
                            pred_phi_2,pred_eta_bb_2,pred_theta_bb_2,pred_chi_2,norm)

        pred_dis_n_ens = pred_dis_n_1 + pred_dis_n_2
        pred_dis_c4_ens = pred_dis_c4_1 + pred_dis_c4_2
        pred_dis_p_ens = pred_dis_p_1 + pred_dis_p_2
        pred_omg_ens = pred_omg_1 + pred_omg_2
        pred_theta_ens = pred_theta_1 + pred_theta_2
        pred_phi_ens = pred_phi_1 + pred_phi_2
        pred_eta_bb_ens = pred_eta_bb_1 + pred_eta_bb_2
        pred_theta_bb_ens = pred_theta_bb_1 + pred_theta_bb_2
        pred_chi_ens = pred_chi_1 + pred_chi_2

        model_out_dir = os.path.join(input_dir,'model_1')
        norm = denominator_1
        save_restraints(model_out_dir,pred_dis_n_ens,pred_dis_c4_ens,pred_dis_p_ens,pred_omg_ens,pred_theta_ens,
                        pred_phi_ens,pred_eta_bb_ens,pred_theta_bb_ens,pred_chi_ens,norm)

        if generate_all_models:
            pred_dis_n_ens = torch.zeros([L,L,40]).to(device)
            pred_dis_c4_ens = torch.zeros([L,L,40]).to(device)
            pred_dis_p_ens = torch.zeros([L,L,40]).to(device)
            pred_omg_ens = torch.zeros([L,L,25]).to(device)
            pred_theta_ens = torch.zeros([L,L,25]).to(device)
            pred_phi_ens = torch.zeros([L,L,13]).to(device)
            pred_chi_ens = torch.zeros([L,24]).to(device)
            pred_eta_bb_ens = torch.zeros([L,24]).to(device)
            pred_theta_bb_ens = torch.zeros([L,24]).to(device)

            pred_dis_n_1 = torch.zeros([L,L,40]).to(device)
            pred_dis_c4_1 = torch.zeros([L,L,40]).to(device)
            pred_dis_p_1 = torch.zeros([L,L,40]).to(device)
            pred_omg_1 = torch.zeros([L,L,25]).to(device)
            pred_theta_1 = torch.zeros([L,L,25]).to(device)
            pred_phi_1 = torch.zeros([L,L,13]).to(device)
            pred_chi_1 = torch.zeros([L,24]).to(device)
            pred_eta_bb_1 = torch.zeros([L,24]).to(device)
            pred_theta_bb_1 = torch.zeros([L,24]).to(device)

            pred_dis_n_2 = torch.zeros([L,L,40]).to(device)
            pred_dis_c4_2 = torch.zeros([L,L,40]).to(device)
            pred_dis_p_2 = torch.zeros([L,L,40]).to(device)
            pred_omg_2 = torch.zeros([L,L,25]).to(device)
            pred_theta_2 = torch.zeros([L,L,25]).to(device)
            pred_phi_2 = torch.zeros([L,L,13]).to(device)
            pred_chi_2 = torch.zeros([L,24]).to(device)
            pred_eta_bb_2 = torch.zeros([L,24]).to(device)
            pred_theta_bb_2 = torch.zeros([L,24]).to(device)

            for parameter_file in weights_model_1_2:
                print('####### Predicting Restraints for Network 3 Parameter File %s #######\n' %parameter_file)
                model = network_1.Pred_Network(device)
                model.eval()
                model.to(device)
                model_dict = torch.load(os.path.join(weight_file_location_1_2,parameter_file), map_location=device)
                model.load_state_dict(model_dict,strict=True)

                output = model(feature_dict)

                pred_dis_n = torch.exp(output["geometry_head"]["pred_dis_n"])
                pred_dis_c4 = torch.exp(output["geometry_head"]["pred_dis_c4"])
                pred_dis_p = torch.exp(output["geometry_head"]["pred_dis_p"])
                pred_omg = torch.exp(output["geometry_head"]["pred_omg"])
                pred_theta = torch.exp(output["geometry_head"]["pred_theta"])
                pred_phi = torch.exp(output["geometry_head"]["pred_phi"])
                pred_eta_bb = torch.exp(output["torsion_head"]["angles_dis"][:,7,:])
                pred_theta_bb = torch.exp(output["torsion_head"]["angles_dis"][:,8,:])
                pred_chi = torch.exp(output["torsion_head"]["angles_dis"][:,6,:])

                pred_dis_n_1+=pred_dis_n
                pred_dis_c4_1+=pred_dis_c4
                pred_dis_p_1+=pred_dis_p
                pred_omg_1+=pred_omg
                pred_theta_1+=pred_theta
                pred_phi_1+=pred_phi
                pred_eta_bb_1+=pred_eta_bb
                pred_theta_bb_1+=pred_theta_bb
                pred_chi_1+=pred_chi

            model_out_dir = os.path.join(input_dir,'model_5')
            norm = len(weights_model_1_2)
            save_restraints(model_out_dir,pred_dis_n_1,pred_dis_c4_1,pred_dis_p_1,pred_omg_1,pred_theta_1,
                            pred_phi_1,pred_eta_bb_1,pred_theta_bb_1,pred_chi_1,norm)

            for parameter_file in weights_model_2_2:
                print('####### Predicting Restraints for Network 4 Parameter File %s #######\n' %parameter_file)
                model = network_2.Pred_Network(device)
                model.eval()
                model.to(device)
                model_dict = torch.load(os.path.join(weight_file_location_2_2,parameter_file), map_location=device)
                model.load_state_dict(model_dict,strict=True)

                output = model(feature_dict)

                pred_dis_n = torch.exp(output["geometry_head"]["pred_dis_n"])
                pred_dis_c4 = torch.exp(output["geometry_head"]["pred_dis_c4"])
                pred_dis_p = torch.exp(output["geometry_head"]["pred_dis_p"])
                pred_omg = torch.exp(output["geometry_head"]["pred_omg"])
                pred_theta = torch.exp(output["geometry_head"]["pred_theta"])
                pred_phi = torch.exp(output["geometry_head"]["pred_phi"])
                pred_eta_bb = torch.exp(output["torsion_head"]["angles_dis"][:,7,:])
                pred_theta_bb = torch.exp(output["torsion_head"]["angles_dis"][:,8,:])
                pred_chi = torch.exp(output["torsion_head"]["angles_dis"][:,6,:])

                pred_dis_n_2+=pred_dis_n
                pred_dis_c4_2+=pred_dis_c4
                pred_dis_p_2+=pred_dis_p
                pred_omg_2+=pred_omg
                pred_theta_2+=pred_theta
                pred_phi_2+=pred_phi
                pred_eta_bb_2+=pred_eta_bb
                pred_theta_bb_2+=pred_theta_bb
                pred_chi_2+=pred_chi

            model_out_dir = os.path.join(input_dir,'model_6')
            norm = len(weights_model_2_2)
            save_restraints(model_out_dir,pred_dis_n_2,pred_dis_c4_2,pred_dis_p_2,pred_omg_2,pred_theta_2,
                            pred_phi_2,pred_eta_bb_2,pred_theta_bb_2,pred_chi_2,norm)

            pred_dis_n_ens = pred_dis_n_1 + pred_dis_n_2
            pred_dis_c4_ens = pred_dis_c4_1 + pred_dis_c4_2
            pred_dis_p_ens = pred_dis_p_1 + pred_dis_p_2
            pred_omg_ens = pred_omg_1 + pred_omg_2
            pred_theta_ens = pred_theta_1 + pred_theta_2
            pred_phi_ens = pred_phi_1 + pred_phi_2
            pred_eta_bb_ens = pred_eta_bb_1 + pred_eta_bb_2
            pred_theta_bb_ens = pred_theta_bb_1 + pred_theta_bb_2
            pred_chi_ens = pred_chi_1 + pred_chi_2

            model_out_dir = os.path.join(input_dir,'model_4')
            norm = denominator_2
            save_restraints(model_out_dir,pred_dis_n_ens,pred_dis_c4_ens,pred_dis_p_ens,pred_omg_ens,pred_theta_ens,
                            pred_phi_ens,pred_eta_bb_ens,pred_theta_bb_ens,pred_chi_ens,norm)

def run_folding(input_dir,generate_all_models):
    num = 3
    if generate_all_models:
        num=6
    for i in range(num):
        name = "model_"+str(i+1)    
        model_dir_name = os.path.join(input_dir,name)
        for f in input_file_list:
            cmd = "cp " + input_dir + "/" + f + " " + model_dir_name
            os.system(cmd)
        
        cmd = "python3 " + bin_dir + "/npz2txt.py " + model_dir_name + "/DeepRNA_40.npz"
        os.system(cmd)
    
        cmd = fold_dir + "/DeepFold " + model_dir_name + "/seq.fasta " + model_dir_name + " " + fold_dir + "/library/"
        os.system(cmd)

def run_refinement(input_dir,num_refinement_steps,generate_all_models):
    num = 3
    if generate_all_models:
        num=6
    for i in range(num):
        name = "model_"+str(i+1)
        model_dir_name = os.path.join(input_dir,name)
        cmd = "perl " + bin_dir + "/refine.pl " + bin_dir + "/SimRNA_64bitIntel_Linux " + bin_dir + "/QRNAS " + model_dir_name + " " + str(num_refinement_steps)
        os.system(cmd)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', type=str, required=True, help="Required input. Path to the directory where the input fasta sequence file is saved (Must be named seq.fasta). The first line of the fasta file should start with a > and the second line must contain the entire upper case RNA sequence (valid characters: A,G,C,U).")

    parser.add_argument('--generate_all_models', type=bool, default=True, help="Optional input. Generate all 6 models. If false, 3 models will be generated. (True/False) Default: True")

    parser.add_argument('--num_refinement_steps', type=int, default=5000, help="Optional input. Number of QRNAS refinement steps to perform. Default value: 5000")
    
    parser.add_argument('--generate_input_only', type=bool, default=False, help="Optional input. Only run the input generation step. Will output the generated MSA (seq.afa), predicted secondary structure (ss.txt), hidden Markov model (seq.cm), and the predicted torsion angles (eta.txt, theta.txt, and chi.txt). (True/False) Default value: False")
    
    parser.add_argument('--predict_restraints_only', type=bool, default=False, help="Optional input. Only run the deep learning restraint prediction step. Note, the following files must be present in the input directory: seq.fasta, seq.afa, seq.cm, and ss.txt. (True/False) Default value: False")
    
    parser.add_argument('--run_folding_only', type=bool, default=False, help="Optional input. Only run the folding pipeline. Note, the following files must be present in the input directory: seq.fasta, DeepRNA_40.npz, eta.txt, theta.txt, and chi.txt. (True/False) Default value: False")

    parser.add_argument('--run_refinement_only', type=bool, default=False, help="Optional input. Only run the refinement pipeline. Note, the following file must be present in the input model directories: unrefined_model.pdb. (True/False) Default value: False")

    args = parser.parse_args()

    run_mode = 1
    if args.generate_input_only:
        run_mode = 2
    elif args.predict_restraints_only:
        run_mode = 3
    elif args.run_folding_only:
        run_mode = 4
    elif args.run_refinement_only:
        run_mode = 5

    if run_mode==1 or run_mode==2:
        print('####### Generating Input #######\n')
        generate_input(args.input_dir)
        print('####### Input Generation Completed #######\n')
    
    if run_mode==1 or run_mode==3:
        print('####### Predicting restraints. This may take some time if running on a CPU. ########\n')
        run_prediction(args.input_dir,args.generate_all_models)
        print('####### Restraint Prediction Completed #######\n')
    
    if run_mode==1 or run_mode==4:
        print('####### Running Folding Simulations #######\n')
        run_folding(args.input_dir,args.generate_all_models)
        print('####### Folding Simulations Completed #######\n')

    if run_mode==1 or run_mode==5:
        print('####### Running Structure Refinement #######\n')
        run_refinement(args.input_dir,args.num_refinement_steps,args.generate_all_models)
        print('####### Refinement Completed #######\n')
