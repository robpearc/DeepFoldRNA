################################################################################

      ______               ______    _     _______ _   _   ___  
      |  _  \              |  ___|  | |   | | ___ \ \ | | / _ \ 
      | | | |___  ___ _ __ | |_ ___ | | __| | |_/ /  \| |/ /_\ \
      | | | / _ \/ _ \ '_ \|  _/ _ \| |/ _` |    /| . ` ||  _  |
      | |/ /  __/  __/ |_) | || (_) | | (_| | |\ \| |\  || | | |
      |___/ \___|\___| .__/\_| \___/|_|\__,_\_| \_\_| \_/\_| |_/
                     | |                                        
                     |_|                                        

                      (Version 1.0, 03/28/2022)

(Copyrighted by the Regents of the University of Michigan, All rights reserved)

DeepFoldRNA is a program for accurate de novo RNA tertiary structure prediction
using Deep Learning

Author: Robin Pearce

To report bugs and questions email: robpearc@umich.edu

If you use this program, please cite:
Pearce, R. Omenn, GS, and Zhang Y. De Novo RNA Tertiary Structure 
Prediction at Atomic Resolution Using Geometric Potentials from Deep Learning, 
Submitted, 2022.

Note, this is the stand-alone program, users may also submit jobs to 
https://zhanggroup.org/DeepFoldRNA/

The source code is freely available to academic/non-profit users under the 
PolyForm Noncommercial License

################################################################################



######################### INSTALLATION INSTRUCTIONS ############################

1.  System requirements: x86_64 machine, Linux kernel OS, Free disk space of more 
	than 500GB. The vast majority of the disk space will be used to store the 
	sequence databases.


2.  From the package directory, first install the third party programs by running the command:
	./scripts/install_third_party.sh

	Note, users are responsible for adhering to the licenses for these programs, 
		which include:
		
		PETfold: used for secondary structure prediction 
		(https://rth.dk/resources/petfold/)
		This program was developed Rolf Backofen's and Søren Brun's Lab
		Citation: Seemann SE, Gorodkin J, Backofen R. Nucleic Acids Res., 
			36(20):6355-62, 2008

		rMSA: used for MSA generation (https://github.com/pylelab/rMSA)
		This program was developed by Chengxin Zhang at Anna Pyle's Lab	
                Citation: Zhang C, Zhang Y, Pyle AM (2021) rMSA: accurate multiple 
			sequence alignment generation to improve RNA structure modeling. 
			ISMB, webinar.

		SimRNA and QRNAS: used for structure refinement 
		(http://genesilico.pl/simrna/, http://genesilico.pl/software/stand-alone/qrnas)
		These programs were developed by Janusz M. Bujnicki's Lab (http://genesilico.pl)
		SimRNA citation: Boniecki et al. Nucleic Acids Res. 2016 Apr 20;44(7):e63. 
			doi: 10.1093/nar/gkv1479.
		QRNAS citation: Stasiewicz et al. BMC Struct Biol 19, 5 (2019). 
			https://doi.org/10.1186/s12900-019-0103-1
		
		SPOT-RNA-1D: used to generate initial conformations for the folding simulations 
		(https://sparks-lab.org/server/spot-rna-1d/)
		This program was developed by Yaoqi Zhou's Lab
		Citation: Singh et al. J Chem Inf Model. 2021 Jun 28;61(6):2610-2622. 
			doi: 10.1021/acs.jcim.1c00153.


3.  Install conda locally by running the command:
	./scripts/install_conda.sh


4.  Install SPOT-RNA-1D dependencies by running the command:
	source scripts/activate_conda_env.sh
	conda activate venv
	pip install -r bin/SPOT-RNA-1D/requirements.txt


5.  To install the sequence databases used for MSA generation, run the command:
	./scripts/install_sequence_database.sh

	Note, the sequence databases require >500 GB of storage space and may take a 
	number of hours to download and process. The sequence databases will be 
	saved at <package_dir>/bin/rMSA/database

################################################################################



####################### RUNNING INSTRUCTIONS ###################################

1.  First activate the conda environment by running the command:
	source scripts/activate_conda_env.sh


2.  If you don't wish to run multiple threads set OMP_NUM_THREADS to 1 by:
	export OMP_NUM_THREADS=1

	Note, if this option isn't set properly, it may cause the modeling to 
	fail, so it is safer to set the number of threads to 1.


3.  Next, run DeepFoldRNA by running the command:

	python3 runDeepFoldRNA.py --input_dir <path to input directory>


	The only file required in the input directory is the input fasta file, which 
	should be named seq.fasta and contain two lines:

	>5v17A
	GGAUCAACCCCAGGUGUGGCACACCAGUCAUACCUUGAUCC

	where the first line is the description of the RNA and the second line is 
	the entire uppercase RNA sequence (valid nucleic acids: A, G, C, U)

	If GPU resources are available, the modeling will automatically run on the GPU, 
	which will be much faster than running the model on a CPU.
	
	Following successful modeling, the unrefined full atom models will be saved as 
	unrefined_model.pdb, while the refined models will be saved as refined_model.pdb 
	in the directories <input_dir>/model_$NUM/final_models. By default, 6 models will 
	be generated using different neural network parameters.

	Note, the unrefined model will typcally have slightly better backbone RMSD, but 
	worse base-pairing than the refined model. This is because the deep learning 
	restraints don't include the base atoms. Also, if there are chain breaks in the 
	refined model, which may be caused during SimRNA refinement due to restraints on 
	the backbone heavy atoms, or clashes, you can run more refinement steps by providing 
	the flag --num_refinement_steps and specifying a value greater than the default 5000 steps.

        Optional Flags:

		--generate_all_models (True/False) Generate all 6 models. If false, 3 models will 
			be generated. True by default.

		--num_refinement_steps Number of QRNAS refinement steps to perform (must be an integer 
			value, Default=5000).

		--generate_input_only (True/False) Only run the input generation step. Will output 
			the generated MSA (seq.afa), predicted secondary structure (ss.txt), hidden 
			Markov model (seq.cm), and the predicted torsion angles (eta.txt, theta.txt, 
			and chi.txt). False by default.

		--predict_restraints_only (True/False) Only run the deep learning restraint prediction 
			step. Note, the following files must be present in the input directory: 
			seq.fasta, seq.afa, seq.cm, and ss.txt. False by default.

		--run_folding_only (True/False) Only run the folding pipeline. Note, the following 
			files must be present in the input model directories: seq.fasta, DeepRNA_40.npz, 
			eta.txt, theta.txt, and chi.txt. False by default.

		--run_refinement_only (True/False) Only run the refinement pipeline. Note, the file 
			unrefined_model.pdb must be present in the input model directories. False by default.


################################################################################
