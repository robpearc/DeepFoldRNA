#!/usr/bin/env perl

########## Get Directories ############
my $bin_dir = "$ARGV[0]";
my $qrna_dir = "$ARGV[1]";
my $datadir = "$ARGV[2]";
my $num_steps = "$ARGV[3]";

chdir "$datadir";

############ Generate Full-Atom Model ###################
`cp -r $bin_dir/data ./`;
`$bin_dir/SimRNA -p unrefined_model.pdb -o unrefined -n 0`;
`$bin_dir/SimRNA_trafl2pdbs unrefined-000001.pdb unrefined.trafl 1:1 AA`;
`cp  unrefined-000001_AA.pdb unrefined_model_full_atom.pdb`;
`rm unrefined-000001_AA.pdb`;
`rm unrefined.trafl`;
`rm unrefined-000001.ss_detected`;
`rm unrefined-000001.pdb`;
`rm unrefined.bonds`;


################ Fix Backbone C3' and P Atoms During Refinement ###################
chomp(my @pdb = `cat unrefined_model_full_atom.pdb | grep -v "TER" | grep -v "END"`);
open(FILE,">pdb_for_refinement.pdb");
foreach my $line (@pdb){
    if(substr($line,13,3) eq "C3\'" || substr($line,13,3) eq "P  " ){
        substr($line,56,4)="0.00";
    }
    else{
        substr($line,56,4)="1.00";
    }
    print FILE "$line\n";
}
close(FILE);


################ Run Refinement  ###################
`$bin_dir/SimRNA -P pdb_for_refinement.pdb -o refine -n 1000`;
chomp(my @lowest_e = `python2 $bin_dir/trafl_extract_lowestE_frame.py refine.trafl`);
my @num = split(/\s+/,$lowest_e[0]);
`$bin_dir/SimRNA_trafl2pdbs refine-000001.pdb refine.trafl $num[0]:$num[0] AA`;

my $length = length($num[0]);

if($length<2){
    `cp refine-00000$num[0]_AA.pdb intermediate_refined_model.pdb`;
    `rm refine-00000$num[0]_AA.pdb`;
    `rm refine-00000$num[0].pdb`;
    `rm refine-00000$num[0].ss_detected`;
}
elsif($length<3){
    `cp refine-0000$num[0]_AA.pdb intermediate_refined_model.pdb`;
    `rm refine-0000$num[0]_AA.pdb`;
    `rm refine-0000$num[0].pdb`;
    `rm refine-0000$num[0].ss_detected`;
}
elsif($length<4){
    `cp refine-000$num[0]_AA.pdb intermediate_refined_model.pdb`;
    `rm refine-000$num[0]_AA.pdb`;
    `rm refine-000$num[0].pdb`;
    `rm refine-000$num[0].ss_detected`;
}
elsif($length<5){
    `cp refine-00$num[0]_AA.pdb intermediate_refined_model.pdb`;
    `rm refine-00$num[0]_AA.pdb`;
    `rm refine-00$num[0].pdb`;
    `rm refine-00$num[0].ss_detected`;
}
`rm refine-000001.pdb`;
`rm refine-000001.ss_detected`;
`rm refine.bonds`;
`rm refine_minE.trafl`;
`rm refine.trafl`;
`rm -r data`;

open(FILE,">configfile.txt");
print FILE "WRITEFREQ  1000\n";
print FILE "NSTEPS     $num_steps\n";
close(FILE);

`$qrna_dir/QRNA -i intermediate_refined_model.pdb -o refined_model.pdb -c configfile.txt`;

####### Save Ouputs ###### 
`mkdir final_models`;
`cp refined_model.pdb final_models/`;
`cp unrefined_model_full_atom.pdb final_models/unrefined_model.pdb`;

####### Clean Up #########
`rm unrefined_model_full_atom.pdb`;
`rm intermediate_refined_model.pdb`;
`rm refined_model.pdb`;
`rm unrefined_model.pdb`;
`rm pdb_for_refinement.pdb`;
