#!/usr/bin/env bash

cd bin/

###### Install MSA generation programs ############
git clone https://github.com/kad-ecoli/rMSA

######## Install programs for structure refinement ############
wget --no-check-certificate https://ftp.users.genesilico.pl/software/simrna/version_3.20/SimRNA_64bitIntel_Linux.tgz
tar -xvf SimRNA_64bitIntel_Linux.tgz
rm SimRNA_64bitIntel_Linux.tgz

wget http://genesilico.pl/QRNAS/QRNAS.tar.gz
gunzip QRNAS.tar.gz
tar -xvf QRNAS.tar
rm QRNAS.tar

###### Install program used to predict torsion angles for initial conformation generation #######
git clone https://github.com/jaswindersingh2/SPOT-RNA-1D.git
cd SPOT-RNA-1D
wget -O checkpoints.tar.xz 'https://www.dropbox.com/s/r9hp20gk30unptf/checkpoints.tar.xz' || wget -O checkpoints.tar.xz 'https://app.nihaocloud.com/f/4d2385c633554ccaa85c/?dl=1'
tar -xvf checkpoints.tar.xz && rm checkpoints.tar.xz
