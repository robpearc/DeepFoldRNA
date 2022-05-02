/*******************************************************************************************
 *   ______               ______    _     _ 
 *   |  _  \              |  ___|  | |   | |
 *   | | | |___  ___ _ __ | |_ ___ | | __| |
 *   | | | / _ \/ _ \ '_ \|  _/ _ \| |/ _` |
 *   | |/ /  __/  __/ |_) | || (_) | | (_| |
 *   |___/ \___|\___| .__/\_| \___/|_|\__,_|
 *                  | |                     
 *                  |_|                      
 *
 *  This program was written by Robin Pearce at the Zhang Lab
 *  Department of Computational Medicine and Bioinformatics 
 *  University of Michigan 
 *           
 *  Please report bugs and questions to robpearc@umich.edu
 *
*******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include "CubicSpline.h"
#include "ParseSeq.h"
#include "Geometry.h"
#include "Energy.h"
#include "LBFGS_optimizer.h" 
#include "CommonParameters.h"
#include "Operations.h"

using namespace std;

Energy energyFunction;


void show_interface(){
    printf(
      "#########################################################################\n"
      "                                                                         \n"
      "          ______               ______    _     _                         \n"
      "          |  _  \\              |  ___|  | |   | |                        \n"
      "          | | | |___  ___ _ __ | |_ ___ | | __| |                        \n"
      "          | | | / _ \\/ _ \\ '_ \\|  _/ _ \\| |/ _` |                        \n"
      "          | |/ /  __/  __/ |_) | || (_) | | (_| |                        \n"
      "          |___/ \\___|\\___| .__/\\_| \\___/|_|\\__,_|                        \n"
      "                         | |                                             \n"
      "                         |_|                                             \n"
      "                                                                         \n"
      "       A program for de novo RNA structure prediction guided by          \n"
      "                 potentials from deep learning                           \n"
      "\n\n"
      "  Written by Robin Pearce at the Yang Zhang Lab\n"
      "  Dept. of Computational Medicine & Bioinformatics\n"
      "  University of Michigan\n"
      "  For questions email robpearc@umich.edu\n"
      "#########################################################################\n");
}


int main(int argc, char** argv){
    string seqfile = argv[1]; 
    string datadir = argv[2]; 
    string libdir  = argv[3];

    show_interface();

    //------- read sequence and secondary structure ---------->
    ParseSeq ps,ps2;
    Geometry geo;
    
    if (!ps.loadseq(seqfile.c_str())) return 1;
    int numseq=ps.seqnum;

    string restraint_path_40 = datadir+"/DeepRNA_40.npz";
    string weight_file1 = libdir+"/weight_stg1.txt";
    string weight_file2 = libdir+"/weight_stg2.txt";
    string weight_file3 = libdir+"/weight_stg3.txt";
    string theta_file=datadir+"/theta.txt";
    string eta_file=datadir+"/eta.txt";
    string chi_file=datadir+"/chi.txt";

    // Set up energy function
    energyFunction.load_files( numseq, datadir, libdir, weight_file1, restraint_path_40 );

    // Generate initial decoy
    point3f *mcdecstr=new point3f[numseq];
    double *vars = new double[2*(numseq-1)];
    
    for ( int i=0; i<numseq; i++ )
    {
        mcdecstr[i].nt=ps.seqdata[i];
        mcdecstr[i].ntind=(unsigned char)(getntid(ps.seqdata[i]));
        if(mcdecstr[i].ntind>4) mcdecstr[i].ntind=4;
    }
    geo.setinitdecoyfromfile( mcdecstr, numseq, theta_file.c_str(), eta_file.c_str(), chi_file.c_str() );

    for ( int i=1; i<numseq; i++ )
    {
        vars[((i-1)*2)]=mcdecstr[i].theta;
        vars[((i-1)*2)+1]=mcdecstr[i].eta;
    }

    // Parameters for LBFGS
    int max_iterations = 2000;
    int scale_factor = 2.0;
    int problem_dim = 2 * ( numseq - 1 );
    int MAX_CYCLE = 10;
    int TRAJECTORIES = 3;
    string converge_test = "absolute";

    // Run LBFGS folding, probably only long proteins will use the full 10 cycles
    for ( int i=0; i<TRAJECTORIES; i++ )
    {
           
        if ( i>0 )
        {
            for ( int k=1; k<numseq; k++ )
            {
                double delta_eta = Random()*10.0;
                double delta_theta = Random()*10.0;
                if(Random()<0.5){
                    delta_eta*=-1.0;
                }
                if(Random()<0.5){
                    delta_theta*=-1.0;
                }
                mcdecstr[k].eta+=delta_eta;
                mcdecstr[k].theta+=delta_theta;

                if(mcdecstr[k].eta>=360.0){
                    mcdecstr[k].eta-=360.0;
                }
                else if(mcdecstr[k].eta<0.0){
                    mcdecstr[k].eta+=360.0;
                }
                if(mcdecstr[k].theta>=360.0){
                    mcdecstr[k].theta-=360.0;
                }
                else if(mcdecstr[k].theta<0.0){
                    mcdecstr[k].theta+=360.0;
                }

                vars[((k-1)*2)]=mcdecstr[k].theta;
                vars[((k-1)*2)+1]=mcdecstr[k].eta;
            }
        }

        for ( int j=0; j<MAX_CYCLE; j++ ) 
        {
            Minimizer fold;
            cout << "Running LBFGS cycle " << j+1 << " of " << MAX_CYCLE << endl; 
            fold.run( vars, mcdecstr, problem_dim, max_iterations, converge_test );
        }
        if(i==0){
            energyFunction.read_energy_weights( weight_file2.c_str() );
        }
        else if(i==1){
            energyFunction.read_energy_weights( weight_file3.c_str() );
        }
    }
    
    // Output final model
    {
        geo.tor2coord_base(mcdecstr,numseq);
        FILE *fp;
        string outfile_name = datadir+"/unrefined_model.pdb";
        fp=fopen(outfile_name.c_str(),"wt");
        int indatom=1;
        for ( int i=0; i<numseq; i++ )
        {
            const char *resn;
            int nt_idx=0;
            for ( int j=0; j<5;j++)
            {
                if (mcdecstr[i].nt==ntid[j])
                {
                    nt_idx=j;
                    resn=aad3[j];
                    break;
                }
            }

            fprintf(fp, "ATOM  %5d %s %s  %4d    %8.3f%8.3f%8.3f\n",			 
                indatom++, " P  ", resn, i+1,
                mcdecstr[i].ptp.x,     mcdecstr[i].ptp.y,     mcdecstr[i].ptp.z);
            fprintf(fp, "ATOM  %5d %s %s  %4d    %8.3f%8.3f%8.3f\n",			 
                indatom++, " C4\'", resn, i+1,
                mcdecstr[i].x,     mcdecstr[i].y,     mcdecstr[i].z);
            
            if(nt_idx==0 || nt_idx==2){
                fprintf(fp, "ATOM  %5d %s %s  %4d    %8.3f%8.3f%8.3f\n",			 
                    indatom++, " N9 ", resn, i+1,
                    mcdecstr[i].ptn.x,     mcdecstr[i].ptn.y,     mcdecstr[i].ptn.z);
                fprintf(fp, "ATOM  %5d %s %s  %4d    %8.3f%8.3f%8.3f\n",			 
                    indatom++, " C2 ", resn, i+1,
                    mcdecstr[i].ptc2.x,     mcdecstr[i].ptc2.y,     mcdecstr[i].ptc2.z);
                fprintf(fp, "ATOM  %5d %s %s  %4d    %8.3f%8.3f%8.3f\n",			 
                    indatom++, " C6 ", resn, i+1,
                    mcdecstr[i].ptc6.x,     mcdecstr[i].ptc6.y,     mcdecstr[i].ptc6.z);
            }
            else{
                fprintf(fp, "ATOM  %5d %s %s  %4d    %8.3f%8.3f%8.3f\n",			 
                    indatom++, " N1 ", resn, i+1,
                    mcdecstr[i].ptn.x,     mcdecstr[i].ptn.y,     mcdecstr[i].ptn.z);
                fprintf(fp, "ATOM  %5d %s %s  %4d    %8.3f%8.3f%8.3f\n",			 
                    indatom++, " C2 ", resn, i+1,
                    mcdecstr[i].ptc2.x,     mcdecstr[i].ptc2.y,     mcdecstr[i].ptc2.z);
                fprintf(fp, "ATOM  %5d %s %s  %4d    %8.3f%8.3f%8.3f\n",			 
                    indatom++, " C4 ", resn, i+1,
                    mcdecstr[i].ptc4.x,     mcdecstr[i].ptc4.y,     mcdecstr[i].ptc4.z);
           }
        }
        fclose(fp);
    }
  
    // Free memory
    delete[]mcdecstr;
    delete[]vars;
    mcdecstr=NULL;
    vars=NULL;

    return 0;
}


