/*******************************************************************************************
**  
**  Functions for calculating energy terms and gradients.
**
**  Please report bugs and questions to robpearc@umich.edu
**
*******************************************************************************************/

#ifndef ENERGY_H
#define ENERGY_H

#include "Operations.h"
#include "CubicSpline.h"
#include "CommonParameters.h"
#include "Derivatives.h"

class Energy  
{
public:
    Energy();

    bool load_dfire( const char *filename );
    bool load_general_contact( const char *filename );
    void load_files( int seqnum, string datadir, string libdir, string weight_file,
                     string restraint_path_40 );

    void torsion_energy(point3f *decstr, double *enelist);
    void torsion_energy(point3f *decstr, double *enelist,double *dE_dvars );


    /* prepare pn compact data structure to host the xyz coordinate. initialize
     * Epair_mat, Epair_tmp, CA_SC_mat, CA_SC_tmp, CB_CT_mat, CB_CT_tmp */
    void genpn();

    void pairwise_orientation( point3f *decstr );

    //double energytorsion2(point3f *decstr,int numseq);//log
    void calcallenergy(point3f *decstr, double *enelist);

    /* Calculates energy and gradients when grad_flag is true. Gradients are stored in dE_dvars */
    double calcrmsdenergy( point3f *decstr, double *vars, double *dE_dvars, bool grad_flag );

    /* Only calculates energy without gradients */
    double calcrmsdenergy( point3f *decstr, double *vars );
    
    /* Calculates the backbone hydrogen bonding energy*/
    void energy_hbond( point3f *decstr, double *enelist );

    void read_eta_bb( string Eta_path );
    void read_theta_bb( string Theta_path );

    bool read_energy_weights( const char *weightfile );
    void read_cb_distance_40( string distance_path );
    void read_p_distance_40( string distance_path );
    void read_ca_distance_40( string distance_path );
    void read_omega_40( string omega_path );
    void read_phi_40( string phi_path );
    void read_theta_40( string theta_path );

    virtual ~Energy();

private:
    Derivatives deriv;
    Geometry geo;
  
    int number_atoms=3;
    int SHORT_RANGE_MIN, SHORT_RANGE_MAX;
    int MED_RANGE_MIN, MED_RANGE_MAX;
    int LONG_RANGE_MIN;
    int numseq;
    bool flag_grad;
  
    int nbin_40 = 39;
    int nbin_20 = 37;
    int nbin_10 = 33;

    int nbin_omega = 25;
    int nbin_theta = nbin_omega;
    int nbin_phi = 13;
    
    int nbin_theta_bb = 24;
    int nbin_eta_bb = 24;
    point3d **atom_derivatives_f1;
    point3d **atom_derivatives_f2;

    double **Pcontact;
    double **PcontactCA;
    double ***dist_energy_new_40;
    double  **cb_dist_prob_new_40;
    double   *cb_dist_bin_value_new_40;
 
    double ***p_dist_energy_40;
    double  **p_dist_prob_40;
    double   *p_dist_bin_value_40; 

    double ***ca_dist_energy_40;
    double  **ca_dist_prob_40;
    double   *ca_dist_bin_value_40; 

    double ***omega_energy_40;
    double  **omega_prob_40;
    double   *omega_bin_value_40;

    double **theta_energy_bb;
    double   *theta_bin_value_bb;

    double **eta_energy_bb;
    double   *eta_bin_value_bb;
 
    double ***theta_energy_40;
    double  **theta_prob_40;
    double   *theta_bin_value_40;
 
    double ***phi_energy_40;
    double  **phi_prob_40;
    double   *phi_bin_value_40; 

    bool flagTheta, flagPhi, flagOmega, flagPdist;
    CubicSpline *spline_eta_bb;
    CubicSpline *spline_theta_bb;

    CubicSpline **spline_dist_new_40;
    CubicSpline **spline_dist_ca_40;
    CubicSpline **spline_dist_p_40;
    CubicSpline **spline_omega_40;
    CubicSpline **spline_theta_40;
    CubicSpline **spline_phi_40;

    bool flag_ca_dist,flag_ca_cont,flag_cb_cont;
    double dist_epsilon;

    /*
    int solnum;
    double *solseq;
    */

    double weights[MAX_ENERGY_TERM_NUM];
    double dist_cut_cb_40; 
    double dist_cut_ca_40; 
    double dist_cut_p_40; 
    double theta_cut_40; 
    double omg_cut_40; 
    double phi_cut_40; 

    // Coordinate arrays
    double **pn;
    double **Epair_mat;
    double **Epair_tmp;
    double **C_N_tmp;
    double **PHI_mat;
    double **PHI_tmp;
    double **THETA_mat;
    double **THETA_tmp;
    double **OMEGA_mat;
    double **OMEGA_tmp;

    double enelist[20];//used for output
     
    bool flagcont;//use contact restraints or not;
};

Energy::Energy()  // initialization of parameters/options
{
    SHORT_RANGE_MIN     =  1; 
    SHORT_RANGE_MAX     = 11;
    MED_RANGE_MIN       = 12; 
    MED_RANGE_MAX       = 23;
    LONG_RANGE_MIN      = 24;
    dist_epsilon = 1e-8;
    for ( int i=0; i<MAX_ENERGY_TERM_NUM; i++ )
    {
        weights[i] = 0.0;
    }   

    flag_grad=true;

    int i,j;
    spline_dist_new_40=NULL;
    spline_dist_ca_40=NULL;
    spline_dist_p_40=NULL;
    spline_eta_bb=NULL;
    spline_theta_bb=NULL;
    spline_omega_40=NULL;
    spline_theta_40=NULL;
    spline_phi_40=NULL;

    pn=NULL; // xyz
    Epair_mat=NULL;
    Epair_tmp=NULL;
    C_N_tmp=NULL;

    PHI_mat=NULL;   // (i,j) is <CAi-CBi-CBj>
    PHI_tmp=NULL;   // (i,j) is <CAi-CBi-CBj>
    THETA_mat=NULL; // (i,j) is <Ni-CAi-CBi-CBj>
    THETA_tmp=NULL; // (i,j) is <Ni-CAi-CBi-CBj>
    OMEGA_mat=NULL; // (i,j) and (j,i) are both <CAi-CBi-CBj-CAj>
    OMEGA_tmp=NULL; // (i,j) and (j,i) are both <CAi-CBi-CBj-CAj>

    flagcont=false;
}

Energy::~Energy() = default;

void Energy::load_files( int seqnum, string datadir, string libdir, string weight_file, 
                         string restraint_path_40 )
{
    read_energy_weights( weight_file.c_str() );
    
    numseq=seqnum;
    genpn();

    // Allocate memory for splines
    spline_dist_new_40 = new2DArrT< CubicSpline >( numseq, numseq );
    spline_dist_ca_40 = new2DArrT< CubicSpline >( numseq, numseq );
    spline_dist_p_40 = new2DArrT< CubicSpline >( numseq, numseq );
    spline_omega_40 = new2DArrT< CubicSpline >( numseq, numseq );
    spline_phi_40 = new2DArrT< CubicSpline >( numseq, numseq );
    spline_theta_40 = new2DArrT< CubicSpline >( numseq, numseq );
    spline_eta_bb = new1DArrT< CubicSpline >( numseq );
    spline_theta_bb = new1DArrT< CubicSpline >( numseq );

    // Allocate memory for atom derivatives
    atom_derivatives_f1 = new point3d*[numseq];
    atom_derivatives_f2 = new point3d*[numseq];
    for(int i=0;i<numseq;i++)
    {
        atom_derivatives_f1[i]=new point3d[number_atoms];
        atom_derivatives_f2[i]=new point3d[number_atoms];
    }

    // Deep learning restraints
    bool flagcon =true;//= checkfile(restraint_path);
    if( flagcon )
    {
        string distance_path_40 = restraint_path_40 + "_dist_n.txt";
        string distance_path_ca_40 = restraint_path_40 + "_dist_c4.txt";
        string distance_path_p_40 = restraint_path_40 + "_dist_p.txt";
        string theta_path_40 = restraint_path_40 + "_theta.txt";
        string omega_path_40 = restraint_path_40 + "_omega.txt";
        string phi_path_40 = restraint_path_40 + "_phi.txt";
        string eta_path_bb = restraint_path_40 + "_eta_bb.txt";
        string theta_path_bb = restraint_path_40 + "_theta_bb.txt";
        read_cb_distance_40( distance_path_40 );
        read_ca_distance_40( distance_path_ca_40 );
        read_p_distance_40( distance_path_p_40 );
        read_omega_40( omega_path_40 );
        read_theta_40( theta_path_40 );
        read_phi_40( phi_path_40 );
        read_eta_bb( eta_path_bb );
        read_theta_bb( theta_path_bb );
    }

    //----------------- Backbon torsion angles ------------------->
    for( int i=0; i<numseq; i++ ){
        CubicSpline SPL;
        SPL.set_points( eta_bin_value_bb, eta_energy_bb[i], nbin_eta_bb );
        spline_eta_bb[i] = SPL;
    }
    for( int i=0; i<numseq; i++ ){
        CubicSpline SPL;
        SPL.set_points( theta_bin_value_bb, theta_energy_bb[i], nbin_theta_bb );
        spline_theta_bb[i] = SPL;
    }
    
    //----------------- N Distance Energy ------------------->
    for( int i=0; i<numseq; i++ ){
        for( int j=i+1; j<numseq; j++ ){
            CubicSpline SPL;
            SPL.set_points( cb_dist_bin_value_new_40, dist_energy_new_40[i][j], nbin_40-1 );
            spline_dist_new_40[i][j] = SPL;
        }
    }

    //----------------- P Distance Energy ------------------->
    for( int i=0; i<numseq; i++ ){
        for( int j=i+1; j<numseq; j++ ){
            CubicSpline SPL;
            SPL.set_points( p_dist_bin_value_40, p_dist_energy_40[i][j], nbin_40-1 );
            spline_dist_p_40[i][j] = SPL;
        }
    }

    //----------------- C4 Distance Energy ------------------->
    for( int i=0; i<numseq; i++ ){
        for( int j=i+1; j<numseq; j++ ){
            CubicSpline SPL;
            SPL.set_points( ca_dist_bin_value_40, ca_dist_energy_40[i][j], nbin_40-1 );
            spline_dist_ca_40[i][j] = SPL;
        }
    }

    //----------------- OMEGA Energy ------------------->
    // Fit spline, note currently we do not institute periodic boundary conditions
    for( int i=0; i<numseq; i++ ){
        for( int j=i+1; j<numseq; j++ ){
            CubicSpline SPL;
            SPL.set_points( omega_bin_value_40, omega_energy_40[i][j], nbin_omega-1 );
            spline_omega_40[i][j] = SPL;
        }
    }

    //----------------- PHI Energy ------------------->
    // Fit spline, note currently we do not institute periodic boundary conditions
    for( int i=0; i<numseq; i++ ){
        for( int j=0; j<numseq; j++ ){
            CubicSpline SPL;
            SPL.set_points( phi_bin_value_40, phi_energy_40[i][j], nbin_phi-1+6 );
            spline_phi_40[i][j] = SPL;
        }
    }

    //----------------- THETA Energy ------------------->
    // Fit spline, note currently we do not institute periodic boundary conditions
    for( int i=0; i<numseq; i++ ){
        for( int j=0; j<numseq; j++ ){
            CubicSpline SPL;
            SPL.set_points( theta_bin_value_40, theta_energy_40[i][j], nbin_theta-1 );
            spline_theta_40[i][j] = SPL;
        }
    }
}

bool Energy::read_energy_weights( const char *weightfile )
{
    for(int i=0;i<MAX_ENERGY_TERM_NUM;i++)
    {
        weights[i]=0.0;
    }
    FILE* pf=fopen(weightfile,"r");
    if(pf!=NULL)
    {
        char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
        while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,pf))
        {
            char term[MAX_LENGTH_ONE_LINE_IN_FILE+1];
            double val=0.0;
            sscanf(line,"%s %lf",term,&val);
            if(!strcmp(term,"Deep_Short_Range_Dist_CB_40"))  weights[ 0]=val;
            else if(!strcmp(term,"Deep_Med_Range_Dist_CB_40"))    weights[ 1]=val;
            else if(!strcmp(term,"Deep_Long_Range_Dist_CB_40"))   weights[ 2]=val;
            else if(!strcmp(term,"Deep_Short_Range_Dist_CA_40"))  weights[ 3]=val;
            else if(!strcmp(term,"Deep_Med_Range_Dist_CA_40"))    weights[ 4]=val;
            else if(!strcmp(term,"Deep_Long_Range_Dist_CA_40"))   weights[ 5]=val;
            else if(!strcmp(term,"Deep_Short_Range_Cont_CB"))  weights[ 6]=val;
            else if(!strcmp(term,"Deep_Med_Range_Cont_CB"))    weights[ 7]=val;
            else if(!strcmp(term,"Deep_Long_Range_Cont_CB"))   weights[ 8]=val;
            else if(!strcmp(term,"Deep_Short_Range_Cont_CA"))  weights[ 9]=val;
            else if(!strcmp(term,"Deep_Med_Range_Cont_CA"))    weights[10]=val;
            else if(!strcmp(term,"Deep_Long_Range_Cont_CA"))   weights[11]=val;
            else if(!strcmp(term,"Deep_Short_Range_Omg_40"))      weights[12]=val;
            else if(!strcmp(term,"Deep_Med_Range_Omg_40"))        weights[13]=val;
            else if(!strcmp(term,"Deep_Long_Range_Omg_40"))       weights[14]=val;
            else if(!strcmp(term,"Deep_Short_Range_Theta_40"))    weights[15]=val;
            else if(!strcmp(term,"Deep_Med_Range_Theta_40"))      weights[16]=val;
            else if(!strcmp(term,"Deep_Long_Range_Theta_40"))     weights[17]=val;
            else if(!strcmp(term,"Deep_Short_Range_Phi_40"))      weights[18]=val;
            else if(!strcmp(term,"Deep_Med_Range_Phi_40"))        weights[19]=val;
            else if(!strcmp(term,"Deep_Long_Range_Phi_40"))       weights[20]=val;
            else if(!strcmp(term,"Deep_HB_AA"))                   weights[21]=val;
            else if(!strcmp(term,"Deep_HB_BB"))                   weights[22]=val;
            else if(!strcmp(term,"Deep_HB_CC"))                   weights[23]=val;
            else if(!strcmp(term,"General_Energy_VDW"))           weights[24]=val;
            else if(!strcmp(term,"General_Energy_HB"))            weights[25]=val;
            else if(!strcmp(term,"General_Energy_Dfire"))         weights[26]=val;
            else if(!strcmp(term,"General_Energy_SG_Contact"))    weights[27]=val;
            else if(!strcmp(term,"General_Energy_Rad_Gyr"))       weights[28]=val;
            else if(!strcmp(term,"Tor_Constr"))                   weights[29]=val;
            else if(!strcmp(term,"Deep_Short_Range_Dist_P_40"))   weights[30]=val;
            else if(!strcmp(term,"Deep_Med_Range_Dist_P_40"))     weights[31]=val;
            else if(!strcmp(term,"Deep_Long_Range_Dist_P_40"))    weights[32]=val;
            else if(!strcmp(term,"PCUT_DIST_P_40"))               dist_cut_p_40=val;
            else if(!strcmp(term,"PCUT_DIST_CB_40"))              dist_cut_cb_40=val;
            else if(!strcmp(term,"PCUT_DIST_CA_40"))              dist_cut_ca_40=val;
            else if(!strcmp(term,"PCUT_OMG_40"))                  omg_cut_40=val;
            else if(!strcmp(term,"PCUT_PHI_40"))                  phi_cut_40=val;
            else if(!strcmp(term,"PCUT_THETA_40"))                theta_cut_40=val;
        }
        fclose(pf);
    }
    else{
	cout << "The energy weights file was not found!" << endl;
	return false;
    }

    for(int i=0;i<MAX_ENERGY_TERM_NUM;i++)
    {
        weights[i]*=1.0;
    }

    return true;
}

void Energy::read_p_distance_40( string distance_path )
{

    double ***Pdist_p_40 = new3DArr( numseq, numseq, nbin_40 );
    p_dist_energy_40 = new3DArr( numseq, numseq, nbin_40-1 );
    p_dist_bin_value_40 = new1DArr( nbin_40 - 1 );
    p_dist_prob_40 = new2DArr( numseq, numseq );

    bool flag = checkfile( distance_path );
    if( !flag )
    {
        //transform_npz2txt( );
        flag = checkfile( distance_path );
        if( !flag )
        {
            cout<<"There is no predicted distance file: "<<distance_path<<endl;
            return;
        }
    }
        
    ifstream in ( distance_path.c_str() );
    string line;
    vector<string > content2;
    for(int i=0; i<nbin_40 + 2; i++){
        string x = "";
        content2.push_back(x);
    }
    while ( getline (in, line) ){
        stringstream word( line );
        for(int i=0; i<content2.size(); i++){
            word >> content2[i];
        }
        int a = atoi(content2[0].c_str());
        int b = atoi(content2[1].c_str());
        if( a > numseq || b > numseq ) continue;
        for(int i=2; i<content2.size(); i++){
            double prob = atof( content2[i].c_str() );
            Pdist_p_40[a-1][b-1][i-2] = prob;
            if(i==2){
                double prob2 = atof( content2[i].c_str() );
                p_dist_prob_40[a-1][b-1]=prob;
            }
        }
    }

    double dist_bin[39][2];
    double bin_width=1.0;
    dist_bin[0][0] = 40.0, dist_bin[0][1] = 999.0;
    dist_bin[1][0] =double(2.0+bin_width), dist_bin[1][1] = double(2.0+bin_width);
    p_dist_bin_value_40[0] = -0.0001;
    p_dist_bin_value_40[1] = double(2.0+bin_width)+double(bin_width/2.0);
    for(int i=3; i<nbin_40; i++ ){
        p_dist_bin_value_40[i-1] = p_dist_bin_value_40[i-2]+bin_width;//0.5*( dist_bin[i][0] + dist_bin[i][1] );
    }

    // Calculate Distance Energy
    for( int i=0; i<numseq-1; i++ ){
        for( int j=i+1; j<numseq; j++ ){
            double Pref = 1E-8;
            for( int k=1; k<nbin_40; k++ )
            {
                double score = -log( ( Pdist_p_40[i][j][k] + Pref ) / ( Pref + Pdist_p_40[i][j][nbin_40-1] ) );
                p_dist_energy_40[i][j][k-1] = score;
                p_dist_energy_40[j][i][k-1] = score;
            }
            double score0 = max( 10.0, p_dist_energy_40[i][j][1] + 4.0 );
            p_dist_energy_40[i][j][0] = score0;
            p_dist_energy_40[j][i][0] = score0;
       }
    }

    //release3DArr( numseq, numseq, Pdist_ca );
    cout<<"The predicted distance file is available."<<endl;
    flagPdist = true;
    in.close();
    //vector<string >().swap( content );
}

void Energy::read_ca_distance_40( string distance_path )
{

    double ***Pdist_ca_40 = new3DArr( numseq, numseq, nbin_40 );
    ca_dist_energy_40 = new3DArr( numseq, numseq, nbin_40-1 );
    ca_dist_bin_value_40 = new1DArr( nbin_40 - 1 );
    ca_dist_prob_40 = new2DArr( numseq, numseq );

    bool flag = checkfile( distance_path );
    if( !flag )
    {
        //transform_npz2txt( );
        flag = checkfile( distance_path );
        if( !flag )
        {
            cout<<"There is no predicted distance file: "<<distance_path<<endl;
            return;
        }
    }
        
    ifstream in ( distance_path.c_str() );
    string line;
    vector<string > content2;
    for(int i=0; i<nbin_40 + 2; i++){
        string x = "";
        content2.push_back(x);
    }
    while ( getline (in, line) ){
        stringstream word( line );
        for(int i=0; i<content2.size(); i++){
            word >> content2[i];
        }
        int a = atoi(content2[0].c_str());
        int b = atoi(content2[1].c_str());
        if( a > numseq || b > numseq ) continue;
        for(int i=2; i<content2.size(); i++){
            double prob = atof( content2[i].c_str() );
            Pdist_ca_40[a-1][b-1][i-2] = prob;
            if(i==2){
                double prob2 = atof( content2[i].c_str() );
                ca_dist_prob_40[a-1][b-1]=prob;
            }
        }
    }

    double dist_bin[39][2];
    double bin_width=1.0;
    dist_bin[0][0] = 40.0, dist_bin[0][1] = 999.0;
    dist_bin[1][0] =double(2.0+bin_width), dist_bin[1][1] = double(2.0+bin_width);
    ca_dist_bin_value_40[0] = -0.0001;
    ca_dist_bin_value_40[1] = double(2.0+bin_width)+double(bin_width/2.0);
    for(int i=3; i<nbin_40; i++ ){
        ca_dist_bin_value_40[i-1] = ca_dist_bin_value_40[i-2]+bin_width;//0.5*( dist_bin[i][0] + dist_bin[i][1] );
    }

    // Calculate Distance Energy
    for( int i=0; i<numseq-1; i++ ){
        for( int j=i+1; j<numseq; j++ ){
            double Pref = 1E-8;
            for( int k=1; k<nbin_40; k++ )
            {
                double score = -log( ( Pdist_ca_40[i][j][k] + Pref ) / ( Pref + Pdist_ca_40[i][j][nbin_40-1] ) );
                ca_dist_energy_40[i][j][k-1] = score;
                ca_dist_energy_40[j][i][k-1] = score;
            }
            double score0 = max( 10.0, ca_dist_energy_40[i][j][1] + 4.0 );
            ca_dist_energy_40[i][j][0] = score0;
            ca_dist_energy_40[j][i][0] = score0;
       }
    }

    //release3DArr( numseq, numseq, Pdist_ca );
    cout<<"The predicted distance file is available."<<endl;
    flagPdist = true;
    in.close();
    //vector<string >().swap( content );
}

void Energy::read_cb_distance_40( string distance_path )
{

    double ***Pdist_cb_new_40 = new3DArr( numseq, numseq, nbin_40 );
    dist_energy_new_40 = new3DArr( numseq, numseq, nbin_40-1 );
    cb_dist_bin_value_new_40 = new1DArr( nbin_40 - 1 );
    cb_dist_prob_new_40 = new2DArr( numseq, numseq );

    bool flag = checkfile( distance_path );
    if( !flag )
    {
        //transform_npz2txt( );
        flag = checkfile( distance_path );
        if( !flag )
        {
            cout<<"There is no predicted distance file: "<<distance_path<<endl;
            return;
        }
    }
        
    ifstream in ( distance_path.c_str() );
    string line;
    vector<string > content;
    for(int i=0; i<nbin_40 + 2; i++){
        string x = "";
        content.push_back(x);
    }
    while ( getline (in, line) ){
        stringstream word( line );
        for(int i=0; i<content.size(); i++){
            word >> content[i];
        }
        int a = atoi(content[0].c_str());
        int b = atoi(content[1].c_str());
        if( a > numseq || b > numseq ) continue;
        for(int i=2; i<content.size(); i++){
            double prob = atof( content[i].c_str() );
            Pdist_cb_new_40[a-1][b-1][i-2] = prob;
            if(i==2){
                double prob2 = atof( content[i].c_str() );
                cb_dist_prob_new_40[a-1][b-1]=prob;
            }
        }
    }

    double dist_bin[39][2];
    double bin_width=1.0;
    dist_bin[0][0] = 40.0, dist_bin[0][1] = 999.0;
    dist_bin[1][0] =double(2.0+bin_width), dist_bin[1][1] = double(2.0+bin_width);
    cb_dist_bin_value_new_40[0] = -0.0001;
    cb_dist_bin_value_new_40[1] = double(2.0+bin_width)+double(bin_width/2.0);
    for(int i=3; i<nbin_40; i++ ){
        cb_dist_bin_value_new_40[i-1] = cb_dist_bin_value_new_40[i-2]+bin_width;
    }

    // Calculate Distance Energy
    for( int i=0; i<numseq-1; i++ ){
        for( int j=i+1; j<numseq; j++ ){
            double Pref = 1E-8;
            for( int k=1; k<nbin_40; k++ )
            {
                double score = -log( ( Pdist_cb_new_40[i][j][k] + Pref ) / ( Pref + Pdist_cb_new_40[i][j][nbin_40-1] ) );
                dist_energy_new_40[i][j][k-1] = score;
                dist_energy_new_40[j][i][k-1] = score;
            }
            double score0 = max( 10.0, dist_energy_new_40[i][j][1] + 4.0 );
            dist_energy_new_40[i][j][0] = score0;
            dist_energy_new_40[j][i][0] = score0;
       }
    }

    //release3DArr( numseq, numseq, Pdist_cb );
    cout<<"The predicted distance file is available."<<endl;
    flagPdist = true;
    in.close();
    //vector<string >().swap( content );
}

void Energy::read_omega_40( string omega_path )
{
    double ***Pomega_40 = new3DArr( numseq, numseq, nbin_omega );
    omega_energy_40 = new3DArr( numseq, numseq, nbin_omega );
    omega_bin_value_40 = new1DArr( nbin_omega-1 );
    omega_prob_40 = new2DArr( numseq, numseq );

    bool flag = checkfile( omega_path );
    if( !flag ){
        //transform_npz2txt( );
        flag = checkfile( omega_path );
        if( !flag ){
            cout<<"There is no predicted omega angle: "<<omega_path<<endl;
            return;
        }
    }

    ifstream in ( omega_path.c_str() );
    string line;
    vector<string > content;
    for(int i=0; i<nbin_omega + 2; i++){
        string x = "";
        content.push_back(x);
    }
    while ( getline (in, line) ){
        stringstream word( line );
        for(int i=0; i<content.size(); i++){
            word >> content[i];
        }
        int a = atoi(content[0].c_str());
        int b = atoi(content[1].c_str());
        if( a > numseq || b > numseq ) continue;
        for(int i=2; i<content.size(); i++){
            Pomega_40[a-1][b-1][i-2] = atof( content[i].c_str() );
	    //cout << "Pomega " <<a <<" " <<b<<" "<<i-2 <<" " <<Pomega[a-1][b-1][i-2] <<endl;
            if(i==2) omega_prob_40[a-1][b-1] = atof( content[i].c_str() );
        }
    }

    double omega_bin[25][2];
    omega_bin[0][0]   = 180.0*raddeg, omega_bin[0][1] = 999.0;
    omega_bin[1][0]   = -180.0*raddeg, omega_bin[1][1] = -165.0*raddeg;
    omega_bin_value_40[0]= -180.0*raddeg+7.5*raddeg;
    for(int i=2; i<nbin_omega; i++ ){
        omega_bin[i][0] = omega_bin[i-1][0] + 15.0*raddeg;
        omega_bin[i][1] = omega_bin[i-1][1] + 15.0*raddeg;
        omega_bin_value_40[i-1] = 0.5*( omega_bin[i][0] + omega_bin[i][1] );
    }

    //calculate omega energy
    for( int i=0; i<numseq; i++ ){
        for( int j=i+1; j<numseq; j++ ){
            for( int k=0; k<nbin_omega; k++ ){
                double prob = Pomega_40[i][j][k+1];
                if( k == nbin_omega-1 ) prob = Pomega_40[i][j][1];
                double score = -1.0*log( ( prob + dist_epsilon ) );
                omega_energy_40[i][j][k] = score;
                omega_energy_40[j][i][k] = score;
            }
        }
    }

    release3DArr( numseq, numseq, Pomega_40 );
    cout<<"The predicted omega angle file is available."<<endl;
    flagOmega = true;
    in.close();
    //vector<string >().swap( content );
}

void Energy::read_eta_bb( string Eta_path )
{
    double **Peta_bb = new2DArr( numseq, nbin_eta_bb );
    eta_energy_bb = new2DArr( numseq, nbin_eta_bb );
    eta_bin_value_bb = new1DArr( nbin_eta_bb );

    bool flag = checkfile( Eta_path );
    if( !flag )
    {
        flag = checkfile( Eta_path );
        if( !flag )
        {
            cout << "There is no predicted eta angle: " << Eta_path << endl;
            return;
        }
    }

    ifstream in ( Eta_path.c_str() );
    string line;
    vector<string > content;
    for(int i=0; i<nbin_eta_bb + 1; i++){
        string x = "";
        content.push_back(x);
    }
    while ( getline (in, line) ){
        stringstream word( line );
        for(int i=0; i<content.size(); i++){
            word >> content[i];
        }
        int a = atoi(content[0].c_str());
        if( a > numseq ) continue;
        for(int i=1; i<content.size(); i++){
            Peta_bb[a-1][i-1] = atof( content[i].c_str() );
        }
    }

    double eta_bin[24][2];
    eta_bin[0][0]    =  180.0*raddeg, eta_bin[0][1] = 999.0;
    eta_bin[1][0]    = -180.0*raddeg, eta_bin[1][1] = -165.0*raddeg;
    eta_bin_value_bb[0] = -180.0*raddeg+7.5*raddeg;
    for(int i=2; i<nbin_eta_bb; i++ ){
        eta_bin[i][0] = eta_bin[i-1][0] + 15.0*raddeg;
        eta_bin[i][1] = eta_bin[i-1][1] + 15.0*raddeg;
        eta_bin_value_bb[i-1] = 0.5*( eta_bin[i][0] + eta_bin[i][1] );
    }
    eta_bin_value_bb[nbin_eta_bb-1] = 180.0*raddeg-7.5*raddeg;

    //calculate eta energy
    for( int i=0; i<numseq; i++ ){
        for( int k=0; k<nbin_eta_bb; k++ ){
            double prob = Peta_bb[i][k];
            double score = -1.0*log( ( prob + dist_epsilon ) );
            eta_energy_bb[i][k] = score;
        }
    }

    cout<<"The predicted eta angle file is available."<<endl;
    flagTheta = true;
    in.close();
    //vector<string >().swap( content );
}

void Energy::read_theta_bb( string Theta_path )
{
    double **Ptheta_bb = new2DArr( numseq, nbin_theta_bb );
    theta_energy_bb = new2DArr( numseq, nbin_theta_bb );
    theta_bin_value_bb = new1DArr( nbin_theta_bb );

    bool flag = checkfile( Theta_path );
    if( !flag )
    {
        flag = checkfile( Theta_path );
        if( !flag )
        {
            cout << "There is no predicted theta angle: " << Theta_path << endl;
            return;
        }
    }

    ifstream in ( Theta_path.c_str() );
    string line;
    vector<string > content;
    for(int i=0; i<nbin_theta_bb + 1; i++){
        string x = "";
        content.push_back(x);
    }
    while ( getline (in, line) ){
        stringstream word( line );
        for(int i=0; i<content.size(); i++){
            word >> content[i];
        }
        int a = atoi(content[0].c_str());
        if( a > numseq ) continue;
        for(int i=1; i<content.size(); i++){
            Ptheta_bb[a-1][i-1] = atof( content[i].c_str() );
        }
    }

    double theta_bin[24][2];
    theta_bin[0][0]    =  180.0*raddeg, theta_bin[0][1] = 999.0;
    theta_bin[1][0]    = -180.0*raddeg, theta_bin[1][1] = -165.0*raddeg;
    theta_bin_value_bb[0] = -180.0*raddeg+7.5*raddeg;
    for(int i=2; i<nbin_theta_bb; i++ ){
        theta_bin[i][0] = theta_bin[i-1][0] + 15.0*raddeg;
        theta_bin[i][1] = theta_bin[i-1][1] + 15.0*raddeg;
        theta_bin_value_bb[i-1] = 0.5*( theta_bin[i][0] + theta_bin[i][1] );
    }
    theta_bin_value_bb[nbin_theta_bb-1] = 180.0*raddeg-7.5*raddeg;


    //calculate theta energy
    for( int i=0; i<numseq; i++ ){
        for( int k=0; k<nbin_theta_bb; k++ ){
            double prob = Ptheta_bb[i][k];
            double score = -1.0*log( ( prob + dist_epsilon ) );
            theta_energy_bb[i][k] = score;
        }
    }

    cout<<"The predicted theta angle file is available."<<endl;
    flagTheta = true;
    in.close();
    //vector<string >().swap( content );
}

void Energy::read_theta_40( string Theta_path )
{
    double ***Ptheta_40 = new3DArr( numseq, numseq, nbin_theta );
    theta_energy_40 = new3DArr( numseq, numseq, nbin_theta );
    theta_prob_40 = new2DArr( numseq, numseq );
    theta_bin_value_40 = new1DArr( nbin_theta-1 );

    bool flag = checkfile( Theta_path );
    if( !flag )
    {
        flag = checkfile( Theta_path );
        if( !flag )
        {
            cout << "There is no predicted theta angle: " << Theta_path << endl;
            return;
        }
    }

    ifstream in ( Theta_path.c_str() );
    string line;
    vector<string > content;
    for(int i=0; i<nbin_theta + 2; i++){
        string x = "";
        content.push_back(x);
    }
    while ( getline (in, line) ){
        stringstream word( line );
        for(int i=0; i<content.size(); i++){
            word >> content[i];
        }
        int a = atoi(content[0].c_str());
        int b = atoi(content[1].c_str());
        if( a > numseq || b > numseq ) continue;
        for(int i=2; i<content.size(); i++){
            Ptheta_40[a-1][b-1][i-2] = atof( content[i].c_str() );
            if(i==2) theta_prob_40[a-1][b-1] = atof( content[i].c_str() );
        }
    }

    double theta_bin[25][2];
    theta_bin[0][0]    =  180.0*raddeg, theta_bin[0][1] = 999.0;
    theta_bin[1][0]    = -180.0*raddeg, theta_bin[1][1] = -165.0*raddeg;
    theta_bin_value_40[0] = -180.0*raddeg+7.5*raddeg;
    for(int i=2; i<nbin_theta; i++ ){
        theta_bin[i][0] = theta_bin[i-1][0] + 15.0*raddeg;
        theta_bin[i][1] = theta_bin[i-1][1] + 15.0*raddeg;
        theta_bin_value_40[i-1] = 0.5*( theta_bin[i][0] + theta_bin[i][1] );
    }

    //calculate theta energy
    for( int i=0; i<numseq; i++ ){
        for( int j=0; j<numseq; j++ ){
            for( int k=0; k<nbin_theta; k++ ){
                double prob = Ptheta_40[i][j][k+1];
                if( k == nbin_theta-1 ) prob = Ptheta_40[i][j][1];
                double score = -1.0*log( ( prob + dist_epsilon ) );
                theta_energy_40[i][j][k] = score;
            }
        }
    }

    cout<<"The predicted theta angle file is available."<<endl;
    flagTheta = true;
    in.close();
    //vector<string >().swap( content );
}

void Energy::read_phi_40( string Phi_path )
{
    double ***Pphi_40 = new3DArr( numseq, numseq, nbin_phi );
    phi_energy_40 = new3DArr( numseq, numseq, nbin_phi -1 +6);
    phi_bin_value_40 = new1DArr( nbin_phi-1+6 );
    phi_prob_40 = new2DArr( numseq, numseq );

    bool flag = checkfile( Phi_path );
    if ( !flag ){
        ///transform_npz2txt( );
        flag = checkfile( Phi_path );
        if( !flag )
        {
            cout<<"There is no predicted phi angle: "<<Phi_path<<endl;
                        return;
        }
    }

    ifstream in ( Phi_path.c_str() );
    string line;
    vector<string > content;
    for ( int i=0; i<nbin_phi + 2; i++ )
    {
        string x = "";
        content.push_back(x);
    }
    while ( getline (in, line) )
    {
        stringstream word( line );
        for( int i=0; i<content.size(); i++ )
        {
            word >> content[i];
        }
        int a = atoi(content[0].c_str());
        int b = atoi(content[1].c_str());
        if( a > numseq || b > numseq ) continue;
        for(int i=2; i<content.size(); i++)
        {
            Pphi_40[a-1][b-1][i-2] = atof( content[i].c_str() );
            if(i==2)
            {
                phi_prob_40[a-1][b-1] = atof( content[i].c_str() );
            }
        }
    }

    double phi_bin[13][2];
    phi_bin[0][0] = PI, phi_bin[0][1] = 999.0;
    phi_bin[1][0] = 0.0, phi_bin[1][1] = 15.0*raddeg;
    phi_bin_value_40[0]=7.5*raddeg-(15.0*3)*raddeg;
    phi_bin_value_40[1]=7.5*raddeg-(15.0*2)*raddeg;
    phi_bin_value_40[2]=7.5*raddeg-(15.0)*raddeg;
    phi_bin_value_40[3]=7.5*raddeg;

    for(int i=2; i<nbin_phi; i++ )
    {
        phi_bin[i][0] = phi_bin[i-1][0] + (15.0)*raddeg;
        phi_bin[i][1] = phi_bin[i-1][1] + (15.0)*raddeg;
        phi_bin_value_40[i-1+3] = 0.5*( phi_bin[i][0] + phi_bin[i][1] );
    }

    phi_bin_value_40[nbin_phi-1+3]=phi_bin_value_40[nbin_phi-2+3]+15.0*raddeg;
    phi_bin_value_40[nbin_phi-1+4]=phi_bin_value_40[nbin_phi-2+3]+15.0*2*raddeg;
    phi_bin_value_40[nbin_phi-1+5]=phi_bin_value_40[nbin_phi-2+3]+15.0*3*raddeg;
        
    //calculate phi energy
    for( int i=0; i<numseq; i++ )
    {
        for( int j=0; j<numseq; j++ )
        {
            if( i==j ) continue;
            for( int k=0; k<nbin_phi-1; k++ )
            {
                double prob = Pphi_40[i][j][k+1];
                double score = -1.0*log( ( prob + dist_epsilon ) );
                phi_energy_40[i][j][k+3] = score;
            }
            phi_energy_40[i][j][0]              = phi_energy_40[i][j][2+3];
            phi_energy_40[i][j][1]              = phi_energy_40[i][j][1+3];
            phi_energy_40[i][j][2]              = phi_energy_40[i][j][0+3];
            phi_energy_40[i][j][nbin_phi-1+3]   = phi_energy_40[i][j][nbin_phi-1-1+3];
            phi_energy_40[i][j][nbin_phi-1+1+3] = phi_energy_40[i][j][nbin_phi-1-2+3];
            phi_energy_40[i][j][nbin_phi-1+2+3] = phi_energy_40[i][j][nbin_phi-1-3+3];
        }
    }

    cout<<"The predicted phi angle file is available."<<endl;
    flagPhi = true;
    in.close();
    //vector<string >().swap( content );
}

void Energy::genpn()
{
    int seqnum=numseq;
    pn=new double*[number_atoms*seqnum];
    int i,j,ii;
    for (i=0;i<seqnum;i++)
    {
        for (ii=0;ii<number_atoms;ii++)
        {
            pn[i*number_atoms+ii]=new double[3];
            pn[i*number_atoms+ii][0]=0;
            pn[i*number_atoms+ii][1]=0;
            pn[i*number_atoms+ii][2]=0;
        }
    }

    Epair_mat=new double*[seqnum];
    Epair_tmp=new double*[seqnum];
    C_N_tmp=new double*[seqnum];
    PHI_mat  =new double*[seqnum];
    PHI_tmp  =new double*[seqnum];
    THETA_mat=new double*[seqnum];
    THETA_tmp=new double*[seqnum];
    OMEGA_mat=new double*[seqnum];
    OMEGA_tmp=new double*[seqnum];
    for (i=0;i<seqnum;i++)
    {
        Epair_mat[i]=new double[seqnum];
        Epair_tmp[i]=new double[seqnum];
        C_N_tmp[i]=new double[seqnum];

        PHI_mat[i]  =new double[seqnum];
        PHI_tmp[i]  =new double[seqnum];
        THETA_mat[i]=new double[seqnum];
        THETA_tmp[i]=new double[seqnum];
        OMEGA_mat[i]=new double[seqnum];
        OMEGA_tmp[i]=new double[seqnum];
        for (j=0;j<seqnum;j++)
        {
            Epair_mat[i][j]=0;
            Epair_tmp[i][j]=0;
            C_N_tmp[i][j]=0;
            PHI_mat[i][j]  =360;
            PHI_tmp[i][j]  =360;
            THETA_mat[i][j]=360;
            THETA_tmp[i][j]=360;
            OMEGA_mat[i][j]=360;
            OMEGA_tmp[i][j]=360;
        }
    }
}


/* calculate orientations of decoy structure */
void Energy::pairwise_orientation( point3f *decstr )
{
    int i,j,r1,r2;
    double THETA,OMEGA;
    point3d pd,pd2;
    for ( i=0; i<numseq; i++ )
    {
        r1=i*number_atoms;
        for ( j=i+1; j<numseq; j++ )
        {
            r2=j*number_atoms;

            if ( theta_prob_40[i][j]<theta_cut_40 )
            {
                // THETAij = <Ni-CAi-CBi-CBj>  , [-pi,pi]
                THETA=phi( pn[r1+1][0], pn[r1+1][1], pn[r1+1][2],
                           pn[r1  ][0], pn[r1  ][1], pn[r1  ][2],
                           pn[r1+2][0], pn[r1+2][1], pn[r1+2][2],
                           pn[r2+2][0], pn[r2+2][1], pn[r2+2][2] );
                if ( THETA>180 ) THETA-=360;
                THETA_tmp[i][j]=THETA;
            }

            if ( theta_prob_40[j][i]<theta_cut_40 )
            {
                // THETAji = <Nj-CAj-CBj-CBi>  , [-pi,pi]
                THETA=phi( pn[r2+1][0], pn[r2+1][1], pn[r2+1][2],
                           pn[r2  ][0], pn[r2  ][1], pn[r2  ][2],
                           pn[r2+2][0], pn[r2+2][1], pn[r2+2][2],
                           pn[r1+2][0], pn[r1+2][1], pn[r1+2][2] );
                if ( THETA>180 ) THETA-=360;
                THETA_tmp[j][i]=THETA;
            }

            if ( phi_prob_40[i][j]<phi_cut_40 )
            {
                // PHIij   = <CAi-CBi-CBj>     , [0,pi]
                pd.x  = pn[r1  ][0]-pn[r1+2][0];
                pd.y  = pn[r1  ][1]-pn[r1+2][1];
                pd.z  = pn[r1  ][2]-pn[r1+2][2];
                pd2.x = pn[r2+2][0]-pn[r1+2][0];
                pd2.y = pn[r2+2][1]-pn[r1+2][1];
                pd2.z = pn[r2+2][2]-pn[r1+2][2];
                PHI_tmp[i][j] = angv( pd, pd2 );
            }

            if ( phi_prob_40[j][i]<phi_cut_40 )
            {
                // PHIji   = <CAj-CBj-CBi>     , [0,pi]
                pd.x  = pn[r2  ][0]-pn[r2+2][0];
                pd.y  = pn[r2  ][1]-pn[r2+2][1];
                pd.z  = pn[r2  ][2]-pn[r2+2][2];
                pd2.x = pn[r1+2][0]-pn[r2+2][0];
                pd2.y = pn[r1+2][1]-pn[r2+2][1];
                pd2.z = pn[r1+2][2]-pn[r2+2][2];
                PHI_tmp[j][i] = angv( pd, pd2 );
            }

            if ( omega_prob_40[i][j]<omg_cut_40 )
            {
                // OMEGAij = <CAi-CBi-CBj-CAj> , [-pi,pi]
                OMEGA=phi( pn[r1  ][0], pn[r1  ][1], pn[r1  ][2],
                           pn[r1+2][0], pn[r1+2][1], pn[r1+2][2],
                           pn[r2+2][0], pn[r2+2][1], pn[r2+2][2],
                           pn[r2  ][0], pn[r2  ][1], pn[r2  ][2] );
                if ( OMEGA>180 ) OMEGA-=360;
                OMEGA_tmp[j][i] = OMEGA_tmp[i][j] = OMEGA;
            }
        }
    }
    return;
}


/* calculate most of the pairwise energy terms in QUARK 
 * enelist   - energy terms */
void Energy::calcallenergy( point3f *decstr, double *enelist )
{
    int i,j,ii,jj; // i,j are residue indexes; ii,jj are atom indexes
    int tbin,ind1,ind2,indm,ind5[5],ind6[5];
    point3d tp;
    point3d atom1,atom2,atom3,atom4;
    double tdist,tdist2,tdist3;
    double tdist_CA,tdist2_CA;
    double tdist_CB,tdist2_CB;
    double scalf=0.00707;
    double sqdist,tdists;
    double weight=0;

    enelist[0]=0;//I-TASSER restraints
    enelist[1]=0;//Excluded volume, ww1
    enelist[4]=0;//E_RWplus_backbone_only, RW part tuned by wRW
    enelist[5]=0;//pairwise side-chain atomic potential, 'sgpolarity5.txt', ww5
    enelist[6]=0;//sequence-specific solvation potential, 'sol.txt', ww6
    enelist[7]=0;//radias gyration, RG part by ww7
    enelist[8]=0;//distant profile from either fragments (QE/QN) or init.dat (QP/QT), wtdp
    enelist[12]=0;//Edist ResTriplet2 distance map
    enelist[14]=0;//CA bond-length, ww14
    enelist[15]=0;//Econtact
    enelist[16]=0;//Econtact
    enelist[17]=0;//Econtact

    int r1,r2;
    int idx; // for accessing distance map
    for ( i=0; i<numseq; i++ )
    {
        r1=i*number_atoms;

        pn[r1  ][0]=decstr[i].x;     //C4
        pn[r1  ][1]=decstr[i].y;
        pn[r1  ][2]=decstr[i].z;
        pn[r1+1][0]=decstr[i].ptp.x; //P
        pn[r1+1][1]=decstr[i].ptp.y;
        pn[r1+1][2]=decstr[i].ptp.z;
        pn[r1+2][0]=decstr[i].ptn.x; //N
        pn[r1+2][1]=decstr[i].ptn.y;
        pn[r1+2][2]=decstr[i].ptn.z;
    }

    /* fill up C_N_tmp and CB_CT_tmp */
    for(i=0;i<numseq;i++)
    {
        r1=i*number_atoms;

        for(j=i+1;j<numseq;j++)
        {
            idx=(j-1)+(2*numseq-3-i)*i/2;
            r2=j*number_atoms;

            /* CA-CA distance calculation */
            tp.x=pn[r1][0]-pn[r2][0];
            tp.y=pn[r1][1]-pn[r2][1];
            tp.z=pn[r1][2]-pn[r2][2];
            tdist2_CA=sqrt(tp.x*tp.x+tp.y*tp.y+tp.z*tp.z);
            C_N_tmp[i][j]=tdist2_CA;

            tp.x=pn[r1+2][0]-pn[r2+2][0];
            tp.y=pn[r1+2][1]-pn[r2+2][1];
            tp.z=pn[r1+2][2]-pn[r2+2][2];
            C_N_tmp[j][i]=sqrt(tp.x*tp.x+tp.y*tp.y+tp.z*tp.z);

        }
    }

    /* pairwise orientation terms
                               Cj
                              /
                     THETAji /
                   CBj----CAj             ---  covalent bond
                  .__/       \            ...  non-covalent interaction
                 .   PHIji    \
              D .              Nj
         OMEGA .                  
              .                   
             .                    
          CBi \                   Symmetric properties:         
           | _/ PHIij                 Dij     = |CBi-CBj|
    THETAij|                          OMEGAij = <CAi-CBi-CBj-CAj> , [-pi,pi]
           |                      Asymmetric properties:
          CAi                         THETAij = <Ni-CAi-CBi-CBj>  , [-pi,pi]
         /   \                        THETAji = <Nj-CAj-CBj-CBi>  , [-pi,pi]
        /     \                       PHIij   = <CAi-CBi-CBj>     , [0,pi]
      Ni       Ci                     PHIji   = <CAj-CBj-CBi>     , [0,pi] */
    
    {
        pairwise_orientation( decstr );

        for ( i=0; i<numseq; i++ )
        {
            int r1=i*number_atoms;
            for (j=i+1;j<numseq;j++)
            {
                int r2=j*number_atoms;

                double weight=0.0;
                if ( abs(i-j)>=SHORT_RANGE_MIN && abs(i-j)<=SHORT_RANGE_MAX )
                {
                    weight=weights[12];	   
                }
                else if ( abs(i-j)>=MED_RANGE_MIN && abs(i-j)<=MED_RANGE_MAX )
                {
                    weight=weights[13];	   
                }
                else if ( abs(i-j)>=LONG_RANGE_MIN )
                {
                    weight=weights[14];	   
                }
                else{
                    weight=0.0;
                }

                if (weight>0.0 && omega_prob_40[i][j]<omg_cut_40 && OMEGA_tmp[i][j]!=180.0)
                {
                    enelist[16]+=weight*spline_omega_40[i][j].cubic_spline( ( OMEGA_tmp[i][j]*raddeg ) );
                    if ( flag_grad )
                    {
                        double dfunc = weight*spline_omega_40[i][j].cubic_dspline( ( OMEGA_tmp[i][j]*raddeg ) );
                        double phi=0.0;
                        atom1 = setv( pn[r1  ][0], pn[r1  ][1], pn[r1  ][2] );
                        atom2 = setv( pn[r1+2][0], pn[r1+2][1], pn[r1+2][2] );
                        atom3 = setv( pn[r2+2][0], pn[r2+2][1], pn[r2+2][2] );
                        atom4 = setv( pn[r2  ][0], pn[r2  ][1], pn[r2  ][2] );
                        
                        point3d scal_f1 = setv( 0.0, 0.0, 0.0 ), scal_f2 = setv( 0.0, 0.0, 0.0 );
                        point3d f1 = setv( 0.0, 0.0, 0.0 ), f2 = setv( 0.0, 0.0, 0.0 );
                        deriv.atom1_dih_deriv( atom1, atom2, atom3, atom4, phi, f1, f2 );
                        scal_f1 = scal( f1, dfunc );
                        scal_f2 = scal( f2, dfunc );
                        atom_derivatives_f1[i][0] = addv( atom_derivatives_f1[i][0], scal_f1 );
                        atom_derivatives_f2[i][0] = addv( atom_derivatives_f2[i][0], scal_f2 ); 

                        scal_f1 = setv( 0.0, 0.0, 0.0 ), scal_f2 = setv( 0.0, 0.0, 0.0 );
                        f1 = setv( 0.0, 0.0, 0.0 ), f2 = setv( 0.0, 0.0, 0.0 );
                        deriv.atom2_dih_deriv( atom1, atom2, atom3, atom4, phi, f1, f2 );
                        scal_f1 = scal( f1, dfunc );
                        scal_f2 = scal( f2, dfunc );
                        atom_derivatives_f1[i][2] = addv( atom_derivatives_f1[i][2], scal_f1 );
                        atom_derivatives_f2[i][2] = addv( atom_derivatives_f2[i][2], scal_f2 ); 
                        
                        scal_f1 = setv( 0.0, 0.0, 0.0 ), scal_f2 = setv( 0.0, 0.0, 0.0 );
                        f1 = setv( 0.0, 0.0, 0.0 ), f2 = setv( 0.0, 0.0, 0.0 );
                        deriv.atom2_dih_deriv( atom4, atom3, atom2, atom1, phi, f1, f2 );
                        scal_f1 = scal( f1, dfunc );
                        scal_f2 = scal( f2, dfunc );
                        atom_derivatives_f1[j][2] = addv( atom_derivatives_f1[j][2], scal_f1 );
                        atom_derivatives_f2[j][2] = addv( atom_derivatives_f2[j][2], scal_f2 ); 

                        scal_f1 = setv( 0.0, 0.0, 0.0 ), scal_f2 = setv( 0.0, 0.0, 0.0 );
                        f1 = setv( 0.0, 0.0, 0.0 ), f2 = setv( 0.0, 0.0, 0.0 );
                        deriv.atom1_dih_deriv( atom4, atom3, atom2, atom1, phi, f1, f2 );
                        scal_f1 = scal( f1, dfunc );
                        scal_f2 = scal( f2, dfunc );
                        atom_derivatives_f1[j][0] = addv( atom_derivatives_f1[j][0], scal_f1 );
                        atom_derivatives_f2[j][0] = addv( atom_derivatives_f2[j][0], scal_f2 ); 
                    } 
                }

                if ( abs(i-j)>=SHORT_RANGE_MIN && abs(i-j)<=SHORT_RANGE_MAX )
                {
                    weight=weights[18];	   
                }
                else if ( abs(i-j)>=MED_RANGE_MIN && abs(i-j)<=MED_RANGE_MAX )
                {
                    weight=weights[19];	   
                }
                else if ( abs(i-j)>=LONG_RANGE_MIN )
                {
                    weight=weights[20];	   
                }
                else{
                    weight=0.0;
                }


                if ( weight>0.0 && phi_prob_40[i][j]<phi_cut_40 )
                {
                    enelist[18]+=weight*spline_phi_40[i][j].cubic_spline( ( PHI_tmp[i][j] ) );

                    if ( flag_grad )
                    {
                        double dfunc=weight*spline_phi_40[i][j].cubic_dspline( ( PHI_tmp[i][j] ) ); 
                        point3d scal_f1=setv(0.0,0.0,0.0),scal_f2=setv(0.0,0.0,0.0);
                        point3d f1_p1, f2_p1, f1_p2, f2_p2, f1_p3, f2_p3;
                        point3d scalf1_p1, scalf2_p1, scalf1_p2, scalf2_p2, scalf1_p3, scalf2_p3;
                        double theta=0.0;
                        atom1=setv(pn[r1  ][0],pn[r1  ][1],pn[r1  ][2]);  //Ca i
                        atom2=setv(pn[r1+2][0],pn[r1+2][1],pn[r1+2][2]);  //Cb i
                        atom3=setv(pn[r2+2][0],pn[r2+2][1],pn[r2+2][2]);  //Cb j

                        deriv.angular_deriv(atom1,atom2,atom3,f1_p1,f2_p1,f1_p2,f2_p2,f1_p3,f2_p3);
                        scalf1_p1=scal(f1_p1,dfunc);
                        scalf2_p1=scal(f2_p1,dfunc);
                        atom_derivatives_f1[i][0]=addv(atom_derivatives_f1[i][0],scalf1_p1);
                        atom_derivatives_f2[i][0]=addv(atom_derivatives_f2[i][0],scalf2_p1);
                        scalf1_p2=scal(f1_p2,dfunc);
                        scalf2_p2=scal(f2_p2,dfunc);
                        atom_derivatives_f1[i][2]=addv(atom_derivatives_f1[i][2],scalf1_p2);
                        atom_derivatives_f2[i][2]=addv(atom_derivatives_f2[i][2],scalf2_p2);
                        scalf1_p3=scal(f1_p3,dfunc);
                        scalf2_p3=scal(f2_p3,dfunc);
                        atom_derivatives_f1[j][2]=addv(atom_derivatives_f1[j][2],scalf1_p3);
                        atom_derivatives_f2[j][2]=addv(atom_derivatives_f2[j][2],scalf2_p3);
                    }
                }
                
                if ( weight>0.0 && phi_prob_40[j][i]<phi_cut_40 )
                {
                    enelist[18] += weight*spline_phi_40[j][i].cubic_spline( ( PHI_tmp[j][i] ) );

                    if ( flag_grad )
                    {
                        double dfunc = weight * spline_phi_40[j][i].cubic_dspline( ( PHI_tmp[j][i] ) ); 
                        point3d scal_f1=setv(0.0,0.0,0.0),scal_f2=setv(0.0,0.0,0.0);
                        point3d f1_p1, f2_p1, f1_p2, f2_p2, f1_p3, f2_p3;
                        point3d scalf1_p1, scalf2_p1, scalf1_p2, scalf2_p2, scalf1_p3, scalf2_p3;
                        double theta=0.0;
                        atom1 = setv( pn[r2  ][0], pn[r2  ][1], pn[r2  ][2] );  //Ca j
                        atom2 = setv( pn[r2+2][0], pn[r2+2][1], pn[r2+2][2] );  //Cb j
                        atom3 = setv( pn[r1+2][0], pn[r1+2][1], pn[r1+2][2] );  //Cb i
                        deriv.angular_deriv( atom1, atom2, atom3, f1_p1, f2_p1, 
                                             f1_p2, f2_p2, f1_p3, f2_p3 );
                        scalf1_p1 = scal( f1_p1, dfunc );
                        scalf2_p1 = scal( f2_p1, dfunc);
                        atom_derivatives_f1[j][0] = addv( atom_derivatives_f1[j][0], scalf1_p1 );
                        atom_derivatives_f2[j][0] = addv( atom_derivatives_f2[j][0], scalf2_p1 );
                        scalf1_p2 = scal( f1_p2, dfunc );
                        scalf2_p2 = scal( f2_p2, dfunc );
                        atom_derivatives_f1[j][2] = addv( atom_derivatives_f1[j][2], scalf1_p2 );
                        atom_derivatives_f2[j][2] = addv( atom_derivatives_f2[j][2], scalf2_p2 );
                        scalf1_p3 = scal( f1_p3, dfunc );
                        scalf2_p3 = scal( f2_p3, dfunc );
                        atom_derivatives_f1[i][2] = addv( atom_derivatives_f1[i][2], scalf1_p3 );
                        atom_derivatives_f2[i][2] = addv( atom_derivatives_f2[i][2], scalf2_p3 );              
                    }
                }

                if ( abs(i-j)>=SHORT_RANGE_MIN && abs(i-j)<=SHORT_RANGE_MAX )
                {
                    weight=weights[15];	   
                }
                else if ( abs(i-j)>=MED_RANGE_MIN && abs(i-j)<=MED_RANGE_MAX )
                {
                    weight=weights[16];	   
                }
                else if ( abs(i-j)>=LONG_RANGE_MIN )
                {
                    weight=weights[17];	   
                }
                else{
                    weight=0.0;
                }

                if( weight>0.0 && theta_prob_40[i][j]<theta_cut_40 && THETA_tmp[i][j]!=180.0 )
                {
                    enelist[17]+=weight*spline_theta_40[i][j].cubic_spline( ( THETA_tmp[i][j]*raddeg ) );

                    if(flag_grad)
                    {
                        double dfunc = weight*spline_theta_40[i][j].cubic_dspline( ( THETA_tmp[i][j]*raddeg ) );
                        double phi=0.0;
                        point3d scal_f1=setv(0.0,0.0,0.0),scal_f2=setv(0.0,0.0,0.0);
                        point3d f1, f2;
                        atom1=setv(pn[r1+1][0],pn[r1+1][1],pn[r1+1][2]);  //N  i
                        atom2=setv(pn[r1  ][0],pn[r1  ][1],pn[r1  ][2]);  //Ca i
                        atom3=setv(pn[r1+2][0],pn[r1+2][1],pn[r1+2][2]);  //Cb i
                        atom4=setv(pn[r2+2][0],pn[r2+2][1],pn[r2+2][2]);  //Cb j

                        deriv.atom1_dih_deriv( atom1, atom2, atom3, atom4, phi, f1, f2 );
                        scal_f1=scal(f1,dfunc);
                        scal_f2=scal(f2,dfunc);
                        atom_derivatives_f1[i][1]=addv(atom_derivatives_f1[i][1],scal_f1);
                        atom_derivatives_f2[i][1]=addv(atom_derivatives_f2[i][1],scal_f2); 
   
                        deriv.atom2_dih_deriv( atom1, atom2, atom3, atom4, phi, f1, f2 );
                        scal_f1=scal(f1,dfunc);
                        scal_f2=scal(f2,dfunc);
                        atom_derivatives_f1[i][0]=addv(atom_derivatives_f1[i][0],scal_f1);
                        atom_derivatives_f2[i][0]=addv(atom_derivatives_f2[i][0],scal_f2); 

                        deriv.atom2_dih_deriv( atom4, atom3, atom2, atom1, phi, f1, f2 );
                        scal_f1=scal(f1,dfunc);
                        scal_f2=scal(f2,dfunc);
                        atom_derivatives_f1[i][2]=addv(atom_derivatives_f1[i][2],scal_f1);
                        atom_derivatives_f2[i][2]=addv(atom_derivatives_f2[i][2],scal_f2); 

                        deriv.atom1_dih_deriv( atom4, atom3, atom2, atom1, phi, f1, f2 );
                        scal_f1=scal(f1,dfunc);
                        scal_f2=scal(f2,dfunc);
                        atom_derivatives_f1[j][2]=addv(atom_derivatives_f1[j][2],scal_f1);
                        atom_derivatives_f2[j][2]=addv(atom_derivatives_f2[j][2],scal_f2);
                    } 
                }
                if( weight>0.0 && theta_prob_40[j][i]<theta_cut_40 && THETA_tmp[j][i]!=180.0 )
                {
                    enelist[17]+=weight*spline_theta_40[j][i].cubic_spline( ( THETA_tmp[j][i]*raddeg ) );
                    if(flag_grad)
                    {
                        double dfunc = weight*spline_theta_40[j][i].cubic_dspline( ( THETA_tmp[j][i]*raddeg ) );
                        double phi=0.0;
                        point3d scal_f1=setv(0.0,0.0,0.0),scal_f2=setv(0.0,0.0,0.0);
                        point3d f1, f2;
                        atom1=setv(pn[r2+1][0],pn[r2+1][1],pn[r2+1][2]);  //N  j
                        atom2=setv(pn[r2  ][0],pn[r2  ][1],pn[r2  ][2]);  //Ca j
                        atom3=setv(pn[r2+2][0],pn[r2+2][1],pn[r2+2][2]);  //Cb j
                        atom4=setv(pn[r1+2][0],pn[r1+2][1],pn[r1+2][2]);  //Cb i
                        deriv.atom1_dih_deriv( atom1, atom2, atom3, atom4, phi, f1, f2 );
                        scal_f1=scal(f1,dfunc);
                        scal_f2=scal(f2,dfunc);
                        atom_derivatives_f1[j][1]=addv(atom_derivatives_f1[j][1],scal_f1);
                        atom_derivatives_f2[j][1]=addv(atom_derivatives_f2[j][1],scal_f2); 
   
                        deriv.atom2_dih_deriv( atom1, atom2, atom3, atom4, phi, f1, f2 );
                        scal_f1=scal(f1,dfunc);
                        scal_f2=scal(f2,dfunc);
                        atom_derivatives_f1[j][0]=addv(atom_derivatives_f1[j][0],scal_f1);
                        atom_derivatives_f2[j][0]=addv(atom_derivatives_f2[j][0],scal_f2); 

                        deriv.atom2_dih_deriv( atom4, atom3, atom2, atom1, phi, f1, f2 );
                        scal_f1=scal(f1,dfunc);
                        scal_f2=scal(f2,dfunc);
                        atom_derivatives_f1[j][2]=addv(atom_derivatives_f1[j][2],scal_f1);
                        atom_derivatives_f2[j][2]=addv(atom_derivatives_f2[j][2],scal_f2); 

                        deriv.atom1_dih_deriv( atom4, atom3, atom2, atom1, phi, f1, f2 );
                        scal_f1=scal(f1,dfunc);
                        scal_f2=scal(f2,dfunc);
                        atom_derivatives_f1[i][2]=addv(atom_derivatives_f1[i][2],scal_f1);
                        atom_derivatives_f2[i][2]=addv(atom_derivatives_f2[i][2],scal_f2); 
                    } 
                }
            }
        }
    }
    
    /* calculate the rest of the terms */
    for(i=0;i<numseq;i++)
    {
        r1=i*number_atoms;
        
	for(j=i+1;j<numseq;j++)
        {
            r2=j*number_atoms;
            indm=i*numseq+j;
            idx=(j-1)+(2*numseq-3-i)*i/2;

            /**** CB-CB specific energy ****/
            tdist2_CB=C_N_tmp[j][i];
            tdist_CB=tdist2_CB*tdist2_CB;

            double weight=0.0;
            if ( abs(i-j)>=SHORT_RANGE_MIN && abs(i-j)<=SHORT_RANGE_MAX )
            {
                weight=weights[0];
            }
            else if ( abs(i-j)>=MED_RANGE_MIN && abs(i-j)<=MED_RANGE_MAX )
            {
                weight=weights[1];
            }
            else if ( abs(i-j)>=LONG_RANGE_MIN )
            {
                weight=weights[2];
            }
            else
            {
                weight=0.0;
            }
            
            if ( tdist2_CB<39.0 && weight>0.0 && ( cb_dist_prob_new_40[i][j]<dist_cut_cb_40 ) )
            {
                enelist[12] += weight * spline_dist_new_40[i][j].cubic_spline( tdist2_CB );
                if(flag_grad)
                {
                    double dfunc = weight * spline_dist_new_40[i][j].cubic_dspline( tdist2_CB );
                    point3d f1=setv(0.0,0.0,0.0),f2=setv(0.0,0.0,0.0);
                    point3d scal_f1_res1=setv(0.0,0.0,0.0),scal_f2_res1=setv(0.0,0.0,0.0);
                    point3d scal_f1_res2=setv(0.0,0.0,0.0),scal_f2_res2=setv(0.0,0.0,0.0);
		    atom1.x=decstr[i].ptn.x,atom1.y=decstr[i].ptn.y,atom1.z=decstr[i].ptn.z;
                    atom2.x=decstr[j].ptn.x,atom2.y=decstr[j].ptn.y,atom2.z=decstr[j].ptn.z;
                    deriv.dist_deriv( atom1, atom2, tdist2_CB, f1, f2 );
                    scal_f1_res1=scal(f1,dfunc);
                    scal_f2_res1=scal(f2,dfunc);
                    scal_f1_res2=scal(f1,-dfunc);
                    scal_f2_res2=scal(f2,-dfunc);
                    atom_derivatives_f1[i][2]=addv(atom_derivatives_f1[i][2],scal_f1_res1);
                    atom_derivatives_f2[i][2]=addv(atom_derivatives_f2[i][2],scal_f2_res1);
                    atom_derivatives_f1[j][2]=addv(atom_derivatives_f1[j][2],scal_f1_res2);
                    atom_derivatives_f2[j][2]=addv(atom_derivatives_f2[j][2],scal_f2_res2);
                }
            }

            /**** CA-CA specific energy ****/
            
            tdist2_CA=C_N_tmp[i][j];
            tdist_CA=tdist2_CA*tdist2_CA;

            weight=0.0;
            
            if ( abs(i-j)>=SHORT_RANGE_MIN && abs(i-j)<=SHORT_RANGE_MAX )
            {
                weight=weights[3];
            }
            else if ( abs(i-j)>=MED_RANGE_MIN && abs(i-j)<=MED_RANGE_MAX )
            {
                weight=weights[4];
            }
            else if ( abs(i-j)>=LONG_RANGE_MIN )
            {
                weight=weights[5];
            }
            else
            {
                weight=0.0;
            }
            
            
            if ( tdist2_CA<39.0 && weight>0.0 && ( ca_dist_prob_40[i][j]<dist_cut_ca_40 ) )
            {
                enelist[8] += weight * spline_dist_ca_40[i][j].cubic_spline( tdist2_CA );
                if(flag_grad)
                {
                    double dfunc = weight * spline_dist_ca_40[i][j].cubic_dspline( tdist2_CA );
                    point3d f1=setv(0.0,0.0,0.0),f2=setv(0.0,0.0,0.0);
                    point3d scal_f1_res1=setv(0.0,0.0,0.0),scal_f2_res1=setv(0.0,0.0,0.0);
                    point3d scal_f1_res2=setv(0.0,0.0,0.0),scal_f2_res2=setv(0.0,0.0,0.0);
                    atom1.x=decstr[i].x,atom1.y=decstr[i].y,atom1.z=decstr[i].z;
                    atom2.x=decstr[j].x,atom2.y=decstr[j].y,atom2.z=decstr[j].z;
                    deriv.dist_deriv( atom1, atom2, tdist2_CA, f1, f2 );
                    scal_f1_res1=scal(f1,dfunc);
                    scal_f2_res1=scal(f2,dfunc);
                    scal_f1_res2=scal(f1,-dfunc);
                    scal_f2_res2=scal(f2,-dfunc);
                    atom_derivatives_f1[i][0]=addv(atom_derivatives_f1[i][0],scal_f1_res1);
                    atom_derivatives_f2[i][0]=addv(atom_derivatives_f2[i][0],scal_f2_res1);
                    atom_derivatives_f1[j][0]=addv(atom_derivatives_f1[j][0],scal_f1_res2);
                    atom_derivatives_f2[j][0]=addv(atom_derivatives_f2[j][0],scal_f2_res2);
                }
            }

            tp.x=pn[r1+1][0]-pn[r2+1][0];
            tp.y=pn[r1+1][1]-pn[r2+1][1];
            tp.z=pn[r1+1][2]-pn[r2+1][2];
            double tdist2_P=sqrt(tp.x*tp.x+tp.y*tp.y+tp.z*tp.z);

            weight=0.0;
            
            if ( abs(i-j)>=SHORT_RANGE_MIN && abs(i-j)<=SHORT_RANGE_MAX )
            {
                weight=weights[30];
            }
            else if ( abs(i-j)>=MED_RANGE_MIN && abs(i-j)<=MED_RANGE_MAX )
            {
                weight=weights[31];
            }
            else if ( abs(i-j)>=LONG_RANGE_MIN )
            {
                weight=weights[32];
            }
            else
            {
                weight=0.0;
            }

            if ( tdist2_P<39.0 && weight>0.0 && ( p_dist_prob_40[i][j]<dist_cut_p_40 ) )
            {
                enelist[8] += weight * spline_dist_p_40[i][j].cubic_spline( tdist2_P );
                if(flag_grad)
                {
                    double dfunc = weight * spline_dist_p_40[i][j].cubic_dspline( tdist2_P );
                    point3d f1=setv(0.0,0.0,0.0),f2=setv(0.0,0.0,0.0);
                    point3d scal_f1_res1=setv(0.0,0.0,0.0),scal_f2_res1=setv(0.0,0.0,0.0);
                    point3d scal_f1_res2=setv(0.0,0.0,0.0),scal_f2_res2=setv(0.0,0.0,0.0);
                    atom1.x=decstr[i].ptp.x,atom1.y=decstr[i].ptp.y,atom1.z=decstr[i].ptp.z;
                    atom2.x=decstr[j].ptp.x,atom2.y=decstr[j].ptp.y,atom2.z=decstr[j].ptp.z;
                    deriv.dist_deriv( atom1, atom2, tdist2_P, f1, f2 );
                    scal_f1_res1=scal(f1,dfunc);
                    scal_f2_res1=scal(f2,dfunc);
                    scal_f1_res2=scal(f1,-dfunc);
                    scal_f2_res2=scal(f2,-dfunc);
                    atom_derivatives_f1[i][1]=addv(atom_derivatives_f1[i][1],scal_f1_res1);
                    atom_derivatives_f2[i][1]=addv(atom_derivatives_f2[i][1],scal_f2_res1);
                    atom_derivatives_f1[j][1]=addv(atom_derivatives_f1[j][1],scal_f1_res2);
                    atom_derivatives_f2[j][1]=addv(atom_derivatives_f2[j][1],scal_f2_res2);
                }
            }
        } //j
    } //i
}

void Energy::torsion_energy(point3f *decstr, double *enelist ){
    
    if(weights[29]==0.0) return;
    
    for(int i=1;i<numseq-1;i++)
    {
        point3d atom1,atom2,atom3,atom4;
        double phi=decstr[i].theta*raddeg;       
        double psi=decstr[i+1].eta*raddeg;    
        enelist[14]+=weights[29]*spline_eta_bb[i].cubic_spline( ( phi ) );
        enelist[14]+=weights[29]*spline_theta_bb[i].cubic_spline( ( psi ) );
    }
}

void Energy::torsion_energy(point3f *decstr, double *enelist,double *dE_dvars ){
   
    if(weights[29]==0.0) return;
    
    for(int i=1;i<numseq-1;i++)
    {
        double phi=decstr[i].theta*raddeg;       
        double psi=decstr[i+1].eta*raddeg;       
        enelist[14]+=weights[29]*spline_eta_bb[i].cubic_spline( ( phi ) );
        enelist[14]+=weights[29]*spline_theta_bb[i].cubic_spline( ( psi ) );
        if(flag_grad)
        {
            dE_dvars[(i-1)*2] +=weights[29]*spline_eta_bb[i].cubic_dspline( ( phi ) );
            dE_dvars[(i-1)*2+1]+=weights[29]*spline_theta_bb[i].cubic_dspline( ( psi ) );
        }
    }
}

double Energy::calcrmsdenergy( point3f *decstr, double *vars, double *dE_dvars, bool calc_gradients )
{
    double trmsd=0.0; // total energy E_total = sum of enelist
    int i,j,ii;
    for(i=0;i<20;i++) enelist[i]=0;
    geo.apply_torsions( decstr, vars, numseq );
    flag_grad=calc_gradients;
    
    // Reset gradients with respect to each variable (torsions)
    if(flag_grad)
    {
        for(i=0;i<(numseq-1)*2;i++)
        {
            dE_dvars[i]=0.0;
        }

        // Initialize f1/f2 as zero
        for(i=0;i<numseq;i++)
        {
            for(ii=0;ii<number_atoms;ii++)
            {
                atom_derivatives_f1[i][ii].x=0.0;
                atom_derivatives_f1[i][ii].y=0.0;
                atom_derivatives_f1[i][ii].z=0.0;
                atom_derivatives_f2[i][ii].x=0.0;
                atom_derivatives_f2[i][ii].y=0.0;
                atom_derivatives_f2[i][ii].z=0.0;
            }
        }
    }

    // Calculate most atomic energy terms and their gradients
    calcallenergy(decstr, enelist);

    // Calculate torsion energy from predicted phi/psi
    torsion_energy(decstr, enelist,dE_dvars);
    
    for(i=0;i<20;i++) trmsd+=enelist[i];
    
    // Sum all f1/f2 contributions from atoms moved by each given torsion angle. This can be done
    // recurrently but the time savings is probably minimal
    if ( flag_grad )
    {
        for ( i=1; i<numseq; i++ )
        {
            decstr[i].f1_theta.x=0.0,decstr[i].f1_theta.y=0.0,decstr[i].f1_theta.z=0.0;
            decstr[i].f2_theta.x=0.0,decstr[i].f2_theta.y=0.0,decstr[i].f2_theta.z=0.0;
            decstr[i].f1_eta.x=0.0,decstr[i].f1_eta.y=0.0,decstr[i].f1_eta.z=0.0;
            decstr[i].f2_eta.x=0.0,decstr[i].f2_eta.y=0.0,decstr[i].f2_eta.z=0.0;
            
            // Add intra-residue f1/f2 contributions
            for ( ii=0; ii<number_atoms; ii++ )
            {
                if(ii==0 || ii==1) continue;
                decstr[i].f1_eta=addv(decstr[i].f1_eta,atom_derivatives_f1[i][ii]);
                decstr[i].f2_eta=addv(decstr[i].f2_eta,atom_derivatives_f2[i][ii]);
            }

            // Add inter-residue f1/f2 contributions
            for ( j=i+1; j<numseq; j++ )
            {
                for ( ii=0; ii<number_atoms; ii++ )
                {
                    decstr[i].f1_eta=addv(decstr[i].f1_eta,atom_derivatives_f1[j][ii]);
                    decstr[i].f2_eta=addv(decstr[i].f2_eta,atom_derivatives_f2[j][ii]);
                }
            }

            // Set psi f1/f2 equal to phi since all atoms moved by phi are moved by psi
            decstr[i].f1_theta.x=decstr[i].f1_eta.x, decstr[i].f1_theta.y=decstr[i].f1_eta.y;
            decstr[i].f1_theta.z=decstr[i].f1_eta.z;
            decstr[i].f2_theta.x=decstr[i].f2_eta.x, decstr[i].f2_theta.y=decstr[i].f2_eta.y, 
            decstr[i].f2_theta.z=decstr[i].f2_eta.z;

            // Add in contributions from the backbone C4 to psi
            decstr[i].f1_theta=addv(decstr[i].f1_theta,atom_derivatives_f1[i][0]);
            decstr[i].f2_theta=addv(decstr[i].f2_theta,atom_derivatives_f2[i][0]);
        }
    }
    
    // Calculate gradients with respect to each torsion angle. 
    if( flag_grad )
    { 
        for( i=1; i<numseq; i++ )
        {
            double dtheta=0.0,deta=0.0;
            point3d axis;
            point3d end_pos;
            end_pos.x = decstr[i].ptp.x;
            end_pos.y = decstr[i].ptp.y;
            end_pos.z = decstr[i].ptp.z;
            axis.x = decstr[i].ptp.x-decstr[i-1].x; 
            axis.y = decstr[i].ptp.y-decstr[i-1].y; 
            axis.z = decstr[i].ptp.z-decstr[i-1].z;
            axis = unit( axis );
            dtheta -= raddeg * ( dotv( axis, decstr[i].f1_theta ) + 
                               dotv( prod( axis, end_pos ), decstr[i].f2_theta ) );
            dE_dvars[(i-1)*2] = dtheta;

            end_pos.x = decstr[i].x;
            end_pos.y = decstr[i].y;
            end_pos.z = decstr[i].z;
            axis.x = decstr[i].x-decstr[i].ptp.x;
            axis.y = decstr[i].y-decstr[i].ptp.y;
            axis.z = decstr[i].z-decstr[i].ptp.z;
            axis = unit( axis );
            deta -= raddeg * ( dotv( axis, decstr[i].f1_eta ) + 
                               dotv( prod( axis, end_pos ), decstr[i].f2_eta ) );
            dE_dvars[(i-1)*2+1] = deta;
        }
    }
    
    return trmsd;
}

double Energy::calcrmsdenergy( point3f *decstr, double *vars )
{
    double trmsd=0.0; // total energy E_total = sum of enelist
    int i,j,ii;
    for(i=0;i<20;i++) enelist[i]=0;
    geo.apply_torsions( decstr, vars, numseq );
    flag_grad=false; // we don't need to calculate gradients here
    
    calcallenergy(decstr, enelist);
    torsion_energy(decstr, enelist);
    for(i=0;i<20;i++) trmsd+=enelist[i];
    return trmsd;
}

#endif
