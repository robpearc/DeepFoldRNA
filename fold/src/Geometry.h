/*******************************************************************************************
**  
**  Functions for converting between torsion/cartesian space and for applying gradients to 
**  the given degrees of freedom
**
**  Please report bugs and questions to robpearc@umich.edu
**
*******************************************************************************************/

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <cstdio>
#include <cassert>
#include <vector>
#include <algorithm>
#include "CommonParameters.h"
#include "Operations.h"

using namespace std;

class Geometry{
public:
    //void loadsgpos2( const char *filename, int ndim );

    void randomdih( double *phi, double *psi );

    bool tor2coord( point3f *decstr, int seqnum );
    bool tor2coordsg( point3f *decstr, int seqnum );
    //bool tor2coordsg_all( point3f *decstr, int seqnum );
    bool tor2coord_base( point3f *decstr, int seqnum );

    bool setinitdecoy( point3f *decstr, int numseq );
    bool setinitdecoyfromfile( point3f *decstr, int numseq, const char *theta_file, 
                               const char *eta_file, const char *chi_file );

    bool apply_torsions( point3f *decstr, double *vars, int numseq );

    virtual ~Geometry();

private:
    //double ****sgposdat;
    bool flagsgpos2;

};

Geometry::~Geometry()
{
}

bool Geometry::tor2coord( point3f *decstr, int seqnum )
{
    int i;
    bool flagts;
    bool flagwhole=true;
    double len_p_c=3.8761;
    double len_c_p=3.8705;
    double len_c_n=3.3912;
    double ang_p_c_p=104.8433;
    double ang_c_p_c=105.6791;

    point3s pfirst;
    decstr[0].ptp.x=double(-0.42318);
    decstr[0].ptp.y=double(3.7758);
    decstr[0].ptp.z=0.0;
    decstr[0].x=0;//decstr[0].len_p_c;
    decstr[0].y=0;
    decstr[0].z=0;
    decstr[1].ptp.x=0-decstr[0].len_c_p*cos(decstr[0].ang_p_c_p*raddeg);
    decstr[1].ptp.y=decstr[0].len_c_p*sin(decstr[0].ang_p_c_p*raddeg);
    decstr[1].ptp.z=0;
    flagts=tor2pos22(decstr[0].ptp.x,decstr[0].ptp.y,decstr[0].ptp.z,
                     decstr[0].x,decstr[0].y,decstr[0].z,
                     decstr[1].ptp.x,decstr[i].ptp.y,decstr[i].ptp.z,
                     decstr[1].theta*raddeg,decstr[1].len_p_c,
                     decstr[1].ang_c_p_c*raddeg,&pfirst.x,&pfirst.y,&pfirst.z);
    if(!flagts)
    {
        flagwhole=false;
    }
    decstr[1].x=pfirst.x;
    decstr[1].y=pfirst.y;
    decstr[1].z=pfirst.z;
   
    for(i=2;i<seqnum;i++)
    {
        point3s pt,pp;
        flagts=tor2pos22(decstr[i-2].x,decstr[i-2].y,decstr[i-2].z,
                         decstr[i-1].ptp.x,decstr[i-1].ptp.y,decstr[i-1].ptp.z,
                         decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
                         decstr[i-1].eta*raddeg,decstr[i].len_c_p,
                         decstr[i].ang_p_c_p*raddeg,&pp.x,&pp.y,&pp.z);
        if(!flagts)
        {
            flagwhole=false;
        }
        decstr[i].ptp.x=pp.x;
        decstr[i].ptp.y=pp.y;
        decstr[i].ptp.z=pp.z;
        flagts=tor2pos22(decstr[i-1].ptp.x,decstr[i-1].ptp.y,decstr[i-1].ptp.z,
                         decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
                         decstr[i].ptp.x,decstr[i].ptp.y,decstr[i].ptp.z,
                         decstr[i].theta*raddeg,decstr[i].len_p_c,
                         decstr[i].ang_c_p_c*raddeg,&pt.x,&pt.y,&pt.z);
        if(!flagts)
        {
            flagwhole=false;
        }
        decstr[i].x=pt.x;
        decstr[i].y=pt.y;
        decstr[i].z=pt.z;
    }
    tor2coordsg( decstr, seqnum );
    return flagwhole;
}

bool Geometry::tor2coordsg( point3f *decstr, int seqnum )
{
    int i;
    int ntind;

    bool flagts;
    bool flagwhole=true;

    for(i=0;i<seqnum-1;i++)
    {
        point3s pn;
        ntind=getntid(decstr[i].nt);
        if(ntind>3) ntind=3;
        double tor=torN[ntind]*raddeg;
        double len=distN[ntind];
        double ang=angN[ntind]*raddeg;
        
        flagts=tor2pos22(decstr[i].ptp.x,decstr[i].ptp.y,decstr[i].ptp.z,
                         decstr[i+1].ptp.x,decstr[i+1].ptp.y,decstr[i+1].ptp.z,
                         decstr[i].x,decstr[i].y,decstr[i].z,
                         tor,len,ang,&pn.x,&pn.y,&pn.z);
        if(!flagts)
        {
            flagwhole=false;
        }
        decstr[i].ptn.x=pn.x;
        decstr[i].ptn.y=pn.y;
        decstr[i].ptn.z=pn.z;
    }  
    
    point3s pvirt,pn;   
    ntind=getntid(decstr[seqnum-1].nt);
    if(ntind>3) ntind=3;
    double tor=torN[ntind]*raddeg;
    double len=distN[ntind];
    double ang=angN[ntind]*raddeg;
    
    flagts=tor2pos22(decstr[seqnum-2].x,decstr[seqnum-2].y,decstr[seqnum-2].z,
                     decstr[seqnum-1].ptp.x,decstr[seqnum-1].ptp.y,decstr[seqnum-1].ptp.z,
                     decstr[seqnum-1].x,decstr[seqnum-1].y,decstr[seqnum-1].z,
                     decstr[seqnum-1].eta*raddeg,decstr[seqnum-1].len_c_p,
                     decstr[seqnum-1].ang_p_c_p*raddeg,&pvirt.x,&pvirt.y,&pvirt.z);
    
    if(!flagts)
    {
        flagwhole=false;
    }
    
    flagts=tor2pos22(decstr[seqnum-1].ptp.x,decstr[seqnum-1].ptp.y,decstr[seqnum-1].ptp.z,
                     pvirt.x,pvirt.y,pvirt.z,
                     decstr[seqnum-1].x,decstr[seqnum-1].y,decstr[seqnum-1].z,
                     tor,len,ang,&pn.x,&pn.y,&pn.z);
    if(!flagts)
    {
        flagwhole=false;
    }
    decstr[seqnum-1].ptn.x=pn.x;
    decstr[seqnum-1].ptn.y=pn.y;
    decstr[seqnum-1].ptn.z=pn.z;
 
    return flagwhole;
}

bool Geometry::tor2coord_base( point3f *decstr, int seqnum )
{
    int i;
    int ntind;

    bool flagts;
    bool flagwhole=true;

    point3s pvirt;   
    flagts=tor2pos22(decstr[seqnum-2].x,decstr[seqnum-2].y,decstr[seqnum-2].z,
                     decstr[seqnum-1].ptp.x,decstr[seqnum-1].ptp.y,decstr[seqnum-1].ptp.z,
                     decstr[seqnum-1].x,decstr[seqnum-1].y,decstr[seqnum-1].z,
                     decstr[seqnum-1].eta*raddeg,decstr[seqnum-1].len_c_p,
                     decstr[seqnum-1].ang_p_c_p*raddeg,&pvirt.x,&pvirt.y,&pvirt.z);
    
    if(!flagts)
    {
        flagwhole=false;
    }

    for(i=0;i<seqnum;i++)
    {
        point3s po4,pc1,pc2,pc4;
        ntind=getntid(decstr[i].nt);
        if(ntind>3) ntind=3;
        if(ntind==0 || ntind==2) continue;
        
        double tor=torO4[ntind]*raddeg;
        double len=distO4[ntind];
        double ang=angO4[ntind]*raddeg;
        
        if(i<seqnum-1)
        {
            flagts=tor2pos22(decstr[i].ptp.x,decstr[i].ptp.y,decstr[i].ptp.z,
                             decstr[i+1].ptp.x,decstr[i+1].ptp.y,decstr[i+1].ptp.z,
                             decstr[i].x,decstr[i].y,decstr[i].z,
                             tor,len,ang,&po4.x,&po4.y,&po4.z);
            if(!flagts)
            {
                flagwhole=false;
            }
        }
        else
        {
            flagts=tor2pos22(decstr[i].ptp.x,decstr[i].ptp.y,decstr[i].ptp.z,
                             pvirt.x,pvirt.y,pvirt.z,
                             decstr[i].x,decstr[i].y,decstr[i].z,
                             tor,len,ang,&po4.x,&po4.y,&po4.z);
            if(!flagts)
            {
                flagwhole=false;
            }
        }
        
        tor=torC1[ntind]*raddeg;
        len=distC1[ntind];
        ang=angC1[ntind]*raddeg;
        flagts=tor2pos22(decstr[i].ptp.x,decstr[i].ptp.y,decstr[i].ptp.z,
                         decstr[i].x,decstr[i].y,decstr[i].z,
                         po4.x,po4.y,po4.z,
                         tor,len,ang,&pc1.x,&pc1.y,&pc1.z);
  
        len=distChi[ntind];
        ang=angChi[ntind]*raddeg;
        flagts=tor2pos22(po4.x,po4.y,po4.z,
                         pc1.x,pc1.y,pc1.z,
                         decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
                         decstr[i].chi,len,ang,&pc2.x,&pc2.y,&pc2.z);
        
        decstr[i].ptc2.x=pc2.x;
        decstr[i].ptc2.y=pc2.y;
        decstr[i].ptc2.z=pc2.z;

        tor=torAfterChi[ntind]*raddeg;
        len=distAfterChi[ntind];
        ang=angAfterChi[ntind]*raddeg;
        flagts=tor2pos22(pc1.x,pc1.y,pc1.z,
                         decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
                         pc2.x,pc2.y,pc2.z,
                         tor,len,ang,&pc4.x,&pc4.y,&pc4.z);

        decstr[i].ptc4.x=pc4.x;
        decstr[i].ptc4.y=pc4.y;
        decstr[i].ptc4.z=pc4.z;
    }

    for(i=0;i<seqnum;i++)
    {
        point3s po4,pc1,pc2,pc4,pc6;
        ntind=getntid(decstr[i].nt);
        if(ntind>3) ntind=3;
        if(ntind==1 || ntind==3) continue;
        
        double tor=torO4[ntind]*raddeg;
        double len=distO4[ntind];
        double ang=angO4[ntind]*raddeg;

        if(i<seqnum-1)
        {
            flagts=tor2pos22(decstr[i].ptp.x,decstr[i].ptp.y,decstr[i].ptp.z,
                             decstr[i+1].ptp.x,decstr[i+1].ptp.y,decstr[i+1].ptp.z,
                             decstr[i].x,decstr[i].y,decstr[i].z,
                             tor,len,ang,&po4.x,&po4.y,&po4.z);
            if(!flagts)
            {
                flagwhole=false;
            }
        }
        else
        {
            flagts=tor2pos22(decstr[i].ptp.x,decstr[i].ptp.y,decstr[i].ptp.z,
                             pvirt.x,pvirt.y,pvirt.z,
                             decstr[i].x,decstr[i].y,decstr[i].z,
                             tor,len,ang,&po4.x,&po4.y,&po4.z);
            if(!flagts)
            {
                flagwhole=false;
            }
        }
        
        tor=torC1[ntind]*raddeg;
        len=distC1[ntind];
        ang=angC1[ntind]*raddeg;
        flagts=tor2pos22(decstr[i].ptp.x,decstr[i].ptp.y,decstr[i].ptp.z,
                         decstr[i].x,decstr[i].y,decstr[i].z,
                         po4.x,po4.y,po4.z,
                         tor,len,ang,&pc1.x,&pc1.y,&pc1.z);
  
        len=distChi[ntind];
        ang=angChi[ntind]*raddeg;
        flagts=tor2pos22(po4.x,po4.y,po4.z,
                         pc1.x,pc1.y,pc1.z,
                         decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
                         decstr[i].chi,len,ang,&pc4.x,&pc4.y,&pc4.z);
        
        decstr[i].ptc4.x=pc4.x;
        decstr[i].ptc4.y=pc4.y;
        decstr[i].ptc4.z=pc4.z;

        tor=torAfterChi[ntind]*raddeg;
        len=distAfterChi[ntind];
        ang=angAfterChi[ntind]*raddeg;
        flagts=tor2pos22(pc1.x,pc1.y,pc1.z,
                         decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
                         pc4.x,pc4.y,pc4.z,
                         tor,len,ang,&pc2.x,&pc2.y,&pc2.z);

        decstr[i].ptc2.x=pc2.x;
        decstr[i].ptc2.y=pc2.y;
        decstr[i].ptc2.z=pc2.z;
        
        tor=torC6[ntind]*raddeg;
        len=distC6[ntind];
        ang=angC6[ntind]*raddeg;
        flagts=tor2pos22(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
                         pc4.x,pc4.y,pc4.z,
                         pc2.x,pc2.y,pc2.z,
                         tor,len,ang,&pc6.x,&pc6.y,&pc6.z);
        
        decstr[i].ptc6.x=pc6.x;
        decstr[i].ptc6.y=pc6.y;
        decstr[i].ptc6.z=pc6.z;
    }
}

bool Geometry::setinitdecoyfromfile( point3f *decstr, int numseq, const char *theta_file, const char *eta_file, const char *chi_file )
{
    double len_p_c=3.8761;
    double len_c_p=3.8705;
    double len_c_n=3.3912;
    double ang_p_c_p=104.8433;
    double ang_c_p_c=105.6791;

    for ( int i=0; i<numseq; i++ )
    {
        decstr[i].len_p_c=len_p_c;
        decstr[i].len_c_p=len_c_p;
        decstr[i].len_c_n=len_c_n;
        decstr[i].ang_p_c_p=ang_p_c_p;
        decstr[i].ang_c_p_c=ang_c_p_c;
    }

    FILE *file;
    if ( ( file=fopen(eta_file,"rt") )!=NULL )
    {
        char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
        while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,file))
        {
            int res;
            double val=0.0;
            sscanf(line,"%d %lf",&res,&val);
	    //if(res==0) continue;
            if(val<0) val+=360.0;
            decstr[res].eta=val;
        }
    }
    fclose(file);
    
    decstr[0].theta=0.0;
    if ( ( file=fopen(theta_file,"rt") )!=NULL )
    {
        char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
        while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,file))
        {
            int res;
            double val=0.0;
            sscanf(line,"%d %lf",&res,&val);
	    if(res==numseq-1) continue;
            if(val<0) val+=360.0;
            decstr[res+1].theta=val;
        }
    }
    fclose(file);

    if ( ( file=fopen(chi_file,"rt") )!=NULL )
    {
        char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
        while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,file))
        {
            int res;
            double val=0.0;
            sscanf(line,"%d %lf",&res,&val);
            if(val<0) val+=360.0;
            decstr[res].chi=val;
        }
    }
    fclose(file);

    return tor2coord( decstr, numseq );
}

bool Geometry::apply_torsions( point3f *decstr, double *vars, int numseq )
{
    double len_p_c=3.8761;
    double len_c_p=3.8705;
    double len_c_n=3.3912;
    double ang_p_c_p=104.8433;
    double ang_c_p_c=105.6791;

    for ( int i=0; i<numseq; i++ )
    {
        decstr[i].len_p_c=len_p_c;
        decstr[i].len_c_p=len_c_p;
        decstr[i].len_c_n=len_c_n;
        decstr[i].ang_p_c_p=ang_p_c_p;
        decstr[i].ang_c_p_c=ang_c_p_c;
    }

    //dummy torsions
    decstr[0].theta = 0.0;
    decstr[0].eta = 0.0;
    for ( int i=1; i<numseq; i++ )
    {
        if ( vars[(i-1)*2]>=360.0 )
        {
            while ( vars[(i-1)*2]>=360.0 )
            {
                vars[(i-1)*2]-=360.0;
            }
        }
        else if ( vars[(i-1)*2]<0.0 )
        {
            while ( vars[(i-1)*2]<0.0 )
            {
                vars[(i-1)*2]+=360.0;
            }
        }

        if ( vars[(i-1)*2+1]>=360.0 )
        {
            while ( vars[(i-1)*2+1]>=360.0 )
            {
                vars[(i-1)*2+1]-=360.0;
            }
        }
        else if ( vars[(i-1)*2+1]<0.0 )
        {
            while ( vars[(i-1)*2+1]<0.0 )
            {
                vars[(i-1)*2+1]+=360.0;
            }
        }

        decstr[i].theta=vars[(i-1)*2];
        decstr[i].eta=vars[(i-1)*2+1];
    }
    return tor2coord( decstr, numseq );
}

#endif
