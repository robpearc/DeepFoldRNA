/*******************************************************************************************
**  
**  Functions for calculating the derivatives/gradients
**
**  Please report bugs and questions to robpearc@umich.edu
**
*******************************************************************************************/

#ifndef DERIVATIVES
#define DERIVATIVES

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "CommonParameters.h"
#include "Operations.h"

class Derivatives
{
public:
    Derivatives();
    
    void dih_func( point3d M, point3d v, point3d w, point3d & F1, point3d & F2 );

    void ang_func( point3d M, point3d w, point3d & F1, point3d & F2 );

    void dist_deriv( point3d atom1, point3d atom2, double distance, point3d & f1, point3d & f2 );

    bool atom1_deriv_dih_first( point3d p1, point3d p2, point3d p3, point3d p4, double & x,
                                point3d & F1, point3d & F2);


    bool atom2_deriv_dih_first( point3d p1, point3d p2, point3d p3, point3d p4, double & x,
                                point3d & F1, point3d & F2 );

    void atom1_dih_deriv( point3d p1, point3d p2, point3d p3, point3d p4, double & theta,
                          point3d & f1, point3d & f2 );

    void atom2_dih_deriv( point3d p1, point3d p2, point3d p3, point3d p4, double & theta,
                          point3d & f1, point3d & f2 );

    void sec_deriv_dih( point3d p1, point3d p2, point3d p3, point3d p4, double x,
                        double & theta, point3d & F1, point3d & F2 );


    void atom1_angle_deriv( point3d atom1, point3d vec1, point3d vec2, double len1, double len2, 
                            double angle, double dtheta_dangle, point3d & F1,point3d & F2 );

    void angular_deriv( point3d p1, point3d p2, point3d p3,
                        point3d & f1_p1, point3d & f2_p1, point3d & f1_p2,
                        point3d & f2_p2, point3d & f1_p3, point3d & f2_p3 );

    virtual ~Derivatives();
};

Derivatives::Derivatives()
{
}

void Derivatives::dih_func( point3d M, point3d v, point3d w, point3d & F1, point3d & F2 )
{
    point3d f2 = prod(v,w);
    F2 = addv(F2,f2);
    F1 = addv(F1,prod(f2,M));
}

void Derivatives::ang_func( point3d M, point3d w, point3d & F1, point3d & F2 )
{
    F2 = addv(F2,w);
    F1 = addv(F1,prod(w,M));
}

void Derivatives::dist_deriv( point3d atom1, point3d atom2, double distance, 
                              point3d &f1, point3d &f2 )
{
    f2 = minu( atom1, atom2 );
    if ( distance != 0.0 ) 
    {
        f1 = prod( atom1, atom2 );
        f1 = scal( f1, 1.0 / distance );
        f2 = scal( f2, 1.0 / distance );
    } 
    else
    {
        f1.x=0.0,f1.y=0.0,f1.z=0.0;
    }
}

bool Derivatives::atom1_deriv_dih_first( point3d atom1, point3d atom2, point3d atom3, point3d atom4, 
                                         double &angle, point3d &F1, point3d &F2)
{
    point3d scalv1, scalv2;
    F1.x = 0.0, F1.y=0.0, F1.z=0.0;
    F2.x = 0.0, F2.y=0.0, F2.z=0.0;

    point3d vec1 = minu(atom1,atom2);
    point3d vec2 = minu(atom2,atom3);
    point3d vec3 = minu(atom3,atom4);
    point3d cross12=prod(vec1,vec2);
    point3d cross23=prod(vec2,vec3);

    double len12=norm(cross12);
    double len23=norm(cross23);
    if ( len12 < 1e-9 || len23 < 1e-9 ) return true;

    angle = dotv( cross12, cross23 ) / ( len12 * len23 );
    
    scalv2 = scal( vec2, 1.0/(len12*len23) );
    dih_func( atom1, scalv2, cross23 , F1, F2);

    scalv2 = scal( vec2, -1.0*angle/(len12*len12) );
    dih_func( atom1, scalv2, cross12, F1, F2 );
    
    return false;
}

bool Derivatives::atom2_deriv_dih_first( point3d atom1, point3d atom2, point3d atom3, point3d atom4, 
                                        double &angle, point3d & F1, point3d & F2 )
{
    point3d scalv1, scalv2, scalv3;
    F1.x=0.0,F1.y=0.0,F1.z=0.0;
    F2.x=0.0,F2.y=0.0,F2.z=0.0;

    point3d vec1=minu(atom1,atom2);
    point3d vec2=minu(atom2,atom3);
    point3d vec3=minu(atom3,atom4);
    point3d cross12 = prod(vec1,vec2);
    point3d cross23 = prod(vec2,vec3);

    double len12 = norm(cross12);
    double len23 = norm(cross23);

    if ( len12 < 1e-9 || len23 < 1e-9 ) return true;

    angle = dotv( cross12, cross23 ) / ( len12 * len23 );

    scalv2 = scal( vec2, -1.0/(len12*len23) );
    dih_func( atom2, scalv2, cross23 , F1, F2);

    scalv1 = scal( vec1, -1.0/(len12*len23) );
    dih_func( atom2, scalv1, cross23 , F1, F2);

    scalv3 = scal( vec3, 1.0/(len12*len23) );
    dih_func( atom2, scalv3, cross12 , F1, F2);

    scalv2 = scal( vec2, angle/(len12*len12) );
    dih_func( atom2, scalv2, cross12, F1, F2 );

    scalv1 = scal( vec1, angle/(len12*len12) );
    dih_func( atom2, scalv1, cross12, F1, F2 );

    scalv3 = scal( vec3, -1.0*angle/(len23*len23) );
    dih_func( atom2, scalv3, cross23, F1, F2 );

    return false;
}

void Derivatives::sec_deriv_dih( point3d atom1, point3d atom2, point3d atom3, point3d atom4, 
                                 double angle, double &theta, point3d & F1, point3d & F2 )
{
    double thetaU = acos( angle );
    double tors=phi(atom1.x,atom1.y,atom1.z,atom2.x,atom2.y,atom2.z,
                    atom3.x,atom3.y,atom3.z,atom4.x,atom4.y,atom4.z);
    if(tors>180.0) tors-=360.0;
    theta = (tors)*raddeg;
    
    if( abs( abs( theta ) - thetaU ) >= 1e-2 ) return;
    angle = min( max( cos(179.9*raddeg), angle ), cos(0.1*raddeg) );
    double dthetaU_dx = -1 / sqrt( 1- angle*angle );
    double dtheta_dthetaU( theta < 0.0 ? -1 : 1 );
    F1 = scal(F1,(dtheta_dthetaU * dthetaU_dx));
    F2 = scal(F2,(dtheta_dthetaU * dthetaU_dx));
}

void Derivatives::atom1_dih_deriv( point3d atom1, point3d atom2, point3d atom3, point3d atom4, 
                                   double &theta, point3d &f1, point3d &f2 )
{
    double angle;
    bool is_colinear = atom1_deriv_dih_first( atom1, atom2, atom3, atom4, angle, f1, f2 );
    if ( is_colinear )
    {
        return;
    }
    sec_deriv_dih( atom1, atom2, atom3, atom4, angle, theta, f1, f2 );
}

void Derivatives::atom2_dih_deriv( point3d atom1, point3d atom2, point3d atom3, point3d atom4, double & theta, 
                                   point3d & f1, point3d & f2)
{
    double angle = 0.0;
    bool is_colinear = atom2_deriv_dih_first( atom1, atom2, atom3, atom4, angle, f1, f2 );
    if ( is_colinear )
    {
        return;
    }
    sec_deriv_dih( atom1, atom2, atom3, atom4, angle, theta, f1, f2 );
}

void Derivatives::atom1_angle_deriv( point3d atom1, point3d vec1, point3d vec2, double len1, double len2,
                                     double angle, double dtheta_dangle, point3d & F1, point3d & F2 )
{
    if ( len1 < 1e-9 || len2 < 1e-9 ) return;

    point3d f1=setv(0.0,0.0,0.0),f2=setv(0.0,0.0,0.0);
    point3d scalv1=setv(0.0,0.0,0.0),scalv2=setv(0.0,0.0,0.0);

    scalv2 = scal( vec2, 1.0 /(len1*len2) );
    ang_func( atom1, scalv2, f1, f2 );

    scalv1 = scal(vec1, -1.0*angle/(len1*len1) );
    ang_func( atom1, scalv1, f1, f2 );

    f1 = scal( f1, dtheta_dangle );
    f2 = scal( f2, dtheta_dangle );
    F1 = addv( F1, f1 );
    F2 = addv( F2, f2 );
}

void Derivatives::angular_deriv( point3d atom1, point3d atom2, point3d atom3,
                                 point3d & f1_atom1, point3d & f2_atom1, 
                                 point3d & f1_atom2, point3d & f2_atom2, 
                                 point3d & f1_atom3, point3d & f2_atom3 )
{
    f1_atom1=setv(0.0,0.0,0.0);
    f2_atom1=setv(0.0,0.0,0.0);
    f1_atom2=setv(0.0,0.0,0.0);
    f2_atom2=setv(0.0,0.0,0.0);
    f1_atom3=setv(0.0,0.0,0.0);
    f2_atom3=setv(0.0,0.0,0.0);

    point3d atom1_atom2 = minu( atom1, atom2 );
    point3d atom3_atom2 = minu( atom3, atom2 );
    double len12 = norm( atom1_atom2 );
    double len23 = norm( atom3_atom2 );
    double length = ( len12 * len23 );
    if ( length < 1e-12 ) return;
    
    point3d atom2_atom1 = minu( atom2, atom1 );
    point3d atom3_atom1 = minu( atom3, atom1 );
    point3d atom2_atom3 = minu( atom2, atom3 );
    point3d atom1_atom3 = minu( atom1, atom3 );

    double angle = dotv( atom1_atom2, atom3_atom2 ) / ( len12 * len23 );
    angle = min( max( cos( 179.9*raddeg ), angle ), cos( 0.1*raddeg ) );
    double dtheta_dangle = -1.0 / sqrt( 1.0 - angle*angle );
    atom1_angle_deriv( atom1, atom1_atom2, atom3_atom2, len12, len23, angle, dtheta_dangle, f1_atom1, f2_atom1 );
    atom1_angle_deriv( atom3, atom3_atom2, atom1_atom2, len23, len12, angle, dtheta_dangle, f1_atom3, f2_atom3 );

    double len13 = norm( atom3_atom1 );
    angle = dotv( atom2_atom1, atom3_atom1 ) / ( len12 * len13 );
    angle = min( max( cos( 179.9*raddeg ), angle ), cos( 0.1*raddeg ) );
    dtheta_dangle = -1.0 / sqrt( 1.0 - angle*angle );
    atom1_angle_deriv( atom2, atom2_atom1, atom3_atom1, len12, len13, angle, dtheta_dangle, f1_atom2, f2_atom2 );

    angle = dotv( atom2_atom3, atom1_atom3 ) / ( len23 * len13 );
    angle = min( max( cos( 179.9*raddeg ), angle ), cos( 0.1*raddeg ) );
    dtheta_dangle = -1.0 / sqrt( 1.0 - angle*angle );
    atom1_angle_deriv( atom2, atom2_atom3, atom1_atom3, len23, len13, angle, dtheta_dangle, f1_atom2, f2_atom2 );

    f1_atom2 = scal( f1_atom2, -1.0 );
    f2_atom2 = scal( f2_atom2, -1.0 );
}

Derivatives::~Derivatives() = default;

#endif 

