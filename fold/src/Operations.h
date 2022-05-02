/*******************************************************************************************
**  
**  Functions for performing various operations used by different classes
**
**  Please report bugs and questions to robpearc@umich.edu
**
*******************************************************************************************/

#ifndef OPERATIONS_H
#define OPERATIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <malloc.h>
#include <limits.h>
#include <unistd.h> //for chdir and getexepath in comet

#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <complex>
#include "CommonParameters.h"

#define eps_margin  1e-7
#define MODULUS    2147483647 /* DON'T CHANGE THIS VALUE                  */
#define MULTIPLIER 48271      /* DON'T CHANGE THIS VALUE                  */
#define STREAMS    256        /* # of streams, DON'T CHANGE THIS VALUE    */
#define DEFAULT_STATE    123456789  /* initial seed, use 0 < DEFAULT < MODULUS  */

static long seed[STREAMS] = {DEFAULT_STATE};  /* current state of each stream   */
static int  stream        = 0;          /* stream index, 0 is the default */
static int  initialized   = 0;          /* test for stream initialization */

/* returns a pseudo-random real number uniformly distributed 
 * between 0.0 and 1.0. */
double Random(void)
{
    const long Q = MODULUS / MULTIPLIER;
    const long R = MODULUS % MULTIPLIER;
    long t;

    t = MULTIPLIER * (seed[stream] % Q) - R * (seed[stream] / Q);
    if (t > 0) seed[stream] = t;
    else seed[stream] = t + MODULUS;
    return ((double) seed[stream] / MODULUS);
}

double* new1DArr( int row )
{
    double *ans = new double[row];
    for(int i=0; i<row; i++)
    {
        ans[i] = 0.0;
    }
    return ans;
}

template<class T>
T* new1DArrT( int row )
{
    T *ans=new T[row];
    return ans;
}


double** new2DArr( int row, int col )
{
    double **ans=new double*[row];
    for( int i=0; i<row; i++ )
    {
        ans[i]=new double[col];
    }
    for( int i=0; i<row; i++ )
    {
        for( int j=0; j<col; j++ )
            ans[i][j] = 0.0;
    }

    return ans;
}

template<class T>
T** new2DArrT( int row, int col )
{
    T **ans=new T*[row];
    for( int i=0; i<row; i++ )
    {
        ans[i]=new T[col];
    }
    return ans;
}

void release2DArr( int n, double ** Arr )
{
    if( Arr != NULL )
    {
        for(int i = 0; i < n; i++){
            delete [] Arr[i];
        }
        delete [] Arr;
        Arr = NULL;
    }
}

template<class T>
void release2DArrT( int n, T** Arr )
{
    if( Arr != NULL ){
        for( int i = 0; i < n; i++ )
        {
            delete [] Arr[i];
        }
        delete [] Arr;
        Arr = NULL;
    }
}

double*** new3DArr(int row, int col, int dim){
    double ***ans;
	ans = new double **[row];
	for(int i=0; i<row; i++){
		ans[i] = new double *[col];
	}
	for(int i=0; i<row; i++){
		for(int j=0; j< col; j++)
			ans[i][j] = new double [dim];
	}
	for(int t=0; t<row; t++){
		for(int i=0; i< col; i++)
			for(int j=0; j< dim; j++)
				ans[t][i][j] = 0.0;
	}
	return ans;
}


void release3DArr(int row, int col, int *** Arr){
	if( Arr != NULL ){
		for(int t=0; t<row; t++){
			for(int i=0; i< col; i++)
				delete [] Arr[t][i];
		}
		for(int t=0; t<row; t++)
			delete [] Arr[t];
		delete [] Arr;
		Arr = NULL;
	}
}


void release3DArr(int row, int col, double *** Arr){
	if( Arr != NULL ){
		for(int t=0; t<row; t++){
			for(int i=0; i< col; i++)
				delete [] Arr[t][i];
		}
		for(int t=0; t<row; t++)
			delete [] Arr[t];
		delete [] Arr;
		Arr = NULL;
	}
}

int*** new3DIntArr(int row, int col, int dim){
	int ***ans;
	ans = new int **[row];
	for(int i=0; i<row; i++){
		ans[i] = new int *[col];
	}
	for(int i=0; i<row; i++){
		for(int j=0; j< col; j++)
			ans[i][j] = new int [dim];
	}
	for(int t=0; t<row; t++){
		for(int i=0; i< col; i++)
			for(int j=0; j< dim; j++)
				ans[t][i][j] = 0;
	}
	return ans;
}

void release3DIntArr(int row, int col, int *** Arr){
	if( Arr != NULL ){
		for(int t=0; t<row; t++){
			for(int i=0; i< col; i++)
				delete [] Arr[t][i];
		}
		for(int t=0; t<row; t++)
			delete [] Arr[t];
		delete [] Arr;
		Arr = NULL;
	}
}

point3d setv( double tx, double ty, double tz )
{
    point3d temp;
    temp.x=tx;
    temp.y=ty;
    temp.z=tz;
    return temp;
}

double norm( point3d p1 )
{
    double t = sqrt(p1.x*p1.x+p1.y*p1.y+p1.z*p1.z);
    return t;
}

point3d scal( point3d p1, double f )
{
    point3d temp;
    temp.x = f*p1.x;
    temp.y = f*p1.y;
    temp.z = f*p1.z;
    return temp;
}

point3d unit( point3d p1 )
{
    point3d temp;
    double t=sqrt(p1.x*p1.x+p1.y*p1.y+p1.z*p1.z);
    if(t<1e-40) 
    {
        temp.x=0;
        temp.y=0;
        temp.z=0;
        return temp;
    }
    temp.x=p1.x/t;
    temp.y=p1.y/t;
    temp.z=p1.z/t;
    return temp;
}

double dotv( point3d p1, point3d p2 )
{
    double temp;
    temp=p1.x*p2.x+p1.y*p2.y+p1.z*p2.z;
    return temp;
}

point3d minu( point3d p1, point3d p2 )
{
    point3d temp;
    temp.x=p1.x-p2.x;
    temp.y=p1.y-p2.y;
    temp.z=p1.z-p2.z;
    return temp;
}

point3d addv( point3d p1, point3d p2 )
{
    point3d temp;
    temp.x=p1.x+p2.x;
    temp.y=p1.y+p2.y;
    temp.z=p1.z+p2.z;
    return temp;
}

int maxnormal( point3d p1 )
{
    int i;
    double tem,temp[3];	
	temp[0]=fabs(p1.x);
	temp[1]=fabs(p1.y);
	temp[2]=fabs(p1.z);
	i=0;
	tem=temp[0];
	if(temp[1]>tem){
		i=1;
		tem=temp[1];
	}
	if(temp[2]>tem){
		i=2;
		tem=temp[2];
	}
	return i;
}

point3d prod( point3d p1, point3d p2 )
{
    point3d temp;
    temp.x = p1.y*p2.z-p1.z*p2.y;
    temp.y = p1.z*p2.x-p1.x*p2.z;
    temp.z = p1.x*p2.y-p1.y*p2.x;
    return temp;
}

double dsqugaussian( double x, double sigma, double miu )
{
    double tup=-2*(x-miu)/(2.0*sigma*sigma);
    return tup;
}

double squgaussian( double x, double sigma, double miu )
{
    double tup=-(x-miu)*(x-miu)/(2.0*sigma*sigma);
    return tup;
}

void q2rot( double *q, double rmax[] )
{
    int i,j;
    for(i=0;i<3;i++)
    {
        for(j=0;j<3;j++)
        {
            rmax[i*3+j]=0;
        }
    }
    rmax[0]=q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3];
    rmax[1]=2*(q[1]*q[2]-q[0]*q[3]);
    rmax[2]=2*(q[1]*q[3]+q[0]*q[2]);
    rmax[3]=2*(q[1]*q[2]+q[0]*q[3]);
    rmax[4]=q[0]*q[0]-q[1]*q[1]+q[2]*q[2]-q[3]*q[3];
    rmax[5]=2*(q[2]*q[3]-q[0]*q[1]);
    rmax[6]=2*(q[1]*q[3]-q[0]*q[2]);
    rmax[7]=2*(q[2]*q[3]+q[0]*q[1]);
    rmax[8]=q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3];
}

point3d mmat( double mat[9], point3d p1 )
{
    point3d temp;
    temp.x=p1.x*mat[0]+p1.y*mat[1]+p1.z*mat[2];
    temp.y=p1.x*mat[3]+p1.y*mat[4]+p1.z*mat[5];
    temp.z=p1.x*mat[6]+p1.y*mat[7]+p1.z*mat[8];
    return temp;
}

double angv( point3d p1, point3d p2 )
{
    double t1=norm(p1);
    double t2=norm(p2);
    if ( t1<eps_margin || t2<eps_margin )
    {
        return 0;
    }
    double t3=dotv(p1,p2)/t1/t2;
    if ( t3<-1.0 )
    {
        t3=-1.0;
    }
    else if ( t3>1.0 )
    {
		t3=1.0;
    }
    return(acos(t3));
}

// get the residue ID
int getntid( char nucleicname )
{
    int i;
    //empty
    if ( nucleicname==' ' )
    {
        return 4;
    }
    // one in 26
    for ( i=0; i<5; i++ )
    {
        if(nucleicname==ntid[i])
        {
            return i;
        }
    }
    // unknown
    //printf("unknown amino name %c%c%c\n",aminoname[0],aminoname[1],aminoname[2]);
    return 4;
}

// torsion angle to coordinates
bool tor2pos22(float xi,float yi,float zi,float xj,float yj,float zj,float xk, float yk,
               float zk,float tang,float tleng,float tinner, float *xl, float *yl, float *zl){
    point3d p12, p23, e1, e2, e3;
    point3d p1, p2, p3;
    double tpangle,tcos,tsin,rmat[9];
    double q[4];
    p12.x=xi-xj; p12.y=yi-yj; p12.z=zi-zj;
    p23.x=xk-xj; p23.y=yk-yj; p23.z=zk-zj;

    if(norm(p12)<eps_margin*0.00001 || norm(p23)<eps_margin*0.00001 || angv(p12, p23)<eps_margin*0.001 || (PI-angv(p12,p23))<eps_margin*0.001){
        int imax=maxnormal(p23);
        if(imax==0){
            yj+=0.01f;
        }
        else if(imax==1){
            zj+=0.01f;
        }
        else{
            xj+=0.01f;
        }
        p12.x=xi-xj; p12.y=yi-yj; p12.z=zi-zj;
        p23.x=xk-xj; p23.y=yk-yj; p23.z=zk-zj;
        cout<<"make adjustment tor2pos22"<<endl;
    }

    e2=unit(p23);
    e3=unit(p12);
    e1=prod(e2,e3);
    e1=unit(e1);
    if(norm(e1)<eps_margin || norm(e2)<eps_margin){
        printf("wrong in tor2pos22 [%f %f %f] [%f %f %f] [%f %f %f]\n",xi,yi,zi,xj,yj,zj,xk,yk,zk);
        *xl = xk; *yl = yk; *zl = zk;
        return false;
    }

    p1 = scal(e2,tleng);
    tpangle = (PI-tinner)/2.0;
    tcos = cos(tpangle);
    tsin = sin(tpangle);
    q[0]=tcos; q[1]=tsin*e1.x; q[2]=tsin*e1.y; q[3]=tsin*e1.z;
    q2rot(q, rmat);
    p2 = mmat(rmat, p1);

    tpangle=tang/2.0;
    tcos=cos(tpangle);
    tsin=sin(tpangle);
    q[0]=tcos; q[1]=tsin*e2.x; q[2]=tsin*e2.y; q[3]=tsin*e2.z;
    q2rot(q,rmat);
    p3=mmat(rmat,p2);

    *xl=p3.x+xk; *yl=p3.y+yk; *zl=p3.z+zk;
    return true;
}

// get the torsion angle from the coordinates
double phi( double xi, double yi, double zi,
            double xj, double yj, double zj,
            double xk, double yk, double zk,
            double xl, double yl, double zl )
{
    double xij,yij,zij,
    xkj,ykj,zkj,  
    xkl,ykl,zkl,
    dxi,dyi,dzi,
    gxi,gyi,gzi,
    bi,bk,ct,
    boi2,boj2,
    z1,z2,ap,s,
    bioj,bjoi;
		
    /* Calculate the vectors C,B,C*/
    xij = xi - xj;
    yij = yi - yj;
    zij = zi - zj;
    xkj = xk - xj;
    ykj = yk - yj;
    zkj = zk - zj;
    xkl = xk - xl;
    ykl = yk - yl;
    zkl = zk - zl;
	
    /* Calculate the normals to the two planes n1 and n2
       this is given as the cross products:
         AB x BC
        --------- = n1
        |AB x BC|
	
         BC x CD
        --------- = n2
        |BC x CD|
    */
    dxi = yij * zkj - zij * ykj;     /* Normal to plane 1                */
    dyi = zij * xkj - xij * zkj;
    dzi = xij * ykj - yij * xkj;
    gxi = zkj * ykl - ykj * zkl;     /* Mormal to plane 2                */
    gyi = xkj * zkl - zkj * xkl;
    gzi = ykj * xkl - xkj * ykl;
	
    /* Calculate the length of the two normals                           */
    bi = dxi * dxi + dyi * dyi + dzi * dzi;
    bk = gxi * gxi + gyi * gyi + gzi * gzi;
    ct = dxi * gxi + dyi * gyi + dzi * gzi;
	
    boi2 = 1./bi;
    boj2 = 1./bk;
    bi   = (double)sqrt((double)bi);
    bk   = (double)sqrt((double)bk);
    if(bi<eps_margin*0.01 || bk<eps_margin*0.01) return 180;
    z1   = 1./bi;
    z2   = 1./bk;
    bioj = bi * z2;
    bjoi = bk * z1;
    ct   = ct * z1 * z2;
    if (ct >  1.0)   ct = 1.0;
    if (ct < (-1.0)) ct = -1.0;
    ap   = acos(ct);
	
    s = xkj * (dzi * gyi - dyi * gzi)
      + ykj * (dxi * gzi - dzi * gxi)
      + zkj * (dyi * gxi - dxi * gyi);
	
    if (s < 0.0) ap = -ap;
	
    ap = (ap > 0.0) ? PI-ap : -(PI+ap);

    // angle
    ap *= (double)180.0/PI;
    if(ap<0) ap+=360;
    return(ap);
}

double maxinthree( float fx, float fy, float fz )
{
    double tmax=fx;
    if(tmax<fy) tmax=fy;
    if(tmax<fz) tmax=fz;
    return tmax;
}

// determine the protein trype according to secondary sctructure
int getprotype( vector<string> SStype ){
    int protype;
    int numstrand = 0; 
    for(int i=0; i<SStype.size(); i++){
        if(SStype[i] == "E"){
            numstrand ++;
            while(SStype[i] == "E" && i<SStype.size()) i ++;
        }
    } 
    if (numstrand <2) protype = 1; //alpha
    else protype = 2;
    return protype;
}

bool checkfile( string path )
{
    bool flag = false;
    std::fstream _file;
    _file.open(path.c_str(), ios::in);
    if(_file)
    {
        flag = true;
        _file.close();
    }
    return flag;
}


#endif
