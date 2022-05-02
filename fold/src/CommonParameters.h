/*******************************************************************************************
**  
**  Common parameters and data tables
**
**  Please report bugs and questions to robpearc@umich.edu
**
*******************************************************************************************/

#ifndef COMMONPARAMETERS_H
#define COMMONPARAMETERS_H

#include <string>// TO use string variables
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <sys/stat.h>

using namespace std;

#define PI 3.14159265358979323846
#define raddeg PI/180.0
#define degrad 180.0/PI
#define lennca 1.460f
#define delnca 0.004f*20
#define lencac 1.525f
#define delcac 0.004f*20
#define lencn 1.338f
#define delcn 0.005f*15
#define lennn 3.131f
#define lencc 3.231f
#define lencaca 3.813f
#define delcaca 0.019f*10
#define lencan1 2.441f
#define delcan1 0.036f*4
#define lencca1 2.446f
#define delcca1 0.036f*4
#define lennc 2.460f
#define delnc 0.012f*10
#define angncac 111.008f
#define angcacn 116.617f
#define angcnca 121.614f
#define MAX_LENGTH_ONE_LINE_IN_FILE  1024
#define MAX_ENERGY_TERM_NUM           100

//-------------------- Common Data Structures ----------------------->
typedef struct ssef
{
    char res;
    char ss;
    float a[3];//coil helix sheet
}ssef;//no change

typedef struct point3s
{
    float x,y,z;
}point3s;

typedef struct point3d{
        double x,y,z;
}point3d;

typedef struct point3f
{
    float x,y,z;
    float eta, theta, chi;
    float len_c_p,len_p_c,len_c_n;
    float ang_p_c_p,ang_c_p_c;
    float angl,leng;
    point3s ptp,ptn,ptc2,ptc4,ptc6;
    point3d f1_eta,f2_eta,f1_theta,f2_theta;

    //char ss2;//psipred in the fragment
    //char ssm;//mycalculation
    char nt; // one letter amino acid code (realseq ACD)
    unsigned char ntind;//realACD index
}point3f;

typedef struct sssegment
{
    int init;
    int term;
    char ss;
}sssegment;

//---------------------- Data Tables -------------------->
static double torN[]={
105.814,102.751,104.225,105.158,
};

static double angN[]={
101.164,99.153,99.483,100.687,
};

static double distN[]={
3.392,3.398,3.383,3.413,
};

static double torO4[]={
114.0196,108.4575,112.4488,112.2701,
};

static double angO4[]={
131.5051,131.0396,130.9602,130.8516,
};

static double distO4[]={
1.4505,1.4512,1.4512,1.4509,
};

static double torC1[]={
128.6302,129.5699,131.9457,127.7324,
};

static double angC1[]={
109.8054,109.7353,109.7436,109.7484,
};

static double distC1[]={
1.4127,1.4141,1.4132,1.4134,
};

static double angChi[]={
126.5274,118.8795,126.4582,117.6994,
};

static double distChi[]={
1.3724,1.3982,1.3735,1.3814,
};

static double torAfterChi[]={
359.3861,179.8174,359.8648,179.7748,
};

static double angAfterChi[]={
159.6609,89.6304,161.8754,88.5634,
};

static double distAfterChi[]={
2.1999,2.3242,2.2144,2.4563,
};

static double torC6[]={
180.2857,0.0,180.2485,0.0,
};

static double angC6[]={
63.4597,0.0,61.8550,0.0,
};

static double distC6[]={
2.3104,0.0,2.4477,0.0,
};

static char ntid[]= {
'A','C','G','U', 'X',
};

static const char aad3[][4] = {
"  A","  C","  G","  U","  X"
};

#endif
