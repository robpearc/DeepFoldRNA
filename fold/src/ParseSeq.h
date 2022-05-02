/*******************************************************************************************
**  
**  Functions for loading the RNA sequence
**
**  Please report bugs and questions to robpearc@umich.edu
**
*******************************************************************************************/

#ifndef PARSESEQ_H
#define PARSESEQ_H

#include "CommonParameters.h"
#include "Operations.h"

class ParseSeq  
{
public:
    ParseSeq();
    int seqnum;
    char *seqdata;
    char *seqheader;
    bool loadseq(const char *seqfile);
    bool loadseq2(const char *seqfile);
    virtual ~ParseSeq();
};

ParseSeq::ParseSeq()
{
    seqnum=0;
    seqheader=NULL;
    seqdata=NULL;
}

ParseSeq::~ParseSeq()
{
    if(seqdata)
    {
        delete[]seqdata;
        seqdata=NULL;
    }
    if(seqheader)
    {
        delete[]seqheader;
        seqheader=NULL;
    }
}

bool ParseSeq::loadseq2( const char *seqfile )
{
    FILE *file;
    char tmpc;
    file= fopen(seqfile, "rt");
    if(file==NULL)
    {
        fprintf(stderr,"Error when loading sequence file %s\n",seqfile);
        return false;
    }
    seqnum=0;
    if(!seqdata) seqdata=new char[65525];
    while(!feof(file))
    {
        tmpc=fgetc(file);
        if(feof(file))
        {
            fclose(file);
            seqdata[seqnum]='\0';
            return true;
        }
        if (!((tmpc>='A' && tmpc<='Z') || tmpc=='\r' || tmpc=='\n'))
            seqdata[seqnum++]='A';
        else if(tmpc>='A' && tmpc<='Z') seqdata[seqnum++]=tmpc;
    }
    fclose(file);
    seqdata[seqnum]='\0';
    return true;
}

/* read fasta or plain text protein sequence */
bool ParseSeq::loadseq(const char *seqfile)
{
    FILE *file;
    char tmpc;
    file= fopen(seqfile, "rt");
    if(file==NULL)
    {
        fprintf(stderr,"Error when loading sequence file %s\n",seqfile);
        return false;
    }
    seqnum=0;
    if(!seqheader) seqheader=new char[65525];
    if(!seqdata) seqdata=new char[65525];
    fgets(seqheader, 65525, file);
    if(seqheader[0]!='>')
    {
        delete[]seqdata;
        seqdata=NULL;
        delete[]seqheader;
        seqheader=NULL;
        fclose(file);
        return loadseq2(seqfile);
    }
    while(!feof(file))
    {
        tmpc=fgetc(file);
        if(feof(file))
        {
            fclose(file);
            seqdata[seqnum]='\0';
            return true;
        }
        if (!((tmpc>='A' && tmpc<='Z') || tmpc=='\r' || tmpc=='\n'))
            seqdata[seqnum++]='A';
        else if(tmpc>='A' && tmpc<='Z')
            seqdata[seqnum++]=tmpc;
    }
    fclose(file);
    seqdata[seqnum]='\0';
    return true;
}

#endif 
