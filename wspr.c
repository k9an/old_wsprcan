/*
 This file is part of k9an-wsprd.
 
 File name: wspr.c
 Description: k9an-wsprd is a detector/demodulator/decoder for K1JT's
 Weak Signal Propagation Reporter (WSPR) mode.
 
 Copyright 2014, Steven Franke, K9AN
 License: GNU GPL v3
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <sndfile.h>
#include <math.h>
#include <strings.h>
#include <fftw3.h>
#include "fano.h"

unsigned char pr3[162]=
{1,1,0,0,0,0,0,0,1,0,0,0,1,1,1,0,0,0,1,0,
    0,1,0,1,1,1,1,0,0,0,0,0,0,0,1,0,0,1,0,1,
    0,0,0,0,0,0,1,0,1,1,0,0,1,1,0,1,0,0,0,1,
    1,0,1,0,0,0,0,1,1,0,1,0,1,0,1,0,1,0,0,1,
    0,0,1,0,1,1,0,0,0,1,1,0,1,0,1,0,0,0,1,0,
    0,0,0,0,1,0,0,1,0,0,1,1,1,0,1,1,0,0,1,1,
    0,1,0,0,0,1,1,1,0,0,0,0,0,1,0,1,0,0,1,1,
    0,0,0,0,0,0,0,1,1,0,1,0,1,1,0,0,0,1,1,0,
    0,0};

unsigned long readc2file(char *ptr_to_infile, double *idat, double *qdat, float *freq)
{
    float buffer[2*65536];
    double dfreq;
    int i,ntrmin;
    char *c2file[15];
    FILE* c2fh;
    
    c2fh=fopen(ptr_to_infile,"r");
    unsigned long nread=fread(c2file,sizeof(char),14,c2fh);
    nread=fread(&ntrmin,sizeof(int),1,c2fh);
    nread=fread(&dfreq,sizeof(double),1,c2fh);
    *freq=(float)dfreq;
    nread=fread(buffer,sizeof(float),2*45000,c2fh);
    
    for(i=0; i<45000; i++) {
        idat[i]=buffer[2*i];
        qdat[i]=-buffer[2*i+1];
    }
    
    if( nread == 2*45000 ) {
        return nread/2;
    } else {
        return 1;
    }
}

unsigned long readwavfile(char *ptr_to_infile, double *idat, double *qdat )
{
    // The part that reads from a wav file is based on work by:
    // Andrew Greensted - Feb 2010
    // http://www.labbookpages.co.uk
    // Version 1

    float *buffer;
    int i, j;
    int nfft1=2*1024*1024;
    int nfft2=nfft1/32; //65536, only the first 45000 points are nonzero
    int nh2=nfft2/2;
    double df=12000.0/nfft1;
    int i0=1500.0/df+0.5;
    fftw_complex *fftin, *fftout;
    fftw_plan MYPLAN;
    
// Open sound file
    SF_INFO sndInfo;
    SNDFILE *sndFile = sf_open(ptr_to_infile, SFM_READ, &sndInfo);
    if (sndFile == NULL) {
        fprintf(stderr, "Error reading source file '%s': %s\n", ptr_to_infile, sf_strerror(sndFile));
    return 1;
    }

// Check format - 16bit PCM
    if (sndInfo.format != (SF_FORMAT_WAV | SF_FORMAT_PCM_16)) {
        fprintf(stderr, "Input should be 16bit Wav\n");
        sf_close(sndFile);
        return 1;
    }

// Check channels - mono
    if (sndInfo.channels != 1) {
        fprintf(stderr, "Wrong number of channels\n");
        sf_close(sndFile);
        return 1;
    }

// Allocate memory
    buffer = malloc(sndInfo.frames * 2 * sizeof(float));
    if (buffer == NULL) {
        fprintf(stderr, "Could not allocate memory for file\n");
        sf_close(sndFile);
        return 1;
    }

// Load data
    unsigned long npoints = sf_readf_float(sndFile, buffer, sndInfo.frames);

// Check correct number of samples loaded
    if (npoints != sndInfo.frames) {
        fprintf(stderr, "Did not read enough frames for source\n");
        sf_close(sndFile);
        free(buffer);
        return 1;
    }

// Output Info
/*    printf("Read %ld frames from %s, Sample rate: %d, Length: %fs\n",
           numFrames, argv[1], sndInfo.samplerate, (float)numFrames/sndInfo.samplerate);
 */

//    printf("Read %ld frames from %s \n",npoints, argv[3]);

    sf_close(sndFile);

    if ( npoints < 12000*110 ) {
        printf("file length is %lu seconds\n",npoints/12000);
        return 1;
    }

//    printf("%d %d %d %d\n",i0, nfft1, nfft2, npoints);
    fftin=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nfft1);
    fftout=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nfft1);
    MYPLAN = fftw_plan_dft_1d(nfft1, fftin, fftout, FFTW_FORWARD, FFTW_ESTIMATE);
    
    for (i=0; i<npoints; i++) {
        fftin[i][0]=buffer[i];
        fftin[i][1]=0.0;
    }
    for (i=npoints; i<nfft1; i++) {
        fftin[i][0]=0.0;
        fftin[i][1]=0.0;
    }
    
    fftw_execute(MYPLAN);
    
    fftw_free(fftin);
    fftw_destroy_plan(MYPLAN);
    
    fftin=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nfft2);
    for (i=0; i<nfft2; i++){
        j=i0+i;
        if( i>nh2 )
        j=j-nfft2;
        fftin[i][0]=fftout[j][0];
        fftin[i][1]=fftout[j][1];
    }

    fftw_free(fftout);
    fftout=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nfft2);
    MYPLAN = fftw_plan_dft_1d(nfft2, fftin, fftout, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(MYPLAN);
    
    for (i=0; i<nfft2; i++) {
        idat[i]=fftout[i][0]/1000.0;
        qdat[i]=fftout[i][1]/1000.0;
    }

    fftw_free(fftin);
    fftw_free(fftout);
    fftw_destroy_plan(MYPLAN);

    return nfft2;
}

void getStats(double *id, double *qd, long np, double *mi, double *mq, double *mi2, double *mq2, double *miq)
{
    double sumi=0.0;
    double sumq=0.0;
    double sumi2=0.0;
    double sumq2=0.0;
    double sumiq=0.0;
    float imax=-1e30, imin=1e30, qmax=-1e30, qmin=1e30;
    
    int i;
    
    for (i=0; i<np; i++) {
        sumi=sumi+id[i];
        sumi2=sumi2+id[i]*id[i];
        sumq=sumq+qd[i];
        sumq2=sumq2+qd[i]*qd[i];
        sumiq=sumiq+id[i]*qd[i];
        if( id[i]>imax ) imax=id[i];
        if( id[i]<imin ) imin=id[i];
        if( qd[i]>qmax ) qmax=qd[i];
        if( qd[i]<qmin ) qmin=qd[i];
    }
    *mi=sumi/np;
    *mq=sumq/np;
    *mi2=sumi2/np;
    *mq2=sumq2/np;
    *miq=sumiq/np;
    
//    printf("imax %f  imin %f    qmax %f  qmin %f\n",imax, imin, qmax, qmin);
}

void sync_and_demodulate(
double *id,
double *qd,
long np,
unsigned char *symbols,
float *f1,
float fstep,
int *shift1,
int lagmin, int lagmax, int lagstep,
float *drift1,
float *sync,
int mode)
{
    // mode is the last argument:
    // 0 no frequency or drift search. find best time lag.
    // 1 no time lag or drift search. find best frequency.
    // 2 no frequency or time lag search. calculate soft-decision symbols
    //   using passed frequency and shift.

    float dt=1.0/375.0, df=375.0/256.0,fbest;
    long int i, j, k;
    double pi=4.*atan(1.0);
    float f0=0.0,fp,ss;
    int lag;
    
    double
    i0[162],q0[162],
    i1[162],q1[162],
    i2[162],q2[162],
    i3[162],q3[162];
    
    double p0,p1,p2,p3,cmet,totp,syncmax,fac;
    double
    c0[256],s0[256],
    c1[256],s1[256],
    c2[256],s2[256],
    c3[256],s3[256];
    double
    dphi0, cdphi0, sdphi0,
    dphi1, cdphi1, sdphi1,
    dphi2, cdphi2, sdphi2,
    dphi3, cdphi3, sdphi3;
    float fsum=0.0, f2sum=0.0, fsymb[162];

    int best_shift = 0, ifreq;
    int ifmin=0, ifmax=0;
    
    syncmax=-1e30;


// mode is the last argument:
// 0 no frequency or drift search. find best time lag.
// 1 no time lag or drift search. find best frequency.
//      2 no frequency or time lag search. calculate soft-decision symbols
//        using passed frequency and shift.
    if( mode == 0 ) {
        ifmin=0;
        ifmax=0;
        fstep=0.0;
        f0=*f1;
    }
    if( mode == 1 ) {
        lagmin=*shift1;
        lagmax=*shift1;
        ifmin=-5;
        ifmax=5;
        f0=*f1;
    }
    if( mode == 2 ) {
        lagmin = *shift1;
        lagmax = *shift1;
        ifmin=0;
        ifmax=0;
        f0=*f1;
    }

    for(ifreq=ifmin; ifreq<=ifmax; ifreq++)
    {
        f0=*f1+ifreq*fstep;
        for(lag=lagmin; lag<=lagmax; lag=lag+lagstep)
        {
            ss=0.0;
            totp=0.0;
            for (i=0; i<162; i++)
            {
                fp = f0 + ((float)*drift1/2.0)*((float)i-81.0)/81.0;
                
                dphi0=2*pi*(fp-1.5*df)*dt;
                cdphi0=cos(dphi0);
                sdphi0=sin(dphi0);
                dphi1=2*pi*(fp-0.5*df)*dt;
                cdphi1=cos(dphi1);
                sdphi1=sin(dphi1);
                dphi2=2*pi*(fp+0.5*df)*dt;
                cdphi2=cos(dphi2);
                sdphi2=sin(dphi2);
                dphi3=2*pi*(fp+1.5*df)*dt;
                cdphi3=cos(dphi3);
                sdphi3=sin(dphi3);
                
                c0[0]=1;
                s0[0]=0;
                c1[0]=1;
                s1[0]=0;
                c2[0]=1;
                s2[0]=0;
                c3[0]=1;
                s3[0]=0;
                
                for (j=1; j<256; j++) {
                    c0[j]=c0[j-1]*cdphi0-s0[j-1]*sdphi0;
                    s0[j]=c0[j-1]*sdphi0+s0[j-1]*cdphi0;
                    c1[j]=c1[j-1]*cdphi1-s1[j-1]*sdphi1;
                    s1[j]=c1[j-1]*sdphi1+s1[j-1]*cdphi1;
                    c2[j]=c2[j-1]*cdphi2-s2[j-1]*sdphi2;
                    s2[j]=c2[j-1]*sdphi2+s2[j-1]*cdphi2;
                    c3[j]=c3[j-1]*cdphi3-s3[j-1]*sdphi3;
                    s3[j]=c3[j-1]*sdphi3+s3[j-1]*cdphi3;
                }
                
                i0[i]=0.0;
                q0[i]=0.0;
                i1[i]=0.0;
                q1[i]=0.0;
                i2[i]=0.0;
                q2[i]=0.0;
                i3[i]=0.0;
                q3[i]=0.0;
 
                for (j=0; j<256; j++)
                {
                    k=lag+i*256+j;
                    if( (k>0) & (k<np) ) {
                    i0[i]=i0[i]+id[k]*c0[j]+qd[k]*s0[j];
                    q0[i]=q0[i]-id[k]*s0[j]+qd[k]*c0[j];
                    i1[i]=i1[i]+id[k]*c1[j]+qd[k]*s1[j];
                    q1[i]=q1[i]-id[k]*s1[j]+qd[k]*c1[j];
                    i2[i]=i2[i]+id[k]*c2[j]+qd[k]*s2[j];
                    q2[i]=q2[i]-id[k]*s2[j]+qd[k]*c2[j];
                    i3[i]=i3[i]+id[k]*c3[j]+qd[k]*s3[j];
                    q3[i]=q3[i]-id[k]*s3[j]+qd[k]*c3[j];
                    }
                }
                p0=i0[i]*i0[i]+q0[i]*q0[i];
                p1=i1[i]*i1[i]+q1[i]*q1[i];
                p2=i2[i]*i2[i]+q2[i]*q2[i];
                p3=i3[i]*i3[i]+q3[i]*q3[i];
            
                totp=totp+p0+p1+p2+p3;
                cmet=(p1+p3)-(p0+p2);
                ss=ss+cmet*(2*pr3[i]-1);
                
                if( mode == 2)
                {
                    fsymb[i]=pr3[i]*(p3-p1)+(1-pr3[i])*(p2-p0);
                }

            }
            
            if( (mode <= 1) && (ss/totp > syncmax) ) {
                syncmax=ss/totp;
                best_shift=lag;
                fbest=f0;
            }
            

        } // lag loop
    } //freq loop

    if( mode <=1 ) {
        *sync=syncmax;
        *shift1=best_shift;
        *f1=fbest;
        return;
    }
    
    if( mode == 2 ) {
        for (i=0; i<162; i++) {
            fsum=fsum+fsymb[i]/162.0;
            f2sum=f2sum+fsymb[i]*fsymb[i]/162.0;
            //            printf("%d %f\n",i,fsymb[i]);
        }
        fac=sqrt(f2sum-fsum*fsum);
        for (i=0; i<162; i++) {
            fsymb[i]=128*fsymb[i]/fac;
            if( fsymb[i] > 127)
                fsymb[i]=127.0;
            if( fsymb[i] < -128 )
                fsymb[i]=-128.0;
            symbols[i]=fsymb[i]+128;
            //            printf("symb: %lu %d\n",i, symbols[i]);
        }
        return;
        }

    return;

}

void unpack50( signed char *dat, int32_t *n1, int32_t *n2 )
{
    int32_t i,i4;
    
    i=dat[0];
    i4=i&255;
    *n1=i4<<20;
    
    i=dat[1];
    i4=i&255;
    *n1=*n1+(i4<<12);
    
    i=dat[2];
    i4=i&255;
    *n1=*n1+(i4<<4);
    
    i=dat[3];
    i4=i&255;
    *n1=*n1+((i4>>4)&15);
    *n2=(i4&15)<<18;
    
    i=dat[4];
    i4=i&255;
    *n2=*n2+(i4<<10);

    i=dat[5];
    i4=i&255;
    *n2=*n2+(i4<<2);
    
    i=dat[6];
    i4=i&255;
    *n2=*n2+((i4>>6)&3);
}

void unpackcall( int32_t ncall, char *call )
{
    char c[]={'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z',' '};
    int32_t n;
    int i,j;
    char tmp[7];
    
    n=ncall;
    
    strcpy(call,"......");
    if (n < 262177560 ) {
        i=n%27+10;
        tmp[5]=c[i];
        n=n/27;
        i=n%27+10;
        tmp[4]=c[i];
        n=n/27;
        i=n%27+10;
        tmp[3]=c[i];
        n=n/27;
        i=n%10;
        tmp[2]=c[i];
        n=n/10;
        i=n%36;
        tmp[1]=c[i];
        n=n/36;
        i=n;
        tmp[0]=c[i];
        j=0;
        tmp[6]='\0';
// remove leading whitespace
        for(i=0; i<5; i++) {
            if( tmp[i] != c[36] )
                break;
        }
        sprintf(call,"%-6s",&tmp[i]);
// remove trailing whitespace
        for(i=0; i<6; i++) {
            if( call[i] == c[36] ) {
                call[i]='\0';
            }
        }
    }
    
}

void unpackgrid( int32_t ngrid, char *grid)
{
    char c[]={'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z',' '};
    int dlat, dlong;
    
    ngrid=ngrid>>7;
    if( ngrid < 32400 ) {
        dlat=(ngrid%180)-90;
        dlong=(ngrid/180)*2 - 180 + 2;
        if( dlong < -180 )
            dlong=dlong+360;
        if( dlong > 180 )
            dlong=dlong+360;
        int nlong = 60.0*(180.0-dlong)/5.0;
        int n1 = nlong/240;
        int n2 = (nlong - 240*n1)/24;
        int n3 = nlong -40*n1 - 24*n2;
        grid[0] = c[10+n1];
        grid[2]=  c[n2];

        int nlat = 60.0*(dlat+90)/2.5;
        n1 = nlat/240;
        n2 = (nlat-240*n1)/24;
        n3 = nlong - 240*n1 - 24*n2;
        grid[1]=c[10+n1];
        grid[3]=c[n2];
    } else {
        strcpy(grid,"XXXX");
    }
}

void unpackpfx( int32_t nprefix, char *call)
{
    char i, nc, pfx[4]="", tmpcall[7]="";
    int32_t n;
    
    strcpy(tmpcall,call);

    if( nprefix < 60000 ) {
// add a prefix of 1 to 3 characters
        n=nprefix;
        for (i=2; i>=0; i--) {
            nc=n%37;
            if( nc >= 0 & nc <= 9 ) {
                pfx[i]=nc+48;
            }
            else if( nc >=10 & nc <= 35 ) {
                pfx[i]=nc+55;
            }
            else {
                pfx[i]=' ';
            }
            n=n/37;
        }

        strcpy(call,pfx);
        strncat(call,"/",1);
        strncat(call,tmpcall,strlen(tmpcall));

    }
    else
// add a suffix of 1 or 2 characters
    {
        nc=nprefix-60000;
        if( nc >= 0 & nc <= 9 ) {
            pfx[0]=nc+48;
            strcpy(call,tmpcall);
            strncat(call,"/",1);
            strncat(call,pfx,1);
        }
        else if( nc >= 10 & nc <= 35) {
            pfx[0]=nc+55;
            strcpy(call,tmpcall);
            strncat(call,"/",1);
            strncat(call,pfx,1);
        }
        else if( nc >= 36 & nc <= 125 ) {
            pfx[0]=(nc-26)/10;
            pfx[1]=(nc-26)%10;
            strcpy(call,tmpcall);
            strncat(call,"/",1);
            strncat(call,pfx,2);
        }
    }

}
void deinterleave(unsigned char *sym)
{
    unsigned char tmp[162];
    unsigned char p, i, j;
    
    p=0;
    i=0;
    while (p<162) {
        j=((i * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
        if (j < 162 ) {
            tmp[p]=sym[j];
            p=p+1;
        }
        i=i+1;
    }
    for (i=0; i<162; i++)
        sym[i]=tmp[i];
}

// used by qsort
int floatcomp(const void* elem1, const void* elem2)
{
    if(*(const float*)elem1 < *(const float*)elem2)
        return -1;
    return *(const float*)elem1 > *(const float*)elem2;
}

void usage(void)
{
    printf("Usage: k9an-wsprd [options...] infile\n");
    printf("       infile must have suffix .wav or .c2\n");
    printf("\n");
    printf("Options:\n");
    printf("       -f x (x is transceiver dial frequency in MHz)\n");
// blanking is not yet implemented. The options are accepted for compatibility
// with development version of wsprd.
//    printf("       -t n (n is blanking duration in milliseconds)\n");
//    printf("       -b n (n is pct of time that is blanked)\n");
    printf("       -H do not use (or update) the hash table\n");
    printf("       -n write noise estimates to file noise.dat\n");
    printf("       -q quick mode - doesn't dig deep for weak signals\n");
    printf("       -v verbose mode\n");
    printf("       -w wideband mode - decode signals within +/- 150 Hz of center\n");
}



int main(int argc, char *argv[])
{
    extern char *optarg;
    extern int optind;
    long int i,j,k;
    unsigned char *symbols, *decdata;
    signed char message[]={-9,13,-35,123,57,-39,64,0,0,0,0};
    char *callsign,*grid,*grid6, *call_loc_pow, *cdbm;
    char *ptr_to_infile,*ptr_to_infile_suffix;
    char uttime[5],date[7];
    int c, delta, nfft2=65536, verbose=0, quickmode=0, writenoise=0, usehashtable=1;
    int shift1, lagmin, lagmax, lagstep, worth_a_try, not_decoded, nadd, ndbm;
    int32_t n1, n2, n3;
    unsigned int nbits;
    unsigned long npoints, metric, maxcycles, cycles;
    float df=375.0/256.0/2;
    float freq0[200],snr0[200],drift0[200],sync0[200];
    int shift0[200];
    float dt=1.0/375.0;
    float dialfreq_cmdline=0.0, dialfreq;
    float fmin=-110, fmax=110;
    float f1, fstep, sync1, drift1, tblank, fblank;

    double *idat, *qdat;
    double mi, mq, mi2, mq2, miq;

    fftw_complex *fftin, *fftout;
    fftw_plan MYPLAN;
#include "./mettab.c"

    idat=malloc(sizeof(double)*nfft2);
    qdat=malloc(sizeof(double)*nfft2);
    
    while ( (c = getopt(argc, argv, "b:f:Hnqt:wv")) !=-1 ) {
        switch (c) {
            case 'b':
                fblank = strtof(optarg,NULL);
                break;
            case 'f':
                dialfreq_cmdline = strtof(optarg,NULL);
                break;
            case 'H':
                usehashtable = 0;
                break;
            case 'n':
                writenoise = 1;
                break;
            case 'q':
                quickmode = 1;
                break;
            case 't':
                tblank = strtof(optarg,NULL);
                break;
            case 'v':
                verbose = 1;
                break;
            case 'w':
                fmin=-150.0;
                fmax=150.0;
                break;
            case '?':
                usage();
                return 1;
        }
    }

    if( optind+1 > argc) {
        usage();
        return 1;
    } else {
        ptr_to_infile=argv[optind];
    }
    
    if( strstr(ptr_to_infile,".wav") ) {
        ptr_to_infile_suffix=strstr(ptr_to_infile,".wav");
        npoints=readwavfile(ptr_to_infile, idat, qdat);
        if( npoints == 1 ) {
            return 1;
        }
        dialfreq=dialfreq_cmdline;
    }    else if ( strstr(ptr_to_infile,".c2") !=0 )  {
        ptr_to_infile_suffix=strstr(ptr_to_infile,".c2");
        npoints=readc2file(ptr_to_infile, idat, qdat, &dialfreq);
        if( npoints == 1 ) {
            return 1;
        }
    }   else {
        printf("Error: infile must have suffix .wav or .c2\n");
        return 1;
    }

// parse the date and time from the given filename
    strncpy(date,ptr_to_infile_suffix-11,6);
    strncpy(uttime,ptr_to_infile_suffix-4,4);
    date[6]='\0';
    uttime[4]='\0';
    if( verbose )
        printf("date %s uttime %s dialfreq %f \n",date, uttime, dialfreq);
    
    getStats(idat, qdat, npoints, &mi, &mq, &mi2, &mq2, &miq);
    if( verbose )
        printf("total power: %5.1f dB\n",10*log10(mi2+mq2));

// Do ffts over 2 symbols, stepped by half symbols
    int nffts=4*floor(npoints/512)-1;
    fftin=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*512);
    fftout=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*512);
    MYPLAN = fftw_plan_dft_1d(512, fftin, fftout, FFTW_FORWARD, FFTW_ESTIMATE);
    
    float ps[512][nffts];
    memset(ps,0.0, sizeof(float)*512*nffts);
    for (i=0; i<nffts; i++) {
        for(j=0; j<512; j++ ){
            k=i*128+j;
            fftin[j][0]=idat[k];
            fftin[j][1]=qdat[k];
        }
        fftw_execute(MYPLAN);
        for (j=0; j<512; j++ ){
            k=j+256;
            if( k>511 )
                k=k-512;
            ps[j][i]=fftout[k][0]*fftout[k][0]+fftout[k][1]*fftout[k][1];
        }
    }
    
    fftw_destroy_plan(MYPLAN);
    fftw_free(fftin);
    fftw_free(fftout);

// get average spectrum
    float psavg[512];
    memset(psavg,0.0, sizeof(float)*512);
    for (i=0; i<nffts; i++) {
        for (j=0; j<512; j++) {
            psavg[j]=psavg[j]+ps[j][i];
        }
    }
    
// smooth with 7-point window and limit the spectrum to +/-150 Hz
    int window[7]={1,1,1,1,1,1,1};
    float smspec[411];
    for (i=0; i<411; i++) {
        smspec[i]=0.0;
        for(j=-3; j<=3; j++) {
            k=256-205+i+j;
            smspec[i]=smspec[i]+window[j+3]*psavg[k];
        }
    }

// sort spectrum values so that we can pick off noise level as a percentile
    float tmpsort[411];
    for (j=0; j<411; j++)
        tmpsort[j]=smspec[j];
    qsort(tmpsort, 411, sizeof(float), floatcomp);

// noise level of spectrum is estimated as 123/411= 30'th percentile
// of the smoothed spectrum. Matched filter sidelobes from very strong signals
// will cause the estimated noise level to be biased high, causing estimated
// snrs to be biased low.
// my SNRs differ from wspr0/wspr-x when there are strong signals in the band.
// This suggests that K1JT's approach and mine have different biases.
    
    float noise_level = tmpsort[122];
    
    if( verbose ) {
        printf("noise level: %5.1f dB\n",10*log10(noise_level)-62.0);
    }
    
    if( writenoise ) {
        FILE *fnoise;
        fnoise=fopen("noise.dat","a");
        fprintf(fnoise,"%s %s %f %5.1f %5.1f\n",date,uttime,dialfreq,10*log10(noise_level)-62.0,10*log10(mi2+mq2));
        fclose(fnoise);
    }

// renormalize spectrum so that (large) peaks represent an estimate of snr
    float min_snr_neg33db = pow(10.0,(-33+26.5)/10.0);
    for (j=0; j<411; j++) {
        smspec[j]=smspec[j]/noise_level - 1.0;
        if( smspec[j] < min_snr_neg33db)
            smspec[j]=0.1;
            continue;
    }

// find all local maxima in smoothed spectrum.
    for (i=0; i<200; i++) {
        freq0[i]=0.0;
        snr0[i]=0.0;
        drift0[i]=0.0;
        shift0[i]=0;
        sync0[i]=0.0;
    }
    int npk=0;
    for(j=1; j<410; j++) {
        if((smspec[j]>smspec[j-1]) && (smspec[j]>smspec[j+1]) && (npk<200)) {
            freq0[npk]=(j-205)*df;
            snr0[npk]=10*log10(smspec[j])-26.5;
            npk++;
        }
    }

// let's not waste time on signals outside of the range [fmin,fmax].
    i=0;
    for( j=0; j<npk; j++) {
        if( freq0[j] >= fmin && freq0[j] <= fmax ) {
            freq0[i]=freq0[j];
            snr0[i]=snr0[j];
            i++;
        }
    }
    npk=i;
    
/* do coarse estimates of freq, drift and shift using k1jt's basic approach, 
more or less.

- look for time offsets of up to +/- 8 symbols relative to the nominal start
 time, which is 2 seconds into the file
- this program calculates shift relative to the beginning of the file
- negative shifts mean that the signal started before the beginning of the 
 file
- The program prints shift-2 seconds to give values consistent with K1JT's
definition
- shifts that cause the sync vector to fall off of an end of the data vector 
 are accommodated by "partial decoding", such that symbols that cannot be 
 decoded due to missing data produce a soft-decision symbol value of 128 
 for the fano decoder.
- the frequency drift model is linear, deviation of +/- drift/2 over the
  span of 162 symbols, with deviation equal to 0 at the center of the signal
  vector (should be consistent with K1JT def). 
*/


    int idrift,ifr,if0,ifd,k0;
    long int kindex;
    float smax,ss,pow,p0,p1,p2,p3;
    for(j=0; j<npk; j++) {
        smax=-1e30;
        if0=freq0[j]/df+256;
        for (ifr=if0-1; ifr<=if0+1; ifr++) {
        for( k0=-10; k0<22; k0++)
        {
            for (idrift=-4; idrift<=4; idrift++)
            {
                ss=0.0;
                pow=0.0;
                for (k=0; k<162; k++)
                {
                    ifd=ifr+((float)k-81.0)/81.0*( (float)idrift )/(2.0*df);
                    kindex=k0+2*k;
                    if( kindex < nffts ) {
                    p0=ps[ifd-3][kindex];
                    p1=ps[ifd-1][kindex];
                    p2=ps[ifd+1][kindex];
                    p3=ps[ifd+3][kindex];
                    ss=ss+(2*pr3[k]-1)*(p3+p1-p0-p2);
                    pow=pow+p0+p1+p2+p3;
                    sync1=ss/pow;
                    }
                }
                if( sync1 > smax ) {
                    smax=sync1;
                    shift0[j]=128*(k0+1);
                    drift0[j]=idrift;
                    freq0[j]=(ifr-256)*df;
                    sync0[j]=sync1;
                }
//                printf("drift %d  k0 %d  sync %f\n",idrift,k0,smax);
            }
        }
        }
        if ( verbose ) {
            printf("npk %2ld snr %6.1f freq %6.1f drift %4.1f shift %5d sync %4.2f\n",
                   j,snr0[j],freq0[j],drift0[j],shift0[j],sync0[j]); }
    }

    nbits=81;
    symbols=malloc(sizeof(char)*nbits*2);
    memset(symbols,0,sizeof(char)*nbits*2);
    decdata=malloc((nbits+7)/8);
    grid=malloc(sizeof(char)*5);
    grid6=malloc(sizeof(char)*7);
    callsign=malloc(sizeof(char)*13);
    call_loc_pow=malloc(sizeof(char)*23);
    cdbm=malloc(sizeof(char)*3);
    float allfreqs[npk];
    memset(allfreqs,0,sizeof(float)*npk);
    char allcalls[npk][13];
    memset(allcalls,0,sizeof(char)*npk*13);
    memset(grid,0,sizeof(char)*5);
    memset(grid6,0,sizeof(char)*6);
    memset(callsign,0,sizeof(char)*13);
    memset(call_loc_pow,0,sizeof(char)*23);
    memset(cdbm,0,sizeof(char)*3);
    char hashtab[32768][13];
    memset(hashtab,0,sizeof(char)*32768*13);
    uint32_t nhash( const void *, size_t, uint32_t);
    int nh;
    
    FILE *fall_wspr, *fwsprd, *fhash;
    
    fall_wspr=fopen("ALL_WSPR.TXT","a");
    fwsprd=fopen("wsprd.out","w");
    
    if( usehashtable ) {
        char line[80], hcall[12];
        if( (fhash=fopen("hashtable.txt","r+")) ) {
           while (fgets(line, sizeof(line), fhash) != NULL) {
              sscanf(line,"%d %s",&nh,hcall);
              strcpy(*hashtab+nh*13,hcall);
           }
        } else {
            fhash=fopen("hashtable.txt","w");
            fclose(fhash);
        }
    }
    
    int uniques=0, noprint=0;
    
    for (j=0; j<npk; j++) {

// now refine the estimates of freq, shift using sync as a metric.
// sync is calculated such that it is a float taking values in the range
// [0.0,1.0].
        
// function sync_and_demodulate has three modes of operation
// mode is the last argument:
//      0 no frequency or drift search. find best time lag.
//      1 no time lag or drift search. find best frequency.
//      2 no frequency or time lag search. calculate soft-decision symbols
//        using passed frequency and shift.
        f1=freq0[j];
        drift1=drift0[j];
        shift1=shift0[j];
        sync1=sync0[j];
        if( verbose ) {
            printf("coarse   : %7.2f %5d %4.1f %6.3f\n",f1,shift1,drift1,sync1); }

// fine search for best sync lag (mode 0)
        fstep=0.0;
        lagmin=shift1-144;
        lagmax=shift1+144;
        lagstep=8;
        sync_and_demodulate(idat, qdat, npoints, symbols, &f1, fstep, &shift1, lagmin, lagmax, lagstep, &drift1, &sync1, 0);
        if( verbose ) {
            printf("fine sync: %7.2f %5d %4.1f %6.3f\n",f1,shift1,drift1,sync1); }

// fine search for frequency peak (mode 1)
        fstep=0.1;
        sync_and_demodulate(idat, qdat, npoints, symbols, &f1, fstep, &shift1, lagmin, lagmax, lagstep, &drift1, &sync1, 1);
        if( verbose ) {
            printf("fine freq: %7.2f %5d %4.1f %6.3f\n",f1,shift1,drift1,sync1); }

        if( sync1 > 0.2 ) {
            worth_a_try = 1;
        } else {
            worth_a_try = 0;
        }

        int idt=0, ii, jiggered_shift;
        uint32_t ihash;
        delta=50;
        maxcycles=10000;
        not_decoded=1;

        while ( worth_a_try && not_decoded && idt<=128) {
            ii=(idt+1)/2;
            if( idt%2 == 1 ) ii=-ii;
            jiggered_shift=shift1+ii;
        
// use mode 2 to get soft-decision symbols
            sync_and_demodulate(idat, qdat, npoints, symbols, &f1, fstep, &jiggered_shift, lagmin, lagmax, lagstep, &drift1, &sync1, 2);
        
            deinterleave(symbols);

            not_decoded = fano(&metric,&cycles,decdata,symbols,nbits,mettab,delta,maxcycles);
            
            idt++;

            if( quickmode ) {
                break;
            }
        }
        
        if( worth_a_try && !not_decoded )
        {
            for(i=0; i<11; i++) {
                if( decdata[i]>127 )
                {
                    message[i]=decdata[i]-256;
                } else {
                    message[i]=decdata[i];
                }
            }

            unpack50(message,&n1,&n2);
            unpackcall(n1,callsign);
            unpackgrid(n2, grid);
            int ntype = (n2&127) - 64;

//            printf("%s\n",callsign);
//            printf("%s\n",grid);
//            printf("%2d\n",ntype);

// Type 1: 6 digit call, grid, power - ntype is positive and is a member of the set {0,3,7,10,13,17,20...60}
// Type 2: extended callsign, power - ntype is positive but not
//         a member of the set of allowed powers
// Type 3: hash, 6 digit grid, power - ntype is negative.

            if( (ntype >= 0) && (ntype <= 62) ) {
                int nu=ntype%10;
                if( nu == 0 || nu == 3 || nu == 7 ) {
                    ndbm=ntype;
                    memset(call_loc_pow,0,sizeof(char)*23);
                    sprintf(cdbm,"%2d",ndbm);
                    strncat(call_loc_pow,callsign,strlen(callsign));
                    strncat(call_loc_pow," ",1);
                    strncat(call_loc_pow,grid,4);
                    strncat(call_loc_pow," ",1);
                    strncat(call_loc_pow,cdbm,2);
                    strncat(call_loc_pow,"\0",1);
                    
                    ihash=nhash(callsign,strlen(callsign),(uint32_t)146);
                    strcpy(*hashtab+ihash*13,callsign);

                    noprint=0;
                } else {
                    nadd=nu;
                    if( nu > 3 ) nadd=nu-3;
                    if( nu > 7 ) nadd=nu-7;
                    n3=n2/128+32768*(nadd-1);
                    unpackpfx(n3,callsign);
                    ndbm=ntype-nadd;

                    memset(call_loc_pow,0,sizeof(char)*23);
                    sprintf(cdbm,"%2d",ndbm);
                    strncat(call_loc_pow,callsign,strlen(callsign));
                    strncat(call_loc_pow," ",1);
                    strncat(call_loc_pow,cdbm,2);
                    strncat(call_loc_pow,"\0",1);

                    noprint=0;
                }
            } else if ( ntype < 0 ) {
                ndbm=-(ntype+1);
                strncat(grid6,callsign+5,1);
                strncat(grid6,callsign,5);
                ihash=(n2-ntype-64)/128;
                if( strncmp(hashtab[ihash],"\0",1) != 0 ) {
                    sprintf(callsign,"<%s>",hashtab[ihash]);
                } else {
                    sprintf(callsign,"%5s","<...>");
                }

                memset(call_loc_pow,0,sizeof(char)*23);
                sprintf(cdbm,"%2d",ndbm);
                strncat(call_loc_pow,callsign,strlen(callsign));
                strncat(call_loc_pow," ",1);
                strncat(call_loc_pow,grid6,strlen(grid6));
                strncat(call_loc_pow," ",1);
                strncat(call_loc_pow,cdbm,2);
                strncat(call_loc_pow,"\0",1);
                
                noprint=0;
                
// I don't know what to do with these... They show up as "A000AA" grids.
                if( ntype == -64 ) {
                    noprint=1;
                }
            }
            
// de-dupe using callsign and freq (only a dupe if freqs are within 1 Hz
            int dupe=0;
            for (i=0; i<npk; i++) {
                if( !strcmp(callsign,allcalls[i]) && (fabs( f1-allfreqs[i] ) < 1.0) )
                    dupe=1;
            }
            if( (verbose || !dupe) && !noprint) {
                uniques++;
                strcpy(allcalls[uniques],callsign);
                allfreqs[uniques]=f1;
                printf("%4s %3.0f %4.1f %10.6f %2d  %-s\n",
                   uttime, snr0[j],(shift1*dt-2.0), dialfreq+(1500+f1)/1e6,
                   (int)drift1, call_loc_pow);
                fprintf(fall_wspr,"%6s %4s %3.0f %3.0f %4.1f %10.6f  %-22s %2d %5lu %4d\n",
                       date,uttime,sync1*10,snr0[j],
                       shift1*dt-2.0, dialfreq+(1500+f1)/1e6,
                       call_loc_pow, (int)drift1, cycles/81, ii);
                fprintf(fwsprd,"%6s %4s %3d %3.0f %4.1f %10.6f  %-22s %2d %5lu %4d\n",
                        date,uttime,(int)(sync1*10),snr0[j],
                        shift1*dt-2.0, dialfreq+(1500+f1)/1e6,
                        call_loc_pow, (int)drift1, cycles/81, ii);
            }
        }
    }
    printf("<DecodeFinished>\n");
    fclose(fall_wspr);
    fclose(fwsprd);

    if( usehashtable ) {
       fclose(fhash);
       fhash=fopen("hashtable.txt","w");
       for (i=0; i<32768; i++) {
           if( strncmp(hashtab[i],"\0",1) != 0 ) {
               fprintf(fhash,"%5ld %s\n",i,*hashtab+i*13);
           }
       }
       fclose(fhash);
    }

    return 0;

}