// Andrew Greensted - Feb 2010
// http://www.labbookpages.co.uk
// Version 1

#include <stdio.h>
#include <stdlib.h>
#include <sndfile.h>
#include <math.h>
#include <strings.h>
#include <fftw3.h>

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
        printf("file length is only %lu seconds\n",npoints/12000);
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
    // 2 no frequency or time lag search. find best drift.

    float dt=1.0/375.0, df=375.0/256.0;
    long int i, j, k;
    double pi=4.*atan(1.0);
    float f0=0.0,fp,ss;
    int lag;
    
    double
    i0[162],q0[162],
    i1[162],q1[162],
    i2[162],q2[162],
    i3[162],q3[162];
    
    double p0,p1,p2,p3,cmet,totp,covmax,pmax, phase, fac;
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
    float fsymb[162];
    
    df=375.0/256.0;
    int best_shift = 0, ifreq;
    
    covmax=-1e30;
    pmax=-1e30;

    int ifmin=0, ifmax=0;

// mode is the last argument:
// 0 no frequency or drift search. find best time lag.
// 1 no time lag or drift search. find best frequency.
// 2 no frequency or time lag search. find best drift.
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
        best_shift = *shift1;
        f0=*f1;
    }

    if( mode != 2 ) {
    for(ifreq=ifmin; ifreq<=ifmax; ifreq++)
    {
        f0=*f1+ifreq*fstep;
// search lag range
        
        for(lag=lagmin; lag<=lagmax; lag=lag+lagstep)
        {
            ss=0;
            totp=0;
            for (i=0; i<162; i++)
            {
                fp = f0 + (*drift1/2.0)*(i-81)/81.0;
                
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
            }
            
            if( ss > covmax ) {
                covmax=ss;
                pmax=totp;

                best_shift=lag;
                
                pmax=totp;
                *f1=f0;
            }
        } // lag loop
        

    } //freq loop
    *sync=covmax/pmax;
    *shift1=best_shift;
        return;
    } //if not mode 2
    
    if( mode == 2 )
    {
//    printf("fbest: %f t0: %f\n", *f1, best_shift*dt);

        for (i=0; i<162; i++)
        {
            i0[i]=0.0;
            q0[i]=0.0;
            i1[i]=0.0;
            q1[i]=0.0;
            i2[i]=0.0;
            q2[i]=0.0;
            i3[i]=0.0;
            q3[i]=0.0;
            
            fp=f0+(*drift1/2.0)*(i-81.0)/81.0;
            for (j=0; j<256; j++)
            {
                k=best_shift+i*256+j;
                if( (k>0) & (k<np) ) {
                phase=2*pi*(fp-1.5*df)*k*dt;
                i0[i]=i0[i]+id[k]*cos(phase)+qd[k]*sin(phase);
                q0[i]=q0[i]-id[k]*sin(phase)+qd[k]*cos(phase);
                phase=2*pi*(fp-0.5*df)*k*dt;
                i1[i]=i1[i]+id[k]*cos(phase)+qd[k]*sin(phase);
                q1[i]=q1[i]-id[k]*sin(phase)+qd[k]*cos(phase);
                phase=2*pi*(fp+0.5*df)*k*dt;
                i2[i]=i2[i]+id[k]*cos(phase)+qd[k]*sin(phase);
                q2[i]=q2[i]-id[k]*sin(phase)+qd[k]*cos(phase);
                phase=2*pi*(fp+1.5*df)*k*dt;
                i3[i]=i3[i]+id[k]*cos(phase)+qd[k]*sin(phase);
                q3[i]=q3[i]-id[k]*sin(phase)+qd[k]*cos(phase);
                }
            }
            
            p0=i0[i]*i0[i]+q0[i]*q0[i];
            p1=i1[i]*i1[i]+q1[i]*q1[i];
            p2=i2[i]*i2[i]+q2[i]*q2[i];
            p3=i3[i]*i3[i]+q3[i]*q3[i];

            
            if( pr3[i] == 1 )
            {
                fsymb[i]=(p3-p1);
            } else {
                fsymb[i]=(p2-p0);
            }
            
        }
        float fsum=0.0, f2sum=0.0;
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
//            printf("symb: %lu %5.1f\n",i, fsymb[i]);
        }
    }
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
    char tmp[6];
    
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
// remove leadig whitespace
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
    printf("K9AN wsprd v0.1\n");
    printf("Usage:\n");
    printf("wsprd -f 14.0956 140710_1822.wav\n");
}



int main(int argc, char *argv[])
{
    long int i,j,k;
    double *idat, *qdat;
    double mi, mq, mi2, mq2, miq;
    unsigned char *symbols, *decdata;
    signed char message[]={-9,13,-35,123,57,-39,64,0,0,0,0};
    int32_t n1,n2;
    char *callsign,*grid;
    int ierr, delta;
    unsigned int nbits;
    unsigned long npoints, metric, maxcycles, cycles;
    float df=375.0/256.0/2;
    float freq0[200],snr0[200],drift0[200],shift0[200];
    int nfft2=65536;
    char *ptr_to_infile,*ptr_to_infile_suffix;
    char uttime[5],date[7];
    float dt=1.0/375.0;
    float dialfreq_cmdline=0.0, dialfreq;

    fftw_complex *fftin, *fftout;
    fftw_plan MYPLAN;
#include "./mettab.c"


    idat=malloc(sizeof(double)*nfft2);
    qdat=malloc(sizeof(double)*nfft2);

    if( argc ==2 && (strstr(argv[1],".wav") || strstr(argv[1],".c2") ) ) {
        ptr_to_infile=argv[1];
    } else if ( argc == 4 && !strcmp(argv[1],"-f") ) {
        ptr_to_infile=argv[3];
        dialfreq_cmdline=strtof(argv[2],NULL);
    } else {
        usage();
        return 1;
    }
    
    if( strstr(ptr_to_infile,".wav") ) {
        ptr_to_infile_suffix=strstr(ptr_to_infile,".wav");
        npoints=readwavfile(ptr_to_infile, idat, qdat);
        dialfreq=dialfreq_cmdline;
    } else {
        ptr_to_infile_suffix=strstr(ptr_to_infile,".c2");
        npoints=readc2file(ptr_to_infile, idat, qdat, &dialfreq);
    }

// parse the date and time from the given filename
    strncpy(date,ptr_to_infile_suffix-11,6);
    strncpy(uttime,ptr_to_infile_suffix-4,4);
    date[6]='\0';
    uttime[4]='\0';
//    printf("date %s uttime %s dialfreq %f \n",date, uttime, dialfreq);
    
    getStats(idat, qdat, npoints, &mi, &mq, &mi2, &mq2, &miq);
//    printf("total power: %4.1f dB\n",10*log10(mi2+mq2));

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

    float psavg[512];
    memset(psavg,0.0, sizeof(float)*512);
    for (i=0; i<nffts; i++) {
        for (j=0; j<512; j++) {
            psavg[j]=psavg[j]+ps[j][i];
        }
    }
    
// smooth with 7-point window
    int window[7]={1,1,1,1,1,1,1};
    float smspec[411];
    for (i=0; i<411; i++) {
        smspec[i]=0.0;
        for(j=-3; j<=3; j++) {
            k=256-205+i+j;
            smspec[i]=smspec[i]+window[j+3]*psavg[k];
        }
    }

    float tmpsort[411];
    for (j=0; j<411; j++)
        tmpsort[j]=smspec[j];
    qsort(tmpsort, 411, sizeof(float), floatcomp);

// noise level of spectrum is estimated as 123/411= 30'th percentile
// of the smoothed spectrum. on a very crowded band, estimated noise level will be
// too high, reducing estimated snr's.
// need to investigate why my SNRs differ from wspr/wspr-x when there are strong
// signals in the band. One of us is biased.

    float noise_level = tmpsort[122];
    
    for (j=0; j<411; j++) {
        smspec[j]=smspec[j]/noise_level - 1.0;
        if( smspec[j] < pow(10.0,(-33+26.5)/10))
            smspec[j]=0.1;
            continue;
    }

/* find all local maxima in smoothed spectrum. this is looser/lazier
 than K1JT's algorithm */
    
    for (i=0; i<200; i++) {
        freq0[i]=0.0;
        snr0[i]=0.0;
        drift0[i]=0.0;
        shift0[i]=0.0;
    }
    int npk=0;
    for(j=1; j<410; j++) {
        if((smspec[j]>smspec[j-1]) & (smspec[j]>smspec[j+1])) {
            freq0[npk]=(j-205)*df;
            snr0[npk]=10*log10(smspec[j])-26.5;
            npk++;
        }
    }

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
    float smax, pmax,ss,pow,p0,p1,p2,p3;
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
                    }
                }
                if( ss > smax ) {
                    smax=ss;
                    pmax=pow;
                    shift0[j]=128*(k0+1)*dt;
                    drift0[j]=idrift;
                    freq0[j]=(ifr-256)*df;
                }
//                printf("drift %d  k0 %d  sync %f\n",idrift,k0,ss/pow);
            }
        }
        }
//        printf("npk %d freq %.1f drift %.1f t0 %.1f sync %.2f\n",j,freq0[j],drift0[j],shift0[j],smax/pmax);
    }

    nbits=81;
    symbols=malloc(sizeof(char)*nbits*2);
    memset(symbols,0,sizeof(char)*nbits*2);
    float f1, fstep, sync, drift1;
    int shift1, lagmin, lagmax, lagstep;
    decdata=malloc((nbits+7)/8);
    grid=malloc(sizeof(char)*5);
    callsign=malloc(sizeof(char)*7);
    char allcalls[npk][7];
    memset(allcalls,0,sizeof(char)*npk*7);
    memset(grid,0,sizeof(char)*5);
    memset(callsign,0,sizeof(char)*7);

    FILE *fall_wspr, *fwsprd;
    fall_wspr=fopen("ALL_WSPR.TXT","a");
    fwsprd=fopen("wsprd.out","w");
    
    int uniques=0, noprint=0;
    
    for (j=0; j<npk; j++) {

// now refine the estimates of freq, shift using sync (0<sync<1) as a metric.
// assume that optimization over freq and shift can be done sequentially
// use function sync_and_demodulate - it has three modes of operation:
// mode is the last argument:
// 0 no frequency or drift search. find best time lag.
// 1 no time lag or drift search. find best frequency.
// 2 no frequency or time lag search. calculate soft-decision symbols using passed frequency and shift.
        f1=freq0[j];
        drift1=drift0[j];

// fine search for best sync lag (mode 0)
        fstep=0.0;
        lagmin=shift0[j]/dt-128;
        lagmax=shift0[j]/dt+128;
        lagstep=8;
        sync_and_demodulate(idat, qdat, npoints, symbols, &f1, fstep, &shift1, lagmin, lagmax, lagstep, &drift1, &sync, 0);
//        printf("after demodulate %f %d %f %f\n",f1,shift1,drift1,sync);

// fine search for frequency peak (mode 1)
        fstep=0.1;
        sync_and_demodulate(idat, qdat, npoints, symbols, &f1, fstep, &shift1, lagmin, lagmax, lagstep, &drift1, &sync, 1);
//        printf("after demodulate %f %d %f %f\n",f1,shift1,drift1,sync);

// use mode 2 to get soft-decision symbols
        sync_and_demodulate(idat, qdat, npoints, symbols, &f1, fstep, &shift1, lagmin, lagmax, lagstep, &drift1, &sync, 2);
        
        deinterleave(symbols);

        delta=50;
        maxcycles=10000;

        ierr = fano(&metric,&cycles,decdata,symbols,nbits,mettab,delta,maxcycles);

//        printf("ierr %d metric %d  cycles %d\n",ierr,metric,cycles/81);

        if( !ierr )
        {
            for(i=0; i<11; i++) {
                if( decdata[i]>128 )
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

// this is where the extended messages would be taken care of
            if( (ntype >= 0) && (ntype <= 62) ) {
                int nu=ntype%10;
                if( nu == 0 || nu == 3 || nu == 7 ) {
                    noprint=0;
                } else {
                    noprint=1;
                }
            } else if ( ntype < 0 ) {
                noprint=1;
            }
            
// de-dupe using callsign
            int dupe=0;
            for (i=0; i<npk; i++) {
                if( !strcmp(callsign,allcalls[i]) )
                    dupe=1;
            }
            if( !dupe && !noprint) {
                uniques++;
                strcpy(allcalls[uniques],callsign);
            printf("%4s %3.0f %4.1f %10.6f %2d  %-s %4s %2d\n",
                   uttime, snr0[j],(shift1*dt-2.0), dialfreq+(1500+f1)/1e6,
                   (int)drift1, callsign, grid, ntype);
            fprintf(fall_wspr,"%6s %4s %3.0f %3.0f %4.1f %10.6f  %-s %4s %2d          %2.0f     %lu\n",
                       date,uttime,sync*10,snr0[j],
                       shift1*dt-2.0, dialfreq+(1500+f1)/1e6,
                       callsign, grid, ntype, drift1, cycles/81);
            fprintf(fwsprd,"%6s %4s %3d %3.0f %4.1f %10.6f  %-s %4s %2d          %2d     %lu\n",
                        date,uttime,(int)(sync*10),snr0[j],
                        shift1*dt-2.0, dialfreq+(1500+f1)/1e6,
                        callsign, grid, ntype, (int)drift1, cycles/81);

            }
        }
    }
    printf("<DecodeFinished>\n");
    fclose(fall_wspr);
    fclose(fwsprd);
	return 0;

}