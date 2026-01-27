/*___________________________________________________
 ¦                                                   ¦
 ¦                   test_lib.c                      ¦
 ¦       by Sergio Frasca - ((( 0 ))) Virgo          ¦
 ¦                    March 2000                     ¦
 ¦___________________________________________________¦*/


#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

#include "pss_serv.h"
#include "pss_snag.h"
#include "pss_math.h"
#include "pss_frame.h"
#include "pss_r87.h"
#include "pss_sfc.h"

#include "FrameL.h"

#define DIMEN 1000

time_t tim0,*tim1;
clock_t cl0,*cl1;

void main()
{
	FILE* fid;
	double data[4];
	int prec[4];

	data[0]=123.434;
	data[1]=123.533;
	data[2]=123.67;
	data[3]=123.34;
	prec[0]=0;
	prec[1]=1;
	prec[2]=3;
	prec[3]=10;

	fid=logfile_open("PROVA");

	logfile_comment(fid,"commento 1");
	logfile_par(fid,"ASD",1234.567,0);
	logfile_input(fid,"FILE_IN","CHANNEL 1");
	logfile_output(fid,"FILE_OUT","scratch file");
	logfile_ev(fid,"PYW",4,data,prec);
	logfile_comment(fid,"commento 2");
	logfile_stop(fid);
	logfile_close(fid);
	
	/* maggio 2005

	char ch;
	SFC_ *sfc_;
	//SFC_ *sfc1;
	double *t0,tt0,xmin,xmax,xmed;
	float *vec,*buffer;
	int *nholes,nper;
	long *nztot,*nzeros,*kzeros;
	ALLPERS_ *allpers;
	HOLES_ *holes;
	long len,i;

	tim0=time(0);
	cl0=clock();
	tim1=&tim0;
	cl1=&cl0;

	sfc_=sfc_open("D:\\Data\\pss\\virgo\\sd\\sds\\c5\\VIR_hrec_20041203_111051_.sds");

	sfc_show(sfc_);

	len=1000000;
	vec=malloc(len*sizeof(float));
	buffer=malloc(len*sfc_->nch*sizeof(float));
	t0=malloc(sizeof(double));
	allpers=malloc(sizeof(ALLPERS_));
	holes=malloc(sizeof(HOLES_));

	vec_from_sds(buffer,vec,t0,holes,sfc_,1,len,allpers);

	printf(" %g ** %g ** %g ** %g ** %g ** %g \n",vec[0],vec[1],vec[2],vec[3],vec[4],vec[5]);

	xmin=vec[0];
	xmax=vec[0];
	xmed=vec[0];

	for(i=1; i<len; i++){
		if(vec[i]<xmin)xmin=vec[i];
		if(vec[i]>xmax)xmax=vec[i];
		xmed+=vec[i];
	}

	printf(" min = %g  max = %g  med = %g \d",xmin,xmax,xmed/len);

	printf(" t0 = %f \n",*t0);

	//sfc_=sfc_openw("prova.sfc",sfc_);
	//fclose(sfc_->fid);

	//sfc_=sfc_open("prova.sfc");
	
	//printf("\n *** Caption  :  %s \n\n",sfc_->capt);

	//sfc_show(sfc_);
	
	wait_inp();
	*/
}