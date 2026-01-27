/**pss_sds.c: pss_frameS software: to read data from SDS files and put into a DS the h-reconstructed data**/
/**Last version: Dic, 04, 2009**/
/**First version: Nov. 2009**/
/**Author: Sabrina e Pia**/
// . /opt/exp_software/virgo/VCS-5.1/System/Env.sh
//FrDump -i  /storage/gpfs_virgo4/virgo/data/VSR2/proc/V-HrecOnline-9354* -d 0 >filehrec_.ffl
//  ./compila
//  ./crea_sfdb.out <input_VSR2_H4096Hz
//
#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>

/****** PSS libraries***************/
#include "../pss_lib/pss_math.h"
#include "../pss_lib/pss_serv.h"
#include "../pss_lib/pss_snag.h"
#include "../pss_lib/pss_sfc.h"
/****** Antenna libraries**************/
#include "pss_ante.h"

/******SFDB libraries**************/
#include "pss_sfdb.h"

/******Pia's SDS libraries***************/
#include "pss_sds.h"
#include "pss_frameS.h"
/******Pss Astro  libraries**************/
#include"../pss_astro/pss_astro.h"


//double frame2hds(DS *d_h,INPUT_PARAM *input_param, DETECTOR_PARAM *detector_param,EVEN_PARAM *even_param,EVF_PARAM *evf_param,HEADER_PARAM *header_param, SCIENCE_SEGMENT* science_segment, FrFile *iFileNoise,double  *StartingTime,int *gps_startDB,int *gps_endDB)
double frame2hds(DS *d_h,INPUT_PARAM *input_param, DETECTOR_PARAM *detector_param,EVEN_PARAM *even_param,EVF_PARAM *evf_param,HEADER_PARAM *header_param, SCIENCE_SEGMENT* science_segment, FrFile *iFileNoise,double *StartingTime,double *gps_startDB,double *gps_endDB)

{
  double itt;
  int verb;
  float sum2;
  int i;
  int res;
  
  //FrameH *frame;  //piapia3
  //int iii;
  res=input_param->resampling;
  header_param->gps_sec=0.; //inizialization
  verb=input_param->verb;
  header_param->wink=input_param->wink;
  
  //iFileNoise = FrFileINew(input_param->filename); //pia03--va sistemato
  
  if(d_h->lcr==0){
     /* Open the frame file and create the structure */
     //iFileNoise = FrFileINew(input_param->filename);
     *gps_startDB = FrFileITStart(iFileNoise);
     *gps_endDB = FrFileITEnd(iFileNoise);
     printf("----------------------------\n");
     printf("Approx. GPS_START = %17.6f\n",*gps_startDB);
     printf("Approx.GPS_END = %17.6f\n",*gps_endDB);
     printf("----------------------------\n\n");
     *StartingTime= (double) *gps_startDB;
    //sfc_data=sds_open(input_param->filename,allpers); //allocate sfc_data 
    strcpy(header_param->sfdbname,input_param->sfdbname);
    printf("Sfdb filename, as written in the header param= %s \n", header_param->sfdbname);
  }

  /*
  iii=0;
  frame = FrameRead(iFileNoise);
	while(frame == NULL & iii < 6){
		iii++;
		printf(" *** Strangely null frame ! \n");
		frame = FrameRead(iFileNoise);
		
	}
	if(frame != NULL) printf("****NON NULL******\n");
       
	*/


  if(d_h->lcr==0)detector_param->itf_var_param->nfft=0;
  if(d_h->lcr==0)detector_param->itf_var_param->kdati_tot=0;
  //2012
  if(d_h->lcr==0){
    detector_param->itf_var_param->n_flagappo=0;
  }
   
  detector_param->itf_var_param->n_flag=0; //number in every fft
  itt=frame2dsS(d_h,input_param,verb,detector_param,iFileNoise,StartingTime, gps_startDB,gps_endDB,header_param, science_segment);
 

  d_h->lcr++;  //another chunk has been served


  if(verb>4)printf("FFT number = %d \n",detector_param->itf_var_param->nfft); 
  /* Definition of the header parameters */
  header_param->nfft=detector_param->itf_var_param->nfft; 
  header_param->n_flag=detector_param->itf_var_param->n_flag; 
  header_param->endian=1.0;
  header_param->gps_sec=(int) detector_param->itf_var_param->gpstime;
  header_param->gps_nsec=(int)((detector_param->itf_var_param->gpstime-(int)detector_param->itf_var_param->gpstime)*1.0e9);
  header_param->tbase=detector_param->itf_var_param->fftlen*detector_param->itf_var_param->tsamplu;
  header_param->firstfreqindex=(int) (detector_param->itf_var_param->frinit/detector_param->itf_var_param->deltanu);
  header_param->nsamples=(int) detector_param->itf_var_param->fftlen/2;
  header_param->mjdtime=detector_param->itf_var_param->mjdtime; //ATT: it can be non-continuous. There are holes !!!
  header_param->frinit=detector_param->itf_var_param->frinit;
  header_param->tsamplu=detector_param->itf_var_param->tsamplu;
  header_param->deltanu=detector_param->itf_var_param->deltanu;
  
  header_param->red=input_param->red;        //reduction factor for the very short FFTs
  header_param->typ=input_param->typ;        //0, 1 = non overlapping data  2=overlapping data

  //printf("AIUTO normd header_param->nsamples=%d header_param->tsamplu =%f\n",header_param->nsamples,header_param->tsamplu);
  //printf("AIUTO normd nsample= %d header_param->nsamples res %d \n",header_param->nsamples,res);
  header_param->normd=res*sqrt(header_param->tsamplu/(2*header_param->nsamples)); //piapia 24 oct 2005s
  if(header_param->wink !=5)header_param->normw=1; //to be added the factor for Hamming,Blackman
  if(header_param->wink==1)header_param->normw=sqrt(3.0); //Hanning
  if(header_param->nfft==1 && header_param->wink==5){ 

    sum2=0;
    for (i=0; i <header_param->nsamples/(2*res) ; i++){  //first 1/4 window
      sum2+=pow((0.5-0.5*cos(i*PIG/(header_param->nsamples/(2*res)))),2);
    }
    for (i=header_param->nsamples/(2*res);i<2*header_param->nsamples/(res)-header_param->nsamples/(2*res); i++){
      sum2+=1.;
    }
    for (i=0; i <header_param->nsamples/(2*res) ; i++){  //last 1/4 window = first 1/4 window
      sum2+=pow((0.5-0.5*cos(i*PIG/(header_param->nsamples/(2*res)))),2);
    }	   
    header_param->normw=sqrt(2*header_param->nsamples/(res*sum2));
    

    //flat + cos edges
    /*    sum2=0;
    for (i=0; i <header_param->nsamples/(2*res) ; i++){  //first 1/4 window
      sum2+=pow((0.5-0.5*cos(i*PIG/(header_param->nsamples/(2*res)))),2);
    }
    for (i=header_param->nsamples/(2*res); i<2*header_param->nsamples/res-header_param->nsamples/(2*res) ; i++){
     sum2+=1.;
    }
    for (i=0; i <header_param->nsamples/(2*res) ; i++){  //last 1/4 window = first 1/4 window
      sum2+=pow((0.5-0.5*cos(i*PIG/(header_param->nsamples/(2*res)))),2);
    }	   
    header_param->normw=sqrt(2*header_param->nsamples/(res*sum2));
    printf("res=%d\n",res);*/
  }
  //printf("AIUTO normw nsample= %d header_param->nsamples\n",header_param->nsamples);
  //printf("header_param->normw =%f normd = %f res=%d \n",header_param->normw,header_param->normd,res); 
  if(verb>4)puts("Exit from sds2hds:time domain data"); 
  if(header_param->nfft==1){
    
     /* Definition of parameters for the Log file */
     logfile_comment(LOG_INFO,"par GEN: general parameters of the run");
    
     logfile_comment(LOG_INFO,"   GEN_BEG is the beginning time (mjd)");
     logfile_comment(LOG_INFO,"   GEN_NSAM the number of samples in 1/2 FFT");
     logfile_comment(LOG_INFO,"   GEN_DELTANU the frequency resolution");
     logfile_comment(LOG_INFO,"   GEN_FRINIT the beginning frequency of the FFT");
     logfile_comment(LOG_INFO,"   GEN_TSAM sampling time [s]");
     logfile_par(LOG_INFO,"GEN_BEG",header_param->mjdtime,3);
     logfile_par(LOG_INFO,"GEN_NSAM",header_param->nsamples/res,1);
     logfile_par(LOG_INFO,"GEN_DELTANU",header_param->deltanu,4);
     logfile_par(LOG_INFO,"GEN_FRINIT",header_param->frinit,0);
      logfile_par(LOG_INFO,"GEN_TSAM",header_param->tsamplu,0);
     logfile_comment(LOG_INFO,"   EVT_CR is the threshold");
     logfile_comment(LOG_INFO,"   EVT_TAU the memory time of the AR estimation");
     logfile_comment(LOG_INFO,"   EVT_DEADT the dead time [s]");  
     logfile_comment(LOG_INFO,"   EVT_EDGE seconds purged around the event");
     //    logfile_par(LOG_INFO,"EVT_CR",even_param->cr,0);
     //logfile_par(LOG_INFO,"EVT_TAU",even_param->tau,0);
     //logfile_par(LOG_INFO,"EVT_DEADT",even_param->deadtime,0);
     //logfile_par(LOG_INFO,"EVT_EDGE",even_param->edge,0);
     //printf("evf param tau %f\n",evf_param->tau_Hz);
     //logfile_par(LOG_INFO,"EVF_THR",evf_param->maxdin,0);
     //logfile_par(LOG_INFO,"EVF_TAU",evf_param->tau_Hz,0);
     //logfile_par(LOG_INFO,"EVF_MAXAGE",evf_param->maxage,0);
     //logfile_par(LOG_INFO,"EVF_FAC",evf_param->factor_write,0);
     //logfile_comment(LOG_INFO,"   EVF_THR is the threshold in amplitude");
     //logfile_comment(LOG_INFO,"   EVF_TAU the memory frequency of the AR estimation");
     //logfile_comment(LOG_INFO,"   EVF_MAXAGE [Hz] the max age of the process. If age>maxage the AR is re-evaluated");
     //logfile_comment(LOG_INFO,"   EVF_FAC is the factor for which the threshold is multiplied, to write less EVF in the log file");

   
  if(header_param->deltanu <= 0.001){
    evf_param->tau_Hz=0.02;
    evf_param->maxage=0.02;
    even_param->tau=600;
  }
  else
    {
      evf_param->tau_Hz= 0.02*header_param->deltanu/0.001; 
      evf_param->maxage=evf_param->tau_Hz;
      even_param->tau=600*0.001/header_param->deltanu;
    }

  //Write in the LOG file

  logfile_par(LOG_INFO,"EVT_CR",even_param->cr,0);
  logfile_par(LOG_INFO,"EVT_TAU",even_param->tau,0);
  logfile_par(LOG_INFO,"EVT_DEADT",even_param->deadtime,0);
  logfile_par(LOG_INFO,"EVT_EDGE",even_param->edge,0);
  //printf("-->evaluated evf param tau_Hz %f\n",evf_param->tau_Hz);
  //printf("-->evaluated evt param tau %f\n",even_param->tau);
  logfile_par(LOG_INFO,"GEN_THR",evf_param->maxdin,3);
  logfile_par(LOG_INFO,"GEN_TAU",evf_param->tau_Hz,3);
  logfile_par(LOG_INFO,"GEN_MAXAGE",evf_param->maxage,3);
  logfile_par(LOG_INFO,"GEN_FAC",evf_param->factor_write,3);

  
     logfile_par(LOG_INFO,"SAT_LEV",even_param->sat_level,0);
     logfile_comment(LOG_INFO,"   SAT_LEV is the saturation level");
     logfile_par(LOG_INFO,"SAT_MAX",even_param->sat_maxn,0);
     logfile_comment(LOG_INFO,"   SAT_MAX is the max allowed number of sat data");
  }
   return itt;
}



double frame2dsS(DS *d,INPUT_PARAM *input_param,int verb, DETECTOR_PARAM *detector_param,FrFile *iFileNoise,double *StartingTime,double *gps_startDB, double *gps_endDB, HEADER_PARAM *header_param, SCIENCE_SEGMENT* science_segment)
{
  double itt;
  long len;
  char fflfile[MAXMAXLINE+1];
  float *vec=NULL;
  float *vecF=NULL;
  double tim0=-1; 
  double exp_diff;
  float norm;  //such that the FFT of time data gives the noise spectral amplitude
               //norm=sqrt(dt/fftlen)
  int i;
  int j=0;
  long len_call;
  int len_flag;
  int nleap;  //gpsleap; NOT USED
  
  //Specific frame
  double sampling;
  struct FrVect *VectNoise = NULL;
  struct FrVect *VectFlag = NULL;
  float nzeroes, sum;
  double stepnull=30.0;
  int appo;
  double utc;
  
  // ------ 03/2013 TEST Alberto write debug info
  long len_call2=0;
  //FILE *outFlag=0;
  //outFlag=fopen("test_flag.out","a");
  //int prevFlag=0;
  //int science=0;

  //FILE *debugfile=0;
  //debugfile=fopen("debugfile.out","a");

  //int iVETATI;
  
  len=input_param->len;
  strcpy(fflfile,input_param->fflname);
  //fflfile=input_param->fflname;
  detector_param->itf_var_param->mjdtime=0.; //4/02/08
  detector_param->itf_var_param->gpstime=0.; //4/02/08

  vec=(float *)malloc((size_t) (d->len)*sizeof(float));
  vecF=(float *)malloc((size_t) (d->len)*sizeof(float));
  for (j=0; j<d->len; j++){
  	vec[j]=0.0;
  	vecF[j]=0.0;
  }
  
  
  if(input_param->iopen==0){
    printf("----------------------------\n");
    printf("Approx. GPS_START = %f\n",*gps_startDB);
    printf("Approx.GPS_END = %f\n",*gps_endDB);
    printf("----------------------------\n\n");
    *StartingTime= (double) *gps_startDB;
    printf("--->>>Starting time of the new chunk = %f\n",*StartingTime);
    VectNoise = FrFileIGetVect(iFileNoise,input_param->framechannel,*StartingTime,(double) 1.0); //to get parameters\\tolto DN
    if(VectNoise==NULL){    //sab
      int inull=0;
      while((VectNoise == NULL)&& *StartingTime < *gps_endDB) {
	//FrVectFree(VectNoise);//sab
	if (inull % 10000 == 0) { printf("******** Channel=%s i=%d  no data found at %f ********  \n",input_param->framechannel , inull, *StartingTime);}
	inull +=1;
	*StartingTime+=stepnull;
	VectNoise = FrFileIGetVect(iFileNoise,input_param->framechannel,*StartingTime,(double) 1.0); //to get parameters\\tolto DN
      }
    }//sab

    if(VectNoise == NULL){
       printf("******** No data found in the whole file! ******** \n");
      input_param->iopen=-1;
 
      free(vec);
      free(vecF);
      itt=0.0;
      return itt;
    
    }
    puts("A new run has started\n");
    sampling=VectNoise->dx[0];
    //printf("Frame file parameters= t0 (gps,s) dt (s) %f %f\n",*StartingTime,sampling);
    FrVectFree(VectNoise);
    VectNoise=NULL;
        
    itt=*StartingTime; //gps time
    if(d->type<=1)exp_diff=sampling*d->len/day2sec;
    if(d->type==2)exp_diff=sampling*(d->len/d->type)/day2sec;
    d->cont=0;
    d->dt=sampling;
    d->tini1=*StartingTime/day2sec;
    d->tini2=*StartingTime/day2sec;        
    detector_param->itf_fixed_param->type=d->type;
    if(verb>3) printf("initial time of the FIRST file  (beginning gps) %f \n",itt);
    detector_param->itf_fixed_param->time=itt/day2sec+dot2gps; /*Beginning time of the run in days from 1900*/
    if(verb>3) printf("initial time of the FIRST filein days from 1900 %f \n",detector_param->itf_fixed_param->time);  
  }// end first time
  sampling= d->dt; //for next times 
  if(input_param->iopen==0)input_param->iopen=1;
  //len=d->len;
  d->cont=0; //control variable.
  
  //modifica 17 luglio 2015. Per usare la function. itt=*StartingTime. Comunque non si usa...LEVO TUTTO !!!
  //if(*StartingTime >=(float)24537832.5)nleap=34;
  //if(*StartingTime >=(float)24537832.5){
  //if(*StartingTime >=(float)2453736.5){ //! Jan 2006     
  //  nleap=33;
  //}
  //else {
  //  nleap=32;
  //}
    
  //}
  //nleap=PSS_leapseconds(itt);
  
  
  //gpsleap=nleap-19;

  detector_param->itf_fixed_param->gpstime=*StartingTime;
  //(detector_param->itf_fixed_param->time-dot2gps)*day2sec+gpsleap;/*run gps time*/
  /* Var parameters */ 
  detector_param->itf_var_param->nfft++;
  detector_param->itf_var_param->tsamplu=sampling;
  detector_param->itf_var_param->fftlen=d->len;
  detector_param->itf_var_param->frinit=0.;
  detector_param->itf_var_param->deltanu=1/(d->len*detector_param->itf_var_param->tsamplu);
  if(verb >4)printf("deltanu,d->len %f %ld \n",detector_param->itf_var_param->deltanu,d->len);
  detector_param->itf_var_param->frinit=0.;
  detector_param->itf_var_param->freq_min=0;  //piapia: to be defined as input parameters
  detector_param->itf_var_param->freq_max=detector_param->itf_var_param->freq_min+1/(2*sampling);
  if(verb>3)printf("fixed time,gpstime= %f %f\n",detector_param->itf_fixed_param->time,detector_param->itf_fixed_param->gpstime);
  if(d->type<=1)exp_diff=sampling*len/day2sec;
  if(d->type==2)exp_diff=sampling*(len/d->type)/day2sec;
  if(verb>3) printf("exp_diff from subsequent ffts, in seconds %f \n", exp_diff*day2sec);
  /*Now read the data in the file:*/
  /*vec_from_frame gives the data and then ds_from_vec a ds with the data*/
  /* The "real" time should be read again here, because vec contains zero padding, to remove "short" holes in the data */
  /*IF the holes are "long" there is no zero padding and we must close the SFDB file, hence put iopen=-1, and start with a new one*/
  
  
  norm=sqrt(sampling/d->len); //such that the FFT of time data gives the noise spectral amplitude. pia:CHECK
  len_call=len;
  if(d->type==2&&detector_param->itf_var_param->nfft>1)len_call=len/2;
  //printf("fft numb.(the next to be done) len_call norm= %d %ld %f\n",detector_param->itf_var_param->nfft,len_call,norm);
  // **************************Read the data**********************
  if(d->type==2&&detector_param->itf_var_param->nfft>1){
    
    //VectNoise = FrFileIGetVect(iFileNoise,input_param->framechannel,*StartingTime+len_call*sampling,(double) len_call*sampling);
    if(verb>3) printf("prima di VectNoise filePoint=%d file=%s startingtime=%f len*sampling=%f nfft=%d\n", iFileNoise, iFileNoise->path,*StartingTime+len_call*sampling,(double) len_call*sampling,detector_param->itf_var_param->nfft);
    VectNoise = FrFileIGetVectFN(iFileNoise,input_param->framechannel,*StartingTime+len_call*sampling,(double) len_call*sampling);
    //VectFlag = FrFileIGetVectFN(iFileNoise,input_param->flagchannel,*StartingTime,(double) len*sampling);
    //printf("prima di VectFlag file=%s startingtime=%f ceil(len*sampling)=%f nfft=%d\n", iFileNoise->path,*StartingTime+len_call*sampling,(double) ceil(len*sampling), detector_param->itf_var_param->nfft);
    VectFlag = FrFileIGetVectFN(iFileNoise,input_param->flagchannel,*StartingTime+len_call*sampling,(double) ceil(len_call*sampling));
  }
  
  else
    {
      //VectNoise = FrFileIGetVect(iFileNoise,input_param->framechannel,*StartingTime,(double) len*sampling);
      if(verb>3) printf("prima di VectNoise ifilePoint=%d file=%s startingtime=%f len*sampling=%f\n", iFileNoise, iFileNoise->path,*StartingTime,(double) len*sampling);
      VectNoise = FrFileIGetVectFN(iFileNoise,input_param->framechannel,*StartingTime,(double) len*sampling);
      //printf("Flag channel %s len %f\n",input_param->flagchannel,(double) len);
      //printf("prima di VectFlag file=%s startingtime=%f ceil(len*sampling)=%f\n", iFileNoise->path,*StartingTime,(double) ceil(len*sampling));
      //VectFlag = FrFileIGetVectFN(iFileNoise,input_param->flagchannel,*StartingTime,(double) len*sampling); // Albe se 131.1 NON FUNGE!
      VectFlag = FrFileIGetVectFN(iFileNoise,input_param->flagchannel,*StartingTime,(double) ceil(len*sampling));
      
      //printf("VectFlag\n");
    }
  if(iFileNoise == NULL) printf("Null ffl file!\n");
  if(verb>3) printf("******* FFT number: %d - channel: %s - StartingTime: %f - Chunk lenght (s): %f ******* \n", detector_param->itf_var_param->nfft, input_param->framechannel, *StartingTime, (double) len_call*sampling);
         

  if(VectNoise != NULL)
    {
      if(verb>3) printf("VectNoise->GTime %f\n", VectNoise->GTime);
      if(verb>3) printf("VectNoise->nData %d\n",(int) VectNoise->nData);
      
      //if((int) VectNoise->nData != len)printf("ATTT");
      if(VectNoise->next==NULL){    // No holes --> take data
	if(verb>3) printf("VectNoise->next==NULL !!!!!!!!!!!!!!!!!!\n");
	sum=0., nzeroes=0;
	switch(VectNoise->type)
	  {
	  case FR_VECT_4R:
	    
	    for(i=0; i<len_call; i++)
	      {
		vec[i] = VectNoise->dataF[i];
		
	      }
	    
	    break;
	    
	  case FR_VECT_8R:
	    
	    for(i=0; i<len_call; i++)
	      {
		vec[i] = VectNoise ->dataD[i];
		
	      }
	    break;
	  }
	
	
      } //IF NEXT NULL	 
     else
	{
	  for(i=0; i<len_call; i++)
	    {
	      if(VectNoise->next->data[i]==0)vec[i] = VectNoise->next->data[i];
	      
	      if(VectNoise->next->data[i]!=0){
		switch(VectNoise->type)
		  {
		  case FR_VECT_4R:
		    vec[i] = VectNoise->dataF[i];
		    break;
		    
		  case FR_VECT_8R:
		    vec[i] = VectNoise ->dataD[i];
		    break;
		  }
	      }
	    }  //
	} //matcheselse
         

      //added for the flags, 18 april 2011
      // 2013-02-26 Alberto: Sposto per testare data=0 se flag=bad >>>>>>>>>>>>
      if(VectFlag != NULL)
	{
	  len_flag=(int )VectFlag->nData;
	  //printf("len_Flag %d\n",len_flag);
	  //printf("VectFlag->GTime %f\n", VectFlag->GTime);
	  //printf("VectFlag->nData %d\n",(int) VectFlag->nData);

	if(VectFlag->next==NULL){    // No holes --> take data

	  switch(VectFlag->type)
	    {
	    case FR_VECT_4R:


	      for(i=0; i<len_flag; i++)
		{
		  vecF[i] = VectFlag->dataF[i];

		}

	      break;

	    case FR_VECT_8R:

	      for(i=0; i<len_flag; i++)
		{
		  vecF[i] = VectFlag->dataD[i];

		}
	      break;
	    }

	  if(detector_param->itf_var_param->nfft==1) {len_call2=0;}
	  else{len_call2=len_call;}

	  //printf("len_flag= %d\n",len_flag);
	  //fprintf(outFlag, ">>> New FFT StartingTime=%f length=%f len_flag=%d\n", *StartingTime+len_call2*sampling,ceil(len_call*sampling), len_flag);
	  
	  for(i=0;i<len_flag;i++){
	  //  
	  //  if (vecF[i] != prevFlag) {
	  //  	
// 		if(vecF[i] < input_param->flagstatus ){
// 			science=0;
// 		}	
// 		else {science=1;}
// 		
// 		fprintf(outFlag, "%f \t %f \t %d\n", *StartingTime+len_call2*sampling+i*(ceil(len_call*sampling)/len_flag), vecF[i], science);
// 	    
// 	       prevFlag=vecF[i];
// 	    } 
	    
	    
	    
	    if(vecF[i] != input_param->flagstatus )detector_param->itf_var_param->n_flag+=1;	
	    if(d->type==2&&(i >=len_flag/2)&&detector_param->itf_var_param->nfft==1){
	      detector_param->itf_var_param->n_flagappo+=1;
	    } //2012			

	  } //end del for

	  //2012: overlapping situation
	  if (d->type==2 && detector_param->itf_var_param->nfft !=1){
	    appo=detector_param->itf_var_param->n_flag+detector_param->itf_var_param->n_flagappo;
	    detector_param->itf_var_param->n_flagappo=detector_param->itf_var_param->n_flag;
	    detector_param->itf_var_param->n_flag=appo;
	  }

	if(len_flag!=0)detector_param->itf_var_param->n_flag/=len_flag;
	if(len_flag!=0 && d->type==2 &&detector_param->itf_var_param->nfft!=1)detector_param->itf_var_param->n_flag/=2;	


	} //IF NEXT NULL 

	else
	  {
	    //printf(" NOT NULL Next Flag -> NO DATA");
	    //2012
	    switch(VectFlag->type)
	     {
	    case FR_VECT_4R:


	      for(i=0; i<len_flag; i++)
		{
		  vecF[i] = VectFlag->dataF[i];

		}

	      break;

	    case FR_VECT_8R:

	      for(i=0; i<len_flag; i++)
		{
		  vecF[i] = VectFlag->dataD[i];

		}
	      break;
	     }

		for(i=0;i<len_flag;i++){
		  if(VectFlag->next->data[i]==0 || vecF[i] != input_param->flagstatus ){
		    detector_param->itf_var_param->n_flag+=1;
		  }
		    if(d->type==2&&(i >=len_flag/2)&&detector_param->itf_var_param->nfft==1){
		      detector_param->itf_var_param->n_flagappo+=1;
		    }


		} //fine DEL FOR pezzo 2012
	      	  //2012: overlapping situation
		if (d->type==2 && detector_param->itf_var_param->nfft !=1){
	           appo=detector_param->itf_var_param->n_flag+detector_param->itf_var_param->n_flagappo;
	           detector_param->itf_var_param->n_flagappo=detector_param->itf_var_param->n_flag;
		   detector_param->itf_var_param->n_flag=appo;
		}
		// NO detector_param->itf_var_param->n_flag=-2;
	      	  if(len_flag!=0){
		    detector_param->itf_var_param->n_flag/=(-len_flag);
		    if(d->type==2 &&detector_param->itf_var_param->nfft!=1)detector_param->itf_var_param->n_flag/=2;	
		  }

	  }

	}
      else{ // IF VectFlag == NULL !!!
	if(verb>4) printf("NULL Flag\n");
	detector_param->itf_var_param->n_flag=-1;
      }
      // 2013-02-26 Alberto: Sposto per testare data=0 se flag=bad <<<<<<<<<<<<<
      

      //  if (science_segment != NULL) {
      //	printf("Applying science segment list selection!\n");
      //}
      //  iVETATI=0;  //for test only
      for(i=0;i<len_call;i++){
	
	
	// --------- 03/2013 TEST alberto set data=0 in non science segments (letti da file)
	/*
	if (science_segment != NULL) {
		int vectime = *StartingTime+len_call2*sampling+i*sampling; 
		jjj=0;
	
		while (jjj<science_segment->nsegm-1 && vectime>science_segment->segtstart[jjj+1]){
			jjj++; 
		}
			
		if (vectime>=science_segment->segtstart[jjj]+science_segment->twin && vectime<science_segment->segtstop[jjj]-science_segment->twin) { 
 			//if (i % 1000000 == 0) fprintf(debugfile,"time = %f segm start+twin = %d segm stop-twin = %d ----> in science mode segment n. %d! vec[i]=%e\n", 
			//	*StartingTime+len_call2*sampling+i*sampling, science_segment->segtstart[jjj]+science_segment->twin, science_segment->segtstop[jjj]-science_segment->twin, science_segment->num[jjj], vec[i]);
		} else if (vectime >= science_segment->segtstart[jjj] && vectime < science_segment->segtstart[jjj]+science_segment->twin ) {
 			vec[i] *= ((float) (vectime-science_segment->segtstart[jjj]))/((float)science_segment->twin);			
			//if (i % 10000 == 0) fprintf(debugfile,"time = %f segm start = %d segm start+twin = %d ----> in edge 1 n. %d! --- vectime-segtstart=%d factor = %f vec=%e\n", 
			//	*StartingTime+len_call2*sampling+i*sampling, science_segment->segtstart[jjj], science_segment->segtstart[jjj]+science_segment->twin, science_segment->num[jjj], vectime-science_segment->segtstart[jjj],((float) (vectime-science_segment->segtstart[jjj]))/((float)science_segment->twin), vec[i]);
		} else if (vectime >= science_segment->segtstop[jjj]-science_segment->twin && vectime < science_segment->segtstop[jjj] ) {
			vec[i] *= ((float) (science_segment->segtstop[jjj]-vectime))/((float)science_segment->twin);			
 			//if (i % 10000 == 0) fprintf(debugfile,"time = %f segm stop-twin = %d segm stop = %d ----> in edge 2 n. %d! --- segtstop-vectime=%d factor = %f vec=%e\n", 
			//	*StartingTime+len_call2*sampling+i*sampling, science_segment->segtstop[jjj]-science_segment->twin, science_segment->segtstop[jjj], science_segment->num[jjj], science_segment->segtstop[jjj]-vectime, ((float) (science_segment->segtstop[jjj]-vectime))/((float)science_segment->twin), vec[i]);
		} else {
			vec[i]=0.;	
		}
	}
	*/
	/*
		if (science_segment != NULL) {
		float vectime = *StartingTime+len_call2*sampling+i*sampling;
		int iALL;
		if (i==0) {
		  iALL=0; //to know if the FFT is all inside science segment
		  jjj=0;  //to run faster
		}
		while (jjj<science_segment->nsegm-1 && vectime>science_segment->segtstart[jjj+1]){
			jjj++; 
		}
		if (i==0) {
		    if (vectime>=science_segment->segtstart[jjj] && (vectime + len_call*sampling)<science_segment->segtstop[jjj]) {
		      printf("===> FFT all inside a science period %f %f\n",vectime,vectime+len_call*sampling);
		      iALL=1;
		    }
		  }
		if (iALL==0) {
		if (vectime>=science_segment->segtstart[jjj] && vectime<science_segment->segtstop[jjj]) { 
 			if (i < 2) printf("OK time = %f segm start = %d segm stop = %d ----> in science mode segment n. %d! vec[i]=%e\n", *StartingTime+len_call2*sampling+i*sampling, science_segment->segtstart[jjj], science_segment->segtstop[jjj], science_segment->num[jjj], vec[i]);
		
		} else {
			vec[i]=0.;
			iVETATI+=1;
			if (iVETATI < 2 ) printf("VETOED time = %f segm start = %d segm stop = %d ----> in science mode segment n. %d! vec[i]=%e\n", *StartingTime+len_call2*sampling+i*sampling, science_segment->segtstart[jjj], science_segment->segtstop[jjj], science_segment->num[jjj], vec[i]);
		}
		}  //iALL==0
		}  //Science segment not null
	*/
	
	vec[i]*=(1.0/header_param->einstein);//output data in the ds will be in EINSTEIN units
	//piapia tolto norm e messo 1.0
	if (vec[i] == 0.) {
	  nzeroes+=1.;  //non giusti se interallacciato. Comunque non usati
	}
	sum+=abs(vec[i]);
	//if(i<=2)printf("i=%d vec=%f\n",i,vec[i]);
	//if(i>=len_call-2)printf("i=%d vec=%f\n",i,vec[i]);
	//if(i>=len_call/2-10 && i<=len_call/2+10)printf("i=%d vec=%f\n",i,vec[i]);
      }  //end del for
      //   printf("Fraction of zeroes: %f - Average vector value %f \n", nzeroes/len_call, sum/len_call);
      
      printf("---tim0 in mjd, VectNoise->GTime,VectNoise->ULeapS, %f %f %d \n", tim0,VectNoise->GTime,VectNoise->ULeapS);
      if(VectNoise->ULeapS==0){
	///PROVAAAAA       tolto messo sotto nleap=PSS_leapseconds(itt);
	utc=detector_param->itf_fixed_param->time+dot2utc;
	nleap=PSS_leapseconds(utc);
	VectNoise->ULeapS=nleap; //34; //ATT: valido dal 1 Jan 2009. Modicficato 17 luglio 2015 per essere generale. 
      }
      tim0= (*StartingTime-(VectNoise->ULeapS-19))/day2sec+44244; //mjd time
      printf("---AFTER tim0 in utc mjd, VectNoise->GTime,VectNoise->ULeapS nleap, %f %f %f %d %d \n", utc,tim0,VectNoise->GTime,VectNoise->ULeapS,nleap);

      ds_from_vec(d, vec,&tim0);
      d->tini1=tim0;
      d->tini2=tim0;
      
      /*
      //sposto qui incremento startingTime 
      //if(input_param->iopen!=0)
      //{
	  
	//   *StartingTime+=d->len*sampling/d->type;
      //}
      */
      
      
      // Write the proper beginning time of the trunck in detector_param->itf_var_param->mjdtime 
      
      if(d->type == 2){
	if (d->lcw%2 == 1){ 
	  //printf("**time from the chunk: d->lcw d->tini1 %ld %f\n",d->lcw,d->tini1);
	  detector_param->itf_var_param->mjdtime=d->tini1;
	}
	else{ 
	  //printf("**time from the chunk: d->lcw d->tini2= %ld %f\n",d->lcw,d->tini2);
	  detector_param->itf_var_param->mjdtime=d->tini2;
	}
      }
      else{
	//printf("**time from the chunk: d->lcw d->tini1 %ld %f\n",d->lcw,d->tini1);
	detector_param->itf_var_param->mjdtime=d->tini1;
      }

      //Prova Pia 2012 . Metto l' if 
      if (d->type < 2) {
	detector_param->itf_var_param->gpstime=VectNoise->GTime; /////PIA2012
      }
      else {
	detector_param->itf_var_param->gpstime=*StartingTime;  /////PIA 2012
	}

      if(verb>4)printf("Actual (var_param) gpstime= %f \n",detector_param->itf_var_param->gpstime);
           
    } // VectNoise!=NULL
  

    else { 
    
    	printf("Frame Noise not Found-->fill vec with zeros!\n");
        ds_from_vec(d, vec,&tim0);
        //nleap=PSS_leapseconds(itt);
        //VectNoise->ULeapS=nleap; //36 dal giugno 2015; //ATT: valido dal 1 Jan 2009. Modicficato 17 luglio 2015 per essere generale.
        //tim0= (*StartingTime-(VectNoise->ULeapS-19))/day2sec+44244; //mjd time
        //era cosi, 19 febbraio 2016  !!
	utc=*StartingTime/day2sec+dot2gps+dot2utc;
	nleap=PSS_leapseconds(utc);
	tim0= (*StartingTime-(nleap-19))/day2sec+44244; //mjd time
	printf(" nleap utc tim0=%d %f %f\n",nleap,utc,tim0);
	d->tini1=tim0;
        d->tini2=tim0;


	//////////%%%%%%%% 2012

   // Write the proper beginning time of the trunck in detector_param->itf_var_param->mjdtime 
      
     if(d->type == 2){
	if (d->lcw%2 == 1){ 
	  //printf("**time from the chunk: d->lcw d->tini1 %ld %f\n",d->lcw,d->tini1);
	  detector_param->itf_var_param->mjdtime=d->tini1;
	}
	else{ 
	  //printf("**time from the chunk: d->lcw d->tini2= %ld %f\n",d->lcw,d->tini2);
	  detector_param->itf_var_param->mjdtime=d->tini2;
	}
      }
      else{
	//printf("**time from the chunk: d->lcw d->tini1 %ld %f\n",d->lcw,d->tini1);
	detector_param->itf_var_param->mjdtime=d->tini1;
      }


      //Prova Pia 2012 . Metto l' if 
      if (d->type < 2) {
	detector_param->itf_var_param->gpstime=*StartingTime; /////PIA2012 --- Albe2014: was -10, corrected to avoid wrong date in SFDB09 files
      }
      else {
	detector_param->itf_var_param->gpstime=*StartingTime;  /////PIA 2012 --- Albe2014: was -10, corrected to avoid wrong date in SFDB09 files
	}


      if(verb>4)printf("Actual (var_param) gpstime= %f \n",detector_param->itf_var_param->gpstime);
           

	//////////%%%%%%%%%% 2012

    } // end frame noise not found
   

    if (VectNoise != NULL) {
      FrVectFree(VectNoise);
      VectNoise=NULL;   
    }
   
 
    
    free(vec);
    vec=NULL;

    
    if(input_param->iopen!=0)
    {
	  
	 *StartingTime+=d->len*sampling/d->type;
    }
    

    //printf("StartingTime,d->len*sampling/d->type,gps_endDB  %f %f %d\n",*StartingTime,d->len*sampling/d->type,*gps_endDB);
    //if((*StartingTime+d->len*sampling/d->type) > (float) *gps_endDB){
    




    //%%%%%
      
      if(d->type==2){
	if( *StartingTime+len_call*sampling >= (float) *gps_endDB){
	  puts("End of the run and SFDB file !!");
	  input_param->iopen=-1;
	}
      }
      //%%%%%%%%%

    if( *StartingTime > (float) *gps_endDB){
	//non e' gia' stato incrementato ? toglierei il +dlen ecc. 
	puts("End of the run and SFDB file !!");
	input_param->iopen=-1;
    }
  
    
      
    if (VectFlag != NULL) {
      FrVectFree(VectFlag);
      VectFlag=NULL;
    }
    

    free(vecF);
    vecF=NULL;
    //fclose(outFlag);
    //fclose(debugfile);

    return itt;
   
}

