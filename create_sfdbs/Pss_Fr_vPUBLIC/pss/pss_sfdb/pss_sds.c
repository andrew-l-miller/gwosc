/**pss_sds.c: pss_sds software: to read data from SDS files and put into a DS the h-reconstructed data**/
/**Last version: August, 11, 2005**/
/**First version: August 2005**/
/**Author: Pia**/

#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>

/*******PSS astro*******/
#include"../pss_astro/jplbin.h"
#include"../pss_astro/novas/novas.h"
#include "../pss_astro/daspostare.h"
#include"../pss_astro/pss_astro.h"
/****** PSS libraries***************/
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



double sds2hds(DS *d_h,INPUT_PARAM *input_param, DETECTOR_PARAM *detector_param,EVEN_PARAM *even_param,EVF_PARAM *evf_param,ALLPERS_ *allpers,HEADER_PARAM *header_param,HOLES_ *holes,SFC_ *sfc_data)

{
  double itt;
  int verb;
  float sum2;
  int i;

  header_param->gps_sec=0.; //inizialization
  verb=input_param->verb;
  header_param->wink=input_param->wink;
  
  if(d_h->lcr==0){
    //sfc_data=sds_open(input_param->filename,allpers); //allocate sfc_data 
    strcpy(header_param->sfdbname,input_param->sfdbname);
    printf("Sfdb filename, as written in the header param= %s \n", header_param->sfdbname);
  }
  if(d_h->lcr==0)detector_param->itf_var_param->nfft=0;
  if(d_h->lcr==0)detector_param->itf_var_param->kdati_tot=0;
  detector_param->itf_var_param->n_flag=0; //number in every fft
  itt=sds2ds(d_h,input_param,verb,detector_param,allpers,holes,sfc_data);
  d_h->lcr++;  //another chunk has been served


  if(verb>0)printf(" fft number = %d \n",detector_param->itf_var_param->nfft); 
  /* Definition of the header parameters */
  header_param->nfft=detector_param->itf_var_param->nfft; 
  header_param->n_flag=detector_param->itf_var_param->n_flag; //how many ??? in that FFT have a flag !=0 
  header_param->endian=1.0;
  header_param->detector=1; /*if 1 sds format (itf); if 2 frame format (itf)*/ 
  header_param->gps_sec=(int) detector_param->itf_var_param->gpstime;
  header_param->gps_nsec=(int)((detector_param->itf_var_param->gpstime-(int)detector_param->itf_var_param->gpstime)*1.0e9);
  header_param->tbase=detector_param->itf_var_param->fftlen*detector_param->itf_var_param->tsamplu;
  header_param->firstfreqindex=(int) (detector_param->itf_var_param->frinit/detector_param->itf_var_param->deltanu);
  header_param->nsamples=(int) detector_param->itf_var_param->fftlen/2;
  header_param->mjdtime=detector_param->itf_var_param->mjdtime; //ATT: it can be non-continuous. There are holes !!!
  //header_param->einstein=EINSTEIN;  //scaling factor
  header_param->frinit=detector_param->itf_var_param->frinit;
  header_param->tsamplu=detector_param->itf_var_param->tsamplu;
  header_param->deltanu=detector_param->itf_var_param->deltanu;
  
  header_param->red=input_param->red;        //reduction factor for the very short FFTs
  header_param->typ=input_param->typ;        //0, 1 = non overlapping data  2=overlapping data

  header_param->normd=sqrt(header_param->tsamplu/(2*header_param->nsamples)); //piapia 24 oct 2005
  if(header_param->wink !=5)header_param->normw=1; //to be added the factor for Hamming,Blackman
  if(header_param->wink==1)header_param->normw=sqrt(3.0); //Hanning
  if(header_param->nfft==1 && header_param->wink==5){     //flat + cos edges
    sum2=0;
    for (i=0; i <header_param->nsamples/2 ; i++){  //first 1/4 window
     sum2+=pow((0.5-0.5*cos(i*PIG/(header_param->nsamples/2))),2);
    }
    for (i=header_param->nsamples/2; i<2*header_param->nsamples-header_param->nsamples/2 ; i++){
     sum2+=1.;
    }
    for (i=0; i <header_param->nsamples/2 ; i++){  //last 1/4 window = first 1/4 window
     sum2+=pow((0.5-0.5*cos(i*PIG/(header_param->nsamples/2))),2);
    }	   
    header_param->normw=sqrt(2*header_param->nsamples/sum2);
	
  }
  printf("header_param->normw = %f \n",header_param->normw); 
  //if(verb>1)puts("Exit from sds2hds:time domain data"); 
  if(header_param->nfft==1){
     /* Definition of parameters for the Log file */
     logfile_comment(LOG_INFO,"par GEN: general parameters of the run");
     logfile_par(LOG_INFO,"GEN_BEG",header_param->mjdtime,3);
     logfile_par(LOG_INFO,"GEN_NSAM",header_param->nsamples,1);
     logfile_par(LOG_INFO,"GEN_DELTANU",header_param->deltanu,0);
     logfile_par(LOG_INFO,"GEN_FRINIT",header_param->frinit,0);
     logfile_comment(LOG_INFO,"   GEN_BEG is the beginning time (mjd)");
     logfile_comment(LOG_INFO,"   GEN_NSAM the number of samples in 1/2 FFT");
     logfile_comment(LOG_INFO,"   GEN_DELTANU the frequency resolution");
     logfile_comment(LOG_INFO,"   GEN_FRINIT the beginning frequency of the FFT");
     //logfile_par(LOG_INFO,"EVT_CR",even_param->cr,0);
     //logfile_par(LOG_INFO,"EVT_TAU",even_param->tau,0);
     //logfile_par(LOG_INFO,"EVT_DEADT",even_param->deadtime,0);
     //logfile_par(LOG_INFO,"EVT_EDGE",even_param->edge,0);
     logfile_comment(LOG_INFO,"   EVT_CR is the threshold");
     logfile_comment(LOG_INFO,"   EVT_TAU the memory time of the AR estimation");
     logfile_comment(LOG_INFO,"   EVT_DEADT the dead time [s]");  
     logfile_comment(LOG_INFO,"   EVT_EDGE seconds purged around the event");
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
  printf("evf param tau_Hz %f\n",evf_param->tau_Hz);
  printf("even param tau %f\n",even_param->tau);
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

double sds2ds(DS *d,INPUT_PARAM *input_param,int verb, DETECTOR_PARAM *detector_param,ALLPERS_ *allpers,HOLES_ *holes,SFC_ *sfc_data)
{
  double itt;
  long len;
 
  char filesds[160];
  float *buffer; 
  float *vec; 
  double tim0; 
  int chn;
  int dim_buffer;
  double diff_tim0,exp_diff;
  float norm;  //such that the FFT of time data gives the noise spectral amplitude
               //norm=sqrt(dt/fftlen)
  int i;
  long len_call;
  double mjdtime;

   detector_param->itf_var_param->mjdtime=0.; //4/02/08
   detector_param->itf_var_param->gpstime=0.; //4/02/08

   if(input_param->iopen==0){
    strcpy(filesds,input_param->filename);
    puts("A new run (hence file in the SFDB) has started\n");
    //sfc_data=sds_open(filesds,allpers);
    //sfc_data=sfc_open(filesds);
    printf("SDS file parameters= t0 (mjd) dt (s) %f %f\n",sfc_data->t0,sfc_data->dt);
    printf("SDS file parameters= name point0 hlen %s %ld %ld\n",sfc_data->filme,sfc_data->point0,sfc_data->hlen);
    printf("SDS file parameters= file %s \n",sfc_data->file);
    printf("SDS file parameters= nch %ld \n",sfc_data->nch);
    input_param->iopen=1;
    itt=sfc_data->t0-dot2mjd-daybegin; //days from 1 Jan 2001
    if(d->type<=1)exp_diff=sfc_data->dt*d->len/day2sec;
    if(d->type==2)exp_diff=sfc_data->dt*(d->len/d->type)/day2sec;
    d->cont=0;
    d->dt=sfc_data->dt;
    d->tini1=sfc_data->t0-exp_diff;
    d->tini2=sfc_data->t0-exp_diff;
    detector_param->itf_fixed_param->type=d->type;
    printf("initial time of the FIRST file  (days from 1 Jan 2001) %f \n",itt);
    detector_param->itf_fixed_param->time=itt+daybegin; /*Beginning time of the run in days from 1900*/
        sfc_show(sfc_data);
   }
  
   len=d->len;
   d->cont=0; //control variable.
   
  
      mjdtime=detector_param->itf_fixed_param->time+dot2mjd;
      detector_param->itf_fixed_param->gpstime=mjd2gps(mjdtime);
   /* Var parameters */ 
  detector_param->itf_var_param->nfft++;
  detector_param->itf_var_param->tsamplu=sfc_data->dt;
  detector_param->itf_var_param->fftlen=d->len;
  detector_param->itf_var_param->frinit=0.;
  detector_param->itf_var_param->deltanu=1/(d->len*detector_param->itf_var_param->tsamplu);
  if(verb >1)printf("deltanu,d->len %f %ld \n",detector_param->itf_var_param->deltanu,d->len);
  detector_param->itf_var_param->frinit=0.;
  detector_param->itf_var_param->freq_min=0;  //piapia: to be defined as input parameters
  detector_param->itf_var_param->freq_max=detector_param->itf_var_param->freq_min+1/(2*sfc_data->dt);

  if(verb>0)printf("fixed time,gpstime= %f %f\n",detector_param->itf_fixed_param->time,detector_param->itf_fixed_param->gpstime);
  if(d->type<=1)exp_diff=sfc_data->dt*len/day2sec;
  if(d->type==2)exp_diff=sfc_data->dt*(len/d->type)/day2sec;
  printf("exp_diff from subsequent ffts, in seconds %f \n", exp_diff*day2sec);
  /*Now read the data in the file:*/
  /*vec_from_sds gives the data and then ds_from_vec a ds with the data*/
  /* The "real" time should be read again here, because vec contains zero padding, to remove "short" holes in the data */
  /*IF the holes are "long" there is no zero padding and we must close the SFDB file, hence put iopen=-1, and start with a new one*/

  vec=(float *)malloc((size_t) (d->len)*sizeof(float));
  dim_buffer=(int) len*sfc_data->nch;
  //printf("len sfc_data->nch dim_buffer= %ld %ld %d\n",len, sfc_data->nch,dim_buffer);
  buffer=(float *)malloc((size_t) (dim_buffer)*sizeof(float));
  chn=1; //piapia: to be defined as input parameter
  norm=sqrt(sfc_data->dt/d->len); //such that the FFT of time data gives the noise spectral amplitude. pia:CHECK
  len_call=len;
  if(d->type==2&&detector_param->itf_var_param->nfft>1)len_call=len/2;
  printf("fft numb.(the next to be done) len_call norm= %d %ld %f\n",detector_param->itf_var_param->nfft,len_call,norm);
  vec_from_sds1(buffer,vec, &tim0, holes, sfc_data, chn, len_call, allpers); //pia version of vec_from_sds
  printf(" *** Actual sds file name file %s \n",sfc_data->filme);
  //printf("**time from vec_from_sds: tim0= %f\n",tim0);
  if(d->type==2&&d->nc1!=0)printf("**beg. time of the chunk: tim0-exp_diff= %f\n",tim0-exp_diff);

  if(d->type == 2){
    if (d->lcw%2 == 1){ 
      //printf("**time from the previous chunk: d->lcw d->tini1 %ld %f\n",d->lcw,d->tini1);
       diff_tim0=tim0-exp_diff-d->tini1; //if overlaps the beg. time is tim0-exp_diff
    }
    else{ 
      //printf("**time from the previous chunk: d->lcw d->tini2= %ld %f\n",d->lcw,d->tini2);
       diff_tim0=tim0-exp_diff-d->tini2; //if overlaps the beg. time is tim0-exp_diff
    }
  }
  else{
    //printf("**time from the previous chunk: d->lcw d->tini1 %ld %f\n",d->lcw,d->tini1);
     diff_tim0=tim0-d->tini1;          //if not overlap the beg. time is tim0
  }
   if(d->nc1==0)diff_tim0=tim0-d->tini1; //in the first chunk, even if overlapped
   if((fabs) ((diff_tim0-exp_diff)*day2sec) >  sfc_data->dt){
       printf("ATT !! diff_tim0 exp_diff  %e %e \n",diff_tim0, exp_diff);
       if(d->type==2)puts("The overlapping should restart. (d->cont put =1) !!");
       d->cont=1; //use d->cont (control variable) in the situation of jumps in the data
   }
    
   //puts("dopo vec_from_sds: print vec (before EINSTEIN normalization) for test");
  //typeVect(vec,512,522,5,1);
  //typeVect(vec,len_call-11,len_call-1,5,1);
  free(buffer);
  for(i=0;i<len_call;i++)vec[i]*=(1.0/EINSTEIN);//output data in the ds will be in EINSTEIN units
 
                                              //piapia tolto norm e messo 1.0
 ds_from_vec(d, vec,&tim0); 
  free(vec);
  /*Write the proper beginning time of the trunck in detector_param->itf_var_param->mjdtime */
  
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

   mjdtime=detector_param->itf_var_param->mjdtime;
   detector_param->itf_var_param->gpstime=mjd2gps(mjdtime);
  
  if(verb>0)printf("Actual (var_param) gpstime= %f \n",detector_param->itf_var_param->gpstime);
   if(verb>0)printf("Actual (var_param) mjdtime= %f \n",detector_param->itf_var_param->mjdtime);
  if(sfc_data->eof==2 || sfc_data->eof==3 || sfc_data->fid==NULL){ 
    puts("End of the run and SFDB file !!");
    input_param->iopen=-1;
  }
  return itt;
}

/*****************to manipulate vectors and ds ************/
void ds_from_vec(DS *d, float *vec, double *tim0)
/*
Puts the data in the vector vec in a DS, overlapping if needed (d->type=2) 
The DS has to be already created, with crea_ds. tim0= Beg. time of the chunk which is going to be written
 */
{
  int i;
  long nc1,nc2;
  long len2;
 
  len2=d->len/2;
  nc1=d->nc1;
  nc2=d->nc2;
  if (d->lcw == 0)
    {
      d->cont=0;
    }

  if(d->type==0){ //y1=actual chunk. y2 not used.
    for(i=0;i<d->len;i++)d->y1[i]=vec[i];
    d->nc1++;
    d->lcw++;
    d->tini1=*tim0;
  } 

  if(d->type==1){
    //y1= actual chunk--y2=always the previous chunk. Not overlapping.
    if(nc1!=0)
    {
      if(d->cont==0)for(i=0;i<d->len;i++)d->y2[i]=d->y1[i];
      if(d->cont==1)for(i=0;i<d->len;i++)d->y2[i]=0.;
    }
    d->nc2++;
    d->tini2= d->tini1;
    for(i=0;i<d->len;i++)d->y1[i]=vec[i];
    d->nc1++;
    d->lcw++;
    d->tini1=*tim0;
  }
  if(d->type==2){
    if(nc1==0){ //first chunk
       for(i=0;i<d->len;i++)d->y1[i]=vec[i];
       d->nc1++;
       d->lcw++;
       d->tini1=*tim0;
       d->cont=0;
       return;
    } 
    if(nc1<=nc2){
      if(nc1>0){ //odd chunks, but the first
	if(d->cont==0)for(i=0;i<len2;i++)d->y1[i]=d->y2[i+len2];
	if(d->cont==1)for(i=0;i<len2;i++)d->y1[i]=0.;
	for(i=0;i<len2;i++)d->y1[i+len2]=vec[i];
	d->nc1++;
	d->lcw++;
	d->tini1=*tim0-(len2*d->dt/day2sec);
      }
    }
    else  /* even chunks :  produces d.y2 */
    {
      if(d->cont==0)for(i=0;i<len2;i++)d->y2[i]=d->y1[i+len2];
      if(d->cont==1)for(i=0;i<len2;i++)d->y2[i]=0;
      for(i=0;i<len2;i++)d->y2[i+len2]=vec[i];
      d->nc2++;
      d->lcw++;
      d->tini2=*tim0-(len2*d->dt/day2sec);
    }
 
  }
  d->cont=0;
  return ;
 }



