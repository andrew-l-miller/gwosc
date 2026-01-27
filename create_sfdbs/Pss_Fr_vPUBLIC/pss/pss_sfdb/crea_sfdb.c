
/*******Program to read and itf files and put data into an output FFT file *********/
/***********by Pia---Last version March 2009**********/
/*** To test woth valgrind:  valgrind --log-file=minnie ./crea_sfdb.out  ***/
/***To compile:  ./compila ***/
// valgrind  --leak-check=full --log-file=minnie ./crea_sfdb.out <inputcorto

#include<stdio.h>

#include<math.h>
#include<string.h>
#include<stdlib.h>
#include <values.h>
#include <malloc.h>
#include "../../framelib/FrameL.h"
#include "../../framelib/FrVect.h"
/****** PSS libraries****************/
#include "../pss_lib/pss_math.h"
#include "../pss_lib/pss_serv.h"
#include "../pss_lib/pss_snag.h"
#include "../pss_lib/pss_sfc.h"
/****** Antenna libraries**************/
#include "pss_ante.h"

/****** SFDB libraries**************/
#include "pss_sfdb.h"


/****** Pia's SDS and frame libraries***************/
#include "pss_sds.h"
#include "pss_frameS.h"



int main(void)
{

  double itt;
  int ifatti;
  int i;
  INPUT_PARAM *input_param;  /*input parameters*/
  DETECTOR_PARAM *detector_param;  /*detector (bar or interferometer) parameters*/
  HEADER_PARAM *header_param; /*Structure with parameters to be written in the sfdb file*/
  SCIENCE_SEGMENT *science_segment=NULL; /*Structure with information about science segments, read from text file*/
  EVEN_PARAM *even_param; /*Structure with the event parameters*/
  EVEN_PARAM *even_paramb; /*Structure with the BAND event parameters*/
  EVF_PARAM *evf_param; /*Structure with the EVF- peaks in frequency domain- parameters*/
  DS *d;                   /*DS with data from files*/
  DS *d_h;                 /*DS with h reconstructed data*/
  GD *gd;                   /* Virtual GD to be used for data in frequency domain*/
  GD *gd_highpass;         /*GD with highpassed data*/
  GD *gd_clean;            /*GD with clean data, without large events*/
  GD *gd_freq;            /*GD with freq. FFT data, alternate: Real, Imag*/
  GD *gd_res;             /*if resampling on the data*/
  int ii;
  int errorcode;

  ALLPERS_ *allpers;  //ALLPERS structure with allowed periods in SDS files
  HOLES_ *holes;  //HOLES structure with holes in SDS files
  SFC_ *sfc_data;  //SFC structure to read parameters and data in sds files. 
  

  struct FrFile *iFileNoise;
  


  double gps_startDB,gps_endDB;  //for frames. ATT:era int modificato LUG 2018
  double StartingTime;  //for frames

  GD *gd_appo;  /*GD with clean data, needed in the function which cleans the subbands*/
  GD *gd_band ;  /*GD needed to  clean in the subbands  Jan 2009*/
  int defaults=1; //0=reads defaults values. Else: asks for input values
  float *ps;
  float *psTOT; //needed for the peakmaps
  int dim;
  int res;
  int casin,casout;
  int nfftfile;

  int iGRID;  //if 1 we are running under GRID. Filenames in the ffl file begin with lfn:
  FILE *FFLLISTA;
  FILE *TESTdump;
 
  //char filelista[MAXMAXLINE+1];
  
  casin=2; // Change for your use !
           // 1 -> input from sds
           // 2 -> input from frames

  casout=1; // will be changed in crea_input_param  for your use !
           // 1 -> output on sfdb09 files
           // 2 -> output on "void" sfdb09 files, containing header+short spectrum+ only 1 datum of FFts..
          // 3 -> output on sfdb09 files and peakmap .p10 files

  
  /*Crea INPUT_PARAM*/
  LOG_INFO=logfile_open("crea_sfdb");
  input_param=crea_input(defaults,casin,&casout);
  res=input_param->resampling;
  dim=(input_param->len/(res*input_param->red)); 
  ps=(float*)malloc(sizeof(float)*2*dim);
  psTOT=(float*)malloc(sizeof(float)*2*input_param->len);
 
  /*Crea detector parameters */
  detector_param=crea_parameters(casin);
  /*Crea header (of the output file) and even, evf parameters */
  header_param=crea_sfdbheader(input_param);

  /*Crea science segment structure  */
  science_segment=crea_science_segment(input_param);
  
  if (science_segment != NULL) {
  	printf("Using science segment file -> Setting header_param->sciseg=1\n");
  	header_param->sciseg=1;
  } else {
  	printf("Not using science segment file -> Setting header_param->sciseg=0\n");
  }
  
  even_param=crea_evenparam(input_param->len/res,1);
  even_paramb=crea_evenparam(input_param->len/res,0);
  puts("You might want to change even_param->absvalue");
  printf("Now it is %d, which means it uses\n",even_param->absvalue);
  if(even_param->absvalue==1)puts("Abs values to look for events, and thrs=6");
  if(even_param->absvalue!=1)puts("Values with sign to look for events, and thrs=5");
  puts("Change, if needed, its value in the main code. 1=uses abs values. Else: uses values with sign");
  evf_param=crea_evfparam(input_param);
  /*Crea DS with the data*/
  d=crea_ds(input_param->len,input_param->typ,input_param->file_extens);
 /*Crea DS with h reconstructed data: needed if raw are not in h reconstructed*/
  d_h=crea_ds(input_param->len,input_param->typ_h,input_param->capt_h);
  gd=crea_gd(input_param->len,0.,0.,input_param->capt_gd);
  gd_highpass=crea_gd(gd->n/res,gd->ini,res*gd->dx,input_param->capt_gd);
  gd_clean=crea_gd(gd->n/res,gd->ini,res*gd->dx,input_param->capt_gd);
  gd_appo=crea_gd(gd->n/res,gd->ini,res*gd->dx,input_param->capt_gd);
  gd_band=crea_gd(gd->n/res,gd->ini,res*gd->dx,input_param->capt_gd);
  gd_freq=crea_gd(gd->n/res,gd->ini,res*gd->dx,input_param->capt_gd);
  gd_res=crea_gd(gd->n/res,gd->ini,res*gd->dx,input_param->capt_gd);
 

 for (ifatti=0;ifatti<input_param->itotal;ifatti++) /*itotal= max number of  chunks */
   {

     switch (casin)
	  {
		  case 1:
		    if(ifatti==0){
		       /*SDS: allpers allocation and creation*/
		      allpers=(ALLPERS_ *)malloc(sizeof(ALLPERS_));
		      allpers->nper=0;
		      holes=(HOLES_ *)malloc(sizeof(HOLES_));
		      sfc_data=sds_open(input_param->filename,allpers);//alloc sfc_data
		    }
		    itt=sds2hds(d_h,input_param,detector_param,even_param,evf_param,allpers,header_param,holes,sfc_data);
		     ii= gd_ds(gd,d_h);    //from a chunk of a DS to a GD   

			  break;
		  case 2:
			  // chiamata ai frame
		    if(ifatti==0){
		      FFLLISTA=fopen(input_param->filename,"r");
		      //TESTdump=fopen("testdump.txt","w");
		      puts("Read the first frame file name. To check if we are running under grid (begins with /storage/gpfs_virgo4/)");
		      printf("input param filename:  %s\n",input_param->filename);
		      fscanf(FFLLISTA,"%s",input_param->filelista);
		      printf("Name of the first frame file in the list:  %s\n",input_param->filelista);
		      fclose(FFLLISTA);
		      iGRID=0; //not using GRID as default
		      if (strncmp(input_param->filelista,"/storage/gpfs_virgo4/",21) == 0) iGRID=1;		     
		      printf("iGRID %d \n",iGRID);
		      printf("READ file from input param filename:  %s\n",input_param->filename);
		      iFileNoise = FrFileINew(input_param->filename);
		      //PROVA
		     
		      //FrFileIClose(iFileNoise);
		      //FrFileIOpen(iFileNoise);  //bohhh vediamo...luglio 2018
		      //FrSetIni(iFileNoise);   
		      //DUMP
		      //FrameH *frame= NULL;
		      //double ftime;
		      //printf("PRIMA del while\n");
		      //while((frame = FrameReadRecycle(iFileNoise,frame)) != NULL){
		      //ftime = frame->GTimeS + 1.e-9 * frame->GTimeN;
		      //printf("===>ftime=%f\n",ftime);
		      //FrameDump(frame, TESTdump, 2);}
		      //fclose(TESTdump);
		       //iFileNoise = FrFileNew(input_param->filename, 0, NULL, 0);
		      //FrFileH *fileH;
		      //fileH=iFileNoise->fileH;
		      printf(" ===>After FrFileINew \n");
		      //printf("====> fileName= %s\n",fileH->fileName);
		      //printf("====> Tstart da fileH = %.0f\n",fileH->tStart);
		      printf("Nel main. Approx. GPS_START = %17.6f\n",FrFileITStart(iFileNoise));
		      printf("Nel main. Approx. GPS_STOP = %17.6f\n",FrFileITEnd(iFileNoise));
		      //printf("==>int startime dal file %d\n",iFileNoise->startTime);
		      //printf("==>int end time dal file %d\n",iFileNoise->closingTime);
		      //end prova
		    }
		    printf("Before frame2hds");
		    itt=frame2hds(d_h,input_param,detector_param,even_param,evf_param,header_param,science_segment, iFileNoise,&StartingTime,&gps_startDB,&gps_endDB);
		    printf("After frame2hds");
		     ii= gd_ds(gd,d_h);    //from a chunk of a DS to a GD   
		     if(input_param->resampling > 1){
		       ii=resampling(gd_res,gd,header_param);
       
		     }
			  break;
	      
		  case 3:
			  // ...
			  break;
	  }
	 
     if(res<=1){
       ii=gd2gdfreq(gd_clean,gd_freq,gd,gd_highpass,gd_band,gd_appo,header_param,even_param,even_paramb,input_param,science_segment);
     }
     else {
       ii=gd2gdfreq(gd_clean,gd_freq,gd_res,gd_highpass,gd_band,gd_appo,header_param,even_param,even_paramb,input_param,science_segment);
	     
     }
     
     if (header_param->nfft>=1)
     {

          switch (casout)
	  {
	    case 1:
	      //puts("gd2sfdbfile called");
	      errorcode=gd2sfdbfile(gd_clean,header_param,evf_param,ps,even_param,input_param,psTOT,iGRID);
	      //	      errorcode=gd2sfdbfile(gd_res,header_param,evf_param,ps,even_param,input_param);
	      //errorcode=gd2sfdbfile(gd,header_param,evf_param,ps,even_param,input_param);
	      if(res>1) {
		header_param->nsamples*=res;
		header_param->tsamplu/=res;
	      }
		    break;
	    case 2:	            
	      //puts("gd2sfdbfilevoid called");
	      errorcode=gd2sfdbfilevoid(gd_clean,header_param,evf_param,ps,even_param,input_param,psTOT);
	      if(res>1) {
		header_param->nsamples*=res;
		header_param->tsamplu/=res;
	      }
	            break;
	    case 3:
	      //puts("gd2sfdbfile called. Does also the peakmaps");
	      errorcode=gd2sfdbfile(gd_clean,header_param,evf_param,ps,even_param,input_param,psTOT,iGRID);
	      //casout 3 also peakmaps are done
              if(input_param->iopen ==-1){
                  //re-write the number of fft done in the pm file.
	
		nfftfile=input_param->nfftONE;
                  fseek(P10,0,SEEK_SET);
                  fwrite((void*)&nfftfile, sizeof(int),1,P10);  //da correggere
                  printf("Re-write nfft in the peakmap file. It is: %d \n", nfftfile);
                  //end of rewrite
                }
	      if(res>1) {
		header_param->nsamples*=res;
		header_param->tsamplu/=res;
	            break;
	  }
	  } 	
    } else
    {
    	printf(" No FFT done!\n");
    
    }	
	     
     if(input_param->iopen==-1){
       
       if(input_param->flagonline !=1)fclose(SFDB);
       if(input_param->flagonline>=3)fclose(P10);
       if(iGRID>=1)ToCopyFiles(input_param,header_param);     //rimettere 1
       puts("******END ALL DONE !!");
       break;}
     
   } //end for ifatti
    
  
 if(casin==1){
  free(allpers);
  free(holes);
  free(sfc_data);
  
 }
 if(casin>=2){
  FrFileFree(iFileNoise);
 }
  free(evf_param);
  free(even_param->ener);
  free(even_param->xamed);
  free(even_param->xastd);
  free(even_param->begin);
  free(even_param->duration);
  free(even_param->imax);
  free(even_param->crmax);
  free(even_paramb->ener);
  free(even_paramb->xamed);
  free(even_paramb->xastd);
  free(even_paramb->begin);
  free(even_paramb->duration);
  free(even_paramb->imax);
  free(even_paramb->crmax);
  free(even_param);
  free(even_paramb);
  free(input_param);
  free(detector_param->itf_fixed_param);
  free(detector_param->itf_var_param);
  free(detector_param);
  free(header_param);
  if (science_segment != NULL) {
  	free(science_segment);
  }
 
  i=logfile_stop(LOG_INFO);
  i=logfile_close(LOG_INFO);
  free(ps);
  free (psTOT);
  free(d->y1);
  free(d->y2);
  free(d_h->y1);
  free(d_h->y2);
  free(d_h);
  free(d);
  free(gd->y);
  free(gd_highpass->y);
  free(gd_clean->y);
  free(gd_freq->y);
  free(gd_appo->y);
  free(gd_band->y);
  free(gd);    
  free(gd_highpass); 
  free(gd_clean);    
  free(gd_freq);
  free(gd_appo);
  free(gd_band);
  puts("End of the data processing on the whole run");
   
  return 0;
}
