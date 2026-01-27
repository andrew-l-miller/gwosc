
/*******Program to read and itf files and put data into an output FFT file *********/
/***********by Pia---Last version March 2009**********/
/*** To test woth valgrind:  valgrind --log-file=minnie ./crea_sfdb.out  ***/
/***To compile:  ./compila ***/

#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>

/****** PSS libraries****************/
#include "../pss_lib/pss_math.h"
#include "../pss_lib/pss_serv.h"
#include "../pss_lib/pss_snag.h"
#include "../pss_lib/pss_sfc.h"
/****** Antenna libraries**************/
#include "pss_ante.h"

/****** SFDB libraries**************/
#include "pss_sfdb.h"


/****** Pia's SDS libraries***************/
#include "pss_sds.h"

int main(void)
{
  double itt;
  int ifatti;
  INPUT_PARAM *input_param;  /*input parameters*/
  DETECTOR_PARAM *detector_param;  /*detector (bar or interferometer) parameters*/
  HEADER_PARAM *header_param; /*Structure with parameters to be written in the sfdb file*/
  EVEN_PARAM *even_param; /*Structure with the event parameters*/
  EVF_PARAM *evf_param; /*Structure with the EVF- peaks in frequency domain- parameters*/
  DS *d;                   /*DS with data from files*/
  DS *d_h;                 /*DS with h reconstructed data*/
  GD *gd;                   /* Virtual GD to be used for data in frequency domain*/
  GD *gd_highpass;         /*GD with highpassed data*/
  GD *gd_clean;            /*GD with clean data, without large events*/
  GD *gd_freq;            /*GD with freq. FFT data, alternate: Real, Imag*/
  int ii;
  int errorcode;

  ALLPERS_ *allpers;  //ALLPERS structure with allowed periods in SDS files
  HOLES_ *holes;  //HOLES structure with holes in SDS files
  SFC_ *sfc_data;  //SFC structure to read parameters and data in sds files. 

  GD *gd_appo;  /*GD with clean data, needed in the function which cleans the subbands*/
  GD *gd_band ;  /*GD needed to  clean in the subbands  Jan 2009*/
  int defaults=1; //0=reads defaults values. Else: asks for input values
  float *ps;
  int dim;

  int casin,casout;

  casin=1; // Change for your use !
           // 1 -> input from sds
           // 2 -> input from frames

  casout=1; // Change for your use !
           // 1 -> output on sfdb09 files
           // 2 -> ...

  
  /*Crea INPUT_PARAM*/
  LOG_INFO=logfile_open("crea_sfdb");
  input_param=crea_input(defaults,casout);
  dim=(input_param->len/input_param->red); 
  ps=(float*)malloc(sizeof(float)*2*dim);  //prova
  
  /*Crea detector parameters */
  detector_param=crea_parameters(casin);
  /*Crea header (of the output file) and even, evf parameters */
  header_param=crea_sfdbheader();
  even_param=crea_evenparam(input_param->len);
  puts("You might want to change even_param->absvalue");
  printf("Now it is %d\n",even_param->absvalue);
  puts("which means it uses");
  if(even_param->absvalue==1)puts("Abs values to look for events, and thrs=6");
  if(even_param->absvalue!=1)puts("Values with sign to look for events, and thrs=5");
  puts("Change, if needed, its value in the main code. 1=uses abs values. Else: uses values with sign");
  evf_param=crea_evfparam();

  /*Crea DS with the data*/
  d=crea_ds(input_param->len,input_param->typ,input_param->file_extens);
 /*Crea DS with h reconstructed data: needed if raw are not in h reconstructed*/
  d_h=crea_ds(input_param->len,input_param->typ_h,input_param->capt_h);
  gd=crea_gd(input_param->len,0.,0.,input_param->capt_gd);
  gd_highpass=crea_gd(gd->n,gd->ini,gd->dx,input_param->capt_gd);
  gd_clean=crea_gd(gd->n,gd->ini,gd->dx,input_param->capt_gd);
  gd_appo=crea_gd(gd->n,gd->ini,gd->dx,input_param->capt_gd);
  gd_band=crea_gd(gd->n,gd->ini,gd->dx,input_param->capt_gd);
  gd_freq=crea_gd(gd->n,gd->ini,gd->dx,input_param->capt_gd);
  
 

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
			  // ...
			  break;
		  case 3:
			  // ...
			  break;
	  }
	 
     
     //      ii=gd2gdfreq(gd_clean,gd_freq,gd,gd_highpass,gd_band,gd_appo,header_param,even_param,input_param);
   

          switch (casout)
	  {
	    case 1:
	      puts("Chiamo gd2sfdbfile");
	      errorcode=gd2sfdbfile(gd,header_param,evf_param,ps,even_param,input_param,sfc_data);

	       
		    break;
	    case 2:
	            //...
	            break;
	    case 3:
	            //...
	            break;
	  }
      	
	     
      if(input_param->iopen==-1)break;
      } //end for ifatti
  
 if(casin==1){
  free(allpers);
  free(holes);
  free(sfc_data);
 }
  logfile_stop(LOG_INFO);
  logfile_close(LOG_INFO);
  puts("End of the data processing on the whole run");
   
  return 0;
}
