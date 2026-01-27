/**pss_sds.h: pss_sds software: to read data from SDS files and put into a DS the h-reconstructed data**/
/**Last version: August, 11, 2005**/
/**Author: Pia**/

FILE *insds; //maybe not needed

#define utc2mjd  2400000.5E0
/***************************To fill the parameter structures*********/
void def_fixed_param_sds(short int k_par,DETECTOR_PARAM *detector_param,double valore);

/*************************To read the channels and put into a DS****************/

double sds2ds(DS *d, INPUT_PARAM *input_param,int verb,DETECTOR_PARAM *detector_param,ALLPERS_ *allpers,HOLES_ *holes,SFC_ *sfc_data);


double sds2hds(DS *d_h,INPUT_PARAM *input_param, DETECTOR_PARAM *detector_param,EVEN_PARAM *even_param,EVF_PARAM *evf_param,ALLPERS_ *allpers,HEADER_PARAM *header_param,HOLES_ *holes,SFC_ *sfc_data);



/*****************To manipulate vectors and ds ************/
void ds_from_vec(DS *d, float *vec,double *tim0);
/*
Puts the data in the vector vec in a DS, overlapping if needed (d->type=2) 
The DS has to be already created, with crea_ds. tim0= Beg. time of the chunk which is going to be written
 */
