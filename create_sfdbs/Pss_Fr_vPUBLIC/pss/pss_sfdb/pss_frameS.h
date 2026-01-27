/**pss_sds.h: pss_sds software: to read data from SDS files and put into a DS the h-reconstructed data**/
/**Last version: August, 11, 2005**/
/**Author: Pia**/

FILE *insds; //maybe not needed


#include <values.h>
#include <stdarg.h>
#include  "../../framelib/FrameL.h"
#include "../../framelib/FrVect.h"
//#include "Frv.h" 
//#include "Fd.h"

#define utc2mjd  2400000.5E0


/*************************To read the channels and put into a DS****************/

double frame2dsS(DS *d, INPUT_PARAM *input_param,int verb,DETECTOR_PARAM *detector_param,FrFile *iFileNoise,double *StartingTime,double *gps_startDB,double *gps_endDB,HEADER_PARAM *header_param, SCIENCE_SEGMENT* science_segment);


//double frame2hds(DS *d_h,INPUT_PARAM *input_param, DETECTOR_PARAM *detector_param,EVEN_PARAM *even_param,EVF_PARAM *evf_param,HEADER_PARAM *header_param, SCIENCE_SEGMENT* science_segment, FrFile *iFileNoise,double *StartingTime,int *gps_startDB,int *gps_endDB);

double frame2hds(DS *d_h,INPUT_PARAM *input_param, DETECTOR_PARAM *detector_param,EVEN_PARAM *even_param,EVF_PARAM *evf_param,HEADER_PARAM *header_param, SCIENCE_SEGMENT* science_segment, FrFile *iFileNoise,double *StartingTime,double *gps_startDB,double *gps_endDB);
