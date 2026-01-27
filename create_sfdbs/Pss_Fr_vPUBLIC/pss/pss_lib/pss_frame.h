/*___________________________________________________
 ¦                                                   ¦
 ¦                   pss_frame.h                     ¦
 ¦       by Sergio Frasca - ((( 0 ))) Virgo          ¦
 ¦                    March 2000                     ¦
 ¦___________________________________________________¦*/

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

//struct r_struct{
//};


struct Frame_ch{
	long frameN;
	long nch;
	char name[1000];
	char comment[2000];
	double samp[50];
	long len[50];
	double tin;
};


int read_fmnl_files1(int nfiles, char** files);


struct FrFile* open_fr(char* fileName);


void inquire_fr_file(char* fileName);


struct Frame_ch inquire_ch(char* fileName);


struct FrAdcData* inquire_adc(char* fileName, char *adcName, int verb);


int get_1ch(struct FrFile *iFile, char *adcName, GD *Dat, int verb);


int get_1chA(struct FrFile *iFile, char *adcName,
			float **data1, long* nData);


int get_1ch_vect(char *fileName, char *chName,
			float *dat, long kframe, long ndat, long lenfrvect, int verb);
/*	dat must be created before
	dat.y must be allocated for a dimension at least equal to
		the value of ndat
	kframe is the sequence number of the frame in the file, starting from 1
	verb = (0,1,...) verbosity
*/


double mjd_frame(struct FrameH *fr);


int frame2ds_single(DS *d, struct FrFile *iFile, char *adcName, RING *r, GD *g, int verb);

int frame2ds(DS *d, struct pss_infile_db *dbf, char *adcName, RING *r, GD *g, int verb);


int read_fmnl_files0(struct pss_infile_db *dbf, int verb);


int check_infile_db(struct pss_infile_db* dbf, char *outfile);


struct tm *GPS2UTC(long GPSsec, int ULeapS);


#ifdef __cplusplus
}
#endif
