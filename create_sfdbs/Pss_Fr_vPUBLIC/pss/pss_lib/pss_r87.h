/*___________________________________________________
 ¦                                                   ¦
 ¦                   pss_r87.h                       ¦
 ¦       by Sergio Frasca - ((( 0 ))) Virgo          ¦
 ¦                    June 2001                      ¦
 ¦___________________________________________________¦*/

/*			R87 DATA FORMAT


   The archived sampled data are collected in direct access files on disc.

   The files are composed of records, all of the same length. There are
   different types of records. The types are:

               1    sampled data
               2    physical parameters
               3    experimenters comments
               4    DAGA system setup

   All records have a header of at least 90 integer*2 words.

   The general structure of the records of the first type is 

                 Header
                 First Data Field
                 Second Data Field (optional)
                 Third Data Field (optional)

   Any data field has a particular sampling time.

   Any data field contains the samples of a certain number of channels, so
   the samples are, in the case of n channels,

                CH1 sample 1
                CH2 sample 1
                ............
                CHn sample 1
                CH1 sample 2
                ............

   The first sample of each channel of each field has the same time, that
   is in the header.

   There are different methods of representing the samples, using one or two
   bytes. In the case of the data of Explorer 1990-91, the samples are
   INTEGER*2 numbers and the sampling delta is 0.004883 V.

   The meaning of the header words is

    1  H_LREC  rec.length (in words)   2  H_LHEA  total header length (in w.)
    3  H_KREC  record number           4  H_ANTG  gravitational antenna
    5  H_KRUN  run number              6  H_TIPR  record type
    7  H_LCO1  length of data field 1  8  H_NCA1  field 1 channels number
       (total number of samples)
    9  H_LCO2  length of data field 2 10  H_NCA2  field 2 channels number
   11  H_LCO3  length of data field 3 12  H_NCA3  field 3 channels number
   13  H_ANNO  year                   14  H_DOY   day of the year
   15  H_ORE   hour                   16  H_MIN   minutes
   17  H_SEC   seconds                18  H_MSEC  milliseconds
   19          sampl. ms for field 1  20          sampling microsec. field 1
   21  H_OFTIM time offset  (ms)      22  H_ADC   A/D code (tab.)
   23  H_LHE1  header field 1 length  24  H_LHE2  header field 2 length
   25  H_OPFL  operation flag         26   -
   27  H_NSIN  number of sintet.      28  Freq. sint. 1 (Hz)
   29  Freq. sint. 1 (10^-4 Hz)       30         "    2
   31        "     2                  32         "    3
   33        "     3                  34         "    4
   35        "     4                  36  H_AMP   electr.ampl. [F16(1000)]
   37  -----  H_CDIR   coeff. V^2 -> K   direct acquisit. [in RHEADER(19)]
   39  -----  H_CRIS1    "        "    resonance nu_minus [in RHEADER(20)]
   41  -----  H_CRIS2    "        "        "     nu_plus  [in RHEADER(21)]
   43  V_0     trasducer voltage (V)  44  V_1
   45  V_2                            46  V_3
   47  V_4                            48  V_5
   49  V_6                            50  MUX_IND
   51  KFFT_K  n. of AAP FFT          52  KFFT_L  lunghezza FFT
   53  KFFT_NP number of AAP pieces   54  KFFT_I1 beginning of piece 1
   55  KFFT_L1 length of piece 1      56  KFFT_I2 beginning of piece 2
   57  KFFT_L2 length of piece 2      58  KFFT_I3 beginning of piece 3
   59  KFFT_L3 length of piece 3      60  KFFT_NN n. of FFT in this rec (1 o 2)
   61  H_TAU1  tau nu(-)[F16(.01 s)]  62  H_TAU2    tau nu(+) [F16(.01 s)]
   63  H_SIG1 sigma  " [F16(.0001 V)] 64  H_SIG2  sigma   "   [F16(.0001 V)]
   65  H_S0 sqr.den.r.elet.[F16(1E-6)]66  H_VART  total var. [F16(.0001 V^2)]
   67  H_TEQR T_eq_rison [F16(1E-5 K)]68  H_TEQWB T_eq_wide_band [F16(1E-5 K)]
   69  H_TEF1 T_eff nu(-)[F16(1E-6 K)]70  H_TEF2  T_eff nu(+) [F16(1E-6 K)]
   71  H_TEFMF T_eff M.F.[F16(1E-6 K)]72  H_KSPET spectrum number of the run
   73  -----  H_CAL      calibration value by lock-in   [in RHEADER(37)]
   75  -----  H_CALSPET       "        "   by spectrum  [in RHEADER(38)]
   77  H_TAUIMP_1  imposed tau_1      78  H_TAUIMP_2  imposed tau_2 
             [F16(.01 s)]                     [F16(.01 s)]
   79  H_FBIMP  imposed F.B.          80
          [F16(.000001 V/sqr(Hz))]
   81  H_CL1_CORR (dms)               82 H_CL2_CORR (dms)
   83  H_D2TIM_CORR (dms)             84
   85                                 86
   87                                 88
   89                                 90
   91                                 92
   93                                 94
   95                                 96  
   97                                 98  
   99                                100 

                     N_MISx = H_LCOx / H_NCAx
                     H_TSAM = H(19)*1000 + H(20)
                     H_TSA1 = H_TSAM
                     H_TSA2 = H_TSAM*N_MIS2/N_MIS1
                     H_TSA3 = H_TSAM*N_MIS3/N_MIS1
                     F_SINx = synthesizer frequency x (real)

   In the case of the data of Explorer 1990-91, the header values were

	 H_LREC      3000             
	 H_LHEA       200             
	 H_LCO1      1200             
	 H_NCA1        12             
	 H_LCO2      1600             
	 H_NCA2         1             
	 H_LCO3         0             
	 H_NCA3         1             
	 N_MIS1       100             
	 N_MIS2      1600             
	 N_MIS3         0             

	 H_KRUN       192       (example)      
	 H_KREC         1           "
	 H_TIPR         2           " 
	 H_ANNO      1991           " 
	 H_DOY        206           " 
	 H_ORE         10           " 
	 H_MIN          4           " 
	 H_SEC         31           " 
	 H_MSEC       850           " 
	 H_TSA1    290816           
	 H_TSA2     18176             
	 H_TSA3         0             
	 H_ADC          0             

   Some of the header variables were not defined at that time.

   There were 12 channels in field 1 and 1 in field 2:

                   Field 1    [sampling time 290816 us]

   1            _X-             X nu -          
   2            _Y-             Y nu -          
   3            _X+             X nu +          
   4            _Y+             Y nu +          
   5            CAL             Calibration     
   6            DIR     Direct Acquisition      
   7            ROB             R out reson     
   8            SLB             Sism l.freq.    
   9            SRF     Sism resonance          
  10            EMM     Electro-magn. mon.      
  11            CLK             Top second      
  12            ALF     Antenna low freq. or Jolly

                  Field 2    [sampling time 18176 us]

   1            AAP     Direct Acquisition after Anti-Aliasing Procedure

*/
/*
%READ_HEADER_R87  builds the header structure for an R87 record
%
%          header=read_header_r87(A)
%
%        A  the R87 record
%
%           Header structure
%
%     header.len         record length (in I16 words) [=0 -> end of file]
%     header.headlen     header length      "
%     header.recnum      record number
%     header.antenna     gravitational antenna code
%     header.runnum      run number
%     header.type        record type (1 sampled data, 2 physical parameters,
%                        3 experimenters comments, 4 DAGA setup (script)
%     header.len1        length of field 1
%     header.nc1         number of channels of field 1
%     header.len2        length of field 2
%     header.nc2         number of channels of field 2
%     header.len3        length of field 3
%     header.nc3         number of channels of field 3 
%     header.time        time (in MJD)
%     header.st          sampling time (in seconds)
%
%     header.adccode     ADC code
%     header.opflag      operational flag
*/
#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>




typedef struct r87_file{
	FILE *fil;
	long reclen;
	int initime[6];
	double samptim;
	short *A;
} r87_file;


typedef struct r87_header{
	long reclen;
	long headlen;
	long recnum;
	long antenna;
	long runnum;
	long type;

	long len1;
	long nc1;
	long len2;
	long nc2;
	long len3;
	long nc3;

	double time;
	double st;

	long adccode;
	long opflag;
	
	long field;
	long chan;
	long len;
	long nc;
	long ndata;
	double tsamp;
} r87_header;


typedef struct r87_inq_ch{
	long nch;
	char name[1000];
	char comment[2000];
	long num[50];
	double samp[50];
	long len[50];
	double tin;
} r87_inq_ch;


typedef struct r87_chdata{
	r87_header header;
	float *data;
} r87_chdata;


r87_file open_r87(char *fileName,int verb);


r87_header read_header_r87(short *A);


r87_inq_ch r87_inquire_ch(char *fileName,int verb);


r87_chdata read_r87rec_ch(FILE *fil,long reclen,long ch);


r87_chdata get_r87_1ch_vect(char *fileName, long chNum,
					 long krec, long ndat, int verb);

char *read_r87_info(char *fileName, long *numchar, int verb);

double mjd_r87(long *it);

#ifdef __cplusplus
}
#endif
