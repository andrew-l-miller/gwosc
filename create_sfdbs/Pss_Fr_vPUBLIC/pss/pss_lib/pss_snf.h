/*___________________________________________________
 ¦                                                   ¦
 ¦                    pss_snf.h                      ¦
 ¦       by Sergio Frasca - ((( 0 ))) Virgo          ¦
 ¦                    March 2000                     ¦
 ¦___________________________________________________¦*/


#ifdef __cplusplus
extern "C" {
#endif

typedef struct SNF_ACCESS{
	short	verb;		/* verbosity level */
	short	type;		/* 0 -> simple sequential 
						   1 -> variable length sequential
						   2 -> simple direct
						   3 -> variable length direct */
	char**	ch_select;	/* char array with the selected channels (for output) */
	long*	ich_select;	/* vector containing sequence numbers of output  channels */
	long	nrec;		/* number of the record */
	long	nrecextra;	/* number of the extra record */
	double	trec;		/* record time (typically in mjd) */
	double	dtrec;		/* time duration of a record (typically in sec) */
	long	inirec;		/* start of the new record */
	long*	selbias;	/* vector containing offsets in the record for the 
							output  channels */
} SNF_ACCESS;



typedef struct SNF_GD{    /* GD R_substructure */
	char	name[16];
	long	n;
	double	ini;
	double	dx;
	short	type;		/* =1 virtual abscissa, =2 real abscissa */
	short	complex;	/* =1 complex, =0 real */
	
	char	nform[8];   /* numerical format ('ascii' 'int8' 'int16' 'log8'
					     'log16' 'float' 'double') */
	char	sform[10];  /* sparse vector format ('sparse' 'bsparse' 'dsparse'
					     'nsparse' 'nosparse') */
} SNF_GD;


typedef struct SNF_DM{    /* DM R_substructure */
	char	name[16];
	long	n;
	double	ini;
	double	dx;
	long	m;
	double	ini2;
	double	dx2;
	short	type;		/* =1 virtual abscissa, =2 real abscissa */
	short	complex;	/* =1 complex, =0 real */
	
	char	nform[8];   /* numerical format ('ascii' 'int8' 'int16' 'log8'
					     'log16' 'float' 'double') */
	char	sform[10];  /* sparse vector format ('sparse' 'bsparse' 'dsparse'
					     'nsparse' 'nosparse') */
} SNF_DM;


typedef struct SNF5{/*  */
	char	name[16];
	char	caption[70];
	long	ivalue;
	double	dvalue;
	char	svalue[16];
} SNF5;


typedef struct SNF_HEADV{
	long	lheadi;
	long	lheadd;
	long	lheads;
	long	lrheadi;
	long	lrheadd;
	long	lrheads;
	long	lfheadi;
	long	lfheadd;
	long	lfheads;

	SNF5	*hi;
	SNF5	*hd;
	SNF5	*hs;
	SNF5	*rhi;
	SNF5	*rhd;
	SNF5	*rhs;
	SNF5	*fhi;
	SNF5	*fhd;
	SNF5	*fhs;
} SNF_HEADV;


typedef struct SNF_FIELD{	  /* e.g. DS channel substructure */
	unsigned long	fieldlen; /* number of bytes */
	char			name[20];
	long			len;      /* number of elements */
	double			dt;
	double			del;
	short			complex;  /* =1 complex, =0 real */
	char			nform[8];
	char			sform[8];
} SNF_FIELD;


typedef struct SNF_DS{		/* DS object substructure - OBSOLETE */
	short nch;

} SNF_DS;


typedef struct SNF_TFM{    /* TFM object substructure */
	char*	name;
	char	type[8];/* "sfdb", "spec", "peak" */
	long	len;	/* length of the band */
	long	inif;	/* starting point in the original fft (0 -> beginning) */
	double	df;		/* frequency step (start frequency inif*df */
} SNF_TFM;


typedef struct SNF_HM{		/* HM object substructure */
	char*	name;
	long	totlen;
	long	lamlen;
} SNF_HM;


typedef struct SNF_STRUCT{	/* Read/Write structure */
	char*	fileName;		/* input or output file */
	FILE*	file;

	long	headerlen;		/* file header length (in bytes) */
	long	reclen;			/* record length (in bytes) */
	long	reclenread;		/* record length as read in record header (in bytes) */
	long	nfields;		/* number of fields */

	char	protocol[10];	/* now is snf02 */
	char	userprot[16];	/* user protocol, for particular application program needs */
	
	char	binasc[6];		/* "bin" or "ascii" */
	char	bform[20];		/* bynary format (ieee-bigendian, ieee-littleendian,
							vaxd, vaxg,...) */

	char	obj[8];			/* object ('gd' 'dm' 'ds' 'tfm' 'tf_sfdb' 
							  'tf_spec' 'tf_peak' 'hm' 'psc' 'ev') */
	
	int		datlen;			/* data length (in byte; 1,2,4,8,...) */

	char	directaccess[15];
	
	char*	caption[122];   /* caption (SNF2 lines) */

	char	xunit[10];		/* abscissa unit (optional; may be 'seconds', 'mjd',
								'Hz') */

	long	nSNF2;			/* number of SNF2 lines */
	long	nSNF3;			/* number of SNF3 lines */
	long	nSNF5;			/* number of SNF5 lines */
	long	nSNF6;			/* number of SNF6 lines */
	char*	strSNF6[129];	/* SNF6 records as strings */

	SNF_ACCESS	*access;	/* read-write access substructure */
	SNF_HEADV	*headv;		/* rfhead substructure */
	
	SNF_FIELD	*field;

	SNF_GD		*gd;		/* GD substructure */
	SNF_DM		*dm;		/* DM substructure */
	SNF_DS		*ds;		/* DS substructure */
	SNF_TFM		*tfm;		/* TFM substructure */
	SNF_HM		*hm;		/* HM substructure */
} SNF_STRUCT;


/* --------------------------------------------- */

typedef struct LogX_STRUCT{
	short	stat;	/* format state variable */
	short	dimen;	/*  X : 1, 2, 4, 8, 16 or 32 */
	short	sign;	/* -1  all negative,  1  all positive,  0  mixed */
	short	linlog;	/*  3  look-up,   2  constant,   1  linear,   0  logarithmic */
	double	m;
	double	b;
	float	*lookup;/* look-up table for the tresholds */
	float	epsval;	/* for two-signs log automatic format */
	float	satur;	/* positive and negative saturation if > 0 */
	float	errmax;
	float	errcode;/* error code: 0 no error,... */
	long	len;	/* length of vec */
	float*	vec;
	long	nchar;  /* total number of bytes in coded data */
	char*	coded;  /* coded data array, with length len */
} LogX_STRUCT;


typedef struct SpVec_STRUCT{
	short	stat;     /* format state variable */
	char	dimen;	  /* 4, 8, 16, 32 (4*2^b) bits */
	char	binsp;    /* 0 non-sparse, 1 non-sparse binary, 2 sparse non-binary, 3 sparse binary */
	char	derive;   /* 0 no,   1 yes */
	char	logx;     /* LogX format for the non-zero elements: 0 no,   1 yes */
	long	lenvec;   /* length of uncoded or decoded vec array */
	float*	vec;	  /* normal vector (uncoded or decoded) */
	long	n0run;	  /* number of 0-runs */
	long	lenspar;  /* length of the sparse vector */ 
	long	nchar;    /* number of bytes of coded data array */
	char*	coded;    /* SpVec coded data array */
	float*	nzelem;   /* non-zero elements */
	long*	i4runs;	  /* should be allocated with enough values */
} SpVec_STRUCT;

/* ------------------------------------------------------------- */

struct LogX_STRUCT* write_LogX(struct LogX_STRUCT* str, short coding, FILE *file);
/* coding = 0  -> use input parameters 
          > 0  optimization level 
		  > 0  optimization on base, min and sign
		  = 2  optimization on dimen 
		  = 3  optimization on log-lin
		  = 4  full optimization */


struct LogX_STRUCT* read_LogX(struct LogX_STRUCT* str, FILE *file);


struct LogX_STRUCT* code_LogX(struct LogX_STRUCT* str, short coding);
/* coding = 0  -> use input parameters 
          > 0  optimization level 
		  > 0  optimization on base, min and sign
		  = 2  optimization on dimen 
		  = 3  optimization on log-lin
		  = 4  full optimization */


struct LogX_STRUCT* decode_LogX(struct LogX_STRUCT* str);


short code_LogX_stat(struct LogX_STRUCT* str);
/* Creates the stat variable with the coding information */


short decode_LogX_stat(struct LogX_STRUCT* str);
/* Reads the LogX stat variable and updates the LogX structure */


long nchar_LogX(struct LogX_STRUCT* str);


int show_LogX_coding(struct LogX_STRUCT* str);
/* shows coding parameters */



struct SpVec_STRUCT* write_SpVec(struct SpVec_STRUCT* str, short coding,
								 struct LogX_STRUCT* lxstr, short lxcoding, FILE *file);
/* coding = 0  -> use input parameters 
          > 0  -> optimize parameters

  lxstr, lxcoding = LogX structure and coding for the non-zero elements
*/

struct SpVec_STRUCT* read_SpVec(struct SpVec_STRUCT* str, FILE *file);


struct SpVec_STRUCT* code_SpVec(struct SpVec_STRUCT* str, short coding);
/* coding = 0  -> use input parameters 
          > 0  -> optimize parameters
*/

struct SpVec_STRUCT* decode_SpVec(struct SpVec_STRUCT* str);


short code_SpVec_stat(struct SpVec_STRUCT* str);
/* Creates the stat variable with the coding information */


short decode_SpVec_stat(struct SpVec_STRUCT* str);
/* Reads the LogX stat variable and updates the LogX structure */


int show_SpVec_coding(struct SpVec_STRUCT* str);
/* shows coding parameters */


SNF_STRUCT *snf_crea_struct(char* type, long mult);


int del_snf_struct(char* type);


/* ------------------- Open (write) -------------------------- */


SNF_STRUCT *more_comments(SNF_STRUCT* snf_str, char** comm, int ncomlin);


SNF_STRUCT* SNF5_lines_i(SNF_STRUCT* snf_str, int n, char** name, char** capt, long* ival);


SNF_STRUCT* SNF5_lines_d(SNF_STRUCT* snf_str, char** name, char** capt, double* val);


SNF_STRUCT* SNF5_lines_s(SNF_STRUCT* snf_str, char** name, char** capt, char** sval);


SNF_STRUCT* SNF5_lines_ri(SNF_STRUCT* snf_str, char** name, char** capt, long* ival);


SNF_STRUCT* SNF5_lines_rd(SNF_STRUCT* snf_str, char** name, char** capt, double* val);


SNF_STRUCT* SNF5_lines_rs(SNF_STRUCT* snf_str, char** name, char** capt, char** sval);


SNF_STRUCT* SNF5_lines_fi(SNF_STRUCT* snf_str, char** name, char** capt, long* ival);


SNF_STRUCT* SNF5_lines_fd(SNF_STRUCT* snf_str, char** name, char** capt, double* val);


SNF_STRUCT* SNF5_lines_fs(SNF_STRUCT* snf_str, char** name, char** capt, char** sval);


SNF_STRUCT* SNF6_lines(SNF_STRUCT* snf_str, char** lin, int nSNF6lin);


SNF_STRUCT *all_header_struct(SNF_STRUCT* snf_str);


SNF_STRUCT *all_snf_struct(SNF_STRUCT* snf_str);


int snf_open_w(SNF_STRUCT *w_str);


SNF_STRUCT *snf_all_ds(char *chname, long *chlen, double *dt, int nds, 
					   char *filename, char *caption, char **nform);
/*  
    Creates the structure for archiving a ds (or a multiple ds).

    chnames is a char array with dimension 20*nds containing, at multiples of 20,
	the names of the channels.
	chlen contains the number of data for each channel.
	These informations can be obtained by the ds objects or by the mch structure.
	
	If the LogX format is used, one or more  LogX structures should be 
	created out of this routine.

    This routine creates the easiest structure, that is with:
		- one caption line
		- the record time (rhd[0])
		- the fields delays (fhd[0])
	More complex features,as more caption lines or (file, record, field) header
	variables and non-default items should be added or modified after the call.
	
	Then the snf_open_w should be called.
*/


SNF_STRUCT *snf_open_tfm(char *filename, char *caption, char *nform, char *userprot,
						 double dt, char *typ, long lfft, long interlac, char *window,
						 int nband, long *blen, long *inif, double df);


SNF_STRUCT *snf_open_hm(SNF_STRUCT *w_str);


SNF_STRUCT *snf_open_psc(SNF_STRUCT *w_str);



/* ----------------- Write --------------------*/


int snf_write_rec(SNF_STRUCT *str, long *oldrecini);
/* writes only the header and complete the header of the previous record 
   at the first record *oldrecini should be set to 0
   the record header variables should be set before the call */ 


int snf_write_recextra(SNF_STRUCT *str, long *oldrecini);
/* writes only the header and complete the header of the previous record 
   at the first record *oldrecini should be set to 0 
   str->access->nrecextra should be 0 at the beginning;
   the record header variables should be set before the call */ 


int snf_write_field(SNF_STRUCT *str, long *prevfieldini);
/* writes only the header and complete the header of the previous field 
   at the beginning *prevfieldini should be set to 0;
   the field header variables should be set before the call */ 



int snf_write_vect(FILE* file, char* binasc, char* nform, char* sform, long n, float* f_v, 
				   LogX_STRUCT *logx, SpVec_STRUCT *spvect);


int snf_write_dvect(FILE* file, char* binasc, char* nform, char* sform, long n, double* d_v,
					   struct LogX_STRUCT *logx, struct SpVec_STRUCT *spvect);


int snf_write_lvect(FILE* file, char* binasc, char* nform, char* sform, long n, long* l_v);


int snf_write_ivect(FILE* file, char* binasc, char* nform, char* sform, long n, int* i_v);


int snf_write_gd(GD* gd, long ncomlin, char *binasc, char* filgd);
/*   
*/


int snf_write_dm(DM* dm, long ndm, char *binasc, char* fildm);
/* Writes data in a snf dm file; the header should be already written.
   The record and field header variables should be set before this call.
   The matrix is entered as a vector, by columns.

   In case of N dms (multiple files), one should do N calls.
*/


int snf_write_ds(SNF_STRUCT *w_str, float **dat, long rec, double tim,
				 struct LogX_STRUCT *logx);
/* writes data in a snf ds file; the header should be already written.
   Only one record is written.
   The record and field header variables should be set before this call.
   
   In case of N channels (multiple files), one should do a single call,
   with all the data for all the channels set serially in y (first all 
   the data of the first channel, then...).
*/


int snf_write_tfm_sfdb(SNF_STRUCT *snf_str, struct LogX_struct *logx);
/* Writes data in a snf sfdb file; the header should be already written.
   The record and field header variables should be set before this call.
   
   The data are entered as full or half ffts, one per call (one fft is 
   a record).
   A record can be also composed by (not adjacent) pieces of the same fft,
   resulting in a multiple file.
*/


int snf_write_tfm_spec(SNF_STRUCT *snf_str, struct LogX_STRUCT *logx);
/* writes data of entire time-frequency map in a snf tfm file; 
   the header should be already written.
   The record and field header variables should be set before this call.

   The data (periodograms), in str->vect, are archived with LogX format.
   
   The data are entered as half periodograms, one per call 
   (one fft is a record).

   A record can be also composed by (not adjacent) pieces of the same fft,
   resulting in a multiple file.
*/


int snf_write_tfm_peak(SNF_STRUCT *snf_str, struct SpVec_STRUCT *spvect);
/* writes data of a time-frequency peak map in a snf tfm file; 
   the header should be already written.
   The record and field header variables should be set before this call.

   The data are archived with SpVect. In this case obviously the 
   records are not all equal. 
   
   The data are entered as half periodograms, one per call 
   (one fft is a record).
   A record can be also composed by (not adjacent) pieces of the same fft,
   resulting in a multiple file.
*/


int snf_write_hm(SNF_STRUCT *snf_str, double* y);
/* writes data in a snf hm file; the header should be already written.
   The record and field header variables should be set before this call.
*/


/* ----------------- Open read --------------------*/



SNF_STRUCT *snf_open_r(char *file, int verb);


int test_nidg(SNF_STRUCT *r_str);


int snf_analyze(char *fileName, short verb, char *typ);
/*	verb	verbosity
	typ		type of particular analysis:
				"time"		time consistency
*/


int snf_show_header(SNF_STRUCT *str);


int snf_show_recs(SNF_STRUCT *str, int verb);


int snf_show_recfields(SNF_STRUCT *str, int verb);


int snf_show_fields(SNF_STRUCT *str, int verb);


/* ----------------- Read --------------------*/


int snf_read_rec(SNF_STRUCT *str);
/* 
	return variable :
	  0		error
	 >0		position of next record, ordinary record
	 <0		-position of next record, extraordinary record
	 -1		last record

  if ASCII, return 1
*/

int snf_read_field(SNF_STRUCT *str);
/* 
	return variable: 
	
	  the position of next field

	 -1		last field of the file

  if ASCII, return 1
*/


int snf_get_rec(SNF_STRUCT str);


float *snf_read_vect(FILE* file,char* binasc, char* nform, char* sform, long n, 
				   struct LogX_STRUCT *logx, struct SpVec_STRUCT *spvect);
/*
	logx, spvect	externally allocated only if necessary (depending on nform and sform)
*/


double *snf_read_dvect(FILE* file,char* binasc, char* nform, char* sform, long n, 
				   struct LogX_STRUCT *logx, struct SpVec_STRUCT *spvect);
/*
	logx, spvect	externally allocated only if necessary (depending on nform and sform)
*/


GD *snf_read_gd(char *file, long kgd, int verb)
/*
	file	file to open
	kgd		take the k-th (1,2,...,n) gd of the file (if multiple file) 
	verb	verbosity
*/;


int snf2ds(DS *ds, short kds, struct pss_infile_db *dbf,
			 RING *r, GD *g, int verb);
/*
	ds		the attached ds (in case of multiple ds, the element of the array of ds)
	kds		sequence number in case of multiple ds (the first is 0), 0 in case of single ds
	dbf		the structure with the infile data-base
	r,g		service structures, should be defined outside with all_ringxds and all_gdxds
	verb	verbosity level
*/


int snf_close(SNF_STRUCT *snf_str);


/*---------------------------- Templates -------------------------------*/


SNF_STRUCT template_write_ds();


SNF_STRUCT template_read_ds();


SNF_STRUCT template_write_tfm_sfdb();


SNF_STRUCT template_read_tfm_sfdb();


SNF_STRUCT template_write_tfm_spec();


SNF_STRUCT template_read_tfm_spec();


SNF_STRUCT template_write_tfm_peak();


SNF_STRUCT template_read_tfm_peak();


int template_snf2snf();


int template_r872snf();


int template_frame2snf();



/*---------------------------- Applications -------------------------------*/


int snf_file_extract(char *filein, char *fileout,
				int nitem, int *fields, long *minind, long *maxind,
				long minrec, long maxrec);


int snf_file_merge();


#ifdef __cplusplus
}
#endif

