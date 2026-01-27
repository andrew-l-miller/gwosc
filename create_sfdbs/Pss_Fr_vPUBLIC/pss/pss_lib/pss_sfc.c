/*___________________________________________________
 ¦                                                   ¦
 ¦                    pss_sfc.c                      ¦
 ¦       by Sergio Frasca - ((( 0 ))) Virgo          ¦
 ¦                    June 2003                      ¦
 ¦___________________________________________________¦*/


/*
                           Simple File Format Collection

 The basic feature of the file formats here collected is the ease of access to the data.
 "ease of access" means:
		- the software to access the data consists in a few lines of basic code
		- the data can be accessed easily by any environment and language
		- the byte level structure is immediately intelligible
		- no unneeded information is present
		- the number of pointers and structures is minimized
		- the structure fits the needs
		- the need for generality is tempered by the need for easyness

  The collection is composed by:
		- sds		simple data stream format, for finite or "infinite" number of equispaced
					samples, in one or more channels, all with the same sampling time
		- 

One or more files can be:
		- collected		logically linked by a collection script file
		                (for example subsequent data files, or to put together different
						sampling time channels)
		- embedded		in a single file, with a toc at the beginning
		- wrapped		by adding one or more external headers
						(for example describing the computer which wrote the file)

--------------------------------------------------------------------------------------------

                        SFC file format

    *** SDS and SBL are some of the file formats of the SFC collection ***

  The format is the following :

   label    - 8 bytes string (1-8)         type label (e.g. #SFC#SDS)
   prot     - 4 bytes int32  (9-12)        protocol (integer; now 1)
   coding   - 4 bytes int32  (13-16)       machine and data coding (normally 0)

   inidat   - 4 bytes int32  (17-20)       pointer to the beginning of data (>= 944+N*128)
   nch      - 4 bytes int32  (21-24)       number of channels (N; integer)
   len      - 8 bytes int64  (25-32)       number of blocks (sbl) or data per channel (sds)
                                             (0 if it is not known)

   t0       - 8 bytes double (33-40)       beginning time (double); may be not meaningful
   dt       - 8 bytes double (41-48)       sampling time (double); may be not meaningful

   capt     - 128 bytes string             caption
   filme    - 128 bytes string             original name of the file
   filmaster- 128 bytes string             master file
   filspre  - 128 bytes string             serial preceding file
   filspost - 128 bytes string             serial subsequent file
   filppre  - 128 bytes string             parallel preceding file
   filppost - 128 bytes string             parallel subsequent file

     Here ends the general part, read by this function; then:

   ch()     - N*128 bytes                  depends on the sfc type
                            *** The channels can be sampled data, single parameters,
                                sets of parameters, matrices, strings, etc.

   user     - (free)                       user file header
                            *** The user header length is (inidat - 944+N*128) bytes
                            
   data     - ...                          the data (depends on the sfc type)
                            *** The data are in N parallel streams (sds) or divided
                                in blocks (sbl)

  To access the user header, fseek(fid,944+N*128,'bof')

                        sfc_ structure

    sfc_.file      file name (with path)
    sfc_.pnam      file path
    sfc_.fid       fid
    sfc_.label     label (e.g. #SFC#SDS)
    sfc_.prot      protocol
    sfc_.coding    machine and data coding
    sfc_.t0        beginning time
    sfc_.dt        sampling time
    sfc_.capt      caption
    sfc_.nch       number of channels
    sfc_.ch(nch)   channel structures (depends on type)
    sfc_.hlen      user header length (bytes)
    sfc_.len       length of a stream (the total number of data is len*nch)
                   or number of blocks
    sfc_.point0    pointer to beginning of data
    sfc_.eof       end of file (1; for chained files, 2 -> end of chain, 3 -> error,
                    -1 -> end chosen period)
    sfc_.acc       access number; 0 at beginning, incremented by user if 
                    not accessed in standard ways
    sfc_filme      original name of the file
    sfc_filmaster  master file
    sfc_filspre    serial preceding file
    sfc_filspost   serial subsequent file
    sfc_filppre    parallel preceding file
    sfc_filppost   parallel subsequent file

*/

#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

//#include "../v6r19/src/FrameL.h"  //path cambiato da pia. Ho v6r19 sotto pss.
                                 //Aggiunta pia 10 Agosto 2005: dava l' errore seguente
/*
../pss/PSS_lib/pss_frame.h:61: warning: `struct FrameH' declared inside parameter list
../pss/PSS_lib/pss_frame.h:61: warning: its scope is only this definition or declaration, which is probably not what you want */

#include "pss_snag.h"
#include "pss_serv.h"
//#include "pss_frame.h"
#include "pss_sfc.h"


SFC_* sfc_open(char *file){

	SFC_* sfc_;
	FILE* fid;
	char* label;
	char* strapp;
	int nnn;

	int i;
	int ii;

	sfc_=malloc(sizeof(struct SFC_));
	label=(char*)malloc(9); 
	strapp=(char*)malloc(129);
	
	
	//pia: fatta str_trim1 (in pss_serv) e sostituita a str_trim  16 Aug 2005
	//pia: fatta path_from_file1 (in pss_serv) e sostituita a path_from_file  17 Aug 2005
	sfc_->acc=0;
	sfc_->eof=0;
	strcpy(sfc_->file,file); 
	ii=path_from_file1(file);
	strapp[0]=0;
	strncpy(strapp,file,ii);  
	strapp[ii]='\0'; //modifica pia e fede 14/05/07
	str_trim1(strapp);
	nnn=strlen(strapp);
	strcpy(sfc_->pname,strapp);
	sfc_->pname[nnn]='\0'; //modifica pia e fede 14/05/07

	//fid=malloc(sizeof(FILE));
	
	fid=fopen(sfc_->file,"rb");
	sfc_->fid=fid;
	if(fid == NULL){
		printf(" *** %s  non-existent file ! \n",file);
		sfc_->eof=3;
		return sfc_;
	}


	fscanf(fid,"%8s",label);
	str_trim1(label);
	strcpy(sfc_->label,label);
	fread(&sfc_->prot,4,1,fid);
	fread(&sfc_->coding,4,1,fid); 

	fread(&sfc_->nch,4,1,fid);
	fread(&sfc_->point0,4,1,fid);
	fread(&sfc_->len,8,1,fid);

	fread(&sfc_->t0,8,1,fid);
	fread(&sfc_->dt,8,1,fid);

	fscanf(fid,"%128s",strapp);
	str_trim1(strapp);
	strcpy(sfc_->capt,strapp);

	fscanf(fid,"%128s",strapp);
	str_trim1(strapp);
	strcpy(sfc_->filme,strapp);
	fscanf(fid,"%128s",strapp);
	str_trim1(strapp);
	strcpy(sfc_->filmaster,strapp);
	fscanf(fid,"%128s",strapp);
	str_trim1(strapp);
	strcpy(sfc_->filspre,strapp);
	fscanf(fid,"%128s",strapp);
	str_trim1(strapp);
	strcpy(sfc_->filspost,strapp);
	fscanf(fid,"%128s",strapp);
	str_trim1(strapp);
	strcpy(sfc_->filppre,strapp);
	fscanf(fid,"%128s",strapp);
	str_trim1(strapp);
	strcpy(sfc_->filppost,strapp);

	sfc_->hlen=sfc_->point0-(944+128*sfc_->nch);

	sfc_->ch=malloc(sfc_->nch*sizeof(char*));
	
	for(i=0;i<sfc_->nch;i++){
	        sfc_->ch[i]=(char*)malloc(129);
	      	fgets(strapp,129,fid); //printf("\n *** %s *** \n",strapp); 
		//fscanf(fid,"%128s",strapp); //NOTE THAT fscanf DOESN'T WORK ! WORKAROUND WITH fgets ; NOTE 129
		str_trim1(strapp);
		strcpy(sfc_->ch[i],strapp);
	}
	

	return sfc_;
}

void  sfc_open1(char *file,SFC_ *sfc_){
/**pia version, 23 Aug. 2005:does not allocate SFC_ again !!**/
        //pia SFC_* sfc_;
	FILE* fid;
	char* label;
	char* strapp;


	int i;
	int ii;
	int nnn; //modifica Pia e fede 15 Maggio 2007

	// it is allocated only using sfc_open, at the beginning.  sfc_=malloc(sizeof(struct SFC_));
	label=(char*)malloc(9); 
	strapp=(char*)malloc(129);
	
	
	//pia: fatta str_trim1 (in pss_serv) e sostituita a str_trim  16 Aug 2005
	//pia: fatta path_from_file1 (in pss_serv) e sostituita a path_from_file  17 Aug 2005
	sfc_->acc=0;
	sfc_->eof=0;
	strcpy(sfc_->file,file); 
	ii=path_from_file1(file);
	strapp[0]=0;
	strncpy(strapp,file,ii);
	strapp[ii]='\0'; //modifica pia e fede 14/05/07 
	str_trim1(strapp);
	nnn=strlen(strapp);
	strcpy(sfc_->pname,strapp);
	sfc_->pname[nnn]='\0'; //modifica pia e fede 14/05/07

	//fid=malloc(sizeof(FILE));
	
	fid=fopen(sfc_->file,"rb");
	sfc_->fid=fid;
	if(fid == NULL){
		printf(" *** %s  non-existent file ! \n",file);
		sfc_->eof=3;
		return; //pia sfc_;
	}


	fscanf(fid,"%8s",label);
	str_trim1(label);
	strcpy(sfc_->label,label);
	fread(&sfc_->prot,4,1,fid);
	fread(&sfc_->coding,4,1,fid); 

	fread(&sfc_->nch,4,1,fid);
	fread(&sfc_->point0,4,1,fid);
	fread(&sfc_->len,8,1,fid);

	fread(&sfc_->t0,8,1,fid);
	fread(&sfc_->dt,8,1,fid);

	fscanf(fid,"%128s",strapp);
	str_trim1(strapp);
	strcpy(sfc_->capt,strapp);

	fscanf(fid,"%128s",strapp);
	str_trim1(strapp);
	strcpy(sfc_->filme,strapp);
	fscanf(fid,"%128s",strapp);
	str_trim1(strapp);
	strcpy(sfc_->filmaster,strapp);
	fscanf(fid,"%128s",strapp);
	str_trim1(strapp);
	strcpy(sfc_->filspre,strapp);
	fscanf(fid,"%128s",strapp);
	str_trim1(strapp);
	strcpy(sfc_->filspost,strapp);
	fscanf(fid,"%128s",strapp);
	str_trim1(strapp);
	strcpy(sfc_->filppre,strapp);
	fscanf(fid,"%128s",strapp);
	str_trim1(strapp);
	strcpy(sfc_->filppost,strapp);

	sfc_->hlen=sfc_->point0-(944+128*sfc_->nch); 
	
	for(i=0;i<sfc_->nch;i++){
	        //pia removed the allocation  sfc_->ch[i]=malloc(129);
	      	fgets(strapp,129,fid); //printf("\n *** %s *** \n",strapp); 
		//fscanf(fid,"%128s",strapp); //NOTE THAT fscanf DOESN'T WORK ! WORKAROUND WITH fgets ; NOTE 129
		str_trim1(strapp);
		strcpy(sfc_->ch[i],strapp);
	}
	

	return; // sfc_;
}



SFC_* sfc_openw(char* file, struct SFC_* sfc_){
/* 
*/
	FILE* fid;
	char* strapp;
	char *nofile="#NOFILE";
	//int prot=1,coding=1;
	int nnn;  //modifica pia e federica 14/05/07

	strapp=(char*)malloc(129);

	sfc_->acc=0;
	strcpy(sfc_->file,file); printf(" *** %s  - %s - %ld \n",file,sfc_->file,sfc_->acc); //pia ld

	fid=fopen(file,"wb");
	sfc_->fid=fid;
	 
	strapp=path_from_file(file); 
	nnn=strlen(str_trim(strapp)); //modifica pia e fede 14/05/07
	strcpy(sfc_->pname,str_trim(strapp));
	sfc_->pname[nnn]='\0'; //modifica pia e fede 14/05/07
	 
	strapp=filename_from_file(file);
	strcpy(sfc_->filme,str_trim(strapp));

	fprintf(fid,"%8s",sfc_->label);printf("\n%8s\n",sfc_->label);
	fwrite(&sfc_->prot,4,1,fid);
	fwrite(&sfc_->coding,4,1,fid);

	fwrite(&sfc_->nch,4,1,fid);
	sfc_->point0=944+128*sfc_->nch;
	fwrite(&sfc_->point0,4,1,fid);
	fwrite(&sfc_->len,8,1,fid);

	fwrite(&sfc_->t0,8,1,fid);
	fwrite(&sfc_->dt,8,1,fid);
	
	fprintf(fid,"%128s",sfc_->capt);

	fprintf(fid,"%128s",sfc_->filme);
	fprintf(fid,"%128s",nofile);
	fprintf(fid,"%128s",nofile);
	fprintf(fid,"%128s",nofile);
	fprintf(fid,"%128s",nofile);
	fprintf(fid,"%128s",nofile);

	strcpy(sfc_->filmaster,nofile);
	strcpy(sfc_->filppre,nofile);
	strcpy(sfc_->filppost,nofile);
	strcpy(sfc_->filspre,nofile);
	strcpy(sfc_->filspost,nofile);

	return sfc_;
}



void sfc_show(SFC_* sfc_){

	char* strapp;
	int i;

	strapp=(char*)malloc(130);

	printf(" Label    :  %s \n",sfc_->label);
	printf(" Protocol :  %ld \n",sfc_->prot); //pia ld
	printf(" Coding   :  %ld \n\n",sfc_->coding);
	printf(" Channels :  %ld \n",sfc_->nch);
	printf(" Point0   :  %ld \n",sfc_->point0);
	printf(" Length   :  %lld \n\n",sfc_->len); //pia lld fino a qui
	printf(" Time     :  %f \n",sfc_->t0);
	printf(" Sampling :  %f \n\n",sfc_->dt);
	printf(" Caption  :  %s \n\n",sfc_->capt);
	printf(" Fileme   :  %s \n",sfc_->filme);
	printf(" FileMast :  %s \n",sfc_->filmaster);
	printf(" FileSpre :  %s \n",sfc_->filspre);
	printf(" FileSpost:  %s \n",sfc_->filspost);
	printf(" FilePpre :  %s \n",sfc_->filppre);
	printf(" FilePpost:  %s \n\n",sfc_->filppost);

	for(i=0;i<sfc_->nch;i++){
		printf("   channel %4d :  %s \n",i+1,sfc_->ch[i]);
	}

	printf("\n Path     :  %s \n",sfc_->pname);

	//pia: fatta filename_from_file1 (in pss_serv) e sostituita a fileneme_from_file  22 Aug 2005
	//strapp=filename_from_file(sfc_->file);
	i=filename_from_file1(sfc_->file);
	strcpy(strapp,sfc_->file+i);
	str_trim1(strapp);
	printf(" File     :  %s \n",str_trim(strapp));
}


/* ------------------------  sds Routines ---------------------------*/

/*
                         SDS file format
%
    *** SDS is one of the file formats of the SFC collection ***
%
  The SDS (Simple DS) format is intended for storing one or more data
  streams, with the same time beginning and the same sampling time.
  
  The format is the following :
%
   label    - 8 bytes string (1-8)         #SFC#SDS
%
  First 944 bytes : see the sfc_open function.
%
   ch()     - N*128 bytes strings          names and captions for the channels
   user     - (free)                       user header
   data     - N*4*len bytes float          the data (l for each channel; float)
%
  To access the user header, fseek(fid,944+N*128,'bof')
%
  The sds_ structure derives from the sfc_ structure.
  The added or modified members are:
%
    sds_.ch(nch)   channel captions
    sds_.len       length of a stream (the total number of data is len*nch)
    sds_.point     pointer for next data

*/



SFC_* sds_open(char *file, struct ALLPERS_ *allpers){
/*  

  allpers is a double array of dimension nper that contains the allowed periods (typically in mjd).
  If nper = 0 , pers not considered. allpers features works if sds_->len is correctly set.

*/
	SFC_ *sds_;
	int ck0;
	double tim0,tim1;
	//FILE *fid;  //pia tolto non e' usato

	ck0=0;

	//sds_=malloc(sizeof(FILE));

	while(ck0==0){
		sds_=sfc_open(file);
		ck0=5;

		if(sds_->fid==NULL){
			return sds_;
		}

		if(allpers->nper>0){
			tim0=sds_->t0;
			tim1=sds_->t0+sds_->len*sds_->dt/86400;
			ck0=sds_check_time(tim0,tim1,allpers);

			if(ck0==0){
				fclose(sds_->fid);

				if(strcmp(sds_->filspost,"#NOFILE")){
					sds_->eof=2;
					printf(" --------> End of files ! \n");
					return sds_;
				}

				file=(char*)malloc(160);

				strcpy(file,sds_->pname);
				strcat(file,sds_->filspost);
			}
		}
	}

	sds_->eof=0;

	return sds_;
}

void sds_open1(char *file, struct ALLPERS_ *allpers,SFC_ *sds_){
/**pia version, 23 Aug. 2005: to call sfc_open1 and not sfc_open:
does not allocate SFC_ again !!**/
/*  

  allpers is a double array of dimension nper that contains the allowed periods (typically in mjd).
  If nper = 0 , pers not considered. allpers features works if sds_->len is correctly set.

*/
        //pia SFC_ *sds_;
	int ck0;
	double tim0,tim1;
	//FILE *fid; //pia tolto non e' usato

	ck0=0;

	//sds_=malloc(sizeof(FILE));

	while(ck0==0){
	        sfc_open1(file,sds_);  //pia: not sfc_open(file);
		ck0=5;

		if(sds_->fid==NULL){
		       return; // pia sds_;
		}

		if(allpers->nper>0){
			tim0=sds_->t0;
			tim1=sds_->t0+sds_->len*sds_->dt/86400;
			ck0=sds_check_time(tim0,tim1,allpers);

			if(ck0==0){
				fclose(sds_->fid);

				if(strcmp(sds_->filspost,"#NOFILE")){
					sds_->eof=2;
					printf(" --------> End of files ! \n");
					return; //pia sds_;
				}

				file=(char*)malloc(160);

				strcpy(file,sds_->pname);
				strcat(file,sds_->filspost);
			}
		}
	}

	sds_->eof=0;

	return; //pia  sds_;
}


SFC_* sds_openw(char* file, struct SFC_* sds_){
/*
*/
  printf(" *** NOT YET IMPLEMENTED ! file= capt= %s %s \n",file,sds_->capt); //aggiunto printf pia
  return sds_; //pia aggiunto, chiedere a Sergio: non ritornava nulla !
}


int sds_check_time(double tim0, double tim1, struct ALLPERS_ *allpers){

        int i; //pia

	printf(" *** NOT YET IMPLEMENTED !"); //pia
	printf(" tim0 tim1 allpers->nper %lf %lf %d\n",tim0,tim1,allpers->nper); //pia
	return i; //pia
}



void vec_from_sds(float *buffer, float *vec, double *tim0, struct HOLES_* holes,
			struct SFC_* sds_, int chn, long len, struct ALLPERS_ *allpers){

/* buffer should be allocated for len*nch 

    chn is the serial number of the channel (1,2,...)

	allpers->nper should be put equal to 0, if allowed periods are not used

*/

	int kskipped;
	long pos,ndat,len1,len2,kvec,count,nread,i,hole;
	double timend_exp;
	char *file;

	sds_->acc++;
	sds_->eof=0;
	holes->nztot=0;
	holes->nholes=0;
	kskipped=0;

	if(sds_->fid==NULL){
		printf(" *** End of files or error ! \n");
		return;
	}

	pos=ftell(sds_->fid);
	ndat=(pos-sds_->point0)/(sds_->nch*4);
	*tim0=sds_->t0+sds_->dt*ndat/86400;

	timend_exp=*tim0+len*sds_->dt/86400;

	len2=0;
	len1=len-len2;
	kvec=0;

	while(len1>0){
	    count=fread(buffer,4,sds_->nch*len1,sds_->fid);
	 	nread=count/sds_->nch;
		len2+=nread;
		len1=len-len2;
	
		if(nread>0){
			for(i=kvec;i<len2;i++)vec[i]=buffer[(i-kvec)*sds_->nch+chn-1];				     
	         	kvec=len2;		      
		}

		if(len1>0){
			sds_->eof=1;
			fclose(sds_->fid);
			if(strcmp(sds_->filspost,"#NOFILE")==0){
				sds_->eof=2;
				printf(" ---------> End of concatenated files \n");
				for(i=kvec;i<len;i++)vec[i]=0;
				len2=len;
				len1=0;
				holes->nholes++;
				holes->nzeros[holes->nholes-1]=len-kvec+1;
				holes->nztot+=holes->nzeros[holes->nholes-1];
				holes->kzeros[holes->nholes-1]=kvec;

				return;
			}
			else{
				file=strcat(sds_->pname,sds_->filspost);
				sds_=sds_open(file,allpers);

				if(sds_->eof==3){
					printf(" ---------> End of files \n");
					for(i=kvec;i<len;i++)vec[i]=0;
					len2=len;
					len1=0;
					holes->nholes++;
					holes->nzeros[holes->nholes-1]=len-kvec+1;
					holes->nztot+=holes->nzeros[holes->nholes-1];
					holes->kzeros[holes->nholes-1]=kvec;

					return;
				}
				printf(" *** open file %s \n",sds_->filme);
				hole=floor((sds_->t0-*tim0)*86400/sds_->dt-len2+0.5);

				if(hole>len1){
					for(i=kvec;i<len;i++)vec[i]=0;
					len2=len;
					len1=0;
					holes->nholes++;
					holes->nzeros[holes->nholes-1]=len-kvec+1;
					holes->nztot+=holes->nzeros[holes->nholes-1];
					holes->kzeros[holes->nholes-1]=kvec;

				}
				else{
					if(hole>0){			
						for(i=kvec;i<kvec+hole;i++)vec[i]=0;
						len2=len2+hole;
						len1=len-len2;
						holes->nholes++;
						holes->nzeros[holes->nholes-1]=hole;
						holes->nztot+=holes->nzeros[holes->nholes-1];
						holes->kzeros[holes->nholes-1]=kvec;
						kvec=len2;
					}
				}
			}
		}

		if(len1<=0){
			if(allpers->nper>0){
				// TO BE DONE WITH FUNCTION ZERO_PERS
			}
		}
	}

	return;
}



void vec_from_sds_0(float *buffer, float *vec, double *tim0, struct HOLES_* holes,
			struct SFC_* sds_, int chn, long len, struct ALLPERS_ *allpers){

/* Version without the short hole embodiment

   buffer should be allocated for len*nch 

    chn is the serial number of the channel (1,2,...)

	allpers->nper should be put equal to 0, if allowed periods are not used

*/

	int kskipped;
	long pos,ndat,len1,len2,kvec,count,nread,i,hole;
	double timend_exp;
	char *file;

	sds_->acc++;
	sds_->eof=0;
	holes->nztot=0;
	holes->nholes=0;
	kskipped=0;

	if(sds_->fid==NULL){
		printf(" *** End of files or error ! \n");
		return;
	}

	pos=ftell(sds_->fid);
	ndat=(pos-sds_->point0)/(sds_->nch*4);
	*tim0=sds_->t0+sds_->dt*ndat/86400;

	timend_exp=*tim0+len*sds_->dt/86400;

	len2=0;
	len1=len-len2;
	kvec=0;

	while(len1>0){
	    count=fread(buffer,4,sds_->nch*len1,sds_->fid);
	 	nread=count/sds_->nch;
		len2+=nread;
		len1=len-len2;
	
		if(nread>0){
			for(i=kvec;i<len2;i++)vec[i]=buffer[(i-kvec)*sds_->nch+chn-1];				     
	         	kvec=len2;		      
		}

		if(len1>0){
			sds_->eof=1;
			fclose(sds_->fid);
			if(strcmp(sds_->filspost,"#NOFILE")==0){
				sds_->eof=2;
				printf(" ---------> End of concatenated files \n");
				for(i=kvec;i<len;i++)vec[i]=0;
				len2=len;
				len1=0;
				holes->nholes++;
				holes->nzeros[holes->nholes-1]=len-kvec+1;
				holes->nztot+=holes->nzeros[holes->nholes-1];
				holes->kzeros[holes->nholes-1]=kvec;

				return;
			}
			else{
				file=strcat(sds_->pname,sds_->filspost);
				sds_=sds_open(file,allpers);

				if(sds_->eof==3){
					printf(" ---------> End of files \n");
					for(i=kvec;i<len;i++)vec[i]=0;
					len2=len;
					len1=0;
					holes->nholes++;
					holes->nzeros[holes->nholes-1]=len-kvec+1;
					holes->nztot+=holes->nzeros[holes->nholes-1];
					holes->kzeros[holes->nholes-1]=kvec;

					return;
				}
				printf(" *** open file %s \n",sds_->filme);
				hole=floor((sds_->t0-*tim0)*86400/sds_->dt-len2+0.5);

				if(hole>len1){
					for(i=kvec;i<len;i++)vec[i]=0;
					len2=len;
					len1=0;
					holes->nholes++;
					holes->nzeros[holes->nholes-1]=len-kvec+1;
					holes->nztot+=holes->nzeros[holes->nholes-1];
					holes->kzeros[holes->nholes-1]=kvec;

				}
				else{
					if(hole>0){		
						for(i=kvec;i<len;i++)vec[i]=0;
						len2=len;
						len1=0;
						holes->nholes++;
						holes->nzeros[holes->nholes-1]=len-kvec+1;
						holes->nztot+=holes->nzeros[holes->nholes-1];
						holes->kzeros[holes->nholes-1]=kvec;
						kvec=len2;
					}
				}
			}
		}

		if(len1<=0){
			if(allpers->nper>0){
				// TO BE DONE WITH FUNCTION ZERO_PERS
			}
		}
	}

	return;
}


void vec_from_sds1(float *buffer, float *vec, double *tim0, struct HOLES_* holes,
			struct SFC_* sds_, int chn, long len, struct ALLPERS_ *allpers){
/**pia version, 23 Aug. 2005: to call sfc_open1 and sds_open1:
so does not allocate SFC_ again !!**/

/* buffer should be allocated for len*nch 

    chn is the serial number of the channel (1,2,...)

    allpers->nper should be put equal to 0, if allowed periods are not used

*/

	int kskipped;
	long pos,ndat,len1,len2,kvec,count,nread,i,hole;
	double timend_exp;
	char *file;

	sds_->acc++;
	sds_->eof=0;
	holes->nztot=0;
	holes->nholes=0;
	kskipped=0;

	if(sds_->fid==NULL){
		printf(" *** End of files or error ! \n");
		return;
	}

	pos=ftell(sds_->fid);
	ndat=(pos-sds_->point0)/(sds_->nch*4);
	*tim0=sds_->t0+sds_->dt*ndat/86400;

	timend_exp=*tim0+len*sds_->dt/86400;

	len2=0;
	len1=len-len2;
	kvec=0;

	while(len1>0){
	    count=fread(buffer,4,sds_->nch*len1,sds_->fid);
	 	nread=count/sds_->nch;
		len2+=nread;
		len1=len-len2;
	
		if(nread>0){
			for(i=kvec;i<len2;i++)vec[i]=buffer[(i-kvec)*sds_->nch+chn-1];				     
	         	kvec=len2;		      
		}

		if(len1>0){
			sds_->eof=1;
			fclose(sds_->fid);
			if(strcmp(sds_->filspost,"#NOFILE")==0){
				sds_->eof=2;
				printf(" ---------> End of concatenated files \n");
				for(i=kvec;i<len;i++)vec[i]=0;
				len2=len;
				len1=0;
				holes->nholes++;
				holes->nzeros[holes->nholes-1]=len-kvec+1;
				holes->nztot+=holes->nzeros[holes->nholes-1];
				holes->kzeros[holes->nholes-1]=kvec;

				return;
			}
			else{
				file=strcat(sds_->pname,sds_->filspost);
				sds_open1(file,allpers,sds_); //pia: not sds_open

				if(sds_->eof==3){
					printf(" ---------> End of files \n");
					for(i=kvec;i<len;i++)vec[i]=0;
					len2=len;
					len1=0;
					holes->nholes++;
					holes->nzeros[holes->nholes-1]=len-kvec+1;
					holes->nztot+=holes->nzeros[holes->nholes-1];
					holes->kzeros[holes->nholes-1]=kvec;

					return;
				}
				printf(" *** open file %s \n",sds_->filme);
				puts(" *** A new file opened. Parameters:");
				sfc_show(sds_);
				hole=floor((sds_->t0-*tim0)*86400/sds_->dt-len2+0.5);
				if(hole>len1){
					for(i=kvec;i<len;i++)vec[i]=0;
					len2=len;
					len1=0;
					holes->nholes++;
					holes->nzeros[holes->nholes-1]=len-kvec+1;
					holes->nztot+=holes->nzeros[holes->nholes-1];
					holes->kzeros[holes->nholes-1]=kvec;

				}
				else{
					if(hole>0){			
						for(i=kvec;i<kvec+hole;i++)vec[i]=0;
						len2=len2+hole;
						len1=len-len2;
						holes->nholes++;
						holes->nzeros[holes->nholes-1]=hole;
						holes->nztot+=holes->nzeros[holes->nholes-1];
						holes->kzeros[holes->nholes-1]=kvec;
						kvec=len2;
					}
				}
			}
		}

		if(len1<=0){
			if(allpers->nper>0){
				// TO BE DONE WITH FUNCTION ZERO_PERS
			}
		}
	}

	return;
}




void vec_from_sds1_0(float *buffer, float *vec, double *tim0, struct HOLES_* holes,
			struct SFC_* sds_, int chn, long len, struct ALLPERS_ *allpers){
/**pia version, 23 Aug. 2005: to call sfc_open1 and sds_open1:
so does not allocate SFC_ again !!  

   doesn't embodies short holes ! ---   **/ 

/* buffer should be allocated for len*nch 

    chn is the serial number of the channel (1,2,...)

    allpers->nper should be put equal to 0, if allowed periods are not used

*/

	int kskipped;
	long pos,ndat,len1,len2,kvec,count,nread,i,hole;
	double timend_exp;
	char *file;

	sds_->acc++;
	sds_->eof=0;
	holes->nztot=0;
	holes->nholes=0;
	kskipped=0;

	if(sds_->fid==NULL){
		printf(" *** End of files or error ! \n");
		return;
	}

	pos=ftell(sds_->fid);
	ndat=(pos-sds_->point0)/(sds_->nch*4);
	*tim0=sds_->t0+sds_->dt*ndat/86400;

	timend_exp=*tim0+len*sds_->dt/86400;

	len2=0;
	len1=len-len2;
	kvec=0;

	while(len1>0){
	    count=fread(buffer,4,sds_->nch*len1,sds_->fid);
	 	nread=count/sds_->nch;
		len2+=nread;
		len1=len-len2;
	
		if(nread>0){
			for(i=kvec;i<len2;i++)vec[i]=buffer[(i-kvec)*sds_->nch+chn-1];				     
	         	kvec=len2;		      
		}

		if(len1>0){
			sds_->eof=1;
			fclose(sds_->fid);
			if(strcmp(sds_->filspost,"#NOFILE")==0){
				sds_->eof=2;
				printf(" ---------> End of concatenated files \n");
				for(i=kvec;i<len;i++)vec[i]=0;
				len2=len;
				len1=0;
				holes->nholes++;
				holes->nzeros[holes->nholes-1]=len-kvec+1;
				holes->nztot+=holes->nzeros[holes->nholes-1];
				holes->kzeros[holes->nholes-1]=kvec;

				return;
			}
			else{
				file=strcat(sds_->pname,sds_->filspost);
				sds_open1(file,allpers,sds_); //pia: not sds_open

				if(sds_->eof==3){
					printf(" ---------> End of files \n");
					for(i=kvec;i<len;i++)vec[i]=0;
					len2=len;
					len1=0;
					holes->nholes++;
					holes->nzeros[holes->nholes-1]=len-kvec+1;
					holes->nztot+=holes->nzeros[holes->nholes-1];
					holes->kzeros[holes->nholes-1]=kvec;

					return;
				}
				printf(" *** open file %s \n",sds_->filme);
				puts(" *** A new file opened. Parameters:");
				sfc_show(sds_);
				hole=floor((sds_->t0-*tim0)*86400/sds_->dt-len2+0.5);
				if(hole>len1){
					for(i=kvec;i<len;i++)vec[i]=0;
					len2=len;
					len1=0;
					holes->nholes++;
					holes->nzeros[holes->nholes-1]=len-kvec+1;
					holes->nztot+=holes->nzeros[holes->nholes-1];
					holes->kzeros[holes->nholes-1]=kvec;

				}
				else{
					if(hole>0){		
						for(i=kvec;i<len;i++)vec[i]=0;
						len2=len;
						len1=0;
						holes->nholes++;
						holes->nzeros[holes->nholes-1]=len-kvec+1;
						holes->nztot+=holes->nzeros[holes->nholes-1];
						holes->kzeros[holes->nholes-1]=kvec;
						kvec=len2;
					}
				}
			}
		}

		if(len1<=0){
			if(allpers->nper>0){
				// TO BE DONE WITH FUNCTION ZERO_PERS
			}
		}
	}

	return;
}




/* ----------------------  sbl Routines -------------------------*/

/*
%
                         SBL file format
%
    *** SBL is one of the file formats of the SFC collection ***
%
  The SBL (Simple block data format) format is intended for storing one or more data
  sets by means of "blocks" composed of sub-blocks, in a variety of different cases.
  
  The format is the following :
%
   label        - 8 bytes string (1-8)        #SFC#SBL
%
  First 944 bytes : see the sfc_open function.
%
   ch()         - N*128 bytes structures       names and captions for the channels
       .dx      - double (1-8)                 sampling first dimension (if any, otherwise 0)
       .dy      - double (9-16)                sampling second dimension (if any, otherwise 0)
       .lenx    - int32  (17-20)               length first dimension (number of rows)
       .leny    - int32  (21-24)               length second dimension (for single dim arrays, 1)
       .type    - int32  (25-28)               type of data:
                                                  1   byte
                                                  2   int16
                                                  3   int32
                                                  4   float
                                                  5   float complex
                                                  6   double
                                                  7   double complex
                                                 ..   compressed formats   
       .name    - 100 bytes string (29-128)    name or caption
%
   user         - (free)                       user header
%
   data blocks  - ....                         the data (l for each channel; float)
%
  To access the user header, fseek(fid,944+N*128,'bof')
%
  The sbl_ structure derives from the sfc_ structure.
  The added or modified members are:
%
    sds_.ch(nch)      channel structures
           .dx
           .dy
           .lenx      number of rows
           .leny      number of columns
           .type      data type of the channel
           .name      channel name or caption
           .len       length of the sub-block (in bytes)
           .inix      initial value of the first abscissa
           .iniy      initial value of the second abscissa
           ....       other (depending by the data types)
           .bias      position of the first data of the block (in bytes)
           .k         number of read data (pointer to the next data)
    sds_.len          number of blocks
    sds_.blen         block length (in bytes)
    sds_.point(nch)   pointer for next data
%
%
                      Structure of the blocks
%
  All blocks have the same length and structure.
  Each block is composed by "channels" (or sub-blocks), containing a short
  header and data, in the following way:
%
       block number (a 16 byte string as "[BLNxxxxxxxxxxx]")
       block time   (double; may be not meaningful)
%
       ch(1).inix (always present; may be not meaningful)
       ch(1).iniy (always present; may be not meaningful)
       ch(1).par1 (data type specific; may be absent)
       ch(1).par2 (data type specific; may be absent)
       ...
       A(lenx,leny)   data of channel 1
%
  then the same for channel 2 and so on.
*/


SFC_* sbl_open(struct SFC_* sbl_){
/*
*/
  printf(" *** NOT YET IMPLEMENTED !"); //pia
  return sbl_;  //pia
}


SFC_* sbl_openw(struct SFC_* sbl_){
/*
*/
  printf(" *** NOT YET IMPLEMENTED !"); //pia
  return sbl_; //pia
}


/* ----------------------  vbl Routines -------------------------*/
/*
%
                         VBL file format
%
    *** VBL is one of the file formats of the SFC collection ***
%
  The vbl (Variable block data format) format is intended for storing one or more data
  sets by means of "blocks" composed of sub-blocks, in a variety of different cases.
  The base structure is similar to the SBL format, the unique difference is
  that every sub-block array is preceded by a label as  [CHxxxx] and two
  int32 with the two dimensions of the array (plus the inix and iniy as in
  SBL). At the end an "index" with the pointers to the blocks may be
  present.
  
  The format is the following :
%
   label        - 8 bytes string (1-8)        #SFC#VBL
%
  First 944 bytes : see the sfc_open function. Then the ch description,
  identical to the sbl format, except the lenx and leny parameters that are
  variable block by block.
%
   ch()         - N*128 bytes structures       names and captions for the channels
       .dx      - double (1-8)                 sampling first dimension (if any, otherwise 0)
       .dy      - double (9-16)                sampling second dimension (if any, otherwise 0)
       .lenx    - int32  (17-20)               0
       .leny    - int32  (21-24)               0
       .type    - int32  (25-28)               type of data:
                                                  1   byte
                                                  2   int16
                                                  3   int32
                                                  4   float
                                                  5   float complex
                                                  6   double
                                                  7   double complex
                                                 ..   compressed formats   
       .name    - 100 bytes string (29-128)    name or caption
%
   user         - (free)                       user header
%
   data blocks  - ....                         the data (l for each channel; float)
%
  To access the user header, fseek(fid,944+N*128,'bof')
%
  The vbl_ structure derives from the sfc_ structure. 
  The added or modified members are:
%
    vbl_.ch(nch)      channel structures
           .dx
           .dy
           .lenx      number of rows
           .leny      number of columns
           .type      data type of the channel
           .name      channel name or caption
           .len       length of the sub-block (in bytes)
           .inix      initial value of the first abscissa
           .iniy      initial value of the second abscissa
           ....       other (depending by the data types)
           .bias      position of the first data of the block (in bytes)
           .k         number of read data (pointer to the next data)
    vbl_.len          number of blocks
    vbl_.point(nch)   pointer for next data
%
%
                      Structure of the blocks
%
  All blocks have the same length and structure.
  Each block is composed by "channels" (or sub-blocks), containing a short
  header and data, in the following way:
%
       block number (a 16 byte string as "[BLNxxxxxxxxxxx]")
       block time   (double; may be not meaningful)
%
       ch number  (a 8 byte string as "[CHxxxx]")
       ch(1).lenx 
       ch(1).leny
       ch(1).inix (always present; may be not meaningful)
       ch(1).iniy (always present; may be not meaningful)
       ch(1).par1 (data type specific; may be absent)
       ch(1).par2 (data type specific; may be absent)
       ...
       A(lenx,leny)   data of channel 1
%
  then the same for channel 2 and so on.
%
  At the end, after all blocks, may be present a block index, with the following
  structure:
%
      pointer to bl 1           8 byte integer
      pointer to bl 2           8 byte integer
          ...                       ...
      total number of blocks    8 byte integer
      [-INDEX]                  8 byte string
%
  The length of the index is (Nbl+2)*8 bytes.
*/


SFC_* vbl_open(struct SFC_* vbl_){
/*
*/
  printf(" *** NOT YET IMPLEMENTED !"); //pia
  return vbl_; //pia
}


SFC_* vbl_openw(struct SFC_* vbl_){
/*
*/
  printf(" *** NOT YET IMPLEMENTED !"); //pia
  return vbl_; //pia
}




/* ----------------------  Service Routines -------------------------*/

void sfc_util(){
	char *title="SFC_Util menu";
	char *items[12]={"Directory choice",
					"File choice",
					"Channel choice",
					"DataType file name block",
					"File Resume",
					" *** prova",
					"Time check",
					"Batch mode (0,1)",
					"...",
					"Dump file standard (on stdumpfile.txt)",
					"Dump file (on dumpfile.txt)",
					"Exit"};
	char filename[200];
	char direc[80];
	char file[80];
	char channel[40];
	char datatype[5];
	char sdsfile[80];
	//char* dumpfile="dumpfile.txt"; //pia: tolta perche' non usata
	//char* stdumpfile="stdumpfile.txt"; //pia: tolta perche' non usata
	char strou[80],inmenu[200];
	int itout,batmode=0;

	SFC_* sfc_;

	sfc_=malloc(sizeof(SFC_));

	strcpy(sfc_->pname,"D:\\Data\\Virgo\\E4\\sds\\");
	strcpy(sfc_->filme,"VIR_20020712_135947_hrec_000.sds");
	//strcpy(channel,"Pr_B1_ACq");
	//strcpy(channel,"Pr_B1_ACq_4kHz");
	//strcpy(channel,"dL_4kHz_NoLines");
	//strcpy(channel,"dL_20kHz");
	//strcpy(channel,"dL_20kHz_NoLines");
	strcpy(channel,"dL_4kHz");
	strcpy(datatype,"xxxx");
	strcpy(sdsfile,"prova.sds");

	while (1)
	{
		strcpy(filename,direc);
		strcat(filename,file);printf(" |%s| |%s| |%s|\n",direc,file,filename);
		printf("Working on %s  -> channel %s\n",filename,channel);
		printf("Output file :  %s \n",sdsfile);
		//pia:tolta da' 2 errori:
		//implicit declaration of function `GPS2UTC'
		//e passing arg 1 of `asctime' makes pointer from integer without a cast
		//printf(" ---> GPS time = 0 is %s \n",asctime(GPS2UTC(0,19)));
		itout=getmenu(title,items,12);

		printf("\n   %s chosen\n",items[itout-1]);
		

		switch(itout){
		case 1:
			strcpy(inmenu,"Directory ?  - now ");
			strcat(inmenu,direc);
			strcat(inmenu," >");
			if(batmode==0)
				get_line(inmenu,strou);
			else
				get_line_bat(inmenu,strou);
			strcpy(direc,strou);
			break;
		case 2:
			strcpy(inmenu,"File ?  - now ");
			strcat(inmenu,file);
			if(batmode==0)
				get_line(inmenu,strou);
			else
				get_line_bat(inmenu,strou);
			strcpy(file,strou);
			break;
		case 3:
			strcpy(inmenu,"Channel ?  - now ");
			strcat(inmenu,channel);
			if(batmode==0)
				get_line(inmenu,strou);
			else
				get_line_bat(inmenu,strou);
			strcpy(channel,strou);
			break;
		case 4:
			strcpy(inmenu,"Output file name DataType ?  - now ");
			strcat(inmenu,datatype);
			if(batmode==0)
				get_line(inmenu,strou);
			else
				get_line_bat(inmenu,strou);
			strcpy(datatype,strou);
			break;
		case 5:
			break;
		case 6:
			sfc_=sfc_open(file);
			printf(" label %s \n",sfc_->label);
		case 7:
		        //piapia tolto, come dice Sergio 16 aug 2005 inquire_fr_file(filename);
			break;
		case 8:
			printf("Now batch mode is %d ; new mode ? (0,1) ",batmode);
			scanf("%d",&batmode);
		case 9:
		        //piapia tolto, come dice Sergio 16 aug 2005  inquire_ch(filename);
			break;
		case 10:
			//print_dump_std(stdumpfile, filename);
			break;
		case 11:
			//print_dump(dumpfile, filename);
			break;
		case 12:
			return;
			break;
//		default:
//			return;
		}
	}
	return;
}





void gw_filename(struct tm *t, char *ant_abbr, char *dat_abbr, char *sub_typ, char *ext, char *out){
/*
	t			time structure
	ant_abbr	antenna abbreviation
	dat_abbr	data type abbreviation
	num			sequence number
	ext			file extension

  The format for the data files is the following:

    AAA_DDDD_dddddddd_tttttt_nnnn.xxx

  (32 bytes), where:
	- AAA			antenna abbreviation (e.g. VIR)
	- DDDD			data type (e.g. hrec)
	- dddddddd		date (e.g. 20030612)
	- tttttt		time (e.g. 120420)
	- nnnn			subtype; for example sequence number (for the same AAA_DDDD_dddddddd_tttttt; e.g. 000)
	- xxx			extension depending on the file format (e.g. sds). Even if the 
					suggested length of the extension is three characters, it can be
					any value (in this case the length of the name is not 32)
*/
	char dat[9],tim[7];

	strcpy(out,ant_abbr);
	strcat(out,"_");
	strcat(out,dat_abbr);
	strcat(out,"_");
	sprintf(dat,"%04d%02d%02d",t->tm_year+1900,t->tm_mon+1,t->tm_mday);
	strcat(out,dat);
	strcat(out,"_");
	sprintf(tim,"%02d%02d%02d",t->tm_hour,t->tm_min,t->tm_sec);
	strcat(out,tim);
	strcat(out,"_");
	strcat(out,sub_typ);
	strcat(out,".");
	strcat(out,ext);

	printf("\n  |%s| \n",out);
	return;
}
