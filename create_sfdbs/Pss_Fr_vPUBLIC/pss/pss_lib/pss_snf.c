/*___________________________________________________
 ¦                                                   ¦
 ¦                    pss_snf.c                      ¦
 ¦       by Sergio Frasca - ((( 0 ))) Virgo          ¦
 ¦                    March 2000                     ¦
 ¦___________________________________________________¦*/


#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "pss_snag.h"
#include "pss_math.h"
#include "pss_serv.h"
#include "pss_snf.h"

#define SNF_PROTOCOL "snf02"

struct SNF_NIDG{
	char	C[6];
	short	I[4];
	long	L[4];
	float	F[6];
	double 	D[6];
} NIDG0={{'#','N','I','D','G','#'},{1,2000,-3000,12345},{1,2000000,-3000000,123456789},
{1,2.e10,-3.e-12,4.e-14,5.e16,123456.},
{1,2.e30,-3.e-40,4.e-50,5.e60,1234567890.}};


/* --------------------------------------------- */


struct LogX_STRUCT* write_LogX(struct LogX_STRUCT* str, short coding, FILE *file)
/* coding = 0  -> use input parameters 
          > 0  optimization level 
		  > 0  optimization on base, min and sign
		  = 2  optimization on dimen 
		  = 3  optimization on log-lin
		  = 4  full optimization */
{
	str=code_LogX(str,coding);

	fwrite(&str->stat,2,1,file);
	fwrite(&str->len,4,1,file);

	if(str->linlog<2){
		fwrite(&str->m,8,1,file);
		fwrite(&str->b,8,1,file);
		fwrite(str->coded,1,str->nchar,file);
	}
	else if(str->linlog==2){
		fwrite(&str->m,4,1,file);
		fwrite(str->coded,4,1,file);
	}
	else{
		fwrite(str->coded,1,str->nchar,file);
	}

	return str;
}


struct LogX_STRUCT* read_LogX(struct LogX_STRUCT* str, FILE *file){
	fread(&str->stat,2,1,file);
	decode_LogX_stat(str);
	fread(&str->len,4,1,file);
	str->nchar=nchar_LogX(str);

	if(str->linlog<2){
		fread(&str->m,8,1,file);
		fread(&str->b,8,1,file);
		fread(str->coded,1,str->nchar,file);
	}
	else if(str->linlog==2){
		fread(&str->m,4,1,file);
		fread(str->coded,4,1,file);
	}
	else{
		fread(str->coded,1,str->nchar,file);
	}
	
	return str;
}


struct LogX_STRUCT* code_LogX(struct LogX_STRUCT* str, short coding)
/* coding = 0  -> use input parameters 
          > 0  optimization level 
		  > 0  optimization on base, min and sign
		  = 2  optimization on dimen 
		  = 3  optimization on log-lin
		  = 4  full optimization */
{
	long i,ii,j;
	long imax,iapp;
	short lenexp,kbuff,maxkbuff;
	float xmin,xmax,amin,amax,buffer[8],epsval,semi;
	float *outf;
	short *outs;
	unsigned short *uouts;
	unsigned char *uoutc,uch;
	char have0,ok;
	double logb,logm;

	str->errcode=0;
	xmin=xmax=str->vec[0];
	amin=1.e38;
	amax=0; printf(" \n");

	ok=1;
	if((str->dimen<8)&&(!str->linlog))ok=0;
	if((!ok) && (coding>0))printf("\n *** Attention ! Erroneous coding ! \n"); 

	if((coding>0)&&ok){            /* computes some coding parameters */
		have0=0;
		for(i=0;i<str->len;i++){
			if(xmin>str->vec[i])xmin=str->vec[i]; printf(" >%f ",str->vec[i]);
			if(xmax<str->vec[i])xmax=str->vec[i];
			if(str->vec[i]!=0){
				if(amin>fabs(str->vec[i]))amin=fabs(str->vec[i]);
				if(amax<fabs(str->vec[i]))amax=fabs(str->vec[i]);
			}
			else have0=1;
		}
		epsval=str->epsval;
		if(str->epsval==0)epsval=amin;

		if(str->satur>0){
			if(xmin>str->satur)xmin=str->satur*xmin/fabs(xmin);
			if(xmax>str->satur)xmax=str->satur*xmax/fabs(xmax);
			if(amin>str->satur)amin=str->satur;
			if(amax>str->satur)amax=str->satur;
		}

//printf("\n @@@ xmax,xmin = %f %f ",xmax,xmin);
		if((coding==2)||(coding==4))str->dimen=16; /* UNDER DEVELOPMENT */

		if((coding>2)&&(str->dimen>7)){            /* linear if small max/min ratio */
			if(amax/amin<=20)str->linlog=1;
			else str->linlog=0;
		}

		if(xmin*xmax<0){           /* mixed signs */
			str->sign=0;

			lenexp=str->dimen-1;
		}
		else{                      /* fixed sign */
			str->sign=-1;
			if(xmax>0)str->sign=1;

			lenexp=str->dimen;

			if(xmin==xmax){
				str->linlog=2;
				*(float *)str->coded=xmin;
				str->nchar=4;

				str->stat=code_LogX_stat(str);

				printf("\n xmin=xmax= %f",xmin); // debug

				return str;
			}
		}
		
/*	choice of m and b  */
		if(str->linlog){         /* linear */
			if(str->dimen<7||str->sign){   /* the code must be unsigned */
				str->m=xmin;
				str->b=(xmax-xmin)/pow(2,str->dimen); //printf("\n @@@ xmax,xmin,b = %f %f %g ",xmax,xmin,str->b);
			}
			else{                          /* the code must be signed */
				str->m=0;
				str->b=amax/(pow(2,str->dimen-1)-0.5);
			}
		}
		else{                    /* logarithmic */
			if(str->dimen>7)
				if(str->sign){
					str->m=amin;
					str->b=pow((amax/amin),1./(pow(2,lenexp)-0.5)); printf("\n @@@ lenexp = %d",lenexp);
				}
				else{
					str->m=epsval;
					str->b=pow((amax/epsval),1./(pow(2,lenexp)-1.5)); printf("\n @@@ lenexp A = %d",lenexp);
				}
			else{
				str->m=0;
				str->b=0.;
			}
		}
	}

	show_LogX_coding(str);

/*	starts coding  */

	maxkbuff=0;
	if(str->dimen<9)maxkbuff=8/str->dimen-1;

	if(str->dimen==32){
		outf=(float *)str->coded; /* gives the same pointer of the character array str->coded */
		for(i=0;i<str->len;i++)outf[i]=str->vec[i];
		str->nchar=4*str->len;

		str->stat=code_LogX_stat(str);

		return str;
	}

	if(str->dimen==16){
		outs=(short *)str->coded;
		uouts=(short *)str->coded;

		if(str->linlog)
			if(str->sign){
				imax=pow(2,16)-1;

				for(i=0;i<str->len;i++){
					iapp=(str->vec[i]-str->m)/str->b;
					if(iapp>imax)iapp=imax;
					uouts[i]=iapp;
				}
			}
			else{
				imax=pow(2,15)-1;
				semi=str->b/2;

				for(i=0;i<str->len;i++){
					iapp=(str->vec[i]-str->m+semi)/str->b;
					if(iapp>imax)iapp=imax;
					if(iapp<-imax)iapp=-imax;
					outs[i]=iapp;
				}
			}
		else{
			logb=log(str->b);
			logm=log(str->m);
			
			if(str->sign)
				for(i=0;i<str->len;i++)
					uouts[i]=(log(fabs(str->vec[i]))-logm)/logb;
			else
				for(i=0;i<str->len;i++){
					if(str->vec[i]==0.)
						uouts[i]=0xffff;
					else{
						uouts[i]=(log(fabs(str->vec[i]))-logm)/logb; //show_bit_s(&uouts[i], "A");
						if(str->vec[i]<0)uouts[i]=uouts[i]|0x8000;
					}
				}
		}
		
		str->nchar=2*str->len;

		str->stat=code_LogX_stat(str);

		return str;
	}

	if(str->dimen==8){
		if(str->linlog){
			uoutc=str->coded;

			if(str->sign){
				imax=pow(2,8)-1;

				for(i=0;i<str->len;i++){
					iapp=(str->vec[i]-str->m)/str->b;
					if(iapp>imax)iapp=imax;
					uoutc[i]=iapp;
				}
			}
			else{
				imax=pow(2,7)-1;
				semi=str->b/2;

				for(i=0;i<str->len;i++){
					iapp=(str->vec[i]-str->m+semi)/str->b;
					if(iapp>imax)iapp=imax;
					str->coded[i]=iapp;
				}
			}
		}
		else{
			logb=log(str->b);
			logm=log(str->m);

			if(str->sign)
				for(i=0;i<str->len;i++)
					str->coded[i]=log(fabs(str->vec[i]/str->m))/logb;
			else{
				for(i=0;i<str->len;i++){
					if(str->vec[i]==0.)
						str->coded[i]=0xff;
					else{
						str->coded[i]=log(fabs(str->vec[i]/str->m))/logb;
						if(str->vec[i]<0)str->coded[i]=str->coded[i]|0x80;
					}
				}
			}
		}
			
		str->nchar=str->len;

		str->stat=code_LogX_stat(str);

		return str;
	}

	/* coding for sub-byte formats */

	ii=0;
	str->coded[0]=0;
	str->nchar=(str->len-1)/(8/str->dimen)+1; //printf("\n  @@@ nchar = %d ",str->nchar);
	uoutc=str->coded;
	imax=pow(2,str->dimen)-1; //printf("\n @@@ imax = %d",imax);
	if(str->linlog==0)str->linlog=3;

	if(str->linlog){
		kbuff=0;
		for(i=0;i<str->len;i++){
			buffer[kbuff]=str->vec[i];
			if(str->satur>0){
				if(buffer[kbuff]>str->satur)buffer[kbuff]=str->satur;
			}
			kbuff++;
			if(kbuff>maxkbuff){
				kbuff=0;
				for(j=0;j<=maxkbuff;j++){
					iapp=(buffer[j]-str->m)/str->b;
					if(iapp>imax)iapp=imax;
					uch=iapp;
					write_bit_c(str->coded+ii,uch,j*str->dimen,str->dimen); printf("\n uch %d coded %d",uch,uoutc[ii]);
				}
				ii++;
				if(ii<=str->nchar)str->coded[ii]=0; show_bit_c(&str->coded[ii-1],"coded");
			}
		}
		if(ii<str->nchar)
			for(j=0;j<kbuff;j++){
				iapp=(buffer[j]-str->m)/str->b;
					if(iapp>imax)iapp=imax;
					uch=iapp;
					write_bit_c(str->coded+ii,uch,j*str->dimen,str->dimen); printf("\n uch %d coded %d",uch,str->coded[ii]);
			}
		show_bit_c(&str->coded[ii],"coded");
	}
	else{
	}

	str->stat=code_LogX_stat(str);

	return str;
}


struct LogX_STRUCT* decode_LogX(struct LogX_STRUCT* str){
	long i,ii,j;
	int maxkbuff;
	float *outf,semi;
	short *outs,outs15,outs1;
	unsigned short *uouts;
	char outc7,outc1;
	unsigned char *uoutc;

	if(str->linlog==2){
		outf=(float *)str->coded;
		for(i=0;i<str->len;i++)str->vec[i]=outf[0];
		return str;
	}

	if(str->dimen==32){        /* 32 bits */
		outf=(float *)str->coded;
		for(i=0;i<str->len;i++)str->vec[i]=outf[i];
	}
	else if(str->dimen==16){   /* 16 bits */
		if(str->linlog){
			if(str->sign){
				uouts=(short *)str->coded;
				semi=str->b/2;
				for(i=0;i<str->len;i++)
					str->vec[i]=str->m+str->b*uouts[i]+semi;
			}
			else{
				outs=(short *)str->coded;
				for(i=0;i<str->len;i++)
					str->vec[i]=str->m+str->b*outs[i];
			}
		}
		else{
			if(str->sign){
				uouts=(unsigned short *)str->coded;
				for(i=0;i<str->len;i++)
					str->vec[i]=str->sign*str->m*pow(str->b,uouts[i]);
			}
			else{
				outs=(short *)str->coded;
				for(i=0;i<str->len;i++){
					if(outs[i]==0xffff)str->vec[i]=0;
					else{
						outs15=read_bit_s(outs[i],0,15);
						outs1=read_bit_s(outs[i],15,1);
						str->vec[i]=str->m*pow(str->b,outs15);
						if(outs1)str->vec[i]=-str->vec[i];
					}
				}
			}
		}
	}
	else if(str->dimen==8){   /* 8 bits */
		if(str->linlog==1){
			if(str->sign){
				uoutc=str->coded;
				semi=str->b/2;
				for(i=0;i<str->len;i++)
					str->vec[i]=str->m+str->b*uoutc[i]+semi;
			}
			else{
				for(i=0;i<str->len;i++)
					str->vec[i]=str->m+str->b*str->coded[i];
			}
		}
		else{
			if(str->sign){
				uoutc=(unsigned char *)str->coded;
				for(i=0;i<str->len;i++)
					str->vec[i]=str->sign*str->m*pow(str->b,uoutc[i]);
			}
			else{
				for(i=0;i<str->len;i++){
					if(str->coded[i]==0xff)str->vec[i]=0;
					else{
						outc7=read_bit_s(str->coded[i],0,7);
						outc1=read_bit_s(str->coded[i],7,1);
						str->vec[i]=str->m*pow(str->b,outc7);
						if(outc1)str->vec[i]=-str->vec[i];
					}
				}
			}
		}
	}
	else{                   /* < 8 bits */
		maxkbuff=8/str->dimen;
		ii=0;
		str->nchar=(str->len-1)/(8/str->dimen)+1;
		semi=str->b/2;

		if(str->linlog==1){
			for(i=0;i<str->nchar;i++){
				for(j=0;j<maxkbuff;j++){
					str->vec[ii]=str->m+str->b*read_bit_c(str->coded[i],j*str->dimen,str->dimen)+semi;
					ii++;
				}
			}
		}
		else{
		}
	}

	return str;
}


short code_LogX_stat(struct LogX_STRUCT* str)
/* Creates the stat variable with the coding information */
{
	short stat=0,expX;

	if(str->sign==1)write_bit_s(&stat,1,0,1);
	if(str->sign==0)write_bit_s(&stat,1,1,1);

	if(str->linlog==2)write_bit_s(&stat,1,2,1);

	expX=iLog2(str->dimen);
	if(expX>5)expX=5;

	write_bit_s(&stat,expX,3,3);// printf("\n expX %d",expX);show_bit_s(&stat,"@@@");
	if(str->linlog==1)write_bit_s(&stat,1,6,1);

	return stat;
}


short decode_LogX_stat(struct LogX_STRUCT* str)
/* Reads the LogX stat variable and updates the LogX structure */
{
	short stat=0,app;

	stat=str->stat;

	str->sign=-1;
	app=read_bit_s(stat,0,1);
	if(app)str->sign=1;

	app=read_bit_s(stat,1,1);
	if(app)str->sign=0;

	app=read_bit_s(stat,6,1);
	str->linlog=app;

	app=read_bit_s(stat,2,1);
	if(app)str->linlog=2;
	
	app=read_bit_s(stat,3,3); //printf("\n @@@ s3 = %d",app);
	str->dimen=pow(2,app); // show_bit_s(&app,"@@@");

	if((str->dimen<8)&(str->linlog==0))str->linlog=3;

	return 1;
}


long nchar_LogX(struct LogX_STRUCT* str){
	long nchar;

	nchar=(str->len*str->dimen)/8;

	if(str->dimen<8)
		if((nchar*8)<(str->len*str->dimen))nchar=nchar+1;

	return nchar;
}


int show_LogX_coding(struct LogX_STRUCT* str)
/* from LogX structure, not from stat */
{
	char *linlog[3]={"logarithmic","linear","constant"};
	char *sign[3]={"negative","mixed","positive"};

	printf("\n  --- LogX coding : %s  %s",linlog[str->linlog],sign[str->sign+1]);
	printf("\n  ---  m , b  :  %g , %.10g",str->m,str->b);
	printf("\n  ---  data dimension (bits) :  %d ",str->dimen);
	printf("\n  ---  number of data, byte number in coded data : %d , %d",str->len,str->nchar);
	show_bit_s(&str->stat,"stat:");

	return 1;
}



struct SpVec_STRUCT* write_SpVec(struct SpVec_STRUCT* str, short coding,
								 struct LogX_STRUCT* lxstr, short lxcoding,FILE *file)
/* coding = 0  -> use input parameters 
          > 0  -> optimize parameters

  lxstr, lxcoding = LogX structure and coding for the non-zero elements
*/
{
	str=code_SpVec(str,coding);

	if(str->binsp>1){       /* sparse */
		fwrite(&str->lenvec,4,1,file);
		fwrite(&str->lenspar,4,1,file);
		fwrite(str->coded,1,str->nchar,file);

		if(str->binsp==2){  /* sparse binary */
			if(str->logx){
			}
			else{
				fwrite(&str->n0run,4,1,file);
				fwrite(str->vec,4,str->n0run,file);
			}
		}
		else{              /* sparse non binary */
			fwrite(&str->lenvec,4,1,file);
			fwrite(str->vec,4,str->lenvec,file);
		}
	}
	else{                  /* non sparse */
		if(str->binsp==0){ /* non sparse non binary */
			fwrite(&str->lenvec,4,1,file);
			fwrite(str->vec,4,str->lenvec,file);
		}
		else {             /* non sparse binary */
			fwrite(&str->lenvec,4,1,file);
			fwrite(str->coded,1,str->nchar,file);
		}
	}

	return str;
}


struct SpVec_STRUCT* read_SpVec(struct SpVec_STRUCT* str, FILE *file){
	fread(&str->stat,2,1,file);
	decode_SpVec_stat(str);

	if(str->binsp>1){       /* sparse */
		fread(&str->lenvec,4,1,file);
		fread(&str->lenspar,4,1,file);
		str->nchar=(((str->lenspar+1)/2)*2*str->dimen)/8;
		fread(str->coded,1,str->nchar,file);

		if(str->binsp==2){  /* sparse binary */
			if(str->logx){
			}
			else{
				fread(&str->n0run,4,1,file);
				fread(str->vec,4,str->n0run,file);
			}
		}
		else{              /* sparse non binary */
			fread(&str->lenvec,4,1,file);
			fread(str->vec,4,str->lenvec,file);
		}
	}
	else{                  /* non sparse */
		if(str->binsp==0){ /* non sparse non binary */
			if(str->logx){
			}
			else{
				fread(&str->lenvec,4,1,file);
				fread(str->vec,4,str->lenvec,file);
			}
		}
		else {             /* non sparse binary */
			fread(&str->lenvec,4,1,file);
			str->nchar=(str->lenvec-1)/8+1;
			fread(str->coded,1,str->nchar,file);
		}
	}

	return str;
}


struct SpVec_STRUCT* code_SpVec(struct SpVec_STRUCT* str, short coding)
/* coding = 0  -> use input parameters 
          > 0  -> optimize parameters
*/
{
	long run[4]={0,0,0,0};
	long drun[4]={0,0,0,0};
	long runlen=0,drunlen=0,i,j,k,jj;
	long nbitmax=0,jmax=0,dnbitmax=0,djmax=0;
	long nbit,dnbit,lenvar,nelem;
	float valnz=0;
	char binary=1,deriv=0,bit,bit1,der;

	float *vec;
	unsigned long *lrun,imax;
	unsigned short *srun;
	unsigned char *crun,crun1[2],crind,cbyte,coded;

	if(coding){
		for(i=0;i<str->lenvec;i++){
			if(str->vec[i]){
				for(j=0;j<4;j++){
					lenvar=4*pow(2,j);
					run[j]=run[j]+(runlen-1)/pow(2,lenvar)+1;
				}
				runlen=0;
				if(valnz)
					if(str->vec[i]!=valnz)binary=0;
			}
			else{
				runlen++;
			}
			if(i>0){
				der=str->vec[i]-str->vec[i-1];
				if(der){
				for(j=0;j<4;j++){
					lenvar=4*pow(2,j);
					drun[j]=drun[j]+(drunlen-1)/pow(2,lenvar)+1;
				}
					drunlen=0;
				}
				else{
					drunlen++;
				}
			}
		}

		for(j=0;j<4;j++){
			lenvar=4*pow(2,j);
			nbit=lenvar*run[j];
			dnbit=lenvar*drun[j];
			if(nbit>nbitmax){
				nbitmax=nbit;
				jmax=j;
			}
			if(dnbit>dnbitmax){
				dnbitmax=dnbit;
				djmax=j;
			}
		}

		if(nbitmax>dnbitmax){
			deriv=1;
			jmax=djmax;
			nbitmax=dnbitmax;
		}

		str->derive=deriv*binary;
		if(binary)
			if(nbitmax>str->lenvec)binary=2;
		str->binsp=binary;
		str->dimen=4*pow(2,jmax);
	}

	if(str->binsp<2){		/* non-sparse coding */
		if(str->binsp==0){	/* 000 non-sparse, normal */
			vec=(float *)str->coded;
			for(i=0;i<str->lenvec;i++)vec[i]=str->vec[i];

			str->nchar=str->lenvec*4;
		}
		else{				/* 001 non-sparse, binary */
			i=0;
			j=0;
			while(i<str->lenvec){
				coded=0;
				cbyte=1;
				k=0;
				while((k<8)&(i<str->lenvec)){
					if(str->vec[i])coded=coded|cbyte;
					else coded=coded&~cbyte;

					cbyte*=2;
					i++;
					k++;
				}
				str->coded[j]=coded;
				j++;
			}
			str->nchar=j;
		}

		return str;
	}

	/* sparse coding; creates the i4runs vector */
	
	runlen=0;

	if(str->derive){	/* 111 derive, sparse, binary */
		j=0;
		bit1=0;
		for(i=0;i<str->lenvec;i++){
			bit=0;
			if(str->vec[i])bit=1;
			der=1;
			if(bit==bit1)der=0; //if((i/10)*10==i)printf("\n %d :  ",i); printf(" %d-%d-%d  ",bit,bit1,der);
			bit1=bit; 
			if(der){
				str->i4runs[j]=runlen;
				j++;
				runlen=0;
				imax=i;
			}
			else{
				runlen++;
			}
		}
	}
	else{				/* 01x sparse */
		for(i=0;i<str->lenvec;i++){
			if(str->vec[i]){
				str->i4runs[j]=runlen;
				if(str->binsp==2)str->nzelem[j]=str->vec[i];	/* 010 */
				j++;
				runlen=0;
				imax=i;
			}
			else{
				runlen++;
			}
		}
	}

	str->i4runs[j]=str->lenvec-imax;
	j++;
	str->n0run=j;

	/* code 0-runs for sparse */

	if(str->dimen==32){
		lrun=(long *)str->coded;
		for(i=0;i<str->n0run;i++)lrun[i]=str->i4runs[i]+1; /* note +1 */
		j=str->n0run;
	}
	else if(str->dimen==16){
		srun=(short *)str->coded;
		j=0;
		for(i=0;i<str->n0run;i++){
			nelem=str->i4runs[i]/65534;
			for(k=0;k<nelem;k++){
				srun[j]=0;
				str->i4runs[i]-=65534;
				j++;
			}
			srun[j]=str->i4runs[i]+1;
			j++;
		}
	}
	else if(str->dimen==8){
		crun=str->coded;
		j=0;
		for(i=0;i<str->n0run;i++){
			nelem=str->i4runs[i]/254;
			for(k=0;k<nelem;k++){
				crun[j]=0;
				str->i4runs[i]-=254;
				j++;
			}
			crun[j]=str->i4runs[i]+1;
			j++;
		}
	}
	else if(str->dimen==4){  
		crun=str->coded;
		crind=0;
		j=0;
		jj=0;
		for(i=0;i<str->n0run;i++){
			nelem=str->i4runs[i]/14;
			for(k=0;k<nelem;k++){
				crun1[crind]=0;
				crind++;
				j++;
				str->i4runs[i]-=14;
				if(crind>=2){
					crind=0;
					crun[jj]=crun1[0]+crun1[1]*16; printf("\n < %d %d %d %d",j,crun1[0],crun1[1],crun[j]);
					jj++; 
				}
			}
			crun1[crind]=str->i4runs[i]+1;
			crind++;
			j++;
			if(crind>=2){
				crind=0;
				crun[jj]=crun1[0]+crun1[1]*16; printf("\n > %d %d %d %d %d",j,str->i4runs[i],crun1[0],crun1[1],crun[j]);
				jj++; 
			}
		}
		if(crind){
			crun[jj]=crun1[0];
			jj++;
		}
	}

	str->lenspar=j;
	str->nchar=(((str->lenspar+1)/2)*2*str->dimen)/8;

	return str;
}


struct SpVec_STRUCT* decode_SpVec(struct SpVec_STRUCT* str){
	unsigned long *lrun;
	unsigned short *srun;
	unsigned char *crun,crun1[2],cr;
	float *vec;

	long i,j,k,lenrun,krun;
	char onoff;

	j=0;

	if(str->binsp<2){
		if(str->binsp==0){		/* 000 non-sparse, normal */
			vec=(float *)str->coded; //printf("\n @@@ ci sono ");
			for(i=0;i<str->lenvec;i++)
				str->vec[i]=vec[i];
		}
		else{					/* 001 non-sparse, binary */
			for(i=0;i<str->nchar;i++){
				cr=str->coded[i]; // show_bit_c(&cr," O-> "); printf("\n");
				for(k=0;k<8;k++){
					str->vec[j]=(cr&0x1)*str->nzelem[0]; //printf(" %g ",str->vec[j]);
					j++;
					cr=cr>>1;
				}
			}
		}
	}
	else{
		if(str->binsp){
			if(str->derive){	/* 111 derive, sparse, binary */
				onoff=0;
				if(str->dimen==32){
					lrun=(long *)str->coded; //printf("\n @@@ 32");
					j=0;

					for(i=0;i<str->lenspar;i++){// if(((i/10)*10)==i)printf("\n --- %d : ",i); printf(" %d ",lrun[i]);
						for(k=0;k<lrun[i]-1;k++){
							str->vec[j]=onoff*str->nzelem[0];
							j++;
						}
						onoff=onoff-1;
						if(onoff)onoff=1;
						str->vec[j]=onoff*str->nzelem[0];
						j++;
					}
				}
				else if(str->dimen==16){
					srun=(short *)str->coded; //printf("\n @@@ 16");
					lenrun=0;
					j=0;

					for(i=0;i<str->lenspar;i++){ //if(((i/10)*10)==i)printf("\n --- %d : ",i); printf(" %d ",srun[i]);
						if(srun[i]){
							lenrun+=srun[i]-1;
							for(k=0;k<lenrun;k++){
								str->vec[j]=onoff*str->nzelem[0];
								j++;
							}
							lenrun=0;
							onoff=onoff-1;
							if(onoff)onoff=1;
							str->vec[j]=onoff*str->nzelem[0];
							j++;
						}
						else lenrun+=65534;
					}
				}
				else if(str->dimen==8){
					crun=str->coded;
					lenrun=0;
					j=0;

					for(i=0;i<str->lenspar;i++){
						if(crun[i]){
							lenrun+=crun[i]-1;
							for(k=0;k<lenrun;k++){
								str->vec[j]=onoff*str->nzelem[0];
								j++;
							}
							lenrun=0;
							onoff=onoff-1;
							if(onoff)onoff=1;
							str->vec[j]=onoff*str->nzelem[0];
							j++;
						}
						else lenrun+=254;
					}
				}
				else if(str->dimen==4){
					crun=str->coded;
					lenrun=0;
					j=0;

					for(i=0;i<str->nchar;i++){
						crun1[0]=crun[i]&0xf;
						crun1[1]=crun[i]>>4; //if((i/10)*10==i)printf("\n %d -> ",i);printf(" %d %d ",crun1[0],crun1[1]);

						if(crun1[0]){
							lenrun+=crun1[0]-1;
							for(k=0;k<lenrun;k++){
								str->vec[j]=onoff*str->nzelem[0];
								j++;
							}//printf(" [%d] ",lenrun);
							lenrun=0;
							onoff=onoff-1;
							if(onoff)onoff=1;
							str->vec[j]=onoff*str->nzelem[0];
							j++;
						}
						else lenrun+=14;

						if(crun1[1]){
							lenrun+=crun1[1]-1;
							for(k=0;k<lenrun;k++){
								str->vec[j]=onoff*str->nzelem[0];
								j++;
							}//printf(" {%d} ",lenrun);
							lenrun=0;
							onoff=onoff-1;
							if(onoff)onoff=1;
							str->vec[j]=onoff*str->nzelem[0];
							j++;
						}
						else lenrun+=14;
					} //printf("\n @@@ nchar,j = %d, %d ",str->nchar,j);
				}
			}
			else{				/* 011 sparse, binary */
				if(str->dimen==32){
					lrun=(long *)str->coded; //printf("\n @@@ 32");
					j=0;

					for(i=0;i<str->lenspar;i++){// if(((i/10)*10)==i)printf("\n --- %d : ",i); printf(" %d ",lrun[i]);
						for(k=0;k<lrun[i]-1;k++){
							str->vec[j]=0;
							j++;
						}
						str->vec[j]=str->nzelem[0];
						j++;
					}
				}
				else if(str->dimen==16){
					srun=(short *)str->coded; //printf("\n @@@ 16");
					lenrun=0;
					j=0;

					for(i=0;i<str->lenspar;i++){ //if(((i/10)*10)==i)printf("\n --- %d : ",i); printf(" %d ",srun[i]);
						if(srun[i]){
							lenrun+=srun[i]-1;
							for(k=0;k<lenrun;k++){
								str->vec[j]=0;
								j++;
							}
							lenrun=0;
							str->vec[j]=str->nzelem[0];
							j++;
						}
						else lenrun+=65534;
					}
				}
				else if(str->dimen==8){
					crun=str->coded;
					lenrun=0;
					j=0;

					for(i=0;i<str->lenspar;i++){
						if(crun[i]){
							lenrun+=crun[i]-1;
							for(k=0;k<lenrun;k++){
								str->vec[j]=0;
								j++;
							}
							lenrun=0;
							str->vec[j]=str->nzelem[0];
							j++;
						}
						else lenrun+=254;
					}
				}
				else if(str->dimen==4){
					crun=str->coded;
					lenrun=0;
					j=0;

					for(i=0;i<str->nchar;i++){
						crun1[0]=crun[i]&0xf;
						crun1[1]=crun[i]>>4; //if((i/10)*10==i)printf("\n %d -> ",i);printf(" %d %d ",crun1[0],crun1[1]);

						if(crun1[0]){
							lenrun+=crun1[0]-1;
							for(k=0;k<lenrun;k++){
								str->vec[j]=0;
								j++;
							}//printf(" [%d] ",lenrun);
							lenrun=0;
							str->vec[j]=str->nzelem[0];
							j++;
						}
						else lenrun+=14;

						if(crun1[1]){
							lenrun+=crun1[1]-1;
							for(k=0;k<lenrun;k++){
								str->vec[j]=0;
								j++;
							}//printf(" {%d} ",lenrun);
							lenrun=0;
							str->vec[j]=str->nzelem[0];
							j++;
						}
						else lenrun+=14;
					} //printf("\n @@@ nchar,j = %d, %d ",str->nchar,j);
				}
			}
		}
		else{					/* x10 sparse, normal */
			if(str->dimen==32){
				lrun=(long *)str->coded; //printf("\n @@@ 32");
				j=0;

				for(i=0;i<str->lenspar;i++){// if(((i/10)*10)==i)printf("\n --- %d : ",i); printf(" %d ",lrun[i]);
					for(k=0;k<lrun[i]-1;k++){
						str->vec[j]=0;
						j++;
					}
					str->vec[j]=str->nzelem[krun];
					j++;
					krun++;
				}
			}
			else if(str->dimen==16){
				srun=(short *)str->coded; //printf("\n @@@ 16");
				lenrun=0;
				j=0;

				for(i=0;i<str->lenspar;i++){ //if(((i/10)*10)==i)printf("\n --- %d : ",i); printf(" %d ",srun[i]);
					if(srun[i]){
						lenrun+=srun[i]-1;
						for(k=0;k<lenrun;k++){
							str->vec[j]=0;
							j++;
						}
						lenrun=0;
						str->vec[j]=str->nzelem[krun];
						j++;
						krun++;
					}
					else lenrun+=65534;
				}
			}
			else if(str->dimen==8){
				crun=str->coded;
				lenrun=0;
				j=0;

				for(i=0;i<str->lenspar;i++){
					if(crun[i]){
						lenrun+=crun[i]-1;
						for(k=0;k<lenrun;k++){
							str->vec[j]=0;
							j++;
						}
						lenrun=0;
						str->vec[j]=str->nzelem[krun];
						j++;
						krun++;
					}
					else lenrun+=254;
				}
			}
			else if(str->dimen==4){
				crun=str->coded;
				lenrun=0;
				j=0;

				for(i=0;i<str->nchar;i++){
					crun1[0]=crun[i]&0xf;
					crun1[1]=crun[i]>>4; //if((i/10)*10==i)printf("\n %d -> ",i);printf(" %d %d ",crun1[0],crun1[1]);

					if(crun1[0]){
						lenrun+=crun1[0]-1;
						for(k=0;k<lenrun;k++){
							str->vec[j]=0;
							j++;
						}//printf(" [%d] ",lenrun);
						lenrun=0;
						str->vec[j]=str->nzelem[krun];
						j++;
						krun++;
					}
					else lenrun+=14;

					if(crun1[1]){
						lenrun+=crun1[1]-1;
						for(k=0;k<lenrun;k++){
							str->vec[j]=0;
							j++;
						}//printf(" {%d} ",lenrun);
						lenrun=0;
						str->vec[j]=str->nzelem[krun];
						j++;
						krun++;
					}
					else lenrun+=14;
				} //printf("\n @@@ nchar,j = %d, %d ",str->nchar,j);
			}
		}
	}

	return str;
}


short code_SpVec_stat(struct SpVec_STRUCT* str)
/* Creates the stat variable with the coding information */
{
	short stat=0,dimen;

	write_bit_s(&stat,str->binsp,0,2);
	if(str->derive)write_bit_s(&stat,1,2,1);

	dimen=iLog2(str->dimen/4);
	write_bit_s(&stat,dimen,3,2);

	if(str->logx)write_bit_s(&stat,1,6,1);

	return stat;
}


short decode_SpVec_stat(struct SpVec_STRUCT* str)
/* Reads the stat variable and updates the LogX structure */
{
	short stat,app;
	float app1;

	stat=str->stat;

	app=read_bit_s(stat,0,2);
	str->binsp=app;

	app=read_bit_s(stat,2,1);
	str->derive=app;

	app1=read_bit_s(stat,3,2);
	str->dimen=4*pow(2.,app1);

	str->logx=read_bit_s(stat,6,1);

	return 1;
}


int show_SpVec_coding(struct SpVec_STRUCT* str)
/* shows coding parameters */
{
	printf("\n  ---  SpVec coding : dimen  %d   binsp %d   derive  %d",
		str->dimen,str->binsp,str->derive);
	printf("\n  ---   lenvec: %d      nchar:   %d",str->lenvec,str->nchar);
	printf("\n  ---   n0run:  %d      lenspar: %d",str->n0run,str->lenspar);
	//show_bit_s(&str->stat,"stat:");

	return 1;
}



/* --------------------------------------------- */



SNF_STRUCT *snf_crea_struct(char* type, long mult){
	SNF_STRUCT *snf;
	long i;

	if(mult<1)mult=1;

	snf=(SNF_STRUCT *)malloc(sizeof(SNF_STRUCT));
	strcpy(snf->binasc,"bin");
	strcpy(snf->userprot," ");
	snf->nSNF2=1;
	snf->nSNF3=mult;
	snf->nSNF5=0;
	snf->nSNF6=0;
	snf->datlen=0;
	snf->nfields=0;
	snf->reclen=0;
	strcpy(snf->protocol,SNF_PROTOCOL);
	strcpy(snf->directaccess,"directaccess");
	snf->caption[0]=(char *)calloc(122,sizeof(char));
	strcpy(snf->caption[0],"no caption");
	snf->fileName=(char *)calloc(32,sizeof(char));
	strcpy(snf->fileName," ");
	sprintf(snf->fileName,"file_%s%03d.snf",type,mult);

	snf->access=(SNF_ACCESS *)malloc(sizeof(SNF_ACCESS));
	snf->access->verb=1;
	snf->access->inirec=0;
	snf->access->nrec=0;
	snf->access->trec=0.;

	snf->headv=(SNF_HEADV *)malloc(sizeof(SNF_HEADV));
	snf->headv->lheadi=0;
	snf->headv->lheadd=0;
	snf->headv->lheads=0;
	snf->headv->lrheadi=0;
	snf->headv->lrheadd=0;
	snf->headv->lrheads=0;
	snf->headv->lfheadi=0;
	snf->headv->lfheadd=0;
	snf->headv->lfheads=0;

	if(!strcmp(type,"gd")){
		strcpy(snf->obj,"gd");
		//strcpy(snf->fileName,"fil_gd.snf");
		snf->gd=(SNF_GD *)calloc(snf->nSNF3,sizeof(SNF_GD));
		for(i=0;i<snf->nSNF3;i++){
			sprintf(snf->gd[i].name,"gd%03d",i+1);
			snf->gd[i].n=0;
			snf->gd[i].type=1;
			snf->gd[i].complex=0;
			snf->gd[i].ini=0.;
			snf->gd[i].dx=1.;
			//snf->gd[i].nform=malloc(8*sizeof(char));
			//snf->gd[i].sform=malloc(10*sizeof(char));
			strcpy(snf->gd[i].nform,"float");
			strcpy(snf->gd[i].sform,"nosparse");
		}
	}
	else if(!strcmp(type,"dm")){
		strcpy(snf->obj,"dm");
		strcpy(snf->fileName,"fil_dm.snf");
		snf->dm=(SNF_DM *)calloc(snf->nSNF3,sizeof(SNF_DM));
		for(i=0;i<snf->nSNF3;i++){
			sprintf(snf->dm[i].name,"dm%03d",i+1);
			snf->dm[i].n=0;
			snf->dm[i].type=1;
			snf->dm[i].complex=0;
			snf->dm[i].ini=0.;
			snf->dm[i].dx=1.;
			snf->dm[i].m=0;
			snf->dm[i].ini2=0.;
			snf->dm[i].dx2=1.;
			//snf->dm[i].nform=malloc(8*sizeof(char));
			//snf->dm[i].sform=malloc(10*sizeof(char));
			strcpy(snf->dm[i].nform,"float");
			strcpy(snf->dm[i].sform,"nosparse");
		}
	}
	else if(!strcmp(type,"ds")){
		strcpy(snf->obj,"ds");
		//strcpy(snf->fileName,"fil_ds.snf");
//OBSOLETE		snf->ds=(SNF_DS *)calloc(snf->nSNF3,sizeof(SNF_DS));
		snf->field=(SNF_FIELD *)calloc(snf->nSNF3,sizeof(SNF_FIELD));
//OBSOLETE		snf->ds->nch=snf->nSNF3;
		for(i=0;i<snf->nSNF3;i++){
			sprintf(snf->field[i].name,"ch%03d",i+1);
			snf->field[i].complex=0;
			snf->field[i].len=0;
			snf->field[i].dt=1.;
			snf->field[i].del=0;
			strcpy(snf->field[i].nform,"float");
			strcpy(snf->field[i].sform,"nosparse");
		}
	}
	else if(!strcmp(type,"tfm")){
		strcpy(snf->obj,"tfm");
		//strcpy(snf->fileName,"fil_tfm.snf");
		snf->tfm=(SNF_TFM *)calloc(snf->nSNF3,sizeof(SNF_TFM));
		for(i=0;i<snf->nSNF3;i++){
			sprintf(snf->tfm[i].name,"tfm%03d",i+1);
			snf->tfm[i].len=0;
			snf->tfm[i].inif=0;
			snf->tfm[i].df=0;
		}
	}
	else if(!strcmp(type,"tf_sfdb")){
		strcpy(snf->obj,"tf_sfdb");
		//strcpy(snf->fileName,"fil_sfdb.snf");
		snf->tfm=(SNF_TFM *)calloc(snf->nSNF3,sizeof(SNF_TFM));
		for(i=0;i<snf->nSNF3;i++){
			sprintf(snf->tfm[i].name,"tf_sfdb%03d",i+1);
			snf->tfm[i].len=0;
			snf->tfm[i].inif=0;
			snf->tfm[i].df=0;
		}
	}
	else if(!strcmp(type,"tf_spec")){
		strcpy(snf->obj,"tf_spec");
		snf->tfm=(SNF_TFM *)calloc(snf->nSNF3,sizeof(SNF_TFM));
		for(i=0;i<snf->nSNF3;i++){
			sprintf(snf->tfm[i].name,"tf_spec%03d",i+1);
			snf->tfm[i].len=0;
			snf->tfm[i].inif=0;
			snf->tfm[i].df=0;
		}
	}
	else if(!strcmp(type,"tf_peak")){
		strcpy(snf->obj,"tf_peak");
		snf->tfm=(SNF_TFM *)calloc(snf->nSNF3,sizeof(SNF_TFM));
		for(i=0;i<snf->nSNF3;i++){
			sprintf(snf->tfm[i].name,"tf_peak%03d",i+1);
			snf->tfm[i].len=0;
			snf->tfm[i].inif=0;
			snf->tfm[i].df=0;
		}
	}
	else if(!strcmp(type,"hm")){
		strcpy(snf->obj,"hm");
		//strcpy(snf->fileName,"fil_hm.snf");
		snf->hm=(SNF_HM *)calloc(snf->nSNF3,sizeof(SNF_HM));
		for(i=0;i<snf->nSNF3;i++){
			sprintf(snf->hm[i].name,"hm%03d",i+1);
			snf->hm[i].totlen=0;
			snf->hm[i].lamlen=0;
		}
	}
/*	else if(!strcmp(type,"ev"))
		snf->ev=(SNF_EV *)calloc(sizeof(snf->nSNF3,SNF_EV)); */

	return snf;
}


int del_snf_struct(char* type){
}


SNF_STRUCT *all_snf_struct(SNF_STRUCT* snf_str){
	/* allocates arrays for snf structure */

	snf_str->caption[0]=(char *)calloc(snf_str->nSNF2*122,sizeof(char));
	strcpy(snf_str->caption[0]," ");

	if(!strcmp(snf_str->obj,"gd"))
		snf_str->gd=(SNF_GD *)calloc(snf_str->nSNF3,sizeof(SNF_GD));
	else if(!strcmp(snf_str->obj,"dm"))
		snf_str->dm=(SNF_DM *)calloc(snf_str->nSNF3,sizeof(SNF_DM));
	else if(!strcmp(snf_str->obj,"ds"))
		snf_str->field=(SNF_FIELD *)calloc(snf_str->nSNF3,sizeof(SNF_FIELD));
	else if(!strcmp(snf_str->obj,"tfm"))
		snf_str->tfm=(SNF_TFM *)calloc(snf_str->nSNF3,sizeof(SNF_TFM));
	else if(!strcmp(snf_str->obj,"tf_sfdb"))
		snf_str->tfm=(SNF_TFM *)calloc(snf_str->nSNF3,sizeof(SNF_TFM));
	else if(!strcmp(snf_str->obj,"tf_spec"))
		snf_str->tfm=(SNF_TFM *)calloc(snf_str->nSNF3,sizeof(SNF_TFM));
	else if(!strcmp(snf_str->obj,"tf_peak"))
		snf_str->tfm=(SNF_TFM *)calloc(snf_str->nSNF3,sizeof(SNF_TFM));
	else if(!strcmp(snf_str->obj,"hm"))
		snf_str->hm=(SNF_HM *)calloc(snf_str->nSNF3,sizeof(SNF_HM));

	if(snf_str->headv->lheadi>0)
		snf_str->headv->hi=(SNF5 *)calloc(snf_str->headv->lheadi,sizeof(SNF5));
	if(snf_str->headv->lheadd>0)
		snf_str->headv->hd=(SNF5 *)calloc(snf_str->headv->lheadd,sizeof(SNF5));
	if(snf_str->headv->lheads>0)
		snf_str->headv->hs=(SNF5 *)calloc(snf_str->headv->lheads,sizeof(SNF5));
	if(snf_str->headv->lrheadi>0)
		snf_str->headv->rhi=(SNF5 *)calloc(snf_str->headv->lrheadi,sizeof(SNF5));
	if(snf_str->headv->lrheadd>0)
		snf_str->headv->rhd=(SNF5 *)calloc(snf_str->headv->lrheadd,sizeof(SNF5));
	if(snf_str->headv->lrheads>0)
		snf_str->headv->rhs=(SNF5 *)calloc(snf_str->headv->lrheads,sizeof(SNF5));
	if(snf_str->headv->lfheadi>0)
		snf_str->headv->fhi=(SNF5 *)calloc(snf_str->headv->lfheadi,sizeof(SNF5));
	if(snf_str->headv->lfheadd>0)
		snf_str->headv->fhd=(SNF5 *)calloc(snf_str->headv->lfheadd,sizeof(SNF5));
	if(snf_str->headv->lfheads>0)
		snf_str->headv->fhs=(SNF5 *)calloc(snf_str->headv->lfheads,sizeof(SNF5));
	
	if(snf_str->nSNF6>0)
		snf_str->strSNF6[0]=(char *)calloc(120*snf_str->nSNF6,sizeof(char));


	return snf_str;
}


SNF_STRUCT *more_comments(SNF_STRUCT* snf_str, char** comm, int ncomlin){
	int i;

	snf_str->nSNF2=ncomlin;

	for(i=0;i<ncomlin;i++){
		//printf(" -> %s\n",comm[i]);
		//snf_str->caption[i]=(char *)calloc(122,sizeof(char));
		snf_str->caption[i]=comm[i];
		//printf(" => %s\n",snf_str->caption[i]);
	}

	return snf_str;
}


SNF_STRUCT *all_header_struct(SNF_STRUCT* snf_str){
	/* allocates arrays for snf header structures */

	if(snf_str->headv->lheadi>0)
		snf_str->headv->hi=(SNF5 *)calloc(snf_str->headv->lheadi,sizeof(SNF5));
	if(snf_str->headv->lheadd>0)
		snf_str->headv->hd=(SNF5 *)calloc(snf_str->headv->lheadd,sizeof(SNF5));
	if(snf_str->headv->lheads>0)
		snf_str->headv->hs=(SNF5 *)calloc(snf_str->headv->lheads,sizeof(SNF5));
	if(snf_str->headv->lrheadi>0)
		snf_str->headv->rhi=(SNF5 *)calloc(snf_str->headv->lrheadi,sizeof(SNF5));
	if(snf_str->headv->lrheadd>0)
		snf_str->headv->rhd=(SNF5 *)calloc(snf_str->headv->lrheadd,sizeof(SNF5));
	if(snf_str->headv->lrheads>0)
		snf_str->headv->rhs=(SNF5 *)calloc(snf_str->headv->lrheads,sizeof(SNF5));
	if(snf_str->headv->lfheadi>0)
		snf_str->headv->fhi=(SNF5 *)calloc(snf_str->headv->lfheadi,sizeof(SNF5));
	if(snf_str->headv->lfheadd>0)
		snf_str->headv->fhd=(SNF5 *)calloc(snf_str->headv->lfheadd,sizeof(SNF5));
	if(snf_str->headv->lfheads>0)
		snf_str->headv->fhs=(SNF5 *)calloc(snf_str->headv->lfheads,sizeof(SNF5));

	return snf_str;
}


SNF_STRUCT* SNF5_lines_i(SNF_STRUCT* snf_str, int n, char** name, char** capt, long* ival){
	int i;

	snf_str->headv->lheadi=n;
	snf_str->headv->hi=(SNF5 *)calloc(snf_str->headv->lheadi,sizeof(SNF5));

	for(i=0;i<snf_str->headv->lheadi;i++){
		strcpy(snf_str->headv->hi[i].name,name[i]);
		strcpy(snf_str->headv->hi[i].caption,capt[i]);
		snf_str->headv->hi[i].ivalue=ival[i];
	}

	return snf_str;
}


SNF_STRUCT* SNF5_lines_d(SNF_STRUCT* snf_str, char** name, char** capt, double* val){
	int i;

	snf_str->headv->hd=(SNF5 *)calloc(snf_str->headv->lheadd,sizeof(SNF5));

	for(i=0;i<snf_str->headv->lheadd;i++){
		strcpy(snf_str->headv->hd[i].name,name[i]);
		strcpy(snf_str->headv->hd[i].caption,capt[i]);
		snf_str->headv->hd[i].dvalue=val[i];
	}

	return snf_str;
}


SNF_STRUCT* SNF5_lines_s(SNF_STRUCT* snf_str, char** name, char** capt, char** sval){
	int i;

	snf_str->headv->hs=(SNF5 *)calloc(snf_str->headv->lheads,sizeof(SNF5));

	for(i=0;i<snf_str->headv->lheads;i++){
		strcpy(snf_str->headv->hs[i].name,name[i]);
		strcpy(snf_str->headv->hs[i].caption,capt[i]);
		strcpy(snf_str->headv->hs[i].svalue,sval[i]);
	}

	return snf_str;
}


SNF_STRUCT* SNF5_lines_ri(SNF_STRUCT* snf_str, char** name, char** capt, long* ival){
	int i;

	snf_str->headv->rhi=(SNF5 *)calloc(snf_str->headv->lrheadi,sizeof(SNF5));

	for(i=0;i<snf_str->headv->lrheadi;i++){
		strcpy(snf_str->headv->rhi[i].name,name[i]);
		strcpy(snf_str->headv->rhi[i].caption,capt[i]);
		snf_str->headv->rhi[i].ivalue=ival[i];
	}

	return snf_str;
}


SNF_STRUCT* SNF5_lines_rd(SNF_STRUCT* snf_str, char** name, char** capt, double* val){
	int i;

	snf_str->headv->rhd=(SNF5 *)calloc(snf_str->headv->lrheadd,sizeof(SNF5));

	for(i=0;i<snf_str->headv->lrheadd;i++){
		strcpy(snf_str->headv->rhd[i].name,name[i]);
		strcpy(snf_str->headv->rhd[i].caption,capt[i]);
		snf_str->headv->rhd[i].dvalue=val[i];
	}

	return snf_str;
}


SNF_STRUCT* SNF5_lines_rs(SNF_STRUCT* snf_str, char** name, char** capt, char** sval){
	int i;

	snf_str->headv->rhs=(SNF5 *)calloc(snf_str->headv->lrheads,sizeof(SNF5));

	for(i=0;i<snf_str->headv->lrheads;i++){
		strcpy(snf_str->headv->rhs[i].name,name[i]);
		strcpy(snf_str->headv->rhs[i].caption,capt[i]);
		strcpy(snf_str->headv->rhs[i].svalue,sval[i]);
	}

	return snf_str;
}


SNF_STRUCT* SNF5_lines_fi(SNF_STRUCT* snf_str, char** name, char** capt, long* ival){
	int i;

	snf_str->headv->fhi=(SNF5 *)calloc(snf_str->headv->lfheadi,sizeof(SNF5));

	for(i=0;i<snf_str->headv->lfheadi;i++){
		strcpy(snf_str->headv->fhi[i].name,name[i]);
		strcpy(snf_str->headv->fhi[i].caption,capt[i]);
		snf_str->headv->fhi[i].ivalue=ival[i];
	}

	return snf_str;
}


SNF_STRUCT* SNF5_lines_fd(SNF_STRUCT* snf_str, char** name, char** capt, double* val){
	int i;

	snf_str->headv->fhd=(SNF5 *)calloc(snf_str->headv->lfheadd,sizeof(SNF5));

	for(i=0;i<snf_str->headv->lheadd;i++){
		strcpy(snf_str->headv->fhd[i].name,name[i]);
		strcpy(snf_str->headv->fhd[i].caption,capt[i]);
		snf_str->headv->fhd[i].dvalue=val[i];
	}

	return snf_str;
}


SNF_STRUCT* SNF5_lines_fs(SNF_STRUCT* snf_str, char** name, char** capt, char** sval){
	int i;

	snf_str->headv->hs=(SNF5 *)calloc(snf_str->headv->lfheads,sizeof(SNF5));

	for(i=0;i<snf_str->headv->lfheads;i++){
		strcpy(snf_str->headv->fhs[i].name,name[i]);
		strcpy(snf_str->headv->fhs[i].caption,capt[i]);
		strcpy(snf_str->headv->fhs[i].svalue,sval[i]);
	}

	return snf_str;
}


SNF_STRUCT* SNF6_lines(SNF_STRUCT* snf_str, char** lin, int nSNF6lin){
	int i;

	snf_str->nSNF6=nSNF6lin;
	snf_str->strSNF6[0]=(char *)calloc(120*snf_str->nSNF6,sizeof(char));

	for(i=0;i<snf_str->nSNF6;i++)
		snf_str->strSNF6[i]=lin[i];
	
	return snf_str;
}


int snf_open_w(SNF_STRUCT *w_str){
	FILE *fil;
	int status,recfilhead=0;
	char *A;
	char *proto;
	long i;

	A=(char *)calloc(128,sizeof(char));

	proto=SNF_PROTOCOL;

	printf("nSNF3 = %d \n",w_str->nSNF3);

	if(!strcmp(w_str->binasc,"ascii"))
		fil=fopen(w_str->fileName, "wt");
	else
		fil=fopen(w_str->fileName, "wb");

	w_str->file=fil;
	w_str->headerlen=0;

	clear_string(A,128);

/* -------------- #SNF#1# ----------------- */

	sprintf(A,"#SNF#1#%s|%s|%s|%s|%s| \r\n", proto, w_str->obj, 
		w_str->binasc, w_str->directaccess, w_str->userprot);

	if(w_str->access->verb>0)
		printf(" %s created as a %s %s file\n",w_str->fileName,w_str->obj,
			w_str->binasc);

	status=fwrite(A,1,128,fil);
	w_str->headerlen+=128;

/* -------------- #SNF#2# ----------------- */

	if(w_str->access->verb>0)
		printf(" %d SNF2 records\n",w_str->nSNF2);


	for(i=0;i<w_str->nSNF2;i++){
		clear_string(A,128);

		sprintf(A,"#SNF#2# %s \r\n",w_str->caption[i],116);
		fseek(w_str->file,w_str->headerlen,SEEK_SET);
		status=fwrite(A,1,128,fil);
		w_str->headerlen+=128;
	}

/* -------------- #SNF#3# ----------------- */

	if(w_str->access->verb>0)
		printf(" %d SNF3 records\n",w_str->nSNF3);

	if(!strcmp(w_str->obj,"gd"))
		for(i=0;i<w_str->nSNF3;i++){
			clear_string(A,128);
			sprintf(A,"#SNF#3# %16s %d %d %d %#g %#g \r\n",w_str->gd[i].name,
				w_str->gd[i].n,w_str->gd[i].type,w_str->gd[i].complex,
				w_str->gd[i].ini,w_str->gd[i].dx);
		
			fseek(w_str->file,w_str->headerlen,SEEK_SET);
			status=fwrite(A,1,128,fil);
			w_str->headerlen+=128;
		}
	else if(!strcmp(w_str->obj,"dm")){
		for(i=0;i<w_str->nSNF3;i++){
			clear_string(A,128);

			sprintf(A,"#SNF#3# %16s %d %d %d %g %g %d %g %g \r\n",
				w_str->dm[i].name, w_str->dm[i].n, w_str->dm[i].type, 
				w_str->dm[i].complex, w_str->dm[i].ini, w_str->dm[i].dx, 
				w_str->dm[i].m, w_str->dm[i].ini2, w_str->dm[i].dx2);

			fseek(w_str->file,w_str->headerlen,SEEK_SET);
			status=fwrite(A,1,128,fil);
			w_str->headerlen+=128;
		}
	}
	else if(!strcmp(w_str->obj,"ds")){
		printf(" *** %d \n",w_str->nSNF3);
		for(i=0;i<w_str->nSNF3;i++){
			clear_string(A,128);

			sprintf(A,"#SNF#3# %16s %8s %d %g \r\n",
				w_str->field[i].name, w_str->field[i].nform, 
				w_str->field[i].len, w_str->field[i].dt);

			fseek(w_str->file,w_str->headerlen,SEEK_SET);
			status=fwrite(A,1,128,fil);
			w_str->headerlen+=128;
		}
		recfilhead=1;
	}
	else if(!strcmp(w_str->obj,"tfm")){
		for(i=0;i<w_str->nSNF3;i++){
			clear_string(A,128);

			sprintf(A,"#SNF#3# %16s %8s %d %d %g\r\n'",
				&w_str->tfm[i].name, &w_str->tfm[i].type, w_str->tfm[i].len,
				w_str->tfm[i].inif, w_str->tfm[i].df);

			fseek(w_str->file,w_str->headerlen,SEEK_SET);
			status=fwrite(A,1,128,fil);
			w_str->headerlen+=128;
		}
		recfilhead=1;
	}
	else if(!strcmp(w_str->obj,"tf_sfdb")){
		for(i=0;i<w_str->nSNF3;i++){
			clear_string(A,128);

			sprintf(A,"#SNF#3# %16s %8s %d %d %g\r\n'",
				&w_str->tfm[i].name, &w_str->tfm[i].type, w_str->tfm[i].len,
				w_str->tfm[i].inif, w_str->tfm[i].df);

			fseek(w_str->file,w_str->headerlen,SEEK_SET);
			status=fwrite(A,1,128,fil);
			w_str->headerlen+=128;
		}
		recfilhead=1;
	}
	else if(!strcmp(w_str->obj,"tf_spec")){
		for(i=0;i<w_str->nSNF3;i++){
			clear_string(A,128);

			sprintf(A,"#SNF#3# %16s %d %g %g\r\n",
				w_str->tfm[i].name, w_str->tfm[i].len, 
				w_str->tfm[i].inif, w_str->tfm[i].df);

			fseek(w_str->file,w_str->headerlen,SEEK_SET);
			status=fwrite(A,1,128,fil);
			w_str->headerlen+=128;
		}
		recfilhead=1;
	}
	else if(!strcmp(w_str->obj,"tf_peak")){
		for(i=0;i<w_str->nSNF3;i++){
			clear_string(A,128);

			sprintf(A,"#SNF#3# %16s %d %g %g\r\n",
				w_str->tfm[i].name, w_str->tfm[i].len, 
				w_str->tfm[i].inif, w_str->tfm[i].df);

			fseek(w_str->file,w_str->headerlen,SEEK_SET);
			status=fwrite(A,1,128,fil);
			w_str->headerlen+=128;
		}
		recfilhead=1;
	}
	else if(!strcmp(w_str->obj,"hm")){
		for(i=0;i<w_str->nSNF3;i++){
			clear_string(A,128);

			sprintf(A,"#SNF#3# %16s %d %d \r\n",
				w_str->hm[i].name, w_str->hm[i].totlen, 
				w_str->hm[i].lamlen);

			fseek(w_str->file,w_str->headerlen,SEEK_SET);
			status=fwrite(A,1,128,fil);
			w_str->headerlen+=128;
		}
		recfilhead=1;
	}
	else if(w_str->obj=="psc"){
		recfilhead=1;
	}
	else if(w_str->obj=="ev"){
		recfilhead=1;
	}

	recfilhead=recfilhead+w_str->headv->lheadi+w_str->headv->lheadd+w_str->headv->lheads+
		w_str->headv->lrheadi+w_str->headv->lrheadd+w_str->headv->lrheads+
		w_str->headv->lfheadi+w_str->headv->lfheadd+w_str->headv->lfheads;

	if(recfilhead>0){
		clear_string(A,128);

/* -------------- #SNF#4# ----------------- */

		sprintf(A,"#SNF#4# %d %d %d %d %d %d %d %d %d",w_str->headv->lheadi, 
			w_str->headv->lheadd,w_str->headv->lheads,w_str->headv->lrheadi, 
			w_str->headv->lrheadd,w_str->headv->lrheads,w_str->headv->lfheadi,
			w_str->headv->lfheadd,w_str->headv->lfheads);
		status=fwrite(A,1,128,fil);
		w_str->headerlen+=128;

/* -------------- #SNF#5# ----------------- */

		for(i=0;i<w_str->headv->lheadi;i++){
			clear_string(A,128);

			sprintf(A,"#SNF#5#%2s%04d %16s %20d %70s \r\n","HI",i,
				w_str->headv->hi[i].name,w_str->headv->hi[i].ivalue,
				w_str->headv->hi[i].caption);
			status=fwrite(A,1,128,fil);
			w_str->headerlen+=128;
		}
		
		for(i=0;i<w_str->headv->lheadd;i++){
			clear_string(A,128);

			sprintf(A,"#SNF#5#%2s%04d %16s %20.12g %70s \r\n","HD",i,
				w_str->headv->hd[i].name,w_str->headv->hd[i].dvalue,
				w_str->headv->hd[i].caption);
			status=fwrite(A,1,128,fil);
			w_str->headerlen+=128;
		}
		
		for(i=0;i<w_str->headv->lheads;i++){
			clear_string(A,128);

			sprintf(A,"#SNF#5#%2s%04d %16s %16s %70s \r\n","HS",i,
				w_str->headv->hs[i].name,w_str->headv->hs[i].svalue,
				w_str->headv->hs[i].caption);
			status=fwrite(A,1,128,fil);
			w_str->headerlen+=128;
		}

		for(i=0;i<w_str->headv->lrheadi;i++){
			clear_string(A,128);

			sprintf(A,"#SNF#5#%2s%04d %16s %20d %70s \r\n","RI",i,
				w_str->headv->rhi[i].name,w_str->headv->rhi[i].ivalue,
				w_str->headv->rhi[i].caption);
			status=fwrite(A,1,128,fil);
			w_str->headerlen+=128;
		}
	
		for(i=0;i<w_str->headv->lrheadd;i++){
			clear_string(A,128);

			sprintf(A,"#SNF#5#%2s%04d %16s %20.12g %70s \r\n","RD",i,
				w_str->headv->rhd[i].name,w_str->headv->rhd[i].dvalue,
				w_str->headv->rhd[i].caption);
			status=fwrite(A,1,128,fil);
			w_str->headerlen+=128;
		}
		
		for(i=0;i<w_str->headv->lrheads;i++){
			clear_string(A,128);

			sprintf(A,"#SNF#5#%2s%04d %16s %16s %70s \r\n","RS",i,
				w_str->headv->rhs[i].name,w_str->headv->rhs[i].svalue,
				w_str->headv->rhs[i].caption);
			status=fwrite(A,1,128,fil);
			w_str->headerlen+=128;
		}
		
		for(i=0;i<w_str->headv->lfheadi;i++){
			clear_string(A,128);

			sprintf(A,"#SNF#5#%2s%04d %16s %20d %70s \r\n","FI",i,
				w_str->headv->fhi[i].name,w_str->headv->fhi[i].ivalue,
				w_str->headv->fhi[i].caption);
			status=fwrite(A,1,128,fil);
			w_str->headerlen+=128;
		}
		
		for(i=0;i<w_str->headv->lfheadd;i++){
			clear_string(A,128);

			sprintf(A,"#SNF#5#%2s%04d %16s %20.12g %70s \r\n","FR",i,
				w_str->headv->fhd[i].name,w_str->headv->fhd[i].dvalue,
				w_str->headv->fhd[i].caption);
			status=fwrite(A,1,128,fil);
			w_str->headerlen+=128;
		}
		
		for(i=0;i<w_str->headv->lfheads;i++){
			clear_string(A,128);

			sprintf(A,"#SNF#5#%2s%04d %16s %16s %70s \r\n","FS",i,
				w_str->headv->fhs[i].name,w_str->headv->fhs[i].svalue,
				w_str->headv->fhs[i].caption);
			status=fwrite(A,1,128,fil);
			w_str->headerlen+=128;
		}
	}

/* -------------- #SNF#6# ----------------- */

	if(w_str->nSNF6>0)for(i=0;i<w_str->nSNF6;i++){
		clear_string(A,128);

		sprintf(A,"#SNF#6#%s\n",w_str->strSNF6[i]);
		status=fwrite(A,1,128,fil);
		w_str->headerlen+=128;
	}

	clear_string(A,128);
	status=fwrite(A,1,128,fil);

	fseek(fil,-128,SEEK_CUR);

	if(strcmp(w_str->binasc,"ascii")){

/* -------------- #NIDG# ----------------- */

		w_str->headerlen+=128;
		fwrite(&NIDG0,sizeof(struct SNF_NIDG),1,fil);
	}
	fseek(fil,w_str->headerlen,SEEK_SET);

	return w_str->headerlen;
}


int snf_derived(SNF_STRUCT *snf_str){ /*
	if(strcmp(snf_str->nform,"int8")){
	}
	else if(strcmp(snf_str->nform,"int16")){
	}
	else if(strcmp(snf_str->nform,"log8")){
	}
	else if(strcmp(snf_str->nform,"log16")){
	}
	else if(strcmp(snf_str->nform,"float")){
	}
	else if(strcmp(snf_str->nform,"double")){
	}
	else if(strcmp(snf_str->nform,"sparse8")){
	}
	else if(strcmp(snf_str->nform,"sparse")){
	} */
}


SNF_STRUCT *snf_all_ds(char *chname, long *chlen, double *dt, int nds, 
					   char *filename, char *caption, char **nform)
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
{
	SNF_STRUCT *w_str;
	int i;

	w_str=malloc(sizeof(SNF_STRUCT));
	strcpy(w_str->protocol,SNF_PROTOCOL);
	strcpy(w_str->obj,"ds");
	strcpy(w_str->userprot," ");
	w_str->fileName=filename;
	strcpy(w_str->binasc,"bin");
	strcpy(w_str->directaccess,"directaccess");
	strcpy(w_str->xunit,"mjd");

	strcpy(w_str->binasc,"bin");
	w_str->file=fopen(w_str->fileName, "wb");

	w_str->nSNF2=1;
	w_str->nSNF3=nds;
	w_str->nSNF5=1;
	w_str->nSNF6=0;
	w_str->caption[0]=(char *)malloc(122*sizeof(char));
	strcpy(w_str->caption[0],caption);

	w_str->access=malloc(sizeof(SNF_ACCESS));
	w_str->access->verb=1;

	w_str->headv=malloc(sizeof(SNF_HEADV));
	w_str->headv->lheadi=0;
	w_str->headv->lheadd=0;
	w_str->headv->lheads=0;
	w_str->headv->lrheadi=0;
	w_str->headv->lrheadd=1;
	w_str->headv->lrheads=0;
	w_str->headv->lfheadi=0;
	w_str->headv->lfheadd=1;
	w_str->headv->lfheads=0;
	
	w_str->headv->rhd=malloc(sizeof(SNF5));
	strcpy(w_str->headv->rhd[0].caption,
		"time of the first sample of every channel");
	strcpy(w_str->headv->rhd[0].name,"t0");
	w_str->headv->rhd[0].dvalue=0;

	w_str->headv->fhd=malloc(sizeof(SNF5));
	strcpy(w_str->headv->fhd[0].caption,
		"field delay");
	strcpy(w_str->headv->fhd[0].name,"dt0");
	w_str->headv->fhd[0].dvalue=0;

//OBSOLETE	w_str->ds=malloc(nds*sizeof(SNF_DS));
//OBSOLETE	w_str->ds->nch=nds;
	w_str->field=malloc(nds*sizeof(SNF_FIELD));
	for(i=0;i<nds;i++){
		w_str->field[i].complex=0;
		w_str->field[i].len=chlen[i];
		w_str->field[i].dt=dt[i];
		strcpy(w_str->field[i].name,chname+20*i);
		w_str->field[i].del=0;
		strcpy(w_str->field[i].nform,nform[i]);
	}
	
	return w_str;
}


SNF_STRUCT *snf_open_tfm(char *filename, char *caption, char *nform, char *userprot,
						 double dt, char *typ, long lfft, long interlac, char *window,
						 int nband, long *blen, long *inif, double df)
/*  
    Creates the structure for archiving a tfm (or a multiple tfm) of type typ.

		filename
		caption		max 121 bytes
		nform		equal for all the bands
		userprot	is the user protocol (for example "sfdb_mar02")
		dt			sampling time
		typ			can be "sfdb", "spec" or "peak"
		lfft		FFT length (total, in samples)
		interlac	interlacing (in samples)
		window		"no", "bartlett", "hanning", ...
		nband		number of bands
		blen[]		the number of data for each band
	
	If the LogX format is used, one or more  LogX structures should be 
	created out of this routine.

    This routine creates the easiest structure, that is with:
		- one caption line
		- the mandatory header fields
	It calls snf_open_w.
	If one needs more complex features,as more caption lines or (file, record, field) header
	variables, snf_close should be called and non-default items should be added or modified 
	after this call; then the snf_open_w should be called.
*/
{
	SNF_STRUCT *w_str;
	int i;
	char sobj[8]="tf_";

	strncat(sobj,typ,4);
	w_str=malloc(sizeof(SNF_STRUCT));
	strcpy(w_str->protocol,SNF_PROTOCOL);
	strcpy(w_str->obj,sobj);
	strcpy(w_str->userprot,userprot);
	w_str->fileName=filename;
	strcpy(w_str->binasc,"bin");
	strcpy(w_str->directaccess,"directaccess");

	strcpy(w_str->binasc,"bin");
	w_str->file=fopen(w_str->fileName, "wb");

	w_str->nSNF2=1;
	w_str->nSNF3=nband;
	//w_str->nSNF5=19;
	w_str->nSNF6=0;
	w_str->caption[0]=(char *)malloc(122*sizeof(char));
	strcpy(w_str->caption[0],caption);

	w_str->access=malloc(sizeof(SNF_ACCESS));
	w_str->access->verb=1;

	w_str->headv=malloc(sizeof(SNF_HEADV));
	w_str->headv->lheadi=2;
	w_str->headv->lheadd=1;
	w_str->headv->lheads=1;
	w_str->headv->lrheadi=1;
	w_str->headv->lrheadd=14;
	w_str->headv->lrheads=0;
	w_str->headv->lfheadi=0;
	w_str->headv->lfheadd=0;
	w_str->headv->lfheads=0;

	w_str->headv->hi=(SNF5 *)malloc(w_str->headv->lheadi*sizeof(SNF5));
	w_str->headv->hd=(SNF5 *)malloc(w_str->headv->lheadd*sizeof(SNF5));
	w_str->headv->hs=(SNF5 *)malloc(w_str->headv->lheads*sizeof(SNF5));

	w_str->headv->rhi=(SNF5 *)malloc(w_str->headv->lrheadi*sizeof(SNF5));
	w_str->headv->rhd=(SNF5 *)malloc(w_str->headv->lrheadd*sizeof(SNF5));
	
	strcpy(w_str->headv->hi[0].caption,
		"Total length of the FFT");
	strcpy(w_str->headv->hi[0].name,"lfft");
	w_str->headv->hi[0].ivalue=lfft;

	strcpy(w_str->headv->hi[1].caption,
		"Number of interlacing samples");
	strcpy(w_str->headv->hi[1].name,"interlac");
	w_str->headv->hi[1].ivalue=interlac;

	strcpy(w_str->headv->hd[0].caption,
		"Sampling time");
	strcpy(w_str->headv->hd[0].name,"dt");
	w_str->headv->hd[0].dvalue=dt;

	strcpy(w_str->headv->hs[0].caption,
		"Applied window");
	strcpy(w_str->headv->hs[0].name,"window");
	strcpy(w_str->headv->hs[0].svalue,window);
	
	strcpy(w_str->headv->rhi[0].caption,
		"number of samples 0-padded");
	strcpy(w_str->headv->rhi[0].name,"n_0pad");
	w_str->headv->rhi[0].ivalue=0;
	
	strcpy(w_str->headv->rhd[0].caption,
		"mjd time of the first datum");
	strcpy(w_str->headv->rhd[0].name,"t_mjd");
	w_str->headv->rhd[0].dvalue=0;
	
	strcpy(w_str->headv->rhd[1].caption,
		"gps time of the first datum");
	strcpy(w_str->headv->rhd[1].name,"t_gps");
	w_str->headv->rhd[1].dvalue=0;

	strcpy(w_str->headv->rhd[2].caption,
		"velocity of the detector at the beginning (fraction of c) - x");
	strcpy(w_str->headv->rhd[2].name,"v_ini_x");
	w_str->headv->rhd[2].dvalue=0;

	strcpy(w_str->headv->rhd[3].caption,
		"velocity of the detector at the beginning (fraction of c) - y");
	strcpy(w_str->headv->rhd[3].name,"v_ini_y");
	w_str->headv->rhd[3].dvalue=0;

	strcpy(w_str->headv->rhd[4].caption,
		"velocity of the detector at the beginning (fraction of c) - z");
	strcpy(w_str->headv->rhd[4].name,"v_ini_z");
	w_str->headv->rhd[4].dvalue=0;

	strcpy(w_str->headv->rhd[5].caption,
		"velocity of the detector at the middle (fraction of c) - x");
	strcpy(w_str->headv->rhd[5].name,"v_mid_x");
	w_str->headv->rhd[5].dvalue=0;

	strcpy(w_str->headv->rhd[6].caption,
		"velocity of the detector at the middle (fraction of c) - y");
	strcpy(w_str->headv->rhd[6].name,"v_mid_y");
	w_str->headv->rhd[6].dvalue=0;

	strcpy(w_str->headv->rhd[7].caption,
		"velocity of the detector at the middle (fraction of c) - z");
	strcpy(w_str->headv->rhd[7].name,"v_mid_z");
	w_str->headv->rhd[7].dvalue=0;

	strcpy(w_str->headv->rhd[8].caption,
		"velocity of the detector at the end (fraction of c) - x");
	strcpy(w_str->headv->rhd[8].name,"v_fin_x");
	w_str->headv->rhd[8].dvalue=0;

	strcpy(w_str->headv->rhd[9].caption,
		"velocity of the detector at the end (fraction of c) - y");
	strcpy(w_str->headv->rhd[9].name,"v_fin_y");
	w_str->headv->rhd[9].dvalue=0;

	strcpy(w_str->headv->rhd[10].caption,
		"velocity of the detector at the end (fraction of c) - z");
	strcpy(w_str->headv->rhd[10].name,"v_fin_z");
	w_str->headv->rhd[10].dvalue=0;

	strcpy(w_str->headv->rhd[11].caption,
		"mean velocity of the detector (fraction of c) - x");
	strcpy(w_str->headv->rhd[11].name,"v_mean_x");
	w_str->headv->rhd[11].dvalue=0;

	strcpy(w_str->headv->rhd[12].caption,
		"mean velocity of the detector (fraction of c) - y");
	strcpy(w_str->headv->rhd[12].name,"v_mean_y");
	w_str->headv->rhd[12].dvalue=0;

	strcpy(w_str->headv->rhd[13].caption,
		"mean velocity of the detector (fraction of c) - z");
	strcpy(w_str->headv->rhd[13].name,"v_mean_z");
	w_str->headv->rhd[13].dvalue=0;

	w_str->field=malloc(nband*sizeof(SNF_FIELD));
	w_str->tfm=malloc(nband*sizeof(SNF_TFM));

	for(i=0;i<nband;i++){
		w_str->field[i].complex=1;
		if(strcmp(typ,"sfdb"))w_str->field[i].complex=0;
		w_str->field[i].len=blen[i];
		w_str->field[i].dt=1./(dt*lfft);
		sprintf(w_str->field[i].name,"band%00d",i+1);
		w_str->field[i].del=0;
		strcpy(w_str->field[i].nform,nform);

		sprintf(&w_str->tfm[i].name,"band%00d",i+1);
		sprintf(w_str->tfm[i].type,typ,i+1);
		w_str->tfm[i].len=blen[i];
		w_str->tfm[i].inif=inif[i];
		w_str->tfm[i].df=df;
	}

	snf_open_w(w_str);
	
	return w_str;
}


SNF_STRUCT *snf_open_hm(SNF_STRUCT *w_str)
{
	SNF_STRUCT *str;

	
	return str;
}


SNF_STRUCT *snf_open_psc(SNF_STRUCT *w_str)
{
	SNF_STRUCT *str;

	
	return str;
}


/* ----------------- Write --------------------*/


int snf_write_rec(SNF_STRUCT *str, long *oldrecini){
/* writes only the header and complete the header of the previous record 
   at the first record *oldrecini should be set to 0
   the record header variables should be set before the call */ 
	long i,lr;
	long ini,lenrec=-1;
	char* label="#SNF#ORD";

	if(!strcmp(str->binasc,"ascii")){
		str->access->nrec++;
		fprintf(str->file,"%s %d\n",label,str->access->nrec);

		for(i=0;i<str->headv->lfheadi;i++)
			fprintf(str->file,"%d\n",str->headv->fhi[i].ivalue);
		for(i=0;i<str->headv->lfheadd;i++)
			fprintf(str->file,"%g\n",str->headv->fhd[i].dvalue);
		for(i=0;i<str->headv->lfheads;i++)
			fprintf(str->file,"%s\n",str->headv->fhs[i].svalue);

		return 0;
	}
	
	ini=ftell(str->file);
	if(*oldrecini>0){
		fseek(str->file,*oldrecini+12,SEEK_SET);
		lr=ini-*oldrecini;
		fwrite(&lr,4,1,str->file);
		fseek(str->file,ini,SEEK_SET);
	}
	else str->access->nrec=0;

	*oldrecini=ini;
	str->access->nrec++;
	fwrite(label,1,8,str->file);
	fwrite(&str->access->nrec,4,1,str->file);
	fwrite(&lenrec,4,1,str->file);

	for(i=0;i<str->headv->lrheadi;i++)
		fwrite(&str->headv->rhi[i].ivalue,4,1,str->file);
	for(i=0;i<str->headv->lrheadd;i++)
		fwrite(&str->headv->rhd[i].dvalue,8,1,str->file);
	for(i=0;i<str->headv->lrheads;i++)
		fwrite(str->headv->rhs[i].svalue,1,16,str->file);

	return 0;
}


int snf_write_recextra(SNF_STRUCT *str, long *oldrecini){
/* writes only the header and complete the header of the previous record 
   at the first record *oldrecini should be set to 0 
   str->access->nrecextra should be 0 at the beginning;
   the record header variables should be set before the call */ 
	long i,lr;
	long ini,lenrec=-1;
	char* label="#SNF#ORD";

	if(!strcmp(str->binasc,"ascii")){
		str->access->nrec++;
		fprintf(str->file,"%s %d\n",label,str->access->nrec);

		return 0;
	}
	
	ini=ftell(str->file);
	if(*oldrecini>0){
		fseek(str->file,*oldrecini+12,SEEK_SET);
		lr=ini-*oldrecini;
		fwrite(&lr,4,1,str->file);
		fseek(str->file,ini,SEEK_SET);
	}

	*oldrecini=ini;
	str->access->nrecextra++;
	fwrite(label,1,8,str->file);
	fwrite(&str->access->nrec,4,1,str->file);
	fwrite(&lenrec,4,1,str->file);

	return 0;
}


int snf_write_field(SNF_STRUCT *str,  long *prevfieldini){
	/* writes only the header and complete the header of the previous field 
       at the beginning *prevfieldini should be set to 0;
	   the field header variables should be set before the call */ 
	long i,lr;
	long ini,lenfiel=-1;

	if(!strcmp(str->binasc,"ascii")){
		for(i=0;i<str->headv->lfheadi;i++)
			fprintf(str->file,"%d\n",str->headv->fhi[i].ivalue);
		for(i=0;i<str->headv->lfheadd;i++)
			fprintf(str->file,"%g\n",str->headv->fhd[i].dvalue);
		for(i=0;i<str->headv->lfheads;i++)
			fprintf(str->file,"%s\n",str->headv->fhs[i].svalue);

		return 0;
	}
	
	ini=ftell(str->file);
	if(*prevfieldini>0){
		fseek(str->file,*prevfieldini,SEEK_SET);
		lr=ini-*prevfieldini;
		fwrite(&lr,4,1,str->file);
		fseek(str->file,ini,SEEK_SET);
	}

	*prevfieldini=ini;
	str->access->nrecextra++;
	fwrite(&lenfiel,4,1,str->file);

	for(i=0;i<str->headv->lfheadi;i++)
		fwrite(&str->headv->fhi[i].ivalue,4,1,str->file);
	for(i=0;i<str->headv->lfheadd;i++)
		fwrite(&str->headv->fhd[i].dvalue,8,1,str->file);
	for(i=0;i<str->headv->lfheads;i++)
		fwrite(str->headv->fhs[i].svalue,1,16,str->file);

	return 0;
}



int snf_write_vect(FILE* file, char* binasc, char* nform, char* sform, long n, float* f_v, 
				   LogX_STRUCT *logx, SpVec_STRUCT *spvect)
/*
	logx, spvect	externally allocated only if necessary (depending on nform and sform)
*/
{
	long nwrit,i;

	if(!strcmp(binasc,"ascii")){
		for(i=0;i<n;i++){
			fprintf(file," %g ",f_v[i]);
			if(((i+1)/10)*10==(i+1))fprintf(file,"\n");
		}
	}
	else{
		if(!strcmp(nform,"float")){
			nwrit=fwrite(f_v,sizeof(float),n,file);
		}
		else if(!strcmp(nform,"lin")){
			nwrit=fwrite(f_v,sizeof(float),n,file);
		}
	}

	return 1;
}


int snf_write_dvect(FILE* file, char* binasc, char* nform, char* sform, long n, double* d_v,
					struct LogX_STRUCT *logx, struct SpVec_STRUCT *spvect)
/*
*/
{
	long nwrit,i;

	
	if(!strcmp(binasc,"ascii")){
		for(i=0;i<n;i++){
			fprintf(file," %g ",d_v[i]);
			if(((i+1)/10)*10==(i+1))fprintf(file,"\n");
		}
	}
	else{
		if(!strcmp(nform,"float")){
			nwrit=fwrite(d_v,sizeof(float),n,file);
		}
	}

	return 1;
}


int snf_write_gd(GD* gd, long ngd, char *binasc, char* filgd){
/*  Stores data in an snf gd file.

	gd		gd array to store
	ngd		number of gds
	binasc	"bin" or "ascii"
	filgd	file name
*/
	SNF_STRUCT *w_str;
	struct LogX_STRUCT *logx;
	struct SpVec_STRUCT *spvect;
	int i,stat;
	unsigned long flen;
	long prevfieldini;
	
	if(ngd<1)ngd=1;
	w_str=snf_crea_struct("gd",ngd);

	if(!strcmp(binasc,"ascii"))strcpy(w_str->binasc,binasc);

	if(strlen(filgd)>1)w_str->fileName=filgd;

	w_str->nSNF2=ngd;
	
	for(i=0;i<ngd;i++){
		w_str->caption[i]=(char *)calloc(122,sizeof(char));
		w_str->caption[i]=gd[i].capt;

		sprintf(w_str->gd[i].name,"gd(%d)",i);
		w_str->gd[i].n=gd[i].n;
		w_str->gd[i].type=gd[i].type;
		w_str->gd[i].complex=gd[i].complex;
		w_str->gd[i].ini=gd[i].ini;
		w_str->gd[i].dx=gd[i].dx;
	}

	printf("open file %s \n",w_str->fileName);
	printf("gd[0] length %d \n",w_str->gd[0].n);
	printf("gd[1] length %d \n",w_str->gd[1].n);
	printf("gd[2] length %d \n",w_str->gd[2].n);

	stat=snf_open_w(w_str);
	prevfieldini=0;

	for(i=0;i<ngd;i++){
		flen=0;
//		ini=ftell(w_str->file);
//		if(strcmp(binasc,"ascii"))fwrite(&flen,4,1,w_str->file);
		snf_write_field(w_str,&prevfieldini);

		if(gd[i].type==2)snf_write_dvect(w_str->file,w_str->binasc,"double","nosparse",
			gd[i].n,gd[i].x,logx,spvect);
		snf_write_vect(w_str->file,w_str->binasc,"float","nosparse",
			gd[i].n*(gd[i].complex+1),gd[i].y,logx,spvect);
		
//		fin=ftell(w_str->file);
//		flen=fin-ini;
//		fseek(w_str->file,ini,SEEK_SET);
//		if(strcmp(binasc,"ascii"))fwrite(&flen,4,1,w_str->file);
//		fseek(w_str->file,fin,SEEK_SET);
	} 

	snf_close(w_str);

	return 1;
}


int snf_write_dm(DM* dm, long ndm, char *binasc, char* fildm)
/*  Stores data in an snf dm file.
	y data matrix is input by columns.
   
	dm		gd array to store
	ndm		number of gds
	binasc	"bin" or "ascii"
	filgd	file name
*/
{
	SNF_STRUCT *w_str;
	struct LogX_STRUCT *logx;
	struct SpVec_STRUCT *spvect;
	int i,stat;
	unsigned long flen;
	long prevfieldini;
	
	if(ndm<1)ndm=1;
	w_str=snf_crea_struct("dm",ndm);

	if(!strcmp(binasc,"ascii"))strcpy(w_str->binasc,binasc);

	if(strlen(fildm)>1)w_str->fileName=fildm;

	w_str->nSNF2=ndm;
	
	for(i=0;i<ndm;i++){
		w_str->caption[i]=(char *)calloc(122,sizeof(char));
		w_str->caption[i]=dm[i].capt;

		sprintf(w_str->dm[i].name,"dm(%d)",i);
		w_str->dm[i].n=dm[i].n;
		w_str->dm[i].type=dm[i].type;
		w_str->dm[i].complex=dm[i].complex;
		w_str->dm[i].ini=dm[i].ini;
		w_str->dm[i].dx=dm[i].dx;
		w_str->dm[i].m=dm[i].m;
		w_str->dm[i].ini2=dm[i].ini2;
		w_str->dm[i].dx2=dm[i].dx2;
	}

	printf("open file %s \n",w_str->fileName);
	printf("dm[0] length %d \n",w_str->dm[0].n);
	printf("dm[1] length %d \n",w_str->dm[1].n);
	printf("dm[2] length %d \n",w_str->dm[2].n);

	stat=snf_open_w(w_str);
	prevfieldini=0;

	for(i=0;i<ndm;i++){
		flen=0;
		snf_write_field(w_str,&prevfieldini);

		if(dm[i].type==2)snf_write_dvect(w_str->file,w_str->binasc,"double","nosparse",
			dm[i].m,dm[i].x,logx,spvect);
		snf_write_vect(w_str->file,w_str->binasc,"float","nosparse",
			dm[i].n*dm[i].m*(dm[i].complex+1),dm[i].y,logx,spvect);
	} 

	snf_close(w_str);

	return 1;
}


int snf_write_ds(SNF_STRUCT *w_str, float **dat, long rec, double tim,
				 struct LogX_STRUCT *logx)
/* writes one data record in a snf ds file; the file header should be already written.
   The record and field header variables should be set before this call, except the 
   mandatory ones (record beginning time and the record number).
   
   In case of N channels (multiple files), one should do a single call,
   with all the data for all the channels set serially in y (first all 
   the data of the first channel, then...).
*/
{
	long i,j;
	unsigned long flen,ini,fin;
	struct SpVec_STRUCT *spvec;

	//w_str->headv->rhi[0].ivalue=rec;
	w_str->headv->rhd[0].dvalue=tim;

	for(i=0;i<w_str->headv->lrheadi;i++)
		fwrite(&w_str->headv->rhi[i].ivalue,4,1,w_str->file);
	for(i=0;i<w_str->headv->lrheadd;i++)
		fwrite(&w_str->headv->rhd[i].dvalue,8,1,w_str->file);
	for(i=0;i<w_str->headv->lrheads;i++)
		fwrite(w_str->headv->rhs[i].svalue,1,16,w_str->file); printf(" @@@ nSNF3 = %d\n",w_str->nSNF3);

//OBSOLETE	for(j=0;j<w_str->ds->nch;j++){
	for(j=0;j<w_str->nSNF3;j++){
		flen=0;
		ini=ftell(w_str->file);
		fwrite(&flen,4,1,w_str->file);

		for(i=0;i<w_str->headv->lfheadi;i++)
			fwrite(&w_str->headv->fhi[i].ivalue,4,1,w_str->file);
		for(i=0;i<w_str->headv->lfheadd;i++)
			fwrite(&w_str->headv->fhd[i].dvalue,8,1,w_str->file);
		for(i=0;i<w_str->headv->lfheads;i++)
			fwrite(w_str->headv->fhs[i].svalue,1,16,w_str->file);

		snf_write_vect(w_str->file,w_str->binasc, w_str->field[j].nform,
			"nosform",w_str->field[j].len,dat[j],logx,spvec);

		fin=ftell(w_str->file);
		flen=fin-ini;
		fseek(w_str->file,ini,SEEK_SET);
		fwrite(&flen,4,1,w_str->file);
		fseek(w_str->file,fin,SEEK_SET);
	}
	
	return 1;
}


int snf_write_tfm_sfdb(SNF_STRUCT *snf_str, struct LogX_struct *logx)
/* Writes data in a snf sfdb file; the header should be already written.
   The record and field header variables should be set before this call.
   
   The data logx->vect are entered as full or half ffts, one per call (one fft is 
   a record). logx should be fully allocated.
   A record can be also composed by (not adjacent) pieces of the same fft,
   resulting in a multiple file.
*/
{

	
	return 1;
}


int snf_write_tfm_spec(SNF_STRUCT *snf_str, struct LogX_STRUCT *logx)
/* writes data of entire time-frequency map in a snf tfm file; 
   the header should be already written.
   The record and field header variables should be set before this call.

   The data (periodograms), in str->vect, are archived with LogX format.
   
   The data are entered as half periodograms, one per call 
   (one fft is a record).

   A record can be also composed by (not adjacent) pieces of the same fft,
   resulting in a multiple file.
*/
{

	
	return 1;
}



int snf_write_tfm_peak(SNF_STRUCT *snf_str, struct SpVec_STRUCT *spvect)
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
{

	
	return 1;
}


int snf_write_hm(SNF_STRUCT *snf_str, double* y)
/* writes data in a snf hm file; the header should be already written.
   The record and field header variables should be set before this call.
*/
{

	
	return 1;
}


/* ----------------- Open --------------------*/




int test_nidg(SNF_STRUCT *r_str){
	int test=0,i;
	long status;

	union U{
		struct	SNF_NIDG nidg;
		unsigned char	C[128];
	} u;

	for(i=0;i<128;i++)u.C[i]=0;

	//status=fread(&u.nidg,sizeof(struct SNF_NIDG),1,r_str->file);
	status=fread(&u.C,1,128,r_str->file);
	if(r_str->access->verb>1)for(i=0;i<100;i++)printf(" %X",u.C[i]);printf("\n  status  %d \n",status);
	
	for(i=0;i<6;i++)
		if(u.nidg.C[i]!=NIDG0.C[i]){
			if(r_str->access->verb>0)
				printf(" *** C %d %d instead of %d\n",i,u.nidg.C[i],NIDG0.C[i]);
			test=1;
		}
	for(i=0;i<4;i++)
		if(u.nidg.I[i]!=NIDG0.I[i]){
			if(r_str->access->verb>0)
				printf(" *** I %d %d instead of %d\n",i,u.nidg.I[i],NIDG0.I[i]);
			test=1;
		}
	for(i=0;i<4;i++)
		if(u.nidg.L[i]!=NIDG0.L[i]){
			if(r_str->access->verb>0)
				printf(" *** L %d %d instead of %d\n",i,u.nidg.L[i],NIDG0.L[i]);
			test=1;
		}
	for(i=0;i<6;i++)
		if(u.nidg.F[i]!=NIDG0.F[i]){
			if(r_str->access->verb>0)
				printf(" *** F %d %g instead of %g\n",i,u.nidg.F[i],NIDG0.F[i]);
			test=1;
		}
	for(i=0;i<6;i++)
		if(u.nidg.D[i]!=NIDG0.D[i]){
			if(r_str->access->verb>0)
				printf(" *** D %d %g instead of %18.12g\n",i,u.nidg.D[i],NIDG0.D[i]);
			test=1;
		}
	return test;
}


SNF_STRUCT *snf_open_r(char *file, int verb){
//int snf_open_r(SNF_STRUCT *r_str){
	FILE *fil;
	int status,recfilhead=0,snf4ok;
	char *A,*B,label[8],*token;
	long startSNF2,startSNF3,startSNF4,startSNF6,startNIDG,nsk=0;
	long i,lf0,lf1;
	SNF_STRUCT *r_str;
	STR_TOK tok;
	float f1,f2,f3,f4; /* used for a bug of MS sscanf */

	r_str=(SNF_STRUCT *)malloc(sizeof(SNF_STRUCT));
	r_str->fileName=file;
	r_str->datlen=0;
	r_str->reclen=0;
	r_str->nfields=0;

	r_str->headv=(SNF_HEADV *)malloc(sizeof(SNF_HEADV));
	r_str->headv->lheadi=0;
	r_str->headv->lheadd=0;
	r_str->headv->lheads=0;
	r_str->headv->lrheadi=0;
	r_str->headv->lrheadd=0;
	r_str->headv->lrheads=0;
	r_str->headv->lfheadi=0;
	r_str->headv->lfheadd=0;
	r_str->headv->lfheads=0;
	
	r_str->access=(SNF_ACCESS *)malloc(sizeof(SNF_ACCESS));
	r_str->access->verb=verb;
	r_str->access->nrec=0;
	r_str->access->nrecextra=0;
	r_str->access->inirec=0;


	A=(char *)calloc(128,sizeof(char));
	B=(char *)calloc(128,sizeof(char));


	printf("Opening file %s\n",r_str->fileName);

	fil=fopen(r_str->fileName,"rb");
	r_str->file=fil;
	r_str->headerlen=0;

/* #SNF#1# */

	clear_string(A,128);
	status=fread(A,1,128,fil);
	r_str->headerlen+=128;

	sscanf(A,"%7s",label);

	if(!strcmp(label,"#SNF#1#")){
		if(r_str->access->verb>0)
			printf("%s -> %s \n",label,A+7);
		nsk++;

		token=strtok(A+7,"|");
		strcpy(r_str->protocol,token);
		if(r_str->access->verb>0)printf("SNF protocol ->%s\n",r_str->protocol);

		token=strtok(NULL,"|");
		strcpy(r_str->obj,token);
		if(r_str->access->verb>0)printf("SNF object ->%s\n",r_str->obj);

		token=strtok(NULL,"|");
		strcpy(r_str->binasc,token);
		if(r_str->access->verb>0)printf("SNF file type ->%s\n",r_str->binasc);

		token=strtok(NULL,"|");
		strcpy(r_str->directaccess,token);
		if(r_str->access->verb>0)printf("SNF access ->%s\n",r_str->directaccess);
		
		token=strtok(NULL,"|");
		strcpy(r_str->userprot,token);
		if(r_str->access->verb>0)printf("SNF user protocol ->%s\n",r_str->userprot);
	}
	else{
		if(r_str->access->verb>0)
			printf(" *** Not SNF file : -> %s",A);
		return 0;
	}

/* #SNF#2# */

	startSNF2=r_str->headerlen;
	clear_string(A,128);
	status=fread(A,1,128,fil);
	r_str->nSNF2=0;

	sscanf(A,"%7s",label);

	while(!strcmp(label,"#SNF#2#")&!feof(fil)){
		if(r_str->access->verb>0)
			printf("%s -> %s \n",label,A+7);

		r_str->nSNF2++;

		clear_string(A,128);
		status=fread(A,1,128,fil);
		r_str->headerlen+=128;
		startSNF3=r_str->headerlen;

		sscanf(A,"%7s",label);
	}

	if(r_str->access->verb>0)printf(" >>> %d #SNF#2# records found\n\n",r_str->nSNF2);

	nsk+=r_str->nSNF2;

	fseek(fil,startSNF2,SEEK_SET);
	
	for(i=0;i<r_str->nSNF2;i++){
		clear_string(A,128);
		status=fread(A,1,128,fil);

		r_str->caption[i]=(char *)malloc(122*sizeof(char));
		strncpy(r_str->caption[i],A+7,121);
		printf("%s\n",r_str->caption[i]);
	}

/* #SNF#3# */

	r_str->nSNF3=0;

	clear_string(A,128);
	status=fread(A,1,128,fil);

	while(!strcmp(label,"#SNF#3#")&!feof(fil)){
		if(r_str->access->verb>0)
			printf("%s -> %s \n",label,A+7);

		r_str->nSNF3++;

		startSNF4=ftell(fil);
		clear_string(A,128);
		status=fread(A,1,128,fil); //printf("%s\n",A);
		r_str->headerlen+=128;

		sscanf(A,"%7s",label);
	}
	
	if(r_str->access->verb>0)printf(" >>> %d #SNF#3# records found\n\n",r_str->nSNF3);
	nsk+=r_str->nSNF3;

	fseek(fil,startSNF3,SEEK_SET);

	r_str->field=(SNF_DS *)calloc(r_str->nSNF3,sizeof(SNF_FIELD));

	if(!strcmp(r_str->obj,"gd")){
		r_str->gd=(SNF_GD *)calloc(r_str->nSNF3,sizeof(SNF_GD));

		for(i=0;i<r_str->nSNF3;i++){
			r_str->field[i].del=0;
			r_str->field[i].fieldlen=0;
			clear_string(A,128);
			status=fread(A,1,128,fil); //printf("@@@ A : %s\n",A);
			sscanf(A+8,"%16s ",r_str->gd[i].name);
			tok=str2numtok(A+24,5); 
			printf("### %d \n",tok.n);
			r_str->gd[i].n=atoi(tok.str[0]);
			r_str->gd[i].type=atoi(tok.str[1]);
			r_str->gd[i].complex=atoi(tok.str[2]);
			r_str->gd[i].ini=atof(tok.str[3]);
			r_str->gd[i].dx=atof(tok.str[4]);

	//		sscanf(A+8,"%16s %d %d %d %g %g ",r_str->gd[i].name,
	//			&r_str->gd[i].n,&r_str->gd[i].type,&r_str->gd[i].complex,
	//			&f1,&f2); 
	//		r_str->gd[i].ini=f1;
	//		r_str->gd[i].dx=f2;
			strcpy(r_str->field[i].name,r_str->gd[i].name);
			r_str->field[i].len=r_str->gd[i].n;
			r_str->field[i].dt=r_str->gd[i].dx;
		}
	}
	else if(!strcmp(r_str->obj,"dm")){
		r_str->dm=(SNF_DM *)calloc(r_str->nSNF3,sizeof(SNF_DM));

		for(i=0;i<r_str->nSNF3;i++){	
			r_str->field[i].del=0;
			r_str->field[i].fieldlen=0;
			clear_string(A,128);
			status=fread(A,1,128,fil);
			sscanf(A+8,"%16s ",r_str->dm[i].name);
			tok=str2numtok(A+24,8); 
			printf("### %d \n",tok.n);
			r_str->dm[i].n=atoi(tok.str[0]);
			r_str->dm[i].type=atoi(tok.str[1]);
			r_str->dm[i].complex=atoi(tok.str[2]);
			r_str->dm[i].ini=atof(tok.str[3]);
			r_str->dm[i].dx=atof(tok.str[4]);
			r_str->dm[i].m=atoi(tok.str[5]);
			r_str->dm[i].ini2=atof(tok.str[6]);
			r_str->dm[i].dx2=atof(tok.str[7]);
	//		sscanf(A+8,"%16s %d %d %d %g %g ",r_str->dm[i].name,
	//			&r_str->dm[i].n,&r_str->dm[i].type,&r_str->dm[i].complex,
	//			&f1,&f2,&r_str->dm[i].m,&f3,&f4);
	//		strcpy(r_str->field[i].name,r_str->dm[i].name);
			r_str->field[i].len=r_str->dm[i].n;
			r_str->field[i].dt=r_str->dm[i].dx;
		}
	}
	else if(!strcmp(r_str->obj,"ds")){
		for(i=0;i<r_str->nSNF3;i++){	
			r_str->field[i].del=0;
			r_str->field[i].fieldlen=0;
			clear_string(A,128);
			status=fread(A,1,128,fil);
			sscanf(A+8,"%16s %8s",&r_str->field[i].name,&r_str->field[i].nform);
			tok=str2numtok(A+32,2);
			r_str->field[i].len=atoi(tok.str[0]);
			r_str->field[i].dt=atof(tok.str[1]);
		}
	}
	else if(r_str->obj=="tfm"){
	}
	else if(r_str->obj=="tf_sfdb"){
	}
	else if(r_str->obj=="tf_spec"){
	}
	else if(r_str->obj=="tf_peak"){
	}
	else if(r_str->obj=="hm"){
	}
	else if(r_str->obj=="psc"){
	}
	else if(r_str->obj=="ev"){
	}

/* #SNF#4# and #SNF#5# */

	fseek(fil,startSNF4,SEEK_SET);
	clear_string(A,128);
	startSNF6=ftell(fil);
	status=fread(A,1,128,fil); 
	snf4ok=0;

	if(!strcmp(label,"#SNF#4#")){
		if(r_str->access->verb>0)
		printf("%s -> %s \n",label,A+7);
		snf4ok=1;

		nsk++;

		sscanf(A+7," %d %d %d %d %d %d %d %d %d",
			&r_str->headv->lheadi,&r_str->headv->lheadd,&r_str->headv->lheads,
			&r_str->headv->lrheadi,&r_str->headv->lrheadd,&r_str->headv->lrheads,
			&r_str->headv->lfheadi,&r_str->headv->lfheadd,&r_str->headv->lfheads);

		if(r_str->headv->lheadi>0)
			r_str->headv->hi=calloc(r_str->headv->lheadi,
				sizeof(SNF5));
		if(r_str->headv->lheadd>0)
			r_str->headv->hd=calloc(r_str->headv->lheadd,
				sizeof(SNF5));
		if(r_str->headv->lheads>0)
			r_str->headv->hs=calloc(r_str->headv->lheads,
				sizeof(SNF5));
		if(r_str->headv->lrheadi>0)
			r_str->headv->rhi=calloc(r_str->headv->lrheadi,
				sizeof(SNF5));
		if(r_str->headv->lrheadd>0)
			r_str->headv->rhd=calloc(r_str->headv->lrheadd,
				sizeof(SNF5));
		if(r_str->headv->lrheads>0)
			r_str->headv->rhs=calloc(r_str->headv->lrheads,
				sizeof(SNF5));
		if(r_str->headv->lfheadi>0)
			r_str->headv->fhi=calloc(r_str->headv->lfheadi,
				sizeof(SNF5));
		if(r_str->headv->lfheadd>0)
			r_str->headv->fhd=calloc(r_str->headv->lfheadd,
				sizeof(SNF5));
		if(r_str->headv->lfheads>0)
			r_str->headv->fhs=calloc(r_str->headv->lfheads,
				sizeof(SNF5));

		for(i=0;i<r_str->headv->lheadi;i++){
			clear_string(A,128);
			status=fread(A,1,128,fil);
			sscanf(A+13,"%16s %20d",r_str->headv->hi[i].name,
				&r_str->headv->hi[i].ivalue);
			B=str_trim(A+51);
			strcpy(r_str->headv->hi[i].caption,B);
			r_str->headerlen+=128;

			if(r_str->access->verb>0)
			printf("%s -> %s \n",label,A+7);
		}
		for(i=0;i<r_str->headv->lheadd;i++){
			clear_string(A,128);
			status=fread(A,1,128,fil);
			sscanf(A+13,"%16s %20g",r_str->headv->hd[i].name,
				&r_str->headv->hd[i].dvalue);
			B=str_trim(A+51);
			strcpy(r_str->headv->hd[i].caption,B);
			r_str->headerlen+=128;

			if(r_str->access->verb>0)
			printf("%s -> %s \n",label,A+7);
		}
		for(i=0;i<r_str->headv->lheads;i++){
			clear_string(A,128);
			status=fread(A,1,128,fil);
			sscanf(A+13,"%16s %16s",r_str->headv->hs[i].name,
				&r_str->headv->hs[i].svalue);
			B=str_trim(A+47);
			strcpy(r_str->headv->hs[i].caption,B);
			r_str->headerlen+=128;

			if(r_str->access->verb>0)
			printf("%s -> %s \n",label,A+7);
		}
		for(i=0;i<r_str->headv->lrheadi;i++){
			clear_string(A,128);
			status=fread(A,1,128,fil);
			sscanf(A+13,"%16s %20d",r_str->headv->rhi[i].name,
				&r_str->headv->rhi[i].ivalue);
			B=str_trim(A+51);
			strcpy(r_str->headv->rhi[i].caption,B);
			r_str->headerlen+=128;

			if(r_str->access->verb>0)
			printf("%s -> %s \n",label,A+7);
		}
		for(i=0;i<r_str->headv->lrheadd;i++){
			clear_string(A,128);
			status=fread(A,1,128,fil);
			sscanf(A+13,"%16s %20g",r_str->headv->rhd[i].name,
				&r_str->headv->rhd[i].dvalue);
			B=str_trim(A+51);
			strcpy(r_str->headv->rhd[i].caption,B);
			r_str->headerlen+=128;

			if(r_str->access->verb>0)
			printf("%s -> %s \n",label,A+7);
		}
		for(i=0;i<r_str->headv->lrheads;i++){
			clear_string(A,128);
			status=fread(A,1,128,fil);
			sscanf(A+13,"%16s %16s",r_str->headv->rhs[i].name,
				&r_str->headv->rhs[i].svalue);
			B=str_trim(A+47);
			strcpy(r_str->headv->rhs[i].caption,B);
			r_str->headerlen+=128;

			if(r_str->access->verb>0)
			printf("%s -> %s \n",label,A+7);
		}
		for(i=0;i<r_str->headv->lfheadi;i++){
			clear_string(A,128);
			status=fread(A,1,128,fil);
			sscanf(A+13,"%16s %20d",r_str->headv->fhi[i].name,
				&r_str->headv->fhi[i].ivalue);
			B=str_trim(A+51);
			strcpy(r_str->headv->fhi[i].caption,B);
			r_str->headerlen+=128;

			if(r_str->access->verb>0)
			printf("%s -> %s \n",label,A+7);
		}
		for(i=0;i<r_str->headv->lfheadd;i++){
			clear_string(A,128);
			status=fread(A,1,128,fil);
			sscanf(A+13,"%16s %20g",r_str->headv->fhd[i].name,
				&r_str->headv->fhd[i].dvalue);
			B=str_trim(A+51);
			strcpy(r_str->headv->fhd[i].caption,B);
			r_str->headerlen+=128;

			if(r_str->access->verb>0)
			printf("%s -> %s \n",label,A+7);
		}
		for(i=0;i<r_str->headv->lfheads;i++){
			clear_string(A,128);
			status=fread(A,1,128,fil);
			sscanf(A+13,"%16s %16s",r_str->headv->fhs[i].name,
				&r_str->headv->fhs[i].svalue);
			B=str_trim(A+47);
			strcpy(r_str->headv->fhs[i].caption,B);
			r_str->headerlen+=128;

			if(r_str->access->verb>0)
			printf("%s -> %s \n",label,A+7);
		}
	}

	if(snf4ok>0)
		r_str->nSNF5=r_str->headv->lheadi+r_str->headv->lheadd+r_str->headv->lheads+
			r_str->headv->lrheadi+r_str->headv->lrheadd+r_str->headv->lrheads+
			r_str->headv->lfheadi+r_str->headv->lfheadd+r_str->headv->lfheads;
	else
		r_str->nSNF5=0;

	nsk+=r_str->nSNF5;


/* #SNF#6# */

	if(snf4ok>0)startSNF6=ftell(fil);
	startNIDG=startSNF6;
	clear_string(A,128);
	status=fread(A,1,128,fil);

	sscanf(A,"%7s",label);

	r_str->nSNF6=0;

	while(!strcmp(label,"#SNF#6#")&!feof(fil)){
		if(r_str->access->verb>0)
			printf("%s -> %s \n",label,A+7);

		r_str->nSNF6++; 

		clear_string(A,128);
		status=fread(A,1,128,fil);
		r_str->headerlen+=128;

		strcpy(label,"     ");
		sscanf(A,"%7s",label);
	}
	
	if(r_str->access->verb>0)
		printf(" >>> %d #SNF#6# records found\n\n",r_str->nSNF6);

	fseek(fil,startSNF6,SEEK_SET);
	for(i=0;i<r_str->nSNF6;i++){
		clear_string(A,128);
		status=fread(A,1,128,fil);
		r_str->strSNF6[i]=malloc(128);
		strcpy(r_str->strSNF6[i],A);
	}


	nsk+=r_str->nSNF6;

	startNIDG=startSNF6+r_str->nSNF6*128;
	fseek(r_str->file,startNIDG,SEEK_SET);

	if(!strcmp(r_str->binasc,"bin")){
		test_nidg(r_str);
		r_str->headerlen+=128;
		fseek(r_str->file,startNIDG+128,SEEK_SET);
		nsk++;
	}
	
	printf(" %d header records; pointing at %d byte \n",nsk,nsk*128);

/* computes access parameters */ 

	r_str->access->inirec=nsk*128;

	r_str->nfields=r_str->nSNF3;

	r_str->reclen=-1;
	if(!strcmp(r_str->directaccess,"directaccess")){
		r_str->reclen=16+r_str->headv->lrheadd*8+r_str->headv->lrheadi*4+
			r_str->headv->lrheads*16; //printf("@@@ r_str->reclen %d \n",r_str->reclen);
		r_str->nfields=r_str->nSNF3;

		for(i=0;i<r_str->nfields;i++){
			lf0=0;//printf("@@@0 %d r_str->field[i].len %d \n",i,r_str->field[0].len);
			lf1=4;
			r_str->field[i].fieldlen=r_str->headv->lfheadd*8+r_str->headv->lfheadi*4+
				r_str->headv->lfheads*16;//printf("@@@ r_str->field[i].fieldlen %d \n",r_str->field[i].fieldlen);
			if(!strcmp(r_str->field[i].nform,"adc8")){
				lf0=4;
				lf1=1;
			}
			else if(!strcmp(r_str->field[i].nform,"adc16")){
				lf0=4;
				lf1=2;
			}
			else if(!strcmp(r_str->field[i].nform,"lin8")){
				lf0=22;
				lf1=1;
			}
			else if(!strcmp(r_str->field[i].nform,"lin16")){
				lf0=22;
				lf1=2;
			}
			else if(!strcmp(r_str->field[i].nform,"log8")){
				lf0=22;
				lf1=1;
			}
			else if(!strcmp(r_str->field[i].nform,"log16")){
				lf0=22;
				lf1=2;
			}
			else if(!strcmp(r_str->field[i].nform,"int8")){
				lf0=0;
				lf1=1;
			}
			else if(!strcmp(r_str->field[i].nform,"int16")){
				lf0=0;
				lf1=2;
			}
			else if(!strcmp(r_str->field[i].nform,"int32")){
				lf0=0;
				lf1=4;
			}
			else if(!strcmp(r_str->field[i].nform,"double")){
				lf0=0;
				lf1=8;
			}
			else{
				printf(" *** no direct access possible ! \n");
				strcpy(r_str->directaccess,"nodirectaccess");
			}


			r_str->field[i].fieldlen=4+r_str->field[i].fieldlen+lf0+
				r_str->field[i].len*lf1;//printf("@@@ %d r_str->field[i].len %d \n",i,r_str->field[i].len);

			r_str->reclen=r_str->reclen+r_str->field[i].fieldlen;
		}

		r_str->access->selbias=malloc(r_str->nfields*sizeof(long));
		r_str->access->selbias[0]=0;
		for(i=1;i<r_str->nfields;i++){
			r_str->access->selbias[i]=r_str->access->selbias[i-1]+
				r_str->field[i-1].fieldlen;
		}
	}
//printf("@@@ r_str->reclen %d \n",r_str->reclen);
	return r_str;
}


int snf_analyze(char *fileName, short verb, char *typ){
/*	verb	verbosity
	typ		type of particular analysis:
				"time"		time consistency
*/
	SNF_STRUCT *str;
	long startdata;
	int stat;

	str=snf_open_r(fileName, verb);
	startdata=ftell(str->file);

	printf("*_____________________________________________________________________________\n|");
	printf("\n|         snf_analyze utility on file %s \n|\n",fileName);

	snf_show_header(str);

	printf("*_____________________________________________________________________________\n\n");

	fseek(str->file,startdata,SEEK_SET);

//	if(str->reclen!=0){
	if(strcmp(str->obj,"gd")*strcmp(str->obj,"dm")){
		printf("\n\n                        Record analysis \n\n");
		stat=snf_show_recs(str,verb);
		
		fseek(str->file,startdata,SEEK_SET);

		printf("\n\n                        Field analysis \n\n");
		stat=snf_show_recfields(str,verb);
	}
	else{
		printf("\n\n                      Data fields analysis \n\n");
		stat=snf_show_fields(str,verb);
	}

	snf_close(str);
	printf("\n\n");
	
	return 1;
}


int snf_show_header(SNF_STRUCT *str){
	int i;

	printf("|   protocol %s  -  %s file  -  object: %s  -  %s\n|\n",
		str->protocol,str->binasc,str->obj,str->directaccess);
	printf("|              user protocol: %s \n|\n",str->userprot);

	printf("|    %d comment lines: \n",str->nSNF2);
	for(i=0;i<str->nSNF2;i++)
	    printf("| %00d * %s \n",i+1,str->caption[i]);

	printf("|\n|    %d data fields: \n",str->nSNF3);
	for(i=0;i<str->nSNF3;i++){
		if(!strcmp(str->obj,"gd"))
			printf("| gd %d: %s %d %d %d %g %g\n",i+1,str->gd[i].name,str->gd[i].n,
			 str->gd[i].type,str->gd[i].complex,str->gd[i].ini,str->gd[i].dx);
		else if(!strcmp(str->obj,"dm"))
			printf("| dm %d: %s %d %d %d %g %g %d %g %g\n",i+1,str->dm[i].name,
			 str->dm[i].n,str->dm[i].type,str->dm[i].complex,str->dm[i].ini,str->dm[i].dx,
			 str->dm[i].m,str->dm[i].ini2,str->dm[i].dx2);
		else if(!strcmp(str->obj,"ds"))
			printf("| ds %d: %s %s %d %g\n",i+1,str->field[i].name,
			 str->field[i].nform,str->field[i].len,str->field[i].dt);
	}

	printf("|\n|  i: %d  d: %d  s: %d  file header variables \n",str->headv->lheadi,
		str->headv->lheadd,str->headv->lheads);
	for(i=0;i<str->headv->lheadi;i++)
		printf("| I%4d: %16s %d  ->%s\n",i+1,str->headv->hi[i].name,str->headv->hi[i].ivalue,
		 str->headv->hi[i].caption);
	for(i=0;i<str->headv->lheadd;i++)
		printf("| D%4d: %16s %g  ->%s\n",i+1,str->headv->hd[i].name,str->headv->hd[i].dvalue,
		 str->headv->hd[i].caption);
	for(i=0;i<str->headv->lheads;i++)
		printf("| S%4d: %16s %16s  ->%s\n",i+1,str->headv->hs[i].name,str->headv->hs[i].svalue,
		 str->headv->hs[i].caption);

	printf("|\n|  i: %d  d: %d  s: %d  record header variables \n",str->headv->lrheadi,
		str->headv->lrheadd,str->headv->lrheads);
	for(i=0;i<str->headv->lrheadi;i++)
		printf("| I%4d: %16s %d  ->%s\n",i+1,str->headv->rhi[i].name,str->headv->rhi[i].ivalue,
		 str->headv->rhi[i].caption);
	for(i=0;i<str->headv->lrheadd;i++)
		printf("| D%4d: %16s %g  ->%s\n",i+1,str->headv->rhd[i].name,str->headv->rhd[i].dvalue,
		 str->headv->rhd[i].caption);
	for(i=0;i<str->headv->lrheads;i++)
		printf("| S%4d: %16s %s  ->%s\n",i+1,str->headv->rhs[i].name,str->headv->rhs[i].svalue,
		 str->headv->rhs[i].caption);

	printf("|\n|  i: %d  d: %d  s: %d  field header variables \n",str->headv->lfheadi,
		str->headv->lfheadd,str->headv->lfheads);
	for(i=0;i<str->headv->lfheadi;i++)
		printf("| I%4d: %16s %d  ->%s\n",i+1,str->headv->fhi[i].name,str->headv->fhi[i].ivalue,
		 str->headv->fhi[i].caption);
	for(i=0;i<str->headv->lfheadd;i++)
		printf("| D%4d: %16s %g  ->%s\n",i+1,str->headv->fhd[i].name,str->headv->fhd[i].dvalue,
		 str->headv->fhd[i].caption);
	for(i=0;i<str->headv->lfheads;i++)
		printf("| S%4d: %16s %16s  ->%s\n",i+1,str->headv->fhs[i].name,str->headv->fhs[i].svalue,
		 str->headv->fhs[i].caption);

	printf("|\n|    %d user defined header lines: \n",str->nSNF6);
	for(i=0;i<str->nSNF6;i++)
		printf("|%4d: %s\n",i+1,str->strSNF6[i]);
	
	printf("|\n|    header length: %d    record length: %d   fields number %d\n",
		str->headerlen,str->reclen,str->nfields);
	for(i=0;i<str->nfields;i++)
		printf("|     field %4d :  length %d   bias %d\n",i+1,str->field[i].fieldlen,
		 str->access->selbias[i]);

	return 1;
}


int snf_show_recs(SNF_STRUCT *str, int verb){
	int i=0,iord=0,iext=0,nextrec;

	while(!feof(str->file)){
		nextrec=snf_read_rec(str);

		if(nextrec>0){
			iord++;
			i++;
			printf("      ordinary record %d  expected length: %d - real length: %d\n",
				iord,str->reclen,str->reclenread);
			fseek(str->file,nextrec,SEEK_SET);
		}
		else if(nextrec==-1){
			iord++;
			i++;
			printf("      ordinary record %d  expected length: %d - real length: %d\n",
				iord,str->reclen,str->reclenread);
			return 1;
		}
		else if(nextrec==-2){
			iext++;
			i++;
			printf("      extra-ordinary record %d  expected length: %d - real length: %d\n",
				iord,str->reclen,str->reclenread);
			return 1;
		}
		else if(nextrec<-2){
			iext++;
			i++;
			printf("      extra-ordinary record %d  expected length: %d - real length: %d\n",
				iord,str->reclen,str->reclenread);
			fseek(str->file,-nextrec,SEEK_SET);
		}
		else{
			printf(" *** error in reading record %d\n",i);
			return 0;
		}
	}

	return 1;
}


int snf_show_recfields(SNF_STRUCT *str, int verb){
	int i=0,iord=0,iext=0,nextrec;
	int j,nextfield;

	while(!feof(str->file)){
		nextrec=snf_read_rec(str);

		if(nextrec>0){
			iord++;
			i++;
			printf("\n rec: %d  ord: %d  fields:",i,iord);
			for(j=0;j<str->nfields;j++){
				nextfield=snf_read_field(str);
				printf(" %d ",j);
				fseek(str->file,nextfield,SEEK_SET);
			}
			printf(" -> nextrec, nextfield: %d, %d",nextrec,nextfield);

			fseek(str->file,nextrec,SEEK_SET);
		}
		else if(nextrec==-1){
			iord++;
			i++;
			printf("\n rec: %d  ord: %d  fields:",i,iord);
			for(j=0;j<str->nfields;j++){
				nextfield=snf_read_field(str);
				printf(" %d ",j);
				fseek(str->file,nextfield,SEEK_SET);
			}
			printf(" -> nextrec, nextfield: %d, %d",nextrec,nextfield);

			return 1;
		}
		else if(nextrec==-2){
			iext++;
			i++;

			return 1;
		}
		else if(nextrec<-2){
			iext++;
			i++;

			fseek(str->file,-nextrec,SEEK_SET);
		}
		else{
			printf(" *** error in reading record %d\n",i);
			return 0;
		}
	}

	return 1;
}


int snf_show_fields(SNF_STRUCT *str, int verb){
	int j,nextfield;

	j=ftell(str->file);
	printf("   data fields start at %d \n",j);

	for(j=0;j<str->nfields;j++){
		nextfield=snf_read_field(str);
		printf("   field %d  next field at %d\n",j,nextfield);
		fseek(str->file,nextfield,SEEK_SET);
	}

	fseek(str->file,0,SEEK_END);
	j=ftell(str->file);
	printf("   data fields end at %d\n\n",j);

	return 1;
}


/* ----------------- Read --------------------*/


int snf_read_rec(SNF_STRUCT *str)
/* 
	return variable :
	  0		error
	 >0		position of next record, ordinary record
	 <0		-position of next record, extraordinary record
	 -1		last record

  if ASCII, return 1
*/
{
	char label[9]=" ";
	long i,ini;

	if(!strcmp(str->binasc,"ascii")){
		return 1;
	}

	label[8]=0;
	ini=ftell(str->file);
	fread(label,1,8,str->file);
	if(!strcmp(label,"#SNF#ORD")){
		fread(&str->access->nrec,4,1,str->file);
		fread(&str->reclenread,4,1,str->file);
		if(str->reclenread>0)
			ini=ini+str->reclenread;
		else ini=-1; /* last record is ord */
		
		for(i=0;i<str->headv->lrheadi;i++)
			fread(&str->headv->rhi[i].ivalue,4,1,str->file);
		for(i=0;i<str->headv->lrheadd;i++)
			fread(&str->headv->rhd[i].dvalue,8,1,str->file);
		for(i=0;i<str->headv->lrheads;i++)
			fread(&str->headv->rhs[i].svalue,1,16,str->file);

		return ini;
	}
	else if(!strcmp(label,"#SNF#EXT")){
		fread(&str->access->nrec,4,1,str->file);
		fread(&str->reclenread,4,1,str->file);
		if(str->reclenread>0)
			ini=-(ini+str->reclenread);
		else ini=-2; /* last record is ext */
	}
	else return 0;
}


int snf_read_field(SNF_STRUCT *str)
/* 
	return variable: 
	
	  the position of next field

	 -1		last field of the file

  if ASCII, return 1
*/
{
	long i,ini,flen;

	if(!strcmp(str->binasc,"ascii")){
		return 1;
	}

	ini=ftell(str->file);

	fread(&flen,4,1,str->file); //printf("@@@ ini,flen %d,%d\n",ini,flen);
	if(flen>0)
			ini=ini+flen;
	else ini=-1;
		
	for(i=0;i<str->headv->lfheadi;i++)
		fread(&str->headv->fhi[i].ivalue,4,1,str->file);
	for(i=0;i<str->headv->lfheadd;i++)
		fread(&str->headv->fhd[i].dvalue,8,1,str->file);
	for(i=0;i<str->headv->lfheads;i++)
		fread(&str->headv->fhs[i].svalue,1,16,str->file);

	return ini;
}


int snf_getrec(SNF_STRUCT *str){
	int stat=0;

	if(strcmp(str->directaccess,"directaccess")){
		str->access->inirec=str->access->inirec+str->reclen;
	}
	else{
	}

	return stat;
}



float *snf_read_vect(FILE* file,char* binasc, char* nform, char* sform, long n, 
				   struct LogX_STRUCT *logx, struct SpVec_STRUCT *spvect)
/*
	logx, spvect	externally allocated only if necessary (depending on nform and sform)
*/
{
	float *vout;
	long nwrit,i;

	vout=(float *)malloc(n*sizeof(float));
	//printf(" ***DBG a %d \n",n);
	
	if(!strcmp(binasc,"ascii")){
		for(i=0;i<n;i++)
			fscanf(file," %g ",vout+i);
	}
	else{
		if(!strcmp(nform,"float")){
			nwrit=fread(vout,sizeof(float),n,file);
	//printf(" ***DBG b %d \n",ftell(str->file));
		}
	}

	return vout;
}


double *snf_read_dvect(FILE* file, char* binasc, char* nform, char* wform, long n, 
				   struct LogX_STRUCT *logx, struct SpVec_STRUCT *spvect)
/*
	logx, spvect	externally allocated only if necessary (depending on nform and sform)
*/
{
	double *vout;
	long nwrit;

	vout=(double *)malloc(n*sizeof(double));
	
	if(!strcmp(nform,"double")){
		nwrit=fread(vout,sizeof(double),n,file);
	}

	return vout;
}


GD *snf_read_gd(char *file, long kgd, int verb)
/*
	file	file to open
	kgd		take the k-th (1,2,...,n) gd of the file (if multiple file) 
	verb	verbosity
*/
{
	GD *gd;
	SNF_STRUCT *r_str;
	long kgd1,njump,i;
	struct LogX_STRUCT *logx;
	struct SpVec_STRUCT *spvect;
	unsigned long flen;

	r_str=snf_open_r(file, verb);
	//printf(" ***DBG n %d \n",r_str->gd[0].n);

	if(strcmp(r_str->obj,"gd")){
		printf(" *** %s file doesn't contain a GD ! \n",file);
		return NULL;
	}

	kgd1=1;
	if(r_str->nSNF3>1)kgd1=kgd;
	if(kgd1<1)kgd1=1;
	kgd1--;
	//printf(" ***DBG kgd1 %d \n",kgd1);

	gd=(GD *)malloc(sizeof(GD));

	gd->capt=r_str->caption[kgd1];
	gd->complex=r_str->gd[kgd1].complex;
	gd->cont=0;
	gd->dx=r_str->gd[kgd1].dx;
	gd->ini=r_str->gd[kgd1].ini;
	gd->n=r_str->gd[kgd1].n;
	gd->type=r_str->gd[kgd1].type;

	if(!strcmp(r_str->binasc,"ascii")&&r_str->nSNF3>1){
		printf(" *** Incompatible attributes: ASCII and multiple files\n");
		return gd;
	}

//	njump=0;
	for(i=0;i<kgd1-1;i++){
		fread(&flen,4,1,r_str->file);
		fseek(r_str->file,flen-4,SEEK_CUR);
//		njump=njump+r_str->gd[i].n*sizeof(float)+4;
//		if(r_str->gd[i].type==2)njump=njump+r_str->gd[i].n*sizeof(double);
	}
	fseek(r_str->file,-flen,SEEK_CUR);

	if(gd->type==2){
		gd->x=(double *)malloc(gd->n*sizeof(double));
		gd->x=snf_read_dvect(r_str->file,r_str->binasc,"double","nosparse",gd->n,logx,spvect);
	}
	
	//printf(" ***DBG %d \n",ftell(r_str->file));
	gd->y=(float *)malloc(gd->n*sizeof(float));
	gd->y=snf_read_vect(r_str->file,r_str->binasc,"float","nosparse",gd->n,logx,spvect);
	//printf(" ***DBG %d \n",ftell(r_str->file));

	if(verb>0)printf(" %d data read",gd->n);

	snf_close(r_str);
	
	return gd;
}


int snf2ds(DS *ds, short kds, struct pss_infile_db *dbf,
			 RING *r, GD *g, int verb)
/*
	ds		the attached ds (in case of multiple ds, the element of the array of ds)
	kds		sequence number in case of multiple ds (the first is 0), 0 in case of single ds
	dbf		the structure with the infile data-base
	r,g		service structures, should be defined outside with all_ringxds and all_gdxds
	verb	verbosity level
*/			 
{
	long len,len2,len4;
	long totin,totout;
	long i,nv,nc1,nc2;
	int status;
	double lastim,tinit;
//	float *v;
//	char *lcw;

	len=ds->len;
	len2=len/2;
	len4=len/4;

	if (ds->lcw == 0)
	{
		ds->cont=0;
		ds->capt="frames by pss_frame";
	}

	totin=r->totin;
	totout=r->totout;

	if (ds->debug > 0) printf(" --> len,totin,totout= %d,%d,%d\n",len,totin,totout);

	nc1=ds->nc1;
	nc2=ds->nc2;

	printf(" --> nc1,nc2= %d,%d\n",nc1,nc2);

	if (ds->type == 2)  /* case interlaced */
	{
		if (nc1 <= nc2) /* beginning and/or odd chunks */
		{
			if (nc1 == 0)
			{
				while (totin-totout < len)  /* if there are few ring data, get more */
				{
					ds->cont++;
//					status=get_1ch(dbf->infile,adcName,g,verb);
					nv=g->n;

					status=write_rg(r,g->y,nv,lastim); /* ATTENZIONE LASTIM */
					totin=r->totin;
					ds->dt=g->dx;
					r->dx=g->dx;
				}

				status=read_rg(r,g->y,len,&tinit);

				//printf(" g->y = %12.4e \n",g->y[0]);

				for (i=0; i < len4; i++)    ds->y1[i]=0.;
				for (i=len4; i < 5*len/4; i++)  ds->y1[i]=g->y[i-len4];

				//
				//
			}

			ds->lcw++;
			ds->nc1++;

			printf("ds chunk -> %d",ds->lcw);
		}
		else    /* even chunks :  produces d.y2 */
		{
			while (totin-totout < len)
			{
				ds->cont++;
//				status=get_1ch(dbf->infile,adcName,g,verb);
				nv=g->n;

				status=write_rg(r,g->y,nv,lastim);
				totin=r->totin;
				ds->dt=g->dx;
				r->dx=g->dx;
			}

			status=read_rg(r,g->y,len,&tinit);

			//printf(" g->y = %12.4e \n",g->y[0]);

			for (i=0; i < 3*len4; i++)
				ds->y2[i]=ds->y1[i+len2];
			for (i=3*len4; i < len; i++)
				ds->y2[i]=g->y[i-3*len4];
			for (i=0; i < len4; i++)
				ds->y1[i]=ds->y2[i+len2];
			for (i=len4; i < 5*len4; i++)
				ds->y1[i]=g->y[i-len4];

			ds->lcw++;
			ds->nc2++;

			printf("ds chunk -> %d",ds->lcw);
		}
	}
	else
	{
		while (totin-totout < len)
		{
			ds->cont++;
//			status=get_1ch(dbf->infile,adcName,g,verb);
			nv=g->n;

			status=write_rg(r,g->y,nv,lastim);
			totin=r->totin;
			ds->dt=g->dx;
			r->dx=g->dx;
		}

		if (ds->type == 1)
			for(i=0; i<ds->len; i++)
				ds->y2[i]=ds->y1[i];

		status=read_rg(r,g->y,len,&tinit);

		//printf(" g->y = %12.4e \n",g->y[0]);

		for (i=0; i < len; i++)
			ds->y1[i]=g->y[i];

		ds->lcw++;
		ds->nc2++;

		printf("ds chunk -> %d",ds->lcw);
	}

	return 1;
}



int snf_close(SNF_STRUCT *snf_str){
	fclose(snf_str->file);

	return 1;
}


/*---------------------------- Templates -------------------------------*/


SNF_STRUCT template_write_ds(){
	SNF_STRUCT *w_str;
	LogX_STRUCT *logx;
	SpVec_STRUCT *spvec;

	char chname[3*20];
	long chlen[3]={100,200,300};
	double chdt[3]={0.3,0.2,0.1};
	char* numform[3]={"float",
						"float",
						"float"
	};

	long i,j,j0,k;
	long precfieldini=0,rec=0;
	float *v1;

	strcpy(chname,"ch1");
	strcpy(chname+20,"ch2");
	strcpy(chname+40,"ch3");
	
	v1=(float *)malloc(1000*sizeof(float));

	w_str=snf_all_ds(chname,chlen,chdt,3,"ds_templ.snf",
		"Template 3 channels",numform);

	snf_open_w(w_str);
	j0=0;

	for(i=0;i<10;i++){
		snf_write_rec(w_str,&rec);
		for(j=0;j<w_str->nSNF3;j++){
			snf_write_field(w_str,&precfieldini);
			if(j==0){
				for(k=0;k<w_str->field[j].len;k++)v1[k]=j0+k;
				j0+=w_str->field[j].len;
			}
			if(j==1)for(k=0;k<w_str->field[j].len;k++)v1[k]=k;
			if(j==2)for(k=0;k<w_str->field[j].len;k++)v1[k]=k*k;

			snf_write_vect(w_str->file,w_str->binasc,w_str->field[j].nform,
				w_str->field[j].sform,w_str->field[j].len*(1+w_str->field[j].complex),
				v1,logx,spvec);
		}
	}

	snf_close(w_str);

	return *w_str;
}


SNF_STRUCT template_read_ds(){
	SNF_STRUCT r_str;

	return r_str;
}


SNF_STRUCT template_write_tfm_sfdb(){
	SNF_STRUCT *w_str;
	LogX_STRUCT *logx;
	SpVec_STRUCT *spvec;

	long blen[2]={512,256};
	long inif[2]={0,300};
	double df=0.01;
	double v_ini[3]={0.1,0.2,0.3};		/* ! */
	double v_mid[3]={0.1,0.2,0.3};		/* ! */
	double v_fin[3]={0.1,0.2,0.3};		/* ! */
	double v_mean[3]={0.1,0.2,0.3};		/* ! */
	double t_mjd,t_gps;					/* ! */

	long i,j,k;
	long precfieldini=0,rec=0;
	float *v1,*v2;
	const long lv1=1024;
	const long lv2=512;
	
	v1=(float *)malloc(lv1*sizeof(float));	
	v2=(float *)malloc(lv2*sizeof(float));

	/* call something to update the arguments of the following function */

	w_str=snf_open_tfm("tf_sfdb_templ.snf","Template tf_sfdb","float",
		"sfdb_mar02",0.001,"sfdb",4096,2048,"no",2,blen,inif,df);			/* !!! */

	for(i=0;i<10;i++){		/* produces 10 records */
		/* call something to update the ts and the vs */
		w_str->headv->rhi[0].ivalue=0;

		w_str->headv->rhd[0].dvalue=t_mjd;
		w_str->headv->rhd[1].dvalue=t_gps;

		for(j=0;j<3;j++)w_str->headv->rhd[2+j].dvalue=v_ini[j];
		for(j=0;j<3;j++)w_str->headv->rhd[5+j].dvalue=v_mid[j];
		for(j=0;j<3;j++)w_str->headv->rhd[8+j].dvalue=v_fin[j];
		for(j=0;j<3;j++)w_str->headv->rhd[11+j].dvalue=v_mean[j];

		snf_write_rec(w_str,&rec);		/* !!! */	

		/* call something to update the vectors (now v1 and v2) */

		/* the following lines are unnecessary */
		for(j=0;j<lv1;j++)v1[j]=k++;
		for(j=0;j<lv2;j++)v2[j]=k++;

		j=0;
		snf_write_field(w_str,&precfieldini);		/* !!! */
		snf_write_vect(w_str->file,w_str->binasc,w_str->field[j].nform,
			w_str->field[j].sform,w_str->field[j].len*(1+w_str->field[j].complex),
			v1,logx,spvec);

		j=1;
		snf_write_field(w_str,&precfieldini);	/* !!! */
		snf_write_vect(w_str->file,w_str->binasc,w_str->field[j].nform,
			w_str->field[j].sform,w_str->field[j].len*(1+w_str->field[j].complex),
			v2,logx,spvec);
	}

	snf_close(w_str); // printf(" @@@ file %s\n",w_str->fileName);

	return *w_str;
}


SNF_STRUCT template_read_tfm_sfdb(){
	SNF_STRUCT *r_str;

	return *r_str;
}


SNF_STRUCT template_write_tfm_spec(){
	SNF_STRUCT *w_str;

	return *w_str;
}


SNF_STRUCT template_read_tfm_spec(){
	SNF_STRUCT *r_str;

	return *r_str;
}


SNF_STRUCT template_write_tfm_peak(){
	SNF_STRUCT *w_str;

	return *w_str;
}


SNF_STRUCT template_read_tfm_peak(){
	SNF_STRUCT *r_str;

	return *r_str;
}
