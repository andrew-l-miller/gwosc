/*___________________________________________________
 ¦                                                   ¦
 ¦                   pss_r87.c                       ¦
 ¦       by Sergio Frasca - ((( 0 ))) Virgo          ¦
 ¦                    June 2001                      ¦
 ¦___________________________________________________¦*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pss_snag.h"
#include "pss_serv.h"
#include "pss_r87.h"

r87_file open_r87(char *fileName, int verb){
	short A[30];
	int count,i;

	//_fmode=_O_BINARY;

	r87_file fil;
	fil.fil=fopen(fileName,"rb");

	if(fil.fil==NULL){
		if(verb>0)printf("Error opening %s\n",fileName);
		return fil;
	}
	else if(verb>0)printf("%s opened\n",fileName);

	count=fread(&A,sizeof(short),30,fil.fil);
	if(verb>0)printf("\n  H  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d\n",
		A[0],A[1],A[2],A[3],A[4],A[5],A[6],A[7],A[8],A[9],A[10],A[11],A[12],A[13],A[14],A[15]);

	fil.reclen=A[0];
	for(i=0;i<6;i++)
		fil.initime[i]=A[12+i];
	fil.samptim=(A[18]+A[19]*0.001)*0.001;

	fclose(fil.fil);

	fil.fil=fopen(fileName,"rb");

	fil.A=A;

	return fil;
}


r87_header read_header_r87(short *A){
	r87_header header;
	int it[7],i;
	double time;

	header.len=A[0];
	header.headlen=A[1];
	header.recnum=A[2];
	header.antenna=A[3];
	header.runnum=A[4];
	header.type=A[5];

	header.len1=A[6];
	header.nc1=A[7];
	header.len2=A[8];
	header.nc2=A[9];
	header.len3=A[10];
	header.nc3=A[11];

	it[0]=A[12];
	it[1]=1;
	for(i=2;i<6;i++)it[i]=A[11+i];

	it[6]=A[17]*1000;
	if(it[0]>1998)it[6]=A[17]*100;


	time=mjd_r87(it);
	header.time=time;

	header.st=(A[18]+A[19]*0.001)*0.001;

	header.adccode=A[21];
	header.opflag=A[24];

	return header;
}


r87_inq_ch r87_inquire_ch(char *fileName,int verb){
	r87_file fil;
	r87_header head;
	r87_inq_ch inq;
	long i,ii=0;

	fil=open_r87(fileName,verb);
	head=read_header_r87(fil.A);

	for(i=0;i<head.nc1;i++){
		inq.num[i]=101+i;
		inq.samp[i]=fil.samptim;
		inq.tin=head.time;
		inq.len[i]=head.len1/head.nc1;
		sprintf(inq.name+ii*20,"%d",inq.num[i]);
		ii++;
	}
	for(;i<head.nc1+head.nc2;i++){
		inq.num[i]=201+i-head.nc1;
		inq.samp[i]=(fil.samptim*head.len1*head.nc2)/(head.len2*head.nc1);
		inq.tin=head.time;
		inq.len[i]=head.len2/head.nc2;
		sprintf(inq.name+ii*20,"%d",inq.num[i]);
		ii++;
	}
	for(;i<head.nc1+head.nc2+head.nc3;i++){
		inq.num[i]=301+i-head.nc1-head.nc2;
		inq.samp[i]=(fil.samptim*head.nc3*head.len1)/(head.len3*head.nc1);
		inq.tin=head.time;
		inq.len[i]=head.len3/head.nc3;
		sprintf(inq.name+ii*20,"%d",inq.num[i]);
		ii++;
	}
	fclose(fil.fil);

	inq.nch=ii;

	return inq;
}


r87_chdata read_r87rec_ch(FILE *fil,long reclen,long ch){
	long chfield,chan,count,i,ii;
	long len,nc,bias;
	double stcor;
	short *A;
	r87_chdata chdata;

	chfield=ch/100;
	chan=ch-chfield*100;

	A=(short*)calloc(reclen,sizeof(short));

	ii=0;
	fseek(fil,ii*reclen*2,0);
	count=fread(A,sizeof(short),reclen,fil);ii+=1;

	while(A[5] != 1){fseek(fil,ii*reclen*2,0);
	count=fread(A,sizeof(short),reclen,fil);ii+=1;}

	chdata.header=read_header_r87(A);
	
	switch(chfield){
	case 1:
	   {len=chdata.header.len1;nc=chdata.header.nc1;bias=chdata.header.headlen+chan;break;}
	case 2:
		{len=chdata.header.len2;nc=chdata.header.nc2;bias=chdata.header.headlen+
			chdata.header.len1+chan;break;}
	case 3:
		{len=chdata.header.len3;nc=chdata.header.nc3;bias=chdata.header.headlen+
			chdata.header.len1+chdata.header.len2+chan;break;}
	}

	stcor=(double)(nc*chdata.header.len1)/(chdata.header.nc1*len);
	chdata.header.tsamp=chdata.header.st*stcor;
	chdata.header.field=chfield;
	chdata.header.chan=chan;
	chdata.header.len=len;
	chdata.header.nc=nc;
	chdata.header.ndata=len/nc;
	
	chdata.data=calloc(len/nc,sizeof(float));
	
	for(i=0;i<len/nc;i++){
		if(chdata.header.adccode==0)chdata.data[i]=(A[bias]-2024)*0.004883;
		else chdata.data[i]=A[bias]*0.004883;
		bias=bias+nc;
	}
	return chdata;
}


r87_chdata get_r87_1ch_vect(char *fileName, long chNum,
					 long krec, long ndat, int verb){
	r87_file rfil;
	long chfield,chan,count,i,ii,jj=0;
	long reclen,recndat,len,nc,bias,bias0;
	double stcor;
	short *A;
	r87_chdata chdata;

	//FILE *fil;
	//fil=fopen("dump.dat","w");

	
	chdata.data=(float *)calloc(ndat,sizeof(float));

	rfil=open_r87(fileName,verb);
	reclen=rfil.reclen;

	chfield=chNum/100;
	chan=chNum-chfield*100;

	A=(short*)calloc(reclen,sizeof(short));

	ii=0;
	fseek(rfil.fil,ii*reclen*2,0);
	count=fread(A,sizeof(short),reclen,rfil.fil);ii++;

	while(A[5] != 1){fseek(rfil.fil,ii*reclen*2,0);
		count=fread(A,sizeof(short),reclen,rfil.fil);ii++;
	}

	chdata.header=read_header_r87(A);
	
	switch(chfield){
	case 1:{
			len=chdata.header.len1;
			nc=chdata.header.nc1;
			bias=chdata.header.headlen+chan-1;
			break;
		}
	case 2:{
			len=chdata.header.len2;
			nc=chdata.header.nc2;
			bias=chdata.header.headlen+chdata.header.len1+chan-1;
			break;
		}
	case 3:{
			len=chdata.header.len3;
			nc=chdata.header.nc3;
			bias=chdata.header.headlen+chdata.header.len1+chdata.header.len2+chan-1;
			break;
		}
	}
	bias0=bias;
	recndat=len/nc;

	stcor=(double)(nc*chdata.header.len1)/(chdata.header.nc1*len);
	chdata.header.tsamp=chdata.header.st*stcor;
	chdata.header.field=chfield;
	chdata.header.chan=chan;
	chdata.header.len=len;
	chdata.header.nc=nc;
	chdata.header.ndata=recndat;

	jj=0;
	//printf("\n I-> %d %d %d %d %d \n",chdata.header.headlen,chdata.header.len1,chdata.header.len2,chdata.header.len3,chan);
	//printf("\n O-> %d %d %d %d \n",len,nc,bias,recndat);

	while(jj<ndat){
		bias=bias0;
		if(chdata.header.adccode==0)for(i=0;i<recndat && jj<ndat;i++){
			chdata.data[jj++]=(A[bias]-2024)*0.004883;
			bias=bias+nc; 
			//fprintf(fil," \t%d \t%d  \n",jj,A[bias]);
		}
		else for(i=0;i<recndat && jj<ndat;i++){
			chdata.data[jj++]=A[bias]*0.004883;
			bias=bias+nc; 
			//fprintf(fil," \t%d \t%d  \n",jj,A[bias]);
		}
		//printf(" rec %d  bias %d\n",A[2],bias);
		do{
			fseek(rfil.fil,ii*reclen*2,0);
			count=fread(A,sizeof(short),reclen,rfil.fil);
			ii++;
		}while(A[5] != 1);
	}

	return chdata;
}


char *read_r87_info(char *fileName, long *numchar, int verb){
	FILE *fil;
	long ii=0,count,reclen,headlen,ninfo=0,i,j=0, k=0, tot=0,rectot=0;
	long appo=0;
	short A[10];
	char *B, *C;

	fil=fopen(fileName,"rb");

	if(fil==NULL){
		if(verb>0)printf("Error opening %s\n",fileName);
		return "error";
	}

	count=fread(A,sizeof(short),9,fil);
	reclen=A[0];
	headlen=A[1];

	while(!feof(fil)){
	//for(i=0;i<10;i++){
		tot=ii*reclen*2;
		fseek(fil,tot,0);
		count=fread(A,2,9,fil);
		ii++;
		if(A[5]!=1)ninfo++;
	}

	rectot=ii;

	if(verb>0)printf("\n %d info records found\n",ninfo);

	ii=0;j=0;
	B=(char *)calloc(reclen*2,sizeof(short));
	C=(char *)calloc(reclen*ninfo*2,sizeof(short));

	fseek(fil,0,0);
	for(k=0;k<rectot;k++){
		tot=ii*reclen*2;
		fseek(fil,tot,0);
		count=fread(B,sizeof(char),reclen*2,fil);
		ii++;
		if(B[10]!=1)for(i=headlen*2;i<reclen*2;i++){
			if(B[i]!=0){
				C[j]=B[i];
				j++;
			}
			else appo++;
		}
	}
	*numchar=j;
	fclose(fil);if(verb>0)printf(" %d non-null character read; null %d \n",j,appo);

	return C;
}


double mjd_r87(long *it){
	/* it[0]   year
	   it[1]   month
	   it[2]   day
	   it[3]   hour
	   it[4]   minute
	   it[5]   second
	   it[6]   microsecond
	*/

	return 0;
}


void r87_dump_rec(char *fileName, long nrec, long inirec){
	r87_file rfil;
	FILE *fil;
	long ii=0,reclen,i,count;
	short *A;

	fil=fopen("dump.dat","w");

	rfil=open_r87(fileName,1);
	reclen=rfil.reclen;
	
	A=(short*)calloc(reclen,sizeof(short));

	for(ii=inirec-1;ii<inirec+nrec&&!feof(rfil.fil);ii++){
		fseek(rfil.fil,ii*reclen*2,0);
		count=fread(A,sizeof(short),reclen,rfil.fil);
		if(ferror(rfil.fil))printf("\n Error record %d",ii+1);
		if(feof(rfil.fil))printf("\n EOF fseek %d count %d",ii*reclen*2,count);
		if(count!=reclen)printf("\n  rec %d  count = %d fseek %d",ii+1,count,ii*reclen*2);
		fprintf(fil,"\n\n --------------- Record %d ------------   count %d\n",ii+1,count);
		for(i=0;i<reclen;i++)fprintf(fil,"\n \t%d   \t%d",i+1,A[i]);
	}
}