/*___________________________________________________
 ¦                                                   ¦
 ¦                   pss_frame.c                     ¦
 ¦       by Sergio Frasca - ((( 0 ))) Virgo          ¦
 ¦                    March 2000                     ¦
 ¦___________________________________________________¦*/



#include <stdio.h>
#include <string.h>

#include "FrameL.h"
//Include externally the last framelib
#include "pss_snag.h"
#include "pss_serv.h"  //pia: invertito l' ordine fra pss_serv.h e pss_frame.h: dave problemi la
                       //struttura pss_infile_db, definita in pss_serv.h e usata in pss_frame.h
#include "pss_frame.h"


struct FrFile* open_fr(char* fileName)
{
	struct FrFile* iFile;

	printf(" %s > \n",fileName);

	iFile = FrFileINew(fileName);
	if(iFile == NULL){
		printf("Error during file opening\n"
			"  Last errors are:\n%s",FrErrorGetHistory());
	}
	return iFile;
};


void inquire_fr_file(char* fileName)
{
	struct FrFile *iFile;
	FrameH *frame;
	FrAdcData *adc;
	FrSerData *sms;
	//FrProcData *proc; //pia tolto non lo usa
	FrSimEvent *simE;
	//FrTrigData *trig; 
	FrSimData *sim;
	FrSummary *sum;

	double mjd, gpssec, secini, secfin;
	long NumFr=0;

	iFile=open_fr(fileName);

	if(iFile != NULL){
		frame=FrameRead(iFile);

		if(frame!=NULL){ 
			secini=frame->GTimeS+frame->GTimeN*0.000000001;
		}

		while(frame != NULL){
			mjd=mjd_frame(frame);
			gpssec=frame->GTimeS;
			gpssec+=frame->GTimeN*1.e-9;
			NumFr+=1;
			printf("\nFrame %d/%d GPSsec=%13.3f mjd=%13.7f - %s\n",frame->run,frame->frame,
				gpssec,mjd,asctime(GPS2UTC(frame->GTimeS,frame->ULeapS)));
			if(frame->rawData != NULL){
				for(adc = frame->rawData->firstAdc; adc != NULL; adc=adc->next)
				        printf("\nAdc: %s : %s  : dt= %f,  n= %lu",
					       adc->name,adc->comment,adc->sampleRate,adc->data->nData); //pia %lu
				for(sms = frame->rawData->firstSer; sms != NULL; sms = sms->next)
					printf("\nSer: %s : dt= %f ",sms->name,sms->sampleRate);
			}
			for(sim = frame->simData; sim != NULL; sim = sim->next)
			                 printf("\nSim: %s : %s  : %f,  n= %lu",sim->name,sim->comment,sim->sampleRate,sim->data->nData); //pia %lu
			for(simE = frame->simEvent; simE != NULL; simE = simE->next)
					printf("\nSimE: %s : %s  ",simE->name,simE->comment);
			//		for(trig = frame->trigData; trig != NULL; trig = trig->next)
			//		printf("\nTrig: %s : %s  ",trig->name,trig->comment);
			 //pia: tolto il for e il printf perche' proc->samplerate non esiste
			//for(proc = frame->procData; proc != NULL; proc = proc->next)
			     //	  printf("\nProc: %s : %s  : %f,  n= %d",proc->name,proc->comment,proc->sampleRate,proc->data->nData);
			for(sum = frame->summaryData; sum != NULL; sum = sum->next)
					printf("\nSum: %s : %s  ",sum->name,sum->comment);
 
			secfin=frame->GTimeS+frame->GTimeN*0.000000001;
			frame=FrameRead(iFile);
		}
		if(NumFr>1)printf("\n\n   %ld frames with mean duration %f s\n",NumFr,(secfin-secini)/(NumFr-1)); //pia ld
	}
	FrFileIEnd(iFile);
}


struct Frame_ch inquire_ch(char* fileName)
{
	struct FrFile *iFile;
	FrameH *frame;
	FrAdcData *adc;
	long ich=0,nch,i;

	struct Frame_ch chs;

	iFile=open_fr(fileName);

	if(iFile != NULL){
		frame=FrameRead(iFile);
		if(frame->rawData != NULL)
			for(adc = frame->rawData->firstAdc; adc != NULL; adc=adc->next){
				ich+=1;
			}
	}
	
	FrFileIEnd(iFile);

	nch=ich;
	chs.nch=nch;
	for(i=0;i<1000;i++)chs.name[i]=0;
	for(i=0;i<2000;i++)chs.comment[i]=0;
	ich=0;
	
	iFile=open_fr(fileName);

	if(iFile != NULL){
		frame=FrameRead(iFile);
		if(frame->rawData != NULL){
			chs.frameN=frame->frame;
			chs.tin=mjd_frame(frame);

			for(adc = frame->rawData->firstAdc; adc != NULL; adc=adc->next){
				vec_strinp(ich, 20, adc->name, chs.name);
				if(adc->comment != NULL)vec_strinp(ich, 40, adc->comment, chs.comment);
				chs.samp[ich]=1./adc->sampleRate;
				chs.len[ich]=adc->data->nData;
				ich+=1;
			}
		}
	}
	
	FrFileIEnd(iFile);

	return chs;
}


struct FrAdcData* inquire_adc(char* fileName, char *adcName, int verb)
{
/*	verb = (0,1,...) verbosity
*/
	struct FrAdcData *adc;
	struct FrFile* iFile;

	struct FrameH *frame;

	iFile = open_fr(fileName);

	frame = FrameRead(iFile);
	if(frame == NULL )printf(" *** Error reading file %s ", fileName);
	
	adc = FrAdcDataFind(frame,adcName);

	if (verb > 0)
    {
		printf(" -> Inquire file %s adc %s \n",fileName,adcName);
		printf(" -> Frame %d run %d ",frame->frame, frame->run);
		printf(" -> Starting GPS time:%s (and %d usec)\n",
			FrStrGTime(frame->GTimeS),frame->GTimeN/1000); 
		printf(" -> type,sample rate,nData,nBytes : %d, %f, %lu, %lu \n\n",
		       adc->data->type,adc->sampleRate,adc->data->nData,adc->data->nBytes); //pia %lu %lu
	}

	FrFileIEnd(iFile);

	return adc;
}



int get_1ch(struct FrFile *iFile, char *adcName,
			GD *Dat, int verb)
/*	Dat must be created
	Dat.y must be allocated for a dimension at least equal to
		the value of adc.data.nData
	verb = (0,1,...) verbosity
*/
{
	struct FrameH *frame;
	struct FrAdcData *adc;
	struct FrVect *vect;

	long i;
	int status,nData;

	frame = FrameRead(iFile);
	if(frame == NULL ){status=1; return status;}
	if (verb > 0)
    {
		printf(" Frame %d run %d ",frame->frame, frame->run);
		printf(" Starting GPS time:%s (and %d usec)\n",
	       FrStrGTime(frame->GTimeS),frame->GTimeN/1000); 
	}

	adc = FrAdcDataFind(frame,adcName);
	if(adc == NULL)
		{printf("No %s channel",adcName);
		status=2;
		return status;}

	nData= adc->data->nData;
	if (nData > Dat->nall)
	{
		status=3;
		printf(" *** Dat gd must be allocated for at least %d ; now is %ld\n",
		       nData,Dat->nall); //pia ld
	}
	Dat->n = nData;
	//printf(" > %d \n",nData);
	//data=(float *)malloc((size_t) nData*sizeof(float));
	Dat->dx = 1./adc->sampleRate;
	Dat->ini = mjd_frame(frame);

	vect = adc->data;
    if(vect == NULL) 
		{printf(" no data for frame %u\n",frame->frame);}  //pia ud
    else if(vect->type == FR_VECT_C)
		{for(i=0; i<nData; i++) {Dat->y[i] = vect->data[i];}}
    else if(vect->type == FR_VECT_2S)
        {for(i=0; i<nData; i++) {Dat->y[i] = vect->dataS[i];}}
    else if(vect->type == FR_VECT_4S)
        {for(i=0; i<nData; i++) {Dat->y[i] = vect->dataI[i];}}
    else if(vect->type == FR_VECT_1U)
        {for(i=0; i<nData; i++) {Dat->y[i] = vect->dataU[i];}}
    else if(vect->type == FR_VECT_2U)
        {for(i=0; i<nData; i++) {Dat->y[i] = vect->dataUS[i];}}
    else if(vect->type == FR_VECT_4U)
        {for(i=0; i<nData; i++) {Dat->y[i] = vect->dataUI[i];}}
    else if(vect->type == FR_VECT_4R)
        {for(i=0; i<nData; i++) {Dat->y[i] = vect->dataF[i];}}
    else if(vect->type == FR_VECT_8R)
        {for(i=0; i<nData; i++) {Dat->y[i] = vect->dataD[i];}}

	FrameFree(frame);
	status=0;
	return status;
}



int get_1chA(struct FrFile *iFile, char *adcName,
			float **data1, long* nData)
{
	struct FrameH *frame;
	struct FrAdcData *adc;
	struct FrVect *vect;
	float *data;

	int status=0;
	double dt;
	long i;

	data=*data1;
	frame = FrameRead(iFile);
	if(frame != NULL)
    {
		printf(" Frame %d run %d ",frame->frame, frame->run);
		printf(" Starting GPS time:%s (and %d usec)\n",
	       FrStrGTime(frame->GTimeS),frame->GTimeN/1000); 
	}
	else {return 1;}

	adc = FrAdcDataFind(frame,adcName);
	if(adc == NULL)
		{printf("No %s channel",adcName);
		return 2;}

	*nData = adc->data->nData;
	printf(" > %ld \n",*nData); //pia ld
	free(data);
	data=(float*)malloc((size_t) nData*sizeof(float));
	dt = 1./adc->sampleRate;

	vect = adc->data;
    if(vect == NULL) 
         {printf(" no data for frame %u\n",frame->frame);} //pia u
    else if(vect->type == FR_VECT_C)
		{for(i=0; i<*nData; i++) {data[i] = vect->data[i];}}
    else if(vect->type == FR_VECT_2S)
        {for(i=0; i<*nData; i++) {data[i] = vect->dataS[i];}}
    else if(vect->type == FR_VECT_4S)
        {for(i=0; i<*nData; i++) {data[i] = vect->dataI[i];}}
    else if(vect->type == FR_VECT_1U)
        {for(i=0; i<*nData; i++) {data[i] = vect->dataU[i];}}
    else if(vect->type == FR_VECT_2U)
        {for(i=0; i<*nData; i++) {data[i] = vect->dataUS[i];}}
    else if(vect->type == FR_VECT_4U)
        {for(i=0; i<*nData; i++) {data[i] = vect->dataUI[i];}}
    else if(vect->type == FR_VECT_4R)
        {for(i=0; i<*nData; i++) {data[i] = vect->dataF[i];}}
    else if(vect->type == FR_VECT_8R)
        {for(i=0; i<*nData; i++) {data[i] = vect->dataD[i];}}

	FrameFree(frame);
	*data1=data;
	return status;
}


int get_1ch_vect(char *fileName, char *chName,
			float *dat, long kframe, long ndat, long lenfrvect, int verb)
/*	dat must be created before
	dat.y must be allocated for a dimension at least equal to
		the value of ndat
	kframe is the sequence number of the frame in the file, starting from 1
	verb = (0,1,...) verbosity
*/
{
	struct FrFile *iFile;
	struct FrameH *frame;
	struct FrVect *vect;

	long i,ifr,ii;
	int status;
	float *dat1;
	double mjd;

	dat1=calloc(lenfrvect,sizeof(float));

	iFile=open_fr(fileName);
	ifr=0;
	ii=0;

	if(iFile != NULL){
		frame=FrameRead(iFile);

		while(frame != NULL){
			mjd=mjd_frame(frame);
			ifr+=1;
			if(ifr >= kframe){
				if(verb >0) printf("\n Frame %d/%d GTimeS=%d  mjd=%f\n",frame->run,frame->frame,frame->GTimeS,mjd);
//				vect=FrameGetV(frame,chName); //pia sostituita con quella sotto
				vect=FrameFindVect(frame,chName);

				if(vect == NULL) 
				       printf(" no data for frame %u\n",frame->frame); //pia %u
				else if(vect->type == FR_VECT_C)
					for(i=0; i<lenfrvect; i++) dat1[i] = vect->data[i];
				else if(vect->type == FR_VECT_2S)
					for(i=0; i<lenfrvect; i++) dat1[i] = vect->dataS[i];
				else if(vect->type == FR_VECT_4S)
					for(i=0; i<lenfrvect; i++) dat1[i] = vect->dataI[i];
				else if(vect->type == FR_VECT_1U)
					for(i=0; i<lenfrvect; i++) dat1[i] = vect->dataU[i];
				else if(vect->type == FR_VECT_2U)
					for(i=0; i<lenfrvect; i++) dat1[i] = vect->dataUS[i];
				else if(vect->type == FR_VECT_4U)
					for(i=0; i<lenfrvect; i++) dat1[i] = vect->dataUI[i];
				else if(vect->type == FR_VECT_4R)
					for(i=0; i<lenfrvect; i++) dat1[i] = vect->dataF[i];
				else if(vect->type == FR_VECT_8R)
					for(i=0; i<lenfrvect; i++) dat1[i] = vect->dataD[i];
			}
			i=0;
			while(ii < ndat && i< lenfrvect){
				dat[ii]=dat1[i];
				ii+=1;
				i+=1;
			}

			if(ii<ndat){
				frame=FrameReadRecycle(iFile,frame);
			}
			else {
				FrFileIEnd(iFile);
				status=0;
				return status;
			}
		}
	}
	
	FrFileIEnd(iFile);

	if(verb>0)printf("\n Not enough data: only %ld ; zero padded",ii); //pia ld

	for(i=ii;i<ndat;i++)dat[i]=0;

	status=1;
	return status;
}



double mjd_frame(struct FrameH *fr)
/* time in mjd */
{
	double t,GS,GN;

	GS=fr->GTimeS;
	GN=fr->GTimeN*1.e-9;

	t=(GS+GN-fr->ULeapS+19)/86400+44244;
	//printf(" %d %d %20.10f\n",fr->GTimeS,fr->GTimeN,t);

	return t;
}



int frame2ds_single(DS *d, struct FrFile *iFile, char *adcName, RING *r, GD *g, int verb)
{
	long len,len2,len4;
	long totin,totout;
	long i,nv,nc1,nc2;
	int status;
	double lastim,tinit;
//	float *v;
//	char *lcw;

	len=d->len;
	len2=len/2;
	len4=len/4;

	if (d->lcw == 0)
	{
		d->cont=0;
		d->capt="frames by pss_frame";
	}

	totin=r->totin;
	totout=r->totout;

	if (d->debug > 0) printf(" --> len,totin,totout= %ld,%ld,%ld\n",len,totin,totout); //pia ld

	nc1=d->nc1;
	nc2=d->nc2;

	printf(" --> nc1,nc2= %ld,%ld\n",nc1,nc2); //pia ld

	if (d->type == 2)  /* case interlaced */
	{
		if (nc1 <= nc2) /* beginning and/or odd chunks */
		{
			if (nc1 == 0)
			{
				while (totin-totout < len)  /* if there are few ring data, get more */
				{
					d->cont++;
					status=get_1ch(iFile,adcName,g,verb);
					nv=g->n;

					status=write_rg(r,g->y,nv,lastim); /* ATTENZIONE LASTIM */
					totin=r->totin;
					d->dt=g->dx;
					r->dx=g->dx;
				}

				status=read_rg(r,g->y,len,&tinit);

				//printf(" g->y = %12.4e \n",g->y[0]);

				for (i=0; i < len4; i++)    d->y1[i]=0.;
				for (i=len4; i < 5*len/4; i++)  d->y1[i]=g->y[i-len4];

				//
				//
			}

			d->lcw++;
			d->nc1++;

			printf("ds chunk -> %ld",d->lcw); //pia ld
		}
		else    /* even chunks :  produces d.y2 */
		{
			while (totin-totout < len)
			{
				d->cont++;
				status=get_1ch(iFile,adcName,g,verb);
				nv=g->n;

				status=write_rg(r,g->y,nv,lastim);
				totin=r->totin;
				d->dt=g->dx;
				r->dx=g->dx;
			}

			status=read_rg(r,g->y,len,&tinit);

			//printf(" g->y = %12.4e \n",g->y[0]);

			for (i=0; i < 3*len4; i++)
				d->y2[i]=d->y1[i+len2];
			for (i=3*len4; i < len; i++)
				d->y2[i]=g->y[i-3*len4];
			for (i=0; i < len4; i++)
				d->y1[i]=d->y2[i+len2];
			for (i=len4; i < 5*len4; i++)
				d->y1[i]=g->y[i-len4];

			d->lcw++;
			d->nc2++;

			printf("ds chunk -> %ld",d->lcw); //pia ld
		}
	}
	else
	{
		while (totin-totout < len)
		{
			d->cont++;
			status=get_1ch(iFile,adcName,g,verb);
			nv=g->n;

			status=write_rg(r,g->y,nv,lastim);
			totin=r->totin;
			d->dt=g->dx;
			r->dx=g->dx;
		}

		if (d->type == 1)
			for(i=0; i<d->len; i++)
				d->y2[i]=d->y1[i];

		status=read_rg(r,g->y,len,&tinit);

		//printf(" g->y = %12.4e \n",g->y[0]);

		for (i=0; i < len; i++)
			d->y1[i]=g->y[i];

		d->lcw++;
		d->nc2++;

		printf("ds chunk -> %ld",d->lcw); //pia ld
	}

	return 1;
}



int frame2ds(DS *d, struct pss_infile_db *dbf, char *adcName,  RING *r, GD *g, int verb)
{
	long len,len2,len4;
	long totin,totout;
	long i,nv,nc1,nc2;
	int status;
	double lastim,tinit;
//	float *v;
//	char *lcw;

	len=d->len;
	len2=len/2;
	len4=len/4;

	if (d->lcw == 0)
	{
		d->cont=0;
		d->capt="frames by pss_frame";
	}

	totin=r->totin;
	totout=r->totout;

	if (d->debug > 0) printf(" --> len,totin,totout= %ld,%ld,%ld\n",len,totin,totout); //pia ld

	nc1=d->nc1;
	nc2=d->nc2;

	printf(" --> nc1,nc2= %ld,%ld\n",nc1,nc2); //pia ld

	if (d->type == 2)  /* case interlaced */
	{
		if (nc1 <= nc2) /* beginning and/or odd chunks */
		{
			if (nc1 == 0)
			{
				while (totin-totout < len)  /* if there are few ring data, get more */
				{
					d->cont++;
					status=get_1ch(dbf->frinfile,adcName,g,verb);
					nv=g->n;

					status=write_rg(r,g->y,nv,lastim); /* ATTENZIONE LASTIM */
					totin=r->totin;
					d->dt=g->dx;
					r->dx=g->dx;
				}

				status=read_rg(r,g->y,len,&tinit);

				//printf(" g->y = %12.4e \n",g->y[0]);

				for (i=0; i < len4; i++)    d->y1[i]=0.;
				for (i=len4; i < 5*len/4; i++)  d->y1[i]=g->y[i-len4];

				//
				//
			}

			d->lcw++;
			d->nc1++;

			printf("ds chunk -> %ld",d->lcw); //pia ld
		}
		else    /* even chunks :  produces d.y2 */
		{
			while (totin-totout < len)
			{
				d->cont++;
				status=get_1ch(dbf->frinfile,adcName,g,verb);
				nv=g->n;

				status=write_rg(r,g->y,nv,lastim);
				totin=r->totin;
				d->dt=g->dx;
				r->dx=g->dx;
			}

			status=read_rg(r,g->y,len,&tinit);

			//printf(" g->y = %12.4e \n",g->y[0]);

			for (i=0; i < 3*len4; i++)
				d->y2[i]=d->y1[i+len2];
			for (i=3*len4; i < len; i++)
				d->y2[i]=g->y[i-3*len4];
			for (i=0; i < len4; i++)
				d->y1[i]=d->y2[i+len2];
			for (i=len4; i < 5*len4; i++)
				d->y1[i]=g->y[i-len4];

			d->lcw++;
			d->nc2++;

			printf("ds chunk -> %ld",d->lcw);  //pia ld
		}
	}
	else
	{
		while (totin-totout < len)
		{
			d->cont++;
			status=get_1ch(dbf->frinfile,adcName,g,verb);
			nv=g->n;

			status=write_rg(r,g->y,nv,lastim);
			totin=r->totin;
			d->dt=g->dx;
			r->dx=g->dx;
		}

		if (d->type == 1)
			for(i=0; i<d->len; i++)
				d->y2[i]=d->y1[i];

		status=read_rg(r,g->y,len,&tinit);

		//printf(" g->y = %12.4e \n",g->y[0]);

		for (i=0; i < len; i++)
			d->y1[i]=g->y[i];

		d->lcw++;
		d->nc2++;

		printf("ds chunk -> %ld",d->lcw); //pia ld
	}

	return 1;
}



int read_fmnl_files0(struct pss_infile_db *dbf, int verb)
/* reads an input files file */
{
	FILE *files0;
	char line[80];
	long i=0;
	int linlen;

	if((files0=fopen(dbf->fil_db,"r"))==NULL)
	{
		printf(" *** Error opening file %s !\n",dbf->fil_db);
		return 3;
	}
	else
		printf(" File-DB file %s opened\n\n",dbf->fil_db);

	while(!feof(files0))
	{
		fscanf(files0,"%s\n",line);
		linlen=strlen(line);
		dbf->files[i]=(char*)malloc(linlen+1);
		strcpy(dbf->files[i],line);
		i++;
		if(verb>0)printf(" %ld %d %s \n",i,linlen,line); //pia ld
	}

	printf("\n %ld files in the db \n",i); //pia ld

	return 0;
}



int check_infile_db(struct pss_infile_db* dbf, char *outfile)
{
	FILE *files1;

	printf("NOT USED: outfile %s\n,",outfile); //pia: aggiunto il print  perche' non usati
	printf("NOT USED: dbf->fil_db %s\n,",dbf->fil_db); //pia
	files1=fopen("d:\\scra\\out.dat","W");
	return 1;
}




struct tm *GPS2UTC(long GPSsec, int ULeapS)
/*	
	GPSsec   gps seconds (integer value)
	ULeapS   leap seconds (TAI-UTC : 19 on 1 Jan 1980, 32 on 1 Jan 1999)

*/
{
	struct tm *t;
	int mdays[12]={31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
	long secutc,ndays,nsec,ndaystot;
	long i=0,nyear=0,yeardays,mday;
	int Y,M,D,DoY,DoW,h,m,s;
	time_t now;

	secutc=GPSsec-ULeapS+19+5*86400;
	ndays=secutc/86400;
	nsec=secutc-ndays*86400;
	ndaystot=ndays;

	nyear=0;
	while(ndays > 0){
		yeardays=365;
		if((nyear/4)*4 == nyear)yeardays=366;
		ndays-=yeardays;
		nyear+=1;
	}
	nyear=nyear-1;
	ndays=ndays+yeardays;
	mday=ndays;

	i=0;
	if((nyear/4)*4 == nyear)mdays[1]=29;
	while(mday >= 0){
		mday=mday-mdays[i]; //printf(" --> %d %d %d \n",ndays,i,mday);
		i++;
	}
	mday=mday+mdays[i-1];
	// t=localtime(&secutc);  Error
	now=time(0);
	t=localtime(&now);

	Y=nyear+1980-1900;
	M=i-1;
	D=mday+1;
	DoY=ndays;
	DoW=ndaystot+2-((ndaystot+2)/7)*7;
	h=nsec/3600;
	m=(nsec-h*3600)/60;
	s=nsec-h*3600-m*60;
	//printf(" %d %d %d  %d %d  %d %d %d \n",Y,M,D,DoY,DoW,h,m,s); 

	t->tm_year=Y;
	t->tm_mon=M;
	t->tm_mday=D;

	t->tm_hour=h;
	t->tm_min=m;
	t->tm_sec=s;
	t->tm_yday=DoY;
	t->tm_wday=DoW;
	t->tm_isdst=0;

	return t;
}
