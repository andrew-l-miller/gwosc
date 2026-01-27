/*___________________________________________________
 ¦                                                   ¦
 ¦                  frame_util.c                     ¦
 ¦       by Sergio Frasca - ((( 0 ))) Virgo          ¦
 ¦                    June 2003                      ¦
 ¦___________________________________________________¦*/

/* New version March 2007 */

#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "frame_util.h"
#include "pss_serv.h"
#include "pss_snag.h"
#include "pss_frame.h"
#include "pss_sfc.h"
#include "FrameL.h"

void frame_util(){
	char *title="Frame_Util menu";
	char *items[12]={"Folder choice",
					"File choice",
					"Channel choice",
					"DataType file name block",
					"Output directory",
					"sds file creation",
					"Time check",
					"Batch mode (0,1)",
					"Antenna acronym",
					"Dump file standard (on stdumpfile.txt)",
					"Dump file (on dumpfile.txt)",
					"Exit"};
	char filename[200];
	char direc[80],direcout[80];
	char file[80];
	char channel[40];
	char datatype[10];
	char antacr[10];
	char subtype[10];
	char sdsfile[80];
	char* dumpfile="dumpfile.txt";
	char* stdumpfile="stdumpfile.txt";
	char strou[80],inmenu[200];
	int itout,batmode=0;

	strcpy(direc,"D:\\Data\\pss\\virgo\\sd\\frame\\");
	strcpy(direcout,direc);
	strcpy(file,"HrecV-806888400-01-Aug-2005-01h40-600F.gwf");
	//strcpy(channel,"Pr_B1_ACq");
	//strcpy(channel,"Pr_B1_ACq_4kHz");
	//strcpy(channel,"h_4kHz_NoLines");
	//strcpy(channel,"h_4kHz");
	//strcpy(channel,"h_20kHz");
	//strcpy(channel,"h_20kHz_NoLines");
	//strcpy(channel,"dL4kHz");
	strcpy(channel,"h_4kHzNo50");
	strcpy(datatype,"xxxx");
	strcpy(antacr,"VIR");
	strcpy(subtype,"");
	strcpy(sdsfile,"prova.sds");

	while (1)
	{
		strcpy(filename,direc);
		strcat(filename,file);printf(" |%s| |%s| |%s|\n",direc,file,filename);
		printf("Working on %s  -> channel %s\n",filename,channel);
		printf("Output file :  %s \n",sdsfile);
		printf(" ---> GPS time = 0 is %s \n",asctime(GPS2UTC(0,19)));
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
			strcpy(direcout,direc);
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
			strcpy(inmenu,"Output directory ?  - (after DirIn) now ");
			strcat(inmenu,direcout);
			strcat(inmenu," >");
			if(batmode==0)
				get_line(inmenu,strou);
			else
				get_line_bat(inmenu,strou);
			strcpy(direcout,strou);
			break;
		case 6:
			to_sds(filename,direcout,channel,antacr,datatype,subtype);
			break;
		case 7:
			inquire_fr_file(filename);
			break;
		case 8:
			printf("Now batch mode is %d ; new mode ? (0,1,2) ",batmode);
			scanf("%d",&batmode);
			break;
		case 9:
			strcpy(inmenu,"Antenna acronym ?  - now ");
			strcat(inmenu,antacr);
			if(batmode==0)
				get_line(inmenu,strou);
			else
				get_line_bat(inmenu,strou);
			strcpy(antacr,strou);
			break;
		case 10:
			print_dump_std(stdumpfile, filename);
			break;
		case 11:
			print_dump(dumpfile, filename);
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



int frame_util_bat(char *filbat)
// Batch file format:
//   folder in         ex.: D:\Data\pss\virgo\sd\frame\
//   folder out        ex.: D:\Data\pss\virgo\sd\sds\
//   channel           ex.: h_4kHz
//   antenna acronym   ex.: VIR
//   data type         ex.: hrec
//   frame file1       ex.: HrecV-806888400-01-Aug-2005-01h40-600F.gwf
//   frame file2
//   frame file3
//    ...
{
	char direc[81],direcout[81];
	char file[81];
	char filename[161];
	char channel[41];
	char datatype[10];
	char antacr[10];
	char subtype[10];
	FILE *fil;
	int ii,i=0;

	strcpy(direc,"D:\\Data\\pss\\virgo\\sd\\frame\\");
	strcpy(direcout,direc);
	strcpy(filename,"HrecV-806888400-01-Aug-2005-01h40-600F.gwf");
	//strcpy(channel,"Pr_B1_ACq");
	//strcpy(channel,"Pr_B1_ACq_4kHz");
	//strcpy(channel,"h_4kHz_NoLines");
	//strcpy(channel,"h_4kHz");
	//strcpy(channel,"h_20kHz");
	//strcpy(channel,"h_20kHz_NoLines");
	//strcpy(channel,"dL4kHz");
	strcpy(channel,"h_4kHzNo50");
	strcpy(datatype,"xxxx");
	strcpy(antacr,"VIR");
	strcpy(subtype,"");

	fil=fopen(filbat,"r");

    getlineKR(fil,direc,80); printf(" %s \n",direc);
    getlineKR(fil,direcout,80); printf(" %s \n",direcout);
    getlineKR(fil,channel,40); printf(" %s \n",channel);
    getlineKR(fil,antacr,6); printf(" %s \n",antacr);
    getlineKR(fil,datatype,9); printf(" %s \n",datatype);

	while(1){
		ii=getlineKR(fil,file,80);
		printf(" %s started \n",file);
		strcpy(filename,direc);
		strcat(filename,file);
		if (ii > 0)
			to_sds(filename,direcout,channel,antacr,datatype,subtype);
		else 
			return ii;
		i++;
		printf(" %s processed \n",filename);
	}
}


void print_dump_std(char *dumpfile, char *filename){
	char *cnow;
	FILE *dfil;

	FrFile *file;         
	FrameH *frame;
	FrAdcData *adc;
	FrSerData *sms;
	FrProcData *proc;
	FrSimEvent *simE;
	FrEvent *evt;
	FrSimData *sim;
	FrSummary *sum;

	int debug=0,i=0;

	cnow=showtim();

	dfil=fopen(dumpfile,"wt");

	fprintf(dfil,"\n    ***  %s file dump done on %s  ***\n\n",filename,cnow);

	           /*---------------------------- Open the input file(s) ------------*/

	FrLibSetLvl(debug);
	file = FrFileINew(filename);
	if(file == NULL)
		{fprintf(stderr,"Cannot open input file %s\n %s",
                    filename, FrErrorGetHistory()); 
		return;}


	FrFileIDump(file, dfil, debug, NULL);
	//FrFileIDump(file, dfil, debug, "*");

	//FrFileIRewind(file);

	frame = FrameRead(file); //printf(" indirizzo %p \n",frame);

           /*----------------Start the Main loop ----------------------------*/

	while(frame != NULL) { 
		i++; //printf(" ciao %d \n",i);
		fprintf(dfil,"\n        --> frame sequence number %d \n",i);
		fprintf(dfil,"Frame %d/%d GTimeS=%d\n",frame->run,frame->frame,frame->GTimeS);
        for(adc = frame->rawData->firstAdc; adc != NULL; adc=adc->next) 
			{FrAdcDataDump( adc,  dfil, debug);}
        
        if(frame->detectProc != NULL)
              {fprintf(dfil,"Detector used for reconstruction:\n");
               FrDetectorDump(frame->detectProc, dfil, debug);}
        if(frame->detectSim != NULL)
              {fprintf(dfil,"Detector used for simulation:\n");
               FrDetectorDump(frame->detectSim, dfil, debug);}

		if(frame->rawData != NULL)
			{for(sms = frame->rawData->firstSer; sms != NULL; sms = sms->next) 
				{fprintf(dfil," SMS: %s Data:%30s\n", sms->name, sms->data);}}
        if(frame->detectProc != NULL)
			{fprintf(dfil,"Detector used for reconstruction:\n");
            FrDetectorDump(frame->detectProc, dfil, debug);}
        if(frame->detectSim != NULL)
            {fprintf(dfil,"Detector used for simulation:\n");
            FrDetectorDump(frame->detectSim, dfil, debug);}
        for(sim = frame->simData; sim != NULL; sim = sim->next) 
            {FrSimDataDump(sim, dfil, debug);}
        for(simE = frame->simEvent; simE != NULL; simE = simE->next)
            {FrSimEventDump(simE, dfil,debug);}
        for(evt = frame->event; evt != NULL; evt = evt->next) 
            {FrEventDump(evt, dfil, debug);}
        for(proc = frame->procData; proc != NULL; proc = proc->next) 
            {fprintf(dfil,"Reconstructed Data: %s ", proc->name);
            FrVectDump( proc->data,  dfil, debug); 
            FrVectDump( proc->aux,   dfil, debug);}
        for(sum = frame->summaryData; sum != NULL; sum = sum->next)
            {fprintf(dfil,"Summary: %s %s \n",sum->name, sum->comment); 
			FrVectDump( sum->moments, dfil, debug);}
		FrRawDataDump(frame->rawData, dfil, debug);

		FrameStat(frame, dfil);

		FrameFree(frame);

           /*---------------------- read the next frame --------------------*/

		frame = FrameRead(file);

           /*------------------End the Main loop ----------------------------*/
	}

	fprintf(dfil,"\n\n           ***  FileIStat  *** \n\n");  
	FrFileIStat(file, dfil);

	FrFileIEnd(file);
	
	fclose(dfil);
}


void print_dump(char *dumpfile, char *filename){
	char *cnow;
	FILE *dfil;

	FrFile *file;         
	FrameH *frame;
	FrAdcData *adc;
	FrProcData *proc;
	FrSimEvent *simE;
	FrEvent *evt;
	FrSimData *sim;

	int debug=0,i=0;
	double tini=0,tini1,terr,dt,t0;
	struct tm *timstr;
	char *asctim;

	cnow=showtim();

	dfil=fopen(dumpfile,"wt");

	fprintf(dfil,"\n    ***  %s file dump done on %s  ***\n\n",filename,cnow);

	           /*---------------------------- Open the input file(s) ------------*/

	FrLibSetLvl(debug);
	file = FrFileINew(filename);
	if(file == NULL)
		{fprintf(stderr,"Cannot open input file %s\n %s",
                    filename, FrErrorGetHistory()); 
		return;}


	FrFileIDump(file, dfil, debug, NULL);

	frame = FrameRead(file);

	t0=mjd_frame(frame);

           /*----------------Start the Main loop ----------------------------*/

	while(frame != NULL) { 
		i++;
		if(i == 1){
			fprintf(dfil,"\n     Project %s   run: %d   mjd %f (leap s: %d)\n\n",
				frame->name,frame->run,t0,frame->ULeapS);
			tini=frame->GTimeS+frame->GTimeN*1.e-9;
			dt=frame->dt;
			if(frame->detectProc != NULL){
//				fprintf(dfil,"  Detector %s (long,lat,h %8.5f,%8.5f,%8.5f) (arms %8.5f,%8.5f)\n",
//					frame->detectProc->longitude,frame->detectProc->latitude,frame->detectProc->elevation,
//					frame->detectProc->armXazimuth,frame->detectProc->armYazimuth);
			}
			if(frame->rawData != NULL){
				fprintf(dfil,"              ADC channels\n");
				for(adc = frame->rawData->firstAdc; adc != NULL; adc=adc->next){
					fprintf(dfil," channel %s  group %d n: %d -> sampl. %12.9f Hz, len %d, offset %7.5f ns\n",
						adc->name,adc->channelGroup,adc->channelNumber,adc->sampleRate,adc->data->nData,adc->timeOffset*1.e9);
				}
			}
			if(frame->procData != NULL){
				fprintf(dfil,"              Proc channels\n");
				for(proc = frame->procData; proc != NULL; proc=proc->next){
					fprintf(dfil," channel %s  type %d subtype %d -> sampl. %12.9f Hz, len %d,offset %7.5f ns\n",
						proc->name,proc->type,proc->subType,1/proc->data->dx[0],proc->data->nData,proc->timeOffset*1.e9);
				}
			}
			if(frame->simData != NULL){
				fprintf(dfil,"              Sim channels\n");
				for(sim = frame->simData; sim != NULL; sim=sim->next){
					fprintf(dfil," channel %s  -> sampl. %12.9f Hz  len %d\n",
						sim->name,1/sim->data->dx[0],sim->data->nData);
				}
			}
		}

		dt=frame->dt;
		timstr=GPS2UTC(frame->GTimeS,frame->ULeapS);
		asctim=asctime(timstr);

		tini1=tini;
		tini=frame->GTimeS+frame->GTimeN*1.e-9;
		terr=tini-tini1-dt;
		if(terr>1.e-6){
			printf("  *** Time error frame %d --> at %s  dt = %f \n",
				frame->frame,asctim,terr);
			fprintf(dfil,"  *** Time error frame %d --> at %s  dt = %f \n",
				frame->frame,asctim,terr);
		}

		fprintf(dfil,"\n\n        --> frame sequence number %d   %s\n",i,asctim);
		fprintf(dfil,"Frame %d/%d GTimeS=%d GTimeN=%d  duration:%15.10f\n",
			frame->run,frame->frame,frame->GTimeS,frame->GTimeN,frame->dt);
		if(frame->rawData != NULL){
			for(adc = frame->rawData->firstAdc; adc != NULL; adc=adc->next) 
			{FrAdcDataDump( adc,  dfil, debug);}}
        
/*        if(frame->rawData != NULL)
			{for(sms = frame->rawData->firstSer; sms != NULL; sms = sms->next) 
				{fprintf(dfil," SMS: %s Data:%30s\n", sms->name, sms->data);}} */
        if(frame->detectProc != NULL)
			{fprintf(dfil,"Detector used for reconstruction:");
            FrDetectorDump(frame->detectProc, dfil, debug);} 
        if(frame->detectSim != NULL)
            {fprintf(dfil,"Detector used for simulation:\n");
            FrDetectorDump(frame->detectSim, dfil, debug);}
        for(sim = frame->simData; sim != NULL; sim = sim->next) 
            {FrSimDataDump(sim, dfil, debug);}
        for(simE = frame->simEvent; simE != NULL; simE = simE->next)
            {FrSimEventDump(simE, dfil,debug);}
        for(evt = frame->event; evt != NULL; evt = evt->next) 
            {FrEventDump(evt, dfil, debug);}
        for(proc = frame->procData; proc != NULL; proc = proc->next) 
            {fprintf(dfil,"\n Reconstructed Data: %s ", proc->name);
//            FrVectDump( proc->data,  dfil, debug); 
//            FrVectDump( proc->aux,   dfil, debug);
			}
/*		for(sum = frame->summaryData; sum != NULL; sum = sum->next)
            {fprintf(dfil,"Summary: %s %s \n",sum->name, sum->comment); 
		FrVectDump( sum->moments, dfil, debug);}
		FrRawDataDump(frame->rawData, dfil, debug);  */

//		FrameStat(frame, dfil);

		FrameFree(frame);

           /*---------------------- read the next frame --------------------*/

		frame = FrameRead(file);

           /*------------------End the Main loop ----------------------------*/
	}

	fprintf(dfil,"\n\n           ***  FileIStat  *** \n\n");  
	FrFileIStat(file, dfil);

	FrFileIEnd(file);
	
	fclose(dfil);
}


void to_sds(char *filename, char *direcout, char *channel, char *antacr, char *datatype, char *sub_typ){
	FILE *fil;

	FrFile *file;         
	FrameH *frame;
	FrAdcData *adc;
	FrProcData *proc;
	FrSimData *sim;

	FrVect *vec;
	float *data;
	int icdata,icstart;

	long ndat=0,ndatold,i,ii=0,iii=0,ic_opennew=0;
	long len=0,len1=0;  /* Attention now len is int 32 ! */
	double t0,t0A,t0gps,t0gpsA=0,dt,errt=0;
	long ULeapS;

	int debug=0;
	double tini;
	
	           /*---------------------------- Open the input file(s) ------------*/

	FrLibSetLvl(debug);
	file = FrFileINew(filename);
	if(file == NULL)
		{fprintf(stderr,"Cannot open input file %s\n %s",
                    filename, FrErrorGetHistory()); 
		return;}


	frame = FrameRead(file);
	while(frame == NULL & iii < 6){
		iii++;
		printf(" *** Strangely null frame ! \n");
		frame = FrameRead(file);
	}

           /*----------------Start the Main loop ----------------------------*/
	
	ic_opennew=0;
	icdata=0;
	icstart=1;
	ndatold=0;

	while(frame != NULL) {
		ii++;
		
		t0gps=frame->GTimeS+frame->GTimeN*1.e-9;
		t0=mjd_frame(frame);
		ULeapS=frame->ULeapS;

		if(ii == 1){
			printf("\n     Project %s   run: %d   (leap s: %d)\n\n",
				frame->name,frame->run,frame->ULeapS);
			tini=frame->GTimeS+frame->GTimeN*1.e-9;
			if(frame->detectProc != NULL){
//				fprintf(dfil,"  Detector %s (long,lat,h %8.5f,%8.5f,%8.5f) (arms %8.5f,%8.5f)\n",
//					frame->detectProc->longitude,frame->detectProc->latitude,frame->detectProc->elevation,
//					frame->detectProc->armXazimuth,frame->detectProc->armYazimuth);
			}
			if(frame->rawData != NULL){
				printf("              ADC channels\n");
				for(adc = frame->rawData->firstAdc; adc != NULL; adc=adc->next){
					printf(" channel %s  -> sampl. %12.9f Hz \n",
						adc->name,adc->sampleRate);
				}
			}
			if(frame->procData != NULL){
				printf("              Proc channels\n");
				for(proc = frame->procData; proc != NULL; proc=proc->next){
					printf(" channel %s  -> sampl. %12.9f Hz \n",
						proc->name,1/proc->data->dx[0]);
				}
			}
			if(frame->simData != NULL){
				printf("              Sim channels\n");
				for(sim = frame->simData; sim != NULL; sim=sim->next){
					printf(" channel %s  -> sampl. %12.9f Hz \n",
						sim->name,1/sim->data->dx[0]);
				}
			}
		}

	//	vec=FrameGetV(frame,channel); sostituita con la seguente
		vec=FrameFindVect(frame,channel);

		if(vec == NULL){ 
			printf(" no data for frame %ld\n",frame->frame);
			goto next;}
		
		ndat=vec->nData;

		if(icstart == 1){
			icstart=0;
			dt=vec->dx[0];

			fil=new_sdsfile(t0, t0gps, ULeapS, direcout, antacr, datatype, channel, sub_typ, dt);
			ic_opennew=1;

			printf(">%d data for frame,  dt = %10.7f,  ini at %16.6f \n\n",
				ndat,dt,t0);
		}

		//t0gps=frame->GTimeS+frame->GTimeN*1.e-9;
//		errt=(t0-t0A)*86400-ndat*dt;
		if(t0gpsA>0)errt=(t0gps-t0gpsA)-ndatold*dt;

		if(abs(errt)>1.e-8)printf(" *** -> time error %f s \n",errt);

		if(abs(errt)>dt){
			if(ic_opennew){
				printf(" *** File prematurely closed: frame found: %d \n",ii);
				fclose(fil);
			}
			else{
				dt=vec->dx[0];
					
				printf(">%d data for frame,  dt = %10.7f,  ini at %16.6f \n\n",
					ndat,dt,t0);
			}
			fil=new_sdsfile(t0, t0gps, ULeapS, direcout, antacr, datatype, channel, sub_typ, dt);
			ic_opennew=1;
		}

		ndatold=ndat;
		if(ndatold!=ndat)icdata=0;
		if(icdata == 0){
			data=malloc(ndat*sizeof(float));
			icdata=1;
		}

		if(vec == NULL) 
			{printf(" no data for the frame %ld\n",frame->frame);
			goto next;}
		else if(vec->type == FR_VECT_C)
			{for(i=0; i<ndat; i++) {data[i] = vec->data[i];}}
		else if(vec->type == FR_VECT_2S)
			{for(i=0; i<ndat; i++) {data[i] = vec->dataS[i];}}
		else if(vec->type == FR_VECT_4S)
			{for(i=0; i<ndat; i++) {data[i] = vec->dataI[i];}}
		else if(vec->type == FR_VECT_1U)
			{for(i=0; i<ndat; i++) {data[i] = vec->dataU[i];}}
		else if(vec->type == FR_VECT_2U)
			{for(i=0; i<ndat; i++) {data[i] = vec->dataUS[i];}}
		else if(vec->type == FR_VECT_4U)
			{for(i=0; i<ndat; i++) {data[i] = vec->dataUI[i];}}
		else if(vec->type == FR_VECT_4R)
			{for(i=0; i<ndat; i++) {data[i] = vec->dataF[i];}}
		else if(vec->type == FR_VECT_8R)
			{for(i=0; i<ndat; i++) {data[i] = (float)vec->dataD[i];}}

		//printf(" %10.4e %10.4e %10.4e %10.4e \n",data[0],data[1],data[2],data[3]);
		//for(i=0; i<ndat; i++){if(data[i] > 1)printf(" >>> %d %10.4e \n",i,data[i]);}
			//printf(" ndat = %d  \n",ndat);
		fwrite(data,4,ndat,fil);
		len+=ndat;

		t0A=t0; // corretto mar-2007
		t0gpsA=t0gps;

		next:;
		FrameFree(frame);

           /*---------------------- read the next frame --------------------*/

//		t0A=t0; corretto mar-2007
		frame = FrameRead(file); 

           /*------------------End the Main loop ----------------------------*/
	}

	FrFileIEnd(file);

	fseek(fil,24,SEEK_SET);
	fwrite(&len,4,1,fil);
	fwrite(&len1,4,1,fil);
	
	printf("File closed: frame found: %d \n",ii);
	fclose(fil);
	free(data);// sabrina 15-06-2009


	return;
}


FILE *new_sdsfile(double t0, double t0gps, long ULeapS, char *direcout, char *antacr, char *datatype, char *channel, char *sub_typ, double dt){
	FILE *fil;

	struct tm *timstr;
	char *asctim;
	char sdsfilename[40],sdsfile[100];
	char *label="#SFC#SDS";
	char *nofile="#NOFILE";
	long prot,coding,nch,pdat,len;

	timstr=GPS2UTC(t0gps,ULeapS);
	asctim=asctime(timstr);
	gw_filename(timstr, antacr, datatype, sub_typ, "sds",sdsfilename);
	strcpy(sdsfile,direcout);
	strcat(sdsfile,sdsfilename);
	printf("\n>Start at %s  on file %s \n",asctim,sdsfile);

	fil=fopen(sdsfile,"wb");
			
	fprintf(fil,"%8s",label);
	prot=1;
	fwrite(&prot,4,1,fil);
	coding=0;
	fwrite(&coding,4,1,fil);

	nch=1;
	fwrite(&nch,4,1,fil);
	pdat=944+128*nch;
	fwrite(&pdat,4,1,fil);
	len=0;
	fwrite(&len,8,1,fil);
			
	fwrite(&t0,8,1,fil);
	fwrite(&dt,8,1,fil);

	fprintf(fil,"%128s",channel);

	fprintf(fil,"%128s",sdsfile);
	fprintf(fil,"%128s",nofile);
	fprintf(fil,"%128s",nofile);
	fprintf(fil,"%128s",nofile);
	fprintf(fil,"%128s",nofile);
	fprintf(fil,"%128s",nofile);
	fprintf(fil,"%128s",channel);

	return fil;
}

