#include <stdio.h>
#include <stdlib.h>
#include<netcdf.h>
#include <dirent.h> 
#include <unistd.h>
#include <fnmatch.h>
#include<string.h>
#include<math.h>
#define PARAL
#ifdef PARAL
#include <mpi.h>
#endif
#define NBLONTS 8064
#define NBLATTS 3584
#define NBLONNDVI 43200
#define NBLATNDVI 21600
#define NBPOINTTERRE  155588764
#define NBT 5
#define PLON -180
#define DLON 179.955365724862
#define DLAT 80
#define PLAT -59.9910714285396
#define PYEAR 2003
#define DYEAR 2021
#define NBLNDVI 4481
#define NBCNDVI 5041
#define NBL 10000000
#define NBCTS 8064
#define NBLTS 3584
#define PLTS 0
#define PCTS 0
#define NBNDVI 255
#define T0 273.13
short *datelst,*datelstmax,*yearlstmax,*duree,*durmax;
float *datatmax,tluprec,*datatmaxmax,*datatmaxave;
int *datamask;
void main(int argc,char *argv[])
{
  int year,month,day,h,t,l,c,premier=1,v,nbpndvi,nblr,err,res;
  int idx=0,ncidin,ncid,nvarid,nlon,nlat,ntimep,nxi,nyi,ndim[3],nti,nvarloni,nvarlati,vartsmax,vartsmean,varyear,varday,status,nvarlst,nvartime,t2,vartsmaxave,vartamean,ncidtamean;
  int vartsmaxmax,vardaymax,varyearmax,i,varvalok,vardur,vardurmax;
  int npi,nline,ncol,nmask,nvarline,nvarcol,nvarmask;
 float datalon[NBPOINTTERRE],datalat[NBPOINTTERRE];
 int dataline[NBPOINTTERRE],datacol[NBPOINTTERRE];
 int myrank,size;
 float scal_fact,off,ndvimax;
  long ncount[4],start[4];
  short range[2],fillval=-8000;
  float fillvalf=-1e20,tlu,tsta,ta,tseuil;
  char str70[70],nfic[200],nfic2[200],fillvalb=-1;
FILE *fic;
 char *datatime,*nbvalok;
 unsigned short *datalst;
 double time[365];
  int datampi;
 #ifdef PARAL
  MPI_Status req;
  MPI_Status statuse;
#endif

#ifdef PARAL
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  printf ("process %d of %d\n",myrank,size);
#else
  size=atoi(argv[2]);
  myrank=atoi(argv[1]);
#endif
  l=myrank*NBL;
  if (l+NBL>NBPOINTTERRE) nblr=NBPOINTTERRE-l; else nblr=NBL;
  printf ("ligne debut=%d nb ligne=%d\n",l,nblr);
 duree=malloc(NBL*sizeof(short));
 durmax=malloc(NBL*sizeof(short));
 nbvalok=malloc(NBL*sizeof(char));
 datalst=malloc(73*NBL*sizeof(short));
 datelst=malloc(NBL*sizeof(short));
 datelstmax=malloc(NBL*sizeof(short));
 yearlstmax=malloc(NBL*sizeof(short));
 datatmax=malloc(NBL*sizeof(float));
 datatmaxmax=malloc(NBL*sizeof(float));
 datatmaxave=malloc(NBL*sizeof(float));
 memset(datatmaxave,0,nblr*sizeof(float));
 datamask=malloc(NBLATNDVI*NBLONNDVI*sizeof(int));
 /*   if (myrank==0)
    {
     ncid=nccreate("data/daytmaxymodis.nc",NC_64BIT_OFFSET);
     nxi = ncdimdef(ncid,"Xdim",NBLONNDVI);
     nyi = ncdimdef(ncid,"Ydim",NBLATNDVI);
     npi = ncdimdef(ncid,"nbpointterre",NBPOINTTERRE);
     nti = ncdimdef(ncid,"time",NC_UNLIMITED);
     nlon = ncvardef(ncid,"lon",NC_FLOAT,1,&npi);
     nlat = ncvardef(ncid,"lat",NC_FLOAT,1,&npi);
     nline = ncvardef(ncid,"ligne",NC_INT,1,&npi);
     ncol = ncvardef(ncid,"col",NC_INT,1,&npi);
     ndim[0]=nyi;
     ndim[1]=nxi;
     nmask = ncvardef(ncid,"mask",NC_INT,2,ndim);
     nvartime = ncvardef(ncid,"time",NC_DOUBLE,1,&nti);  
     sprintf (str70,"days since %04d-%02d-%02d %02d:%02d:%02d",PYEAR,1,1,0,0,0);
     ncattput(ncid,nvartime, "units",NC_CHAR, strlen(str70), str70);
     ncattput(ncid,nvartime, "axis",NC_CHAR, 1, "T");
     ndim[0]=nti;
     ndim[1]=npi;
     vartsmax = ncvardef(ncid,"tsmax",NC_FLOAT,2,ndim);
     ncattput(ncid, vartsmax, "_FillValue",NC_FLOAT,1,&fillvalf);
     status=ncattput (ncid,vartsmax,"units",NC_CHAR,7,"Celsius");
     status=ncattput (ncid,vartsmax,"cell_methods",NC_CHAR,20,"area:mean where land");
     varvalok = ncvardef(ncid,"nbvalok",NC_CHAR,2,ndim);
     //ncattput(ncid, vartsmax, "_FillValue",NC_CHAR,1,&fillvalb);
     status=ncattput (ncid,vartsmax,"units",NC_CHAR,1,"-");
     status=ncattput (ncid,vartsmax,"cell_methods",NC_CHAR,20,"area:mean where land");
     vartsmaxmax = ncvardef(ncid,"tsmaxmax",NC_FLOAT,1,ndim+1);
     status=ncattput (ncid,vartsmaxmax,"units",NC_CHAR,7,"Celsius");
     status=ncattput (ncid,vartsmaxmax,"cell_methods",NC_CHAR,20,"area:mean where land");
     ncattput(ncid, vartsmaxmax, "_FillValue",NC_FLOAT,1,&fillvalf);
     vartsmaxave = ncvardef(ncid,"tsmaxave",NC_FLOAT,1,ndim+1);
     status=ncattput (ncid,vartsmaxave,"units",NC_CHAR,7,"Celsius");
     status=ncattput (ncid,vartsmaxave,"cell_methods",NC_CHAR,20,"area:mean where land");
     ncattput(ncid, vartsmaxave, "_FillValue",NC_FLOAT,1,&fillvalf);
     vardaymax = ncvardef(ncid,"day_maxmax",NC_SHORT,1,ndim+1);
     varday = ncvardef(ncid,"day_max",NC_SHORT,2,ndim);
     status=ncattput (ncid,varday,"units",NC_CHAR,3,"doy");
     ncattput(ncid, varday, "_FillValue",NC_SHORT,1,&fillval);
     vardur = ncvardef(ncid,"duration",NC_SHORT,2,ndim);
     status=ncattput (ncid,vardur,"units",NC_CHAR,4,"days");
     ncattput(ncid, vardur, "_FillValue",NC_SHORT,1,&fillval);
     vardurmax = ncvardef(ncid,"durtmaxmax",NC_SHORT,1,ndim+1);
     status=ncattput (ncid,vardurmax,"units",NC_CHAR,4,"days");
     ncattput(ncid, vardur, "_FillValue",NC_SHORT,1,&fillval);
     varyearmax = ncvardef(ncid,"year_maxmax",NC_SHORT,1,ndim+1);
     status=ncattput (ncid,nlon,"units",NC_CHAR,12,"degrees_east");
     status=ncattput (ncid,nlon,"axis",NC_CHAR,1,"X");
     status=ncattput (ncid,nlat,"units",NC_CHAR,13,"degrees_north");
     status=ncattput (ncid,nlat,"axis",NC_CHAR,1,"Y");
     ncattput(ncid, vardaymax, "_FillValue",NC_SHORT,1,&fillval);
     ncattput(ncid, varyearmax, "_FillValue",NC_SHORT,1,&fillval);
     ncattput(ncid, vartsmaxmax, "_FillValue",NC_FLOAT,1,&fillvalf);
     status=ncattput (ncid,vardaymax,"units",NC_CHAR,3,"doy");
     status=ncattput (ncid,varyearmax,"units",NC_CHAR,4,"year");
     status=ncattput (ncid,varday,"cell_methods",NC_CHAR,20,"area:mean where land");
     status=ncattput (ncid,vardaymax,"cell_methods",NC_CHAR,20,"area:mean where land");
     status=ncattput (ncid,varyearmax,"cell_methods",NC_CHAR,20,"area:mean where land");
     status=ncattput (ncid,vartsmax,"units",NC_CHAR,7,"Celsius");
     status=ncattput (ncid,vartsmax,"cell_methods",NC_CHAR,20,"area:mean where land");
     ncendef(ncid);
     printf ("fin def\n");
   }
   sprintf(nfic,"data/modisglst1km_%d.nc",PYEAR);
     ncidin=ncopen(nfic,NC_NOWRITE);
 	 nc_inq_varid(ncidin,"LST_DAY",&nvarlst);
	 printf ("nvarlst=%d\n",nvarlst);
	 nc_inq_varid(ncidin,"lon",&nvarloni);
	 printf ("nvarloni=%d\n",nvarloni);
	 nc_inq_varid(ncidin,"lat",&nvarlati);
	 nc_inq_varid(ncidin,"line",&nvarline);
	 nc_inq_varid(ncidin,"col",&nvarcol);
	 nc_inq_varid(ncidin,"mask",&nvarmask);
	 start[0]=0;
	 ncount[0]=NBPOINTTERRE;
	 ncvarget(ncidin,nvarloni,start,ncount,datalon);
	 printf ("1\n");
	 err=ncvarput(ncid,nlon,start,ncount,datalon);
	 nc_strerror(err);
	 ncvarget(ncidin,nvarlati,start,ncount,datalat);
	 printf ("2\n");
	 ncvarput(ncid,nlat,start,ncount,datalat);
	 ncvarget(ncidin,nvarline,start,ncount,dataline);
	 printf ("3\n");
	 ncvarput(ncid,nline,start,ncount,dataline);
	 ncvarget(ncidin,nvarcol,start,ncount,datacol);
	 printf ("4\n"); 
	 ncvarput(ncid,ncol,start,ncount,datacol);
	 start[0]=start[1]=0;
	 ncount[0]=NBLATNDVI;
	 ncount[1]=NBLONNDVI;
	 ncvarget(ncidin,nvarmask,start,ncount,datamask);
	 printf ("4\n");
	 ncvarput(ncid,nmask,start,ncount,datamask);	     
	 printf ("5\n");
	 ncclose(ncid);
	 ncclose(ncidin);
	 exit(0);*/
	 /*  else
   {
     printf ("sleep\n");
     sleep(30);
     printf ("attente message");
#ifdef PARAL
     MPI_Recv(&datampi,1,MPI_INT,myrank-1,0,MPI_COMM_WORLD,&req);
     if (datampi!=myrank-1) {printf ("erreur message %f\n",datampi);exit(1);} else printf ("message recu\n");
     #endif
 */
   ncid=ncopen("data/daytmaxymodis.nc",NC_64BIT_OFFSET|NC_WRITE);
      printf ("ouverture fichier\n");
     nc_inq_varid(ncid,"tsmax",&vartsmax);
     nc_inq_varid(ncid,"nbvalok",&varvalok);
     nc_inq_varid(ncid,"tsmaxmax",&vartsmaxmax);
     nc_inq_varid(ncid,"tsmaxave",&vartsmaxave);
     nc_inq_varid(ncid,"day_maxmax",&vardaymax);     
     nc_inq_varid(ncid,"day_max",&varday);
     nc_inq_varid(ncid,"year_maxmax",&varyearmax);
     nc_inq_varid(ncid,"duration",&vardur);
     nc_inq_varid(ncid,"durtmaxmax",&vardurmax);
    /*}
     #ifdef PARAL
    datampi=myrank;
    printf ("envoi message proc suivant\n");
   if (myrank<size-1)  MPI_Ssend(&datampi,1,MPI_INT,myrank+1,0,MPI_COMM_WORLD);
    #endif
    */
      /*     sprintf(nfic,"data/modisglst1km_%d.nc",PYEAR);
     ncidin=ncopen(nfic,NC_NOWRITE);
 	 nc_inq_varid(ncidin,"LST_DAY",&nvarlst);
	 printf ("nvarlst=%d\n",nvarlst);
	 nc_inq_varid(ncidin,"lon",&nvarloni);
	 printf ("nvarloni=%d\n",nvarloni);
	 nc_inq_varid(ncidin,"lat",&nvarlati);
	 nc_inq_varid(ncidin,"line",&nvarline);
	 nc_inq_varid(ncidin,"col",&nvarcol);
	 nc_inq_varid(ncidin,"mask",&nvarmask);
	 start[0]=0;
	 ncount[0]=NBPOINTTERRE;
	 ncvarget(ncidin,nvarloni,start,ncount,datalon);
	 printf ("1\n");
	 err=ncvarput(ncid,nlon,start,ncount,datalon);
	 nc_strerror(err);
	 ncvarget(ncidin,nvarlati,start,ncount,datalat);
	 printf ("2\n");
	 ncvarput(ncid,nlat,start,ncount,datalat);
	 ncvarget(ncidin,nvarline,start,ncount,dataline);
	 printf ("3\n");
	 ncvarput(ncid,nline,start,ncount,dataline);
	 ncvarget(ncidin,nvarcol,start,ncount,datacol);
	 printf ("4\n"); 
	 ncvarput(ncid,ncol,start,ncount,datacol);
	 start[0]=start[1]=0;
	 ncount[0]=NBLATNDVI;
	 ncount[1]=NBLONNDVI;
	 ncvarget(ncidin,nvarmask,start,ncount,datamask);
	 printf ("4\n");
	 ncvarput(ncid,nmask,start,ncount,datamask);	     
	 printf ("5\n");
	 ncclose(ncid);
	 ncclose(ncidin);
	 exit(0);*/
   
for (t=0;t<DYEAR-PYEAR+1;t++) time[t]=t*365+182;
  for (year=PYEAR;year<=DYEAR;year++)
   {
      printf ("year %d\n",year);
     sprintf(nfic,"data/modisglst1km_%d.nc",year);
     ncidin=ncopen(nfic,NC_NOWRITE);
	 nc_inq_varid(ncidin,"LST_DAY",&nvarlst);
	 printf ("nvarlst=%d\n",nvarlst);
     /*   if ((year==PYEAR))
       {
	 nc_inq_varid(ncidin,"LST_DAY",&nvarlst);
	 printf ("nvarlst=%d\n",nvarlst);
	 nc_inq_varid(ncidin,"lon",&nvarloni);
	 printf ("nvarloni=%d\n",nvarloni);
	 nc_inq_varid(ncidin,"lat",&nvarlati);
	 nc_inq_varid(ncidin,"line",&nvarline);
	 nc_inq_varid(ncidin,"col",&nvarcol);
	 nc_inq_varid(ncidin,"mask",&nvarmask);
	 start[0]=0;
	 ncount[0]=NBPOINTTERRE;
	 ncvarget(ncidin,nvarloni,start,ncount,datalon);
	 printf ("1\n");
	 err=ncvarput(ncid,nlon,start,ncount,datalon);
	 nc_strerror(err);
	 ncvarget(ncidin,nvarlati,start,ncount,datalat);
	 printf ("2\n");
	 ncvarput(ncid,nlat,start,ncount,datalat);
	 ncvarget(ncidin,nvarline,start,ncount,dataline);
	 printf ("3\n");
	 ncvarput(ncid,nline,start,ncount,dataline);
	 ncvarget(ncidin,nvarcol,start,ncount,datacol);
	 printf ("4\n"); 
	 ncvarput(ncid,ncol,start,ncount,datacol);
	 start[0]=start[1]=0;
	 ncount[0]=NBLATNDVI;
	 ncount[1]=NBLONNDVI;
	 ncvarget(ncidin,nvarmask,start,ncount,datamask);
	 printf ("4\n");
	 ncvarput(ncid,nmask,start,ncount,datamask);	     
	 printf ("5\n");
	 }*/
	 //     for (l=0;l<NBPOINTTERRE;l+=NBL)
       {
	 if (l+NBL>NBPOINTTERRE) nblr=NBPOINTTERRE-l; else nblr=NBL;
	 printf ("nblr=%d\n",nblr);
 	 start[0]=0;
	 start[1]=l;
	 ncount[0]=73;
	 ncount[1]=nblr;
	 ncvarget(ncidin,nvarlst,start,ncount,datalst);
	 memset(datatmax,0,nblr*sizeof(float));
	 memset(duree,0,nblr*sizeof(short));
	 for (i=0;i<nblr;i++)
	   {
	     for (t=nbvalok[i]=0;t<73;t++)
	       {
		 if (datalst[t*nblr+i]>0) nbvalok[i]++;
		 tlu=(float)datalst[t*nblr+i]*0.02-273.13;
		 if ((tlu>datatmax[i])&&(tlu<80)&&((t>0)&&(fabs(tlu-tluprec)<20)))
		   {		 		
		     datatmax[i]=tlu;
		     datelst[i]=t*5;
		   }
		 tluprec=tlu;
	       }
	     if (datatmax[i]>=43)
	       {
		 for (t=datelst[i]/5;(t>=0)&&((datalst[t*nblr+i]==0)||((float)(datalst[t*nblr+i]*0.02-273.13>=43)));t--) duree[i]++;
		 for (t=datelst[i]/5;(t<73)&&((datalst[t*nblr+i]==0)||((float)(datalst[t*nblr+i]*0.02-273.13>=43)));t++) duree[i]++;
	       }
	     else duree[i]=-1;
	    
	   }
	 	 for (i=0;i<nblr;i++)
	   {
	     datatmaxave[i]+=datatmax[i]/(DYEAR-PYEAR+1); 
	     if (datatmax[i]>datatmaxmax[i])
	       {
		 datatmaxmax[i]=datatmax[i];
		 datelstmax[i]=datelst[i];
		 durmax[i]=duree[i];
		 yearlstmax[i]=year;
	       }
	       }
#ifdef PARAL
	 if (myrank>0)
	   {
	     MPI_Recv(&datampi,1,MPI_INT,myrank-1,0,MPI_COMM_WORLD,&req);
	     if (datampi!=myrank-1) {printf ("erreur message %f\n",datampi);exit(1);} else printf ("message recu\n");
	   }
#endif
	 start[0]=year-PYEAR;
	 start[1]=l;
	 ncount[0]=1;
	 ncount[1]=nblr;
	 printf ("1b\n");
	 res=ncvarput(ncid,vardur,start,ncount,duree);
	 ncvarput(ncid,varday,start,ncount,datelst);
	 printf ("2b\n");
	 res=ncvarput(ncid,vartsmax,start,ncount,datatmax);
	 printf ("3b\n");
	 ncvarput(ncid,varvalok,start,ncount,nbvalok);
	 printf ("4b\n");
	 //if (res!=0) nc_strerror(res);
	 printf ("5b\n");
#ifdef PARAL
    datampi=myrank;
    printf ("envoi message proc suivant\n");
   if (myrank<size-1)  MPI_Ssend(&datampi,1,MPI_INT,myrank+1,0,MPI_COMM_WORLD);
#endif
       }
       ncclose(ncidin);
   }
  start[0]=l;
  ncount[0]=nblr;
      ncvarput(ncid,vardaymax,start,ncount,datelstmax);
  ncvarput(ncid,vartsmaxmax,start,ncount,datatmaxmax);
  ncvarput(ncid,vartsmaxave,start,ncount,datatmaxave);
  ncvarput(ncid,varyearmax,start,ncount,yearlstmax); 	 
  ncvarput(ncid,vardurmax,start,ncount,durmax); 	 
  ncclose(ncid);
#ifdef PARAL 
 MPI_Barrier (MPI_COMM_WORLD); 
 MPI_Finalize();
#endif

}





