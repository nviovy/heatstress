#include <stdio.h>
#include <stdlib.h>
#include<netcdf.h>
#include <dirent.h> 
#include <unistd.h>
#include <fnmatch.h>
#include<string.h>
#include<math.h>
#include <stdbool.h>
#define NBLONTS 8064
#define NBLATTS 3584
#define NBT 5
#define PLON -180
#define DLON 179.955365724862
#define DLAT 80
#define PLAT -59.9910714285396
#define PYEAR 2011
#define DYEAR 2019
#define NBLNDVI 4481
#define NBCNDVI 5041
#define NBL 100
#define NBCTS 8064
#define NBLTS 3584
#define PLTS 0
#define PCTS 0
#define NBNDVI 255
#define T0 273.13
short datelst[NBLTS][NBCTS],datelstmax[NBLTS][NBCTS],yearlstmax[NBLTS][NBCTS],datanbj[NBLTS][NBCTS];
float datatmax[NBLTS][NBCTS],datatstam[NBLTS][NBCTS],tluprec[NBLTS][NBCTS],datatmaxmax[NBLTS][NBCTS],datatmaxave[NBLTS][NBCTS],datatamean[NBLTS][NBCTS],datavi[7][NBLTS][NBCTS];
short duree[NBLTS][NBCTS],dureem[NBLTS][NBCTS];
bool stop[NBLTS][NBCTS];
void main(int argc,char *argv[])
{
  int year,month,day,h,t,l,c,premier=1,v,nbpndvi;
  int idx=0,ncidin,ncid,nvarid,nlon,nlat,ntimep,nxi,nyi,ndim[3],nti,nvarloni,nvarlati,vartsmax,vartsmean,varyear,varday,vartsta,status,nvarlst,nvartime,t2,nvartsta,vartsmaxave,varnbj,vartamean,ncidtamean;
  int vartsmaxmax,vardaymax,varyearmax,nvarvi,ncidvi,vardur;
 double datalon[NBCTS],datalat[NBLTS];
 float scal_fact,off,ndvimax;
  long ncount[4],start[4];
  short range[2],fillval=-8000;
  float fillvalf=-1e20,tlu,tsta,ta,tseuil;
  char str70[70],nfic[200],nfic2[200];
FILE *fic;
 char *datatime;
 short *datalst,*datatsta;
 double time[365];
 datalst=malloc(NBT*NBLTS*NBCTS*sizeof(short));
 datatsta=malloc(NBT*NBLTS*NBCTS*sizeof(short));
 ncid=nccreate("data/daytmaxyglob2.nc",NC_CLOBBER);
 nxi = ncdimdef(ncid,"lon",NBCTS);
  nyi = ncdimdef(ncid,"lat",NBLTS);
     nti = ncdimdef(ncid,"time",NC_UNLIMITED);
 nlon = ncvardef(ncid,"lon",NC_DOUBLE,1,&nxi);
  nlat = ncvardef(ncid,"lat",NC_DOUBLE,1,&nyi);
    ntimep = ncvardef(ncid,"time",NC_DOUBLE,1,&nti);
    sprintf (str70,"days since %04d-%02d-%02d %02d:%02d:%02d",PYEAR,1,1,0,0,0);
     ncattput(ncid,ntimep, "units",NC_CHAR, strlen(str70), str70);
     ncattput(ncid,ntimep, "axis",NC_CHAR, 1, "T");
   ndim[0]=nti;
 ndim[1]=nyi;
  ndim[2]=nxi;
 vartsmax = ncvardef(ncid,"tsmax",NC_FLOAT,3,ndim);
 ncattput(ncid, vartsmax, "_FillValue",NC_FLOAT,1,&fillvalf);
 status=ncattput (ncid,vartsmax,"units",NC_CHAR,7,"Celsius");
 status=ncattput (ncid,vartsmax,"cell_methods",NC_CHAR,20,"area:mean where land");
 varnbj = ncvardef(ncid,"nbjcrit",NC_SHORT,3,ndim);
 ncattput(ncid, varnbj, "_FillValue",NC_SHORT,1,&fillval);
 status=ncattput (ncid,varnbj,"units",NC_CHAR,4,"Days");
 status=ncattput (ncid,varnbj,"cell_methods",NC_CHAR,20,"area:mean where land");
 vardur = ncvardef(ncid,"duration",NC_SHORT,3,ndim);
 ncattput(ncid, vardur, "_FillValue",NC_SHORT,1,&fillval);
 status=ncattput (ncid,vardur,"units",NC_CHAR,4,"Days");
 status=ncattput (ncid,vardur,"cell_methods",NC_CHAR,20,"area:mean where land");
 vartsmaxmax = ncvardef(ncid,"tsmaxmax",NC_FLOAT,2,ndim+1);
 status=ncattput (ncid,vartsmaxmax,"units",NC_CHAR,7,"Celsius");
 status=ncattput (ncid,vartsmaxmax,"cell_methods",NC_CHAR,20,"area:mean where land");
 ncattput(ncid, vartsmaxmax, "_FillValue",NC_FLOAT,1,&fillvalf);
 vartsmaxave = ncvardef(ncid,"tsmaxave",NC_FLOAT,2,ndim+1);
 status=ncattput (ncid,vartsmaxave,"units",NC_CHAR,7,"Celsius");
 status=ncattput (ncid,vartsmaxave,"cell_methods",NC_CHAR,20,"area:mean where land");
 ncattput(ncid, vartsmaxave, "_FillValue",NC_FLOAT,1,&fillvalf);
 vartsta = ncvardef(ncid,"tsta_max",NC_FLOAT,3,ndim);
 vardaymax = ncvardef(ncid,"day_maxmax",NC_SHORT,2,ndim+1);
 varday = ncvardef(ncid,"day_max",NC_SHORT,3,ndim);
status=ncattput (ncid,varday,"units",NC_CHAR,3,"doy");
 ncattput(ncid, varday, "_FillValue",NC_SHORT,1,&fillval);
 varyearmax = ncvardef(ncid,"year_maxmax",NC_SHORT,2,ndim+1);
 
 status=ncattput (ncid,nlon,"units",NC_CHAR,12,"degrees_east");
 status=ncattput (ncid,nlon,"axis",NC_CHAR,1,"X");
 status=ncattput (ncid,nlat,"units",NC_CHAR,13,"degrees_north");
 status=ncattput (ncid,nlat,"axis",NC_CHAR,1,"Y");
 ncattput(ncid, vardaymax, "_FillValue",NC_SHORT,1,&fillval);
 ncattput(ncid, varyearmax, "_FillValue",NC_SHORT,1,&fillval);
 ncattput(ncid, vartsta, "_FillValue",NC_FLOAT,1,&fillvalf);
 ncattput(ncid, vartsmaxmax, "_FillValue",NC_FLOAT,1,&fillvalf);
status=ncattput (ncid,vardaymax,"units",NC_CHAR,3,"doy");
status=ncattput (ncid,varyearmax,"units",NC_CHAR,4,"year");
 status=ncattput (ncid,varday,"cell_methods",NC_CHAR,20,"area:mean where land");
 status=ncattput (ncid,vardaymax,"cell_methods",NC_CHAR,20,"area:mean where land");
 status=ncattput (ncid,varyearmax,"cell_methods",NC_CHAR,20,"area:mean where land");
 status=ncattput (ncid,vartsmax,"units",NC_CHAR,7,"Celsius");
 status=ncattput (ncid,vartsmax,"cell_methods",NC_CHAR,20,"area:mean where land");
 ncendef(ncid);
  ncidtamean=ncopen("data/meantemp.nc",NC_NOWRITE);
  status=nc_inq_varid(ncidtamean,"tmean",&vartamean);
  ncount[0]=NBLTS;
  ncount[1]=NBCTS;
  start[0]=0;
  start[1]=0;
  status=ncvarget(ncidtamean,vartamean,start,ncount,datatamean);      
for (t=0;t<DYEAR-PYEAR+1;t++) time[t]=t*365+182;
    for (l=0;l<NBLTS;l++)
      for (c=0;c<NBCTS;c++)   {datatmaxmax[l][c]=fillvalf;datatmaxave[l][c]=0;datelstmax[l][c]=yearlstmax[l][c]=fillval;}
  for (year=PYEAR;year<=DYEAR;year++)
   {
     for (l=0;l<NBLTS;l++)
       for (c=0;c<NBCTS;c++)   { datatmax[l][c]= datatstam[l][c]=fillvalf;datelst[l][c]=fillval;datanbj[l][c]=duree[l][c]=dureem[l][c]=0;}
     
     printf ("year %d\n",year);
     sprintf(nfic,"lst/tsdfb_%d.nc",year);
     ncidin=ncopen(nfic,NC_NOWRITE);
     nc_inq_varid(ncidin,"LST_MAX",&nvarlst);
     nc_inq_varid(ncidin,"TSTA",&nvartsta);
     /*   sprintf(nfic,"data/ndvipfttscopfglob_%d.nc",year);
     ncidvi=ncopen(nfic,NC_NOWRITE);
     nc_inq_varid(ncidvi,"NDVI",&nvarvi);*/
     for (t=0;t<73;t++)
       {
	 printf ("%d.",t);fflush(stdout);
	 start[0]=t*NBT;
	 start[1]=PLTS;
	 start[2]=PCTS;
	 ncount[0]=NBT;
	 ncount[1]=NBLTS;
	 ncount[2]=NBCTS;
         ncvarget(ncidin,nvarlst,start,ncount,datalst);
	 ncvarget(ncidin,nvartsta,start,ncount,datatsta);
	 /*	 if (t%2==0)
	   {
	     if (t>=70)  start[0]=35; else start[0]=t/2;
	     start[1]=0;
	     start[2]=PLTS;
	     start[3]=PCTS;
	     ncount[0]=1;
	     ncount[1]=7;
	     ncount[2]=NBLTS;
	     ncount[3]=NBCTS;
	     ncvarget(ncidvi,nvarvi,start,ncount,datavi);
	     }*/
	 if (premier)
	   {
	     nc_inq_varid(ncidin,"lon",&nvarloni);
	     nc_inq_varid(ncidin,"lat",&nvarlati);
	     start[0]=PCTS;
	     ncount[0]=NBCTS;
	     ncvarget(ncidin,nvarloni,start,ncount,datalon);
	     start[0]=0;
	     ncvarput(ncid,nlon,start,ncount,datalon);
	     start[0]=PLTS;
	     ncount[0]=NBLTS;
	     ncvarget(ncidin,nvarlati,start,ncount,datalat);
	     start[0]=0;
	     ncvarput(ncid,nlat,start,ncount,datalat);
	     premier=0;
	   }
 	 for (t2=0;t2<NBT;t2++)    
	   for (l=0;l<NBLTS;l++)
	     for (c=0;c<NBCTS;c++)
	       	       if ((datalst[t2*NBLTS*NBCTS+l*NBCTS+c]>=-7000)&&((year!=2013)||(c<2640)||(c>2750)||(l<1237)||(l>2000)||(t*NBT+t2<40)||(t*NBT+t2>190)))
	       if ((datalst[t2*NBLTS*NBCTS+l*NBCTS+c]>=-7000))
		 {
		   tlu=(float)datalst[t2*NBLTS*NBCTS+l*NBCTS+c]*0.01;
		   tsta=(float)datatsta[t2*NBLTS*NBCTS+l*NBCTS+c]*0.01;
		   if ((year==2011)&&(t==0)&&(t2==0)) tluprec[l][c]=tlu;
		   //		   if ((tlu>datatmax[l][c])&&(tlu<80)&&(fabs(tlu-tluprec[l][c])<20))
		   if ((tlu>datatmax[l][c]))
		     {
		       datatmax[l][c]=tlu;
		        datatstam[l][c]=tsta;
		       datelst[l][c]=t*NBT+t2+1;
		     }
		   ta=datatamean[l][c];
		   tseuil=(0.6*(ta)+38.1);
		   if (tseuil<38) tseuil=38;
		   if (tlu>tseuil)
		     {
		       tluprec[l][c]=tlu;
		     }
		 }
       }
     for (l=0;l<NBLTS;l++)
       for (c=0;c<NBCTS;c++)
	 {
	   duree[l][c]=0;
	   stop[l][c]=false;
	 }
     
     /*     for (t=0;t<73;t++)
       {
	 printf ("%d.",t);fflush(stdout);
	 start[0]=t*NBT;
	 start[1]=PLTS;
	 start[2]=PCTS;
	 ncount[0]=NBT;
	 ncount[1]=NBLTS;
	 ncount[2]=NBCTS;
         ncvarget(ncidin,nvarlst,start,ncount,datalst);
 	 for (t2=0;t2<NBT;t2++)    
	   for (l=0;l<NBLTS;l++)
	     for (c=0;c<NBCTS;c++)
	       {
		 if (datatmax[l][c]>=42.5)
		   {
		     tlu=(float)datalst[t2*NBLTS*NBCTS+l*NBCTS+c]*0.01;
		     if (datatmax[l][c]>=t*NBT+t2)
		       {
			 if (tlu>=42.5) duree[l][c]++; else  duree[l][c]=0;
		       }
		     else
		       {
			 if  ((tlu >=42.5)&&(!stop[l][c])) duree[l][c]++;else stop[l][c]=true;	 
		       }
		   }
	       }
	       }*/
     ncclose(ncidin);
     //    ncclose(ncidvi);
     for (l=0;l<NBLTS;l++)
       for (c=0;c<NBCTS;c++)
	 {
	   if (datatmax[l][c]>fillvalf) datatmaxave[l][c]+=datatmax[l][c]/(DYEAR-PYEAR+1); else {datatmaxave[l][c]=fillvalf;datanbj[l][c]=dureem[l][c]=fillval;}
	   if (datatmax[l][c]>datatmaxmax[l][c])
	   {
	     datatmaxmax[l][c]=datatmax[l][c];
	     datelstmax[l][c]=datelst[l][c];
	     yearlstmax[l][c]=year;
	   }
	 }
     start[0]=year-PYEAR;
     start[1]=start[2]=0;
     ncount[0]=1;
     ncount[1]=NBLTS;
     ncount[2]=NBCTS;
     ncvarput(ncid,varday,start,ncount,datelst);
     ncvarput(ncid,vartsmax,start,ncount,datatmax);
     ncvarput(ncid,vartsta,start,ncount,datatstam);
     ncvarput(ncid,varnbj,start,ncount,datanbj);
     ncvarput(ncid,vardur,start,ncount,duree);
     ncvarput(ncid,ntimep,start,ncount,time);
   }
  start[0]=start[1]=0;
  ncount[0]=NBLTS;
  ncount[1]=NBCTS;
  ncvarput(ncid,vardaymax,start,ncount,datelstmax);
  ncvarput(ncid,vartsmaxmax,start,ncount,datatmaxmax);
  ncvarput(ncid,vartsmaxave,start,ncount,datatmaxave);
  ncvarput(ncid,varyearmax,start,ncount,yearlstmax);
 	 
ncclose(ncid);
}





