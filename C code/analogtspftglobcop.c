#include <stdio.h>
#include <stdlib.h>
#include<netcdf.h>
#include <dirent.h> 
#include <unistd.h>
#include <fnmatch.h>
#include<string.h>
#include<math.h>
#include <mpi.h>
#define NBLONTS 8064
#define NBLATTS 3584
#define NBLONNDVI 40320
#define NBLATNDVI 15680
#define PLON -180
#define DLON 179.955365724862
#define DLAT 80
#define PLAT -59.9910714285396
#define PYEAR 1999
#define DYEAR 2019
#define NBLNDVI 4481
#define NBCNDVI 5041
#define NBL 512
#define NBCTS 8064
#define NBLTS 3584
#define PLNDVI 561
#define PCNDVI 18480
#define PLTS 0
#define PCTS 0
#define NBNDVI 255
#define NVEG 7
#define ERRHAND(e) if ((e) !=NC_NOERR) printf ("%s\n",nc_strerror((e)))


double time[365],datalonts[NBLONTS],datalatts[NBLATTS];
//char datavi[NBLATNDVI][NBCNDVI];
float  *datavi;
float  *datavim;
float  *dataveg;
/*short datats[NBL][NBCTS];
short datatsta[NBLTS][NBCTS];
float datatstam[NBLTS][NBCTS];
float datatsm[NBLTS][NBCTS];*/
float datapc[NBL][NBCTS];
float datapcm[NBL][NBCTS];
short datatstad[NBL][NBCTS];
short datatstadm[NBL][NBCTS];
float datafitx[4000], datafity[4000],coefa[NBL][NBCTS],coefb[NBL][NBCTS],abdev[NBL][NBCTS],diffano[DYEAR-PYEAR+1][NBL][NBCTS];
int nbvalprec[NBL][NBCTS],nbvv[NVEG][NBL][NBCTS],nbvalprecavant[NBL][NBCTS],nbvalprecapres[NBL][NBCTS],nbvalts[NBL][NBCTS];
int nbvalndvi[NVEG][NBL][NBCTS],nbvalndviavant[NVEG][NBL][NBCTS],nbvalndviapres[NVEG][NBL][NBCTS],l2,c2,idx,nbvmin[NBL][NBCTS],nbvmints[NBL][NBCTS];
float *anondvi,*anondviavant,*anondviapres,*ndvim,*ndvimin,*ndvimax;
//  float predano[DYEAR-PYEAR+1][NBL][NBCTS];
float *anoprec,*anoprecavant,*anoprecapres,*anots,*anotsta;
float *anondvimin,*anondvivar,*anondviminavant,*anondviminapres,minanoprec[NVEG][NBL][NBCTS],minanoprecr[NBL][NBCTS],*yearmin,*anondvimean,anoprecmean[NBL][NBCTS],mints[NBL][NBCTS],mintsta[NBL][NBCTS],*ndviap;
void main(int argc,char *argv[])
{
  long ncount[4],start[4],varh,varder,var2Dder,var2Dderm,var2Dsec,var2Dsecm;
  int nvartime,nvarts,nvarvi,vartest,nvard2m,nvarpc,nvarpcm,nv,nvarveg,vveg;
  int ncidndvi,ncidndvi2,ncidts,ncidpc,ncidtsm,nvarlonvi,nvarlatvi,nvarlonts,nvarlatts,ncid,nxi,nyi,nxxi,nyyi,ndim[4],status,vartestvi,ncidmean,ncidtsmean,nvartsm,nvarvim,varanoprec,varanondvi,varpredanondvi,varanots,varanota,varts,varta,varanandvi,varanandvivar,varminp,varyearmin,varmeanp,varanandvim;
  int varanandviav,varanandviap,varanatstaav,varanatstaap,varanats,varanatsta,varanotsta,varanoprecav,varanoprecap,varanotstad,varanotstadap,varanotstadav,varantstaav,varantstaap,varndvim,varndvimin,varndvimax,varndviap;
  int ncidparam,vardifano,ncidpcmean,nvartstam,nvartsta,nvardmax,nvartsmax,nvartstad,nvartsd,ncidtsd;
  int i,j,pcndvi,ncndvi,plndvi,nlndvi,pcts,ncts,plts,nlts,year,l,c,satus,t,nxip,nyip,ntip,nlonp,nlatp,ntimep;
  int cur,curm1,curp1,jj,ti,yearref,v,y;
  char nfic[100],nficts[100],str70[70];
  float *datatsmax,*datatsta;
  short *datadmax,dmax;
  double val,val2,val3,val4,val5,val6,dts;
  float fillval=-1e34;
  MPI_Status req;
  MPI_Status statuse;
  int datampi;
  //float anota[DYEAR-PYEAR+1][NBLTS][NBCTS],anots[DYEAR-PYEAR+1][NBLTS][NBCTS], ta[DYEAR-PYEAR+1][NBLTS][NBCTS],ts[DYEAR-PYEAR+1][NBLTS][NBCTS];
 int blk,nbl;
 int myrank=0,size=1,plproc,atraiter,nbproc;
    MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  nbproc=NBLTS/NBL;
 yearref=2011+myrank/nbproc;
 blk=myrank%nbproc;
  anoprec=calloc((DYEAR-PYEAR+1)*NBL*NBCTS,sizeof(float));
 anoprecavant=calloc((DYEAR-PYEAR+1)*NBL*NBCTS,sizeof(float));
 anoprecapres=calloc((DYEAR-PYEAR+1)*NBL*NBCTS,sizeof(float));
 anots=calloc((DYEAR-PYEAR+1)*NBL*NBCTS,sizeof(float));
 anotsta=calloc((DYEAR-PYEAR+1)*NBL*NBCTS,sizeof(float));
 anondvi=calloc((DYEAR-PYEAR+1)*NBL*NBCTS*NVEG,sizeof(float));
 yearmin=calloc(NBL*NBCTS*NVEG,sizeof(float));

 anondviavant=calloc((DYEAR-PYEAR+1)*NBL*NBCTS*NVEG,sizeof(float));
 anondviapres=calloc((DYEAR-PYEAR+1)*NBL*NBCTS*NVEG,sizeof(float));
 anondvimean=calloc(NBL*NBCTS*NVEG,sizeof(float));
anondvimin=calloc (NBL*NBCTS*NVEG,sizeof(float));
anondvivar=calloc (NBL*NBCTS*NVEG,sizeof(float));
anondviminavant=calloc (NBL*NBCTS*NVEG,sizeof(float));
anondviminapres=calloc (NBL*NBCTS*NVEG,sizeof(float));
ndvim=calloc (NBL*NBCTS*NVEG,sizeof(float));
ndvimin=calloc (NBL*NBCTS*NVEG,sizeof(float));
 for (l=0;l<NBL;l++) for (c=0;c<NBCTS;c++) for (v=0;v<NVEG;v++) ndvimin[v*NBL*NBCTS+l*NBCTS+c]=1000;
ndvimax=calloc (NBL*NBCTS*NVEG,sizeof(float));

ndviap=calloc (NBL*NBCTS*NVEG,sizeof(float));
 dataveg=malloc (NBL*NBCTS*NVEG*sizeof(float));
 datavi=malloc (NBL*NBCTS*NVEG*sizeof(float));

 datavim=malloc (NBL*NBCTS*NVEG*sizeof(float));
 datadmax=calloc(9*NBL*NBCTS,sizeof(short));
 datatsmax=calloc(9*NBL*NBCTS,sizeof(float));
 datatsta=calloc(9*NBL*NBCTS,sizeof(float));
    sprintf (nfic,"data/ndvipfttscopfglob_%d.nc",yearref);
ncidmean=ncopen(nfic,NC_NOWRITE);
 sprintf (nfic,"era5/rain_%d.nc",yearref);
 ncidpcmean=ncopen(nfic,NC_NOWRITE);
 status=nc_inq_varid(ncidmean,"NDVI",&nvarvim);
 printf ("ncidmean=%d nvar=%d\n",ncidmean,nvarvim);
 status=nc_inq_varid(ncidpcmean,"rain",&nvarpcm);
 if (blk==0)
   {
     printf("creat\n");
    sprintf (nfic,"data/analogtpftglob_%d.nc",yearref);
     ncidparam=nccreate(nfic,NC_64BIT_OFFSET);
     ntip = ncdimdef(ncidparam,"time",NC_UNLIMITED);
     nxip = ncdimdef(ncidparam,"lon",NBCTS);
     nyip = ncdimdef(ncidparam,"lat",NBLTS);
     nv  = ncdimdef(ncidparam,"veg",NVEG);
     nlonp = ncvardef(ncidparam,"lon",NC_DOUBLE,1,&nxip);
     nlatp = ncvardef(ncidparam,"lat",NC_DOUBLE,1,&nyip);
     ntimep = ncvardef(ncidparam,"time",NC_DOUBLE,1,&ntip);
     ndim[0]=ntip;
     ndim[1]=nv;
     ndim[2]=nyip;
     ndim[3]=nxip;
     varanondvi = ncvardef(ncidparam,"anondvi",NC_FLOAT,4,ndim);
     varanandvi = ncvardef(ncidparam,"anandvi",NC_FLOAT,3,ndim+1);
     varanandvivar = ncvardef(ncidparam,"anandvivar",NC_FLOAT,3,ndim+1);
     varanandviav = ncvardef(ncidparam,"anandviav",NC_FLOAT,3,ndim+1);
     varanandviap = ncvardef(ncidparam,"anandviap",NC_FLOAT,3,ndim+1);
    varanandvim = ncvardef(ncidparam,"anondvim",NC_FLOAT,3,ndim+1);
     varndvim=ncvardef(ncidparam,"ndvim",NC_FLOAT,3,ndim+1);
     varndvimin=ncvardef(ncidparam,"ndvimin",NC_FLOAT,3,ndim+1);
     varndvimax=ncvardef(ncidparam,"ndvimax",NC_FLOAT,3,ndim+1);
   varndviap=ncvardef(ncidparam,"ndviap",NC_FLOAT,3,ndim+1);
     vveg = ncvardef(ncidparam,"veg",NC_FLOAT,3,ndim+1);
   varyearmin= ncvardef(ncidparam,"yearmin",NC_FLOAT,3,ndim+1);
      ndim[0]=ntip;
     ndim[1]=nyip;
     ndim[2]=nxip;
     varanoprec = ncvardef(ncidparam,"anoprec",NC_FLOAT,3,ndim);
     varanots = ncvardef(ncidparam,"anots",NC_FLOAT,3,ndim);
     varanotsta = ncvardef(ncidparam,"anotsta",NC_FLOAT,3,ndim);
     varanats = ncvardef(ncidparam,"anats",NC_FLOAT,2,ndim+1);
      varanatsta = ncvardef(ncidparam,"anatsta",NC_FLOAT,2,ndim+1);
 
     varminp= ncvardef(ncidparam,"mindifprec",NC_FLOAT,2,ndim+1);
     varanoprecav= ncvardef(ncidparam,"anoprecav",NC_FLOAT,2,ndim+1);
     varanoprecap= ncvardef(ncidparam,"anoprecap",NC_FLOAT,2,ndim+1);
     varmeanp= ncvardef(ncidparam,"meandifprec",NC_FLOAT,2,ndim+1);
     sprintf (str70,"days since %04d-%02d-%02d %02d:%02d:%02d",1999,1,1,0,0,0);
     ncattput(ncidparam,ntimep, "units",NC_CHAR, strlen(str70), str70);
     ncattput(ncidparam,ntimep, "axis",NC_CHAR, 1, "T");
     status=ncattput (ncidparam,nlonp,"units",NC_CHAR,12,"degrees_east");
     status=ncattput (ncidparam,nlonp,"axis",NC_CHAR,1,"X");
     status=ncattput (ncidparam,nlatp,"units",NC_CHAR,13,"degrees_north");
     status=ncattput (ncidparam,nlatp,"axis",NC_CHAR,1,"Y");
     ncattput(ncidparam, varanoprec, "_FillValue",NC_FLOAT,1,&fillval);
     ncattput(ncidparam, varanondvi, "_FillValue",NC_FLOAT,1,&fillval);
      ncattput(ncidparam, varanots, "_FillValue",NC_FLOAT,1,&fillval);
      ncattput(ncidparam, varanotsta, "_FillValue",NC_FLOAT,1,&fillval);
     ncattput(ncidparam, varanats, "_FillValue",NC_FLOAT,1,&fillval);
     printf ("ncat 1\n");
     ncattput(ncidparam, varanatsta, "_FillValue",NC_FLOAT,1,&fillval);
    ncattput(ncidparam, varanandvi, "_FillValue",NC_FLOAT,1,&fillval);
   ncattput(ncidparam, vveg, "_FillValue",NC_FLOAT,1,&fillval);
   ncattput(ncidparam, varanandvivar, "_FillValue",NC_FLOAT,1,&fillval);
      ncattput(ncidparam, varanandviav, "_FillValue",NC_FLOAT,1,&fillval);
      ncattput(ncidparam, varanandviap, "_FillValue",NC_FLOAT,1,&fillval);
     ncattput(ncidparam, varanandvim, "_FillValue",NC_FLOAT,1,&fillval);
     printf ("ncat 2\n");
     ncattput(ncidparam, varminp, "_FillValue",NC_FLOAT,1,&fillval);
     ncattput(ncidparam, varmeanp, "_FillValue",NC_FLOAT,1,&fillval);
     ncattput(ncidparam, varyearmin, "_FillValue",NC_FLOAT,1,&fillval);
    ncattput(ncidparam, varndvim, "_FillValue",NC_FLOAT,1,&fillval);
    ncattput(ncidparam, varndvimin, "_FillValue",NC_FLOAT,1,&fillval);
    ncattput(ncidparam, varndvimax, "_FillValue",NC_FLOAT,1,&fillval);
    ncattput(ncidparam, varndviap, "_FillValue",NC_FLOAT,1,&fillval);
     ncattput(ncidparam, varndviap, "_FillValue",NC_FLOAT,1,&fillval);
   printf ("ncat 3\n");
     printf ("ncat 4\n");
     ncendef(ncidparam);
     ncclose(ncidparam);
   }
 printf ("block %d year %d\n",blk,yearref);
 if (blk*NBL+NBL>NBLTS) nbl=NBLTS-blk*NBL; else nbl=NBL;
 printf ("nbl=%d\n",nbl);
 for (v=0;v<NVEG;v++)
   for (l=0;l<nbl;l++)
     for (c=0;c<NBCTS;c++)
       {
	 minanoprec[v][l][c]=1e10;
       }
 ncidts=ncopen("data/daytmaxyglob.nc",NC_NOWRITE);
 nc_inq_varid(ncidts,"day_max",&nvardmax);
 nc_inq_varid(ncidts,"tsmax",&nvartsmax);
 nc_inq_varid(ncidts,"tsta_max",&nvartsta);
 start[0]=0;
 ncount[0]=9;
 start[1]=(blk*NBL);
 ncount[1]=nbl;
 start[2]=0;
 ncount[2]=NBCTS;
 // printf ("get 1 %d %d\n",ncidts,nvardmax);
 status=ncvarget(ncidts,nvardmax,start,ncount,datadmax);
   ERRHAND(status);
   //printf ("get 2\n");
 status=ncvarget(ncidts,nvartsmax,start,ncount,datatsmax);
  ERRHAND(status);
  // printf ("get 3\n");
status=ncvarget(ncidts,nvartsta,start,ncount,datatsta);
   ERRHAND(status);
 for (year=PYEAR;year<=DYEAR;year++)
   {
     printf (" **************** %d **********************\n",year);
     sprintf (nfic,"data/ndvipfttscopfglob_%d.nc",year);
     ncidndvi=ncopen(nfic,NC_NOWRITE);
     sprintf (nfic,"era5/rain_%d.nc",year);
     ncidpc=ncopen(nfic,NC_NOWRITE);
     if (year<2011)  sprintf (nfic,"lst/tsdfb_2011.nc");
     else sprintf (nfic,"lst/tsdfb_%d.nc",year);
     ncidtsd=ncopen(nfic,NC_NOWRITE);
     nc_inq_varid(ncidndvi,"NDVI",&nvarvi);
    nc_inq_varid(ncidndvi,"VEG",&nvarveg);
      nc_inq_varid(ncidpc,"rain",&nvarpc);
     nc_inq_varid(ncidts,"LST_MAX",&nvarts);
     nc_inq_varid(ncidts,"TSTA",&nvartsta);
     memset(nbvalndvi,0,NVEG*nbl*NBCTS*sizeof(int));
     memset(nbvalndviavant,0,NVEG*nbl*NBCTS*sizeof(int));
     memset(nbvalndviapres,0,NVEG*nbl*NBCTS*sizeof(int));
     memset(nbvalprec,0,nbl*NBCTS*sizeof(int));
     memset(nbvalts,0,nbl*NBCTS*sizeof(int));
     memset(&anoprec[(year-PYEAR)*nbl*NBCTS],0,nbl*NBCTS*sizeof(float));
     start[0]=start[2]=0;
     start[1]=blk*NBL;
     ncount[0]=NVEG;
     ncount[1]=nbl;
     ncount[2]=NBCTS;
    ncvarget(ncidndvi,nvarveg,start,ncount,dataveg);
      for (t=60;t<244;t++)
       {
	 start[0]=t;
	 ncount[0]=1;
	 start[1]=PLTS+blk*NBL;
	 ncount[1]=nbl;
	 start[2]=PCTS;
	 ncount[2]=NBCTS;
	     status=ncvarget(ncidpc,nvarpc,start,ncount,datapc);
 	  status=ncvarget(ncidpcmean,nvarpcm,start,ncount,datapcm);
  	  //printf ("get4\n");
	 ERRHAND(status);
	 if (t%10==0)
	   {	   
	     printf ("%d\n",t);
	     start[0]=(t/10+1)%36;
	     ncount[0]=1;
	     start[1]=0;
	     ncount[1]=NVEG;
	     start[2]=(blk*NBL);
	     ncount[2]=nbl;
	     start[3]=0;
	     ncount[3]=NBCTS;
	     printf ("get5 %d %d %d %d %d %d %d %d\n",start[0],ncount[0],start[1],ncount[1],start[2],ncount[2],start[3],ncount[3]); 
	     ncvarget(ncidndvi,nvarvi,start,ncount,datavi);
	     printf ("get6 %d %d %d %d %d %d %d %d\n",start[0],ncount[0],start[1],ncount[1],start[2],ncount[2],start[3],ncount[3]); 
	     ncvarget(ncidmean,nvarvim,start,ncount,datavim);
	      //printf ("get6\n");
	   }
	 for (l=0;l<nbl;l++)
	   for (c=0;c<NBCTS;c++)
	     if (datadmax[(yearref-2011)*nbl*NBCTS+(l)*NBCTS+c]>=0)
	     {
	       dmax=datadmax[(yearref-2011)*nbl*NBCTS+(l)*NBCTS+c];
	       ti=t/10;
	       if (ti>35) ti=35;
	       for (v=0;v<NVEG;v++)
		 if ((t%10==0)&&(datavi[v*nbl*NBCTS+l*NBCTS+c]>0)/*&&(dmax-30>=60)&&(dmax+30<244)*/)
		 {
		   anondvi[(year-PYEAR)*NVEG*nbl*NBCTS+v*nbl*NBCTS+l*NBCTS+c]+=(double)datavim[v*nbl*NBCTS+l*NBCTS+c]-(double)datavi[v*nbl*NBCTS+l*NBCTS+c];
		   if (year==yearref)
		     {
		       ndvim[v*nbl*NBCTS+l*NBCTS+c]+=datavim[v*nbl*NBCTS+l*NBCTS+c];
		       ndvimin[v*nbl*NBCTS+l*NBCTS+c]=(datavim[v*nbl*NBCTS+l*NBCTS+c]<ndvimin[v*nbl*NBCTS+l*NBCTS+c])?datavim[v*nbl*NBCTS+l*NBCTS+c]:ndvimin[v*nbl*NBCTS+l*NBCTS+c];
		       ndvimax[v*nbl*NBCTS+l*NBCTS+c]=(datavim[v*nbl*NBCTS+l*NBCTS+c]>ndvimax[v*nbl*NBCTS+l*NBCTS+c])?datavim[v*nbl*NBCTS+l*NBCTS+c]:ndvimax[v*nbl*NBCTS+l*NBCTS+c];
		     }
		   nbvalndvi[v][l][c]++;
		   if ((t<dmax)&&(t>=dmax-30))
		     {
		       anondviavant[(year-PYEAR)*NVEG*nbl*NBCTS+v*nbl*NBCTS+l*NBCTS+c]=(double)datavim[v*nbl*NBCTS+l*NBCTS+c]-(double)datavi[v*nbl*NBCTS+l*NBCTS+c];
		       nbvalndviavant[v][l][c]++;
		     }
		   else
		     if ((t>dmax)&&(t<=dmax+30))
		     {
		       anondviapres[(year-PYEAR)*NVEG*nbl*NBCTS+v*nbl*NBCTS+l*NBCTS+c]+=(double)datavim[v*nbl*NBCTS+l*NBCTS+c]-(double)datavi[v*nbl*NBCTS+l*NBCTS+c];
		       nbvalndviapres[v][l][c]++;
		       if (year==yearref) ndviap[v*nbl*NBCTS+l*NBCTS+c]+=datavim[v*nbl*NBCTS+l*NBCTS+c];
		     }															 
		 }
	     }
	 for (l=0;l<nbl;l++)
	   for (c=0;c<NBCTS;c++)
	     {
	       if ((datapc[l][c]>=0)&&(datapcm[l][c]>=0))
		 {
		   anoprec[(year-PYEAR)*nbl*NBCTS+l*NBCTS+c]+=(datapcm[l][c]-datapc[l][c]);
		   nbvalprec[l][c]++;
		 }
	       if (t<datadmax[(yearref-2011)*nbl*NBCTS+(l)*NBCTS+c])
		 {
		   anoprecavant[(year-PYEAR)*nbl*NBCTS+l*NBCTS+c]+=(datapcm[l][c]-datapc[l][c]);
		   nbvalprecavant[l][c]++;
		 }
	       else
		 {
		   anoprecapres[(year-PYEAR)*nbl*NBCTS+l*NBCTS+c]+=(datapcm[l][c]-datapc[l][c]);
		   nbvalprecapres[l][c]++;
		 }	       
	     }
       }
     
     for (l=0;l<nbl;l++)
       for (c=0;c<NBCTS;c++)
	 for (v=0;v<NVEG;v++)
	     {
	       if (nbvalndvi[v][l][c]>0) 
		 {
		   anondvi[(year-PYEAR)*NVEG*nbl*NBCTS+v*nbl*NBCTS+l*NBCTS+c]/=nbvalndvi[v][l][c]; 
		    if (year==yearref) ndvim[v*nbl*NBCTS+l*NBCTS+c]/=nbvalndvi[v][l][c];
		   if ( anondvi[(year-PYEAR)*NVEG*nbl*NBCTS+v*nbl*NBCTS+l*NBCTS+c]<-1) anondvi[(year-PYEAR)*NVEG*nbl*NBCTS+v*nbl*NBCTS+l*NBCTS+c]=fillval;
		   if ( anondvi[(year-PYEAR)*NVEG*nbl*NBCTS+v*nbl*NBCTS+l*NBCTS+c]>1) anondvi[(year-PYEAR)*NVEG*nbl*NBCTS+v*nbl*NBCTS+l*NBCTS+c]=fillval;
		 }
	       else 
		 {
		   anondvi[(year-PYEAR)*NVEG*nbl*NBCTS+v*nbl*NBCTS+l*NBCTS+c]=fillval;
		      if (year==yearref) ndvim[v*nbl*NBCTS+l*NBCTS+c]= ndvimin[v*nbl*NBCTS+l*NBCTS+c]=ndvimax[v*nbl*NBCTS+l*NBCTS+c]=fillval;
		 }
	       if (nbvalndviavant[v][l][c]>0) 
		 {
		   anondviavant[(year-PYEAR)*NVEG*nbl*NBCTS+v*nbl*NBCTS+l*NBCTS+c]/=nbvalndviavant[v][l][c]; 
		   if ( anondviavant[(year-PYEAR)*NVEG*nbl*NBCTS+v*nbl*NBCTS+l*NBCTS+c]<-1) anondviavant[(year-PYEAR)*NVEG*nbl*NBCTS+v*nbl*NBCTS+l*NBCTS+c]=fillval;
		   if ( anondviavant[(year-PYEAR)*NVEG*nbl*NBCTS+v*nbl*NBCTS+l*NBCTS+c]>1) anondviavant[(year-PYEAR)*NVEG*nbl*NBCTS+v*nbl*NBCTS+l*NBCTS+c]=fillval;
		 }
	       else  anondviavant[(year-PYEAR)*NVEG*nbl*NBCTS+v*nbl*NBCTS+l*NBCTS+c]=fillval;
	       if (nbvalndviapres[v][l][c]>0) 
		 {
		   anondviapres[(year-PYEAR)*NVEG*nbl*NBCTS+v*nbl*NBCTS+l*NBCTS+c]/=nbvalndviapres[v][l][c]; 
		   if (year==yearref) ndviap[v*nbl*NBCTS+l*NBCTS+c]/=nbvalndviapres[v][l][c];
		   if ( anondviapres[(year-PYEAR)*NVEG*nbl*NBCTS+v*nbl*NBCTS+l*NBCTS+c]<-1) anondviapres[(year-PYEAR)*NVEG*nbl*NBCTS+v*nbl*NBCTS+l*NBCTS+c]=fillval;
		   if ( anondviapres[(year-PYEAR)*NVEG*nbl*NBCTS+v*nbl*NBCTS+l*NBCTS+c]>1) anondviapres[(year-PYEAR)*NVEG*nbl*NBCTS+v*nbl*NBCTS+l*NBCTS+c]=fillval;
		 }
	       else 
		 {
		   anondviapres[(year-PYEAR)*NVEG*nbl*NBCTS+v*nbl*NBCTS+l*NBCTS+c]=fillval;
		    if (year==yearref) ndviap[v*nbl*NBCTS+l*NBCTS+c]=fillval;
		 }
	     }
     for (l=0;l<nbl;l++)
       for (c=0;c<NBCTS;c++)
	 {
	   if (nbvalprec[l][c]>0) anoprec[(year-PYEAR)*nbl*NBCTS+l*NBCTS+c]/=nbvalprec[l][c];
	   else anoprec[(year-PYEAR)*nbl*NBCTS+l*NBCTS+c]=fillval;
	   if (nbvalprecavant[l][c]>0) anoprecavant[(year-PYEAR)*nbl*NBCTS+l*NBCTS+c]/=nbvalprecavant[l][c];
	   else anoprecavant[(year-PYEAR)*nbl*NBCTS+l*NBCTS+c]=fillval;
	   if (nbvalprecapres[l][c]>0) anoprecapres[(year-PYEAR)*nbl*NBCTS+l*NBCTS+c]/=nbvalprecapres[l][c];
	   else anoprecapres[(year-PYEAR)*nbl*NBCTS+l*NBCTS+c]=fillval;
	 }
   }

 for (l=0;l<nbl;l++)
   for (c=0;c<NBCTS;c++) 
	 {
	   anoprecmean[l][c]=nbvmin[l][c]=nbvmints[l][c]=0;
	   mints[l][c]= mintsta[l][c]= minanoprecr[l][c]=0;
	   for (v=0;v<NVEG;v++)	 anondviminavant[v*nbl*NBCTS+l*NBCTS+c]=anondviminapres[v*nbl*NBCTS+l*NBCTS+c]=anondvimin[v*nbl*NBCTS+l*NBCTS+c]=fillval;
	   for (year=PYEAR;year<=DYEAR;year++)
	     if (year!=yearref)
	       {
		 if ((anoprec[(year-PYEAR)*nbl*NBCTS+l*NBCTS+c]>fillval))
		   {
		     nbvmin[l][c]=1;
		     if (datadmax[(yearref-2011)*nbl*NBCTS+l*NBCTS+c]>=0)
		       {
			 if ((year>=2011)&&(yearref>=2011)) anots[(year-PYEAR)*nbl*NBCTS+l*NBCTS+c]=datatsmax[(yearref-2011)*nbl*NBCTS+l*NBCTS+c]-datatsmax[(year-2011)*nbl*NBCTS+l*NBCTS+c];
			 if ((year>=2011)&&(yearref>=2011)) anotsta[(year-PYEAR)*nbl*NBCTS+l*NBCTS+c]=datatsta[(yearref-2011)*nbl*NBCTS+l*NBCTS+c]-datatsta[(year-2011)*nbl*NBCTS+l*NBCTS+c];
			 if ((year<2011)||((yearref<2011)||anots[(year-PYEAR)*nbl*NBCTS+l*NBCTS+c]>100)) anotsta[(year-PYEAR)*nbl*NBCTS+l*NBCTS+c]=anots[(year-PYEAR)*nbl*NBCTS+l*NBCTS+c]=fillval;
		       }
		     else  anotsta[(year-PYEAR)*nbl*NBCTS+l*NBCTS+c]=anots[(year-PYEAR)*nbl*NBCTS+l*NBCTS+c]=fillval;
		     //		     if ((fabs(anoprec[(year-PYEAR)*nbl*NBCTS+l*NBCTS+c])<minanoprec[l][c])&&(fabs(anoprec[(year-PYEAR)*nbl*NBCTS+l*NBCTS+c]*1000<0.1)))
		     for (v=0;v<NVEG;v++)
		       {
		
			 if (((fabs(anondvi[(year-PYEAR)*NVEG*nbl*NBCTS+v*nbl*NBCTS+l*NBCTS+c])<minanoprec[v][l][c]))&&(anondvi[(year-PYEAR)*NVEG*nbl*NBCTS+v*nbl*NBCTS+l*NBCTS+c]>fillval))
			   {
			     minanoprec[v][l][c]=fabs(anondvi[(year-PYEAR)*NVEG*nbl*NBCTS+v*nbl*NBCTS+l*NBCTS+c]);
			     minanoprecr[l][c]=anoprec[(year-PYEAR)*nbl*NBCTS+l*NBCTS+c];
			     if (year>=2011)
			       {
				 mints[l][c]=anots[(year-PYEAR)*nbl*NBCTS+l*NBCTS+c];
				 mintsta[l][c]=anotsta[(year-PYEAR)*nbl*NBCTS+l*NBCTS+c]; 
				 //   nbvmints[l][c]++;
				 nbvmints[l][c]=1;
			       }
			     yearmin[v*nbl*NBCTS+l*NBCTS+c]=year;
			     anondvimin[v*nbl*NBCTS+l*NBCTS+c]=anondvi[(year-PYEAR)*NVEG*nbl*NBCTS+v*nbl*NBCTS+l*NBCTS+c]; 
			     anondviminavant[v*nbl*NBCTS+l*NBCTS+c]=anondviavant[(year-PYEAR)*NVEG*nbl*NBCTS+v*nbl*NBCTS+l*NBCTS+c]; 
			     anondviminapres[v*nbl*NBCTS+l*NBCTS+c]=anondviapres[(year-PYEAR)*NVEG*nbl*NBCTS+v*nbl*NBCTS+l*NBCTS+c]; 
			   }
		       }
			 anoprecmean[l][c]+=anoprec[(year-PYEAR)*nbl*NBCTS+l*NBCTS+c]/(float)(DYEAR-PYEAR);
		   }	
	 //else anoprecmean[l][c]=minanoprecr[l][c]=yearmin[l][c]=anondvimin[l][c]= anondviminavant[l][c]=anondviminapres[l][c]= mints[l][c]=mintsta[l][c]=fillval;
		 v=0;
		 for (nbvv[v][l][c]=0;v<NVEG;v++) 
		   if (anondvi[(year-PYEAR)*NVEG*nbl*NBCTS+v*nbl*NBCTS+l*NBCTS+c]>fillval) 
		     {
		         anondvimean[v*nbl*NBCTS+l*NBCTS+c]+=anondvi[(year-PYEAR)*NVEG*nbl*NBCTS+v*nbl*NBCTS+l*NBCTS+c];
		        nbvv[v][l][c]++;
		       }

	       }
	 }
 for (l=0;l<nbl;l++)
   for (c=0;c<NBCTS;c++) 
     {
       if (nbvmints[l][c]>0)
	 {
	   mints[l][c]/=nbvmints[l][c];
	   mintsta[l][c]/=nbvmints[l][c];
	   if (( mints[l][c]<-20)||(mints[l][c]>100)) mints[l][c]=mintsta[l][c]=fillval;
	 }
       else
	 mints[l][c]=mintsta[l][c]=fillval;
       for (v=0;v<NVEG;v++) if (nbvv[v][l][c]>0) {anondvimean[v*nbl*NBCTS+l*NBCTS+c]/=nbvv[v][l][c];} else anondvimean[v*nbl*NBCTS+l*NBCTS+c]=fillval;
     }
 for (l=0;l<nbl;l++)
   for (c=0;c<NBCTS;c++) 
     {
       if (dataveg[l*NBCTS+c]<0)
	 {
	   for (y=0;y<(DYEAR-PYEAR+1);y++) anoprec[y*nbl*NBCTS+l*NBCTS+c]=anots[y*nbl*NBCTS+l*NBCTS+c]=anoprecavant[y*nbl*NBCTS+l*NBCTS+c]=anoprecapres[y*nbl*NBCTS+l*NBCTS+c]=anotsta[y*nbl*NBCTS+l*NBCTS+c]=fillval;
	 mints[l][c]=mintsta[l][c]=minanoprecr[l][c]=anoprecmean[l][c]=fillval;
	 }
       for (v=0;v<NVEG;v++)
	 if (dataveg[v*nbl*NBCTS+l*NBCTS+c]<0.1)
	   {
	     anondvimean[v*nbl*NBCTS+l*NBCTS+c]=anondvimin[v*nbl*NBCTS+l*NBCTS+c]=anondviminavant[v*nbl*NBCTS+l*NBCTS+c]=anondviminapres[v*nbl*NBCTS+l*NBCTS+c]=yearmin[v*nbl*NBCTS+l*NBCTS+c]=fillval;
	     ndvim[v*nbl*NBCTS+l*NBCTS+c]=ndviap[v*nbl*NBCTS+l*NBCTS+c]=fillval;
	     for (y=0;y<(DYEAR-PYEAR+1);y++) anondvi[(year-PYEAR)*NVEG*nbl*NBCTS+v*nbl*NBCTS+l*NBCTS+c]=fillval;
	   }
     }
 if (blk==0) for (t=0;t<DYEAR-PYEAR+1;t++) time[t]=t*365+182;
 if (myrank>0) 
   {
     MPI_Recv(&datampi,1,MPI_INT,myrank-1,0,MPI_COMM_WORLD,&req);
     if (datampi!=myrank-1) {printf ("erreur message %f\n",datampi);exit(1);}
   }
 sprintf (nfic,"data/analogtpftglob_%d.nc",yearref);
 ncidparam=ncopen(nfic,NC_WRITE);
 printf ("inq1\n");
 status=nc_inq_varid(ncidparam,"anoprec",&varanoprec);
 ERRHAND(status);
 printf ("inq2\n");
 status=nc_inq_varid(ncidparam,"anondvi",&varanondvi);
 ERRHAND(status);
 printf ("inq3\n");
 status=nc_inq_varid(ncidparam,"anots",&varanots);
 ERRHAND(status);
 printf ("inq4\n");
 status=nc_inq_varid(ncidparam,"anotsta",&varanotsta);
 ERRHAND(status);
 printf ("inq5\n");
 status=nc_inq_varid(ncidparam,"anandvi",&varanandvi);
 ERRHAND(status);
 status=nc_inq_varid(ncidparam,"anandvivar",&varanandvivar);
 ERRHAND(status);
     printf ("inq6\n");
    status=nc_inq_varid(ncidparam,"anandviav",&varanandviav);
   ERRHAND(status);
     printf ("inq7\n");
   status=nc_inq_varid(ncidparam,"anandviap",&varanandviap);
   ERRHAND(status);
     printf ("inq8\n");
     status=nc_inq_varid(ncidparam,"anondvim",&varanandvim);
   ERRHAND(status);
     printf ("inq9\n");
     status=nc_inq_varid(ncidparam,"anondvi",&varanondvi);
   ERRHAND(status);
     printf ("inq10\n");
   status= nc_inq_varid(ncidparam,"mindifprec",&varminp);
   ERRHAND(status);
  status= nc_inq_varid(ncidparam,"anoprecav",&varanoprecav);
   ERRHAND(status);
  status= nc_inq_varid(ncidparam,"anoprecap",&varanoprecap);
   ERRHAND(status);
     printf ("inq11\n");
   status= nc_inq_varid(ncidparam,"meandifprec",&varmeanp);
   ERRHAND(status);
     printf ("inq12\n");
   status=nc_inq_varid(ncidparam,"yearmin",&varyearmin);
   ERRHAND(status);
     printf ("inq13\n");
   status=nc_inq_varid(ncidparam,"anats",&varanats);
   ERRHAND(status);
     printf ("inq14\n");
   status=nc_inq_varid(ncidparam,"anatsta",&varanatsta);
	       ERRHAND(status);
     printf ("inq15\n");
   status=nc_inq_varid(ncidparam,"ndvim",&varndvim);
   status=nc_inq_varid(ncidparam,"ndvimin",&varndvimin);
 status=nc_inq_varid(ncidparam,"ndvimax",&varndvimax);
   status=nc_inq_varid(ncidparam,"ndviap",&varndviap);
   ERRHAND(status);
  status=nc_inq_varid(ncidparam,"veg",&vveg);
   ERRHAND(status);
     printf ("inq15\n");

  if (blk==0)
     {
       sprintf (nfic,"lst/tsdfb_%d.nc",yearref);
       ncidts=ncopen(nfic,NC_NOWRITE);
       nc_inq_varid(ncidmean,"lon",&nvarlonvi);
       nc_inq_varid(ncidmean,"lat",&nvarlatvi);
       nc_inq_varid(ncidts,"lon",&nvarlonts);
       nc_inq_varid(ncidts,"lat",&nvarlatts);
       start[0]=0;
       ncount[0]=NBLONTS;
       ncvarget(ncidts,nvarlonts,start,ncount,datalonts);
       ncount[0]=NBLATTS;
       ncvarget(ncidts,nvarlatts,start,ncount,datalatts);
       for (i=0;(i<NBLONTS)&&(datalonts[i]<PLON);i++);
       pcts=i;      
       for (;(i<NBLONTS)&&(datalonts[i]<=DLON);i++);
       ncts=i-pcts+1;
       for (j=0;(j<NBLATTS)&&(datalatts[j]>DLAT);j++);
       plts=j;
       for (;(j<NBLATTS)&&(datalatts[j]>PLAT);j++);
       nlts=j-plts+1;
       printf ("plts=%d nblts=%d pcts=%d nbcts=%d\n",plts,nlts,pcts,ncts);
       start[0]=0;
       ncount[0]=NBLTS;
       printf ("write1\n");
       ncvarput(ncidparam,nlatp,start,ncount,&datalatts[plts]);	  
       start[0]=0;
       ncount[0]=NBCTS;
       ncvarput(ncidparam,nlonp,start,ncount,&datalonts[pcts]);
       printf("write2\n");
     }
   start[0]=0;
   start[2]=0;
   start[1]=blk*NBL;
   ncount[0]=DYEAR-PYEAR+1;
   ncount[1]=nbl;
   ncount[2]=NBCTS; 
   printf ("%d %d %d %d %d %d\n",start[0],start[1],start[2],ncount[0],ncount[1],ncount[2]);
   ncvarput(ncidparam,varanoprec,start,ncount,anoprec);
   printf ("put 1\n");
   printf ("put 2n");
   ncvarput(ncidparam,varanots,start,ncount,anots);
   printf ("put 3\n");
   ncvarput(ncidparam,varanotsta,start,ncount,anotsta);
   if (blk==0) ncvarput(ncidparam,ntimep,start,ncount,time);
   ncount[0]=nbl;
   ncount[1]=NBCTS;
   start[0]=blk*NBL;
   start[1]=0;
  printf ("put 4\n");
   printf ("put 3\n");
   status=ncvarput(ncidparam,varanats,start,ncount,mints);
   printf ("put 4\n");
   status=ncvarput(ncidparam,varanatsta,start,ncount,mintsta);
  printf ("put 5\n");
    ERRHAND(status);
  printf ("put 6\n");
    status=ncvarput(ncidparam,varminp,start,ncount,minanoprecr);
   ERRHAND(status);
    status=ncvarput(ncidparam,varanoprecav,start,ncount,anoprecavant);
   ERRHAND(status);
    status=ncvarput(ncidparam,varanoprecap,start,ncount,anoprecapres);
   ERRHAND(status);
  printf ("put 7\n");
    status=ncvarput(ncidparam,varmeanp,start,ncount,anoprecmean);
   ERRHAND(status);
  printf ("put 8\n");
   ncount[0]=NVEG;
   ncount[1]=nbl;
   ncount[2]=NBCTS;
   start[1]=blk*NBL;
   start[2]=start[0]=0;
 printf ("put 10\n");
    status=ncvarput(ncidparam,varanandvi,start,ncount,anondvimin);
  printf ("put 11\n");
    status=ncvarput(ncidparam,varanandviav,start,ncount,anondviminavant);
 printf ("put 12\n");
  status=ncvarput(ncidparam,varanandviap,start,ncount,anondviminapres);
 printf ("put 13\n");
   status= ncvarput(ncidparam,varanandvim,start,ncount,anondvimean);
   printf ("put 9 %d \n",vveg);
   status= ncvarput(ncidparam,vveg,start,ncount,dataveg);
   status= ncvarput(ncidparam,varndvim,start,ncount,ndvim);
   status= ncvarput(ncidparam,varndvimin,start,ncount,ndvimin);
   status= ncvarput(ncidparam,varndvimax,start,ncount,ndvimax);
    ERRHAND(status);
   status=ncvarput(ncidparam,varndviap,start,ncount,ndviap); 
   ERRHAND(status);
    status=ncvarput(ncidparam,varyearmin,start,ncount,yearmin); 
   ERRHAND(status);
  printf ("put 9\n");
  ncount[0]=DYEAR-PYEAR+1;
   ncount[1]=NVEG;
   ncount[2]=nbl;
   ncount[3]=NBCTS;
   start[2]=blk*NBL;
   start[3]=start[1]=start[0]=0;
   ncvarput(ncidparam,varanondvi,start,ncount,anondvi);
 printf ("put 9\n");
   status=ncclose(ncidparam);
   printf ("fin ecriture, on envoie au suivant\n");
   datampi=myrank;
   if (myrank<size-1)  MPI_Ssend(&datampi,1,MPI_INT,myrank+1,0,MPI_COMM_WORLD);
   printf ("message reçu on a fini\n");
  MPI_Barrier (MPI_COMM_WORLD);
   MPI_Finalize();  
}
     

      
