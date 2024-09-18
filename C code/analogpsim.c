#include <stdio.h>
#include <stdlib.h>
#include<netcdf.h>
#include <dirent.h> 
#include <unistd.h>
#include <fnmatch.h>
#include<string.h>
#include<math.h>
//#include <mpi.h>
#define NBLTS 360
#define NBCTS 720
#define PYEAR 2001
#define DYEAR 2022
#define TP0 273.16
#define NBPOINTTERRE 67420
#define NVEG 13
#define ERRHAND(e) if ((e) !=NC_NOERR) printf ("%s\n",nc_strerror((e)))


float time[365],datalonvi[NBCTS],datalatvi[NBLTS],datalonts[NBCTS],datalatts[NBLTS];
float *datandvi;
float *datavim; 
float  *dataveg;
float *datapcm;
float *datapc;
float *datatsm;
float *datats;
//float datatam[NBPOINTTERRE];
int nbvalprec[NBPOINTTERRE],nbvalprecavant[NBPOINTTERRE],nbvalprecapres[NBPOINTTERRE],l2,c2,idx,nbvmin[NVEG][NBPOINTTERRE],nbvmints[NBPOINTTERRE];
int nbvalndvi[NVEG][NBPOINTTERRE],nbvalndviavant[NVEG][NBPOINTTERRE],nbvalndviapres[NVEG][NBPOINTTERRE],l2,c2,idx,nbvv[NVEG][NBPOINTTERRE];
float *anondvi,*anondviavant,*anondviapres,*ndvim,*ndvimin,*ndvimax;
float *anondvi,*anoprec,*anoprecavant,*anoprecapres,*anondviavant,*anondviapres,*anots,*anotsta;
float *anondvimin,*anondviminavant,*anondviminapres,minanoprec[NVEG][NBPOINTTERRE],minanoprecr[NBPOINTTERRE],*yearmin,*anondvimean,anoprecmean[NBPOINTTERRE],mints[NBPOINTTERRE],mintsta[NBPOINTTERRE];
float *ndvim,*ndviap,vtmax[NBPOINTTERRE];
int dataidligne[NBPOINTTERRE],dataidcol[NBPOINTTERRE],datamask[NBLTS][NBCTS];
float fpar(float lai)
{
  // return 1.-exp(-0.5*lai);
  return lai;
}
void main(int argc,char *argv[])
{
  long ncount[3],start[3],varh,varder,var2Dder,var2Dderm,var2Dsec,var2Dsecm;
  int nvartime,nvarts,nvarvi,vartest,nvard2m,nvarpc,nvarpcm;
  int ncidndvi,ncidndvi2,ncidts,ncidpc,ncidtsm,nvarlonvi,nvarlatvi,nvarlonts,nvarlatts,ncid,nxi,nyi,nxxi,nyyi,ndim[3],status,vartestvi,ncidmean,ncidtsmean,nvartsm,nvarvim,varanoprec,varanondvi,varpredanondvi,varanots,varanota,varts,varta,varanandvi,varminp,varyearmin,varmeanp,varanandvim,vartsmax;
  int varanandviav,varanandviap,varanats,varanatsta,varanotsta,varanoprecav,varanoprecap,varndvim,varndvimin,varndvimax,nvarmask,nvaridligne,nvaridcol,nmaskp,nvarveg,nv,varndviap,vveg;
  int ncidparam,vardifano,ncidpcmean,nvartstam,nvartsta,nvardmax,nvartsmax,nidlignep,nidcolp;
  int i,j,pcndvi,ncndvi,plndvi,nlndvi,pcts,ncts,plts,nlts,year,l,c,satus,t,nxip,nyip,ntip,npip,nlonp,nlatp,ntimep,v,vardmax;
  int cur,curm1,curp1,jj,ti,yearref,t1,t2,t3;
  char nfic[100],nficts[100],str70[70];
  float *datatsmax,*datatsta;
  float *datadmax,lat;
  double val,val2,val3,val4,val5,val6,dts,y;
  float fillval=-1e34;
  int datampi;
  //float anota[DYEAR-PYEAR+1][NBLTS][NBCTS],anots[DYEAR-PYEAR+1][NBLTS][NBCTS], ta[DYEAR-PYEAR+1][NBLTS][NBCTS],ts[DYEAR-PYEAR+1][NBLTS][NBCTS];
 int blk=0;
 int myrank=0,size=1,plproc,atraiter;
 /*    MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Barrier (MPI_COMM_WORLD);
  yearref=atoi(argv[1])+myrank;*/
 yearref=atoi(argv[1]);
 anoprec=calloc((DYEAR-PYEAR+1)*NBPOINTTERRE,sizeof(float));
 anoprecavant=calloc((DYEAR-PYEAR+1)*NBPOINTTERRE,sizeof(float));
 anoprecapres=calloc((DYEAR-PYEAR+1)*NBPOINTTERRE,sizeof(float));
 anots=calloc((DYEAR-PYEAR+1)*NBPOINTTERRE,sizeof(float));
 anotsta=calloc((DYEAR-PYEAR+1)*NBPOINTTERRE,sizeof(float));
 anondvi=calloc((DYEAR-PYEAR+1)*NBPOINTTERRE*NVEG,sizeof(float));
 yearmin=calloc(NBPOINTTERRE*NVEG,sizeof(float));

 anondviavant=calloc((DYEAR-PYEAR+1)*NBPOINTTERRE*NVEG,sizeof(float));
 anondviapres=calloc((DYEAR-PYEAR+1)*NBPOINTTERRE*NVEG,sizeof(float));
 anondvimean=calloc(NBPOINTTERRE*NVEG,sizeof(float));
anondvimin=calloc (NBPOINTTERRE*NVEG,sizeof(float));
anondviminavant=calloc (NBPOINTTERRE*NVEG,sizeof(float));
anondviminapres=calloc (NBPOINTTERRE*NVEG,sizeof(float));
ndvim=calloc (NBPOINTTERRE*NVEG,sizeof(float));
ndvimin=calloc (NBPOINTTERRE*NVEG,sizeof(float));
 for (i=0;i<NBPOINTTERRE;i++) for (v=0;v<NVEG;v++) ndvimin[v*NBPOINTTERRE+i]=1000;
ndvimax=calloc (NBPOINTTERRE*NVEG,sizeof(float));

ndviap=calloc (NBPOINTTERRE*NVEG,sizeof(float));
 dataveg=malloc (NBPOINTTERRE*NVEG*sizeof(float));
 datandvi=malloc (365*NBPOINTTERRE*NVEG*sizeof(float));
 datavim=malloc (365*NBPOINTTERRE*NVEG*sizeof(float));
 datapcm=malloc(365*NBPOINTTERRE*sizeof(float));
 if (datapcm==NULL) {printf ("error alloc datapcm\n");exit(1);}
 datapc=malloc(365*NBPOINTTERRE*sizeof(float));
 if (datapc==NULL) {printf ("error alloc datapc\n");exit(1);}
 datatsm=malloc(365*NBPOINTTERRE*sizeof(float));
 if (datatsm==NULL) {printf ("error alloc datatsm\n");exit(1);}
datats=malloc(365*NBPOINTTERRE*sizeof(float));
 datadmax=calloc(NBPOINTTERRE,sizeof(float));
 datatsmax=calloc((DYEAR-PYEAR+1)*NBPOINTTERRE,sizeof(float));

   //  if (myrank==0)
   {
     sprintf (nfic,"output/analogpftorchths_%d.nc",yearref);
     ncidparam=nccreate(nfic,NC_CLOBBER);
     printf ("creation \n");
     nxip = ncdimdef(ncidparam,"longitude",NBCTS);
     nyip = ncdimdef(ncidparam,"latitude",NBLTS);
     ntip = ncdimdef(ncidparam,"time_counter",NC_UNLIMITED);
     npip = ncdimdef(ncidparam,"points_terre",NBPOINTTERRE);
     nv = ncdimdef(ncidparam,"veget",NVEG);
     nlonp = ncvardef(ncidparam,"longitude",NC_FLOAT,1,&nxip);
     nlatp = ncvardef(ncidparam,"latitude",NC_FLOAT,1,&nyip);
     ntimep = ncvardef(ncidparam,"time",NC_FLOAT,1,&ntip);
     nidlignep=ncvardef(ncidparam,"indice_ligne",NC_INT,1,&npip);
     nidcolp=ncvardef(ncidparam,"indice_col",NC_INT,1,&npip);
     ndim[0]=nyip;
     ndim[1]=nxip;
     nmaskp=ncvardef(ncidparam,"ref_point",NC_INT,2,ndim);   
     ndim[0]=ntip;
     ndim[1]=nv;
     ndim[2]=npip;
     varanondvi = ncvardef(ncidparam,"anondvi",NC_FLOAT,3,ndim);
     varanandvi = ncvardef(ncidparam,"anandvi",NC_FLOAT,2,ndim+1);
     varanandviav = ncvardef(ncidparam,"anandviav",NC_FLOAT,2,ndim+1);
     varanandviap = ncvardef(ncidparam,"anandviap",NC_FLOAT,2,ndim+1);
     varanandvim = ncvardef(ncidparam,"anondvim",NC_FLOAT,2,ndim+1);
     varndvim=ncvardef(ncidparam,"ndvim",NC_FLOAT,2,ndim+1);
     varndvimin=ncvardef(ncidparam,"ndvimin",NC_FLOAT,2,ndim+1);
     varndvimax=ncvardef(ncidparam,"ndvimax",NC_FLOAT,2,ndim+1);
     varndviap=ncvardef(ncidparam,"ndviap",NC_FLOAT,2,ndim+1);
     vveg = ncvardef(ncidparam,"veg",NC_FLOAT,2,ndim+1);
     varyearmin= ncvardef(ncidparam,"yearmin",NC_FLOAT,2,ndim+1);
      ndim[0]=ntip;
     ndim[1]=npip;
     varanoprec = ncvardef(ncidparam,"anoprec",NC_FLOAT,2,ndim);
     varanots = ncvardef(ncidparam,"anots",NC_FLOAT,2,ndim);
      vartsmax = ncvardef(ncidparam,"tsmax",NC_FLOAT,2,ndim);
     varanotsta = ncvardef(ncidparam,"anotsta",NC_FLOAT,2,ndim);
     varanats = ncvardef(ncidparam,"anats",NC_FLOAT,1,ndim+1);
      varanatsta = ncvardef(ncidparam,"anatsta",NC_FLOAT,1,ndim+1);
 
     varminp= ncvardef(ncidparam,"mindifprec",NC_FLOAT,1,ndim+1);
     varanoprecav= ncvardef(ncidparam,"anoprecav",NC_FLOAT,1,ndim+1);
     varanoprecap= ncvardef(ncidparam,"anoprecap",NC_FLOAT,1,ndim+1);
     varmeanp= ncvardef(ncidparam,"meandifprec",NC_FLOAT,1,ndim+1);
     vardmax= ncvardef(ncidparam,"dmax",NC_FLOAT,1,ndim+1);

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
     ncattput(ncidparam, vartsmax, "_FillValue",NC_FLOAT,1,&fillval);
      ncattput(ncidparam, varanotsta, "_FillValue",NC_FLOAT,1,&fillval);
     ncattput(ncidparam, varanats, "_FillValue",NC_FLOAT,1,&fillval);
     printf ("ncat 1\n");
     ncattput(ncidparam, varanatsta, "_FillValue",NC_FLOAT,1,&fillval);
    ncattput(ncidparam, varanandvi, "_FillValue",NC_FLOAT,1,&fillval);
   ncattput(ncidparam, vveg, "_FillValue",NC_FLOAT,1,&fillval);
       ncattput(ncidparam, varanandviav, "_FillValue",NC_FLOAT,1,&fillval);
      ncattput(ncidparam, varanandviap, "_FillValue",NC_FLOAT,1,&fillval);
     ncattput(ncidparam, varanandvim, "_FillValue",NC_FLOAT,1,&fillval);
     printf ("ncat 2\n");
     ncattput(ncidparam, varminp, "_FillValue",NC_FLOAT,1,&fillval);
     ncattput(ncidparam, varmeanp, "_FillValue",NC_FLOAT,1,&fillval);
     ncattput(ncidparam, vardmax, "_FillValue",NC_FLOAT,1,&fillval);
     ncattput(ncidparam, varyearmin, "_FillValue",NC_FLOAT,1,&fillval);
    ncattput(ncidparam, varndvim, "_FillValue",NC_FLOAT,1,&fillval);
    ncattput(ncidparam, varndvimin, "_FillValue",NC_FLOAT,1,&fillval);
    ncattput(ncidparam, varndvimax, "_FillValue",NC_FLOAT,1,&fillval);
    ncattput(ncidparam, varndviap, "_FillValue",NC_FLOAT,1,&fillval);
     ncattput(ncidparam, varndviap, "_FillValue",NC_FLOAT,1,&fillval);
   printf ("ncat 3\n");
     printf ("ncat 4\n");
     ncendef(ncidparam);
    }
  sprintf (nfic,"output/crujrad_%d.nc",yearref);
 ncidmean=ncopen(nfic,NC_NOWRITE);
     nc_inq_varid(ncidmean,"LAI",&nvarvim);
     nc_inq_varid(ncidmean,"Raint",&nvarpcm);
     nc_inq_varid(ncidmean,"tsol_max",&nvartsm);
     nc_inq_varid(ncidmean,"tair",&nvartstam);
    nc_inq_varid(ncidmean,"maxvegetfrac",&nvarveg);
  nc_inq_varid(ncidmean,"longitude",&nvarlonvi);
  nc_inq_varid(ncidmean,"latitude",&nvarlatvi);
  nc_inq_varid(ncidmean,"ref_point",&nvarmask);
  nc_inq_varid(ncidmean,"indice_ligne",&nvaridligne);
  nc_inq_varid(ncidmean,"indice_col",&nvaridcol);
  start[0]=0;
  ncount[0]=NBCTS;
  printf("varp0\n");
  ncvarget(ncidmean,nvarlonvi,start,ncount,datalonvi);
  printf("varp00\n");
 ncvarput(ncidparam,nlonp,start,ncount,datalonvi);	  
  printf ("varp1\n");
  ncount[0]=NBLTS;
  ncvarget(ncidmean,nvarlatvi,start,ncount,datalatvi);
  ncvarput(ncidparam,nlatp,start,ncount,datalatvi);	  
   printf ("varp1\n");
 ncount[0]=NBPOINTTERRE;
  ncvarget(ncidmean,nvaridligne,start,ncount,dataidligne);
  ncvarget(ncidmean,nvaridcol,start,ncount,dataidcol);
  ncvarput(ncidparam,nidlignep,start,ncount,dataidligne);	  
  printf ("varp1\n");
  ncvarput(ncidparam,nidcolp,start,ncount,dataidcol);	  
  printf ("varp1\n");
  start[0]=start[1]=0;
  ncount[0]=NBLTS;
  ncount[1]=NBCTS;
  ncvarget(ncidmean,nvarmask,start,ncount,datamask);
  ncvarput(ncidparam,nmaskp,start,ncount,datamask);	  
  printf ("varp1\n");
     start[0]=0;
     ncount[0]=365;
     start[1]=0;
     ncount[1]=NVEG;
     start[2]=0;
     ncount[2]=NBPOINTTERRE;
        printf ("varvim=%d\n",nvarvim);
     status=ncvarget(ncidmean,nvarvim,start,ncount,datavim);
     // ncount[0]=DYEAR-PYEAR+1;
     printf ("lec ndvim\n");
     //printf ("%d\n",nvarpcm);
     start[0]=0;
   ncount[0]=365;
     start[1]=0;
     ncount[1]=NBPOINTTERRE;
      status=ncvarget(ncidmean,nvarpcm,start,ncount,datapcm);  
     printf ("lec datapcm\n");
     // printf ("%d\n",nvartsm);
      status=ncvarget(ncidmean,nvartsm,start,ncount,datatsm);
    printf ("lec datatsm\n");
     start[0]=0;
     ncount[0]=1;
     ncount[1]=NVEG;
     start[1]=0;
     ncount[2]=NBPOINTTERRE;
     start[2]=0;
     status=ncvarget(ncidmean,nvarveg,start,ncount,dataveg);
    for (i=0;i<NBPOINTTERRE;i++)
       {
	 for (t=7,vtmax[i]=-1000;t<365;t++) 
	   {
	     if (datatsm[t*NBPOINTTERRE+i]-TP0>vtmax[i]) 
	       {
		 datadmax[i]=t;
		 vtmax[i]=datatsm[t*NBPOINTTERRE+i]-TP0;
	       
	       }
	   }
       }
   for (i=0;i<NBPOINTTERRE;i++)for(v=0;v<NVEG;v++)
     {
       minanoprec[v][i]=1e10;
       ndvim[v*NBPOINTTERRE+i]=0;
     }
   for (i=0;i<NBPOINTTERRE;i++) for (year=0;year<DYEAR-PYEAR+1;year++) datatsmax[(year)*NBPOINTTERRE+i]=-1000;

 for (year=PYEAR;year<=DYEAR;year++)
   {
     printf (" **************** %d **********************\n",year);
     sprintf (nfic,"output/crujrad_%d.nc",year);
     printf ("fichier %s\n",nfic);
     ncidndvi=ncopen(nfic,NC_NOWRITE);
     printf ("ncidndvi=%d\n",ncidndvi);
     nc_inq_varid(ncidndvi,"LAI",&nvarvi);
     nc_inq_varid(ncidndvi,"Raint",&nvarpc);
     nc_inq_varid(ncidndvi,"tsol_max",&nvarts);
     nc_inq_varid(ncidndvi,"tair",&nvartsta);
     memset(nbvalndvi,0,NVEG*NBPOINTTERRE*sizeof(int));
     memset(nbvalndviavant,0,NVEG*NBPOINTTERRE*sizeof(int));
     memset(nbvalndviapres,0,NVEG*NBPOINTTERRE*sizeof(int));
     memset(nbvalprec,0,NBPOINTTERRE*sizeof(int));
     start[0]=0;
     ncount[0]=365;
     start[1]=0;
     ncount[1]=NBPOINTTERRE;
     status=ncvarget(ncidndvi,nvarpc,start,ncount,datapc);  
     status=ncvarget(ncidndvi,nvarts,start,ncount,datats);
     start[0]=0;
     ncount[0]=365;
     start[1]=0;
     ncount[1]=NVEG;
     start[2]=0;
     ncount[2]=NBPOINTTERRE;
         status=ncvarget(ncidndvi,nvarvi,start,ncount,datandvi);
       //     status=ncvarget(ncidndvi,nvatsta,start,ncount,datata);  
    for (i=0;i<NBPOINTTERRE;i++) 
        {
	  lat=datalatvi[dataidligne[i]];
	  if (lat>0) {t1=60;t2=240;}
	  else {t1=240;t2=424;}
	  for (t=2;t<365;t++) 
	    {
	     for (v=0;v<NVEG;v++)
	       if (dataveg[v*NBPOINTTERRE+i]>0)
	       {
		 if (((t1==60)&&(t>t1)&&(t<t2))||(t1==240)&&((t<60)||(t>240))) 
		   {
		     anondvi[(year-PYEAR)*NBPOINTTERRE*NVEG+v*NBPOINTTERRE+i]+=(fpar(datavim[t*NBPOINTTERRE*NVEG+v*NBPOINTTERRE+i])-fpar(datandvi[t*NBPOINTTERRE*NVEG+v*NBPOINTTERRE+i]));
		     //    printf ("( %f %e)",dataveg[v*NBPOINTTERRE+i],datandvi[t*NBPOINTTERRE*NVEG+v*NBPOINTTERRE+i]); 
		     if (year==yearref)
		       {
			 ndvim[v*NBPOINTTERRE+i]+=fpar(datandvi[t*NBPOINTTERRE*NVEG+v*NBPOINTTERRE+i]);
			 ndvimin[v*NBPOINTTERRE+i]=(datavim[v*NBPOINTTERRE+i]<ndvimin[v*NBPOINTTERRE+i])?fpar(datavim[v*NBPOINTTERRE+i]):ndvimin[v*NBPOINTTERRE+i];
			 ndvimax[v*NBPOINTTERRE+i]=(datavim[v*NBPOINTTERRE+i]>ndvimax[v*NBPOINTTERRE+i])?fpar(datavim[v*NBPOINTTERRE+i]):ndvimax[v*NBPOINTTERRE+i];
		       }
		     nbvalndvi[v][i]++;
		   }
		 if ((t<datadmax[i])&&(t>=datadmax[i]-30))
		 {
		   anondviavant[(year-PYEAR)*NBPOINTTERRE*NVEG+v*NBPOINTTERRE+i]+=(fpar(datavim[t*NBPOINTTERRE*NVEG+v*NBPOINTTERRE+i])-fpar(datandvi[t*NBPOINTTERRE*NVEG+v*NBPOINTTERRE+i]));
		   nbvalndviavant[v][i]++;
		 }
		 else
		   if ((t>datadmax[i])&&(t<=datadmax[i]+30))
		     {
		       anondviapres[(year-PYEAR)*NBPOINTTERRE*NVEG+v*NBPOINTTERRE+i]+=(fpar(datavim[t*NBPOINTTERRE*NVEG+v*NBPOINTTERRE+i])-fpar(datandvi[t*NBPOINTTERRE*NVEG+v*NBPOINTTERRE+i]));
		       nbvalndviapres[v][i]++;
		       if (year==yearref) ndviap[v*NBPOINTTERRE+i]+=fpar(datavim[t*NBPOINTTERRE*NVEG+v*NBPOINTTERRE+i]);
		     }	       
	       }
	     if ((datats[t*NBPOINTTERRE+i]-TP0>datatsmax[(year-PYEAR)*NBPOINTTERRE+i])) 
	       {
		 datatsmax[(year-PYEAR)*NBPOINTTERRE+i]=datats[t*NBPOINTTERRE+i]-TP0;
	       }
	    }
	 /*	 for (i=0;i<NBPOINTTERRE;i++)
	     {
	       anoprec[(year-PYEAR)*NBPOINTTERRE+i]+=(datapcm[t*NBPOINTTERRE+i]-datapc[t*NBPOINTTERRE+i]);
	       nbvalprec[i]++;
	       if (t<datadmax[i])
		 {
		   anoprecavant[(year-PYEAR)*nbl*NBCTS+l*NBCTS+c]+=(datapcm[t*NBPOINTTERRE+i]-datapc[t*NBPOINTTERRE+i]);
		   nbvalprecavant[i]++;
		 }
	       else
		 {
		   anoprecapres[(year-PYEAR)*nbl*NBCTS+l*NBCTS+c]+=(datapcm[t*NBPOINTTERRE+i]-datapc[t*NBPOINTTERRE+i]);
		   nbvalprecapres[i]++;
		 }	                
		 }*/
	}
     for (i=0;i<NBPOINTTERRE;i++) if ( datatsmax[(year-PYEAR)*NBPOINTTERRE+i]<-100)  datatsmax[(year-PYEAR)*NBPOINTTERRE+i]=0;
     for (i=0;i<NBPOINTTERRE;i++)
	 for (v=0;v<NVEG;v++)
	   {
	     if (nbvalndvi[v][i]>0) 
	       {
		 anondvi[(year-PYEAR)*NVEG*NBPOINTTERRE+v*NBPOINTTERRE+i]/=nbvalndvi[v][i];
		 if (year==yearref) ndvim[v*NBPOINTTERRE+i]/=nbvalndvi[v][i];
	       }
	     if (nbvalndviavant[v][i]>0) 
	       {
		 anondviavant[(year-PYEAR)*NVEG*NBPOINTTERRE+v*NBPOINTTERRE+i]/=nbvalndviavant[v][i];
	       }
	     if (nbvalndviapres[v][i]>0) 
	       {
		 anondviapres[(year-PYEAR)*NVEG*NBPOINTTERRE+v*NBPOINTTERRE+i]/=nbvalndviapres[v][i];
		 if (year==yearref) ndviap[v*NBPOINTTERRE+i]/=nbvalndviapres[v][i];

	       }
	   }
     for (i=0;i<NBPOINTTERRE;i++)
       {
	 if (nbvalprec[i]>0) anoprec[(year-PYEAR)*NBPOINTTERRE+i]/=nbvalprec[i];
	 else anoprec[(year-PYEAR)]=fillval;
	 if (nbvalprecavant[i]>0) anoprecavant[(year-PYEAR)*NBPOINTTERRE+i]/=nbvalprecavant[i];
	 else anoprecavant[(year-PYEAR)*NBPOINTTERRE+i]=fillval;
	 if (nbvalprecapres[i]>0) anoprecapres[(year-PYEAR)*NBPOINTTERRE+i]/=nbvalprecapres[i];
	 else anoprecapres[(year-PYEAR)*NBPOINTTERRE+i]=fillval;
       }
   }

 for (i=0;i<NBPOINTTERRE;i++) 
   {
     anondvimean[i]=anoprecmean[i]=nbvmints[i]=0;
     for (v=0;v<NVEG;v++) nbvmin[v][i]=0;
     mints[i]= mintsta[i]= minanoprecr[i]= anondvimin[i]= anondviminavant[i]=anondviminapres[i]=0;
     //     ndvim[i]/=nbvalprec[i];
     for (year=PYEAR;year<=DYEAR;year++)
       {
       if (year!=yearref)
	 {
	   if (datadmax[i]>=0)
	     {
	       //	       if (year>=2011)
		 {
		   anots[(year-PYEAR)*NBPOINTTERRE+i]=datatsmax[(yearref-PYEAR)*NBPOINTTERRE+i]-datatsmax[(year-PYEAR)*NBPOINTTERRE+i];
		   if (anots[(year-PYEAR)*NBPOINTTERRE+i]<-50) anots[(year-PYEAR)*NBPOINTTERRE+i]=0;
		 }
	       for (v=0;v<NVEG;v++)
		 {
		   //		   if (fabs(anondvi[(year-PYEAR)*NVEG*NBPOINTTERRE+v*NBPOINTTERRE+i]<minanoprec[v][i])&&(anondvi[(year-PYEAR)*NVEG*NBPOINTTERRE+v]>fillval)&&(fabs(anondvi[(year-PYEAR)*NVEG*NBPOINTTERRE+v*NBPOINTTERRE+i])<0.1))
		   //		   if ((anondvi[(year-PYEAR)*NVEG*NBPOINTTERRE+v]>fillval)&&(fabs(anondvi[(year-PYEAR)*NVEG*NBPOINTTERRE+v*NBPOINTTERRE+i])<0.1))
		   if (fabs(anondvi[(year-PYEAR)*NVEG*NBPOINTTERRE+v*NBPOINTTERRE+i]<minanoprec[v][i])&&(anondvi[(year-PYEAR)*NVEG*NBPOINTTERRE+v]>fillval))
		     {
		       minanoprec[v][i]=fabs(anondvi[(year-PYEAR)*NBPOINTTERRE*NVEG+v*NBPOINTTERRE+i]);
		       minanoprecr[i]=anoprec[(year-PYEAR)*NBPOINTTERRE+i];
		       mints[i]+=anots[(year-PYEAR)*NBPOINTTERRE+i];
		       mintsta[i]+=anotsta[(year-PYEAR)*NBPOINTTERRE+i]; 
		       nbvmin[v][i]++;
		       /*		       anondvimin[v*NBPOINTTERRE+i]=+anondvi[(year-PYEAR)*NBPOINTTERRE*NVEG+v*NBPOINTTERRE+i]; 
		       anondviminavant[v*NBPOINTTERRE+i]=+anondviavant[(year-PYEAR)*NBPOINTTERRE*NVEG+v*NBPOINTTERRE+i]; 
		       anondviminapres[v*NBPOINTTERRE+i]=+anondviapres[(year-PYEAR)*NBPOINTTERRE*NVEG+v*NBPOINTTERRE+i]; */
		       anondvimin[v*NBPOINTTERRE+i]=anondvi[(year-PYEAR)*NBPOINTTERRE*NVEG+v*NBPOINTTERRE+i]; 
		       anondviminavant[v*NBPOINTTERRE+i]=anondviavant[(year-PYEAR)*NBPOINTTERRE*NVEG+v*NBPOINTTERRE+i]; 
		       anondviminapres[v*NBPOINTTERRE+i]=anondviapres[(year-PYEAR)*NBPOINTTERRE*NVEG+v*NBPOINTTERRE+i]; 
		       
		     }
		 }
	       anoprecmean[i]+=anoprec[(year-PYEAR)*NBPOINTTERRE+i]/(float)(DYEAR-PYEAR);
	     }
	   v=0;
	   for (nbvv[v][i]=0;v<NVEG;v++) 
	     if (anondvi[(year-PYEAR)*NVEG*NBPOINTTERRE+v*NBPOINTTERRE+i]>fillval)
	       {
		 anondvimean[v*NBPOINTTERRE+i]+=anondvi[(year-PYEAR)*NVEG*NBPOINTTERRE+v*NBPOINTTERRE+i];
		 nbvv[v][i]++;
	       }
	 }
       }
   }
   for (i=0;i<NBPOINTTERRE;i++) 
     {
       if (nbvmints[i]>0)
	 {
	   mints[i]/=nbvmints[i];
	   mintsta[i]/=nbvmints[i];
	 }
       for (v=0;v<NVEG;v++)
	 if (nbvmin[v][i]>0)
	   {
	     anondvimin[v*NBPOINTTERRE+i]/=nbvmin[v][i];
	     anondviminavant[v*NBPOINTTERRE+i]/=nbvmin[v][i];
	     anondviminapres[v*NBPOINTTERRE+i]/=nbvmin[v][i];
	   }
    for (v=0;v<NVEG;v++) if (nbvv[v][i]>0) {anondvimean[v*NBPOINTTERRE+i]/=nbvv[v][i];} else anondvimean[v*NBPOINTTERRE+i]=fillval;
      }
   /*   for (i=0;i<NBPOINTTERRE;i++) 
     {
       for (v=0;v<NVEG;v++)
	 if (dataveg[v*NBPOINTTERRE+i]<0.1)
	   {
	     anondvimean[v*NBPOINTTERRE+i]=anondvimin[v*NBPOINTTERRE+i]=anondviminavant[v*NBPOINTTERRE+i]=anondviminapres[v*NBPOINTTERRE+i]=yearmin[v*NBPOINTTERRE+i]=fillval;
	     ndvim[v*NBPOINTTERRE+i]=ndviap[v*NBPOINTTERRE+i]=fillval;
	     for (y=0;y<(DYEAR-PYEAR+1);y++) anondvi[(year-PYEAR)*NVEG*NBPOINTTERRE+v*NBPOINTTERRE+i]=fillval;
	   }
	   }*/
   for (t=0;t<DYEAR-PYEAR+1;t++) time[t]=t*365+182;


   start[0]=0;
   start[1]=0;
   ncount[0]=DYEAR-PYEAR+1;
   ncount[1]=NBPOINTTERRE;
   printf ("%d %d %d %d %d %d\n",start[0],start[1],start[2],ncount[0],ncount[1],ncount[2]);
   ncvarput(ncidparam,varanoprec,start,ncount,anoprec);
   printf ("put 1\n");
   printf ("put 2n");
   ncvarput(ncidparam,varanots,start,ncount,anots);
   printf ("put 3\n");
  ncvarput(ncidparam,vartsmax,start,ncount,datatsmax);
   ncvarput(ncidparam,varanotsta,start,ncount,anotsta);
  ncvarput(ncidparam,ntimep,start,ncount,time);
    start[0]=0;
  ncount[0]=NBPOINTTERRE;
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
    status=ncvarput(ncidparam,vardmax,start,ncount,datadmax);
   ERRHAND(status);
  printf ("put 8\n");
   ncount[0]=NVEG;
   ncount[1]=NBPOINTTERRE;
   start[1]=start[0]=0;
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
   ncount[2]=NBPOINTTERRE;
   start[2]=start[1]=start[0]=0;
   ncvarput(ncidparam,varanondvi,start,ncount,anondvi);
  status=ncclose(ncidparam);
  printf ("fin ecriture, on envoie au suivant\n");
  // MPI_Finalize();
}
     

      
