#include <stdio.h>
#include <stdlib.h>
#include<netcdf.h>
#include <dirent.h> 
#include <unistd.h>
#include <fnmatch.h>
#include<string.h>
#include<math.h>
#include "/ccc/cont003/home/drf/p529viov/nr/nr.h"
#include "/ccc/cont003/home/drf/p529viov/nr/nrutil.h"

#define PYEAR 2011
#define DYEAR 2019
#define NBLTS 3584
#define NBCTS 8064
#define NBLVEG 7200
#define NBCVEG 14400
#define NBVAL 160
#define NBVAL2 180
#define NBVAL3 300
#define NBVAL4 190
#define NBVAL5 300
#define NBVALTS 120
#define NVEG 7
#define NPFT 15
#define NBYEAR 9
#define NBL  128
#define FACVEG 1.785714328
#define DTS 22.399999466
#define DIVWC 120
#define NBVALSEUIL 697805
#define ERRHAND(e) if ((e) !=NC_NOERR) printf ("%s\n",nc_strerror((e)))
float danandvi[NVEG][NBL][NBCTS], danats[NBL][NBCTS],danandviav[NVEG][NBL][NBCTS],danandviap[NVEG][NBL][NBCTS],dndvim[NVEG][NBL][NBCTS],dndviap[NVEG][NBL][NBCTS],datatsmax[DYEAR-PYEAR+1][NBL][NBCTS],datatstamax[NBL][NBCTS],datatstamax[NBL][NBCTS],meanmaxts[NBL][NBCTS],tsmseuil[NPFT][NBL][NBCTS],datavego[NPFT][NBL][NBCTS],nbvalseuil[NPFT][NBL][NBCTS],year1[NPFT][NBL][NBCTS],year2[NPFT][NBL][NBCTS],diffout[NPFT][NBL][NBCTS],anaav[NPFT][NBL][NBCTS],nbvaldiff[NPFT][NBL][NBCTS],datatamean[NBL][NBCTS];
float nvalmm[NBL][NBCTS],nbptot;
float hist2D[NPFT][NBVAL2][NBVAL2];
float histseuil[NBVAL4][NBVAL5];
int nval2D[NPFT][NBVAL2][NBVAL2];
float xval[NBVAL5],yval[NBVAL4];
float *xvseuil,*yvseuil;
float histtamseuil[NPFT][NBVAL3], histtaeseuil[NPFT][NBVAL3],nbvaltamseuil[NPFT][NBVAL3];
void main(int argc,char *argv[])
{
  FILE *fout[NPFT],*foutts[NPFT],*foutseuil,*fouttt[NPFT],*fouttt2[NPFT];
  int ncidparam[NBYEAR],ncidts,ncidveg,status,nxip,nyip,ntip,nlonp,nlap,ntimep,varanondvi,ncido,nvarlat,nvarlon,nlatp,ndim[4],v2,vardiff,varveget,varvego,vartamean,ncidtamean;
  int varanandvi[NBYEAR],varanandviav[NBYEAR],varanandviap[NBYEAR],vartsmax,vartstamax,varndvim[NBYEAR],varanats[NBYEAR],var,nvarlonts,nvarlatts,varmtmax,nxii,nyii,nxiii,nyiii,nh2D,nxv,nyv,vartseuil,varnbyear,varyear1,varyear2,nveg,varndviap[NBYEAR],varanaav;
  float diffvi,tsm,histom[NPFT][NBVAL2],histts[NPFT][NBVALTS],varm[NPFT][NBVAL2],histom2[NPFT][NBVAL2],histom3[NPFT][NBVAL2];
  float *dataveg;
  int nvalts[NPFT][NBVALTS],year,i,l,c,t,v,l2,nblveg,nbvalv[NPFT],l3,c2,lveg,cveg;
  float nval[NPFT][NBVAL2];
  long ncount[4],start[4];
  char nfic[100],str70[70];
  float fillval=-1e32,ats,fillval2=-1e20,tsmm,lat,hist[NPFT];
  double meantseuil[NPFT],etseuil[NPFT],nbvaltseuil[NPFT];
  double datalon[NBCTS],datalat[NBLTS],time[10];
  char ficout[20];
  int nbvalxy=0;
  float a,b,siga,sigb,chi2,q,r,prob,z;
  int corre[NPFT]={6,2,3,4,2,3,4,2,5,2,2,2,1,1,1};
  memset (histseuil,0,NBVAL4*NBVAL5*sizeof(float));
  memset (hist2D,0,NBVAL2*NBVAL2*NPFT*sizeof(float));
  memset (nval2D,0,NBVAL2*NBVAL2*NPFT*sizeof(int));
  dataveg=malloc(NPFT*NBL*FACVEG*NBCVEG*sizeof(float));
  for (v=0;v<NPFT;v++) meantseuil[NPFT]=etseuil[NPFT]=nbvaltseuil[v]=0;
 ncidts=ncopen("data/daytmaxyglob.nc",NC_NOWRITE);
  status=nc_inq_varid(ncidts,"tsta_max",&vartstamax);
  status=nc_inq_varid(ncidts,"tsmax",&vartsmax);
  ncidveg=ncopen("data/PFTmap_2015.nc",NC_NOWRITE);
  status=nc_inq_varid(ncidveg,"maxvegetfrac",&varveget);
  ncidtamean=ncopen("data/meantemp.nc",NC_NOWRITE);
  status=nc_inq_varid(ncidtamean,"tmean",&vartamean);
  xvseuil=malloc(NBVALSEUIL*sizeof(float));
  yvseuil=malloc(NBVALSEUIL*sizeof(float));
 ERRHAND(status);
for (v=0;v<NPFT;v++)
   {
     sprintf (ficout,"histtndvig_%d",v+1);
     fout[v]=fopen(ficout,"w");
     sprintf (ficout,"histtstndvig_%d",v+1);
     foutts[v]=fopen(ficout,"w");
     sprintf (ficout,"teuilvst_%d",v+1);
     fouttt[v]=fopen(ficout,"w");
     sprintf (ficout,"tseuilmt_%d",v+1);
     fouttt2[v]=fopen(ficout,"w");

     for (i=0;i<NBVAL2;i++) histom[v][i]=histom2[v][i]=histom3[v][i]=varm[v][i]=nval[v][i]=0;
     for (i=0;i<NBVALTS;i++) histts[v][i]=nvalts[v][i]=0;
   } 
 foutseuil=fopen("histoseuil","w");
 ncido=nccreate("data/repndvipftglob.nc",NC_64BIT_OFFSET);
 nxip = ncdimdef(ncido,"lon",NBCTS);
 nyip = ncdimdef(ncido,"lat",NBLTS);
 ntip = ncdimdef(ncido,"time",NC_UNLIMITED);
 nveg = ncdimdef(ncido,"veg",NPFT);
 nlonp = ncvardef(ncido,"lon",NC_DOUBLE,1,&nxip);
 nlatp = ncvardef(ncido,"lat",NC_DOUBLE,1,&nyip);
 ntimep = ncvardef(ncido,"time",NC_DOUBLE,1,&ntip);
 nxii= ncdimdef(ncido,"TSmax",NBVAL2);
 nyii= ncdimdef(ncido,"TSmax_mean",NBVAL2);
 nxiii= ncdimdef(ncido,"Ta_mean",NBVAL5);
 nyiii= ncdimdef(ncido,"Tcrit",NBVAL4);
nxv=ncvardef(ncido,"Ta_mean",NC_FLOAT,1,&nxiii);
nyv=ncvardef(ncido,"Tcrit",NC_FLOAT,1,&nyiii);

ndim[0]=nyiii;
 ndim[1]=nxiii;
 nh2D=ncvardef(ncido,"histoseuil",NC_FLOAT,2,ndim);
 ndim[0]=ntip;
 ndim[1]=nveg;
 ndim[2]=nyip;
 ndim[3]=nxip;
 varanondvi = ncvardef(ncido,"anondvi",NC_FLOAT,4,ndim);
 varmtmax = ncvardef(ncido,"mtmax",NC_FLOAT,2,ndim+2);
 vartseuil = ncvardef(ncido,"tseuil",NC_FLOAT,3,ndim+1);
 varvego = ncvardef(ncido,"veget",NC_FLOAT,3,ndim+1);
 varnbyear = ncvardef(ncido,"nbyear",NC_FLOAT,3,ndim+1);
 vardiff = ncvardef(ncido,"diff",NC_FLOAT,3,ndim+1);
 varanaav = ncvardef(ncido,"anav",NC_FLOAT,3,ndim+1);
 varyear1 = ncvardef(ncido,"year1",NC_FLOAT,2,ndim+2);
 varyear2 = ncvardef(ncido,"year2",NC_FLOAT,2,ndim+2);
 ncattput(ncido, varanondvi, "_FillValue",NC_FLOAT,1,&fillval);
 ncattput(ncido, varmtmax, "_FillValue",NC_FLOAT,1,&fillval);
 ncattput(ncido, vartseuil, "_FillValue",NC_FLOAT,1,&fillval);
 ncattput(ncido, varnbyear, "_FillValue",NC_FLOAT,1,&fillval);
 ncattput(ncido, varyear1, "_FillValue",NC_FLOAT,1,&fillval);
 ncattput(ncido, vardiff, "_FillValue",NC_FLOAT,1,&fillval);
 ncattput(ncido, varanaav, "_FillValue",NC_FLOAT,1,&fillval);
  ncattput(ncido, varyear2, "_FillValue",NC_FLOAT,1,&fillval);
  ncattput(ncido, varvego, "_FillValue",NC_FLOAT,1,&fillval);
 sprintf (str70,"days since %04d-%02d-%02d %02d:%02d:%02d",2011,1,1,0,0,0);
 ncattput(ncido,ntimep, "units",NC_CHAR, strlen(str70), str70);
 ncattput(ncido,ntimep, "axis",NC_CHAR, 1, "T");
 status=ncattput (ncido,nlonp,"units",NC_CHAR,12,"degrees_east");
 status=ncattput (ncido,nlonp,"axis",NC_CHAR,1,"X");
 status=ncattput (ncido,nlatp,"units",NC_CHAR,13,"degrees_north");
 status=ncattput (ncido,nlatp,"axis",NC_CHAR,1,"Y");
 ncendef(ncido);
  for (i=0;i<NBVAL4;i++) yval[i]=i/4.+20.;
     start[0]=0;
 ncount[0]=NBVAL4;
  ncvarput(ncido,nyv,start,ncount,yval);
 ncount[0]=NBVAL5;
  for (i=0;i<NBVAL5;i++) xval[i]=((float)i)/10.+5.;
  ncvarput(ncido,nxv,start,ncount,xval);
 nc_inq_varid(ncidts,"lon",&nvarlonts);
 nc_inq_varid(ncidts,"lat",&nvarlatts);
 start[0]=0;
 ncount[0]=NBLTS;
 ncvarget(ncidts,nvarlatts,start,ncount,datalat);
 ncvarput(ncido,nlatp,start,ncount,datalat);	  
 start[0]=0;
 ncount[0]=NBCTS;
 ncvarget(ncidts,nvarlonts,start,ncount,datalon);
 ncvarput(ncido,nlonp,start,ncount,&datalon);
 for (year=PYEAR;year<=DYEAR;year++)
   {
     sprintf (nfic,"data/analogtpftglob_%d.nc",year);
     ncidparam[year-PYEAR]=ncopen(nfic,NC_NOWRITE);
     status=nc_inq_varid(ncidparam[year-PYEAR],"anats",&varanats[year-PYEAR]);
     status=nc_inq_varid(ncidparam[year-PYEAR],"anandvi",&varanandvi[year-PYEAR]);
     status=nc_inq_varid(ncidparam[year-PYEAR],"anandviav",&varanandviav[year-PYEAR]);
     status=nc_inq_varid(ncidparam[year-PYEAR],"anandviap",&varanandviap[year-PYEAR]);
     status=nc_inq_varid(ncidparam[year-PYEAR],"ndvim",&varndvim[year-PYEAR]);
     status=nc_inq_varid(ncidparam[year-PYEAR],"ndviap",&varndviap[year-PYEAR]);
   }
 l=0;
 for (l=0;l<NBLTS;l+=NBL)
   {
     printf ("l=%d\n",l);
     memset (tsmseuil,0,NBL*NBCTS*NPFT*sizeof(int));
     memset (nbvalseuil,0,NBL*NBCTS*NPFT*sizeof(int));
     memset (year1,0,NBL*NBCTS*NPFT*sizeof(int));
     memset (year2,0,NBL*NBCTS*NPFT*sizeof(int));
     for (v=0;v<NPFT;v++) for(l2=0;l2<NBL;l2++) for (c=0;c<NBCTS;c++) {diffout[v][l2][c]=anaav[v][l2][c]=0;nbvaldiff[v][l2][c]=0;}
     start[0]=PYEAR-2011;
     start[1]=l;
       start[2]=0;
     ncount[0]=DYEAR-PYEAR+1;
     ncount[1]=NBL;
     ncount[2]=NBCTS;
     status=ncvarget(ncidts,vartsmax,start,ncount,datatsmax); 
     start[0]=start[1]=0;
     start[2]=((float)l+10.*DTS)*FACVEG;
     start[3]=0;
     ncount[0]=1;
     ncount[1]=15;
     ncount[2]=NBL*FACVEG;
     if (start[2]+ncount[2]>NBLVEG) ncount[2]=NBLVEG-start[2];
     nblveg=ncount[2];
     ncount[3]=NBCVEG;
     ERRHAND(status);
     status=ncvarget(ncidveg,varveget,start,ncount,dataveg);      
      ncount[0]=NBL;
     ncount[1]=NBCTS;
    start[0]=l;
       start[1]=0;
     status=ncvarget(ncidtamean,vartamean,start,ncount,datatamean);      

    for (l2=0;l2<NBL;l2++)
       for (c=0;c<NBCTS;c++) meanmaxts[l2][c]=nvalmm[l2][c]=0;
     for (year=PYEAR;year<=DYEAR;year++)
       for (l2=0;l2<NBL;l2++)
	 for (c=0;c<NBCTS;c++)
	   if (datatsmax[year-PYEAR][l2][c]>fillval2) 
	     {
	       //	   printf ("(%d %e) ",year-PYEAR,datatsmax[year-PYEAR][l][c]);
	       meanmaxts[l2][c]+=datatsmax[year-PYEAR][l2][c];
	       nvalmm[l2][c]++;
	     }
     for (l2=0;l2<NBL;l2++)
       for (c=0;c<NBCTS;c++)
	 if (nvalmm[l2][c]>0) meanmaxts[l2][c]/= nvalmm[l2][c]; else  meanmaxts[l2][c]=fillval;
     for (year=PYEAR;year<=DYEAR;year++)
       {
	 start[0]=l;
	 ncount[0]=NBL;
	 start[1]=0;
	 ncount[1]=NBCTS;
	 status=ncvarget(ncidparam[year-PYEAR],varanats[year-PYEAR],start,ncount,danats);
	 ERRHAND(status);
	 start[0]=0;
	 ncount[0]=NVEG;
	 start[1]=l;
	 ncount[1]=NBL;
	 start[2]=0;
	 ncount[2]=NBCTS;
	 status=ncvarget(ncidparam[year-PYEAR],varanandvi[year-PYEAR],start,ncount,danandvi);
	 ERRHAND(status);
	 ERRHAND(status);
	 status=ncvarget(ncidparam[year-PYEAR],varanandviav[year-PYEAR],start,ncount,danandviav);
	 ERRHAND(status);
	 ERRHAND(status);
	 status=ncvarget(ncidparam[year-PYEAR],varanandviap[year-PYEAR],start,ncount,danandviap);
	 ERRHAND(status);
	 ERRHAND(status);
	 status=ncvarget(ncidparam[year-PYEAR],varndvim[year-PYEAR],start,ncount,dndvim);
	 ERRHAND(status);
	 status=ncvarget(ncidparam[year-PYEAR],varndviap[year-PYEAR],start,ncount,dndviap);
	 start[0]=year-2011;
	 start[1]=l;
	 start[2]=0;
	 ncount[0]=1;
	 ncount[1]=NBL;
	 ncount[2]=NBCTS;
         status=ncvarget(ncidts,vartstamax,start,ncount,datatstamax);
	 ERRHAND(status);
	   for (l2=0;l2<NBL;l2++)
	     for (c=0;c<NBCTS;c++)
	       {
		 for (v=0;v<NPFT;v++)
		   {
		     hist[v]=nbvalv[v]=0;
		     datavego[v][l2][c]=fillval;
		     cveg=c*FACVEG;
		     lveg=l2*FACVEG;
		     for (l3=0;l3<2;l3++) 
		       for (c2=0;c2<2;c2++)
			 if ((lveg+l3<NBL*FACVEG)&&(cveg+c2<NBCVEG))
			   if ((dataveg[v*nblveg*NBCVEG+(lveg+l3)*NBCVEG+cveg+c2]<16)) 
			     {nbvalv[v]++;hist[v]+=dataveg[v*nblveg*NBCVEG+(lveg+l3)*NBCVEG+cveg+c2];}
		     if (nbvalv[v]>0) hist[v]/=nbvalv[v];
		   }
		 for (v2=0;v2<NVEG;v2++) 
		   {
		     switch (v2)
		       {
		       case 0:
			 v=(hist[14]>hist[13])? 14:13;
			 break;
		       case 1:
			 if ((hist[9]>hist[10])&&(hist[9]>hist[11])&&(hist[9]>hist[12])) v=9;
			 else 
			   {
			     if ((hist[11]>hist[10])&&(hist[11]>hist[12])) v=11;
			     else 
			       if (hist[12]>hist[10]) v=12;
				   else v=10;
			   }
			 break;
		       case 2:
			 v=(hist[4]>hist[1])? 4:1;
			 break;
		       case 3:
			 if ((hist[2]>hist[5])&&(hist[2]>hist[7])) v=2;
			 else 
			   {
			     if (hist[5]>hist[7]) v=5;
			     else v=7;
			   }
			 break;
		       case 4:
			 v=(hist[6]>hist[3])? 6:3;
			 break;
		       case 5:
			 v=8;
			 break;
		       case 6:
			 v=0;
			 break;
		       }		       		       
		     datavego[v][l2][c]=hist[v];
		     tsm=datatsmax[year-PYEAR][l2][c];
		     tsmm=meanmaxts[l2][c];
		     if ((hist[v]>0.1)&&(danandviav[v2][l2][c]>-1)&&(danandviav[v2][l2][c]<=1)&&(danandviap[v2][l2][c]>-1)&&(danandviap[v2][l2][c]<=1))
		       {
			 // if ((v==0)) printf ("(%e %e %e)",danandviav[v][l][c],dndvim[l][c],danandviap[v][l][c]);		
			 diffvi=(danandviap[v2][l2][c]-danandviav[v2][l2][c]);
			 if ((diffvi<-0.01) && (danats[l2][c]>4)&&(dndvim[v2][l2][c]>0.45)&&(dndviap[v2][l2][c]>0.4))
			   //if ((danats[l][c]>4))
			   {
			     //	     if (v>7) printf ("ok tseuil\n");
			     if(datatamean[l2][c]>25) tsmseuil[v][l2][c]+=datatsmax[year-PYEAR][l2][c];
			     if (year1[v][l2][c]==0) year1[v][l2][c]=year; else  year2[v][l2][c]=year;
			     if(datatamean[l2][c]>25) nbvalseuil[v][l2][c]++;
			      if ((datatamean[l2][c]>=5)&&(datatamean[l2][c]<35))
				{
				  histseuil[(int)((tsm-20.)*4.)][(int)((datatamean[l2][c]-5)*10)]++;
				  xvseuil[nbvalxy]=datatamean[l2][c];
				  yvseuil[nbvalxy]=tsm;
				  nbvalxy++;
				}
			   }
		
			 if ((tsm>40))
			   {
			     diffout[v][l2][c]+=diffvi;
			     anaav[v][l2][c]+=danandviav[v2][l2][c];
			     nbvaldiff[v][l2][c]++;
			   }
			 if ((tsm>20)&&(tsm<65)&&(tsmm>20)&&(tsmm<65)&&(dndvim[v2][l2][c]>0.45)&&( meanmaxts[l2][c]<45)&&(dndviap[v2][l2][c]>0.4))
			   {
			     /*  hist2D[v][(int)((tsmm-20.)*4.)][(int)((tsm-20.)*4.)]+=diffvi;
				 nval2D[v][(int)((tsmm-20.)*4.)][(int)((tsm-20.)*4.)]++;*/
			   }
			 if ((tsm>20)&&(tsm<65)&&(dndvim[v2][l2][c]>0.4)&&( meanmaxts[l2][c]<45)&&(dndviap[v2][l2][c]>0.3))
			   {
			     histom[v][(int)((tsm-20.)*4.)]+=diffvi;
			     histom2[v][(int)((tsm-20.)*4.)]+=danandvi[v2][l2][c];
			     histom3[v][(int)((tsm-20.)*4.)]+=datatstamax[l2][c];
			   nval[v][(int)((tsm-20.)*4.)]++;
			   varm[v][(int)((tsm-20.)*4.)]+=diffvi*diffvi;
			 }
		     }
		   if ((hist[v]>0.1)&&(danandviav[v2][l2][c]>-1)&&(danandviav[v2][l2][c]<=1)&&(danandviap[v2][l2][c]>-1)&&(danandviap[v2][l2][c]<=1))
		   {
		     ats=danats[l2][c];
		     if ((ats>-15)&&(ats<15)&&(dndvim[v2][l2][c]>0.3))
		       {
			 histts[v][(int)((ats+15.)*4.)]+=diffvi;
			 nvalts[v][(int)((ats+15.)*4.)]++;
		       }
		   }
		 }
	       }
	 ncount[0]=1;
	 ncount[1]=NPFT;
	 ncount[2]=NBL;
	 ncount[3]=NBCTS; 
	 start[0]=year-PYEAR;
	 start[2]=l;
	 start[1]=start[3]=0;
	 printf ("put1 %d %d %d %d %d %d %d %d\n",start[0],ncount[0],start[1],ncount[1],start[2],ncount[2],start[3],ncount[3]);
	 status=ncvarput(ncido,varanondvi,start,ncount,danandvi);
	 printf ("apres put\n");
      //	 status=ncvarput(ncido,varanondvi,start,ncount,dndvim);
       }
     for (l2=0;l2<NBL;l2++)
       for (c=0;c<NBCTS;c++)
	 for (v=0;v<NPFT;v++)
	   {
	     if (nbvalseuil[v][l2][c]>0)  
	       {
		 tsmseuil[v][l2][c]/=nbvalseuil[v][l2][c]; 
		 if ((datatamean[l2][c]>0)&&(tsmseuil[v][l2][c]>30)) 
		   {
		     fprintf (fouttt[v],"%f %f %f\n",datatamean[l2][c],meanmaxts[l2][c], tsmseuil[v][l2][c]);
		     if ((datatamean[l2][c]>5)&&(datatamean[l2][c]<35))
		       {
			 histtamseuil[0][(int)((datatamean[l2][c]-5)*10)]+=tsmseuil[v][l2][c]; 
			 nbvaltamseuil[0][(int)((datatamean[l2][c]-5)*10)]++;
			 histtaeseuil[0][(int)((datatamean[l2][c]-5)*10)]+=(tsmseuil[v][l2][c]*tsmseuil[v][l2][c]);
		       } 
		       }
	       }
	        else tsmseuil[v][l2][c]=nbvalseuil[v][l2][c]=year1[v][l2][c]=year2[v][l2][c]=fillval;
	     if (nbvaldiff[v][l2][c]>0)  {diffout[v][l2][c]/=nbvaldiff[v][l2][c]; anaav[v][l2][c]/=nbvaldiff[v][l2][c];} else diffout[v][l2][c]=anaav[v][l2][c]=fillval;
	     if ((tsmseuil[v][l2][c]>30)&&(tsmseuil[v][l2][c]<60))
	       {
		 meantseuil[v]+= tsmseuil[v][l2][c];
		 etseuil[v]+= (tsmseuil[v][l2][c]*tsmseuil[v][l2][c]);
		 nbvaltseuil[v]++;
	       }
	   }
     ncount[0]=NPFT;
     ncount[1]=NBL;
     ncount[2]=NBCTS; 
     start[1]=l;start[0]=start[2]=0;
     printf ("put1\n");
     status=ncvarput(ncido,varmtmax,start+1,ncount+1,meanmaxts);
     printf ("put2\n");
     status=ncvarput(ncido,vartseuil,start,ncount,tsmseuil);
     status=ncvarput(ncido,varvego,start,ncount,datavego);
     printf ("put3\n");
     status=ncvarput(ncido,varnbyear,start,ncount,nbvalseuil);
     printf ("put4\n");
     status=ncvarput(ncido,varyear1,start,ncount,year1);
     printf ("put5\n");
     status=ncvarput(ncido,varyear2,start,ncount,year2);
     status=ncvarput(ncido,vardiff,start,ncount,diffout);
    status=ncvarput(ncido,varanaav,start,ncount,anaav);
   }
 for(v=0;v<NPFT;v++)
   for (l=0;l<NBVAL2;l++)
     for (c=0;c<NBVAL2;c++)
       if (nval2D[v][l][c]>0) hist2D[v][l][c]/=nval2D[v][l][c];
 for (t=0;t<DYEAR-PYEAR+1;t++) time[t]=t*365+182;
 start[0]=0;
 ncount[0]=DYEAR-PYEAR+1;
 ncvarput(ncido,ntimep,start,ncount,time);
 start[0]=start[1]=start[2]=0;
 ncount[0]=NBVAL4;
 ncount[1]=NBVAL5;
 printf ("put6\n");
 ncvarput(ncido,nh2D,start,ncount,histseuil);
 for (v=0;v<NPFT;v++)
   {
     for (i=nbptot=0;i<NBVAL3;i++) nbptot+=nbvaltamseuil[v][i];
     for (i=0;i<NBVAL3;i++)
	 if (nbvaltamseuil[v][i])
	   {
	     histtamseuil[v][i]/=nbvaltamseuil[v][i];
	     histtaeseuil[v][i]=sqrt(histtaeseuil[v][i]/nbvaltamseuil[v][i]-histtamseuil[v][i]*histtamseuil[v][i]);
	   }
	 else histtamseuil[v][i]=-999;
      for (i=0;i<NBVAL3;i++)
	if (histtamseuil[v][i]>-999) fprintf(fouttt2[v],"%f %f %f %f\n",(float)i/10.+5.,histtamseuil[v][i],histtamseuil[v][i]-histtaeseuil[v][i],histtamseuil[v][i]+histtaeseuil[v][i]);
     for (i=nbptot=0;i<NBVAL2;i++) nbptot+=nval[v][i];
     for (i=0;i<NBVAL2;i++) 
       if ((nbptot>0)&&(nval[v][i]/nbptot>0.0005))
	 {
	   histom[v][i]/=nval[v][i];
	   histom2[v][i]/=nval[v][i];
	   histom3[v][i]/=nval[v][i];
	   varm[v][i]=1.862*sqrt(varm[v][i]/nval[v][i]-histom[v][i]*histom[v][i])/sqrt(nval[v][i]);
	   fprintf (fout[v],"%f %f %f %f %f %f %f\n",(float)i/4.+20,histom[v][i],histom2[v][i],histom3[v][i],histom[v][i]-varm[v][i],histom[v][i]+varm[v][i],nval[v][i]);
	 }
       else 
	 fprintf (fout[v],"%f 0 0 0 0 0 0\n",(float)i/4.+20);  
     if ( nbvaltseuil[v]>0) 
       {
	 meantseuil[v]/= nbvaltseuil[v];
	 fprintf (foutseuil,"%f %f %d\n",meantseuil[v],1.862*sqrt(etseuil[v]/nbvaltseuil[v]-meantseuil[v]*meantseuil[v]/sqrt(nval[v][i])),nbvaltseuil[v]);
       }
     else fprintf (foutseuil,"0 0 0\n");
   }
     for (v=0;v<NPFT;v++)
       for (i=0;i<NBVALTS;i++) 
	 if (nvalts[v][i]>100)
	   {
	     histts[v][i]/=nvalts[v][i];
	     fprintf (foutts[v],"%f %d\n",histts[v][i],nvalts[v][i]);
	   }
	 else 
	   fprintf (foutts[v],"0 0\n");  
     printf ("nbval seuil %d\n",nbvalxy);
          fit(xvseuil,yvseuil,nbvalxy,NULL,0,&a,&b,&siga,&sigb,&chi2,&q);
     printf ("a=%f b=%f,siga=%f sigb=%f chi2=%f,q=%f",a,b,siga,sigb,chi2,q);
     pearsn(xvseuil,yvseuil,nbvalxy,&r,&prob,&z);
     printf ("r=%f prob=%f z=%f\n",r,prob,z);
     
     
 status=ncclose(ncido);
}	
