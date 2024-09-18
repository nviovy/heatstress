#include <stdio.h>
#include <stdlib.h>
#include<netcdf.h>
#include <dirent.h> 
#include <unistd.h>
#include <fnmatch.h>
#include<string.h>
#include<math.h>
#define PYEAR 2001
#define DYEAR 2022
#define NBLTS 360
#define NBCTS 720
#define NBLVEG 7200
#define NBCVEG 14400
#define NBVAL 160
#define NBVAL2 180
#define NBVAL3 300
#define NBVALTS 120
#define NVEG 13
#define NVEG2 7
#define NPFT 13
#define NBYEAR 24
#define NBL  128
#define FACVEG 1.785714328
#define DTS 22.399999466
#define DIVWC 120
#define NBPOINTTERRE 67420
#define ERRHAND(e) if ((e) !=NC_NOERR) printf ("%s\n",nc_strerror((e)))
float danandvi[NVEG][NBPOINTTERRE], yearmin[NVEG][NBPOINTTERRE], danandvii[NVEG][NBPOINTTERRE], dyearmin[NVEG][NBPOINTTERRE], danandvio[NVEG][NBLTS][NBCTS],danaprec[NBPOINTTERRE],danaprecap[NBPOINTTERRE],/*danaevapav[NBPOINTTERRE],danaevapap[NBPOINTTERRE],danamrsoav[NBPOINTTERRE],danamrsoap[NBPOINTTERRE],*/danats[NBPOINTTERRE],danandviav[NVEG][NBPOINTTERRE],danandviap[NVEG][NBPOINTTERRE],dndvim[NVEG][NBPOINTTERRE],dndvimin[NVEG][NBPOINTTERRE],dndvimax[NVEG][NBPOINTTERRE],dndviap[NVEG][NBPOINTTERRE],datatsmax[DYEAR-PYEAR+1][NBPOINTTERRE],datatstamax[NBPOINTTERRE],datatstamax[NBPOINTTERRE],meanmaxts[NBLTS][NBCTS],nbval1[NBLTS][NBCTS],tsmseuil[NPFT][NBLTS][NBCTS],datavego[NPFT][NBLTS][NBCTS],nbvalseuil[NPFT][NBLTS][NBCTS],year1[NPFT][NBLTS][NBCTS],year2[NPFT][NBLTS][NBCTS],diffout[NPFT][NBLTS][NBCTS],anaav[NPFT][NBLTS][NBCTS],anaap[NPFT][NBLTS][NBCTS],ana[NPFT][NBLTS][NBCTS],nbvaldiff[NPFT][NBLTS][NBCTS],datatamean[NBPOINTTERRE];
float nvalmm[NBLTS][NBCTS],nbptot,nbptotv[NVEG2],prec[NBLTS][NBCTS],precap[NBLTS][NBCTS],evapav[NBLTS][NBCTS],evapap[NBLTS][NBCTS]/*,mrsoav[NBLTS][NBCTS],mrsoap[NBLTS][NBCTS]*/;
float hist2D[NPFT][NBVAL2][NBVAL2];
int nval2D[NPFT][NBVAL2][NBVAL2];
float xval[NBVAL2],yval[NBVAL2];
float histtamseuil[NPFT][NBVAL3], histtaeseuil[NPFT][NBVAL3],nbvaltamseuil[NPFT][NBVAL3];
int indiceligne[NBPOINTTERRE],indicecol[NBPOINTTERRE];
  float dataveg[NPFT][NBPOINTTERRE];
void main(int argc,char *argv[])
{
  FILE *fout[NPFT],*foutts[NPFT],*foutseuil,*fouttt[NPFT],*fouttt2[NPFT],*foutaf;
  int ncidparam[NBYEAR],ncidts,ncidveg,status,nxip,nyip,ntip,nlonp,nlap,ntimep,varanondvi,ncido,nvarlat,nvarlon,nlatp,ndim[4],vardiff,varveget,varvego,vartamean,ncidtamean,ncidmean,nvaridligne,nvaridcol;
  int varanoprec[NBYEAR],varanoprecap[NBYEAR],varanoevapav[NBYEAR],varanoevapap[NBYEAR],varanomrsoav[NBYEAR],varanomrsoap[NBYEAR],varanandvi[NBYEAR],varanandvii[NBYEAR],varanandviav[NBYEAR],varyearmin[NBYEAR],varanandviap[NBYEAR],vartsmax,vartstamax,varndvim[NBYEAR],varndvimin[NBYEAR],varndvimax[NBYEAR],varanats[NBYEAR],var,nvarlonts,nvarlatts,varmtmax,nxii,nyii,nh2D,nxv,nyv,vartseuil,varnbyear,varyear1,varyear2,nveg,varndviap[NBYEAR],varanaav,varanaap,varana,nvarmask,varpreco,varprecapo,varevapavo,varevapapo,varmrsoavo,varmrsoapo;
  float diffvi,tsm,histom[NPFT][NBVAL2],histts[NPFT][NBVALTS],varm[NPFT][NBVAL2],histom2[NPFT][NBVAL2],histom3[NPFT][NBVAL2];
  float hist1[NPFT][NBVAL2],hist2[NPFT][NBVAL2],hist3[NPFT][NBVAL2],varm1[NPFT][NBVAL2];
  int nvalts[NPFT][NBVALTS],year,i,l,c,t,v,l2,nblveg,nbvalv[NPFT],l3,c2,lveg,cveg,v2;
  float nval[NPFT][NBVAL2],ampli;
  long ncount[4],start[4];
  char nfic[100],str70[70];
  float fillval=-1e32,ats,fillval2=-1e20,tsmm,lat,hist[NPFT];
  double meantseuil[NPFT],etseuil[NPFT],nbvaltseuil[NPFT];
  float datalon[NBCTS],datalat[NBLTS],time[20];
  char ficout[20];
  int corre[NPFT]={7,2,3,4,2,3,4,3,5,1,1,0,0},i2;
  float nbpixaf[NBYEAR];
  memset (hist2D,0,NBVAL2*NBVAL2*NPFT*sizeof(float));
  memset (nval2D,0,NBVAL2*NBVAL2*NPFT*sizeof(int));
  memset (hist1,0,NBVAL2*NPFT*sizeof(float));
  memset (hist2,0,NBVAL2*NPFT*sizeof(float));
  memset (hist3,0,NBVAL2*NPFT*sizeof(float));
  memset (varm1,0,NBVAL2*NPFT*sizeof(float));
   for (i=0;i<NBVAL2;i++) xval[i]=yval[i]=i/4.+20.;
  for (v=0;v<NPFT;v++) meantseuil[NPFT]=etseuil[NPFT]=nbvaltseuil[v]=0;
  /* ncidts=ncopen("data/daytmaxyglob.nc",NC_NOWRITE);
  status=nc_inq_varid(ncidts,"tsta_max",&vartstamax);
  status=nc_inq_varid(ncidts,"tsmax",&vartsmax);
  ncidveg=ncopen("data/PFTmap_2015.nc",NC_NOWRITE);
  status=nc_inq_varid(ncidveg,"maxvegetfrac",&varveget);
  ncidtamean=ncopen("data/meantemp.nc",NC_NOWRITE);
  status=nc_inq_varid(ncidtamean,"tmean",&vartamean);*/

  foutaf=fopen("fpixaf","w");

for (v=0;v<NPFT;v++)
   {
     sprintf (ficout,"histtndvighs_%d",v+1);
     fout[v]=fopen(ficout,"w");
     sprintf (ficout,"histtstndvighs_%d",v+1);
     foutts[v]=fopen(ficout,"w");
     sprintf (ficout,"teuilvsths_%d",v+1);
     fouttt[v]=fopen(ficout,"w");
     sprintf (ficout,"tseuilmths_%d",v+1);
     fouttt2[v]=fopen(ficout,"w");

     for (i=0;i<NBVAL2;i++) histom[v][i]=histom2[v][i]=histom3[v][i]=varm[v][i]=nval[v][i]=0;
     for (i=0;i<NBVALTS;i++) histts[v][i]=nvalts[v][i]=0;
   } 
 foutseuil=fopen("histoseuilhs","w");
 ncido=nccreate("output/repndvipftglobhs.nc",NC_64BIT_OFFSET);
 nxip = ncdimdef(ncido,"lon",NBCTS);
 nyip = ncdimdef(ncido,"lat",NBLTS);
 ntip = ncdimdef(ncido,"time",NC_UNLIMITED);
 nveg = ncdimdef(ncido,"veg",NPFT);
 nlonp = ncvardef(ncido,"lon",NC_FLOAT,1,&nxip);
 nlatp = ncvardef(ncido,"lat",NC_FLOAT,1,&nyip);
 ntimep = ncvardef(ncido,"time",NC_DOUBLE,1,&ntip);
 nxii= ncdimdef(ncido,"TSmax",NBVAL2);
 nyii= ncdimdef(ncido,"TSmax_mean",NBVAL2);
nxv=ncvardef(ncido,"TSmax",NC_FLOAT,1,&nxii);
nyv=ncvardef(ncido,"TSmax_mean",NC_FLOAT,1,&nyii);
 ndim[0]=nveg;
 ndim[1]=nyii;
 ndim[2]=nxii;
 nh2D=ncvardef(ncido,"histo",NC_FLOAT,3,ndim);
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
 varana = ncvardef(ncido,"ana",NC_FLOAT,3,ndim+1);
 varanaav = ncvardef(ncido,"anav",NC_FLOAT,3,ndim+1);
 varanaap = ncvardef(ncido,"anap",NC_FLOAT,3,ndim+1);
 varyear1 = ncvardef(ncido,"year1",NC_FLOAT,3,ndim+1);
 varyear2 = ncvardef(ncido,"year2",NC_FLOAT,3,ndim+1);
 varpreco = ncvardef(ncido,"dprec",NC_FLOAT,2,ndim+2);
 varprecapo = ncvardef(ncido,"dprecap",NC_FLOAT,2,ndim+2);
 /* varevapavo = ncvardef(ncido,"evapav",NC_FLOAT,2,ndim+2);
 varevapapo = ncvardef(ncido,"evapap",NC_FLOAT,2,ndim+2);
 varmrsoavo = ncvardef(ncido,"mrsoav",NC_FLOAT,2,ndim+2);
 varmrsoapo = ncvardef(ncido,"mrsoap",NC_FLOAT,2,ndim+2);*/

 ncattput(ncido, varanondvi, "_FillValue",NC_FLOAT,1,&fillval);
 ncattput(ncido, varmtmax, "_FillValue",NC_FLOAT,1,&fillval);
 ncattput(ncido, vartseuil, "_FillValue",NC_FLOAT,1,&fillval);
 ncattput(ncido, varnbyear, "_FillValue",NC_FLOAT,1,&fillval);
 ncattput(ncido, varyear1, "_FillValue",NC_FLOAT,1,&fillval);
 ncattput(ncido, vardiff, "_FillValue",NC_FLOAT,1,&fillval);
 ncattput(ncido, varanaav, "_FillValue",NC_FLOAT,1,&fillval);
 ncattput(ncido, varana, "_FillValue",NC_FLOAT,1,&fillval);
 ncattput(ncido, varanaap, "_FillValue",NC_FLOAT,1,&fillval);
  ncattput(ncido, varyear2, "_FillValue",NC_FLOAT,1,&fillval);
  ncattput(ncido, varpreco, "_FillValue",NC_FLOAT,1,&fillval);
  ncattput(ncido, varprecapo, "_FillValue",NC_FLOAT,1,&fillval);
  ncattput(ncido, varvego, "_FillValue",NC_FLOAT,1,&fillval);
  /* ncattput(ncido, varevapavo, "_FillValue",NC_FLOAT,1,&fillval);
  ncattput(ncido, varevapapo, "_FillValue",NC_FLOAT,1,&fillval);
  ncattput(ncido, varmrsoavo, "_FillValue",NC_FLOAT,1,&fillval);
  ncattput(ncido, varmrsoapo, "_FillValue",NC_FLOAT,1,&fillval);*/
 sprintf (str70,"days since %04d-%02d-%02d %02d:%02d:%02d",2012,1,1,0,0,0);
 ncattput(ncido,ntimep, "units",NC_CHAR, strlen(str70), str70);
 ncattput(ncido,ntimep, "axis",NC_CHAR, 1, "T");
 status=ncattput (ncido,nlonp,"units",NC_CHAR,12,"degrees_east");
 status=ncattput (ncido,nlonp,"axis",NC_CHAR,1,"X");
 status=ncattput (ncido,nlatp,"units",NC_CHAR,13,"degrees_north");
 status=ncattput (ncido,nlatp,"axis",NC_CHAR,1,"Y");
 ncendef(ncido);
 ncidts=ncopen("output/analogpftorchths_2018.nc",NC_NOWRITE);
 nc_inq_varid(ncidts,"longitude",&nvarlonts);
 nc_inq_varid(ncidts,"latitude",&nvarlatts);
 nc_inq_varid(ncidts,"ref_point",&nvarmask);
  nc_inq_varid(ncidts,"indice_ligne",&nvaridligne);
  nc_inq_varid(ncidts,"indice_col",&nvaridcol);
 nc_inq_varid(ncidts,"tsmax",&vartsmax);
 start[0]=0;
 ncount[0]=NBLTS;
 ncvarget(ncidts,nvarlatts,start,ncount,datalat);
ncvarput(ncido,nlatp,start,ncount,datalat);	  
 start[0]=0;
 ncount[0]=NBCTS;
 ncvarget(ncidts,nvarlonts,start,ncount,datalon);
 ncvarput(ncido,nlonp,start,ncount,&datalon);
 start[0]=0;
 ncount[0]=NBPOINTTERRE;
 ncvarget(ncidts,nvaridligne,start,ncount,indiceligne);
ncvarget(ncidts,nvaridcol,start,ncount,indicecol);
 start[0]=0;
 start[1]=start[2]=0;
 ncount[0]=DYEAR-PYEAR+1;
 ncount[1]=NBPOINTTERRE;
	 printf ("varget\n");
 status=ncvarget(ncidts,vartsmax,start,ncount,datatsmax);
	 printf ("varget\n");
 for (l=0;l<NBLTS;l++) for(c=0;c<NBCTS;c++)  meanmaxts[l][c]=nvalmm[l][c]=0;
 for (year=PYEAR;year<=DYEAR;year++)
     for (i=0;i<NBPOINTTERRE;i++)
	 {
	   l=indiceligne[i];
	   c=indicecol[i];
	   meanmaxts[l][c]+=datatsmax[year-PYEAR][i];
	   nvalmm[l][c]++;
	 }
  for (l=0;l<NBLTS;l++) for(c=0;c<NBCTS;c++)   
     if (nvalmm[l][c]>0) meanmaxts[l][c]/= nvalmm[l][c]; else  meanmaxts[l][c]=fillval;


 for (year=PYEAR;year<=DYEAR;year++)
   {
     sprintf (nfic,"output/analogpftorchths_%d.nc",year);
     ncidparam[year-PYEAR]=ncopen(nfic,NC_NOWRITE);
     status=nc_inq_varid(ncidparam[year-PYEAR],"anoprecav",&varanoprec[year-PYEAR]);
     status=nc_inq_varid(ncidparam[year-PYEAR],"anoprecap",&varanoprecap[year-PYEAR]);
     status=nc_inq_varid(ncidparam[year-PYEAR],"anats",&varanats[year-PYEAR]);
     status=nc_inq_varid(ncidparam[year-PYEAR],"anandvi",&varanandvi[year-PYEAR]);
     status=nc_inq_varid(ncidparam[year-PYEAR],"anandvii",&varanandvii[year-PYEAR]);
     status=nc_inq_varid(ncidparam[year-PYEAR],"anandviav",&varanandviav[year-PYEAR]);
     status=nc_inq_varid(ncidparam[year-PYEAR],"anandviap",&varanandviap[year-PYEAR]);
     status=nc_inq_varid(ncidparam[year-PYEAR],"ndvim",&varndvim[year-PYEAR]);
      status=nc_inq_varid(ncidparam[year-PYEAR],"ndvimin",&varndvimin[year-PYEAR]);
     status=nc_inq_varid(ncidparam[year-PYEAR],"ndvimax",&varndvimax[year-PYEAR]);
   status=nc_inq_varid(ncidparam[year-PYEAR],"ndviap",&varndviap[year-PYEAR]);
    status=nc_inq_varid(ncidparam[year-PYEAR],"yearmin",&varyearmin[year-PYEAR]);
    /*    status=nc_inq_varid(ncidparam[year-PYEAR],"anoevapav",&varanoevapav[year-PYEAR]);
    status=nc_inq_varid(ncidparam[year-PYEAR],"anoevapap",&varanoevapap[year-PYEAR]);
    status=nc_inq_varid(ncidparam[year-PYEAR],"anomrsoav",&varanomrsoav[year-PYEAR]);
    status=nc_inq_varid(ncidparam[year-PYEAR],"anomrsoap",&varanomrsoap[year-PYEAR]);*/
    if (year==PYEAR) status=nc_inq_varid(ncidparam[year-PYEAR],"veg",&varveget);
    nbpixaf[year-PYEAR]=0;
   }
 for (l=0;l<NBLTS;l+=NBL)
   {
     printf ("l=%d\n",l);
     memset (tsmseuil,0,NBLTS*NBCTS*NPFT*sizeof(int));
     memset (nbvalseuil,0,NBLTS*NBCTS*NPFT*sizeof(int));
     memset (nbvaldiff,0,NBLTS*NBCTS*NPFT*sizeof(int));
     memset (year1,0,NBLTS*NBCTS*NPFT*sizeof(int));
     memset (year2,0,NBLTS*NBCTS*NPFT*sizeof(int));
    memset (prec,0,NBLTS*NBCTS*sizeof(float));
     memset (precap,0,NBLTS*NBCTS*sizeof(float));
     /*    memset (evapav,0,NBLTS*NBCTS*sizeof(float));
    memset (evapap,0,NBLTS*NBCTS*sizeof(float));
    memset (mrsoav,0,NBLTS*NBCTS*sizeof(float));
    memset (mrsoap,0,NBLTS*NBCTS*sizeof(float));*/
    memset (nbval1,0,NBLTS*NBCTS*sizeof(float));

     for (v=0;v<NPFT;v++) for(l2=0;l2<NBLTS;l2++) for (c=0;c<NBCTS;c++) {diffout[v][l2][c]=anaav[v][l2][c]=anaap[v][l2][c]=ana[v][l2][c]=0;}
     start[0]=0;
     start[1]=l;
       start[2]=0;
     ncount[0]=DYEAR-PYEAR+1;
     ncount[1]=NBL;
     ncount[2]=NBCTS;
     status=ncvarget(ncidts,vartsmax,start,ncount,datatsmax); 
      printf ("get tsmax\n");
    start[0]=0;
     start[1]=0;
     ncount[0]=NPFT;
     ncount[1]=NBPOINTTERRE;
     status=ncvarget(ncidparam[0],varveget,start,ncount,dataveg); 
     printf ("get dataveg\n");
      for (year=PYEAR;year<=DYEAR;year++)
       {
	 start[0]=0;
	 ncount[0]=NBPOINTTERRE;
	 status=ncvarget(ncidparam[year-PYEAR],varanats[year-PYEAR],start,ncount,danats);
	 status=ncvarget(ncidparam[year-PYEAR],varanoprec[year-PYEAR],start,ncount,danaprec);
	 status=ncvarget(ncidparam[year-PYEAR],varanoprecap[year-PYEAR],start,ncount,danaprecap);
     printf ("get dataprec\n");
	 /*	 status=ncvarget(ncidparam[year-PYEAR],varanoevapav[year-PYEAR],start,ncount,danaevapav);
	 status=ncvarget(ncidparam[year-PYEAR],varanoevapap[year-PYEAR],start,ncount,danaevapap);
	 status=ncvarget(ncidparam[year-PYEAR],varanomrsoav[year-PYEAR],start,ncount,danamrsoav);
	 status=ncvarget(ncidparam[year-PYEAR],varanomrsoap[year-PYEAR],start,ncount,danamrsoap);*/
	 ERRHAND(status);
	 start[0]=0;
	 ncount[0]=NVEG;
	 start[1]=0;
	 ncount[1]=NBPOINTTERRE;
	 status=ncvarget(ncidparam[year-PYEAR],varanandvi[year-PYEAR],start,ncount,danandvi);
	 ERRHAND(status);
     printf ("get dataana\n");
     //	 status=ncvarget(ncidparam[year-PYEAR],varanandvii[year-PYEAR],start,ncount,danandvii);
	 ERRHAND(status);
     printf ("get dataana\n");
	 status=ncvarget(ncidparam[year-PYEAR],varanandviav[year-PYEAR],start,ncount,danandviav);
     printf ("get dataana\n");
	 ERRHAND(status);
	 status=ncvarget(ncidparam[year-PYEAR],varanandviap[year-PYEAR],start,ncount,danandviap);
     printf ("get dataana\n");
	 ERRHAND(status);
	 status=ncvarget(ncidparam[year-PYEAR],varndvim[year-PYEAR],start,ncount,dndvim);
     printf ("get dataana\n");
	 status=ncvarget(ncidparam[year-PYEAR],varndvimin[year-PYEAR],start,ncount,dndvimin);
     printf ("get dataana\n");
	 status=ncvarget(ncidparam[year-PYEAR],varndvimax[year-PYEAR],start,ncount,dndvimax);
     printf ("get dataana\n");
	 ERRHAND(status);
	 status=ncvarget(ncidparam[year-PYEAR],varndviap[year-PYEAR],start,ncount,dndviap);
     printf ("get dataana\n");
	 ERRHAND(status);
	 status=ncvarget(ncidparam[year-PYEAR],varyearmin[year-PYEAR],start,ncount,dyearmin);
	 ERRHAND(status);
     printf ("get dataana\n");
	   for (i=0;i<NBPOINTTERRE;i++)
	       {
		 l2=indiceligne[i];
		 c=indicecol[i];
		 for (v=0;v<NPFT;v++)
		   {
		     v2=corre[v];
		     hist[v]=nbvalv[v]=danandvio[v][l2][c]=0;
		     datavego[v][l2][c]=dataveg[v][i];
		     //   danandvio[v][l2][c]=danandvi[v][i];
		     if ((dataveg[v][i]>0.1)&&(danandviav[v][i]>-1)&&(danandviav[v][i]<=1)&&(danandviap[v][i]>-1)&&(danandviap[v][i]<=1))
		       {
			 // if ((v==0)) printf ("(%e %e %e)",danandviav[v][l][c],dndvim[l][c],danandviap[v][l][c]);		
			 ampli=dndvimax[v][i]-dndvimin[v][i];
			   if (ampli>0.2) diffvi=(danandviap[v][i]-danandviav[v][i])/ampli; else diffvi=(danandviap[v][i]-danandviav[v][i]);
			   //diffvi=(danandviap[v][i]-danandviav[v][i]);
			 // printf ("(%f %f %f)",danandviav[v][i],danandviap[v][i],danats[i]); 
			 if ((diffvi<-0.01) && (danats[i]>4))
			   //if ((danats[l][c]>4))
			   {
			     //	     if (v>7) printf ("ok tseuil\n");use tes
			     if (datatsmax[year-PYEAR][i]>100) printf ("%d %d %f\n",year,i,datatsmax[year-PYEAR][i]);
			     tsmseuil[v][l2][c]+=datatsmax[year-PYEAR][i];
			     nbvalseuil[v][l2][c]++;
			     if (tsm>35)  nbpixaf[year-PYEAR]+=3086.358025*cosf(fabs(90.-l2*0.5)*M_PI/180.);
			   }
			 tsm=datatsmax[year-PYEAR][i];
			 tsmm=meanmaxts[l2][c];
			 if ((tsm>50)&&(danats[i]>4)&&(fabs(danandvi[v][i])<0.05))
			   {
			     diffout[v][l2][c]+=diffvi;
			     if (v==10)
			       {
				  prec[l2][c]+=danaprec[i];
				 precap[l2][c]+=danaprecap[i];
				 /*				 evapav[l2][c]+=danaevapav[i];
				 evapap[l2][c]+=danaevapap[i];
				 mrsoav[l2][c]+=danamrsoav[i];
				 mrsoap[l2][c]+=danamrsoap[i];*/
				 nbval1[l2][c]++;
			       }
			     danandvio[v][l2][c]=diffvi;
			     anaav[v][l2][c]+=danandviav[v][i];
			     ana[v][l2][c]+=danandvi[v][i];
			     anaap[v][l2][c]+=danandviap[v][i];
			     nbvaldiff[v][l2][c]++;
			     year1[v][l2][c]=year; 
			     year2[v][l2][c]=dyearmin[v][i];
			     //year2[v][l2][c]= nbvaldiff[v][l2][c];

			   }
			 if ((tsm>20)&&(tsm<65)&&(tsmm>20)&&(tsmm<65)&&(dndvim[v][i]>0.45)&&(dndviap[v][i]>0.4)&&(danats[i]>4)&&(fabs(danandvii[v][i])<0.05))
			   {
			     hist2D[v][(int)((tsmm-20.)*4.)][(int)((tsm-20.)*4.)]+=diffvi;
			     nval2D[v][(int)((tsmm-20.)*4.)][(int)((tsm-20.)*4.)]++;

			   }
			 if ((tsm>20)&&(tsm<65)&&(dndvim[v][i]>0.4)&&( meanmaxts[l2][c]<45)&&(dndviap[v][i]>0.3)&&(danats[i]>4)&&(fabs(danandvii[v][i])<0.05))
			   {
			     histom[v2][(int)((tsm-20.)*4.)]+=diffvi;
			     //  if ((v2==1)&&((int)((tsm-20.)*4.)==40)) printf ("(%f %f: %f %f %f)\n",danandviav[v][i], danandviap[v][i],histom[v2][(int)((tsm-20.)*4.)],diffvi,nval[v2][(int)((tsm-20.)*4.)]+1);
			     histom2[v2][(int)((tsm-20.)*4.)]+=danandvi[v][i];

			     //    histom3[v][(int)((tsm-20.)*4.)]+=datatstamax[l2][c];
			   nval[v2][(int)((tsm-20.)*4.)]++;
			   varm[v2][(int)((tsm-20.)*4.)]+=diffvi*diffvi;
			 }
		     }
		   if ((dataveg[v][i]>0.1)&&(danandviav[v][i]>-1)&&(danandviav[v][i]<=1)&&(danandviap[v][i]>-1)&&(danandviap[v][i]<=1))
		   {
		     ats=danats[i];
		     if ((ats>-15)&&(ats<15)&&(dndvim[v][i]>0.3))
		       {
			 histts[v][(int)((ats+15.)*4.)]+=diffvi;
			 nvalts[v][(int)((ats+15.)*4.)]++;
		       }
		   }
		 }
	       }
	 ncount[0]=1;
	 ncount[1]=NPFT;
	 ncount[2]=NBLTS;
	 ncount[3]=NBCTS; 
	 start[0]=year-PYEAR;
	 start[2]=0;
	 start[1]=start[3]=0;
	 printf ("put1 %d %d %d %d %d %d %d %d\n",start[0],ncount[0],start[1],ncount[1],start[2],ncount[2],start[3],ncount[3]);
	 status=ncvarput(ncido,varanondvi,start,ncount,danandvio);
	 printf ("apres put\n");
      //	 status=ncvarput(ncido,varanondvi,start,ncount,dndvim);
       }

    for (l2=0;l2<NBLTS;l2++)
       for (c=0;c<NBCTS;c++)
	 {
	   for (v=0;v<NPFT;v++)
	     {
	       if (nbvalseuil[v][l2][c]>0)  
		 {
		   tsmseuil[v][l2][c]/=nbvalseuil[v][l2][c]; 
		 }
	       if (nbvaldiff[v][l2][c]>0)
		 {
		   diffout[v][l2][c]/=nbvaldiff[v][l2][c]; 
		   anaav[v][l2][c]/=nbvaldiff[v][l2][c];
		   ana[v][l2][c]/=nbvaldiff[v][l2][c];
		   anaap[v][l2][c]/=nbvaldiff[v][l2][c];
		  } 
	     }
	   if (nbval1[l2][c]>0)
	     {
	       prec[l2][c]/=nbval1[l2][c];
	       precap[l2][c]/=nbval1[l2][c];
	       /*	       evapav[l2][c]/=nbval1[l2][c];
	       evapap[l2][c]/=nbval1[l2][c];
	       mrsoav[l2][c]/=nbval1[l2][c];
	       mrsoap[l2][c]/=nbval1[l2][c];*/
	     }
	 }
    ncount[0]=NPFT;
     ncount[1]=NBLTS;
     ncount[2]=NBCTS; 
     start[1]=0;start[0]=start[2]=0;
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
    status=ncvarput(ncido,varana,start,ncount,ana);
   status=ncvarput(ncido,varanaav,start,ncount,anaav);
    status=ncvarput(ncido,varpreco,start,ncount+1,prec);
    printf ("putprec\n");
    status=ncvarput(ncido,varprecapo,start,ncount+1,precap);
    /*    status=ncvarput(ncido,varevapavo,start,ncount+1,evapav);
    status=ncvarput(ncido,varevapapo,start,ncount+1,evapap);
    status=ncvarput(ncido,varmrsoavo,start,ncount+1,mrsoav);
    status=ncvarput(ncido,varmrsoapo,start,ncount+1,mrsoap);*/
   printf ("putprecap\n");
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
 ncount[0]=NPFT;
 ncount[1]=NBVAL2;
 ncount[2]=NBVAL2;
 printf ("put6\n");
 ncvarput(ncido,nh2D,start,ncount,hist2D);
 printf ("put7\n");

 for (v=0;v<NPFT;v++)
   {
     v2=corre[v];
   for (i=nbptot=0;i<NBVAL3;i++) nbptot+=nbvaltamseuil[v][i];
     for (i=0;i<NBVAL3;i++)
	 if (nbvaltamseuil[v][i])
	   {
	     histtamseuil[v][i]/=nbvaltamseuil[v][i];
	     histtaeseuil[v][i]=sqrt(histtaeseuil[v][i]/nbvaltamseuil[v][i]-histtamseuil[v][i]*histtamseuil[v][i]);
	   }
	 else histtamseuil[v][i]=-999;
      for (i=0;i<NBVAL3;i++)
	fprintf(fouttt2[v],"%f %f %f\n",(float)i/10.+5.,histtamseuil[v][i],histtaeseuil[v][i]);
   }
 for (v2=0;v2<NVEG2;v2++)
   {
     for (i=nbptotv[v2]=0;i<NBVAL2;i++) nbptotv[v2]+=nval[v2][i];
     for (i=0;i<NBVAL2;i++) 
       if ((nbptotv[v2]>0)&&(nval[v2][i]/nbptot>0.0005))
	 {
	   printf ("%d %d: %f %f %f \n", v2,i,histom[v2][i], nval[v2][i],histom[v2][i]/nval[v2][i]);
	   	   histom[v2][i]/=nval[v2][i];
	   histom2[v2][i]/=nval[v2][i];
	   histom3[v2][i]/=nval[v2][i];
	   varm[v2][i]=1.862*sqrt(varm[v2][i]/nval[v2][i]-histom[v2][i]*histom[v2][i])/sqrt(nval[v2][i]);
	 }
     for (i=0;i<NBVAL2;i++) 	   
      if ((nbptotv[v2]>0)&&(nval[v2][i]/nbptot>0.0005))
	 {
	   if ((i>=2)&&(i<NBVAL2-2)) 
	     {
	       for (i2=-2;i2<=2;i2++)
		 {
		   hist1[v2][i]+=histom[v2][i+i2]/5;
		   hist2[v2][i]+=histom2[v2][i+i2]/5;
		   hist3[v2][i]+=histom3[v2][i+i2]/5;
		   varm1[v2][i]+=varm[v2][i+i2]/5;
		 }
	     }
	   else
	     {
	       hist1[v2][i]=histom[v2][i];
	       hist2[v2][i]=histom2[v2][i];
	       hist3[v2][i]=histom3[v2][i];
	       varm1[v2][i]=varm[v2][i];
	     }
	 }
   }
 for (v=0;v<NPFT;v++)
   if ( nbvaltseuil[v]>0) 
       {
	 meantseuil[v]/= nbvaltseuil[v];
	 fprintf (foutseuil,"%f %f %d\n",meantseuil[v],sqrt(etseuil[v]/nbvaltseuil[v]-meantseuil[v]*meantseuil[v]),nbvaltseuil[v]);
       }
     else fprintf (foutseuil,"0 0 0\n");
 for (v2=0;v2<NVEG2;v2++) 
     for (i=0;i<NBVAL2;i++) 
      if ((nbptotv[v2]>0)&&(nval[v2][i]/nbptot>0.0005))
		   fprintf (fout[v2],"%f %e %e %e %e %e %e\n",(float)i/4.+20,hist1[v2][i],hist2[v2][i],hist3[v2][i],hist1[v2][i]-varm1[v2][i],hist1[v2][i]+varm1[v2][i],nval[v2][i]);
       else 
	 	 fprintf (fout[v2],"%f 0 0 0 0 0 0\n",(float)i/4.+20);  
 /*
     if ( nbvaltseuil[v]>0) 
     for (v=0;v<NPFT;v++)
       for (i=0;i<NBVALTS;i++) 
	 if (nvalts[v][i]>100)
	   {
	     histts[v][i]/=nvalts[v][i];
	     fprintf (foutts[v],"%f %d\n",histts[v][i],nvalts[v][i]);
	   }
	 else 
	 fprintf (foutts[v],"0 0\n");  */
     start[0]=0;
 ncount[0]=NBVAL2;
  ncvarput(ncido,nxv,start,ncount,xval);
  ncvarput(ncido,nyv,start,ncount,yval);
 status=ncclose(ncido);
  for (year=PYEAR;year<=DYEAR;year++)  fprintf (foutaf,"%d %f\n",year,nbpixaf[year-PYEAR]);
 fclose (foutaf);

}	
