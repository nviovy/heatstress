#include <stdio.h>
#include <stdlib.h>
#include<netcdf.h>
#include <dirent.h> 
#include <unistd.h>
#include <fnmatch.h>
#include<string.h>
#include<math.h>
#define PYEAR 2011
#define DYEAR 2019
#define NBLTS 3584
#define NBCTS 8064
#define NBVAL 160
#define NBVAL2 180
#define NBVALTS 120
#define NVEG 7
#define NVEG2 7
#define NBVAL3 300
#define NBVAL4 300
#define NBYEAR 9
#define NBL  128
#define ERRHAND(e) if ((e) !=NC_NOERR) printf ("%s\n",nc_strerror((e)))
float danandvi[NVEG][NBL][NBCTS], danats[NBL][NBCTS],danandviav[NVEG][NBL][NBCTS],danandviap[NVEG][NBL][NBCTS],dndvim[NVEG][NBL][NBCTS],dndvimin[NVEG][NBL][NBCTS],dndvimax[NVEG][NBL][NBCTS],dndviap[NVEG][NBL][NBCTS],datatsmax[DYEAR-PYEAR+1][NBL][NBCTS],datatstamax[NBL][NBCTS],map[NVEG][NBL][NBCTS],datatstamax[NBL][NBCTS],meanmaxts[NBL][NBCTS],tsmseuil[NVEG][NBL][NBCTS],nbvalseuil[NVEG][NBL][NBCTS],year1[NVEG][NBL][NBCTS],year2[NVEG][NBL][NBCTS],diffout[NVEG][NBL][NBCTS],diffout2[NVEG][NBL][NBCTS],anaav[NVEG][NBL][NBCTS],nbvaldiff[NVEG][NBL][NBCTS];
float nvalmm[NBL][NBCTS],nbptot;
float hist2D[NVEG2][NBVAL2][NBVAL2];
int nval2D[NVEG2][NBVAL2][NBVAL2];
float xval[NBVAL2],yval[NBVAL2];
float histtamseuil[NVEG][NBVAL3], histtaeseuil[NVEG][NBVAL3],nbvaltamseuil[NVEG][NBVAL3], histtamseuilnv[NBVAL3], histtaeseuilnv[NBVAL3],nbvaltamseuilnv[NBVAL3];
void main(int argc,char *argv[])
{
  FILE *fout[NVEG2],*foutts[NVEG2];
  int ncidparam[NBYEAR],ncidts,status,nxip,nyip,ntip,nlonp,nlap,ntimep,varanondvi,ncido,nvarlat,nvarlon,nlatp,ndim[4],v2,vardiff;
  int varanandvi[NBYEAR],varanandviav[NBYEAR],varanandviap[NBYEAR],vartsmax,vartstamax,varndvim[NBYEAR],varndvimin[NBYEAR],varndvimax[NBYEAR],varanats[NBYEAR],var,nvarlonts,nvarlatts,varmtmax,nxii,nyii,nh2D,nxv,nyv,vartseuil,varnbyear,varyear1,varyear2,nveg,varndviap[NBYEAR],varanaav;
  float diffvi,tsm,histom[NVEG2][NBVAL2],histts[NVEG2][NBVALTS],varm[NVEG2][NBVAL2],histom2[NVEG2][NBVAL2],histom3[NVEG2][NBVAL2],ampli;
  int nvalts[NVEG2][NBVALTS],year,i,l,c,t,v,l2;
  float nval[NVEG2][NBVAL2];
  long ncount[4],start[4];
  char nfic[100],str70[70];
  float fillval=-1e32,ats,fillval2=-1e20,tsmm,lat,vdiff;
  double datalon[NBCTS],datalat[NBLTS],time[10],hist1,hist2,hist3,varm1;
  char ficout[20];
  int cr,lr,c2,l3,passe,pointtraite,rayon,i2;
  float valr,nbvalr;
  memset (hist2D,0,NBVAL2*NBVAL2*NVEG*sizeof(float));
  memset (nval2D,0,NBVAL2*NBVAL2*NVEG*sizeof(int));

 ncidts=ncopen("data/daytmaxyglob.nc",NC_NOWRITE);
  ERRHAND(status);
 status=nc_inq_varid(ncidts,"tsta_max",&vartstamax);
 	status=nc_inq_varid(ncidts,"tsmax",&vartsmax);
 ERRHAND(status);
for (v=0;v<NVEG2;v++)
   {
     sprintf (ficout,"histtlaimod_%d",v+1);
     fout[v]=fopen(ficout,"w");
     sprintf (ficout,"histtslaimod_%d",v+1);
     foutts[v]=fopen(ficout,"w");
     for (i=0;i<NBVAL2;i++) histom[v][i]=histom2[v][i]=histom3[v][i]=varm[v][i]=nval[v][i]=0;
     for (i=0;i<NBVALTS;i++) histts[v][i]=nvalts[v][i]=0;
   } 
 ncido=nccreate("data/repndvipftglobmodlai.nc",NC_64BIT_OFFSET);
 nxip = ncdimdef(ncido,"lon",NBCTS);
 nyip = ncdimdef(ncido,"lat",NBLTS);
 ntip = ncdimdef(ncido,"time",NC_UNLIMITED);
 nveg = ncdimdef(ncido,"veg",NVEG2);
 nlonp = ncvardef(ncido,"lon",NC_DOUBLE,1,&nxip);
 nlatp = ncvardef(ncido,"lat",NC_DOUBLE,1,&nyip);
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
 sprintf (str70,"days since %04d-%02d-%02d %02d:%02d:%02d",2011,1,1,0,0,0);
 ncattput(ncido,ntimep, "units",NC_CHAR, strlen(str70), str70);
 ncattput(ncido,ntimep, "axis",NC_CHAR, 1, "T");
 status=ncattput (ncido,nlonp,"units",NC_CHAR,12,"degrees_east");
 status=ncattput (ncido,nlonp,"axis",NC_CHAR,1,"X");
 status=ncattput (ncido,nlatp,"units",NC_CHAR,13,"degrees_north");
 status=ncattput (ncido,nlatp,"axis",NC_CHAR,1,"Y");
 ncendef(ncido);
  for (i=0;i<NBVAL2;i++) xval[i]=yval[i]=i/4.+20.;
  for (i=0;i<NBVAL2;i++) yval[i]=((float)i-80.)/100.;
  for (i=0;i<NBVAL2;i++) printf ("%f ",yval[i]);
  printf ("\n");
     start[0]=0;
 ncount[0]=NBVAL2;
  ncvarput(ncido,nxv,start,ncount,xval);
  ncvarput(ncido,nyv,start,ncount,yval);

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
     sprintf (nfic,"data/analogtpftglobmodlai_%d.nc",year);
     ncidparam[year-PYEAR]=ncopen(nfic,NC_NOWRITE);
     status=nc_inq_varid(ncidparam[year-PYEAR],"anats",&varanats[year-PYEAR]);
     status=nc_inq_varid(ncidparam[year-PYEAR],"anandvi",&varanandvi[year-PYEAR]);
     status=nc_inq_varid(ncidparam[year-PYEAR],"anandviav",&varanandviav[year-PYEAR]);
     status=nc_inq_varid(ncidparam[year-PYEAR],"anandviap",&varanandviap[year-PYEAR]);
     status=nc_inq_varid(ncidparam[year-PYEAR],"ndvim",&varndvim[year-PYEAR]);
     status=nc_inq_varid(ncidparam[year-PYEAR],"ndvimin",&varndvimin[year-PYEAR]);
     status=nc_inq_varid(ncidparam[year-PYEAR],"ndvimax",&varndvimax[year-PYEAR]);
     printf ("ndvimax=%d %d %d \n",year-PYEAR,ncidparam[year-PYEAR],varndvimax[year-PYEAR]);
     status=nc_inq_varid(ncidparam[year-PYEAR],"ndviap",&varndviap[year-PYEAR]);
   }
 for (l=0;l<NBLTS;l+=NBL)
   {
     printf ("l=%d\n",l);
     memset (tsmseuil,0,NBL*NBCTS*NVEG*sizeof(int));
     memset (nbvalseuil,0,NBL*NBCTS*NVEG*sizeof(int));
     memset (year1,0,NBL*NBCTS*NVEG*sizeof(int));
     memset (year2,0,NBL*NBCTS*NVEG*sizeof(int));
     for (v=0;v<NVEG;v++) for(l2=0;l2<NBL;l2++) for (c=0;c<NBCTS;c++) {diffout[v][l2][c]=anaav[v][l2][c]=0;nbvaldiff[v][l2][c]=0;}
     start[0]=PYEAR-2011;
     start[1]=l;
       start[2]=0;
     ncount[0]=DYEAR-PYEAR+1;
     ncount[1]=NBL;
     ncount[2]=NBCTS;
     status=ncvarget(ncidts,vartsmax,start,ncount,datatsmax); 
      ERRHAND(status);
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
	 for (v=0;v<NVEG;v++)
	   for (l2=0;l2<NBL;l2++)
	     for (c=0;c<NBCTS;c++) map[v][l2][c]=fillval;
	 start[0]=l;
	 ncount[0]=NBL;
	 start[1]=0;
	 ncount[1]=NBCTS;
	 status=ncvarget(ncidparam[year-PYEAR],varanats[year-PYEAR],start,ncount,danats);
   printf ("get1\n");
 	 ERRHAND(status);
	 start[0]=0;
	 ncount[0]=NVEG;
	 start[1]=l;
	 ncount[1]=NBL;
	 start[2]=0;
	 ncount[2]=NBCTS;
	 	 status=ncvarget(ncidparam[year-PYEAR],varanandvi[year-PYEAR],start,ncount,danandvi);
   printf ("get1\n");
 	 ERRHAND(status);
	 ERRHAND(status);
	  status=ncvarget(ncidparam[year-PYEAR],varanandviav[year-PYEAR],start,ncount,danandviav);
   printf ("get1\n");
 	 ERRHAND(status);
	 ERRHAND(status);
	status=ncvarget(ncidparam[year-PYEAR],varanandviap[year-PYEAR],start,ncount,danandviap);
   printf ("get1\n");
 	 ERRHAND(status);
	 ERRHAND(status);
	 status=ncvarget(ncidparam[year-PYEAR],varndvim[year-PYEAR],start,ncount,dndvim);
   printf ("get1\n");
   status=ncvarget(ncidparam[year-PYEAR],varndvimin[year-PYEAR],start,ncount,dndvimin);
	 printf ("get1 %d\n",varndvimax[year-PYEAR]);
	 status=ncvarget(ncidparam[year-PYEAR],varndvimax[year-PYEAR],start,ncount,dndvimax);
   printf ("get1\n");
 	 ERRHAND(status);
	 status=ncvarget(ncidparam[year-PYEAR],varndviap[year-PYEAR],start,ncount,dndviap);
   printf ("get1\n");
 	 start[0]=year-2011;
	 start[1]=l;
	 start[2]=0;
	 ncount[0]=1;
	 ncount[1]=NBL;
	 ncount[2]=NBCTS;
         status=ncvarget(ncidts,vartstamax,start,ncount,datatstamax);
	 ERRHAND(status);
	 for (v=0;v<NVEG;v++) 
	   for (l2=0;l2<NBL;l2++)
	     for (c=0;c<NBCTS;c++)
	       {
		 lat=80.-(l+l2)/10.;
		 if ((lat>65)||(lat<-65)) v2=v;
		 else  if ((lat>-35)&&(lat<35)) v2=12+v;
		 else v2=6+v;
		 v2=v;
		 tsm=datatsmax[year-PYEAR][l2][c];
		 tsmm=meanmaxts[l2][c];
		 if ((danandviav[v][l2][c]>-5)&&(danandviav[v][l2][c]<=5)&&(danandviap[v][l2][c]>-5)&&(danandviap[v][l2][c]<=5)&&(danats[l2][c]>0)&&(dndviap[v][l2][c]>0.2)&&(tsm>20)&&(tsm<65)&&(tsmm>20)&&(tsmm<45)&&(dndvimax[v][l2][c]>0.3)&&(meanmaxts[l2][c]<45)&&(dndvim[v][l2][c]>0.4))
		   {
		     // if ((v==0)) printf ("(%e %e %e)",danandviav[v][l][c],dndvim[l][c],danandviap[v][l][c]);	
		     ampli=dndvimax[v][l2][c];
		     if (dndvimax[v][l2][c]>0.2) diffvi=(danandviap[v][l2][c]-danandviav[v][l2][c])/dndvimax[v][l2][c]; else diffvi=(danandviap[v][l2][c]-danandviav[v][l2][c]);
		     if (diffvi<-0.2)
		       {
			 tsmseuil[v2][l2][c]+=datatsmax[year-PYEAR][l2][c];
			 if (year1[v2][l2][c]==0) year1[v2][l2][c]=year; else  year2[v2][l2][c]=year;
			 nbvalseuil[v2][l2][c]++;
		       }
		    
			 if ((tsm>35)&&(tsm<45)&&(tsmm>38)&&(tsmm<52)) map[v2][l2][c]=danandvi[v][l2][c];
			 // map[l][c]=tsm-tsmm;
			 /*			 hist2D[v2][(int)((tsmm-20.)*4.)][(int)((tsm-20.)*4.)]+=diffvi;
						 nval2D[v2][(int)((tsmm-20.)*4.)][(int)((tsm-20.)*4.)]++;*/
			 vdiff=(diffvi)*100+80;
			 if (vdiff<0) vdiff=0;
			 if (vdiff>159) vdiff=159;
			 hist2D[v2][(int)vdiff][(int)((tsm-20.)*4.)]++;
			 histom[v2][(int)((tsm-20.)*4.)]+=diffvi;
			 histom2[v2][(int)((tsm-20.)*4.)]+=danandvi[v][l2][c];
			 histom3[v2][(int)((tsm-20.)*4.)]+=datatstamax[l2][c];
			 nval[v2][(int)((tsm-20.)*4.)]++;
			 varm[v2][(int)((tsm-20.)*4.)]+=diffvi*diffvi;
		     if (diffvi<-0.2)
		       {
			 diffout[v][l2][c]+=diffvi;
			 anaav[v][l2][c]+=danandviav[v][l2][c];
			 nbvaldiff[v][l2][c]++;
		       }

		   }
		 if ((danandviav[v][l2][c]>-1)&&(danandviav[v][l2][c]<=1)&&(danandviap[v][l2][c]>-1)&&(danandviap[v][l2][c]<=1))
		   {
		     ats=danats[l2][c];
		     if ((ats>-15)&&(ats<15)&&(dndvim[v][l2][c]>0.3))
		       {
			 histts[v2][(int)((ats+15.)*4.)]+=diffvi;
			 nvalts[v2][(int)((ats+15.)*4.)]++;
		       }
		   }
	       }
	 ncount[0]=1;
	 ncount[1]=NVEG;
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
	 for (v=0;v<NVEG2;v++)
	   {
	     if (nbvalseuil[v][l2][c]>0)  tsmseuil[v][l2][c]/=nbvalseuil[v][l2][c]; else tsmseuil[v][l2][c]=nbvalseuil[v][l2][c]=year1[v][l2][c]=year2[v][l2][c]=fillval;
	     if (nbvaldiff[v][l2][c]>0)  {diffout[v][l2][c]/=nbvaldiff[v][l2][c]; anaav[v][l2][c]/=nbvaldiff[v][l2][c];} else diffout[v][l2][c]=anaav[v][l2][c]=fillval;
	   }
     //    for (passe=0;passe<5;passe++)
       {
	 pointtraite=0;
	 for (l2=0;l2<NBL;l2++)
	   for (c=0;c<NBCTS;c++)
	     for (v=0;v<NVEG2;v++)
	       {
		 if (datatsmax[0][l2][c]>-50) 
		   {
		     valr=nbvalr=0;
		     //for (rayon=5;(rayon<50)&&(diffout[v][l2][c]<0);rayon+=5);
		     rayon=5;
		     {
		       for (l3=-rayon;l3<rayon;l3++)
			 for (c2=-rayon;c2<rayon;c2++)
			   {
			     lr=l2+l3;
			     if (lr<0) lr=0;
			     if (lr>=NBL) lr=NBL-1;
			     cr=c+c2;
			     if (cr<0) cr=0;
			     if (cr>=NBCTS) cr=NBCTS-1;
			     if (diffout[v][lr][cr]>fillval)
			       {
				 if (diffout[v][lr][cr]>fillval) {valr+=diffout[v][lr][cr];nbvalr++;}  
			       }
			   }
		     }
		     if (nbvalr>20) diffout2[v][l2][c]=valr/nbvalr; else diffout2[v][l2][c]=fillval;
		   }
		 else diffout2[v][l2][c]=fillval;
	       }
       }
     ncount[0]=NVEG;
     ncount[1]=NBL;
     ncount[2]=NBCTS; 
     start[1]=l;start[0]=start[2]=0;
     printf ("put1\n");
     status=ncvarput(ncido,varmtmax,start+1,ncount+1,meanmaxts);
     printf ("put2\n");
     status=ncvarput(ncido,vartseuil,start,ncount,tsmseuil);
     printf ("put3\n");
     status=ncvarput(ncido,varnbyear,start,ncount,nbvalseuil);
     printf ("put4\n");
     status=ncvarput(ncido,varyear1,start+1,ncount+1,datatsmax[0]);
     printf ("put5\n");
     status=ncvarput(ncido,varyear2,start+1,ncount+1,year2);
     status=ncvarput(ncido,vardiff,start,ncount,diffout2);
    status=ncvarput(ncido,varanaav,start,ncount,anaav);
   }
  for(v=0;v<NVEG2;v++)
   for (l=0;l<NBVAL2;l++)
     for (c=0;c<NBVAL2;c++)
     if (nval[v][c]>0) hist2D[v][l][c]/=nval[v][c];
 for (t=0;t<DYEAR-PYEAR+1;t++) time[t]=t*365+182;
 start[0]=0;
 ncount[0]=DYEAR-PYEAR+1;
 ncvarput(ncido,ntimep,start,ncount,time);
 start[0]=start[1]=start[2]=0;
 ncount[0]=NVEG;
 ncount[1]=NBVAL2;
 ncount[2]=NBVAL2;
 printf ("put6\n");
 ncvarput(ncido,nh2D,start,ncount,hist2D);
 for (v=0;v<NVEG2;v++)
   {
     for (i=nbptot=0;i<NBVAL2;i++) nbptot+=nval[v][i];
     for (i=0;i<NBVAL2;i++) 
       if (nbptot>0)
	 {
	   histom[v][i]/=nval[v][i];
	   histom2[v][i]/=nval[v][i];
	   histom3[v][i]/=nval[v][i];
	   varm[v][i]=1.862*sqrt(varm[v][i]/nval[v][i]-histom[v][i]*histom[v][i])/sqrt(nval[v][i]);
	 }
    for (i=0;i<NBVAL2;i++) 	   
       if ((nbptot>0)&&(nval[v][i]/nbptot>0.0005))
	 {
	   hist1=hist2=hist3=varm1=0;
	   if ((i>=2)&&(i<NBVAL2-2)) 
	     {
	       for (i2=-2;i2<=2;i2++)
		 {
		   hist1+=histom[v][i+i2]/5;
		   hist2+=histom2[v][i+i2]/5;
		   hist3+=histom3[v][i+i2]/5;
		   varm1+=varm[v][i+i2]/5;
		 }
	     }
	   else
	     {
	       hist1=histom[v][i];
	       hist2=histom2[v][i];
	       hist3=histom3[v][i];
	       varm1=varm[v][i];
	     }
	   fprintf (fout[v],"%f %lf %lf %lf %lf %lf %f\n",(float)i/4.+20,hist1,hist2,hist3,hist1-varm1,hist1+varm1,nval[v][i]);
	 }
       else 
	 fprintf (fout[v],"%f 0 0 0 0 0 0\n",(float)i/4.+20);  
   }
     for (v=0;v<NVEG2;v++)
       for (i=0;i<NBVALTS;i++) 
	 if (nvalts[v][i]>100)
	   {
	     histts[v][i]/=nvalts[v][i];
	     fprintf (foutts[v],"%f %d\n",histts[v][i],nvalts[v][i]);
	   }
	 else 
	   fprintf (foutts[v],"0 0\n");  
 status=ncclose(ncido);
}	
