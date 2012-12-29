/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Main program
author: Taposh Dutta Roy - dutttap@iit.edu
creation date: Wednesday, August 25, 2004
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include <stdio.h>
#include <stdlib.h>
#include "parse.h"
#include "simu.h"
#include "last.h"
#include <string.h>
#define REG "_regevol.raw"
#define PLT "_regevol.plt"
#define DAT "_regevol"
#include <unistd.h>
#define _GNU_SOURCE
#include <getopt.h> 
#include "misce.h"
#include "last.h"     
#define MAXSTRLEN (80)
#define TIME 0
#define PART 1
#define DISP 2
#define TOTAL 3
#define ASCIIC 48
#define CURRENTNO 3
//#define PI        (3.14159265359)
//#define ABS(a)    ((a)<0 ? -(a) : (a))
//#define ps        (1e-12)
//#define GHz       (1e9)




/*------------------------------------------------------------------------*/
static struct option outtype[]=                                        /* graphics programs type list */
{ 
  {"tecplot", no_argument,0,TECPLOT},
 {"generic", no_argument,0,GENERIC}, 
 {"help",    optional_argument,0,HELPO}, 
 {"file",    required_argument,0,FILEN}, 
 {"output",  required_argument,0,OUTN},  
 {"list",    0,0,LISTV},  
 {0,0,0,0} 
};  

typedef  struct {
  Toplot *tpl;
  double t1;
  int numcontreg;
  int curfilenum;
  char *filename;
 }simpts;

typedef struct {
  double dt;
  double t1;
  double t2;
  char *terminal;
  double vs1;
  double vs2;
  double vg1;
  double vg2;
  double vd1;
  double vd2;
  double deltaV;
  int totfilenum;
  int totfilereq; 
}params;

typedef struct
{
  MYFLOAT re;
  MYFLOAT im;
}
COMPLEX;

/*------------------------------------------------------------------------*/



/*-----------------------------------------------------------------------*/
static Toplot *getvar     (Toplot*, UINT1*, int, char **,int);
static void usage (char**);
static void scan_file(char*);
static MYFLOAT** convert_plt(simpts*);
static void convert_dat(Toplot*,char *,int);
static int test_compatability(Toplot*);
static MYFLOAT** getcurrents(Toplot*,int,char*);
static void writefile_dat(char*,float**,int,Point*,int);
void check1plot(Toplot*,const char * ,UINT1,int);
static simpts *read_sim(const char *,const char *,int ,char **,params* );
static simpts *checkplot(simpts*,int,char**,int);
static MYFLOAT***  finalcurrents(simpts*,int);
static MYFLOAT autocorrelation(MYFLOAT*,int);
static int stepanalysis(MYFLOAT ***,params *);
static int  myfileplot(MYFLOAT ***,simpts *,int );
     


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Main program
author: Taposh Dutta Roy - dutttap@iit.edu
creation date: Wednesday, August 25, 2004
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int main (int argc, char **argv) 
{  									
  const char simrul[]="simulation.rul";
  char *simdef;
  MYFLOAT ***currents=NULL,***currents0=NULL,***currents1=NULL;
  int ii=0,numfiles=0;
  
  


  simpts *dataout=NULL,*datain=NULL;
  params *genparams = NULL;
 

  if (argc < 3 ) {
     simdef=(char*)(calloc(strlen(*(argv+1)),sizeof(char)));
     strcpy(simdef,*(argv+1));// No flag specified
  }
   else {
     simdef=(char*)(calloc(strlen(*(argv+2)),sizeof(char)));
     strcpy(simdef,*(argv+2));
   }

 
  if (strlen(simdef)==0) {
    IERRS("Please Enter the File name");
  }
 
  genparams = (params*)calloc(1,sizeof(params));

  datain=read_sim(simdef,simrul,argc,argv,&(*genparams));
  if(datain ==0 ){
    IERRS("FreqAnalysis:datain cannot be zero.\n");
  }
  dataout = checkplot(datain,argc,argv,(genparams->totfilereq));
  if(dataout ==0 ){
    IERRS("FreqAnalysis:dataout cannot be zero.\n");
  }
  currents =finalcurrents(dataout,(genparams->totfilereq));
  if(currents ==0 ){
    IERRS("FreqAnalysis:currents cannot be zero.\n");
  }
  /* Plot the structure */
  numfiles = (genparams->totfilereq);
  
  for (ii=0;ii<numfiles;ii++)
    {
      myfileplot(currents,dataout,ii);
    }
  stepanalysis(currents,genparams);
  
  free(genparams);
  exit(0);
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
author: Taposh Dutta Roy  - dutttap@iit.edu
creation date: Thursday, Sept 11,2004
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

static int stepanalysis(MYFLOAT ***currents,params *genparams)
{
  MYFLOAT ISS1;
  MYFLOAT ISS2;
  static COMPLEX     *cur_four0=NULL,   *cur_four=NULL,   *cur_four1=NULL;
  static COMPLEX     *four_smooth0=NULL,*four_smooth1=NULL;
  int N1 = 0;
  int N2 =0,ii=0,totfiles=0;
  MYFLOAT **cur0=NULL,**cur1=NULL,**cur=NULL;
  

  /*Parameters */ 
  N1 =(int) ((genparams->t1)/(genparams->dt));
  N2 = (int)((genparams->t2)/(genparams->dt));
  
 /* compute autocorrelation to get ISS1 and ISS2 */
  printf("computing autocorrelation: ");
    
  if (genparams->totfilereq ==1){
    WERRS("\n Two files needed in PLATEAU mode");
    // ISS1 = autocorrelation(*(*(currents+0)+1),N1);
    exit(0);
  }
    /* Passing Particle Currents in All Files */
  cur0 = *(currents+0);
  ISS1 = autocorrelation(*(cur0+1),N1); 
  cur1 = *(currents+1);
  ISS2 = autocorrelation(*(cur1+1),N2); 
  printf("ISS1=%e\tISS2=%e\tdeltaI=%e\n",ISS1,ISS2,ISS2-ISS1);
  
 }
 
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
author: Taposh Dutta Roy  - dutttap@iit.edu
creation date: Thursday, Sept 11,2004
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

static int  myfileplot(MYFLOAT ***currents,simpts *dataout,int filereq)
{
    FILE *myoutput;
    int ii=0,jj=0;
    int numcontreg,numpts,totcols;
    char *outname;
    
    numcontreg = (dataout+filereq)->numcontreg;
    numpts= (dataout+filereq)->tpl->nptx;
    
    outname = (char*) calloc(strlen((dataout+filereq)->filename),sizeof(char));
    strcpy(outname,((dataout+filereq)->filename));
    myoutput = fopen(outname,"ab");
   
    totcols =((numcontreg * 3)+1);
    
    for (jj=0;jj<numpts;jj++)
      {
      for (ii=0;ii<totcols;ii++){
	fprintf( myoutput,MYFORMAT,*(*(*(currents+filereq)+ii)+jj));
	 fprintf(myoutput,"\t");
      }
      fprintf(myoutput,"\n");
    }
    
    fclose(myoutput);
    free(outname);
}



/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
author: Taposh Dutta Roy - dutttap@iit.edu
creation date: 
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*---------------------------------------------------------------------------------------------------*/
/* take the Fourrier transform of the current, the discrete frequencies are stored in four0, the     */
/* the Fourier transform of the total current variation  is stored in four0+1                        */
/*---------------------------------------------------------------------------------------------------*/
/*void fourier(COMPLEX *four0, MYFLOAT *fl0, MYFLOAT ISS1, MYFLOAT ISS2)
{
  int kk,nn;
  MYFLOAT *fl1=NULL,   *fl=NULL;
  COMPLEX *four1=NULL, *four=NULL;
  MYFLOAT omegak;

  fl1=fl0+N;
  four1=four0+N/2;
  
  for (four=four0,kk=0;four<four1;four++,kk++) {
    omegak=2*PI*kk/(T*ps);
    four->re=0.0;
    four->im=0.0;
    for (fl=fl0,nn=0;fl<fl1;fl++,nn++) {
      four->re+=(*fl-ISS1)*cos(omegak*nn*DT*ps);
      four->im-=(*fl-ISS1)*sin(omegak*nn*DT*ps);
    }
    four->re*=(DT*ps);
    four->im*=(DT*ps);

    if (kk!=0){      
      four->im-=(ISS2-ISS1)/omegak;
    }
    else{
      four->re-=T*ps*(ISS2-ISS1);
    }

  }
}*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
author: Taposh Dutta Roy - dutttap@iit.edu
creation date: 
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


static MYFLOAT autocorrelation(MYFLOAT *fl0, int Ntot)
{
  MYFLOAT ISS;
  MYFLOAT sum=0.0; 
  MYFLOAT *fl=NULL;
  int nn, nprime;

  nprime=((int)(Ntot/3));
  for (nn=0;nn<(Ntot-nprime);nn++){
    fl=fl0+nn;
    sum+=(*(fl))*(*(fl+nprime));
    // Taposh-Debugging 
    if (sum >= 1e+4){
      printf("Look for error in sum %d",(nprime));
    }
    //----------------
  }
  if(Ntot-nprime<=0) {
    printf("autocorrelation: error, division by zero\n");
    exit(0);
  }
  sum/=((int)(Ntot-nprime));
  if(sum<0) printf("autocorrelation: error, square root of a negative number\n");
  ISS=sqrt(ABS(sum));
  return(ISS);
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
author: Taposh Dutta Roy  - dutttap@iit.edu
creation date: Thursday, Sept 9,2004
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

static MYFLOAT***  finalcurrents(simpts *dataout,int totfilereq)
{
  int nvar =0,numpt=0,ii=0,jj=0;
  int numcontreg =0,totcols =0;
  MYFLOAT   ***cur,***cur0=NULL,***cur1=NULL,***cur2=NULL;
    

  cur2=(MYFLOAT***)(calloc(totfilereq+1,sizeof(MYFLOAT**)));
  if(cur2==NULL){
    IERRS("finalcurrents: Memory Allocation Failed for the pointer \n");
  }
  for (ii=0;ii<(totfilereq);ii++)
    {
      numcontreg = ((dataout+ii)->numcontreg);
      nvar =((dataout+ii)->tpl->numvars);
      numpt =((dataout+ii)->tpl->nptx);
      totcols = (((nvar/2)*CURRENTNO)+1);
      
      if( nvar ==0 ){
	IERRS("convert_plt: Number of Variables cannot be Zero.\n");
      }
      if (totcols < ((numcontreg * CURRENTNO) +1)) {
	WERRS("getcurrents: Raw file has less number of Regions than device file\n");
      }
      *(cur2+ii)=(MYFLOAT**)(calloc(totcols,sizeof(MYFLOAT*)));
      if((cur2+ii)==NULL){
	IERRS("getcurrents: Memory Allocation Failed for the pointer \n");
      }
      cur1=(cur2+ii)+totcols;
      for(cur=cur2;cur<cur1;cur++){
	*(cur+ii)=(MYFLOAT**)(calloc(numpt,sizeof(MYFLOAT*)));
	if(cur==NULL){
	  IERRS("getcurrents: Memory Allocation Failed for the pointer\n");
	}
      }
    }
  
  
  for (ii=0;ii<totfilereq;ii++){

 *(cur+ii) = convert_plt(dataout+ii);

  }
  
  return cur;
  
}




/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
author: Taposh Dutta Roy  - dutttap@iit.edu
creation date: Thursday, Sept 2,2004
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

static simpts *checkplot(simpts *sp,int argc,char **argv,int numfiles)
{
  int ii=0;
  UINT1 nvout=0; 
  Toplot *tpout=NULL;
  simpts *newsp=NULL;
  tpout=(Toplot*)(calloc(1,sizeof(Toplot)));
  newsp=(simpts*)(calloc(numfiles,sizeof(simpts)));
  

  
  for(ii=0;ii<=(numfiles-1);ii++)
    {
      tpout =getvar(((sp+ii)->tpl),&nvout,argc,argv,((sp+ii)->curfilenum));
       (newsp+ii)->tpl=(Toplot*)(calloc(1,sizeof(Toplot)));
      memcpy(((Toplot*)(newsp+ii)->tpl),(Toplot*)tpout,sizeof(Toplot));
      (newsp+ii)->t1 = (sp+ii)->t1;
      (newsp+ii)->numcontreg = (sp+ii)->numcontreg;
      (newsp+ii)->curfilenum = (sp+ii)->curfilenum;
      (newsp+ii)->filename = (sp+ii)->filename;
      
    }
  
  return newsp;
}  


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Returns a list of selected variables codes to insert in structure.
Parameter numcalsout returns the number of class of output variables.
Parameter frame returns code of included frame.
author: claudio berselli (berselli@neumann.ece.iit.edu)
creation date: Wensday, September 6, 2000
first revision: Monday, January 1, 2001 (marco saraniti)
second revision: Taposh Dutta Roy Current Saturday, Jul24 5:33
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

static Toplot *getvar(Toplot *tp,     
		  UINT1 *nvout, 
		  int   argc,
		  char  **argv,int fnum)           
{
  struct option *opt=NULL, *opt0=NULL, *opt1=NULL;
  Toplot      *tpout=NULL;
  Point       *po=NULL, *pi=NULL, *po1=NULL;
  Pline         *pli=NULL, *plo=NULL, *plo1=NULL;
  MYFLOAT       *valo=NULL, *vali=NULL;
  UINT2         namel=0;
  UINT4         tout=0, numpt=0;
  FLAG          *thisone=NULL, gefunden=NO;
  char          **vout=NULL, **lout=NULL, *sout=NULL, *sopt0=NULL, *sopt=NULL;
  char          **vname=NULL, **lopt=NULL;
  char          **vlopt=NULL, *vsopt=NULL;
  int           lo=0, cout=0, ch=0;

  opt0=opt=(struct option*)(calloc(tp->numvars+1,sizeof(struct option)));
  opt1=opt0+tp->numvars;
  thisone=(FLAG*)(calloc(tp->numvars,sizeof(FLAG)));
  sopt=sopt0=(char*)(calloc(tp->numvars+1,sizeof(char)));
  vlopt=tp->varlopt;
  vsopt=tp->varsopt;
  sopt=sopt0;
  for(opt=opt0; opt<opt1; opt++, sopt++, vlopt++, vsopt++) {
    namel=strlen(*vlopt);
    opt->name=(char*)(calloc(namel+1,sizeof(char)));
    sprintf((char*)(opt->name),"%s",*vlopt);
    *sopt=opt->val=(int)(*vsopt);      
    thisone[opt-opt0]=NO;
  }
  
  fnum = ASCIIC +fnum;/* Adding 48 for char value of a number*/
  optind =0;
   while(optind<argc) {

     ch=getopt(argc,argv,"1234567890");        
       if (ch == fnum){
	 fprintf(OUTPTR,"\t%d\t\n",(fnum-48));
      while(optind<argc) {          
	cout=getopt_long(argc,argv,sopt0,opt0,&lo); 

	  if (cout!='?') {
	    opt=opt0; gefunden=NO;
	    while((opt<opt1)&&(gefunden==NO)) {
	      if((opt->val==cout) && (thisone[opt-opt0]==NO)) {
		thisone[opt-opt0]=YES;
		gefunden=YES;
		tout++;
	      }
	      opt++;
	    }
	  }
	  else {
	    break;
	  }
	      	    
	}
      }
   }
 
 
 
  if(tout !=0) {
    /* for(tout=0; tout<tp->numvars; tout++) thisone[tout]=YES;*/
   
  
  tpout=(Toplot*)(calloc(1,sizeof(Toplot))); 
  tpout->dimension=tp->dimension;
  tpout->nptx=tp->nptx;
  tpout->npty=tp->npty;
  tpout->nptz=tp->nptz;
  tpout->numvars=tout;
  tpout->title=(char*)(calloc(strlen(tp->title)+1,sizeof(char)));
  sprintf(tpout->title,"%s",tp->title);
  tpout->xlabel=(char*)(calloc(strlen(tp->xlabel)+1,sizeof(char)));
  sprintf(tpout->xlabel,"%s",tp->xlabel);
  tpout->ylabel=(char*)(calloc(strlen(tp->ylabel)+1,sizeof(char)));
  sprintf(tpout->ylabel,"%s",tp->ylabel);
  tpout->zlabel=(char*)(calloc(strlen(tp->zlabel)+1,sizeof(char)));
  sprintf(tpout->zlabel,"%s",tp->zlabel);
  vout=tpout->varname=(char**)(calloc(tpout->numvars,sizeof(char*)));
  lout=tpout->varlopt=(char**)(calloc(tpout->numvars,sizeof(char*)));
  sout=tpout->varsopt=(char*)(calloc(tpout->numvars,sizeof(char)));
  vname=tp->varname; lopt=tp->varlopt; sopt=tp->varsopt;
  for(vname=tp->varname; vname < (tp->varname+tp->numvars); vname++,lopt++,sopt++) { 
    if(thisone[vname-tp->varname]==YES) {
      namel=strlen(*vname);        
      *vout=(char*)(calloc(namel+1,sizeof(char)));
      sprintf(*vout,"%s",*vname);
      namel=strlen(*lopt);
      *lout=(char*)(calloc(namel+1,sizeof(char)));
      sprintf(*lout,"%s",*lopt);
      *sout=*sopt;
#     if(VERBOSE)
      DOVE; fprintf(OUTPTR,"getvar: extracting \"%s\"for \"%d\"\n",*vout,fnum);
#     endif
      vout++; lout++; sout++;
    }
  }
  switch(tp->dimension) {
  case 1: 
    if(tp->nptx==0) {IERRS("getvar: nptx=0 in 1-D plot data.");}
    if(tp->npty!=0) {IERRS("getvar: npty!=0 in 1-D plot data.");}
    if(tp->nptz!=0) {IERRS("getvar: nptz!=0 in 1-D plot data.");}
    numpt=tp->nptx; 
    break;
  case 2: 
    if(tp->nptx==0) {IERRS("getvar: nptx=0 in 1-D plot data.");}
    if(tp->npty==0) {IERRS("getvar: npty=0 in 2-D plot data.");}
    if(tp->nptz!=0) {IERRS("getvar: nptz!=0 in 2-D plot data.");}
    numpt=tp->nptx*tp->npty; 
    break;
  case 3: 
    if(tp->nptx==0) {IERRS("getvar: nptx=0 in 3-D plot data.");}
    if(tp->npty==0) {IERRS("getvar: npty=0 in 3-D plot data.");}
    if(tp->nptz==0) {IERRS("getvar: nptz=0 in 3-D plot data.");}
    numpt=tp->nptx*tp->npty*tp->nptz; 
    break;
  default:
    IERRS("getvar: wrong data dimension.");
    break;
  }
  tpout->pt0=po=(Point*)(calloc(numpt,sizeof(Point)));
  po1=tpout->pt0+numpt;
  for(pi=tp->pt0, po=tpout->pt0; po<po1; po++, pi++) {
    po->x=pi->x;
    po->y=pi->y;
    po->z=pi->z;
    valo=po->value=(MYFLOAT*)(calloc(tpout->numvars,sizeof(MYFLOAT)));
    for(vali=pi->value; vali<(pi->value+tp->numvars); vali++) {
      if(thisone[vali-pi->value]==YES) {
        if(valo==(po->value+tpout->numvars)) {
	  IERRS("getvar: pointer out of range.");
	}
	*valo++=*vali;
      }
    }
    if(valo!=(po->value+tpout->numvars)) {
      IERRS("getvar: pointer out of range.");
    }
  }
  tpout->numpl=tp->numpl;
  if(tpout->numpl) {
    tpout->pl0=(Pline*)(calloc(tpout->numpl,sizeof(Pline))); 
    plo1=tpout->pl0+tpout->numpl;
    for(pli=tp->pl0, plo=tpout->pl0; plo<plo1; plo++, pli++) { 
      memcpy((Pline*)(plo),(Pline*)(pli),sizeof(Pline));
      po=plo->pt0=(Point*)(calloc(plo->numpt,sizeof(Point)));
      memcpy((Point*)(po),(Point*)(pli->pt0),sizeof(Point)*plo->numpt);
    }
  }
  for(opt=opt0; opt<opt1; free((char*)(opt->name)), opt++);
  FREE(sopt0);
  FREE(thisone);
  relax_plot(tp);
  return(tpout);
  }
  return(0);
}




/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
author: Taposh Dutta Roy  - dutttap@iit.edu
creation date: Thursday, Sept 2,2004
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

static simpts *read_sim(const char *def,const char *rul,int argc,char **argv,params *genparams)
{
  FILE         *myoutput;
  unsigned int numentry=0;
  unsigned int ii=0; 
  double       *rtmp=NULL;
  double       *r1tmp=NULL;
  double       *r2tmp=NULL;
  char         **stmp=NULL;
  char         *s2tmp=NULL,*s1tmp=NULL,*argchar=NULL;
  double       dt,vd1,vd2,deltaV;
  double       t1,t2,t10; 
  char         *str2,*str1;
  char         *terminal;
  double       vs1,vs2,vg1,vg2; 
  int          fnamelen, numcontreg =0;
  char         *outname=NULL, *filename=NULL,*fname=NULL,*fname1=NULL,*outname1=NULL;
  Toplot        *datain,*dataout=NULL;
  simpts       *data;
  // UINT1        nvout=0; 
  int          cout =0,filereq=0,ch,lo=0,fnum=1,mynum=0,jj=0;
  UINT1        outprog=0;
  FLAG         gefunden =NO;
  simpts *totdata=NULL;/*totdata0=NULL, *totdata1=NULL;*/
 
  
 
   optarg = 0;
  
    stmp=sval(rul,def,15,1,&numentry); /* get the device.rul file */
    if(numentry==0) {
      SERRS("read_sim: .");
    }
    s1tmp=(char*)(calloc(numentry,sizeof(char)));
    if(s1tmp==NULL){
      IERRS("read_sim: error allocating string pointer memory.");
    }
    strcpy(s1tmp,*stmp);
    if(s1tmp==NULL){
      IERRS("read_sim: error allocating string pointer memory ");
    }
    

    stmp=sval(rul,def,15,2,&numentry); /* get the device.def */
    if(numentry==0) {
      SERRS("read_sim: .");
    }
    s2tmp=(char*)(calloc(numentry,sizeof(char)));
    if(s2tmp==NULL){
      IERRS("read_sim: error allocating string pointer memory.");
    }
    strcpy(s2tmp,*stmp);
    if(s2tmp==NULL){
      IERRS("read_sim: error allocating string pointer memory ");
    }

 
    stmp=sval(s1tmp,s2tmp,38,1,&numentry);
    if(numentry==0) {
      SERRS("read_sim: NO '2d Contact Regions' entry found.");
    }
    else {
      numcontreg = numentry;
    }
    free(s1tmp);
    free(s2tmp);
    
    /* reading Poisson Time */
    rtmp=rval(rul,def,7,1,&numentry);
    if(numentry==0)  {
      SERRS("read_sim: no 'poisson time' entry found.");
    }
    dt = *rtmp; 
    
    /* Reading Duration of Steps T1 and T2*/        
    r1tmp=rval(rul,def,9,7,&numentry); 
    if(numentry == 0)  {
      SERRS("read_sim:NO  'duration time' entry found.");
    }
    t1 = *r1tmp;
    t2 = *(r1tmp+1);

    /* reading Drain Potential*/
    stmp=sval(rul,def,9,10,&numentry); 
    if(numentry==0) {
      SERRS("read_sim: NO 'Drain Potential' entry found.");
    }
    s2tmp=(char*)(calloc(numentry,sizeof(char)));
    if(s2tmp==NULL) {
      IERRS("read_sim: error allocating string pointer memory.");
    }
    strcpy(s2tmp,*stmp);
    if(s2tmp==NULL) {
      IERRS("read_sim: error allocating string pointer memory ");
    }
    vd1 = atof(s2tmp);/*vd1*/
    free(s2tmp);
    
    s1tmp=(char*)(calloc(numentry,sizeof(char)));
    strcpy(s1tmp,*(stmp+1));
    if(s1tmp==NULL) {
      IERRS("read_sim: error allocating string pointer memory ");
    }
    vd2 = atof(s1tmp);/*vd2*/
    free(s1tmp);
    
    stmp=sval(rul,def,9,9,&numentry);
    if(numentry==0) {
      SERRS("read_sim: NO 'Drain Potential' entry found.");
    }
    s2tmp=(char*)(calloc(numentry,sizeof(char)));
    if(s2tmp==NULL){
      IERRS("read_sim: error allocating string pointer memory.");
    }
    strcpy(s2tmp,*stmp);
    if(s2tmp==NULL) {
      IERRS("read_sim: error allocating string pointer memory ");
      }
    vg1 = atof(s2tmp);/*vg1*/
    free(s2tmp);
    
    s1tmp=(char*)(calloc(numentry,sizeof(char)));
    strcpy(s1tmp,*(stmp+1));
    if(s1tmp==NULL) {
      IERRS("read_sim: error allocating string pointer memory ");
    }
    vg2 = atof(s1tmp);/*vg2*/
    free(s1tmp);
    
    stmp=sval(rul,def,9,8,&numentry);
    if(numentry==0){
      SERRS("read_sim: NO 'Drain Potential' entry found.");
    }
    s2tmp=(char*)(calloc(numentry,sizeof(char)));
    if(s2tmp==NULL){
      IERRS("read_sim: error allocating string pointer memory.");
    }
    strcpy(s2tmp,*stmp);
    if(s2tmp==NULL){
      IERRS("read_sim: error allocating string pointer memory ");
    }
    vs1 = atof(s2tmp);/*vs1*/
    free(s2tmp);
    
    s1tmp=(char*)(calloc(numentry,sizeof(char)));
    strcpy(s1tmp,*(stmp+1));
    if(s1tmp==NULL) {
      IERRS("read_sim: error allocating string pointer memory ");
      }
    vs2 = atof(s1tmp);/*vs2*/
    free(s1tmp);
    
    if (vd1-vd2 != 0.0) {  /*  Determine the Type of Simulation */
      VERBO("read_sim: Drain Type Simualtion ."); 
      deltaV =(vd1-vd2);
      terminal = "Drain";
    } else {
      if (vg1-vg2 !=0) 
	{
	  VERBO("read_sim:Gate Type Simualtion ");
	  deltaV =(vg1-vg2);
	  terminal = "Gate"; 
	}
      else
	{
	  if (vs1-vs2 !=0) 
	    {
	        VERBO("read_sim:Source Type Simualtion ."); 
	        terminal = "Source";
        	deltaV = (vs1-vs2);
	    }
	  else
	    {
	      IERRS("read_sim:No difference Potential in any terminal");  
	    }			    
	}   
    }
    
  
     
    stmp=sval(rul,def,9,1,&numentry);     /* Read the label or base-name of the output Files*/
    
    if(numentry==0) 
      {
	SERRS("read_sim: no Label or BASE-NAME  entry found.");
      }
 
    fprintf(OUTPTR,"*******************************************\n");
    fprintf(OUTPTR,"\tList of RAW  Files in Definition File\n");
    fprintf(OUTPTR,"*******************************************\n");
    for(ii = 0; ii<(numentry);ii++)
      {
	
	fnamelen = (strlen(*(stmp+ii)) + strlen(REG) +1);
	fname1 = (char*)(calloc(fnamelen,sizeof(char)));
	outname1 = (char*)(calloc(fnamelen,sizeof(char)));
	
	if(fname1==NULL) {
	  IERRS("read_sim: error allocating string pointer memory.");
	}
	strcpy (fname1,*(stmp+ii));
	strcat(fname1,REG);
   	strcpy(outname1,*(stmp+ii));
	
	fprintf(OUTPTR,"\t%d\t%s\n",(ii+1),fname1);
	
      }
    

    free(fname1);
    free(outname1);
    
    
    filereq=calnumfiles(argc,argv,numentry);    
    fprintf(OUTPTR,"Number of Files Requested: \t%d\t\n",filereq);
	
      ((genparams)->dt) =dt;
      ((genparams)->t1)= t1;
      ((genparams)->t2)= t2;
      ((genparams)->terminal)= (char*)(calloc(strlen(terminal),sizeof(char)));
      ((genparams)->terminal)= terminal;
      ((genparams)->deltaV)= deltaV;
      ((genparams)->vs1)= vs1;
      ((genparams)->vs1)= vs1;
      ((genparams)->vs2)= vs2;
      ((genparams)->vg1)= vg1;
      ((genparams)->vg2)= vg2;
      ((genparams)->vd1)= vd1;
      ((genparams)->vd2)= vd2;
      ((genparams)->totfilenum)= numentry;
       ((genparams)->totfilereq)= filereq;
       
       outprog =0;
       totdata=(simpts*)calloc(filereq+1,sizeof(simpts));
       
       
       for(ii = 0; ii<(numentry);ii++)
	 {
	   fnamelen = (strlen(*(stmp+ii)) + strlen(REG) +1);
	   fname = (char*)(calloc(fnamelen,sizeof(char)));
	   outname = (char*)(calloc(fnamelen,sizeof(char)));
	
	   if(fname==NULL) {
	     IERRS("read_sim: error allocating string pointer memory.");
	   }
	   strcpy (fname,*(stmp+ii));
	   strcat(fname,REG);
	   strcpy(outname,*(stmp+ii));
	   
	   
	   if(argc==1) { usage(argv); EXIT; }
	   cout=0,opterr=0;
	   
	   while(optind<argc) {
	     
	     cout=getopt_long(argc,argv,"o:h:f:l",outtype,&lo);
	     
	     switch(cout) {
	     case 'h':
	       usage(argv); 
	    if(optarg) scan_file(optarg); EXIT;
	    break;
	     case 'f':
	       fname = fname;  break;
	     case TECPLOT:
	       outprog=TECPLOT; break;
	     case GENERIC:
	       outprog=GENERIC; break;
	     case '?':
	       break;
	     }
	   }
	   optind =0;
	   /*if (outprog == 5){
	     strcat(outname,PLT);
	}
	   else if (outprog ==1) {
	     strcat(outname,DAT);
	     }*/
	   
	
	fnum =ii+1;
	fnum = ASCIIC +fnum;
	while(optind<argc) {
	  ch=getopt(argc,argv,"1234567890");        
	  if (ch == fnum){
	    //datain = (Toplot *)(calloc(1,sizeof(Toplot)));
	    datain = convRAW4d(fname);
	    ((totdata+mynum)->tpl)=(calloc(1,sizeof(Toplot)));
	    ((totdata+mynum)->tpl) =datain;
	    ((totdata+mynum)->t1)= t1;
	    ((totdata+mynum)->numcontreg)=numcontreg;
	    ((totdata+mynum)->curfilenum)= (fnum - ASCIIC);
	    ((totdata+mynum)->filename)= (char*)(calloc(fnamelen,sizeof(char)));
	    ((totdata+mynum)->filename)= strcat(outname,PLT); 
	    mynum =mynum +1;
	    break;
	  }
	}

      }
    free(fname);
    // free(outname);
    return totdata;
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Calculates the number of files requested by the user
author: Taposh Dutta Roy  - dutttap@iit.edu
creation date: Thursday,August  25, 2004
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int calnumfiles(int argc,char **argv,int num)
{
  int ii=0,filereq=0,fnum=0,ch=0;
  opterr=0;
  
  for (ii=1;ii<=num;ii++)
    {
      fnum =ii;
      fnum = ASCIIC +fnum;
      optind =0;
      while(optind<argc) {
	ch=getopt(argc,argv,"1234567890");        
	if (ch == fnum){
	  filereq =filereq+1;
	  break;
	}
	
      }
      
    }
  
 return filereq;
} 

 
  

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
author: Taposh Dutta Roy  - dutttap@iit.edu
creation date: PENDING
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

 static void usage(char **argv)
  { 
  fprintf(OUTPTR,"NAME\n\n");
  fprintf(OUTPTR,"\t\t\t %s - 2D Frequency Analysis.\n\n",argv[0]);
  fprintf(OUTPTR,"SYNOPSIS\n\n");
  fprintf(OUTPTR,"\t\t\t %s [OPTION#1]...[FILE]..[OPTION#2]..[OPTION#3]...\n\n",argv[0]);
  fprintf(OUTPTR,"DESCRIPTION\n\n");
  fprintf(OUTPTR,"\t\t\t Does Frequency Analysis of 2D Structures.Should be used with all 3 \n");
  fprintf(OUTPTR,"\t\t\t options  \n");
   fprintf(OUTPTR,"OPTION#1 \n");
  fprintf(OUTPTR,"     -h [FILE], --help[=FILE]\n");
  fprintf(OUTPTR,"   -           produces this message and the list of Options \n\n");
  fprintf(OUTPTR,"         -f FILE, --file=FILE\n");
  fprintf(OUTPTR,"              sets the file to be processed.\n\n");
  fprintf(OUTPTR,"              These options are listed by -h [FILE] or -f FILE -l .\n\n");
  return;
  }





/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Finds out the region,direction of current and gets all 3 currents in a pointer 
author: Taposh Dutta Roy  - dutttap@iit.edu
creation date: Thursday,August  25, 2004
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

static int test_compatability(Toplot *tpl)
{
 
  MYFLOAT *val=NULL;
  Point   *pt=NULL,*pt0=NULL;
  UINT4   numpt=0;
  char    **vname=NULL,*mylvar=NULL,*mycontact=NULL,**vname1=NULL; 
  float   disdr=0,tmp0=0,tmp1=0,tmp2=0,tmp3=0;
  int     nvar=0,ncontact=0,timeflag =0,compatible =0,leng=0 ;
  char    *myreg1=NULL,*mydirect1=NULL,*mydirect2=NULL,*myreg2=NULL,*mytmp1=NULL,*mytmp2=NULL;
  int     ii =0,jj=0;
  
  
  pt0=tpl->pt0;
  numpt=tpl->nptx;  
                                                                           
  nvar = tpl->numvars;
  if( nvar ==0 ){
  IERRS("test_compatability: Number of Variables cannot be Zero.");
  }
  
  else{
  
   
    //  for(vname= tpl->varlopt; vname < (tpl->varlopt + nvar);vname++) {
    for(ii=0 ; ii <  nvar;ii++) {
      leng = strlen(*(tpl->varlopt+ii));
      mylvar = (char*)calloc(leng,sizeof(char));
      mydirect1 = (char*)calloc(leng,sizeof(char));
      strcpy(mylvar,*(tpl->varlopt+ii));
      strcpy(mydirect1,*(tpl->varlopt+ii));
      
      for(jj=(ii+1); jj <nvar;jj++) {
	// if(ii==(nvar-1)){
	//vname1 = (tpl->varlopt+0);
	//}else{
	  vname1 = (tpl->varlopt+jj);
	  //}
	 
	
	  mycontact = (char*)calloc(strlen(*vname1),sizeof(char));
	  strcpy(mycontact,*vname1);
	  mydirect2 = (char*)calloc(strlen(*vname1),sizeof(char));
	  strcpy(mydirect2,*vname1);
		
	if(strcmp(mylvar,mycontact)!=0){
	  myreg1 = (char*)calloc(leng,sizeof(char));
	  strcpy(myreg1,strtok(mydirect1,"_"));
	  strcpy(mydirect1,mylvar);
	  myreg2 = (char*)calloc(strlen(*vname1),sizeof(char));
	  strcpy(myreg2,strtok(mydirect2,"_"));
	  strcpy(mydirect2,*vname1);
	  if (strcmp(myreg1,myreg2)==0){
	    if ((strchr(mydirect1,'x') != NULL ) && (strchr(mydirect2,'x') !=NULL)) {
	      fprintf(OUTPTR,"test_compatability: contact region is %s and direction x \n",myreg2);
	      compatible = compatible +1;
	      break;
	      	
	    }
	     if ((strchr(mydirect1,'y') != NULL ) && (strchr(mydirect2,'y') !=NULL)) {
	       fprintf(OUTPTR,"test_compatability: contact region is %s and direction y \n ",myreg2);
	       break;
	       
	     }
	     if ((strchr(mydirect1,'z') != NULL ) && (strchr(mydirect2,'z') !=NULL)) {
		fprintf(OUTPTR,"test_compatability: contact region is %s and direction z \n",myreg2);
		break;
		
	     }
	  }
	}
      }
  }
  }
  free(mylvar);
  free(vname);
  free(mydirect1);
  free(mydirect2);
  free(myreg2);
  free(myreg1);
  return compatible;
}




/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
author: Taposh Dutta Roy  - dutttap@iit.edu
creation date: Thursday,August  25, 2004
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  
static MYFLOAT** getcurrents(Toplot *tpl,int numcontreg,char* filename )
{
  MYFLOAT *val=NULL;
  Point   *pt=NULL,*pt0=NULL;
  UINT4   numpt=0;
  static FLAG firstime=YES;
  float   **cur,**cur0=NULL,**cur1=NULL;
  int nptx=0,nvar,ii,thiscurrent=0,whichcurrent=0;
  int totcols =0,regno=0 ;
  
  
  pt0=tpl->pt0;
  numpt=tpl->nptx;  

  nvar = tpl->numvars;
  if( nvar ==0 ){
    IERRS("convert_plt: Number of Variables cannot be Zero.\n");
  }
  totcols = ((nvar/2)*CURRENTNO)+1;
  if (totcols < ((numcontreg * CURRENTNO) +1)) {
    WERRS("getcurrents: Raw file has less number of Regions than device file\n");
  }
  if(firstime==YES) {
    firstime==NO;
    cur0=(MYFLOAT**)(calloc(totcols,sizeof(MYFLOAT*)));
    if(cur0==NULL){
      IERRS("getcurrents: Memory Allocation Failed for the pointer \n");
    }
    cur1=cur0+totcols;
    for(cur=cur0;cur<cur1;cur++){
      *cur=(MYFLOAT*)(calloc(numpt,sizeof(MYFLOAT)));
      if(cur==NULL){
	IERRS("getcurrents: Memory Allocation Failed for the pointer\n");
      }
    }
  } else {
    for(cur=cur0;cur<cur1;cur++){
      if(cur==NULL){
	WERRS("getcurrents: Release of Non-Existent Memory\n");
      } else {
	free(cur);
      }
    }
    if(cur0==NULL){
      WERRS("getcurrents: Release of Non-Existent Memory\n");
    } else {
      free(cur0);
    }
    firstime==YES;
    
    return;
  }
  for(pt=pt0; pt<pt0+numpt;pt++) 
    {
      ii=(int)(pt-pt0);
      regno =0;
      
      for(thiscurrent=0;thiscurrent < totcols  ; thiscurrent++){
	whichcurrent = (thiscurrent -( CURRENTNO * regno));
	switch(whichcurrent){
	case( TIME):
	  *((*cur0)+ii)=pt->x;
	  break;
	case(PART):
	  *(*(cur0+thiscurrent)+ii)=((*(pt->value+(2*regno))));
	  break;
	case(DISP): 
	  *(*(cur0+thiscurrent)+ii)=(((*(pt->value+1+(2*regno)))-(*(pt->value+(2*regno)))));	
	  break;
	case(TOTAL):
	  *(*(cur0+thiscurrent)+ii)=(*(pt->value+1+(2*regno)));
	  regno = regno +1;
	  break;
	}
      }
    }  
  // writefile_dat(filename,cur0,numpt,pt0,totcols);
  return(cur0);
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
author: Taposh Dutta Roy  - dutttap@iit.edu
creation date: Thursday,August  28, 2004
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


static MYFLOAT** convert_plt(simpts *sp)
{
  Toplot *tp;
  char *filename;
  int numcontreg;
  
  tp = (Toplot *)calloc(1,sizeof(Toplot));
  tp = sp->tpl;
  
  filename = (char*)calloc(strlen(sp->filename),sizeof(char));
  filename = sp->filename;
  
  numcontreg = sp->numcontreg;
  
  MYFLOAT   **currents=NULL;
  int     compatible;
  compatible = test_compatability(tp );
  if (compatible !=0){
    currents = getcurrents(tp,numcontreg,filename);
    //getcurrents(tpl,numcontreg,filename);
  } else {
    IERRS("convert_plt: Valid parameters not supplied");
  }
  return currents;
  
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
author: Taposh Dutta Roy  - dutttap@iit.edu
creation date: Thursday,August  28, 2004
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

static void convert_dat(Toplot *tpl,char* filename,int numcontreg)
{
  float   **currents=NULL;
  int     compatible;
  compatible = test_compatability(tpl );
  if (compatible !=0){
    currents = getcurrents(tpl,numcontreg,filename);
    getcurrents(tpl,numcontreg,filename);
  } else {
    IERRS("convert_dat: Valid parameters not supplied");
  }
 
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
author: Taposh Dutta Roy  - dutttap@iit.edu
creation date: Thursday,August  28, 2004
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


static void writefile_dat(char* filename,float **cur0,int numpt, Point *pt0,int totcols)
{
  FILE    *outfile=NULL;
  Point   *pt=NULL;
  int     thiscurrent =0,ii=0;
  

  if((outfile=fopen(filename,"w"))==NULL) {
    IERRS("writefile_dat: open file failed.");
  }
  
  for(pt=pt0; pt<pt0+numpt;pt++) 
    {
      ii=(int)(pt-pt0);
      for(thiscurrent=0;thiscurrent < totcols  ; thiscurrent++){
	fprintf(outfile,MYFORMAT,*(*(cur0+thiscurrent)+ii));
	fprintf(outfile,"\t");
	
      }
      fprintf(outfile,"\n");
    }
  if(fclose(outfile)==EOF) {
    IERRS("writefile_dat: close file failed.");
  }
  
}
 
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

author: Taposh Dutta Roy - dutttap@iit.edu
creation date: Saturday, January 20, 2004
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
static void scan_file(char *filename)
{
  UINT1        compression=GZIP;
  char         *title, **vname0, **vname, **vlopt0, **vlopt;
  char         **vsopt, **vsopt0;
  FILE         *file=NULL; 
  UINT1        dimen=0;
  UINT2        namel=0, maxn=0, maxl=0, nal=0;
  UINT4        numpt=0, numvars=0;

  if((file=stream_open(filename,"r",compression))==NULL) {
    fprintf(ERRPTR,"outgra: file \"%s\" not found.\n",optarg);
    EXIT;
  }
  stream_read(file,&dimen,sizeof(UINT1),compression); 
  stream_read(file,&namel,sizeof(UINT2),compression);
  title=(char*)(calloc(namel+1,sizeof(char))); 
  stream_read(file,title,namel,compression);
  stream_read(file,&numpt,sizeof(UINT4),compression);
  stream_read(file,&numpt,sizeof(UINT4),compression);
  stream_read(file,&numpt,sizeof(UINT4),compression);
  stream_read(file,&numvars,sizeof(UINT4),compression);
  fprintf(OUTPTR,"%d-D data set: %s, ",dimen,title);
  fprintf(OUTPTR,"recorded %d quantities:\n\n",numvars); 
  vname=vname0=(char**)(calloc(numvars+1,sizeof(char*)));
  vlopt=vlopt0=(char**)(calloc(numvars+1,sizeof(char*)));
  vsopt=vsopt0=(char**)(calloc(numvars+1,sizeof(char*)));
  *vname=(char*)(calloc(90,sizeof(char))); 
  sprintf(*vname,"DATA");
  *vlopt=(char*)(calloc(90,sizeof(char))); 
  sprintf(*vlopt,"LONG OPTION");
  *vsopt=(char*)(calloc(90,sizeof(char))); 
  sprintf(*vsopt,"SHORT OPTION");
  maxn=strlen(*vname);
  maxl=strlen(*vlopt);
  vlopt++;vsopt++;vname++;
  for(; vsopt<(vsopt0+numvars+1); vsopt++, vlopt++,vname++) {
    stream_read(file,&namel,sizeof(UINT2),compression);
    maxn=(maxn<namel ? namel: maxn);
    *vname=(char*)(calloc(namel+1,sizeof(char))); 
    stream_read(file,*vname,namel,compression);
    stream_read(file,&namel,sizeof(UINT2),compression);
    maxl=(maxl<namel ? namel: maxl);
    *vlopt=(char*)(calloc(namel+1,sizeof(char))); 
    stream_read(file,*vlopt,namel,compression);
    *vsopt=(char*)(calloc(2,sizeof(char))); 
    stream_read(file,*vsopt,sizeof(char),compression);        
  }
  stream_close(file,compression);
  maxn+=2; maxl+=2;
  vsopt=vsopt0;
  vlopt=vlopt0;
  vname=vname0;
  for(; vsopt<(vsopt0+numvars+1); vsopt++, vlopt++,vname++) {
    nal=strlen(*vname);
    *vname=(char*)(realloc(*vname,maxn+1));
    memset(((*vname)+nal),(int)(' '),maxn-nal);
    memset((*vname+maxn),0,1);
    nal=strlen(*vlopt);
    *vlopt=(char*)(realloc(*vlopt,maxl+1));
    memset((*vlopt+nal),(int)(' '),maxl-nal);
    memset((*vlopt+maxl),0,1);
    fprintf(OUTPTR,"%s  %s  %s\n", *vname, *vlopt, *vsopt);  
  }
  fprintf(OUTPTR,"\n\n");
  return;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
author :Taposh Dutta Roy
Creation :Thursday Aug 12 2004
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*
 void check1plot(Toplot  *data,const char *file, UINT1 mode,int numcontreg)                                            
{
  Point  *pt=NULL,*pt0=NULL,*pt1=NULL;
  UINT1    compression=GZIP;
  UINT2    namel=0;
  UINT4    numpt;
  Pline    *pl=NULL, *pl1=NULL;
  FILE     *outfile=NULL;
  char     filename[MAXNAMELEN], **vname=NULL;
  char     *sopt, **lopt, *dotpos=NULL;
  
  if((file==NULL) || strlen(file)==0) {   
    IERRS("outplot: wrong namefile in input.");
  }
  if(data==NULL) {data=convRAW4d(file);}
  switch(data->dimension) {
  case 1: numpt=data->nptx; break;
  default: IERRS("outplot: wrong data dimension."); break;
  }  
  if(numpt<=0) {SERRS("outplot: illegal number of points.");}
  strcpy(filename,file);
  if((dotpos=strstr(filename,".raw\0"))!=NULL) *dotpos=0;
  switch(mode) {
  case GENERIC:
    if((dotpos=strstr(filename,".dat\0"))!=NULL) *dotpos=0;
    strcat(filename,".dat");
    switch(data->dimension) {
    case 1:convert_dat(data,filename,numcontreg); break;
    default: IERRS("outplot: wrong data dimension."); break;
    }
    break;
  case TECPLOT:
    if((dotpos=strstr(filename,".plt\0"))!=NULL) *dotpos=0;
    strcat(filename,".plt");
    switch(data->dimension) {
    case 1: convert_plt(data,filename,numcontreg); break;
    default: IERRS("outplot: wrong data dimension."); break;
    }
    break;
  }
}
 */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Returns a list of selected variables codes to insert in structure.
Parameter numcalsout returns the number of class of output variables.
Parameter frame returns code of included frame.
author: claudio berselli (berselli@neumann.ece.iit.edu)
creation date: Wensday, September 6, 2000
first revision: Monday, January 1, 2001 (marco saraniti)
second revision: Taposh Dutta Roy Current Saturday, Jul24 5:33
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*
static Toplot *findvar(Toplot *tp,     
		  UINT1 *nvout, 
		  int   argc,
		  char  **argv,char *fname,int fnum)           
{
  struct option *opt=NULL, *opt0=NULL, *opt1=NULL;
  Toplot      *tpout=NULL;
  Point       *po=NULL, *pi=NULL, *po1=NULL;
  Pline         *pli=NULL, *plo=NULL, *plo1=NULL;
  MYFLOAT       *valo=NULL, *vali=NULL;
  UINT2         namel=0;
  UINT4         tout=0, numpt=0;
  FLAG          *thisone=NULL, gefunden=NO;
  char          **vout=NULL, **lout=NULL, *sout=NULL, *sopt0=NULL, *sopt=NULL;
  char          **vname=NULL, **lopt=NULL;
  char          **vlopt=NULL, *vsopt=NULL;
  int           lo=0, cout=0, ch=0;

  opt0=opt=(struct option*)(calloc(tp->numvars+1,sizeof(struct option)));
  opt1=opt0+tp->numvars;
  thisone=(FLAG*)(calloc(tp->numvars,sizeof(FLAG)));
  sopt=sopt0=(char*)(calloc(tp->numvars+1,sizeof(char)));
  vlopt=tp->varlopt;
  vsopt=tp->varsopt;
  sopt=sopt0;
  for(opt=opt0; opt<opt1; opt++, sopt++, vlopt++, vsopt++) {
    namel=strlen(*vlopt);
    opt->name=(char*)(calloc(namel+1,sizeof(char)));
    sprintf((char*)(opt->name),"%s",*vlopt);
    *sopt=opt->val=(int)(*vsopt);      
    thisone[opt-opt0]=NO;
  }
  
  fnum = ASCIIC +fnum;

   while(optind<argc) {

     ch=getopt(argc,argv,"1234567890");        
       if (ch == fnum){
	 fprintf(OUTPTR,"\t%d\t%s\n",(fnum-48),fname);
      while(optind<argc) {          
	cout=getopt_long(argc,argv,sopt0,opt0,&lo); 

	  if (cout!='?') {
	    opt=opt0; gefunden=NO;
	    while((opt<opt1)&&(gefunden==NO)) {
	      if((opt->val==cout) && (thisone[opt-opt0]==NO)) {
		thisone[opt-opt0]=YES;
		gefunden=YES;
		tout++;
	      }
	      opt++;
	    }
	  }
	  else {
	    break;
	  }
	      	    
	}
      }
   }
 
 
 
 if(tout !=0) {
*/ 
   /* for(tout=0; tout<tp->numvars; tout++) thisone[tout]=YES;*/
   
/* 
  tpout=(Toplot*)(calloc(1,sizeof(Toplot))); 
  tpout->dimension=tp->dimension;
  tpout->nptx=tp->nptx;
  tpout->npty=tp->npty;
  tpout->nptz=tp->nptz;
  tpout->numvars=tout;
  tpout->title=(char*)(calloc(strlen(tp->title)+1,sizeof(char)));
  sprintf(tpout->title,"%s",tp->title);
  tpout->xlabel=(char*)(calloc(strlen(tp->xlabel)+1,sizeof(char)));
  sprintf(tpout->xlabel,"%s",tp->xlabel);
  tpout->ylabel=(char*)(calloc(strlen(tp->ylabel)+1,sizeof(char)));
  sprintf(tpout->ylabel,"%s",tp->ylabel);
  tpout->zlabel=(char*)(calloc(strlen(tp->zlabel)+1,sizeof(char)));
  sprintf(tpout->zlabel,"%s",tp->zlabel);
  vout=tpout->varname=(char**)(calloc(tpout->numvars,sizeof(char*)));
  lout=tpout->varlopt=(char**)(calloc(tpout->numvars,sizeof(char*)));
  sout=tpout->varsopt=(char*)(calloc(tpout->numvars,sizeof(char)));
  vname=tp->varname; lopt=tp->varlopt; sopt=tp->varsopt;
  for(vname=tp->varname; vname < (tp->varname+tp->numvars); vname++,lopt++,sopt++) { 
    if(thisone[vname-tp->varname]==YES) {
      namel=strlen(*vname);        
      *vout=(char*)(calloc(namel+1,sizeof(char)));
      sprintf(*vout,"%s",*vname);
      namel=strlen(*lopt);
      *lout=(char*)(calloc(namel+1,sizeof(char)));
      sprintf(*lout,"%s",*lopt);
      *sout=*sopt;
#     if(VERBOSE)
      DOVE; fprintf(OUTPTR,"findvar: extracting \"%s\"for \"%s\"\n",*vout,fname);
#     endif
      vout++; lout++; sout++;
    }
  }
  switch(tp->dimension) {
  case 1: 
    if(tp->nptx==0) {IERRS("findvar: nptx=0 in 1-D plot data.");}
    if(tp->npty!=0) {IERRS("findvar: npty!=0 in 1-D plot data.");}
    if(tp->nptz!=0) {IERRS("findvar: nptz!=0 in 1-D plot data.");}
    numpt=tp->nptx; 
    break;
  case 2: 
    if(tp->nptx==0) {IERRS("findvar: nptx=0 in 1-D plot data.");}
    if(tp->npty==0) {IERRS("findvar: npty=0 in 2-D plot data.");}
    if(tp->nptz!=0) {IERRS("findvar: nptz!=0 in 2-D plot data.");}
    numpt=tp->nptx*tp->npty; 
    break;
  case 3: 
    if(tp->nptx==0) {IERRS("findvar: nptx=0 in 3-D plot data.");}
    if(tp->npty==0) {IERRS("findvar: npty=0 in 3-D plot data.");}
    if(tp->nptz==0) {IERRS("findvar: nptz=0 in 3-D plot data.");}
    numpt=tp->nptx*tp->npty*tp->nptz; 
    break;
  default:
    IERRS("findvar: wrong data dimension.");
    break;
  }
  tpout->pt0=po=(Point*)(calloc(numpt,sizeof(Point)));
  po1=tpout->pt0+numpt;
  for(pi=tp->pt0, po=tpout->pt0; po<po1; po++, pi++) {
    po->x=pi->x;
    po->y=pi->y;
    po->z=pi->z;
    valo=po->value=(MYFLOAT*)(calloc(tpout->numvars,sizeof(MYFLOAT)));
    for(vali=pi->value; vali<(pi->value+tp->numvars); vali++) {
      if(thisone[vali-pi->value]==YES) {
        if(valo==(po->value+tpout->numvars)) {
	  IERRS("findvar: pointer out of range.");
	}
	*valo++=*vali;
      }
    }
    if(valo!=(po->value+tpout->numvars)) {
      IERRS("findvar: pointer out of range.");
    }
  }
  tpout->numpl=tp->numpl;
  if(tpout->numpl) {
    tpout->pl0=(Pline*)(calloc(tpout->numpl,sizeof(Pline))); 
    plo1=tpout->pl0+tpout->numpl;
    for(pli=tp->pl0, plo=tpout->pl0; plo<plo1; plo++, pli++) { 
      memcpy((Pline*)(plo),(Pline*)(pli),sizeof(Pline));
      po=plo->pt0=(Point*)(calloc(plo->numpt,sizeof(Point)));
      memcpy((Point*)(po),(Point*)(pli->pt0),sizeof(Point)*plo->numpt);
    }
  }
  for(opt=opt0; opt<opt1; free((char*)(opt->name)), opt++);
  FREE(sopt0);
  FREE(thisone);
  relax_plot(tp);
  return(tpout);
  }
  return(0);
}
*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
author: Taposh Dutta Roy  - dutttap@iit.edu
creation date: Thursday, Sept 9,2004
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*
static MYFLOAT**  allocatemem(simpts *dataout,int totfilereq)
{
  int nvar =0,numpt=0,ii=0;
  int numcontreg =0;
  MYFLOAT   **cur,**cur0=NULL,**cur1=NULL;
    
  for (ii=0;ii<(totfilereq);ii++)
    {
      numcontreg = ((dataout+ii)->numcontreg);
      nvar =((dataout+ii)->tpl->numvars);
      numpt =((dataout+ii)->tpl->nptx);
      
      if( nvar ==0 ){
	IERRS("convert_plt: Number of Variables cannot be Zero.\n");
      }
      if ((((nvar/2)*CURRENTNO)+1) < ((numcontreg * CURRENTNO) +1)) {
	WERRS("getcurrents: Raw file has less number of Regions than device file\n");
      }
      cur0=(MYFLOAT**)(calloc(((nvar/2)*CURRENTNO)+1,sizeof(MYFLOAT*)));
      if(cur0==NULL){
	IERRS("getcurrents: Memory Allocation Failed for the pointer \n");
      }
      cur1=cur0+(((nvar/2)*CURRENTNO)+1);
      for(cur=cur0;cur<cur1;cur++){
	*cur=(MYFLOAT*)(calloc(numpt,sizeof(MYFLOAT)));
	if(cur==NULL){
	  IERRS("getcurrents: Memory Allocation Failed for the pointer\n");
	}
      }
    }
}
*/
