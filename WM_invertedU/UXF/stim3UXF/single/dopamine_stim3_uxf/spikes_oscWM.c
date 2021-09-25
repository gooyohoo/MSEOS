#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <signal.h>

#include "randdev.h" 
#define ROUND(x) ((int)(x+0.5))

/* ------------------------------------------------------------------------*/
/* ------------------------------- PARAMETRI ----------------------------- */
/* ------------------------------------------------------------------------*/
#define Dopamine 2


#define N 10000 
#define Pe 0.8
#define STIM_DELAY_E 0.0
#define STIM_DELAY_I 0.0
int PE=(int)(Pe*100);  
int D1=(int)(Dopamine*100);
float LIFE=0,De,Di;

float pop_part=1.0;
int pops=5; 
unsigned int Ne=N*Pe, Ni=N-N*Pe;
float pconn=0.20, Cp=.10;

float f=.10f, f_noise=.15f, x=0.00f;  

float Jdep , Jpot, Jie , Jei , Jii;
/*
float Jdep=0.10 , Jpot=0.45;
float Jie=-.25f , Jei=.135 , Jii=-.20f; */

float TauNMDA=100.0, xNMDA_EE=0.0, xNMDA_IE=0.0;
float Taue=15, Taui=10;	      
float Tarpge=2, Tarpgi=2;     
float ThrE=20, ThrI=20; 
float He=.8, Hi=.65;	       

float width_stim=0.25;
float width_popout=0.00;

//Initialization Parameters
float nue=.001, nui=.005;
float mu_e_mf, s2_e_mf, mu_i_mf, s2_i_mf;
float avg_x, avg_x2, avg_u;
float muV_e, muV_i, s2V_e, s2V_i;

float contrastE_sel=1.1f, contrastE_nonsel=1.0f;
float contrastI_sel=1.0, contrastI_nonsel=1.00;
float contrastE_sel_noise=1.0f, contrastE_nonsel_noise=1.0f;
float contrastI_sel_noise=1.0, contrastI_nonsel_noise=1.00;

float muEext=23.80, muIext=21.00, sEext=1.00, sIext=1.00;

float Tprestim, Tcue, Tdelay, Ttest, Tonset, Toffset;	 
float Tprestim_SA=15000, Tcue_SA=0, Tdelay_SA=0, Ttest_SA=0;	 
float Tprestim_PT=1000, Tcue_PT=200, Tdelay_PT=2800;	 
float T_popout=2000, DT_popout=0, g_popout=1.00;//1.05 

float dt = 0.1f;     /* step minimo in ms                      */
float TCamp = 10.0;      /* Tempo[ms] di campionamento attivita' neurale*/

float MinDelay  = 0.1f; /* Ritardo sinaptico minimo in ms (per eccitatori) */
float StepDelay = 0.1f; /* Larghezza step delay  */
int MaxStep = 10.0;

int num_pres;

/* Parametri per la dinamica sinaptica */
float tau_D=200.0, tau_F=1500.0, U=0.2;
//float tau_D=700.0, tau_F=20.0, U=0.3;
//float tau_D=350.0, tau_F=350.0, U=0.3;
float qp=1.0, qm=.1;

/**** FLAGS ****/
int flagstruct=1, flag_fac=1, flag_SA, flag_popout, flag_fs;

/* ------------------------------------------------------------------------*/
/* ---------------------- INIZIO DEFINIZIONE DATI ------------------------ */
/* ------------------------------------------------------------------------*/

unsigned int Ntot;            /* Neuroni totali */

float *Inmda, *Iext;       

char *gg, **P, **D, **sel_old;

int **nofs_conn, **nofs_conn_inv, *nconn;
int *delay_stim, *duration_stim, *duration_popout;
int *spikes, spikesI;
float *curr, currI;

float theta; 

typedef struct _neuron Neuron;
typedef struct _population Population;
typedef struct _efficacy Efficacy; 
typedef struct _arrayroot EffRoot; 

struct _neuron
{
  float V;
  int status;
  int tspike;
  float Hr;
  float xtm;
  float utm;
  float uxtm;
};
  
struct _population
{
  unsigned int num;
  int *index;
};
  
struct _efficacy
{
  int dist; /* distanza del prossimo neurone */ 
  float efficacy; /* efficacia sinaptica */
  char delay; 
};

struct _arrayroot
{
  unsigned int num;
  Efficacy *eff;
};


Population *Pop;
Neuron *neur;
EffRoot *J;         
float **Events;     
float **EventsNMDA;

FILE *rates, *curr_out, *ratepops, *spikeout, *synvar_u, *synvar_x,*num_para;


unsigned int Time;

/*************************************************************************/

void shifttspike()
{
  int life;
  int i;
  life= (int)ROUND((Tprestim+Tcue+Tdelay+Ttest)/dt);
  for (i=0; i<Ntot; i++)
    neur[i].tspike-=life;
}

void STP(int i, int t1)
{
  int t0;
  float x0, x1, delt, u0, u1;
  x0=neur[i].xtm;
  u0=neur[i].utm;
  t0=neur[i].tspike;
  delt=(t1-t0)*dt;

  if (flag_fac==0)
    u0=U;
  x1=1.0-(1.0-(1.0-u0)*x0)*exp(-delt/tau_D);
  u0=U+u0*(1-U);
  u1=U-(U-u0)*exp(-delt/tau_F);
  neur[i].uxtm=u1*x1;
  neur[i].xtm=x1;
  if (flag_fac==0)
    u1=U;
  neur[i].utm=u1;
}

float readout(int i,int t1,char flag_ux)
{
  int t0;
  float x0,x1,delt, u0, u1;
  float ris;
  x0=neur[i].xtm;
  u0=neur[i].utm;
  t0=neur[i].tspike;
  delt=(t1-t0)*dt;
  if (flag_fac==0)
    u0=U;
  if (flag_ux==0)
    ris=1.0-(1.0-x0)*exp(-delt/tau_D);
  else
    ris=U-(U-u0)*exp(-delt/tau_F);
  return(ris);
}

/** Populations' routines - BEGIN **/

void makePop(char **sel)
{
  int i, p, count;
  
  Pop=(Population*)malloc(pops*sizeof(Population));
  for (p=0; p<pops; p++)
    {
      
      count=0;
      for (i=0; i<Ne; i++)
	count+=sel[i][p];
      Pop[p].num=count;
      Pop[p].index=(int*)malloc(count*sizeof(int));
      
      count=0;
      for (i=0; i<Ne; i++)
	{
	  if (sel[i][p]==1)
	    {
	      Pop[p].index[count]=i;
	      count++;
	    }
	}
    }
}

char** createPop(char flag_fix)
{
  char **sel;
  int i,p, mul, count, Nfix_E=ROUND(f*Ne), Nfix_I=ROUND(f*Ni);
  
  // Allocation and Initialization of matrix selectivity

  sel=(char **)malloc((Ne+Ni)*sizeof(char *));
  for (i=0; i<(Ne+Ni); i++)
    {
      sel[i]=(char *)malloc((pops+1)*sizeof(char));
      for (p=0; p<=pops; p++)
	sel[i][p]=0;
    }

  // Generation of Population

  for (p=0; p<pops; p++)
    {
      if (flag_fix==0)
	{
	  for (i=0; i<(Ne+Ni); i++)
	    if (Random()<f)
	      sel[i][p]=1;
	}
      else
	{
	  count=0;
	  while (count<Nfix_E)
	    {
	      i=((int)(Random()*Ne));
	      if (sel[i][p]==0)
		{
		  sel[i][p]=1;
		  count++;
		}
	    }
	  count=0;
	  while (count<Nfix_I)
	    {
	      i=Ne+((int)(Random()*Ni));
	      if (sel[i][p]==0)
		{
		  sel[i][p]=1;
		  count++;
		}
	    }
	}
    }
  
  // Generation of Multiplicity

  for (i=0; i<(Ne+Ni); i++)
    {
      mul=0;
      for (p=0; p<pops; p++)
	mul+=sel[i][p];
      sel[i][pops]=mul;
    }
  
  return(sel);

}

char** createPop_noov()
{
  char **sel;
  int i,p, mul, count, Nfix_E=ROUND(f*Ne), Nfix_I=ROUND(f*Ni);
  int N1, N2;
  // Allocation and Initialization of matrix selectivity

  sel=(char **)malloc((Ne+Ni)*sizeof(char *));
  for (i=0; i<(Ne+Ni); i++)
    {
      sel[i]=(char *)malloc((pops+1)*sizeof(char));
      for (p=0; p<=pops; p++)
	sel[i][p]=0;
    }

  // Generation of Population
  for (p=0; p<pops; p++)
    {
      N1=p*Nfix_E; N2=(p+1)*Nfix_E;
      for (i=N1; i<N2; i++)
	sel[i][p]=1;
      
      N1=p*Nfix_I; N2=(p+1)*Nfix_I;
      for (i=Ne+N1; i<Ne+N2; i++)
	sel[i][p]=1;
    }
  
  
  // Generation of Multiplicity

  for (i=0; i<(Ne+Ni); i++)
    {
      mul=0;
      for (p=0; p<pops; p++)
	mul+=sel[i][p];
      sel[i][pops]=mul;
    }
  
  return(sel);

}


void computePD(char **sel)
{
  int p, i, j;
  int running_P, running_D;

  //printf("> computePD()\n");

  P=(char **)malloc(Ne*sizeof(char *));
  D=(char **)malloc(Ne*sizeof(char *));

  for (i=0; i<Ne; i++)
    {
      P[i]=(char *)malloc(Ne*sizeof(char));
      D[i]=(char *)malloc(Ne*sizeof(char));
    }

  for (i=0; i<Ne; i++)
    {
      for (j=0; j<Ne; j++)
	{
	  running_P=running_D=0;
	  for (p=0; p<pops; p++)
	    {
	      running_P+=((sel[i][p]==1)&&(sel[j][p]==1));
	      running_D+=((sel[i][p]==1)&&(sel[j][p]==0));
	    }
	  P[i][j]=running_P;
	  D[i][j]=running_D;
	}
    }


  //printf("> end computePD()\n");
}

/** Populations' routines - END **/



/** Synapses' routines - BEGIN **/

void initSynapses_nofs()
{
  //June 23 2007
  //it works only for non-overlapping memories

  char theend, theend_mem;
  char *excit, *inhib;
  int mem, ne, ni, Nmem, Nrest;
  int C_mem, C_rest, C_totE, C_totI;
  int post, pre, count_mem, count_others, count_I, count;
  int k, run_pre;

  //nofs_conn[N][cN] - to be declared outside
  
  nofs_conn=(int **)calloc((Ne+Ni),sizeof(int *));
  for (k=0; k<(Ne+Ni); k++)
    nofs_conn[k]=(int *)calloc(ROUND(pconn*(Ne+Ni)),sizeof(int));
  
  //excitatory cells
  Nmem=ROUND(f*Ne); Nrest=Ne-(pops*Nmem);
  C_totE=ROUND(pconn*Ne);  C_totI=ROUND(pconn*Ni);
  C_mem=ROUND(pconn*f*Ne);  C_rest=C_totE-C_mem;
  
  excit=(char *)calloc(Ne,sizeof(char));
  inhib=(char *)calloc(Ni,sizeof(char));

  for (post=0; post<(Ne+Ni); post++)
    {
      for (k=0; k<Ne; k++)
	excit[k]=0;
      for (k=0; k<Ni; k++)
	inhib[k]=0;
      
      theend=0; count=0; mem=0;  
      
      while(theend==0)
	{
	  theend_mem=0; count_mem=0; 
	  while (theend_mem==0)
	    {
	      pre=((int)(Random()*Nmem));
	      run_pre=mem*Nmem+pre;
	      
	      if (excit[run_pre]==0)
		{
		  excit[run_pre]=1;
		  nofs_conn[post][count]=run_pre;
		  count++; count_mem++;
		  
		  if (count_mem==C_mem)
		    {
		      theend_mem=1;
		      mem++;
		    }
		}
	      
	      if (mem==pops)
		theend=1;
	    }
	}
	  
      theend=0;
      while(theend==0)
	{
	  pre=((int)(Random()*Nrest));
	  run_pre=mem*Nmem+pre;
	  if (excit[run_pre]==0)
	    {
	      excit[run_pre]=1;
	      nofs_conn[post][count]=run_pre;
	      count++;
	      if (count==(C_mem+C_rest))
		theend=1;
	    }
	}
      
      //select inhibitory
      theend=0;
      while(theend==0)
	{
	  pre=((int)(Random()*Ni));
	  run_pre=Ne+pre;
	  if (inhib[pre]==0)
	    {
	      inhib[pre]=1;
	      nofs_conn[post][count]=run_pre;
	      count++;
	      if (count==(C_mem+C_rest+C_totI))
		theend=1;
	    }
	} 
    }
  free(excit); free(inhib);
}

void invert_nofs_udo()
{
  int post, k, max_conn, run_i;
  int C_tot=ROUND(pconn*(Ne+Ni));
  int *count;

  nconn=(int *)calloc((Ne+Ni),sizeof(int));
  
  for (post=0; post<(Ne+Ni); post++)
    for (k=0; k<C_tot; k++)
      nconn[nofs_conn[post][k]]++;
    
  max_conn=0;
  for (k=0; k<(Ne+Ni); k++)
    if (max_conn<nconn[k])
      max_conn=nconn[k];
  
  count=(int *)calloc((Ne+Ni),sizeof(int));

  nofs_conn_inv=(int **)calloc((Ne+Ni),sizeof(int*));
  for (k=0; k<(Ne+Ni); k++)
    nofs_conn_inv[k]=(int *)calloc(nconn[k],sizeof(int));
  
  for (post=0; post<(Ne+Ni); post++)
    {
      for (k=0; k<C_tot; k++)
	{
	  run_i=nofs_conn[post][k];
	  nofs_conn_inv[run_i][count[run_i]]=post;
	  count[run_i]++;
	}
    }

  free(count); 

  for (post=0; post<(Ne+Ni); post++)
    free(nofs_conn[post]);
  free(nofs_conn);

} 


void InitSynapses()
{
  
  int i1, i2, i3, kk, count, oldi2, dummy;
  int* numarray=NULL;
  float ppot;
  float Jr;
  FILE *fpo_chk;
  
  Ntot=Ne+Ni;
  /* alloca la matrice sinaptica */
  //printf("Drawing synaptic matrix: please, wait...\n");
  J=(EffRoot*)malloc(Ntot*sizeof(EffRoot));
  
  if (flagstruct==0)
    printf(">> no structuring!!\n");
  else
    computePD(sel_old);
  
  fpo_chk=fopen("chk_syndist.log","w");

  for(i1=0; i1<Ntot; i1++)
    { 
      if(!numarray)
	numarray=(int*)malloc(sizeof(int)*((int)(5*pconn*Ntot)));
      
      if (flag_fs==0)
	{
	  count=oldi2=0;
	  for(i2=0; i2<Ntot; i2++)
	    {
	      if(Random()<pconn)
		{
		  numarray[count]=i2-oldi2 ;
		  if( ((count+1)>(5*pconn*Ntot-1))) 
		    break;
		  count++;
		  oldi2=i2;
		}
	    } 
	}
      else
	{
	  numarray[0]=nofs_conn_inv[i1][0];
	  oldi2=0;
	  for (i2=0; i2<nconn[i1]; i2++)
	    {
	      numarray[i2]=nofs_conn_inv[i1][i2]-oldi2;
	      oldi2=nofs_conn_inv[i1][i2];
	    }
	  count=nconn[i1];
	}
	 
      fprintf(fpo_chk,"%d\n",count);
      
      J[i1].num=count;
      J[i1].eff=(Efficacy*)malloc(sizeof(Efficacy)*count);
      i3=0;
      for(i2=0; i2<count; i2++)
	{
	  J[i1].eff[i2].dist=numarray[i2];
	  i3+=numarray[i2];
	  if (i1<Ne)
	    (i3<Ne) ? (Jr=Jdep) : (Jr=Jei);
	  else
	    (i3<Ne) ? (Jr=Jie) : (Jr=Jii);
	  J[i1].eff[i2].efficacy=Jr;
	  if ((i1<Ne)&&(i3<Ne))  
	    {
	
	      if (flagstruct==1)
		{
		  if ((P[i1][i3]==0)&&(D[i1][i3]==0))
		    ppot=Cp;
		  else
		    ppot=qp*P[i1][i3]/(qp*P[i1][i3]+qm*D[i1][i3]);
		}
	      else
		ppot=Cp;
	      if (Random()<ppot){
		J[i1].eff[i2].efficacy=Jpot;      

	      }
	    }
	  
	  if (i1>Ne)
	    J[i1].eff[i2].delay=(char)(0);
	  else
	    {
	      kk=(int)(Random()*MaxStep);
	      J[i1].eff[i2].delay=(char)(kk); 
	    } 
	}
      free(numarray);
      numarray=NULL;      
    }
  
  fclose(fpo_chk);

  for (i1=0; i1<Ne; i1++)
    {
      free(P[i1]);
      free(D[i1]);
    }
  free(P); free(D);
  
  if (flag_fs==1)
    {
      for (i1=0; i1<(Ne+Ni); i1++)
	free(nofs_conn_inv[i1]);
      free(nofs_conn_inv);
    }
}


/** Initialization routines - BEGIN **/

void InitEventMatrix()
{
  int i, t, tbound, deltamin;
  float run_mu;

  tbound=(int)((MinDelay + MaxStep*StepDelay)/dt+0.5)+1;
  deltamin=(int)((MinDelay)/dt+0.5)+1;
  Events=(float **)malloc(tbound*sizeof(*Events));
  EventsNMDA=(float **)malloc(tbound*sizeof(*EventsNMDA));
  for(t=0; t<tbound; t++)
    {
      Events[t]=(float *)malloc((Ne+Ni)*sizeof(**Events));
      EventsNMDA[t]=(float *)malloc((Ne+Ni)*sizeof(**EventsNMDA));
      for(i=0; i<(Ne+Ni); i++)
	{
	  Events[t][i]=0.0;
	  EventsNMDA[t][i]=0.0; 
	}	  
    }
}

void InitV()
{
  int count;
  int tauarp;
  float run_V, w;
  float gu=1.0;

  Ntot=Ne+Ni;

  neur=(Neuron *)malloc(Ntot*sizeof(Neuron));

  gg=(char *)malloc(Ntot*sizeof(char));
  Inmda=(float*)malloc(Ntot*sizeof(float));
  Iext=(float*)malloc(Ntot*sizeof(float));
  delay_stim = (int*) malloc(Ntot*sizeof(int));
  duration_stim = (int*) malloc(Ntot*sizeof(int));
  duration_popout = (int*) malloc(Ntot*sizeof(int));
  spikes=(int*) malloc( Ne * sizeof(int));
  curr=(float*) malloc( Ne * sizeof(float));
  

  for(count=0; count<Ntot; count++)
    {
      if (count<Ne)
	{
	  //neur[count].Hr=.5+.4*Random();
	  neur[count].Hr=He;
	  if (Random()<nue*Tarpge)
	    {
	      tauarp=ROUND(Tarpge/dt); 
	      neur[count].status=-(int)(tauarp*Random());
	      neur[count].V=0;
	      neur[count].tspike=0.0;
	      neur[count].xtm=0.1;
	      neur[count].utm=U*gu;
	    }
	  else
	    {
	      run_V=muV_e+sqrt(s2V_e)*NormDev();
	      if (run_V>=ThrE)
		neur[count].V=muV_e;
	      else
		neur[count].V=run_V;
	      neur[count].status=0;
	      neur[count].xtm=0.1;
	      neur[count].utm=U*gu;
	      w=1200.0;
	      neur[count].tspike=-((int)(w/dt));
	    }
	}
      else
	{
	  neur[count].Hr=Hi;
	  if (Random()<nui*Tarpgi)
	    {
	      tauarp=ROUND(Tarpgi/dt); 
	      neur[count].status=-(int)(tauarp*Random());
	      neur[count].V=0;
	      neur[count].tspike=0.0;
	    }
	  else
	    {
	      run_V=muV_i+sqrt(s2V_i)*NormDev();
	      if (run_V>=ThrI)
		neur[count].V=muV_i;
	      else
		neur[count].V=run_V;
	      neur[count].status=0;
	    }
	  
	  neur[count].xtm=0;
	  neur[count].utm=0;
	  neur[count].tspike=0.0;
	}
      
      Inmda[count]=0.0;
      Iext[count]=0.0;
      delay_stim[count]=0;
      gg[count]=0;
      
      if (count<Ne)
	spikes[count]=0;
    }
}

/** Initialization routines - END **/

// simple pseudo-code for running a single trial

void prepare_popout()
{
  int i;
  for (i=0; i<Ntot; i++)
    duration_popout[i]=ROUND((DT_popout*(1+width_popout*(.5-Random())))/dt);
}

void select_VR(int cue)
{
  int i;
  for (i=0; i<Ntot; i++)
    {
      gg[i]=0;
      delay_stim[i]=0;
      duration_stim[i]=0;
    }
  
  if (cue<pops)
    {
      //select one of the stored items

      /* select excitatory visual responsive */
      
      for (i=0; i<Ne; i++)
	{
	  if ((sel_old[i][cue]==1)&&(Random()<(1-x*(1-f))))
	    {
	      gg[i]=1;
	      delay_stim[i]=-ROUND(Random()*STIM_DELAY_E/dt);
	      duration_stim[i]=ROUND((Tcue_PT*(1+width_stim*(.5-Random())))/dt);
	    }
	  if ((sel_old[i][cue]==0)&&(Random()<f*x))
	    {
	      gg[i]=1;
	      delay_stim[i]=-ROUND(Random()*STIM_DELAY_E/dt);
	      duration_stim[i]=ROUND((Tcue_PT*(1+width_stim*(.5-Random())))/dt);
	    }
	}
      
      /* select inhibitory visual responsive */
      
      for (i=Ne; i<Ne+Ni; i++)
	{
	  if ((sel_old[i][cue]==1)&&(Random()<(1-x*(1-f))))
	    {
	      gg[i]=1;
	      delay_stim[i]=-ROUND(Random()*STIM_DELAY_I/dt);
	      duration_stim[i]=ROUND((Tcue_PT*(1+width_stim*(.5-Random())))/dt);
	    }
	  if ((sel_old[i][cue]==0)&&(Random()<f*x))
	    {
	      gg[i]=1;
	      delay_stim[i]=-ROUND(Random()*STIM_DELAY_I/dt);
	      duration_stim[i]=ROUND((Tcue_PT*(1+width_stim*(.5-Random())))/dt);
	    }
	}
    }
  else
    {
      //generate a noise pattern
      for (i=0; i<Ne; i++)
	{
	  if (Random()<f_noise)
	    {
	      duration_stim[i]=ROUND((Tcue_PT*(1+width_stim*(.5-Random())))/dt);
	      delay_stim[i]=-ROUND(Random()*STIM_DELAY_E/dt);
	      gg[i]=1;  
	    }
	}
    }
}


void single_trial(int cue)
{
  int i, i2, n;
  
  char flag_swap;
  float mext, vext, tau, h;
  float cE_sel, cE_nsel, cI_sel, cI_nsel;
  int tauarp;	
  
  int p;

  unsigned int life;
  unsigned int totsp_cue, totsp_test, totsp_others, totsp_bk, totsp_total;
  
  int popstim, flagstim;
  int tbound, deltatmin, deltatstep, deltat;
  int Ncue, Ntest, Nn, Ntotal;

  unsigned int tcamp, tsyn, tcampV;
  int stimrate[pops], i_sel, i_nonsel;
  
  char buf[30], flag_i_sel=0, flag_i_nonsel=0;
  float nuI_s, nuI_n;
  int count_swap;

  float tot_currE;
  float u_s, x_s, ux_s, u_n, x_n, ux_n, run_u, run_x;
  
  if (flag_SA==1)
    {
      Tprestim=Tprestim_SA;
      Tcue=Tcue_SA;
      Tdelay=Tdelay_SA;
    }
  else
    {
      Tprestim=Tprestim_PT;
      Tcue=Tcue_PT;
      Tdelay=Tdelay_PT;
    }

  Time=0;
  life=(unsigned int)ROUND((Tprestim+(num_pres+1)*Tcue+Tdelay+Ttest)/dt);
  tbound=ROUND((MinDelay + MaxStep*StepDelay)/dt)+1;
  deltatmin=ROUND((MinDelay)/dt);
  deltatstep=ROUND((StepDelay)/dt);
  
  tcamp=(unsigned int)ROUND(TCamp/dt);
  
  (void)select_VR(cue);
  Tonset=Tprestim; Toffset=Tprestim+Tcue;
  (void)prepare_popout();
  if (cue<pops)
    {
      cE_sel=contrastE_sel; cE_nsel=contrastE_nonsel;
      cI_sel=contrastI_sel; cI_nsel=contrastI_nonsel;
    }
  else
    {
      cE_sel=contrastE_sel_noise; cE_nsel=contrastE_nonsel_noise;
      cI_sel=contrastI_sel_noise; cI_nsel=contrastI_nonsel_noise;
    }

  currI=0; spikesI=0;
  while(Time<life)
    { 
      tauarp=ROUND(Tarpge/dt); theta=ThrE;
      tau=Taue;
      
      /* sampling rate */
      
      if (Time%tcamp==tcamp-1)
	{
	  fprintf(ratepops,"%.2f ",(Time+1)*dt);
	  fprintf(curr_out,"%.2f ",(Time+1)*dt);
	  totsp_cue=0; totsp_others=0; totsp_total=0; 
	  for (p=0; p<pops; p++)
	    {
	      totsp_cue=0; tot_currE=0; Nn=0;
	      for (i=0; i<Ne; i++)
		if (sel_old[i][p]==1)
		  {
		    totsp_cue+=spikes[i];
		    tot_currE+=curr[i];
		    Nn++;
		  }
	      fprintf(ratepops,"%.2f ",(float)totsp_cue/TCamp/Nn*1000);
	      fprintf(curr_out,"%.4f ",tot_currE/TCamp/Nn);
	    }
	  
	  fprintf(ratepops,"%.2f\n",(float)spikesI/TCamp/Ni*1000);
	  fprintf(curr_out,"%.4f\n",currI/TCamp/Ni);
	  
	  fflush(ratepops); fflush(curr_out);
	  
	  spikesI=0; currI=0;
	  for (i=0; i<Ne; i++)
	    { spikes[i]=0; curr[i]=0; }
	  
	  /* end sampling rate */
	}
    
      /* Gaussian Currents  - Prepare external currents */
      
      flag_popout=0;
      if ((flag_SA==0)&&((Time*dt)>=T_popout))
	flag_popout=1;
      
      for (i=0; i<Ntot; i++)
	{
	  mext=muEext*dt/Taue; vext=sEext*dt/Taue;
	  if (i>Ne)
	    { mext=muIext*dt/Taui; vext=sIext*dt/Taui; }
	  
	  
	  if ((flag_popout==1)&&(i<Ne))
	    {
	      duration_popout[i]--;
	      if (duration_popout[i]>=0)
	      {if((gg[i]==1)&&(Random()> pop_part))
	          mext*=1;
   		   else
      	      mext*=g_popout;
		  }
		
	      else
		duration_popout[i]=-1;
	    }
	  
	  flagstim=0;
	  if (((Time*dt)>=Tonset)&&((Time*dt)<Toffset)&&(flag_SA==0))
	    flagstim=1;
	  
	  if ((flagstim==1))
	    {
	      delay_stim[i]++;
	      if (delay_stim[i]>=0)
		{
		  delay_stim[i]=1;
		  duration_stim[i]--;
		  if (duration_stim[i]>0)
		    {
		      if (i<Ne)
			{
			  if (gg[i]==1)
			    { mext*=cE_sel; vext*=cE_sel; }
			  else
			    { mext*=cE_nsel; vext*=cE_nsel; }
			}
		      else
			{
			  if (gg[i]==1)
			    { mext*=cI_sel; vext*=cI_sel; }
			  else
			    { mext*=cI_nsel; vext*=cI_nsel; }
			}
		    }
		  else
		    duration_stim[i]=-1;
		}
	    }
	  
	  Iext[i]=mext+NormDev()*sqrt(vext);
	} 
      
      /* Cycle over neurons - i = pre-synaptic */
      
      for (i=0; i<Ntot; i++)   
	{
	  tau=Taue;
	  if (i>Ne)
	    {
	      tauarp=ROUND(Tarpgi/dt); 
	      theta=ThrI;
	      tau=Taui;
	    }
	  
	  /* manage NMDA */
	  
	  Inmda[i]*=exp(-dt/TauNMDA);
	  Inmda[i]+=EventsNMDA[Time%tbound][i];
	  EventsNMDA[Time%tbound][i]=0;
	  
	  if (i<Ne)
	    curr[i]+=Events[Time%tbound][i];
	  else
	    currI+=Events[Time%tbound][i];
	  
	  if(neur[i].status>=0) /* ci troviamo in fase refrattaria? */
	    {
	      neur[i].status++;
	      neur[i].V*=exp(-dt/tau);
	      neur[i].V+=Inmda[i]*dt/tau;
	      neur[i].V+=Iext[i];
	      neur[i].V+=Events[Time%tbound][i]; /* incoming spikes */
	      Events[Time%tbound][i]=0;
	      
	      if(neur[i].V>=theta)
		{  
		  neur[i].V=-1.0;
		  neur[i].status=-tauarp;

		  if (i<Ne)
		    (void)STP(i,Time); 
		  
		  //write out spike rasters
		  //if ((i%10==0))
		   // {
		      fprintf(spikeout,"%d\t%.2f\t",i,(Time*dt));
		      if (i<Ne)
			{
			  fprintf(spikeout,"%d\t",sel_old[i][pops]);
			  if (gg[i]==1)
			    fprintf(spikeout,"1\n");
			  else
			    fprintf(spikeout,"0\n");
			}
		      else
			fprintf(spikeout,"-1\t-1\n");
		      
		      fflush(spikeout);
		  //  }
		  
		  /* propagate spike */
		  i2=0;
		  for(n=0;n<J[i].num;n++)
		    {
		      i2+=J[i].eff[n].dist;
		      deltat=(int)(deltatmin+deltatstep*(J[i].eff[n].delay));
		      if (i<Ne)
			{
			  if(i2<Ne)
			    {
			      Events[(Time+deltat)%tbound][i2]+=((1-xNMDA_EE)*neur[i].uxtm*J[i].eff[n].efficacy);
			      EventsNMDA[(Time+deltat)%tbound][i2]+=(Taue*xNMDA_EE*neur[i].uxtm*J[i].eff[n].efficacy/TauNMDA);				
			    }
			  else
			    {
			      Events[(Time+deltat)%tbound][i2]+=(1-xNMDA_IE)*J[i].eff[n].efficacy;
			      EventsNMDA[(Time+deltat)%tbound][i2]+=Taui*xNMDA_IE*(J[i].eff[n].efficacy)/TauNMDA;
			    } 
			}
		      else
			Events[(Time+deltat)%tbound][i2]+=J[i].eff[n].efficacy;
		    }
		  //update tspike 
		  neur[i].tspike=Time;
		  
		  if (i<Ne)
		    spikes[i]++;
		  else
		    spikesI++;
		}
	    }
	  else
	    {
	      /* fase refrattaria */
	      neur[i].status++;
	      Events[Time%tbound][i]=0;
	      
	      if(neur[i].status==0)  /* Finito il tarp? */
		neur[i].V=(float)(theta*neur[i].Hr); 
	    }
	}
      
      if (Time%tcamp==tcamp-1)
	{
	  /* start sampling synaptic variables */
	  fprintf(synvar_u,"%.2f\t",(Time+1)*dt);
	  fprintf(synvar_x,"%.2f\t",(Time+1)*dt);
	  for (p=0; p<pops; p++)
	    {
	      u_s=0; x_s=0; Nn=0;
	      for (i=0; i<Ne; i++)
		if (sel_old[i][p]==1)
		  {
		    u_s+=readout(i,Time,1);
		    x_s+=readout(i,Time,0);
		    Nn++;
		  }
	      fprintf(synvar_u,"%.5f\t",u_s/Nn);
	      fprintf(synvar_x,"%.5f\t",x_s/Nn);
	    }
	  fprintf(synvar_u,"\n"); fprintf(synvar_x,"\n");
	  fflush(synvar_u); fflush(synvar_x);
	}
      Time++;  /* e via si riparte */
    }
}

void naive_mf()
{
  float Ce, Ci, Jee;
  
  avg_u=U*(1+tau_F*nue)/(1+U*tau_F*nue);
  avg_x=1/(1+(1-avg_u)*tau_D*nue);
  avg_x2=avg_x/(1+avg_u*tau_D*nue*(1-.5*avg_u));
  
  Ce=pconn*Ne; Ci=pconn*Ni;
  Jee=Cp*Jpot+(1-Cp)*Jdep;
  
  mu_e_mf=Ce*Jee*avg_u*avg_x*nue+Ci*Jie*nui;
  mu_i_mf=Ce*Jei*nue+Ci*Jii*nui;
  s2_e_mf=Ce*Jee*Jee*avg_u*avg_u*avg_x2*nue+Ci*Jie*Jie*nui;
  s2_i_mf=Ce*Jei*Jei*nue+Ci*Jii*Jii*nui;

  muV_e=Taue*mu_e_mf+muEext; s2V_e=Taue*s2_e_mf+sEext;
  muV_i=Taui*mu_i_mf+muIext; s2V_i=Taui*s2_i_mf+sIext;
  
}
  
int main()
{
	
	float Cmax,Cmin,d0,kc,xv,d0e,d0i,Cv,Co,xs,ce,ci,d; 
	Cmax=1.6;Cmin=1.2;//强度因子最大和最小值
	d0=1;//性能最佳的浓度
	kc=0.12;//确定多巴胺调节的剧烈程度（即：对应倒U形的变化快慢）
	xv=0.185;//同步调节x0e和x0i(即调xv确保（x0e+x0i）/2=x0);
	
	d0e=d0-xv;d0i=d0+xv;//使得在x0两曲线的值等于Jei和Jee在最佳浓度的值；
	Cv=Cmax-Cmin;//强度最大偏移量
	Co=1-Cmin;//偏1因子
	
	xs=0;//归一化 
	ce=1/(1+Cv/(1+exp((d0e-xs)/kc)));
	ci=1/(1+Cv/(1+exp((d0i-xs)/kc)));

while(D1>=140) 
{
		
  d=(float)D1/100;			
  De=(1+Cv/(1+exp((d0e-d)/kc)))*ce-Co;
  Di=(1+Cv/(1+exp((d0i-d)/kc)))*ci-Co;
  De=(int)(De*100);Di=(int)(Di*100);//让De和Di只取两位有效精度； 
  //printf("%d,%.16f,%.16f; ",D1,De,Di);	
  printf("%d; ",D1);
  
  Jdep=0.10*De/100; Jpot=0.45*De/100;Jei=0.135*Di/100;
  Jie=-0.25f; Jii=-.20f; 
  

  int p, cue, test, nt, ciclo, count, stim_to_be[pops], running_stim;
  int i,j, stim_lastprevious;
  float *structuring, old_gss, old_gsn;
  FILE *fp;
  char name[50];
  
  (void)naive_mf();
  
  Randomize();

  InitV();
  InitEventMatrix();
 
  flag_fs=1; 
  
  if (flag_fs==1)
    {
      sel_old=createPop_noov(); makePop(sel_old); 
      initSynapses_nofs(); 
      invert_nofs_udo();
      InitSynapses();
    }
  else
    {
      sel_old=createPop(1); makePop(sel_old); 
      InitSynapses();
    }
      
  //SA trial
  flag_SA=1;
  ratepops=fopen("SA_rates_pops.log","w");
  spikeout=fopen("SA_rasters.log","w");
  synvar_u=fopen("SA_stp_u.log","w");
  synvar_x=fopen("SA_stp_x.log","w");
  curr_out=fopen("SA_currents.log","w");
  single_trial(0);
  shifttspike();
  
  fflush(ratepops); fflush(curr_out);//释放内存磁盘内存文件， 
  fflush(spikeout);fflush(synvar_u);fflush(synvar_x);
  
  fclose(ratepops); fclose(curr_out);
  fclose(spikeout); fclose(synvar_u); fclose(synvar_x);
  //end SA trial
  
  flag_SA=0;
  for (nt=0; nt<1; nt++)
    {
   	if ((D1>=100)&&(D1<1000))
     {sprintf(name,"rates_pops_%1d_%3d.log",nt,D1);
      ratepops=fopen(name,"w");
      sprintf(name,"rasters_%1d_%3d.log",nt,D1);
      spikeout=fopen(name,"w");
      sprintf(name,"stp_u_%1d_%3d.log",nt,D1);
      synvar_u=fopen(name,"w");
      sprintf(name,"stp_x_%1d_%3d.log",nt,D1);
      synvar_x=fopen(name,"w");
      sprintf(name,"currents_%1d_%3d.log",nt,D1);
      curr_out=fopen(name,"w");
     } 
     else 
     {
      if(D1<10)
	  {sprintf(name,"rates_pops_%1d_%1d.log",nt,D1);
      ratepops=fopen(name,"w");
      sprintf(name,"rasters_%1d_%1d.log",nt,D1);
      spikeout=fopen(name,"w");
      sprintf(name,"stp_u_%1d_%1d.log",nt,D1);
      synvar_u=fopen(name,"w");
      sprintf(name,"stp_x_%1d_%1d.log",nt,D1);
      synvar_x=fopen(name,"w");
      sprintf(name,"currents_%1d_%1d.log",nt,D1);
      curr_out=fopen(name,"w");
	  }
      else
      {sprintf(name,"rates_pops_%1d_%2d.log",nt,D1);
      ratepops=fopen(name,"w");
      sprintf(name,"rasters_%1d_%2d.log",nt,D1);
      spikeout=fopen(name,"w");
      sprintf(name,"stp_u_%1d_%2d.log",nt,D1);
      synvar_u=fopen(name,"w");
      sprintf(name,"stp_x_%1d_%2d.log",nt,D1);
      synvar_x=fopen(name,"w");
      sprintf(name,"currents_%1d_%2d.log",nt,D1);
      curr_out=fopen(name,"w");
      	
      }
     }
     
      
      cue=nt;
      single_trial(cue);
      shifttspike();

      fflush(ratepops); fflush(curr_out);//释放内存磁盘内存文件， 
      fflush(spikeout);fflush(synvar_u);fflush(synvar_x);

      fclose(ratepops); fclose(curr_out); 
      fclose(spikeout); fclose(synvar_u);  fclose(synvar_x);
    }
    
    LIFE=Tprestim_PT+Tdelay_PT+Tcue_PT;
    if ((D1>=100)&&(D1<1000))
      sprintf(name,"num_parameter_%1d_%3d.log",0,D1);
    else
      {
	   if(D1<10)
	    sprintf(name,"num_parameter_%1d_%1d.log",0,D1);
	    else
	    sprintf(name,"num_parameter_%1d_%2d.log",0,D1);
	  }
    
    num_para=fopen(name,"w");
    fprintf(num_para,"%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t",(Ni+Ne),PE,D1,f,dt,LIFE,muEext,muIext,Tprestim_PT,TCamp,Tcue_PT);
    fflush(num_para);
    fclose(num_para);
    
	D1=D1-10;   
}

/*
 free(Pop);free(P);free(D); free(J); free(Events); //释放内存空间 
 free(EventsNMDA); free(neur); free(gg); free(Inmda); free(Iext);free(Iext);
 free(duration_stim);free(duration_popout);free(spikes);free(curr);free(delay_stim);
 int i;
 for (i=0;i<Ne+Ni;i++)
 {
 if (i<Ne)
 {free(P[i]);free(D[i]);}
 	
 }
 int p;
 for (p=0; p<pops; p++)
    {   
     free(Pop[p].index);
    }
 for(i=0; i<(int)(MinDelay + MaxStep*StepDelay/dt+0.5)+1;i++)
    {
      free(Events[i]);
  	  free(EventsNMDA[i]);      
    }    */
}
