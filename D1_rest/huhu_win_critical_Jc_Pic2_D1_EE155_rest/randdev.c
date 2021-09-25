
/*----------------------------------------------*
 *                                              *
 *   RANDDEV.C                                  *
 *                                              *
 *   Questa e' una libreria per la generazione  *
 *   di numeri pseudocasuali con varie distri-  *
 *   buzioni (RANDom DEViations).               *
 *                                              *
 *   Realizzato da Mattia Maurizio.             *
 *   Iniziazto il  21 novembre 1993.            *
 *                                              *
 *----------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "randdev.h"



static int    TimeSeed = 0;      /* Seme funzione del tempo del generatore. */
static int rand49_idum = -77531; /* Varibile di appoggio per gli algoritmi. */



#ifdef RAND2

/*------------------------------------------------*
 *                                                *
 *   rand2 e Random                               *
 *                                                *
 *   Riporta  un float casuale  nell'intervallo   *
 *   [0,1[ con DISTRIBUZIONE UNIFORME. L'algo-    *
 *   ritmo usato e' da rand2 nel Numerical Reci-  *
 *   pice.                                        *
 *   NOTA: gli int devono essere a 32bit (Watcom  *
 *         e Unix, no Borland).                   *
 *------------------------------------------------*/

#define M 714025
#define IA 1366
#define IC 150889

float rand2 (int * idum) {

	static long iy,ir[98];
	static int iff=0;
	int j;

	if (*idum < 0  || iff == 0) {
		iff=1;
		if ((*idum=(IC-(*idum)) % M) < 0) *idum = -(*idum);
		for (j=1;j<=97;j++) {
			*idum=(IA*(*idum)+IC) % M;
			ir[j]=(*idum);
		}
		*idum=(IA*(*idum)+IC) % M;
		iy=(*idum);
	}
	j=(int)(1 + 97.0*iy/M);
	if (j > 97 || j < 1) printf("RAN2: This cannot happen.");
	iy=ir[j];
	*idum=(IA*(*idum)+IC) % M;
	ir[j]=(*idum);
	return (float) iy/M;
}

#undef M
#undef IA
#undef IC

float Random (void)  { return rand2(&rand49_idum); }

#endif



#ifdef RAND3

/*------------------------------------------------------*
 *                                                      *
 *   rand3 e Random                                     *
 *                                                      *
 *   Riporta  un float casuale  nell'intervallo         *
 *   [0,1[ con DISTRIBUZIONE UNIFORME.                  *
 *   E' l'algoritmo preso da numerical recipice dovuto  *
 *   a knuth e nominato ran3. E' molto veloce ma poco   *
 *   studiato (non sanno bene che tipo di correlazioni  *
 *   puo' indurre).                                     *
 *   NOTA: gli int devono essere a 32bit (Watcom        *
 *         e Unix, no Borland). Altrimenti dove e' com- *
 *         mentato 'long' inserire 'long'.              *
 *------------------------------------------------------*/

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)


double rand3(int * idum)
{
   static int  inext, inextp;
   static int  ma[56];        /*** long ***/
   static int  iff = 0;
          int  mj, mk;        /*** long ***/
          int  i, ii, k;
/*** Da togliere in fase di debug.
   double out;
***/

   /*** Percorso di inizializzazione. ***/
   if (*idum < 0 || iff == 0) {
      iff=1;
      mj=MSEED-(*idum < 0 ? -*idum : *idum);

      /*** Modifica mia. ***/
      mj = (mj < 0 ? -mj : mj); /*** Questa aggiunta l'ho fatta io, perche' ***
                                 *** se non ci fosse non vengono numeri tra ***
                                 *** 0.0 e 1.0. Io penso che in NR ci sia   ***
                                 *** un errore di interpretazione sul si-   ***
                                 *** gnificato di '%': pensano che riporti  ***
                                 *** sempre un numero positivo.             ***
                                 *** Se non ci credi prova con questo seed  ***
                                 *** 878104453 (da passare a setRanomSeed). ***/
      mj %= MBIG;
      ma[55]=mj;
      mk=1;
      for (i=1;i<=54;i++) {
         ii=(21*i) % 55;
         ma[ii]=mk;
         mk=mj-mk;
         if (mk < MZ) mk += MBIG;
         mj=ma[ii];
      }
      for (k=1;k<=4;k++)
         for (i=1;i<=55;i++) {
            ma[i] -= ma[1+(i+30) % 55];
            if (ma[i] < MZ) ma[i] += MBIG;
         }
      inext=0;
      inextp=31;
      *idum=1;
   }

   if (++inext == 56) inext=1;
   if (++inextp == 56) inextp=1;
   mj=ma[inext]-ma[inextp];
   if (mj < MZ) mj += MBIG;
   ma[inext]=mj;

/*** Da scommentare in fase di debug.
   out = mj*FAC;
   if (out < 0 || out >= 1) {
      fprintf(stderr, "ERROR: Bad random number.");
      abort();
   }
   return out;
***/

   return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

double Random (void)  { return rand3(&rand49_idum); }

#endif /* of RAND3 */



#ifdef RAND1

/*------------------------------------------------*
 *                                                *
 *   UniDev.                                      *
 *                                                *
 *   Riporta  un  int  casuale  nell'intervallo   *
 *   [0,RAND_MAX] con DISTRIBUZIONE UNIFORME.     *
 *   Per inizializzare  il seme del  generatore   *
 *   utilizzare la seguente sintassi:             *
 *      UniDev( &Seed )                           *
 *   dove Seed deve essere una variabile intera   *
 *   contenente il seme iniziale. Per mantenere   *
 *   intatto il seme basta richiamare la funzio-  *
 *   ne con la seguente sintassi:                 *
 *      UniDev( NULL )                            *
 *   Questa procedura si basa sull'algoritmo di   *
 *   BAYS e DURHAM.                               *
 *                                                *
 *   NOTE: Il periodo del  generatore e' pratica- *
 *         infinito, ma i numeri forniti non so-  *
 *         no piu' di RAND_MAX+1.                 *
 *         Deve essere usata con i piedi di piom- *
 *         bo quando si utilizzano i bit meno si- *
 *         gnificativi.                           *
 *                                                *
 *------------------------------------------------*/

int UniDev (int * SeedPntr) /* Puntatore al seme del generatore. Se NULL il generatore *
                             * non viene inizializzato.                                */

{
   /*** Dichiarazione delle variabili locali. ***/
   #define SEQUENCE_DIM 98.0
   static int y, Sequence[(int)SEQUENCE_DIM];
   static int FirstCall = 1;         /* Prima chiamata della funzione? */
   static int j;                     /* Contatore.                     */

   /*** Inizializzazione della sequenza. ***/
   if (SeedPntr != NULL || FirstCall)
     {
      FirstCall = 0;
      if (SeedPntr != NULL) srand(*SeedPntr);
      for (j=0; j<(int)SEQUENCE_DIM; j++) rand();
      for (j=0; j<(int)SEQUENCE_DIM; j++) Sequence[j] = rand();
      y = rand();
     }

   /*** Generazione del numero ps eudo-casuale. ***/
   j           = (SEQUENCE_DIM*y)/(RAND_MAX+1);
   y           = Sequence[j];
   Sequence[j] = rand();
   return y;
}

#endif /* of RAND1 */



/*----------------------------------------------*
 *                                              *
 *   SetRandomSeed (int Seed)                   *
 *                                              *
 *   Assegna un valore al seme del generatore   *
 *   numeri pseudo-casuali.                     *
 *                                              *
 *----------------------------------------------*/

void SetRandomSeed (int Seed)

{
   TimeSeed = Seed;
   
#ifdef RAND2
   rand49_idum=-Seed;
   rand2(&rand49_idum);
#endif
#ifdef RAND3
   rand49_idum=-Seed;
   rand3(&rand49_idum);
#endif
#ifdef RAND1
   UniDev(&TimeSeed);
#endif
}



/*------------------------------------------------------*
 *   GetRandomSeed ()                                   *
 *                                                      *
 *   Riporta il seme del generatore di numeri pseudo-   *
 *   casuali TimeSeed.                                  *
 *------------------------------------------------------*/

int GetRandomSeed (void)

{
   return TimeSeed;
}



/*----------------------------------------------*
 *                                              *
 *   Randomize.                                 *
 *                                              *
 *   Inizializza in modo casuale (attraverso    *
 *   la funzione time) il generatore di nume-   *
 *   ri pseudo-casuali UniDev.                  *
 *                                              *
 *----------------------------------------------*/

void Randomize (void)

{
   SetRandomSeed((unsigned)time(NULL));
}



/*------------------------------------------------*
 *                                                *
 *   ExpDev.                                      *
 *                                                *
 *   Riporta un numero  casuale reale positivo    *
 *   (float) con  distribuzione esponenziale e    *
 *   media  unitaria, usando come  sorgente di    *
 *   numeri casuali con distribuzione uniforme    *
 *   Random.                                      *
 *     Fun. di distribuzione = e^-x con x>=0.0    *
 *   Questa procedura si basa sull'algoritmo di   *
 *   Montecarlo.                                  *
 *                                                *
 *------------------------------------------------*/

float ExpDev ( void )

{
   /*** Dichiarazione delle variabili locali. ***/
   float Dum;      /* Variabile di appoggio. */

   /*** Ciclo di generazione del singolo numero. ***/
   do
      Dum = (float)Random();
   while (Dum == 0.0);
   return (float)-log(Dum);
}



/*------------------------------------------------*
 *                                                *
 *   NormDev.                                     *
 *                                                *
 *   Riporta un  numero  casuale reale (float)    *
 *   con distribuzione gaussiana,  media nulla    *
 *   e varianza unitaria. Come sorgente di nu-    *
 *   meri  casuali con  distribuzione uniforme    *
 *   viene usata Random().                        *
 *     Fun. di distr. = 1/(sqrt(2pi))*e^-x^2      *
 *   Questa procedura si basa sull'algoritmo di   *
 *   BOX-MULLER.                                  *
 *                                                *
 *------------------------------------------------*/

float NormDev ( void )

{
   /*** Dichiarazione delle variabili locali. ***/
   static int   Set = 0;
   static float IIran;
   float Fac, r, v1, v2;

   /*** Decisione su quale set prendere il singolo numero. ***/
   if (Set == 0)
     {do
        {v1 = (float)(2.0 * Random() - 1.0);
         v2 = (float)(2.0 * Random() - 1.0);
         r  = v1 * v1 + v2 * v2;}
      while (r >= 1.0 || r==0.0);
      Fac   = (float)sqrt(-2.0 * log(r) / r);
      IIran = v1 * Fac;
      Set   = 1;
      return v2 * Fac;}
   else
     {Set = 0;
      return IIran;}
}

