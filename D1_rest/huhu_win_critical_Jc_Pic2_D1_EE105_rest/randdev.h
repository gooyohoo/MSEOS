/*------------------------------------------------*
 *                                                *
 *   RANDDEV.H                                    *
 *                                                *
 *   Questa e' l'intestazione della libreria per  *
 *   la generazione di numeri pseudocasuali con   *
 *   varie distribuzioni.                         *
 *                                                *
 *   Realizzato da Mattia Maurizio.               *
 *   Iniziazto il  21 novembre 1993.              *
 *                                                *
 *------------------------------------------------*/


/*** Solo una delle tre seguenti define deve esistere. ***/
#define _RAND1       /* Se e' definita questa variabile un numero pseudo *
                     * casuale viene generato con l'algoritmo di ran1   *
                     * di numerical recipice.                           */
#define _RAND2       /* Se e' definita questa variabile un numero pseudo *
                     * casuale viene generato con l'algoritmo di ran2.  *
                     * di numerical recipice.                           */
#define RAND3      /* Se e' definita questa variabile un numero pseudo *
                     * casuale viene generato con l'algoritmo di ran3.  *
                     * di numerical recipice.                           */



#ifdef RAND2

/*------------------------------------------------*
 *                                                *
 *   Random.                                      *
 *                                                *
 *   Riporta  un float casuale  nell'intervallo   *
 *   [0,1[ con DISTRIBUZIONE UNIFORME. L'algo-    *
 *   ritmo usato e' da rand2 nel Numerical Reci-  *
 *   pice.                                        *
 *   NOTA: gli int devono essere a 32bit (Watcom  *
 *         e Unix, no Borland).                   *
 *------------------------------------------------*/

float rand2 (int * idum);

#endif /* of RAND2 */



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
 *         e Unix, no Borland).                         *
 *------------------------------------------------------*/

double rand3(int * idum);

#endif



#ifdef RAND1

/*------------------------------------------------*
 *                                                *
 *   UniDev.                                      *
 *                                                *
 *   Riporta  un  int  casuale  nell'intervallo   *
 *   [0,RAND_MAX] con distrubuzione uniforme.     *
 *   Per inizializzare  il seme del  generatore   *
 *   utilizzare la seguente sintassi:             *
 *      RandBD( &Seed )                           *
 *   dove Seed deve essere una variabile intera   *
 *   contenente il seme iniziale. Per mantenere   *
 *   intatto il seme basta richiamare la funzio-  *
 *   ne con la seguente sintassi:                 *
 *      RansBD( NULL )                            *
 *   Questa procedura si basa sull'algoritmo di   *
 *   BAYS e DURHAM.                               *
 *                                                *
 *   NOTE: Il periodo del  generatore Š pratica-  *
 *         infinito, ma i numeri forniti non so-  *
 *         no pi— di RAND_MAX+1.                  *
 *         Deve essere usata con i piedi di piom- *
 *         bo quando si utilizzano i bit meno si- *
 *         gnificativi.                           *
 *                                                *
 *------------------------------------------------*/

int UniDev (int * SeedPntr);

#endif



/*----------------------------------------------*
 *                                              *
 *   SetRandomSeed (int Seed)                   *
 *                                              *
 *   Assegna un valore al seme del generatore   *
 *   numeri pseudo-casuali.                     *
 *                                              *
 *----------------------------------------------*/

void SetRandomSeed (int Seed);



/*------------------------------------------------------*
 *   GetRandomSeed ()                                   *
 *                                                      *
 *   Riporta il seme del generatore di numeri pseudo-   *
 *   casuali TimeSeed.                                  *
 *------------------------------------------------------*/

int GetRandomSeed (void);



/*----------------------------------------------*
 *                                              *
 *   Randomize.                                 *
 *                                              *
 *   Inizializza in modo casuale (attraverso    *
 *   la funzione time) il generatore di nume-   *
 *   ri pseudo-casuali UniDev.                  *
 *                                              *
 *----------------------------------------------*/

void Randomize (void);



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

float ExpDev ( void );



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

float NormDev ( void );



/*----------------------------------------------------*
 *   Macro Utili:                                     *
 *                                                    *
 *     Random    - genera un numero casuale con di-   *
 *                 stribuzione uniforme nell'inter-   *
 *                 vallo [0,1[.                       *
 *----------------------------------------------------*/

#ifdef RAND2
float Random (void);
#endif


#ifdef RAND3
double Random (void);
#endif


#ifdef RAND1
#define Random()   UniDev(NULL)/(RAND_MAX+1.0)
#endif
