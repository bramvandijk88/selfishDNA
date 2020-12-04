#include <vector>
#include <string>
#include <list>
#include "../PoaS/Genome.hh" // Includes the header needed for PoaS-genomes
#include "../PoaS/eDNA.hh"   // Includes the header needed for eDNA objects (a vector of StringOfPearls)
extern "C"{
#include "mersenne.h"
}

/******************************************************
 ************Here you define your own type*************
 ******************************************************/

#ifndef CASH2_HPP
#define CASH2_HPP
#define NUMABS 20

typedef struct __type2{
/*****************************************
 *** Val = State (0 = dead, >0 living) ***
 ******** Fval5 = Fitness storage ********
 *****************************************/

  int val, val2,val3,val4,val5;
  float fval, fval2, fval3, fval4, fval5;
  float fitness;
  float ab_effect;
  Genome * G;							// The class Genome is defined in the OOP files in PoaS
  eDNA * DNA;						// The class eDNA is defined in the OOP files, its just a list of small genomes


} TYPE2;


/*****************************************************
 *************Do not change from here*****************
 *****************************************************/

/* This structure is used in asynchronous updating to determine
   the order of the cells to update. Look at the functions in
   chash2.c which use "struct updorder" for details and the usage
   of this structure. Look at UpdOrdReset() in cash2.c  */
struct updorder {
  int row;
  int col;
};

/* This structure is used in neighborhoods retrieving. Look
   at NeighSet() in cash2.c */
struct point {
  int row;
  int col;
};

/* This structure is used in Margolus neighborhood
   retrieving. Read MarGolusNeigh(). If you say "MARGOLUS
   something_you_defined[*][i][j]" where [*] denote even-or-odd
   phase, then "----[*][i][j].m[CW].row" and
   "----[*][i][j].m[CW].col" give you the coordinate of the CW
   neighbor (see below) of Margolus neighborhood of the [i][j]
   cell in the phase [*] (even/odd, 0 or 1).

   -----------
   |HERE| CW |
   -----------
   | CCW| OPP|
   -----------

   Note that the exact place of "HERE" does not matter in that
   "CW", "OPP" and "CCW" notations are invariant for the position
   of "HERE".

   What is even and what is odd? Let "SYD" denote
   "something_you_defined". Then, SYD[0] and SYD[1] gives you tow
   alternative Margolus partitioning of the plane. Then, SYD[0]'s
   partition will be such that [0][0] cell of your plane will be
   the upper-left corner (UL) of the 2x2 square. SYD[1]'s
   partition will be such that [0][0] cell  of your plane will be
   the lower-right corner (LR) of the 2x2 square. See the
   followings too.

   UL,UR,LL,LR are like,
   -----------
   | UL | UR |
   -----------
   | LL | LR |
   -----------


   in even phase (SYD[0]), the plane is divided as follows:
   -----------------
   | 1,1 | 1,2 | etc.
   -------------
   | 2,1 | 2,2 | etc.
   -------------
   | etc. etc.   etc.

   in odd phase (SYD[1]) and boundary=WRAP, the plane is divided as follows:
   ------------------
   |nr,nc| nr,1| etc.
   -------------
   | 1,nc|  1,1| etc.
   -------------
   | etc. etc.   etc.

   where "nr" is the number of row, and "nc" is the number of
   column. If boundary==FIXED, then read the avobe table such
   that "nr=nc=0".
*/
#define CW (0)
#define OPP (1)
#define CCW (2)
typedef struct __margolus{
  struct point m[3];/* 0:clockwise 1:opposite 2:counter-clockwise */
} MARGOLUS;


/*********************************************basic*/
TYPE2 **NewP2(void);
TYPE2 **New2(void);
double **NewDB(void);
int PlaneFree2(TYPE2**);
void UpdOrdReset(struct updorder**,bool reinit);
TYPE2 **Copy2(TYPE2 **,TYPE2 **);
TYPE2 **Fill2(TYPE2 **,TYPE2);
void DiffusionDB(double **,double,struct point**[]);
int InitialPlaneSet(TYPE2**,int,TYPE2,...);

/*********************************************shift*/
TYPE2 **Boundaries2(TYPE2**);
double **BoundariesDB(double **,double);
/*****************************************neighbors*/
struct point **NeighSet(int);
void NeighFree(struct point**);
TYPE2 ***NeighSetP(TYPE2**,int);
void NeighFreeP(TYPE2***);

/******************************************margolus*/
MARGOLUS ***MargolusNeigh(void);
void MargolusFree(MARGOLUS***);
void MargolusDiffusion(TYPE2**,MARGOLUS***,int);

/*********************************************noise*/
void PerfectMix(TYPE2**); /* Shake() in CASH */

#endif
