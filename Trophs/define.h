/*********************************************************************
**********************************************************************
****      <<< iRBC adhesion Code >>>                       	  ****
****                                                              ****
****        Anil Kumar Dasanna                                    ****
****        Prof. Dr. Ulrich Schwarz Group                        ****
****        Institute of Theoretical Physics                      ****
****        Heidelberg-69115, Germany                             ****
****        contact:anil.dasanna@bioquant.uni-heidelberg.de       ****
**********************************************************************
*********************************************************************/
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <stddef.h>
#include <iostream>
#include <new>
#include <vector>
#include <algorithm>
#include <fstream>
//Random number generation files
#include "Include/solfuncs.h"
#include "Include/memfuncs.h"
#include "Include/randomc.h"
#include "Include/stocc.h"
#include "Include/poly34.h"





#define Pi 3.141592653
#define randno 12345678

//1000 nm resolution
//viscosity = 21.5

//for solvent
#define a 1.0
#define Lx 40
#define Ly 30
#define Lz 30
#define dens 10
#define Ns (Lx*Ly*Lz*dens)
#define navg (1.0*dens)
#define dtcd 0.02	//time step
#define mass 1.0    //mass
#define Theta 0.1
//alpha = 135
#define cosalpha -0.70710678118
#define sinalpha 0.70710678118




// for membrane
#define MAXLINE 400   //used for # of lines in file
//tuning parameters//
#define Nt 1000  //no of vertices
#define shear (15.0*30.0)
#define kappab (140.0*Theta)
////////////////////
#define pmass (navg*mass) 	// in units of mass
#define dtmd 0.002	//in units of sqrt( m*a**2/K_bT_0)
#define mfactor 2.2
#define ks 400.0
#define dthr 0.35
#define FD 200.0
#define Fg 100.0
#define el0 0.05
#define memdis 0.0

#define PRmass (navg*mass)
#define paraR 2.0




using namespace std ;

typedef pair<double,int> Pair ;

 struct particle
 {

	double xpos ;
	double ypos ;
	double zpos ;
 	double vx ;
	double vy ;
	double vz ;
 	int pid ;
 	double mas ;
 	int indio ; //index that tells the particle is inside or outside
 	bool operator < (const particle& rhs) const { return indio < rhs.indio; }
 };


 struct cell
 {
	double xmax ;
 	double xmin ;
 	double ymax ;
 	double ymin ;
 	double zmax ;
 	double zmin ;
	double vcm[3] ;
	int nc ;
	double ncc ;
	double rand1,rand2 ;
	double dvv ;  //to compute scaling factor
	double ThE ;
	int flgThE ;
	// vertex list for bounce-back collisions
 	vector <int> llist ;
	vector <int> liglist ;
	vector <int> sol ;
	vector <int> mem ;
 };


struct vert
{
	double x ;
	double y ;
	double z ;
	double vx ;
	double vy ;
	double vz ;
	int pid ;
	int ix ;
	int iy ;
	int bondi ; //bond index
	int bondl ; // index of ligand
	double bondlen ; //initial bond length
	int knob ;
	double poff ;
	double pon ;
	double ti ;
	double forc ;
	vector < int > liglist ;
	vector < Pair > Pairs ;
	vector < int > nnlist ;
	vector < int > nlist ;
};

struct verto
{
	double x ;
	double y ;
	double z ;
};


struct linkk
{
	int l1 ;
	int l2 ;
	double l0 ;
	double kappap ;
	double elap ;
};

struct linkkk
{
	int l1 ;
	int l2 ;
};


struct triangle
{
	int l1 ;
	int l2 ;
	int l3 ;
	double cmx ;
	double cmy ;
	double cmz ;
};

struct force
{
	double fx ;
	double fy ;
	double fz ;
};


#ifndef MAIN
extern vector < CRandomMersenne > RanGen ;
#else
vector < CRandomMersenne > RanGen ;
#endif


#ifndef MAIN
extern vector < StochasticLib1 > sto ;
#else
vector < StochasticLib1 > sto ;
#endif


 //1351.351351351 : conversion factor
#ifndef MAIN
extern double onrate ;
#else
double onrate ;
#endif


#ifndef MAIN
extern double offrate ;
#else
double offrate ;
#endif


#ifndef MAIN
extern vector < vert > vts ;
#else
vector < vert > vts ;
#endif


#ifndef MAIN
extern vector < linkk > edges ;
#else
vector < linkk > edges ;
#endif

#ifndef MAIN
extern vector < linkkk > edgest ;
#else
vector < linkkk > edgest ;
#endif

#ifndef MAIN
extern vector < int > iedge ;
#else
vector < int > iedge ;
#endif

#ifndef MAIN
extern vector < triangle > tries ;
#else
vector < triangle > tries ;
#endif


#ifndef MAIN
extern vector < double > tarea ;
#else
vector < double > tarea ;
#endif



#ifndef MAIN
extern vector < force > fmapdt ;
#else
vector < force > fmapdt ;
#endif



#ifndef MAIN
extern double Sarea ;
#else
double Sarea ;
#endif



#ifndef MAIN
extern double Volume ;
#else
double Volume ;
#endif


#ifndef MAIN
extern double center[3] ;
#else
double center[3] ;
#endif


#ifndef MAIN
extern double lmax ;
#else
double lmax ;
#endif


#ifndef MAIN
extern double ll0 ;
#else
double ll0 ;
#endif


#ifndef MAIN
extern double Cq ;
#else
double Cq ;
#endif


#ifndef MAIN
extern double facT,factorT ;
#else
double facT,factorT ;
#endif


#ifndef MAIN
extern vector <particle> solvent ; //particles
#else
vector <particle> solvent ; //particles
#endif



#ifndef MAIN
extern vector <cell> Cells ; //Cells
#else
vector <cell> Cells ; //Cells
#endif


#ifndef MAIN
extern int cindex[(int)(Lx+2)][(int)(Ly+2)][(int)(Lz+2)];
#else
int cindex[(int)(Lx+2)][(int)(Ly+2)][(int)(Lz+2)];
#endif


#ifndef MAIN
extern double Sarea0 ;
#else
double Sarea0 ;
#endif


#ifndef MAIN
extern double Volume0 ;
#else
double Volume0 ;
#endif


#ifndef MAIN
extern int nsolb ;
#else
int nsolb ;
#endif


#ifndef MAIN
extern vector < double > xlig ;
#else
vector < double > xlig ;
#endif


#ifndef MAIN
extern vector < double > ylig ;
#else
vector < double > ylig ;
#endif

#ifndef MAIN
extern vector < int > ligi ;
#else
vector < int > ligi ;
#endif


#ifndef MAIN
extern double para[3] ;
#else
double para[3] ;
#endif

#ifndef MAIN
extern double paraV[3] ;
#else
double paraV[3] ;
#endif


#ifndef MAIN
extern double paraft[3] ;
#else
double paraft[3] ;
#endif


#ifndef MAIN
extern double parafdt[3] ;
#else
double parafdt[3] ;
#endif


#ifndef MAIN
extern int paraid ;
#else
int paraid ;
#endif

#ifndef MAIN
extern int paradx ;
#else
int paradx ;
#endif

#ifndef MAIN
extern int parady ;
#else
int parady ;
#endif


#ifndef MAIN
extern double theta0,avg_edge ;
#else
double theta0,avg_edge ;
#endif


#ifndef MAIN
extern double kappasa ;
#else
double kappasa ;
#endif


#ifndef MAIN
extern double kappad ;
#else
double kappad ;
#endif

#ifndef MAIN
extern double kappav,wallspeed ;
#else
double kappav,wallspeed,walls ;
#endif

#ifndef MAIN
extern int knob_size,sold ;
#else
int knob_size,sold ;
#endif
