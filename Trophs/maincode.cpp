/*********************************************************************
**********************************************************************
****      <<< iRBC deformation adhesion Code >>>                  ****
****                                                              ****
****        Anil Kumar Dasanna                                    ****
****        Prof. Dr. Ulrich Schwarz Group                        ****
****        Institute of Theoretical Physics                      ****
****        Heidelberg-69115, Germany                             ****
****        contact:anil.dasanna@bioquant.uni-heidelberg.de       ****
**********************************************************************
*********************************************************************/
#define MAIN
#include "define.h"


void Exclusion()
{
	double cutoff = 1.5*avg_edge ;
	double epsilon = 1.0 ;
	double sigma = 0.890898718*cutoff ;
	#pragma omp parallel default(shared)
	{
		double A[3],X[3] ;
		X[0] = X[1] = X[2] = 0.0 ;
		#pragma omp for
		for ( int i = 0 ; i < vts.size() ; ++i )
		{
			int fx = (int)floor(vts[i].x+1.0) ;
			int fy = (int)floor(vts[i].y+1.0) ;
			int fz = (int)floor(vts[i].z+1.0) ;
			int pid = fx*(int)(Ly+2)*(int)(Lz+2) + fy*(int)(Lz+2) + fz ;
			for ( unsigned int j = 0 ; j < Cells[pid].llist.size() ; ++j )
			{
				int jj = Cells[pid].llist[j] ;
				int p1 = tries[jj].l1 ;
				int p2 = tries[jj].l2 ;
				int p3 = tries[jj].l3 ;
				if(find(vts[i].nlist.begin(), vts[i].nlist.end(), p1) == vts[i].nlist.end())
				{
					A[0] = (vts[i].x + vts[i].ix*Lx) - (vts[p1].x + vts[p1].ix*Lx) ;
					A[1] = (vts[i].y + vts[i].iy*Ly) - (vts[p1].y + vts[p1].iy*Ly) ;
					A[2] = vts[i].z - vts[p1].z ;
					double aa =  A[0]*A[0] + A[1]*A[1] + A[2]*A[2]  ;
					if ( aa <= cutoff*cutoff)
					{
						double ex = (sigma * sigma * sigma * sigma * sigma * sigma)/(aa * aa * aa);
						X[0] = (24.0 * epsilon/aa) * ((2.0 * ex * ex) - ex ) * A[0];
						X[1] = (24.0 * epsilon/aa) * ((2.0 * ex * ex) - ex ) * A[1];
						X[2] = (24.0 * epsilon/aa) * ((2.0 * ex * ex) - ex ) * A[2];
						fmapdt[i].fx += X[0] ;
						fmapdt[i].fy += X[1] ;
						fmapdt[i].fz += X[2] ;
					}
				}
				if(find(vts[i].nlist.begin(), vts[i].nlist.end(), p2) == vts[i].nlist.end())
				{
					A[0] = (vts[i].x + vts[i].ix*Lx) - (vts[p2].x + vts[p2].ix*Lx) ;
					A[1] = (vts[i].y + vts[i].iy*Ly) - (vts[p2].y + vts[p2].iy*Ly) ;
					A[2] = vts[i].z - vts[p2].z ;
					double aa =  A[0]*A[0] + A[1]*A[1] + A[2]*A[2]  ;
					if ( aa <= cutoff*cutoff)
					{
						double ex = (sigma * sigma * sigma * sigma * sigma * sigma)/(aa * aa * aa);
						X[0] = (24.0 * epsilon/aa) * ((2.0 * ex * ex) - ex ) * A[0];
						X[1] = (24.0 * epsilon/aa) * ((2.0 * ex * ex) - ex ) * A[1];
						X[2] = (24.0 * epsilon/aa) * ((2.0 * ex * ex) - ex ) * A[2];
						fmapdt[i].fx += X[0] ;
						fmapdt[i].fy += X[1] ;
						fmapdt[i].fz += X[2] ;
					}
				}
				if(find(vts[i].nlist.begin(), vts[i].nlist.end(), p3) == vts[i].nlist.end())
				{
					A[0] = (vts[i].x + vts[i].ix*Lx) - (vts[p3].x + vts[p3].ix*Lx) ;
					A[1] = (vts[i].y + vts[i].iy*Ly) - (vts[p3].y + vts[p3].iy*Ly) ;
					A[2] = vts[i].z - vts[p3].z ;
					double aa =  A[0]*A[0] + A[1]*A[1] + A[2]*A[2]  ;
					if ( aa <= cutoff*cutoff)
					{
						double ex = (sigma * sigma * sigma * sigma * sigma * sigma)/(aa * aa * aa);
						X[0] = (24.0 * epsilon/aa) * ((2.0 * ex * ex) - ex ) * A[0];
						X[1] = (24.0 * epsilon/aa) * ((2.0 * ex * ex) - ex ) * A[1];
						X[2] = (24.0 * epsilon/aa) * ((2.0 * ex * ex) - ex ) * A[2];
						fmapdt[i].fx += X[0] ;
						fmapdt[i].fy += X[1] ;
						fmapdt[i].fz += X[2] ;
					}
				}
			}
		}
	}
}



void bond_dynamics(int mpcds,FILE *fp)
{
	double pon ;
	double B[3],bb ;
	double koff ;

	for ( int ii = 0 ; ii < vts.size() ; ++ii )
	{
		if (vts[ii].bondi == 1 )   //already bonded
		{
			double randoff = RanGen[0].Random() ;
			int jj = vts[ii].bondl ;
			int p1modx = vts[ii].ix ;
			int p1mody = vts[ii].iy ;
			if (vts[ii].poff > randoff)
			{
  	 			vts[ii].bondi = 0 ;
  	 			vts[ii].bondl = randno ;
  	 			ligi[jj] = 0 ; //ligand is freed
  	 			//For foot printing
				fprintf( fp,"%lf\t%lf\t%lf\t%lf\t%lf\n",(double)mpcds*dtmd,xlig[jj] + (double)p1modx*Lx,ylig[jj] + (double)p1mody*Ly,((double)mpcds*dtmd)-vts[ii].ti,vts[ii].forc) ;
				fflush(fp) ;
			}
		}
		//for vacant ligands
		else if ( ( vts[ii].z < dthr ) && (vts[ii].bondi == 0) && (vts[ii].knob == 1) )
		{
			for ( unsigned int jkl = 0 ; jkl < vts[ii].Pairs.size() ; ++jkl )
			{
					int j = vts[ii].Pairs[jkl].second ;
					double pon = vts[ii].pon ;
					double randon = RanGen[0].Random() ;
					if (( pon > randon ) && (ligi[j] == 0))
					{
						ligi[j] = 1 ; //ligand is filled
						vts[ii].bondi = 1 ;
						vts[ii].bondl = j ;
						vts[ii].ti = (double)mpcds*dtmd ;
						break ;
					}
			 }
		}
	}
}


void write_lammpsfile (FILE *flamp,int filec)
{
	fprintf (flamp,"ITEM: TIMESTEP\n");
	fprintf (flamp,"%d\n",filec) ;
	fprintf (flamp,"ITEM: NUMBER OF ATOMS\n") ;
	fprintf (flamp,"%d\n",vts.size()) ;
	fprintf (flamp,"ITEM: BOX BOUNDS pp pp ff\n") ;
	fprintf (flamp,"%lf %lf\n",0.0,(double)Lx) ;
	fprintf (flamp,"%lf %lf\n",0.0,(double)Ly) ;
	fprintf (flamp,"%lf %lf\n",0.0,(double)Lz) ;
	fprintf (flamp,"ITEM: ATOMS id type x y z ix iy iz\n");
	int type = 1 ;
	for ( int i = 0 ; i < vts.size() ; ++i )
	{
		type = 1 ;
		if ( vts[i].knob == 1) { type = 2 ;}
		fprintf(flamp, "%d %d %lf %lf %lf %d %d %d\n",i+1,type,vts[i].x,vts[i].y,vts[i].z,0,0,0) ;
	}
}


double set_kp(double eel0)
{
	double kp ;
	double lmax = mfactor*eel0 ;
	double x0 = eel0/lmax ;
	double onemx0 = 1.0 - x0 ;
	double shear1 = 1.0/(4.0*onemx0*onemx0) ;
	double shear2 = x0/(2.0*onemx0*onemx0*onemx0) ;
	double alph = shear1 - 0.25 + x0 ;
	double fact1 = ( shear2 - shear1 + 0.25 )*1.732050/(4.0*eel0*eel0*eel0*alph) ;
	double fact2 = 3.0*1.73205/(4.0*eel0*eel0*eel0) ;
	double fact = fact1+fact2 ;
	kp = shear/fact ;
	double pp = eel0*eel0*Theta*alph/kp ;

	return kp ;
}

double set_lp(double eel0)
{
	double lp ;
	double lmax = mfactor*eel0 ;
	double x0 = eel0/lmax ;
	double onemx0 = 1.0 - x0 ;
	double shear1 = 1.0/(4.0*onemx0*onemx0) ;
	double shear2 = x0/(2.0*onemx0*onemx0*onemx0) ;
	double alph = shear1 - 0.25 + x0 ;
	double fact1 = ( shear2 - shear1 + 0.25 )*1.732050/(4.0*eel0*eel0*eel0*alph) ;
	double fact2 = 3.0*1.73205/(4.0*eel0*eel0*eel0) ;
	double fact = fact1+fact2 ;
	double kpp = shear/fact ;
	lp = eel0*eel0*Theta*alph/kpp ;

	return lp ;
}



int main(int argc, char *argv[])
{
   int seed = (int) time(0) ;
   int nthreads = omp_get_max_threads();
   for ( int i = 0 ; i < nthreads ; ++i )
   {
	   int seedn = seed + i*20 ;
	   RanGen.push_back(CRandomMersenne(seedn)) ;
	   sto.push_back(StochasticLib1(seedn)) ;
   }
   char filename[256];
   int mm ;
   int fx,fy,fz,pi,nn ;
   double A[3],B[3],aa ;

   cout <<"List of arguments: v_wall_z file_index onrate offrate knob_size" << endl ;

   walls = atof(argv[1]) ;   //speed of the wall
   onrate = (walls/Lz)*0.5 ;
   offrate = (walls/Lz)*0.5 ;
   //onrate = atof(argv[3]) ;
   //offrate = atof(argv[4]) ;


  //ligands positions
  double xdelta = 0.5 ;
  double ydelta = 0.5 ;
  double xmin = 0.0 ;
  double ymin = 0.0 ;

  while ( xmin < Lx )
  {
  	while ( ymin < Ly )
  	{
  		xlig.push_back(xmin) ;
  		ylig.push_back(ymin) ;
  		ligi.push_back(0) ;
  		ymin += ydelta ;
  	}
  	xmin += xdelta ;
  	ymin = 0 ;
  }
   //Solvent initialization
	double vmax = sqrt(Theta/mass) ;
	int nrt = omp_get_thread_num() ;
 	for ( int j = 0 ; j < Ns ; ++j )
 	{
 	 	solvent.push_back(particle()) ;
		solvent[j].xpos = (Lx*RanGen[nrt].Random()) ;
 	 	solvent[j].ypos = (Ly*RanGen[nrt].Random()) ;
		solvent[j].zpos = (Lz*RanGen[nrt].Random()) ;
		solvent[j].vx = vmax*sto[nrt].Normal(0,1) ;
		solvent[j].vy = vmax*sto[nrt].Normal(0,1) ;
		solvent[j].vz = vmax*sto[nrt].Normal(0,1) ;
		solvent[j].mas = mass ;
	}



 	 ///////////////////Particles are placed////////////////
 	int count = -1 ;
	int lxi = (int)(Lx/a) + 2 ;
	int lyi = (int)(Ly/a) + 2 ;
	int lzi = (int)(Lz/a) + 2 ;
	for ( int i = 0 ; i < lxi ; ++i )
	{
		for ( int j = 0 ; j < lyi ; ++j )
		{
			 for (int k = 0 ; k < lzi ; ++k)
			 {
				   Cells.push_back(cell()) ;
				   count += 1 ;
				   Cells[count].xmin =  -1.0 + (double)i*a ;
				   Cells[count].xmax =  -1.0 + (double)i*a + a ;
				   Cells[count].ymin =  -1.0 + (double)j*a ;
				   Cells[count].ymax =  -1.0 + (double)j*a + a ;
				   Cells[count].zmin =  -1.0 + (double)k*a ;
				   Cells[count].zmax =  -1.0 + (double)k*a + a ;
				   fx = (int)floor(Cells[count].xmin + 1.0) ;
				   fy = (int)floor(Cells[count].ymin + 1.0) ;
				   fz = (int)floor(Cells[count].zmin + 1.0) ;
				   cindex[fx][fy][fz] = count ;
			 }
		 }
	 }

	//FILE *fle = fopen("liglist.dat","w");
	 for ( int i = 0 ; i < Cells.size() ; ++i )
 	 {
		 if ( Cells[i].zmax <= 2.0)
		 {
			 for ( int j = 0 ; j < xlig.size() ; ++j )
			 {
				A[0] = (Cells[i].xmin+Cells[i].xmax)*0.5 - xlig[j] ;
				A[1] = (Cells[i].ymin+Cells[i].ymax)*0.5 - ylig[j] ;
				if ( xlig[j] > Lx-a )
				{
					A[0] += ((Cells[i].xmin+Cells[i].xmax)*0.5 < a)?Lx:0 ;
				}
				if ( xlig[j] < a )
				{
					A[0] -= ((Cells[i].xmin+Cells[i].xmax)*0.5 > Lx-a)?Lx:0 ;
				}
				if ( ylig[j] > Ly-a )
				{
					A[1] += ((Cells[i].ymin+Cells[i].ymax)*0.5 < a)?Ly:0 ;
				}
				if ( ylig[j] < a )
				{
					A[1] -= ((Cells[i].ymin+Cells[i].ymax)*0.5 > Ly-a)?Ly:0 ;
				}
				double xx = sqrt ( A[0]*A[0] + A[1]*A[1] ) ;
				 if ( xx < 1.5*a )
				 {
					 Cells[i].liglist.push_back(j) ;
				 }
			 }
			// fprintf(fle,"%d %d\n",i,Cells[i].liglist.size() ) ;
		 }
	 }
	//fclose(fle) ;

   for ( int i = 0 ; i <  Nt ; ++i )
   {
	   vts.push_back(vert()) ;
	   fmapdt.push_back(force()) ;
   }


   //Reading of Data from file...
   	FILE *fpv,*fp ;  //file pointer for membrane files
	double xx, yy, zz ;
	double xxx,yyy,zzz ;
	double Center[3] ;
	int mapi ;
	count = -1 ;
	vector < double > zcrds ;
	int mm1 = sprintf(filename,"V%d.dat",Nt) ;
	fpv = fopen(filename,"r" ) ;
	for ( int i = 0 ; i < Nt ; ++i )
	{
		int mmm = fscanf(fpv,"%lf\t%lf\t%lf\n",&xxx,&yyy,&zzz ) ;
		zcrds.push_back(zzz) ;
		vts[i].x = xxx ;
		vts[i].y = yyy ;
		vts[i].z = zzz ;
		vts[i].vx =  0 ;
		vts[i].vy =  0 ;
		vts[i].vz =  0 ;
		fmapdt[i].fx = 0 ;
		fmapdt[i].fy = 0 ;
		fmapdt[i].fz = 0 ;
		vts[i].ix = 0 ;
		vts[i].iy = 0 ;
		vts[i].bondi = 0 ;
		vts[i].bondl = randno ;
		vts[i].bondlen = 0 ;
		vts[i].knob = 0 ;
	}
	fclose(fpv) ; //closing of vertices file
	sort(zcrds.begin(),zcrds.end()) ;
	Center[2] = -zcrds.front() + 0.6 ; //initial height
	Center[0] = Lx*0.5 ;
	Center[1] = Ly*0.5 ;

	for ( int i = 0 ; i < vts.size() ; ++i )
	{
		vts[i].x += Center[0] ;
		vts[i].y += Center[1] ;
		vts[i].z += Center[2] ;
	}

	double numer = 1.73205080757*((double)Nt - 2.0) - 5.0*Pi ;
	double denom = 1.73205080757*((double)Nt - 2.0) - 3.0*Pi ;
	theta0 = 180.0*acos(numer/denom)/Pi ;
	theta0 = 0.0 ;

	//parasite
  	para[0] = Center[0] ;
  	para[1] = Center[1] ;
  	para[2] = Center[2] ;
	vmax = sqrt(Theta/PRmass) ;
	paraV[0] =  -vmax + 2.0*vmax*RanGen[nrt].Random() ;
	paraV[1] =  -vmax + 2.0*vmax*RanGen[nrt].Random() ;
	paraV[2] =  -vmax + 2.0*vmax*RanGen[nrt].Random() ;
	paradx = 0 ;
	parady = 0 ;
	for ( int i = 0 ; i <  3 ; ++i )
	{
		paraft[i] = 0.0 ;
		parafdt[i] = 0.0 ;
	}

	//We now read links(edges) and triangles..
	double l00 ;
	int mm2 = sprintf(filename,"E%d.dat",Nt) ;
	fpv = fopen(filename,"r" ) ;
	int nol = linecount(fpv) ;
	fseek(fpv, 0, SEEK_SET);
	int l1,l2 ;
	double l0 ;
	double meanedge = 0 ;
	double meanedge2 = 0 ;
	//FILE *fuk = fopen("edges.dat","w") ;
	for ( int i = 0 ; i < nol ; ++i )
	{
		int mmm = fscanf(fpv,"%d\t%d\n",&l1,&l2) ;
		//Add this link to edges vector
		edges.push_back(linkk()) ;
		edges[i].l1 = l1-1 ;
		edges[i].l2 = l2-1 ;
		edges[i].kappap = 0 ;
		edges[i].elap = 0 ;
		double distx =  (vts[edges[i].l1].x - vts[edges[i].l2].x)*(vts[edges[i].l1].x - vts[edges[i].l2].x) ;
		double disty =  (vts[edges[i].l1].y - vts[edges[i].l2].y)*(vts[edges[i].l1].y - vts[edges[i].l2].y) ;
		double distz =  (vts[edges[i].l1].z - vts[edges[i].l2].z)*(vts[edges[i].l1].z - vts[edges[i].l2].z) ;
		edges[i].l0 = sqrt( distx + disty + distz ) ;
		//fprintf(fuk,"%lf\n",edges[i].l0 ) ;
		meanedge += edges[i].l0 ;
		meanedge2 += edges[i].l0 * edges[i].l0 ;
	}
	//fclose(fuk) ;
	cout << "Mean edge is:\t"<< meanedge/(double)nol << endl ;
	avg_edge = meanedge/(double)nol ;
	cout << "Standard deviation is:\t"<< sqrt((meanedge2/(double)nol) - ((meanedge/(double)nol)*(meanedge/(double)nol))) << endl ;
	fclose(fpv) ;  //closing of links file..

	int mm3 = sprintf(filename,"T%d.dat",Nt) ;
	fpv = fopen(filename,"r" ) ;
	nol = linecount(fpv) ;
	fseek(fpv, 0, SEEK_SET);
	int t1,t2,t3 ;
	for ( int i = 0 ; i < nol ; ++i )
	{
			int mmm  = fscanf(fpv,"%d\t%d\t%d\n",&t1,&t2,&t3) ;
			tries.push_back(triangle()) ;
			tries[i].l1 = t1-1 ;
			tries[i].l2 = t2-1 ;
			tries[i].l3 = t3-1 ;
	}
	fclose(fpv) ;
	/////////////////////////////////////////
	//Now we make other list "edgest" for bending energy
	for ( int i = 0 ; i <  edges.size() ; ++i )
	{
		edgest.push_back(linkkk()) ;
	}
	int flg1 = 0 ;
	int flg2 = 0 ;
	for ( int i = 0 ; i <  edges.size() ; ++i )
	{
		int l1 = edges[i].l1 ;
		int l2 = edges[i].l2 ;
		flg1 = 0 ;
		flg2 = 0 ;
		for ( int j = 0 ; j < tries.size() ; ++j )
		{
			flg1 = 0 ;
			flg2 = 0 ;
			flg1 += (l1 == tries[j].l1 ) ? 1 : 0 ;
			flg1 += (l1 == tries[j].l2 ) ? 1 : 0 ;
			flg1 += (l1 == tries[j].l3 ) ? 1 : 0 ;
			flg2 += (l2 == tries[j].l1 ) ? 1 : 0 ;
			flg2 += (l2 == tries[j].l2 ) ? 1 : 0 ;
			flg2 += (l2 == tries[j].l3 ) ? 1 : 0 ;
			if (( flg1 > 0 ) && (flg2 > 0 ) )
			{
				if (( tries[j].l1 != l1 ) && (tries[j].l1 != l2 ) )
				{
					edgest[i].l1 = tries[j].l1 ;
				}
				if (( tries[j].l2 != l1 ) && (tries[j].l2 != l2 ) )
				{
					edgest[i].l1 = tries[j].l2 ;
				}
				if (( tries[j].l3 != l1 ) && (tries[j].l3 != l2 ) )
				{
					edgest[i].l1 = tries[j].l3 ;
				}
				break ;
			}
		}
		for ( int j = tries.size()-1 ; j >= 0 ; j-- )
		{
			flg1 = 0 ;
			flg2 = 0 ;
			flg1 += (l1 == tries[j].l1 ) ? 1 : 0 ;
			flg1 += (l1 == tries[j].l2 ) ? 1 : 0 ;
			flg1 += (l1 == tries[j].l3 ) ? 1 : 0 ;
			flg2 += (l2 == tries[j].l1 ) ? 1 : 0 ;
			flg2 += (l2 == tries[j].l2 ) ? 1 : 0 ;
			flg2 += (l2 == tries[j].l3 ) ? 1 : 0 ;
			if (( flg1 > 0 ) && (flg2 > 0 ) )
			{
				if (( tries[j].l1 != l1 ) && (tries[j].l1 != l2 ) )
				{
					edgest[i].l2 = tries[j].l1 ;
				}
				if (( tries[j].l2 != l1 ) && (tries[j].l2 != l2 ) )
				{
					edgest[i].l2 = tries[j].l2 ;
				}
				if (( tries[j].l3 != l1 ) && (tries[j].l3 != l2 ) )
				{
					edgest[i].l2 = tries[j].l3 ;
				}
				break ;
			}
		}
	}
	//////////////////////////////////////////////
	//spring energy
	TotalVA0() ;  //initial S.area and volume
	cout << "Sarea:\t" << Sarea0 << endl ;
	cout << "Volume:\t" <<  Volume0 << endl ;
	//setting surface area higher to make it more rigid
	Volume0 = Volume0  ;
	kappasa = shear*40.0 ;
	kappad = shear ;
	kappav = 4000.0 ;
	double kpp,elapp ;


	for ( int i = 0 ; i <  edges.size() ; ++i )
	{
		ll0 = edges[i].l0 ;
		edges[i].kappap = set_kp(ll0) ;
		edges[i].elap = set_lp(ll0) ;
	}
	cout << "Kappa_p and p are:\t"<< edges[0].kappap <<"\t"<< edges[0].elap << endl ;

	for ( int i = 0 ; i < tries.size() ; ++i )
	{
		double are = Areatri(tries[i].l1,tries[i].l2,tries[i].l3) ;
		tarea.push_back(are) ;
	}
	////////////////////////////////////////////////
	//Indexing of edgest
	for ( int i = 0 ; i <  edges.size() ; ++i )
	{
		iedge.push_back(int()) ;
	}

	centerofmass() ;
	makeiedge() ;
	TotalVA() ;
	//pfdt() ;
	force_calculate_dt() ;
	//////////////////////////////////////////////////

	tricm() ;
	sortsolvent_initial1() ;
	llistupdate() ;

	/////////////////////////////////////////////////////////
	//knob selection
	knob_size = atoi(argv[3]) ;
	int sri = atoi(argv[2]) ;
	
	// setting knob positions
	Choose_Knob_Locations(knob_size,Center) ;

	int cnt = 0 ;
	for ( int i = 0 ; i < vts.size() ; ++i )
	{
		if ( vts[i].knob == 1 )
		{
			cnt += 1 ;
		}
	}
	cout << "no of vertices selected:\t" << cnt << endl ;
   /////////////////////////////////////////////////////////////////
	//We start iterations here ...
	int filec = 0 ;

	int filecc = 0 ;
	int bbi = 0 ;
	int upd = 0 ;
	double rpos[3] ; //Random vector for grid shifting
	int csteps = (int)(dtcd/dtmd) ;
	int colsteps = csteps-1 ;


	char filename1[256], filename2[256],filename3[256];

	mm1 = sprintf(filename1,"data_%d.dat",sri) ;
	FILE *fdat = fopen(filename1,"w") ;
	mm2 = sprintf(filename2,"fp_%d.dat",sri) ;
	FILE *foot = fopen(filename2,"w") ;
	mm3 = sprintf(filename3,"data_%d.dump",sri) ;
	FILE *fprs = fopen(filename3,"w");

	//Now we create shear flow...

	int uplist = 25 ;

	int peel_id = 0 ;

	for ( int mpcd_steps = 0 ; mpcd_steps < 10000000 ; ++mpcd_steps )
	{
		if ( mpcd_steps < 500000 )
		{
			wallspeed = 0.0 ;
		}
		else
		{
			wallspeed = walls ;
		}
		colsteps += 1 ;
		bbi += 1 ;
		upd += 1 ;
		filecc += 1 ;
		/////////////////////////////////
		if (colsteps == csteps )
		{
			collision_step() ;
			colsteps = 0;
		}
		////////////////////////////////


	#pragma omp  parallel for
	for ( unsigned int i = 0 ; i < vts.size() ; ++i )
	{
		vts[i].vx += (fmapdt[i].fx*dtmd/(2.0*pmass) ) ;
		vts[i].vy += (fmapdt[i].fy*dtmd/(2.0*pmass) ) ;
		vts[i].vz += (fmapdt[i].fz*dtmd/(2.0*pmass) ) ;
		fmapdt[i].fx = 0.0 ;
		fmapdt[i].fy = 0.0 ;
		fmapdt[i].fz = 0.0 ;
	}

	//Membrane part
	//centerofmass() ;
	if ( upd == uplist)
	{
		llistupdate() ;
		tricm() ;
		upd = 0 ;
	}

	//Solvent
	//Streaming + Bounce back .. For Poisellue flow
	//Change Bounce Back boundary conditions
	 #pragma omp parallel
	 {
		 #pragma omp for private(A,B)
		for ( int i = 0 ; i < Ns ; ++i )
		{
			double fullts = reflectionFed(i,0.0) ;
			solvent[i].xpos += solvent[i].vx*(dtmd-fullts) ;
			solvent[i].ypos += solvent[i].vy*(dtmd-fullts) ;
			solvent[i].zpos += solvent[i].vz*(dtmd-fullts) ;
			double ul = Lz ;
			double ll = 0 ;
			if (solvent[i].zpos > ul)
			{
				A[0] = solvent[i].xpos - solvent[i].vx*(dtmd-fullts) ;
				A[1] = solvent[i].ypos - solvent[i].vy*(dtmd-fullts) ;
				A[2] = solvent[i].zpos - solvent[i].vz*(dtmd-fullts) ;
				double s = (ul-A[2])/(solvent[i].zpos-A[2]) ; //parameter of line equation
				B[0] = A[0] + s*(solvent[i].xpos-A[0]) ; // point of crossing
				B[1] = A[1] + s*(solvent[i].ypos-A[1]) ;
				B[2] = A[2] + s*(solvent[i].zpos-A[2]) ;
				double dist = sqrt( ((B[0]-A[0])*(B[0]-A[0])) + ((B[1]-A[1])*(B[1]-A[1])) + ((B[2]-A[2])*(B[2]-A[2])) ) ;
				double dtt = dist/sqrt( (solvent[i].vx*solvent[i].vx) + (solvent[i].vy*solvent[i].vy) + (solvent[i].vz*solvent[i].vz)) ;
				solvent[i].vx *= -1.0 ;
				solvent[i].vy *= -1.0 ;
				solvent[i].vz *= -1.0 ;
				solvent[i].vx += wallspeed ;

				solvent[i].xpos = B[0] + solvent[i].vx*((dtmd-fullts)-dtt) ;
				solvent[i].ypos = B[1] + solvent[i].vy*((dtmd-fullts)-dtt) ;
				solvent[i].zpos = B[2] + solvent[i].vz*((dtmd-fullts)-dtt) ;

			}
			else if(solvent[i].zpos < ll)
			{
				A[0] = solvent[i].xpos - solvent[i].vx*(dtmd-fullts) ;
				A[1] = solvent[i].ypos - solvent[i].vy*(dtmd-fullts) ;
				A[2] = solvent[i].zpos - solvent[i].vz*(dtmd-fullts) ;
				double s = (ll-A[2])/(solvent[i].zpos-A[2]) ; //parameter of line equation
				B[0] = A[0] + s*(solvent[i].xpos-A[0]) ; // point of crossing
				B[1] = A[1] + s*(solvent[i].ypos-A[1]) ;
				B[2] = A[2] + s*(solvent[i].zpos-A[2]) ;
				double dist = sqrt( ((B[0]-A[0])*(B[0]-A[0])) + ((B[1]-A[1])*(B[1]-A[1])) + ((B[2]-A[2])*(B[2]-A[2])) ) ;
				double dtt = dist/sqrt( (solvent[i].vx*solvent[i].vx) + (solvent[i].vy*solvent[i].vy) + (solvent[i].vz*solvent[i].vz)) ;
				solvent[i].vx *= -1.0 ;
				solvent[i].vy *= -1.0 ;
				solvent[i].vz *= -1.0 ;
				solvent[i].xpos = B[0] + solvent[i].vx*((dtmd-fullts)-dtt) ;
				solvent[i].ypos = B[1] + solvent[i].vy*((dtmd-fullts)-dtt) ;
				solvent[i].zpos = B[2] + solvent[i].vz*((dtmd-fullts)-dtt) ;
			}
			else
			{
				solvent[i].vx += 0.0 ;
			}
			solvent[i].xpos = bbound(solvent[i].xpos,Lx) ;
			solvent[i].ypos = bbound(solvent[i].ypos,Ly) ;
		}
	 }



		#pragma omp parallel for
		for ( unsigned int i = 0 ; i < vts.size() ; ++i )
		{
			vts[i].x += (vts[i].vx*dtmd)  ;
			vts[i].y += (vts[i].vy*dtmd)  ;
			vts[i].z += (vts[i].vz*dtmd)  ;
			if ( vts[i].z < 0.0 )
			{
			     double dtreamin =  (vts[i].z-0.0)/vts[i].vz ;
			     vts[i].x -= 2.0*dtreamin*vts[i].vx ;
			     vts[i].y -= 2.0*dtreamin*vts[i].vy ;
			     vts[i].z -= 2.0*dtreamin*vts[i].vz ;
			     vts[i].vx *= -1.0 ;
			     vts[i].vy *= -1.0 ;
			     vts[i].vz *= -1.0 ;
			}
			else if ( vts[i].z > Lz )
			{
			     double dtreamin =  (vts[i].z-Lz)/vts[i].vz ;
				 vts[i].x -= 2.0*dtreamin*vts[i].vx ;
				 vts[i].y -= 2.0*dtreamin*vts[i].vy ;
				 vts[i].z -= 2.0*dtreamin*vts[i].vz ;
				 vts[i].vx *= -1.0 ;
				 vts[i].vy *= -1.0 ;
				 vts[i].vz *= -1.0 ;
			}
			// seperate than solvent function because of updating ix and iy
			vts[i].x = bboundmx(i,vts[i].x,Lx) ;
			vts[i].y = bboundmy(i,vts[i].y,Ly) ;
		}


 	//Computation of force of dt and updating velocity
	TotalVA() ;
	force_calculate_dt() ;
	Exclusion() ;


   // Gravitational force
   #pragma omp parallel for
   for ( int i = 0 ; i < vts.size() ; ++i )
   {
   	fmapdt[i].fz += -Fg/(double)vts.size() ;
   }
   //Bond forces
   int bonds = 0 ;
   double B1[3] ;
   #pragma omp parallel for default(shared) private(A,B1) reduction(+:bonds)
   for ( int i = 0 ; i < vts.size() ; ++i )
   {
	if ( vts[i].bondi == 1 )
	{
			bonds += 1 ;
			int jj = vts[i].bondl ;
			A[0] = vts[i].x - xlig[jj] ;
			A[1] = vts[i].y - ylig[jj];
			A[2] = vts[i].z ;
			A[0] -= Lx*round(A[0]/Lx) ;
			A[1] -= Ly*round(A[1]/Ly) ;
			double aa = sqrt( A[0]*A[0] + A[1]*A[1] + A[2]*A[2] ) ;
			A[0] /= aa ;
			A[1] /= aa ;
			A[2] /= aa ;
			B1[0] = -ks * (aa - el0) * A[0] ;
			B1[1] = -ks * (aa - el0) * A[1] ;
			B1[2] = -ks * (aa - el0) * A[2] ;
			fmapdt[i].fx += B1[0] ;
			fmapdt[i].fy += B1[1] ;
			fmapdt[i].fz += B1[2] ;
		}
	}
	//////////////////////////////////////
	//Preparation to bond dynamics
	#pragma omp parallel for private(A,B)
	for ( unsigned int ii = 0 ; ii < vts.size() ; ++ii )
	{
		if ( vts[ii].bondi == 1)
		{
			int jj = vts[ii].bondl ;
			A[0] = (vts[ii].x ) - (xlig[jj]) ;
			A[1] = (vts[ii].y ) - (ylig[jj]) ;
			A[2] = vts[ii].z - 0.0 ;
			if ( xlig[jj] > Lx-a )
			{
				A[0] += (vts[ii].x < a)?Lx:0 ;
			}
			if ( xlig[jj] < a )
			{
				A[0] -= (vts[ii].x > Lx-a)?Lx:0 ;
			}
			if ( ylig[jj] > Ly-a )
			{
				A[1] += (vts[ii].y < a)?Ly:0 ;
			}
			if ( ylig[jj] < a )
			{
				A[1] -= (vts[ii].y > Ly-a)?Ly:0 ;
			}
			double xx = sqrt ( A[0]*A[0] + A[1]*A[1] + A[2]*A[2] ) ;
			vts[ii].forc = ks*fabs(xx-el0)/FD ;
			double koff = offrate*exp(ks*fabs(xx-el0)/FD) ;
			vts[ii].poff = 1.0 - exp(-koff*dtmd) ;
			if ( xx > 1.0*dthr ) { vts[ii].poff = 1.0 ; }
		}

		else if ( ( vts[ii].z < dthr ) && (vts[ii].bondi == 0) && (vts[ii].knob == 1) )
		{
			vts[ii].Pairs.clear() ;
			int pid = vts[ii].pid ;
			double kon = onrate ;
			vts[ii].pon = 1.0 - exp(-kon*dtmd) ;
			for ( unsigned int jkl = 0 ; jkl < Cells[pid].liglist.size() ; ++jkl )
			{
					int j = Cells[pid].liglist[jkl] ;
					Pair apair ;
					A[0] = (vts[ii].x ) - (xlig[j]) ;
					A[1] = (vts[ii].y ) - (ylig[j]) ;
					A[2] = vts[ii].z - 0.0 ;
					if ( xlig[j] > Lx-a )
					{
						A[0] += (vts[ii].x < a)?Lx:0 ;
					}
					if ( xlig[j] < a )
					{
						A[0] -= (vts[ii].x > Lx-a)?Lx:0 ;
					}
					if ( ylig[j] > Ly-a )
					{
						A[1] += (vts[ii].y < a)?Ly:0 ;
					}
					if ( ylig[j] < a )
					{
						A[1] -= (vts[ii].y > Ly-a)?Ly:0 ;
					}

					double aaa =  ( A[0]*A[0] + A[1]*A[1] + A[2]*A[2] ) ;
					if (aaa <= dthr*dthr )
					{
						apair.first = aaa ;
						apair.second = j ;
						vts[ii].Pairs.push_back(apair) ;
					}
			 }
			 sort(vts[ii].Pairs.begin(),vts[ii].Pairs.end()) ;
		}
	}
	//Bond dynamics
	bond_dynamics(mpcd_steps,foot) ;


	 #pragma omp parallel for
	for ( int i = 0 ; i < vts.size() ; ++i )
	{
		vts[i].vx += (fmapdt[i].fx*dtmd/(2.0*pmass) ) ;
		vts[i].vy += (fmapdt[i].fy*dtmd/(2.0*pmass) ) ;
		vts[i].vz += (fmapdt[i].fz*dtmd/(2.0*pmass) ) ;
	}



	/////////////////////////////////////////////////////////
	if (( filecc == 15000 )||(mpcd_steps == 0 ) )
	{
		filec += 1 ;
		cout << mpcd_steps*dtmd << "\t" << Sarea << "\t" << Volume << "\t" << bonds << endl ;
		fprintf(fdat,"#%d\n",filec) ;
		for ( int i = 0 ; i < Nt ; ++i )
		{
			fprintf(fdat,"%lf\t%lf\t%lf\t%d\t%d\t%d\n",vts[i].x,vts[i].y,vts[i].z,vts[i].ix,vts[i].iy,vts[i].bondi) ;
		}
		fflush(fdat) ;
	  	write_lammpsfile(fprs,filec) ;
		fflush(fprs) ;
		if ( bonds == 0 )
		{
			peel_id += 1 ;
		}
		if ( bonds != 0 )
		{
			peel_id = 0 ;
		}
		filecc = 0 ;
	}
	if ( peel_id == 100 )
	{
		break ;
	}
}

fclose(fprs) ;
fclose(foot) ;
fclose(fdat) ;





return 0 ;

}
