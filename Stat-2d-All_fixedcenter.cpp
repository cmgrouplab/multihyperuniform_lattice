//CPP file for obtaining statistics of the packing, including:
//     (1) pair correlation functions
//     (2) number variance
//     (3) structure factor  


//version: 04/10/2020
//author: duyu chen

//modified by Yang Jiao on 07.20.22 
//implement node and edge-based sampling to study effects of correlation on variance scaling

//modified 10/14/22
//read in a hyperuniform system and apply random but tiny pertrubtations ...

using namespace std;

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
//#include <random>

//#define N 2500 
#define N 167
//#define N 421
#define SD 2
#define Bin 1.5
//define Bin 0.008
double Center[N][SD];
double Lx;  //side length of box, these need to be consistent with the chocie of bin and bin_num for both real and Fourier space statistics 
double Ly;
#define Nbin 80 //for g(r)

double G[N][Nbin];//this is the local un-normalized g2, to compute node-based variance 
double BG[3*N][Nbin]; //this is for bond-based statistics 


//need to check the change these numbers for SK calculation
int Nk = 25; //this is the number of mesh point along each direction in K space, considering the symmetry the unit of mesh is 2pi/L
double Kbin = 0.025; //this should be on the magnitude of 2pi/L , need to know an estimate of L
#define Kbin_num 60 //this is for the radial average, roughly 2 times of Nk 

/*
int Nk = 300;
double Kbin = 0.15;
#define Kbin_num 100
*/

double pi = 3.14159265358979;
#define Rbin 2  //for number variance
#define Nr 10000 //sampling points for number variance at each r
#define MAXY 20000 //the maximal number is 32767 

ofstream fout;
ifstream fin;


void read_config()
{
	FILE * fp;
//if ((fp = fopen("../../10000/0.04/test2/vertex_graphene_process.txt", "r")) == NULL)
	if ((fp = fopen("vertex_relax.txt", "r")) == NULL)
	{
		perror("Cannot open file!\n");
		exit(1);
	}
	else
	{
		double temp_t;
		fscanf(fp, "%lf", &temp_t);
		fscanf(fp, "%lf", &temp_t);
		fscanf(fp, "%lf", &Lx);
		fscanf(fp, "%lf", &temp_t);
		fscanf(fp, "%lf", &temp_t);
		fscanf(fp, "%lf", &Ly);
		cout << "Lx: " << Lx << endl;
		cout << "Ly: " << Ly << endl;
		//L = 1029.0 / L;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < SD; j++)
			{
				fscanf(fp, "%lf", &Center[i][j]);
			}
		//	for(int i = 0; i < N; i++)
		//		std::cout << Center[i][0] << " " << Center[i][1] << std::endl;
		//Nbin = (int)floor(Ly / 2.0 / Bin);
		fclose(fp);
	}
	
    fp = fopen("centers.xls","w");
	//fprintf(fp, "%lf\n", L);
	for(int i = 0; i < N; i++)
		{
		for(int j = 0; j < SD; j++)
			fprintf(fp, "%lf\t", Center[i][j]);
		fprintf(fp, "\n");
		}	
	fclose(fp);

}

/*
void print_config()
{
	cout<<"printing out the pertrubed configuration ..."<<endl;
	
	fout.open("config_perturbed.txt");
	
	fout<<N<<endl;
	fout<<Lx<<"\t"<<0<<endl;
	fout<<"0"<<"\t"<<Ly<<endl;
	
	for(int i=0; i<N; i++)
		fout<<Center[i][0]<<"\t"<<Center[i][1]<<endl;
		
	fout.close();
	
	
}
*/

/*
void perturb_config()
{
	//apply very tiny random perturbation to the hyperunfiorm configurations...
	
	double rand_mod = 0.0005; //this is with respect to the box length Lx
	
	cout<<"rand_mod = "<<rand_mod<<endl;
	
	cout<<"petrubing the configuration ..."<<endl;
	
	double dx, dy;
	
	for(int i=0; i<N; i++)
	{
		dx = (((double)(rand()%100)/100.0)-0.5)*rand_mod*Lx;
		dy = (((double)(rand()%100)/100.0)-0.5)*rand_mod*Ly;
		
		Center[i][0] = Center[i][0] + dx;
		if(Center[i][0]>=Lx) Center[i][0] = Center[i][0] - Lx;
		else if(Center[i][0]<0) Center[i][0] = Center[i][0] + Lx;
		
		Center[i][1] = Center[i][1] + dy;
		if(Center[i][1]>=Ly) Center[i][1] = Center[i][1] - Ly;
		else if(Center[i][1]<0) Center[i][1] = Center[i][1] + Ly;
	}
}
*/


double MinDis(int m, int n)
{
  //find the minimal distance between the centers of two polyhedra in Ecludean space...
  //by checking all the images of Poly[m], while keeping Poly[n] in the central box
  //record the index of the box which Poly[m] is in...
  //the order (m,n) is important, especially when getting the NNL.... 
  double dx = Center[m][0] - Center[n][0];
  double dy = Center[m][1] - Center[n][1];
 

  double dist = 1000000000.0; //just a large number...


  //loop over all possible images of Point m, keep n fixed in the center simulation box....
  for(int i = -1; i <= 1; i++)
    for(int j = -1; j <= 1; j++)
	{
	  double tempd[SD]; 
	  tempd[0] = dx + (double)i*Lx;
	  tempd[1] = dy + (double)j*Ly;
	  
	  double tempdist = tempd[0]*tempd[0]+tempd[1]*tempd[1]; 
	  
	  //printf("TempDist[%d][%d][%d] = %f\n", i, j, k, tempdist);
	  
	  if(tempdist < dist) // store the smallest distance...
	    {
			dist = tempdist; 
	    }
	  
	}
  
  
  return sqrt(dist); //this is center-to-center distance...
}

double MinDis(double x1, double y1, double x2, double y2)
{
  //this is a simplified version, for a rectanglar box, no need to check all images, just the orthogonal minimal distance should be OK
  
  double dx = fabs(x1 - x2);
  double dy = fabs(y1 - y2);
  
  if(dx >= Lx/2.0) dx = Lx - dx;
  if(dy >= Ly/2.0) dy = Ly - dy;
  	
  return sqrt(dx*dx + dy*dy); 
  
 
 
  /*
  double dx = (x1 - x2);
  double dy = (y1 - y2);
 
  double dist = 1000000000.0; //just a large number...


  //loop over all possible images of Point m, keep n fixed in the center simulation box....
  for(int i = -1; i <= 1; i++)
    for(int j = -1; j <= 1; j++)
	{
	  double tempd[SD]; 
	  tempd[0] = dx + (double)i*Lx;
	  tempd[1] = dy + (double)j*Ly;
	  
	  double tempdist = tempd[0]*tempd[0]+tempd[1]*tempd[1]; 
	  
	  //printf("TempDist[%d][%d][%d] = %f\n", i, j, k, tempdist);
	  
	  if(tempdist < dist) // store the smallest distance...
	    {
			dist = tempdist; 
	    }
	  
	}
  
  
  return sqrt(dist); //this is center-to-center distance...
  */
  
}



double Get_NDensity() // the number density, for computing g2...
{
  double VLambda;
  VLambda = Lx * Ly;

  return (double)N/VLambda;
}


void Get_PairCorr()
{

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //compute the minimum length of the lattice vectors
  double * g;
  
  g = new double [Nbin];
  for(int i = 0; i < Nbin; i++)
  {
		g[i] = 0.0;
   }	
  printf("Computing G2 and g2 now....\n");

  
  //loop over all particles
  for(int i = 0; i < N; i++)
      for(int j = 0; j < N; j++)
	{
	  
	  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  //for the radial g2....
	  double temp_dis = MinDis(j, i)/Bin;
	  int int_dis = (int)floor(temp_dis);
	  //the upper limit should be NLcounter/2.0, here relax this to MaxNL
	  // && temp_dis<(MaxNL-1)
	  if(j != i && temp_dis < Nbin)
	    {
	      
	      //double delta_dis = temp_dis - int_dis;
	      //if(delta_dis>0.5) int_dis++;
	      
	      g[int_dis] = g[int_dis] + 1.0;
	    }
	  }

	 
  double rho_n = Get_NDensity();

  for(int r = 0; r < Nbin; r++)
    {
      g[r] = g[r] / (N * rho_n * pi * ((r+1.0) * (r+1.0) - r * r) * Bin * Bin);
    }

 
  
  FILE* fp = fopen("./g2_relax.xls","w");
  for(int r = 0; r < Nbin; r++)
    fprintf(fp, "%lf\t%lf\n", (r + 0.5) * Bin, g[r]);
  fclose(fp);
  delete [] g;
}

/*
//using similar method used by Reka and Eli for the real space system
void num_var_node()
{
	cout<<"computing node-based variance now ...."<<endl;
	
	for(int i=0; i<N; i++)
	 	for(int j=0; j<Nbin; j++)
	 	{
	 		G[i][j] = 0;
		 }
	
	
	//loop over all particles
  for(int i = 0; i < N; i++)
      for(int j = 0; j < N; j++)
	{
	  
	  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  //for the radial g2....
	  double temp_dis = MinDis(j, i)/Bin;
	  int int_dis = (int)floor(temp_dis);
	  //the upper limit should be NLcounter/2.0, here relax this to MaxNL
	  // && temp_dis<(MaxNL-1)
	  if(j != i && temp_dis < Nbin)
	    {
	      
	      //double delta_dis = temp_dis - int_dis;
	      //if(delta_dis>0.5) int_dis++;
	     
		  for(int k=int_dis; k<Nbin; k++)
	         G[i][k] = G[i][k] + 1.0;
	    }
	  }
	  
	  
	//now need to compute the average of variance
	
	double temp_sum[Nbin];
	double temp_var2[Nbin];
	double temp_ave[Nbin];
	
	for(int i=0; i<Nbin; i++)
		{
			temp_sum[i] = 0;
			temp_var2[i] = 0;
			temp_ave[i] = 0;
		}
	
		
	for(int i=0; i<Nbin; i++)
		for(int j=0; j<N; j++)
			temp_sum[i] += G[j][i];
			
	//compute the average		
	for(int i=0; i<Nbin; i++)
		temp_ave[i] = temp_sum[i]/(double)N;
		
	//now compute the variance
	for(int i=0; i<Nbin; i++)
	{
		for(int j=0; j<N; j++)
		{
			temp_var2[i] = temp_var2[i]+(G[j][i]-temp_ave[i])*(G[j][i]-temp_ave[i]);
		}
		
		temp_var2[i] = temp_var2[i]/(double)N;
	}
	
	fout.open("var2_node.xls");
	for(int i=0; i<Nbin; i++)
		fout<<Bin*i<<'\t'<<temp_var2[i]<<endl;
	fout.close();
	
	fout.open("Nave_node.xls");
	for(int i=0; i<Nbin; i++)
		fout<<Bin*i<<'\t'<<temp_ave[i]<<endl;
	fout.close();
	
	//int a;
	//cin>>a;
}
*/


/*
//using similar method used by Reka and Eli for the real space system
void num_var_point()
{
	cout<<"computing point-based variance now ...."<<endl;
	
	for(int i=0; i<N; i++)
	 	for(int j=0; j<Nbin; j++)
	 	{
	 		G[i][j] = 0;
		 }
	
	
	//i is the index for the random sampling point
	//j is the particle center on honeycomb
  for(int i = 0; i < N; i++)
  {
  
     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  //first generate the coords for the random points ...
	  double temp_x = ((double)(rand()%20000)/20000.0)*Lx;
	  double temp_y = ((double)(rand()%20000)/20000.0)*Ly;
   
      for(int j = 0; j < N; j++)
	{
	  
	  
	  double temp_dis = MinDis(temp_x, temp_y, Center[j][0], Center[j][1])/Bin;
	  int int_dis = (int)floor(temp_dis);
	  //the upper limit should be NLcounter/2.0, here relax this to MaxNL
	  // && temp_dis<(MaxNL-1)
	  if(temp_dis < Nbin)
	    {
	      
	      //double delta_dis = temp_dis - int_dis;
	      //if(delta_dis>0.5) int_dis++;
	     
		  for(int k=int_dis; k<Nbin; k++)
	         G[i][k] = G[i][k] + 1.0;
	    }
	  }
    }
	  
	//now need to compute the average of variance
	
	double temp_sum[Nbin];
	double temp_var2[Nbin];
	double temp_ave[Nbin];
	
	for(int i=0; i<Nbin; i++)
		{
			temp_sum[i] = 0;
			temp_var2[i] = 0;
			temp_ave[i] = 0;
		}
	
		
	for(int i=0; i<Nbin; i++)
		for(int j=0; j<N; j++)
			temp_sum[i] += G[j][i];
			
	//compute the average		
	for(int i=0; i<Nbin; i++)
		temp_ave[i] = temp_sum[i]/(double)N;
		
	//now compute the variance
	for(int i=0; i<Nbin; i++)
	{
		for(int j=0; j<N; j++)
		{
			temp_var2[i] = temp_var2[i]+(G[j][i]-temp_ave[i])*(G[j][i]-temp_ave[i]);
			
			//temp_var2[i] = temp_var2[i]+(G[j][i]*G[j][i]-temp_ave[i]*temp_ave[i]);
		}
		
		temp_var2[i] = temp_var2[i]/(double)N;
	}
	
	fout.open("var2_point.xls");
	for(int i=0; i<Nbin; i++)
		fout<<Bin*i<<'\t'<<temp_var2[i]<<endl;
	fout.close();
	
	fout.open("Nave_point.xls");
	for(int i=0; i<Nbin; i++)
		fout<<Bin*i<<'\t'<<temp_ave[i]<<endl;
	fout.close();
	
	//int a;
	//cin>>a;
}
*/

/*
//this use the full bond configuration
void num_var_bond_II()
{
	cout<<"computing bond-based variance now ...."<<endl;
	
	for(int i=0; i<3*N; i++)
	 	for(int j=0; j<Nbin; j++)
	 	{
	 		BG[i][j] = 0;
		 }
	
	//still using A*N samples, read in the connectivity matrix, using the first A*N bonds for sampling .... (where A is the coordination number)

	fin.open("connectivity_matrix.txt");

	int temp_int; //temp int for reading in the connectivity matrix
	
	int temp_ct = 0; //this is the counter, the maximum is N
	
	cout<<"checking the connectivity matrix now, this can be slow ...."<<endl;
	
//m, n are the index for the points in the lattice, for connectivity matrix
  for(int m = 0; m < N; m++)
    for(int n=0; n<N; n++)
  {
      fin>>temp_int;
      if(temp_int == 1) //this is a bond...
      {
	     cout<<"for bound "<<m<<" "<<n<<endl;
  
     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  //first generate the coords for the random points ...
	  double ratio = ((double)(rand()%20000)/20000.0);
	  //double temp_y = ((double)(rand()%20000)/20000.0)*Ly;
   	 double temp_x = Center[m][0]*ratio + Center[n][0]*(1-ratio);
   	 double temp_y = Center[m][1]*ratio + Center[n][1]*(1-ratio);
   
   //now compute the distance to all lattice sites
      for(int j = 0; j < N; j++)
    	{
	  
	  
	  double temp_dis = MinDis(temp_x, temp_y, Center[j][0], Center[j][1])/Bin;
	  int int_dis = (int)floor(temp_dis);
	  //the upper limit should be NLcounter/2.0, here relax this to MaxNL
	  // && temp_dis<(MaxNL-1)
	  if(temp_dis < Nbin)
	    {
	      
	      //double delta_dis = temp_dis - int_dis;
	      //if(delta_dis>0.5) int_dis++;
	     
		  for(int k=int_dis; k<Nbin; k++)
	         BG[temp_ct][k] = BG[temp_ct][k] + 1.0;
	    }
	   }
	   
	   temp_ct ++;
	   
      }
      
      if(temp_ct==3*N) goto L1;
    }
	  
	//now need to compute the average of variance
	
	L1: fin.close();
	
	double temp_sum[Nbin];
	double temp_var2[Nbin];
	double temp_ave[Nbin];
	
	for(int i=0; i<Nbin; i++)
		{
			temp_sum[i] = 0;
			temp_var2[i] = 0;
			temp_ave[i] = 0;
		}
	
		
	for(int i=0; i<Nbin; i++)
		for(int j=0; j<3*N; j++)
			temp_sum[i] += BG[j][i];
			
	//compute the average		
	for(int i=0; i<Nbin; i++)
		temp_ave[i] = temp_sum[i]/(double)(3*N);
		
	//now compute the variance
	for(int i=0; i<Nbin; i++)
	{
		for(int j=0; j<3*N; j++)
		{
			temp_var2[i] = temp_var2[i]+(BG[j][i]-temp_ave[i])*(BG[j][i]-temp_ave[i]);
			
			//temp_var2[i] = temp_var2[i]+(G[j][i]*G[j][i]-temp_ave[i]*temp_ave[i]);
		}
		
		temp_var2[i] = temp_var2[i]/(double)(3*N);
	}
	

	
	fout.open("var2_bond.xls");
	for(int i=0; i<Nbin; i++)
		fout<<Bin*i<<'\t'<<temp_var2[i]<<endl;
	fout.close();
	
	fout.open("Nave_bond.xls");
	for(int i=0; i<Nbin; i++)
		fout<<Bin*i<<'\t'<<temp_ave[i]<<endl;
	fout.close();
	
	//int a;
	//cin>>a;
}
*/

/*
//this version only uses first N bonds
void num_var_bond()
{
	cout<<"computing bond-based variance now ...."<<endl;
	
	for(int i=0; i<N; i++)
	 	for(int j=0; j<Nbin; j++)
	 	{
	 		G[i][j] = 0;
		 }
	
	//still using A*N samples, read in the connectivity matrix, using the first A*N bonds for sampling .... (where A is the coordination number)

	fin.open("connectivity_matrix.txt");

	int temp_int; //temp int for reading in the connectivity matrix
	
	int temp_ct = 0; //this is the counter, the maximum is N
	
	cout<<"checking the connectivity matrix now, this can be slow ...."<<endl;
	
//m, n are the index for the points in the lattice, for connectivity matrix
  for(int m = 0; m < N; m++)
    for(int n=0; n<N; n++)
  {
      fin>>temp_int;
      if(temp_int == 1) //this is a bond...
      {
	     cout<<"for bound "<<m<<" "<<n<<endl;
  
     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  //first generate the coords for the random points ...
	  double ratio = ((double)(rand()%20000)/20000.0);
	  //double temp_y = ((double)(rand()%20000)/20000.0)*Ly;
   	 double temp_x = Center[m][0]*ratio + Center[n][0]*(1-ratio);
   	 double temp_y = Center[m][1]*ratio + Center[n][1]*(1-ratio);
   
   //now compute the distance to all lattice sites
      for(int j = 0; j < N; j++)
    	{
	  
	  
	  double temp_dis = MinDis(temp_x, temp_y, Center[j][0], Center[j][1])/Bin;
	  int int_dis = (int)floor(temp_dis);
	  //the upper limit should be NLcounter/2.0, here relax this to MaxNL
	  // && temp_dis<(MaxNL-1)
	  if(temp_dis < Nbin)
	    {
	      
	      //double delta_dis = temp_dis - int_dis;
	      //if(delta_dis>0.5) int_dis++;
	     
		  for(int k=int_dis; k<Nbin; k++)
	         G[temp_ct][k] = G[temp_ct][k] + 1.0;
	    }
	   }
	   
	   temp_ct ++;
	   
      }
      
      if(temp_ct==N) goto L1;
    }
	  
	//now need to compute the average of variance
	
	L1: fin.close();
	
	double temp_sum[Nbin];
	double temp_var2[Nbin];
	double temp_ave[Nbin];
	
	for(int i=0; i<Nbin; i++)
		{
			temp_sum[i] = 0;
			temp_var2[i] = 0;
			temp_ave[i] = 0;
		}
	
		
	for(int i=0; i<Nbin; i++)
		for(int j=0; j<N; j++)
			temp_sum[i] += G[j][i];
			
	//compute the average		
	for(int i=0; i<Nbin; i++)
		temp_ave[i] = temp_sum[i]/(double)(N);
		
	//now compute the variance
	for(int i=0; i<Nbin; i++)
	{
		for(int j=0; j<N; j++)
		{
			temp_var2[i] = temp_var2[i]+(G[j][i]-temp_ave[i])*(G[j][i]-temp_ave[i]);
			
			//temp_var2[i] = temp_var2[i]+(G[j][i]*G[j][i]-temp_ave[i]*temp_ave[i]);
		}
		
		temp_var2[i] = temp_var2[i]/(double)(N);
	}
	

	
	fout.open("var2_bond.xls");
	for(int i=0; i<Nbin; i++)
		fout<<Bin*i<<'\t'<<temp_var2[i]<<endl;
	fout.close();
	
	fout.open("Nave_bond.xls");
	for(int i=0; i<Nbin; i++)
		fout<<Bin*i<<'\t'<<temp_ave[i]<<endl;
	fout.close();
	
	//int a;
	//cin>>a;
}
*/

void num_var()
{
	printf("Computing number variance now....\n");
	int Ns = (int)floor(Ly / 4.0 / Rbin) - 1;
	cout << "Ns = " << Ns << endl;
	double* sigma = new double [Ns];
	double ave, square_ave;
	int Nc;
	double r = Rbin;
	
	/*
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<double> distribution(0.0, 1.0);
	*/
	
	for(int t = 0; t < Ns; t++)
	{
		cout << "r = " << r << endl;
		ave = 0.0;
		square_ave = 0.0;
		for(int n = 0; n < Nr; n++)
		{
			Nc = 0;
			//double cx = distribution(gen) * Lx;
			//double cy = distribution(gen) * Ly;
			//cout << cx << " " << cy << endl;
			double cx = (double)(rand() % MAXY) / (double)MAXY * Lx;
			double cy = (double)(rand() % MAXY) / (double)MAXY * Ly;
			
			for(int i = 0; i < N; i++)
			{
				double dx = fabs(Center[i][0] - cx);
				if (dx > Lx / 2)
				{
					dx = Lx - dx;
				}
				double dy = fabs(Center[i][1] - cy);
				if (dy > Ly / 2)
				{
					dy = Ly - dy;
				}
				double dist = sqrt(dx * dx + dy * dy);
				if(dist < r)
				{
					Nc++; 
				}
			}
//			cout << "Nc = " << Nc << endl;
			ave = ave + (double)Nc;
			square_ave = square_ave + (double)(Nc * Nc);
		}
		ave = ave / Nr;
		square_ave = square_ave / Nr;
		cout << "ave = " << ave << endl;
		cout << "square_ave = " << square_ave << endl;
		sigma[t] = square_ave - ave * ave;
		r = r + Rbin;
	}
	FILE* fp = fopen("num_var_relax.xls","w");
	for(int r = 0; r < Ns; r++)
		fprintf(fp, "%lf\t%lf\n", (r + 0.5)*Rbin, sigma[r]);
	fclose(fp);
	delete [] sigma;
}


double GetInnerProduct(double Vector1[SD], double Vector2[SD])
{
	double sum = 0;

	for(int i = 0; i < SD; i++)
		sum += Vector1[i] * Vector2[i];
	return sum;
}

double Get_Sk(double Vector_k[SD])
{
	double Sk;
	double sum_cos = 0.0;
	double sum_sin = 0.0;
	double temp;
	for(int i = 0; i < N; i++)
	{
			temp = GetInnerProduct(Center[i], Vector_k); 
			sum_cos += cos(temp);
			sum_sin += sin(temp);
	}
	Sk = (sum_cos * sum_cos + sum_sin * sum_sin) / (double)N;
	return Sk;
}

void Print_Sk(double K_Histo[Kbin_num], int K_Counter[Kbin_num])
{
	ofstream histo_out;
	histo_out.open("./Sk_relax.xls");
	for(int t = 0; t < Kbin_num; t++)
		if(K_Counter[t] > 0)
			histo_out << (t + 0.5) * Kbin << "\t" << K_Histo[t] << endl;
			//histo_out << (t + 0.5) * Kbin * mean_nnb_dist / 2.0 / pi << " " << K_Histo[t] << endl;
	histo_out.close();
}

void Get_KHistogram()
{
	printf("Computing Sk now....\n");
	double Sk, k_dis;
	int t;
	double KPoint[SD];
	double K_Histo[Kbin_num];
	int K_Counter[Kbin_num];
	for(int i = 0; i < Kbin_num; i++)
	{
		K_Histo[i] = 0.0;
		K_Counter[i] = 0;
	}
	for(int i = 1; i <= Nk; i++)
		for(int j = - Nk; j <= Nk; j++)
		{
			KPoint[0] = i * 2.0 * pi / Lx;
			KPoint[1] = j * 2.0 * pi / Ly;
			Sk = Get_Sk(KPoint);
			k_dis = 0.0;
			for(int d = 0; d < SD; d++)
				k_dis += KPoint[d] * KPoint[d];
			k_dis = sqrt(k_dis);
			t = floor(k_dis / Kbin);
			if(t < Kbin_num)
			{
				K_Histo[t] += Sk;
				K_Counter[t] ++;
			}
/*
			if(t == 0)
			{
				cout << "Sk = " << Sk << endl;
				cout << KPoint[0] << " " << KPoint[1] << endl;
			}
*/
		}
	for(int i = 1; i <= Nk; i++)
	{
		KPoint[0] = 0.0;
		KPoint[1] = i * 2.0 * pi / Ly;
		Sk = Get_Sk(KPoint);
		k_dis = 0.0;
		k_dis = KPoint[1];
		t = floor(k_dis / Kbin);
		if(t < Kbin_num)
		{
			K_Histo[t] += Sk;
			K_Counter[t] ++;
		}
/*		if(t == 0)
			{
				cout << "Sk = " << Sk << endl;
				cout << KPoint[0] << " " << KPoint[1] << endl;
			}

*/
	}
	for(t = 0; t < Kbin_num; t++)
		if(K_Counter[t] != 0)		
			K_Histo[t] = K_Histo[t] / K_Counter[t];
	//cout << "first bin: " << K_Counter[0] << " " << K_Histo[0] << endl;

	Print_Sk(K_Histo, K_Counter);
/*	
	for(int i = -Nk; i <= Nk; i++)
		for(int j = -Nk; j <= Nk; j++)
		{
			KPoint[0] = i * Kspace;
			KPoint[1] = j * Kspace;
			Sk_2d[i+Nk][j+Nk] = Get_Sk(KPoint);

		}
	FILE * fp = fopen("sk_2d.txt", "w");
	for(int i = 0; i < 2 * Nk + 1; i++)
	{
		for(int j = 0; j < 2 * Nk + 1; j++)
			fprintf(fp, "%lf\t", Sk_2d[i][j]);
		fprintf(fp, "\n");
	}
	fclose(fp);
*/
}



int main()
{
	read_config();
	
	Get_PairCorr();
	Get_KHistogram();
	//srand(time(NULL));
	
	//perturb_config();
	
	//print_config();
	
	num_var();
	
    //num_var_node();
    
    //num_var_point();
    
    //num_var_bond();
    
    //num_var();

	return 1;
}

