/* 
	Description: The ANSI C code of the ECO approach
	Programmer:  Rafael Stubs Parpinelli
	E-mail:      rsparpin@gmail.com
	Date:	     20/11/2012
	Lisence:     Free
	Note:        The system was developed using Linux. If you use this code, pls cite one of the following works:
	             1. PARPINELLI, R.S.; LOPES, H.S. "Biological Plausibility in Optimization: An Ecosystemic View". 
			International Journal of Bio-Inspired Computation (Online), v. 4, p. 345-358, 2012.
		     2. PARPINELLI, R.S.; Lopes, H.S. "A Hierarchical Clustering Strategy to Improve the Biological Plausibility 
			of an Ecology-based Evolutionary Algorithm". In: Ibero-American Conference on Artificial Intelligence
			(IBERAMIA), 2012, Cartajena. Lecture Notes in Computer Science/Lecture Notes in Artificial Intelligence. 
			Berlin, Heidelberg: Springer-Verlag, 2012. v. 7637. p. 310-319.
	To compile:  1. Go to ECO directory
		     2. Type: make

	To debug using gdb: gcc mersenne.o -o eco-thread eco-thread.c -lm -lpthread -g -O0 -v -da -Q

	To run: ./eco-cluster input_eco.in
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include"mersenne.h"
#include <pthread.h>
#include <stdbool.h>
#include "./src-clstr/cluster.h" /* The C Clustering Library */
#include "data.h" //data information for shifted functions

#define abss(a)     (a<0 ? (-a) : a)

#define FAIL    0

//for 3D-AB
#define SEQ_MAX_LENGTH 100  //Protein max length
#define PI 3.14159265
char SEQUENCE[SEQ_MAX_LENGTH]; //protein sequence (max length = SEQ_MAX_LENGTH)
int ProteinSize;
typedef struct 
{
	int Tipo; 		//tipo A (1) ou B (-1)
	double x, y, z; //posição no plano cartesiano
}amino;
amino *Aminoacido;

//This variable sets up the report to files containing: 
// .the convergence of each specie
// .the centroid of all species in each eco-step
// .the whole population in each eco-step, including the initial population
int REPORT;
//Chooses the search strategy
int STRATEGY;

int cur_POP_SIZE[1]; //current value of POP_SIZE.

//This variable disable the habitat feature of the ecossystem meaning that the species evolve without interacting to each other.
int NOHABITAT;

//Control parameters for eco
int POP_SIZE; //INITIAL number of candidate solutions in each specie. For ABC, POP_SIZE represents the number of food sources. 
int ECO_STEP; //number of ecological succession steps
int EVO_STEP; //number of function evaluations in each ECO_STEP.
		      //For ABC: The number of cycles for foraging {an evolutionary period stopping criteria}.
int SP_NUMBER; //number of species. 

int RUN;	//Number of times that eco runs. Gives the Avg and StdDev.
int T_SIZE;     //tourney size for solution selection
//Problem definitions
int DIM;    //number of problem variables
float lb; //lower bound of the variables
float ub;  //upper bound of the variables
int FUNCTION;
int STRATEGY_NUMBER; //is the number of different strategies being employed.

//Single-link CLUSTER
//Variables used to linearly scalonate the distances present in the single-link tree. These distances are used as probabilities to build the clusters.
//This gives more biological plausibility to the system.
#define MINVAL 0.01 //1% of chance to NOT link two populations. Gives 1% of chance to the closest populations not belong to the same habitat. 
#define MAXVAL 0.99 //1% of chance to link. This gives 1% of chance to the farest populations belong to the same habitat.

//ABC Control parameters
#define limit 100  /*A food source which could not be improved through "limit" trials is abandoned by its employed bee*/
//Total number of function evaluations for each specie: ECO_STEP * EVO_STEP * (POP_SIZE*2)
//For ABC, employed bees perform POP_SIZE evaluations and onlooker bees perform another POP_SIZE evaluations.

//DE Control parameters
#define F 0.9             //---weight factor----------------------/
#define CR 1.            //---crossing over factor---------------/

//PSO Control Paramaters
#define W 1 / ( 2 * log( 2 ) ) // 0.721 <Inertia weight>
#define C 0.5 + log( 2 ) // 1.193 <confidence coef>

//DE-BBO Control Parameters
#define URAND			((double)rand()/((double)RAND_MAX + 1.0))
#define INF				1e99

//In ProblemFunc.c
double objfunc(double sol[], int sp_index);
void prepararObjFunc(void);
void Rotation(double x,double y,double z,double a,double b,double c,double u,double v,double w,double t,double *rx,double *ry,double *rz);//3D-AB
double TempValue(double x,int a,int k,int m);
double sqr( double x );

/*Write your own objective function name*/
double (*function)(double sol[], int sp_index) = objfunc; //<<<<<<<<<<<-------
//==========================

//Global variables for ECO
//Dynamically allocated
struct specie { 
	double **pop; //pop[POP_SIZE*2][DIM];  //population of candidate solutions
	double *fo;  //fo[POP_SIZE*2];      //objective function value
	double *fit; //fit[POP_SIZE*2];     //fitness value
	double *best; //best[DIM];           //best solution found
	double bestfo[1];         //best fo value
	int best_index[1];     //index for the best solution

	double *Ci; //Ci[DIM];             //represents where species are concentrated. The centroid.
	double CiFO[1];		  //represents the FO value for the centroid.
	int * trial;

	// only for jDE and jDE-DE
	double *CRR; //CRR[POP_SIZE]
	double *FF; //FF[POP_SIZE]

}*Eco; //Eco[SP_NUMBER];

struct habitat {
	int h_sp_count;		//Contains how many species are in each habitat.
	int *h_sp;    //This vector contains which species are in each habitat. h_sp[SP_NUMBER]
} *H; 		    //struct that contains the habitats H[SP_NUMBER]
double *mean_euclid_dist; //mean_euclid_dist[ECO_STEP]; //contain the mean euclidean distance of each ecological succession step.
double *ElapsedTime; //armazena o tempo de execucao de cada rodada. ElapsedTime[RUN]
int *tourney_index; //keeps the selected indexes for tournament. tourney_index[T_SIZE]
int *sp_int;		//This vector contains how many species each specie can interact to.
int **sp_adj; //This matrix contains for each species (lines) what are the other species that they can interact to (column). sp_adj[SP_NUMBER][SP_NUMBER]
				  //Represents the adjacency matrix for each specie.
int *tabu_list; //tabu list for "visited" species. Each specie must be placed in only one habitat. tabu_list[SP_NUMBER]
int *next_sp; //this vector is used in the procedure to create habitats. next_sp[SP_NUMBER]
double **distance_matrix; //contain the distances between centroids. Diagonal matrix. Euclidean distance. distance_matrix[SP_NUMBER][SP_NUMBER]
double **normalized_distance_matrix; //contain the normalized distances between centroids. Diagonal matrix. Euclidean distance. 
					//normalized_distance_matrix[SP_NUMBER][SP_NUMBER]
int *func_eval; //func_eval[SP_NUMBER]
double *Cc; //Cc[DIM] //Centroid of centroids
pthread_t *threads; //threads[SP_NUMBER]

double CcFO;    //Centroid FO
int h_count; //total number of habitats in each eco step.
int suc[1]; //counter for ecological succession loop
pthread_attr_t  pthread_custom_attr;
double Dlenght; //diagonal lenght

//aux functions
int GetParameters(char **argv);
void showParameters();
void freeArrays();
void AllocArrays();

//In ABC.c
int ABC(int sp_index);

//In DE.c
int DE(int sp_index, int strategy);
void  assignd(int D, double a[], double b[]);

//In PSO.c
int PSO(int sp_index);
double alea( double a, double b );
int alea_integer( int a, int b );
//Global variables for PSO
struct velocity
{
   int size;
   double **v; //v[POP_SIZE*2][DIM]
} *Vel; // Velocities; Vel[SP_NUMBER]

struct pbest_position
{
   int size;
   double **x; //x[POP_SIZE*2][DIM];
   double *f;  //f[POP_SIZE*2];
}*Pbest; //Pbest[SP_NUMBER]; // Best positions found by each particle;
 
//In DE_BBO.c
void Run_jDE_BBO(int sp_index);
int	rndint(int);
double	rndreal(double, double);
//Global variables for DE_BBO

//VNS PARAMETERS START
int KMAX; //the number of neighbourhood
double Q; //the geometric progresion rate
double RADII; //the radio aplied
int RADII_T_FLAG; //the flag for radii type test
int P; //the metric used
int G_MAX; //Max generation of pop
int T_MAX; //Max iteration
double *RAIO; //raio
int *k_index; //k radi for each individual
double PC = 0.9; //default value from DE
int KMAX; //default 5

//VNS PARAMETERS END

//in VNS.c
void prepare_prvns();
void PRVNS(int sp_index);
