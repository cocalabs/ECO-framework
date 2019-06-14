/* Standard PSO version 2006

   Motivation
   Quite often some authors say they compare their PSO versions
   to the "standard one" ... which is never the same!
   So the idea is to define a real standard at least for one year, validated by some
   researchers of the field, in particular James Kennedy and Maurice Clerc.
   This PSO version does not intend to be the best one on the market (in particular there is
   no adaptation of the swarm size nor of the coefficients) but simply very near of the
   original version (1995) with just a few improvements based on some recent works.

   So referring to "standard PSO 2006" would mean exactly this version with the default values
   detailed below as,for example, referring to "standard PSO 2006 (w=0.8)" would mean almost
   this version but with a non standard first cognitive/confidence coefficient.

   Parameters
S := swarm size
K := maximum number of particles _informed_ by a given one
T := topology of the information links
w := first cognitive/confidence coefficient
c := second cognitive/confidence coefficient
R := random distribution of c
B := rule "to keep the particle in the box"

Equations
For each particle and each dimension
v(t+1) = w*v(t) + R(c)*(p(t)-x(t)) + R(c)*(g(t)-x(t))
x(t+1) = x(t) + v(t+1)
where
v(t) := velocity at time t
x(t) := position at time t
p(t) := best previous position of the particle
g(t) := best previous position of the informants of the particle

Default values
S = 10+2*sqrt(D) where D is the dimension of the search space
K = 3
T := randomly modified after each step if there has been no improvement
w = 1/(2*ln(2))
c = 0.5 + ln(2)
R = U(0..c), i.e. uniform distribution on [0, c]
B := set the position to the min. (max.) value and the velocity to zero

About information links topology
A lot of works have been done about this topic. The main result is there is no
"best" topology. Hence the random approach used here. Note that is does not mean that
each particle is informed by K ones: the number of particles that informs a given one
may be any value between 1 (for each particle informs itself) and S.

About initialisation
Initial positions are chosen at random inside the search space (which is
supposed to be a hyperparallelepid, and even often a hypercube), according to an uniform
distribution. This is not the best way, but the one of the original PSO.

Each initial velocity is simply defined as the difference of two random
positions. It is simple, and needs no additional parameter.
However, again, it is not the best approach. The resulting distribution is not even
uniform, as for any method that uses an uniform distribution independently for each
component. The mathematically correct approach needs to use an uniform
distribution inside a hypersphere. It is not that difficult, and indeed used in some PSO
versions, but quite different from the original one.

Some results with the standard values
You may want to recode this program in another language. Also you may want to modify it
for your own purposes. Here are some results on classical test functions to help you to
check your code.
Dimension D=30
Acceptable error eps=0
Objective value f_min=0
Number of runs n_exec_max=50
Number of evaluations for each run eval_max=30000

Problem                            Mean best value    Standard deviation
Parabola/Sphere on [-100, 100]^D        0                  0
Griewank on [-600, 600]^D               0.018              0.024
Rosenbrock/Banana on [-30, 30]^D       50.16              36.9
Rastrigin on [-5.12, 5.12]^D           48.35              14.43
Ackley on [-32, 32]^D                   1.12               0.85

Last updates
2006-02-27 Fixed a bug about minimal best value over several runs
2006-02-16 Fixed a bug (S_max for P, V, X, instead of D_max), thanks to Manfred Stickel
2006-02-16 replaced k by i x by xs (in perf()), because of possible confusion with K and X
2006-02-13  added De Jong's f4
*/
// =================================================
int PSO(int sp_index)
{
	struct position
	{
		int size;
		double *x;//[DIM];
		double f;
	}*X; //[cur_POP_SIZE[0]]; // Positions;

	int best; // Best of the best position (rank in the swarm)
	int D; // Search space dimension
	double E; // exp(1). Useful for some test functions
	double f_min; // Objective(s) to reach
	int **LINKS; //[POP_SIZE] [POP_SIZE]; // Information links
	int nb_eval; // Total number of evaluations
	double pi; // Useful for some test functions
	//struct position P[S_max]; // Best positions found by each particle
	int S; // Swarm size
	double *xmin, // [DIM]
	       *xmax; //[DIM]; // Intervals defining the search space

	// File(s)

	double c; // Second onfidence coefficient
	int d; // Current dimension
	double eps; // Admissible error
	double eps_mean; // Average error
	double error; // Error for a given position
	double error_prev; // Error after previous iteration
	int eval_max; // Max number of evaluations
	double eval_mean; // Mean number of evaluations
	//  int function; // Code of the objective function
	int g; // Rank of the best informant
	int init_links; // Flag to (re)init or not the information links
	int i;
	int K; // Max number of particles informed by a given one
	int m;
	double *mean_best; //[RUN];
	double min; // Best result through several runs
	int n_exec, n_exec_max; // Nbs of executions
	int n_failure; // Number of failures
	int s; // Rank of the current particle
	double t1, t2;
	double variance;
	double w; // First confidence coefficient
	int k,j,
	    iter; //iterations
	int local_func_eval;
	double sum; //to calculate the centroid

	//  f_run = fopen( "f_run.txt", "w" );
	E = exp( 1 );
	pi = acos( -1 );

	D = DIM; // Search space dimension

	X = (struct position*) malloc(cur_POP_SIZE[0] * sizeof(struct position));
	for (i = 0; i < cur_POP_SIZE[0]; i++)
	{
		X[i].x = (double*)malloc (DIM * sizeof(double));
	}//for i

	xmin = (double*)malloc (DIM * sizeof(double));
	xmax = (double*)malloc (DIM * sizeof(double));

	mean_best = (double*)malloc (RUN * sizeof(double));

	// D-cube data
	for ( d = 0; d < D; d++ )
	{
		xmin[d] = lb; xmax[d] = ub;
	}

	eps = 0.9999; // Acceptable error
	f_min = 0; // Objective value
	n_exec_max = RUN; // Numbers of runs
	//  eval_max = 1120; // Max number of evaluations for each run

	//-----------------------------------------------------  PARAMETERS
	S = cur_POP_SIZE[0];
	K = 3;

	//printf("\n Swarm size %i", S);
	//printf("\n coefficients %f %f \n",W,C);

	//----------------------------------------------------- INITIALISATION
	//t1 = clock(); // Init time
	// Initialisation of information variables
	n_exec = 0; eval_mean = 0; eps_mean = 0; n_failure = 0;

init:
	n_exec = n_exec + 1;
	local_func_eval = 0;

	for ( s = 0; s < S; s++ )  // Positions
	{
		X[s].size = D; Vel[sp_index].size = D;
		X[s].f = Eco[sp_index].fo[s];
		for ( d = 0; d < D; d++ )
		{
			X[s].x[d] = Eco[sp_index].pop[s][d]; //load population from specie 'sp_index'.
		}
	}

	LINKS = malloc (cur_POP_SIZE[0] * sizeof(int*));
	for (j = 0;j < cur_POP_SIZE[0]; j++)
		LINKS[j] = (int*)malloc (cur_POP_SIZE[0] * sizeof(int));

	if (suc[0] == 0)//initialize velocities. Olnly for the first ecological succession.
	{
		for ( s = 0; s < S; s++ )  // Positions and velocities
		{
			for ( d = 0; d < D; d++ )
			{
				Pbest[sp_index].x[s][d] = X[s].x[d]; // Best position = current one
				Pbest[sp_index].f[s] = X[s].f;
				Vel[sp_index].v[s][d] = (alea( xmin[d], xmax[d] ) - X[s].x[d])/2.0; // Non uniform
				// V[s].v[d] = ( xmin[d]-xmax[d] )*(0.5-alea(0,1)); //Uniform. 2006-02-24
			}
		}
	}//END IF 

	best = Eco[sp_index].best_index[0];
	Pbest[sp_index].f[best] = Eco[sp_index].bestfo[0];
	for ( d = 0; d < D; d++ )
	{
		Pbest[sp_index].x[best][d] = Eco[sp_index].best[d]; // Best position = current one
	}
	// First evaluations
	nb_eval = 0;
	/*    for ( s = 0; s < S; s++ )
	      {
	      X[s].f = function(X[s].x); //fabs( perf( s, function ) - f_min );
	      P[s] = X[s]; // Best position = current one
	      }
	      */
	// Find the best
	/*    best = 0;
	      for ( s = 1; s < S; s++ )
	      if ( Pbest[sp_index].f[s] < Pbest[sp_index].f[best] ) best = s;
	      */
	error =  Pbest[sp_index].f[best] ; // Current min error
	//printf("BEst FO: %.2f \n", Pbest[sp_index].f[best]);
	//printf("Error: %.2f \n", error);
	if(n_exec==1) min=error;
	error_prev=error; // Previous min error

	init_links = 1; // So that information links will be initialized

	//---------------------------------------------- ITERATIONS
	iter = -1;
loop:
	//printf("Error: %.2f \n", error);
	if ( init_links == 1 )
	{
		// Who informs who, at random
		for ( s = 0; s < S; s++ )
		{
			for ( m = 0; m < S; m++ ) 
			{
				LINKS[0][0] = 0;
				LINKS[m][s] = 0; // Init to "no link"
			}
			LINKS[s][s] = 1; // Each particle informs itself
		}
		for ( m = 0; m < S; m++ ) // Other links
		{
			for ( i = 0; i < K; i++ )
			{
				s = alea_integer( 0, S - 1 );
				LINKS[m] [s] = 1;
			}
		}
	}

	// The swarm MOVES
	for ( s = 0; s < S; s++ ) // For each particle ...
	{
		// .. find the best informant
		g=s;
		for ( m = 0; m < S; m++ )
		{
			if ( LINKS[m] [s] == 1 && Pbest[sp_index].f[m]<Pbest[sp_index].f[g]  ) g = m;
		}

		// ... compute the new velocity, and move
		for ( d = 0; d < D; d++ )
		{

			Vel[sp_index].v[s][d] = W * Vel[sp_index].v[s][d] + alea( 0, C ) * ( Pbest[sp_index].x[s][d] - X[s].x[d] );
			Vel[sp_index].v[s][d] = Vel[sp_index].v[s][d] + alea( 0, C ) * ( Pbest[sp_index].x[g][d] - X[s].x[d] );
			X[s].x[d] = X[s].x[d] + Vel[sp_index].v[s][d] ;
		}

		// ... interval confinement (keep in the box)
		for ( d = 0; d < D; d++ )
		{
			if ( X[s].x[d] < xmin[d] )
			{
				X[s].x[d] = xmin[d]; Vel[sp_index].v[s][d] = 0;
			}
			if ( X[s].x[d] > xmax[d] )
			{
				X[s].x[d] = xmax[d]; Vel[sp_index].v[s][d] = 0;
			}
		}


		// ... evaluate the new position
		X[s].f = function(X[s].x, sp_index);//fabs( perf( s, function ) - f_min );
		local_func_eval++;


		// ... update the best previous position
		if ( X[s].f<Pbest[sp_index].f[s] )
		{
			for ( d = 0; d < D; d++ )
			{
				Pbest[sp_index].x[s][d] = X[s].x[d]; // update position
			}
			Pbest[sp_index].f[s] = X[s].f;
			// ... update the best of the bests
			if (  Pbest[sp_index].f[s]<Pbest[sp_index].f[best] ) best = s;
		}
		//printf("func_eval %d\n", local_func_eval);

		//verity stop condition
		if ( local_func_eval > EVO_STEP-1 ) goto stat;

	}//End FOR Swarm Moves
	//printf("TESTE DEPOIS... \n");
	// Check if finished
	// If no improvement, information links will be reinitialized
	error=Pbest[sp_index].f[best];
	if ( error >= error_prev ) init_links = 1;
	else init_links = 0;
	error_prev = error;

	//    if ( error > eps && func_eval < eval_max ) goto loop;
	iter++;
	if ( local_func_eval < EVO_STEP ) goto loop;

stat:
	if ( error > eps ) n_failure = n_failure + 1;

	//printf("\nLocal Func Eval: %d\n",local_func_eval);

	// Result display
	//    printf( "\nExec %i Eval %i. Error %f ", n_exec, func_eval, error );
	//    printf( "\n Position :\n" );
	//    for ( d = 0; d < D; d++ ) printf( " %f", Pbest[sp_index].x[best][d] );

	// Save result
	//    fprintf( f_run, "\n%i %i %f ", n_exec, func_eval,error );
	//    fprintf( f_run, " Position: " );
	//    for ( d = 0; d < D; d++ ) fprintf( f_run, " %f", Pbest[sp_index].x[best][d] );

	// Compute some statistical information
	if ( error < min ) min = error;
	eval_mean = eval_mean + func_eval[sp_index];
	eps_mean = eps_mean + error;
	mean_best[n_exec - 1] = error;

	//    if ( n_exec < n_exec_max ) goto init;

	// END. Display some statistical information
	//    t2 = clock();
	//    printf( "\n\n Total clocks %.0f", t2 - t1 );
	eval_mean = eval_mean / ( double )n_exec;
	eps_mean = eps_mean / ( double )n_exec;
	//    printf( "\n\n Eval. (mean)= %.2f", eval_mean );
	//    printf( "\n Error (mean) = %.2f", eps_mean );

	// Variance
	variance = 0;
	for ( d = 0; d < n_exec_max; d++ ) variance = variance + ( mean_best[d] - eps_mean ) * ( mean_best[d] - eps_mean );
	variance = sqrt( variance / n_exec_max );
	//    printf( "\n Std. dev. %f", variance );

	// Success rate and minimum value
	//    printf( "\n Success rate = %.2f%%", 100 * ( 1 - n_failure / ( double )n_exec ) );
	//    if ( n_exec > 0 ) printf( "\n Best min value = %f", min );

end:;

    //update ecosystem variable
    for ( s = 0; s < S; s++ )  // Positions and velocities
    {
	    Eco[sp_index].fo[s] = X[s].f;
	    for ( d = 0; d < D; d++ )
	    {
		    Eco[sp_index].pop[s][d] = X[s].x[d];
	    }
    }

    Eco[sp_index].bestfo[0]=Pbest[sp_index].f[best];
    Eco[sp_index].best_index[0]=best;

    for (d=0;d<DIM;d++)
	    Eco[sp_index].best[d] = Pbest[sp_index].x[best][d]; 


    //DEFINE Ci value for each specie. Place where a specie is concentrated.
    //In this implementation, the Ci value is the population centroid.
    for (k=0; k<DIM;k++)//each variable
    {
	    sum = 0.0;
	    for (j=0;j<cur_POP_SIZE[0];j++) //each solution
	    {
		    sum += Eco[sp_index].pop[j][k];
	    }
	    Eco[sp_index].Ci[k] = sum/(double)cur_POP_SIZE[0];
    }
    //			printf("Teste1: %.2f, %.2f, %d\n",Eco[i].Ci[0],Eco[i].Ci[1], DIM);
    Eco[sp_index].CiFO[0] = function(Eco[sp_index].Ci, sp_index);
    local_func_eval++;


    for (j = 0;j < cur_POP_SIZE[0]; j++)
	    free(X[j].x);
    free(X);

    for (j = 0;j < cur_POP_SIZE[0]; j++)
	    free(LINKS[j]);
    free(LINKS);

    free(xmin);
    free(xmax);
    free(mean_best);


    //    return 0;
}

//===========================================================
double alea( double a, double b )
{ // random number (uniform distribution) in [a b]
	double r;
	r=(double)rand(); r=r/RAND_MAX;
	return a + r * ( b - a );
}
//===========================================================
int alea_integer( int a, int b )
{ // Integer random number in [a b]
	int ir;
	double r;
	r = alea( 0, 1 ); ir = ( int )( a + r * ( b + 1 - a ) );
	if ( ir > b ) ir = b;
	return ir;
}

