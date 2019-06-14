/* 
	Description: The C code of the DE algorithms with MS VC++ 6.0
	Optimizer:   Basic DE; jDE; jDE-BBO 
	References:  1. J. of Global Optimization 1997, Basic DE
				 2. IEEE TEVC 2006, jDE
				 3. ..., jDE-BBO
	Programmer:  Wenyin Gong
	E-mail:      cug11100304@yahoo.com.cn
	Date:	     27/5/2009
	Lisence:     Free
	Note:		 If you use this code, pls cite the following paper:
	             W. Gong, Z. Cai, and C.X. Ling, DE/BBO: A Hybrid Differential 
				 Evolution with Biogeography-Based Optimization for Global Numerical 
				 Optimization, Soft Computing, In Press.
*/



/* Random number generator defined by URAND should return
   double-precision floating-point values uniformly distributed
   over the interval [0.0, 1.0)	*/

#define max_iteration	1500	/* maximal number of iteration */

void Run_jDE_BBO(int sp_index)
{

typedef struct 
{
	double *xreal; //xreal[DIM];	/* decision variables */
	double fo;			/* fo of the solution */
	double CRR;				/* for jDE and jDE-DE */
	double FF;				/* for jDE and jDE-DE */

	int	   SpeciesCount;	/* species count */
}individual;

individual *parent_pop; //parent_pop[POP_SIZE];/* save the parent population */
individual *child_pop;  //child_pop[POP_SIZE];	/* save the child population */
int		bestIndex;				/* index of the best solution in the current population */
int		worstIndex;				/* index of the worst solution in the current population */

int		gen_no;		
FILE	*fp;

int		evaluations;

individual best_individual;

/* for BBO */
double lambdaLower = 0.0;	/* lower bound for immigration probability per gene */
double lambdaUpper = 1.0;	/* upper bound for immigration probability per gene */
double II = 1.0;			/* max immigration rate for each island */
double EE = 1.0;			/* max emigration rate, for each island */
int    *sort_index; //sort_index[POP_SIZE];/* only for sorting the population */
double *lambda; //lambda[POP_SIZE];	/* the immigration rate */
double *mu; //mu[POP_SIZE];		/* the emigration rate */

	int i, j, k;
	FILE *fp1, *fp2;
	clock_t start, finish;
	double time_consuming;
	double low, up, sum;

	srand(time(NULL));

	//Alloc memory
	parent_pop = (individual*) malloc(cur_POP_SIZE[0] * sizeof(individual));
	for (i = 0; i < cur_POP_SIZE[0]; i++)
		parent_pop[i].xreal = (double*)malloc (sizeof(double)*DIM );

	child_pop = (individual*) malloc(cur_POP_SIZE[0] * sizeof(individual));
	for (i = 0; i < cur_POP_SIZE[0]; i++)
		child_pop[i].xreal = (double*)malloc (sizeof(double)*DIM );

	best_individual.xreal = (double*)malloc (sizeof(double)*DIM );

	lambda = (double*) malloc (cur_POP_SIZE[0] * sizeof(double));
	mu = (double*) malloc (cur_POP_SIZE[0] * sizeof(double));
	sort_index = (int*) malloc (cur_POP_SIZE[0] * sizeof(int));

	start = clock();

//	LOAD population
	if (suc[0] == 0)
	{
		for(i=0;i<cur_POP_SIZE[0];i++)
		{
			parent_pop[i].CRR = rndreal(0.0, 1.0);
			parent_pop[i].FF = rndreal(0.1, 1.0);
			Eco[sp_index].CRR[i] = parent_pop[i].CRR;
			Eco[sp_index].FF[i] = parent_pop[i].FF;
		}
	}

	for(i=0;i<cur_POP_SIZE[0];i++)
	{
		for(j=0;j<DIM;j++)
		{
			parent_pop[i].xreal[j] = Eco[sp_index].pop[i][j];
		}
		parent_pop[i].fo = Eco[sp_index].fo[i];
		if (suc[0] != 0)
		{
			parent_pop[i].CRR = Eco[sp_index].CRR[i];
			parent_pop[i].FF = Eco[sp_index].FF[i];
		}
	}

	gen_no = 1;
	
	double CRR = 0.9;
	double FF = 0.5;
	int    j_rnd;
	int    r1, r2, r3;
	int	   flag_eval = 0;
	evaluations = 0;
	while (evaluations < EVO_STEP)
//	while (gen_no < EVO_STEP)
	{
//		find_best_and_worst();
		bestIndex  = 0;
		worstIndex = 0;
		for (i=1;i<cur_POP_SIZE[0];i++)
		{
			if (parent_pop[i].fo <parent_pop[bestIndex].fo)
			{
				bestIndex = i;
			}
			if (parent_pop[i].fo > parent_pop[worstIndex].fo)
			{
				worstIndex = i;
			}
		}

//		copyIndividuals(&parent_pop[bestIndex], &best_individual);
//			       (individual *source, individual *target)
		for (i=0;i<DIM;i++)
		{
			best_individual.xreal[i] = parent_pop[bestIndex].xreal[i];
		}
		best_individual.fo      = parent_pop[bestIndex].fo;
		best_individual.CRR           = parent_pop[bestIndex].CRR;
		best_individual.FF           = parent_pop[bestIndex].FF;
		best_individual.SpeciesCount = parent_pop[bestIndex].SpeciesCount;


		if (flag_eval == 0 )
		{
			flag_eval = 1;
//			fprintf(fp2, "%d\n", evaluations);
		}
		
		finish = clock();
		time_consuming = (finish-start)/1000.0;

/*		if (gen_no%20 == 0)
		{
			fprintf(fp, "%6d %8d %e\n", gen_no, evaluations, best_individual.fo);
		}
*/
		// only for jDE-BBO
		// sort the population
		for (i=0;i<cur_POP_SIZE[0];i++)
		{
			sort_index[i] = i;
		}
//		shell_sort_pop(parent_pop, sort_index, POP_SIZE);
//		(individual *pop, int *list, int size)
		int done;
		int step, bound;
		int temp;

		step = cur_POP_SIZE[0];  // array length
		while (step > 1) 
		{
			step /= 2;	//halve the step size
			do 
			{
				done   = 1;
				bound  = cur_POP_SIZE[0] - step;
				for (j = 0; j < bound; j++) 
				{
					i = j + step + 1;
					if (parent_pop[sort_index[j]].fo > parent_pop[sort_index[i-1]].fo) 	
					{
						temp      = sort_index[i-1];
						sort_index[i-1] = sort_index[j];
						sort_index[j]   = temp;
						done = 0; // if a swap has been made we are not finished yet
					}  // if
				}  // for
			} while (done == 0);   // while
		} //while (step > 1)
		
		// calculate mu and lambda of each solution
		// Map cost values to species counts.
//		GetSpeciesCounts(parent_pop, POP_SIZE);
		for (i=0;i<cur_POP_SIZE[0];i++)
		{
			if (parent_pop[sort_index[i]].fo < INF)
			{
				parent_pop[sort_index[i]].SpeciesCount = cur_POP_SIZE[0]-i-1;
			}
			else
			{
				parent_pop[sort_index[i]].SpeciesCount = 0;
			}
		}

		// Compute immigration rate and emigration rate for each species count.
//		GetLambdaMu(parent_pop, POP_SIZE);
		for (i=0;i<cur_POP_SIZE[0];i++)
		{
			lambda[i] = II*(1.0-(double)parent_pop[sort_index[i]].SpeciesCount/cur_POP_SIZE[0]);
			mu[i] = EE*((double)parent_pop[sort_index[i]].SpeciesCount/cur_POP_SIZE[0]);
		}

		// Now use lambda and mu to decide how much information to share between habitats.
		double lambdaMin = lambda[0];
		double lambdaMax = lambda[0];
		for (j=1;j<cur_POP_SIZE[0];j++)
		{
			if (lambda[j] < lambdaMin)
			{
				lambdaMin = lambda[j];
			}
			if (lambda[j] > lambdaMax)
			{
				lambdaMax = lambda[j];
			}
		}

		/* DE/rand/1/bin operator */
		for (i=0;i<cur_POP_SIZE[0];i++)
		{
			do {
				r1 = rndint(cur_POP_SIZE[0]);
			} while(r1 == i);
			do {
				r2 = rndint(cur_POP_SIZE[0]);
			} while(r2 == i || r2 == r1);
			do {
				r3 = rndint(cur_POP_SIZE[0]);
			} while(r3 == i || r3 == r2 || r3 == r1);

			int mate1 = sort_index[r1];
			int mate2 = sort_index[r2];
			int mate3 = sort_index[r3]; 

			/* parameter self-adaptation technique proposed by J. Brest */
			child_pop[i].FF = parent_pop[sort_index[i]].FF;
			child_pop[i].CRR = parent_pop[sort_index[i]].CRR;
			if (rndreal(0.0, 1.0) < 0.1)
			{
				child_pop[i].FF = rndreal(0.1, 1.0);
			}
			if (rndreal(0.0, 1.0) < 0.1)
			{
				child_pop[i].CRR = rndreal(0.0, 1.0);
			}
			CRR = child_pop[i].CRR;
			FF = child_pop[i].FF;

			// Normalize the immigration rate.
			double lambdaScale=lambda[i];
			lambdaScale = lambdaLower +
				(lambdaUpper-lambdaLower)*(lambda[i]-lambdaMin)/(lambdaMax-lambdaMin);

			j_rnd = rndint(DIM);

			for (j=0;j<DIM;j++)
			{
				//if (rndreal(0,1) < lambdaScale)
				if (1) // especially for the rotated functions
				{
					if ( (rndreal(0,1) < CRR) || (j == j_rnd) )
					{// differential operator
						double low = lb;
						double up  = ub;
						
						child_pop[i].xreal[j] = parent_pop[mate1].xreal[j]
							+ FF*(parent_pop[mate2].xreal[j]-parent_pop[mate3].xreal[j]);

						if (child_pop[i].xreal[j] < low || 
							child_pop[i].xreal[j] > up)
						{
							child_pop[i].xreal[j] = rndreal(low, up);
						}
					}
					else
					{// BBO migration	
						int SelectIndex;
						SelectIndex = rndint(cur_POP_SIZE[0]);
						if (rndreal(0,1) < mu[SelectIndex])
						{
							child_pop[i].xreal[j] = 
								parent_pop[sort_index[SelectIndex]].xreal[j];
						}
						else
						{
							child_pop[i].xreal[j] = parent_pop[sort_index[i]].xreal[j];
						}
					}
				}
				else
				{
					child_pop[i].xreal[j] = parent_pop[sort_index[i]].xreal[j];
				}				
			}
		}
		
		/* evaluate the children population */
//		evaluate_pop(child_pop, POP_SIZE);
		for (i=0;i<cur_POP_SIZE[0];i++)
		{
			child_pop[i].fo = function(child_pop[i].xreal, sp_index);
			evaluations++;
		}
		
		/* selection */
		for (i=0;i<cur_POP_SIZE[0];i++)
		{
			if (child_pop[i].fo <= parent_pop[sort_index[i]].fo)
			{
//				copyIndividuals(&child_pop[i], &parent_pop[sort_index[i]]);
				for (j=0;j<DIM;j++)
				{
					parent_pop[sort_index[i]].xreal[j] = child_pop[i].xreal[j];
				}
				parent_pop[sort_index[i]].fo      = child_pop[i].fo;
				parent_pop[sort_index[i]].CRR           = child_pop[i].CRR;
				parent_pop[sort_index[i]].FF           = child_pop[i].FF;
				
				parent_pop[sort_index[i]].SpeciesCount = child_pop[i].SpeciesCount;

			}
		}

		gen_no++;
	}

//		find_best_and_worst();
		bestIndex  = 0;
		worstIndex = 0;
		for (i=1;i<cur_POP_SIZE[0];i++)
		{
			if (parent_pop[i].fo <parent_pop[bestIndex].fo)
			{
				bestIndex = i;
			}
			if (parent_pop[i].fo > parent_pop[worstIndex].fo)
			{
				worstIndex = i;
			}
		}

//		copyIndividuals(&parent_pop[bestIndex], &best_individual);
//			       (individual *source, individual *target)
		for (i=0;i<DIM;i++)
		{
			best_individual.xreal[i] = parent_pop[bestIndex].xreal[i];
		}
		best_individual.fo      = parent_pop[bestIndex].fo;
		best_individual.CRR           = parent_pop[bestIndex].CRR;
		best_individual.FF           = parent_pop[bestIndex].FF;
		best_individual.SpeciesCount = parent_pop[bestIndex].SpeciesCount;


  //update ecosystem variable
  for ( i = 0; i < cur_POP_SIZE[0]; i++ )  // Positions
  {
    Eco[sp_index].fo[i] = parent_pop[i].fo;
    Eco[sp_index].CRR[i] = parent_pop[i].CRR;
    Eco[sp_index].FF[i] = parent_pop[i].FF;
    for ( j = 0; j < DIM; j++ )
    {
	Eco[sp_index].pop[i][j] = parent_pop[i].xreal[j];
    }
  }

  Eco[sp_index].bestfo[0]=best_individual.fo;
  Eco[sp_index].best_index[0]=bestIndex;
  for (i=0;i<DIM;i++)
  	Eco[sp_index].best[i] = best_individual.xreal[i]; 


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
			evaluations++;




/*
	printf("gen=%d\teval=%d\tfo=%e\t%f s.\n", 
		gen_no, evaluations, best_individual.fo, time_consuming);
	fprintf(fp1, "%d\t%d\t%e\t\n", gen_no, evaluations, 
		best_individual.fo, time_consuming);

	fclose(fp);
	fclose(fp1);
	fclose(fp2);
*/

   for (j = 0;j < cur_POP_SIZE[0]; j++)
   {
		free(parent_pop[j].xreal);
		free(child_pop[j].xreal);
   }
   free(parent_pop);
   free(child_pop);
   free(lambda);
   free(mu);
   free(sort_index);
   free(best_individual.xreal);
}

/* random number generator */
double rndreal(double low,double high)
{
	double  val;
	val = URAND*(high-low)+low;
	return (val);
}

int rndint(int Nmem)
{
	int newval;
	newval = (int)(URAND*Nmem);
	return(newval);
}

