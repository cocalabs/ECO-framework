//Differential Evolution Code
int DE(int sp_index, int   strategy)
{
	char  chr;             /* y/n choice variable                */
	char  *strat[] =       // strategy-indicator                
	{
		"",
		"DE/best/1/exp",
		"DE/rand/1/exp",
		"DE/rand-to-best/1/exp",
		"DE/best/2/exp",
		"DE/rand/2/exp",
		"DE/best/1/bin",
		"DE/rand/1/bin",
		"DE/rand-to-best/1/bin",
		"DE/best/2/bin",
		"DE/rand/2/bin"
	};

	int   i, j, L, n, k;      /* counting variables                 */
	int   r1, r2, r3, r4;  /* placeholders for random indexes    */
	int   r5;              /* placeholders for random indexes    */
	int   refresh;         /* refresh rate of screen output      */
	//   int   strategy;        /* choice parameter for screen output */
	int   gen;   

	long  nfeval;          /* number of function evaluations     */

	double trial_cost;      /* buffer variable                    */
	double cvar;            /* computes the cost variance         */
	double cmean;           /* mean cost                          */
	//   double F,CR;            /* control variables of DE            */
	double cmin;            /* help variables                     */
	double sum;
	double *tmp; //tmp[DIM];
	double cost[cur_POP_SIZE[0]];    /* obj. funct. values                 */
	double *solution; //solution [DIM];

	/*------Initializations----------------------------*/
	//strategy = sp_index+1
	// strategy = 7;          //---choice of strategy---------------should be ex {1,2,3,4,5,6,7,8,9,10}

	/*-----Initialize random number generator-----------------------------*/

	// cost = (double*)malloc(cur_POP_SIZE[0] * sizeof(double));
	tmp = (double*)malloc(DIM * sizeof(double));
	solution = (double*)malloc(DIM * sizeof(double));

	srand(time(NULL));
	//MT_seed();

	/*------Initialization------------------------------------------------*/
	/*------Right now this part is kept fairly simple and just generates--*/
	/*------random numbers in the range [-initfac, +initfac]. You might---*/
	/*------want to extend the init part such that you can initialize-----*/
	/*------each parameter separately.------------------------------------*/
	/*=======================================================================*/
	/*=========Iteration loop================================================*/
	/*=======================================================================*/
	nfeval       =  0;  /* reset number of function evaluations */
	gen = 0;                          /* generation counter reset */

	//printf("ANTES: POP %d nfeval: %ld\n", cur_POP_SIZE[0], nfeval);

	while (nfeval < EVO_STEP)
	{                                            /* is accepted by compiler    */
		gen++;

		for (i=0; i<cur_POP_SIZE[0]; i++)         /* Start of loop through ensemble  */
		{
			do                        /* Pick a random population member */
			{                         /* Endless loop for POP_SIZE < 2 !!!     */
				r1=(int)( ((double)rand() / ((double)(RAND_MAX)+(double)(1)) ) * cur_POP_SIZE[0]);
			}while(r1==i);            

			do                        /* Pick a random population member */
			{                         /* Endless loop for POP_SIZE < 3 !!!     */
				r2 = (int)( ((double)rand() / ((double)(RAND_MAX)+(double)(1)) ) * cur_POP_SIZE[0]);
			}while((r2==i) || (r2==r1));

			do                        /* Pick a random population member */
			{                         /* Endless loop for POP_SIZE < 4 !!!     */
				r3 = (int)( ((double)rand() / ((double)(RAND_MAX)+(double)(1)) ) * cur_POP_SIZE[0]);
			}while((r3==i) || (r3==r1) || (r3==r2));

			do                        /* Pick a random population member */
			{                         /* Endless loop for POP_SIZE < 5 !!!     */
				r4 = (int)( ((double)rand() / ((double)(RAND_MAX)+(double)(1)) ) * POP_SIZE);
			}while((r4==i) || (r4==r1) || (r4==r2) || (r4==r3));

			do                        /* Pick a random population member */
			{                         /* Endless loop for POP_SIZE < 6 !!!     */
				r5 = (int)( ((double)rand() / ((double)(RAND_MAX)+(double)(1)) ) * cur_POP_SIZE[0]);
			}while((r5==i) || (r5==r1) || (r5==r2) || (r5==r3) || (r5==r4));
			if (r1 == cur_POP_SIZE[0])
				r1 = cur_POP_SIZE[0]-1;
			if (r2 == cur_POP_SIZE[0])
				r2 = cur_POP_SIZE[0]-1;
			if (r3 == cur_POP_SIZE[0])
				r3 = cur_POP_SIZE[0]-1;
			if (r4 == cur_POP_SIZE[0])
				r4 = cur_POP_SIZE[0]-1;
			if (r5 == cur_POP_SIZE[0])
				r5 = cur_POP_SIZE[0]-1;
			//printf("\nTESTE!!! sp%d %d %d %d %d %d \n",sp_index, r1, r2, r3, r4, r5);

			/*=======Choice of strategy===============================================================*/
			/*=======We have tried to come up with a sensible naming-convention: DE/x/y/z=============*/
			/*=======DE :  stands for Differential Evolution==========================================*/
			/*=======x  :  a string which denotes the vector to be perturbed==========================*/
			/*=======y  :  number of difference vectors taken for perturbation of x===================*/
			/*=======z  :  crossover method (exp = exponential, bin = binomial)=======================*/
			/*                                                                                        */
			/*=======There are some simple rules which are worth following:===========================*/
			/*=======1)  F is usually between 0.5 and 1 (in rare cases > 1)===========================*/
			/*=======2)  CR is between 0 and 1 with 0., 0.3, 0.7 and 1. being worth to be tried first=*/
			/*=======3)  To start off POP_SIZE = 10*D is a reasonable choice. Increase POP_SIZE if misconvergence=*/
			/*           happens.                                                                     */
			/*=======4)  If you increase POP_SIZE, F usually has to be decreased============================*/
			/*=======5)  When the DE/best... schemes fail DE/rand... usually works and vice versa=====*/


			/*=======EXPONENTIAL CROSSOVER============================================================*/

			/*-------DE/best/1/exp--------------------------------------------------------------------*/
			/*-------Our oldest strategy but still not bad. However, we have found several------------*/
			/*-------optimization problems where misconvergence occurs.-------------------------------*/
			if (strategy == 1) /* strategy DE0 (not in our paper) */
			{
				assignd(DIM,tmp,Eco[sp_index].pop[i]);
				n = (int)( ((double)rand()/((double)(RAND_MAX)+(double)(1)) ) * DIM);
				L = 0;
				do
				{                       
					tmp[n] = Eco[sp_index].best[n] + F*(Eco[sp_index].pop[r2][n]-Eco[sp_index].pop[r3][n]);
					n = (n+1)%DIM;
					L++;
				}while(( ( (double)rand()/((double)(RAND_MAX)+(double)(1)) < CR ) && (L < DIM) ));
			}
			/*-------DE/rand/1/exp-------------------------------------------------------------------*/
			/*-------This is one of my favourite strategies. It works especially well when the-------*/
			/*-------"bestit[]"-schemes experience misconvergence. Try e.g. F=0.7 and CR=0.5---------*/
			/*-------as a first guess.---------------------------------------------------------------*/
			else if (strategy == 2) /* strategy DE1 in the techreport */
			{
				//	   assignd(DIM,tmp,Eco[sp_index].pop[i]);
				for (j=0;j<DIM;j++)
					tmp[j] = Eco[sp_index].pop[i][j];
				n = (int)( ((double)rand()/((double)(RAND_MAX)+(double)(1)) ) * DIM);
				L = 0;
				do
				{                       
					tmp[n] = Eco[sp_index].pop[r1][n] + F*(Eco[sp_index].pop[r2][n]-Eco[sp_index].pop[r3][n]);
					n = (n+1)%DIM;
					L++;
				}while(( ( (double)rand()/((double)(RAND_MAX)+(double)(1)) < CR ) && (L < DIM) ));
			}
			/*-------DE/rand-to-best/1/exp-----------------------------------------------------------*/
			/*-------This strategy seems to be one of the best strategies. Try F=0.85 and CR=1.------*/
			/*-------If you get misconvergence try to increase POP_SIZE. If this doesn't help you----------*/
			/*-------should play around with all three control variables.----------------------------*/
			else if (strategy == 3) /* similiar to DE2 but generally better */
			{ 
				assignd(DIM,tmp,Eco[sp_index].pop[i]);
				n = (int)( ((double)rand()/((double)(RAND_MAX)+(double)(1)) ) * DIM);
				L = 0;
				do
				{                       
					tmp[n] = tmp[n] + F*(Eco[sp_index].best[n] - tmp[n]) + F*(Eco[sp_index].pop[r1][n]-Eco[sp_index].pop[r2][n]);
					n = (n+1)%DIM;
					L++;
				}while(( ( (double)rand()/((double)(RAND_MAX)+(double)(1)) < CR ) && (L < DIM) ));
			}
			/*-------DE/best/2/exp is another powerful strategy worth trying--------------------------*/
			else if (strategy == 4)
			{ 
				assignd(DIM,tmp,Eco[sp_index].pop[i]);
				n = (int)( ((double)rand()/((double)(RAND_MAX)+(double)(1)) ) * DIM);
				L = 0;
				do
				{                           
					tmp[n] = Eco[sp_index].best[n] + 
						(Eco[sp_index].pop[r1][n]+Eco[sp_index].pop[r2][n]-Eco[sp_index].pop[r3][n]-Eco[sp_index].pop[r4][n])*F;
					n = (n+1)%DIM;
					L++;
				}while(( ( (double)rand()/((double)(RAND_MAX)+(double)(1)) < CR ) && (L < DIM) ));
			}
			/*-------DE/rand/2/exp seems to be a robust optimizer for many functions-------------------*/
			else if (strategy == 5)
			{ 
				assignd(DIM,tmp,Eco[sp_index].pop[i]);
				n = (int)( ((double)rand()/((double)(RAND_MAX)+(double)(1)) ) * DIM);
				L = 0;
				do
				{                           
					tmp[n] = Eco[sp_index].pop[r5][n] + 
						(Eco[sp_index].pop[r1][n]+Eco[sp_index].pop[r2][n]-Eco[sp_index].pop[r3][n]-Eco[sp_index].pop[r4][n])*F;
					n = (n+1)%DIM;
					L++;
				}while(( ( (double)rand()/((double)(RAND_MAX)+(double)(1)) < CR ) && (L < DIM) ));
			}

			/*=======Essentially same strategies but BINOMIAL CROSSOVER===============================*/

			/*-------DE/best/1/bin--------------------------------------------------------------------*/
			else if (strategy == 6) 
			{
				assignd(DIM,tmp,Eco[sp_index].pop[i]);
				n = (int)( ((double)rand()/((double)(RAND_MAX)+(double)(1)) ) * DIM);
				for (L=0; L<DIM; L++) /* perform D binomial trials */
				{
					if (((double)rand()/((double)(RAND_MAX)+(double)(1)) < CR) || L == (DIM-1)) /* change at least one parameter */
					{                       
						tmp[n] = Eco[sp_index].best[n] + F*(Eco[sp_index].pop[r2][n]-Eco[sp_index].pop[r3][n]);
					}
					n = (n+1)%DIM;
				}
			}
			/*-------DE/rand/1/bin-------------------------------------------------------------------*/
			else if (strategy == 7) 
			{
				assignd(DIM,tmp,Eco[sp_index].pop[i]);
				n = (int)( ((double)rand()/((double)(RAND_MAX)+(double)(1)) ) * DIM);
				for (L=0; L<DIM; L++) /* perform D binomial trials */
				{
					if (((double)rand()/((double)(RAND_MAX)+(double)(1)) < CR) || L == (DIM-1)) /* change at least one parameter */
					{                       
						tmp[n] = Eco[sp_index].pop[r1][n] + F*(Eco[sp_index].pop[r2][n]-Eco[sp_index].pop[r3][n]);
					}
					n = (n+1)%DIM;
				}
			}
			/*-------DE/rand-to-best/1/bin-----------------------------------------------------------*/
			else if (strategy == 8) 
			{ 
				assignd(DIM,tmp,Eco[sp_index].pop[i]);
				n = (int)( ((double)rand()/((double)(RAND_MAX)+(double)(1)) ) * DIM);
				for (L=0; L<DIM; L++) /* perform D binomial trials */
				{
					if (((double)rand()/((double)(RAND_MAX)+(double)(1)) < CR) || L == (DIM-1)) /* change at least one parameter */
					{                       
						tmp[n] = tmp[n] + F*(Eco[sp_index].best[n] - tmp[n]) + F*(Eco[sp_index].pop[r1][n]-Eco[sp_index].pop[r2][n]);
					}
					n = (n+1)%DIM;
				}
			}
			/*-------DE/best/2/bin--------------------------------------------------------------------*/
			else if (strategy == 9)
			{ 
				assignd(DIM,tmp,Eco[0].pop[i]);
				n = (int)( ((double)rand()/((double)(RAND_MAX)+(double)(1)) ) * DIM);
				for (L=0; L<DIM; L++) /* perform D binomial trials */
				{
					if (((double)rand()/((double)(RAND_MAX)+(double)(1)) < CR) || L == (DIM-1)) /* change at least one parameter */
					{                       
						tmp[n] = Eco[sp_index].best[n] + 
							(Eco[sp_index].pop[r1][n]+Eco[sp_index].pop[r2][n]-Eco[sp_index].pop[r3][n]-Eco[sp_index].pop[r4][n])*F;
					}
					n = (n+1)%DIM;
				}
			}
			/*-------DE/rand/2/bin--------------------------------------------------------------------*/
			else if (strategy == 0)
			{ 
				assignd(DIM,tmp,Eco[sp_index].pop[i]);
				n = (int)( ((double)rand()/((double)(RAND_MAX)+(double)(1)) ) * DIM);
				for (L=0; L<DIM; L++) /* perform D binomial trials */
				{
					if (((double)rand()/((double)(RAND_MAX)+(double)(1)) < CR) || L == (DIM-1)) /* change at least one parameter */
					{                       
						tmp[n] = Eco[sp_index].pop[r5][n] + 
							(Eco[sp_index].pop[r1][n]+Eco[sp_index].pop[r2][n]-Eco[sp_index].pop[r3][n]-Eco[sp_index].pop[r4][n])*F;
					}
					n = (n+1)%DIM;
				}
			}


			/*=======Trial mutation now in tmp[]. Test how good this choice really was.==================*/
			for (L=0; L<DIM; L++)
			{
				if (tmp[L] < lb || tmp[L] > ub)
					tmp[L] = randon(lb,ub);
			}

			trial_cost = function(tmp, sp_index);  /* Evaluate new vector in tmp[] */
			//printf("trial: %.4f fo-%d:%.4f\n", trial_cost, i, Eco[sp_index].fo[i]);
			nfeval++;
			if (trial_cost <= Eco[sp_index].fo[i])   /* improved objective function value ? */
			{                                  
				Eco[sp_index].fo[i]=trial_cost;       

				for (j=0;j<DIM;j++)
					Eco[sp_index].pop[i][j] = tmp[j];  
				//	    assignd(DIM,Eco[sp_index].pop[i],tmp);

				if (trial_cost<=Eco[sp_index].bestfo[0])          /* Was this a new minimum? */
				{                               /* if so...*/
					Eco[sp_index].bestfo[0]=trial_cost;           /* reset cmin to new low...*/
					Eco[sp_index].best_index[0]=i;

					for (j=0;j<DIM;j++)
						Eco[sp_index].best[j] = tmp[j];  
					//	       assignd(DIM,Eco[sp_index].best,tmp);           
				}                           
			}                            
			else
			{
				//DO NOTHING!!!
			}

			//verity stop condition
			if ( nfeval > EVO_STEP-1 ) goto stop;

		}//END FOR    End mutation loop through pop.

		/*----Compute the energy variance (just for monitoring purposes)-----------*/

stop:

		//printf("DEPOIS POP:%d nfeval: %ld \n", cur_POP_SIZE[0], nfeval);

		cmean = 0.;          /* compute the mean value first */
		for (j=0; j<cur_POP_SIZE[0]; j++)
		{
			cmean += Eco[sp_index].fo[j];
		}
		cmean = cmean/(double)cur_POP_SIZE[0];

		cvar = 0.;           /* now the variance              */
		for (j=0; j<cur_POP_SIZE[0]; j++)
		{
			cvar += (Eco[sp_index].fo[j] - cmean)*(Eco[sp_index].fo[j] - cmean);
		}
		cvar = cvar/(double)(cur_POP_SIZE[0]-1);
	}//end WHILE

	//find best
	for(i=0;i<cur_POP_SIZE[0];i++)
	{
		if (Eco[sp_index].fo[i]<Eco[sp_index].bestfo[0])
		{
			Eco[sp_index].bestfo[0]=Eco[sp_index].fo[i];
			for(j=0;j<DIM;j++)
				Eco[sp_index].best[j]=Eco[sp_index].pop[i][j];

			Eco[sp_index].best_index[0] = i;
		}
	}


	/*----Output part----------------------------------------------------------*/
	/*	printf("\n\n                         PRESS ANY KEY TO ABORT"); 
		printf("\n\n\n Best-so-far cost funct. value=%-15.10g\n",Eco[sp_index].bestfo[0]);

		for (j=0;j<DIM;j++)
		{
		printf("\n best[%d]=%-15.10g",j,Eco[0].best[j]);
		}
		printf("\n\n Generation=%d  NFEs=%ld   Strategy: %s    ",gen,nfeval,strat[strategy]);
		printf("\n POP_SIZE=%d    F=%-4.2g    CR=%-4.2g   cost-variance=%-10.10g\n",
		POP_SIZE,F,CR,cvar);
		*/
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
	nfeval++;
	//			printf("\nCi from sp %d: %.2f\n",sp_index,Eco[sp_index].CiFO);

	free(solution);
	free(tmp);
	//   return(0);
}
/*-----------End of DE()------------------------------------------*/

void  assignd(int D, double a[], double b[])
	/**C*F****************************************************************
	 **                                                                  **
	 ** Assigns D-dimensional vector b to vector a.                      **
	 ** You might encounter problems with the macro ASSIGND on some      **
	 ** machines. If yes, better use this function although it's slower. **
	 **                                                                  **
	 ***C*F*E*************************************************************/
{
	int j;
	for (j=0; j<D; j++)
	{
		a[j] = b[j];
	}
}
