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
#include "eco.h"

//Functions declarations
void AvgStdDev(double *Avg,double *StdDev,double Var[]);
double randon( double inferior, double superior);
void print_eco(int suc);
char* my_itoa(int n, char str[]);
double box_muller(double m, double s);

void *threadExecute(void *threadid);

//Functions

double box_muller(double m, double s)	/* normal random variate generator */
{				        /* mean m, standard deviation s */
	double x1, x2, w, y1;
	static double y2;
	static int use_last = 0;

	if (use_last)		        /* use value from previous call */
	{
		y1 = y2;
		use_last = 0;
	}
	else
	{
		do {
			x1 = 2.0 * randon(0.0,1.0) - 1.0;
			x2 = 2.0 * randon(0.0,1.0) - 1.0;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );

		w = sqrt( (-2.0 * log( w ) ) / w );
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}

	return( m + y1 * s );
}

void AvgStdDev(double *Avg,double *StdDev,double Var[])
{
	int i;

	*Avg = 0;
	*StdDev = 0;
	for (i=0;i<RUN;i++)
		*Avg += Var[i];
	*Avg /= RUN;

	for (i=0;i<RUN;i++)
	{
		*StdDev += pow((*Avg-Var[i]),2);
	}
	*StdDev /= RUN;
   *StdDev = sqrt(*StdDev);
}

double randon( double inferior, double superior)
{
  double aux = (float)inferior + ((superior - inferior)*rand()/(RAND_MAX+1.0));
//  double aux = (float)inferior + ((superior - inferior)*MT_randInt(RAND_MAX)/(RAND_MAX+1.0));
 
  return aux;
}

char* my_itoa(int n, char str[])
{
  int i = 0;               /* Loop counter              */
  int negative = 0;        /* Indicate negative integer */
  int length = 0;          /* Length of string          */
  int temp = 0;            /* Temporary storage         */
 
  if(negative = (n<0))     /* Is it negative?  */
    n = -n;                /* make it positive */

  /* Generate digit characters in reverse order */
  do
  {
    str[i++] = '0'+n%10;    /* Create a rightmost digit        */
    n /= 10;                /* Remove the digit                */              
  }while(n>0);              /* Go again if there's more digits */

  if(negative)              /* If it was negative */
    str[i++] = '-';         /* Append minus       */
  str[i] = '\0';            /* Append terminator  */
  length = i;               /* Save the length    */

  /* Now reverse the string in place */
  /* by switching first and last,    */
  /* second and last but one, etc    */
  for(i = 0 ; i<length/2 ;i++)
  {
    temp = str[i];
    str[i] = str[length-i-1];
    str[length-i-1] = temp;
  }
  return str;                /* Return the string */
}

void init_eco()
{
	int i, j, k;
	double *mean; //mean[DIM];
	//init eco

	mean = (double*)malloc(DIM * sizeof(double));

	for (i=0;i<SP_NUMBER;i++)//each specie
	{
		for (k=0;k<DIM;k++)		
			mean[k] = randon(lb,ub);
		for (j=0;j<POP_SIZE*2;j++)//each individual
		{
			Eco[i].fo[j]  = 0.0;
			Eco[i].fit[j] = 0.0;
			for (k=0; k<DIM;k++) //each dimension
			{
					Eco[i].pop[j][k] = randon(lb,ub); //box_muller(mean[k],(abs(lb)+abs(ub))/10.0);

				if (Eco[i].pop[j][k] < lb)
					Eco[i].pop[j][k] = lb;
				if (Eco[i].pop[j][k] > ub)
					Eco[i].pop[j][k] = ub;
			}
		}
		for (k=0; k<DIM;k++) //each dimension
		{
			Eco[i].best[k] = 0.0;
		}
		Eco[i].bestfo[0] = 0.0;
		Eco[i].CiFO[0] = 0.0;
		func_eval[i] = 0.0;
	}
	free(mean);
}

void create_habitats_cluster(int r)
{
//********************************Single-link Clustering
		int i,j;
  		const int nrows = SP_NUMBER;
		const int ncols = DIM;
		double** data = malloc(nrows*sizeof(double*) );
		int** mask = malloc(nrows*sizeof(int*));
		double** distMatrix;
		double* weight = malloc(ncols*sizeof(double));
  		int nnodes = nrows-1;
		Node* tree;
//Variables used to perform scalonation of single-link distances
		double maxrange, minrange;
		char *file_name, *file_name_aux;
		file_name = (char *) malloc(100);
		file_name_aux = (char *) malloc(100);

		for (i = 0; i < nrows; i++)
		{ data[i] = malloc(ncols*sizeof(double));
		  mask[i] = malloc(ncols*sizeof(int));
  		}

		for (i=0;i<SP_NUMBER;i++)
  		{
		  	for (j=0;j<DIM;j++)
  			{
				data[i][j]=Eco[i].Ci[j];
				mask[i][j]=1;//there is no missing data
			}
  		}

  		distMatrix = malloc(nrows*sizeof(double*));
		if(distMatrix==NULL) return;// NULL; /* Not enough memory available */
		distMatrix[0] = NULL;
		  /* The zeroth row has zero columns. We allocate it anyway for convenience.*/
		for (i = 1; i < nrows; i++)
		{ 	distMatrix[i] = malloc(i*sizeof(double));
    			if (distMatrix[i]==NULL) break; /* Not enough memory available */
  		}

		if (i < nrows) /* break condition encountered */
  		{ 
			j = i;
    			for (i = 1; i < j; i++) free(distMatrix[i]);
  		}

		for (i = 0; i < ncols; i++) weight[i] = 1.0;

		for (i = 1; i < nrows; i++)
			for (j = 0; j < i; j++)
				distMatrix[i][j] = normalized_distance_matrix[i][j];

		tree = treecluster(nrows, ncols, 0, 0, 0, 0, 'e', 's', distMatrix); //<<<<<<<<< generates the single-link TREE!!!

		  /* The distance matrix was modified by treecluster, so we cannot use it any
		   * more. But we still need to deallocate it here.
		   * The first row of distmatrix is a single null pointer; no need to free it.
		   */

		for (i = 1; i < nrows; i++) free(distMatrix[i]);
		  free(distMatrix);
		if (!tree)
		{ /* Indication that the treecluster routine failed */
		    printf ("treecluster routine failed due to insufficient memory\n");
		    free(weight);
		    return;
		}
		//============== Scalonate the distances present in the single-link tree. These distances are used as probabilities.
		maxrange = -1;
		for (i=0;i<SP_NUMBER-1;i++)
			if (maxrange < tree[i].distance)
				maxrange = tree[i].distance;
		minrange = maxrange;
		for (i=0;i<SP_NUMBER-1;i++)
			if (minrange > tree[i].distance)
				minrange = tree[i].distance;
		for (i=0;i<SP_NUMBER-1;i++)
			tree[i].distance=(MAXVAL-MINVAL)/(double)(maxrange-minrange)*tree[i].distance-(MAXVAL-MINVAL)/(double)(maxrange-minrange)*minrange+MINVAL;	
		//============== 
	if (REPORT != 0)
	{
		FILE* fp;
  		if (r == 0 && (suc[0] == 0 || suc[0] == ECO_STEP-1) )
  		{
			strcpy(file_name,"");
			my_itoa( FUNCTION ,file_name_aux);
			strcat(file_name,"./report/f");
			strcat(file_name,file_name_aux);
			strcat(file_name,"-s");
			my_itoa ( STRATEGY ,file_name_aux);
			strcat(file_name,file_name_aux);
			strcat(file_name,"-run1-single-linkage-clustering.txt");
	   	
			fp = fopen( file_name, "a" ); //append to the file.
	  		fprintf(fp,"\n");
 	 		fprintf(fp,"MATLAB: ==== Pairwise single linkage clustering using NORMALIZED DIST MATRIX.\n");

   			fprintf(fp,"Suc: %d   Item 1   Item 2    Distance\n", suc[0]);
   			for(i=0; i<nnodes; i++)
   			{
				if (tree[i].left >= 0 && tree[i].right >= 0)
				    fprintf(fp,"%9d%9d      %g\n",tree[i].left+1, tree[i].right+1, tree[i].distance);
				if (tree[i].left >= 0 && tree[i].right < 0)
				    fprintf(fp,"%9d%9d      %g\n",tree[i].left+1, (-1*tree[i].right)+SP_NUMBER, tree[i].distance);
				if (tree[i].left < 0 && tree[i].right >= 0)
				    fprintf(fp,"%9d%9d      %g\n",(-1*tree[i].left)+SP_NUMBER, tree[i].right+1, tree[i].distance);
				if (tree[i].left < 0 && tree[i].right < 0)
				    fprintf(fp,"%9d%9d      %g\n",(-1*tree[i].left)+SP_NUMBER, (-1*tree[i].right)+SP_NUMBER, tree[i].distance);
			   }
			fprintf(fp,"\n");
		   	fclose(fp);
		}//end IF
	}//#endif REPORT

//********************NEW: CREATE HABITATS

		int aux_i, aux_j;
		double sum_next, aux_rand;
		double sum_left, sum_right;
		int cur_node, cur_node_index;
		int cur_habitat;

		int sum_tabu = 0;
		aux_i = aux_j = 0;
		sum_left = sum_right = sum_next = aux_rand = 0.0;
		h_count = 0;
		cur_node = SP_NUMBER-1;//starts with the highest node. Top-down.
		cur_node_index = 0;
		cur_habitat = h_count; //where h_count == 0
		while (sum_tabu < SP_NUMBER-1)
		{
			aux_rand = randon(0.0,(double)1.0);
			if ( aux_rand >= tree[cur_node-1].distance) //link those two itens
			{
				if (cur_node == SP_NUMBER-1) //first iteration.
				{
//printf("\nJOIN: First iteration...\n");
					H[cur_habitat].h_sp[ H[cur_habitat].h_sp_count ] = tree[cur_node-1].left;
					H[cur_habitat].h_sp[ H[cur_habitat].h_sp_count+1 ] = tree[cur_node-1].right;
					H[cur_habitat].h_sp_count = 2;
					sum_tabu++;
				}else
				{
					if (cur_node == 1)
					{
//printf("\nJOIN...\n");
						H[cur_habitat].h_sp[ cur_node_index ] = tree[0].left;
						H[cur_habitat].h_sp[ H[cur_habitat].h_sp_count ] = tree[0].right;
						H[cur_habitat].h_sp_count += 1;
						sum_tabu++;
					}
					else
					{
//printf("\nJOIN...\n");
						H[cur_habitat].h_sp[ cur_node_index ] = tree[cur_node-1].left;
						H[cur_habitat].h_sp[ H[cur_habitat].h_sp_count ] = tree[cur_node-1].right;
						H[cur_habitat].h_sp_count += 1;
						sum_tabu++;
					}
				}
			}else//split in two habitats
			{	//Then decides which item stays. The cosest item stays and replaces the cur_node value. 
				//The other creates a new habitat!
				if (cur_node == SP_NUMBER-1) //first iteration. One item stays and other item creates a new habitat.
				{
//printf("\nSPLIT: First iteration...\n");
					H[cur_habitat].h_sp[ H[cur_habitat].h_sp_count ] = tree[cur_node-1].left;
					h_count++;
					H[h_count].h_sp[ H[h_count].h_sp_count ] = tree[cur_node-1].right;
					H[cur_habitat].h_sp_count += 1;
					H[h_count].h_sp_count += 1;
					sum_tabu++;
				}else
				{
					if (H[cur_habitat].h_sp_count == 1 && H[cur_habitat].h_sp[0] < 0)
					//means that there is only one negative item inside the habitat.
					{
//printf("\nSPLIT: Habitat with only one negative...\n");
						H[cur_habitat].h_sp[ cur_node_index ] = tree[cur_node-1].left;
						h_count++;
						H[h_count].h_sp[ H[h_count].h_sp_count ] = tree[cur_node-1].right;
						//H[cur_habitat].h_sp_count += 1;
						H[h_count].h_sp_count += 1;
						sum_tabu++;
					}
					if (tree[cur_node-1].left >= 0 && tree[cur_node-1].right >= 0 && H[cur_habitat].h_sp_count > 1) 
					//means that both are positives
					//Here we use the normalized distance matrix.
					{
//printf("\nSPLIT: Both positives...\n");
						//calcula a distância acumulada dos dois itens para todos os itens positivos dentro do habitat.
						sum_left = sum_right = 0.0;
						for (i = 0; i < H[cur_habitat].h_sp_count; i++)
						{
						   if (H[cur_habitat].h_sp[ i ] >= 0)
						   {
							sum_left  += normalized_distance_matrix[ H[cur_habitat].h_sp[ i ] ][ tree[cur_node-1].left ];
							sum_right += normalized_distance_matrix[ H[cur_habitat].h_sp[ i ] ][ tree[cur_node-1].right ];
						   }
						}
						if (sum_left <= sum_right) //left stays and right goes out
						{
							H[cur_habitat].h_sp[ cur_node_index ] = tree[cur_node-1].left;
							h_count++;
							H[h_count].h_sp[ H[h_count].h_sp_count ] = tree[cur_node-1].right;
						}else//left goes out and right stays
						{
							H[cur_habitat].h_sp[ cur_node_index ] = tree[cur_node-1].right;
							h_count++;
							H[h_count].h_sp[ H[h_count].h_sp_count ] = tree[cur_node-1].left;
						}
						//H[cur_habitat].h_sp_count += 1; Doesn't update because the chosen item replaces cur_node_index.
						H[h_count].h_sp_count += 1;
						sum_tabu++;
					}else if (tree[cur_node-1].left < 0 && tree[cur_node-1].right < 0 && H[cur_habitat].h_sp_count > 1)
					//means that both are negative. The closest item is the lowest item.
						{
//printf("\nSPLIT: Both negatives...\n");
						if ( tree[cur_node-1].left < tree[cur_node-1].right )//left stays and right goes out. 
						//Means that it is closest to the cur_node.
						{
							H[cur_habitat].h_sp[ cur_node_index ] = tree[cur_node-1].left;
							h_count++;
							H[h_count].h_sp[ H[h_count].h_sp_count ] = tree[cur_node-1].right;
						}else//right stays and left goes out.
						{
							H[cur_habitat].h_sp[ cur_node_index ] = tree[cur_node-1].right;
							h_count++;
							H[h_count].h_sp[ H[h_count].h_sp_count ] = tree[cur_node-1].left;
						}					
						//H[cur_habitat].h_sp_count += 1; Doesn't update because the chosen item replaces cur_node_index.
						H[h_count].h_sp_count += 1;
						sum_tabu++;
					}else if (H[cur_habitat].h_sp_count > 1) //(tree[cur_node-1].left < 0 || tree[cur_node-1].right < 0). 
					//Means that at least one of the itens are negative.
					{
//printf("\nSPLIT: Positive and negative...\n");
						//The closest item is the non-negative item.
						if (tree[cur_node-1].left >= 0) //left stays and right goes out
						{
							H[cur_habitat].h_sp[ cur_node_index ] = tree[cur_node-1].left;
							h_count++;
							H[h_count].h_sp[ H[h_count].h_sp_count ] = tree[cur_node-1].right;
						}else//left goes out and left stays
						{
							H[cur_habitat].h_sp[ cur_node_index ] = tree[cur_node-1].right;
							h_count++;
							H[h_count].h_sp[ H[h_count].h_sp_count ] = tree[cur_node-1].left;
						}
						//H[cur_habitat].h_sp_count += 1; Doesn't update because the chosen item replaces cur_node_index.
						H[h_count].h_sp_count += 1;
						sum_tabu++;
					}//end ELSE
				}//end ELSE
			}//end ELSE

			//update cur_habitat
			//selects the first habitat with a negative item inside.

			i = 0;
			cur_habitat = -1;
			j = 0;
			while (i <= h_count && cur_habitat < 0)
			{
				while (j < H[i].h_sp_count && cur_habitat < 0)
				{
					if (H[i].h_sp[j] < 0)
					{
						cur_habitat = i;
					}
					else j++;				
				}
				i++;
				j = 0;
			}

			//update 'cur_node' and 'cur_node_index'. Find the first negative item inside cur_habitat.
		        if (cur_habitat >= 0)
			{
				i = 0;
				while (H[cur_habitat].h_sp[i] >= 0)
					i++;
				cur_node = -1*H[cur_habitat].h_sp[i];
				cur_node_index = i;
			}
		}//end WHILE

		free(tree);

		for (i = 0; i < nrows; i++)
		{ free(data[i]);
		  free(mask[i]);
  		}
		free(data);
		free(mask);
		free(weight);
		free(file_name_aux);
		free(file_name);

///end*******************NEW: CREATE HABITATS
}

void print_eco(int suc)
{
	int i, j, k;
	//print eco
////Population, variables, centroids.
	for (i=0;i<SP_NUMBER;i++)
	{
		printf("\nPop %d \n",i);
		for (j=0;j<cur_POP_SIZE[0];j++)
		{
			printf("Variables: ");		
			for (k=0; k<DIM;k++) //variables
			{
				printf("%7.4f ",Eco[i].pop[j][k]);
			}
			printf(" Fo e Fit:");
			printf("%8.4f %7.4f \n",Eco[i].fo[j],Eco[i].fit[j]);
		}
		printf("\nBest:\n");
		for(k=0;k<DIM;k++)
		{
			printf("GlobalParam[%d]: %7.4f\n",k+1,Eco[i].best[k]);
		}
		printf("Fo: ");
		printf("%7.4f \n",Eco[i].bestfo[0]);
		printf("\nCentroid:\n");
		for(k=0;k<DIM;k++)
		{
			printf("Ci[%d]: %7.4f \n",k+1,Eco[i].Ci[k]);
		}
		printf("CiFO: ");
		printf("%7.4f \n",Eco[i].CiFO[0]);
	}
}




/*Main program of the ECO algorithm*/
int main(int argc, char **argv)
{

	FILE *fp, *frun;
	char *file_name, *file_name_aux;
	time_t t;
	double tempo_execucao;
	int i, j, k, r,
	    int_sp, //number of interations in species phase
	    int_hab,//number of interations in habitats phase
	    out_solution, //variable that contains a random solution index to be sent.
	    in_solution, //variable that contains a random solution index to be replaced.
	    out_specie, //variable that contains a random specie index to sent a solution.
	    in_specie, //variable that contains a random specie index to receive a solution.
	    in_habitat; //variable that contains a random habitat index to receive a solution.

	int returnInfo; // info from threads

	int best_t;
	double best_valor;
	double *solution; //solution[DIM];
	double aux1, aux2;
	int tag;

	double sum, *euclid_diff, //euclid_diff[DIM]
	       euclid_sum, euclid_max, euclid_mean, euclid_min, *sol, //sol[DIM],
	       avg, std, *MeanVector, //MeanVector[RUN], 
	       *BestVector, //BestVector[RUN],
	       *eco_stat; //eco_stat[ECO_STEP];

	//Statistics
	double AvgBest, 	//contains the avg value of the best result obtained by each specie.
	       GlobalBest,	//contains the global best value.
	       **GlobalParams; //GlobalParams[RUN][DIM]; //contains the global best parameters.

 	double dist_aux; //auxiliar variable to form adjacent list

	srand(time(NULL));
	MT_seed();

	if (GetParameters(argv) == -1)	//read input file
		return 0;
	showParameters();

	
	//the number of search strategies employed
	STRATEGY_NUMBER = 1; 
	if (STRATEGY == 3) //ABC-PSO
		STRATEGY_NUMBER = 2;
	if (STRATEGY == 5) //All DE
		STRATEGY_NUMBER = 10;
	if (STRATEGY == 7) //ABC-PSO-jDE/BBO
		STRATEGY_NUMBER = 3;
	if (STRATEGY == 8) //ABC-PSO-DE-jDE/BBO
		STRATEGY_NUMBER = 4;

	AllocArrays();			//alloc arrays
	solution = (double*)malloc(DIM * sizeof(double));
	euclid_diff = (double*)malloc(DIM * sizeof(double));
	sol = (double*)malloc(DIM * sizeof(double));
	MeanVector = (double*)malloc(RUN * sizeof(double));
	BestVector = (double*)malloc(RUN * sizeof(double));
	GlobalParams = malloc (RUN * sizeof(double*));
	for (i = 0; i < RUN; i++)
		GlobalParams[i] = (double*)malloc (DIM * sizeof(double));

	eco_stat = (double*)malloc(ECO_STEP * sizeof(double));

	if (FUNCTION == 21) //3D-AB
	{
		if (ProteinSize < 1)
		{
			printf("Info: SEQUENCE is needed for 3D-AB.\nInfo: For 3D-AB, DIM = 2*ProteinSize-5 i.e.: 13 aminoacids --> DIM = 21\n");
			return 0;
		}
		for(i=0;i<ProteinSize;i++)
	   	{
			if (SEQUENCE[i] == 'A' || SEQUENCE[i] == 'a')
				Aminoacido[i].Tipo = 1;
			if (SEQUENCE[i] == 'B' || SEQUENCE[i] == 'b')
				Aminoacido[i].Tipo = -1;
	   	}
	}

	if (FUNCTION == 23) //2D-AB
	{
		if (ProteinSize < 1)
		{
			printf("Info: SEQUENCE is needed for 2D-AB.\nInfo: For 2D-AB, DIM ProteinSize-2 i.e.: 13 aminoacids --> DIM = 11\n");
			return 0;
		}
		for(i=0;i<ProteinSize;i++)
	   	{
			if (SEQUENCE[i] == 'A' || SEQUENCE[i] == 'a')
				Aminoacido[i].Tipo = 1;
			if (SEQUENCE[i] == 'B' || SEQUENCE[i] == 'b')
				Aminoacido[i].Tipo = -1;
	   	}
	}

	pthread_attr_init(&pthread_custom_attr);

	file_name = (char *) malloc(100);
	file_name_aux = (char *) malloc(100);

	prepararObjFunc();

	strcpy(file_name,"");
	my_itoa ( FUNCTION ,file_name_aux);
	strcat(file_name,"./report/f");
	strcat(file_name,file_name_aux);
	strcat(file_name,"-s");
	my_itoa ( STRATEGY ,file_name_aux);
	strcat(file_name,file_name_aux);
	strcat(file_name,"-optimization-result.txt");

   frun = fopen( file_name, "w" ); //create the file.
   fprintf(frun,"Runs: %d \nSpecies: %d \nCandidate solutions in each specie: %d \nEcological succession steps: %d \nEvolutionary steps (gen/cycles): %d\n", RUN,SP_NUMBER, POP_SIZE, ECO_STEP, EVO_STEP);
   fclose(frun);

	strcpy(file_name,"");
	my_itoa ( FUNCTION ,file_name_aux);
	strcat(file_name,"./report/f");
	strcat(file_name,file_name_aux);
	strcat(file_name,"-s");
	my_itoa ( STRATEGY ,file_name_aux);
	strcat(file_name,file_name_aux);
	strcat(file_name,"-best-result.txt");

   frun = fopen(file_name, "w" ); //create the file.
   fprintf(frun,"<Run> <GlobalFo>");
   fclose(frun);

   for (i=0;i<ECO_STEP;i++)
	eco_stat[i] = 0.0;
   printf("<Run> <GlobalFo>");
   for (r=0;r<RUN;r++)
   {
	strcpy(file_name,"");
	my_itoa ( FUNCTION ,file_name_aux);
	strcat(file_name,"./report/f");
	strcat(file_name,file_name_aux);
	strcat(file_name,"-s");
	my_itoa ( STRATEGY ,file_name_aux);
	strcat(file_name,file_name_aux);
	strcat(file_name,"-optimization-result.txt");

	frun = fopen( file_name, "a" ); //create the file.
   	fprintf(frun,"\nRun: %d \n", r);
   	fclose(frun);
	
	tempo_execucao = (unsigned int) time(&t);
	ElapsedTime[r] = 0.0;
	suc[0] = 0;
	init_eco();

	cur_POP_SIZE[0] = POP_SIZE;

//******************* REPORT
if (REPORT != 0)
{
 //FILE	//CREATE FILES
 if (r == 0) //record data for the first run
 {
	strcpy(file_name,"");
	my_itoa ( FUNCTION ,file_name_aux);
	strcat(file_name,"./report/f");
	strcat(file_name,file_name_aux);
	strcat(file_name,"-s");
	my_itoa ( STRATEGY ,file_name_aux);
	strcat(file_name,file_name_aux);
	strcat(file_name,"-run1-single-linkage-clustering.txt");

	fp = fopen( file_name, "w" ); //create the file.
	fclose(fp);

		//Format the name of a file to record the best individual of each specie in each ITERATION.
		strcpy(file_name,"");
		my_itoa ( FUNCTION ,file_name_aux);
		strcat(file_name,"./report/f");
		strcat(file_name,file_name_aux);
		strcat(file_name,"-s");
		my_itoa ( STRATEGY ,file_name_aux);
		strcat(file_name,file_name_aux);
//		my_itoa ( i ,file_name_aux);
		strcat(file_name,"-data-best-sp");
//		strcat(file_name,file_name_aux);
		strcat(file_name,".txt");
		fp = fopen( file_name, "w" ); //create the file.
		fprintf( fp, "# iteration | Best \n" );
		fclose(fp);

	for (i=0;i<SP_NUMBER;i++)
	{
	    if (REPORT == 12)
	    {
		//Format the name of a file to record the centroid of each specie in each ECO_STEP.
		strcpy(file_name,"");
		my_itoa ( FUNCTION ,file_name_aux);
		strcat(file_name,"./report/f");
		strcat(file_name,file_name_aux);
		strcat(file_name,"-s");
		my_itoa ( STRATEGY ,file_name_aux);
		strcat(file_name,file_name_aux);
		my_itoa ( i ,file_name_aux);
		strcat(file_name,"-data-centroid-sp-");
		strcat(file_name,file_name_aux);
		strcat(file_name,".txt");
		fp = fopen( file_name, "w" ); //create the file.
		fprintf( fp, "# eco-step | centroid: x0 .. xd fo \n" );
		fclose(fp);
		//Format the name of a file to record the population of each specie in each ECO_STEP.
		for (j=0;j<ECO_STEP;j++)
		{
		   if (j < 3 || (j == ECO_STEP-2) || (j == ECO_STEP-1))
	   	   {
			strcpy(file_name,"");
			my_itoa ( FUNCTION ,file_name_aux);
			strcat(file_name,"./report/f");
			strcat(file_name,file_name_aux);
			strcat(file_name,"-s");
			my_itoa ( STRATEGY ,file_name_aux);
			strcat(file_name,file_name_aux);
			my_itoa ( i ,file_name_aux);			
			strcat(file_name,"-data-pop-sp-");
			strcat(file_name,file_name_aux);
			my_itoa ( j ,file_name_aux);
			strcat(file_name,"-eco-");
			strcat(file_name,file_name_aux);
			strcat(file_name,".txt");
			strcpy(file_name_aux,"");
			fp = fopen( file_name, "w" ); //create the file.
			fprintf( fp, "# x0 .. xd fo \n" );
			fclose(fp);
		   }
		}
	    }//end REPORT == 12
	}//END FOR

			for (i=0;i<SP_NUMBER;i++)
			{
				for (j=0;j<cur_POP_SIZE[0];j++)
				{
					for(k=0;k<DIM;k++){
						if (Eco[i].pop[j][k]<lb)
					           Eco[i].pop[j][k]=lb;
					        if (Eco[i].pop[j][k]>ub)
					           Eco[i].pop[j][k]=ub;
						sol[k] = Eco[i].pop[j][k];
					}
					Eco[i].fo[j] = function(sol,i);
				}
			}//END FOR

		//CENTROID POP INICIAL	//DEFINE Ci value for each specie. Place where a specie is concentrated.
		//In this implementation, the Ci value is the population centroid.
		for (i=0;i<SP_NUMBER;i++) 
		{	
			for (k=0; k<DIM;k++)//each variable
			{
				sum = 0.0;
				for (j=0;j<cur_POP_SIZE[0];j++) //each solution
				{
					sum += Eco[i].pop[j][k];
				}
				Eco[i].Ci[k] = sum/(double)cur_POP_SIZE[0];
			}
			Eco[i].CiFO[0] = function(Eco[i].Ci,i);
		}

     if (REPORT == 12)
     {
	for (j=0;j<ECO_STEP;j++)
	{
	   if (j < 3 || (j == ECO_STEP-2) || (j == ECO_STEP-1))
	   {
		strcpy(file_name,"");
		my_itoa ( FUNCTION ,file_name_aux);
		strcat(file_name,"./report/f");
		strcat(file_name,file_name_aux);
		strcat(file_name,"-s");
		my_itoa ( STRATEGY ,file_name_aux);
		strcat(file_name,file_name_aux);
		strcat(file_name,"-data-centroid-all-sp-eco-");
		my_itoa ( j ,file_name_aux);
		strcat(file_name,file_name_aux);
		strcat(file_name,".txt");
		strcpy(file_name_aux,"");
		fp = fopen( file_name, "w" ); //create the file.
		fprintf( fp, "# sp | centroid: x0 .. xd fo \n" );
		fclose(fp);
  	   }
	}

	strcpy(file_name,"");
	my_itoa ( FUNCTION ,file_name_aux);
	strcat(file_name,"./report/f");
	strcat(file_name,file_name_aux);
	strcat(file_name,"-s");
	my_itoa ( STRATEGY ,file_name_aux);
	strcat(file_name,file_name_aux);
	strcat(file_name,"-data-centroid-all-sp-eco-ini.txt");
	strcpy(file_name_aux,"");
	fp = fopen( file_name, "w" ); //create the file.
	fprintf( fp, "# sp | centroid: x0 .. xd fo \n" );
	fclose(fp);


//FILE Pop initial
			for (i=0;i<SP_NUMBER;i++)
			{
				strcpy(file_name,"");
				my_itoa ( FUNCTION ,file_name_aux);
				strcat(file_name,"./report/f");
				strcat(file_name,file_name_aux);
				strcat(file_name,"-s");
				my_itoa ( STRATEGY ,file_name_aux);
				strcat(file_name,file_name_aux);
				my_itoa ( i ,file_name_aux); //FILE
				strcat(file_name,"-data-pop-sp-"); //FILE
				strcat(file_name,file_name_aux); //FILE
				my_itoa ( (int)suc[0] ,file_name_aux); //FILE
				strcat(file_name,"-00"); //FILE
				strcat(file_name,".txt"); //FILE
				strcpy(file_name_aux,""); //FILE
				fp = fopen( file_name, "w" ); //create the file.
				for (j=0;j<cur_POP_SIZE[0];j++)
				{
					for(k=0;k<DIM;k++){
						fprintf( fp, "%10.4f\'",Eco[i].pop[j][k]);
					}
					fprintf( fp, "%10.4f\n",Eco[i].fo[j]); //FILE
				}
				fclose(fp); //FILE
			}//END FOR

		//CENTROID INITIAL POP	//DEFINE Ci value for each specie. Place where a specie is concentrated.
		//In this implementation, the Ci value is the population centroid.
		for (i=0;i<SP_NUMBER;i++) 
		{	
				//Format the name of a file to record the CENTROID of each specie in each ECO_STEP.
				strcpy(file_name,"");
				my_itoa ( FUNCTION ,file_name_aux);
				strcat(file_name,"./report/f");
				strcat(file_name,file_name_aux);
				strcat(file_name,"-s");
				my_itoa ( STRATEGY ,file_name_aux);
				strcat(file_name,file_name_aux);
				my_itoa ( i ,file_name_aux); //FILE
				strcat(file_name,"-data-centroid-sp-"); //FILE
				strcat(file_name,file_name_aux); //FILE
				strcat(file_name,".txt"); //FILE
				fp = fopen( file_name, "a" ); //open the file.	
				fprintf( fp, "ini\'"); //FILE
			for (k=0; k<DIM;k++)//each variable
			{
				fprintf( fp, "%10.4f\'",Eco[i].Ci[k]); //FILE
			}
			fprintf( fp, "%10.4f\n",Eco[i].CiFO[0]); //FILE			
			fclose(fp); //close file.

			strcpy(file_name,"");
			my_itoa ( FUNCTION ,file_name_aux);
				strcat(file_name,"./report/f");
			strcat(file_name,file_name_aux);
			strcat(file_name,"-s");
			my_itoa ( STRATEGY ,file_name_aux);
			strcat(file_name,file_name_aux);
			strcat(file_name,"-data-centroid-all-sp-eco-ini.txt"); //FILE
			strcpy(file_name_aux,""); //FILE
			fp = fopen( file_name, "w" ); //open the file.
			fprintf( fp, "%d\'",i); //FILE
			for (k=0; k<DIM;k++)//each variable
			{
				fprintf( fp, "%10.4f\'",Eco[i].Ci[k]); //FILE
			}
			fprintf( fp, "%10.4f\n",Eco[i].CiFO[0]); //FILE			
			fclose(fp); //close file.
			
			strcpy(file_name,"");
			my_itoa ( FUNCTION ,file_name_aux);
			strcat(file_name,"./report/f");
			strcat(file_name,file_name_aux);
			strcat(file_name,"-s");
			my_itoa ( STRATEGY ,file_name_aux);
			strcat(file_name,file_name_aux);
			my_itoa ( i ,file_name_aux); //FILE
			strcat(file_name,"-data-centroid-all-sp-eco-ini-"); //FILE
			strcat(file_name,file_name_aux); //FILE
			strcat(file_name,".txt"); //FILE
			strcpy(file_name_aux,""); //FILE
			fp = fopen( file_name, "w" ); //open the file.
			fprintf( fp, "%d\'",i); //FILE
			for (k=0; k<DIM;k++)//each variable
			{
				fprintf( fp, "%10.4f\'",Eco[i].Ci[k]); //FILE
			}
			fprintf( fp, "%10.4f\n",Eco[i].CiFO[0]); //FILE			
			fclose(fp); //close file.

		}
     } //end REPORT == 12

	strcpy(file_name,"");
	my_itoa ( FUNCTION ,file_name_aux);
	strcat(file_name,"./report/f");
	strcat(file_name,file_name_aux);
	strcat(file_name,"-s");
	my_itoa ( STRATEGY ,file_name_aux);
	strcat(file_name,file_name_aux);
	strcat(file_name,"-run1-data-habitats-eco.txt");
	strcpy(file_name_aux,"");
	fp = fopen( file_name, "w" ); //create the file.
	fprintf( fp, "# eco-step | number of habitats \n" );
	fclose(fp);

		//Euclidan Distance among centroids.
		suc[0] = 0;
		h_count = 0;
		euclid_max = 0.0;
		euclid_mean = 0.0;

		for (i=0;i<SP_NUMBER;i++)
		{

			for (j=0;j<SP_NUMBER;j++)
			{
				for (k=0; k<DIM;k++)//each variable
				{
					euclid_diff[k] = (Eco[i].Ci[k] - Eco[j].Ci[k])*(Eco[i].Ci[k] - Eco[j].Ci[k]);
				}	
				for (k=0; k<DIM;k++)//each variable
				{
					euclid_sum += euclid_diff[k];
					euclid_diff[k] = 0.0;
				}
				distance_matrix[i][j] = sqrt(euclid_sum);

				euclid_mean += distance_matrix[i][j];

				if (euclid_max < distance_matrix[i][j])
					euclid_max = distance_matrix[i][j];

				euclid_sum = 0.0;
			}
		}
		mean_euclid_dist[suc[0]] = euclid_mean/(2*SP_NUMBER);

		for (i=0;i<SP_NUMBER;i++)
		{
			for (j=0;j<SP_NUMBER;j++)
			{
				normalized_distance_matrix[i][j] = distance_matrix[i][j]/euclid_max;
				euclid_mean += normalized_distance_matrix[i][j];
			}
		}

		//INIT ADJACENCY LIST
		for (i=0;i<SP_NUMBER;i++)
		{
			sp_int[i] = 0;
			tabu_list[i] = 0; //initiate the tabu list for species.
			next_sp[i] = 0; //initiate the variable that aux the habitats creation.
		}
		for (i=0;i<SP_NUMBER;i++)
		{
			for (j=0;j<SP_NUMBER;j++)
			{
				sp_adj[i][j] = -1;
				H[i].h_sp[j] = 0; //INIT HABITAT
			}
			H[i].h_sp_count = 0; //INIT HABITAT COUNT
		}
		h_count = 0;

		create_habitats_cluster(r);

		h_count++;
///begin*******************NEW: ADJACENCT LIST
		for (i=0;i<h_count;i++)
		{
			if (H[i].h_sp_count == 1)//means that there is only one population inside habitat 'i'
			{
				//update nothing
			}else if ((H[i].h_sp_count == 2))//means that there are two populations inside habitat 'i'. 
							 //One population comunicates with each other.
			{
					sp_adj[ H[i].h_sp[0] ][0] = H[i].h_sp[1];
					sp_int[H[i].h_sp[0]]++;

					sp_adj[ H[i].h_sp[1] ][0] = H[i].h_sp[0];
					sp_int[H[i].h_sp[1]]++;
			}else// means that there are more than two populations.
			     //each population choses probabilisticaly according to the normalized distance matrix one population to communicate.
			{
				for (j=0;j<H[i].h_sp_count;j++)
				{
					k = (int)randon(0,H[i].h_sp_count);
					while (sp_int[H[i].h_sp[j]] == 0)
					{
						if (normalized_distance_matrix[ H[i].h_sp[j] ][ H[i].h_sp[k] ] == 1.0)//give a small chance to choose
							dist_aux = 0.99;
						else dist_aux = normalized_distance_matrix[ H[i].h_sp[j] ][ H[i].h_sp[k] ];

						if ( (j != k) && (randon(0.0,(double)1.0) >= dist_aux) )
//						if ( (j != k) && (randon(0.0,(double)1.0) >= normalized_distance_matrix[ H[i].h_sp[j] ][ H[i].h_sp[k] ]) )
						{
							sp_adj[ H[i].h_sp[j] ][sp_int[H[i].h_sp[j]]] = H[i].h_sp[k];
							sp_int[H[i].h_sp[j]]++;
							k = H[i].h_sp_count;
						}else k = (int)randon(0,H[i].h_sp_count);
					}
				}
			}//end ELSE
		}
  ///end*******************NEW: ADJACENCY LIST


//FILE Habiatas and AdjList
			strcpy(file_name,"");
			my_itoa ( FUNCTION ,file_name_aux);
			strcat(file_name,"./report/f");
			strcat(file_name,file_name_aux);
			strcat(file_name,"-s");
			my_itoa ( STRATEGY ,file_name_aux);
			strcat(file_name,file_name_aux);
			strcat(file_name,"-run1-data-habitat-adj-eco-ini");
			strcpy(file_name_aux,"");
			fp = fopen( file_name, "w" ); //create the file.
			////Print the habitats at each ecological step
			fprintf(fp,"Habitats:\n");
			for (i=0;i<h_count;i++)
			{	
				fprintf(fp,"H(%d): ",i);
				for (j=0;j<H[i].h_sp_count;j++)
					fprintf(fp,"%4d ",H[i].h_sp[j]);
				fprintf(fp,"\n");
			}		
			////Adjacency list
			fprintf(fp,"\nAdjacency list:\n");
			for (i=0;i<SP_NUMBER;i++)
			{
				fprintf(fp,"Sp-%d: ",i);
				for (j=0;j<sp_int[i];j++)
				{
					fprintf(fp,"%3d ",sp_adj[i][j]);
				}
				fprintf(fp,"\n");
			}
			fclose(fp);

			strcpy(file_name,"");
			my_itoa ( FUNCTION ,file_name_aux);
			strcat(file_name,"./report/f");
			strcat(file_name,file_name_aux);
			strcat(file_name,"-s");
			my_itoa ( STRATEGY ,file_name_aux);
			strcat(file_name,file_name_aux);
			strcat(file_name,"-run1-data-habitats-eco.txt");
			strcpy(file_name_aux,"");
			fp = fopen( file_name, "a" ); //create the file.
			fprintf( fp, "ini' %d\n",h_count);
			fclose(fp);


 }//end IF record for r==0.
}//end //report
//*******************REPORT

	if(STRATEGY == 50)//PRVNS
	{
		prepare_prvns();
	}

	//ECOLOGICAL SUCCESSION loop
	h_count = 0;
	printf("\n%d ",r); //printf RUN

	strcpy(file_name,"");
	my_itoa ( FUNCTION ,file_name_aux);
	strcat(file_name,"./report/f");
	strcat(file_name,file_name_aux);
	strcat(file_name,"-s");
	my_itoa ( STRATEGY ,file_name_aux);
	strcat(file_name,file_name_aux);
	strcat(file_name,"-best-result.txt");
	frun = fopen(file_name,"a");
	fprintf(frun,"\n%d ",r); //printf RUN
	fclose(frun);

	/*All initial solutions are evaluated */
	for (k=0;k<SP_NUMBER;k++)
	{
		for(i=0;i<cur_POP_SIZE[0];i++)
		{
			for (j=0;j<DIM;j++)
   			{
				if (Eco[k].pop[i][j] < lb)
					Eco[k].pop[i][j] = lb;
				if (Eco[k].pop[i][j] > ub)
					Eco[k].pop[i][j] = ub;
			   	solution[j]=Eco[k].pop[i][j];
   			}
   			Eco[k].fo[i]=function(solution, k);
			//fitness
		 	if(Eco[k].fo[i]>=0)
		 	{
				 Eco[k].fit[i]=1/(double)(Eco[k].fo[i]+1.0);
		 	}
		 	else
		 	{
			 Eco[k].fit[i]=1+fabs(Eco[k].fo[i]);
		 	}
		}
		/*The best solution is memorized*/
		Eco[k].bestfo[0]=Eco[k].fo[0];
        	for(i=0;i<DIM;i++)
		        Eco[k].best[i]=Eco[k].pop[0][i];
		for(i=0;i<cur_POP_SIZE[0];i++)
		{
			if (Eco[k].fo[i]<Eco[k].bestfo[0])
			{
        			Eco[k].bestfo[0]=Eco[k].fo[i];
        			for(j=0;j<DIM;j++)
        			   Eco[k].best[j]=Eco[k].pop[i][j];
				Eco[k].best_index[0] = i;
        		}
		}
	}//END FOR

		strcpy(file_name,"");
		my_itoa ( FUNCTION ,file_name_aux);
		strcat(file_name,"./report/f");
		strcat(file_name,file_name_aux);
		strcat(file_name,"-s");
		my_itoa ( STRATEGY ,file_name_aux);
		strcat(file_name,file_name_aux);
		strcat(file_name,"-optimization-result.txt");
		frun = fopen( file_name, "a" ); //create the file.
		if (frun == NULL) {
			perror("fopen failed");  // print a meaningful error message telling why fopen faile
    			return 1;  // return 1 from main, signifying error and stopping your program
		}

		fprintf(frun,"Eco-Step: ");
		fclose(frun);


	while (suc[0] < ECO_STEP) 
	{
		strcpy(file_name,"");
		my_itoa ( FUNCTION ,file_name_aux);
		strcat(file_name,"./report/f");
		strcat(file_name,file_name_aux);
		strcat(file_name,"-s");
		my_itoa ( STRATEGY ,file_name_aux);
		strcat(file_name,file_name_aux);
		strcat(file_name,"-optimization-result.txt");
		frun = fopen( file_name, "a" ); //create the file.
		if (frun == NULL) {
			perror("fopen failed");  // print a meaningful error message telling why fopen faile
    			return 1;  // return 1 from main, signifying error and stopping your program
		}
		fprintf(frun," %d ",suc[0]); //printf eco_step
		fclose(frun);
	
		int_hab = 0;
		int_sp = 0;
		//INIT ADJACENCY LIST
		for (i=0;i<SP_NUMBER;i++)
		{
			sp_int[i] = 0;
			tabu_list[i] = 0; //initiate the tabu list for species.
			next_sp[i] = 0; //initiate the variable that aux the habitats creation.
		}
		for (i=0;i<SP_NUMBER;i++)
		{
			for (j=0;j<SP_NUMBER;j++)
			{
				sp_adj[i][j] = -1;
				H[i].h_sp[j] = 0; //INIT HABITAT
			}
			H[i].h_sp_count = 0; //INIT HABITAT COUNT
		}
			

		//EVOLUTIVE PERIOD for each specie. Number of generations/cycles.
		for (i=0;i<SP_NUMBER;i++) 
		{
			returnInfo = pthread_create(&threads[i], &pthread_custom_attr, threadExecute, ((void *) ((long) i)));
		        if (returnInfo)
		        {
		            printf("ERROR; return code from pthread_create() is %d\ninitialize matrices\n", returnInfo);
		            return -1;
		        }
		}
		for (i = 0; i < SP_NUMBER; i++)
        	{	
            		pthread_join(threads[i], NULL);
        	}


		//Compute the FO mean of the Centroids.
		CcFO = 0.0;
		sum = 0.0;
		for (j=0;j<SP_NUMBER;j++) //each solution
		{
if (REPORT == 12)
{
 if (r == 0) //record data for the first run
 {
				//Format the name of a file to record the CENTROID of each specie in each ECO_STEP.
				strcpy(file_name,"");
				my_itoa ( FUNCTION ,file_name_aux);
				strcat(file_name,"./report/f");
				strcat(file_name,file_name_aux);
				strcat(file_name,"-s");
				my_itoa ( STRATEGY ,file_name_aux);
				strcat(file_name,file_name_aux);
				my_itoa ( j ,file_name_aux); //FILE
				strcat(file_name,"-data-centroid-sp-"); //FILE
				strcat(file_name,file_name_aux); //FILE
				strcat(file_name,".txt"); //FILE
				fp = fopen( file_name, "a" ); //open the file.	
				fprintf( fp, "%d\'",suc[0]); //FILE
			for (k=0; k<DIM;k++)//each variable
			{
				fprintf( fp, "%10.4f\'",Eco[j].Ci[k]); //FILE
			}
			fprintf( fp, "%10.4f\n",Eco[j].CiFO[0]); //FILE			
			fclose(fp); //close file.

		   if (suc[0] < 3 || (suc[0] == ECO_STEP - 2) || (suc[0] == ECO_STEP-1))
		   {
			strcpy(file_name,"");
			my_itoa ( FUNCTION ,file_name_aux);
			strcat(file_name,"./report/f");
			strcat(file_name,file_name_aux);
			strcat(file_name,"-s");
			my_itoa ( STRATEGY ,file_name_aux);
			strcat(file_name,file_name_aux);
			strcat(file_name,"-data-centroid-all-sp-eco-"); //FILE
			my_itoa ( suc[0] ,file_name_aux); //FILE
			strcat(file_name,file_name_aux); //FILE
			strcat(file_name,".txt"); //FILE
			strcpy(file_name_aux,""); //FILE
			fp = fopen( file_name, "a" ); //open the file.
			fprintf( fp, "%d\'",j); //FILE
			for (k=0; k<DIM;k++)//each variable
			{
				fprintf( fp, "%10.4f\'",Eco[j].Ci[k]); //FILE

			}
			fprintf( fp, "%10.4f\n",Eco[j].CiFO[0]); //FILE			
			fclose(fp); //close file.
		   }
 }//end IF r==0.
}//end //report
			sum += Eco[j].CiFO[0];
		}
		CcFO = sum/(double)SP_NUMBER;

		//Euclidan Distance among centroids.
		h_count = 0;
		euclid_max = 0.0;
		euclid_mean = 0.0;
		euclid_min = 0.0;
		euclid_sum = 0.0;

if (NOHABITAT)
{
// NO interation among species!
}
else 
{
		for (i=0;i<SP_NUMBER;i++)
		{
			for (j=0;j<SP_NUMBER;j++)
			{
				for (k=0; k<DIM;k++)//each variable
				{
					euclid_diff[k] = pow (fabs(Eco[i].Ci[k] - Eco[j].Ci[k]), 2);
				}	
				for (k=0; k<DIM;k++)//each variable
				{
					euclid_sum += euclid_diff[k];
					euclid_diff[k] = 0.0;
				}
				distance_matrix[i][j] = (sqrt(euclid_sum));

				if (euclid_max < distance_matrix[i][j])
					euclid_max = distance_matrix[i][j];

				euclid_sum = 0.0;
			}
		}

		euclid_min = euclid_max;
		for (i=0;i<SP_NUMBER;i++)
		{
			for (j=0;j<SP_NUMBER;j++)
			{
				if (euclid_min > distance_matrix[i][j] && i!=j)
					euclid_min = distance_matrix[i][j];
			}
		}

		for (i=0;i<SP_NUMBER;i++)
		{
			for (j=0;j<SP_NUMBER;j++)
			{
				euclid_mean += distance_matrix[i][j];
			}
		}
		euclid_max = 0.0;
		for (i=0;i<SP_NUMBER;i++)
		{
			for (j=0;j<SP_NUMBER;j++)
			{
				if (euclid_max < distance_matrix[i][j] && i!=j)
					euclid_max = distance_matrix[i][j];
			}
		}

		mean_euclid_dist[suc[0]] = euclid_mean/SP_NUMBER;


		for (i=0;i<SP_NUMBER;i++)
		{
			for (j=0;j<SP_NUMBER;j++)
			{
				normalized_distance_matrix[i][j] = distance_matrix[i][j]/euclid_max;
				euclid_mean += normalized_distance_matrix[i][j];
			}
		}

		//INIT ADJACENCY LIST
		for (i=0;i<SP_NUMBER;i++)
		{
			sp_int[i] = 0;
			tabu_list[i] = 0; //initiate the tabu list for species.
			next_sp[i] = 0; //initiate the variable that aux the habitats creation.
		}
		for (i=0;i<SP_NUMBER;i++)
		{
			for (j=0;j<SP_NUMBER;j++)
			{
				sp_adj[i][j] = -1;
			}
			H[i].h_sp_count = 0; //INIT HABITAT COUNT
		}
		h_count = 0;

		create_habitats_cluster(r);

		h_count++;
///begin*******************NEW: ADJACENCT LIST
		for (i=0;i<h_count;i++)
		{
			if (H[i].h_sp_count == 1)//means that there is only one population inside habitat 'i'
			{
				//update nothing
			}else if ((H[i].h_sp_count == 2))//means that there are two populations inside habitat 'i'. 
							 //One population comunicates with each other.
			{
					sp_adj[ H[i].h_sp[0] ][0] = H[i].h_sp[1];
					sp_int[H[i].h_sp[0]]++;

					sp_adj[ H[i].h_sp[1] ][0] = H[i].h_sp[0];
					sp_int[H[i].h_sp[1]]++;
			}else// means that there are more than two populations.
			     //each population choses probabilisticaly according to the normalized distance matrix one population to communicate.
			{
				for (j=0;j<H[i].h_sp_count;j++)
				{
					k = (int)randon(0,H[i].h_sp_count);
					while (sp_int[H[i].h_sp[j]] == 0)
					{
						if (normalized_distance_matrix[ H[i].h_sp[j] ][ H[i].h_sp[k] ] == 1.0)//give a small chance to choose
							dist_aux = 0.99;
						else dist_aux = normalized_distance_matrix[ H[i].h_sp[j] ][ H[i].h_sp[k] ];

						if ( (j != k) && (randon(0.0,(double)1.0) >= dist_aux) )
						{
							sp_adj[ H[i].h_sp[j] ][sp_int[H[i].h_sp[j]]] = H[i].h_sp[k];
							sp_int[H[i].h_sp[j]]++;
							k = H[i].h_sp_count;
						}else k = (int)randon(0,H[i].h_sp_count);
					}
				}
			}//end ELSE
		}
///end*******************NEW: ADJACENCY LIST

if (REPORT != 0)
{
 if (r == 0) //record data for the first run
 {
		strcpy(file_name,"");
		my_itoa ( FUNCTION ,file_name_aux);
		strcat(file_name,"./report/f");
		strcat(file_name,file_name_aux);
		strcat(file_name,"-s");
		my_itoa ( STRATEGY ,file_name_aux);
		strcat(file_name,file_name_aux);
		strcat(file_name,"-run1-data-habitats-eco.txt");
		strcpy(file_name_aux,"");
		fp = fopen( file_name, "a" ); //create the file.
		fprintf( fp, "%d\'%d\n",suc[0],h_count);
		fclose(fp);
 }//end IF r==0.
}//end report

		eco_stat[suc[0]] += h_count;

		//PERFORM species communications for each topology TQj.
		//According to the ADJACENCY LIST, each specie communicates with each other.
		for (i=0;i<SP_NUMBER;i++)
		{
			if (sp_int[i] == 0) //there is no specie to communicate. This means that specie i is an habitat by itself.
			{ } //do nothing
			else
			{   
				//'out_': origin 'in_':destiny
				//Choose a solution from 'i' to reproduce based on a tournament selection.
				//TOURNEY
			        for (j = 0; j < T_SIZE; j++)
				{
					k = (int)randon(0.0,(double)cur_POP_SIZE[0]); //representa uma outra solução escolhida eleatoriamente
					tourney_index[j] = k;
				}
				//seleciona a melhor solução da população torneio tendo como base a FObj
				best_valor = Eco[ i ].fo[tourney_index[0]];
 			        best_t = 0;
				for(j=0;j<T_SIZE;j++)
				{
					if(Eco[ i ].fo[tourney_index[j]] < best_valor)
				  	{
					       best_t = j;
					       best_valor = Eco[ i ].fo[tourney_index[j]];
					}
				}
				out_solution = tourney_index[best_t];

				//Choose a random adjacent specie from sp i.
				in_specie = (int)randon(0,sp_int[i]);

				//Choose a solution from 'in_specie' to reproduce based on a tournament selection.
				//TOURNEY
			        for (j = 0; j < T_SIZE; j++)
				{
					k = (int)randon(0.0,(double)cur_POP_SIZE[0]); //representa uma outra solução escolhida eleatoriamente
					tourney_index[j] = k;
				}
				//seleciona a melhor solução da população torneio tendo como base a FObj
				best_valor = Eco[ sp_adj[i][in_specie] ].fo[tourney_index[0]];
 			        best_t = 0;
				for(j=0;j<T_SIZE;j++)
				{
					if(Eco[ sp_adj[i][in_specie] ].fo[tourney_index[j]] < best_valor)
				  	{
					       best_t = j;
					       best_valor = Eco[ sp_adj[i][in_specie] ].fo[tourney_index[j]];
					}
				}
				in_solution = tourney_index[best_t];

				//PERFORM species communications for each topology TQj.
				//Generate a new solution.
				for(k=0;k<DIM;k++){
					if (randon(0.0,(double)1.0) <= .5)
					{
						sol[k] = Eco[i].pop[out_solution][k];
					}					
					else
					{
						sol[k] = Eco[ sp_adj[i][in_specie] ].pop[in_solution][k];
					}
				}
				
				//Choose a random solution to be replaced.
				in_solution = (int)randon(0,cur_POP_SIZE[0]);
				while (in_solution == Eco[ sp_adj[i][in_specie] ].best_index[0]) //The best solution can't be replaced. 
					in_solution = (int)randon(0,cur_POP_SIZE[0]);

				for (k = 0;k< DIM;k++)
				{
					Eco[ sp_adj[i][in_specie] ].pop[in_solution][k] = sol[k];
					if (Eco[ sp_adj[i][in_specie] ].pop[in_solution][k] < lb)
						Eco[ sp_adj[i][in_specie] ].pop[in_solution][k] = lb;
					if (Eco[ sp_adj[i][in_specie] ].pop[in_solution][k] > ub)
						Eco[ sp_adj[i][in_specie] ].pop[in_solution][k] = ub;
				//For PSO only. Pbest receive current new solution.
					Pbest[ sp_adj[i][in_specie] ].x[in_solution][k] = sol[k];
				//For PSO only. The velocities are renewed at random.
					Vel[sp_adj[i][in_specie]].v[in_solution][k] = (alea(lb,ub)-Pbest[ sp_adj[i][in_specie] ].x[in_solution][k])/2.0;
				}

				Eco[ sp_adj[i][in_specie] ].fo[in_solution] = function(sol, sp_adj[i][in_specie] );

				if(Eco[sp_adj[i][in_specie]].fo[in_solution]>=0)
			 	{
					Eco[sp_adj[i][in_specie]].fit[in_solution]=1/(double)(Eco[sp_adj[i][in_specie]].fo[in_solution]+1.0);
			 	}
			 	else
			 	{
					Eco[sp_adj[i][in_specie]].fit[in_solution]=1+fabs(Eco[sp_adj[i][in_specie]].fo[in_solution]);
			 	}

			//For PSO. Pbest fo is updated.
				Pbest[ sp_adj[i][in_specie] ].f[in_solution] = Eco[ sp_adj[i][in_specie] ].fo[in_solution];
				int_sp++;
			}//end else

		}//END FOR

 		//DEFINE communication topology TH among habitats. << Communication inter habitats!
		if (h_count == 1) //There is only one habitat.
		{  } //nothing to be done
		else
		{
		   for (i=0;i<h_count;i++)
	           {
			if (H[i].h_sp_count == 1) //there is only one specie inside the habitat.
			{  
				//Choose a specie inside the habitat i.
				out_specie = H[i].h_sp[0];
				//The best solution of 'out_specie' is selected to migrate.
				out_solution = Eco[ out_specie ].best_index[0];
				//Choose a random habitat to receive a solution. 
				in_habitat = (int)randon(0,h_count);
				while (in_habitat == i) //The selected habitat MUST be different from itself. 'in_habitat != i'
					in_habitat = (int)randon(0,h_count);
				//Choose a random specie inside 'in_habitat' to receive a solution. 
				in_specie = H[in_habitat].h_sp[ (int)randon(0,H[in_habitat].h_sp_count) ];
				//Choose a random solution inside 'in_specie' to be replaced. 
				in_solution = (int)randon(0,cur_POP_SIZE[0]);
				while (in_solution == Eco[ in_specie ].best_index[0]) //The best solution can't be replaced. 
					in_solution = (int)randon(0,cur_POP_SIZE[0]);

				//PERFORM habitats communications for topology TH.
				for (k = 0;k<DIM;k++)
				{
					Eco[in_specie].pop[in_solution][k] = Eco[out_specie].pop[out_solution][k];
				//For PSO. Pbest migrated with the particle.
					Pbest[ in_specie ].x[in_solution][k] = Pbest[out_specie].x[out_solution][k];
				//For PSO. The velocities are renewed at random.
					Vel[in_specie].v[in_solution][k] = Vel[out_specie].v[out_solution][k];//(alea(lb,ub)-Pbest[ in_specie ].x[in_solution][k])/2.0;
				}

				Eco[ in_specie ].fo[in_solution] = Eco[ out_specie ].fo[out_solution];
				Eco[ in_specie ].fit[in_solution] = Eco[ out_specie ].fit[out_solution];

			//For PSO. Pbest fo is updated.
				Pbest[ in_specie ].f[in_solution] = Pbest[ out_specie ].f[out_solution];


				int_hab++;

 			}//end if
			else //there are more than one species inside the habitat.
			{	
					//Choose a random specie inside the habitat i.
					out_specie = H[i].h_sp[ (int)randon(0,H[i].h_sp_count) ];
					//The best solution of 'out_specie' is selected to migrate.
					out_solution = Eco[ out_specie ].best_index[0];
					//Choose a random habitat to receive a solution. 
					in_habitat = (int)randon(0,h_count);
					while (in_habitat == i) //The selected habitat MUST be different from itself. 'in_habitat != i'
						in_habitat = (int)randon(0,h_count);
					//Choose a random specie inside 'in_habitat' to receive a solution. 
					in_specie = H[in_habitat].h_sp[ (int)randon(0,H[in_habitat].h_sp_count) ];
					//Choose a random solution inside 'in_specie' to be replaced. 
					in_solution = (int)randon(0,cur_POP_SIZE[0]);
					while (in_solution == Eco[ in_specie ].best_index[0]) //The best solution can't be replaced. 
						in_solution = (int)randon(0,cur_POP_SIZE[0]);
	
					//PERFORM habitats communications for topology TH.
					for (k = 0;k<DIM;k++)
					{
						Eco[in_specie].pop[in_solution][k] = Eco[out_specie].pop[out_solution][k];
					//For PSO. Pbest migrated with the particle.
						Pbest[ in_specie ].x[in_solution][k] = Pbest[out_specie].x[out_solution][k];
					//For PSO. The velocities are renewed at random.
						Vel[in_specie].v[in_solution][k] = Vel[out_specie].v[out_solution][k];//(alea(lb,ub)-Pbest[ in_specie ].x[in_solution][k])/2.0;
					}
	
					Eco[ in_specie ].fo[in_solution] = Eco[ out_specie ].fo[out_solution];
					Eco[ in_specie ].fit[in_solution] = Eco[ out_specie ].fit[out_solution];

				//For PSO. Pbest fo is updated.
					Pbest[ in_specie ].f[in_solution] = Pbest[ out_specie ].f[out_solution];


					int_hab++;
			}//end else
		   }//end for
		}//end else
} //end else //NOHABITAT		

		/*The current best solution is updated*/
		for (k=0;k<SP_NUMBER;k++)
		{
			Eco[k].bestfo[0]=Eco[k].fo[0];
        		for(i=0;i<DIM;i++)
			        Eco[k].best[i]=Eco[k].pop[0][i];
			for(i=0;i<cur_POP_SIZE[0];i++)
			{
				if (Eco[k].fo[i]<Eco[k].bestfo[0]) // <<== minimization problem!!!
				{
        				Eco[k].bestfo[0]=Eco[k].fo[i];
        				for(j=0;j<DIM;j++)
        				   Eco[k].best[j]=Eco[k].pop[i][j];
					Eco[k].best_index[0] = i;
        			}
			}
		}//END FOR 
		suc[0]++;

if (REPORT != 0)
{
 if (r == 0) //record data for the first run
 {
		if (suc[0] < 3 || (suc[0] == ECO_STEP-2) || (suc[0] == ECO_STEP-1))
	   	{
		   if (REPORT == 12)
		   {
			//FILE Pop
			for (i=0;i<SP_NUMBER;i++)
			{
				strcpy(file_name,"");
				my_itoa ( FUNCTION ,file_name_aux);
				strcat(file_name,"./report/f");
				strcat(file_name,file_name_aux);
				strcat(file_name,"-s");
				my_itoa ( STRATEGY ,file_name_aux);
				strcat(file_name,file_name_aux);
				my_itoa ( i ,file_name_aux);		 //FILE	
				strcat(file_name,"-data-pop-sp-"); //FILE
				strcat(file_name,file_name_aux); //FILE
				my_itoa ( suc[0] ,file_name_aux); //FILE
				strcat(file_name,"-eco-"); //FILE
				strcat(file_name,file_name_aux); //FILE
				strcat(file_name,".txt"); //FILE
				strcpy(file_name_aux,""); //FILE
				fp = fopen( file_name, "w" ); //create the file.
				for (j=0;j<cur_POP_SIZE[0];j++)
				{
					for(k=0;k<DIM;k++){
						fprintf( fp, "%10.4f\'",Eco[i].pop[j][k]); //FILE
					}
					fprintf( fp, "%10.4f\n",Eco[i].fo[j]); //FILE
				}
				fclose(fp); //FILE
			}
		  }//end REPORT == 12
			//FILE Habiatas and AdjList
			strcpy(file_name,"");
			my_itoa ( FUNCTION ,file_name_aux);
			strcat(file_name,"./report/f");
			strcat(file_name,file_name_aux);
			strcat(file_name,"-s");
			my_itoa ( STRATEGY ,file_name_aux);
			strcat(file_name,file_name_aux);
			strcat(file_name,"-data-habitat-adj-eco-");
			my_itoa ( suc[0] ,file_name_aux);
			strcat(file_name,file_name_aux);
			strcat(file_name,".txt");
			strcpy(file_name_aux,"");
			fp = fopen( file_name, "w" ); //create the file.
			////Print the habitats at each ecological step
			fprintf(fp,"Habitats:\n");
			for (i=0;i<h_count;i++)
			{	
				fprintf(fp,"H(%d): ",i);
				for (j=0;j<H[i].h_sp_count;j++)
					fprintf(fp,"%4d ",H[i].h_sp[j]);
				fprintf(fp,"\n");
			}		
			////Adjacency list
			fprintf(fp,"\nAdjacency list:\n");
			for (i=0;i<SP_NUMBER;i++)
			{
				fprintf(fp,"Sp-%d: ",i);
				for (j=0;j<sp_int[i];j++)
				{
					fprintf(fp,"%3d ",sp_adj[i][j]);
				}
				fprintf(fp,"\n");
			}
			fclose(fp);
		}//end IF record
				//FILE best of each sp in each eco loop
				strcpy(file_name,"");
				my_itoa ( FUNCTION ,file_name_aux);
				strcat(file_name,"./report/f");
				strcat(file_name,file_name_aux);
				strcat(file_name,"-s");
				my_itoa ( STRATEGY ,file_name_aux);
				strcat(file_name,file_name_aux);
				strcat(file_name,"-data-best-sp");
				strcat(file_name,".txt");
				strcpy(file_name_aux,"");
			fp = fopen( file_name, "a" ); //append to the file.
			fprintf(fp,"%d \t",suc[0]);
			for (i=0;i<SP_NUMBER;i++)
			{
				////Print the habitats at each ecological step
				fprintf(fp,"%.5g \t",Eco[i].bestfo[0]);
			}
			fprintf(fp,"\n");
			fclose(fp);
 }//end IF r==0.
} //end //REPORT

	}//END WHILE ECO_STEP

	//FIND GLOBAL BEST
	GlobalBest = Eco[0].bestfo[0];
	AvgBest = 0.0;
	for(i=0;i<SP_NUMBER;i++)
	{
		AvgBest += Eco[i].bestfo[0];
		if (Eco[i].bestfo[0]<=GlobalBest) // <<== minimization problem!!!
		{
        		GlobalBest=Eco[i].bestfo[0];
        		for(j=0;j<DIM;j++)
        		   GlobalParams[r][j]=Eco[i].best[j];
        	}
	}
	AvgBest /= (float)SP_NUMBER;

	printf("%g",GlobalBest);

 	strcpy(file_name,"");
	my_itoa ( FUNCTION ,file_name_aux);
	strcat(file_name,"./report/f");
	strcat(file_name,file_name_aux);
	strcat(file_name,"-s");
	my_itoa ( STRATEGY ,file_name_aux);
	strcat(file_name,file_name_aux);
	strcat(file_name,"-best-result.txt"); //FILE
	frun = fopen(file_name, "a" ); //create the file.
   	fprintf(frun,"%g ",GlobalBest);

	MeanVector[r] = AvgBest;
	BestVector[r] = GlobalBest;

	if ( (FUNCTION == 21) || (FUNCTION == 23) )
	{
		for(k=0;k<DIM;k++)
		{
			fprintf(frun,"%g ",GlobalParams[r][k]);
		}
		fprintf(frun,"Fo: %g \n",BestVector[r]);
	}
	fclose(frun);

	tempo_execucao = time(&t) - tempo_execucao;

	ElapsedTime[r] = tempo_execucao;

   }//END FOR r

	//FIND GLOBAL BEST
	int aux;
	GlobalBest = BestVector[0];
	for(i=0;i<RUN;i++)
	{
		if (BestVector[i]<=GlobalBest) // <<== minimization problem!!!
		{
        		GlobalBest=BestVector[i];
			aux = i;
        	}
	}
	strcpy(file_name,"");
	my_itoa ( FUNCTION ,file_name_aux);
	strcat(file_name,"./report/f");
	strcat(file_name,file_name_aux);
	strcat(file_name,"-s");
	my_itoa ( STRATEGY ,file_name_aux);
	strcat(file_name,file_name_aux);
	strcat(file_name,"-overall-best-solution.txt"); //FILE
	frun = fopen( file_name, "w" ); //create the file.
	for(k=0;k<DIM;k++)
	{
		fprintf(frun,"%g ",GlobalParams[aux][k]);
	}
	fprintf(frun,"Fo: %g \n",GlobalBest);
   	fclose(frun);

if (REPORT != 0)
{
 	strcpy(file_name,"");
	my_itoa ( FUNCTION ,file_name_aux);
	strcat(file_name,"./report/f");
	strcat(file_name,file_name_aux);
	strcat(file_name,"-s");
	my_itoa ( STRATEGY ,file_name_aux);
	strcat(file_name,file_name_aux);
	strcat(file_name,"-data-habitats-run.txt"); //FILE
   fp = fopen( file_name, "w" ); //create the file.
	for (i=0;i<ECO_STEP;i++)
		fprintf(fp,"%d\'%10.4f\n",i,eco_stat[i]/(float)RUN);
   fclose(fp);
}

	printf("\nRuns: %d \nSpecies: %d \nCandidate solutions in each specie: %d \nEcological succession steps: %d \nEvolutionary steps (eval/ECO_STEP): %d\n", RUN,SP_NUMBER, POP_SIZE, ECO_STEP, EVO_STEP);

   sum = 0.0;
   for (i=0;i<SP_NUMBER;i++)
	sum += func_eval[i];
   printf("Mean Function Evaluations for each Specie: %.2f\n",sum/(double)SP_NUMBER);

 	strcpy(file_name,"");
	my_itoa ( FUNCTION ,file_name_aux);
	strcat(file_name,"./report/f");
	strcat(file_name,file_name_aux);
	strcat(file_name,"-s");
	my_itoa ( STRATEGY ,file_name_aux);
	strcat(file_name,file_name_aux);
	strcat(file_name,"-optimization-result.txt"); //FILE
   frun = fopen( file_name, "a" ); //append to file.

   fprintf(frun,"\n===============\n");
   printf("===============\n");

   AvgStdDev(&avg,&std,BestVector);
   fprintf(frun," Avg of All Species Global Best: %g +/- %7.4f\n===============\n",avg,std);
   printf(" Avg of All Species Global Best: %g +/- %7.4f\n===============\n",avg,std);

   if ( (REPORT != 0) && (RUN > 1) )
   {
	   avg = 0;
	   std = 0;
	   for (i=1;i<RUN;i++)
		avg += ElapsedTime[i];
	   avg /= (RUN-1);

	   for (i=1;i<RUN;i++)
	   {
		std += pow((avg-ElapsedTime[i]),2);
	   }
	   std /= (RUN-1);
	   std = sqrt(std);
   }
   else {
	   AvgStdDev(&avg,&std,ElapsedTime);
   }

   printf("Mean Time (sec): %.2f +/- %.2f\n", avg, std);
   fprintf(frun,"Mean Time (sec): %.2f +/- %.2f\n", avg, std);

if (!NOHABITAT)
{
	printf("CLUSTER: Single Link Information.\n");
	fprintf(frun,"CLUSTER: Single Link Information.\n");
}

if (NOHABITAT)
{
	printf("NOHABITAT: enabled.\n");
	fprintf(frun,"NOHABITAT: enabled.\n");
}
else
{
	printf("NOHABITAT: NOT enabled.\n");
	fprintf(frun,"NOHABITAT: NOT enabled.\n");
}

if (REPORT != 0)
{
	printf("REPORT: enabled.\n");
	fprintf(frun,"REPORT: enabled.\n");
}
else
{
	printf("REPORT: NOT enabled.\n");
	fprintf(frun,"REPORT: NOT enabled.\n");
}


   fclose(frun);
   return 0;
}

void *threadExecute(void *threadid)
{
    int i = (int) ((long) threadid);

    	switch (STRATEGY)
	{
		case 1: //ABC
			ABC(i);
			break;
		case 2: //PSO
			PSO(i);
			break;
		case 3: //ABC-PSO
			if( ( i % STRATEGY_NUMBER ) == 0) 
			      ABC(i); //printf("Even");
			if( ( i % STRATEGY_NUMBER ) ==1) 
			      PSO(i); //printf("Odd");
			break;
		case 41: //DE
			DE(i, 1);
			break;
		case 42: //DE
			DE(i, 2);
			break;
		case 43: //DE
			DE(i, 3);
			break;
		case 44: //DE
			DE(i, 4);
			break;
		case 45: //DE
			DE(i, 5);
			break;
		case 46: //DE
			DE(i, 6);
			break;
		case 47: //DE
			DE(i, 7);
			break;
		case 48: //DE
			DE(i, 8);
			break;
		case 49: //DE
			DE(i, 9);
			break;
		case 40: //DE
			DE(i, 0);
			break;
		case 5: //All DE
			DE(i, i % STRATEGY_NUMBER);
			break;
		case 6: //jDE-BBO
			Run_jDE_BBO(i);
			break;
		case 7: //ABC-PSO-jDE/BBO
			if( ( i % STRATEGY_NUMBER ) == 0) 
			      ABC(i); //printf("Even");
			if( ( i % STRATEGY_NUMBER ) ==1) 
			      PSO(i); //printf("Odd");
			if( ( i % STRATEGY_NUMBER ) ==2) 
			      Run_jDE_BBO(i); //printf("Odd");
			break;
		case 8: //ABC-PSO-DE-jDE/BBO
			if( ( i % STRATEGY_NUMBER ) == 0) 
			      ABC(i);
			if( ( i % STRATEGY_NUMBER ) ==1) 
			      PSO(i);
			if( ( i % STRATEGY_NUMBER ) ==2) 
			      Run_jDE_BBO(i);
			if( ( i % STRATEGY_NUMBER ) ==3) 
			      DE(i, 2);
			break;
		case 50:
			PRVNS(i);
			break;

	}
	pthread_exit(NULL);
}

/* Dynamic array allocation */
/* mol and sites */
void AllocArrays()	
{
	int i,j;

	ElapsedTime = (double*) malloc (RUN * sizeof(double));
	tourney_index = (int*) malloc (T_SIZE * sizeof(int));
 	sp_int = (int*) malloc (SP_NUMBER * sizeof(int));

	H = (struct habitat*) malloc(SP_NUMBER * sizeof(struct habitat));
	for (i = 0; i < SP_NUMBER; i++)
		H[i].h_sp = (int*)malloc (sizeof(int)*SP_NUMBER );

	sp_adj = malloc (SP_NUMBER * sizeof(int*));
	for (i = 0; i < SP_NUMBER; i++)
		sp_adj[i] = (int*)malloc (SP_NUMBER * sizeof(int));

	tabu_list = (int*) malloc (SP_NUMBER * sizeof(int));
	threads = (pthread_t*) malloc (SP_NUMBER * sizeof(pthread_t));
	next_sp = (int*) malloc (SP_NUMBER * sizeof(int));

	distance_matrix = malloc (SP_NUMBER * sizeof(double*));
	for (i = 0; i < SP_NUMBER; i++)
		distance_matrix[i] = (double*)malloc (SP_NUMBER * sizeof(double));

	normalized_distance_matrix = malloc (SP_NUMBER * sizeof(double*));
	for (i = 0; i < SP_NUMBER; i++)
		normalized_distance_matrix[i] = (double*)malloc (SP_NUMBER * sizeof(double));

	func_eval = (int*) malloc (SP_NUMBER * sizeof(int));

	//only for PSO
	Vel = (struct velocity*) malloc(SP_NUMBER * sizeof(struct velocity));
	for (i = 0; i < SP_NUMBER; i++)
	{
	   Vel[i].v = malloc (POP_SIZE*2 * sizeof(double*));
	   for (j = 0;j < POP_SIZE*2; j++)
		Vel[i].v[j] = (double*)malloc (DIM * sizeof(double));
	}//for i
	//only for PSO
	Pbest = (struct pbest_position*) malloc(SP_NUMBER * sizeof(struct pbest_position));
	for (i = 0; i < SP_NUMBER; i++)
	{
	   Pbest[i].x = malloc (POP_SIZE*2 * sizeof(double*));
	   for (j = 0;j < POP_SIZE*2; j++)
		Pbest[i].x[j] = (double*)malloc (DIM * sizeof(double));
	   Pbest[i].f = (double*)malloc (POP_SIZE*2 * sizeof(double));
	}//for i

	Cc = (double*) malloc (DIM * sizeof(double));
	mean_euclid_dist = (double*) malloc (ECO_STEP * sizeof(double));

	Eco = (struct specie*) malloc(SP_NUMBER * sizeof(struct specie));
	for (i = 0; i < SP_NUMBER; i++)
	{
	   Eco[i].pop = malloc (POP_SIZE*2 * sizeof(double*));
	   for (j = 0;j < POP_SIZE*2; j++)
		Eco[i].pop[j] = (double*)malloc (DIM * sizeof(double));
	   Eco[i].fo = (double*)malloc (POP_SIZE*2 * sizeof(double));
	   Eco[i].fit = (double*)malloc (POP_SIZE*2 * sizeof(double));
	   Eco[i].best = (double*)malloc (DIM * sizeof(double));
	   Eco[i].Ci = (double*)malloc (DIM * sizeof(double));
	   Eco[i].trial = (int*)malloc (POP_SIZE * sizeof(int));

	   Eco[i].CRR = (double*)malloc (POP_SIZE * sizeof(double));
	   Eco[i].FF = (double*)malloc (POP_SIZE * sizeof(double));
	}//for i

	if ( (FUNCTION == 21) || (FUNCTION == 23) ) //3D-AB or 2D-AB
		Aminoacido = (amino *) malloc ((ProteinSize) * sizeof (amino));

}


/* Free arrays */
void freeArrays()
{
	int i,j;
	free(mean_euclid_dist);
	free(Cc);
	free(next_sp);
	free(func_eval);
	free(ElapsedTime);
	free(tourney_index);
	free(sp_int);
	free(tabu_list);
	for (i = 0; i < SP_NUMBER; i++)
		free(sp_adj[i]);
	free(sp_adj);
	for (i = 0; i < SP_NUMBER; i++)
		free(H[i].h_sp);
	free(H);
	free(threads);
	for (i = 0; i < SP_NUMBER; i++)
		free(distance_matrix[i]);
	free(distance_matrix);
	for (i = 0; i < SP_NUMBER; i++)
		free(normalized_distance_matrix[i]);
	free(normalized_distance_matrix);

	for (i = 0; i < SP_NUMBER; i++)
	{
	   for (j = 0;j < POP_SIZE*2; j++)
		free(Vel[i].v[j]);
	   free(Vel[i].v);
	}//for i
	free(Vel);

	for (i = 0; i < SP_NUMBER; i++)
	{
	   for (j = 0;j < POP_SIZE*2; j++)
		free(Pbest[i].x[j]);
	   free(Pbest[i].x);
	   free(Pbest[i].f);
	}//for i
	free(Pbest);

	for (i = 0; i < SP_NUMBER; i++)
	{
	   for (j = 0;j < POP_SIZE*2; j++)
		free(Eco[i].pop[j]);
	   free(Eco[i].pop[j]);
	   free(Eco[i].fo);
	   free(Eco[i].fit);
	   free(Eco[i].best);
	   free(Eco[i].Ci);
	   free(Eco[i].trial);
	}//for i
	free(Eco);

	if ( (FUNCTION == 21) || (FUNCTION == 23) ) //3D-AB or 2D-AB
		free(Aminoacido);

}


//
/*Input file reading*/
int GetParameters(char **argv)
{
	FILE *file = fopen( argv[1], "r" );
    
    if (file == 0)
    {
    	printf( "Could not open ini file! Usage ./<exec> <file.in>\n" );
	return -1;
    }
    else 
    {
		ffscanf("RUN", file, "%d", &RUN);
		ffscanf("SP_NUMBER", file, "%d", &SP_NUMBER);
		ffscanf("POP_SIZE", file, "%d", &POP_SIZE);
		ffscanf("ECO_STEP", file, "%d", &ECO_STEP);
		ffscanf("EVO_STEP", file, "%d", &EVO_STEP);
		ffscanf("T_SIZE", file, "%d", &T_SIZE);
		ffscanf("DIM", file, "%d", &DIM); 
		ffscanf("FUNCTION",file, "%d", &FUNCTION);
		ffscanf("NOHABITAT", file, "%d", &NOHABITAT);
		ffscanf("REPORT", file, "%d", &REPORT);
		ffscanf("STRATEGY", file, "%d", &STRATEGY);
		if ( (FUNCTION == 21) || (FUNCTION == 23) ) //3D-AB or 2D-AB
		{
			ffscanf("SEQUENCE", file, "%s", &SEQUENCE);
			ProteinSize = strlen(SEQUENCE);
		}
	return 1;
    }
    fclose(file);
}


void showParameters()
{
	printf("***PARAMETERS***\n");
	printf("RUNS = %d\n", RUN);
	printf("SP_NUMBER = %d\n", SP_NUMBER);
	printf("POP_SIZE = %d\n", POP_SIZE);
	printf("ECO_STEP = %d\n", ECO_STEP);
	printf("EVO_STEP = %d\n", EVO_STEP);
	printf("T_SIZE = %d\n", T_SIZE);
	printf("DIM = %d\n", DIM);
	if (NOHABITAT == 0 && SP_NUMBER <2)
	{
		printf("Info: For hierarchical clustering, SP_NUMBER must be greater than 2.\n") ;
		exit(0);
	}
	printf("NOHABITAT = %s\n",(NOHABITAT)?"YES":"NO");
	printf("REPORT = %s\n",(REPORT)?"YES":"NO");
	switch (STRATEGY)
	{
		case 1:
			printf("STRATEGY = %s\n","ABC");
			break;
		case 2:
			printf("STRATEGY = %s\n","PSO");
			break;
		case 3:
			printf("STRATEGY = %s\n","ABC-PSO");
			break;
		case 41:
			printf("STRATEGY = %s\n","DE/best/1/exp");
			break;
		case 42:
			printf("STRATEGY = %s\n","DE/rand/1/exp");
			break;
		case 43:
			printf("STRATEGY = %s\n","DE/rand-to-best/1/exp");
			break;
		case 44:
			printf("STRATEGY = %s\n","DE/best/2/exp");
			break;
		case 45:
			printf("STRATEGY = %s\n","DE/rand/2/exp");
			break;
		case 46:
			printf("STRATEGY = %s\n","DE/best/1/bin");
			break;
		case 47:
			printf("STRATEGY = %s\n","DE/rand/1/bin");
			break;
		case 48:
			printf("STRATEGY = %s\n","DE/rand-to-best/1/bin");
			break;
		case 49:
			printf("STRATEGY = %s\n","DE/best/2/bin");
			break;
		case 40:
			printf("STRATEGY = %s\n","DE/rand/2/bin");
			break;
		case 5:
			printf("STRATEGY = %s\n","DE/All");
			break;
		case 6:
			printf("STRATEGY = %s\n","jDE/BBO");
			break;
		case 7:
			printf("STRATEGY = %s\n","ABC-PSO-jDE/BBO");
			break;
		case 8:
			printf("STRATEGY = %s\n","ABC-PSO-DE/rand/1/exp-jDE/BBO");
			break;
		case 50:
			printf("STRATEGY = %s\n", "PRVNS");
			break;
			
		default:
		        printf("Info: Invalid strategy.\n") ;
		        exit(0);
	}
	switch (FUNCTION)
	{
		case 0:
			printf("FUNCTION = %s\n","Rastrigin");
			break;
		case 1:
			printf("FUNCTION = %s\n","Schaffer f7");
			break;
		case 2:
			printf("FUNCTION = %s\n","Griewank");
			break;
		case 3:
			printf("FUNCTION = %s\n","Ackley");
			break;
		case 4:
			printf("FUNCTION = %s\n","Rosenbrock");
			break;
		case 5:
			printf("FUNCTION = %s\n","Sphere");
			break;
		case 6:
			printf("FUNCTION = %s\n","StretchedV");
			break;
		case 7:
			printf("FUNCTION = %s\n","Schwefel's function 2.22");
			break;
		case 8:
			printf("FUNCTION = %s\n","Step function");
			break;
		case 9:
			printf("FUNCTION = %s\n","Generalized Schwefel's function");
			break;
		case 10:
			printf("FUNCTION = %s\n","Generalized Penalized #1");
			break;
		case 11:
			printf("FUNCTION = %s\n","Generalized Penalized #2");
			break;
		case 12:
			printf("FUNCTION = %s\n","Levy");
			break;
		case 13:
			printf("FUNCTION = %s\n","Zakharov");
			break;
		case 14:
			printf("FUNCTION = %s\n","Egg Holder");
			break;
		case 15:
			printf("FUNCTION = %s\n","Generalized Holzman");
			break;
		case 16:
			printf("FUNCTION = %s\n","Michalewitz");
			break;
		case 17:
			printf("FUNCTION = %s\n","Multimod");
			break;
		case 18:
			printf("FUNCTION = %s\n","Powell");
			break;
		case 19:
			printf("FUNCTION = %s\n","Rana");
			break;
		case 20:
			printf("FUNCTION = %s\n","Shubert");
			break;
		case 21:
			printf("FUNCTION = %s\n","3D-AB");
			printf("PROTEIN SIZE: %d\n", strlen(SEQUENCE));
			break;
		case 22:
			printf("FUNCTION = %s\n","10-bar-truss");
			break;
		case 23:
			printf("FUNCTION = %s\n","2D-AB");
			printf("PROTEIN SIZE: %d\n", strlen(SEQUENCE));
			break;
		case 24:
			printf("FUNCTION = %s\n","Schaffer f6");
			break;
		case 25:
			printf("FUNCTION = %s\n","Shifted Sphere");
			break;
		case 26:
			printf("FUNCTION = %s\n","Shifted Schwefel 2.21");
			break;
		case 27:
			printf("FUNCTION = %s\n","Shifted Rosenbrock");
			break;
		case 28:
			printf("FUNCTION = %s\n","Shifted Rastrigin");
			break;
		case 29:
			printf("FUNCTION = %s\n","Shifted Griewank");
			break;
		case 30:
			printf("FUNCTION = %s\n","Shifted Ackley");
			break;
		default:
		        printf("Info: Invalid function.\n") ;
		        exit(0);
	}
	printf("****************\n");
}



int ffscanf(char *fieldname, FILE *fp, char *format, void *inbuffer)
{
	char buffer[100];
	int len;
    int commentflag = 0;
	char *pch;
	char *pch2;
	
	do
	{
		if(fgets(buffer, 99, fp) == 0) return FAIL;
		buffer[99] = '\0';
		len = strlen(buffer);
	    if (buffer[len - 1] == '\n') buffer[len - 1] = '\0';

	    switch (commentflag) 
		{
		    case 0:
				if (strstr(buffer, "/*") != 0) 
				{
			    	commentflag = 1;
			    	if (strstr(buffer, "*/") != 0)
						commentflag = 2;
				}
				break;

		    case 1:
				if (strstr(buffer, "*/") != 0)
				    commentflag = 2;
				break;

			    case 2:
				if (strstr(buffer, "/*") != 0) 
				{
				    commentflag = 1;
				    if (strstr(buffer, "*/") != 0)
					commentflag = 2;
				}
				else
				    commentflag = 0;
				break;
			    }	
	}while(commentflag != 0);	

	//separate field name: token = "="
	if (strstr (buffer, fieldname) != 0)
	{
		pch = strtok (buffer,"=");

		while (pch != NULL)
		{
			pch2 = pch;
    		pch = strtok (NULL, "= ");
		}
		sscanf(pch2, format, inbuffer);
		return 1;//ok
	}
	else return 0; 
}


//==============
#include "ProblemFunc.c"
#include "./algorithms/ABC.c"
#include "./algorithms/DE.c"
#include "./algorithms/PSO.c"
#include "./algorithms/DE_BBO.c"
#include "./algorithms/PRVNS.c"
