typedef struct _VNS_SOLUTION
{
	double bestfo;         //best fo value
	int c_aval;		//avaliation count total
	int c_aval_best;	//alaviation count when found best			
	time_t stime;//start time
	time_t etime;//end time
	double t_total;//total time


}pVNS_SOLUTION;

typedef struct _VNS
{
	int P; /*The metric used, l1,l2,l_inf */
	int FUNCTION;
	//Problem definitions
	int DIM;    //number of problem variables
	int RUN;  /*Algorithm can run many times in order to see its robustness*/
	double *best; //[DIM];           //best solution found
	double bestfo;         //best fo value
	double LB; //lower bound of the variables
	double UB; //upper bound of the variables
	int KMAX,TMAX,AVAL_MAX;
	int RADII_T_FLAG;
	int know;
	double Q,RADII,RHO,EPSILON;
	double *r;
	int METHOD;
	pVNS_SOLUTION solv;	
	int ECO_STEP;
	int EVO_STEP;
	int VNS_POP;
	int POP_INDEX;
	double DELTA_MUTATION;
	int G_MAX;
	double PC; //probabilidade de crossover
}pVNS;

void prepare_prvns(){

	PC = 0.9; //default value from DE
	KMAX = 5; //default

	//definir raio
	RAIO = malloc(sizeof(double) * KMAX);
	RAIO[0] = 0.1f;
	RAIO[1] = 0.3f;
	RAIO[2] = 0.5f;
	RAIO[3] = 0.7f;
	RAIO[4] = 0.9f;

	k_index = malloc(sizeof(int) * 50);
	int i;
	for(i = 0; i<50; i++){
		k_index[i] = 0;
	}
}

//Populational Reduced VNS
void PRVNS(int sp_index){

	int i,j, t=0,iter=0;
	int r1,r2,p, best_aval, best_index;
	int a;
	double *y = malloc(sizeof(double) * DIM);

	double fy;
	double bestfo;
	/* 	while((iter + 1 < G_MAX) && ((t + cur_POP_SIZE[0]) <= EVO_STEP)){ */
	while((t + cur_POP_SIZE[0]) <= EVO_STEP){

		//etapa populacional
		for(i =0; i < cur_POP_SIZE[0];i++){

			//usar os whiles separados para tentar melhorar o tempo
			//seleciona um por vez pode ser melhor que selecionar todos
			//e testar no mesmo while
			do{
				r1 = randon(0,cur_POP_SIZE[0]);
			}while(r1 == i);

			do{
				r2 = randon(0,cur_POP_SIZE[0]);
			}while( r1 == r2 || r2 == i );


			//pos na dim para sempre perturbar
			p = randon(0,DIM);

			//etapa da evolucao diferencial
			for(j=0;j<DIM;j++){

				if( (j == p) || (randon(0,1) < PC)){

					y[j] = Eco[sp_index].pop[r2][j] + ( randon(-1.0f * RAIO[k_index[i]], RAIO[k_index[i]]) * Eco[sp_index].pop[r1][j] ); 
					/* y[j] = x[i]->best[j] + ( randon(-1.0f * vns->r[x[i]->know], vns->r[x[i]->know]) * x[r1]->best[j] );  */
				}else{
					//nao perturbar
					y[j] = Eco[sp_index].pop[i][j];
				}	


				//trata se esta no intervalo permitido da FO
				if(ub < y[j]){
					y[j] = ub;
				}else if(lb > y[j]){
					y[j] = lb;
				}
			}

			/* printf("PC = %lf x = [ ", PC); */
			/* for(a = 0; a < DIM; a++) */
			/* { */
			/* 	printf("%lf, ",Eco[sp_index].pop[i][a]); */
			/* } */
			/* printf("]\n"); */

			/* printf("y = [ "); */
			/* for(a = 0; a < DIM; a++) */
			/* { */
			/* 	printf("%lf, ",y[a]); */
			/* } */
			/* printf("]\n"); */
			//etapa de avaliacao e troca de vizinhanca
			fy = function(y,sp_index);
			t++;

			if(fy < Eco[sp_index].fo[i]){

				/* memcpy(Eco[sp_index].pop[i],y,sizeof(double) * DIM); */
				/* memcpy(Eco[sp_index].best[i],y,sizeof(double) * DIM); */
				Eco[sp_index].fo[i] = fy;
				for(a = 0; a < DIM; a++)
				{
					Eco[sp_index].pop[i][a] = y[a];
				}

				k_index[i] = 0;

			}else{
				if(k_index[i] <= KMAX -1){
					k_index[i] += 1;
				}else{
					k_index[i] = 0; //reinicia k
				}
			}
		}
		iter++;
	}
	//==============================================================
	// Localizar o melhor
	//==============================================================
	int best = 0;
	for(a = 1; a < cur_POP_SIZE[0]; a++)
	{
		if(Eco[sp_index].fo[a] < Eco[sp_index].fo[best])
		{
			best = a;
		}
	}
	//==============================================================
	// Atualizar o Best
	//==============================================================
	for(a = 0; a < DIM; a++)
	{
		Eco[sp_index].best[a] = Eco[sp_index].pop[best][a];
	}
	Eco[sp_index].bestfo[0] = Eco[sp_index].fo[best];
	Eco[sp_index].best_index[0] = best;
	//==============================================================
	// Calcular a Centroide
	//==============================================================
	int b;
	double obj;
	for (a = 0; a < DIM; a++)//each variable
	{
		obj = 0.0;
		for (b = 0; b < cur_POP_SIZE[0]; b++) //each solution
		{
			obj += Eco[sp_index].pop[b][a];
		}
		Eco[sp_index].Ci[a] = obj / (double)cur_POP_SIZE[0];
	}

	/* free(f_name); */
	/* free(vns_plot); */

	free(y);
	/* free(vns); */
	/* for(i=0;i<POP_SIZE;i++){ */
	/* 	free(x[i]->best); */
	/* 	free(x[i]); */
	/* } */
	/* free(x); */
	/* free(x); */
	/* free(y); */
	/* free(r); */

	//when return, return the best of the best
}



