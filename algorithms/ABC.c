//Functions ABC
/*Main program of the ABC algorithm*/
int ABC(int sp_index)
{
	int iter,run,j,k,i,t,
	    file_iter_count; //counter for the global iteration number to put into file
	double mean, maxfit, result;

	double prob[cur_POP_SIZE[0]]; /*prob is a vector holding probabilities of food sources (solutions) to be chosen*/
	double *solution; //solution[DIM]; /*New solution (neighbour) produced by v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) j is a randomly chosen parameter and k is a randomlu chosen solution different from i*/
	double ObjValSol; /*Objective function value of new solution*/
	double FitnessSol; /*Fitness value of new solution*/
	int neighbour, param2change; /*param2change corrresponds to j, neighbour corresponds to k in equation v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij})*/
	double GlobalMins[1]; /*GlobalMins holds the GlobalMin of each run in multiple runs*/
	double r; /*a random number in the range [0,1)*/
	double sum; //for centroid calculation

	int local_func_eval;


	FILE *fp;

	solution = (double*)malloc(DIM * sizeof(double));

	mean=0;

	for(run=0;run<1;run++)
	{
		//printf("SUC DENTRO: %d \n", suc[0]);
		if (suc[0] == 0)
		{
			for (i = 0;i< cur_POP_SIZE[0]; i++)
				Eco[sp_index].trial[i] = 0;
		}

		//each evolutionary period runs EVO_STEP functions evaluations
		local_func_eval = 0;
		while (local_func_eval < EVO_STEP)
			// for (iter=0;iter<EVO_STEP;iter++)
		{
			//    CalculateProbabilities(fitness);
			/* A food source is chosen with the probability which is proportioal to its quality*/
			/*Different schemes can be used to calculate the probability values*/
			/*For example prob(i)=fitness(i)/sum(fitness)*/
			/*or in a way used in the metot below prob(i)=a*fitness(i)/max(fitness)+b*/
			/*probability values are calculated by using fitness values and normalized by dividing maximum fitness value*/
			maxfit=Eco[sp_index].fit[0];
			for (i=1;i<cur_POP_SIZE[0];i++)
			{
				if (Eco[sp_index].fit[i]>maxfit)
					maxfit=Eco[sp_index].fit[i];
			}
			for (i=0;i<cur_POP_SIZE[0];i++)
			{
				prob[i]=(0.9*(Eco[sp_index].fit[i]/(double)maxfit))+0.1;
			}

			/*Employed Bee Phase*/
			for (i=0;i<cur_POP_SIZE[0];i++)
			{
				/*The parameter to be changed is determined randomly*/
				r = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );
				param2change=(int)(r*DIM);

				/*A randomly chosen solution is used in producing a mutant solution of the solution i*/
				r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
				neighbour=(int)(r*cur_POP_SIZE[0]);

				/*Randomly selected solution must be different from the solution i*/        
				while(neighbour==i)
				{
					r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
					neighbour=(int)(r*cur_POP_SIZE[0]);
				}

				for(j=0;j<DIM;j++)
					solution[j]=Eco[sp_index].pop[i][j];

				/*v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
				r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
				solution[param2change]=Eco[sp_index].pop[i][param2change]+(Eco[sp_index].pop[i][param2change]-Eco[sp_index].pop[neighbour][param2change])*(r-0.5)*2;

				/*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
				if (solution[param2change]<lb)
					solution[param2change]=lb;
				if (solution[param2change]>ub)
					solution[param2change]=ub;
				ObjValSol=function(solution, sp_index);
				local_func_eval++;

				//fitness
				result=0;
				if(ObjValSol>=0)
				{
					FitnessSol=1/(ObjValSol+1);
				}
				else
				{
					FitnessSol=1+fabs(ObjValSol);
				}

				/*a greedy selection is applied between the current solution i and its mutant*/
				if (FitnessSol>Eco[sp_index].fit[i])
				{
					/*If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i*/
					Eco[sp_index].trial[i]=0;
					for(j=0;j<DIM;j++)
						Eco[sp_index].pop[i][j]=solution[j];
					Eco[sp_index].fo[i]=ObjValSol;
					Eco[sp_index].fit[i]=FitnessSol;
				}
				else
				{   /*if the solution i can not be improved, increase its trial counter*/
					Eco[sp_index].trial[i]=Eco[sp_index].trial[i]+1;
				}

				//verity stop condition
				if ( local_func_eval > EVO_STEP-1 ) goto stop;

			}//END FOR
			/*end of employed bee phase*/

			//    CalculateProbabilities(fitness);
			/* A food source is chosen with the probability which is proportioal to its quality*/
			/*Different schemes can be used to calculate the probability values*/
			/*For example prob(i)=fitness(i)/sum(fitness)*/
			/*or in a way used in the metot below prob(i)=a*fitness(i)/max(fitness)+b*/
			/*probability values are calculated by using fitness values and normalized by dividing maximum fitness value*/
			maxfit=Eco[sp_index].fit[0];
			for (i=1;i<cur_POP_SIZE[0];i++)
			{
				if (Eco[sp_index].fit[i]>maxfit)
					maxfit=Eco[sp_index].fit[i];
			}
			for (i=0;i<cur_POP_SIZE[0];i++)
			{
				prob[i]=(0.9*(Eco[sp_index].fit[i]/(double)maxfit))+0.1;
			}

			//    SendOnlookerBees(Foods,f,fitness);
			i=0;
			t=0;

			/*onlooker Bee Phase*/
			while(t<cur_POP_SIZE[0])
			{
				r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
				if(r<prob[i]) /*choose a food source depending on its probability to be chosen*/
				{        
					t++;

					/*The parameter to be changed is determined randomly*/
					r = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );
					param2change=(int)(r*DIM);

					/*A randomly chosen solution is used in producing a mutant solution of the solution i*/
					r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
					neighbour=(int)(r*cur_POP_SIZE[0]);

					/*Randomly selected solution must be different from the solution i*/        
					while(neighbour==i)
					{
						r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
						neighbour=(int)(r*cur_POP_SIZE[0]);
					}
					for(j=0;j<DIM;j++)
						solution[j]=Eco[sp_index].pop[i][j];

					/*v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
					r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
					solution[param2change]=Eco[sp_index].pop[i][param2change]+(Eco[sp_index].pop[i][param2change]-Eco[sp_index].pop[neighbour][param2change])*(r-0.5)*2;

					/*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
					if (solution[param2change]<lb)
						solution[param2change]=lb;
					if (solution[param2change]>ub)
						solution[param2change]=ub;
					ObjValSol=function(solution, sp_index);
					local_func_eval++;
					//fitness
					result=0;
					if(ObjValSol>=0)
					{
						FitnessSol=1/(ObjValSol+1);
					}
					else
					{
						FitnessSol=1+fabs(ObjValSol);
					}

					/*a greedy selection is applied between the current solution i and its mutant*/
					if (FitnessSol>Eco[sp_index].fit[i])
					{
						/*If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i*/
						Eco[sp_index].trial[i]=0;
						for(j=0;j<DIM;j++)
							Eco[sp_index].pop[i][j]=solution[j];
						Eco[sp_index].fo[i]=ObjValSol;
						Eco[sp_index].fit[i]=FitnessSol;
					}
					else
					{   /*if the solution i can not be improved, increase its trial counter*/
						Eco[sp_index].trial[i]=Eco[sp_index].trial[i]+1;
					}
				} /*if */
				i++;
				if (i==cur_POP_SIZE[0]-1)
					i=0;

				//verity stop condition
				if ( local_func_eval > EVO_STEP-1 ) goto stop;

			}/*END while*/
			/*end of onlooker bee phase     */

stop:

			//    MemorizeBestSource(Foods,f,GlobalParams,GlobalMin,best_index);
			/*The best food source is memorized*/
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

			//    SendScoutBees();
			int maxtrialindex;
			maxtrialindex=0;
			for (i=1;i<cur_POP_SIZE[0];i++)
			{
				if (Eco[sp_index].trial[i]>Eco[sp_index].trial[maxtrialindex])
					maxtrialindex=i;
			}
			if(Eco[sp_index].trial[maxtrialindex]>=limit)
			{

				for (j=0;j<DIM;j++)
				{
					Eco[sp_index].pop[maxtrialindex][j]= randon(lb,ub);
					solution[j]=Eco[sp_index].pop[maxtrialindex][j];
				}
				Eco[sp_index].fo[maxtrialindex]=function(solution, sp_index);
				local_func_eval++;

				//fitness
				result=0;
				if(Eco[sp_index].fo[maxtrialindex]>=0)
				{
					Eco[sp_index].fit[maxtrialindex]=1/(double)(Eco[sp_index].fo[maxtrialindex]+1.0);
				}
				else
				{
					Eco[sp_index].fit[maxtrialindex]=1+fabs(Eco[sp_index].fo[maxtrialindex]);
				}
				Eco[sp_index].trial[maxtrialindex]=0;
			}//end if
			/*
			   for (j=0;j<cur_POP_SIZE;j++)
			   {
			   printf("Variaveis: ");			
			   for (k=0; k<DIM;k++) //variables
			   {
			   printf("%.2f ",Foods[j][k]);
			   }
			   printf("\n Fo e Fit:");
			   printf("%.2f %.2f \n",f[j],fitness[j]);
			   }
			   printf("\n Melhores Param e FO min: ");
			   for (k=0; k<DIM;k++) //variables
			   {
			   printf("%.2f ",GlobalParams[k]);
			   }
			   printf("%.2f \n",GlobalMin[0]);


			   for(j=0;j<DIM;j++)
			   {
			   printf("GlobalParam[%d]: %f\n",j+1,GlobalParams[j]);
			   }
			   printf("%d. run: %e \n",run+1,GlobalMin[0]);
			   */
			GlobalMins[run]=Eco[sp_index].bestfo[0];
			mean=mean+Eco[sp_index].bestfo[0];

			// if (sp_index == 1) //record the best solution of sp 1 of each loop of ABC
			// {
			//	fp = fopen( "./report/sp1.txt", "a" ); //append to the file.
			//	fprintf(fp,"%.4f\n", Eco[sp_index].bestfo[0]);
			//	fclose(fp);
			// }//end if

		}// END for iter while
		//printf("\nLocal Iter: %d\n",local_func_eval);
	}//END RUN
	//mean=mean/runtime;
	//printf("Means of %d runs: %e\n",runtime,mean);

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
	//			printf("\nCi from sp %d: %.2f\n",sp_index,Eco[sp_index].CiFO);


	free(solution);

} //end MAIN

