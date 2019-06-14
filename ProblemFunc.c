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
static double pi;
static double e;

double objfunc(double sol[], int sp_index)
{
    int j, i;
    double top = 0.00 , top1 = 0.00, top2 = 0.00;
    double aux = 0.0;
    double aux1 = 0.0;

    double z = 0.0;
    double zx[DIM];
    double Fx = 0.0;

	pi = acos(-1.0);
	e = exp(1.0);

	//for Multimod
	double t, s; //p
	//for Rana
    	double t1, t2, sum;

    double *y; //for Penalized function #1
    
    //for 10-bar-truss
	double *area;
	double stress[10] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
	double displX[6] = {0., 0., 0., 0., 0., 0.};
	double displY[6] = {0., 0., 0., 0., 0., 0.};
	double somF = 0.;
	double somg1 = 0.;
	double somg2 = 0.;

double ro = 0.1;
double P  = 10000.;

double E  = 10000.;

int num_node = 0, 
    num_elem = 0;
int kode[6] = {0, 0, 0, 0, 0, 0};
double coordx[6] = {0., 0., 0., 0., 0., 0.},
       coordy[6] = {0., 0., 0., 0., 0., 0.},
       fX[6] = {0., 0., 0., 0., 0., 0.},
       fY[6] = {0., 0., 0., 0., 0., 0.};
int lm[4] = {0,0,0,0}, 
    neq = 0;
int member[10][2];
double ymod[10] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
double a[12][12], 
       b[12] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, 
       stiffmatrix[4][4];

int n = 0, nn = 0;

double dx = 0.,dy = 0.,xl = 0.;
double du = 0.,dv = 0.,p = 0.,dl = 0.;

int nl = 0;
int n1 = 0,m1 = 0;
int neq1 = 0;

int m = 0, k = 0, g = 0;
double cosa = 0.,sina = 0.,comm = 0.;

double length[10] = {
	360.,
	360.,
	360.,
	360.,
	360.,
	360.,
	509.116882454,
	509.116882454,
	509.116882454,
	509.116882454
};

    //for 3D-AB
    double w1, w2, LJ, LJ_aux;
    double r, rx, ry, rz, u1,u2,u3,v3;
    double nx,ny,nz,vx,vy,vz,tx,ty,tz,ux,uy,uz; 
    int dim_aux = ProteinSize - 2;
    amino *amino_pos; //
    //for 2D-AB
    double a_ab,b_ab,c_ab,d_ab;
    double v1,v2;

    switch (FUNCTION) {
    case 0: //Rastrigin

        for(j=0;j<DIM;j++)
        {
            top=top+(pow(sol[j],(double)2)-10*cos(2*M_PI*sol[j])+10);
        }
	func_eval[sp_index] += 1;
        return top;

    case 1: //Schaffer F7

        top1 = 0;
        for(j=0;j<DIM;j++)
        {
        top=top+(pow(sol[j],(double)2));
        }
        top = pow(top,(double)0.25);
        for(j=0;j<DIM;j++)
        {
        top1=top1+(pow(sol[j],(double)2));
        }
        top1=pow(top1,(double)0.1);
        top1 = pow(sin(50*top1),(double)2) +1.0;

	func_eval[sp_index] += 1;

        return top*top1;

    case 2: //Griewank

        top=0;
        top1=0;
        top2=1;
        for(j=0;j<DIM;j++)
        {
        top1=top1+pow((sol[j]),(double)2);
        top2=top2*cos((((sol[j])/sqrt((double)(j+1)))*M_PI)/180);
        }
        top=(1/(double)4000)*top1-top2+1;

	func_eval[sp_index] += 1;

        return top;

    case 3: //Ackley
/*
-	Dimension: n arbitrary
-       Domain:   -32 <= | x_i | <= 32.0
-       Minimum 0 at x_i = 0.0
*/
        for (i = 0; i < DIM; i++)
        {
        aux += sol[i]*sol[i];
        }
        for (i = 0; i < DIM; i++)
        {
        aux1 += cos(2.0*M_PI*sol[i]);
        }

	func_eval[sp_index] += 1;

        return (-20.0*(exp(-0.2*sqrt(1.0/(float)DIM*aux)))-exp(1.0/(float)DIM*aux1)+20.0+exp(1));

    case 4: //Rosenbrock

        for (i = 0; i < DIM-1; i++)
        {
            top=top+100.*pow((sol[i+1] - pow(sol[i],2.)),2) + pow((1. - sol[i]),2);
        }

	func_eval[sp_index] += 1;

        return top;

    case 5: //Sphere
	for(j=0;j<DIM;j++)
	{
		top=top+sol[j]*sol[j];
	}

        func_eval[sp_index] += 1;

	return top;

    case 6: //StretchedV
/*
- Domain:  | x_i | <= 10.0
Global minimum is 0.0 at x_i = 0.00
*/
	  sum = 0.0;
	  for (i = 0; i < DIM-1; i++) {
	    aux = sol[i+1]*sol[i+1] + sol[i]*sol[i];
	    sum += pow(aux, 0.25) * (pow(sin(50.0 * pow(aux, 0.1)), 2.0)+1.0);
	  }
          func_eval[sp_index] += 1;
	  return sum;

    case 7: // Schwefel's function 2.22
/*
- Domain:  | x_i | <= 10.0
Global minimum is 0.0 at x_i = 0.00
*/
	for (i=0;i<DIM;i++)
	{
		 aux += fabs(sol[i]);
		 aux1 *= fabs(sol[i]);
	}

	func_eval[sp_index] += 1;
	return (aux+aux1);

    case 8: // Step function
/*
-   Domain: | x_i  | < 100.0
-   Global minimum is 0.0 at x_i = 0.5
*/
	for (i=0;i<DIM;i++)
	{
		aux = (sol[i]+0.5);
		aux1 += aux*aux; 
	}
	func_eval[sp_index] += 1;
	return (aux1);

    case 9: //Generalized Schwefel's function 2.26
        //known_optimal = -418.982887272433 at  sol(i)=420.9687
	for (i=0;i<DIM;i++)
	{
		aux += sol[i]*sin(sqrt(fabs(sol[i]))); 
	}
	func_eval[sp_index] += 1;
	return(-1*aux/DIM);

    case 10: //Generalized Penalized function #1
	// -500 <= xi <= 500
	//known_optimal = 1.57044103551786e-032;
	y = (double*) malloc (DIM * sizeof(double));
	for (i=0;i<DIM;i++)
	{
		y[i]=0.0;
	}

	for (i=0;i<DIM;i++)
	{
		y[i]=1+(sol[i]+1)/4.0;
	}

	for (i=0;i<DIM-1;i++)
	{
		aux += pow(y[i]-1,2.0)*(1.0+10.0*pow(sin(PI*y[i+1]),2.0)); 
	}
	for (i=0;i<DIM;i++)
	{
		aux1 += TempValue(sol[i],10,100,4);
	}
	func_eval[sp_index] += 1;

	aux = (10.0*pow(sin(PI*y[0]),2.0)+aux+pow(y[DIM-1]-1,2))*PI/DIM+aux1;

        free(y);

	return ( aux );
	
     case 11: //Generalized Penalized function #2
	//known_optimal = 1.34969464963992e-032;
	for (i=0;i<DIM-1;i++)
	{
		aux += pow(sol[i]-1,2.0)*(1.0+10.0*pow(sin(3*PI*sol[i+1]),2.0)); 
	}
	for (i=0;i<DIM;i++)
	{
		aux1 += TempValue(sol[i],5,100,4);
	}
	 
	func_eval[sp_index] += 1;
	return ( (pow(sin(3.0*PI*sol[0]),2.0)+aux+pow(sol[DIM-1]-1,2.0)
		*(1.0+pow(sin(2.0*PI*sol[DIM-1]),2.0)))/10.0+aux1 );

    case 12: //Levy Function
	//x[i] = 1 f(x[i])=0
	y = (double*) malloc (DIM * sizeof(double));
	for (i = 0; i< DIM; i++)
		y[i] = 1+(sol[i]-1)/4.0;
	aux = pow(sin(PI*y[0]),2.0);
	for (i = 0; i<DIM-1;i++)
	    aux = aux + pow(y[i]-1,2.0)*(1+10*pow(sin(PI*y[i]+1),2.0));

	aux = aux+pow(y[DIM-1]-1,2.0)*(1+pow(sin(2*PI*y[DIM-1]),2.0) );

	func_eval[sp_index] += 1;
	free (y);
	return ( aux );

    case 13: // Zakharov function  //x[i] = 0 f(x[i])=0
	aux = aux1 = 0.;
	for (j = 0; j< DIM; j++)
	{
	    aux = aux + pow(sol[j],2.0);
	    aux1 = aux1+0.5*j*sol[j];
	}

	func_eval[sp_index] += 1;

	return ( aux+pow(aux1,2.0)+pow(aux1,4.0) );

    case 14: //Egg holder
/*
- Dimension: n
- Domain:  -512 < | x_i | < 512
- Minimum for n=2 fmin(512, 404.2319) = -959.641
*/
	  aux = 0.0;
	  for (i = 0; i < DIM-1; i++) 
	  {
	    aux += -(sol[i+1] + 47.0) * sin(sqrt(fabs(sol[i+1] + sol[i] * 0.5 + 47.0))) + sin(sqrt(fabs(sol[i] - (sol[i+1] + 47.0)))) * (-sol[i]);
	  }
	  func_eval[sp_index] += 1;
	  return (aux);

    case 15: //Generalized Holzman
/*
Dimension: n
Domain: | x[i] | <= 10
Global minimum: 0 at x[i] = 0
*/
  	for (i = 0; i < DIM; i++) 
  	{
  	  aux += i * pow(sol[i] , 4);
  	}
        func_eval[sp_index] += 1;
  	return aux;

    case 16://Michalewitz
/*
Dimension: n
Domain: 0< | x[i] | <= PI
Global minimum: x[] = -0.966*n
*/
	aux=0;
	for (i=0;i<DIM;i++) {
		aux = aux + sin(sol[i])
			* pow(sin((i+1)*sol[i]*sol[i]/(float)PI), 2.0 * 10.0);
	}
        func_eval[sp_index] += 1;
	return(-1*aux);

    case 17: //Multimod
/*
Dimension: n
Domain: -10<= x[i] <= 10
Global minimum: x[0] = 0
*/
  	s = p = fabs(sol[0]);
	for (i = 1; i < DIM; i++) {
     		t = fabs(sol[i]);
     		s += t;
     		p *= t;
  	}
        func_eval[sp_index] += 1;
  	return s + p;

    case 18: //Powell
/*
Dimension: n > 4
Domain: -4<= x[i] <= 5
Global minimum: at (3, -1, 0, 1, ..., 3, -1, 0, 1) with fmin = 0
*/
	  aux = 0.0;
	  for (j = 1; j <= (int)DIM/4; j++) {
		aux +=  pow(sol[4*j-4] + 10 * sol[4*j-3],2.0)
        	  + 5 * pow(sol[4*j-2] - sol[4*j-1],2.0)
        	  + pow(sol[4*j-3] - 2 * sol[4*j-2], 4.0)
		  + 10 * pow(sol[4*j - 4] - sol[4*j-1], 4.0);
	  }
          func_eval[sp_index] += 1;
	  return aux;

    case 19: //Rana
/*
Dimension: n
Domain: -520<= x[i] <= 520
Global minimum: ???
*/
	  sum = 0.0;
	  for (i = 0; i < DIM-1; i++) {
	    t1 = sqrt(fabsf(sol[i+1] + sol[i] + 1.0));
	    t2 = sqrt(fabsf(sol[i+1] - sol[i] + 1.0));
	    sum += (sol[i+1] + 1.0) * cos(t2) * sin(t1) + cos(t1) * sin(t2) * sol[i];
	  }
          func_eval[sp_index] += 1;
	  return sum/(double)(DIM-1);

    case 20: //Shubert
/*
-   Domain  |x| <= 10.0
-   Number of local minimum = 400
-   Global minimum fmin = -24.062499 at the ff. points
-    (-6.774576, -6.774576), ..., (5.791794, 5.791794)
*/
	  sum = 0.0;
	  for (i = 0; i < DIM; i++) {
		sum += -sin(2.0*sol[i]+1.0)
	          -2.0*sin(3.0*sol[i]+2.0)
        	  -3.0*sin(4.0*sol[i]+3.0)
        	  -4.0*sin(5.0*sol[i]+4.0)
        	  -5.0*sin(6.0*sol[i]+5.0);
	  }


          func_eval[sp_index] += 1;
   	  return sum/(DIM/2.0);

    case 21: 
//3D-AB
    amino_pos = (amino *) malloc ((ProteinSize) * sizeof (amino));
	//polar --> cartesian
    amino_pos[0].x = amino_pos[0].y = amino_pos[0].z = 0;
    amino_pos[1].x = amino_pos[1].z = 0;  amino_pos[1].y = 1;

	nx=0;
	ny=0; 
	nz=1;

    for (i = 1; i < ProteinSize-1; i++)
    {  

		if(sol[i-1]>PI)
			sol[i-1] = randon(0.0,PI);
		if(sol[i-1]<-PI)
			sol[i-1] = randon(-PI,0.0);
		if(sol[i+dim_aux-1]>PI)
			sol[i+dim_aux-1] = randon(0.0,PI);
		if(sol[i+dim_aux-1]<-PI)
			sol[i+dim_aux-1] = randon(-PI,0.0);	

		tx=0;ty=0;tz=0;
		vx = amino_pos[i].x - amino_pos[i-1].x;
   		vy = amino_pos[i].y - amino_pos[i-1].y;
		vz = amino_pos[i].z - amino_pos[i-1].z;

     Rotation(amino_pos[i].x+vx, amino_pos[i].y+vy, amino_pos[i].z+vz, amino_pos[i].x, amino_pos[i].y, amino_pos[i].z,nx,ny,nz,sol[i-1],&tx,&ty,&tz);

	 if (i == 1)
	 	 Rotation(tx,ty,tz,amino_pos[i].x,amino_pos[i].y,amino_pos[i].z,vx,vy,vz,0,&amino_pos[i+1].x,&amino_pos[i+1].y,&amino_pos[i+1].z);		
	 else
	 	Rotation(tx,ty,tz,amino_pos[i].x,amino_pos[i].y,amino_pos[i].z,vx,vy,vz,sol[(i-2)+dim_aux],&amino_pos[i+1].x,&amino_pos[i+1].y,&amino_pos[i+1].z);


		ux = amino_pos[i+1].x - amino_pos[i].x;
	   	uy = amino_pos[i+1].y - amino_pos[i].y;
   		uz = amino_pos[i+1].z - amino_pos[i].z;

  		nx = vy*uz-vz*uy;
  		ny = vz*ux-vx*uz;
  		nz = vx*uy-vy*ux;

    }

	//bond angles
    w1 = 0;
    for (i = 1; i < ProteinSize-1; i++) 
	{
     	u1 = amino_pos[i].x - amino_pos[i-1].x;
     	u2 = amino_pos[i].y - amino_pos[i-1].y;
     	u3 = amino_pos[i].z - amino_pos[i-1].z;
     	v1 = amino_pos[i+1].x - amino_pos[i].x;
    	v2 = amino_pos[i+1].y - amino_pos[i].y;
     	v3 = amino_pos[i+1].z - amino_pos[i].z;
     	w1 += (u1*v1 + u2*v2 + u3*v3)/(sqrt(pow(u1,2) + pow(u2,2) + pow(u3,2)) * sqrt(pow(v1,2) + pow(v2,2) + pow(v3,2)));
	}


	//torsion angles
    w2 = 0;
    for (i = 1; i < ProteinSize-2; i++)

    {   
     	u1 = amino_pos[i].x - amino_pos[i-1].x;
    	u2 = amino_pos[i].y - amino_pos[i-1].y;
     	u3 = amino_pos[i].z - amino_pos[i-1].z;
     	v1 = amino_pos[i+2].x - amino_pos[i+1].x;
     	v2 = amino_pos[i+2].y - amino_pos[i+1].y;
     	v3 = amino_pos[i+2].z - amino_pos[i+1].z;
     	w2 = w2 - ((u1*v1 + u2*v2 + u3*v3)/(sqrt(pow(u1,2) + pow(u2,2) + pow(u3,2))*sqrt(pow(v1,2) + pow(v2,2)+pow(v3,2))))/2.; 
    }

	//Lennard-Jones
	LJ = 0;
    for (i = 0; i < ProteinSize-2; i++)
    {   for (j = i + 2; j < ProteinSize; j++)
        {   
			rx = amino_pos[i].x - amino_pos[j].x;
			ry = amino_pos[i].y - amino_pos[j].y;
			rz = amino_pos[i].z - amino_pos[j].z;
			r = sqrt(rx * rx + ry * ry + rz * rz);

			 LJ_aux = 4 * (1/pow(r, 12) - 1/pow(r, 6));

			if ((Aminoacido[i].Tipo != Aminoacido[j].Tipo) || (Aminoacido[i].Tipo == -1 && Aminoacido[j].Tipo == -1))
			{	
				 LJ_aux *= 0.5;
			}
			
			LJ +=  LJ_aux;
        }
    }
	    func_eval[sp_index] += 1;
	    free(amino_pos);
	    return(w1 + w2 + LJ);
/*end 3D-AB function*/
    case 22: //10-bar-truss


/*	truss(area, stress, displX, displY);         call truss FE to get lengths,
						 stresses and displacements

init */
	for (i = 0; i < 12; i++)
		for (j = 0; j < 12; j++)
			a[i][j] = 0.;

	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
		       stiffmatrix[i][j] = 0.;

	area = (double*) malloc (DIM * sizeof(double));
	for (i=0; i < DIM; i++)
	{
		area[i] = sol[i];
	}

	num_node = 6;
	num_elem = 10;


	kode[4] = 3; /* fixo X e Y*/
	coordx[4] = 0;
	coordy[4] = 0;
	fX[4] = 0;
	fY[4] = 0;
	kode[2] = 0;
	coordx[2] = 360;
	coordy[2] = 0;
	fX[2] = 0;
	fY[2] = 50;
	kode[0] = 0;
	coordx[0] = 720;
	coordy[0] = 0;
	fX[0] = 0;
	fY[0] = 50;
	kode[5] = 3;
	coordx[5] = 0;
	coordy[5] = -360;
	fX[5] = 0;
	fY[5] = 0;
	kode[3] = 0;
	coordx[3] = 360;
	coordy[3] = -360;
	fX[3] = 0;
	fY[3] = -150;
	kode[1] = 0;
	coordx[1] = 720;
	coordy[1] = -360;
	fX[1] = 0;
	fY[1] = -150;


	member[0][0] = 5; member[0][1] = 3;
	member[1][0] = 3; member[1][1] = 1;
	member[2][0] = 6; member[2][1] = 4;
	member[3][0] = 4; member[3][1] = 2;
	member[4][0] = 3; member[4][1] = 4;
	member[5][0] = 1; member[5][1] = 2;
	member[6][0] = 5; member[6][1] = 4;
	member[7][0] = 6; member[7][1] = 3;
	member[8][0] = 3; member[8][1] = 2;
	member[9][0] = 1; member[9][1] = 4;

	ymod[0] = E;
	ymod[1] = E;
	ymod[2] = E;
	ymod[3] = E;
	ymod[4] = E;
	ymod[5] = E;
	ymod[6] = E;
	ymod[7] = E;
	ymod[8] = E;
	ymod[9] = E;
/*end init

//stiff */
	/* Initialize Stiffness Matrix and Load Vector */
	neq = 2 * num_node;
	for (m=0; m<neq;m++)
	{
		b[m] = 0.;
		for(k=0;k<neq;k++) a[m][k]=0.;
	}
	/*  >>>>>>> loop over all elements  <<<<<<< */
	for(m=0;m<num_elem;m++)
	{
		/* Determination of Length, Direction */
		i=member[m][0] - 1;
		j=member[m][1] - 1;
		dx = coordx[j]-coordx[i];
		dy = coordy[j]-coordy[i];
		xl = sqrt((dx * dx) + (dy * dy));
		cosa=dx/xl;
		sina=dy/xl;
		comm=area[m]*ymod[m]/xl;
		/* Construction of the 4X4 Element Stiffness Matrix */
		stiffmatrix[0][0] = cosa * cosa * comm;
		stiffmatrix[0][1] = cosa * sina * comm;
		stiffmatrix[0][2] = -stiffmatrix[0][0];
		stiffmatrix[0][3] = -stiffmatrix[0][1];
		stiffmatrix[1][0] =  stiffmatrix[0][1];
		stiffmatrix[1][1] = sina * sina * comm;
		stiffmatrix[1][2] = -stiffmatrix[0][1];
		stiffmatrix[1][3] = -stiffmatrix[1][1];
		stiffmatrix[2][0] =  stiffmatrix[0][2];
		stiffmatrix[2][1] =  stiffmatrix[1][2];
		stiffmatrix[2][2] =  stiffmatrix[0][0];
		stiffmatrix[2][3] =  stiffmatrix[0][1];
		stiffmatrix[3][0] =  stiffmatrix[0][3];
		stiffmatrix[3][1] =  stiffmatrix[1][3];
		stiffmatrix[3][2] =  stiffmatrix[2][3];
		stiffmatrix[3][3] =  stiffmatrix[1][1];

		/* Assembly of Element Stiffness to Overall
		   Stiffness   */
		lm[1] = 2 * member[m][0] - 1;
		lm[0] = lm[1] - 1;
		lm[3] = 2 * member[m][1] - 1;
		lm[2] = lm[3] - 1;
		for (k = 0; k < 4; k++)
		{
			for(g=0; g<4; g++)
			{
				i = lm [k];
				j = lm [g];
				a[i][j] = a[i][j] + stiffmatrix[k][g];
			}
		}

	}
/*end stiff*/

/*displ*/
	for(n = 0; n < num_node; n++)
	{
		nn = (2 * n) + 1;
		b[nn] = fY[n];
		b[nn-1] = fX[n];
	}

	for(n = 0; n < num_node; n++)
	{
		nn = (2 * n) + 1;
		switch(kode[n])
		{
			case 1:
				for(i=0;i<neq;i++)
				{
					b[i]=b[i]-(a[i][nn-1]*fX[n]);
					a[i][nn-1]=0.;
					a[nn-1][i]=0.;
				}
				a[nn-1][nn-1]=1.00;
				b[nn-1]=fX[n];
				break;
			case 2:
				for(i=0;i<neq;i++)
				{
					b[i]=b[i]-(a[i][nn]*fY[n]);
					a[i][nn]=0.;
					a[nn][i]=0.;
				}
				a[nn][nn]=1.00;
				b[nn]=fY[n];
				break;
			case 3:
				for(i=0;i<neq;i++)
				{
					b[i]=b[i]-a[i][nn-1]*fX[n];
					a[i][nn-1]=0.;
					a[nn-1][i]=0.;
				}
				a[nn-1][nn-1]=1.00;
				b[nn-1]=fX[n];
				for(i=0;i<neq;i++)
				{
					b[i]=b[i]-a[i][nn]*fY[n];
					a[i][nn]=0.;
					a[nn][i]=0.;
				}
				a[nn][nn]=1.00;
				b[nn]=fY[n];
				break;
		}
	};


/*	symsol();*/
	neq1 = nl = 0;
	neq1=neq-1;
	nl = neq - 1;

	for (n=0;n<nl;n++)
	{
		if (a[n][n]<=0.)
		{
			printf("zero or negative main-diagonal %d\n",n);
			exit(2);
		}

		n1 = n +1 ;
		for (j=n1;j<neq;j++) {
			/*printf("l279: %f %f\n",a[n][j],a[n][n]);*/
			a[n][j]=a[n][j]/a[n][n];
		}
		for(i=n1;i<neq;i++)
		{
			if(a[n][i]==0.0)
				b[i]=b[i]-a[n][i]*b[n];
			else
			{
				for(j=i;j<neq;j++)
				{
					/*printf("l288: %f %f %f\n",a[i][j],a[i][n],a[n][j]);*/
					a[i][j]=a[i][j]-a[i][n]*a[n][j];
					a[j][i]=a[i][j];
				}
			}
			b[i]=b[i]-a[n][i]*b[n];
		}
		b[n]=b[n]/a[n][n];
	};

	m=neq1;
	b[m]=b[m]/a[m][m];
	for(n=0;n<nl;n++)
	{
		m1=m;
		m=m-1;
		for(j=m1;j<neq;j++) {
			b[m]=b[m]-b[j]*a[m][j];
			/*printf("l306: %f %f %f\n",b[m],b[j],a[m][j]);*/
		}
	}
/*end symbol*/

	for(i=0;i<num_node;i++)
	{
		displX[i] = b[2*i];
		displY[i] = b[2*i+1];
	}

	/*printf("\n NODE    VEL-X     VEL-Y\n");
	for (m=0;m<num_node;m++)
		printf("%3d  %10.7f %10.7f\n", (m+1),displX[m],displY[m]);*/
/*end displ*/


/*stress*/
	/*printf("\nELEM I-NOD J-NOD    MEM-FOR     STRESS\n");*/
	for(m=0;m<num_elem;m++)
	{
		i   =  member[m][0] - 1;
		j   =  member[m][1] - 1;
		dx  =  coordx[j]-coordx[i];
		dy  =  coordy[j]-coordy[i];
		xl  =  sqrt((dx*dx) + (dy * dy));
		du  =  displX[j] -displX[i];
		dv  =  displY[j] -displY[i];
		dl  =  dv*dy/xl+du*dx/xl;
		stress[m] =  dl*ymod[m]/xl;
		p   =  area[m]*stress[m];

		/*printf("%d    %d    %d  %10.4f  %10.4f\n", (m+1), (i+1), (j+1), p,strss[m]);*/
	}
/*end stress*/

	for (i = 0; i < 10; i++) {
		/*weight of the truss*/
		somF += ro*area[i]*length[i];
	}
	for(i = 0; i < 10; i++) {
		float aux = fabs(stress[i]);
		somg1 += (aux > 25) ? aux - 25. : 0;          /*stress constraint (25 ksi)*/
//		if (print_fitness == 1)
//			printf("G1[%d] = %f\n", i, aux);
	}
	for(i = 0; i < 6; i++) {
		/*displacement constraint (2 in.)*/
		float aux = fabs(displX[i])/* * 1000.*/;
		somg2 += (aux > 2) ? aux - 2. : 0;
//		if (print_fitness == 1)
//			printf("G2X[%d] = %f\n", i, aux);
	}
	for(i = 0; i < 6; i++) {
		/*displacement constraint (2 in.)*/
		float aux = fabs(displY[i])/* * 1000.*/;
		somg2 += (aux > 2) ? aux - 2. : 0;
//		if (print_fitness == 1)
//			printf("G2Y[%d] = %f\n", i, aux);
	}
//	if (print_fitness == 1) {
//		printf("F %f G1 %f G2 %f\n", somF, somg1, somg2);
//	}

        func_eval[sp_index] += 1;

	return (somF+P*(somg1+somg2));                /*Transformed Function*/
/////END 10-bar-truss

    case 23:
//2D-AB
    amino_pos = (amino *) malloc ((ProteinSize) * sizeof (amino));

    amino_pos[0].x = amino_pos[0].y = 0;
    amino_pos[1].x = 1;
    amino_pos[1].y = 0;
    for (i = 1; (i < (ProteinSize-1)); i++)
    {   
	if(sol[i]>PI)
		sol[i] = randon(0.0,PI);
	if(sol[i]<-PI)
		sol[i] = randon(-PI,0.0);

	a_ab = amino_pos[i].x-amino_pos[i-1].x;
        b_ab = amino_pos[i].y-amino_pos[i-1].y;
        amino_pos[i+1].x = amino_pos[i].x+a_ab*cos(sol[i-1])-b_ab*sin(sol[i-1]);
        amino_pos[i+1].y = amino_pos[i].y+b_ab*cos(sol[i-1])+a_ab*sin(sol[i-1]);
    }

    v1 = 0;
    for (i = 1; (i < (ProteinSize-1) ); i++) 
	v1 += (1.0-cos(sol[i-1]))/4.0;

    v2 = 0;
    for (i = 0; (i < (ProteinSize-2)); i++)
    {   for (j = (i+2); (j < ProteinSize); j++)
        {
	   //c_ab = (1.0+Aminoacido[i].Tipo+Aminoacido[j].Tipo+5.0*Aminoacido[i].Tipo*Aminoacido[j].Tipo)/8.0; //Energy for the Bonds AA - BB - AB
	    if (Aminoacido[i].Tipo == 1 && Aminoacido[j].Tipo == 1) //AA bond
		c_ab = 1;
	    else if (Aminoacido[i].Tipo == -1 && Aminoacido[j].Tipo == -1) //BB bond
		c_ab = 0.5;
	    else
		c_ab = -0.5; //AB or BA bond

            d_ab = sqrt(((amino_pos[i].x-amino_pos[j].x)*(amino_pos[i].x-amino_pos[j].x))+((amino_pos[i].y-amino_pos[j].y)*(amino_pos[i].y-amino_pos[j].y))); //Distance for Lennard-Jones        
	    v2 += 4.0*(1/pow(d_ab,12)-c_ab/pow(d_ab,6));
        }
    }

       func_eval[sp_index] += 1;

	free(amino_pos);

       return(v1 + v2);

	case 24: //Schaffer F6
	top=0;
        top1=0;
        top2=0;

	for(j=0;j<DIM-1;j++)
        {
	    top1 = pow( sin( sqrt( pow(sol[j],(double)2) + pow(sol[j+1],(double)2) ) ), (double)2) - 0.5;
	    top2 = pow( 1.0 + 0.001* ( pow(sol[j],(double)2) + pow(sol[j+1],(double)2) ), (double)2);
	    top = top + (0.5 + top1/top2);
        }
	func_eval[sp_index] += 1;
        return top;
        
	case 25: //Shifted_Sphere
	// x* = o , F(x*) = f_bias1 = - 450
	    for(i=0;i<DIM;i++){
	        z = sol[i] - sphere[i];
	        Fx += z*z;
	    }
	func_eval[sp_index] += 1;
	return Fx + f_bias[0];

	case 26: //Shifted Schwefel Problem 2.21
	// x* = o , F(x*) = f_bias1 = - 450
		Fx = abss(sol[0]);
		for(i=1;i<DIM;i++){
	    	  z = sol[i] - schwefel[i];
	          Fx = max(Fx , abss(z));
    		}
	func_eval[sp_index] += 1;
	return Fx + f_bias[1]; 

	case 27: //Shifted_Rosenbrock
	// x* = o , F(x* ) = f_bias3 = 390

	    for(i=0;i<DIM;i++) zx[i] = sol[i] - rosenbrock[i] + 1;   

	    for(i=0;i<DIM-1;i++){    
	        Fx = Fx + 100*( pow((pow(zx[i],2)-zx[i+1]) , 2) ) + pow((zx[i]-1) , 2);
	    }
	    func_eval[sp_index] += 1;
	return Fx + f_bias[2]; 

	case 28: //Shifted Rastrigin
	//x* = o , F( x * ) = f_bias4 = - 330
	    for(i=0;i<DIM;i++){  
	        z = sol[i] - rastrigin[i];
	        Fx = Fx + ( pow(z,2) - 10*cos(2*pi*z) + 10);
	    }
	    func_eval[sp_index] += 1;
	    return Fx + f_bias[3];

	case 29: //Shifted Griewank
	//x* = o , F(x* ) = f_bias5 = -180
	    top1 = 0;
	    top2 = 1;
	    for(i=0;i<DIM;i++){       
	        z = sol[i] - griewank[i];
	        top1 = top1 + ( pow(z,2) / 4000 );
	        top2 = top2 * ( cos(z/sqrt(i+1)));
	
	    }
	    func_eval[sp_index] += 1;
	    return (top1 - top2 + 1 + f_bias[4]);

	case 30://Shifted Ackley
	// x* = o , F(x* ) = f_bias6 = - 140
	    top1 = 0;
	    top2 = 0;
	    Fx = 0;
	    for(i=0;i<DIM;i++){   
	        z = sol[i] - ackley[i];
	        top1 = top1 + pow(z , 2 );
	        top2 = top2 + cos(2*pi*z);
	    }
	    Fx = -20*exp(-0.2*sqrt(top1/DIM)) -exp(top2/DIM) + 20 + e + f_bias[5];
	    func_eval[sp_index] += 1;
	    return Fx;
        

    default:
        printf("Info: Invalid function.\n") ;
        exit(0);
    }
}

double sqr( double x )
{
  return x*x;
};

double TempValue(double x,int a,int k,int m)
{
	double temp = 0.0;
	if( x > a)
	{
		temp = k*pow(x-a,m);
	}
	else if( x <= a && x >= -a)
	{
		temp = 0.0;
	}
	else
	{
		temp = k*pow(-x-a,m);
	}
	return temp;
}

void prepararObjFunc(void)//domain for each variable of each problem
{
    switch (FUNCTION) {
    case 0: //Rastrigin
        lb = -5.12;
        ub = 5.12;
        break;
    case 1: //Schaffer F7
        lb = -100.00;
        ub = 100.00;
        break;
    case 2: //Griewank
        lb = -600.00;
        ub = 600.00;
        break;
    case 3: //Ackley
        lb = -32.00;
        ub = 32.00;
        break;
    case 4: //Rosenbrock
        lb = -30.00;
        ub = 30.00;
        break;
    case 5: //Sphere
	lb = -100.00;
	ub = 100.00;
	break;
    case 6: //StretchedV
	lb = -10.00;
	ub = 10.00;
	break;
    case 7: // Schwefel's function 2.22
	lb = -10.00;
	ub = 10.00;
	break;
    case 8: // Step function
	lb = -100.00;
	ub = 100.00;
	break;
    case 9: // Generalized Schwefel's function 2.26
	lb = -500.00;
	ub = 500.00;
	break;
    case 10: // Generalized penalized function #1
	lb = -50.00;
	ub = 50.00;
	break;
    case 11: // Generalized penalized function #2
	lb = -50.00;
	ub = 50.00;
	break;
    case 12: // Levy Function
	lb = -10.00;
	ub = 10.00;
	break;
    case 13: // Zakharov
	lb = -5.00;
	ub = 10.00;
	break;
    case 14: // EggHolder
	lb = -512.00;
	ub = 512.00;
	break;
    case 15: // Holzman
	lb = -10.00;
	ub = 10.00;
	break;
    case 16: // Michalewitz
	lb = 0.00;
	ub = PI;
	break;
    case 17: // Multimod
	lb = -10.00;
	ub = 10.00;
	break;
    case 18: // Powell
	lb = -4.00;
	ub = 5.00;
	break;
    case 19: // Rana
	lb = -512.00;
	ub = 512.00;
	break;
    case 20: // Shubert
	lb = -10.00;
	ub = 10.00;
	break;
    case 21: //3D-AB
	lb = -1*PI;
	ub = PI;
	break;
    case 22: //10-bar-truss
	lb = 0.10;
	ub = 35.00;
	break;
    case 23: //2D-AB
	lb = -1*PI;
	ub = PI;
	break;
    case 24: //Schaffer F6
        lb = -100.00;
        ub = 100.00;
        break;
    case 25: //Shifted Sphere
        lb = -100.00;
        ub = 100.00;
        break;
    case 26: //Shifted Schwefel's Problem 2.21
        lb = -100.00;
        ub = 100.00;
        break;
    case 27: //Shifted Rosenbrock
        lb = -100.00;
        ub = 100.00;
        break;
    case 28: //Shifted Rastrigin
        lb = -5.12;
        ub = 5.12;
        break;
    case 29: //Griewank
        lb = -600.00;
        ub = 600.00;
        break;
    case 30: //Ackley
        lb = -32.00;
        ub = 32.00;
        break;
    default:
        printf("Info: Invalid function.\n") ;
        exit(0);
    }
}

//3D-AB
void Rotation(double x,double y,double z,double a,double b,double c,double u,double v,double w,double t,double *rx,double *ry,double *rz)
{
   double u2 = u*u;
   double v2 = v*v;
   double w2 = w*w;
   double cost = cos(t);
   double sint = sin(t);
   double l2 = u2+v2+w2;
   double l = sqrt(l2);

   *rx = (a*(v2+w2)+u*(-b*v-c*w+u*x+v*y+w*z)+(-a*(v2+w2)+u*(b*v+c*w-v*y-w*z)+(v2+w2)*x)*cost+l*(-c*v+b*w-w*y+v*z)*sint)/l2;
   *ry = (b*(u2+w2)+v*(-a*u-c*w+u*x+v*y+w*z)+(-b*(u2+w2)+v*(a*u+c*w-u*x-w*z)+(u2+w2)*y)*cost+l*(c*u-a*w+w*x-u*z)*sint)/l2;
   *rz = (c*(u2+v2)+w*(-a*u-b*v+u*x+v*y+w*z)+(-c*(u2+v2)+w*(a*u+b*v-u*x-v*y)+(u2+v2)*z)*cost+l*(-b*u+a*v-v*x+u*y)*sint)/l2;
}

