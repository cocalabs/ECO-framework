# ECO Framework

![ECO Framework](https://github.com/cocalabs/ECO-framework/blob/master/eco-view.jpg "ECO Framework View")


	Description: The ANSI C code of the ECO approach
	Programmer:  Rafael Stubs Parpinelli
	E-mail:      rsparpin@gmail.com
	Date:	     20/11/2012
	Lisence:     Free
	Note:        The system was developed using Linux. If you use this code, pls cite one of the following works:
	             1. PARPINELLI, R.S.; LOPES, H.S. "Biological Plausibility in Optimization: An Ecosystemic View". 
			International Journal of Bio-Inspired Computation (Online), v. 4, p. 345-358, 2012.
		     2. PARPINELLI, R.S.; LOPES, H. S. "A computational ecosystem for optimization: review and 
		     	perspectives for future research". Memetic Computing, v. 7, p. 29-41, 2014.
	To compile:  1. Go to ECO directory
		     2. Type: make

	To debug using gdb: gcc mersenne.o -o eco-thread eco-thread.c -lm -lpthread -g -O0 -v -da -Q

	To run: ./eco-cluster input_eco.in

	Files:	     1. 'eco.h': headers and global variables
		     2. 'eco-cluster.c': the main system itself
		     3. 'input_eco.in': set up the prarameters
		     4. 'ProblemFunc.c': contains the objective functions

	Directories: 1. 'algorithms': contains the meta-heuristics that can be used
		     2. 'report': contains the output of the system
		     3. 'src-clstr': contains the source code of the clustering algorithm
