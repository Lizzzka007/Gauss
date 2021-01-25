#include <iostream>
#include "functions.h"
#include <sched.h>
#include <sys/sysinfo.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

using namespace std;




int main(int argc, char* argv[])
{
	cpu_set_t cpu;
	int nprocs = get_nprocs();
	CPU_ZERO(&cpu);
	CPU_SET(nprocs - 1, &cpu);
	sched_setaffinity(getpid(), sizeof(cpu), &cpu);
	
	if( argc < 5 || argc > 6) 
	{
		std::cout<<"argc is less or greater then needed"<<endl;

		return 1;
	}

	struct timeval tv1;
	struct timeval tv2;
 	struct timezone tz;

	int r, n, m, l, k, s; 

	if( !( sscanf(argv[1], "%d", &n) == 1 && n > 0 
		&& sscanf(argv[2], "%d", &m) == 1 && m > 0  
		&& sscanf(argv[3], "%d", &r) == 1 && r >= 0 
		&& sscanf(argv[4], "%d", &s) == 1 && s >= 0 && s < 5) )
	{

		if( sscanf(argv[1], "%d", &n) == 0 )
		{
			printf("n should be an integer\n");
			
			return 1;
		}

		if( sscanf(argv[2], "%d", &m) == 0 )
		{
			printf("m should be an integer\n");
			
			return 1;
		}

		if( sscanf(argv[3], "%d", &r) == 0 )
		{
			printf("r should be an integer\n");
			
			return 1;
		}

		if( sscanf(argv[4], "%d", &s) == 0 )
		{
			printf("s should be an integer\n");
			
			return 1;
		}

		if( n == 0) 
		{
			printf("As n = 0 there are no equations to solve\n");
			
			return 1;
		}

		if( m == 0) printf("m, which equals 0, is not an acceptable parameter for the algorithm\n");
		return 1;
	}

	if( (atoi(argv[2]) > atoi(argv[1])) )
	{
		if( atoi(argv[2]) > atoi(argv[1]) ) std::cout<<"m is greater then n, m = "<<m<<", n = "<<n<<endl;

		return 1;
	}

	if( atoi(argv[3]) > atoi(argv[1]) ) 
	{
		if(n <= 8) r = n;
		else r = 8;
	}

	double * b, *a, *x, *v_1, *v_2;

	k = n/m;
	l = n - m*k;
	a = (double*) malloc(n*n * sizeof(double));
	b = (double*) malloc(n * sizeof(double));
	x = (double*) malloc(n * sizeof(double));
	v_1 = (double*) malloc(m*m * sizeof(double));
	v_2 = (double*) malloc(m*m * sizeof(double));
	

	if(argc == 5) //enter_matrix(n, m, k, s, l, a);   //matrix initialization
	{
		if( s == 0)
		{
			std::cout<<" Wrong s: if you want to input matrix using the function, you should enter s such that 0 < s < 5 "<<endl;
			printf(" Usage: %s n m r s [name].txt\n", argv[0]);

			free(a);
			free(b);
			free(x);
			free(v_1);
			free(v_2);

			return 1;
		}

		enter_matrix(n, m, k, s, l, a); 
	}

	if(argc == 6) 
	{
		const char* filename;
		if( s != 0)
		{
			std::cout<<" Wrong s: if you want to input matrix from the txt file, you should enter s as 0 "<<endl;
			printf(" Usage: %s n m r s [name].txt\n", argv[0]);

			free(a);
			free(b);
			free(x);
			free(v_1);
			free(v_2);

			return 1;
		}

		filename = argv[5];

		if(enter_matrix_from_file(a, n, filename) == 1)
		{
			free(a);
			free(b);
			free(x);
			free(v_1);
			free(v_2);
			
			return 1;
		}
	}

	printf("\n\n");

	if ( (print_matrix(n, n, r, a)) == 1) 
	{
		free(a);
		free(b);
		free(x);
		free(v_1);
		free(v_2);

		return 1;
	}

	vector(n, a, b); // enter vector b

  	gettimeofday(&tv1, &tz);

	if(solution(n, m, a, b, x, v_1, v_2) == 1)
	{
		std::cout<<"Algorithm is not acceptable"<<endl;
		free(a);
		free(b);
		free(x);
		free(v_1);
		free(v_2);

		return 1;
	}

	gettimeofday(&tv2, &tz);

	if(argc == 5) enter_matrix(n, m, k, s, l, a);
	if(argc == 6) 
	{
		const char* filename;

		filename = argv[5];

		enter_matrix_from_file(a, n, filename);
	}

	vector(n, a, b); // enter vector b

	std::cout<<" Vector b:"<<endl;                    // vector b
	for(int i = 0; i < r; i++) printf(" %10.3e\n", b[i]);
	printf("\n");

	std::cout<<" Vector x:"<<endl;                    // vector x
	for(int i = 0; i < r; i++) printf(" %10.3e\n", x[i]);
	printf("\n");    

	print_discrepancy(s, n, m, a, b, x, v_1, v_2, (float)(tv2.tv_sec+tv2.tv_usec/1000000.0 - tv1.tv_sec-tv1.tv_usec/1000000.0), argv[0]); // norma nevyazki

	if(argc == 5) enter_matrix(n, m, k, s, l, a);
	if(argc == 6) 
	{
		const char* filename;

		filename = argv[5];

		enter_matrix_from_file(a, n, filename);
	}

	free(a);
	free(b);
	free(x);
	free(v_1);
	free(v_2);

	return 0;
}


