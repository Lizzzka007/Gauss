#include "functions.h"
#include <iostream> 
#include <math.h>
#include <sys/time.h>

using namespace std;



void vector(int n, double * a, double * b)  // make vector b
{
	for(int i = 0; i < n; i++)
	{
		b[i] = 0;

		for(int j = 0; j < (n+1)/2  ; j++)
		{
			b[i] += a[i*n + 2*j];
		}
	}
}


int print_matrix(int l, int n, int r, double * a) // print matrix l * n
{
	if ( r > n )
	{
		r = n;
	}

	if ( r > l )
	{
		r = l;
	}

	for(int i = 0; i < r; i ++)
	{
		for(int j = 0; j < r; j ++) printf("%10.3e ", a[i*n + j]);
		cout<<endl;
	}
	cout<<endl;

	return 0;
}



void print_discrepancy(int s5, int n, int m, double * a, double * b, double * x, double *v_1, double *v_2, float elapsed, const char * name) // norma nevyazki
{
	int k, l, last, min, s, q;
	double sum = 0.0, norm_b = 0.0, margin_error = 0.0, sum1 = 0.0;
	int last1, s1;
	int  a_v, a_h, t = 0;
	double *p_a;

	struct timeval tv1;
	struct timeval tv2;
 	struct timezone tz;
  	gettimeofday(&tv1, &tz);

	k = n/m;
	l = n - m*k;

	if( l == 0) 
	{
		last = k;
		min = m;
	}
	else 
	{
		last = k + 1;
		min = l;
	}

	for(int i = 0; i < last; i ++)
	{
		//return_i_vector_part(i, m, k, l, x, v_2);// put X_i into V_2

		if( m % 2 == 0)
		{
			for( q = 0; q < m; q ++)
			{
				if( q % 2 == 0) v_2[q] = 1.0;
				else v_2[q] = 0.0;
			}
		}
		else
		{
			if( i % 2 == 0)
			{
				for(q = 0; q < m; q ++)
				{
					if( q % 2 == 0) v_2[q] = 1.0;
					else v_2[q] = 0.0;
				}
			}
			else
			{
				for(q = 0; q < m; q ++)
				{
					if( q % 2 == 0) v_2[q] = 0.0;
					else v_2[q] = 1.0;
				}
			}
		}


		if( i == last - 1)	for(s = 0; s < min; s++)	margin_error += (x[i*m + s] - v_2[s])*(x[i*m + s] - v_2[s]);
		else	for(s = 0; s < m; s++)	margin_error += (x[i*m + s] - v_2[s])*(x[i*m + s] - v_2[s]);

		//return_i_vector_part(0, m, k, l, x, v_2);// put X_i into V_2
		//return_i_vector_part(i, m, k, l, b, v_3);// put B_i into V_3

		if( l != 0)
		{
			if( i == k) last1 = l;
			else last1 = m;
		}
		else last1 = m;

		for(s1 = 0; s1 < last1; s1++) norm_b += b[i*m + s1] * b[i*m + s1];
			
		for(int j = 0; j < last; j ++)
		{
			

			//return_i_vector_part(j, m, k, l, x, v_2);// put X_j into V_2
			//return_i_j_matrix_block(i, j, n, m, k, l, a, v_1);// put A_ij into V_1

			t = 0;
			if( i == k) a_v = l;
			else a_v = m;
			if( j == k) a_h = l;
			else a_h = m;

			p_a = a + (i*n+j)*m;

			for( s1 = 0; s1 < a_v; s1++)
			{
				for(q = 0; q < a_h; q++) 
				{
					v_1[t] = p_a[ s1*n + q ];
					p_a[ s1*n + q ] = 0.0;
					t++;
				}
			}
			t = 0;

			//null_matrix(min, m, v_3);

			//put_back_i_j_matrix_block(i, j, n, m, k, l, a, v_3);

			if( i == last - 1)
			{
				if( j == last - 1)//	matrix_multiplication(min, min, 1, v_1, v_2, v_3);
				{
					for(s1 = 0; s1 < min; s1 ++)
					{
						for(t = 0; t < 1; t ++ )
						{
							sum1 = 0.0;

							for(q = 0; q < min; q ++ ) sum1 += v_1[s1*min + q] * x[j*m + q*1 + t];

							v_2[s1*1 + t] = sum1;	
						}
					}
				}
				else//	matrix_multiplication(min, m, 1, v_1, v_2, v_3);
				{
					for(s1 = 0; s1 < min; s1 ++)
					{
						for(t = 0; t < 1; t ++ )
						{
							sum1 = 0.0;

							for(q = 0; q < m; q ++ ) sum1 += v_1[s1*m + q] * x[j*m + q*1 + t];

							v_2[s1*1 + t] = sum1;	
						}
					}
				}

				for( s = 0; s < min; s ++){  a[(i*n)*m + s] += v_2[s];}

			}
			else
			{
				if( j == last - 1)//	matrix_multiplication(m, min, 1, v_1, v_2, v_3);
				{
					for(s1 = 0; s1 < m; s1 ++)
					{
						for(t = 0; t < 1; t ++ )
						{
							sum1 = 0.0;

							for(q = 0; q < min; q ++ ) sum1 += v_1[s1*min + q] * x[j*m + q*1 + t];

							v_2[s1*1 + t] = sum1;	
						}
					}
				}

				else//	matrix_multiplication(m, m, 1, v_1, v_2, v_3);
				{
					for(s1 = 0; s1 < m; s1 ++)
					{
						for(t = 0; t < 1; t ++ )
						{
							sum1 = 0.0;

							for(q = 0; q < m; q ++ ) sum1 += v_1[s1*m + q] * x[j*m + q*1 + t];

							v_2[s1*1 + t] = sum1;	
						}
					}
				}

				for( s = 0; s < m; s ++) { a[(i*n)*m + s] += v_2[s];}

			}
		}
		
		//return_i_vector_part(i, m, k, l, b, v_1);// put X_i into V_1
					
		if( i == last - 1) for(s = 0; s < min; s++)	sum += (a[(i*n)*m + s] - b[i*m + s])*(a[(i*n)*m + s] - b[i*m + s]);
		else	for(s = 0; s < m; s++)	sum += (a[(i*n )*m + s] - b[i*m + s])*(a[(i*n)*m + s] - b[i*m + s]);
	}

	sum = sqrt(sum) / sqrt(norm_b);

	gettimeofday(&tv2, &tz);

	//printf(" Residual norm is %10.3Le\n", sum);
	//std::cout<<" Time of searching of the residual norm and margin error = "<<tv2.tv_sec+tv2.tv_usec/1000000.0 - tv1.tv_sec-tv1.tv_usec/1000000.0<<" sec"<<endl;
	//printf(" Margin error is %10.3Le\n", sqrt(margin_error));

	std::cout<<" Time of searching of the residual norm and margin error = "<<tv2.tv_sec+tv2.tv_usec/1000000.0 - tv1.tv_sec-tv1.tv_usec/1000000.0<<" sec"<<endl;
	printf(" Margin error is %10.3e\n", sqrt(margin_error));
	printf ("%s : residual = %e elapsed = %.2f for s = %d n = %d m = %d\n", name, sum, elapsed, s5, n, m);
	//printf (" Residual = %e elapsed = %.2f for s = %d n = %d m = %d\n", sum, tv_seconds, s5, n, m);
}














